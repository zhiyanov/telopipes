# Environment findings (Phase 0)

Recorded 2026-08-29 on the development workstation. These measurements decide the container
architecture for the whole project, so they are kept under version control rather than in a
chat log.

## Host

| | |
|---|---|
| OS | Fedora 44, Asahi Linux |
| Kernel | `7.1.6-400.asahi.fc44.aarch64+16k` |
| Arch | `aarch64` (Apple Silicon) |
| **Page size** | **16384 bytes (16 KB)** |
| CPU / RAM / free disk | 8 cores / 15 GB / ~85 GB |
| Container engine | podman 5.8.4, rootless, `overlay` driver |
| SELinux | **Enforcing** |
| Filesystem | btrfs (`/home` on `/dev/nvme0n1p6`) |
| podman graphroot | `~/.local/share/containers/storage` |
| `podman.socket` | present but disabled — we shell out to the CLI, so it stays disabled |

## Finding 1 — x86_64 emulation is impossible on this machine

The original plan was to run amd64 images under FEX-Emu, since `binfmt_misc` has an
`FEX-x86_64` entry registered with flags `POCF` — and the `F` (fix-binary) flag normally means
the interpreter is reachable from inside a container mount namespace. That reasoning was
correct but irrelevant: FEX cannot run at all here.

```
$ podman run --rm --platform linux/amd64 docker.io/library/debian:bookworm-slim uname -m
exec container process (missing dynamic library?) `/usr/bin/uname`: No such file or directory

$ /usr/bin/FEXInterpreter -R ./amd64rootfs ./amd64rootfs/bin/uname -m
<jemalloc>: Unsupported system page size
terminate called without an active exception
(core dumped)
```

**Root cause: the 16 KB page size.** x86_64 architecturally guarantees 4 KB pages, and mmap
granularity cannot be faithfully emulated on a host with larger pages. This is a property of
the running kernel, not a misconfiguration — no amount of FEX RootFS setup fixes it.

`qemu-user-static` is packaged for Fedora 44 (`qemu-user-static-x86`) but hits the same
limitation: user-mode QEMU also cannot emulate a smaller-page guest on a larger-page host.

Asahi Linux does ship a 4 KB-page kernel variant. Booting it would re-enable emulation, at the
cost of a system-wide performance regression on Apple Silicon and a reboot. **Rejected** in
favour of native builds, which are both faster and less invasive.

## Finding 2 — native arm64 works, and every dependency is available

```
$ podman run --rm --platform linux/arm64 docker.io/library/debian:bookworm-slim \
    sh -c 'uname -m; getconf PAGESIZE'
aarch64
16384
```

| Dependency | arm64 availability |
|---|---|
| `minimap2` | Debian bookworm `2.24+dfsg-3+b1`, Ubuntu noble newer |
| `samtools` | Debian bookworm `1.16.1-1` |
| `rocker/r-ver` | multi-arch (arm64 in manifest) |
| `python:3.12-slim` | multi-arch |
| `bioconductor/bioconductor_docker` | **amd64 only — unusable here** |

## Consequences for the design

1. **Both images are built natively for the host architecture.** The `--platform linux/amd64`
   decision is void.
2. **Pipeline A changes base image.** `bioconductor/bioconductor_docker:RELEASE_3_18` is amd64
   only, so it moves to `rocker/r-ver:4.4.1` with `BiocManager` installing the same three
   packages (`Biostrings`, `GenomicAlignments`, `S4Vectors`). This shifts R 4.3 → 4.4 and
   Bioconductor 3.18 → 3.19. Those three are core, long-stable packages and the pipeline uses
   only their basic APIs, but the change is a real deviation from the reference and is why the
   smoke test against the published numbers matters.
3. **On arm64, Bioconductor compiles from source.** rocker points R at Posit Package Manager,
   which serves binaries for amd64 Linux only; the Dockerfile overrides `repos` to CRAN so the
   behaviour is explicit rather than a silent fallback. This is a one-time build cost.
4. **Pipeline B is unaffected** and is the bigger win: it is a tight pure-Python loop over every
   base of every read, which is precisely where emulation would have hurt most. It now runs at
   native speed.

## Risks this introduces, to be closed by the smoke test

- **Tool versions differ** from those that produced `reference/results/`. Debian/Ubuntu ship
  minimap2 2.24+ and samtools 1.16+; the reference image installed whatever its Debian base had
  at the time. Different aligner versions can shift individual alignments.
- **Floating-point scoring.** The TeloTag barcode step uses `Biostrings::pairwiseAlignment`
  with a `PhredQuality`-derived substitution matrix, producing fractional scores compared
  against a threshold of 20. A read sitting exactly on that boundary could in principle be
  classified differently on a different architecture or library version.

Both are detected by reproducing the reference figures: pipeline A `n=4362, median=4589.5`,
pipeline B `n=5049, median=4378.0`.

## Finding 3 — the vendored TeloBP matches upstream exactly

The reference project vendors GreiderLab TeloBP as a local editable install. Installing it from
GitHub at a pinned commit is only safe if that commit is what actually produced
`reference/results/`. It is:

```
$ git clone https://github.com/GreiderLab/TeloBP.git && git checkout 05bd9d8
05bd9d8bb27b9d52e69e9c1c4b6607744365ea9f  2025-12-26  "Fixing dependancy issue"

$ diff -rq TeloBP/ reference/pipeline/TeloBP/TeloBP/     # → no differences
$ diff -rq Scripts/ reference/pipeline/TeloBP/Scripts/   # → no differences
```

So `pipelines/telonp/requirements.txt` pins
`TeloBP @ git+https://github.com/GreiderLab/TeloBP@05bd9d8bb27b9d52e69e9c1c4b6607744365ea9f`
rather than vendoring a second copy into this repo.

## Images as built

| Image | Base | Build time | Size |
|---|---|---|---|
| `telomers/pipeline-a:0.1.0` | `rocker/r-ver:4.4.1` | 5 min 06 s | 1.32 GB |
| `telomers/pipeline-b:0.1.0` | `python:3.12-slim` | 30 s | 634 MB |

Resulting toolchain, for comparison against whatever produced `reference/results/`:

```
minimap2 2.24-r1122          samtools 1.13
R 4.4.1                      Bioconductor 3.19
Biostrings 2.72.1            GenomicAlignments 1.40.0     S4Vectors 0.42.1
RColorBrewer 1.1.3           optparse 1.8.2               Snakemake 9.26.1
```

### Three build problems worth recording

1. **`short-name-mode = "enforcing"`.** `/etc/containers/registries.conf` refuses to resolve
   unqualified image names without a TTY prompt, so every `FROM` must be fully qualified
   (`docker.io/rocker/r-ver:4.4.1`, not `rocker/r-ver:4.4.1`).
2. **`rocker/r-ver` ships runtime libraries but no headers.** With no arm64 Bioconductor
   binaries everything compiles from source, so the build died on `zlib.h`, `curl/curl.h` and
   the OpenSSL headers, cascading into ~10 misleading "dependency not available" errors.
   `littler` was collateral damage: `BiocManager::install(version=)` tries to refresh rocker's
   own packages, so it now passes `update=FALSE` as well.
3. **`rocker/r-ver` ships no Python at all**, and the base's Ubuntu release decides what apt
   would give — `r-ver:4.4.x` is Ubuntu 22.04, i.e. Python 3.10, below Snakemake 9's floor.
   Pinning Snakemake 7 for this image would have split the two pipelines across Snakemake
   majors and dragged in `datrie`, which has its own build problems. Instead the image gets a
   standalone Python 3.12 via `uv`, so the orchestrator is independent of the base image and
   both pipelines stay on Snakemake 9.26.1. This is also why R stays at 4.4 / Bioc 3.19 rather
   than jumping to a newer rocker tag purely to chase a newer system Python.
