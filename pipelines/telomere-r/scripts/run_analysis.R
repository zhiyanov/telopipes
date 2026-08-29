## Method 1 driver: mapping-based telomere length.
## Sources the GreiderLab pipeline functions and runs them with the parameters for one sample.
## Runs inside the telomers/pipeline-a image.
##
## Differences from the reference driver (see docs/PATCHES.md):
##   * named arguments instead of positional ones, so the Snakefile and the web form can pass
##     parameters that were previously hardcoded here;
##   * --barcode-seq accepts a literal sequence. Upstream resolved a *name* with get(), which
##     works only for the 18 TeloTags defined in telobp_functions.R and fails with a bare
##     "object 'NB99uq' not found" for anything else;
##   * stable output filenames, so Snakemake can declare them without globbing;
##   * --emit-all writes the per-read table for every primary read, including those the filters
##     reject, with the reason. The filter flags already exist inside writeTelTableFn; upstream
##     just never wrote them out, leaving read accounting to be scraped from stdout.

suppressPackageStartupMessages(library(optparse))

opts <- list(
    make_option("--bam",           type = "character", help = "coordinate-sorted, indexed BAM"),
    make_option("--outdir",        type = "character", help = "output directory"),
    make_option("--exp-name",      type = "character", dest = "exp_name",
                help = "experiment/sample label used in plot titles"),
    make_option("--chr-arm-ln",    type = "integer", default = 500000L, dest = "chr_arm_ln",
                help = "per-record length of the cut reference [%default]"),
    make_option("--barcode-name",  type = "character", default = "NB65uq", dest = "barcode_name",
                help = "TeloTag name, used for column labels [%default]"),
    make_option("--barcode-seq",   type = "character", default = NULL, dest = "barcode_seq",
                help = "TeloTag sequence; if omitted, --barcode-name is looked up"),
    make_option("--read-len-min",  type = "integer", default = 3000L, dest = "read_len_min",
                help = "drop reads shorter than this [%default]"),
    make_option("--aln-sc-th",     type = "double", default = 20, dest = "aln_sc_th",
                help = "barcode alignment score threshold [%default]"),
    make_option("--map-start-th",  type = "integer", default = 1000L, dest = "map_start_th",
                help = "max distance from the subtelomere boundary [%default]"),
    make_option("--seq-ln",        type = "integer", default = 300L, dest = "seq_ln",
                help = "bases from each read end searched for the barcode [%default]"),
    make_option("--threads",       type = "integer", default = 4L,
                help = "threads passed to scanBam [%default]"),
    make_option("--emit-all",      action = "store_true", default = FALSE, dest = "emit_all",
                help = "also write reads.table.all.txt including filtered-out reads")
)
args <- parse_args(OptionParser(option_list = opts))

for (req in c("bam", "outdir", "exp_name")) {
    if (is.null(args[[req]])) stop("missing required argument: --", gsub("_", "-", req))
}

bamFN <- normalizePath(args$bam, mustWork = TRUE)

scriptDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(scriptDir) == 0 || scriptDir == "") scriptDir <- "."
source(file.path(scriptDir, "telobp_functions.R"))

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
outdir <- normalizePath(args$outdir)

## The pipeline functions write their PNG/TXT outputs into the working directory.
oldwd <- getwd(); setwd(outdir); on.exit(setwd(oldwd))

## Resolve the barcode. A literal sequence wins; otherwise fall back to the named DNAStrings
## defined in telobp_functions.R, but fail with something readable rather than R's default.
if (!is.null(args$barcode_seq)) {
    seqUp <- toupper(args$barcode_seq)
    if (!grepl("^[ACGT]+$", seqUp)) {
        stop("--barcode-seq must contain only A, C, G and T; got: ", args$barcode_seq)
    }
    tagSeq <- DNAString(seqUp)
} else if (exists(args$barcode_name)) {
    tagSeq <- get(args$barcode_name)
} else {
    stop("unknown --barcode-name '", args$barcode_name,
         "'. Known TeloTags: ", paste(grep("^NB[0-9]+uq$", ls(), value = TRUE), collapse = ", "),
         ". Pass an arbitrary barcode with --barcode-seq instead.")
}
tagName <- paste0("Samp.", args$barcode_name)
tags <- setNames(list(tagSeq), tagName)

cat("== Method 1: mapping-based telomere length ==\n")
cat("bam:      ", bamFN, "\n")
cat("outdir:   ", outdir, "\n")
cat("expName:  ", args$exp_name, "\n")
cat("chrArmLn: ", args$chr_arm_ln, "\n")
cat("barcode:  ", args$barcode_name, " (", as.character(tagSeq), ")\n\n", sep = "")

## Step 1: read BAM, compute per-read telomere length from softclip vs mapping position.
readsPrim <- processBamFn(bamFN = bamFN, expN = args$exp_name, chrArmLn = args$chr_arm_ln,
                          readLnMin = args$read_len_min, nThreads = args$threads)

## Step 2: align the TeloTag to read ends, filter, write table + plots.
readsPrim.dm <- alignWriteTableFn(reads = readsPrim, expN = args$exp_name, tags = tags,
                                  seqLn = args$seq_ln, alnScTh = args$aln_sc_th,
                                  mapStartTh = args$map_start_th, emitAll = args$emit_all)

## ---- stable filenames ----
## The upstream functions bake the experiment name, tag name and threshold into the filename.
## Rename to fixed names so the Snakefile can declare its outputs without globbing.
upstreamTbl <- paste0("reads.table_", args$exp_name, ".prim.", tagName,
                      "_alnScTh.", args$aln_sc_th, ".txt")
if (file.exists(upstreamTbl)) file.rename(upstreamTbl, "reads.table.txt")
upstreamAll <- paste0("reads.table.all_", args$exp_name, ".prim.", tagName,
                      "_alnScTh.", args$aln_sc_th, ".txt")
if (file.exists(upstreamAll)) file.rename(upstreamAll, "reads.table.all.txt")

## ---- consolidated summary across chromosome ends ----
if (file.exists("reads.table.txt")) {
    tbl <- read.delim("reads.table.txt", stringsAsFactors = FALSE)
    telCol <- grep("^Telomere", names(tbl), value = TRUE)[1]
    tbl$telomere_bp <- tbl[[telCol]]
    tbl$chrEnd <- gsub("_M$|_P$", "", as.character(tbl$chr))

    overall <- data.frame(
        metric = c("tagged_reads", "mean_telomere_bp", "median_telomere_bp",
                   "sd_telomere_bp", "min_telomere_bp", "max_telomere_bp"),
        value  = c(nrow(tbl),
                   round(mean(tbl$telomere_bp, na.rm = TRUE), 1),
                   round(median(tbl$telomere_bp, na.rm = TRUE), 1),
                   round(sd(tbl$telomere_bp, na.rm = TRUE), 1),
                   min(tbl$telomere_bp, na.rm = TRUE),
                   max(tbl$telomere_bp, na.rm = TRUE)))
    write.csv(overall, file = "summary.overall.csv", row.names = FALSE)

    perChr <- do.call(rbind, lapply(split(tbl, tbl$chrEnd), function(d) {
        data.frame(chrEnd = d$chrEnd[1], n = nrow(d),
                   median_telomere_bp = round(median(d$telomere_bp, na.rm = TRUE), 1),
                   mean_telomere_bp   = round(mean(d$telomere_bp, na.rm = TRUE), 1),
                   sd_telomere_bp     = round(sd(d$telomere_bp, na.rm = TRUE), 1))
    }))
    perChr <- perChr[order(gsub("chr", "", perChr$chrEnd)), ]
    write.csv(perChr, file = "summary.by_chr.csv", row.names = FALSE)

    cat("\n== Overall (Method 1) ==\n"); print(overall)
    cat("\n== Per chromosome end (Method 1) ==\n"); print(perChr, row.names = FALSE)
} else {
    cat("\nWARNING: no per-read table was produced",
        "\n(no reads passed the barcode + position filters?)\n")
}

cat("\n[run_analysis.R] done. Outputs in", outdir, "\n")
