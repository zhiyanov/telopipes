## Method 1 driver: mapping-based telomere length for HG002 / NB65uq.
## Sources the GreiderLab pipeline functions and runs them with the parameters
## for this dataset. Runs inside the telomere-r Docker image.
##
## Usage:
##   Rscript run_analysis.R <sorted_indexed.bam> <outdir> <expName> [chrArmLn] [barcode]
## Defaults: chrArmLn=500000, barcode=NB65uq

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript run_analysis.R <bam> <outdir> <expName> [chrArmLn] [barcode]")
}
bamFN    <- normalizePath(args[[1]])
outdir   <- args[[2]]
expName  <- args[[3]]
chrArmLn <- if (length(args) >= 4) as.integer(args[[4]]) else 500000L
barcodeN <- if (length(args) >= 5) args[[5]] else "NB65uq"

scriptDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(scriptDir) == 0 || scriptDir == "") scriptDir <- "."
source(file.path(scriptDir, "telobp_functions.R"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outdir <- normalizePath(outdir)

## The pipeline functions write their PNG/TXT outputs into the working dir.
oldwd <- getwd(); setwd(outdir); on.exit(setwd(oldwd))

## Barcode for this sample. The functions file defines NB65uq (and the rest) as
## DNAString objects; pick the requested one by name.
tagSeq <- get(barcodeN)
tags <- setNames(list(tagSeq), paste0("Samp.", barcodeN))

cat("== Method 1: mapping-based telomere length ==\n")
cat("bam:      ", bamFN, "\n")
cat("outdir:   ", outdir, "\n")
cat("expName:  ", expName, "\n")
cat("chrArmLn: ", chrArmLn, "\n")
cat("barcode:  ", barcodeN, " (", as.character(tagSeq), ")\n\n", sep = "")

## Step 1: read BAM, compute per-read telomere length from softclip vs mapping
## position, keep primary mappings >= 3 kb.
readsPrim <- processBamFn(bamFN = bamFN, expN = expName, chrArmLn = chrArmLn,
                          readLnMin = 3000, nThreads = 4)

## Step 2: align the TeloTag barcode to read ends, filter (mapStart <=1kb of the
## subtelomere boundary, alnScore>=20, positive telomere), write table + plots.
readsPrim.dm <- alignWriteTableFn(reads = readsPrim, expN = expName, tags = tags,
                                  seqLn = 300L, alnScTh = 20, mapStartTh = 1000)

## ---- consolidated summary across chromosome ends ----
## Reload the per-read table the pipeline just wrote for a tidy summary.
tagName <- names(tags)[1]
tblFN <- paste0("reads.table_", expName, ".prim.", tagName, "_alnScTh.20.txt")
if (file.exists(tblFN)) {
    tbl <- read.delim(tblFN, stringsAsFactors = FALSE)
    telCol <- grep("^Telomere", names(tbl), value = TRUE)[1]
    tbl$telomere_bp <- tbl[[telCol]]
    ## chromosome-end label without haplotype/strand suffix noise
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
    write.csv(overall, file = paste0("summary.overall_", expName, ".csv"), row.names = FALSE)

    perChr <- do.call(rbind, lapply(split(tbl, tbl$chrEnd), function(d) {
        data.frame(chrEnd = d$chrEnd[1], n = nrow(d),
                   median_telomere_bp = round(median(d$telomere_bp, na.rm = TRUE), 1),
                   mean_telomere_bp   = round(mean(d$telomere_bp, na.rm = TRUE), 1),
                   sd_telomere_bp     = round(sd(d$telomere_bp, na.rm = TRUE), 1))
    }))
    perChr <- perChr[order(gsub("chr", "", perChr$chrEnd)), ]
    write.csv(perChr, file = paste0("summary.by_chr_", expName, ".csv"), row.names = FALSE)

    cat("\n== Overall (Method 1) ==\n"); print(overall)
    cat("\n== Per chromosome end (Method 1) ==\n"); print(perChr, row.names = FALSE)
} else {
    cat("\nWARNING: expected per-read table not found:", tblFN,
        "\n(no reads passed the barcode + position filters?)\n")
}

cat("\n[run_analysis.R] done. Outputs in", outdir, "\n")
