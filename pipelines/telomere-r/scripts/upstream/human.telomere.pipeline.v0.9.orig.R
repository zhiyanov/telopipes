## == ## 3/29/2024 version 0.9
##       Andreas Rechtsteiner

## this script reads in a bam file
## bam file should be sorted and indexed
## map fastq file with telomere reads to a reference that contains the chr ends of the genome.
## reference genome needs to have telomere repeats cut at chosen subtelomere boundary
## all chr ends in the cut reference genome need to have the same length
## length of chr end needs to be passed to functions below

## make sure following packages are installed and load them:

library(Biostrings)
library(GenomicAlignments)
library(S4Vectors)
library(writexl)
library(RColorBrewer)

## == ## 24 bp unique Nanopore barcodes for demultiplexing and filtering for reads that have a barcode at the end
NB10uq <- DNAString("GAGAGGACAAAGGTTTCAACGCTT")
NB12uq <- DNAString("TCCGATTCTGCTTCTTTCTACCTG")
NB13uq <- DNAString("AGAACGACTTCCATACTCGTGTGA")
NB16uq <- DNAString("CGTCAACTGACAGTGGTTCGTACT")
NB50uq <- DNAString("ATGGACTTTGGTAACTTCCTGCGT")
NB65uq <- DNAString("TTCTCAGTCTTCCTCCAGACAAGG")
NB68uq <- DNAString("GAATCTAAGCAAACACGAAGGTGG")
NB88uq <- DNAString("TCTTCTACTACCGATCCGAAGCAG")
NB01uq <- DNAString("CACAAAGACACCGACAACTTTCTT")
NB02uq <- DNAString("ACAGACGACTACAAACGGAATCGA")
NB15uq <- DNAString("AGGTCTACCTCGCTAACACCACTG")
NB19uq <- DNAString("GTTCCTCGTGCAGTGTCAAGAGAT")
NB20uq <- DNAString("TTGCGTCCTGTTACGAGAACTCAT")
NB66uq <- DNAString("CCGATCCTTGTGGCTTCTAACTTC")
NB67uq <- DNAString("GTTTGTCATACTCGTGTGCTCACC")
NB69uq <- DNAString("TACAGTCCGAGCCTCATGTGATCT")
NB70uq <- DNAString("ACCGAGATCCTACGAATGGAGTGT")
NB72uq <- DNAString("TAGCTGACTGTCTTCCATACCGAC")

## == ##
## load all of the functions below:

## === === === === ##
## two functions are being called directly:
## 1. processBamFn() and
## 2. alignWriteTableFn():
## === === === === ##
## input:  bamFN     bam File Name (also needs index .bai file)
##         expN      experiment Name
##         chrArmLn  chr arm elngth in reference genome reads were mapped to
##         readLnMin filter out reads shorter than readLnMin
##         nThreads  number of threads for scanBam function (might not do anything for scanBam)
## output: object with telomere lengths based on mapping and reference

processBamFn <- function(bamFN, expN,  chrArmLn=500000, readLnMin=3000, nThreads=4){

    print("print bam input file name: "); print(bamFN)

    cat("\nRead in bam file\n")
    reads <- scanBam(bamFN, nThreads=nThreads)[[1]]

    ## adjust positions of q end chromosome arms (use negative coordinates for q arms)
    reads <- data.frame(posEndPQchrFn(reads,  chrArmLn))

    ## reorder mappings by chromosome name and mapStart position
    ind <- order(gsub("chr","", reads$chr), reads$mapStart)
    reads <- reads[ind, ]

    ## get indices for primary (secondary ... ) mappings:
    mapStats <- mapStatsFn(reads)

    cat("\nGet telomere length based on mapping of read relative to subtelomere boundary\n")
    reads <- getSclipTelPQchrFn(reads)

    ## histogram of read lengths of primary and unmapped reads:
    try({
        png(paste0("hist.readLn.prim.unmap_", gsub(" ", ".", expN), ".png"), w=640, h=480)
        ## par(mfrow=c(2,1), cex=1.4)
        par(mfrow=c(1,1), cex=1.4)

        xlim.max <- quantile(reads$qwidth[mapStats[["primary"]]], 0.99, na.rm=T)
        (hist.step <- 10*round(xlim.max/500))

        hist(tmp  <- reads$qwidth[mapStats[["primary"]]], seq(0,max(tmp, na.rm=TRUE)+hist.step, hist.step), main=paste0("read length primary map. ", expN), xlab="read length bp", xlim=c(0,xlim.max))
        abline(v=median(tmp, na.rm=TRUE), col="green3", lwd=2)
        axis(1, at=median(tmp, na.rm=TRUE), label=round(median(tmp, na.rm=TRUE)), col="green3", line=-1)

        ## hist(reads$qwidth[mapStats[["unmapped"]]], seq(0,max(tmp, na.rm=TRUE)+hist.step, hist.step), main=paste0("read length unmapped ", expN), xlab="read length bp", xlim=c(0,xlim.max))
        ## abline(v=median(tmp, na.rm=TRUE), col="red3", lwd=2)
        ## axis(1, at=median(tmp, na.rm=TRUE), label=round(median(tmp, na.rm=TRUE)), col="red3", line=-1)

        dev.off()
    })

    ## histogram mapq scores primary:
    try({
        png(paste0("hist.mapq_", gsub(" ", ".", expN), ".primary.map.png"), w=640, h=640)
        par(mfrow=c(1,1), cex=1.3)
        hist(reads[["mapq"]][mapStats[["primary"]]], seq(-0.5,61,1), main=paste0(expN, ", primary mappings"), xlab="mapq")
        dev.off()
    })

    ## plot hist of mapStart position relative to our subtelomere boundary, overlay all chromsomes:
    mapqCutL <- list("all_mapq"=list("index"= TRUE, "name"= "_all.mapq", "main" = "all mapq "),
                     "mapq>1"=list("index"=reads$mapq>1, "name"="_filt.mapq>1", "main"="mapq > 1 "),
                     "mapq>40"=list("index"=reads$mapq>40, "name"="_filt.mapq>40", "main"="mapq > 40 "))

    hist.step <- 1000
    try({
        for (mapqc in mapqCutL){
            png(paste0("hist.mapStart.pos.rel_prim_", gsub(" ", ".", expN), mapqc[["name"]], ".png"), w=640, h=480)
            par(mfrow=c(1,1), cex=1.3)
            hist(reads[["mapStart"]][mapqc[["index"]] & mapStats[["primary"]]], seq(-chrArmLn,chrArmLn,hist.step), main=paste0("primary map start pos ", expN, " ", mapqc[["main"]]), xlab="map start relative subtelomere (bp)", xlim=c(0, chrArmLn))
            dev.off()
        }
    })

    ## only look at primary mapped reads:
    reads.prim <- lapply(reads, function(x){x[mapStats$primary]})
    ## filter short reads:
    reads.prim <- lapply(reads.prim, function(x){x[reads.prim$qwidth >= readLnMin]})
    cat(paste0("\nNr mappings returned in reads.prim object: ", length(reads.prim$mapStart), "\n\n"))

    ## return primary mapped reads
    return(reads.prim)

}


## == ## function that calls routines to align to tags and write tables
##  input:
##       reads      all mapping info with telomere lengths
##       expN       name of experiment
##       tags       list(tagName=DNAStringsObject(tagSequence))
##       seqLn      end sequence of read to look for tag/barcode
##       alnScTh    alignment score threshold for barcode
##       mapStartTh mapping needs to start within this distance from subtelomere boundary
## outputs:
##       plots of mapping and telomere statistics
##       input object with added analysis results

alignWriteTableFn <- function(reads, expN, tags, seqLn=300L, alnScTh=20, mapStartTh=1000){

    for (tagN in names(tags)){

        cat("\nAlign reads from expN: ", expN, "\n\t to tag: ", tagN, "\n\t with alnSc thresh: ", alnScTh, "\n")
        ## align tag/barcode to mapped reads and return alignment results
        reads <- align2tagPQchrFn(reads=reads, expN=expN, tagN = tagN, tagSeq = tags[[tagN]], seqLn=seqLn, alnScTh=alnScTh)

        ## write table to txt file
        writeTelTableFn(reads=reads, expN=expN, tagN=tagN, fileNTagged=paste0("reads.table_", expN, ".prim.", tagN, "_alnScTh.", alnScTh, ".txt"), posTh = mapStartTh, alnScTh=alnScTh)

    }

    return(reads)

}

## === ## all functions being called in above two functions:

## == ## read in reads object
##       output a list of logical vectors indicating primary mappings, secondary mappings, unmapped reads:
mapStatsFn <- function(reads){

    ## identify primary alignments:
    primf <- reads$flag==0
    primr <- reads$flag==16
    sum(prim <- primf | primr)
    length(unique(reads[["qname"]][prim]))
    ## secondary alignments:
    sum(secf <- reads$flag==256)
    sum(secr <- reads$flag==272)
    sum(sec <- secf | secr)
    length(unique(reads[["qname"]][sec]))
    length(unique(reads[["qname"]][prim | sec]))
    ## supplementary alignments
    sum(supf <- reads$flag==2048)
    sum(supr <- reads$flag==2064)
    sum(sup <- supf | supr)
    length(unique(reads[["qname"]][sup]))
    length(unique(reads[["qname"]][prim | sup]))
    ## unmapped alignments
    sum(unm <- reads$flag==4)
    sum(unm | prim)
    ## either a read is unmapped or it has a primary mapping
    length(unique(reads[["qname"]][unm]))
    length(unique(reads[["qname"]][prim | unm]))
    (uniq.reads <- length(unique(reads$qname)))

    return(list("primary"=prim, "secondary"=sec, "suppplemental" = sup, "unmapped"=unm, "primary or unmapped"=prim | unm, "nr unique reads: " = uniq.reads))

}


## =============== ##
## input: reads object
##        chrArmLn: chr arm length from reference genome
##    extracts chromosome from bam rname object
##    calculates end position of mapping in reference coordinates
##    for mappings to q arm, convert to negative reference coordinates
##    Note that reference genome needs chr names in correct format used in script, e.g. chr11p, chr11q, ....
## output: reads object with above additions

posEndPQchrFn <- function(reads, chrArmLn=500000){

    ## add chr as list element:
    reads$chr <- gsub("_CP.*","", reads$rname)        ## adjust T2T_CHM13 genome chr naming
    reads$chr <- gsub("chr([0-9][q|p])$", "chr0\\1", reads$chr) ## if  chr1, chr2,  ... found make it into chr01, chr02, chr03, chr04, ....
    reads$chr <- gsub("ATERNAL.*","", reads$rname)    ## HG002 has PATERNAL/MATERNAL naming on chr ends, change to just "M" and "P"
    reads$chr <- as.factor(reads$chr)

    ## add posEnd for mapping from pos and cigar string:
    reads$posEnd <- reads[["pos"]] + cigarWidthAlongReferenceSpace(reads$cigar)

    qchr <- grepl("chr([0-9]*|X|Y)q", reads$chr)
    pchr <- grepl("chr([0-9]*|X|Y)p", reads$chr)

    ## make for right chr end mapped reads the reference coordinates negative:
    reads$pos[qchr]  <- reads$pos[qchr] - chrArmLn - 2
    reads$posEnd[qchr]  <- reads$posEnd[qchr] - chrArmLn - 2

    ## output closest mapping position relative to subtelomere boundary (in positive coordinates)
    reads$mapStart[qchr] <- abs(reads$posEnd)[qchr]
    reads$mapStart[pchr] <- abs(reads$pos)[pchr]

    return(reads)
}


## == ##
##  calculate telomere length from mapping and sofclipped bases
##       this can be repalced with TeloNP output which is based on the read sequence
##       input: reads object
##       output: reads object with telomere lengths

getSclipTelPQchrFn <- function(reads){

    ## cigar string:
    reads.cig <- reads[["cigar"]]
    length(reads.cig)
    head(reads.cig)
    ## Softclipped reads left
    head(reads$SoftClipp <- as.numeric(gsub("^([0-9]+)S.*", "\\1", reads.cig)))
    ## Softclipped reads right
    head(reads$SoftClipq <- as.numeric(gsub(".*[MIDNSHPX=]([0-9]+)S$", "\\1", reads.cig)))

    ## get q-chr and p-chr arms vectors
    sum(qchr <- grepl("chr([0-9]*|X|Y)q", reads$chr))
    sum(pchr <- grepl("chr([0-9]*|X|Y)p", reads$chr))

    ## to get telomere length from subtelomere I substract starting position of mapping:
    reads$Telomere[qchr] <- reads$SoftClipq[qchr] - abs(reads$posEnd)[qchr] + 1
    reads$Telomere[pchr] <- reads$SoftClipp[pchr] - abs(reads$pos)[pchr] + 1
    ## telomere length could be negative in that approach if start position of mapping is more internal than length of softclipped bases; likely wrong alignment. remove negtaive telomere lengths:
    reads$Telomere[reads$Telomere < 0] <- NA

    return(reads)

}


## == ## if you have a reads object with tags/barcodes,
##       this function can get the alignment scores for the barcodes
## inputs:
##    reads    object
##    expN     experiment name
##    tagSeq   tag/barcode sequences as DNAString object
##    tagN     name of tag/barcode
##    seqLn    how far in from chr end should be searched for tag
##    alnScTh  threshold on aligbment score.
##             For 24 bp barcodes and the current pairwise alignment parameters I mostly used 20
##    alnScMax for figure layouts, what is the expected maximum alignmetn score
##             (48 for 24bp barcodes and gapOpening and gapExtension of 2)
##    nrReadsNeeded  plot figures and safe alignment info if there are more than nrReadsNeeded mapped reads for that tag/barcode
## output:
##    return input object with added alignment paremters

align2tagPQchrFn <- function(reads, expN, tagSeq, tagN, seqLn=300L, alnScTh=20, alnScMax = 50, nrReadsNeeded = 1){

    ## pairwise alignment parameters:
    gapOpening = 2; gapExtension = 2
    mat <- NULL     ## not sure what Tom uses, might check his package

    print(paste0("tagSeq: ", tagSeq))

    ## pick out the first and last seqLn bp of reads on both left and right ends
    ## left end
    readLEnd <- DNAStringSet(reads$seq, start=1, end=seqLn)
    names(readLEnd) <- reads$qname
    ## right end
    readLn <- nchar(reads$seq)
    readREnd <- DNAStringSet(reads$seq, start=readLn-(seqLn-1), end=readLn)
    names(readREnd) <- reads$qname

    ## pairwise alignment to left side of read
    alnLEndRC <- pairwiseAlignment(pattern = readLEnd, subject = reverseComplement(tagSeq), type = "overlap", substitutionMatrix = mat, gapOpening = gapOpening, gapExtension = gapExtension)
    ## pairwise alignment to right side of read
    alnREnd <- pairwiseAlignment(pattern = readREnd, subject = tagSeq, type = "overlap", substitutionMatrix = mat, gapOpening = gapOpening, gapExtension = gapExtension)

    ## index vector for p and q arm of chromosmes
    ## format of chr name needs to be something like chr01q (chr01q_M for HG002 with maternal and paternal haplotypes)
    qchr <- grepl("chr([0-9]*|X|Y)q", reads$chr)
    pchr <- grepl("chr([0-9]*|X|Y)p", reads$chr)

    ## alignment scores of tag to left end of reads
    alnLEndRCLog <- score(alnLEndRC) >= alnScTh
    ## alignment scores of tag to right end of reads
    alnREndLog <- score(alnREnd) >= alnScTh

    ## are there more reads than nrReadsNeeded that have an alignemnt to the tag, otherwise we do not plot them
    if (sum(alnLEndRCLog | alnREndLog, na.rm=T) > nrReadsNeeded) {

        try({

            ## create hist of alignment scores:
            png(paste0("hist.aln.score_", expN, "_endln.", seqLn, "_aln2.", tagN, "_alnScTh.", alnScTh, ".png"), w=640, h=800)
            par(mfrow=c(2,1), cex=1.3)
            h=hist(score(alnLEndRC)[pchr], seq(0,alnScMax,2), xlab="align. score reads' left end", main=paste0(expN, " aligned to revC ", tagN)); abline(v=alnScTh,lt=2,col=2,lwd=2)
            nr = sum((left.Al <- score(alnLEndRC) >= alnScTh)[pchr]); perc = round(nr/sum(pchr), 3)*100L
            text(0.5*alnScMax, 0.6*max(h$counts),  paste0(nr, ", ", perc, "%"), cex=1.3, pos=4, col="red")

            h=hist(score(alnREnd)[qchr], seq(0,alnScMax,2), xlab="align. score reads' right end", main=paste0(expN, " aligned to ", tagN)); abline(v=alnScTh,lt=2,col=2,lwd=2)
            nr = sum((right.Al <- score(alnREnd) >= alnScTh)[qchr]); perc = round(nr/sum(qchr), 3)*100L
            text(0.5*alnScMax, 0.6*max(h$counts),  paste0(nr, ", ", perc, "%"), cex=1.3, pos=4, col="red")
            dev.off()

            ## assign alignment scores to reads object:
            reads[[tagN]][qchr] <- (score(alnREnd))[qchr]
            reads[[tagN]][pchr] <- (score(alnLEndRC))[pchr]
            head(reads[[tagN]])

            ## position in read coordinates for barcode alignment
            reads[[paste0(tagN,".alignStart")]][pchr] <- start(pattern(alnLEndRC)@range)[pchr]
            reads[[paste0(tagN,".alignStart")]][qchr] <- (abs(seqLn - end(pattern(alnREnd)@range))+1)[qchr]
            reads[[paste0(tagN,".alignEnd")]][pchr] <- end(pattern(alnLEndRC)@range)[pchr]
            reads[[paste0(tagN,".alignEnd")]][qchr] <- (abs(seqLn - start(pattern(alnREnd)@range))+1)[qchr]

        })

    }

    return(reads)

}


## ======== ##
## write out telomere length table and map statistics and make figures
## ======== ##

## input:
##   reads         object with all the data, created in previous analysis steps
##   expN          experiment Name added to created file names
##   tagN          name of tag or barcode used to retrieve info in reads object
##   fileNTagged:  name of table with map data and telomere length
##   posTh         reads need to map withing posTh of subtelomere boundary as indicated by end of chromosomes in reference genome
##   alnScTh       tag/barcode considered found in read if pairwise alignment score is above alnScTh
##   plChrStats    logical. If TRUE, some statistics of mappings by chr are plotted like reads by chr, telomere length by chr, sub-telomere length by chr
##
## output:
##   various plots, like telomere and fragment length histograms, telomere length by chr
##   txt file with table of mapping statistics and tag alignments incluiding telomere lengths
##   returns object containing same data written out to file

writeTelTableFn <- function(reads, expN, tagN, fileNTagged, posTh = 1000, alnScTh=20, plChrStats=TRUE){

    cat(paste0("\nNr mapped reads with alignment score above threshold ", alnScTh, ": ", sum(indAlnScTh <- reads[[tagN]]>=alnScTh, na.rm=TRUE)), sep="\n")
    indAlnScTh[is.na(indAlnScTh)] <- F

    ## mapping within posTh from subtelomere boundary
    indMapPosTh <- reads$mapStart <= posTh

    ## positive telomere lengths (should all be positive already):
    indPosTel <- reads$Tel > 0; indPosTel[is.na(indPosTel)] <- F;

    indTagged <- indAlnScTh & indMapPosTh & indPosTel
    indNonTagged <- !(indAlnScTh) & indMapPosTh & indPosTel

    ## Build table to write out and make plots
    try({
        telTaggedDf <- data.frame("qname"=reads$qname, "qwidth"=reads$qwidth, "chr"=reads$chr, "strand"=reads$strand, "posStart"=reads$pos, "posEnd"=reads$posEnd, "mapPosFromTel"= reads$mapStart, "mapq"=reads$mapq, "Telomere+Adap"= reads$Tel)[indTagged, ]
        telTaggedDf[paste0("alnSc.",tagN)] = reads[[tagN]][indTagged]
        telTaggedDf[paste0(tagN,".alignStart")] = reads[[paste0(tagN,".alignStart")]][indTagged]
        telTaggedDf[paste0(tagN,".alignEnd")] = reads[[paste0(tagN,".alignEnd")]][indTagged]
    })

    ## only write out if telTaggedDf has elements
    if (nrow(telTaggedDf) > 0){
        ## == ## just tagged and primary mapped frag lengths and telomeres
        try({

            write.table(telTaggedDf, file=fileNTagged, row.names = F, sep = "\t", quote=F)

            cat("\nPlot histogram of fragment lengths and genome-wide telomere lengths\n");

            png(paste0("hist.Frag.ln.and.Telo.ln.", expN, ".", tagN, "_alnScTh.", alnScTh, ".png"), w=640, h=800)

            par(mfrow=c(2,1), cex=1.3)
            xlim.max <- max(15000, quantile(telTaggedDf$qwidth, 0.99, na.rm=TRUE));
            (hist.step <- 10*round(xlim.max/500))   ## plot 50 steps for hstogram;
            ## Plot read length
            hist(tmp  <- telTaggedDf$qwidth, seq(0, max(tmp, na.rm=TRUE)+hist.step, hist.step), main=paste0("frag ln ", expN, " ", tagN), xlab="read frag length bp", xlim=c(0,xlim.max))
            abline(v=median(tmp, na.rm=TRUE), col="green3", lwd=2)
            axis(1, at=median(tmp, na.rm=TRUE), label=round(median(tmp, na.rm=TRUE)), col="green3", line=-1)
            ## plot Telomere length
            hist(tmp  <- telTaggedDf$Tel, seq(0, max(tmp, na.rm=T)+hist.step, hist.step), main=paste0("Telomere ln ", expN, " ", tagN), xlab="telomere length bp", xlim=c(0,xlim.max))
            abline(v=median(tmp, na.rm=TRUE), col="blue3", lwd=2)
            axis(1, at=median(tmp, na.rm=TRUE), label=round(median(tmp, na.rm=TRUE)), col="blue3", line=-1)

            dev.off()

            cat(paste0("\nWrite table and plot telomere length histograms for ", expN, " & tag ", tagN,  "\n"))

        })

        ## plot statistics by chromosome end:
        if (plChrStats){

            try({

                ## remove "chr" from chr names:
                chrs <- as.factor(gsub("chr", "", telTaggedDf$chr))
                ## if more than 46 chr ends, reads are likely mapped to both haplotypes
                ## in that case, plot two chr ends with same color, assuming
                ## maternal and paternal ends follow each other in the plot
                if (length(chrs) > 46) {
                    colVec <- rep(brewer.pal(12, "Set3"), each=2)
                } else {
                    colVec <- brewer.pal(12, "Set3")
                }

                cat("\nMake chrosomome specific Telomere plots.\n")
                ## various mapq cutoffs
                mapqCutL <- list("all_mapq"=list("index"= TRUE, "name"= "_all.mapq", "main" = "all mapq "),
                                 "mapq>1"=list("index"=telTaggedDf$mapq>1, "name"="_filt.mapq>1", "main"="mapq > 1 "),
                                 "mapq>5"=list("index"=telTaggedDf$mapq>5, "name"="_filt.mapq>5", "main"="mapq > 5 "),
                                 "mapq>40"=list("index"=telTaggedDf$mapq>40, "name"="_filt.mapq>40", "main"="mapq > 40 "))

                for (mapqc in mapqCutL){

                    ## nr reads by chr end
                    png(paste0("barpl.chr.tagged.", gsub(" ", ".", expN), "_aln2.", tagN, mapqc[["name"]], ".v3.png"), w=length(unique(reads$chr))*14, h=640)
                    barplot((table(c(chrs[mapqc[["index"]]], unique(chrs)))-1), main=paste("Nr tagged reads", expN, tagN, mapqc[["main"]]), col=colVec, las=2)
                    dev.off()

                    ## telomere length by boxplot
                    png(paste0("boxpl.Tel.by.chr.tagged.", gsub(" ", ".", expN), "_aln2.", tagN, mapqc[["name"]], ".v2notch.png"), w=length(unique(reads$chr))*14, h=640)
                    boxplot(split(telTaggedDf$Tel[mapqc[["index"]]], chrs[mapqc[["index"]]]), main=paste("Tel", expN, tagN, mapqc[["main"]]), col=colVec, las=2, notch=T, ylab="Telomere length in bp")
                    dev.off()

                    ## non-telomere length boxplot by chr
                    png(paste0("boxpl.non-Tel.by.chr.tagged.", gsub(" ", ".", expN), "_aln2.", tagN, mapqc[["name"]], ".v2notch.png"), w=length(unique(reads$chr))*14, h=640)
                    nonTelLn <- telTaggedDf$qwidth - telTaggedDf$Tel
                    boxplot(split(nonTelLn[mapqc[["index"]]], chrs[mapqc[["index"]]]), main=paste("Non-Tel ", expN, tagN, mapqc[["main"]]), col=colVec, las=2, notch=T, ylab="Non-Telomere length in bp")

                    dev.off()

                }
            })
        }

    } else {

        cat(paste0(expN, " for tag ", tagN, " did not have reads, nothing to write out \n\n"))

    }

    return(telTaggedDf)

}



## ## ======================== ##
## ## how to run the functions:
## ## ======================== ##

chrArmLn <- 500000L       ## uploaded CHM13 and HG002 genomes were cut with TeloBP and both 500,000 bp chr ends
## barcodes name and sequence:
tags <- list("Samp1.NB50uq"=NB50uq, "Samp2.NB88uq"=NB88uq)
## bam file name:
bamFN <- "file.sort.bam"
## experiment name used in names of output files:
expName <- "test.data"

## read bam file and process.
## return primary mapped reads with telomere length and various statistics
readsPrim <- processBamFn(bamFN=bamFN, expN=expName,  chrArmLn=chrArmLn)
## returned fields in list object:
str(readsPrim)

## align reads to tag or barcode
## demultiplex and write out a table for each barcode
## make figures
readsPrim.dm <- alignWriteTableFn(reads=readsPrim, expN=expName, tags=tags, alnScTh=20, mapStartTh=1000)
## see objects in list:
str(readsPrim.dm)

