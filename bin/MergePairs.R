#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("--errFor"), type="character", default=NULL, help="Forward RDS file"),
    make_option(c("--errRev"), type="character", default=NULL, help="Reverse RDS file"),
    make_option(c("--cpus"), type="numeric", , default=1, help="cpus"),
    make_option(c("--pool"), type="character", help="Pooling"),

    make_option(c("--minOverlap"), type="numeric", help="Min overlap"),
    make_option(c("--maxMismatch"), type="numeric", help="Max mismatch"),
    make_option(c("--trimOverhang"), default=FALSE, action="store_true", help="Trim overhanging sequence if overlapping"),
    make_option(c("--justConcatenate"), default=FALSE, action="store_true", help="Just concatenate sequences")
)

opt <- parse_args(OptionParser(option_list=option_list))

environment(mergePairsRescue) <- asNamespace('dada2')

filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)

errF <- readRDS(opt$errFor)
errR <- readRDS(opt$errRev)
cat("Processing all samples\n")

#Variable selection from CLI input flag --pool
pool <- opt$pool
if(pool == "T" || pool == "TRUE"){
  pool <- as.logical(pool)
}

derepFs <- derepFastq(filtFs)

ddFs <- dada(derepFs, err=errF, multithread=opt$cpus, pool=pool)

derepRs <- derepFastq(filtRs)

ddRs <- dada(derepRs, err=errR, multithread=opt$cpus, pool=pool)

mergers <- mergePairs(ddFs, derepFs, ddRs, derepRs,
    returnRejects = TRUE,
    minOverlap = opt$minOverlap,
    maxMismatch = opt$maxMismatch,
    trimOverhang = as.logical(opt$trimOverhang),
    justConcatenate = as.logical(opt$justConcatenate)
    )

# TODO: make this a single item list with ID as the name, this is lost
# further on
saveRDS(mergers, "all.mergers.RDS")

saveRDS(ddFs, "all.ddF.RDS")
saveRDS(derepFs, "all.derepFs.RDS")

saveRDS(ddRs, "all.ddR.RDS")
saveRDS(derepRs, "all.derepRs.RDS")
