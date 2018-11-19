#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(optparse))


option_list = list(
    make_option(c("--id"), type="character", default=NULL, help="Sample ID"),
    make_option(c("--fwd"), type="character", default=NULL, help="Forward/Read 1"),
    make_option(c("--rev"), type="character", default=NULL, help="Rev/Read 2"),
    make_option(c("--cpus"), type="numeric", , default=1, help="cpus"),
    make_option(c("--trimFor"), type="numeric", , help="trim from 5' for R1"),
    make_option(c("--trimRev"), type="numeric", , help="trim from 5' for R2"),
    make_option(c("--truncFor"), type="numeric", , help="truncate R1 to this"),
    make_option(c("--truncRev"), type="numeric", , help="trucate R2 to this"),
    make_option(c("--maxEEFor"), type="numeric", , help="maxEE for R1"),
    make_option(c("--maxEERev"), type="numeric", , help="maxEE for R2"),
    make_option(c("--truncQ"), type="numeric", , help="Quality score to truncate to"),
    make_option(c("--maxN"), type="numeric", help="Maximum Ns allowed"),
    make_option(c("--maxLen"), type="numeric", help="Max length"),
    make_option(c("--minLen"), type="numeric", help="Min length"),
    make_option(c("--rmPhiX"), default=TRUE, action="store_true", help="Remove possible PhiX sequences")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

tmpFunc <- function(...) {
    print(as.list(match.call()))
    filterAndTrim(...)
}

out <- filterAndTrim(fwd        = opt$fwd,
                    filt        = paste0(opt$id, ".R1.filtered.fastq.gz"),
                    rev         = opt$rev,
                    filt.rev    = paste0(opt$id, ".R2.filtered.fastq.gz"),
                    trimLeft    = c(opt$trimFor,  opt$trimRev),
                    truncLen    = c(opt$truncFor, opt$truncRev),
                    maxEE       = c(opt$maxEEFor, opt$maxEERev),
                    truncQ      = opt$truncQ,
                    maxN        = opt$maxN,
                    rm.phix     = opt$rmPhiX,
                    maxLen      = opt$maxLen,
                    minLen      = opt$minLen,
                    compress    = TRUE,
                    verbose     = TRUE,
                    multithread = opt$cpus
                    )
write.csv(out, paste0(opt$id, ".trimmed.txt"))
