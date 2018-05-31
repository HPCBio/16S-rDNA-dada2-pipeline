library(dada2); packageVersion("dada2")
    out <- filterAndTrim(fwd = "${reads[0]}",
                        filt = paste0("${pairId}", ".R1.filtered.fastq.gz"),
                        rev = "${reads[1]}",
                        filt.rev = paste0("${pairId}", ".R2.filtered.fastq.gz"),
                        trimLeft = c(${params.trimFor},${params.trimRev}),
                        truncLen = c(${params.truncFor},${params.truncRev}),
                        maxEE = c(${params.maxEEFor},${params.maxEERev}),
                        truncQ = ${params.truncQ},
                        maxN = ${params.maxN},
                        rm.phix = ${params.rmPhiX},
                        maxLen = ${params.maxLen},
                        minLen = ${params.minLen},
                        compress = TRUE,
                        verbose = TRUE,
                        multithread = ${task.cpus})
    write.csv(out, paste0("${pairId}", ".trimmed.txt"))
