library(dada2); packageVersion("dada2")
    # Forward Reads
    pdf("R1.pdf")
    fnFs <- list.files('.', pattern="_R1_*.fastq*", full.names = TRUE)
    plotQualityProfile(fnFs, aggregate = TRUE)
    dev.off()
    # Reverse Reads
    pdf("R2.pdf")
    fnRs <- list.files('.', pattern="_R2_*.fastq*", full.names = TRUE)
    plotQualityProfile(fnRs, aggregate = TRUE)
    dev.off()
