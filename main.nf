#!/usr/bin/env nextflow

// Some help with time stamps
import java.text.SimpleDateFormat

// version
version = 0.4

timestamp = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())

// Pass this in to run
params.reads = "./raw-seq/*_R{1,2}.fastq.gz"
params.outdir = "./" + timestamp + "-dada2"
params.ticket = 0

// Trimming
params.trimFor = false
params.trimRev = false
params.truncFor = 0
params.truncRev = 0
params.maxEEFor = 2
params.maxEERev = 2
params.truncQ = 2
params.maxN = 0
params.maxLen = "Inf"
params.minLen = 20
params.rmPhiX = "FALSE"

// Merging
params.minOverlap = 20
params.maxMismatch = 0
params.trimOverhang = "FALSE" // This should be true with some sequences (V4)

params.reference = false
params.species = false

// NYI, for dada sample inference pooling (requires all samples)
params.pool = false

if ( params.trimFor == false ) {
    exit 1, "Must set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)"
}

if ( params.trimRev == false ) {
    exit 1, "Must set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)"
}

if ( params.reference == false ) {
    exit 1, "Must set reference database using --reference"
}

// if (params.ticket == 0) exit 1, "Must set Redmine ticket for pipeline summary to be sent"

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { dada2ReadPairsToQual; dada2ReadPairs }

refFile = file(params.reference)

// Test R library
dada2Mod = 'R-lib/3.4.2'

myQueue = 'normal'

// TODO: maybe have a way to check the params and fill this out automatically?
runInfo = """

====================================
 dada2 : Paired-end workflow
====================================
Pipeline ver : ${version}
Reads        : ${params.reads}
Ticket       : ${params.ticket}
dada2        : ${dada2Mod}
Trim-For     : ${params.trimFor}
Trim-Rev     : ${params.trimRev}
Trunc-For    : ${params.truncFor}
Trunc-Rev    : ${params.truncRev}
Reference    : ${params.reference}
Species      : ${params.species}
Current home : $HOME
Current user : $USER
Current path : $PWD
Script dir   : $baseDir
Working dir  : $workDir
Output dir   : ${params.outdir}
====================================
"""

/*
 *
 * Step 1: Filter and trim (run per sample?)
 *
 */

// TODO: Note we need to hard trim reads to remove the primers at the 5' end,
// these mess with dada2 (overpredict chimeras)

process plotQual {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "12 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "link"

    input:
    file allReads from dada2ReadPairsToQual.flatMap({ it[1] }).collect()

    output:
    file "R1.pdf" into forQualPDF
    file "R2.pdf" into revQualPDF

    script:
    """
    #!/usr/bin/env Rscript
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
    """
}

process filterAndTrim {
    cpus 4
    executor 'slurm'
    queue myQueue
    memory "12 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "link"

    input:
    set pairId, file(reads) from dada2ReadPairs

    output:
    set val(pairId), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" into filteredReads
    file "*.R1.filtered.fastq.gz" into forReads
    file "*.R2.filtered.fastq.gz" into revReads
    file "*.trimmed.txt" into trimTracking

    script:
    """
    #!/usr/bin/env Rscript
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
    """
}

process mergeTrimmedTable {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "link"

    input:
    file trimData from trimTracking.collect()

    output:
    file "all.trimmed.csv" into trimmedReadTracking

    script:
    """
    #!/usr/bin/env Rscript
    trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
    sample.names <- sub('.trimmed.txt', '', trimmedFiles)
    trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
    colnames(trimmed)[1] <- "Sequence"
    trimmed\$SampleID <- sample.names
    write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)
    """
}

/*
 *
 * Step 2: Learn error rates (run on all samples)
 *
 */

// TODO: combine For and Rev process to reduce code duplication?

process LearnErrorsFor {
    cpus 8
    executor 'slurm'
    queue myQueue
    memory "12 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "link"

    input:
    file fReads from forReads.collect()

    output:
    file "errorsF.RDS" into errorsFor

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2);
    packageVersion("dada2")

    # File parsing
    filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
    sample.namesF <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn forward error rates
    errF <- learnErrors(filtFs, nread=1e6, multithread=${task.cpus})
    pdf("R1.err.pdf")
    plotErrors(errF, nominalQ=TRUE)
    dev.off()
    saveRDS(errF, "errorsF.RDS")
    """
}

process LearnErrorsRev {
    cpus 8
    executor 'slurm'
    queue myQueue
    memory "12 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "link"

    input:
    file rReads from revReads.collect()

    output:
    file "errorsR.RDS" into errorsRev

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2);
    packageVersion("dada2")

    # load error profiles

    # File parsing
    filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)
    sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn forward error rates
    errR <- learnErrors(filtRs, nread=1e6, multithread=${task.cpus})
    pdf("R2.err.pdf")
    plotErrors(errR, nominalQ=TRUE)
    dev.off()
    saveRDS(errR, "errorsR.RDS")
    """
}

/*
 *
 * Step 3: Dereplication, Sample Inference, Merge Pairs
 *
 */

// TODO: allow serial processing of this step?

process SampleInferDerepAndMerge {
    cpus 4
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-Derep", mode: "link"

    input:
    set val(pairId), file(filtFor), file(filtRev) from filteredReads
    file errFor from errorsFor
    file errRev from errorsRev

    output:
    file "*.merged.RDS" into mergedReads
    file "*.ddF.RDS" into dadaFor
    file "*.ddR.RDS" into dadaRev

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    packageVersion("dada2")

    errF <- readRDS("${errFor}")
    errR <- readRDS("${errRev}")
    cat("Processing:", "${pairId}", "\\n")

    derepF <- derepFastq("${filtFor}")
    ddF <- dada(derepF, err=errF, multithread=${task.cpus})

    derepR <- derepFastq("${filtRev}")
    ddR <- dada(derepR, err=errR, multithread=${task.cpus})

    merger <- mergePairs(ddF, derepF, ddR, derepR,
        minOverlap = ${params.minOverlap},
        maxMismatch = ${params.maxMismatch},
        trimOverhang = TRUE
        )

    # TODO: make this a single item list with ID as the name, this is lost
    # further on
    saveRDS(merger, paste("${pairId}", "merged", "RDS", sep="."))
    saveRDS(ddF, paste("${pairId}", "ddF", "RDS", sep="."))
    saveRDS(ddR, paste("${pairId}", "ddR", "RDS", sep="."))
    """
}

// TODO: step may be obsolete if we run the above serially

process mergeDadaRDS {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-Inference", mode: "link"

    input:
    file ddFs from dadaFor.collect()
    file ddRs from dadaRev.collect()

    output:
    file "all.ddF.RDS" into dadaForReadTracking
    file "all.ddR.RDS" into dadaRevReadTracking

    script:
    '''
    #!/usr/bin/env Rscript
    library(dada2)
    packageVersion("dada2")

    dadaFs <- lapply(list.files(path = '.', pattern = '.ddF.RDS$'), function (x) readRDS(x))
    names(dadaFs) <- sub('.ddF.RDS', '', list.files('.', pattern = '.ddF.RDS'))
    dadaRs <- lapply(list.files(path = '.', pattern = '.ddR.RDS$'), function (x) readRDS(x))
    names(dadaRs) <- sub('.ddR.RDS', '', list.files('.', pattern = '.ddR.RDS'))
    saveRDS(dadaFs, "all.ddF.RDS")
    saveRDS(dadaRs, "all.ddR.RDS")
    '''
}

/*
 *
 * Step 4: Construct sequence table
 *
 */

process SequenceTable {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-SeqTable", mode: "link"

    input:
    file mr from mergedReads.collect()

    output:
    file "seqtab.RDS" into seqTable
    file "mergers.RDS" into mergerTracking

    script:
    '''
    #!/usr/bin/env Rscript
    library(dada2)
    packageVersion("dada2")

    mergerFiles <- list.files(path = '.', pattern = '.*.RDS$')
    pairIds <- sub('.merged.RDS', '', mergerFiles)
    mergers <- lapply(mergerFiles, function (x) readRDS(x))
    names(mergers) <- pairIds
    seqtab <- makeSequenceTable(mergers)

    saveRDS(seqtab, "seqtab.RDS")
    saveRDS(mergers, "mergers.RDS")
    '''
}

/*
 *
 * Step 8: Remove chimeras
 *
 */

if (params.species) {

    speciesFile = file(params.species)
    process ChimeraTaxonomySpecies {
        cpus 24
        executor 'slurm'
        queue myQueue
        memory "48 GB"
        module dada2Mod
        publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "link"

        input:
        file st from seqTable
        file ref from refFile
        file sp from speciesFile

        output:
        file "seqtab_final.RDS" into seqTableFinal,seqTableFinalTree,seqTableFinalTracking
        file "tax_final.RDS" into taxFinal

        script:
        """
        #!/usr/bin/env Rscript
        library(dada2)
        packageVersion("dada2")

        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=${task.cpus})

        # Assign taxonomy
        tax <- assignTaxonomy(seqtab, "${ref}", multithread=${task.cpus})
        tax <- addSpecies(tax, "${sp}")

        # Write to disk
        saveRDS(seqtab, "seqtab_final.RDS")
        saveRDS(tax, "tax_final.RDS")
        """
    }

} else {

    process ChimeraTaxonomy {
        cpus 24
        executor 'slurm'
        queue myQueue
        memory "48 GB"
        module dada2Mod
        publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "link"

        input:
        file st from seqTable
        file ref from refFile

        output:
        file "seqtab_final.RDS" into seqTableFinal,seqTableFinalTree,seqTableFinalTracking
        file "tax_final.RDS" into taxFinal

        script:
        """
        #!/usr/bin/env Rscript
        library(dada2)
        packageVersion("dada2")

        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=${task.cpus})

        # Assign taxonomy
        tax <- assignTaxonomy(seqtab, "${ref}", multithread=${task.cpus})

        # Write to disk
        saveRDS(seqtab, "seqtab_final.RDS")
        saveRDS(tax, "tax_final.RDS")
        """
    }
}

/*
 *
 * Step 9: Construct phylogenetic tree
 *
 */

// TODO: break into more steps?  phangorn takes a long time...

process AlignAndGenerateTree {
    cpus 12
    executor 'slurm'
    queue myQueue
    memory "12 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-Alignment", mode: "link"

    input:
    file sTable from seqTableFinalTree

    output:
    file "aligned_seqs.fasta" into alnFile
    file "phangorn.tree.RDS" into treeRDS
    file "tree.newick" into treeFile
    file "tree.GTR.newick" into treeGTRFile

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(DECIPHER)
    library(phangorn)

    seqs <- getSequences(readRDS("${sTable}"))
    names(seqs) <- seqs # This propagates to the tip labels of the tree
    alignment <- AlignSeqs(DNAStringSet(seqs),
                           anchor=NA,
                           processors = ${task.cpus})
    writeXStringSet(alignment, "aligned_seqs.fasta")

    # TODO: optimize this, or maybe split into a second step?
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit = pml(treeNJ, data=phang.align)
    write.tree(fit\$tree, file = "tree.newick")

    ## negative edges length changed to 0!
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))
    saveRDS(fitGTR, "phangorn.tree.RDS")
    write.tree(fitGTR\$tree, file = "tree.GTR.newick")
    """
}

process BiomFile {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-BIOM", mode: "link"

    input:
    file sTable from seqTableFinal
    file tTable from taxFinal

    output:
    file "dada2.biom" into biomFile

    script:
    """
    #!/usr/bin/env Rscript
    library(biomformat)
    packageVersion("biomformat")
    seqtab <- readRDS("${sTable}")
    taxtab <- readRDS("${tTable}")
    st.biom <- make_biom(t(seqtab), observation_metadata = taxtab)
    write_biom(st.biom, "dada2.biom")
    """
}

/*
 *
 * Step 10: Track reads
 *
 */

// Broken: needs a left-join on the initial table

process ReadTracking {
    cpus 2
    executor 'slurm'
    queue myQueue
    memory "8 GB"
    module dada2Mod
    publishDir "${params.outdir}/dada2-ReadTracking", mode: "link"

    input:
    file trimmedTable from trimmedReadTracking
    file sTable from seqTableFinalTracking
    file mergers from mergerTracking
    file ddFs from dadaForReadTracking
    file ddRs from dadaRevReadTracking

    output:
    file "all.readtracking.txt"

    script:
    """
    #!/usr/bin/env Rscript

    library(dada2)
    packageVersion("dada2")
    library(dplyr)

    getN <- function(x) sum(getUniques(x))

    dadaFs <- as.data.frame(sapply(readRDS("${ddFs}"), getN))
    dadaFs\$SampleID <- rownames(dadaFs)

    dadaRs <- as.data.frame(sapply(readRDS("${ddRs}"), getN))
    dadaRs\$SampleID <- rownames(dadaRs)

    mergers <- as.data.frame(sapply(readRDS("${mergers}"), getN))
    mergers\$SampleID <- rownames(mergers)

    seqtab.nochim <- as.data.frame(rowSums(readRDS("${sTable}")))
    seqtab.nochim\$SampleID <- rownames(seqtab.nochim)

    trimmed <- read.csv("${trimmedTable}")

    track <- Reduce(function(...) merge(..., by = "SampleID"),  list(trimmed, dadaFs, dadaRs, mergers, seqtab.nochim))
    colnames(track) <- c("SampleID", "SequenceR1", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    write.table(track, "all.readtracking.txt", sep = "\t", row.names = FALSE)
    """
}

workflow.onComplete {
    def subject = "[Task #${params.ticket}] BS-Seq aligment and methylaton calling pipeline"
    def recipient = 'hpcbiohelp@igb.illinois.edu'

    finalLog = """
Pipeline execution summary
---------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}
---------------------------

"""

    // ['mail', '-s', subject, recipient].execute() << "${finalLog}\n${runInfo}"
}
