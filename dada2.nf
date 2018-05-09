#!/usr/bin/env nextflow
/*
* USAGE: nextflow run foo.nf -qs 8
* Note that "-qs" is similar to "#PBS -t" and will only run a specified # of jobs at a time.
* This script creates hard links to data that exists in nextflow's work directory for everything
*/

/*
 * Set parameter values here. Only change the values within quotations
 * Many of the values already present will be okay to use, but make sure the pathes
 * are correct for your dataset.
*/

/* Path to base project results directory */
/* TODO: Need exception check? */

import java.text.SimpleDateFormat

// version
version = 0.3

timestamp = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())

// Pass this in to run
params.reads = "./raw-seq/*_R{1,2}.fastq.gz"
params.outdir = "./" + timestamp + "-dada2"
params.ticket = 0

// removing primers
params.trimFor = false
params.trimRev = false

// truncating final length
params.truncFor = 250
params.truncRev = 250

params.reference = false
params.species = false

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

    # TODO: add forward and reverse hard-trimming (see trimLeft)
    out <- filterAndTrim(fwd="${reads[0]}", filt=paste0("${pairId}", ".R1.filtered.fastq.gz"),
                  rev="${reads[1]}", filt.rev=paste0("${pairId}", ".R2.filtered.fastq.gz"),
                  trimLeft = c(${params.trimFor},${params.trimRev}),
                  truncLen=c(${params.truncFor},${params.truncRev}),
                  maxEE=2,
                  truncQ=11,
                  maxN=0,
                  compress=TRUE,
                  verbose=TRUE,
                  multithread=${task.cpus})

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
    file "*.trimmed.txt" from trimTracking.collect()

    output:
    file "all.trimmed.txt"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr); packageVersion("dplyr")

    # TODO: add forward and reverse hard-trimming (see trimLeft)
    trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
    trimmed <- lapply(trimmedFiles, function (x) read.csv(x))
    all.trimmed <- bind_rows(trimmed)
    colnames(all.trimmed)[1] <- "SampleID"

    write.table(all.trimmed, "all.trimmed.txt",
        row.names = FALSE,
        quote = FALSE,
        sep="\t")
    """
}

// TODO: process to merge individual tables

/*
 *
 * Step 2: Learn error rates (run on all samples)
 *
 */

// TODO: combine For and Rev process to reduce code duplication

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

    merger <- mergePairs(ddF, derepF, ddR, derepR)

    # TODO: make this a single item list with ID as the name, this is lost
    # further on
    saveRDS(merger, paste("${pairId}", "merged", "RDS", sep="."))

    # Construct sequence table and remove chimeras
    # seqtab <- makeSequenceTable(mergers)
    # saveRDS(seqtab, "/path/to/run1/output/seqtab.rds") # CHANGE ME to where you want sequence table saved
    """
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
    file "mergers.RDS" into mergers

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
        file "seqtab_final.RDS" into seqTableFinal
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
        file "seqtab_final.RDS" into seqTableFinal
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
 * Step 9: Track reads
 *
 */

// TODO: This step doesn't seem to be included in the example dada2 PE Big Data template


/*
 *
 * Step 10: Assign taxonomy
 *
 */

// TODO: Note this is in the above step #8, can it be split out?

/*
 *
 * Step 11: Evaluate accuracy
 *
 */

// TODO: This step doesn't seem to be included in the example dada2 PE Big Data template


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
