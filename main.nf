#!/usr/bin/env nextflow
/*
========================================================================================
               D A D A 2   P I P E L I N E
========================================================================================
 DADA2 NEXTFLOW PIPELINE FOR UCT CBIO
 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ===================================
     uct-cbio/16S-rDNA-dada2-pipeline  ~  version ${params.version}
    ===================================
    Usage:
    nextflow run uct-cbio/16S-rDNA-dada2-pipeline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex
    The typical command for running the pipeline with command line flags is as follows:
    nextflow run -c <dada2.conf>  <dada2.nf> -profile uct_hex
    where: 
    dada2.conf is the configuration file
    dada2.nf   is the pipeline script
    
    To override existing values from the command line, please type these parameters:
    
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' - create your own if necessary
                                    NB -profile should always be specified on the command line, not in the config file
      --trimFor                     integer. headcrop of read1 (set 0 if no trimming is needed)
      --trimRev                     integer. headcrop of read2 (set 0 if no trimming is needed)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz)
    
    All available read preparation parameters:
      --trimFor                     integer. headcrop of read1
      --trimRev                     integer. headcrop of read2
      --truncFor                    integer. tailcrop of read1. enforced before trimming
      --truncRev                    integer. tailcrop of read2. enforced before trimming
      --maxEEFor                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --maxEERev                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --truncQ                      integer. Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2
      --maxN                        integer. Discard reads with more than maxN number of Ns in read; default=0
      --maxLen                      integer. maximum length of sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum)
      --minLen                      integer. minLen is enforced after trimming and truncation; default=50
      --rmPhiX                      {"T","F"}. remove PhiX from read              
      --minOverlap                  integer. minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
      --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. 
                                    "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction                                               primer region. Default="F" (false)
  
    Other arguments:
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are                                         
                                    pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run                                     
                                    sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    
     Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.email = false
params.plaintext_email = false

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//Validate inputs
if ( params.trimFor == false ) {
    exit 1, "Must set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)"
}

if ( params.trimRev == false ) {
    exit 1, "Must set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)"
}

if ( params.reference == false ) {
    exit 1, "Must set reference database using --reference"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { dada2ReadPairsToQual; dada2ReadPairs }

refFile = file(params.reference)

// Header log info
log.info "==================================="
log.info " uct-cbio/16S-rDNA-dada2-pipeline  ~  version ${params.version}"
log.info "==================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['trimFor'] = params.trimFor
summary['trimRev'] = params.trimRev
summary['truncFor'] = params.truncFor
summary['truncRev'] = params.truncRev
summary['truncQ'] = params.truncQ
summary['maxEEFor'] = params.maxEEFor
summary['maxEERev'] = params.maxEERev
summary['maxN'] = params.maxN
summary['maxLen'] = params.maxLen
summary['minLen'] = params.minLen
summary['rmPhiX'] = params.rmPhiX
summary['minOverlap'] = params.minOverlap
summary['maxMismatch'] = params.maxMismatch
summary['trimOverhang'] = params.trimOverhang
summary['species'] = params.species
summary['pool'] = params.pool
summary['Reference'] = params.reference
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 *
 * Step 1: Filter and trim (run per sample?)
 *
 */

process runFastQC {
    tag { "rFQC.${pairId}" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: false

    input:
        set pairId, file(in_fastq) from dada2ReadPairsToQual

    output:
        file("${pairId}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${pairId}_fastqc
    fastqc --outdir ${pairId}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

process runMultiQC{
    tag { "rMQC" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

process filterAndTrim {
    tag { "filterAndTrim" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: false
  
    input:
    set pairId, file(reads) from dada2ReadPairs

    output:
    set val(pairId), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" into filteredReadsforQC, filteredReads
    file "*.R1.filtered.fastq.gz" into forReads
    file "*.R2.filtered.fastq.gz" into revReads
    file "*.trimmed.txt" into trimTracking

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    
    #Variable selection from CLI input flag --rmPhix
    if("${params.rmPhiX}"=="F"){
      rm.phix <- FALSE
    } else if("${params.rmPhiX}"=="T"){
      rm.phix <- TRUE
    }
    print(rm.phix)
    out <- filterAndTrim(fwd = "${reads[0]}",
                        filt = paste0("${pairId}", ".R1.filtered.fastq.gz"),
                        rev = "${reads[1]}",
                        filt.rev = paste0("${pairId}", ".R2.filtered.fastq.gz"),
                        trimLeft = c(${params.trimFor},${params.trimRev}),
                        truncLen = c(${params.truncFor},${params.truncRev}),
                        maxEE = c(${params.maxEEFor},${params.maxEERev}),
                        truncQ = ${params.truncQ},
                        maxN = ${params.maxN},
                        rm.phix = rm.phix,
                        maxLen = ${params.maxLen},
                        minLen = ${params.minLen},
                        compress = TRUE,
                        verbose = TRUE,
                        multithread = ${task.cpus})
    write.csv(out, paste0("${pairId}", ".trimmed.txt"))
    """
}

process runFastQC_postfilterandtrim {
    tag { "rFQC_post_FT.${pairId}" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: "copy", overwrite: false

    input:
    set val(pairId), file(filtFor), file(filtRev) from filteredReadsforQC
    
    output:
        file("${pairId}_fastqc_postfiltertrim/*.zip") into fastqc_files_2

    """
    mkdir ${pairId}_fastqc_postfiltertrim
    fastqc --outdir ${pairId}_fastqc_postfiltertrim \
    ${filtFor} \
    ${filtRev}
    """
}

process runMultiQC_postfilterandtrim {
    tag { "rMQC_post_FT" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files_2.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

process mergeTrimmedTable {
    tag { "mergTrimmedTable" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: false
  
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
    tag { "LearnErrorsFor" }
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: false
  
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
    errF <- learnErrors(filtFs, multithread=${task.cpus})
    pdf("R1.err.pdf")
    plotErrors(errF, nominalQ=TRUE)
    dev.off()
    saveRDS(errF, "errorsF.RDS")
    """
}

process LearnErrorsRev {
    tag { "LearnErrorsRev" }
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: false
  
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
    errR <- learnErrors(filtRs, multithread=${task.cpus})
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
    tag { "SampleInferDerepAndMerge" }
    publishDir "${params.outdir}/dada2-Derep", mode: "copy", overwrite: false
  
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
    
    #Variable selection from CLI input flag --pool
   
    if("${params.pool}"=="pseudo"){
      pool <- "pseudo"
    } else if("${params.pool}"=="F"){
      pool <- FALSE
    } else if("${params.pool}"=="T"){
      pool <- TRUE 
    }
    print(pool)
    derepF <- derepFastq("${filtFor}")
    
    ddF <- dada(derepF, err=errF, multithread=${task.cpus}, pool=pool)

    derepR <- derepFastq("${filtRev}")
    ddR <- dada(derepR, err=errR, multithread=${task.cpus},pool=pool)

    #Variable selection from CLI input flag --trimOverhang
    if("${params.trimOverhang}"=="F"){
      trimOverhang <- FALSE
    } else if("${params.trimOverhang}"=="T"){
      trimOverhang <- TRUE
    }  
    print(trimOverhang)
    
    merger <- mergePairs(ddF, derepF, ddR, derepR,
        minOverlap = ${params.minOverlap},
        maxMismatch = ${params.maxMismatch},
        trimOverhang = trimOverhang
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
    tag { "mergeDadaRDS" }
    publishDir "${params.outdir}/dada2-Inference", mode: "copy", overwrite: false
  
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
    tag { "SequenceTable" }
    publishDir "${params.outdir}/dada2-SeqTable", mode: "copy", overwrite: false
  
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
        tag { "ChimeraTaxonomySpecies" }
        publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: false
      
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
        tag { "ChimeraTaxonomy" }
        publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: false
      
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
    tag { "AlignAndGenerateTree" }
    publishDir "${params.outdir}/dada2-Alignment", mode: "copy", overwrite: false
  
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
    tag { "BiomFile" }
    publishDir "${params.outdir}/dada2-BIOM", mode: "copy", overwrite: false
  
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
    tag { "ReadTracking" }
    publishDir "${params.outdir}/dada2-ReadTracking", mode: "copy", overwrite: false
  
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

/*
 * Completion e-mail notification
 */
workflow.onComplete {
  
    def subject = "[uct-cbio/16S-rDNA-dada2-pipeline] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[uct-cbio/16S-rDNA-dada2-pipeline] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[uct-cbio/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[uct-cbio/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (mail)"
        }
    }
}
