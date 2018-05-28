# ![kviljoen/16S-rDNA-dada2-pipeline](/assets/cbio_logo.png)

# dada2 Nextflow workflow

(BETA)

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow for implementation on the UCT high-performance compute cluster

Usage:
    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/16S-rDNA-dada2-pipeline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. local / uct_hex
      --trimFor                     Set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)
      --trimRev                     Set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz)
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.



