# ![kviljoen/16S-rDNA-dada2-pipeline](/assets/cbio_logo.png)

# 16S rRNA gene amplicon sequencing pipeline using DADA2, implemented in Nextflow

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow for implementation on the UCT high-performance compute cluster

## Basic usage:

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
      Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline --help

## Prerequisites

Nextflow, dada2 (>= 1.8), R (>= 3.2.0), Rcpp (>= 0.11.2), methods (>= 3.2.0), DECIPHER, phangorn, biomformat

## Documentation
The uct-cbio/16S-rDNA-dada2-pipeline pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)

## Built With

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://singularity.lbl.gov/)


## Credits

The initial implementation of the DADA2 pipeline as a Nextflow workflow (https://github.com/HPCBio/dada2-Nextflow) was done by members of the high performance computational biology unit at the University of Illinois (http://www.hpcbio.illinois.edu). Please remember to cite the authors of [DADA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) when using this pipeline. Further development to the Nextflow workflow and containerisation in Docker and Singularity for implementation on UCT's HPC was done by Dr Katie Lennard and Gerrit Botha, with inspiration and code snippets from Phil Ewels http://nf-co.re/

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


