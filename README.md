# dada2 Nextflow workflow

(BETA)

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow for implementation on the UCT high-performance compute cluster

# Caveats

* This currently is hard-coded to utilize the Biocluster module system, in particular a single module containing all the necessary R libraries:
  * `dada2`
  * `DECIPHER`
  * `phangorn`


