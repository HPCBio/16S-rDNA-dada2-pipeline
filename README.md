# dada2 workflow

(BETA)

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. 

# Caveats

* This currently is hard-coded to utilize the Biocluster module system, in particular a single module containing all the necessary R libraries:
  * `dada2`
  * `DECIPHER`
  * `phangorn`


