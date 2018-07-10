# uct-cbio/16S-rDNA-dada2-pipeline Installation

To start using the uct-cbio/16S-rDNA-dada2-pipeline, follow the steps below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
sudo mv nextflow /usr/local/bin
```

###For Univeristy of Cape Town users working on HPC (hex):
```
#From your home directory on hex install nextflow
curl -fsSL get.nextflow.io | bash

#Add the following to ~/.bashrc:
export JAVA_HOME=/opt/exp_soft/java/jdk1.8.0_31/

#Do not run nextflow from the headnode, it requires substantial memory to run java. Please therefore first start an interactive job as follows: 
qsub -I -q UCTlong -l nodes=1:series600:ppn=1 -d `pwd`
```

**You need NextFlow version >= 0.24 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `uct-cbio/16S-rDNA-dada2-pipeline` is specified as the pipeline name when executing `nextflow run uct-cbio/16S-rDNA-dada2-pipeline`. If for some reason you need to use the development branch, this can be specified as `nextflow run uct-cbio/16S-rDNA-dada2-pipeline -r dev`

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/uct-cbio/16S-rDNA-dada2-pipeline/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/16S-rDNA-dada2-pipeline-master
```

---

[![UCT Computational Biology](/assets/cbio_logo.png)](http://www.cbio.uct.ac.za/)

---
