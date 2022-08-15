# Sourmash analysis of long read assemblies from a subset of run accessions in PRJNA853785

PRJNA853785 used nanopore sequencing to profile time series cheese production.

## Getting started with this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```
mamba env create -n prjna853785 --file environment.yml
conda activate prjna853785
```

To start the pipeline, run:
```
snakemake k31 --use-conda -j 1
```

### Running this repository on AWS

While this pipeline can be run on any computational infrastructure, it was designed to be run on AWS.
These results can be replicated by following the instructions and commands below.

Launch an AWS EC2 instance (i3.4xlarge machine running AMI `ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-20220609` with 100GB of hard drive space).
100GB is sufficient to run two k-mer size end points (e.g. `k21`, `k31`), but ~150GB is needed to run all three end points.
After launching your EC2 instance, run the following commands to install and configure miniconda.

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda install mamba # install mamba for faster software installation.
```

Once minicond is configured, clone the repository, `cd` into the repository, and follow the "Getting started" instructions above.
```
git clone https://github.com/Arcadia-Science/2022-prjna853785-sourmash.git
cd 2022-prjna853785-sourmash
```


