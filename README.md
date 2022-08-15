# Sourmash analysis of long read assemblies from a subset of run accessions in PRJNA853785

PRJNA853785 used nanopore sequencing to profile time series cheese production.
The data analyzed here are a subset of Nanopore raw reads assembled with Flye 2.9-b1768 using the metag flag.

## Purpose of this repository

This repository implements snakemake pipeline to run sourmash commands on metagenome assemblies and visualize their outputs.
While the results of this pipeline were needed for a biology project, this pipeline alpha tests the generalization of a sourmash commands and visualizations with the idea that these could be agnostically applied to any new metagenome assembly.
In keeping with this goal the visualizations incorporate minimal metadata and use only the sample names, as these would be available in all use cases.
This strategy my change as this repository matures.
This pipeline may also eventually be translated from snakemake into Nextflow.

## Background information on sourmash

[Sourmash](sourmash.readthedocs.io) is a tool that facilitates lightweight and rapid comparisons between potentially very large sets of sequencing data (see [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/) and [here](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract)).
Sourmash uses FracMinHash sketching to generate compressed representations of sequencing data.
Two core parameters determine which sequences are represented in these sketches: k-mers size and scaled value.
The k-mer size dictates the length of sequences represented in the sketch; defaults are *k* = 21, *k* = 31, and *k* = 51.
While these sizes are somewhat empiric, overlap between related genomic sequences for a given k-mer size roughly approximates taxonomic relationships, with substantial overlap at *k* = 21 roughly approximating genus-level relatedness, *k* = 31 species-level, and *k* = 51 strain-level.
The scaled value dictates the fraction of k-mers in the original sequence that get included in the sketch.
Approximately 1/scaled value of k-mers in the original sequence get included in the sketch.
The sourmash defaults are 1000, 2000, and 10,000.
The scaled value takes advantage of the fact that a subset of data can be used to accurately estimate things like similarity and containment and dramatically sub-samples the original sequences thus facilitating rapid comparisons even between very large data sets.

This repository uses the sourmash commands `sketch`, `compare`, `gather`, and `taxonomy`.
+ `sketch`: Produces sketches for each metagenome assembly using parameters *k* = 21, *k* = 31, and *k* = 51, scaled = 1000, and abundance tracking of k-mers turned on.
+ `compare`: Compares sketches of each metagenome to produce cosine distance estimates (defaults to including abundance information).
+ `plot`: Plots the output of the compare command.
+ `gather`: Compares each metagenome assembly against all microbial (bacterial, archaeal, viral, fungal, protozoa) in GenBank (March 2022) to identify the minimum set of genomes that contain all of the known k-mers in the metagenome.
+ `taxonomy`: Assigns taxonomic lineages to the sourmash gather matches.

To aid in the interpretation of these commands, this repository also has notebooks that visualize the outputs of sourmash.
TBD.

## Navigating this repository

TBD

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

The pipeline has four encoded endpoints: `k21`, `k31`, `k51`, and `all`. 
These endpoints are used to execute a subset of rules associated with a single k-mer size at a given time.
If the pipeline is called without any endpoints specified, `all` endpoints will be run by default.

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

Once miniconda is configured, clone the repository, `cd` into the repository, and follow the "Getting started" instructions above.
```
git clone https://github.com/Arcadia-Science/2022-prjna853785-sourmash.git
cd 2022-prjna853785-sourmash
```

