# this snakefile runs checkV on the PRJNA853785 assemblies to use results as a pseudo-validation set for the sourmash viral results.

import pandas as pd
import os

metadata = pd.read_csv("inputs/metadata.csv")
SAMPLES = metadata['run_accession'].unique().tolist()


rule all:
    input: expand("outputs/checkv/{sample}/completeness.tsv", sample = SAMPLES)
        
rule checkv_download_and_build_database:
    output: "inputs/checkv-database-built/checkv-db-v1.4/genome_db/checkv_reps.dmnd"
    params: outdir = "inputs/checkv-database-built"
    conda: "envs/checkv.yml"
    shell:'''
    checkv download_database {params.outdir}
    '''

rule checkv:
    input:
        db = "inputs/checkv-database-built/checkv-db-v1.4/genome_db/checkv_reps.dmnd",
        fasta="inputs/raw/{sample}ass.fasta"
    output: "outputs/checkv/{sample}/completeness.tsv"
    params:
        outdir = lambda wildcards: "outputs/checkv/" + wildcards.sample,
        dbdir  = "inputs/checkv-database-built/checkv-db-v1.4"
    threads: 16
    conda: "envs/checkv.yml"
    shell:'''
    checkv end_to_end {input.fasta} {params.outdir} -t {threads} -d {params.dbdir}
    '''

# calculate num nucleotides/percent of sample confidently identified as viral
# try and figure out what the confident calls match (I think best match is in completeness.tsv)
