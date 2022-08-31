import pandas as pd
import os

metadata = pd.read_csv("inputs/metadata.csv")
SAMPLES = metadata['run_accession'].unique().tolist()
PROTEIN_KSIZES = [7, 10]
DNA_KSIZES = [21, 31, 51]


rule protein:
    input: expand("outputs/sourmash_gather/{sample}ass-vs-checkv-db-v1.4-protein-k{protein_ksize}.csv", sample = SAMPLES, protein_ksize = PROTEIN_KSIZES)

rule dna:
    input: expand("outputs/sourmash_gather/{sample}ass-vs-checkv-db-v1.4-dna-k{dna_ksize}.csv", sample = SAMPLES, dna_ksize = DNA_KSIZES)

##############################################################
## Download CheckV database
##############################################################

rule download_checkv_database:
    output: "inputs/checkv-database/checkv-db-v1.4.tar.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://portal.nersc.gov/CheckV/checkv-db-v1.4.tar.gz
    '''

rule decompress_checkv_database:
    input: "inputs/checkv-database/checkv-db-v1.4.tar.gz"
    output: 
        "inputs/checkv-database/checkv-db-v1.4/genome_db/checkv_reps.fna",
        "inputs/checkv-database/checkv-db-v1.4/genome_db/checkv_reps.faa"
    params: outdir = "inputs/checkv-database"
    shell:'''
    tar xf {input} -C {params.outdir} 
    '''

###############################################################
## Create sourmash database for CheckV nucleotide sequences
###############################################################

rule sourmash_sketch_checkv_fna:
    """
    Sketch the checkV nucleotide database file.
    Each viral genome is a single contig in the file, so the --singleton parameter will produce a sketch for each of these contigs.
    The --name-from-first command makes sure every sketch has an appropriate name pulled from the fasta name.
    This results in a single file that contains all of the sketches for every contig and every parameter.
    """
    input: "inputs/checkv-database/checkv-db-v1.4/genome_db/checkv_reps.fna"
    output: "outputs/checkv-db-v1.4-sourmash-sketch/checkv-db-v1.4-reps-dna.sig"
    benchmark: "benchmarks/sourmash-sketch-dna-checkv-db-v1.4.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=100,abund --name-from-first --singleton -o {output} {input}
    '''

rule sourmash_sig_cat_checkv_fna:
    """
    Select a specific k-mer size and create a zip file that will be used as a gather database.
    """
    input: "outputs/checkv-db-v1.4-sourmash-sketch/checkv-db-v1.4-reps-dna.sig"
    output: "outputs/sourmash-databases/checkv-db-v1.4-reps-dna-k{dna_ksize}.zip"
    benchmark: "benchmarks/sourmash-sig-cat-checkv-db-v1.4-dna-k{dna_ksize}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig cat -k {wildcards.dna_ksize} -o {output} {input}
    '''

################################################################
## Create sourmash database for CheckV protein sequences
################################################################

checkpoint separate_checkv_protein_sequences_into_genomes:
    input: checkv="inputs/checkv-database/checkv-db-v1.4/genome_db/checkv_reps.faa"
    output: directory("outputs/checkv-db-v1.4-protein-sequences-by-genome/")
    benchmark: "benchmarks/separate-checkv-db-v1.4-protein-sequences-into-genomes.txt"
    conda: "envs/biostrings.yml"
    script: "scripts/separate-checkv-protein-sequences-into-genomes.R"

rule sourmash_sketch_checkv_faa:
    input: "outputs/checkv-db-v1.4-protein-sequences-by-genome/{viral_genome}.faa.gz"
    output: "outputs/checkv-db-v1.4-sourmash-sketch/{viral_genome}-protein.sig"
    benchmark:"benchmarks/sourmash-sketch-protein-checkv-db-v1.4/{viral_genome}-protein.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch protein -p k=7,k=10,scaled=200,abund --name {wildcards.viral_genome} -o {output} {input}
    '''

def checkpoint_separate_checkv_protein_sequences_into_genomes(wildcards):
    checkpoint_output = checkpoints.separate_checkv_protein_sequences_into_genomes.get(**wildcards).output[0]    
    file_names = expand("outputs/checkv-db-v1.4-sourmash-sketch/{viral_genome}-protein.sig", 
                        viral_genome = glob_wildcards(os.path.join(checkpoint_output, "{viral_genome}.faa.gz")).viral_genome)
    return file_names


rule sourmash_sig_cat_checkv_faa:
    input: checkpoint_separate_checkv_protein_sequences_into_genomes
    output: "outputs/sourmash-databases/checkv-db-v1.4-reps-protein-k{protein_ksize}.zip"
    benchmark: "benchmarks/sourmash-sig-cat-checkv-db-v1.4-protein-k{protein_ksize}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig cat -k {wildcards.protein_ksize} {input}
    '''


###################################################################
## Sourmash sketch input sequences in nucleotide and protein space
###################################################################

rule sourmash_sketch_input_files_dna_for_viral:
    input: "inputs/raw/{sample}ass.fasta"
    output: "outputs/sourmash-sketch/{sample}ass-dna-scaled100.sig"
    benchmark: "benchmarks/sourmash-sketch/{sample}-dna-scaled100.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=100,abund --name {wildcards.sample} -o {output} {input}
    '''

rule sourmash_sketch_input_files_protein_for_viral:
    input: "inputs/raw/{sample}ass.fasta"
    output: "outputs/sourmash-sketch/{sample}ass-protein-scaled200.sig"
    benchmark: "benchmarks/sourmash-sketch/{sample}-protein-scaled200.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch protein -p k=7,k=10,scaled=200,abund --name {wildcards.sample} -o {output} {input}
    '''

###################################################################
## Sourmash gather against checkV viral databases
###################################################################


rule sourmash_gather_dna:
    input:
        sig="outputs/sourmash-sketch/{sample}ass-dna-scaled100.sig",
        databases="outputs/sourmash-databases/checkv-db-v1.4-reps-dna-k{dna_ksize}.zip"
    output: csv="outputs/sourmash_gather/{sample}ass-vs-checkv-db-v1.4-dna-k{dna_ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.dna_ksize} --scaled 100 --threshold-bp 0 -o {output.csv} {input.sig} {input.databases}
    '''

rule sourmash_gather_protein:
    input:
        sig="outputs/sourmash-sketch/{sample}ass-protein-scaled200.sig",
        databases="outputs/sourmash-databases/checkv-db-v1.4-reps-protein-k{protein_ksize}.zip"
    output: csv="outputs/sourmash_gather/{sample}ass-vs-checkv-db-v1.4-protein-k{protein_ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.protein_ksize} --scaled 200 --threshold-bp 0 -o {output.csv} {input.sig} {input.databases}
    '''

