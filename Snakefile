import pandas as pd

metadata = pd.read_csv("inputs/metadata.csv")
SAMPLES = metadata['run_accession'].unique().tolist()
KSIZES = [21, 31, 51]
LINEAGES=['bacteria', 'viral', 'archaea', 'fungi', 'protozoa']

rule all:
    input: expand("outputs/sourmash_taxonomy/{samples}-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = KSIZES, samples = SAMPLES)

rule k21:
    input: expand("outputs/sourmash_taxonomy/{samples}-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 21, samples = SAMPLES)

rule k31:
    input: 
        expand("outputs/sourmash_taxonomy/{samples}-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 31, samples = SAMPLES),
        expand("outputs/sourmash_compare/comp_k{ksize}.csv", ksize = 31),
        expand("outputs/sourmash_sketch_csv/{samples}_k{ksize}.csv", ksize = 31, samples = SAMPLES),
        expand("outputs/sourmash_sketch_csv_abund/{samples}_k{ksize}.csv", ksize = 31, samples = SAMPLES)

rule k51:
    input: expand("outputs/sourmash_taxonomy/{samples}-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 51, samples = SAMPLES)

##########################################################
## Download sourmash databases & taxonomy files
##########################################################

rule download_genbank_bacteria_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/6qxfp/download
    '''

rule download_genbank_bacteria_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/9ue5g/download
    '''

rule download_genbank_bacteria_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/5gvbw/download
    '''

rule download_genbank_bacteria_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/4agsp/download
    '''

rule download_genbank_fungi_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/fy82q/download
    '''
rule download_genbank_fungi_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/4pdbj/download
    '''
rule download_genbank_fungi_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/b9a86/download
    '''

rule download_genbank_fungi_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/s4b85/download
    '''

rule download_genbank_archaea_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k21.zip"
    conda: "envs/wget.yml"
    shell:''' 
    wget -O {output} https://osf.io/g94n5/download
    '''

rule download_genbank_archaea_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/hfybv/download
    '''

rule download_genbank_archaea_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/dehrc/download
    '''

rule download_genbank_archaea_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/kcbpn/download
    '''

rule download_genbank_viral_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/updvc/download
    '''

rule download_genbank_viral_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/egkt2/download
    '''

rule download_genbank_viral_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/z8scg/download
    '''

rule download_genbank_viral_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-viral.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/j4tsu/download
    '''

rule download_genbank_protozoa_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/m23r6/download 
    '''

rule download_genbank_protozoa_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k31.zip"
    conda: "envs/wget.yml"
    shell:''' 
    wget -O {output} https://osf.io/zm5vg/download
    '''

rule download_genbank_protozoa_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/32y98/download
    '''

rule download_genbank_protist_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/2x8u4/download
    '''

##########################################################
## Sketch files
##########################################################

rule sourmash_sketch_input_files:
    input: "inputs/raw/{samples}_flyeass.fasta"
    output: "outputs/sourmash_sketch/{samples}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.samples} -o {output} {input}
    '''

##########################################################
## Sketch compare
##########################################################

rule sourmash_compare:
    input: expand("outputs/sourmash_sketch/{samples}.sig", samples = SAMPLES)
    output: 
        comp = "outputs/sourmash_compare/comp_k{ksize}", 
        csv = "outputs/sourmash_compare/comp_k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare -k {wildcards.ksize} --csv {output.csv} -o {output.comp} {input}
    '''

##########################################################
## Default sourmash plots
##########################################################

rule sourmash_plot:
   input: "outputs/sourmash_compare/comp_k{ksize}" 
   output: "outputs/sourmash_plot/comp_k{ksize}.hist.pdf"
   params: outdir = "outputs/sourmash_plot"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash plot --pdf --labels --output-dir {params.outdir} 
   '''

##########################################################
## Sourmash gather
##########################################################

rule sourmash_gather:
    input:
        sig="outputs/sourmash_sketch/{samples}.sig",
        databases=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{{ksize}}.zip", lineage = LINEAGES)
    output: csv="outputs/sourmash_gather/{samples}-vs-genbank-2022.03-k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} --scaled 1000 --threshold-bp 0 -o {output.csv} {input.sig} {input.databases}
    '''
   
##########################################################
## Sourmash taxonomy
##########################################################

rule gunzip_lineage_csvs:
    input: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv.gz"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule sourmash_taxonomy_prepare:
    input: expand("inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv", lineage = LINEAGES),
    output: "outputs/sourmash_taxonomy/genbank-2022.03-prepared-lineages.sqldb"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax prepare --taxonomy-csv {input} -o {output}
    '''

rule sourmash_taxonomy_annotate:
   input:
       lin_prepared="outputs/sourmash_taxonomy/genbank-2022.03-prepared-lineages.sqldb",
       gather="outputs/sourmash_gather/{samples}-vs-genbank-2022.03-k{ksize}.csv"
   output: "outputs/sourmash_taxonomy/{samples}-vs-genbank-2022.03-k{ksize}.with-lineages.csv"
   params: outdir = "outputs/sourmash_taxonomy/"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash tax annotate -g {input.gather} -t {input.lin_prepared} -o {params.outdir}
   '''

# visualizations:
# upset plot on hashes
# upset plot on taxonomy?
# taxonomy bar charts as fraction of sample
# ...

####################################################################
## Visualization: upset plot for the intersection of shared 
##                hashes (subsampled k-mers) between samples
####################################################################

rule sourmash_sketch_convert_to_csv:
    input: "outputs/sourmash_sketch/{samples}.sig"
    output: "outputs/sourmash_sketch_csv/{samples}_k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    python scripts/sig_to_csv.py {wildcards.ksize} {input} {output}
    '''

# TODO: update download link to main when script is merged into mtx repo
rule download_sig_to_csv_abund_script:
    output: "scripts/sig_to_csv_abund.py"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://raw.githubusercontent.com/Arcadia-Science/2022-mtx-not-in-mgx-pairs/ter/specaccum/scripts/sig_to_csv_abund.py
    '''

rule sourmash_sketch_convert_to_csv_abund:
    input:
        py = "scripts/sig_to_csv_abund.py", 
        sig="outputs/sourmash_sketch/{samples}.sig"
    output: "outputs/sourmash_sketch_csv_abund/{samples}_k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    python scripts/sig_to_csv_abund.py {wildcards.ksize} {input.sig} {output}
    '''
