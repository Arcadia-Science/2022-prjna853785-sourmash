import pandas as pd

metadata = pd.read_csv("inputs/metadata.csv")
SAMPLES = metadata['run_accession'].unique().tolist()
KSIZES = [21, 31, 51]
LINEAGES=['bacteria', 'viral', 'archaea', 'fungi', 'protozoa']

rule all:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = KSIZES, samples = SAMPLES)

rule k21:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 21, samples = SAMPLES)

rule k31:
    input: 
        expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 31, samples = SAMPLES),
        expand("outputs/sourmash_compare/comp_k{ksize}.csv", ksize = 31),
        expand("outputs/sourmash_sketch_csv/{samples}_k{ksize}.csv", ksize = 31, samples = SAMPLES),
        expand("outputs/sourmash_sketch_csv_abund/{samples}_k{ksize}.csv", ksize = 31, samples = SAMPLES)

rule k51:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 51, samples = SAMPLES)

##########################################################
## Download sourmash databases & taxonomy files
##########################################################

rule download_sourmash_databases_genbank:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}-k{ksize}.zip"
    run:
        sourmash_database_info = pd.read_csv(input)
        lineage_df = sourmash_database_info.loc[sourmash_database_info['lineage'] == wildcards.lineage]
        db_df = lineage_df.loc[lineage_df['ksize'] == wildcards.ksize]
        osf_hash = db_df['osf_hash'].values[0] 
        shell("wget -O {output} https://osf.io/{osf_hash}/download")

rule download_sourmash_lineages_genbank:
    input: "inputs/sourmash_databases/sourmash-lineage-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv.gz"
    run:
        sourmash_lineage_info = pd.read_csv(input)
        lineage_df = sourmash_lineage_info.loc[sourmash_lineage_info['lineage'] == wildcards.lineage]
        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("wget -O {output} https://osf.io/{osf_hash}/download")

##########################################################
## Sketch files
##########################################################

rule sourmash_sketch_input_files:
    input: "inputs/raw/{samples}ass.fasta"
    output: "outputs/sourmash_sketch/{samples}ass.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.samples} -o {output} {input}
    '''

##########################################################
## Sketch compare
##########################################################

rule sourmash_compare:
    input: expand("outputs/sourmash_sketch/{samples}ass.sig", samples = SAMPLES)
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
        sig="outputs/sourmash_sketch/{samples}ass.sig",
        databases=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{{ksize}}.zip", lineage = LINEAGES)
    output: csv="outputs/sourmash_gather/{samples}ass-vs-genbank-2022.03-k{ksize}.csv"
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
       gather="outputs/sourmash_gather/{samples}ass-vs-genbank-2022.03-k{ksize}.csv"
   output: "outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv"
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
    input: "outputs/sourmash_sketch/{samples}ass.sig"
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
        sig="outputs/sourmash_sketch/{samples}ass.sig"
    output: "outputs/sourmash_sketch_csv_abund/{samples}_k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    python scripts/sig_to_csv_abund.py {wildcards.ksize} {input.sig} {output}
    '''
