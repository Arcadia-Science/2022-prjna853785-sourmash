import pandas as pd

metadata = pd.read_csv("inputs/metadata.csv")
SAMPLES = metadata['run_accession'].unique().tolist()
KSIZES = [21, 31, 51]

rule all:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = KSIZES, samples = SAMPLES)

rule k21:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 21, samples = SAMPLES)

rule k31:
    input: 
        expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 31, samples = SAMPLES),
        expand("outputs/sourmash_compare/comp_k{ksize}.csv", ksize = 31)

rule k51:
    input: expand("outputs/sourmash_taxonomy/{samples}ass-vs-genbank-2022.03-k{ksize}.with-lineages.csv", ksize = 51, samples = SAMPLES)

##########################################################
## Download sourmash databases & taxonomy files
##########################################################

rule download_genbank_bacteria_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeif2hdztfrevkngnfqk3bsoyajxxf67o57u4dezbz647jwcf6gnwoy
    '''

rule download_genbank_bacteria_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeigkcvizvhe3xzxsuzv3ryf3ogvgvcmms2e5nfk7epl5egts22jyue
    '''

rule download_genbank_bacteria_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeie3eyyectnh5xqxz44oa3qj5vura3bffqdwfqk6jjuzzadkh7e2sq
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
    wget -O {output} https://dweb.link/ipfs/bafybeibrirvek4lxn6hh3wgsmtsd5vz5gtmewpzeg364bix3hojghwmygq
    '''
rule download_genbank_fungi_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeidhhwvwujkteno5ugwgjy4brhrv5dff2aumifcuew73qolfktdndq
    '''
rule download_genbank_fungi_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeibnrtt45f7wez2xb3fy5rxhatpeevc3rilm2gs65u5h6gc4u72fam
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
    wget -O {output} https://dweb.link/ipfs/bafybeiepywe7c6zjzgh3rksqiwpo5zyb7uuefbgvbn5nkgiq77iaavpzl4
    '''

rule download_genbank_archaea_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeidn6epju7yrdxrktq5wjko2yiwp6nrx3mq37htiuwecm7lffrbcdi
    '''

rule download_genbank_archaea_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeifyrwbx5dnay4mflboc5zai2de3xrvcxtgiu4j7adzj6qrxhb3zva
    '''

rule download_genbank_archaea_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/kcbpn/download
    '''

rule download_genbank_virus_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeicjyx6qkhdtw6q4cxs6fyl46gqfhd4q5eqje5lkswf2npljnyytzi
    '''

rule download_genbank_virus_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeibqsldwsztjf66rwvwnb6hamjtsfkmdk5bmfqbzwrod6wwwkqz2ya
    '''

rule download_genbank_virus_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeibgifuv4q3mihfubnhjhwm2esjnoseudpnzwahlkp3hvlbmtd4s2q
    '''

rule download_genbank_virus_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-viral.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/j4tsu/download
    '''

rule download_genbank_protozoa_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeicfh4xl4wuxd4xy2tf73hfamxlqa3s2higghnsjay4t5wtlmsdo5y
    '''

rule download_genbank_protozoa_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeicpxjhfrzem7f34eghbbwm3vglz2njxo72vpqcw7foilfomexsghi
    '''

rule download_genbank_protozoa_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://dweb.link/ipfs/bafybeigfpxkmzyq6sdkob53l6ztiy5ro44dzkad7dxakhuaao6cw4gp4eu
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
        db1="inputs/sourmash_databases/genbank-2022.03-bacteria-k{ksize}.zip",
        db2="inputs/sourmash_databases/genbank-2022.03-viral-k{ksize}.zip",
        db3="inputs/sourmash_databases/genbank-2022.03-protozoa-k{ksize}.zip",
        db4="inputs/sourmash_databases/genbank-2022.03-archaea-k{ksize}.zip",
        db5="inputs/sourmash_databases/genbank-2022.03-fungi-k{ksize}.zip" 
    output: csv="outputs/sourmash_gather/{samples}ass-vs-genbank-2022.03-k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} --scaled 1000 --threshold-bp 0 -o {output.csv} {input.sig} {input.db1} {input.db2} {input.db3} {input.db4} {input.db5}
    '''
   
##########################################################
## Sourmash taxonomy
##########################################################

rule sourmash_taxonomy_prepare:
    input:
        lin1="inputs/sourmash_databases/genbank-2022.03-bacteria.lineages.csv.gz",
        lin2="inputs/sourmash_databases/genbank-2022.03-viral.lineages.csv.gz",
        lin3="inputs/sourmash_databases/genbank-2022.03-protozoa.lineages.csv.gz", 
        lin4="inputs/sourmash_databases/genbank-2022.03-archaea.lineages.csv.gz",
        lin5="inputs/sourmash_databases/genbank-2022.03-fungi.lineages.csv.gz",
    output: "outputs/sourmash_taxonomy/genbank-2022.03-prepared-lineages.sqldb"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax prepare --taxonomy-csv {input.lin1} {input.lin2} {input.lin3} {input.lin4} {input.lin5} -o {output}
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
