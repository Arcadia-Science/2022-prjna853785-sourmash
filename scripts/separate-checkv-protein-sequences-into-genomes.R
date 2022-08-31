library(Biostrings)

viral_protein_sequences <- readAAStringSet(snakemake@input[['checkv']], format = "fasta")

# make the output directory for the protein sequences
dir.create(file.path("outputs/checkv-db-v1.4-protein-sequences-by-genome/"), showWarnings = FALSE)
# remove protein sequence increment from end of fasta name and save distinct names
viral_genome_names <- unique(gsub("_[^_]+$", "", viral_protein_sequences@ranges@NAMES))
for(genome in viral_genome_names){
  # subset to genomes with the same genome name
  genome_sequences <- viral_protein_sequences[grep(pattern = genome, x = viral_protein_sequences@ranges@NAMES)]
  # write protein sequences from the same genome to a file with that genomes name
  filepath <- paste0("outputs/checkv-db-v1.4-protein-sequences-by-genome/", genome, ".faa.gz")
  writeXStringSet(genome_sequences, 
                  filepath = filepath,
                  format = "fasta", 
                  compress = "gzip")
}
