library(Biostrings)

viral_protein_sequences <- readAAStringSet(snakemake@input[['checkv']], format = "fasta")

# remove protein sequence increment from end of fasta name and save distinct names
viral_genome_names <- unique(gsub("_[^_]+$", "", viral_protein_sequences@ranges@NAMES))
for(genome in viral_genome_names){
  # subset to genomes with the same genome name
  genome_sequences <- viral_protein_sequences[grep(pattern = genome, x = viral_protein_sequences@ranges@NAMES)]
  # write protein sequences from the same genome to a file with that genomes name
  writeXStringSet(genome_sequences, 
                  filepath = paste0("outputs/checkv_db_v1.4_protein_sequences_by_genome", genome, ".faa.gz"), 
                  format = "fasta", 
                  compress = "gzip")
}

