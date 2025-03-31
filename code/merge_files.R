library(dplyr)
library(tidyr)

NapA_Alpha <- read.csv("SI_TS_NapA_alpha_diversiy.csv")

## Abundance

merge_depths <- function(file_list) {
  results <- data.frame(stringsAsFactors = FALSE)
  
  for (file in file_list) {
    # Read the TSV file
    data <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Extract depth from filename (remove 'm' before converting to numeric)
    depth <- as.numeric(sub("m.*", "", sub("_.*", "", basename(file))))
    
    # Filter for NapA marker and add depth column
    napa_data <- subset(data, Marker == "NapA")
    napa_data$Depth_m <- depth
    
    # Append to results
    results <- rbind(results, napa_data)
  }
  
  return(results)
}

extract_abundance_by_depth <- function(file_list) {
  results <- data.frame(Depth_m = numeric(), Abundance = numeric(), stringsAsFactors = FALSE)
  
  for (file in file_list) {
    # Extract depth from filename (remove 'm' before converting to numeric)
    depth <- as.numeric(sub("m.*", "", sub("_.*", "", basename(file))))
    
    # Get NapA abundance
    abundance <- extract_napa_abundance(file)
    
    # Append to results
    results <- rbind(results, data.frame(Depth_m = depth, Abundance = abundance))
  }
  
  return(results)
}

DNA_file_list = c("dna/10m_classifications.tsv", "dna/100m_classifications.tsv", "dna/120m_classifications.tsv", "dna/135m_classifications.tsv", "dna/150m_classifications.tsv", "dna/165m_classifications.tsv", "dna/200m_classifications.tsv")

RNA_file_list = c("rna/10m_classifications_rna.tsv", "rna/100m_classifications_rna.tsv", "rna/120m_classifications_rna.tsv", "rna/135m_classifications_rna.tsv", "rna/150m_classifications_rna.tsv", "rna/165m_classifications_rna.tsv", "rna/200m_classifications_rna.tsv")


#Save merged data as file

  write.table(merge_depths(DNA_file_list), file = "merged_depths_dna.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(merge_depths(RNA_file_list), file = "merged_depths_rna.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
