library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

dna_merged_data = read.delim("merged_depths_dna.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rna_merged_data = read.delim("merged_depths_rna.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define comparisons
depth_comparisons <- list(
  c("10", "100"), 
  c("100", "120"), 
  c("120", "135"),
  c("135", "150"),
  c("150", "165"),
  c("165", "200")
)

## DNA

# Count number of data points per depth
depth_counts_dna <- dna_merged_data %>%
    filter(Abundance > 0) %>%  # Remove zero values
    group_by(Depth_m) %>%
    summarize(n = n())

# Merge the counts back into the main dataset
dna_merged_data <- dna_merged_data %>%
    left_join(depth_counts_dna, by = "Depth_m")

#Plot abundance by depth  
ggplot(dna_merged_data, aes(x = as.factor(Depth_m), y = Abundance, fill = n)) +
    geom_violin(alpha = 0.8) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Make plots more narrow
    scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient based on density
    scale_y_log10() +  
    labs(x = "Depth (m)", y = "NapA Abundance (Log10)", 
         title = "NapA Abundance by Depth", fill = "Observed Samples") +
    theme_test() +
    geom_signif(comparisons = depth_comparisons, 
                test = "wilcox.test", 
                map_signif_level = TRUE)

##RNA

# Count number of data points per depth
depth_counts_rna <- rna_merged_data %>%
    filter(Abundance > 0) %>%  # Remove zero values
    group_by(Depth_m) %>%
    summarize(n = n())

rna_merged_data <- rna_merged_data %>%
    filter(Abundance > 0) %>%  # Remove zero values
    left_join(depth_counts_rna, by = "Depth_m")

#Plot abundance by depth    
ggplot(rna_merged_data, aes(x = as.factor(Depth_m), y = Abundance, fill = n)) +
    geom_violin(alpha = 0.8) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Make plots more narrow
    scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient based on density
    scale_y_log10() +  
    labs(x = "Depth (m)", y = "NapA Abundance (Log10)", 
         title = "NapA Abundance by Depth", fill = "Observed Samples") +
    theme_test() +
    geom_signif(comparisons = depth_comparisons, 
                test = "wilcox.test", 
                map_signif_level = TRUE)

