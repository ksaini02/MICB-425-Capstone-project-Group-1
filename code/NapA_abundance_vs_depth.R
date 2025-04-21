library(ggplot2)
library(tidyverse)
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

depths_sum_dna_df <- dna_merged_data %>%
    group_by(Depth_m) %>%
    summarize(total_abundance = sum(Abundance))

#Plot abundance by depth  
dna_abundance_vs_depth = ggplot(dna_merged_data, aes(x = as.factor(Depth_m), y = Abundance, fill = n)) +
    geom_violin(alpha = 0.8) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Make plots more narrow
    scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient based on density
    scale_y_log10() +  
    labs(x = "Depth (m)", y = expression(italic("NapA")~" Abundance (Log10)"), 
         title = "", fill = "Classified Sequences") +
    theme_test() +
    geom_signif(comparisons = depth_comparisons, 
                test = "wilcox.test", 
                map_signif_level = TRUE) +
    geom_text(data = depths_sum_dna_df, 
              aes(x = as.factor(Depth_m), y = 0.08,  # adjust y position as needed
                  label = round(total_abundance, 0)),
              inherit.aes = FALSE,
              size = 5, color = "black")
a
##RNA

rna_merged_data <- rna_merged_data %>%
    filter(Abundance > 0)  # Remove zero values


# Count number of data points per depth
depth_counts_rna <- rna_merged_data %>%
    group_by(Depth_m) %>%
    summarize(n = n())

rna_merged_data <- rna_merged_data %>%
    left_join(depth_counts_rna, by = "Depth_m")

depths_sum_df <- rna_merged_data %>%
    group_by(Depth_m) %>%
    summarize(total_abundance = sum(Abundance))

#Plot abundance by depth    
rna_abundance_vs_depth = ggplot(rna_merged_data, aes(x = as.factor(Depth_m), y = Abundance, fill = n)) +
    geom_violin(alpha = 0.8) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Make plots more narrow
    scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient based on density
    scale_y_log10() +  
    labs(x = "Depth (m)", y = expression(italic("NapA")~"Transcript Abundance (Log10)"), 
         title = "", fill = "Classified Sequences") +
    theme_test() +
    geom_signif(comparisons = depth_comparisons, 
                test = "wilcox.test", 
                map_signif_level = TRUE) +
    geom_text(data = depths_sum_df, 
              aes(x = as.factor(Depth_m), y = 1,  # adjust y position as needed
                  label = round(total_abundance, 0)),
              inherit.aes = FALSE,
              size = 5, color = "black")

ggsave(filename = "plots/rna_abundance_vs_depth.png", plot = rna_abundance_vs_depth)

ggsave(filename = "plots/dna_abundance_vs_depth.png", plot = dna_abundance_vs_depth)
