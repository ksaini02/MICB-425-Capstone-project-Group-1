library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(ggsignif)

### Question 1
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



### Question 2

##Figure 2: Plotting Depth against Nutrient Concentrations 

#Filter saanich dataset
saanich_data <- read.csv("Saanich_data.csv") |>
  filter(Cruise== 72) |>
  filter(Depth == 0.010 | Depth == 0.100 | Depth == 0.12 | Depth == .135 | Depth==.150 | Depth==.165 | Depth==.200)|> 
  mutate("Depth" =Depth*1000)

merged_data<- cbind(saanich_data, alpha_div)
long_data <- merged_data %>%
  pivot_longer(cols = c("WS_NO3", "Mean_NO2", "Mean_NH4", "WS_O2", "WS_H2S"), 
               names_to = "Nutrient", 
               values_to = "Concentration")
# Create the plot with facet wrap for each variable
plot<- ggplot(long_data, aes(x = Concentration, y = Depth, color= Nutrient)) +
  geom_point(size=3)+  # Increase point size for clarity
  scale_y_reverse() +  
  labs(title = "Depth vs Nutrient Concentrations", 
       x = "Concentration (uM)",
       y = "Depth (m)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(. ~Nutrient, scales= "free_x")


ggsave("Depth_v_Concentration.png", plot, height = 6, width = 8, units = "in", dpi = 300)

##Figure 3: Alpha diversity metrics of NapA DNA in Saanich Inlet.
# plotting DNA alpha diversity against depth
alpha_napA <- read.csv("capstone/alpha_diversity/SI_TS_NapA_alpha_diversiy.csv")
alpha_napA$placerun <- as.numeric(gsub("SI072_(\\d+)m.*", "\\1", alpha_napA$placerun))
alpha_napA_longer <- alpha_napA |>
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")
order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")
alpha_napA_longer$metric <- factor(alpha_napA_longer$metric, levels = order)
ggplot(alpha_napA_longer, aes(y=placerun, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Depth",
       x = "Value",
       y = "Depth") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x")

# reading in and filtering Saanich data to include desired columns
saanich <- read.csv("Saanich_Data.csv") |>
  mutate(depth = Depth*1000) |>
  select(depth, WS_O2, WS_NO3, WS_H2S, Mean_CH4)

#calculating mean of each concentration at each depth
saanich_mod <- aggregate(cbind(WS_O2, WS_NO3, WS_H2S, Mean_CH4) ~ depth, data = saanich, FUN = mean)

# combining DNA alpha diversity with Saanich data
DNA_combined <- merge(alpha_napA, saanich_mod, by.x = "placerun", by.y = "depth")

#plotting DNA alpha diversity (rooted PD) against oxygen, sulfide, and nitrate concentrations
DNA_combined_longer <- DNA_combined |>
  pivot_longer(cols = c("WS_O2", "WS_NO3", "WS_H2S", "placerun", "Mean_CH4"),
               names_to = "measure",
               values_to = "value")
combined_order <- c("depth", "WS_O2", "WS_NO3", "WS_H2S", "Mean_CH4")
DNA_combined_longer$measure <- factor(DNA_combined_longer$measure, levels = combined_order)
rootedpd <- ggplot(DNA_combined_longer, aes(y=value, x=rooted_pd)) +
  geom_point(aes(colour=measure, shape=measure), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity",
       x = "Rooted PD",
       y = "Measure") +
  theme_minimal() +
  facet_wrap(. ~measure, scales = "free_y")
rootedpd

#plotting DNA alpha diversity
DNA_conc_longer <- DNA_combined |>
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value") 
DNA_conc_longer$metric <- factor(DNA_conc_longer$metric, levels = order)

redox <- ggplot(DNA_conc_longer, aes(y=WS_O2, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Redox Level",
       x = "Value",
       y = "log O2 (uM)") +
  facet_grid(. ~metric, scales = "free_x") +
  scale_y_log10()
redox

depth <- ggplot(DNA_conc_longer, aes(y=placerun, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Depth",
       x = "Value",
       y = "Depth (m)") +
  facet_grid(. ~metric, scales = "free_x")
depth

sulfide <- ggplot(DNA_conc_longer, aes(y=WS_H2S, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Hydrogen Sulfide",
       x = "Value",
       y = "log H2S (uM)") +
  facet_grid(. ~metric, scales = "free_x")  +
  scale_y_log10()
sulfide

nitrate <- ggplot(DNA_conc_longer, aes(y=WS_NO3, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  labs(title = "DNA Alpha Diversity Metric against Nitrate",
       x = "Value",
       y = "NO3 (uM)") +
  facet_grid(. ~metric, scales = "free_x") 
nitrate

methane <- ggplot(DNA_conc_longer, aes(y=Mean_CH4, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  labs(title = "DNA Alpha Diversity Metric against Methane",
       x = "Value",
       y = "CH4 (uM)") +
  facet_grid(. ~metric, scales = "free_x") 
methane

##Figure S1: Alpha diversity metrics of NapA RNA in Saanich Inlet.
# plotting RNA alpha diversity against depth
alpha_napA <- read.csv("metatranscriptomes/SI_TS_NapA_alpha_diversiy_RNA.csv")
alpha_napA$placerun <- as.numeric(gsub("SI072_(\\d+)m.*", "\\1", alpha_napA$placerun))
alpha_napA_longer <- alpha_napA |>
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")
order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")
alpha_napA_longer$metric <- factor(alpha_napA_longer$metric, levels = order)
alpha_v_depth<- ggplot(alpha_napA_longer, aes(y=placerun, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "RNA Alpha Diversity Metric against Depth",
       x = "Value",
       y = "Depth") +
  facet_grid(. ~metric, scales = "free_x")
ggsave("alpha_div_vs_depth.png", alpha_v_depth, height = 6, width = 8, units = "in", dpi = 300)

# reading in and filtering Saanich data to include desired columns
saanich <- read.csv("Saanich_Data.csv") |>
  mutate(depth = Depth*1000) |>
  select(depth, WS_O2, WS_NO3, WS_H2S, Mean_CH4)

#calculating mean of each concentration at each depth
saanich_mod <- aggregate(cbind(WS_O2, WS_NO3, WS_H2S, Mean_CH4) ~ depth, data = saanich, FUN = mean)

# combining RNA alpha diversity with Saanich data
RNA_combined <- merge(alpha_napA, saanich_mod, by.x = "placerun", by.y = "depth")

#plotting RNA alpha diversity (rooted PD) against oxygen, sulfide, and nitrate concentrations
RNA_combined_longer <- RNA_combined |>
  pivot_longer(cols = c("WS_O2", "WS_NO3", "WS_H2S", "Mean_CH4", "placerun"),
               names_to = "measure",
               values_to = "value")
combined_order <- c("depth", "WS_O2", "WS_NO3", "WS_H2S", "Mean_CH4")
RNA_combined_longer$measure <- factor(RNA_combined_longer$measure, levels = combined_order)
rootedpd <- ggplot(RNA_combined_longer, aes(y=value, x=rooted_pd)) +
  geom_point(aes(colour=measure, shape=measure), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity",
       x = "Rooted PD",
       y = "Measure") +
  theme_minimal() +
  facet_wrap(. ~measure, scales = "free_y")
rootedpd

#plotting RNA alpha diversity
RNA_conc_longer <- RNA_combined |>
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value") 
RNA_conc_longer$metric <- factor(RNA_conc_longer$metric, levels = order)

redox <- ggplot(RNA_conc_longer, aes(y=WS_O2, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  labs(title = "RNA Alpha Diversity Metric against Redox Level",
       x = "Value",
       y = "log O2 (uM)") +
  scale_y_reverse()+
  facet_grid(. ~metric, scales = "free_x") +
  scale_y_log10()
redox

ggsave("redox.png", redox, height = 6, width = 8, units = "in", dpi = 300)


depth <- ggplot(RNA_conc_longer, aes(y=placerun, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "RNA Alpha Diversity Metric against Depth",
       x = "Value",
       y = "Depth (m)") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x")
depth

sulfide <- ggplot(RNA_conc_longer, aes(y=WS_H2S, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "RNA Alpha Diversity Metric against Hydrogen Sulfide",
       x = "Value",
       y = "log H2S (uM)") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x") +
  scale_y_log10()
sulfide

nitrate <- ggplot(RNA_conc_longer, aes(y=WS_NO3, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  labs(title = "RNA Alpha Diversity Metric against Nitrate",
       x = "Value",
       y = "NO3 (uM)") +
  facet_grid(. ~metric, scales = "free_x") 
nitrate
ggsave("nitrate.png", nitrate, height = 6, width = 8, units = "in", dpi = 300)

methane <- ggplot(RNA_conc_longer, aes(y=Mean_CH4, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  labs(title = "RNA Alpha Diversity Metric against Methane",
       x = "Value",
       y = "CH4 (uM)") +
  facet_grid(. ~metric, scales = "free_x") 
methane
ggsave("methane.png", methane, height = 6, width = 8, units = "in", dpi = 300)



### Question 3
# Load DNA and RNA datasets
dna <- read.delim("merged_depths_dna.tsv")
rna <- read.delim("merged_depths_rna.tsv")

# Add type column
dna$type <- "DNA"
rna$type <- "RNA"

# Combine into a single dataframe
merged <- bind_rows(dna, rna)

# Split the taxonomy into multiple columns
merged <- merged %>%
  separate(Taxonomy,
           into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ")

## Genus side-by-side stacked bar plot using top 10 most abundant overall genera
top_10_genus <- merged %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

bar_data <- merged %>%
  filter(Genus %in% top_10_genus) %>%
  group_by(Depth_m, type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Free Y-axis. DNA & RNA at different scales
ggplot(bar_data, aes(x = factor(Depth_m), y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Top 10 NapA-Carrying Genera Across Depths (DNA vs RNA)",
    x = "Depth (m)",
    y = "Total Abundance",
    fill = "Genus"
  )




### Question 4