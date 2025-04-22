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
saanich_data <- read.csv("Saanich_Data.csv")
dna <- read_tsv("merged_depths_dna.tsv")
rna <- read_tsv("merged_depths_rna.tsv")

view(saanich_data)
head(dna)

saanich_select <- saanich_data %>%
  select(Cruise, Date, Depth, WS_NO3, Mean_NH4, Std_NH4, Mean_NO2, Std_NO2, Mean_N2, Std_n2, Mean_N2O, Std_N2O) 

view(saanich_select) #check

#dna_abundance = dna$Abundance #testing 
#saanich_select$dna_abundance = dna_abundance #testing

# sum the dna abundance from each depth 
abundance_10m <- dna %>%
  filter(Sample == "SI072_10m_contig") 
sum_abundance_10m <- sum(abundance_10m$Abundance)


abundance_100m <- dna %>%
  filter(Sample == "SI072_100m_contig") 
sum_abundance_100m <- sum(abundance_100m$Abundance)

abundance_120m <- dna %>%
  filter(Sample == "SI072_120m_contig") 
sum_abundance_120m <- sum(abundance_120m$Abundance)

abundance_135m <- dna %>%
  filter(Sample == "SI072_135m_contig") 
sum_abundance_135m <- sum(abundance_135m$Abundance)

abundance_150m <- dna %>%
  filter(Sample == "SI072_150m_contig") 
sum_abundance_150m <- sum(abundance_150m$Abundance)

abundance_165m <- dna %>%
  filter(Sample == "SI072_165m_contig") 
sum_abundance_165m <- sum(abundance_165m$Abundance)

abundance_200m <- dna %>%
  filter(Sample == "SI072_200m_contig") 
sum_abundance_200m <- sum(abundance_200m$Abundance)

view(abundance_10m)


#WS_NO3 in each depth (dna)
WS_NO3_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(WS_NO3) 
count_WS_NO3_10m_dna <- sum(!is.na(WS_NO3_10m))
sum_WS_NO3_10m <- sum(WS_NO3_10m$WS_NO3, na.rm = TRUE) / count_WS_NO3_10m_dna


WS_NO3_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(WS_NO3) 
count_WS_NO3_100m_dna <- sum(!is.na(WS_NO3_100m))
sum_WS_NO3_100m <- sum(WS_NO3_100m$WS_NO3, na.rm = TRUE) / count_WS_NO3_100m_dna


WS_NO3_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
count_WS_NO3_120m_dna <- sum(!is.na(WS_NO3_120m))
sum_WS_NO3_120m <- sum(WS_NO3_120m$WS_NO3, na.rm = TRUE) / count_WS_NO3_120m_dna

WS_NO3_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(WS_NO3) 
count_WS_NO3_135m_dna <- sum(!is.na(WS_NO3_135m))
sum_WS_NO3_135m <- sum(WS_NO3_135m$WS_NO3, na.rm = TRUE) / count_WS_NO3_135m_dna

WS_NO3_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(WS_NO3) 
count_WS_NO3_150m_dna <- sum(!is.na(WS_NO3_150m))
sum_WS_NO3_150m <- sum(WS_NO3_150m$WS_NO3, na.rm = TRUE) / count_WS_NO3_150m_dna

WS_NO3_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
count_WS_NO3_165m_dna <- sum(!is.na(WS_NO3_165m))
sum_WS_NO3_165m <- sum(WS_NO3_165m$WS_NO3, na.rm = TRUE) / count_WS_NO3_165m_dna


WS_NO3_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
count_WS_NO3_200m_dna <- sum(!is.na(WS_NO3_200m))
sum_WS_NO3_200m <- sum(WS_NO3_200m$WS_NO3, na.rm = TRUE) / count_WS_NO3_200m_dna

sum_WS_NO3_200m

sum_WS_NO3_200m

#WS_NO3_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

WS_NO3_abundance <- data.frame(
  Nitrogen = c(sum_WS_NO3_10m, sum_WS_NO3_100m, sum_WS_NO3_120m, sum_WS_NO3_135m, sum_WS_NO3_150m, sum_WS_NO3_165m, sum_WS_NO3_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)

WS_NO3_abundance$Nitrogen

head(WS_NO3_abundance)
WS_NO3_abundance


WS_NO3_dna <- ggplot(WS_NO3_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst WS_NO3 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "WS_NO3 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(aes(label = Depth), hjust=0, vjust=1)

WS_NO3_dna #check





################################################



#Mean_NH4 in each depth (dna)
Mean_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NH4) 
count_Mean_NH4_10m_dna <- sum(!is.na(Mean_NH4_10m))
sum_Mean_NH4_10m <- sum(Mean_NH4_10m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_10m_dna


Mean_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
count_Mean_NH4_100m_dna <- sum(!is.na(Mean_NH4_100m))
sum_Mean_NH4_100m <- sum(Mean_NH4_100m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_100m_dna

Mean_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
count_Mean_NH4_120m_dna <- sum(!is.na(Mean_NH4_120m))
sum_Mean_NH4_120m <- sum(Mean_NH4_120m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_120m_dna

Mean_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NH4) 
count_Mean_NH4_135m_dna <- sum(!is.na(Mean_NH4_135m))
sum_Mean_NH4_135m <- sum(Mean_NH4_135m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_135m_dna

Mean_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NH4) 
count_Mean_NH4_150m_dna <- sum(!is.na(Mean_NH4_150m))
sum_Mean_NH4_150m <- sum(Mean_NH4_150m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_150m_dna

Mean_NH4_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
count_Mean_NH4_165m_dna <- sum(!is.na(Mean_NH4_165m))
sum_Mean_NH4_165m <- sum(Mean_NH4_165m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_165m_dna

Mean_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
count_Mean_NH4_200m_dna <- sum(!is.na(Mean_NH4_200m))
sum_Mean_NH4_200m <- sum(Mean_NH4_200m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_200m_dna


#Mean_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NH4_abundance <- data.frame(
  Nitrogen = c(sum_Mean_NH4_10m, sum_Mean_NH4_100m, sum_Mean_NH4_120m, sum_Mean_NH4_135m, sum_Mean_NH4_150m, sum_Mean_NH4_165m, sum_Mean_NH4_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Mean_NH4_dna <- ggplot(Mean_NH4_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_NH4 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Mean_NH4 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(aes(label = Depth), hjust=0, vjust=1)

Mean_NH4_dna #check

#######################################################

#Std_NH4 in each depth (dna)
Std_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NH4) 
count_Std_NH4_10m_dna <- sum(!is.na(Std_NH4_10m))
sum_Std_NH4_10m <- sum(Std_NH4_10m$Std_NH4, na.rm = TRUE) / count_Std_NH4_10m_dna


Std_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
count_Std_NH4_100m_dna <- sum(!is.na(Std_NH4_100m))
sum_Std_NH4_100m <- sum(Std_NH4_100m$Std_NH4, na.rm = TRUE) / count_Std_NH4_100m_dna

Std_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NH4) 
count_Std_NH4_120m_dna <- sum(!is.na(Std_NH4_120m))
sum_Std_NH4_120m <- sum(Std_NH4_120m$Std_NH4, na.rm = TRUE) / count_Std_NH4_120m_dna

Std_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NH4) 
count_Std_NH4_135m_dna <- sum(!is.na(Std_NH4_135m))
sum_Std_NH4_135m <- sum(Std_NH4_135m$Std_NH4, na.rm = TRUE) / count_Std_NH4_135m_dna

Std_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NH4) 
count_Std_NH4_150m_dna <- sum(!is.na(Std_NH4_150m))
sum_Std_NH4_150m <- sum(Std_NH4_150m$Std_NH4, na.rm = TRUE) / count_Std_NH4_150m_dna

Std_NH4_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
count_Std_NH4_165m_dna <- sum(!is.na(Std_NH4_165m))
sum_Std_NH4_165m <- sum(Std_NH4_165m$Std_NH4, na.rm = TRUE) / count_Std_NH4_165m_dna

Std_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
count_Std_NH4_200m_dna <- sum(!is.na(Std_NH4_200m))
sum_Std_NH4_200m <- sum(Std_NH4_200m$Std_NH4, na.rm = TRUE) / count_Std_NH4_200m_dna


#Std_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NH4_abundance <- data.frame(
  Nitrogen = c(sum_Std_NH4_10m, sum_Std_NH4_100m, sum_Std_NH4_120m, sum_Std_NH4_135m, sum_Std_NH4_150m, sum_Std_NH4_165m, sum_Std_NH4_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Std_NH4_dna <-ggplot(Std_NH4_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_NH4 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Std_NH4 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1, aes(label = Depth))

Std_NH4_dna #check

###############################

#Mean_NO2 in each depth (dna)
Mean_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NO2) 
count_Mean_NO2_10m_dna <- sum(!is.na(Mean_NO2_10m))
sum_Mean_NO2_10m <- sum(Mean_NO2_10m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_10m_dna


Mean_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NO2) 
count_Mean_NO2_100m_dna <- sum(!is.na(Mean_NO2_100m))
sum_Mean_NO2_100m <- sum(Mean_NO2_100m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_100m_dna

Mean_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_NO2) 
count_Mean_NO2_120m_dna <- sum(!is.na(Mean_NO2_120m))
sum_Mean_NO2_120m <- sum(Mean_NO2_120m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_120m_dna

Mean_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NO2) 
count_Mean_NO2_135m_dna <- sum(!is.na(Mean_NO2_135m))
sum_Mean_NO2_135m <- sum(Mean_NO2_135m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_135m_dna

Mean_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NO2) 
count_Mean_NO2_150m_dna <- sum(!is.na(Mean_NO2_150m))
sum_Mean_NO2_150m <- sum(Mean_NO2_150m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_150m_dna

Mean_NO2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
count_Mean_NO2_165m_dna <- sum(!is.na(Mean_NO2_165m))
sum_Mean_NO2_165m <- sum(Mean_NO2_165m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_165m_dna

Mean_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
count_Mean_NO2_200m_dna <- sum(!is.na(Mean_NO2_200m))
sum_Mean_NO2_200m <- sum(Mean_NO2_200m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_200m_dna


#Mean_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NO2_abundance <- data.frame(
  Nitrogen = c(sum_Mean_NO2_10m, sum_Mean_NO2_100m, sum_Mean_NO2_120m, sum_Mean_NO2_135m, sum_Mean_NO2_150m, sum_Mean_NO2_165m, sum_Mean_NO2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)

Mean_NO2_abundance

Mean_NO2_dna <- ggplot(Mean_NO2_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_NO2 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Mean_NO2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label = Depth))

Mean_NO2_dna
############################################

#Std_NO2 in each depth (dna)
Std_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NO2) 
count_Std_NO2_10m_dna <- sum(!is.na(Std_NO2_10m))
sum_Std_NO2_10m <- sum(Std_NO2_10m$Std_NO2, na.rm = TRUE) / count_Std_NO2_10m_dna


Std_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_NO2) 
count_Std_NO2_100m_dna <- sum(!is.na(Std_NO2_100m))
sum_Std_NO2_100m <- sum(Std_NO2_100m$Std_NO2, na.rm = TRUE) / count_Std_NO2_100m_dna


Std_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NO2) 
count_Std_NO2_120m_dna <- sum(!is.na(Std_NO2_120m))
sum_Std_NO2_120m <- sum(Std_NO2_120m$Std_NO2, na.rm = TRUE) / count_Std_NO2_120m_dna


Std_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NO2) 
count_Std_NO2_135m_dna <- sum(!is.na(Std_NO2_135m))
sum_Std_NO2_135m <- sum(Std_NO2_135m$Std_NO2, na.rm = TRUE) / count_Std_NO2_135m_dna


Std_NO2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
count_Std_NO2_165m_dna <- sum(!is.na(Std_NO2_165m))
sum_Std_NO2_165m <- sum(Std_NO2_165m$Std_NO2, na.rm = TRUE) / count_Std_NO2_165m_dna

Std_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NO2) 
count_Std_NO2_150m_dna <- sum(!is.na(Std_NO2_150m))
sum_Std_NO2_150m <- sum(Std_NO2_150m$Std_NO2, na.rm = TRUE) / count_Std_NO2_150m_dna


Std_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
count_Std_NO2_200m_dna <- sum(!is.na(Std_NO2_200m))
sum_Std_NO2_200m <- sum(Std_NO2_200m$Std_NO2, na.rm = TRUE) / count_Std_NO2_200m_dna



#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NO2_abundance <- data.frame(
  Nitrogen = c(sum_Std_NO2_10m, sum_Std_NO2_100m, sum_Std_NO2_120m, sum_Std_NO2_135m, sum_Std_NO2_150m, sum_Std_NO2_165m, sum_Std_NO2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)

Std_NO2_abundance

Std_NO2_dna <- ggplot(Std_NO2_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_NO2 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Std_NO2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label = Depth))

Std_NO2_dna
##################################################

#Mean_N2 in each depth (dna)
Mean_N2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2) 
count_Mean_N2_10m_dna <- sum(!is.na(Mean_N2_10m))
sum_Mean_N2_10m <- sum(Mean_N2_10m$Mean_N2, na.rm = TRUE) / count_Mean_N2_10m_dna


Mean_N2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2) 
count_Mean_N2_100m_dna <- sum(!is.na(Mean_N2_100m))
sum_Mean_N2_100m <- sum(Mean_N2_100m$Mean_N2, na.rm = TRUE) / count_Mean_N2_100m_dna

Mean_N2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2) 
count_Mean_N2_120m_dna <- sum(!is.na(Mean_N2_120m))
sum_Mean_N2_120m <- sum(Mean_N2_120m$Mean_N2, na.rm = TRUE) / count_Mean_N2_120m_dna

Mean_N2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2) 
count_Mean_N2_135m_dna <- sum(!is.na(Mean_N2_135m))
sum_Mean_N2_135m <- sum(Mean_N2_135m$Mean_N2, na.rm = TRUE) / count_Mean_N2_135m_dna

Mean_N2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2) 
count_Mean_N2_150m_dna <- sum(!is.na(Mean_N2_150m))
sum_Mean_N2_150m <- sum(Mean_N2_150m$Mean_N2, na.rm = TRUE) / count_Mean_N2_150m_dna

Mean_N2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
count_Mean_N2_165m_dna <- sum(!is.na(Mean_N2_165m))
sum_Mean_N2_165m <- sum(Mean_N2_165m$Mean_N2, na.rm = TRUE) / count_Mean_N2_165m_dna

Mean_N2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
count_Mean_N2_200m_dna <- sum(!is.na(Mean_N2_200m))
sum_Mean_N2_200m <- sum(Mean_N2_200m$Mean_N2, na.rm = TRUE) / count_Mean_N2_200m_dna


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2_abundance <- data.frame(
  Nitrogen = c(sum_Mean_N2_10m, sum_Mean_N2_100m, sum_Mean_N2_120m, sum_Mean_N2_135m, sum_Mean_N2_150m, sum_Mean_N2_165m, sum_Mean_N2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Mean_N2_dna <- ggplot(Mean_N2_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_N2 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Mean_N2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_N2_dna
#####################################################

#Std_n2 in each depth (dna)
Std_n2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_n2) 
count_Std_n2_10m_dna <- sum(!is.na(Std_n2_10m))
sum_Std_n2_10m <- sum(Std_n2_10m$Std_n2, na.rm = TRUE) / count_Std_n2_10m_dna


Std_n2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_n2) 
count_Std_n2_100m_dna <- sum(!is.na(Std_n2_100m))
sum_Std_n2_100m <- sum(Std_n2_100m$Std_n2, na.rm = TRUE) / count_Std_n2_100m_dna

Std_n2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_n2) 
count_Std_n2_120m_dna <- sum(!is.na(Std_n2_120m))
sum_Std_n2_120m <- sum(Std_n2_120m$Std_n2, na.rm = TRUE) / count_Std_n2_120m_dna

Std_n2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_n2) 
count_Std_n2_135m_dna <- sum(!is.na(Std_n2_135m))
sum_Std_n2_135m <- sum(Std_n2_135m$Std_n2, na.rm = TRUE) / count_Std_n2_135m_dna

Std_n2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_n2) 
count_Std_n2_150m_dna <- sum(!is.na(Std_n2_150m))
sum_Std_n2_150m <- sum(Std_n2_150m$Std_n2, na.rm = TRUE) / count_Std_n2_150m_dna

Std_n2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
count_Std_n2_165m_dna <- sum(!is.na(Std_n2_165m))
sum_Std_n2_165m <- sum(Std_n2_165m$Std_n2, na.rm = TRUE) / count_Std_n2_165m_dna

Std_n2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
count_Std_n2_200m_dna <- sum(!is.na(Std_n2_200m))
sum_Std_n2_200m <- sum(Std_n2_200m$Std_n2, na.rm = TRUE) / count_Std_n2_200m_dna


#Std_n2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_n2_abundance <- data.frame(
  Nitrogen = c(sum_Std_n2_10m, sum_Std_n2_100m, sum_Std_n2_120m, sum_Std_n2_135m, sum_Std_n2_150m, sum_Std_n2_165m, sum_Std_n2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Std_n2_dna <- ggplot(Std_n2_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_n2 concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Std_n2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_n2_dna

#######################################################

#Mean_N2O in each depth (dna)
Mean_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
count_Mean_N2O_10m_dna <- sum(!is.na(Mean_N2O_10m))
sum_Mean_N2O_10m <- sum(Mean_N2O_10m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_10m_dna


Mean_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2O) 
count_Mean_N2O_100m_dna <- sum(!is.na(Mean_N2O_100m))
sum_Mean_N2O_100m <- sum(Mean_N2O_100m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_100m_dna

Mean_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2O) 
count_Mean_N2O_120m_dna <- sum(!is.na(Mean_N2O_120m))
sum_Mean_N2O_120m <- sum(Mean_N2O_120m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_120m_dna

Mean_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2O) 
count_Mean_N2O_135m_dna <- sum(!is.na(Mean_N2O_135m))
sum_Mean_N2O_135m <- sum(Mean_N2O_135m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_135m_dna

Mean_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2O) 
count_Mean_N2O_150m_dna <- sum(!is.na(Mean_N2O_150m))
sum_Mean_N2O_150m <- sum(Mean_N2O_150m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_150m_dna

Mean_N2O_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
count_Mean_N2O_165m_dna <- sum(!is.na(Mean_N2O_165m))
sum_Mean_N2O_165m <- sum(Mean_N2O_165m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_165m_dna

Mean_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
count_Mean_N2O_200m_dna <- sum(!is.na(Mean_N2O_200m))
sum_Mean_N2O_200m <- sum(Mean_N2O_200m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_200m_dna


#SMean_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2O_abundance <- data.frame(
  Nitrogen = c(sum_Mean_N2O_10m, sum_Mean_N2O_100m, sum_Mean_N2O_120m, sum_Mean_N2O_135m, sum_Mean_N2O_150m, sum_Mean_N2O_165m, sum_Mean_N2O_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Mean_N2O_dna <- ggplot(Mean_N2O_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_N2O concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Mean_N2O concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_N2O_dna #check

###############################################

#Std_N2O in each depth (dna)
Std_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
count_Std_N2O_10m_dna <- sum(!is.na(Std_N2O_10m))
sum_Std_N2O_10m <- sum(Std_N2O_10m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_10m_dna


Std_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_N2O) 
count_Std_N2O_100m_dna <- sum(!is.na(Std_N2O_100m))
sum_Std_N2O_100m <- sum(Std_N2O_100m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_100m_dna

Std_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_N2O) 
count_Std_N2O_120m_dna <- sum(!is.na(Std_N2O_120m))
sum_Std_N2O_120m <- sum(Std_N2O_120m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_120m_dna

Std_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_N2O) 
count_Std_N2O_135m_dna <- sum(!is.na(Std_N2O_135m))
sum_Std_N2O_135m <- sum(Std_N2O_135m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_135m_dna

Std_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_N2O) 
count_Std_N2O_150m_dna <- sum(!is.na(Std_N2O_150m))
sum_Std_N2O_150m <- sum(Std_N2O_150m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_150m_dna

Std_N2O_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
count_Std_N2O_165m_dna <- sum(!is.na(Std_N2O_165m))
sum_Std_N2O_165m <- sum(Std_N2O_165m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_165m_dna

Std_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
count_Std_N2O_200m_dna <- sum(!is.na(Std_N2O_200m))
sum_Std_N2O_200m <- sum(Std_N2O_200m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_200m_dna

Std_N2O_abundance
#Std_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_N2O_abundance <- data.frame(
  Nitrogen = c(sum_Std_N2O_10m, sum_Std_N2O_100m, sum_Std_N2O_120m, sum_Std_N2O_135m, sum_Std_N2O_150m, sum_Std_N2O_165m, sum_Std_N2O_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_165m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Std_N2O_dna <- ggplot(Std_N2O_abundance, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_N2O concentration (μM)", 
       x = "NapA DNA Abundance",
       y = "Std_N2O concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_N2O_dna #check

###################################################

dna_graph <- ggarrange(Std_NO2_dna, Mean_NH4_dna, Std_NH4_dna, Mean_NO2_dna, Mean_N2_dna, Std_n2_dna, Mean_N2O_dna, Std_N2O_dna, WS_NO3_dna,
                       ncol = 3, nrow = 3,
                       heights = c(0.1, 0.1, 0.1))


dna_graph_selected_nitrogen_species <- ggarrange(Mean_NO2_dna, WS_NO3_dna,  
                                                    ncol = 1, nrow = 2,
                                                    heights = c(0.1, 0.1, 0.1))

dna_graph_mean_WS <- ggarrange(Mean_NH4_dna, Mean_NO2_dna, Mean_N2_dna, Mean_N2O_dna, WS_NO3_dna,
                       ncol = 2, nrow = 3,
                       heights = c(0.1, 0.1, 0.1))

dna_graph_mean_WS

#decided not to include other nitrogen species that do not directly associated with NapA plus the Std data cannot be calculated and combine 
#with the mean data for each nitrogen species without further clarification of the original data set

dna_graph_selected_nitrogen_species #check

#dna_abundance <- data.frame(
#  "Species" = c("WS_NO3", "Mean_NH4", "Std_NH4", "Mean_NO2", "Std_NO2", "Mean_N2", "Std_n2", "Mean_N2O", "Std_N2O"),
#  "Nitrogen Concentration (μM)" = c(sum_WS_NO3_10m, sum_Mean_NH4_10m, sum_Std_NH4_10m, sum_Mean_NO2_10m, sum_Std_NO2_10m, sum_Mean_N2_10m, sum_Std_n2_10m, sum_Mean_N2O_10m, sum_Std_N2O_10m)
  
#)



##############################################################################################################################################################################
#RNA 

# sum the dna abundance from each depth 
abundance_10m_rna <- rna %>%
  filter(Sample == "SI072_10m_contig_classified") 
sum_10m_rna <- sum(abundance_10m_rna$Abundance)


abundance_100m_rna <- rna %>%
  filter(Sample == "SI072_100m_contig_classified") 
sum_100m_rna <- sum(abundance_100m_rna$Abundance)

view(abundance_100m_rna)

abundance_120m_rna <- rna %>%
  filter(Sample == "SI072_120m_contig_classified") 
sum_120m_rna <- sum(abundance_120m_rna$Abundance)

abundance_135m_rna <- rna %>%
  filter(Sample == "SI072_135m_contig_classified") 
sum_135m_rna <- sum(abundance_135m_rna$Abundance)

abundance_150m_rna <- rna %>%
  filter(Sample == "SI072_150m_contig_classified") 
sum_150m_rna <- sum(abundance_150m_rna$Abundance)

abundance_165m_rna <- rna %>%
  filter(Sample == "SI072_165m_contig_classified") 
sum_165m_rna <- sum(abundance_165m_rna$Abundance)

abundance_200m_rna <- rna %>%
  filter(Sample == "SI072_200m_contig_classified") 
sum_200m_rna <- sum(abundance_200m_rna$Abundance)


#WS_NO3 in each depth (rna)
WS_NO3_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(WS_NO3) 
count_WS_NO3_10m_rna <- sum(!is.na(WS_NO3_10m))
sum_WS_NO3_10m_rna <- sum(WS_NO3_10m$WS_NO3, na.rm = TRUE) / count_WS_NO3_10m_rna


WS_NO3_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(WS_NO3) 
count_WS_NO3_100m_rna <- sum(!is.na(WS_NO3_100m))
sum_WS_NO3_100m_rna <- sum(WS_NO3_100m$WS_NO3, na.rm = TRUE) / count_WS_NO3_100m_rna

WS_NO3_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
count_WS_NO3_120m_rna <- sum(!is.na(WS_NO3_120m))
sum_WS_NO3_120m_rna <- sum(WS_NO3_120m$WS_NO3, na.rm = TRUE) / count_WS_NO3_120m_rna

WS_NO3_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(WS_NO3) 
count_WS_NO3_135m_rna <- sum(!is.na(WS_NO3_135m))
sum_WS_NO3_135m_rna <- sum(WS_NO3_135m$WS_NO3, na.rm = TRUE) / count_WS_NO3_135m_rna

WS_NO3_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(WS_NO3) 
count_WS_NO3_150m_rna <- sum(!is.na(WS_NO3_150m))
sum_WS_NO3_150m_rna <- sum(WS_NO3_150m$WS_NO3, na.rm = TRUE) / count_WS_NO3_150m_rna

WS_NO3_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
count_WS_NO3_165m_rna <- sum(!is.na(WS_NO3_165m))
sum_WS_NO3_165m_rna <- sum(WS_NO3_165m$WS_NO3, na.rm = TRUE) / count_WS_NO3_165m_rna

WS_NO3_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
count_WS_NO3_200m_rna <- sum(!is.na(WS_NO3_200m))
sum_WS_NO3_200m_rna <- sum(WS_NO3_200m$WS_NO3, na.rm = TRUE) / count_WS_NO3_200m_rna



#WS_NO3_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

WS_NO3_abundance_rna <- data.frame(
  Nitrogen = c(sum_WS_NO3_10m_rna, sum_WS_NO3_100m_rna, sum_WS_NO3_120m_rna, sum_WS_NO3_135m_rna, sum_WS_NO3_150m_rna, sum_WS_NO3_165m_rna, sum_WS_NO3_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)



WS_NO3_rna <- ggplot(WS_NO3_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst WS_NO3 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "WS_NO3 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

WS_NO3_rna #check



################################################



#Mean_NH4 in each depth (rna)
Mean_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NH4) 
count_Mean_NH4_10m_rna <- sum(!is.na(Mean_NH4_10m))
sum_Mean_NH4_10m_rna <- sum(Mean_NH4_10m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_10m_rna


Mean_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
count_Mean_NH4_100m_rna <- sum(!is.na(Mean_NH4_100m))
sum_Mean_NH4_100m_rna <- sum(Mean_NH4_100m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_100m_rna

Mean_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
count_Mean_NH4_120m_rna <- sum(!is.na(Mean_NH4_120m))
sum_Mean_NH4_120m_rna <- sum(Mean_NH4_120m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_120m_rna

Mean_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NH4) 
count_Mean_NH4_135m_rna <- sum(!is.na(Mean_NH4_135m))
sum_Mean_NH4_135m_rna <- sum(Mean_NH4_135m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_135m_rna

Mean_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NH4) 
count_Mean_NH4_150m_rna <- sum(!is.na(Mean_NH4_150m))
sum_Mean_NH4_150m_rna <- sum(Mean_NH4_150m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_150m_rna

Mean_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
count_Mean_NH4_165m_rna <- sum(!is.na(Mean_NH4_165m))
sum_Mean_NH4_165m_rna <- sum(Mean_NH4_165m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_165m_rna

Mean_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
count_Mean_NH4_200m_rna <- sum(!is.na(Mean_NH4_200m))
sum_Mean_NH4_200m_rna <- sum(Mean_NH4_200m$Mean_NH4, na.rm = TRUE) / count_Mean_NH4_200m_rna


#Mean_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NH4_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_NH4_10m_rna, sum_Mean_NH4_100m_rna, sum_Mean_NH4_120m_rna, sum_Mean_NH4_135m_rna, sum_Mean_NH4_150m_rna, sum_Mean_NH4_165m_rna, sum_Mean_NH4_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Mean_NH4_rna <- ggplot(Mean_NH4_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_NH4 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Mean_NH4 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_NH4_rna #check

#######################################################

#Std_NH4 in each depth (rna)
Std_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NH4) 
count_Std_NH4_10m_rna <- sum(!is.na(Std_NH4_10m))
sum_Std_NH4_10m_rna <- sum(Std_NH4_10m$Std_NH4, na.rm = TRUE) / count_Std_NH4_10m_rna


Std_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
count_Std_NH4_100m_rna <- sum(!is.na(Std_NH4_100m))
sum_Std_NH4_100m_rna <- sum(Std_NH4_100m$Std_NH4, na.rm = TRUE) / count_Std_NH4_100m_rna

Std_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NH4) 
count_Std_NH4_120m_rna <- sum(!is.na(Std_NH4_120m))
sum_Std_NH4_120m_rna <- sum(Std_NH4_120m$Std_NH4, na.rm = TRUE) / count_Std_NH4_120m_rna

Std_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NH4) 
count_Std_NH4_135m_rna <- sum(!is.na(Std_NH4_135m))
sum_Std_NH4_135m_rna <- sum(Std_NH4_135m$Std_NH4, na.rm = TRUE) / count_Std_NH4_135m_rna

Std_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NH4) 
count_Std_NH4_150m_rna <- sum(!is.na(Std_NH4_150m))
sum_Std_NH4_150m_rna <- sum(Std_NH4_150m$Std_NH4, na.rm = TRUE) / count_Std_NH4_150m_rna

Std_NH4_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
count_Std_NH4_165m_rna <- sum(!is.na(Std_NH4_165m))
sum_Std_NH4_165m_rna <- sum(Std_NH4_165m$Std_NH4, na.rm = TRUE) / count_Std_NH4_165m_rna

Std_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
count_Std_NH4_200m_rna <- sum(!is.na(Std_NH4_200m))
sum_Std_NH4_200m_rna <- sum(Std_NH4_200m$Std_NH4, na.rm = TRUE) / count_Std_NH4_200m_rna


#Std_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NH4_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_NH4_10m_rna, sum_Std_NH4_100m_rna, sum_Std_NH4_120m_rna, sum_Std_NH4_135m_rna, sum_Std_NH4_150m_rna, sum_Std_NH4_165m_rna, sum_Std_NH4_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna,sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)



Std_NH4_rna <- ggplot(Std_NH4_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_NH4 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Std_NH4 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_NH4_rna #check

###############################

#Mean_NO2 in each depth (rna)
Mean_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NO2) 
count_Mean_NO2_10m_rna <- sum(!is.na(Mean_NO2_10m))
sum_Mean_NO2_10m_rna <- sum(Mean_NO2_10m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_10m_rna

Mean_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NO2) 
count_Mean_NO2_100m_rna <- sum(!is.na(Mean_NO2_100m))
sum_Mean_NO2_100m_rna <- sum(Mean_NO2_100m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_100m_rna

Mean_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_NO2) 
count_Mean_NO2_120m_rna <- sum(!is.na(Mean_NO2_120m))
sum_Mean_NO2_120m_rna <- sum(Mean_NO2_120m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_120m_rna

Mean_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NO2) 
count_Mean_NO2_135m_rna <- sum(!is.na(Mean_NO2_135m))
sum_Mean_NO2_135m_rna <- sum(Mean_NO2_135m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_135m_rna

Mean_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NO2) 
count_Mean_NO2_150m_rna <- sum(!is.na(Mean_NO2_150m))
sum_Mean_NO2_150m_rna <- sum(Mean_NO2_150m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_150m_rna

Mean_NO2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
count_Mean_NO2_165m_rna <- sum(!is.na(Mean_NO2_165m))
sum_Mean_NO2_165m_rna <- sum(Mean_NO2_165m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_165m_rna

Mean_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
count_Mean_NO2_200m_rna <- sum(!is.na(Mean_NO2_200m))
sum_Mean_NO2_200m_rna <- sum(Mean_NO2_200m$Mean_NO2, na.rm = TRUE) / count_Mean_NO2_200m_rna


#Mean_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NO2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_NO2_10m_rna, sum_Mean_NO2_100m_rna, sum_Mean_NO2_120m_rna, sum_Mean_NO2_135m_rna, sum_Mean_NO2_150m_rna, sum_Mean_NO2_165m_rna, sum_Mean_NO2_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Mean_NO2_rna <- ggplot(Mean_NO2_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_NO2 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Mean_NO2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_NO2_rna #check

############################################

#Std_NO2 in each depth (rna)
Std_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NO2) 
count_Std_NO2_10m_rna <- sum(!is.na(Std_NO2_10m))
sum_Std_NO2_10m_rna <- sum(Std_NO2_10m$Std_NO2, na.rm = TRUE) / count_Std_NO2_10m_rna


Std_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_NO2) 
count_Std_NO2_100m_rna <- sum(!is.na(Std_NO2_100m))
sum_Std_NO2_100m_rna <- sum(Std_NO2_100m$Std_NO2, na.rm = TRUE) / count_Std_NO2_100m_rna

Std_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NO2) 
count_Std_NO2_120m_rna <- sum(!is.na(Std_NO2_120m))
sum_Std_NO2_120m_rna <- sum(Std_NO2_120m$Std_NO2, na.rm = TRUE) / count_Std_NO2_120m_rna

Std_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NO2) 
count_Std_NO2_135m_rna <- sum(!is.na(Std_NO2_135m))
sum_Std_NO2_135m_rna <- sum(Std_NO2_135m$Std_NO2, na.rm = TRUE) / count_Std_NO2_135m_rna

Std_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NO2) 
count_Std_NO2_150m_rna <- sum(!is.na(Std_NO2_150m))
sum_Std_NO2_150m_rna <- sum(Std_NO2_150m$Std_NO2, na.rm = TRUE) / count_Std_NO2_150m_rna

Std_NO2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
count_Std_NO2_165m_rna <- sum(!is.na(Std_NO2_165m))
sum_Std_NO2_165m_rna <- sum(Std_NO2_165m$Std_NO2, na.rm = TRUE) / count_Std_NO2_165m_rna

Std_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
count_Std_NO2_200m_rna <- sum(!is.na(Std_NO2_200m))
sum_Std_NO2_200m_rna <- sum(Std_NO2_200m$Std_NO2, na.rm = TRUE) / count_Std_NO2_200m_rna


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NO2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_NO2_10m_rna, sum_Std_NO2_100m_rna, sum_Std_NO2_120m_rna, sum_Std_NO2_135m_rna, sum_Std_NO2_150m_rna, sum_Std_NO2_165m_rna, sum_Std_NO2_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Std_NO2_rna <- ggplot(Std_NO2_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_NO2 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Std_NO2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_NO2_rna #check

##################################################

#Mean_N2 in each depth (dna)
Mean_N2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2) 
count_Mean_N2_10m_rna <- sum(!is.na(Mean_N2_10m))
sum_Mean_N2_10m_rna <- sum(Mean_N2_10m$Mean_N2, na.rm = TRUE) / count_Mean_N2_10m_rna


Mean_N2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2) 
count_Mean_N2_100m_rna <- sum(!is.na(Mean_N2_100m))
sum_Mean_N2_100m_rna <- sum(Mean_N2_100m$Mean_N2, na.rm = TRUE) / count_Mean_N2_100m_rna

Mean_N2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2) 
count_Mean_N2_120m_rna <- sum(!is.na(Mean_N2_120m))
sum_Mean_N2_120m_rna <- sum(Mean_N2_120m$Mean_N2, na.rm = TRUE) / count_Mean_N2_120m_rna

Mean_N2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2) 
count_Mean_N2_135m_rna <- sum(!is.na(Mean_N2_135m))
sum_Mean_N2_135m_rna <- sum(Mean_N2_135m$Mean_N2, na.rm = TRUE) / count_Mean_N2_135m_rna

Mean_N2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2) 
count_Mean_N2_150m_rna <- sum(!is.na(Mean_N2_150m))
sum_Mean_N2_150m_rna <- sum(Mean_N2_150m$Mean_N2, na.rm = TRUE) / count_Mean_N2_150m_rna

Mean_N2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
count_Mean_N2_165m_rna <- sum(!is.na(Mean_N2_165m))
sum_Mean_N2_165m_rna <- sum(Mean_N2_165m$Mean_N2, na.rm = TRUE) / count_Mean_N2_165m_rna

Mean_N2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
count_Mean_N2_200m_rna <- sum(!is.na(Mean_N2_200m))
sum_Mean_N2_200m_rna <- sum(Mean_N2_200m$Mean_N2, na.rm = TRUE) / count_Mean_N2_200m_rna


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_N2_10m_rna, sum_Mean_N2_100m_rna, sum_Mean_N2_120m_rna, sum_Mean_N2_135m_rna, sum_Mean_N2_150m_rna, sum_Mean_N2_165m_rna, sum_Mean_N2_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)
Mean_N2_abundance_rna

Mean_N2_rna <- ggplot(Mean_N2_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_N2 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Mean_N2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_N2_rna #check

#####################################################

#Std_n2 in each depth (rna)
Std_n2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_n2) 
count_Std_n2_10m_rna <- sum(!is.na(Std_n2_10m))
sum_Std_n2_10m_rna <- sum(Std_n2_10m$Std_n2, na.rm = TRUE) / count_Std_n2_10m_rna


Std_n2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_n2) 
count_Std_n2_100m_rna <- sum(!is.na(Std_n2_100m))
sum_Std_n2_100m_rna <- sum(Std_n2_100m$Std_n2, na.rm = TRUE) / count_Std_n2_100m_rna

Std_n2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_n2) 
count_Std_n2_120m_rna <- sum(!is.na(Std_n2_120m))
sum_Std_n2_120m_rna <- sum(Std_n2_120m$Std_n2, na.rm = TRUE) / count_Std_n2_120m_rna

Std_n2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_n2) 
count_Std_n2_135m_rna <- sum(!is.na(Std_n2_135m))
sum_Std_n2_135m_rna <- sum(Std_n2_135m$Std_n2, na.rm = TRUE) / count_Std_n2_135m_rna

Std_n2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_n2) 
count_Std_n2_150m_rna <- sum(!is.na(Std_n2_150m))
sum_Std_n2_150m_rna <- sum(Std_n2_150m$Std_n2, na.rm = TRUE) / count_Std_n2_150m_rna

Std_n2_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
count_Std_n2_165m_rna <- sum(!is.na(Std_n2_165m))
sum_Std_n2_165m_rna <- sum(Std_n2_165m$Std_n2, na.rm = TRUE) / count_Std_n2_165m_rna

Std_n2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
count_Std_n2_200m_rna <- sum(!is.na(Std_n2_200m))
sum_Std_n2_200m_rna <- sum(Std_n2_200m$Std_n2, na.rm = TRUE) / count_Std_n2_200m_rna


#Std_n2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_n2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_n2_10m_rna, sum_Std_n2_100m_rna, sum_Std_n2_120m_rna, sum_Std_n2_135m_rna, sum_Std_n2_150m_rna, sum_Std_n2_165m_rna, sum_Std_n2_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna), 
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)


Std_n2_rna <- ggplot(Std_n2_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_n2 concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Std_n2 concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_n2_rna #check

#######################################################

#Mean_N2O in each depth (rna)
Mean_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
count_Mean_N2O_10m_rna <- sum(!is.na(Mean_N2O_10m))
sum_Mean_N2O_10m_rna <- sum(Mean_N2O_10m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_10m_rna


Mean_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2O) 
count_Mean_N2O_100m_rna <- sum(!is.na(Mean_N2O_100m))
sum_Mean_N2O_100m_rna <- sum(Mean_N2O_100m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_100m_rna

Mean_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2O) 
count_Mean_N2O_120m_rna <- sum(!is.na(Mean_N2O_120m))
sum_Mean_N2O_120m_rna <- sum(Mean_N2O_120m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_120m_rna

Mean_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2O) 
count_Mean_N2O_135m_rna <- sum(!is.na(Mean_N2O_135m))
sum_Mean_N2O_135m_rna <- sum(Mean_N2O_135m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_135m_rna

Mean_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2O) 
count_Mean_N2O_150m_rna <- sum(!is.na(Mean_N2O_150m))
sum_Mean_N2O_150m_rna <- sum(Mean_N2O_150m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_150m_rna

Mean_N2O_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
count_Mean_N2O_165m_rna <- sum(!is.na(Mean_N2O_165m))
sum_Mean_N2O_165m_rna <- sum(Mean_N2O_165m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_165m_rna

Mean_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
count_Mean_N2O_200m_rna <- sum(!is.na(Mean_N2O_200m))
sum_Mean_N2O_200m_rna <- sum(Mean_N2O_200m$Mean_N2O, na.rm = TRUE) / count_Mean_N2O_200m_rna


#SMean_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2O_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_N2O_10m_rna, sum_Mean_N2O_100m_rna, sum_Mean_N2O_120m_rna, sum_Mean_N2O_135m_rna, sum_Mean_N2O_150m_rna, sum_Mean_N2O_165m_rna, sum_Mean_N2O_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)
Mean_N2O_abundance_rna

Mean_N2O_rna <- ggplot(Mean_N2O_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_N2O concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Mean_N2O concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Mean_N2O_rna #check

###############################################

#Std_N2O in each depth (rna)
Std_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
count_Std_N2O_10m_rna <- sum(!is.na(Std_N2O_10m))
sum_Std_N2O_10m_rna <- sum(Std_N2O_10m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_10m_rna


Std_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_N2O) 
count_Std_N2O_100m_rna <- sum(!is.na(Std_N2O_100m))
sum_Std_N2O_100m_rna <- sum(Std_N2O_100m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_100m_rna

Std_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_N2O) 
count_Std_N2O_120m_rna <- sum(!is.na(Std_N2O_120m))
sum_Std_N2O_120m_rna <- sum(Std_N2O_120m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_120m_rna

Std_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_N2O) 
count_Std_N2O_135m_rna <- sum(!is.na(Std_N2O_135m))
sum_Std_N2O_135m_rna <- sum(Std_N2O_135m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_135m_rna

Std_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_N2O) 
count_Std_N2O_150m_rna <- sum(!is.na(Std_N2O_150m))
sum_Std_N2O_150m_rna <- sum(Std_N2O_150m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_150m_rna

Std_N2O_165m <- saanich_select %>%
  filter(Depth <= 0.165) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
count_Std_N2O_165m_rna <- sum(!is.na(Std_N2O_165m))
sum_Std_N2O_165m_rna <- sum(Std_N2O_165m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_165m_rna

Std_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
count_Std_N2O_200m_rna <- sum(!is.na(Std_N2O_200m))
sum_Std_N2O_200m_rna <- sum(Std_N2O_200m$Mean_N2O, na.rm = TRUE) / count_Std_N2O_200m_rna


#Std_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 rna abundance data, therefore excluding this depth concentration

Std_N2O_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_N2O_10m_rna, sum_Std_N2O_100m_rna, sum_Std_N2O_120m_rna, sum_Std_N2O_135m_rna, sum_Std_N2O_150m_rna, sum_Std_N2O_165m_rna, sum_Std_N2O_200m_rna), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_165m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "165m", "200m")
)

Std_N2O_abundance_rna

Std_N2O_rna <- ggplot(Std_N2O_abundance_rna, aes(y = Nitrogen, x = Abundance)) +
  geom_point(aes(color = Depth), size = 4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black", linetype = "solid", linewidth = 0.7) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_N2O concentration (μM)", 
       x = "NapA RNA Abundance",
       y = "Std_N2O concentration (μM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2, aes(label =  Depth))

Std_N2O_rna #check

###################################################

rna_graph <- ggarrange(Std_NO2_rna, Mean_NH4_rna, Std_NH4_rna, Mean_NO2_rna, Mean_N2_rna, Std_n2_rna, Mean_N2O_rna, Std_N2O_rna, WS_NO3_rna,
                       ncol = 3, nrow = 3,
                       heights = c(0.1, 0.1, 0.1))

rna_graph_selected_nitrogen_species <- ggarrange(Mean_NO2_rna, WS_NO3_rna,
                                                 ncol = 1, nrow = 2,
                                                 heights = c(0.1, 0.1, 0.1))
rna_graph_selected_nitrogen_species #testing

rna_graph #testing
dna_graph #testing

rna_graph_mean_WS <- ggarrange(Mean_NH4_rna, Std_NH4_rna, Mean_NO2_rna, Mean_N2_rna, Std_n2_rna, Mean_N2O_rna, Std_N2O_rna, WS_NO3_rna,
                       ncol = 3, nrow = 3,
                       heights = c(0.1, 0.1, 0.1))

###################################################
# merging graphs 
merge <- ggarrange(dna_graph, rna_graph, 
                   ncol = 1, nrow = 2,
                   heights = c(0.1, 0.1, 0.1))

#install.package("patchwork") #testing
#library(patchwork) #testing


merge_selected_nitrogen_species <- ggarrange(dna_graph_selected_nitrogen_species, rna_graph_selected_nitrogen_species, 
                                             ncol = 2, nrow = 1,
                                             heights = c(0.1, 0.1, 0.1))
merge_selected_nitrogen_species #testing



merge_other <- ggarrange(Mean_NH4_dna, Mean_NH4_rna, Mean_NO2_dna, Mean_NO2_rna, Mean_N2_dna, Mean_N2_rna, Mean_N2O_dna, Mean_N2O_rna, WS_NO3_dna, WS_NO3_rna,
                               ncol = 2, nrow = 5,
                               heights = c(0.1, 0.1, 0.1))

merge_other
