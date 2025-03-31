# Load libraries
library(dplyr) # Wrangling
library(tidyr) # Reshape and pivot
library(ggplot2) # Visualize
library(pheatmap) # Heatmaps
library(vegan)  # Bray-Curtis and PERMANOVA
library(tibble)

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

## Genus
# Just to see
filtered_table_genus <- merged %>%
  filter(!is.na(Genus)) %>%
  select(Depth_m, type, Genus, Abundance)

# Group by Depth, Type, and Genus. Aggreate abundance
agg_abundance_genus <- merged %>%
  filter(!is.na(Genus)) %>%
  group_by(Depth_m, type, Genus) %>%
  summarise(TotalAbundance = sum(Abundance))

# Filter for top abundant genera for each depth
top_genus <- agg_abundance_genus %>%
  group_by(Depth_m, type) %>%
  slice_max(order_by = TotalAbundance, n = 1)

top_3_genus <- agg_abundance_genus %>%
  group_by(Depth_m, type) %>%
  slice_max(order_by = TotalAbundance, n = 3)


## Genus side-by-side stacked bar plot using top 3 abundant genera at each depth
ggplot(top_3_genus, aes(x = factor(Depth_m), y = TotalAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type) +
  theme_minimal() +
  labs(
    title = "Top 3 Abundant Genera at Each Depth (DNA vs RNA)",
    x = "Depth (m)",
    y = "Total Abundance",
    fill = "Genus"
  )

# Free Y-axis. DNA & RNA at different scales
ggplot(top_3_genus, aes(x = factor(Depth_m), y = TotalAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Top 3 Abundant Genera at Each Depth (DNA vs RNA)",
    x = "Depth (m)",
    y = "Total Abundance",
    fill = "Genus"
  )

# Relative abundance
top_3_genus_rel <- top_3_genus %>%
  group_by(Depth_m, type) %>%
  mutate(RelAbundance = TotalAbundance / sum(TotalAbundance)) %>%
  ungroup()

ggplot(top_3_genus_rel, aes(x = factor(Depth_m), y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type) +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Top 3 Genera by Depth (DNA vs RNA)",
    x = "Depth (m)",
    y = "Proportion",
    fill = "Genus"
  )


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

# To scale, hard to compare
ggplot(bar_data, aes(x = factor(Depth_m), y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type) +
  theme_minimal() +
  labs(
    title = "Top 10 NapA-Carrying Genera Across Depths (DNA vs RNA)",
    x = "Depth (m)",
    y = "Total Abundance",
    fill = "Genus"
  )

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

# Relative abundance. Compare composition with percentage
bar_data_rel <- bar_data %>%
  group_by(Depth_m, type) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

ggplot(bar_data_rel, aes(x = factor(Depth_m), y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type) +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Top 10 NapA-Carrying Genera (DNA vs RNA)",
    x = "Depth (m)",
    y = "Proportion",
    fill = "Genus"
  )


## Family
# Group by Depth, Type, and Family. Aggreate abundance
agg_abundance_family <- merged %>%
  filter(!is.na(Family)) %>%
  group_by(Depth_m, type, Family) %>%
  summarise(TotalAbundance = sum(Abundance))

# Filter for top abundant family for each depth
top_family <- agg_abundance_family %>%
  group_by(Depth_m, type) %>%
  slice_max(order_by = TotalAbundance, n = 1)

top_3_family <- agg_abundance_family %>%
  group_by(Depth_m, type) %>%
  slice_max(order_by = TotalAbundance, n = 3)

## Family side-by-side stacked bar plot using top 3 abundant faamily at each depth
# Free Y-axis. DNA & RNA at different scales
ggplot(top_3_family, aes(x = factor(Depth_m), y = TotalAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Top 3 Abundant Family at Each Depth (DNA vs RNA)",
    x = "Depth (m)",
    y = "Total Abundance",
    fill = "Genus"
  )

# Relative abundance
top_3_family_rel <- top_3_family %>%
  group_by(Depth_m, type) %>%
  mutate(RelAbundance = TotalAbundance / sum(TotalAbundance)) %>%
  ungroup()

ggplot(top_3_family_rel, aes(x = factor(Depth_m), y = RelAbundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~type) +
  theme_minimal() +
  labs(
    title = "Relative Abundance of Top 3 Family by Depth (DNA vs RNA)",
    x = "Depth (m)",
    y = "Proportion",
    fill = "Genus"
  )

## Code that didn't work well

# Pivot to compare DNA and RNA
abund_compare_genus <- top_genus %>%
  pivot_wider(names_from = type, values_from = TotalAbundance, values_fill = 0) %>%
  mutate(ExpressionRatio = ifelse(DNA > 0, RNA / DNA, NA))

# Plot scatter of DNA vs RNA abundance
ggplot(abund_compare_genus, aes(x = DNA, y = RNA, color = Genus)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  theme_minimal() +
  labs(title = "Top Genera: DNA vs RNA Abundance", x = "DNA Abundance", y = "RNA Abundance")
