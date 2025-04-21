library(tidyverse)
library(ggplot2)

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
