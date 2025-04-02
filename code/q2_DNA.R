library(tidyverse)
library(ggplot2)

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
  select(depth, WS_O2, WS_NO3, WS_H2S)

#calculating mean of each concentration at each depth
saanich_mod <- aggregate(cbind(WS_O2, WS_NO3, WS_H2S) ~ depth, data = saanich, FUN = mean)

# combining DNA alpha diversity with Saanich data
DNA_combined <- merge(alpha_napA, saanich_mod, by.x = "placerun", by.y = "depth")

#plotting DNA alpha diversity (rooted PD) against oxygen, sulfide, and nitrate concentrations
DNA_combined_longer <- DNA_combined |>
  pivot_longer(cols = c("WS_O2", "WS_NO3", "WS_H2S", "placerun"),
               names_to = "measure",
               values_to = "value")
combined_order <- c("depth", "WS_O2", "WS_NO3", "WS_H2S")
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
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x") +
  scale_y_log10()
redox

depth <- ggplot(DNA_conc_longer, aes(y=placerun, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Depth",
       x = "Value",
       y = "Depth (m)") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x")
depth

sulfide <- ggplot(DNA_conc_longer, aes(y=WS_H2S, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Hydrogen Sulfide",
       x = "Value",
       y = "log H2S (uM)") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x") +
  scale_y_log10()
sulfide

nitrate <- ggplot(DNA_conc_longer, aes(y=WS_NO3, x=value)) +
  geom_point(aes(colour=metric, shape=metric), size=4) +
  scale_y_reverse() +
  labs(title = "DNA Alpha Diversity Metric against Nitrate",
       x = "Value",
       y = "NO3 (uM)") +
  theme_minimal() +
  facet_grid(. ~metric, scales = "free_x") 
nitrate


# visualizing beta diversity
beta_napA <- read.csv("capstone/beta_diversity/SI_TS_NapA_beta_diversity.csv")





