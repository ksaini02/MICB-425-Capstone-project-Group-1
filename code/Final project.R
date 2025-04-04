install.packages("tidyverse")
install.packages("pheatmap")
install.packages("ggpubr")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggpubr)

saanich_data <- read.csv("Saanich_Data.csv")
dna <- read_tsv("merged_depths_dna.tsv")
rna <- read_tsv("merged_depths_rna.tsv")

view(saanich_data)
head(dna)

saanich_select <- saanich_data %>%
  select(Cruise, Date, Depth, WS_NO3, Mean_NH4, Std_NH4, Mean_NO2, Std_NO2, Mean_N2, Std_n2, Mean_N2O, Std_N2O) 

view(saanich_select)

dna_abundance = dna$Abundance
saanich_select$dna_abundance = dna_abundance

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
sum_WS_NO3_10m <- sum(WS_NO3_10m, na.rm = TRUE)

sum_WS_NO3_10m

WS_NO3_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(WS_NO3) 
sum_WS_NO3_100m <- sum(WS_NO3_100m, na.rm = TRUE)

WS_NO3_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
sum_WS_NO3_120m <- sum(WS_NO3_120m, na.rm = TRUE)

WS_NO3_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(WS_NO3) 
sum_WS_NO3_135m <- sum(WS_NO3_135m, na.rm = TRUE)

WS_NO3_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(WS_NO3) 
sum_WS_NO3_150m <- sum(WS_NO3_150m, na.rm = TRUE)

WS_NO3_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
sum_WS_NO3_200m <- sum(WS_NO3_200m, na.rm = TRUE)



#WS_NO3_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

WS_NO3_abundance <- data.frame(
  Nitrogen = c(sum_WS_NO3_10m, sum_WS_NO3_100m, sum_WS_NO3_120m, sum_WS_NO3_135m, sum_WS_NO3_150m, sum_WS_NO3_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)

head(WS_NO3_abundance)

WS_NO3_dna <- ggplot(WS_NO3_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
#  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst WS_NO3 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "WS_NO3 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)



################################################



#Mean_NH4 in each depth (dna)
Mean_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NH4) 
sum_Mean_NH4_10m <- sum(Mean_NH4_10m, na.rm = TRUE)


Mean_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
sum_Mean_NH4_100m <- sum(Mean_NH4_100m, na.rm = TRUE)

Mean_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
sum_Mean_NH4_120m <- sum(Mean_NH4_120m, na.rm = TRUE)

Mean_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NH4) 
sum_Mean_NH4_135m <- sum(Mean_NH4_135m, na.rm = TRUE)

Mean_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NH4) 
sum_Mean_NH4_150m <- sum(Mean_NH4_150m, na.rm = TRUE)

Mean_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
sum_Mean_NH4_200m <- sum(Mean_NH4_200m, na.rm = TRUE)


#Mean_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NH4_abundance <- data.frame(
  Nitrogen = c(sum_Mean_NH4_10m, sum_Mean_NH4_100m, sum_Mean_NH4_120m, sum_Mean_NH4_135m, sum_Mean_NH4_150m, sum_Mean_NH4_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_NH4_dna <- ggplot(Mean_NH4_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
#  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_NH4 concentration (uM)", 
      x = "NapA DNA Abundance",
      y = "Mean_NH4 concentration (uM)", 
      color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)


#######################################################

#Std_NH4 in each depth (dna)
Std_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NH4) 
sum_Std_NH4_10m <- sum(Std_NH4_10m, na.rm = TRUE)


Std_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
sum_Std_NH4_100m <- sum(Std_NH4_100m, na.rm = TRUE)

Std_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NH4) 
sum_Std_NH4_120m <- sum(Std_NH4_120m, na.rm = TRUE)

Std_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NH4) 
sum_Std_NH4_135m <- sum(Std_NH4_135m, na.rm = TRUE)

Std_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NH4) 
sum_Std_NH4_150m <- sum(Std_NH4_150m, na.rm = TRUE)

Std_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
sum_Std_NH4_200m <- sum(Std_NH4_200m, na.rm = TRUE)


#Std_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NH4_abundance <- data.frame(
  Nitrogen = c(sum_Std_NH4_10m, sum_Std_NH4_100m, sum_Std_NH4_120m, sum_Std_NH4_135m, sum_Std_NH4_150m, sum_Std_NH4_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_NH4_dna <-ggplot(Std_NH4_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
#  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_NH4 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Std_NH4 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)


###############################

#Mean_NO2 in each depth (dna)
Mean_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NO2) 
sum_Mean_NO2_10m <- sum(Mean_NO2_10m, na.rm = TRUE)


Mean_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NO2) 
sum_Mean_NO2_100m <- sum(Mean_NO2_100m, na.rm = TRUE)

Mean_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_NO2) 
sum_Mean_NO2_120m <- sum(Mean_NO2_120m, na.rm = TRUE)

Mean_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NO2) 
sum_Mean_NO2_135m <- sum(Mean_NO2_135m, na.rm = TRUE)

Mean_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NO2) 
sum_Mean_NO2_150m <- sum(Mean_NO2_150m, na.rm = TRUE)

Mean_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
sum_Mean_NO2_200m <- sum(Mean_NO2_200m, na.rm = TRUE)


#Mean_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NO2_abundance <- data.frame(
  Nitrogen = c(sum_Mean_NO2_10m, sum_Mean_NO2_100m, sum_Mean_NO2_120m, sum_Mean_NO2_135m, sum_Mean_NO2_150m, sum_Mean_NO2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_NO2_dna <- ggplot(Mean_NO2_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
#  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_NO2 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Mean_NO2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)

############################################

#Std_NO2 in each depth (dna)
Std_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NO2) 
sum_Std_NO2_10m <- sum(Std_NO2_10m, na.rm = TRUE)


Std_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_NO2) 
sum_Std_NO2_100m <- sum(Std_NO2_100m, na.rm = TRUE)

Std_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NO2) 
sum_Std_NO2_120m <- sum(Std_NO2_120m, na.rm = TRUE)

Std_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NO2) 
sum_Std_NO2_135m <- sum(Std_NO2_135m, na.rm = TRUE)

Std_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NO2) 
sum_Std_NO2_150m <- sum(Std_NO2_150m, na.rm = TRUE)

Std_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
sum_Std_NO2_200m <- sum(Std_NO2_200m, na.rm = TRUE)


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NO2_abundance <- data.frame(
  Nitrogen = c(sum_Std_NO2_10m, sum_Std_NO2_100m, sum_Std_NO2_120m, sum_Std_NO2_135m, sum_Std_NO2_150m, sum_Std_NO2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_NO2_dna <- ggplot(Std_NO2_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_NO2 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Std_NO2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)

##################################################

#Mean_N2 in each depth (dna)
Mean_N2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2) 
sum_Mean_N2_10m <- sum(Mean_N2_10m, na.rm = TRUE)


Mean_N2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2) 
sum_Mean_N2_100m <- sum(Mean_N2_100m, na.rm = TRUE)

Mean_N2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2) 
sum_Mean_N2_120m <- sum(Mean_N2_120m, na.rm = TRUE)

Mean_N2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2) 
sum_Mean_N2_135m <- sum(Mean_N2_135m, na.rm = TRUE)

Mean_N2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2) 
sum_Mean_N2_150m <- sum(Mean_N2_150m, na.rm = TRUE)

Mean_N2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
sum_Mean_N2_200m <- sum(Mean_N2_200m, na.rm = TRUE)


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2_abundance <- data.frame(
  Nitrogen = c(sum_Mean_N2_10m, sum_Mean_N2_100m, sum_Mean_N2_120m, sum_Mean_N2_135m, sum_Mean_N2_150m, sum_Mean_N2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_N2_dna <- ggplot(Mean_N2_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_N2 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Mean_N2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)

#####################################################

#Std_n2 in each depth (dna)
Std_n2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_n2) 
sum_Std_n2_10m <- sum(Std_n2_10m, na.rm = TRUE)


Std_n2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_n2) 
sum_Std_n2_100m <- sum(Std_n2_100m, na.rm = TRUE)

Std_n2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_n2) 
sum_Std_n2_120m <- sum(Std_n2_120m, na.rm = TRUE)

Std_n2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_n2) 
sum_Std_n2_135m <- sum(Std_n2_135m, na.rm = TRUE)

Std_n2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_n2) 
sum_Std_n2_150m <- sum(Std_n2_150m, na.rm = TRUE)

Std_n2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
sum_Std_n2_200m <- sum(Std_n2_200m, na.rm = TRUE)


#Std_n2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_n2_abundance <- data.frame(
  Nitrogen = c(sum_Std_n2_10m, sum_Std_n2_100m, sum_Std_n2_120m, sum_Std_n2_135m, sum_Std_n2_150m, sum_Std_n2_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_n2_dna <- ggplot(Std_n2_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_n2 concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Std_n2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)


#######################################################

#Mean_N2O in each depth (dna)
Mean_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
sum_Mean_N2O_10m <- sum(Mean_N2O_10m, na.rm = TRUE)


Mean_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2O) 
sum_Mean_N2O_100m <- sum(Mean_N2O_100m, na.rm = TRUE)

Mean_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2O) 
sum_Mean_N2O_120m <- sum(Mean_N2O_120m, na.rm = TRUE)

Mean_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2O) 
sum_Mean_N2O_135m <- sum(Mean_N2O_135m, na.rm = TRUE)

Mean_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2O) 
sum_Mean_N2O_150m <- sum(Mean_N2O_150m, na.rm = TRUE)

Mean_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
sum_Mean_N2O_200m <- sum(Mean_N2O_200m, na.rm = TRUE)


#SMean_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2O_abundance <- data.frame(
  Nitrogen = c(sum_Mean_N2O_10m, sum_Mean_N2O_100m, sum_Mean_N2O_120m, sum_Mean_N2O_135m, sum_Mean_N2O_150m, sum_Mean_N2O_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_N2O_dna <- ggplot(Mean_N2O_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Mean_N2O concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Mean_N2O concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) +
  geom_text(hjust=0, vjust=-2)


###############################################

#Std_N2O in each depth (dna)
Std_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
sum_Std_N2O_10m <- sum(Std_N2O_10m, na.rm = TRUE)


Std_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_N2O) 
sum_Std_N2O_100m <- sum(Std_N2O_100m, na.rm = TRUE)

Std_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_N2O) 
sum_Std_N2O_120m <- sum(Std_N2O_120m, na.rm = TRUE)

Std_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_N2O) 
sum_Std_N2O_135m <- sum(Std_N2O_135m, na.rm = TRUE)

Std_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_N2O) 
sum_Std_N2O_150m <- sum(Std_N2O_150m, na.rm = TRUE)

Std_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
sum_Std_N2O_200m <- sum(Std_N2O_200m, na.rm = TRUE)


#Std_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_N2O_abundance <- data.frame(
  Nitrogen = c(sum_Std_N2O_10m, sum_Std_N2O_100m, sum_Std_N2O_120m, sum_Std_N2O_135m, sum_Std_N2O_150m, sum_Std_N2O_200m), 
  Abundance = c(sum_abundance_10m, sum_abundance_100m, sum_abundance_120m, sum_abundance_135m, sum_abundance_150m, sum_abundance_200m),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_N2O_dna <- ggplot(Std_N2O_abundance, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA DNA abundance aginst Std_N2O concentration (uM)", 
       x = "NapA DNA Abundance",
       y = "Std_N2O concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)


###################################################

ggarrange(Std_NO2_dna, Mean_NH4_dna, Std_NH4_dna, Mean_NO2_dna, Mean_N2_dna, Std_n2_dna, Mean_N2O_dna, Std_N2O_dna, WS_NO3_dna,
          ncol = 3, nrow = 3,
          heights = c(0.1, 0.1, 0.1))





dna_abundance <- data.frame(
  "Species" = c("WS_NO3", "Mean_NH4", "Std_NH4", "Mean_NO2", "Std_NO2", "Mean_N2", "Std_n2", "Mean_N2O", "Std_N2O"),
  "Nitrogen Concentration (uM)" = c(sum_WS_NO3_10m, sum_Mean_NH4_10m, sum_Std_NH4_10m, sum_Mean_NO2_10m, sum_Std_NO2_10m, sum_Mean_N2_10m, sum_Std_n2_10m, sum_Mean_N2O_10m, sum_Std_N2O_10m)
  
)



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
sum_WS_NO3_10m <- sum(WS_NO3_10m, na.rm = TRUE)

sum_WS_NO3_10m

WS_NO3_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(WS_NO3) 
sum_WS_NO3_100m <- sum(WS_NO3_100m, na.rm = TRUE)

WS_NO3_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
sum_WS_NO3_120m <- sum(WS_NO3_120m, na.rm = TRUE)

WS_NO3_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(WS_NO3) 
sum_WS_NO3_135m <- sum(WS_NO3_135m, na.rm = TRUE)

WS_NO3_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(WS_NO3) 
sum_WS_NO3_150m <- sum(WS_NO3_150m, na.rm = TRUE)

WS_NO3_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(WS_NO3) 
sum_WS_NO3_200m <- sum(WS_NO3_200m, na.rm = TRUE)



#WS_NO3_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

WS_NO3_abundance_rna <- data.frame(
  Nitrogen = c(sum_WS_NO3_10m, sum_WS_NO3_100m, sum_WS_NO3_120m, sum_WS_NO3_135m, sum_WS_NO3_150m, sum_WS_NO3_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)



WS_NO3_rna <- ggplot(WS_NO3_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst WS_NO3 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "WS_NO3 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)

WS_NO3_rna
WS_NO3_abundance_rna


################################################



#Mean_NH4 in each depth (rna)
Mean_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NH4) 
sum_Mean_NH4_10m <- sum(Mean_NH4_10m, na.rm = TRUE)


Mean_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
sum_Mean_NH4_100m <- sum(Mean_NH4_100m, na.rm = TRUE)

Mean_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(WS_NO3) 
sum_Mean_NH4_120m <- sum(Mean_NH4_120m, na.rm = TRUE)

Mean_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NH4) 
sum_Mean_NH4_135m <- sum(Mean_NH4_135m, na.rm = TRUE)

Mean_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NH4) 
sum_Mean_NH4_150m <- sum(Mean_NH4_150m, na.rm = TRUE)

Mean_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NH4) 
sum_Mean_NH4_200m <- sum(Mean_NH4_200m, na.rm = TRUE)


#Mean_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NH4_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_NH4_10m, sum_Mean_NH4_100m, sum_Mean_NH4_120m, sum_Mean_NH4_135m, sum_Mean_NH4_150m, sum_Mean_NH4_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_NH4_rna <- ggplot(Mean_NH4_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_NH4 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Mean_NH4 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)

Mean_NH4_rna

#######################################################

#Std_NH4 in each depth (rna)
Std_NH4_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NH4) 
sum_Std_NH4_10m <- sum(Std_NH4_10m, na.rm = TRUE)


Std_NH4_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NH4) 
sum_Std_NH4_100m <- sum(Std_NH4_100m, na.rm = TRUE)

Std_NH4_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NH4) 
sum_Std_NH4_120m <- sum(Std_NH4_120m, na.rm = TRUE)

Std_NH4_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NH4) 
sum_Std_NH4_135m <- sum(Std_NH4_135m, na.rm = TRUE)

Std_NH4_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NH4) 
sum_Std_NH4_150m <- sum(Std_NH4_150m, na.rm = TRUE)

Std_NH4_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NH4) 
sum_Std_NH4_200m <- sum(Std_NH4_200m, na.rm = TRUE)


#Std_NH4_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NH4_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_NH4_10m, sum_Std_NH4_100m, sum_Std_NH4_120m, sum_Std_NH4_135m, sum_Std_NH4_150m, sum_Std_NH4_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_NH4_rna <-ggplot(Std_NH4_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_NH4 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Std_NH4 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)


###############################

#Mean_NO2 in each depth (rna)
Mean_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_NO2) 
sum_Mean_NO2_10m <- sum(Mean_NO2_10m, na.rm = TRUE)


Mean_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_NO2) 
sum_Mean_NO2_100m <- sum(Mean_NO2_100m, na.rm = TRUE)

Mean_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_NO2) 
sum_Mean_NO2_120m <- sum(Mean_NO2_120m, na.rm = TRUE)

Mean_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_NO2) 
sum_Mean_NO2_135m <- sum(Mean_NO2_135m, na.rm = TRUE)

Mean_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_NO2) 
sum_Mean_NO2_150m <- sum(Mean_NO2_150m, na.rm = TRUE)

Mean_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_NO2) 
sum_Mean_NO2_200m <- sum(Mean_NO2_200m, na.rm = TRUE)


#Mean_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_NO2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_NO2_10m, sum_Mean_NO2_100m, sum_Mean_NO2_120m, sum_Mean_NO2_135m, sum_Mean_NO2_150m, sum_Mean_NO2_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_NO2_rna <- ggplot(Mean_NO2_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_NO2 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Mean_NO2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)

############################################

#Std_NO2 in each depth (rna)
Std_NO2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_NO2) 
sum_Std_NO2_10m <- sum(Std_NO2_10m, na.rm = TRUE)


Std_NO2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_NO2) 
sum_Std_NO2_100m <- sum(Std_NO2_100m, na.rm = TRUE)

Std_NO2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_NO2) 
sum_Std_NO2_120m <- sum(Std_NO2_120m, na.rm = TRUE)

Std_NO2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_NO2) 
sum_Std_NO2_135m <- sum(Std_NO2_135m, na.rm = TRUE)

Std_NO2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_NO2) 
sum_Std_NO2_150m <- sum(Std_NO2_150m, na.rm = TRUE)

Std_NO2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_NO2) 
sum_Std_NO2_200m <- sum(Std_NO2_200m, na.rm = TRUE)


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_NO2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_NO2_10m, sum_Std_NO2_100m, sum_Std_NO2_120m, sum_Std_NO2_135m, sum_Std_NO2_150m, sum_Std_NO2_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_NO2_rna <- ggplot(Std_NO2_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_NO2 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Std_NO2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=1)

##################################################

#Mean_N2 in each depth (dna)
Mean_N2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2) 
sum_Mean_N2_10m <- sum(Mean_N2_10m, na.rm = TRUE)


Mean_N2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2) 
sum_Mean_N2_100m <- sum(Mean_N2_100m, na.rm = TRUE)

Mean_N2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2) 
sum_Mean_N2_120m <- sum(Mean_N2_120m, na.rm = TRUE)

Mean_N2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2) 
sum_Mean_N2_135m <- sum(Mean_N2_135m, na.rm = TRUE)

Mean_N2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2) 
sum_Mean_N2_150m <- sum(Mean_N2_150m, na.rm = TRUE)

Mean_N2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2) 
sum_Mean_N2_200m <- sum(Mean_N2_200m, na.rm = TRUE)


#Std_NO2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_N2_10m, sum_Mean_N2_100m, sum_Mean_N2_120m, sum_Mean_N2_135m, sum_Mean_N2_150m, sum_Mean_N2_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_N2_rna <- ggplot(Mean_N2_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_N2 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Mean_N2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)

#####################################################

#Std_n2 in each depth (rna)
Std_n2_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Std_n2) 
sum_Std_n2_10m <- sum(Std_n2_10m, na.rm = TRUE)


Std_n2_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_n2) 
sum_Std_n2_100m <- sum(Std_n2_100m, na.rm = TRUE)

Std_n2_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_n2) 
sum_Std_n2_120m <- sum(Std_n2_120m, na.rm = TRUE)

Std_n2_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_n2) 
sum_Std_n2_135m <- sum(Std_n2_135m, na.rm = TRUE)

Std_n2_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_n2) 
sum_Std_n2_150m <- sum(Std_n2_150m, na.rm = TRUE)

Std_n2_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_n2) 
sum_Std_n2_200m <- sum(Std_n2_200m, na.rm = TRUE)


#Std_n2_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Std_n2_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_n2_10m, sum_Std_n2_100m, sum_Std_n2_120m, sum_Std_n2_135m, sum_Std_n2_150m, sum_Std_n2_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_n2_rna <- ggplot(Std_n2_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_n2 concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Std_n2 concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)


#######################################################

#Mean_N2O in each depth (rna)
Mean_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
sum_Mean_N2O_10m <- sum(Mean_N2O_10m, na.rm = TRUE)


Mean_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Mean_N2O) 
sum_Mean_N2O_100m <- sum(Mean_N2O_100m, na.rm = TRUE)

Mean_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Mean_N2O) 
sum_Mean_N2O_120m <- sum(Mean_N2O_120m, na.rm = TRUE)

Mean_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Mean_N2O) 
sum_Mean_N2O_135m <- sum(Mean_N2O_135m, na.rm = TRUE)

Mean_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Mean_N2O) 
sum_Mean_N2O_150m <- sum(Mean_N2O_150m, na.rm = TRUE)

Mean_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Mean_N2O) 
sum_Mean_N2O_200m <- sum(Mean_N2O_200m, na.rm = TRUE)


#SMean_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 dna abundance data, therefore excluding this depth concentration

Mean_N2O_abundance_rna <- data.frame(
  Nitrogen = c(sum_Mean_N2O_10m, sum_Mean_N2O_100m, sum_Mean_N2O_120m, sum_Mean_N2O_135m, sum_Mean_N2O_150m, sum_Mean_N2O_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Mean_N2O_rna <- ggplot(Mean_N2O_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Mean_N2O concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Mean_N2O concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)


###############################################

#Std_N2O in each depth (rna)
Std_N2O_10m <- saanich_select %>%
  filter(Depth <= 0.010) %>%
  select(Mean_N2O) 
sum_Std_N2O_10m <- sum(Std_N2O_10m, na.rm = TRUE)


Std_N2O_100m <- saanich_select %>%
  filter(Depth <= 0.1) %>%
  filter(Depth > 0.010) %>%
  select(Std_N2O) 
sum_Std_N2O_100m <- sum(Std_N2O_100m, na.rm = TRUE)

Std_N2O_120m <- saanich_select %>%
  filter(Depth <= 0.12) %>%
  filter(Depth > 0.1)%>%
  select(Std_N2O) 
sum_Std_N2O_120m <- sum(Std_N2O_120m, na.rm = TRUE)

Std_N2O_135m <- saanich_select %>%
  filter(Depth <= 0.135) %>%
  filter(Depth > 0.12)%>%
  select(Std_N2O) 
sum_Std_N2O_135m <- sum(Std_N2O_135m, na.rm = TRUE)

Std_N2O_150m <- saanich_select %>%
  filter(Depth <= 0.15) %>%
  filter(Depth > 0.135)%>%
  select(Std_N2O) 
sum_Std_N2O_150m <- sum(Std_N2O_150m, na.rm = TRUE)

Std_N2O_200m <- saanich_select %>%
  filter(Depth <= 0.20) %>%
  filter(Depth > 0.15)%>%
  select(Std_N2O) 
sum_Std_N2O_200m <- sum(Std_N2O_200m, na.rm = TRUE)


#Std_N2O_above_200m <- saanich_select %>%
#  filter(Depth > 0.2)
#no above 200 rna abundance data, therefore excluding this depth concentration

Std_N2O_abundance_rna <- data.frame(
  Nitrogen = c(sum_Std_N2O_10m, sum_Std_N2O_100m, sum_Std_N2O_120m, sum_Std_N2O_135m, sum_Std_N2O_150m, sum_Std_N2O_200m), 
  Abundance = c(sum_10m_rna, sum_100m_rna, sum_120m_rna, sum_135m_rna, sum_150m_rna, sum_200m_rna),
  Depth = c("10m", "100m", "120m", "135m", "150m", "200m")
)


Std_N2O_rna <- ggplot(Std_N2O_abundance_rna, aes(y = Nitrogen, x = Abundance, color = Depth, label = Depth)) +
  geom_point(size = 4) +
  #  scale_y_reverse() +
  labs(title = "NapA RNA abundance aginst Std_N2O concentration (uM)", 
       x = "NapA RNA Abundance",
       y = "Std_N2O concentration (uM)", 
       color = "Depth") + 
  theme_minimal(base_size = 10) + 
  geom_text(hjust=0, vjust=-2)


###################################################

ggarrange(Std_NO2_rna, Mean_NH4_rna, Std_NH4_rna, Mean_NO2_rna, Mean_N2_rna, Std_n2_rna, Mean_N2O_rna, Std_N2O_rna, WS_NO3_rna,
          ncol = 3, nrow = 3,
          heights = c(0.1, 0.1, 0.1))


