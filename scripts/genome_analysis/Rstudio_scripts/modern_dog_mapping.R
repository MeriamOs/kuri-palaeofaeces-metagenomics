#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridisLite)
library(ggnewscale)

setwd("~/Documents/Rstudio/r-tidyverse")

##########################
### MODERN DOG MAPPING ###
##########################

modern_dog <- readr::read_csv("modern_dog_info.csv")
head(modern_dog)
# A tibble: 6 × 16
# Sample Country read_count `Collapsed reads` host_mapped host_DNA

sorted <- modern_dog[order(modern_dog$host_DNA), ]
sorted$Sample <- factor(sorted$Sample, levels = sorted$Sample)

### Bin into 5% bins
modern_dog$host_DNA_num_bin <- cut(modern_dog$host_DNA_num, 
 breaks = seq(0, 100, by = 5), right = FALSE, include.lowest = TRUE)

# Count the number of samples in each bin
count_df <- modern_dog %>%
  group_by(host_DNA_num_bin) %>%
  summarise(count = n())

ggplot(count_df, aes(x = host_DNA_num_bin, y = count)) +
  geom_bar(stat = "identity", fill = "darkmagenta") +
  labs(x = "host_DNA_num (Binned)", y = "Sample Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Endogenous DNA percentage (binned)", y = "Sample count")

### Count samples per bin and per country
count_df <- modern_dog %>%
  group_by(Country, host_DNA_num_bin) %>%
  summarise(count = n(), .groups = "drop")

ggplot(count_df, aes(x = host_DNA_num_bin, y = count, fill = Country)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Country) +
  labs(
    x = "Endogenous DNA percentage (binned)",
    y = "Sample count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave("Binned_endogenous_DNA_percentage_per_country.png", 
#              width = 12, height = 6)

ad_levels <- unique(sorted$AD)
country_levels <- unique(sorted$Country)

ad_colors <- setNames(brewer.pal(length(ad_levels), "Paired"), ad_levels)
country_colors <- setNames(brewer.pal(length(country_levels), "Accent"), country_levels)

# Make sure Sample is a factor (helps with x-axis alignment)
sorted$Sample <- factor(sorted$Sample, levels = sorted$Sample)

ggplot(sorted, aes(x = Sample)) +
  geom_bar(aes(y = host_DNA, fill = AD), stat = "identity") +
  scale_fill_manual(name = "AD", values = ad_colors) +
  new_scale_fill() +
  geom_tile(
    aes(y = -max(sorted$host_DNA) * 0.05, fill = Country),
    height = max(sorted$host_DNA) * 0.03) +
  scale_fill_manual(name = "Country", values = country_colors) +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Host DNA")

ggsave("Endogenous_DNA_percentage_barplot_coloured_lifestyle_country.png", 
       width = 8, height = 5)

### Sort per country AND lifestyle
sorted$AD <- factor(sorted$AD, levels = c(
                    "Diet_experiment", "Shelter", "Street", "Laos Dog",
                    "Farm", "Suburban", "Urban", "missing"), 
                    labels = c("USA Dog","India, Shelter",
                    "India, Street", "Laos Dog","SA, Farm",
                    "SA, Suburban", "SA, Urban", "Not available"))

ggplot(sorted, aes(x = AD, y = host_DNA, fill = AD)) +
  geom_boxplot() +
  geom_point() +
  scale_fill_brewer(palette = "Set3",
                    labels = c("USA Dog (n = 10)",
                               "India, Shelter (n = 20)",
                               "India, Street (n = 14)",
                               "Laos Dog (n = 33)",
                               "South Africa, Farm (n = 3)",
                               "South Africa, Suburban (n = 4)",
                               "South Africa, Urban (n = 12)",
                               "Not available (n = 11)")) +
  labs(x = "Dog category", 
       y = "Endogenous DNA percentage", 
       title = "Boxplots of endogenous DNA percentage per dog lifestyle category", size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18),
    legend.position = "bottom",  
    legend.direction = "horizontal") 

#ggsave("Endogenous_DNA_percentage_boxplots_per_lifestyle_ppp.png", 
#       width = 14, height = 8)

####################################
### COMBINE MODERN AND KURI DATA ###
####################################

modern <- dplyr::select(modern_dog, Sample, host_DNA, AD)

kuri <- read_csv("MR_kuri.csv")
head(kuri)
# A tibble: 6 × 3
# Sample  host_DNA AD 

merged <- rbind(modern, kuri)

merged$AD <- factor(merged$AD, 
                    levels = c(
                      "Diet_experiment", "Shelter", "Street", "Laos Dog",
                      "Farm", "Suburban", "Urban", "missing", "palaeofaeces"), 
                    labels = c(
                      "USA Dog", "India, Shelter", "India, Street",
                      "Laos Dog", "SA, Farm", "SA, Suburban",
                      "SA, Urban", "Not available", "Palaeofaeces"))

ggplot(merged, aes(x = AD, y = host_DNA, fill = AD)) +
  geom_boxplot() +
  geom_point() +
  scale_fill_brewer(palette = "Set3",
                    labels = c("USA Dog (n = 10)",
                               "India, Shelter (n = 20)",
                               "India, Street (n = 14)",
                               "Laos Dog (n = 33)",
                               "South Africa, Farm (n = 3)",
                               "South Africa, Suburban (n = 4)",
                               "South Africa, Urban (n = 12)",
                               "Not available (n = 11)",
                               "Palaeofaeces (n = 30)"),
                    guide = guide_legend(nrow = 3)) +
  labs(x = "Dog category", 
       y = "Endogenous DNA percentage", 
       title = "Boxplots of endogenous DNA percentage per dog lifestyle category", size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size = 22),
    axis.title.x = element_text(size = 22, face = "bold"),  
    axis.title.y = element_text(size = 22, face = "bold"),
    legend.position = "bottom",  
    legend.direction = "horizontal")  

#ggsave("Endogenous_DNA_percentage_boxplots_modern_and_kuri_ppp.png", 
#       width = 14, height = 10)

#anova test
anova_result <- aov(`host_DNA` ~ AD, data = merged)
summary(anova_result)

TukeyHSD(anova_result) # Honestly Significant Difference

######################################
### HOST MAPPING Cui et al. (2024) ###
######################################

cui_mapping <- readr::read_csv("Host_mapping_Cui_et_al.csv")
head(cui_mapping)
# A tibble: 6 × 14
# Sample Species `Common name` Diet `CE (ng/ul)` `RR (×100 million)`
# `Fd (%)` `MR (%)` `Mitochondrial genome Dm` `Mitochondrial genome Cm (%)` 
# `Nuclear genome Dn` `Nuclear genome Cn (%)` 
# `Mitochondrial genome Accession ID` `Nuclear genome Accession ID` 

### Group per Dietary class
cui_mapping$Diet <- factor(cui_mapping$Diet, levels = c("Herbivore", "Omnivore", "Carnivore"))

ggplot(cui_mapping, aes(x = `Common name`, y = `MR (%)`, fill = `Common name`)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ Diet, scales = "free_x") +
  scale_fill_brewer(palette = "Set3",
                      labels = c("Amur leopard (n = 1)",
                                 "Amur tiger (n = 14)",
                                 "Artic fox (n = 3)",
                                 "Cow (n = 4)",
                                 "Leopard cat (n = 1)",
                                 "Moose (n = 2)",
                                 "Red deer (n = 15)",
                                 "Roe deer (n = 1)",
                                 "Sika deer (n = 4)",
                                 "Wild Boar (n = 3)")) +
  guides(fill = guide_legend(nrow = 3)) +
  labs(x = "Species", 
       y = "Endogenous DNA percentage", 
       title = "Boxplots of endogenous DNA percentage per species", size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),  
    axis.title.y = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "bottom",  
    legend.direction = "horizontal")  

#ggsave("Endogenous_DNA_percentage_boxplots_per_species_Cui_et_al.png", 
#       width = 12, height = 8)

count_cui <- cui_mapping %>%
  group_by(Diet) %>%
  summarise(count = n(), .groups = "drop")

ggplot(cui_mapping, aes(x = Diet, y = `MR (%)`, fill = Diet)) +
  geom_boxplot() +
  geom_point() +
#  facet_wrap(~ Diet, scales = "free_x") +
  scale_fill_brewer(palette = "Set3",
                    labels = c("Herbivore (n = 26)",
                               "Omnivore (n = 3)",
                               "Carnivore (n = 19)")) +
  labs(x = "Species", 
       y = "Endogenous DNA percentage", 
       title = "Boxplots of endogenous DNA percentage per diet", size = 16) +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Change x-axis title size
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",  
    legend.direction = "horizontal") 

#ggsave("Endogenous_DNA_percentage_boxplots_per_diet_Cui_et_al.png", 
#       width = 8, height = 6)

#################################
### COMBINE CUI AND KURI DATA ###
#################################

cui <- dplyr::select(cui_mapping, Sample, Diet, `MR (%)`)

kuri <- read_csv("MR_kuri.csv")

colnames(kuri)[2] <- "MR (%)"
colnames(kuri)[3] <- "Diet"

kuri <- kuri[, c(1, 3, 2)]
merged <- rbind(cui, kuri)

merged %>%
  group_by(Diet) %>%
  summarise(
    mean_MR = mean(`MR (%)`, na.rm = TRUE),
    median_MR = median(`MR (%)`, na.rm = TRUE),
    sd_MR = sd(`MR (%)`, na.rm = TRUE),
    min_MR = min(`MR (%)`, na.rm = TRUE),
    max_MR = max(`MR (%)`, na.rm = TRUE),
    n = n())

merged$Diet <- factor(merged$Diet, levels = c("Herbivore", "Omnivore", "Carnivore", "palaeofaeces"))

ggplot(merged, aes(x = Diet, y = `MR (%)`, fill = Diet)) +
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(
    values = c("Herbivore" = "#66c2a5",
               "Omnivore" = "#fc8d62",
               "Carnivore" = "#8da0cb",
               "palaeofaeces" = "gray90"),
    labels = c("Herbivore (n = 26)",
               "Omnivore (n = 3)",
               "Carnivore (n = 19)",
               "Palaeofaeces (n = 30)")) +
  labs(x = "Dietary class", 
       y = "Endogenous DNA percentage", 
       title = "Boxplots of endogenous DNA percentage per diet classification", size = 16) +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Change x-axis title size
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",  
    legend.direction = "horizontal") 

ggsave("Endogenous_DNA_percentage_boxplots_per_diet_Cui_and_palaeofaeces_ppp.png", 
       width = 10, height = 6)

######################
#### ANOVA TESTING ###
######################

# data = cui_mapping OR data = merged
anova_result <- aov(`MR (%)` ~ Diet, data = merged) 
summary(anova_result)

TukeyHSD(anova_result) # Honestly Significant Difference


