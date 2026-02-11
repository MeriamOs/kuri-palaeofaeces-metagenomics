#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

###########################
### ENDOGENOUS DNA DATA ###
###########################

classified <- readr::read_csv(file = 'classified.csv')
head(classified)
# A tibble: 6 × 6
# Sample  Site `Host DNA %` `Microbial classified %` `Unclassified %`

classified_long <- classified %>%
  pivot_longer(cols = c("Host DNA %", "Microbial classified %",
                        "Unclassified %"),
               names_to = "fraction",
               values_to = "Value")

classified_long$fraction <- factor(
  classified_long$fraction, levels = 
    c("Unclassified %", "Microbial classified %", 
      "Host DNA %"))

desired_order <- c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                   "Whenua_Hou_post")

last_three_colors <- rev(tail(brewer.pal(12, "Set3"), 3))

ggplot(classified_long, aes(x = Sample, y = Value,
                            fill = fraction)) +
  geom_col() +
  facet_grid(~factor(Site, levels = desired_order), 
             scales = 'free') +
  labs(title = "Fractions of reads classified", 
       y = "Fraction") +
  scale_fill_manual(values = last_three_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggplot(classified, aes(x = Site, y = `Host DNA %`,
                       fill = Site)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Endogenous DNA percentage per site", 
       y = "ENdogenous DNA %") +
  scale_fill_brewer(palette = "Set1",
                    limits = c("Long_Bay", "Kahukura", 
                               "Whenua_Hou_pre", 
                               "Whenua_Hou_post"),
                    labels = c("Long Bay pre-contact (n=4)", 
                               "Kahukura pre-contact (n=4)", 
                               "Whenua Hou pre-contact (n=5)", 
                               "Whenua Hou post-contact (n=3)")) +
  theme(plot.title = element_text(size = 16, face = "bold", , hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom",  
        legend.direction = "horizontal") +
  # guides(shape = guide_legend(ncol = 1)) +
  guides(fill = guide_legend(ncol = 2)) 

#ggsave("Endogenous_DNA_percentage_per_site.png", width = 10, height = 6)

#############################
### MOLECULAR SEXING DATA ###
#############################

kuri_sex <- readr::read_csv(file = "kuri_sexing.csv")
head(kuri_sex)
# A tibble: 6 × 6
# Sample Site x_cov a_cov x_a_ratio Sex  

kuri_sex$Site <- factor(kuri_sex$Site, 
                        levels = c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                                   "Whenua_Hou_post"))

# Arrange by site only (not by x_a_ratio)
kuri_sex <- kuri_sex %>% 
  arrange(Site)

# Reverse sample order to display top-down in plot
kuri_sex$Sample <- factor(kuri_sex$Sample, 
                          levels = rev(kuri_sex$Sample))


ggplot(data = kuri_sex) + 
  geom_point(aes(x = x_a_ratio, 
                 y = Sample, 
                 color = as.character(Site),
                 shape = Sex),
             size = 6) + 
  labs(x = "X/A ratio", y = "Sample", 
       title = "Ratio of X chromosome / autosome coverage")  +
  labs(color = "Site", shape = "Sex") +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact", 
                                 "Kahukura pre-contact", 
                                 "Whenua Hou pre-contact", 
                                 "Whenua Hou post-contact")) +
  scale_shape_manual(values = c("female" = 17, "male" = 16),
                     labels = c("Female", "Male")) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16))

#ggsave("kuri_sexing_plot.png", width = 10, height = 7)

####################
### COMBINE DATA ###
####################

classified$Site <- factor(classified$Site, 
                        levels = c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                                   "Whenua_Hou_post"))

# Arrange by site only (not by x_a_ratio)
classified <- classified %>% 
  arrange(Site)

# Reverse sample order to display top-down in plot
classified$Sample <- factor(classified$Sample, 
                          levels = rev(classified$Sample))

classified_long <- classified %>%
  pivot_longer(cols = c("Host DNA %", "Microbial classified %",
                        "Unclassified %"),
               names_to = "fraction",
               values_to = "Value")

classified_long$fraction <- factor(
  classified_long$fraction, levels = 
    c("Unclassified %", "Microbial classified %", 
      "Host DNA %"))

p1 <- ggplot(classified_long, aes(y = Sample, x = Value, fill = fraction)) +
  geom_col() +
  facet_wrap(~factor(Site, levels = desired_order), ncol = 1, scales = "free_y") +
  labs(title = NULL, x = "Proportion of reads", 
       y = "Sample") +
  labs(fill = "Reads classification") +
  scale_fill_manual(values = last_three_colors) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    strip.text = element_blank(),
    strip.background = element_blank())

p2 <- ggplot(data = kuri_sex) + 
  geom_col(aes(x = a_cov, 
                 y = Sample, 
                 fill = as.character(Site)),
             size = 6) + 
  facet_wrap(~factor(Site, levels = desired_order), ncol = 1, scales = "free_y") +
  labs(x = "Genome coverage", y = "Sample")  +
#  labs(fill = "Site") +
  scale_fill_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact", 
                                 "Kahukura pre-contact", 
                                 "Whenua Hou pre-contact", 
                                 "Whenua Hou post-contact")) +
  guides(fill = "none") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank())

p3 <- ggplot(data = kuri_sex) + 
  geom_point(aes(x = x_a_ratio, 
                 y = Sample, 
                 color = as.character(Site),
                 shape = Sex),
             size = 6) + 
  facet_wrap(~factor(Site, levels = desired_order), ncol = 1, scales = "free_y") +
  labs(x = "X/A coverage ratio", y = "Sample", 
       title = "Combined plot of endogenous DNA and sex")  +
  labs(color = "Site", shape = "Sex") +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact", 
                                 "Kahukura pre-contact", 
                                 "Whenua Hou pre-contact", 
                                 "Whenua Hou post-contact")) +
  scale_shape_manual(values = c("female" = 17, "male" = 16),
                     labels = c("Female", "Male")) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank())

library(patchwork)

combined <- p1 + p2 + p3 +
  plot_layout(ncol = 3, widths = c(1, 1), guides = "collect") & 
  theme(legend.position = "right")
combined

ggsave("Combined_endogenous_DNA_and_sex.png", width = 15, height = 6)
