## compare <50bp and >50bp metaphlan profiles ##
#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025
#Edited in July 2026 M. van Os

library(phyloseq)
library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(plotly)
library(dplyr)
library(mixOmics)
library(microbiome)
library(forcats)
library(scales)

setwd("~/Documents/Rstudio/r-tidyverse")

###########################
### LOAD METAPHLAN DATA ### 
###########################

mpa_pal <- read.table("merged_palaeofaeces_metaphlan_table_species.txt", header = TRUE) %>%
   dplyr::rename_with(~str_remove(., '_jun23')) %>% 
   dplyr::rename_with(~str_remove(., '.metaphlan')) %>% 
   dplyr::rename_with(~str_remove(., '_nreads')) %>% 
   rename("ID" = colnames(.)[1]) %>%
   mutate(ID = gsub("^.*s__|","", ID))

### Combine datasets
mpa_all <- mpa_pal

### transpose dataset
mpaT <- t(mpa_all[-1]) %>%
  as.data.frame() 
colnames(mpaT) <- mpa_all$ID 
#mpaT <- rownames_to_column(mpaT, var = "Sample") 

mpaT[is.na(mpaT)] <- 0.0 #set NAs to 0

##################################################
### filter dataset to remove low abundant taxa ###
##################################################

count_cut_off <- 0.5 #percentage in metaphlan, vs fraction (out of 1.0) in kraken2

mpa_filter <- mpaT[, apply(mpaT, 2, function(col) any(col >= count_cut_off))]

### Centered Log Ratio transformation
CLR <- microbiome::transform(mpa_filter, "clr", pseudocount = 1)

dim(CLR) #check is dimensions are correct
#[1] 16 (samples) X (taxa)

##################
### Scree plot ###
##################

tune_pca <- tune.pca(CLR, ncomp = 10, scale = TRUE)
variation <- as.data.frame(tune_pca$prop_expl_var) %>%
  rownames_to_column()
variation$PC <- rownames(variation)
variation$cum <- tune_pca$cum.var

ggplot(variation, aes(x = fct_inorder(PC), y = X)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal() 

## Check for elbow in explained proportions of principal components
## > elbow occurs at 3/4 components

ggsave("Scree plot - variance explained principal components.png", 
       width = 12, height = 8)

###############################################
### PCA analysis on CLR transformed dataset ###
###############################################

## set number of components (ncomp) to findings from scree plot
pca_result <- mixOmics::pca(CLR, ncomp = 2, scale = TRUE, center = TRUE)
plotIndiv(pca_result, comp = c(1, 2))

pca_scores <- as.data.frame(pca_result$variates$X)  # Get scores for the components
pca_scores$Sample <- rownames(pca_scores)  # Add sample names as a column

site_info <- read.csv("sample_site_info.csv")
site_info <- site_info[103:118,]

## Merge PCA scores with site information
pca_plot_data <- merge(pca_scores, site_info, by = "Sample")

## Create the PCA plot
#scatter <- 
ggplot(data = pca_plot_data, 
       aes(x = PC1, y = PC2, colour = Site, 
           shape = Period)) +
  geom_point(size = 5, alpha = 0.7) + #size = 5, alpha = 0.7
  #  stat_ellipse(aes(group = Island), level = 0.95) +
  geom_line(data = subset(pca_plot_data, grepl("^MS", sample)), 
              aes(group = sample), alpha = 0.5) +
  labs(x = "PC1 (18%)", #NOTE: put in PC values
       y = "PC2 (16%)", #NOTE: put in PC values
       title = "PCA with CLR transformation", size = 16) +
   scale_color_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"
                                 #"modern_dog", "arch_bone", 
                                 #"blank", "modern_human"
                                 ),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)"
                                 # "Modern dog", "Archaeological bone", 
                                 # "Blank", "Modern Human"
                                 ),
                      guide = guide_legend(nrow = 2)) +
  scale_shape_manual(values = c(17, 16), 
                     labels = c("post-contact", "pre-contact"),
                     guide = guide_legend(nrow = 2)) +
    # scale_shape_manual(values = c(16, 17, 15), 
    #                   labels = c("All", "Over 50 bp", "Under 50 bp"),
    #                   guide = guide_legend(nrow = 3)) +
  theme(
    plot.title = element_blank(),
    #    plot.title = element_text(size = 18, 
    #                              face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"), 
    legend.text = element_text(size = 16),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    axis.title.x = element_text(size = 20, face = "bold"),  # Change x-axis title size
    axis.title.y = element_text(size = 20, face = "bold")) 

#ggplotly(scatter)

#ggsave("PCA_CLR_kuri_sp_cutoff_0.1percent_metaphlan_period.png", 
#       width = 10, height = 8)  
  
ggsave("pub_PCA_CLR_palaeofaeces_cutoff_0.5percent_metaphlan.png", 
       width = 5.5, height = 5)

##########################
### PERMANOVA ANALYSIS ###
##########################

site_info <- read.csv("sample_site_info.csv", header=TRUE)
site_info <- site_info[1:16,]

rownames(mpa_filter) <- stringr::str_remove(rownames(mpa_filter), "_under50bp")

CLR <- microbiome::transform(mpa_filter, "clr", pseudocount = 1)

### test for the different variables by changing the variable after ~
## Variables tested: Island, Site, Period, Sex, Relatedness
permanova <- vegan::adonis2(CLR ~ Island,
                            data = site_info,
                            method = "euclidean",
                            permutations = 9999)

permanova # view results

#################################
######## DESeq2 ANALYSIS ########
#################################

## DESeq2 only works on count data, so tested otu_num and otu_frac turned into "count data"
## Analysis in thesis based on the otu_frac approach
mpa_deseq <- mpa_filter * 10000
#otu_deseq <- otu_num

dds <- DESeqDataSetFromMatrix(countData = round(t(mpa_deseq)), 
                              colData = site_info, 
                              design = ~ Island)

dds <- DESeq(dds)

# Test for various variables by changing the data within c()
res <- results(dds, contrast = c("Island", 
                                 "North_Island", 
                                 "South_Island")) %>%
  as.data.frame()

# View the results for the comparison
head(res)

# Change taxa name to column
res_tax <- rownames_to_column(res) %>%
  as.data.frame()
# colnames(res_tax)[1] <- "taxonomy_id"

## Create interactive vulcano plot to visualise taxa drivers
vulcano <- ggplot(res_tax, aes(x = log2FoldChange, y = -log10(padj),
                               text = paste("Taxid:", rowname))) +
  geom_point(aes(color = padj < 0.01), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot of Differentially Abundant Taxa North vs South", # NOTE change title
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value")

vulcano_interactive <- ggplotly(vulcano)
vulcano_interactive

#htmlwidgets::saveWidget(vulcano_interactive, "interactive_vulcano_plot_t0.75_c0.001_pre_vs_post-contact.html")

## Extract significant taxa and output as CSV file
significant_taxa <- subset(res, padj < 0.05) %>%
  as.data.frame()

head(significant_taxa) # View the significant taxa

significant_taxid <- rownames_to_column(significant_taxa) %>%
  as.data.frame()
colnames(significant_taxid)[1] <- "taxonomy_id"

write.csv(significant_taxid, file = "significant_taxa_metaphlan_north_vs_south_cutoff0.5.csv", row.names = FALSE)

###################################
######## MaAslin2 ANALYSIS ########
###################################

#colnames(CLR) <- samples$name[match(colnames(CLR), samples$taxonomy_id)]

maaslin_metadata <- column_to_rownames(site_info, var = "sample")

## Test for various variables by changing the data 
## Normalisation and transformation set to NONE as using CLR data
fit_data = Maaslin2(input_data     = CLR, 
                    input_metadata = maaslin_metadata, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    transform      = "NONE",
                    output         = "Maaslin_Sex_CLR_output", 
                    fixed_effects  = c("Period"),
#                    reference      = c("Sex,male")
                    heatmap_first_n = 100,
                    max_pngs        = 150)
