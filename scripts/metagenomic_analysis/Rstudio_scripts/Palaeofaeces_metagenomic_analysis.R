#Script written by Meriam van Os
#Used for downstream analysis of metagenomic shotgun data from kuri (dog) palaeofaeces
#Uploaded 18/09/2025

library(decontam)
library(vegan)
library(DESeq2)
library(readr)
library(dplyr)
library(glue)
library(stringr)
library(phyloseq)
library(microbiome)
library(mixOmics)
library(ggplot2)
library(plotly)
library(scales)
library(forcats)
library(ecodist)
library(Maaslin2)
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

#setwd("~/Documents/Rstudio/r-tidyverse")

############################
### GENERAL KRAKEN STATS ###
############################

stats <- readr::read_csv(file = 'kraken-standard-2024-overall-stats.csv')

site_info <- read.csv("sample_site_info.csv", header=TRUE)
site_info <- site_info[1:16,]

stats_kraken <- merge(stats, site_info, by = "Sample")

colnames(stats_kraken) <- gsub("_", " ", colnames(stats_kraken))

stats_kraken$"Classified reads" <- stats_kraken$"Classified reads" * 100

stats_long <- stats_kraken %>%
  pivot_longer(cols = c("Number of raw reads", "Classified reads"), 
               names_to = "Reads", 
               values_to = "Value")

stats_long$Site <- factor(stats_long$Site, 
                          levels = c("Long_Bay", "Kahukura", 
                                     "Whenua_Hou_pre", 
                                     "Whenua_Hou_post"))

stats_long$Reads <- factor(stats_long$Reads, 
                           levels = c("Number of raw reads", "Classified reads"))

combined_plot <- ggplot(stats_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  scale_y_continuous(labels = comma) +
  labs(title = "Boxplots of number of non-host reads per site", 
       x = "Reads", 
       y = "Value") +
  facet_wrap(~Reads, scales = "free",
             labeller = as_labeller(c(
               "Number of raw reads" = "A. Number of raw reads",
               "Classified reads" = "B. Classified reads"
             ))) +
  scale_fill_brewer(palette = "Set1",
                    limits = c("Long_Bay", "Kahukura", 
                               "Whenua_Hou_pre", 
                               "Whenua_Hou_post"),
                    labels = c("Long Bay pre-contact (n=4)", 
                               "Kahukura pre-contact (n=4)", 
                               "Whenua Hou pre-contact (n=5)", 
                               "Whenua Hou post-contact (n=3)")) +
  theme(axis.text.x = element_blank(), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 20,, face = "bold"),  # Change facet label size
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),  # Change x-axis title size
        axis.title.y = element_text(size = 20),  # Change y-axis title size
        axis.text.y = element_text(size = 18),
        legend.position = "bottom",  # Move the legend to the bottom
        legend.direction = "horizontal") +
  guides(fill = guide_legend(ncol = 2, title = NULL))

combined_plot
ggsave("boxplots_no_non-host_reads_&_percentage_of_kraken2_classified_reads.png", 
       combined_plot, width = 12, height = 6)

##########################
### RAREFACTION CURVES ###
##########################

samples <- readr::read_csv(file = 'palaeofaeces_combined_2024-standard_minimizer_bracken_species.csv')

#Decontaminate count dataframe
samples_num <- dplyr::select(samples,
                             name, taxonomy_id,
                             contains("num")) %>%
  dplyr::rename_with(~str_remove(., '_filtered.bracken_num')) %>%
  dplyr::rename_with(~str_remove(., '_minimizer'))

samples_numT <- t(samples_num[ ,-c(1:2)]) %>%
  as.data.frame()
colnames(samples_numT) <- samples$taxonomy_id

sample_size <- rowSums(samples_numT)
raremin <- min(rowSums(samples_numT)) #smallest sample size across samples
Srare <- rarefy(samples_numT, raremin)

site_levels <- c("Long_Bay", "Kahukura", "Whenua_Hou_pre", "Whenua_Hou_post")
site_cols   <- setNames(brewer.pal(length(site_levels), "Set1"), site_levels)

col_vector  <- site_cols[ site_info$Site[ match(rownames(samples_numT), site_info$Sample) ] ]

png("rarefaction_curves_kraken2_standard.png", width = 8, height = 6, units = "in", res = 300) 
rarecurve(samples_numT, step = 1000, col = col_vector, cex = 0.8, 
          sample= raremin, xlab = "Sample Size", ylab = "Species", 
          cex.lab = 1.5, cex.axis = 1.3,
          label = TRUE, lwd = 2)

legend("bottomright", 
       legend = c("Long Bay pre-contact (n=4)", 
                  "Kahukura pre-contact (n=4)", 
                  "Whenua Hou pre-contact (n=5)", 
                  "Whenua Hou post-contact (n=3)"), 
       col = brewer.pal(4, "Set1"), 
       lty = 1, lwd = 2, cex = 1.1, bty = "n", 
       text.font = 2, text.siz)

dev.off()


###########################################
### IDENTIFY CONTAMINANTS WITH DECONTAM ###
###########################################

#read in combined sample bracken table and select the fraction columns (rather than read columns)
samples <- readr::read_csv(file = 'palaeofaeces_combined_2024-standard_minimizer_bracken_species.csv')
samples_frac <- dplyr::select(samples, 
                        name, taxonomy_id, 
                        contains("frac")) %>% 
  dplyr::rename_with(~str_remove(., '_filtered.bracken_frac')) %>%
  dplyr::rename_with(~str_remove(., '_minimizer'))

#read in combined environmental control bracken table and select the fraction columns (rather than read columns)
blanks <- readr::read_csv(file = 'negatives_combined_2024-standard_minimizer_bracken_species.csv')
blanks <- dplyr::select(blanks, 
                         name, taxonomy_id, 
                         contains("frac")) %>% 
  dplyr::rename_with(~str_remove(., '_filtered.bracken_frac')) %>%
  dplyr::rename_with(~str_remove(., '_minimizer'))

#merge the two data frames
merged <- merge(samples_frac, blanks, by = "taxonomy_id", all = TRUE)

merged[is.na(merged)] <- 0 # Replace NA with 0

otu <- merged[, !(colnames(merged) %in% c("name.x", "name.y"))]

mergedT <- t(otu[, -1])
colnames(mergedT) <- merged$taxonomy_id %>%
  as.matrix()

#turn data into a phyloseq data frame
otu_table_ps <- otu_table(mergedT, taxa_are_rows=FALSE)  

decontam_info <- readr::read_csv(file = 'decontam_info.csv')
metadata <- sample_data(decontam_info)
rownames(metadata) <- decontam_info$Sample 

ps <- phyloseq(otu_table_ps, metadata)

#Identify contaminants with decontam using the prevalance method, code adjusted from 
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#introduction ###

sample_data(ps)$is.neg <- sample_data(ps)$decontam == "TRUE"

#Data is already in fractions, so ormalisation is turned off
#Tried with normalisation = TRUE as well, but same results
contamdf.prev <- isContaminant(ps, 
            method="prevalence", neg="is.neg", threshold=0.75, normalize = FALSE)
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$decontam == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$decontam == "FALSE", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, 
  color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#interactive_plot <- ggplotly(p = ggplot2::last_plot(), tooltip = "text")
#interactive_plot

row_decontam <- which(contamdf.prev$contaminant)

contaminant_info <- merged[row_decontam, c("taxonomy_id", "name.y")]

write.csv(contaminant_info, file = "contaminant_species_t0.75_kraken_2024-standard.csv", row.names = FALSE)

#################################################
### REMOVE CONTAMINANT AND LOW ABUNDANCE TAXA ###
#################################################

contaminant_info <- readr::read_csv(file = 'contaminant_species_t0.75_kraken_2024-standard.csv')
#contaminant_info <- readr::read_csv(file = 'contaminant_species_t0.75_kraken_2024-standard_extra_added.csv')
#contaminant_info <- readr::read_csv(file = 'contaminant_species_t0.75_kraken_2024-standard_extra_fish_added.csv')

taxids <- contaminant_info$taxonomy_id

#Decontaminate fraction dataframe
samples_fracT <- t(samples_frac[ ,-c(1:2)]) %>%
  as.data.frame()
colnames(samples_fracT) <- samples$taxonomy_id

otu_frac_clean <- samples_fracT[ , !(colnames(samples_fracT) %in% taxids)] %>%
  as.data.frame()

#colnames(otu_frac_clean) <- samples$name[match(colnames(otu_frac_clean), samples$taxonomy_id)]
#otu_frac_cleanT <- t(otu_frac_clean)

frac_cut_off <- 0.001
otu_frac <- otu_frac_clean[, apply(otu_frac_clean, 2, function(col) any(col >= frac_cut_off))]

#Decontaminate count dataframe
samples_num <- dplyr::select(samples,
                         name, taxonomy_id,
                         contains("num")) %>%
  dplyr::rename_with(~str_remove(., '_filtered.bracken_num')) %>%
  dplyr::rename_with(~str_remove(., '_minimizer'))

samples_numT <- t(samples_num[ ,-c(1:2)]) %>%
  as.data.frame()
colnames(samples_numT) <- samples$taxonomy_id

otu_num_clean <- samples_numT[ , !(colnames(samples_numT) %in% taxids)] %>%
  as.data.frame()

#colnames(otu_num_clean) <- samples$name[match(colnames(otu_num_clean), samples$taxonomy_id)]
#otu_num_cleanT <- t(otu_num_clean)

count_cut_off <- 1000
otu_num <- otu_num_clean[, apply(otu_num_clean, 2, function(col) any(col >= count_cut_off))]

#otu_frac <- sweep(otu_num, 1, rowSums(otu_num), FUN = "/")

#################################
######## ALPHA DIVERSITY ########
#################################

otu_alpha <- otu_frac
#otu_alpha <- otu_num

mean(colMeans(otu_alpha))

### Alpha diversity per sample and site ###
shannon <- data.frame(vegan::diversity(otu_alpha, "shannon")) 
names(shannon)[1] <- "shannon" 
simpson <- data.frame(vegan::diversity(otu_alpha, "simpson"))
names(simpson)[1] <- "simpson" 
invsimp <- data.frame(vegan::diversity(otu_alpha, "inv"))
names(invsimp)[1] <- "inv_simpson" 
#fisher <- data.frame(fisher.alpha(otu_filter))
#names(fisher)[1] <- "fisher" 
readcount <- data.frame(rowSums(otu_alpha))
names(readcount)[1] <- "readcount"

## Species richness and Pielou's evenness ##
richness <- data.frame(specnumber(otu_alpha))
names(richness)[1] <- "Richness" 
evenness <- data.frame(shannon/log(specnumber(otu_alpha)))
names(evenness)[1] <- "evenness"


Diversity <- cbind(shannon, simpson, invsimp, 
                   richness, evenness, 
                   readcount)
Diversity$Sample <- row.names(Diversity)
names(Diversity) <- c("Shannon", "Simpson", "InvSimpson", 
                      "Richness", "Evenness", "Proportion of reads", "Sample")
Diversity <- Diversity[,c(7,1,2,3,4,5,6)]

# Add site info to dataframe 
Diversity_site_data <- merge(Diversity, site_info, by = "Sample")

## Create boxplots alpha diversity ##
div_long <- Diversity_site_data %>%
  pivot_longer(cols = c("Shannon", "Simpson", 
                        "Proportion of reads", "InvSimpson", 
                        "Richness", "Evenness"), 
               names_to = "Index", 
               values_to = "Value")

# Reorder the Site and index variable
div_long$Site <- factor(div_long$Site, 
                        levels = c("Long_Bay", "Kahukura", "Whenua_Hou_pre",
                                   "Whenua_Hou_post"))

div_long$Index <- factor(div_long$Index, 
                         levels = c("Shannon", "Simpson", "InvSimpson", 
                                    "Evenness", "Richness", "Proportion of reads"))  

ggplot(div_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  labs(title = glue("Boxplots of alpha index values per site after removing taxa <{frac_cut_off}"), 
       x = "Index", 
       y = "Value") +
  facet_wrap('~Index', scales = 'free') +
  theme(axis.text.x = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 16),  # Change facet label size
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),  # Change x-axis title size
        axis.title.y = element_text(size = 16),  # Change y-axis title size
        axis.text.y = element_text(size = 14),
        legend.position = "bottom",  # Move the legend to the bottom
        legend.direction = "horizontal") +  # Change y-axis tick labels size) 
  scale_fill_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)")) +
  guides(fill = guide_legend(ncol = 2, title = NULL))
#  scale_y_continuous(labels = scales::number_format())

ggsave("boxplots_alpha_diversity_t0.75_0.00001fracfilter_per_site.png", 
       width = 12, height = 8)

# Create scatterplot matrix
# png("scatterplot_alpha_diversity.png", width = 1200, height = 800)
# pairs(cbind(simpson, shannon, richness, 
#             evenness, readcount), pch="+", col="blue")
# dev.off()

# Add site info to dataframe 
stats_kraken$"Classified_reads" <- stats_kraken$"Classified_reads" * 100

stats_kraken$Site <- factor(stats_kraken$Site, 
                            levels = c("Long_Bay", "Kahukura", 
                                       "Whenua_Hou_pre", 
                                       "Whenua_Hou_post"))


Diversity_site_data$Number_of_raw_reads <- stats_kraken$Number_of_raw_reads
#Diversity$Site <- stats_kraken$Site

Species_long <- Diversity_site_data %>%
  pivot_longer(cols = c("Shannon", "Simpson", "InvSimpson", 
                        "Evenness", "Richness"),
               names_to = "metric") %>%
  as.data.frame()
#Species_long <- merge(Species_long, site_info, by = "Sample")

ggplot(Species_long, aes(x = Number_of_raw_reads, y = value)) +
  geom_point(aes(colour = Site)) +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)")) +
  geom_smooth() +
  facet_wrap(~metric, nrow = 6, scales = "free_y") +
  scale_x_continuous(labels = scales::comma)

ggsave("relation_alpha_diversity_number_of_non-host_reads_nofilter_kraken2024-standard.png", 
       width = 9, height = 12)

### per richness plots ###
ggplot(data = div_long, aes(x=Sample, y=Value, fill = Index)) +
  geom_col() +
  facet_wrap(~Index, scales = "free_y") +
  scale_fill_viridis_d(option = "D") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("alpha_diversity_0.001readfilt_per_sample.png", 
       width = 16, height = 12)


#########################################
######## BETA DIVERSITY ANALYSIS ########
#########################################

# CLR transformation with pseudocount added as data contains zeros
CLR <- microbiome::transform(otu_frac, "clr", pseudocount = 0)

dim(CLR) #[1] 16 (samples) X (taxa)

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
#elbow occurs at 3-4

#ggsave("Scree plot - variance explained principal components.png", 
#       width = 12, height = 8)

# PCA analysis on CLR transformed data
pca_result <- mixOmics::pca(CLR, ncomp = 3, scale = TRUE, center = TRUE)
plotIndiv(pca_result, comp = c(1, 2))
#plotIndiv(pca_result, style = '3d') 


pca_scores <- as.data.frame(pca_result$variates$X)  # Get scores for the components
pca_scores$Sample <- rownames(pca_scores)  # Add sample names as a column

# Merge PCA scores with site information
pca_plot_data <- merge(pca_scores, site_info, by = "Sample")

# Creates the plot
ggplot(data = pca_plot_data, 
       aes(x = PC1, y = PC2, colour = Site, shape = Period)) + #, shape = Relatedness
  geom_point(size = 7, alpha = 0.7) + 
#  stat_ellipse(aes(group = Site), level = 0.90) +
  labs(x = "PC1 (18%)", #NOTE: put in PC values
       y = "PC2 (17%)", #NOTE: put in PC values
       title = "PCA with CLR transformation", size = 16) +
  scale_colour_brewer(palette = "Set1",
                     limits = c("Long_Bay", "Kahukura", 
                                "Whenua_Hou_pre", 
                                "Whenua_Hou_post"),
                     labels = c("Long Bay pre-contact (n=4)", 
                                "Kahukura pre-contact (n=4)", 
                                "Whenua Hou pre-contact (n=5)", 
                                "Whenua Hou post-contact (n=3)"),
                     guide = guide_legend(nrow = 2)) +
  scale_shape_manual(values = c(17, 16), 
                     labels = c("Post-contact", "Pre-contact"),
                     guide = guide_legend(nrow = 2)) +
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

ggsave("PCA_CLR_decontam_t0.75_fraction0.01_2024-standard_minimizer_extra_env_removed.png", 
       width = 10, height = 8)

# Calculate Euclidean distances between samples and create cluster dendrogram and heatmap
euclidean_dist <- vegdist(CLR, method = "euclidean")

fit <- hclust(euclidean_dist, method="complete")
plot(fit) 
heatmap(as.matrix(euclidean_dist))

##########################
### PERMANOVA ANALYSIS ###
##########################

#test for the different variables by changing the variable after ~
permanova <- vegan::adonis2(CLR ~ Period,
                                   data = site_info,
                                   method = "euclidean",
                                   permutations = 9999)

permanova

#################################
######## DESeq2 ANALYSIS ########
#################################

# DESeq2 only works on count data, so tested otu_num and otu_frac turned into "count data"
# Analysis in thesis based on the otu_frac approach
otu_deseq <- otu_frac * 100000
#otu_deseq <- otu_num

dds <- DESeqDataSetFromMatrix(countData = round(t(otu_deseq)), 
                              colData = site_info, 
                              design = ~ Period)

dds <- DESeq(dds)

# Test for various variables by changing the data within c()
res <- results(dds, contrast = c("Period", 
                                 "pre-contact", 
                                 "post-contact")) %>%
  as.data.frame()

# View the results for the comparison
head(res)

# Add taxonomy names to taxids
res_tax <- rownames_to_column(res) %>%
  as.data.frame()
colnames(res_tax)[1] <- "taxonomy_id"

res_info <- samples %>%
  filter(taxonomy_id %in% res_tax$taxonomy_id) %>%
  dplyr::select(taxonomy_id, name) 

res_all <- merge(res_info, 
                          res_tax, by = "taxonomy_id", all = TRUE)

vulcano <- ggplot(res_all, aes(x = log2FoldChange, y = -log10(padj),
                            text = paste("Taxid:", name))) +
  geom_point(aes(color = padj < 0.01), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot of Differentially Abundant Taxa North vs South", # NOTE change title
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value")

vulcano <- ggplot(res_all, aes(x = log2FoldChange, y = -log10(padj),
                               text = paste("Taxid:", name))) +
  geom_point(aes(color = abs(log2FoldChange) > 5 & padj < 0.01), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot of Differentially Abundant Taxa North vs South",
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value")


vulcano_interactive <- ggplotly(vulcano)

vulcano_interactive

#htmlwidgets::saveWidget(vulcano_interactive, "interactive_vulcano_plot_t0.75_c0.001_pre_vs_post-contact.html")

# Extract significant taxa and output as CSV file
significant_taxa <- subset(res, padj < 0.01) %>%
  as.data.frame()

significant_taxa <- subset(significant_taxa, log2FoldChange < -5 | log2FoldChange > 5)

head(significant_taxa) # View the significant taxa

significant_taxid <- rownames_to_column(significant_taxa) %>%
  as.data.frame()
colnames(significant_taxid)[1] <- "taxonomy_id"

significant_info <- samples %>%
  filter(taxonomy_id %in% significant_taxid$taxonomy_id) %>%
  dplyr::select(taxonomy_id, name) 

significant_info <- merge(significant_info, 
  significant_taxid, by = "taxonomy_id", all = TRUE)

#write.csv(significant_info, file = "significant_taxa_t0.75_c0.001_pre_vs_post-contact.csv", row.names = FALSE)

###################################
######## MaAslin2 ANALYSIS ########
###################################

colnames(CLR) <- samples$name[match(colnames(CLR), samples$taxonomy_id)]

maaslin_metadata <- column_to_rownames(site_info, var = "sample")

# Test for various variables by changing the data 
# Normalisation and transformation set to NONE as using CLR data
fit_data = Maaslin2(input_data     = CLR, 
                    input_metadata = maaslin_metadata, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    transform      = "NONE",
                    output         = "Maaslin_Sex_CLR_output", 
                    fixed_effects  = c("Sex"),
                    heatmap_first_n = 100,
                    max_pngs        = 150,
                    reference      = c("Sex,male")
                    )

################################
### TOP X SPECIES COPROLITES ###
################################

samples_frac <- as.data.frame(samples_frac)
rownames(samples_frac) <- samples_frac[[1]]
samples_frac <- samples_frac[, -c(1,2)]

# Sum across all samples for each taxonomy_id and get the top X
top_taxa <- samples_frac %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>%
  head(21) %>%
  names()

# Subset the OTU table to only include the top X
otu_top20 <- samples_frac[top_taxa, ]

################################
### TOP X SPECIES PER SAMPLE ###
################################

top_n <- 10

sample_columns <- c(3:18)  # Replace with the actual column indices for your samples

sample_names <- colnames(samples_frac[sample_columns])

results_list <- list()

# Loop through each sample column by its index and get top 10 names and values
for (sample_col in sample_names) {
  top_data <- samples_frac %>%
    arrange(desc(!!sym(sample_col))) %>%
    slice_head(n = top_n) %>%
    dplyr::select(name, all_of(sample_col))  # Select only the name and the current sample column
  colnames(top_data) <- c(paste0("species"), "fraction")
  results_list[[sample_col]] <- top_data
}

### Combine all the top sample data into a single data frame
final_result <- do.call(cbind, results_list)

### Write the combined data frame to a CSV file
write.csv(final_result, "top_species_standard2024.csv", row.names = FALSE)

### print only the species names ##
#species_only <- dplyr::select(final_result, 
#                              contains(".species")) %>% 
#  dplyr::rename_with(~str_remove(., '.species')) #removed suffix
#write.csv(species_only, "top_species_names.csv", row.names = FALSE)

frac_long <- final_result %>%
  pivot_longer(cols = -starts_with("Sample"),  # Select all columns except for sample names
               names_to = c("Sample", ".value"),
               names_sep = "\\.") 

unique_species_count <- frac_long %>%
  summarise(unique_species = n_distinct(species)) %>%
  print()

p <- ggplot(frac_long, aes(x = Sample, y = fraction, fill = species, text = paste("Species:", species, "<br>Fraction:", fraction))) +
  geom_bar(stat = "identity") +  
  labs(x = "Samples", y = "Fraction", fill = "Species") +  
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14))

# Convert ggplot to an interactive plotly plot
interactive_plot <- ggplotly(p, tooltip = "text")
interactive_plot
htmlwidgets::saveWidget(interactive_plot, 
                        "interactive_barplot_top10_species_standard2024_per_sample.html")

#####################
### FISH SPOILERS ###
#####################

otu_fish <- samples_frac[ , -(2)]
otu_fish <- column_to_rownames(otu_fish, var = "name")

fish_CLR <- microbiome::transform(otu_fish, "clr", pseudocount = 0)
dim(fish_CLR)

species_list <- c("Lactococcus cremoris", #"Lactococcus petauri", "Lactococcus garvieae", 
                  "Brochothrix thermosphacta", "Carnobacterium maltaromaticum",
                  "Photobacterium phosphoreum", "Photobacterium toruni",
                  "Rhodococcus sp. 008", "Corynebacterium stationis",
                  "Jeotgalicoccus sp. WY2", "Chryseobacterium balustinum",
                  "Pseudoalteromonas nigrifaciens", "Pseudomonas bubulae",
                  "Pseudomonas fragi", "Pseudomonas lundensis", "Pseudomonas psychrophile",
                  "Psychrobacter alimentarius", "Psychrobacter sp. CLB018", 
                  "Psychrobacter sp. WY6", "Psychrobacter sp. P11F6",
                  "Shewanella baltica", "Shewanella sp. Pdp11", "Shewanella putrefaciens"
                  )

fish_CLR <- fish_CLR[rownames(fish_CLR) %in% species_list, ]

column_order <- c("MS11679", "MS11683", "MS11684", "MS11686", "MS11770", "MS11771", 
                  "MS11774", "MS11775", "MS11674", "MS11675", "MS11676", "MS11677", 
                  "MS11678", "MS11669", "MS11670", "MS11673")

fish_CLR <- fish_CLR[ , column_order]

site_info <- read.csv("sample_site_info.csv", header=TRUE)
meta_data <- site_info[1:16,]

rownames(meta_data) <- meta_data$sample
annotation_col <- meta_data["Site"]

source_colors <- setNames(RColorBrewer::brewer.pal(length(unique(meta_data$Site)), "Set1"), 
                           unique(meta_data$Site))

site_names <- c("Long_Bay", "Kahukura", 
                "Whenua_Hou_pre", 
                "Whenua_Hou_post")

# Extract colors from "Set1" palette (ensure enough colors are chosen)
source_colors <- setNames(brewer.pal(n = length(site_names), "Set1"), site_names)


p <- pheatmap(fish_CLR, 
              annotation_col = annotation_col, 
              annotation_colors = list(Site = source_colors),
              show_colnames = TRUE, 
              show_rownames = TRUE,
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              legend = TRUE,
              annotation_legend = FALSE)

legend_plot <- ggplot(annotation_col, aes(x = 1, fill = Site)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)"),
                      guide = guide_legend(nrow = 2)) +
  theme(legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16))

# Extract only the legend
legend <- get_legend(legend_plot)

# Combine the heatmap and legend
combined_plot <- plot_grid(p$gtable, legend, ncol = 1, rel_heights = c(1, 0.15))

# Save to file (PDF or PNG)
ggsave("heatmap_fish_spoilers_per_sample_no_clustering.png", 
       combined_plot, width = 8, height = 10, bg = "white")


ggsave("heatmap_fish_spoilers_per_sample_no_clustering.png", plot = p, width = 8, height = 8)  

#######################
### LOSS OF SPECIES ###
#######################

otu_period <- samples_frac[ , -(2)]
otu_period <- column_to_rownames(otu_period, var = "name")

period_CLR <- microbiome::transform(otu_period, "clr", pseudocount = 0)

dim(period_CLR)

species_list <- c(#"Streptococcus suis", 
                  #"Streptococcus salivarius",
                  #"Streptococcus thermophilus",
                  "Bacteroides sp. A1C1",  
                  "Bacteroides humanifaecis"#, 
                  #"Bacteroides sp. HF-162",
                  #"Fusobacterium animalis"
)

period_CLR <- period_CLR[rownames(period_CLR) %in% species_list, ]
period_CLR <- period_CLR[grepl("^Streptococcus\\b", rownames(period_CLR), ignore.case = TRUE), ]
period_CLR <- period_CLR[grepl("^Bacteroides\\b", rownames(period_CLR), ignore.case = TRUE), ]

site_info <- read.csv("sample_site_info.csv", header=TRUE)
meta_data <- site_info[1:16,]

period_CLR <- period_CLR[ , column_order]

#annotation_col <- annotation_col %>% 
#    rownames_to_column("Sample")

period_long <- period_CLR %>% 
  as.data.frame() %>%
  rownames_to_column("Species") %>% 
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance") %>% 
  left_join(annotation_col, by = "Sample")

period_long <- period_long %>% 
  mutate(Species = factor(Species, levels = rownames(period_CLR)))

period_long$Site <- factor(period_long$Site, 
                          levels = c("Long_Bay", "Kahukura", 
                                     "Whenua_Hou_pre", 
                                     "Whenua_Hou_post"))

p <- pheatmap(period_CLR, 
              annotation_col = annotation_col, 
              annotation_colors = list(Site = source_colors),
              show_colnames = TRUE, 
              show_rownames = TRUE,
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              legend = TRUE,
              annotation_legend = FALSE)

legend_plot <- ggplot(annotation_col, aes(x = 1, fill = Site)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set1",
                    limits = c("Long_Bay", "Kahukura", 
                               "Whenua_Hou_pre", 
                               "Whenua_Hou_post"),
                    labels = c("Long Bay pre-contact (n=4)", 
                               "Kahukura pre-contact (n=4)", 
                               "Whenua Hou pre-contact (n=5)", 
                               "Whenua Hou post-contact (n=3)"),
                    guide = guide_legend(nrow = 2)) +
  theme(legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16))

# Extract only the legend
legend <- get_legend(legend_plot)

# Combine the heatmap and legend
combined_plot <- plot_grid(p$gtable, legend, ncol = 1, rel_heights = c(1, 0.15))

ggsave("heatmap_temporal difference_no_clustering.png", 
       combined_plot, width = 8, height = 8, bg = "white")


ggplot(period_long, aes(x = Site, y = Abundance, fill = Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +
  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Site comparison",
       x = "Source", y = "CLR abundance",) +
  scale_fill_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)"),
                      guide = guide_legend(nrow = 2)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 20),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    strip.text = element_text(size = 20,, face = "bold"), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20)
  )

ggsave("boxplots_pre-post_strepto_bacter.png", 
       width = 10, height = 7)


ggsave("heatmap_pre-post_per_sample_no_clustering.png", plot = p, width = 8, height = 8)  


###############
### VIRUSES ###
###############

otu_virus <- samples %>% 
  filter(grepl("virus", !!sym(names(samples)[1]), ignore.case = TRUE))

otu_virus <- otu_virus[ , -(2)]
otu_virus <- column_to_rownames(otu_virus, var = "name")

virusT <- t(otu_virus) %>%
  as.data.frame()

count_cut_off <- 100
virusT <- virusT[, apply(virusT, 2, function(col) any(col >= count_cut_off))]

virus_CLR <- microbiome::transform(virusT, "clr", pseudocount = 0)

dim(virus_CLR)

p <- pheatmap(t(virus_CLR), 
         show_colnames = TRUE, 
         show_rownames = TRUE,
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         legend = TRUE)

# Save the plot to a file (e.g., PNG)
ggsave("heatmap_viruses_per_sample_clustering.png", plot = p, width = 8, height = 12)  

tune_pca <- tune.pca(virus_CLR, ncomp = 10, scale = TRUE)
variation <- as.data.frame(tune_pca$prop_expl_var) %>%
  rownames_to_column()
variation$PC <- rownames(variation)
variation$cum <- tune_pca$cum.var

ggplot(variation, aes(x = fct_inorder(PC), y = X)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  #  geom_point(data = variation, aes(x = fct_inorder(PC), y = cum), 
  #             color = "red") +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal() 
#elbow occurs at 2

#ggsave("Scree plot - variance explained principal components.png", 
#       width = 12, height = 8)

# PCA analysis on otu
pca_result <- mixOmics::pca(virus_CLR, ncomp = 2, scale = TRUE, center = TRUE)
plotIndiv(pca_result, comp = c(1, 2))
#plotIndiv(pca_result, style = '3d') 


pca_scores <- as.data.frame(pca_result$variates$X)  # Get scores for the components
pca_scores$Sample <- rownames(pca_scores)  # Add sample names as a column

# Merge PCA scores with site information
pca_plot_data <- merge(pca_scores, site_info, by = "Sample")

# Creates the plot
ggplot(data = pca_plot_data, 
       aes(x = PC1, y = PC2, colour = Site, shape = Period)) + #, shape = Relatedness
  geom_point(size = 5, alpha = 0.7) + #size = 5, alpha = 0.7
  #  stat_ellipse(aes(group = Site), level = 0.90) +
  labs(x = "PC1 (29%)",
       y = "PC2 (21%)",
       title = "PCA viruses with CLR transformation") +
  scale_colour_brewer(palette = "Set1",
                      limits = c("Long_Bay", "Kahukura", 
                                 "Whenua_Hou_pre", 
                                 "Whenua_Hou_post"),
                      labels = c("Long Bay pre-contact (n=4)", 
                                 "Kahukura pre-contact (n=4)", 
                                 "Whenua Hou pre-contact (n=5)", 
                                 "Whenua Hou post-contact (n=3)"),
                      guide = guide_legend(nrow = 2)) +
  scale_shape_manual(values = c(17, 16), 
                     labels = c("Post-contact", "Pre-contact"),
                     guide = guide_legend(nrow = 2)) +
  theme(
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"), 
    legend.text = element_text(size = 16),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    axis.title.x = element_text(size = 16),  # Change x-axis title size
    axis.title.y = element_text(size = 16)) 

######################################################################
#### COMPARISON MODERN AND SEDIMENT SAMPLES TOP X SPECIES OVERALL ####
######################################################################

merged <- read.csv("combined_kuri_sources_bracken_species_march2025.csv", header=TRUE)

merged <- dplyr::select(merged,
                        name, taxonomy_id,
                        contains("num")) %>%
  dplyr::rename_with(~str_remove(., '_kraken_standard2024_conf0.50.txt_filtered.bracken_num')) %>%
  dplyr::rename_with(~str_remove(., '_filtered.bracken_num')) %>%
  dplyr::rename_with(~str_remove(., '.bracken_num')) %>%
  dplyr::rename_with(~str_remove(., '_minimizer')) %>%
  column_to_rownames(var = "name")

merged <- merged[ , -(1)]

# Sum across all samples for each taxonomy_id and get the top X
top_taxa <- merged %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>%
  head(100) %>%
  names()

# Clr transform for better visualization
otu_top20_clr <- microbiome::transform(merged, "clr", pseudocount = 0)

# Subset the OTU table to only include the top X
otu_top20_clr <- otu_top20_clr[grepl("^Helicobacter\\b", rownames(otu_top20_clr), ignore.case = TRUE), ]
otu_top20_clr <- otu_top20_clr[top_taxa, ]

#remove sample names and give colours
meta_data <- read.csv("Merged_ids.csv", header = TRUE)

rownames(meta_data) <- meta_data$sample
annotation_col <- meta_data["Source"]

source_colors <- setNames(RColorBrewer::brewer.pal(length(unique(meta_data$Source)), "Set1"), 
                          unique(meta_data$Source))

otu_top20_clr <- otu_top20_clr[, meta_data$sample]

otu_top20_clr <- otu_top20_clr[, !colnames(otu_top20_clr) %in% c("ERR3761400", "ERR3761401",
                                                                 "ERR3761402", "ERR3761404",
                                                                 "ERR3761405", "ERR3761406",
                                                                 "ERR10114881", "SRR7774469",
                                                                 "SRR7774472", "SRR7774473",
                                                                 "SRR7774471", "SRR7774474",
                                                                 "SRR7774476", "SRR7774477",
                                                                 "Blank1_WH", "Blank2_WH",
                                                                 "KH_blank_1", "KH_blank_2", 
                                                                 "LB_blank_1", "LB_blank_2")]

p <- pheatmap(otu_top20_clr, 
              annotation_col = annotation_col, 
              annotation_colors = list(Source = source_colors),
              show_colnames = FALSE, 
              show_rownames = TRUE,
              cluster_rows = TRUE, 
              cluster_cols = FALSE) 

Human        Dog_India         Dog_Laos Dog_South_Africa 
"#E41A1C"        "#377EB8"        "#4DAF4A"        "#984EA3" 
Dog_USA             Kuri         Sediment          NZ_bone 
"#FF7F00"        "#FFFF33"        "#A65628"        "#F781BF" 

annotation_col$Source <- factor(annotation_col$Source, levels = names(source_colors))

p <- pheatmap(otu_top20_clr, 
              annotation_col = annotation_col, 
              annotation_colors = list(Source = source_colors),
              show_colnames = FALSE, 
              show_rownames = TRUE,
              cluster_rows = TRUE, 
              cluster_cols = FALSE,
              legend = TRUE,
              annotation_legend = FALSE)

legend_plot <- ggplot(annotation_col, aes(x = 1, fill = Source)) +
  geom_bar() +
  scale_fill_manual(values = source_colors) +
  guides(fill = guide_legend(ncol = 4)) +
  theme(legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16))

legend <- get_legend(legend_plot)

# Combine the heatmap and legend
combined_plot <- plot_grid(p$gtable, legend, ncol = 1, rel_heights = c(1, 0.15))

# Assuming otu_top20_clr is a matrix (taxa as rows, samples as columns)
df <- as.data.frame(otu_top20_clr)

# Add taxonomy as a column
df$Taxon <- rownames(df)

# Reshape to long format
df_long <- melt(df, id.vars = "Taxon",
                variable.name = "Sample", value.name = "Abundance")

# Add Source info from annotation_col
df_long <- df_long %>%
  left_join(annotation_col %>% tibble::rownames_to_column("Sample"),
            by = "Sample")

# Make boxplots per Source for each Taxon
p <- ggplot(df_long, aes(x = Source, y = Abundance, fill = Source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +
  facet_wrap(~ Taxon, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Top 20 bacterial species in coprolites",
       x = "Source", y = "CLR abundance",) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"), 
    legend.text = element_text(size = 18),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18)
  )

ggsave("Boxplot_modern_helicobacter.png", plot = p, 
       width = 12, height = 14,bg = "white")  
#ggsave("heatmap_top100_species_all_samples_overall_clustering_filter_march2025.png", plot = p, width = 12, height = 12)  

#########################################################################
#### COMPARISON MODERN AND SEDIMENT SAMPLES TOP X SPECIES PER SAMPLE ####
#########################################################################

# Identify the top 10 species per sample
top_species_per_sample <- apply(merged, 2, function(x) {
  names(sort(x, decreasing = TRUE)[1:3])  # Get top X species for each sample
})

# Get unique species from all top 10 lists
unique_species <- unique(unlist(top_species_per_sample))

# Subset OTU table to include only these species
otu_top <- merged[unique_species, , drop = FALSE]

# Convert to data frame for processing
otu_top_df <- as.data.frame(otu_top)
otu_top_df$Species <- rownames(otu_top_df)  # Keep original row names

# Clean species names (remove .1, .2, etc.)
otu_top_df$Species <- str_replace(otu_top_df$Species, "\\.\\d+$", "")

# Merge duplicate species by summing them
otu_top_df <- otu_top_df %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Convert back to matrix with unique row names
otu_top_matrix <- as.matrix(column_to_rownames(otu_top_df, var = "Species"))

# Log transform for better visualization
otu_top_clr <- microbiome::transform(otu_top_matrix, "clr", pseudocount = 0)

# Create heatmap 
otu_top_clr <- otu_top_clr[, !colnames(otu_top_clr) %in% c("ERR3761400", "ERR3761401",
                                                           "ERR3761402", "ERR3761404",
                                                           "ERR3761405", "ERR3761406",
                                                           "ERR10114881", "SRR7774469",
                                                           "SRR7774472", "SRR7774473",
                                                           "SRR7774471", "SRR7774474",
                                                           "SRR7774476", "SRR7774477",
                                                           "Blank1_WH", "Blank2_WH",
                                                           "KH_blank_1", "KH_blank_2", 
                                                           "LB_blank_1", "LB_blank_2")]

p <- pheatmap(otu_top_clr, 
              annotation_col = annotation_col, 
              annotation_colors = list(source = source_colors),
              show_colnames = TRUE, 
              show_rownames = TRUE,
              cluster_rows = TRUE, 
              cluster_cols = TRUE) 

#ggsave("heatmap_top3_species_all_samples_clustering.png", plot = p, width = 24, height = 36)  

#########################################
#### TOP COPROLITE SPECIES COMPARISON ###
#########################################

merged_clr <- microbiome::transform(merged, "clr", pseudocount = 0)

merged_clr <- merged_clr[, meta_data$sample]

merged_clr <- merged_clr[, !colnames(merged_clr) %in% c("ERR3761400", "ERR3761401",
                                                           "ERR3761402", "ERR3761404",
                                                           "ERR3761405", "ERR3761406",
                                                           "ERR10114881", "SRR7774469",
                                                           "SRR7774472", "SRR7774473",
                                                           "SRR7774471", "SRR7774474",
                                                           "SRR7774476", "SRR7774477",
                                                           "Blank1_WH", "Blank2_WH",
                                                           "KH_blank_1", "KH_blank_2", 
                                                           "LB_blank_1", "LB_blank_2")]

merged_clr <- merged_clr[rownames(merged_clr) %in% c(rownames(otu_top20)), ]

annotation_col <- annotation_col %>% 
  rownames_to_column("Sample")
merged_long <- merged_clr %>% 
  as.data.frame() %>%
  rownames_to_column("Species") %>% 
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance") %>% 
  left_join(annotation_col, by = "Sample")


merged_long <- merged_long %>% 
  mutate(Species = factor(Species, levels = rownames(otu_top20)))

ggplot(merged_long, aes(x = Source, y = Abundance, fill = Source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +
  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "Top 20 bacterial species in coprolites",
       x = "Source", y = "CLR abundance",) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18, 
                              face = "bold", , hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"), 
    legend.text = element_text(size = 18),
    legend.position = "bottom",  
    legend.direction = "horizontal",
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 18)
  )

ggsave("boxplots_top21_species_in_coprolites_vs_modern_sources.png", 
              width = 12, height = 16)
       
#########################################
### MAASLIN ANLYSIS MODERN COMPARISON ###
#########################################

fit_data = Maaslin2(input_data      = otu_top20_clr, 
                    input_metadata  = annotation_col, 
                    min_prevalence  = 0,
                    normalization   = "NONE",
                    transform       = "NONE",
                    output          = "Maaslin_Laos_CLR_output", 
                    fixed_effects   = c("Source"),
                    heatmap_first_n = 100,
                    max_pngs        = 100,
                    reference       = c("Source,Dog_Laos")
)


##################################
### TOP X SPECIES BLANK & BONE ###
##################################

top_n <- 10

sample_columns <- c(3:16)  # Replace with the actual column indices for your samples

sample_names <- colnames(blanks[sample_columns])

results_list <- list()

# Loop through each sample column by its index and get top 10 names and values
for (sample_col in sample_names) {
  top_data <- blanks %>%
    arrange(desc(!!sym(sample_col))) %>%
    slice_head(n = top_n) %>%
    dplyr::select(name, all_of(sample_col))  # Select only the name and the current sample column
  colnames(top_data) <- c(paste0("species"), "fraction")
  results_list[[sample_col]] <- top_data
}

### Combine all the top sample data into a single data frame
final_result <- do.call(cbind, results_list)

### Write the combined data frame to a CSV file
write.csv(final_result, "top_species_blanks_and_bone_standard2024.csv", row.names = FALSE)

### print only the species names ##
#species_only <- dplyr::select(final_result, 
#                              contains(".species")) %>% 
#  dplyr::rename_with(~str_remove(., '.species')) #removed suffix
#write.csv(species_only, "top_species_names.csv", row.names = FALSE)

frac_long <- final_result %>%
  pivot_longer(cols = -starts_with("Sample"),  # Select all columns except for sample names
               names_to = c("Sample", ".value"),
               names_sep = "\\.") 

unique_species_count <- frac_long %>%
  summarise(unique_species = n_distinct(species)) %>%
  print()

p <- ggplot(frac_long, aes(x = Sample, y = fraction, fill = species, text = paste("Species:", species, "<br>Fraction:", fraction))) +
  geom_bar(stat = "identity") +  
  labs(x = "Samples", y = "Fraction", fill = "Species") +  
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14))

# Convert ggplot to an interactive plotly plot
interactive_plot <- ggplotly(p, tooltip = "text")
interactive_plot
htmlwidgets::saveWidget(interactive_plot, 
                        "interactive_barplot_top10_species_standard2024_per_sample.html")