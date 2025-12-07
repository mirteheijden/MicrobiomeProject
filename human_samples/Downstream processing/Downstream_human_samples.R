# Load required packages
library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(ggplot2)

#files

##TSE from GG2
tse_GG2 <- loadFromQIIME2("path_to/RESULTS/GG2/counts.qza", taxonomy = "path_to/RESULTS/GG2/taxonomy.qza", sampleMetaFile="path_to/RESULTS/GG2/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))


##TSE from Metaphlan
tse_Metaphlan <- loadFromMetaphlan("path_to/RESULTS/Metaphlan/metaphlan_db_meta4_combined_reports.txt", sample_meta="path_to/RESULTS/Metaphlan/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))



# Define tse
tse <- tse_GG2
#tse <- tse_Metaphlan



### Beta diversity PCOA    -- DC vs No_DC

# Species level
tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]

# Metaphlan (check availability altExpNames(tse))
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Extract day 1 samples
tse_day1 <- tse_rel[, rownames(subset(colData(tse_rel), Sample_point == "Day_1"))]

#print(assay_names)
tse_day1_rel_abund_assay <- assays(tse_day1)$relabundance
bray_curtis_dist <- vegan::vegdist(t(tse_day1_rel_abund_assay), method = "bray")
#install.packages("ecodist")
#library(ecodist)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
#PCoA1 and PCoA2
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
# Create a plot
# Adds the variable you want to use for coloring to the data frame
bray_curtis_RECIST_pcoa_df <- cbind(bray_curtis_pcoa_df,
                                    Response = colData(tse_day1)$Response.day78.,
                                    label = colData(tse_day1)$Patient)


# For DC vs No_DC
bray_curtis_plot <- ggplot(data = bray_curtis_RECIST_pcoa_df, 
                           aes(x = pcoa1, y = pcoa2, color = Response)) +
  geom_point() +  # Points for each sample
  geom_text(aes(label = label), size = 3, vjust = 1.5, show.legend = FALSE) +  # Labels for each point (Patient names)
  
  # Add circles around the groups with colored outlines and no fill
  geom_mark_ellipse(aes(group = Response), 
                    fill = NA,  # No fill
                    linetype = "solid",  # Outline type
                    show.legend = FALSE) +  # Hide legend for the ellipses
  
  labs(x = "PC1",
       y = "PC2", 
       title = "Beta diversity between Response groups") +
  theme_bw() + #white background
  theme(
    title = element_text(size = 10),
    plot.title = element_text(face = "plain", size = 18),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    text = element_text(size = 12),
    
    # Background and border tweaks
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # Add axis lines with thickness
    axis.line = element_line(color = "black", size = 0.8)
  ) +
  scale_color_manual(values = c("DC" = "blue", "No_DC" = "#C21E56")) + # Define custom colors
  # Zoom out a bit: extend axes limits
  expand_limits(x = c(min(bray_curtis_RECIST_pcoa_df$pcoa1) - 0.1, 
                      max(bray_curtis_RECIST_pcoa_df$pcoa1) + 0.22),
                y = c(min(bray_curtis_RECIST_pcoa_df$pcoa2) - 0.1, 
                      max(bray_curtis_RECIST_pcoa_df$pcoa2) + 0.1))


# Display the plot
bray_curtis_plot





### Beta diversity PCOA    -- ABX

# Species level
tse <- tse_GG2
tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)
tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]
tse <- tse[, which(colData(tse)$Antibiotics %in% c("Yes", "No"))]

# Metaphlan (check availability altExpNames(tse))
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Extract day 1 samples
tse_day1 <- tse_rel[, rownames(subset(colData(tse_rel), Sample_point == "Day_1"))]

#print(assay_names)
tse_day1_rel_abund_assay <- assays(tse_day1)$relabundance
bray_curtis_dist <- vegan::vegdist(t(tse_day1_rel_abund_assay), method = "bray")
#library(ecodist)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
#PCoA1 and PCoA2
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
# Create a plot
# Adds the variable you want to use for coloring to the data frame
bray_curtis_RECIST_pcoa_df <- cbind(bray_curtis_pcoa_df,
                                    Response = colData(tse_day1)$Antibiotics,
                                    label = colData(tse_day1)$Patient)


# For No vs Yes
bray_curtis_plot <- ggplot(data = bray_curtis_RECIST_pcoa_df, 
                           aes(x = pcoa1, y = pcoa2, color = Response)) +
  geom_point() +  # Points for each sample
  geom_text(aes(label = label), size = 3, vjust = 1.5, show.legend = FALSE) +  # Labels for each point (Patient names)
  
  # Add circles around the groups with colored outlines and no fill
  geom_mark_ellipse(aes(group = Response), 
                    fill = NA,  # No fill
                    linetype = "solid",  # Outline type
                    show.legend = FALSE) +  # Hide legend for the ellipses
  
  labs(x = "PC1",
       y = "PC2", 
       title = "Beta diversity between Antibiotic groups",
       color = "Antibiotic use") +
  theme_bw() + #white background
  theme(
    title = element_text(size = 10),
    plot.title = element_text(face = "plain", size = 18),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    text = element_text(size = 12),
    
    # Background and border tweaks
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # Add axis lines with thickness
    axis.line = element_line(color = "black", size = 0.8)
  ) +
  theme(title = element_text(size = 10)) +  # Makes titles smaller
  scale_color_manual(values = c("No" = "purple", "Yes" = "black")) + # Define custom colors
  # Zoom out a bit: extend axes limits
  expand_limits(x = c(min(bray_curtis_RECIST_pcoa_df$pcoa1) - 0.1, 
                      max(bray_curtis_RECIST_pcoa_df$pcoa1) + 0.22),
                y = c(min(bray_curtis_RECIST_pcoa_df$pcoa2) - 0.1, 
                      max(bray_curtis_RECIST_pcoa_df$pcoa2) + 0.1))


# Display the plot
bray_curtis_plot



#######################################################################################
#AlphaDiversity - Response

# Load required packages
library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(ggplot2)
library(scater)

#files

##TSE from GG2
tse_GG2 <- loadFromQIIME2("path_to/RESULTS/GG2/counts.qza", taxonomy = "path_to/RESULTS/GG2/taxonomy.qza", sampleMetaFile="path_to/RESULTS/GG2/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))

##TSE from Metaphlan
tse_Metaphlan <- loadFromMetaphlan("path_to/RESULTS/Metaphlan/metaphlan_db_meta4_combined_reports.txt", sample_meta="path_to/RESULTS/Metaphlan/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))

#Define tse
tse <- tse_GG2
#tse <- tse_Metaphlan

tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]
tse <- tse[, which(colData(tse)$Antibiotics %in% c("Yes", "No"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Estimate (observed) richness
tse_rel <- addAlpha(
  tse_rel, assay.type = "relabundance", index = "shannon", name = "observed",
  detection = 10)

# Check some of the first values in colData
head(tse_rel$observed)

# Set the specific order for the Patient factor levels
colData(tse_rel)$Patient <- factor(colData(tse_rel)$Patient, 
                                   levels = c("20206", "20103", "20212", "30207", 
                                              "20108", "20211", "20104", "30209"))

# Plot
plotColData(
  tse_rel,
  y = "observed",
  x = "Patient_Response",
  colour_by = "Antibiotics",
  size_by = "observed"
) +
  theme(
    title = element_text(size = 10),
    plot.title = element_text(face = "plain", size = 18),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
    axis.text.y = element_text(size = 12, face = "plain"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    text = element_text(size = 12),
    
    # Background and border tweaks
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # Add axis lines with thickness
    axis.line = element_line(color = "black", size = 0.8)
  ) +
  labs(
    x = "Patient",
    y = expression(Richness[Observed]),
    colour = "Antibiotics"
  ) +
  scale_color_manual(values = c("Yes" = "black", "No" = "purple"))



###################################################################################################
## Taxonomic barplots

# Load required packages
library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(RColorBrewer)

#files

##TSE from GG2
tse_GG2 <- importQIIME2("path_to/RESULTS/GG2/counts.qza", taxonomy = "path_to/RESULTS/GG2/taxonomy.qza", sampleMetaFile="path_to/RESULTS/GG2/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))

##TSE from Metaphlan
tse_Metaphlan <- importMetaPhlAn("path_to/RESULTS/Metaphlan/metaphlan_db_meta4_combined_reports.txt", sample_meta="path_to/RESULTS/Metaphlan/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))



# Barplots

# phylum - GG2
tse <- tse_GG2

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]

# Extract phylum level
tse_spec <- agglomerateByRank(tse, rank ="phylum", onRankOnly=TRUE)

tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

#No_DC
No_DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "No_DC"))]

tse_No_DC_phylum <- agglomerateByRank(No_DC_rel, rank ="phylum", onRankOnly=TRUE)
top_No_DC_taxa <- getTop(tse_No_DC_phylum, top = 6)

#DC
DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "DC"))]

tse_DC_phylum <- agglomerateByRank(DC_rel, rank ="phylum", onRankOnly=TRUE)
top_DC_taxa <- getTop(tse_DC_phylum, top = 6)

# Combine the top 10 taxa (5 from No_DC and 5 from DC)
top_taxa <- union(top_No_DC_taxa, top_DC_taxa)

# Renaming the "phylum" rank to keep only top taxa and label others as "Other"
phylum_renamed <- sapply(rowData(tse_rel)$phylum,
                         function(x) { if (x %in% top_taxa) {x} else {"Other"} })

rowData(tse_rel)$phylum <- as.character(phylum_renamed)

# Now aggregate all "Other" taxa into a single row
tse_combined <- mergeFeaturesByRank(tse_rel, rank = "phylum", onRankOnly = TRUE)

# Change rownames of the colData to the Patient ID
col_data_tse <- colData(tse_combined)
rownames(col_data_tse) <- col_data_tse$Patient_Response
colData(tse_combined) <- col_data_tse

# Filter out the "Other" phylum from rowData
tse_filtered <- tse_combined[rowData(tse_combined)$phylum != "Other", ]

# Remove the "p__" prefix from the 'phylum' column in rowData
rowData(tse_filtered)$phylum <- gsub("^p__", "", rowData(tse_filtered)$phylum)

# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

custom_colors <- c(
  "Actinobacteriota" = color_palette[1],       # Specify your taxa and their corresponding colors
  "Bacteroidota" = color_palette[2],
  "Firmicutes_D" = color_palette[3],
  "Firmicutes_A" = color_palette[4],
  "Proteobacteria" = color_palette[5],
  "Synergistetes" = color_palette[6],
  "Verrucomicrobiota" = color_palette[7],
  "Other" = color_palette[8]
)

# Visualize the composition barplot without the "Other" category
plotAbundance(tse_filtered, rank = "phylum", assay.type = "relabundance", 
              order_rank_by = "abund", order_sample_by = "Patient_Response",
              decreasing = FALSE, add_x_text = TRUE) +
  theme(
    plot.title = element_text(face = "plain", size = 20),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(angle = 90, size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14, face = "plain"),
    text = element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")  
  ) +
  theme(title = element_text(size = 10)) +
  theme(legend.key.height = unit(0.5, "cm")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +  # forces start at 0, no padding
  labs(title = "Most abundant phylum per patient", x = "", y = "Rel. Abundance (%)") +
  scale_fill_manual(values = custom_colors) +
  guides(fill = guide_legend(title = "Phylum"),
         color = "none")





# genus -- GG2
tse <- tse_GG2
#tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]

# Transform to relative abundance
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# No_DC
No_DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "No_DC"))]
tse_No_DC_genus <- mergeFeaturesByRank(No_DC_rel, rank = "genus", onRankOnly = TRUE)
top_No_DC_taxa <- getTopFeatures(tse_No_DC_genus, top = 12)

# DC
DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "DC"))]
tse_DC_genus <- mergeFeaturesByRank(DC_rel, rank = "genus", onRankOnly = TRUE)
top_DC_taxa <- getTopFeatures(tse_DC_genus, top = 12)

# Combine the top taxa 
top_taxa <- union(top_No_DC_taxa, top_DC_taxa)

# Renaming the "genus" rank to keep only top taxa and label others as "Other"
genus_renamed <- sapply(rowData(tse_rel)$genus,
                        function(x) { if (x %in% top_taxa) {x} else {"Other"} })

rowData(tse_rel)$genus <- as.character(genus_renamed)

# Now aggregate all "Other" taxa into a single row
tse_combined <- mergeFeaturesByRank(tse_rel, rank = "genus", onRankOnly = TRUE)

# Change rownames of the colData to the Patient ID
col_data_tse <- colData(tse_combined)
rownames(col_data_tse) <- col_data_tse$Patient_Response
colData(tse_combined) <- col_data_tse

# Filter out the "Other" genus from rowData
tse_filtered <- tse_combined[rowData(tse_combined)$genus != "Other", ]

# Remove the "g__" prefix from the 'genus' column in rowData
rowData(tse_filtered)$genus <- gsub("^g__", "", rowData(tse_filtered)$genus)

# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

custom_colors <- c(
  "Akkermansia" = color_palette[1],       # Specify your taxa and their corresponding colors
  "Bacteroides_H" = color_palette[2],
  "Bifidobacterium_388775" = color_palette[3],
  "Citrobacter_A_692098" = color_palette[4],
  "Eggerthella" = color_palette[5],
  "Blautia_A_141781" = color_palette[6],
  "Gemmiger_A_73129" = color_palette[7],
  "Gemmiger_A_73129" = color_palette[8],
  "Agathobacter_164117" = color_palette[9],
  "Agathobaculum" = color_palette[10],
  "Dysosmobacter" = color_palette[11],
  #"" = color_palette[12],
  "Phocaeicola_A_858004" = color_palette[13],
  #"" = color_palette[14],
  "Alistipes_A_871400" = color_palette[15],
  "Scatavimonas" = color_palette[16],
  "Sellimonas" = color_palette[17],
  "Collinsella" = color_palette[18],
  #"" = color_palette[19],
  "Faecalibacterium" = color_palette[20],
  "Eisenbergiella" = color_palette[21],
  "Escherichia_710834" = color_palette[22],
  "Lawsonibacter" = color_palette[23],
  "Gordonibacter" = color_palette[24],
  "Ruthenibacterium" = color_palette[25],
  #"" = color_palette[26],
  "Ruminococcus_B" = color_palette[27],
  "Streptococcus" = color_palette[28],
  "taxa29" = color_palette[29],
  "taxa30" = color_palette[30]
)

# Visualize the composition barplot without the "Other" category
plotAbundance(tse_filtered, rank = "genus", assay.type = "relabundance", 
              order_rank_by = "abund", order_sample_by = "Eggerthella",
              decreasing = FALSE, add_x_text = TRUE) +
  theme(
    plot.title = element_text(face = "plain", size = 20),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(angle = 90, size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14, face = "plain"),
    text = element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black") 
    ) +
  theme(legend.key.height = unit(0.5, "cm")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +  # Set y-axis from 0 to 100%
  labs(title = "Most abundant genus per patient", x = "", y = "Rel. Abundance (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom color palette
  guides(fill = guide_legend(title = "Genus"),  
         color = "none")  # Remove the legend for the outline color



# species
tse <- tse_GG2
#tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$Response.day78. %in% c("No_DC", "DC"))]

# Extract species level
tse_spec <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

#No_DC
No_DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "No_DC"))]

tse_No_DC_species <- mergeFeaturesByRank(No_DC_rel, rank ="species", onRankOnly=TRUE)
top_No_DC_taxa <- getTopFeatures(tse_No_DC_species, top = 12)

#DC
DC_rel <- tse_rel[, rownames(subset(colData(tse_rel), Response.day78. == "DC"))]

tse_DC_species <- mergeFeaturesByRank(DC_rel, rank ="species", onRankOnly=TRUE)
top_DC_taxa <- getTopFeatures(tse_DC_species, top = 12)

# Combine the top taxa 
top_taxa <- union(top_No_DC_taxa, top_DC_taxa)

# Renaming the "species" rank to keep only top taxa and label others as "Other"
species_renamed <- sapply(rowData(tse_rel)$species,
                          function(x) { if (x %in% top_taxa) {x} else {"Other"} })

rowData(tse_rel)$species <- as.character(species_renamed)

# Now aggregate all "Other" taxa into a single row
tse_combined <- mergeFeaturesByRank(tse_rel, rank = "species", onRankOnly = TRUE)

# Change rownames of the colData to the Patient ID
col_data_tse <- colData(tse_combined)
rownames(col_data_tse) <- col_data_tse$Patient_Response
colData(tse_combined) <- col_data_tse

# Filter out the "Other" species from rowData
tse_filtered <- tse_combined[rowData(tse_combined)$species != "Other" & rowData(tse_combined)$species != "s__", ]

# Remove the "s__" prefix from the 'species' column in rowData
rowData(tse_filtered)$species <- gsub("^s__", "", rowData(tse_filtered)$species)

# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

custom_colors <- c(
  "Akkermansia muciniphila_D_776786" = color_palette[1],       # Specify your taxa and their corresponding colors
  "Collinsella sp000434535" = color_palette[2],
  "Bifidobacterium longum" = color_palette[3],
  "Citrobacter_A_692098 werkmanii" = color_palette[4],
  "Eggerthella lenta" = color_palette[5],
  "Bifidobacterium bifidum" = color_palette[6],
  "Gemmiger_A_73129 formicilis_72474" = color_palette[7],
  "Dysosmobacter welbionis" = color_palette[8],
  "Agathobacter rectalis" = color_palette[9],
  "Alistipes_A_871400 putredinis" = color_palette[10],
  "Flavonifractor plautii" = color_palette[11],
  "Gemmiger_A_73129 qucibialis" = color_palette[12],
  "Phocaeicola_vulgatus" = color_palette[13],
  "Sellimonas intestinalis" = color_palette[14],
  "Alistipes_putredinis" = color_palette[15],
  "Bacteroides_uniformis" = color_palette[16],
  "Collinsella sp000434535" = color_palette[17],
  "Collinsella aerofaciens_E" = color_palette[18],
  "Collinsella aerofaciens_F" = color_palette[19],
  "Faecalibacterium prausnitzii_C_71351" = color_palette[20],
  "Gordonibacter pamelaeae" = color_palette[21],
  "Scatavimonas merdipullorum" = color_palette[22],
  "Ruminococcus_B gnavus" = color_palette[23],
  "Gordonibacter_urolithinfaciens" = color_palette[24],
  "Blautia_A_141781 wexlerae" = color_palette[25],
  "Phocaeicola_A_858004 vulgatus" = color_palette[26],
  "Ruthenibacterium lactatiformans" = color_palette[27],
  "Streptococcus gallolyticus" = color_palette[28],
  "taxa29" = color_palette[29],
  "taxa30" = color_palette[30]
)

# Visualize the composition barplot without the "Other" category
plotAbundance(tse_filtered, rank = "species", assay.type = "relabundance", 
              order_rank_by = "abund", order_sample_by = "Eggerthella lenta",
              decreasing = FALSE, add_x_text = TRUE) +
  theme(
    plot.title = element_text(face = "plain", size = 20),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(angle = 90, size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14, face = "plain"),
    text = element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black") 
    ) +
  theme(legend.key.height = unit(0.5, "cm")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +  # Set y-axis from 0 to 100%
  labs(title = "Most abundant species per patient", x = "", y = "Rel. Abundance (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom color palette
  guides(fill = guide_legend(title = "Species"),  
         color = "none")  # Remove the legend for the outline color





#####################################################################################################





#####################################################################################################################


################ Wilcoxon, Log FC plot -  Genus
##TSE from GG2
tse_GG2 <- loadFromQIIME2("path_to/RESULTS/GG2/counts.qza", taxonomy = "path_to/RESULTS/GG2/taxonomy.qza", sampleMetaFile="path_to/RESULTS/GG2/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))

##TSE from Metaphlan
tse_Metaphlan <- loadFromMetaphlan("path_to/RESULTS/Metaphlan/metaphlan_db_meta4_combined_reports.txt", sample_meta="path_to/RESULTS/Metaphlan/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))
library(phyloseq)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Load and preprocess data
tse <- tse_GG2
#tse <- tse_Metaphlan
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]
rownames(tse) <- gsub(" ", "_", rownames(tse))

# Transform to relative abundance
tse_rel <- transformAssay(tse,
                          assay.type = "counts",
                          method = "relabundance")

# Convert to phyloseq and agglomerate at genus level
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse_rel)
pseq_genus <- tax_glom(pseq, taxrank = "genus")

# Extract abundance and metadata
otu_mat <- as.data.frame(otu_table(pseq_genus))
taxa <- tax_table(pseq_genus)[, "genus"]
sample_metadata <- sample_data(pseq_genus)

# Ensure groupings are factors
group <- factor(sample_metadata$Response.day78., levels = c("DC", "No_DC"))

# Wilcoxon test per genus
wilcox_results <- apply(otu_mat, 1, function(x) {
  test <- wilcox.test(x ~ group, exact = FALSE)
  logFC <- log2((mean(x[group == "DC"]) + 0.01) / (mean(x[group == "No_DC"]) + 0.01))
  c(p_value = test$p.value, logFC = logFC)
})

# Convert to data frame
wilcox_df <- as.data.frame(t(wilcox_results))
wilcox_df$taxon <- rownames(wilcox_df)
wilcox_df$genus <- tax_table(pseq_genus)[wilcox_df$taxon, "genus"]

# Adjust p-values
wilcox_df$adj_p_value <- p.adjust(wilcox_df$p_value, method = "fdr")

# Filter results
filtered_results <- wilcox_df %>%
  filter(p_value < 0.05 & (logFC > 2 | logFC < -2)) %>%
  mutate(response = ifelse(logFC > 0, "DC", "No_DC"),
         highlight = ifelse(logFC > 0, "darkturquoise", "pink"),
         hjust_pos = ifelse(logFC > 0, -0.1, 1.1),
         genus = gsub("^g__", "", genus)) %>%
  filter(!is.na(genus) & genus != "")

# Bar plot
bar_plot <- ggplot(filtered_results, aes(x = reorder(genus, logFC), 
                                         y = logFC, 
                                         fill = response)) +
  geom_bar(stat = "identity", show.legend = TRUE) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +  
  coord_flip() +
  scale_fill_manual(values = c("DC" = "#FADFA0", "No_DC" = "#DCDCDC")) +  
  labs(title = "Differential abundant genus between DC and No_DC",
       x = "", y = "Log Fold Change") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "plain", size = 34, hjust = 0.5),
    axis.title.x = element_text(size = 24, margin = margin(t = 10)),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 22),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(title = "", labels = c("DC", "No_DC"))) +
  geom_text(aes(label = genus, hjust = hjust_pos), 
            size = 4.5, 
            nudge_y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0)))

# Show plot
print(bar_plot)

# Export results to Excel
write.xlsx(filtered_results[, c("genus", "logFC", "p_value")], 
           file = "path_to/RESULTS/GG2/wilcoxon_significant_genus.xlsx", 
           row.names = FALSE)





################ Wilcoxon, Log FC plot -  species
##TSE from GG2
tse_GG2 <- loadFromQIIME2("path_to/RESULTS/GG2/counts.qza", taxonomy = "path_to/RESULTS/GG2/taxonomy.qza", sampleMetaFile="path_to/RESULTS/GG2/metadata.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))


library(phyloseq)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Load and preprocess data
tse <- tse_GG2
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]
rownames(tse) <- gsub(" ", "_", rownames(tse))

# Transform to relative abundance
tse_rel <- transformAssay(tse,
                          assay.type = "counts",
                          method = "relabundance")

# Convert to phyloseq and agglomerate at species level
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse_rel)
pseq_species <- tax_glom(pseq, taxrank = "species")

# Extract abundance and metadata
otu_mat <- as.data.frame(otu_table(pseq_species))
taxa <- tax_table(pseq_species)[, "species"]
sample_metadata <- sample_data(pseq_species)

# Ensure groupings are factors
group <- factor(sample_metadata$Response.day78., levels = c("DC", "No_DC"))

# Wilcoxon test per species
wilcox_results <- apply(otu_mat, 1, function(x) {
  test <- wilcox.test(x ~ group, exact = FALSE)
  logFC <- log2((mean(x[group == "DC"]) + 0.01) / (mean(x[group == "No_DC"]) + 0.01))
  c(p_value = test$p.value, logFC = logFC)
})

# Convert to data frame
wilcox_df <- as.data.frame(t(wilcox_results))
wilcox_df$taxon <- rownames(wilcox_df)
wilcox_df$species <- tax_table(pseq_species)[wilcox_df$taxon, "species"]

# Adjust p-values
wilcox_df$adj_p_value <- p.adjust(wilcox_df$p_value, method = "fdr")

# Filter results
filtered_results <- wilcox_df %>%
  filter(p_value < 0.05 & (logFC > 2 | logFC < -2)) %>%
  mutate(response = ifelse(logFC > 0, "DC", "No_DC"),
         highlight = ifelse(logFC > 0, "darkturquoise", "pink"),
         hjust_pos = ifelse(logFC > 0, -0.1, 1.1),
         species = gsub("^s__", "", species)) %>%
  filter(!is.na(species) & species != "")

# Bar plot
bar_plot <- ggplot(filtered_results, aes(x = reorder(species, logFC), 
                                         y = logFC, 
                                         fill = response)) +
  geom_bar(stat = "identity", show.legend = TRUE) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +  
  coord_flip() +
  scale_fill_manual(values = c("DC" = "#F6B9A9", "No_DC" = "#DCDCDC")) +  
  labs(title = "Differential abundant species between DC and No_DC",
       x = "", y = "Log Fold Change") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "plain", size = 34, hjust = 0.5),
    axis.title.x = element_text(size = 24, margin = margin(t = 10)),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 22),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(title = "", labels = c("DC", "No_DC"))) +
  geom_text(aes(label = species, hjust = hjust_pos), 
            size = 4.5, 
            nudge_y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0)))

# Show plot
print(bar_plot)

# Export results to Excel
write.xlsx(filtered_results[, c("species", "logFC", "p_value")], 
           file = "path_to/RESULTS/GG2/wilcoxon_significant_species.xlsx", 
           row.names = FALSE)





######## Wilcoxon bar plots

library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(gridExtra)  
library(grid)       

### genus

# Isolate Day 1 samples with a known disease outcome
tse <- tse_GG2
#tse <- tse_Metaphlan

tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Extract genus level
tse_genus <- mergeFeaturesByRank(tse, rank ="genus", onRankOnly=TRUE)

# Replace all spaces with underscores in row names of tse_genus
rownames(tse_genus) <- gsub(" ", "_", rownames(tse_genus))
rownames(tse_genus) <- gsub("[^a-zA-Z0-9_]", "_", rownames(tse_genus))

# Manually transform counts to relative abundance
counts_data <- assay(tse_genus)

# Calculate the total counts for each sample (sum across rows)
total_counts_per_sample <- rowSums(counts_data)

# Calculate the relative abundance for each genus by dividing by the total counts for each sample
rel_abundance_data <- counts_data / total_counts_per_sample

# Check if the relative abundance sums to 1 for each sample
row_sums <- rowSums(rel_abundance_data)
summary(row_sums)  

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_genus))

# Merge metadata with the transposed relative abundance data
merged_data <- merge(metadata_df, as.data.frame(t(rel_abundance_data)), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# DC first in boxplot
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Filter genus columns that start with "g__Alistipes"
Alistipes_columns <- grep("^g__Alistipes", names(merged_data), value = TRUE)

# Loop over only Alistipes genus columns and plot in RStudio
for (genus in Alistipes_columns) {
  # Perform Wilcoxon rank sum test
  wilcox_result <- wilcox.test(as.formula(paste(genus, "~ Response.day78.")),
                               data = merged_data,
                               subset = merged_data$Response.day78. %in% c("DC", "No_DC"),
                               exact = FALSE)
  
  # Plot using ggplot2
  boxplot <- ggplot(merged_data, aes_string(x = "Response.day78.", y = genus)) +
    geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
    geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1, height = 0), alpha = 0.7, size = 2) +
    geom_text(
      data = subset(merged_data, Patient == "20103(PR)"),
      aes(x = Response.day78., y = .data[[genus]], label = Patient),
      hjust = -0.2,
      vjust = 1.5, 
      size = 5,
      show.legend = FALSE
    ) +
    labs(x = NULL, y = "Relative abundance (%)", title = expression(italic("Alistipes") * " increased in DC")) +
    scale_fill_manual(values = c("DC" = "#FADFA0", "No_DC" = "#DCDCDC")) +
    scale_color_manual(values = c("DC" = "orange", "No_DC" = "darkgrey")) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme(
      plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5) , color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6) 
    )
  
  # Add statistics to the boxplot
  boxplot <- boxplot +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      #map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 4
    )
   
  # Print the plot in RStudio
  print(boxplot)
}

## Patchwork 3 plots
library(ggplot2)
library(ggsignif)
library(patchwork)

# Isolate Day 1 samples with a known disease outcome
tse <- tse_GG2
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Extract genus level
tse_genus <- mergeFeaturesByRank(tse, rank ="genus", onRankOnly=TRUE)

# Replace all spaces with underscores in row names of tse_genus
rownames(tse_genus) <- gsub(" ", "_", rownames(tse_genus))
rownames(tse_genus) <- gsub("[^a-zA-Z0-9_]", "_", rownames(tse_genus))

# Manually transform counts to relative abundance
counts_data <- assay(tse_genus)

# Calculate the total counts for each sample (sum across rows)
total_counts_per_sample <- rowSums(counts_data)

# Calculate the relative abundance for each genus by dividing by the total counts for each sample
rel_abundance_data <- counts_data / total_counts_per_sample

# Check if the relative abundance sums to 1 for each sample
row_sums <- rowSums(rel_abundance_data)
summary(row_sums)  

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_genus))

# Merge metadata with the transposed relative abundance data
merged_data <- merge(metadata_df, as.data.frame(t(rel_abundance_data)), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# DC first in boxplot
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Combine Alistipes and Eggerthella columns
#genus_columns <- grep("^g__Alistipes|^g__Eggerthella", names(merged_data), value = TRUE)
genus_columns <- grep("g__Alistipes_A_871404|g__Faecalibacterium|^g__Odoribacter|^g__Coprococcus_A_187", names(merged_data), value = TRUE)

# Store individual plots
plot_list <- list()

# Loop over genus columns
for (genus in genus_columns) {
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(
    as.formula(paste(genus, "~ Response.day78.")),
    data = merged_data,
    subset = merged_data$Response.day78. %in% c("DC", "No_DC"),
    exact = FALSE
  )
  
  # Build plot
  boxplot <- ggplot(merged_data, aes_string(x = "Response.day78.", y = genus)) +
    geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
    geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1), alpha = 0.7, size = 2) +
    geom_text(
      data = subset(merged_data, Patient == "20103(PR)"),
      aes(x = Response.day78., y = .data[[genus]], label = Patient),
      hjust = -0.2,
      vjust = 1.5,
      size = 5,
      show.legend = FALSE
    ) +
    labs(
      x = NULL,
      y = "Relative abundance (%)",
      title = bquote(italic(.(gsub("^g__", "", genus))) ~ "")
    ) +
    scale_fill_manual(values = c("DC" = "#FADFA0", "No_DC" = "#DCDCDC")) +
    scale_color_manual(values = c("DC" = "orange", "No_DC" = "darkgrey")) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme(
      plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)
    ) +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 4
    )
  
  # Store each plot
  plot_list[[genus]] <- boxplot
}

# Combine plots using patchwork (adjust spacing if needed)
combined_plot <- wrap_plots(plot_list, nrow = 1)

# Show final plot
print(combined_plot)




## Alistipes sum of both
# Combine Alistipes columns into one summed column
merged_data$Alistipes_total <- rowSums(merged_data[, Alistipes_columns], na.rm = TRUE)

# Wilcoxon test on combined Alistipes total
wilcox_result <- wilcox.test(Alistipes_total ~ Response.day78.,
                             data = merged_data,
                             exact = FALSE)

# Create the boxplot
Alistipes_plot <- ggplot(merged_data, aes(x = Response.day78., y = Alistipes_total)) +
  geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
  geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1), alpha = 0.7, size = 2) +
  geom_text(
    data = subset(merged_data, Patient == "20103(PR)"),
    aes(x = Response.day78., y = Alistipes_total, label = Patient),
    hjust = -0.2, vjust = 1.5, size = 5, show.legend = FALSE
  ) +
  labs(x = NULL, y = "Relative abundance (%)", title = expression(italic("Alistipes") * " increased in DC")) +
  scale_fill_manual(values = c("DC" = "#FADFA0", "No_DC" = "#DCDCDC")) +
  scale_color_manual(values = c("DC" = "orange", "No_DC" = "darkgrey")) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(
    plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
    axis.text.y = element_text(size = 16),
    legend.position = "none",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.6)
  ) +
  geom_signif(
    test = "wilcox.test",
    comparisons = list(c("DC", "No_DC")),
    map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
    step_increase = 0.1,
    color = "black",
    size = 0.7,
    textsize = 4
  )

# Print the plot
print(Alistipes_plot)



### GG2 species

# Isolate Day 1 samples with a known disease outcome
tse <- tse_GG2
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

# Replace all spaces with underscores in row names of tse_species
rownames(tse_species) <- gsub(" ", "_", rownames(tse_species))
rownames(tse_species) <- gsub("[^a-zA-Z0-9_]", "_", rownames(tse_species))

# Manually transform counts to relative abundance
counts_data <- assay(tse_species)

# Calculate the total counts for each sample (sum across rows)
total_counts_per_sample <- rowSums(counts_data)

# Calculate the relative abundance for each species by dividing by the total counts for each sample
rel_abundance_data <- counts_data / total_counts_per_sample

# Check if the relative abundance sums to 1 for each sample
row_sums <- rowSums(rel_abundance_data)
summary(row_sums) 

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_species))

# Merge metadata with the transposed relative abundance data
merged_data <- merge(metadata_df, as.data.frame(t(rel_abundance_data)), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# DC first in boxplot
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Filter species columns that start with "s__Eggerthella_lenta"
# Filter species columns that start with "s__Alistipes_A_871400_sp002362235"
alistipes_columns <- grep("^s__Eggerthella_lenta", names(merged_data), value = TRUE)

# Loop over only Alistipes species columns and plot in RStudio
for (species in alistipes_columns) {
  # Perform Wilcoxon rank sum test
  wilcox_result <- wilcox.test(as.formula(paste(species, "~ Response.day78.")),
                               data = merged_data,
                               subset = merged_data$Response.day78. %in% c("DC", "No_DC"),
                               exact = FALSE)
  
  # Plot using ggplot2
  boxplot <- ggplot(merged_data, aes_string(x = "Response.day78.", y = species)) +
    geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
    geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1, height = 0), alpha = 0.7, size = 2) +
    geom_text(
      data = subset(merged_data, Patient == "20103(PR)"),
      aes(x = Response.day78., y = .data[[species]], label = Patient),
      hjust = -0.2, 
      vjust = 1.5, 
      size = 5,
      show.legend = FALSE
    ) +
    labs(
      x = NULL,
      y = "Relative abundance (%)",
      title = expression(italic("Eggerthella lenta"))
    ) +  
    scale_fill_manual(values = c("DC" = "#F6B9A9", "No_DC" = "#DCDCDC")) +
    scale_color_manual(values = c("DC" = "red", "No_DC" = "darkgrey")) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme(
      plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5) , color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)   
    )
  
  # Add statistics to the boxplot
  boxplot <- boxplot +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      #map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 4
    )
  
  # Print the plot in RStudio
  print(boxplot)
}


# Patchwork

# Isolate Day 1 samples with a known disease outcome
## Patchwork, barplots species
# Isolate Day 1 samples with a known disease outcome
tse <- tse_GG2
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

# Replace all spaces with underscores in row names of tse_species
rownames(tse_species) <- gsub(" ", "_", rownames(tse_species))
rownames(tse_species) <- gsub("[^a-zA-Z0-9_]", "_", rownames(tse_species))

# Manually transform counts to relative abundance
counts_data <- assay(tse_species)

# Calculate the total counts for each sample (sum across rows)
total_counts_per_sample <- rowSums(counts_data)

# Calculate the relative abundance for each species by dividing by the total counts for each sample
rel_abundance_data <- counts_data / total_counts_per_sample

# Check if the relative abundance sums to 1 for each sample
row_sums <- rowSums(rel_abundance_data)
summary(row_sums)  

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_species))

# Merge metadata with the transposed relative abundance data
merged_data <- merge(metadata_df, as.data.frame(t(rel_abundance_data)), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# DC first in boxplot
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Define the desired order of species
ordered_species_columns <- c(
  "s__Alistipes_A_871400_excrementavium",
  "s__Alistipes_A_871400_shahii",
  "s__Eggerthella_lenta"
)

# Define corresponding pretty titles
species_map <- list(
  "s__Alistipes_A_871400_excrementavium" = expression(italic("Alistipes excrementavium")),
  "s__Alistipes_A_871400_shahii" = expression(italic("Alistipes shahii")),
  "s__Eggerthella_lenta" = expression(italic("Eggerthella lenta"))
)

# Generate and collect plots
plot_list <- list()
for (species_col in ordered_species_columns) {
  boxplot <- ggplot(merged_data, aes_string(x = "Response.day78.", y = species_col)) +
    geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
    geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1), alpha = 0.7, size = 2) +
    geom_text(
      data = subset(merged_data, Patient == "20103(PR)"),
      aes(x = Response.day78., y = .data[[species_col]], label = Patient),
      hjust = -0.2,
      vjust = 1.5,
      size = 5,
      show.legend = FALSE
    ) +
    labs(
      x = NULL,
      y = "Relative abundance (%)",
      title = species_map[[species_col]]
    ) +
    scale_fill_manual(values = c("DC" = "#F6B9A9", "No_DC" = "#DCDCDC")) +
    scale_color_manual(values = c("DC" = "red", "No_DC" = "darkgrey")) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme(
      plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)
    ) +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 4
    )
  
  plot_list[[species_col]] <- boxplot
}

# Combine with patchwork
library(patchwork)
combined_patchwork <- wrap_plots(plot_list, nrow = 1)

# Print
print(combined_patchwork)



#### LinDA differentially abundance analysis

# Define tse - GG2 - genus
tse <- tse_GG2


tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Extract genus level
tse_spec <- mergeFeaturesByRank(tse, rank ="genus", onRankOnly=TRUE)

# Replace all spaces with underscores in row names of tse_spec
rownames(tse_spec) <- gsub(" ", "_", rownames(tse_spec))

# Transform count assay to relative abundances
tse_rel <- transformAssay(tse_spec,
                          assay.type = "counts",
                          method = "relabundance")

# Run LinDA
linda_out <- linda(feature.dat = as.data.frame(assay(tse_rel)),
                   meta.dat = as.data.frame(colData(tse_rel)),
                   formula = "~ Response.day78.",
                   alpha = 0.05,
                   prev.filter = 0,
                   mean.abund.filter = 0)

head(linda_out)

# Extract the results data frame from linda_out
results <- linda_out$output$Response.day78.No_DC

# Convert row names to a column named 'taxon'
results$taxon <- rownames(results)

# Remove the "g__" prefix from the 'taxon' column in results
results$taxon <- gsub("^g__", "", results$taxon)

# Add a column to indicate if the taxon should be highlighted
results$highlight <- ifelse(
  results$taxon %in% c("Alistipes_A_871400", "Eggerthella"),
  "highlight",
  "no_highlight"
)

head(results)

# Filter for significant results based on adjusted p-value threshold (adjust as needed)
#significant_results <- results[results$reject == TRUE, ]

#significant_results

# Filter for significant results based on adjusted p-value threshold (adjust as needed)
significant_results <- results[results$pvalue < 0.05, ]
# Filter for significant results based on log2FoldChange thresholds
significant_results <- significant_results[significant_results$log2FoldChange < 0 | significant_results$log2FoldChange > 3,]


# Create the bar plot
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log2FoldChange), y = log2FoldChange, fill = highlight)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("highlight" = "blue", "no_highlight" = "grey")) +
  coord_flip() +  # Flip coordinates to make taxa names readable
  labs(title = "Differentially expressed genus (LinDA analysis)", x = "Taxon", y = "Log Fold Change") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "plain", size = 18, hjust = 0.5, margin = margin(b = 10)),
    plot.title.position = "plot",
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
    axis.text.y = element_text(size = 16),
    legend.position = "none",
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.6)
  ) 

# Print the plot
print(bar_plot)

# Export results to Excel
write.xlsx(significant_results[, c("genus", "logFC", "p_value")], 
           file = "path_to/RESULTS/GG2/LinDA_significant_genus.xlsx", 
           row.names = FALSE)



##################################################################################################################
################## Network analysis

# Network modelling using mia (https://microbiome.github.io/OMA/docs/devel/pages/network_learning.html)

#::install_github("microbiome/OMA", dependencies = TRUE, upgrade = TRUE)
install.packages("parallelly")

# Load required packages
library(devtools)
library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)


#install.packages("devtools")
library(devtools)
if(!require(SpiecEasi)){
  devtools::install_github("zdk123/SpiecEasi")
}

if(!require(SPRING)){
  devtools::install_github("GraceYoon/SPRING")
}


# Required packages
#install.packages("devtools")
library(devtools)
#install.packages("BiocManager")
library(BiocManager)

# Install NetCoMi
#devtools::install_github("stefpeschel/NetCoMi", 
#                         dependencies = c("Depends", "Imports", "LinkingTo"),
#                         repos = c("https://cloud.r-project.org/",
#                                   BiocManager::repositories()))


##TSE from GG2
# Data prep
tse <- tse_GG2
tse <- tse[, which(colData(tse)$Response.day78. %in% c("DC", "No_DC"))]
tse <- tse[, which(colData(tse)$Sample_point %in% c("Day_1"))]

# Agglomerate to genus level
tse <- agglomerateByRank(tse, rank = "genus")

# Replace all spaces with underscores in row names of tse
rownames(tse) <- gsub(" ", "_", rownames(tse))
rownames(tse) <- gsub("[^a-zA-Z0-9_]", "_", rownames(tse))

# Manually transform counts to relative abundance
counts_data <- assay(tse)

# Calculate the total counts for each sample (sum across rows)
total_counts_per_sample <- rowSums(counts_data)

# Calculate the relative abundance for each genus by dividing by the total counts for each sample
rel_abundance_data <- counts_data / total_counts_per_sample

# Check if the relative abundance sums to 1 for each sample
row_sums <- rowSums(rel_abundance_data)
summary(row_sums) 

assay(tse, "relabundance") <- rel_abundance_data

# Include only genus present in at least 4 or more patients
relabund_data <- assay(tse, "relabundance")
# Make sure there are no NAs
relabund_data[is.na(relabund_data)] <- 0
genus_to_keep <- rowSums(relabund_data > 0) >= 4
tse <- tse[genus_to_keep, ]
dim(tse)

#DC
DC <- tse[, rownames(subset(colData(tse), Response.day78. == "DC"))]
tse_DC_genus <- mergeFeaturesByRank(DC, rank ="genus", onRankOnly=TRUE)
top_DC <- getTopFeatures(tse_DC_genus, top = 50)
head(top_DC)

#No_DC
No_DC <- tse[, rownames(subset(colData(tse), Response.day78. == "No_DC"))]
tse_No_DC_genus <- mergeFeaturesByRank(No_DC, rank ="genus", onRankOnly=TRUE)
top_No_DC <- getTopFeatures(tse_No_DC_genus, top = 50)
head(top_No_DC)

# Renaming the "genus" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(tse)$genus,
                         function(x){if (x %in% top_DC || x %in% top_No_DC) {x} else {"Other"}})

rowData(tse)$genus <- as.character(genus_renamed)

# Filter out rows where taxonomic Response.day78. is "other"
tse <- tse[rowData(tse)$genus != "Other", ]
dim(tse)


# Filter by prevalence (on relative abundance)
tse <- subsetByPrevalent(
  tse,
  prevalence = 0.2,
  detection = 0,
  assay.type = "relabundance")

# Add log10-transformed abundances (on relative abundance)
tse <- transformAssay(tse, method = "log10", pseudocount = 1, assay.type = "relabundance")

# Add clr-transformed abundances (on relative abundance)
tse <- transformAssay(tse, method = "clr", pseudocount = 1, assay.type = "relabundance")

dim(tse)

# Remove the "g__" from the genus names
data <- rowData(tse)
data$genus <- gsub("^g__", "", data$genus)

# Assign the modified rowData back to the tse object
rowData(tse) <- data
rowData(tse)

# SPRING
library(SPRING)

# Using relative abundance
set.seed(13075)
spring_est <- SPRING(
  t(assay(tse, "relabundance")),  
  Rmethod = "approx",
  thresh = 0.5,
  lambdaseq = "data-specific"
)

# Get index of the optimal lambda selected by StARS
opt.K <- spring_est$output$stars$opt.index

# Store partial correlation matrix belonging to the optimal lambda as matrix
spring_cor <- SpiecEasi::symBeta(as.matrix(spring_est$output$est$beta[[opt.K]]))
spring_cor <- as.matrix(spring_cor)
rownames(spring_cor) <- colnames(spring_cor) <- rownames(tse)
diag(spring_cor) <- 1

# Network analysis with relative abundance
transform_asso <- function(assoMat, thresh = NULL, dissTrans = "unsigned") {
  # Sparsification
  if (!is.null(thresh)) {
    assoMat[abs(assoMat) < thresh] <- 0
  }
  
  # Compute dissimilarity matrix
  if (dissTrans == "signed") {
    dissMat <- sqrt(0.5 * (1 - assoMat))
  } else {
    dissMat <- sqrt(1 - assoMat^2)
  }
  
  # Dissimilarity between nodes with zero correlation is set to 1
  # (these nodes are unconnected and thus should have maximum dissimilarity)
  dissMat[assoMat == 0] <- 1
  
  # Compute similarity matrix
  simMat <- 1 - dissMat
  
  # Turn into igraph object
  graphObj <- SpiecEasi::adj2igraph(simMat)
  
  return(list(graph = graphObj, adja = simMat, asso = assoMat, diss = dissMat))
}

# Create graph object using relative abundance
spring_graph <- transform_asso(spring_cor)$graph

# NetCoMi
library(SpiecEasi)
library(NetCoMi)
library(SPRING)

# Remove the "g__" prefix from the row names in the 'relabundance' assay data
# Extract the relative abundance assay data
assay_data <- assay(tse, "relabundance")  

# Remove the "g__" prefix from the genus names
new_rownames <- gsub("^g__", "", rownames(assay_data))

# Update the row names in the assay data
rownames(assay_data) <- new_rownames
rownames(assay_data)

# Assign the modified assay data back to the tse object (without changing the dimensions)
assay(tse, "relabundance", withDimnames = FALSE) <- assay_data  

# Update the row names of the 'tse' object itself
rownames(tse) <- new_rownames

# Check if the update was successful
assay(tse, "relabundance")
column_sums <- colSums(assay(tse, "relabundance"))

# Print the summary to ensure sums are approximately 1
summary(column_sums)  

netcomi_net <- netConstruct(
  t(assay(tse, "relabundance")),  
  taxRank = "genus",
  filtTax = "numbSamp",
  filtTaxPar = list(numbSamp = 0.4),
  measure = "spring",
  measurePar = list(thresh = 0.9, Rmethod = "approx"),
  sparsMethod = "threshold",
  dissFunc = "signed",
  seed = 13075
)

edgelist <- netcomi_net$edgelist1[
  order(netcomi_net$edgelist1$adja, decreasing = TRUE), ]
head(edgelist)
tail(edgelist)

# Recreate adjacency matrix with signs
signed_adja <- netcomi_net$assoMat1  

# Threshold: Keep only strong positive or strong negative edges
positive_thresh <- 0.2
negative_thresh <- -0.2
signed_adja[signed_adja < positive_thresh & signed_adja > negative_thresh] <- 0

# Update NetCoMi network object with thresholded data
netcomi_net$assoMat1 <- signed_adja
netcomi_net$adjaMat1 <- abs(signed_adja)  # Required internally for plotting and clustering

# Convert to igraph object if you want to use outside of NetCoMi
netcomi_graph <- SpiecEasi::adj2igraph(signed_adja)

# Analyze the updated network
netcomi_netprops <- netAnalyze(
  netcomi_net,
  clustMethod = "cluster_fast_greedy",
  hubPar = "eigenvector",
  normDeg = FALSE
)

summary(netcomi_netprops, numbNodes = 5)

# check: Largest connected component (LCC), and #singletons
# check: # clusters and # nodes
# check: # hub nodes (highest eigenvector centrality)

dev.off()

# include edge filter
plot(netcomi_netprops,
     repulsion = 0.98,
     rmSingles = TRUE,
     shortenLabels = "none",
     labelScale = FALSE,
     nodeSize = "eigenvector",
     nodeSizeSpread = 3,
     nodeColor = "cluster",
     hubBorderCol = "gray40",
     cexNodes = 1,
     cexLabels = 0.9,
     labelRepel = TRUE, 
     edgeTranspHigh = 20,
     title1 = "Network of top genus in all patients",
     showTitle = TRUE,
     cexTitle = 2.1,
     mar = c(1, 3, 3, 8))

# Add legends
legend("topright",                 # Legend position
       inset = c(0.05, 0.15),       # Move legend higher (increase vertical inset)
       cex = 0.8,                  # Legend text size
       title = "Estimated correlation:",
       title.cex = 1.1,
       legend = c("+", "-"), 
       lty = 1, 
       lwd = 4, 
       col = c("#009900", "red"),
       bty = "n", 
       horiz = TRUE)


## Write Excel file
library(openxlsx)  # make sure this is loaded for createWorkbook etc.

# Get adjacency matrix (signed, thresholded) from your updated network
signed_mat <- netcomi_net$assoMat1

# Extract edges from adjacency matrix (non-zero edges)
edges <- which(signed_mat != 0, arr.ind = TRUE)

# Create edge list with bacteria pairs and correlation signs
edge_export <- data.frame(
  Genus_1 = rownames(signed_mat)[edges[,1]],
  Genus_2 = colnames(signed_mat)[edges[,2]],
  Correlation_sign = ifelse(signed_mat[edges] > 0, "+", "-")
)

# Remove duplicate edges (adjacency matrix is symmetric)
edge_export <- edge_export[edge_export$Genus_1 < edge_export$Genus_2, ]

# Check the first few rows
head(edge_export)

# Write to Excel
wb <- createWorkbook()
addWorksheet(wb, "Edges")
writeData(wb, "Edges", edge_export)
saveWorkbook(wb, "path_to/RESULTS/GG2/network_edge_list.xlsx", overwrite = TRUE)


#####################################################################################################

# Dot plot - KEGG
# Gene family data
gene_family_data <- read.delim("path_to/RESULTS/humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)

#Set column 1 as rownames
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Function to extract the desired part of the column names
extract_desired_part <- function(column_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(column_name, "000[^_]*")
}

# Apply the function to all column names of the dataframe
new_colnames <- sapply(colnames(gene_family_data), extract_desired_part)

# Set the new column names to the dataframe
colnames(gene_family_data) <- new_colnames

# Print the new column names to verify
print(colnames(gene_family_data))


## Metadata
metadata <- read.table("path_to/RESULTS/humann3/metadata.TSV", header = TRUE, sep = "\t", row.names = 1, check.names=FALSE)

# Function to extract the desired part of the row names
extract_desired_part <- function(row_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(row_name, "000[^_]*")
}

# Apply the function to all row names of the dataframe
new_rownames <- sapply(rownames(metadata), extract_desired_part)

# Set the new row names to the dataframe
rownames(metadata) <- new_rownames

# Print the new row names to verify
print(rownames(metadata))


# Ensure that both identifiers are sorted and aligned
# Get common identifiers
common_ids <- intersect(colnames(gene_family_data), rownames(metadata))

# Filter and reorder columns and rows to match the common identifiers
gene_family_data <- gene_family_data[, common_ids,drop=FALSE]
metadata <- metadata[common_ids, , drop = FALSE]

# Ensure that the row names of the metadata match the column names of the count data
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# Create a TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(
  assays = list(counts = (gene_family_data)),
  colData = metadata
)

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$'Response(day78)' %in% c("DC", "No_DC"))]


# check the TSE object
print(tse)


# Subset TSE object (already done above)
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$'Response(day78)' %in% c("DC", "No_DC"))]

# Extract count data and metadata
counts <- assay(tse)
meta <- colData(tse)

# Ensure group labels
group <- factor(meta$`Response(day78)`, levels = c("DC", "No_DC"))

# Perform Wilcoxon test for each feature
wilcox_results <- apply(counts, 1, function(feature_counts) {
  test <- wilcox.test(feature_counts ~ group)
  return(c(statistic = test$statistic, pval = test$p.value))
})

# Convert to data frame
wilcox_df <- as.data.frame(t(wilcox_results))
wilcox_df$Feature <- rownames(wilcox_df)
wilcox_df$padj <- p.adjust(wilcox_df$pval, method = "fdr")

# Filter significant features
significant_results <- wilcox_df %>% filter(pval < 0.05)

# Calculate Gene Ratio
abundance_matrix_subset <- counts[significant_results$Feature, , drop = FALSE]
significant_results$GeneRatio <- rowSums(abundance_matrix_subset > 0) / ncol(abundance_matrix_subset)

# Calculate logFC (log2 of median ratio between groups)
logFC <- apply(counts[significant_results$Feature, , drop = FALSE], 1, function(feature_counts) {
  group1_median <- median(feature_counts[group == "DC"])
  group2_median <- median(feature_counts[group == "No_DC"])
  log2((group1_median + 1e-6) / (group2_median + 1e-6))  # Add small constant to avoid log(0)
})

# Add logFC to the results table
significant_results$logFC <- logFC

# Format for plotting
plot_data <- significant_results %>%
  mutate(Feature = gsub("\\.", " ", Feature)) %>%
  arrange(pval) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

# Plot
ggplot(plot_data, aes(x = GeneRatio, y = Feature)) +
  geom_point(aes(color = pval), size = 4) +
  scale_color_gradient(low = "#2E8B57", high = "skyblue") +
  theme_bw(base_size = 14) +
  labs(
    title = "Functional Enrichment of KEGG Features",
    x = "Gene Ratio",
    y = NULL,
    color = "p-value"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold")
  )

# Rename pval to p_value and write Excel
significant_results <- significant_results %>%
  rename(p_value = pval)

write.xlsx(significant_results[, c("Feature", "p_value", "logFC")],
           file = "path_to/RESULTS/GG2/wilcoxon_kegg_results.xlsx",
           row.names = FALSE)




## Heatmap Arginine
library(pheatmap)
library(stringr)

# KO
gene_family_data <- read.delim("path_to/RESULTS/humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)

#Set column 1 as rownames
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]


# Function to extract the desired part of the column names
extract_desired_part <- function(column_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(column_name, "000[^_]*")
}

# Apply the function to all column names of the dataframe
new_colnames <- sapply(colnames(gene_family_data), extract_desired_part)

# Set the new column names to the dataframe
colnames(gene_family_data) <- new_colnames

# Print the new column names to verify
print(colnames(gene_family_data))


## Metadata
metadata <- read.table("path_to/RESULTS/humann3/metadata.TSV", header = TRUE, sep = "\t", row.names = 1, check.names=FALSE)

# Function to extract the desired part of the row names
extract_desired_part <- function(row_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(row_name, "000[^_]*")
}

# Apply the function to all row names of the dataframe
new_rownames <- sapply(rownames(metadata), extract_desired_part)

# Set the new row names to the dataframe
rownames(metadata) <- new_rownames

# Print the new row names to verify
print(rownames(metadata))


# Ensure that both identifiers are sorted and aligned
# Get common identifiers
common_ids <- intersect(colnames(gene_family_data), rownames(metadata))

# Filter and reorder columns and rows to match the common identifiers
gene_family_data <- gene_family_data[, common_ids,drop=FALSE]
metadata <- metadata[common_ids, , drop = FALSE]

# Ensure that the row names of the metadata match the column names of the count data
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# Create a TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(
  assays = list(counts = (gene_family_data)),
  colData = metadata
)

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$'Response(day78)' %in% c("DC", "No_DC"))]
rownames(tse) <- gsub("[^[:alnum:]]", "_", rownames(tse))


# Check the TSE object
print(tse)

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse))

# Transpose assay data
tse_assay <- t(assay(tse))

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# DC first in boxplot
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Extract only the KO columns from merged_data
columns_of_interest <- grep("^K03758|^K01478", names(merged_data), value = TRUE)

# Confirm they exist
print(columns_of_interest)

# Extract those columns and force it into a data frame
heatmap_matrix <- as.data.frame(merged_data[, columns_of_interest, drop = FALSE])

# Assign rownames using Patient info
heatmap_matrix$Patient <- merged_data$Patient
heatmap_matrix <- column_to_rownames(heatmap_matrix, "Patient")

# Rename patients
name_mapping <- c(
  "20211" = "No_DC_20211",
  "20104" = "No_DC_20104",
  "30209" = "No_DC_30209",
  "20212" = "DC_20212",
  "20206" = "DC_20206",
  "20108" = "DC_20108",
  "30207" = "DC_30207",
  "20103(PR)" = "DC_20103(PR)"
)

# Apply renaming
matched_names <- intersect(rownames(heatmap_matrix), names(name_mapping))
rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% matched_names] <- name_mapping[matched_names]

# Reorder rows
ordered_patients <- c("No_DC_20211", "No_DC_20104", "No_DC_30209",
                      "DC_20212", "DC_20206", "DC_20108", "DC_30207", "DC_20103(PR)")

# Filter only available patients
ordered_patients <- ordered_patients[ordered_patients %in% rownames(heatmap_matrix)]
heatmap_matrix <- heatmap_matrix[ordered_patients, , drop = FALSE]

# Rename KO column names directly
colnames(heatmap_matrix) <- gsub(
  pattern = "K01478.*",
  replacement = "Arginine deiminase
  (K01478)",
  x = colnames(heatmap_matrix)
)

colnames(heatmap_matrix) <- gsub(
  pattern = "K03758.*",
  replacement = "Arginine ornithine antiporter 
  / lysine permease
  (K03758)",
  x = colnames(heatmap_matrix)
)

# Color palette
color_palette <- colorRampPalette(c("skyblue", "#2E8B57"))(100)

# Your heatmap plot
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = color_palette,
         fontsize_row = 12,
         fontsize_col = 20,
         angle_col = 45,
         main = "Arginine pathway abundance",
         border_color = "black",
         fontsize = 18,
         legend = TRUE)  


################ barplots 
# KO
gene_family_data <- read.delim("path_to/RESULTS/humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)

#Set column 1 as rownames
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Function to extract the desired part of the column names
extract_desired_part <- function(column_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(column_name, "000[^_]*")
}

# Apply the function to all column names of the dataframe
new_colnames <- sapply(colnames(gene_family_data), extract_desired_part)

# Set the new column names to the dataframe
colnames(gene_family_data) <- new_colnames

# Print the new column names to verify
print(colnames(gene_family_data))


## Metadata
metadata <- read.table("path_to/RESULTS/humann3/metadata.TSV", header = TRUE, sep = "\t", row.names = 1, check.names=FALSE)

# Function to extract the desired part of the row names
extract_desired_part <- function(row_name) {
  # Use regular expression to extract part starting with 000 and ending before _
  str_extract(row_name, "000[^_]*")
}

# Apply the function to all row names of the dataframe
new_rownames <- sapply(rownames(metadata), extract_desired_part)

# Set the new row names to the dataframe
rownames(metadata) <- new_rownames

# Print the new row names to verify
print(rownames(metadata))


# Ensure that both identifiers are sorted and aligned
# Get common identifiers
common_ids <- intersect(colnames(gene_family_data), rownames(metadata))

# Filter and reorder columns and rows to match the common identifiers
gene_family_data <- gene_family_data[, common_ids,drop=FALSE]
metadata <- metadata[common_ids, , drop = FALSE]

# Ensure that the row names of the metadata match the column names of the count data
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# Create a TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(
  assays = list(counts = (gene_family_data)),
  colData = metadata
)

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Sample_point == "Day_1")]
tse <- tse[, which(colData(tse)$'Response(day78)' %in% c("DC", "No_DC"))]
rownames(tse) <- gsub("[^[:alnum:]]", "_", rownames(tse))


# Check the TSE object
print(tse)

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse))

# Transpose assay data
tse_assay <- t(assay(tse))

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
merged_data <- merged_data %>%
  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# Ensure the correct factor order
merged_data$Response.day78. <- factor(merged_data$Response.day78., levels = c("DC", "No_DC"))

# Extract first matching column names for the KOs
K02106_col <- grep("^K02106", colnames(merged_data), value = TRUE)[1]
K01035_col <- grep("^K01035", colnames(merged_data), value = TRUE)[1]
K00929_col <- grep("^K00929", colnames(merged_data), value = TRUE)[1]

# Define colors
fill_colors <- c("DC" = "skyblue", "No_DC" = "#9DC183")
point_colors <- c("DC" = "darkblue", "No_DC" = "darkgreen")

# Function to create each boxplot
library(ggsignif)

make_gene_plot <- function(gene_col, plot_title) {
  # Prepare a clean dataframe with consistent column naming
  sub_df <- merged_data[, c("Patient", "Response.day78.", gene_col)]
  colnames(sub_df)[3] <- "Value"
  
  # Build the plot using the renamed 'Value' column
  boxplot <- ggplot(sub_df, aes(x = Response.day78., y = Value)) +
    geom_boxplot(aes(fill = Response.day78.), color = "black", width = 0.3, outlier.color = NA) +
    geom_point(aes(color = Response.day78.), position = position_jitter(width = 0.1), alpha = 0.7, size = 2) +
    geom_text(
      data = subset(sub_df, Patient == "20103(PR)"),
      aes(label = Patient),
      hjust = -0.2,
      vjust = 1.5,
      size = 5,
      show.legend = FALSE
    ) +
    labs(x = NULL, y = "Gene count", title = plot_title) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = point_colors) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)
    ) +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      y_position = max(sub_df$Value, na.rm = TRUE) * 1.1,  # Position slightly above max value
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 5
    )
  
  return(boxplot)
}


# Create the three individual plots
p1 <- make_gene_plot(K02106_col, "Short chain fatty acids\ntransporter (K02106)")
p2 <- make_gene_plot(K01035_col, "Acetate CoA (K01035)")
p3 <- make_gene_plot(K00929_col, "Butyrate kinase (K00929)")

# Combine using patchwork
combined_plot <- (p1 | plot_spacer() | p2 | plot_spacer() | p3) +
  plot_layout(widths = c(1, 0.1, 1, 0.1, 1))

# Print the combined plot
print(combined_plot)





## SCFA Serum analysis violin plots
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(patchwork)

# Load data
serum <- read_excel("path_to/SCFA/SCFA_serum")

# Define patient columns
patient_columns <- c("20211", "20104", "30209", "20212", "20206", "20108", "20103(PR)", "30209")

# Tidy the data
serum_clean <- serum %>%
  select(Metabolites, all_of(patient_columns)) %>%
  filter(Metabolites %in% c("Butyric-isobutyric Acid", "2-Methylvaleric Acid", "Propionic Acid")) %>%
  pivot_longer(-Metabolites, names_to = "Patient", values_to = "Value") %>%
  mutate(Response = ifelse(Patient %in% c("20211", "20104", "30209"), "No_DC", "DC"))

# Function to generate a violin plot
make_violin_plot <- function(met_name) {
  sub_df <- filter(serum_clean, Metabolites == met_name)
  
  ggplot(sub_df, aes(x = Response, y = Value)) +
    geom_violin(aes(fill = Response), color = "black", width = 0.3) +
    geom_point(aes(color = Response), position = position_jitter(width = 0.1), size = 2, alpha = 0.7) +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 5
    ) +
    geom_text(
      data = subset(sub_df, Patient == "20103(PR)"),
      aes(x = Response, y = Value, label = Patient),
      hjust = -0.2,
      vjust = 1.5,
      size = 5,
      show.legend = FALSE
    ) +
    labs(title = met_name, x = NULL, y = "Concentration (µM)") +
    scale_fill_manual(values = c("DC" = "skyblue", "No_DC" = "#9DC183")) +
    scale_color_manual(values = c("DC" = "darkblue", "No_DC" = "darkgreen")) +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)
    )
}

# Create the 3 plots
p1 <- make_violin_plot("2-Methylvaleric Acid")
p2 <- make_violin_plot("Butyric-isobutyric Acid")
p3 <- make_violin_plot("Propionic Acid")

# Combine using patchwork
combined_plot <- (p1 | plot_spacer() | p2 | plot_spacer() | p3) +
  plot_layout(widths = c(1, 0.1, 1, 0.1, 1))

# Show the plot
print(combined_plot)



#################### Eggerthella and Alistipes over time (30209, 30207)
##  Genus
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(tibble)
library(patchwork)


tse <- tse_GG2
tse <- tse[, colData(tse)$Patient %in% c("30209", "30207")]

# Convert to relative abundance
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Agglomerate at genus level
tse_genus <- agglomerateByRank(tse_rel, rank = "genus", onRankOnly = TRUE)

# Filter to specific genera
keep_genera <- c("g__Eggerthella", "g__Alistipes_A_871400")
tse_filtered <- tse_genus[rowData(tse_genus)$genus %in% keep_genera, ]

# Convert to data frame
df <- assay(tse_filtered, "relabundance") %>%
  as.data.frame() %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "sample", values_to = "abundance")

# add metadata
df <- df %>%
  dplyr::left_join(
    as.data.frame(colData(tse_filtered)) %>% 
      rownames_to_column("sample"),
    by = "sample"
  )

# Clean Sample_point
df$Sample_point <- gsub("_", " ", df$Sample_point)

# Ensure factors have the correct order
df$Sample_point <- factor(df$Sample_point, levels = c("Day 1", "Day 36"))

# Rename Patient values for plotting
df$Patient <- recode(df$Patient, "30207" = "30207(DC)", "30209" = "30209(No_DC)")

plot_genus <- function(genus_name) {
  df_sub <- df %>% filter(genus == genus_name)
  
  # Define x-axis positions for custom spacing
  df_sub <- df_sub %>%
    mutate(
      Sample_point = factor(Sample_point, levels = c("Day 1", "Day 36")),
      x_pos = case_when(
        Patient == "30207(DC)" & Sample_point == "Day 1" ~ 1,
        Patient == "30207(DC)" & Sample_point == "Day 36" ~ 2,
        Patient == "30209(No_DC)" & Sample_point == "Day 1" ~ 4,
        Patient == "30209(No_DC)" & Sample_point == "Day 36" ~ 5
      )
    )
  
  genus_label <- gsub("^g__", "", genus_name)
  
  ggplot(df_sub, aes(x = x_pos, y = abundance)) +
    geom_line(aes(group = Patient), color = "gray60", linewidth = 1) +
    geom_point(
      aes(shape = Sample_point, color = Sample_point),
      size = 3, stroke = 1.2
    ) +
    scale_color_manual(
      values = c("Day 1" = "black", "Day 36" = "#580F8F"),
      guide = guide_legend(title = "Sample point")
    ) +
    scale_shape_manual(
      values = c("Day 1" = 1, "Day 36" = 16),
      guide = guide_legend(title = "Sample point")
    ) +
    scale_x_continuous(
      breaks = c(1.5, 4.5),
      labels = c("30207\n(DC)", "30209\n(No_DC)"),
      expand = expansion(add = 0.5)
    ) +
    labs(title = genus_label, x = NULL, y = "Relative abundance (%)") +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.8)
    )
}

p1 <- plot_genus("g__Alistipes_A_871400") +
  labs(title = "Alistipes_A_871400\nincrease after therapy")

p2 <- plot_genus("g__Eggerthella") +
  labs(title = "Eggerthella\ndecrease after therapy")

# Combine plots with individual legends, no global title
(p1 | plot_spacer() | p2) +
  plot_layout(widths = c(1, 0.1, 1))


# species - Pie chart
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(mia)

# Filter and normalize
tse <- tse_GG2[, colData(tse_GG2)$Patient %in% c("30207")]
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")
tse_species <- agglomerateByRank(tse_rel, rank = "species", onRankOnly = TRUE)

# Filter for Alistipes species only
Alistipes_species <- rowData(tse_species)$species[
  grepl("^s__Alistipes", rowData(tse_species)$species)
]
tse_filtered <- tse_species[rowData(tse_species)$species %in% Alistipes_species, ]

# Extract metadata
meta_df <- as.data.frame(colData(tse_filtered))
meta_df$sample <- rownames(meta_df)

# Melt assay data and group "Other"
df_long <- assay(tse_filtered, "relabundance") %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance") %>%
  dplyr::left_join(meta_df, by = "sample")

# Clean up metadata
df_long$Sample_point <- gsub("_", " ", df_long$Sample_point)
df_long$Sample_point <- factor(df_long$Sample_point, levels = c("Day 1", "Day 36"))
df_long$Patient <- as.factor(df_long$Patient)

# Define color palette (keep your custom colors)
custom_colors <- c(
  "#deebf7", "#9ecae1", "#6baed6", "#3182bd",
  "#8856a7", "#6a51a3", "grey70", "#756bb1"
)

# Identify all species with total abundance > 0%
species_used <- df_long %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  filter(total_abundance > 0) %>%
  pull(species)

# Assign custom colors to species with >0% abundance
if (length(species_used) > length(custom_colors)) {
  warning("More species than colors! Repeating colors. You may want to expand 'custom_colors'.")
  custom_colors <- rep(custom_colors, length.out = length(species_used))
}
species_colors <- setNames(custom_colors[1:length(species_used)], species_used)

# Pie chart data prep function
prepare_pie_data <- function(df, patient, day) {
  df %>%
    filter(Patient == patient, Sample_point == day) %>%
    group_by(species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    mutate(
      fraction = abundance / sum(abundance),
      percent_label = paste0(round(fraction * 100, 1), "%")
    ) %>%
    arrange(desc(species)) %>%
    mutate(
      csum = rev(cumsum(rev(fraction))),
      pos = fraction / 2 + lead(csum, 1),
      pos = if_else(is.na(pos), fraction / 2, pos)
    )
}

# Clean species names and build labels
pie1_data <- prepare_pie_data(df_long, "30207", "Day 1") %>%
  filter(fraction > 0) %>%
  mutate(clean_species = gsub("^s__", "", species))

pie1_labels <- setNames(
  paste0(pie1_data$clean_species, " — ", pie1_data$percent_label),
  pie1_data$species
)

species_colors1 <- species_colors[names(species_colors) %in% pie1_data$species]

pie1 <- ggplot(pie1_data, aes(x = "", y = fraction, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(
    values = species_colors1,
    labels = pie1_labels,
    drop = FALSE
  ) +
  labs(title = "Day 1", fill = "Species") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

pie2_data <- prepare_pie_data(df_long, "30207", "Day 36") %>%
  filter(fraction > 0) %>%
  mutate(clean_species = gsub("^s__", "", species))

pie2_labels <- setNames(
  paste0(pie2_data$clean_species, " — ", pie2_data$percent_label),
  pie2_data$species
)

species_colors2 <- species_colors[names(species_colors) %in% pie2_data$species]

pie2 <- ggplot(pie2_data, aes(x = "", y = fraction, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(
    values = species_colors2,
    labels = pie2_labels,
    drop = FALSE
  ) +
  labs(title = "Day 36", fill = "Species") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# Combine with overall title
(pie1 + pie2) +
  plot_annotation(
    title = "Alistipes species in 30207 (DC)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )


# Plot without %
pie1_nolegend <- ggplot(pie1_data, aes(x = "", y = fraction, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(
    values = species_colors1,
    drop = FALSE  
  ) +
  labs(title = "Day 1") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Clean species names again (no % in legend)
pie2_labels_names_only <- setNames(
  pie2_data$clean_species,
  pie2_data$species
)

pie2_namesonly <- ggplot(pie2_data, aes(x = "", y = fraction, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(
    values = species_colors2,
    labels = pie2_labels_names_only,
    drop = FALSE
  ) +
  labs(title = "Day 36", fill = "Species") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

(pie1_nolegend + pie2_namesonly) +
  plot_annotation(
    title = "Alistipes species in 30207 (DC)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )


#####################################################################################################
############ Functional analysis over time

# Load required libraries
library(TreeSummarizedExperiment)
library(tidyverse)
library(stringr)

# Load gene family data
gene_family_data <- read.delim(
  "path_to/RESULTS/humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt",
  header = TRUE
)
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Extract column IDs
extract_desired_part <- function(name) {
  str_extract(name, "000[^_]*")
}
new_colnames <- sapply(colnames(gene_family_data), extract_desired_part)
if (any(is.na(new_colnames))) {
  warning("Some column names could not be parsed correctly. Check regular expression.")
}
colnames(gene_family_data) <- new_colnames

# Load metadata
metadata <- read.table(
  "path_to/RESULTS/humann3/metadata.TSV",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
new_rownames <- sapply(rownames(metadata), extract_desired_part)
if (any(is.na(new_rownames))) {
  warning("Some row names in metadata could not be parsed correctly. Check regular expression.")
}
rownames(metadata) <- new_rownames

# Align data
common_ids <- intersect(colnames(gene_family_data), rownames(metadata))
gene_family_data <- gene_family_data[, common_ids, drop = FALSE]
metadata <- metadata[common_ids, , drop = FALSE]
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# Create TSE object
tse <- TreeSummarizedExperiment(
  assays = list(counts = gene_family_data),
  colData = metadata
)
tse_selected <- tse[, colData(tse)$Patient %in% c("30207", "30209")]
tse_rel <- transformAssay(tse_selected, assay.type = "counts", method = "relabundance")

# Define KO IDs and titles
ko_list <- c("K02106", "K01035", "K00929")
titles <- c(
  "SCFA transporter\ndecreased after therapy",
  "Acetate CoA\ndecreased after therapy",
  "Butyrate kinase\nincreased after therapy"
)

# Plotting function
make_plot <- function(tse_rel, ko_id, title_text) {
  rows <- grep(paste0("^", ko_id), rownames(tse_rel), value = TRUE)
  if (length(rows) == 0) stop(paste("No rows starting with", ko_id, "found."))
  
  tse_filtered <- tse_rel[rows, ]
  
  # Build metadata table
  meta_df <- data.frame(sample = colnames(tse_filtered), as.data.frame(colData(tse_filtered)))
  
  # Build long format data
  df_long <- assay(tse_filtered, "relabundance") %>%
    as.data.frame() %>%
    rownames_to_column("species") %>%
    pivot_longer(-species, names_to = "sample", values_to = "abundance") %>%
    dplyr::left_join(meta_df, by = "sample")
  
  # Clean and order factors
  df_long$Sample_point <- gsub("_", " ", df_long$Sample_point)
  df_long$Sample_point <- factor(df_long$Sample_point, levels = c("Day 1", "Day 36"))
  df_long$Patient <- factor(df_long$Patient, levels = c("30207", "30209"))
  
  # Define x positions for grouped bars
  df_long <- df_long %>%
    mutate(x_pos = case_when(
      Patient == "30207" & Sample_point == "Day 1" ~ 1,
      Patient == "30207" & Sample_point == "Day 36" ~ 2,
      Patient == "30209" & Sample_point == "Day 1" ~ 4,
      Patient == "30209" & Sample_point == "Day 36" ~ 5
    ))
  
  # Set upper y-limit
  upper_limit <- max(df_long$abundance, na.rm = TRUE) * 1.1
  
  # Plot
  ggplot(df_long, aes(x = x_pos, y = abundance, fill = Sample_point)) +
    geom_col(width = 0.9, color = "black") +
    geom_point(shape = 21, size = 3.5, stroke = 1, color = "black") +
    scale_fill_manual(values = c("Day 1" = "#C3B1E1", "Day 36" = "#3f007d")) +
    scale_x_continuous(
      breaks = c(1.5, 4.5),
      labels = c("30207\n(DC)", "30209\n(No_DC)"),
      expand = expansion(add = 0.5)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, upper_limit)) +
    labs(
      title = title_text,
      x = "Patient",
      y = "Relative abundance (%)",
      fill = "Sample point"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8)
    )
}


plots <- list()
for (i in seq_along(ko_list)) {
  plots[[i]] <- make_plot(tse_rel, ko_list[i], titles[i])
}

# Combine with patchwork
(plots[[1]] | plot_spacer() | plots[[2]] | plot_spacer() | plots[[3]]) +
  plot_layout(widths = c(1, 0.1, 1, 0.1, 1), guides = "collect")


############################### Clinical correlations, BL ###################

library(ggplot2)
library(patchwork)
library(readxl)

# Load the data
clinic_fecal <- read_excel("path_to/Correlation data.xlsx", "Correlation_fecal")
clinic_blood <- read_excel("path_to/Correlation data.xlsx", "Correlation_blood_levels")

# Convert first column to rownames
clinic_fecal[[1]] <- as.character(clinic_fecal[[1]])
clinic_blood[[1]] <- as.character(clinic_blood[[1]])
clinic_fecal <- clinic_fecal[!duplicated(clinic_fecal[[1]]), ]
clinic_blood <- clinic_blood[!duplicated(clinic_blood[[1]]), ]
rownames(clinic_fecal) <- clinic_fecal[[1]]
rownames(clinic_blood) <- clinic_blood[[1]]
clinic_fecal <- clinic_fecal[, -1]
clinic_blood <- clinic_blood[, -1]

# Ensure common patients
common_patients <- intersect(rownames(clinic_fecal), rownames(clinic_blood))
clinic_fecal <- clinic_fecal[common_patients, ]
clinic_blood <- clinic_blood[common_patients, ]

# Variable mappings
plot_pairs <- list(
  "2-Methylvaleric Acid" = "Neutrophils",
  "Propionic Acid" = "Neutrophils"
)

rename_vars <- list(
  "2-Methylvaleric Acid (Day1)" = "2-Methylvaleric Acid",
  "Propionic Acid (Day1)" = "Propionic Acid",
  "Neutrophils_day1" = "Neutrophils"
)

# Plot storage
plot_list <- list()

# Generate plots
for (fecal_var_pretty in names(plot_pairs)) {
  blood_var_pretty <- plot_pairs[[fecal_var_pretty]]
  fecal_var <- names(rename_vars)[which(rename_vars == fecal_var_pretty)]
  blood_var <- names(rename_vars)[which(rename_vars == blood_var_pretty)]
  
  x <- as.numeric(clinic_fecal[[fecal_var]])
  y <- as.numeric(clinic_blood[[blood_var]])
  df <- data.frame(x = x, y = y)
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 3) next
  
  norm_x <- shapiro.test(df$x)$p.value > 0.05
  norm_y <- shapiro.test(df$y)$p.value > 0.05
  method <- ifelse(norm_x && norm_y, "pearson", "spearman")
  
  cor_result <- cor.test(df$x, df$y, method = method)
  cor_val <- signif(cor_result$estimate, 3)
  p_val <- signif(cor_result$p.value, 3)
  
  cor_symbol <- if (method == "pearson") "r" else "ρ"
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "darkgreen", alpha = 0.8, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", fill = "lightgrey") +
    labs(
      x = fecal_var_pretty,
      y = blood_var_pretty,
      title = paste(fecal_var_pretty, "vs", blood_var_pretty),
      subtitle = paste("Method:", method, "|", cor_symbol, "=", cor_val, "| p =", ifelse(p_val < 0.001, "< 0.001", p_val))
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.6)
    )
  
  plot_list[[fecal_var_pretty]] <- p
}

# Display side-by-side with patchwork
if (length(plot_list) == 2) {
  combined_plot <- plot_list[[1]] | plot_list[[2]]
  print(combined_plot)
}




### Clinical relations, Day1 vs Day36, paired barplot ###################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Define renaming and order
metric_rename <- c(
  "Neutrophils" = "Neutrophils",
  "Lymphocyte_count" = "Lymphocytes",
  "Eosinophils" = "Eosinophils",
  "Monocyte_count" = "Monocytes"
)
metric_order <- c("Neutrophils", "Lymphocytes", "Eosinophils", "Monocytes")

# Function to generate log2 FC plot for one patient
make_fc_plot <- function(patient_id, title_text, fill_color) {
  clinic_blood <- read_excel(
    "path_to/Correlation data.xlsx", 
    sheet = "Correlation_blood_levels"
  )
  
  blood_day1 <- c("Neutrophils_day1", "Eosinophils_day1", 
                  "Lymphocyte_count_day1", "Monocyte_count_day1")
  blood_day36 <- gsub("day1", "day36", blood_day1)
  
  clinic_blood <- clinic_blood %>%
    rename(Patient = 1) %>%
    filter(Patient == patient_id) %>%
    select(Patient, all_of(blood_day1), all_of(blood_day36)) %>%
    mutate(across(-Patient, ~ as.numeric(gsub(",", ".", .))))
  
  fold_change_data <- clinic_blood %>%
    rowwise() %>%
    mutate(across(all_of(blood_day1), .names = "log2FC_{.col}", 
                  ~ {
                    metric_name <- gsub("_day1", "", cur_column())
                    day36_col <- paste0(metric_name, "_day36")
                    val_day1 <- .x
                    val_day36 <- get(day36_col)
                    ifelse(val_day1 > 0, log2(val_day36 / val_day1), NA)
                  })) %>%
    select(Patient, starts_with("log2FC_")) %>%
    pivot_longer(cols = -Patient, names_to = "Metric", values_to = "log2_FC") %>%
    mutate(Metric = gsub("log2FC_|_day1", "", Metric)) %>%
    mutate(Metric = recode(Metric, !!!metric_rename)) %>%
    mutate(Metric = factor(Metric, levels = metric_order))
  
  ggplot(fold_change_data, aes(x = Metric, y = log2_FC, fill = Patient)) +
    geom_bar(stat = "identity", fill = fill_color, width = 0.7) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "", 
      y = "log2 Fold Change (Day36 / Day1)",
      title = title_text
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 5)),
      plot.title.position = "plot",
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(margin = margin(t = 5), color = "black", size = 14),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      text = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )
}

# Create both plots
p_dc <- make_fc_plot("30207", "Increased immune\nresponses after therapy in DC", "#3182bd")
p_no_dc <- make_fc_plot("30209", "Decreased immune\nresponses after therapy in No_DC", "#9ecae1")

# Combine plots using patchwork
(p_dc | plot_spacer() | p_no_dc) +
  plot_layout(widths = c(1, 0.01, 1)) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 15))
    )
  )
