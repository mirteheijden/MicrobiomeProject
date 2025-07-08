########## Mice study, taxonomic bar plots
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
tse_GG2 <- importQIIME2("path_to/RESULTS/counts.qza", taxonomy = "path_to/RESULTS/taxonomy.qza", sampleMetaFile="path_to/config/metadata_GG2.tsv")
# change MI_ID rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Mirte_ID
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))


# genus - GG2
tse <- tse_GG2

# Choose samples with sample_points Day_7
tse <- tse[, which(colData(tse)$Sample_time %in% c("Day_7"))]
tse <- tse[, which(colData(tse)$Group %in% c("MOCK", "Virus", "A_shahii", "A_shahii+Virus"))]

tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

assay(tse_rel, "counts") <- sweep(assay(tse, "counts"), 2, colSums(assay(tse, "counts")), "/")

summary(assay(tse_rel, "counts"))


#MOCK
MOCK <- tse_rel[, rownames(subset(colData(tse_rel), Group == "MOCK"))]

tse_MOCK_genus <- mergeFeaturesByRank(MOCK, rank ="genus", onRankOnly=TRUE)
top_MOCK <- getTopFeatures(tse_MOCK_genus, top = 15)

#Virus
Virus <- tse_rel[, rownames(subset(colData(tse_rel), Group == "Virus"))]

tse_Virus_genus <- mergeFeaturesByRank(Virus, rank ="genus", onRankOnly=TRUE)
top_Virus <- getTopFeatures(tse_Virus_genus, top = 15)

#A_shahii
A_shahii <- tse_rel[, rownames(subset(colData(tse_rel), Group == "A_shahii"))]

tse_A_shahii_genus <- mergeFeaturesByRank(A_shahii, rank ="genus", onRankOnly=TRUE)
top_A_shahii <- getTopFeatures(tse_A_shahii_genus, top = 15)

#A_shahii_and_virus_and_virus
A_shahii_and_virus <- tse_rel[, rownames(subset(colData(tse_rel), Group == "A_shahii+virus"))]

tse_A_shahii_and_virus_genus <- mergeFeaturesByRank(A_shahii_and_virus, rank ="genus", onRankOnly=TRUE)
top_A_shahii_and_virus <- getTopFeatures(tse_A_shahii_and_virus_genus, top = 15)


# Renaming the "genus" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(tse_rel)$genus,
                        function(x){if (x %in% top_MOCK || x %in% top_Virus || x %in% top_A_shahii || x %in% top_A_shahii_and_virus ) {x} else {"Other"}})

rowData(tse_rel)$genus <- as.character(genus_renamed)

# Filter out rows where taxonomic group is "other"
tse_filtered <- tse_rel[rowData(tse_rel)$genus != "Other", ]

rowData(tse_filtered)$genus <- gsub("^g__", "", rowData(tse_filtered)$genus)
rownames(colData(tse_filtered)) <- colData(tse_filtered)$Group

# Get the current colData
coldata <- colData(tse_filtered)

# Update the rownames
rownames(coldata) <- gsub("^A_shahii$", "Alistipes", rownames(coldata))
rownames(coldata) <- gsub("^A_shahii\\+Virus$", "Virus + Alistipes", rownames(coldata))

# Save it back
colData(tse_filtered) <- coldata

## Different color scale:
# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

custom_colors <- c(
  "Akkermansia" = color_palette[1],       
  "Adlercreutzia_404257" = color_palette[2],
  "Acutalibacter" = color_palette[3],
  "Bacteroides_H" = color_palette[4],
  "CAG-873" = color_palette[5],
  "CAG-485" = color_palette[6],
  "COE1" = color_palette[7],
  "Duncaniella" = color_palette[8],
  "Dubosiella" = color_palette[9],
  #"Dubosiella" = color_palette[10],
  "Faecalibaculum" = color_palette[11],
  "Limosilactobacillus" = color_palette[12],
  "Mucispirillum" = color_palette[13],
  "Muribaculum" = color_palette[14],
  "NM07-P-09" = color_palette[15],
  "Paramuribaculum" = color_palette[16],
  "Parabacteroides_B_862066" = color_palette[17],
  "QWKK01" = color_palette[18],
  "UBA7173" = color_palette[19],
  "14-2" = color_palette[20],
  "Schaedlerella" = color_palette[21]
  #"Escherichia_710834" = color_palette[22],
  #"Lawsonibacter" = color_palette[23],
  #"Gordonibacter" = color_palette[24],
  #"Ruthenibacterium" = color_palette[25],
  #"" = color_palette[26],
  #"Ruminococcus_B" = color_palette[27],
  #"Streptococcus" = color_palette[28],
  #"taxa29" = color_palette[29],
  #"taxa30" = color_palette[30]
)


plotAbundance(tse_filtered, rank = "genus", order_rank_by = "abund", decreasing = TRUE, add_x_text = TRUE) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(angle = 90, size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14, face = "plain"),
    text = element_text(size = 12),
    panel.border = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +  # Set y-axis from 0 to 100%
  labs(title = "Day 7", x = "", y = "Rel. Abundance (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom color palette
  guides(fill = guide_legend(title = "Genus"),  # Rename the legend for fill to "species"
         color = "none")  # Remove the legend for the outline color (species)





# species - GG2
tse <- tse_GG2

# Choose samples with sample_points Day_7
tse <- tse[, which(colData(tse)$Sample_time %in% c("Day_7"))]
tse <- tse[, which(colData(tse)$Group %in% c("MOCK", "Virus", "A_shahii", "A_shahii+Virus"))]

tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")

assay(tse_rel, "counts") <- sweep(assay(tse, "counts"), 2, colSums(assay(tse, "counts")), "/")

summary(assay(tse_rel, "counts"))


#MOCK
MOCK <- tse_rel[, rownames(subset(colData(tse_rel), Group == "MOCK"))]

tse_MOCK_species <- mergeFeaturesByRank(MOCK, rank ="species", onRankOnly=TRUE)
top_MOCK <- getTopFeatures(tse_MOCK_species, top = 12)

#Virus
Virus <- tse_rel[, rownames(subset(colData(tse_rel), Group == "Virus"))]

tse_Virus_species <- mergeFeaturesByRank(Virus, rank ="species", onRankOnly=TRUE)
top_Virus <- getTopFeatures(tse_Virus_species, top = 12)

#A_shahii
A_shahii <- tse_rel[, rownames(subset(colData(tse_rel), Group == "A_shahii"))]

tse_A_shahii_species <- mergeFeaturesByRank(A_shahii, rank ="species", onRankOnly=TRUE)
top_A_shahii <- getTopFeatures(tse_A_shahii_species, top = 12)

#A_shahii_and_virus_and_virus
A_shahii_and_virus <- tse_rel[, rownames(subset(colData(tse_rel), Group == "A_shahii+virus"))]

tse_A_shahii_and_virus_species <- mergeFeaturesByRank(A_shahii_and_virus, rank ="species", onRankOnly=TRUE)
top_A_shahii_and_virus <- getTopFeatures(tse_A_shahii_and_virus_species, top = 12)


# Renaming the "species" rank to keep only top taxa and the rest to "Other"
species_renamed <- lapply(rowData(tse_rel)$species,
                          function(x){if (x %in% top_MOCK || x %in% top_Virus || x %in% top_A_shahii || x %in% top_A_shahii_and_virus ) {x} else {"Other"}})

rowData(tse_rel)$species <- as.character(species_renamed)

# Filter out rows where taxonomic group is "other"
tse_filtered <- tse_rel[rowData(tse_rel)$species != "Other", ]

rowData(tse_filtered)$species <- gsub("^s__", "", rowData(tse_filtered)$species)
rownames(colData(tse_filtered)) <- colData(tse_filtered)$Group

# Get the current colData
coldata <- colData(tse_filtered)

# Update the rownames
rownames(coldata) <- gsub("^A_shahii$", "Alistipes", rownames(coldata))
rownames(coldata) <- gsub("^A_shahii\\+Virus$", "Virus + Alistipes", rownames(coldata))

# Save it back
colData(tse_filtered) <- coldata

## Different color scale:
# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

## Different color scale:
# Generate a color palette with at least 30 colors from Set3
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(30)

custom_colors <- c(
  "Akkermansia muciniphila_D_776786" = color_palette[1],       
  "Adlercreutzia caecimuris" = color_palette[2],
  "Acutalibacter sp009936055" = color_palette[3],
  "Bacteroides_H xylanisolvens" = color_palette[4],
  "CAG-873 sp011959565" = color_palette[5],
  "CAG-CAG-485 sp002362485" = color_palette[6],
  "CAG-485 sp002493045" = color_palette[7],
  "Duncaniella dubosii" = color_palette[8],
  "Dubosiella newyorkensis" = color_palette[9],
  #"Dubosiella" = color_palette[10],
  "Faecalibaculum rodentium" = color_palette[11],
  "Limosilactobacillus reuteri" = color_palette[12],
  "Mucispirillum schaedleri" = color_palette[13],
  "Muribaculum gordoncarteri" = color_palette[14],
  "NM07-P-09 sp004793665" = color_palette[15],
  "Paramuribaculum intestinale" = color_palette[16],
  "Parabacteroides_B_862066 goldsteinii" = color_palette[17],
  "UBA7173 sp002491305" = color_palette[19],
  "14-2 sp000403315" = color_palette[20],
  "Schaedlerella sp000403295" = color_palette[21],
  "CAG-485 sp002493045" = color_palette[22],
  "UBA7173 sp002491305" = color_palette[23]
  #"Gordonibacter" = color_palette[24],
  #"Ruthenibacterium" = color_palette[25],
  #"" = color_palette[26],
  #"Ruminococcus_B" = color_palette[27],
  #"Streptococcus" = color_palette[28],
  #"taxa29" = color_palette[29],
  #"taxa30" = color_palette[30]
)


plotAbundance(tse_filtered, rank = "species", order_rank_by = "abund", decreasing = TRUE, add_x_text = TRUE) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(angle = 90, size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14, face = "plain"),
    text = element_text(size = 12),
    panel.border = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +  # Set y-axis from 0 to 100%
  labs(title = "Day 7", x = "", y = "Rel. Abundance (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom color palette
  guides(fill = guide_legend(title = "Species"),  # Rename the legend for fill to "species"
         color = "none")  # Remove the legend for the outline color (species)
