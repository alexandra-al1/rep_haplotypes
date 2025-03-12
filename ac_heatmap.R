library("vcfR")
library("viridis")
library("ggplot2")
#Compare aggregate (from BOTH CMT1A-REPs) alt allele count in PSVs of AFR vs EUR inds
data.t <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psvs_per_superpop/heat_plots/AFR_EUR_psvs.vcf.gz")
data.t <- extract.gt(data.t)

names <- rownames(data.t)

split.names <- strsplit(names, "_") #Select first column elements from this list (Sample_names)

vars_in_correct_format <- c();

for(i in 1:length(split.names)) { #For each of the IDs 
  vars_in_correct_format <- c(vars_in_correct_format, paste(split.names[[i]][2], sep = "_"))
  
}

rownames(data.t) <- vars_in_correct_format

#Remove rows with at least one NA (no GT) -> leaves 284 PSVs

data.t <- na.omit(data.t)
dim(data.t)

#gt_string <- "0/0/1/1"
sum_alt_counts <- function(gt_string) {
  
  split_string <- strsplit(gt_string, "/")
  gt_vector <- as.numeric(unlist(split_string))
  sum(gt_vector)
  
}

#Create new matrix of sum of alt allele counts
counts_matrix <- apply(data.t, c(1,2), sum_alt_counts)
dim(counts_matrix)

pops <- read.table(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/1kgp_files/metadata/Code_Continent_2504.txt",header=T, sep="\t")
dim(pops)

colnames(pops) <- c("Sample_name", "Population_code", "Superpopulation_code")
pops <- pops[pops$Superpopulation_code %in% c("AFR", "EUR"), ]
pops <- pops[pops$Population_code %in% c("CEU", "YRI"), ]
#Remove ID with dup?
pops <- pops[pops$Sample_name != "HG01896", ]
dim(pops)

pops_included <- pops[order(pops$Superpopulation_code),]

counts_matrix <- counts_matrix[, pops$Sample_name]
dim(counts_matrix)

pops_of_ids <- pops$Superpopulation_code
rownames(pops) <- pops$Sample_name
counts_matrix_t <- t(counts_matrix)
class(counts_matrix_t)
#Only for PSVs in hotspot zone 1
#14,184,471..14,187,746
counts_matrix <- counts_matrix[which(as.numeric(rownames(counts_matrix)) > 14184471 & as.numeric(rownames(counts_matrix)) < 14187746), ]
dim(counts_matrix)

heatmap(counts_matrix, scale="column", col = viridis(5), labCol = pops_of_ids)
min(counts_matrix)


library("pheatmap")
pheatmap(counts_matrix_t, color = viridis(5), cluster_rows = T, cluster_cols = T, clustering_distance_rows = "euclidean", clustering_method = "centroid", cutree_rows = 2, labels_row = pops_of_ids)

row_annotations <- data.frame(Pop = pops$Population_code)
rownames(row_annotations) <- rownames(pops)

row_annotations  <- row_annotations[order(row_annotations$Superpop),]

unique_superpops <- unique(row_annotations$Pop)

library("ggsci")
?ggsci
blues_palette <- colorRampPalette(c("lightblue", "blue"))
eur_colors <- blues_palette(1)
reds_palette <- colorRampPalette(c("#D62341", "darkred"))
afr_colors <- reds_palette(1)
print(colors)
annotation_colors <- list(Pop = setNames(c(afr_colors, eur_colors), unique_superpops))

pheatmap(
  mat = counts_matrix,
  scale = "none",
  color = viridis(5), # Use viridis for a 5-color palette
  labels_row = rep("", nrow(counts_matrix)),  # Remove column labels
  labels_col = rep("", ncol(counts_matrix)),  # Remove row labels
  cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete",
  annotation_col = row_annotations, # Column annotations
  annotation_colors = annotation_colors, # Annotation colors
  annotation_names_col = TRUE # Keep annotation names visible
)

pheatmap(
  mat = counts_matrix_t,
  scale = "none",
  color = viridis(5, direction = -1), # Use viridis for a 5-color palette
  labels_row = rep("", nrow(counts_matrix_t)),  # Remove column labels
  labels_col = rep("", ncol(counts_matrix_t)),  # Remove row labels
  cluster_rows = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
  cutree_rows = 2,
  annotation_row = row_annotations, # Column annotations
  annotation_colors = annotation_colors, # Annotation colors
  annotation_names_row = TRUE # Keep annotation names visible
)
