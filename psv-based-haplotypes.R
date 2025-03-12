library("haplotypes")
library('vcfR')
library("ape")
rm(list=ls());
psv_vcf <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psv_analysis/quality_psvs_EUR_final.vcf.gz")

test <- as.data.frame(psv_vcf@fix)
test$POS <- as.numeric(test$POS)
nrow(test)
 #Split vcf into 2, distal and proximal
distal_vcf <- psv_vcf[as.numeric(psv_vcf@fix[, "POS"]) <= 15000000, ]
nrow(distal_vcf@fix) #155

proximal_vcf <- psv_vcf[as.numeric(psv_vcf@fix[, "POS"]) >= 15000000, ]
nrow(proximal_vcf@fix) #155

distal_dna <- vcfR2DNAbin(distal_vcf, extract.indels = F, verbose = T, extract.haps = F, unphased_as_NA = F)

proximal_dna <- vcfR2DNAbin(proximal_vcf, extract.indels = F, verbose = T, extract.haps = F, unphased_as_NA = F)
#as.character(distal_dna)

#write.dna(distal_dna, file = "/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psv_analysis/distal.fasta", format = "fasta", nbcol = -1, colsep = "")
#write.dna(proximal_dna, file = "/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psv_analysis/proximal.fasta", format = "fasta", nbcol = -1, colsep = "")

distal_seq <- as.dna(distal_dna)
class(distal_seq)
proximal_seq <- as.dna(proximal_dna)

#Distance calculations
distal_dist <- dist.gene(distal_dna, method = "percentage", pairwise.deletion = F) #distance d between two individuals is the number of loci for which they differ, and the associated variance is d(L − d)/L, where L is the number of loci
distal_dist <- as.matrix(distal_dist)
which(distal_dist == max(distal_dist), arr.ind = TRUE)

proximal_dist <- dist.gene(proximal_dna, method = "percentage", pairwise.deletion = F) #distance d between two individuals is the number of loci for which they differ, and the associated variance is d(L − d)/L, where L is the number of loci
proximal_dist <- as.matrix(proximal_dist)

pops <- read.table(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/1kgp_files/metadata/Code_Continent_2504.txt",header=T, sep="\t")
dim(pops)
colnames(pops) <- c("Sample_name", "Population_code", "Superpopulation_code")

#Remove ID with dup?
pops <- pops[pops$Sample_name != "HG01896", ]

afr_ids <- pops[pops$Superpopulation_code %in% "AFR", "Sample_name"]
eur_ids <- pops[pops$Superpopulation_code %in% "EUR", "Sample_name"]

hap_ids_afr <- c()
hap_ids_eur <- c()


for (i in 1:length(afr_ids)){
  hap_ids_afr <- c(hap_ids_afr, paste(afr_ids[i], "_0", sep=""), paste(afr_ids[i], "_1", sep=""))
}
for (i in 1:length(eur_ids)){
  hap_ids_eur <- c(hap_ids_eur, paste(eur_ids[i], "_0", sep=""), paste(eur_ids[i], "_1", sep=""))
}

distal_eur <- distal_dist[hap_ids_eur, hap_ids_eur]

distal_afr <- distal_dist[hap_ids_afr, hap_ids_afr]

dim(distal_afr)
#distal_afr <- distal_afr[upper.tri(distal_afr)]
#mean(distal_afr)

proximal_eur <- proximal_dist[hap_ids_eur, hap_ids_eur]

proximal_afr <- proximal_dist[hap_ids_afr, hap_ids_afr]

distal_seq <- as.dna(distal_dna)
proximal_seq <- as.dna(proximal_dna)

#Get haplotypes, indel = "5th" where gaps are a fifth character
distal_haps <- haplotypes::haplotype(distal_seq, indel = "5th")
proximal_haps <- haplotypes::haplotype(proximal_seq, indel = "5th")


nams <- labels(distal_seq)


sep <- strsplit(nams,"_");



continent <- c();



for(i in 1:length(sep)){
  
  sample_name <- c(continent,sep[[i]][1]);
  continent <- c(continent, pops$Superpopulation_code[pops$Sample_name %in% sample_name]) 
  
}

distal_haps_grouped <- grouping(distal_haps, factor = continent)


nams <- labels(proximal_seq)


sep <- strsplit(nams,"_");



continent <- c();



for(i in 1:length(sep)){
  
  sample_name <- c(continent,sep[[i]][1]);
  continent <- c(continent, pops$Superpopulation_code[pops$Sample_name %in% sample_name]) 
  
}

proximal_haps_grouped <- grouping(proximal_haps, factor = continent)

#Distal
haps_sum <- rowSums(distal_haps_grouped$hapmat[, 1:2])

distal_haps_grouped$hapmat <- cbind(distal_haps_grouped$hapmat, haps_sum)
D_group <- as.data.frame(distal_haps_grouped$hapmat)

D_group <- arrange(D_group, desc(haps_sum))




D_group <- D_group[(D_group$haps_sum >3), ]

sum(D_group$haps_sum) #2257 out of 2326 total, 97,0% inds
#2257/2326
dim(D_group)

#Proximal
haps_sum <- rowSums(proximal_haps_grouped$hapmat[, 1:2])

proximal_haps_grouped$hapmat <- cbind(proximal_haps_grouped$hapmat, haps_sum)
P_group <- as.data.frame(proximal_haps_grouped$hapmat)

P_group <- arrange(P_group, desc(haps_sum))
dim(P_group)

P_group <- P_group[(P_group$haps_sum >3), ]

sum(P_group$haps_sum) #2216 out of 2326 total, 95,27 % inds

dim(P_group)
#2216/2326
library(tidyr)

P_group <- P_group[,-ncol(P_group)]
D_group <- D_group[,-ncol(D_group)]
####Get haplotypes/cont plot #####
normalize_columns <- function(df) {
  # Calculate the sum of each column 
  col_sums <- colSums(df)
  
  # Divide each value by its column sum
  df <- sweep(df, 2, col_sums, FUN = "/")
  
  return(df)
}

normalized_df <- normalize_columns(D_group)
df <- as.data.frame(normalized_df)

# Add the row names as a column (assuming they are haplotype names)
df$haplotype <- rownames(D_group)
df$haplotype <- factor(df$haplotype, levels = unique(df$haplotype))


library("tidyr")
# Reorder columns to have the haplotype names as the first column
df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
long_df <- pivot_longer(df, cols = -haplotype, names_to = "Population", values_to = "Count")

library(scales)
library(ggsci)
library(pals)
library(ggplot2)
library("randomcoloR")

# Generate 47 distinct colors
palette <- distinctColorPalette(nrow(df))
palette = sample(color, nrow(df))
pal_ig
cool_v3 <- sample(pal_igv()(nrow(df)))
#cool_af <- cool
ggplot(long_df, aes(x = Population, y = -Count, fill = as.factor(haplotype))) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "", title = "") +
  scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0), labels = c("25%", "50%", "75%", "100%")) +
  theme_void() +
  scale_fill_manual(values = cool_v3) +  # Use custom color palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle=90, hjust = 1)) 


#Difference between distal and proximal PSV haps within inds
#Rename seqs to distinguish between distal and proximal in the joint matrix
nams_distal <- labels(distal_dna)

new_nams_distal <- c();

for (i in 1:length(nams_distal)){
  new_nams_distal <- c(new_nams_distal, paste("D", nams_distal[i], sep = "_"))
}

rownames(distal_dna) <- new_nams_distal

nams_proximal <- labels(proximal_dna)

new_nams_proximal <- c();

for (i in 1:length(nams_proximal)){
  new_nams_proximal <- c(new_nams_proximal, paste("P", nams_proximal[i], sep = "_"))
}

rownames(proximal_dna) <- new_nams_proximal

#Bind the two DNA matrices
combined_dna <- rbind(distal_dna, proximal_dna)
dim(combined_dna) # 4652 rows (2326 * 2)

distal_proximal_dist <- dist.gene(combined_dna, method = "percentage", pairwise.deletion = F)
dim(distal_proximal_dist)

#Make the matrix symmetric so that distal seqs are rows and proximal seqs are columns
distal_proximal_dist <- as.matrix(distal_proximal_dist)
distal_proximal_dist <- distal_proximal_dist[new_nams_distal,new_nams_proximal]
dim(distal_proximal_dist)
max(distal_proximal_dist)
#write.dna(combined_dna, file = "/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psv_analysis/combined_seqs.fasta", format = "fasta", nbcol = -1, colsep = "")

#Loop through each id, get cols and rows that only contain the id in the name, get mean(dist) from the 4 seqs/ID
length(eur_ids)
length(afr_ids)

avg_dist_eur_ids <- c();

for (i in 1:length(eur_ids)){
  sample_dist <- distal_proximal_dist[grep(paste(eur_ids[i], collapse = "|"), rownames(distal_proximal_dist), value = TRUE), grep(paste(eur_ids[i], collapse = "|"), colnames(distal_proximal_dist), value = TRUE)]
  avg_dist_sample <- mean(sample_dist)
  avg_dist_eur_ids <- c(avg_dist_eur_ids, avg_dist_sample)

  }

avg_dist_eur_ids <- setNames(avg_dist_eur_ids, eur_ids)

avg_dist_afr_ids <- c();

for (i in 1:length(afr_ids)){
  sample_dist <- distal_proximal_dist[grep(paste(afr_ids[i], collapse = "|"), rownames(distal_proximal_dist), value = TRUE), grep(paste(afr_ids[i], collapse = "|"), colnames(distal_proximal_dist), value = TRUE)]
  avg_dist_sample <- mean(sample_dist)
  avg_dist_afr_ids <- c(avg_dist_afr_ids, avg_dist_sample)
  
}
length(avg_dist_afr_ids)
avg_dist_afr_ids <- setNames(avg_dist_afr_ids, afr_ids)

combined_df <- as.data.frame(c(avg_dist_afr_ids,avg_dist_eur_ids))
combined_df$Sample_name <- rownames(combined_df)
colnames(combined_df)[1] <- "dist"

combined_df <- left_join(combined_df, pops, by = "Sample_name")

ggplot(combined_df, aes(x=Superpopulation_code, y=log(dist))) + 
  geom_boxplot()

#######################################################################################
####Haplotype networks in AFR, EUR separately using distal, proximal as clustering factor
#taking distal_seq from beginning of script

distal_seq <- as.matrix(distal_seq)

distal_row_names <- c();

for (i in 1:length(rownames(distal_seq))){
  distal_row_names <- c(distal_row_names, paste("D_", rownames(distal_seq)[i], sep = ""))
}

rownames(distal_seq) <- distal_row_names
proximal_seq <- as.matrix(proximal_seq)

proximal_row_names <- c();

for (i in 1:length(rownames(proximal_seq))){
  proximal_row_names <- c(proximal_row_names, paste("P_", rownames(proximal_seq)[i], sep = ""))
}
rownames(proximal_seq) <- proximal_row_names
dim(distal_seq) == dim(proximal_seq) #Check for equal dims, equal cols (var pos) AND equal nrow (subsetting by pop so same n of haps)

combined_distal_proximal <- data.frame(matrix(NA, nrow = nrow(distal_seq) * 2, ncol = ncol(distal_seq)))

for (i in 1:ncol(distal_seq)) {
  combined_distal_proximal[, i] <- c(distal_seq[, i], proximal_seq[, i])
}
rownames(combined_distal_proximal) <- c(rownames(distal_seq), rownames(proximal_seq))

combined_distal_proximal <- as.dna(combined_distal_proximal)

class(combined_distal_proximal)

afr_haplotypes <- haplotypes::haplotype(combined_distal_proximal)
#Repeat above process for EUR vcf, create combined DNA object and get haplotypes
eur_haplotypes <- haplotypes::haplotype(combined_distal_proximal_EUR)

#nams from the Dna object we used for haplotypes
nams <- labels(combined_distal_proximal_EUR)
sep <- strsplit(nams,"_");
rep_category <- c();

for(i in 1:length(sep)){
  
  category <- sep[[i]][1]
  rep_category <- c(rep_category, category) 
  
}

afr_haps_grouped <- grouping(afr_haplotypes, factor = rep_category)
#Repeat rep_category for EUR Dna object and get grouping
eur_haps_grouped <- grouping(eur_haplotypes, factor = rep_category)
class(afr_haplotypes)

#To check hap distribution in distal, proximal
eur_haps_grouped$hapmat
afr_haps_grouped$hapmat

#Distance matrix of haplotypes but row,colnames are first id that has the haplotype
afr_haplotypes@d

id_from_each_hap <- sapply(names(eur_haplotypes@haplist), function(hap_name) {
    hap_members <- as.character(eur_haplotypes@haplist[[hap_name]])
    if (length(hap_members) > 0) {
      sample(hap_members, 1) # Sample one ID randomly
    } else {
      NA # Handle empty haplotypes
    }
  })

length(id_from_each_hap)

hap_numbers <- as.numeric(gsub("haplotype", "", names(id_from_each_hap)))

names(id_from_each_hap) <- hap_numbers

#id_from_each_hap_distal <- id_from_each_hap[1:68] # in AFR
#id_from_each_hap_distal <- id_from_each_hap[1:40] #[1:40] in EUR

id_from_each_hap_proximal <- id_from_each_hap[41:length(id_from_each_hap)]

proximal_cols <- grep("^P_", colnames(eur_haplotypes@d), value = TRUE)

distal_rows <- grep("^D_", rownames(eur_haplotypes@d), value = TRUE)

dist_matrix_eur <- eur_haplotypes@d[distal_rows, proximal_cols, drop = FALSE] 
mean(dist_matrix_eur) #178.0616
mean(dist_matrix_afr) #176.851

library(ggplot2)
library(reshape2)

long_df <- melt(dist_matrix_afr)
long_df$cont <- "AFR"

dim(long_df)
sampled_long_df <- long_df[sample(nrow(long_df), nrow(long_dfa)), ]


combined <- rbind(long_dfa, sampled_long_df)

boxplot_plot <- ggplot(combined, aes(x = cont, y = value)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels


##########################################################################################
##Using consensus sequences obtained with vcfx
#For hotspot zone 1 

library("MSA2dist")
library(Biostrings)
library("msa")



distal_eur <- readDNAStringSet("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/fasta/distal_eur_seqs.fa", format = "fasta") 
distal_eur <- msaClustalW(distal_eur)
distal_eur_matrix <- as.matrix(distal_eur$distSTRING)


all_seqs <- read.FASTA("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/fasta/full_distal_proximal_eur_afr_seqs_aligned.fa", type = "DNA")

all_dist <- dist.dna(all_seqs ,model = "TN93", pairwise.deletion = F, as.matrix = T)

range(all_dist)

avg_dist_eur_ids <- c();

for (i in 1:length(eur_ids)){
  sample_dist <- all_dist[grep(paste(eur_ids[i], collapse = "|"), rownames(all_dist), value = TRUE), grep(paste(eur_ids[i], collapse = "|"), colnames(all_dist), value = TRUE)]
  avg_dist_sample <- mean(sample_dist)
  avg_dist_eur_ids <- c(avg_dist_eur_ids, avg_dist_sample)
  
}

min(avg_dist_eur_ids)


avg_dist_afr_ids <- c();

for (i in 1:length(afr_ids)){
  sample_dist <- all_dist[grep(paste(afr_ids[i], collapse = "|"), rownames(all_dist), value = TRUE), grep(paste(afr_ids[i], collapse = "|"), colnames(all_dist), value = TRUE)]
  avg_dist_sample <- mean(sample_dist)
  avg_dist_afr_ids <- c(avg_dist_afr_ids, avg_dist_sample)
  
}

min(avg_dist_afr_ids)

avg_dist_eur_ids <- setNames(avg_dist_eur_ids, eur_ids)

avg_dist_afr_ids <- setNames(avg_dist_afr_ids, afr_ids)

library("dplyr")
library("ggplot2")
combined_df <- as.data.frame(c(avg_dist_afr_ids,avg_dist_eur_ids))
combined_df$Sample_name <- rownames(combined_df)
colnames(combined_df)[1] <- "dist"

combined_df <- left_join(combined_df, pops, by = "Sample_name")

ggplot(combined_df, aes(x=Superpopulation_code, y=log(dist))) + 
  geom_boxplot()


#Empirical distribution test for EUR mean vs random AFR samples 
n_repetitions <- 1000;



ssim <- rep(NA,n_repetitions);#Simulation stats vector


set.seed(123)
count = 1;


while(count <= n_repetitions){
  
  sample_ids <- sample(nrow(distal_afr), 1006)
  sampled_rows <- distal_afr[sample_ids, sample_ids]
  #A symmetric matrix
  #sampled_gwd <- sampled_rows[,rownames(sampled_rows)]
  #Calculate mean for sampled GWD
  sampled_dist_tri <- sampled_rows[upper.tri(sampled_rows)]
  sampled_mean_dist <- mean(sampled_dist_tri)
  
  ssim[count] <- sampled_mean_dist;
  
  
  
  count = count + 1;
  
}

class(ssim)
max(ssim)
min(distal_afr)


pvalue <- mean(ssim <= mean(distal_eur))
ssim_df <- as.data.frame(ssim)

