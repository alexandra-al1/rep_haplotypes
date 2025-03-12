##########################################################################################

#

# Haplotype analysis

#

#


##########################################################################################

rm(list=ls());

library("pegas");

library("stringr")

library("colorRamps")
library("ape")
library("geneHapR")
library("adegenet")
library("haplotypes")
library(dplyr)
library("ggrepel")
library("ggplot2")

library("FactoMineR")

library("factoextra")
library("vcfR")



x <- read.fas(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/merged/common_3/proximal_zone_1-.fasta")
class(x)
x <- fasta2DNAbin(x)
hp <- pegas::haplotype(x, decreasing = T, what = "frequencies") #With maf > 0.05: distal has 7 haps, proximal has 32
diffHaplo(h,a = 1, b=2)
hap.div(x)
write.dna(h, "C:\\Users\\al-al\\Desktop\\Genetics and Genomics\\TFM\\Datasets\\1000G\\Haplotype network\\fasta\\maf_0.05\\proximal_haplotypes.fasta", format= "fasta")

write.nexus.data(h, file= "C:\\Users\\al-al\\Desktop\\Genetics and Genomics\\TFM\\Datasets\\1000G\\Haplotype network\\fasta\\maf_0.05\\distal_haplotypes.nex", format = "dna", datablock = F,
                 interleaved = TRUE, charsperline = NULL,
                 gap = NULL, missing = NULL)


####Grouping, pkg haplotypes

pops <- read.table(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/1kgp_files/metadata/Code_Continent_2504.txt",header=T, sep="\t")
dim(pops)
colnames(pops) <- c("Sample_name", "Population_code", "Superpopulation_code")

x <- read.vcfR(file="/Users/alexandra_al/Downloads/PLOD1/v2/isec_snps/samples_snps.vcf.gz")

x <- vcfR2DNAbin(x, extract.indels = F, verbose = T, extract.haps = T, unphased_as_NA = F)

samples_dist <- dist.gene(x, method = "percentage", pairwise.deletion = F) #distance d between two individuals is the number of loci for which they differ, and the associated variance is d(L − d)/L, where L is the number of loci
#dist.gene from ape: This function is meant to be very general and accepts different kinds of data (alleles, haplotypes, SNP, DNA sequences, ...).
#class(samples_dist)
samples_dist <- as.matrix(samples_dist)
dim(samples_dist) #18 x 18

#Remove one of the ids from the same family from rows and cols
samples_dist <- samples_dist[!grepl("^SEN_HFN_0030_0076", rownames(samples_dist)), ]
samples_dist <- samples_dist[, !grepl("^SEN_HFN_0030_0076", colnames(samples_dist))]

samples_dist <- samples_dist[upper.tri(samples_dist)]

samples_mean_dist <- mean(samples_dist)


y <- read.vcfR(file="/Users/alexandra_al/Downloads/PLOD1/v2/isec_snps/1kgp_snps_GWD.vcf.gz")

y <- vcfR2DNAbin(y, extract.indels = F, verbose = T, extract.haps = T, unphased_as_NA = F)


#haps_gwd <-  haplotypes::haplotype(y, indels = "sic")
gwd_dist <- dist.gene(y, method = "percentage", pairwise.deletion = F) #distance d between two individuals is the number of loci for which they differ, and the associated variance is d(L − d)/L, where L is the number of loci

gwd_dist <- as.matrix(gwd_dist)
dim(gwd_dist) # 226 x 226
gwd_mean_dist <- mean(gwd_dist)


#For gambians
n_repetitions <- 10000;



ssim <- rep(NA,n_repetitions);#Simulation stats vector


#set.seed(123)
count = 1;


while(count <= n_repetitions){
  
  sampled_ids <- sample(pops$Sample_name[pops$Population_code == "GWD"], 8) #8 because our unrelated patient samples are 8
  ids_string <- paste(sampled_ids, collapse = "|")
  #Extract both haplotypes from the sampled IDs
  extracted_haps <- rownames(gwd_dist)[grep(ids_string, rownames(gwd_dist))]
  sampled_gwd <- gwd_dist[extracted_haps,extracted_haps] 
  #Calculate mean for sampled GWD
  gwd_dist_tri <- sampled_gwd[upper.tri(sampled_gwd)]
  gwd_mean_dist <- mean(gwd_dist_tri)
    
  ssim[count] <- gwd_mean_dist;
    
    
    
    count = count + 1;
    
  }
  
class(ssim)

mean(ssim)

pvalue <- mean(ssim <= samples_mean_dist)
ssim_df <- as.data.frame(ssim)
dim(ssim_df)

test <- data.frame(
  pop = factor(c(rep("GWD", length(ssim)),
                   rep("Samples", length(samples_dist)))),
  dist = c(ssim, samples_dist)
)

ggplot(test, aes(x = pop, y = dist, fill = pop)) +
  geom_boxplot() +
  labs(title = "Genetic distances in GWD vs SENEGENE samples", x = "Group", y = "Distance") +
  theme_minimal()

density_plot <- ggplot(ssim_df, aes(x = ssim)) +
  geom_histogram(aes(y = after_stat(density)),
                 colour = "black", 
                 fill = "lightblue", 
                 binwidth = 0.01) +
  geom_density(aes(y = after_stat(density))) +
  theme_classic() +
  scale_y_continuous(expand = c(0.0, 0.02)) +
  geom_vline(aes(xintercept=samples_mean_dist), 
             colour="black", linetype="dashed", size=0.5) + 
#annotate("text", x = 0.222, y = 8,size = 3.5, label = "pvalue of higher distance between GWD samples\n compared to SENEGENE samples:\n 1e-04")+
  labs(y = "Density", x = "Distance between haplotype pairs in GWD individuals")



dpb <- ggplot_build(density_plot)

x1 <- min(which(dpb$data[[1]]$x >=.25))
x2 <- max(which(dpb$data[[1]]$x <= samples_mean_dist))

density_plot + geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                            y=dpb$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey", color= "black", size =0.4 , alpha= 0.7)  +
  geom_text(aes(x = 27, label = paste0("pvalue of higher distance between GWD samples\n compared to SENEGENE samples:\n 0.006"), y = 0.07),
            check_overlap = TRUE, size = 3)



which(h_p@d == 9)
details <- as.list(h_p)

write(details, "", sep = "\t")
lapply(details, write, "~/Desktop/prox_haps_check.txt", append=TRUE)
tab <- details$sequence

#Proximal
#HG00099_EUR_2 -> 4
#Distal
#HG00097_EUR_1 -> 3 1nt 6
#HG00099_EUR_1 -> 4 2nt 1,6
#HG00103_EUR_2 -> 6 4,6
d <- distance(x)

nams <- labels(y)


sep <- strsplit(nams,"_");



continent <- c();



for(i in 1:length(sep))
  
{
  
  continent <- c(continent,sep[[i]][2]);
  
}

sum(continent == "AFR")
continent <- as.vector(continent);
unique(continent)
library(dplyr)

g <- grouping(h_p, factors= continent)
nrow(g$hapmat)
haps_sum <- rowSums(g$hapmat[, 1:7])

g$hapmat <- cbind(g$hapmat, haps_sum)
P_p <- as.data.frame(g$hapmat)
P_p <- arrange(P_p, desc(haps_sum))

nrow(P_p)


P_p<- P_p[(P_p$haps_sum >4), ]
sum(P_p$haps_sum) 
library(tidyr)

P_d
####Get haplotypes/cont plot #####
normalize_columns <- function(df) {
  # Calculate the sum of each column (excluding the first column which contains haplotype names)
  col_sums <- colSums(df)
  
  # Divide each value by its column sum
  df <- sweep(df, 2, col_sums, FUN = "/")
  
  return(df)
}
P_d_df
# Apply the function to normalize the data frame
normalized_df <- normalize_columns(P_p_df)
df <- as.data.frame(normalized_df)

# Add the row names as a column (assuming they are haplotype names)
df$or_hap_proxi <- rownames(P_p_df)
df$or_hap_proxi <- factor(df$or_hap_proxi, levels = unique(df$or_hap_proxi))


library("tidyr")
# Reorder columns to have the haplotype names as the first column
df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
long_df <- pivot_longer(df, cols = -or_hap_proxi, names_to = "Population", values_to = "Count")

library(scales)
library(ggsci)
library(pals)
library(ggplot2)

cool_v3 <- sample(pal_igv()(50))
#cool_af <- cool
ggplot(long_df, aes(x = Population, y = -Count, fill = as.factor(or_hap_proxi))) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "", title = "") +
  scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0), labels = c("25%", "50%", "75%", "100%")) +
  theme_void() +
  scale_fill_manual(values = cool_af) +  # Use custom color palette
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle=90, hjust = 1)) 


#library(Polychrome)
#P50 = createPalette(47,c("#FFC20AFF","#003399FF", cool_v3[-(1:3)]), range = c(30, 60), target = "normal", M= 10000)
#P50 <- as.vector(P50)

P_d <- P_d[,-8]



nrow(P_d) #N filtered haps

head(g)

res.ca <- CA(P_d, graph = FALSE)
res.ca$eig

P_d <- P_d[!(row.names(P_d) %in% outlier_row), ]




get_ca(res.ca, "row")

haps <- get_ca_row(res.ca) #AFR always separated
haps_coords <- haps$coord
haps_2dim <- haps_coords[,1:2]
haps_sum <- rowSums(P_d[, 1:7])
afr_part <- P_d[,1]/sum(haps_sum)
haps_2dim <- cbind(haps_2dim, haps_sum, afr_part)
haps_2dim <- as.data.frame(haps_2dim)
P_d$or_hap <- rownames(P_d)

pie_df <- P_d[,c(9, 4)]
pie_df <- as.data.frame(pie_df)
#colourCount = length(unique(pie_df[,2]))
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#ggplot(pie_df, aes(x="", y=pie_df[,2], fill=pie_df[,1])) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +  theme_void() +
  scale_fill_manual(values=getPalette(colourCount))


  haps_2dim$original_hap <- rownames(haps_2dim) 
rownames(haps_2dim) <- seq(1:nrow(haps_2dim))

superpops <- get_ca_col(res.ca)
superpop_coords <- superpops$coord
superpop_2dim <- superpop_coords[,1:2]

new_min <- 1.5 # Minimum value in the new range
new_max <- 60  # Maximum value in the new range

# Define the current minimum and maximum values in the input data
current_min <- min(haps_2dim$haps_sum)
current_max <- max(haps_2dim$haps_sum)

library(colourvalues)
# Perform linear transformation
transformed_values <- (((haps_2dim$haps_sum - current_min) / (current_max - current_min)) * (new_max - new_min) + new_min)

test <- rank(haps_2dim[,4], ties.method = "first")
shades <- color_values(test, palette = "diverge_hcl")
show_colours("pubu")
#superpop_2dim <- superpop_2dim[-nrow(superpop_2dim),]
custom_ca <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_jitter(inherit.aes = FALSE, data = haps_2dim, aes(x = haps_2dim[, 1], y = haps_2dim[, 2]), shape = 16, size = 2, alpha = 0.8, color = "#832232", position = position_jitter(seed = 7)) +
  #geom_text_repel(data = haps_2dim, aes(x = haps_2dim[, 1], y = haps_2dim[, 2]), label= haps_2dim[,5], max.overlaps = 50) +
  geom_point(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2]), colour = "black", shape = "triangle") +
  geom_text_repel(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2], label = rownames(superpop_2dim)), direction = "both", color = "black") +
  labs(title = "Correspondence analysis of distal zone 1 \nCMT1A-REP haplotypes", x = paste("Dimension 1:", round(res.ca$eig[1,2]), "%"), y = paste("Dimension 2:", round(res.ca$eig[2,2]), "%")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))


library("colourvalues")






fviz_ca_biplot(res.ca, repel = TRUE, title ="Correspondence analysis for proximal CMT1A-REP haplotypes (MAF > 0.1)")



##############################
#####Previous version (does not work with indels, sequence length is not equal)
##### Alternative to get previous method to work: add Ns to alleles of indels
#################

#library("vcfR")
y <- read.vcfR(file="/Users/alexandra_al/Downloads/proximal_pos_patient_no_indel.vcf")
y <- vcfR2DNAbin(y, extract.indels = F, verbose = T)
library("ape")
x <- read.dna(file="/Users/alexandra_al/Desktop/CMT1A/all_vars/for_haps/proximal_haps_snps_zone1.fasta",format="fasta")

h <- pegas::haplotype(y)

h <- subset(h, minfreq= 1)

d <- dist.dna(h);



net <- mst(d);



nams <- labels(x)



sep <- strsplit(nams,"_");

# Pops contains 1kGP metadata

ids <- read.table("~/Downloads/pangenome_hifi.txt", sep = "\t")
ids_vec<- as.vector(ids$V1)
pops_3202 <- read.table("/Users/alexandra_al/Downloads/20130606_g1k_3202_samples_ped_population.txt", sep = " ", header = T)
pops_3202 <- pops_3202[, c(2,6,7)]


# Define a function to match sample names to superpop_names
match_superpop <- function(sample_name) {
  # Find the row where sample_name matches in other_data
  match_row <- pops$Sample_name == sample_name
  # Extract the corresponding superpop_name
  superpop <- pops$Superpopulation_code[match_row]
  # If no match found, return NA
  if (length(superpop) == 0) {
    return(NA)
  } else {
    return(superpop)
  }
}

# Initialize an empty vector to store superpop_names
continents <- character(length(ids_vec))

# Loop through each element in sep
for (i in 1:length(ids_vec)) {
  # Get the sample name (equivalent to sep[[i]][1])
  sample_name <- ids_vec[[i]][1]
  # Match sample name to superpop_name
  continent <- match_superpop(sample_name)
  # Store the superpop_name in the vector
  continents[i] <- continent
}


continent <- c();


for(i in 1:length(sep)){
  
  continent <- c(continent,sep[[i]][2]);
  
}



continents <- as.factor(continents);



categories <- unique(continents)



P <- pegas::haploFreq(x, fac = continents, haplo = h)


write.table(P, file = "~/Desktop/continent_haplotype.txt")

##http://127.0.0.1:26189/graphics/plot_zoom_png?width=1920&height=1017
net.labs <- attr(net, "labels");



rownames(P) <- net.labs

P <- P[net.labs, ]

sz <- summary(h);

sz <- sz[net.labs];

sz <- sz[!is.na(sz)];



heatcols <- rainbow(length(categories))

p <- pegas::haploNet(h)

plot(p, size = sz, show.mutation = 20, scale.ratio = 1, pie = P, bg=heatcols, cex = 0.1, main = "Minimum spanning tree in proximal", legend = FALSE, label =F);

replot()

legend(1200, 100, categories, col=heatcols, legend= categories, fill=heatcols, cex=0.5, title = "Continental groups") 
###


###CA from haploFreq pegas
res.ca <- CA(df, graph = FALSE)

get_ca(res.ca, "row")



haps <- get_ca_row(res.ca) #AFR always separated
haps_coords <- haps$coord
haps_2dim <- haps_coords[,1:2]
haps_sum <- rowSums(P_d[, 1:7])
afr_part <- haps_sum/sum(haps_sum)
haps_2dim <- cbind(haps_2dim, haps_sum, df, afr_part)

P_d

superpops <- get_ca_col(res.ca)
superpop_coords <- superpops$coord
superpop_2dim <- superpop_coords[,1:2]

new_min <- 2  # Minimum value in the new range
new_max <- 60  # Maximum value in the new range

# Define the current minimum and maximum values in the input data
current_min <- min(haps_sum)
current_max <- max(haps_sum)

# Perform linear transformation
transformed_values <- (((haps_sum - current_min) / (current_max - current_min)) * (new_max - new_min) + new_min)*0.8

# Print the transformed values
print(transformed_values)



custom_ca <- ggplot() + 
  geom_jitter(inherit.aes = F, data = haps_2dim, aes(x = haps_2dim[, 1], y = haps_2dim[, 2]), shape = "circle", size = transformed_values, alpha = 0.7,  height = 0.3, width = 0.25, color = shades) +
  geom_point(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2]), colour = "black", shape = "triangle") +
  geom_text(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2], label = rownames(superpop_2dim)), vjust = -0.8, hjust = 0.5) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkblue") +
  geom_hline(yintercept=0, linetype="dashed", color= "darkblue")

shades <- color_values(test, palette = "blues")
library("colourvalues")
install.packages("colourvalues")
shades <- color_values(1:11)
test <- rank(haps_2dim[,9], ties.method = "first")



####Function for plotting boxplots of genetic distance between distal and proximal CMT1A-REP for all individuals within a superpopulation
#"1" = "#7d2430", "AMR" = "#ce8147", "EAS" = "#3e7a4a", "EUR" = "#2e86ab","SAS" = "#dbd053", "MEA" = "#904082" , "OCN" = "#D65A90"
install.packages("ggpubr")
library(ggpubr)
kk <- read.table("/Users/alexandra_al/Downloads/Difference_between_distal_and_proximal_2024.txt", header=T, sep= " ")

kk_cont <- kk[kk$CONTINENT == "EUR", ]
summary(kk_cont$DISTANCE)[1]

conts <- unique(kk$CONTINENT)
summary_stats <- matrix(NA, nrow = length(conts), ncol = 6)
rownames(summary_stats) <- conts

for (i in 1:length(conts)) {
  cont <- conts[i]
  kk_cont <- kk[kk$CONTINENT == cont, ]
  summary_stats[i,] <- summary(kk_cont$DISTANCE)
  
}
colnames(summary_stats) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")

write.table(summary_stats, "/Users/alexandra_al/Desktop/CMT1A/Plots/dist_distal_proximal_superpops/summary_stats.txt", quote = F, sep = "\t", row.names = T, col.names = T)
ggplot(kk, aes(x=CONTINENT, y=DISTANCE, fill = CONTINENT)) + 
  geom_boxplot(alpha= 1) + theme_minimal() +
  #geom_jitter(shape=16, position=position_jitter(width = 0.4, height = 0.05), size = 0.4, color = "black") +
  labs(x= "Superpopulation", y="Distance between distal and proximal \nCMT1A-REP hotspot zone 1") +
  scale_fill_manual(values = c("AFR" = "#915058", "AMR" = "#ce8147", "EAS" = "#3e7a4a", "EUR" = "#2e86ab","SAS" = "#dbd053", "MEA" = "#904082" , "OCN" = "#CE749B" ))+
  guides(fill= "none") +
  stat_compare_means(method = "wilcox.test", ref.group = "AFR", aes(label = after_stat(p.signif)))


####################
###K-means clustering of haplotype sequences
####################

library("ape")
library("geneHapR")
library("adegenet")
library("haplotypes")
library(dplyr)
library("ggrepel")
library("ggplot2")

library("FactoMineR")

library("factoextra")

x <- read.vcfR(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/merged/common_3/distal_zone_1_maf_0.05.vcf.gz")
z <- vcfR2DNAbin(x, extract.indels = F, verbose = T) ##Class DNAbin
z <- as.dna(z) ## Class Dna
#ydf <- as.data.frame(y) rows are ids, those we will replace with the respective pop

merged_pops <- read.table("/Users/alexandra_al/Desktop/CMT1A/Datasets/merged_metadata/1kgp_hgdp_cab_moz.txt", sep = "\t", header = F)
colnames(merged_pops) <- c("Sample_name", "Population_code", "Superpopulation_code")

#In haplotype(), indels can be "sic", "5th", "missing")
h_d <- haplotypes::haplotype(z, indels = "sic") #Get distance matrix from haplotypes, not pairs of individuals sequences. 54 x 54 matrix
h_p@d #The distance matrix
h_p@sequence


dist_d <- h_d@d

# Assuming h_p@hapind is a named list
haplotype_lookup <- unlist(h_d@hapind)
haplotype_lookup <- names(haplotype_lookup)

split_names <- strsplit(haplotype_lookup, "\\.", fixed = FALSE)

split_df <- do.call(rbind, lapply(split_names, function(x) x[1:2]))

colnames(split_df) <- c("Haplotype", "ID")

# Convert to data frame
split_df <- as.data.frame(split_df, stringsAsFactors = FALSE)

common_ids <- split_df[split_df$ID %in% rownames(dist_p), ]

id_to_haplotype <- setNames(common_ids$Haplotype, common_ids$ID)

#Replace the rownames and colnames of dist_p using the mapping
rownames(dist_d) <- id_to_haplotype[rownames(dist_d)]
colnames(dist_d) <- id_to_haplotype[colnames(dist_d)]

#Get haps/cont by matching IDs to Superpop_name and using that to group ind per haplotype
nams <- labels(z)


sep <- strsplit(nams,"_");

# Define a function to match sample names to superpop_names
match_superpop <- function(sample_name) {
  # Find the row where sample_name matches in other_data
  match_row <- merged_pops$Sample_name == sample_name
  # Extract the corresponding superpop_name
  superpop <- merged_pops$Superpopulation_code[match_row]
  # If no match found, return NA
  if (length(superpop) == 0) {
    return(NA)
  } else {
    return(superpop)
  }
}

# Initialize an empty vector to store superpop_names
continents <- character(length(sep))

# Loop through each element in sep
for (i in 1:length(sep)) {
  # Get the sample name (equivalent to sep[[i]][1])
  sample_name <- sep[[i]][1]
  # Match sample name to superpop_name
  continent <- match_superpop(sample_name)
  # Store the superpop_name in the vector
  continents[i] <- continent
}

g <- grouping(h_d, factors= continents)

nrow(g$hapmat)
haps_sum <- rowSums(g$hapmat[, 1:7])

g$hapmat <- cbind(g$hapmat, haps_sum)
P_d <- as.data.frame(g$hapmat)
P_d <- arrange(P_d, desc(haps_sum))

P_d <- P_d[,-8]
colnames(P_d)
ncol(P_p)

P_p<- P_p[(P_p$haps_sum >4), ]
sum(P_p$haps_sum) 
library(tidyr)

res.ca <- CA(P_d, graph = FALSE)
res.ca$eig

get_ca(res.ca, "row")

haps <- get_ca_row(res.ca) #AFR always separated
haps_coords <- haps$coord
haps_2dim <- haps_coords[,1:2]
#haps_sum <- rowSums(P_p[, 1:7])
#afr_part <- P_d[,1]/sum(haps_sum)
#haps_2dim <- cbind(haps_2dim, haps_sum, afr_part)
haps_2dim <- as.data.frame(haps_2dim)
P_d$or_hap <- rownames(P_d)

pie_df <- P_d[,c(9, 4)]
pie_df <- as.data.frame(pie_df)

superpops <- get_ca_col(res.ca)
superpop_coords <- superpops$coord
superpop_2dim <- superpop_coords[,1:2]

custom_ca <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_jitter(inherit.aes = FALSE, data = haps_2dim, aes(x = haps_2dim[, 1], y = haps_2dim[, 2]), shape = 16, size = 2, alpha = 0.8, color = "#832232", position = position_jitter(seed = 7)) +
  #geom_text_repel(data = haps_2dim, aes(x = haps_2dim[, 1], y = haps_2dim[, 2]), label= haps_2dim[,5], max.overlaps = 50) +
  geom_point(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2]), colour = "black", shape = "triangle") +
  geom_text_repel(data = superpop_2dim, aes(x = superpop_2dim[, 1], y = superpop_2dim[, 2], label = rownames(superpop_2dim)), direction = "both", color = "black") +
  labs(title = "Correspondence analysis of distal zone 1 \nCMT1A-REP haplotypes", x = paste("Dimension 1:", round(res.ca$eig[1,2]), "%"), y = paste("Dimension 2:", round(res.ca$eig[2,2]), "%")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#For clustering
max(dist_d) #Normalize distance matrix before clustering, divide all cells in dist() by the total n of variants in the vcf
dist_d <- dist_d/6

kmeans(x = dist_p, centers = 4)

library(factoextra) # fviz_nbclust()
library(NbClust) # NbClust

### Create database with Z-scores from EFA ###

### Check adequacy ###
cortest.bartlett(z_scores_EFA)
cortest.bartlett(z_scores_EFA, n=5831, diag=TRUE, R = ) # Same result

### Optimal cluster number ###
fviz_nbclust(dist_d, kmeans, method = "wss") # Elbow method -> 3 clusters
fviz_nbclust(dist_d, kmeans, method = "silhouette") # Silhouette method -> 3 clusters
fviz_nbclust(dist_d, kmeans, method = "gap_stat") # Gap statistic method -> 3 clusters

clusters <- cutree(hc, k = 3) 

# Plot dendrogram
plot(hc)

## NbClust calculates 30 indices and recommends number of clusters
NbClust(diss, distance = "NULL", min.nc = 2, max.nc = 10, method = "kmeans") # Takes a while

# *** : The Hubert index is a graphical method of determining the number of clusters.
# In the plot of Hubert index, we seek a significant knee that corresponds to a 
# significant increase of the value of the measure i.e the significant peak in Hubert
# index second differences plot. 
# 
# *** : The D index is a graphical method of determining the number of clusters. 
# In the plot of D index, we seek a significant knee (the significant peak in Dindex
#                                                     second differences plot) that corresponds to a significant increase of the value of
# the measure. 
# 


### Hierarchical Clustering ###
hckmo2 <- hclust(dist(dist_d), "ward.D2") #  ward.D2 implements Ward's criterion,  the dissimilarities are squared before cluster updating
hckmo2
plot(hckmo2) # Obtain dendrogram

### Cut point to create the clusters ###
cluster <- cutree(hckmo2, h=40, k=3)

### Create data frame combining Z-scores and cluster ###
mixed2 <- as.data.frame(cbind(z_scores_EFA, cluster))
mixed2
summary(as.factor(mixed2$cluster))

### Differences between clusters in factors ### 
boxplot(PA1~cluster, data=mixed2) 
boxplot(PA2~cluster, data=mixed2) 


#########################################################



### k-means() #####

## K-means clustering, 3 or 5 ####
k_means <- kmeans(dist_p, 3, iter.max = 20, nstart = 25)
k_means

## Means of z-score of each factor in each cluster ##
aggregate(dist_p, by=list(cluster=k_means$cluster), mean)

## Visualize k-means clustering ##
fviz_cluster(k_means, dist_p, geom = "text")

rownames(P_p) <- paste("haplotype", rownames(P_p), sep = "")
cluster <- k_means$cluster 
combined <- cbind(P_p, cluster)
combined <- combined[combined$cluster == 3, ]
colSums(combined)

#####PLOD1 paper


####Grouping, pkg haplotypes
x <- read.vcfR(file="/Users/alexandra_al/Downloads/PLOD1/v2/isec_common/samples_common.vcf")

x <- vcfR2DNAbin(x, extract.indels = F, verbose = T, extract.haps = T, unphased_as_NA = F)

x <- as.dna(x)

haps_samples <- haplotypes::haplotype(x, indels = "sic")
dist_samples <- haps_samples@d

#Remove one of the ids from the same family from rows and cols
dist_samples <- dist_samples[!grepl("^SEN_HFN_0030_0075", rownames(dist_samples)), ]
dist_samples <- dist_samples[, !grepl("^SEN_HFN_0030_0075", colnames(dist_samples))]

y <- read.vcfR(file="/Users/alexandra_al/Downloads/PLOD1/v2/isec_common/1kgp_common.vcf")

y <- vcfR2DNAbin(y, extract.indels = F, verbose = T, extract.haps = T, unphased_as_NA = F)

#y <- as.dna(y)

#haps_1kgp <-  haplotypes::haplotype(y, indels = "sic")
#dist_1kgp <- haps_1kgp@d
#head(dist_1kgp)

h <- pegas::haplotype(y)

pops <- read.table(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/1kgp_files/metadata/Code_Continent_2504.txt",header=T, sep="\t")
dim(pops)


colnames(pops) <- c("Sample_name", "Population_code", "Superpopulation_code")



h <- subset(h, minfreq= 2)

d_1kgp <- dist.dna(y);



net <- mst(d_1kgp);



nams <- labels(y)



sep <- strsplit(nams,"_");

ids_vec <- c();


for(i in 1:length(sep)){
  
  ids_vec <- c(ids_vec,sep[[i]][1]);
  
}


# Define a function to match sample names to superpop_names
match_superpop <- function(sample_name) {
  # Find the row where sample_name matches in other_data
  match_row <- pops$Sample_name == sample_name
  # Extract the corresponding superpop_name
  superpop <- pops$Superpopulation_code[match_row]
  # If no match found, return NA
  if (length(superpop) == 0) {
    return(NA)
  } else {
    return(superpop)
  }
}

# Initialize an empty vector to store superpop_names
continents <- character(length(ids_vec))

length(ids_vec)

# Loop through each element in sep
for (i in 1:length(ids_vec)) {
  # Get the sample name (equivalent to sep[[i]][1])
  sample_name <- ids_vec[i]
  # Match sample name to superpop_name
  continent <- match_superpop(sample_name)
  # Store the superpop_name in the vector
  continents[i] <- continent
}

length(continents)

continents <- as.factor(continents);



categories <- unique(continents)



P <- pegas::haploFreq(y, fac = continents, haplo = h)


write.table(P, file = "~/Desktop/continent_haplotype.txt")

##http://127.0.0.1:26189/graphics/plot_zoom_png?width=1920&height=1017
net.labs <- attr(net, "labels");



plod1_net <- haploNet(h, getProb = T)

rownames(P) <- net.labs

P <- P[net.labs, ]

sz <- summary(h);

sz <- sz[net.labs];

sz <- sz[!is.na(sz)];


plot(plod1_net);

replot()

legend(1200, 100, categories, col=heatcols, legend= categories, fill=heatcols, cex=0.5, title = "Continental groups") 

