rm(list=ls());
library(seqinr)
library(vegan)

library(ggplot2)
library(vcfR)
pops <- read.table(file="/Users/alexandra_al/Desktop/CMT1A/Datasets/1kgp_files/metadata/Code_Continent_2504.txt",header=T, sep="\t")
dim(pops)
colnames(pops) <- c("Sample_name", "Population_code", "Superpopulation_code")
pops <- pops[pops$Sample_name != "HG01896", ]


afr_vcf <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/variants_AFR_phased_biallelic_snps.vcf.gz")

afr_vcf <- extract.gt(afr_vcf, return.alleles = T)
afr_vcf <- as.data.frame(afr_vcf)


#distal_afr_vcf <- afr_vcf[as.numeric(afr_vcf@fix[, "POS"]) <= 15000000, ]

#proximal_afr_vcf <- afr_vcf[as.numeric(afr_vcf@fix[, "POS"]) >= 15000000, ]


eur_vcf <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/variants_EUR_phased_biallelic_snps.vcf.gz")


eur_vcf <- extract.gt(eur_vcf, return.alleles = T)
eur_vcf <- as.data.frame(eur_vcf)

psv_vcf <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psv_analysis/quality_psvs_AFR_EUR_final_biallelic_snps.vcf.gz")
psv_vcf <- extract.gt(psv_vcf, return.alleles = T)
psv_vcf <- as.data.frame(psv_vcf)

eur_psvs <- psv_vcf[, pops[pops$Superpopulation_code %in% "EUR", "Sample_name"]]

afr_psvs <- psv_vcf[, pops[pops$Superpopulation_code %in% "AFR", "Sample_name"]]

dim(eur_psvs)
dim(psv_vcf)

#distal_eur_vcf <- eur_vcf[as.numeric(eur_vcf@fix[, "POS"]) <= 15000000, ]

#proximal_eur_vcf <- eur_vcf[as.numeric(eur_vcf@fix[, "POS"]) >= 15000000, ]

get_split_columns <- function(df) {
  # Create an empty list to store the new columns
  new_df <- list()
  
  # Loop through each column and split it
  for (col in colnames(df)) {
    # Split the values by "|"
    haplotypes <- strsplit(df[[col]], "\\|")
    
    # Create new columns for each haplotype
    new_df[[paste0(col, "_1")]] <- sapply(haplotypes, function(x) x[1])
    new_df[[paste0(col, "_2")]] <- sapply(haplotypes, function(x) x[2])
  }
  
  # Combine all new columns into a new data frame
  new_df <- as.data.frame(new_df)
  original_positions <- c();
  
  for (i in 1:length(rownames(df))) {
    # Extract the position and variant number
    rowname_parts <- unlist(strsplit(rownames(df)[i], "_"))
    position <- as.numeric(rowname_parts[2])
    original_positions <- c(original_positions, position)
  }
  
  rownames(new_df) <- original_positions
  return(new_df)
}


#Get df from whichever vcfR GT df
split_alleles_df <- get_split_columns(eur_psvs)

#Counts per base for diversity()
split_alleles_df$A <- rowSums(split_alleles_df == "A")
split_alleles_df$C <- rowSums(split_alleles_df == "C")
split_alleles_df$G <- rowSums(split_alleles_df == "G")
split_alleles_df$T <- rowSums(split_alleles_df == "T")

base_counts_df <- split_alleles_df[,(ncol(split_alleles_df) - 3):ncol(split_alleles_df)]

#Calculate entropy per pos
entropy_pos <- data.frame(row.names = rownames(base_counts_df))

for (i in 1:nrow(base_counts_df)) {
  entropy_pos[i,1] <- diversity(base_counts_df[i,])
}

entropy_psvs_afr <- entropy_pos
entropy_psvs_eur <- entropy_pos

dim(entropy_pos)

entropy_psvs_afr$Pos <- rownames(entropy_psvs_afr)
entropy_psvs_eur$Pos <- rownames(entropy_psvs_eur)

combined <- full_join(entropy_psvs_afr, entropy_psvs_eur, by = "Pos") %>%
  mutate(across(everything(), ~replace_na(., 0)))

colnames(combined) <- c("AFR", "Pos", "EUR")

ggplot(combined, aes(x = AFR, y = EUR)) +
  geom_point(color = "black", size = 3) +
  labs( x = "AFR", y = "EUR") +
  theme_minimal()


summary(lm(combined$AFR~combined$EUR))
summary(lm(combined$EUR~combined$AFR))
hist((combined$AFR-combined$EUR))
