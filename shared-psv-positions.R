library("vcfR")
library("viridis")
library("ggplot2")
library("dplyr")
data.t <- read.vcfR("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/psvs_per_superpop/filtered_fval/AFR_EUR_psvs_v2.vcf.gz")

proximal_pos <- extract.info(data.t, "pos2")


split.names <- strsplit(proximal_pos, ":") #Select first column elements from this list (Sample_names)

vars_in_correct_format <- c();

vars_in_correct_format <- do.call(rbind, lapply(split.names, function(x) c(x[1], x[2])))
vars_in_correct_format <- as.data.frame(vars_in_correct_format, stringsAsFactors = FALSE)
colnames(vars_in_correct_format) <- c("CHROM", "POS")

distal_pos <- data.t@fix[, c("CHROM", "POS")]
distal_pos <- as.data.frame(distal_pos, stringsAsFactors = FALSE)

all_pos_pairs <- cbind(distal_pos, vars_in_correct_format)
dim(all_pos_pairs)

all_pos_pairs <- all_pos_pairs[,-3]
colnames(all_pos_pairs) <- c("CHROM", "DISTAL_POS", "PROXIMAL_POS")
write.table(all_pos, "~/Desktop/CMT1A/parascopy_calls/quality_psv_positions.txt", sep = "\t", quote = F, col.names = F, row.names = F)
dim(all_pos)

rownames(data.t) <- vars_in_correct_format

class(all_pos)
all_pos$POS <- as.numeric(all_pos$POS)
head(pos_in_common)

pos_in_common <- read.table("/Users/alexandra_al/Downloads/parascopy_AFR_EUR/isec_all_vars/sites.txt",sep = "\t",header = F)
colnames(pos_in_common) <- c("CHROM", "POS", "REF", "ALT", "T")
shared_pos <- inner_join(pos_in_common, all_pos, by = c("CHROM", "POS"))
dim(shared_pos)

#matched_pairs <- all_pos_pairs %>%
#  rowwise() %>%
#  filter(DISTAL_POS %in% pos_in_common$POS & PROXIMAL_POS %in% pos_in_common$POS) %>%
#  ungroup()

write.table(shared_pos[c(1,2)], "/Users/alexandra_al/Downloads/parascopy_AFR_EUR/isec_all_vars/shared_psv_sites.txt", sep = "\t", quote = F, col.names = F, row.names = F)


common_psvs <- read.vcfR("/Users/alexandra_al/Downloads/parascopy_AFR_EUR/isec_all_vars/common_psvs_AFR_EUR.vcf.gz")
counts <- gt.to.popsum(common_psvs)
