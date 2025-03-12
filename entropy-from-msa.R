library(seqinr)
library(vegan)

library(ggplot2)

positionInfo <- function(data, verbose=TRUE){
  info <- list()
  if(verbose){cat("Analyzing positions \n")}
  for (i in 1:(nchar(as.character(datadf$sequence[1])))){
    ss <- substr(datadf[,2],1,1)
    ss <- toupper(ss)
    df <- data.frame(base=unique(ss))
    tmp <- c()
    for ( j in 1:length(unique(ss))){
      tmp[1] <- length(grep(df[1,1], ss))
    }
    df$rep <- tmp
    rownames(df) <- as.character(df[,1])
    df[,1] <- NULL
    
    dff <- subset(df, rownames(df)!="-")
    
    info$shan[1] <- vegan::diversity(t(df))                                    # Shannon entropy of position
    info$shanc[i] <- ifelse(nrow(dff)==0, NA, vegan::diversity(t(dff)))        # Shannon entropy of position removing gaps ("-")
    info$rich[i] <- length(unique(ss))                                  # Position richness
    info$richc[i] <- ifelse(nrow(dff)==0, NA, length(rownames(dff)))    # Position richness removing gaps ("-")
    info$uniq[i] <- list(unique(ss))                                    # Unique bases in position
    info$repe[i] <- list(tmp)                                           # Repetitions of the unique bases in position
    
    if(verbose){
      if((i %% (round(nchar(as.character(data$sequence[1]))*0.1,0))) == 0){message("  ", round(i/nchar(as.character(data$sequence[1]))*100, 0), "%")}
    }
  }
  
  if(verbose){cat("Creating the data frame with the information \n")}
  
  df <- data.frame(posi=c(1:length(info$shan)),
                   shan=info$shan,
                   shanc=info$shanc,
                   rich=info$rich,
                   richc=info$richc,
                   uniq=unlist(lapply(info$uniq, function(x) paste(x, collapse="|"))),
                   repe=unlist(lapply(info$repe, function(x) paste(x, collapse="|"))))
  if(verbose){cat("Finished \n")}
  return(df)
  
}


datadf <-  seqinr::read.fasta("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/fasta/proximal_eur_seqs_zone1_aligned.fa", seqtype= "DNA", as.string = T)
datadf <- data.frame(name=paste(getAnnot(datadf)), sequence=paste0(datadf))
datadf$name <- gsub(">", "",datadf$name)

names_in_correct_format <- c();


#Rename sequences with D/P_ AFR/EUR_ idname
for (i in 1:length(datadf$name)){
  sample_name <- sub("^[^_]+_([^_]+).*", "\\1", datadf$name[i])
  superpop_code <- pops[pops$Sample_name %in% sample_name, "Superpopulation_code"]
  new_name <- paste(superpop_code,"_",datadf$name[i], sep ="")
  names_in_correct_format <- c(names_in_correct_format, new_name)
  
}

datadf$name <- names_in_correct_format

if(!all(duplicated(nchar(as.character(datadf$sequence)))[-1])){stop("Sequences in fasta file are not aligned.")}

#Analyze all seqs together 
df <- positionInfo(datadf)

# Analyzing the alignment by samples
# This understands you have in the same aligned fasta different groups of sequences with the following name: identifier_sequenceName
sample <- unique(gsub("(^[^_]+_[^_]+).*", "\\1", datadf$name))

# And now looping the function through the different subsetted files
df <- data.frame()
for(i in sample){
  cat("Reading sequences starting by '", i, "'. ", grep(paste0("^", i, "$"), sample), "/", length(sample), "\n", sep="")
  ss <- datadf[grepl(paste0("^", i, "_"), datadf$name),]
  dfs <- positionInfo(ss)
  dfs$sample <- i
  df <- rbind(df, dfs)
  cat("_____________________________", "\n")
}; rm(i, ss, dfs)


#________________________________________________________________________
# If you just want to get the entropy of the full alignment file instead
df <- positionInfo(datadf)
#________________________________________________________________________


# Export the table with the information for every position. Note there is a new column for the sample (identifier for the group of sequence)
write.table(df, paste0(gsub("\\.[^\\.]+$", "", file), "_positionInfo.tsv"), quote=FALSE, row.names=FALSE, sep="\t")


# Plotting the entropy with ggplot2
dive <- ggplot(df, aes(x=posi, y=shan, colour=sample))+
  geom_point(aes())+
  geom_smooth(level=.99)+
  #facet_wrap(~sample, scales="free", nrow=3)+
  ggtitle("Shannon's entropy in AFR proximal CMT1A-REP")+
  theme_bw() + labs(x = "Position", y = "Entropy")
dive

distal_eur <- dive

# Exporting the PDF of the plot
pdf(paste(gsub("\\.[^\\.]+$", "", file), "_alignment_entropy.pdf", sep=""), width=11.69, height=8.27, paper='special')
plot(dive)
dev.off()



####Import aligned zone 1
datadf <-  read.FASTA("/Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/fasta//Users/alexandra_al/Desktop/CMT1A/parascopy_calls/all_vars_analysis/fasta/distal_proximal_eur_afr_seqs.fa")
class(datadf)

all_vs_all <- dist.dna(datadf)
