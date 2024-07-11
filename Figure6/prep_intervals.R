options(scipen = 999)
library(ggplot2)
library(dplyr)
library(reshape2)
library(purrr)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "hapname", action = "store")
)

opt <- parse_args(OptionParser(option_list=option_list, add_help_option = TRUE), positional_arguments = TRUE)

chr_labels <- c(
  chr1 = "1",
  chr2 = "2",
  chr3 = "3",
  chr4 = "4",
  chr5 = "5",
  chr6 = "6",
  chr7 = "7",
  chr8 = "8",
  chr9 = "9",
  chr10 = "10",
  chr11 = "11",
  chr12 = "12",
  chr13 = "13",
  chr14 = "14",
  chr15 = "15",
  chr16 = "16",
  chr17 = "17",
  chr18 = "18"
)

haps_chrs <- read.delim("Figure2_chromosome_lengths.txt", sep = "\t", header = TRUE)
colnames(haps_chrs) <- c("HAP","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")
chrlen_melt <-melt(haps_chrs, id.vars = "HAP", variable.name = "CHROM", value.name = "len")
chrlen_melt$CHROM <- factor(chrlen_melt$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
chrlen_melt <- subset(chrlen_melt, HAP %in% c(opt$options$input))
## chr-based plots


lookup <- read.delim("all_varpos_index.txt", sep = "\t", header = FALSE)
colnames(lookup) <- c("HAP","CHROM","POS","ID")
lookup$CHROM <- factor(lookup$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

haplotypes <- c(opt$options$input)

lapply(haplotypes, function(j) {
  print(paste("Doing ref ",j," now...", sep = ""), quote = FALSE)
  cactus_dat <- read.delim(paste("/path/",j,"_setup.txt",sep =""), sep = "\t", header = TRUE)
  cactus_dat[[j]] <- rep("0", length(cactus_dat$ID))
  
  melted <- melt(cactus_dat, id.vars = "ID", variable.name = "HAP", value.name = "GT")
  
  combined <- left_join(melted, lookup, by = c("ID","HAP"))
  final <- na.omit(combined)
  final$POS <- as.numeric(final$POS)
  
  final2 <- left_join(final, chrlen_melt, by = c("HAP","CHROM"))
  final3 <- final2[final2$GT != "0",]
  final3 <- final3[final3$GT != ".",]
  final3 <- data.frame(final3)
  rm(final,combined,melted,cactus_dat)
  haps <- unique(final3$HAP)
  chrs <- unique(final3$CHROM)
  combo  <- lapply(haps, function(i) {
    print(paste("Doing ",i," now...", sep = ""), quote = FALSE)
    subset <- final3[final3$HAP == i,]
    subset2 <- final2[final2$HAP == i,]
    interval_metadat <- lapply(chrs, function(z) {
      print(paste("Doing ",z," now...", sep = ""), quote = FALSE)
      chr_sub <- subset[subset$CHROM == z,]
      chr_sub2 <- subset2[subset2$CHROM == z,]
      
      plot <- ggplot(data = chr_sub, aes(x = POS)) +
        geom_histogram(binwidth = 100000) +
        ggtitle(paste(i,"_",z, sep = ""))
      
      plotdat <- ggplot_build(plot)
      histdat <- plotdat$data[[1]]
      
      breaks <- as.numeric(c(histdat$x))
      intervals <- t(as.data.frame(map2(breaks[-length(breaks)], breaks[-1] -1, c)))
      colnames(intervals) <- c("start","end")
      counts <- apply(intervals, 1, function(x) nrow(subset(chr_sub, POS >= x[1] & POS <= x[2])))
      varcounts <- apply(intervals, 1, function(x) nrow(subset(chr_sub2, POS >= x[1] & POS <= x[2])))
      dat <- data.frame(intervals, rep(z, nrow(intervals)), rep(i, nrow(intervals)))
      cbind(dat, counts, varcounts)
    })
    
    do.call(rbind, interval_metadat)
  })
  
  final_combo <- do.call(rbind, combo)
  rownames(final_combo) <- seq(1:length(final_combo$start))
  colnames(final_combo) <- c("start","end","chr","hap","count","varcounts")
  
  bedtools_test <- final_combo[,c(3,1,2,4,5,6)]
  
  bedtools_test <- bedtools_test[bedtools_test$varcounts > 0,]
  bedtools_test <- bedtools_test %>%
    mutate(strand = if_else(count <= 50, "+", "-"))
  
  bedtools_test <- bedtools_test[,c(1,2,3,4,5,7)]
  bedtools_test$start <- as.integer(bedtools_test$start)
  bedtools_test$end <- as.integer(bedtools_test$end)
  bedtools_test <- na.omit(bedtools_test)
  
  print(paste("Writing to file...", sep = ""), quote = FALSE)
  write.table(bedtools_test, paste("/outpath/",j,"_intervals.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
})
