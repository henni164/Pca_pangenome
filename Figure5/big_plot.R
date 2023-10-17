options(scipen = 999)
library(ggplot2)
library(dplyr)
library(reshape2)
library(purrr)


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

haps_chrs <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/chr_test.txt", sep = "\t", header = TRUE)
colnames(haps_chrs) <- c("HAP","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")
chrlen_melt <-melt(haps_chrs, id.vars = "HAP", variable.name = "CHROM", value.name = "len")
chrlen_melt$CHROM <- factor(chrlen_melt$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
chrlen_melt <- subset(chrlen_melt, HAP %in% c("hap3","hap4","hap5","hap6","hap7","hap8","hap21","hap22","hap23","hap24","hap25","hap26","hap27",'hap28'))


## chr-based plots

lookup <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/Varpos_plotting/all_varpos_index.txt", sep = "\t", header = FALSE)
colnames(lookup) <- c("HAP","CHROM","POS","ID")
lookup$CHROM <- factor(lookup$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

haplotypes <- c("hap3","hap4","hap5","hap6","hap7","hap8","hap23","hap26","hap27","hap28")
#j <- "hap3"
thingy <- lapply(haplotypes, function(j) {
  cactus_dat <- read.delim(paste("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/Varpos_plotting/",j,"_setup.txt",sep =""), sep = "\t", header = TRUE)
  cactus_dat[[j]] <- rep("0", length(cactus_dat$ID))
  
  melted <- melt(cactus_dat, id.vars = "ID", variable.name = "HAP", value.name = "GT")
  
  combined <- left_join(melted, lookup, by = c("ID","HAP"))
  final <- na.omit(combined)
  final$POS <- as.numeric(final$POS)
  
  final2 <- left_join(final, chrlen_melt, by = c("HAP","CHROM"))
  final3 <- final2[final2$GT != "0",]
  final3 <- final3[final3$GT != ".",]
  final3 <- data.frame(final3)
  rm(final2,final,combined,melted,cactus_dat)
  haps <- unique(final3$HAP)
  chrs <- unique(final3$CHROM)
  #i <- "hap7"
  combo  <- lapply(haps, function(i) {
    
    subset <- final3[final3$HAP == i,]
    #z <- "chr1"
    interval_metadat <- lapply(chrs, function(z) {
      
      chr_sub <- subset[subset$CHROM == z,]
      
      plot <- ggplot(data = chr_sub, aes(x = POS)) +
        geom_histogram(binwidth = 10000) +
        ggtitle(paste(i,"_",z, sep = ""))
      
      plotdat <- ggplot_build(plot)
      histdat <- plotdat$data[[1]]
      
      breaks <- as.numeric(c(histdat$x))
      intervals <- t(as.data.frame(map2(breaks[-length(breaks)], breaks[-1] -1, c)))
      colnames(intervals) <- c("start","end")
      counts <- apply(intervals, 1, function(x) nrow(subset(chr_sub, POS >= x[1] & POS <= x[2])))
      dat <- data.frame(intervals, rep(z, nrow(intervals)), rep(i, nrow(intervals)))
      cbind(dat, counts)
    })
    
    do.call(rbind, interval_metadat)
  })
  
  final_combo <- do.call(rbind, combo)
  rownames(final_combo) <- seq(1:length(final_combo$start))
  colnames(final_combo) <- c("start","end","chr","hap","count")
  
  bedtools_test <- final_combo[,c(3,1,2,4,5)]
  bedtools_test <- bedtools_test %>%
    mutate(strand = if_else(count <= 5, "+", "-"))
  
  
  bedtools_test$start <- as.integer(bedtools_test$start)
  bedtools_test$end <- as.integer(bedtools_test$end)
  bedtools_test <- na.omit(bedtools_test)
  
  write.table(bedtools_test, paste("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/Varpos_plotting/",j,"_intervals.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  bedtools_test
})

#combined <- do.call(rbind, thingy)
#combined$length <- (combined$end - combined$start)

#ggplot(data = combined) + 
  #geom_point(aes(x = length, y = density)) +
  #coord_cartesian(ylim = c(0,0.004))

#hist(combined$density, breaks = 500, xlim = c(0,0.004))

#####################
#####################

merged_intervals <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/Varpos_plotting/all_intervals.txt", sep = "\t", header = FALSE)
colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
keepme <- final_joined
keepme <- keepme[keepme$HAP != "hap21" & keepme$HAP != "hap22" & keepme$HAP != "hap24",]

keepme$REFHAP <- factor(keepme$REFHAP, levels = c("hap6","hap23","hap28","hap27","hap8","hap7","hap26","hap5","hap4","hap3"))
keepme$HAP <- factor(keepme$HAP, levels = c("hap3","hap4","hap25","hap5","hap26","hap7","hap8","hap27","hap28","hap23","hap6"))
keepme$CHROM <- factor(keepme$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

plot <- ggplot(data = keepme) +
  geom_segment(aes(x = 1, xend = len, y = HAP, yend = HAP, color = HAP), linewidth = 3) +
  geom_segment(aes(x = start, xend = end, y = HAP, yend = HAP, color = REFHAP), linewidth = 3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free", space = "free", switch = "both", margins = FALSE) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_discrete(limits=rev, expand = c(0.001,0.001)) +
  scale_color_manual(breaks = c("hap3","hap4","hap26","hap7","hap8","hap5","hap27","hap28","hap23","hap6","hap25"), values = c("pink","purple","skyblue","orange","black","hotpink","#004BCF","#B0EE2C","#436200","#DCC60B","grey"), name = "Haplotype") +
  labs(x = "Position on chromosome (Mbp)", y = NULL) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.key = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "bottom",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_pangenome/Pangenome_analysis/Varpos_plotting/pangenome_blocks.tiff", plot, device = "tiff", width = 10, height = 7, units = "in", dpi = 600)



plot <- ggplot(data = final3) +
  geom_segment(aes(x = 1, xend = len, y = HAP, yend = HAP), color = "black", linewidth = 3, lineend = "round") +
  geom_segment(aes(x = POS-300, xend = POS+300, y = HAP, yend = HAP), color = "white", linewidth = 3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free", space = "free", switch = "y") +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(labels = function(x)x/1000000) +
  labs(x = "Position on chromosome (Mbp)", y = NULL) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8),
        strip.background.y = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.key = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        panel.spacing = unit(0.1, units = "lines"),
        legend.position = "none")