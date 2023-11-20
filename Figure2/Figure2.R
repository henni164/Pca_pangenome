options(scipen = 999)
library(ggplot2)
library(reshape2)
library(egg)

## Figure 2a chromosome lengths
chr_dat <- read.delim("Figure2_chromosome_lengths.txt", sep = "\t", header = TRUE)
colnames(chr_dat) <- c("Haplotype", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18")
chr_dat$Haplotype <- as.factor(chr_dat$Haplotype)
chr_dat$Haplotype <- factor(chr_dat$Haplotype, levels = c("hap1","hap2","hap3", "hap4", "hap5", "hap6", "hap7", "hap8", "hap9", "hap10", "hap11", "hap12", "hap13", "hap14", "hap15", "hap16", "hap17", "hap18", "hap19", "hap20","hap21","hap22","hap23","hap24","hap25","hap26","hap27","hap28","hap29","hap30","hap31","hap32"))

chr_melted <- melt(chr_dat, value.name = "Length")
colnames(chr_melted)[2] <- "Chr"
chr_melted$Length <- chr_melted$Length/1000000

chr_lengths_plot <- ggplot(chr_melted) +
  geom_boxplot(aes(x = Chr, y = Length), color = "grey60") +
  geom_point(aes(x = Chr, y = Length), color = "black", size = 0.5) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  labs(x = "Chromosome", y = "Length (Mbp)") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey90"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))

chr_lengths_plot

ggsave("Figure2a.tiff", chr_lengths_plot, device = "tiff", height = 4.5, width = 5, units = "in", dpi = 600)


## Figure 2b chr9 region dotplots

hap8vshap21 <- read.delim("hap25_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
hap8vshap21 <- hap8vshap21[,-(13:17)]
colnames(hap8vshap21) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap8vshap21 <- hap8vshap21[grep("chr", hap8vshap21$qname),]

hap6vshap23 <- read.delim("hap26_vs_hap27_chr9.paf", sep = "\t", header = FALSE)
hap6vshap23 <- hap6vshap23[,-(13:17)]
colnames(hap6vshap23) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")


hap5vshap24 <- read.delim("hap26_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
hap5vshap24 <- hap5vshap24[,-(13:17)]
colnames(hap5vshap24) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap5vshap24 <- hap5vshap24[grep("chr", hap5vshap24$qname),]

hap5vshap24_plot <- ggplot(data = hap5vshap24) +
  geom_segment(data = hap5vshap24, aes(x = qstart, xend = qend, y = tstart, yend = tend), linewidth = 0.1) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "%identity", labels = c("0","25","50","75","100")) +
  coord_cartesian(xlim = c(1600000,3000000), ylim = c(1500000,2700000)) +
  labs(x = "hap13", y = "hap26") +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(0.5,0.3,0.3,0.3), units = "lines"))

hap8vshap21_plot <- ggplot(data = hap8vshap21) +
  geom_segment(data = hap8vshap21, aes(x = qstart, xend = qend, y = tstart, yend = tend), linewidth = 0.1) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "% identity", limits = c(0,100), breaks = c(0,25,50,75,100), labels = c("0","25","50","75","100")) +
  labs(x = "hap13", y = "hap25") +
  coord_cartesian(xlim = c(2000000,2850000), ylim = c(1850000,2800000)) +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.4, units = "cm"),
        legend.key.width = unit(0.2, units = "cm"),
        plot.margin = unit(c(0.3,0.3,0.3,0.3), units = "lines"),
        legend.box.margin=margin(-5,-5,-5,-5))

hap6vshap23_plot <- ggplot(data = hap6vshap23) +
  geom_segment(data = hap6vshap23, aes(x = qstart, xend = qend, y = tstart, yend = tend), linewidth = 0.1) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "% identity", limits = c(0,100), breaks = c(0,25,50,75,100), labels = c("0","25","50","75","100")) +
  coord_cartesian(xlim = c(1600000,2200000), ylim = c(1600000,2400000)) +
  labs(x = "hap27", y = "hap26") +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 8),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        legend.title = element_text(color = "black", size = 6),
        legend.text = element_text(color = "black", size = 6),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(0.5,0.3,0.3,0.3), units = "lines"),
        legend.box.margin=margin(-5,-5,-5,-5))

combo_plot <- ggarrange(hap8vshap21_plot, hap6vshap23_plot, hap5vshap24_plot, ncol = 1, nrow = 3, heights = c(1,1,1))

ggsave("Figure2b.tiff", combo_plot, device = "tiff", width = 1.5, height = 4, units = "in", dpi = 600)
