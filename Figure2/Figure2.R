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

hap25vshap13chr9 <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_population_study/Mating_loci/hap25_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
hap25vshap13chr9 <- hap25vshap13chr9[,-(13:17)]
colnames(hap25vshap13chr9) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap25vshap13chr9 <- hap25vshap13chr9[grep("chr", hap25vshap13chr9$qname),]

hap26vshap27chr9 <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_population_study/Mating_loci/hap26_vs_hap27_chr9.paf", sep = "\t", header = FALSE)
hap26vshap27chr9 <- hap26vshap27chr9[,-(13:17)]
colnames(hap26vshap27chr9) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")


hap26vshap13chr9 <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/Research/Current_projects/Pca_population_study/Mating_loci/hap26_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
hap26vshap13chr9 <- hap26vshap13chr9[,-(13:17)]
colnames(hap26vshap13chr9) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap26vshap13chr9 <- hap26vshap13chr9[grep("chr", hap26vshap13chr9$qname),]

hap25vshap13chr9_plot <- ggplot(data = hap25vshap13chr9) +
  geom_segment(data = hap25vshap13chr9, aes(x = qstart, xend = qend, y = tstart, yend = tend), color = "black", linewidth = 0.3) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000, breaks = c(2000000,2300000,2600000), sec.axis = dup_axis()) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000, breaks = c(2100000,2400000,2700000), sec.axis = dup_axis()) +
  labs(x = "hap13", y = "hap25") +
  coord_cartesian(xlim = c(2000000,2850000), ylim = c(1850000,2800000)) +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title.y.right = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.y.left = element_blank(),
        axis.title.x.top = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.x.bottom = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3, linetype = "dashed"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(b = 4, r = 0.1, l = 0.1))

hap26vshap27chr9_plot <- ggplot(data = hap26vshap27chr9) +
  geom_segment(data = hap26vshap27chr9, aes(x = qstart, xend = qend, y = tstart, yend = tend), color = "black", linewidth = 0.3) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000, breaks = c(1700000,2000000,2300000), sec.axis = dup_axis()) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000, breaks = c(1700000,1900000,2100000), sec.axis = dup_axis()) +
  coord_cartesian(xlim = c(1600000,2200000), ylim = c(1600000,2400000)) +
  labs(x = "hap27", y = "hap26") +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black", size = 8),
        axis.title.y.right = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.y.left = element_blank(),
        axis.title.x.top = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.x.bottom = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3, linetype = "dashed"),
        plot.margin = margin(b = 2, t = 2, r = 1, l = 1))

hap26vshap13chr9_plot <- ggplot(data = hap26vshap13chr9) +
  geom_segment(data = hap26vshap13chr9, aes(x = qstart, xend = qend, y = tstart, yend = tend), color = "black", linewidth = 0.3) +
  scale_y_continuous(position = "left", expand = c(0,0), labels = function(x)x/1000000, breaks = c(1700000,2100000,2500000), sec.axis = dup_axis()) +
  scale_x_continuous(position = "bottom", expand = c(0,0), labels = function(x)x/1000000, breaks = c(1800000,2300000,2800000), sec.axis = dup_axis()) +
  coord_cartesian(xlim = c(1600000,3000000), ylim = c(1500000,2700000)) +
  labs(x = "hap13", y = "hap26") +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title.y.right = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.y.left = element_blank(),
        axis.title.x.top = element_text(color = "black", size = 10, margin = margin(l = 0.01)),
        axis.title.x.bottom = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3, linetype = "dashed"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(t = 4, r = 1, l = 1))

combo_plot <- ggarrange(hap25vshap13chr9_plot, hap26vshap27chr9_plot, hap26vshap13chr9_plot, ncol = 1, nrow = 3, heights = c(1,1,1))

ggsave("Figure2b.tiff", combo_plot, device = "tiff", width = 1.8, height = 4.75, units = "in", dpi = 600)
