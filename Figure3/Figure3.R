options(scipen = 999)
library(ggplot2)
library(reshape2)
library(egg)

## Figure 3b chromosome lengths
chr_dat <- read.delim("Figure3_chromosome_lengths.txt", sep = "\t", header = TRUE)
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

ggsave("Figure3b.tiff", chr_lengths_plot, device = "tiff", height = 4.5, width = 5, units = "in", dpi = 600)


## Figure 3c chr9 region dotplots

hap25vshap13chr9 <- read.delim("hap25_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
hap25vshap13chr9 <- hap25vshap13chr9[,-(13:17)]
colnames(hap25vshap13chr9) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap25vshap13chr9 <- hap25vshap13chr9[grep("chr", hap25vshap13chr9$qname),]

hap26vshap27chr9 <- read.delim("hap26_vs_hap27_chr9.paf", sep = "\t", header = FALSE)
hap26vshap27chr9 <- hap26vshap27chr9[,-(13:17)]
colnames(hap26vshap27chr9) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")


hap26vshap13chr9 <- read.delim("hap26_vs_hap13_chr9.paf", sep = "\t", header = FALSE)
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

ggsave("Figure3c.tiff", combo_plot, device = "tiff", width = 1.8, height = 4.75, units = "in", dpi = 600)


## Figure 3d

hap8vshap21 <- read.delim("map_hap21_to_hap8.paf", sep = "\t", header = FALSE)
hap8vshap21 <- hap8vshap21[,-(13:17)]
colnames(hap8vshap21) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap8vshap21 <- hap8vshap21[grep("chr", hap8vshap21$qname),]

hap6vshap23 <- read.delim("map_hap23_to_hap6.paf", sep = "\t", header = FALSE)
hap6vshap23 <- hap6vshap23[,-(13:17)]
colnames(hap6vshap23) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")

hap8vshap21$qname <- factor(hap8vshap21$qname, levels = c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A"))
hap8vshap21$tname <- factor(hap8vshap21$tname, levels = rev(c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B")))
hap8vshap21$mquality <- as.numeric(hap8vshap21$mquality)
levels(hap8vshap21$qname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
levels(hap8vshap21$tname) <- rev(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))

hap6vshap23$qname <- factor(hap6vshap23$qname, levels = (c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A")))
hap6vshap23$tname <- factor(hap6vshap23$tname, levels = rev(c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B")))
hap6vshap23$mquality <- as.numeric(hap6vshap23$mquality)

levels(hap6vshap23$qname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
levels(hap6vshap23$tname) <- rev(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))



hap5vshap24 <- read.delim("map_hap24_to_hap5.paf", sep = "\t", header = FALSE)
hap5vshap24 <- hap5vshap24[,-(13:17)]
colnames(hap5vshap24) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap5vshap24 <- hap5vshap24[grep("chr", hap5vshap24$qname),]

hap5vshap24$tname <- factor(hap5vshap24$tname, levels = rev(c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A")))
hap5vshap24$qname <- factor(hap5vshap24$qname, levels = c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B"))
hap5vshap24$mquality <- as.numeric(hap5vshap24$mquality)
levels(hap5vshap24$qname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
levels(hap5vshap24$tname) <- rev(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))


hap5vshap24_plot <- ggplot(data = hap5vshap24) +
  geom_point(data = hap5vshap24[hap5vshap24$alnlength <= 100000,], aes(x = qstart, y = tstart, color = mquality), size = 0.25) +
  geom_segment(data = hap5vshap24[hap5vshap24$alnlength > 100000,], aes(x = qstart, xend = qend, y = tstart, yend = tend, color = ((nmatches/alnlength)*100)), linewidth = 1) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_continuous(position = "top", expand = c(0,0)) +
  facet_grid(tname ~ qname, scales = "free", drop = TRUE) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "%identity", labels = c("0","25","50","75","100")) +
  labs(x = "hap24", y = "hap5") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_blank(),
        #strip.text.y = element_text(angle = 0, color = "black", size = 6),
        strip.text.x = element_text(angle = 0, color = "black", size = 6),
        strip.background = element_rect(fill = "white", linewidth = 0.5, colour = "black"),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.2),
        plot.margin = unit(c(0,0.2,0,0), units = "lines"))

hap8vshap21_plot <- ggplot(data = hap8vshap21) +
  geom_point(data = hap8vshap21[hap8vshap21$alnlength <= 100000,], aes(x = qstart, y = tstart, color = mquality), size = 0.25) +
  geom_segment(data = hap8vshap21[hap8vshap21$alnlength > 100000,], aes(x = qstart, xend = qend, y = tstart, yend = tend, color = ((nmatches/alnlength)*100)), linewidth = 1) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_continuous(position = "top", expand = c(0,0)) +
  facet_grid(tname ~ qname, scales = "free", drop = TRUE) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "% identity", limits = c(0,100), breaks = c(0,25,50,75,100), labels = c("0","25","50","75","100")) +
  labs(x = "hap21", y = "hap8") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_blank(),
        #strip.text.y = element_text(angle = 0, color = "black", size = 6),
        strip.text.x = element_text(angle = 0, color = "black", size = 6),
        strip.background = element_rect(fill = "white", linewidth = 0.5, colour = "black"),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "nonw",
        #legend.direction = "horizontal",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.2),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.4, units = "cm"),
        legend.key.width = unit(0.6, units = "cm"),
        plot.margin = unit(c(0,0.2,0,0), units = "lines"),
        legend.box.margin=margin(-5,-5,-5,-5))

hap6vshap23_plot <- ggplot(data = hap6vshap23) +
  geom_point(data = hap6vshap23[hap6vshap23$alnlength <= 100000,], aes(x = qstart, y = tstart, color = mquality), size = 0.25) +
  geom_segment(data = hap6vshap23[hap6vshap23$alnlength > 100000,], aes(x = qstart, xend = qend, y = tstart, yend = tend, color = ((nmatches/alnlength)*100)), linewidth = 1) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_continuous(position = "top", expand = c(0,0)) +
  facet_grid(tname ~ qname, scales = "free", drop = TRUE) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "% identity", limits = c(0,100), labels = c("0","25","50","75","100")) +
  labs(x = "hap23", y = "hap6") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, color = "black", size = 6),
        strip.text.x = element_text(angle = 0, color = "black", size = 6),
        axis.title = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", linewidth = 0.5, colour = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.7, units = "cm"),
        legend.key.width = unit(0.6, units = "cm"),
        panel.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.2),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        legend.box.margin=margin(-5,-5,-5,-5))

combo_plot <- ggarrange(hap5vshap24_plot, hap8vshap21_plot, hap6vshap23_plot, ncol = 3, nrow = 1, widths = c(1,1,1))

ggsave("Figure3_d.tiff", combo_plot, device = "tiff", width = 7.5, height = 2.5, units = "in", dpi = 600)
