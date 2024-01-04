library(ggplot2)
library(egg)

options(scipen = 999)

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
        strip.text.y = element_text(angle = 0, color = "black", size = 8),
        strip.text.x = element_text(angle = 0, color = "black", size = 8),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(0,0,0,0), units = "lines"))

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
        strip.text.y = element_text(angle = 0, color = "black", size = 8),
        strip.text.x = element_blank(),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.4, units = "cm"),
        legend.key.width = unit(0.2, units = "cm"),
        plot.margin = unit(c(1,0,0,0), units = "lines"),
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
        strip.text.y = element_text(angle = 0, color = "black", size = 8),
        strip.text.x = element_blank(),
        axis.title = element_text(color = "black", size = 8),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        legend.title = element_text(color = "black", size = 6),
        legend.text = element_text(color = "black", size = 6),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(1,0,0,0), units = "lines"))

combo_plot <- ggarrange(hap5vshap24_plot, hap8vshap21_plot, hap6vshap23_plot, ncol = 1, nrow = 3, heights = c(1,1,1))

ggsave("Supplemental_figure6_hybrid_alignments.tiff", combo_plot, device = "tiff", width = 3.5, height = 8, units = "in", dpi = 600)
