library(ggplot2)
library(dplyr)
library(ggh4x)

options(scipen = 999)

## Supplemental figure 6a

recombinant_paf <- read.delim("map_hap25_to_20NSW19.paf", sep = "\t", header = FALSE)
recombinant_paf <- recombinant_paf[,-(13:17)]
colnames(recombinant_paf) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
recombinant_paf <- recombinant_paf[grep("chr", recombinant_paf$qname),]

recombinant_paf$qname <- factor(recombinant_paf$qname, levels = rev(c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A")))
recombinant_paf$tname <- factor(recombinant_paf$tname, levels = c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A", "chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B"))
recombinant_paf$nmatches <- as.numeric(recombinant_paf$nmatches)
recombinant_paf$alnlength <- as.numeric(recombinant_paf$alnlength)
recombinant_paf$qlength <- as.numeric(recombinant_paf$qlength)
recombinant_paf$qstart <- as.numeric(recombinant_paf$qstart)
recombinant_paf$qend <- as.numeric(recombinant_paf$qend)
recombinant_paf$tstart <- as.numeric(recombinant_paf$tstart)
recombinant_paf$tend <- as.numeric(recombinant_paf$tend)
recombinant_paf$tlength <- as.numeric(recombinant_paf$tlength)

levels(recombinant_paf$qname) <- rev(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))
levels(recombinant_paf$tname) <- c("1A","2A","3A","4A","5A","6A","7A","8A","9A","10A","11A","12A","13A","14A","15A","16A","17A","18A","1B","2B","3B","4B","5B","6B","7B","8B","9B","10B","11B","12B","13B","14B","15B","16B","17B","18B")
recombinant_paf$perc <- (recombinant_paf$nmatches/recombinant_paf$alnlength*100)
recombinant_paf2 <- recombinant_paf %>% mutate(hap = case_when(grepl("A", tname) ~ "hap3", grepl("B", tname) ~ "hap4"))
levels(recombinant_paf2$tname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")

recombinant_plot <- ggplot(data = recombinant_paf2) +
  geom_point(data = recombinant_paf2[recombinant_paf2$alnlength <= 100000,], aes(x = qstart, y = tstart, color = perc), size = 0.25) +
  geom_point(data = recombinant_paf2, aes(x = 1, y = 1), color = "white", size = 0.01) +
  geom_point(data = recombinant_paf2, aes(x = qlength, y = tlength), color = "white", size = 0.01) +
  geom_segment(data = recombinant_paf2[recombinant_paf2$alnlength > 100000,], aes(x = qstart, xend = qend, y = tstart, yend = tend, color = perc), linewidth = 1) +
  facet_nested(qname ~ hap + tname, scales = "free", drop = TRUE) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_continuous(position = "top", expand = c(0,0)) +
  scale_color_gradient(low = "#D5D5D5", high = "#000000", aesthetics = "color", name = "% identity", limits = c(0,100), breaks = c(0,25,50,75,100), labels = c("0","25","50","75","100")) +
  labs(x = "20NSW19", y = "20WA94 hap25") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, color = "black", size = 6),
        strip.text.x = element_text(angle = 0, color = "black", size = 6, margin = margin(0.2,0.2,0.2,0.2)),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        axis.title.x.top = element_text(color = "black", size = 6),
        axis.title.y.right = element_text(color = "black", size = 6),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(0,0,0,0), units = "lines"))

ggsave("Supplemental_figure6a.tiff", recombinant_plot, device = "tiff", width = 7, height = 4, units = "in", dpi = 600)

## Supplemental figure 6b

hap5_vs_hap26_paf <- read.delim("map_hap5_to_20WA94.paf", sep = "\t", header = FALSE)
hap5_vs_hap26_paf <- hap5_vs_hap26_paf[,-(13:17)]
hap5_vs_hap26_paf <- dplyr::filter(hap5_vs_hap26_paf, grepl("_B", V6))

colnames(hap5_vs_hap26_paf) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap5_vs_hap26_paf <- hap5_vs_hap26_paf[grep("chr", hap5_vs_hap26_paf$qname),]

hap5_vs_hap26_paf$tname <- factor(hap5_vs_hap26_paf$tname, levels = rev(c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B")))
hap5_vs_hap26_paf$qname <- factor(hap5_vs_hap26_paf$qname, levels = c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A"))
hap5_vs_hap26_paf$nmatches <- as.numeric(hap5_vs_hap26_paf$nmatches)
hap5_vs_hap26_paf$alnlength <- as.numeric(hap5_vs_hap26_paf$alnlength)
hap5_vs_hap26_paf$qlength <- as.numeric(hap5_vs_hap26_paf$qlength)
hap5_vs_hap26_paf$qstart <- as.numeric(hap5_vs_hap26_paf$qstart)
hap5_vs_hap26_paf$qend <- as.numeric(hap5_vs_hap26_paf$qend)
hap5_vs_hap26_paf$tstart <- as.numeric(hap5_vs_hap26_paf$tstart)
hap5_vs_hap26_paf$tend <- as.numeric(hap5_vs_hap26_paf$tend)
hap5_vs_hap26_paf$tlength <- as.numeric(hap5_vs_hap26_paf$tlength)

levels(hap5_vs_hap26_paf$qname) <- rev(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))
levels(hap5_vs_hap26_paf$tname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
hap5_vs_hap26_paf$perc <- (hap5_vs_hap26_paf$nmatches/hap5_vs_hap26_paf$alnlength*100)

hap5_vs_hap26_paf2 <- hap5_vs_hap26_paf[hap5_vs_hap26_paf$perc >= 95,]

panelbplot <- ggplot(data = hap5_vs_hap26_paf2) +
  geom_segment(data = hap5_vs_hap26_paf2, aes(x = 1, xend = qlength, y = 1, yend = tlength), color = "white") +
  geom_point(data = hap5_vs_hap26_paf2[hap5_vs_hap26_paf2$alnlength <= 100000,], aes(x = qstart, y = tstart), size = 0.25) +
  geom_segment(data = hap5_vs_hap26_paf2[hap5_vs_hap26_paf2$alnlength > 100000,], aes(x = qstart, xend = qend, y = tstart, yend = tend), linewidth = 1) +
  facet_nested(qname ~ tname, scales = "free", drop = TRUE) +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_x_continuous(position = "top", expand = c(0,0)) +
  labs(x = "20QLD86 hap5", y = "20WA94 hap26") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, color = "black", size = 6),
        strip.text.x = element_text(angle = 0, color = "black", size = 6, margin = margin(0.2,0.2,0.2,0.2)),
        strip.background = element_rect(fill = "white", linewidth = 0.1, colour = "black"),
        axis.title.x.top = element_text(color = "black", size = 6),
        axis.title.y.right = element_text(color = "black", size = 6),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed", linewidth = 0.1),
        plot.margin = unit(c(0,0,0,0), units = "lines"))

ggsave("Supplemental_figure6b.tiff", panelbplot, device = "tiff", width = 4, height = 4, units = "in", dpi = 1200)

