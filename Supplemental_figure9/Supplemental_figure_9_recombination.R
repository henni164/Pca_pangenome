library(ggplot2)
library(dplyr)
library(ggh4x)
library(ggpubr)

options(scipen = 999)

## Supplemental figure 9e

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
recombinant_paf2$qname <- factor(recombinant_paf2$qname, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))

recombinant_paf2<- rbind(recombinant_paf2, data.frame(qname = "3", qlength = 6426316, qstart = 0, qend = 0, strand = "+", tname = "3", tlength = 6518986, tstart = 0, tend = 0, nmatches = 1, alnlength = 0, mquality = 60, perc = 100, hap = "hap4"))

p1 <- ggplot(data = recombinant_paf2) +
  geom_point(aes(x = tlength+50000), y = -0.5, color = "white") +
  geom_point(aes(x = tlength+50000), y = 1.5, color = "white") +
  geom_point(x = -50000, y = -0.5, color = "white") +
  geom_point(x = -50000, y = 1.5, color = "white") +
  geom_segment(aes(x = tstart, xend = tend, color = hap), y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= recombinant_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_grid(hap ~ qname, scales = "free_x", space = "free_x", switch = "y", margins = FALSE, drop = FALSE) +
  scale_x_continuous(expand = c(0.07,0.07)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4"),
                     values = c("#FA0348","#C581FC"), 
                     name = "Color key", drop = FALSE, limits = c("hap3","hap4")) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = NULL) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_text(size = 10, vjust = 1, hjust = 0.5, margin = margin(t = 1, b = 0, l = 0, r = 0)),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank())

p2 <- ggplot(data = recombinant_paf2) +
  geom_point(aes(x = qlength+50000), y = -0.5, color = "white") +
  geom_point(aes(x = qlength+50000), y = 1.5, color = "white") +
  geom_point(x = -50000, y = -0.5, color = "white") +
  geom_point(x = -50000, y = 1.5, color = "white") +
  geom_segment(aes(x = qstart, xend = qend, color = hap), y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= recombinant_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_grid( ~ tname, scales = "free_x", space = "free_x", switch = "y", margins = FALSE, drop = FALSE) +
  scale_x_continuous(expand = c(0.07,0.07)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4"),
                     values = c("#FA0348","#C581FC"), 
                     name = "Color key", drop = FALSE, limits = c("hap3","hap4")) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "hap25") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_blank(),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 1, vjust = 0.5))

recombinant_plot <- ggarrange(p1, p2, nrow = 2, ncol = 1, heights = c(2.2,1))

ggsave("Supplemental_figure9e.tiff", recombinant_plot, device = "tiff", width = 7.5, height = 2, units = "in", dpi = 600)