library(ggplot2)
library(dplyr)
library(ggh4x)
library(cowplot)

options(scipen = 999)

## Supplemental figure 10a

hap5_vs_hap26_paf <- read.delim("map_hap5_to_20WA94.paf", sep = "\t", header = FALSE)
hap5_vs_hap26_paf <- hap5_vs_hap26_paf[,-(13:17)]
hap5_vs_hap26_paf <- dplyr::filter(hap5_vs_hap26_paf, grepl("_B", V6))

colnames(hap5_vs_hap26_paf) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
hap5_vs_hap26_paf <- hap5_vs_hap26_paf[grep("chr", hap5_vs_hap26_paf$qname),]
hap5_vs_hap26_paf$tname <- as.factor(hap5_vs_hap26_paf$tname)
hap5_vs_hap26_paf$qname <- as.factor(hap5_vs_hap26_paf$qname)
hap5_vs_hap26_paf$tname <- factor(hap5_vs_hap26_paf$tname, levels = c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B", "chr7_B", "chr8_B", "chr9_B", "chr10_B", "chr11_B", "chr12_B", "chr13_B", "chr14_B", "chr15_B", "chr16_B", "chr17_B", "chr18_B"))
hap5_vs_hap26_paf$qname <- factor(hap5_vs_hap26_paf$qname, levels = c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A"))
hap5_vs_hap26_paf$nmatches <- as.numeric(hap5_vs_hap26_paf$nmatches)
hap5_vs_hap26_paf$alnlength <- as.numeric(hap5_vs_hap26_paf$alnlength)
hap5_vs_hap26_paf$qlength <- as.numeric(hap5_vs_hap26_paf$qlength)
hap5_vs_hap26_paf$qstart <- as.numeric(hap5_vs_hap26_paf$qstart)
hap5_vs_hap26_paf$qend <- as.numeric(hap5_vs_hap26_paf$qend)
hap5_vs_hap26_paf$tstart <- as.numeric(hap5_vs_hap26_paf$tstart)
hap5_vs_hap26_paf$tend <- as.numeric(hap5_vs_hap26_paf$tend)
hap5_vs_hap26_paf$tlength <- as.numeric(hap5_vs_hap26_paf$tlength)

levels(hap5_vs_hap26_paf$qname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
levels(hap5_vs_hap26_paf$tname) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
hap5_vs_hap26_paf$perc <- (hap5_vs_hap26_paf$nmatches/hap5_vs_hap26_paf$alnlength*100)

hap5_vs_hap26_paf2 <- hap5_vs_hap26_paf[hap5_vs_hap26_paf$perc >= 95,]


p1 <- ggplot(data = hap5_vs_hap26_paf2) +
  geom_segment(aes(x = tstart, xend = tend), color = "#0027D9", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = "#0027D9", size = 0.35) +
  facet_wrap(~ tname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "hap5") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5))

p1v2 <- ggplot(data = hap5_vs_hap26_paf2) +
  geom_segment(aes(x = tstart, xend = tend), color = "#0027D9", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_wrap(~ tname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "hap5") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5))


p2 <- ggplot(data = hap5_vs_hap26_paf2) +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = NA, fill = NA, size = 0.35) +
  geom_segment(aes(x = qstart, xend = qend), color = "#0027D9", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_wrap( ~ qname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "hap26") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5))

recombinant_plot <- cowplot::plot_grid(p1v2 ,p2, align = "hv", nrow = 2, ncol = 1, rel_heights = c(1,1))
ggsave("Supplemental_figure10a.tiff", recombinant_plot, device = "tiff", width = 7.5, height = 1.25, units = "in", dpi = 600)

#panel b

p3 <- ggplot(data = hap5_vs_hap26_paf2) +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = NA, fill = "#6CAD00", size = 0.35) +
  geom_segment(aes(x = qstart, xend = qend), color = "#0027D9", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_wrap(~ qname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "hap26") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 1, vjust = 0.5))

recombinant_plot2 <- cowplot::plot_grid(p1,p3, align = "hv", nrow = 2, ncol = 1, rel_heights = c(1,1))
ggsave("Supplemental_figure10b.tiff", recombinant_plot2, device = "tiff", width = 7.5, height = 1.25, units = "in", dpi = 600)


p4 <- ggplot(data = hap5_vs_hap26_paf2) +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = NA, fill = "#6CAD00", size = 0.35) +
  geom_segment(aes(x = qstart, xend = qend), color = "grey80", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_wrap(~ qname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "?") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 1, vjust = 0.5))

recombinant_plot3 <- cowplot::plot_grid(p1,p4, align = "hv", nrow = 2, ncol = 1, rel_heights = c(1,1))
ggsave("Supplemental_figure10c.tiff", recombinant_plot3, device = "tiff", width = 7.5, height = 1.25, units = "in", dpi = 600)

p5 <- ggplot(data = hap5_vs_hap26_paf2) +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = NA, fill = "#0027D9", size = 0.35) +
  geom_segment(aes(x = qstart, xend = qend), color = "grey80", y = 0.5, yend = 0.5, linewidth = 9.5, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= hap5_vs_hap26_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_wrap(~ tname, scales = "fixed", ncol = 18, drop = FALSE) +
  scale_x_continuous(expand = c(0.05,0.05)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = "?") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x.top = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.left = element_text(angle = 0, size = 10, hjust = 1, vjust = 0.5))

recombinant_plot4 <- cowplot::plot_grid(p5,p3, align = "hv", nrow = 2, ncol = 1, rel_heights = c(1,1))
ggsave("Supplemental_figure10d.tiff", recombinant_plot4, device = "tiff", width = 7.5, height = 1.25, units = "in", dpi = 600)
