library(ggplot2)
library(reshape2)
library(tidyr)
library(ggpubr)


hapAvars <- read.delim("hapA_HD_vars.vcf")

colnames(hapAvars) <- c("CHROM","POS","REF","ALT","Pca2016-15-1","90WI131-1","Pca77-46","17TX62-3","17MNBT-23","16MN105-1","17MNBT-14","17MN150-2","17LA50-3","PCA16-15-1-1","16MN85-1","17FL16-1","17AR69-3","17IA106-3","15OH12-3","15SD11-2","17LA50-2","15ND20-3","15MS7-1","17TX44-3","17TX15-3","18MNBT-39","17TX59-2","17MNBT-11","16MN84-1","17SD132-1","16MN94-1","16MN99-1","18MNBT-56","17MNBT-7","17TX9-2","16MN98-2","15MN23-1","15MN16-3","16MN107-1","16MN102-3","17MNBT-20","16MN107-3","17MN78-2","16SD74-2","18MNBT-37","18MNBT-41","18MNBT-46","17SD126-3","90MN5B-1","18MNBT-55","17MN150-1","17MNBT-6","17TX12-3","17SD133-2")

hapA_sepped <- separate_wider_delim(hapAvars, c("Pca2016-15-1","90WI131-1","Pca77-46","17TX62-3","17MNBT-23","16MN105-1","17MNBT-14","17MN150-2","17LA50-3","PCA16-15-1-1","16MN85-1","17FL16-1","17AR69-3","17IA106-3","15OH12-3","15SD11-2","17LA50-2","15ND20-3","15MS7-1","17TX44-3","17TX15-3","18MNBT-39","17TX59-2","17MNBT-11","16MN84-1","17SD132-1","16MN94-1","16MN99-1","18MNBT-56","17MNBT-7","17TX9-2","16MN98-2","15MN23-1","15MN16-3","16MN107-1","16MN102-3","17MNBT-20","16MN107-3","17MN78-2","16SD74-2","18MNBT-37","18MNBT-41","18MNBT-46","17SD126-3","90MN5B-1","18MNBT-55","17MN150-1","17MNBT-6","17TX12-3","17SD133-2"), delim = "/", names_sep = "_")


hapA_melted <- melt(hapA_sepped, id.vars = c("CHROM", "POS", "REF", "ALT"), value.name = "genotype")
hapA_final <- separate_wider_delim(hapA_melted, variable, delim = "_", names = c("Sample", "Allele"))

hapA_final$Sample <- factor(hapA_final$Sample, levels = c("PCA16-15-1-1", "Pca2016-15-1", "Pca77-46", "17MNBT-14", "17SD133-2",
                                                          "16SD74-2","18MNBT-56",
                                                          "90MN5B-1","15MN23-1","16MN105-1","17MNBT-6","17TX9-2","18MNBT-41",
                                                          "16MN102-3","16MN99-1","17SD126-3","18MNBT-55","15MN16-3","16MN107-3",
                                                          "15OH12-3","16MN85-1","16MN98-2","17TX62-3",
                                                          "17AR69-3","17MNBT-11","17SD132-1","90WI131-1",
                                                          "15SD11-2","17MN150-1","17MN150-2","17MNBT-23","17MNBT-7","18MNBT-46","18MNBT-37",
                                                          "15MS7-1","17IA106-3",
                                                          "17LA50-2", "17LA50-3", "17FL16-1","17TX12-3", "17TX15-3","17TX44-3","17TX59-2","15ND20-3","17MN78-2",
                                                          "16MN107-1","16MN94-1","16MN84-1","17MNBT-20","18MNBT-39"))

hapA_final$genotype <- factor(hapA_final$genotype, levels = c(".","0","1","2","3","4","5"))

hapA_final_sub <- subset(hapA_final, POS >= 4480818 & POS <= 4484474)


hapBvars <- read.delim("hapB_HD_vars.vcf")

colnames(hapBvars) <- c("CHROM","POS","REF","ALT","Pca2016-15-1","90WI131-1","Pca77-46","17TX62-3","17MNBT-23","16MN105-1","17MNBT-14","17MN150-2","17LA50-3","PCA16-15-1-1","16MN85-1","17FL16-1","17AR69-3","17IA106-3","15OH12-3","15SD11-2","17LA50-2","15ND20-3","15MS7-1","17TX44-3","17TX15-3","18MNBT-39","17TX59-2","17MNBT-11","16MN84-1","17SD132-1","16MN94-1","16MN99-1","18MNBT-56","17MNBT-7","17TX9-2","16MN98-2","15MN23-1","15MN16-3","16MN107-1","16MN102-3","17MNBT-20","16MN107-3","17MN78-2","16SD74-2","18MNBT-37","18MNBT-41","18MNBT-46","17SD126-3","90MN5B-1","18MNBT-55","17MN150-1","17MNBT-6","17TX12-3","17SD133-2")

hapB_sepped <- separate_wider_delim(hapBvars, c("Pca2016-15-1","90WI131-1","Pca77-46","17TX62-3","17MNBT-23","16MN105-1","17MNBT-14","17MN150-2","17LA50-3","PCA16-15-1-1","16MN85-1","17FL16-1","17AR69-3","17IA106-3","15OH12-3","15SD11-2","17LA50-2","15ND20-3","15MS7-1","17TX44-3","17TX15-3","18MNBT-39","17TX59-2","17MNBT-11","16MN84-1","17SD132-1","16MN94-1","16MN99-1","18MNBT-56","17MNBT-7","17TX9-2","16MN98-2","15MN23-1","15MN16-3","16MN107-1","16MN102-3","17MNBT-20","16MN107-3","17MN78-2","16SD74-2","18MNBT-37","18MNBT-41","18MNBT-46","17SD126-3","90MN5B-1","18MNBT-55","17MN150-1","17MNBT-6","17TX12-3","17SD133-2"), delim = "/", names_sep = "_")


hapB_melted <- melt(hapB_sepped, id.vars = c("CHROM", "POS", "REF", "ALT"), value.name = "genotype")
hapB_final <- separate_wider_delim(hapB_melted, variable, delim = "_", names = c("Sample", "Allele"))

hapB_final$Sample <- factor(hapB_final$Sample, levels = c("PCA16-15-1-1", "Pca2016-15-1", "Pca77-46", "17MNBT-14", "17SD133-2",
                                                          "16SD74-2","18MNBT-56",
                                                          "90MN5B-1","15MN23-1","16MN105-1","17MNBT-6","17TX9-2","18MNBT-41",
                                                          "16MN102-3","16MN99-1","17SD126-3","18MNBT-55","15MN16-3","16MN107-3",
                                                          "15OH12-3","16MN85-1","16MN98-2","17TX62-3",
                                                          "17AR69-3","17MNBT-11","17SD132-1","90WI131-1",
                                                          "15SD11-2","17MN150-1","17MN150-2","17MNBT-23","17MNBT-7","18MNBT-46","18MNBT-37",
                                                          "15MS7-1","17IA106-3",
                                                          "17LA50-2", "17LA50-3", "17FL16-1","17TX12-3", "17TX15-3","17TX44-3","17TX59-2","15ND20-3","17MN78-2",
                                                          "16MN107-1","16MN94-1","16MN84-1","17MNBT-20","18MNBT-39"))

hapB_final$genotype <- factor(hapB_final$genotype, levels = c(".","0","1","2","3","4","5"))

hapB_final_sub <- subset(hapB_final, POS >= 5295260 & POS <= 5298876)

hapA <- ggplot(data = hapA_final_sub) + 
  geom_segment(x = 4480820, xend = 4484472, y = 0, yend = 0, lwd = 1) +
  geom_rect(xmin = 4480820, xmax = 4482137, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") + #lineend parameter needed to fix a weird graphics processing bug on windows
  geom_rect(xmin = 4482200, xmax = 4482612, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 4483196, xmax = 4483834, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 4483916, xmax = 4484047, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 4484118, xmax = 4484472, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_segment(data = subset(hapA_final_sub, Allele == "1" & genotype == "0"), aes(x = POS, xend = POS, color = genotype), y = -0.3, yend = 0, lwd = 0.3, alpha = 0) +
  geom_segment(data = subset(hapA_final_sub, Allele == "1" & genotype != "0"), aes(x = POS, xend = POS, color = genotype), y = -0.3, yend = 0, lwd = 0.3) +
  geom_segment(data = subset(hapA_final_sub, Allele == "2" & genotype == "0"), aes(x = POS, xend = POS, color = genotype), y = 0, yend = 0.3, lwd = 0.3, alpha = 0) +
  geom_segment(data = subset(hapA_final_sub, Allele == "2" & genotype != "0"), aes(x = POS, xend = POS, color = genotype), y = 0, yend = 0.3, lwd = 0.3) +
  scale_color_manual(breaks = c(".","0","1","2","3","4","5"), values = c("white","black","#FF00D1","#47D45A","#00B0F0","#EA753A","#FBFF00")) +
  #scale_color_manual(breaks = c(".","0","1","2","3","4","5"), values = c("white","black","#FBFF00","#FBFF00","#FBFF00","#FBFF00","#FBFF00")) +
  scale_y_continuous(limits = c(-0.4,0.4)) +
  labs(x = NULL, y = NULL) +
  facet_grid(rows = vars(Sample), switch = "y", margins = FALSE) +
  theme(panel.background = element_rect(fill ="white"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, margin = margin(0,0,0,0)),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0.1,0,0.1,0))


hapB <- ggplot(data = hapB_final_sub) + 
  geom_segment(x = 5295262, xend = 5298874, y = 0, yend = 0, lwd = 1) +
  geom_rect(xmin = 5295262, xmax = 5296585, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 5296644, xmax = 5297050, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 5297604, xmax = 5298236, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 5298318, xmax = 5298449, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_rect(xmin = 5298521, xmax = 5298874, ymin = -0.3, ymax = 0.3, fill = "black", color = "black", lwd = 1, lineend = "square") +
  geom_segment(data = subset(hapB_final_sub, Allele == "1" & genotype == "0"), aes(x = POS, xend = POS, color = genotype), y = -0.3, yend = 0, lwd = 0.3, alpha = 0) +
  geom_segment(data = subset(hapB_final_sub, Allele == "1" & genotype != "0"), aes(x = POS, xend = POS, color = genotype), y = -0.3, yend = 0, lwd = 0.3) +
  geom_segment(data = subset(hapB_final_sub, Allele == "2" & genotype == "0"), aes(x = POS, xend = POS, color = genotype), y = 0, yend = 0.3, lwd = 0.3, alpha = 0) +
  geom_segment(data = subset(hapB_final_sub, Allele == "2" & genotype != "0"), aes(x = POS, xend = POS, color = genotype), y = 0, yend = 0.3, lwd = 0.3) +
  scale_color_manual(breaks = c(".","0","1","2","3","4","5"), values = c("white","black","#FF00D1","#47D45A","#00B0F0","#EA753A","#FBFF00")) +
  #scale_color_manual(breaks = c(".","0","1","2","3","4","5"), values = c("white","black","#FBFF00","#FBFF00","#FBFF00","#FBFF00","#FBFF00")) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  facet_grid(rows = vars(Sample), switch = "y", margins = FALSE) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(0.1,0,0.1,0),
        legend.margin = margin(0,0,0,0.1),
        legend.key = element_rect(fill = "black"),
        legend.text = element_text(size = 8))


combo <- ggarrange(hapA, hapB, ncol = 2, nrow = 1, widths = c(1.3,1))

ggsave("Figure_S7_HD_locus_snps_binary.tiff", combo, device = "tiff", width = 6.25, height = 9.5, units = "in", dpi = 600)

