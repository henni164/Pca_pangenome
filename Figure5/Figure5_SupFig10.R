options(scipen = 999)
library(ggplot2)
library(dplyr)
library(reshape2)
library(purrr)
library(ggchicklet)

## Figure 5 & Supplemental Figure 10

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

merged_intervals <- read.delim("all_intervals.txt", sep = "\t", header = FALSE)
colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
keepme <- final_joined

keepme <- subset(keepme, !(HAP %in% c("hap21","hap22","hap24")))
keepme <- subset(keepme, !(REFHAP %in% c("hap21","hap22","hap24")))

keepme$REFHAP <- factor(keepme$REFHAP)

keepme$HAP <- factor(keepme$HAP, levels = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap13","hap11","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap10","hap12","hap23","hap9"))
keepme$CHROM <- factor(keepme$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

chrlen_melt <- subset(chrlen_melt, !(HAP %in% c("hap21","hap22","hap24")))
chrlen_melt$HAP <- factor(chrlen_melt$HAP, levels = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap13","hap11","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap10","hap12","hap23","hap9"))


## Figure 5

chrlen_aus <- subset(chrlen_melt, !(HAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24")))

aus_intervals <- subset(merged_intervals, !(HAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24")))
aus_intervals <- subset(aus_intervals, !(REFHAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24")))

final_aus <- full_join(aus_intervals, chrlen_aus, by = c("HAP","CHROM"))
final_aus$HAP <- factor(final_aus$HAP)
just_aus <- final_aus

just_aus$HAP <- factor(just_aus$HAP, levels = c("hap26","hap5","hap31","hap32","hap7","hap8","hap29","hap30","hap27","hap28","hap25","hap3","hap4","hap6","hap23"))
just_aus$CHROM <- factor(just_aus$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
chrlen_aus$HAP <- factor(chrlen_aus$HAP, levels = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap6","hap23"))
chrlen_aus$HAP <- factor(chrlen_aus$HAP, levels = c("hap26","hap5","hap31","hap32","hap7","hap8","hap29","hap30","hap27","hap28","hap25","hap3","hap4","hap6","hap23"))
just_aus$REFHAP <- factor(just_aus$REFHAP)


plot2 <- ggplot(data = just_aus, aes(y = HAP)) +
  geom_point(aes(x = len+15000), y = -0.5, color = "white") +
  geom_point(aes(x = len+15000), y = 1.5, color = "white") +
  geom_point(x = -15000, y = -0.5, color = "white") +
  geom_point(x = -15000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = chrlen_aus, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.3, 'npc')) +
  geom_segment(aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 3.7, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= chrlen_aus, aes(xmax = len+15000), xmin = -14999, ymin = -0.1, ymax = 1.1, r = unit(0.3, 'npc'), color = "white", fill = NA, size = 0.3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free_x", space = "free_x", switch = "both", margins = FALSE) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8",
                                "hap27","hap28","hap29","hap30","hap31","hap32","hap6","hap23"),
                     values = c("#FF0000","#FFFB00","black","#FF8000","#30D500","#00FFE8","#002EFF",
                                "#AE00FF","black","black","black","black","black","black","black"), name = "Color key", drop = TRUE, limits = c("hap3","hap4","hap26","hap5","hap7","hap8"), na.value = "black") +
  scale_fill_manual(breaks = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8",
                               "hap27","hap28","hap29","hap30","hap31","hap32","hap6","hap23"),
                    values = c("#FF0000","#FFFB00","black","#FF8000","#30D500","#00FFE8","#002EFF",
                               "#AE00FF","black","black","black","black","black","black","black"), name = "Color key", drop = TRUE, guide = "none") +
  labs(x = "Chromosome", y = NULL) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7, hjust = 1),
        strip.text.x.bottom = element_text(size = 7, vjust = 1),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.15), units = "in"),
        plot.margin = unit(c(0,0,2.25,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(-0.07,-0.1),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 1, fill = NA, linetype = 0)))

ggsave("Figure5_pangenome_blocks.tiff", plot2, device = "tiff", width = 7, height = 6, units = "in", dpi = 600)


## Supplemental Figure 10

plot <- ggplot(data = keepme, aes(y = HAP)) +
  geom_point(aes(x = len+50000), y = -0.5, color = "white") +
  geom_point(aes(x = len+50000), y = 1.5, color = "white") +
  geom_point(x = -50000, y = -0.5, color = "white") +
  geom_point(x = -50000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = chrlen_melt, aes(xmax = len, color = HAP, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.3, 'npc')) +
  geom_segment(aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 3.3, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= chrlen_melt, aes(xmax = len+50000), xmin = -49999, ymin = -0.1, ymax = 1.1, r = unit(0.3, 'npc'), color = "white", fill = NA, size = 0.35) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free_x", space = "free_x", switch = "both", margins = FALSE) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8",
                                "hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap13",
                                "hap11","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap10","hap12","hap23","hap9"),
                     values = c("#FF0000","#FFFB00","black","#FF8000","#30D500","#00FFE8","#002EFF",
                                "#AE00FF","black","black","black","black","black","#B8E9FF","#FF00F7","#006023","#FFCD97",
                                "#DDB8FF","#CFFFB8","#45C48C","#AEB431","#00B5D2","black","black","black","black","black","black","black"), name = "Color key", drop = TRUE) +
  scale_fill_manual(breaks = c("hap3","hap4","hap25","hap26","hap5","hap7","hap8",
                               "hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap13",
                               "hap11","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap10","hap12","hap23","hap9"),
                    values = c("#FF0000","#FFFB00","black","#FF8000","#30D500","#00FFE8","#002EFF",
                               "#AE00FF","black","black","black","black","black","#B8E9FF","#FF00F7","#006023","#FFCD97",
                               "#DDB8FF","#CFFFB8","#45C48C","#AEB431","#00B5D2","black","black","black","black","black","black","black"), name = "Color key", drop = TRUE, guide = "none") +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = "Chromosome", y = NULL) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7, hjust = 1),
        strip.text.x.bottom = element_text(size = 7, vjust = 1),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.15), units = "in"),
        plot.margin = unit(c(0,0,5,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(-0.06,-0.11),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_blank()) +
  guides(color = guide_legend(nrow = 5, override.aes = list(size = 1, fill = NA, linetype = 0)))

ggsave("Supplemental_figure10_pangenome_blocks.tiff", plot, device = "tiff", width = 7, height = 9, units = "in", dpi = 600)
