options(scipen = 999)
library(ggplot2)
library(dplyr)
library(reshape2)
library(purrr)
library(ggchicklet)

## Fig6A, Fig6B, FigS10

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

haps_chrs <- read.delim("FigureS4A_chromosome_lengths.txt", sep = "\t", header = TRUE)
colnames(haps_chrs) <- c("HAP","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")
chrlen_melt <-melt(haps_chrs, id.vars = "HAP", variable.name = "CHROM", value.name = "len")
chrlen_melt$CHROM <- factor(chrlen_melt$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

merged_intervals <- read.delim("all_intervals_fig6A.txt", sep = "\t", header = FALSE)

colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
keepme <- final_joined

keepme <- subset(keepme, !(HAP %in% c("hap21","hap22","hap24")))
keepme <- subset(keepme, !(REFHAP %in% c("hap21","hap22","hap24")))

keepme$REFHAP <- factor(keepme$REFHAP)

keepme$HAP <- factor(keepme$HAP, levels = c("hap11","hap12","hap15","hap16","hap17","hap18","hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap1","hap2","hap10","hap19","hap20","hap13","hap14","hap9"))
keepme$CHROM <- factor(keepme$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

chrlen_melt <- subset(chrlen_melt, !(HAP %in% c("hap21","hap22","hap24")))
chrlen_melt$HAP <- factor(chrlen_melt$HAP, levels = c("hap11","hap12","hap15","hap16","hap17","hap18","hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap1","hap2","hap10","hap19","hap20","hap13","hap14","hap9"))


## Figure 6

## Figure 6A

fig6A_chr_sub <- subset(chrlen_melt, HAP %in% c("hap11","hap12","hap15","hap16","hap17","hap18"))

fig6A_sub <- subset(keepme, HAP %in% c("hap11","hap12","hap15","hap16","hap17","hap18"))
fig6A_sub <- subset(fig6A_sub, REFHAP %in% c("hap11","hap12","hap15","hap16","hap17","hap18"))


plot <- ggplot(data = fig6A_sub, aes(y = HAP)) +
  geom_point(aes(x = len+15000), y = -0.5, color = "white") +
  geom_point(aes(x = len+15000), y = 1.5, color = "white") +
  geom_point(x = -15000, y = -0.5, color = "white") +
  geom_point(x = -15000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = fig6A_chr_sub, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.1, 'npc')) +
  geom_segment(data = na.omit(fig6A_sub), aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 5.7, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= fig6A_chr_sub, aes(xmax = len+15000), xmin = -14999, ymin = -0.01, ymax = 1.01, r = unit(0.1, 'npc'), color = "black", fill = NA, size = 0.3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free_x", space = "free_x", switch = "y", margins = FALSE) +
  scale_x_continuous(expand = c(0.12,0.12)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap11","hap12","hap15","hap16","hap17","hap18"),
                     values = c("#DDB8FF","grey30","#CFFFB8","#71E2ED","#AEB431","#6BA9B4"), name = "Color key", drop = TRUE, limits = c("hap11","hap12","hap15","hap16","hap17","hap18")) +
  scale_fill_manual(breaks = c("hap11","hap12","hap15","hap16","hap17","hap18"),
                    values = c("#DDB8FF","grey30","#CFFFB8","#71E2ED","#AEB431","#6BA9B4"), name = "Color key", drop = TRUE, guide = "none") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_text(size = 10, vjust = 1, hjust = 0.5, margin = margin(t = 1, b = 0, l = 0, r = 0)),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.1), units = "in"),
        plot.margin = unit(c(0,0,1.1,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(0,-0.01),
        legend.background = element_blank(),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1, fill = NA, linetype = 0)))

ggsave("Figure6A.tiff", plot, device = "tiff", width = 7.5, height = 2.8, units = "in", dpi = 600)


## Figure 6B

merged_intervals <- read.delim("all_intervals_fig6B.txt", sep = "\t", header = FALSE)

colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
keepme <- final_joined

chrlen_aus <- subset(chrlen_melt, !(HAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24")))

aus_intervals <- subset(merged_intervals, !(HAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24")))
aus_intervals <- subset(aus_intervals, !(REFHAP %in% c("hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap1","hap2","hap21","hap22","hap24",
                                                       "hap6","hap32","hap25","hap31","hap7","hap8","hap29","hap30","hap27","hap28")))

final_aus <- full_join(aus_intervals, chrlen_aus, by = c("HAP","CHROM"))
final_aus$HAP <- factor(final_aus$HAP)
just_aus <- final_aus

just_aus$HAP <- factor(just_aus$HAP, levels = c("hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32"))
just_aus$CHROM <- factor(just_aus$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
chrlen_aus$HAP <- factor(chrlen_aus$HAP, levels = c("hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32"))
just_aus$REFHAP <- factor(just_aus$REFHAP)

plot2 <- ggplot(data = just_aus, aes(y = HAP)) +
  geom_point(aes(x = len+15000), y = -0.5, color = "white") +
  geom_point(aes(x = len+15000), y = 1.5, color = "white") +
  geom_point(x = -15000, y = -0.5, color = "white") +
  geom_point(x = -15000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = chrlen_aus, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.1, 'npc')) +
  geom_segment(data = na.omit(just_aus), aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 5.7, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= chrlen_aus, aes(xmax = len+15000), xmin = -14999, ymin = -0.01, ymax = 1.01, r = unit(0.1, 'npc'), color = "black", fill = NA, size = 0.3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free_x", space = "free_x", switch = "y", margins = FALSE) +
  scale_x_continuous(expand = c(0.03,0.03)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8",
                                "hap27","hap28","hap29","hap30","hap31","hap32"),
                     values = c("#F45587","#FFA900","#0027D9","#6CAD00","black","grey","white","white","white",
                                "white","white","white","white","white","white"), name = "Color key", drop = TRUE, limits = c("hap3","hap4","hap5","hap26","hap6","hap23")) +
  scale_fill_manual(breaks = c("hap3","hap4","hap5","hap26","hap6","hap23","hap25","hap7","hap8",
                               "hap27","hap28","hap29","hap30","hap31","hap32"),
                    values = c("#F45587","#FFA900","#0027D9","#6CAD00","black","grey","white","white","white",
                               "white","white","white","white","white","white"), name = "Color key", drop = TRUE, guide = "none") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_text(size = 10, vjust = 1, hjust = 0.5, margin = margin(t = 1, b = 0, l = 0, r = 0)),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.1), units = "in"),
        plot.margin = unit(c(0,0,1.1,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(0.7,-0.01),
        legend.background = element_blank(),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1, fill = NA, linetype = 0)))

ggsave("Figure6B.tiff", plot2, device = "tiff", width = 7.5, height = 6, units = "in", dpi = 600)


## Figure S10
merged_intervals <- read.delim("FigS10_intervals.txt", sep = "\t", header = FALSE)

colnames(merged_intervals) <- c("CHROM","start","end","HAP","dens","REFHAP")
merged_intervals <- merged_intervals[,c(4,1,2,3,6)]

final_joined <- full_join(merged_intervals, chrlen_melt, by = c("HAP","CHROM"))
final_joined$HAP <- factor(final_joined$HAP)
figureS10_revised <- final_joined


figS10_chr_sub <- subset(chrlen_melt, HAP %in% c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"))

figureS10_revised$HAP <- factor(figureS10_revised$HAP, levels =c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"))
figureS10_revised$REFHAP <- factor(figureS10_revised$REFHAP, levels =c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"))
figureS10_revised$CHROM <- factor(figureS10_revised$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
figS10_chr_sub$HAP <- factor(figS10_chr_sub$HAP, levels = c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14","hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"))


figureS10_revised <- figureS10_revised[!is.na(figureS10_revised$REFHAP),]


## Supplemental Figure 10

plot3 <- ggplot(data = figureS10_revised, aes(y = HAP)) +
  geom_point(aes(x = len+50000), y = -0.5, color = "white") +
  geom_point(aes(x = len+50000), y = 1.5, color = "white") +
  geom_point(x = -50000, y = -0.5, color = "white") +
  geom_point(x = -50000, y = 1.5, color = "white") +
  ggchicklet:::geom_rrect(data = figS10_chr_sub, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.1, 'npc')) +
  geom_segment(aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 5.6, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= figS10_chr_sub, aes(xmax = len+50000), xmin = -49999, ymin = -0.01, ymax = 1.01, r = unit(0.1, 'npc'), color = "black", fill = NA, size = 0.3) +
  facet_grid(HAP ~ CHROM, labeller = labeller(CHROM = chr_labels), scales = "free_x", space = "free_x", switch = "y", margins = FALSE) +
  scale_x_continuous(expand = c(0.03,0.01)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14",
                                "hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"),
                     values = c("#FF00F7","#00E095","black", "#ffc1f8", "#FFCD97", "#b7a2ff", "grey30", "#da8a16", "#006023",
                                "#CFFFB8","#71E2ED","#AEB431","#6BA9B4","grey70","#fff63b","white","white","white"), 
                     name = "Color key", drop = TRUE, limits = c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14",
                                                                 "hap15","hap16","hap17","hap18","hap19","hap20")) +
  scale_fill_manual(breaks = c("hap1","hap2","hap9","hap10","hap11","hap12","hap13","hap14",
                               "hap15","hap16","hap17","hap18","hap19","hap20","hap5","hap6","hap23"),
                    values = c("#FF00F7","#00E095","black", "#ffc1f8", "#FFCD97", "#b7a2ff", "grey30","#da8a16", "#006023",
                               "#CFFFB8","#71E2ED","#AEB431","#6BA9B4","grey70","#fff63b","white","white","white"), 
                    name = "Color key", drop = TRUE, guide = "none") +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = NULL) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
        strip.text.x.top = element_text(size = 10, vjust = 1, hjust = 0.5, margin = margin(t = 1, b = 0, l = 0, r = 0)),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.1), units = "in"),
        plot.margin = unit(c(0,0,2,0), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = c(-0.01,-0.03),
        legend.background = element_blank(),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 1, fill = NA, linetype = 0)))


ggsave("FigureS10.tiff", plot3, device = "tiff", width = 7.5, height = 7, units = "in", dpi = 600)
