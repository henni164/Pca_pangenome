library(reshape2)
library(purrr)
library(ggchicklet)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(egg)

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
  geom_segment(aes(x = tstart, xend = tend, color = hap), y = 0.5, yend = 0.5, linewidth = 6.8, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= recombinant_paf2, aes(xmax = tlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_grid(hap ~ qname, scales = "fixed", space = "free_x", switch = "y", margins = FALSE, drop = FALSE) +
  scale_x_continuous(expand = c(0.07,0.07)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4"),
                     values = c("#F45587","#FFA900"), 
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
  geom_segment(aes(x = qstart, xend = qend, color = hap), y = 0.5, yend = 0.5, linewidth = 6.8, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= recombinant_paf2, aes(xmax = qlength+30000), xmin = -29999, ymin = -0.05, ymax = 1.05, r = unit(0.05, 'npc'), color = "black", fill = NA, size = 0.35) +
  facet_grid( ~ tname, scales = "fixed", space = "free_x", switch = "y", margins = FALSE, drop = FALSE) +
  scale_x_continuous(expand = c(0.07,0.07)) +
  scale_y_continuous(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4"),
                     values = c("#F45587","#FFA900"), 
                     name = "Color key", drop = FALSE, limits = c("hap3","hap4")) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  labs(x = NULL, y = NULL) +
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

ggsave("Supplemental_figure9e.tiff", recombinant_plot, device = "tiff", width = 7.5, height = 1.5, units = "in", dpi = 600)


## pangenome graph approach
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

haps_chrs <- read.delim("Figure3_chromosome_lengths.txt", sep = "\t", header = TRUE)
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


keepme$HAP <- factor(keepme$HAP, levels = c("hap3","hap4","hap25","hap5","hap26","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap10","hap11","hap12","hap13","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap23","hap9"))
keepme$CHROM <- factor(keepme$CHROM, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))

chrlen_melt <- subset(chrlen_melt, !(HAP %in% c("hap21","hap22","hap24")))
chrlen_melt$HAP <- factor(chrlen_melt$HAP, levels = c("hap3","hap4","hap25","hap5","hap26","hap7","hap8","hap27","hap28","hap29","hap30","hap31","hap32","hap14","hap1","hap2","hap10","hap11","hap12","hap13","hap15","hap16","hap17","hap18","hap19","hap20","hap6","hap23","hap9"))


## Figure 5

## Figure 5
## panel a subset

fig5a_chr_sub <- subset(chrlen_melt, HAP %in% c("hap25"))

fig5a_sub <- subset(keepme, HAP %in% c("hap25"))
fig5a_sub <- subset(fig5a_sub, REFHAP %in% c("hap3","hap4"))

levels(fig5a_chr_sub$CHROM) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
levels(fig5a_sub$CHROM) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")

hap25_bar <- read.delim("C:/Users/HEN294/Downloads/hap25_intervals.txt", header = FALSE)

hap25_subset <- subset(hap25_bar, V4 %in% c("hap3","hap4","hap5","hap6","hap7","hap8","hap21","hap22","hap23","hap24","hap26","hap27","hap28","hap29","hap30","hap31","hap32"))
#hap25_subset <- subset(hap25_subset, V6 %in% c("+"))

hap25_subset$V1 <- factor(hap25_subset$V1, levels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
levels(hap25_subset$V1) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
hap25_subset$V4 <- factor(hap25_subset$V4, levels = c("hap3","hap4","hap5","hap6","hap7","hap8","hap21","hap22","hap23","hap24","hap26","hap27","hap28","hap29","hap30","hap31","hap32"))

hap25_forplot <- subset(hap25_subset, V4 %in% c("hap3","hap4"))


plot <- ggplot(data = fig5a_sub, aes(y = HAP)) +
  ggchicklet:::geom_rrect(data = fig5a_chr_sub, aes(xmax = len, fill = HAP), xmin = 1, ymin = 0, ymax = 1, r = unit(0.1, 'npc')) +
  geom_segment(data = na.omit(fig5a_sub), aes(x = start, xend = end, color = REFHAP, group = HAP), y = 0.5, yend = 0.5, linewidth = 11.2, key_glyph = "pointrange") +
  ggchicklet:::geom_rrect(data= fig5a_chr_sub, aes(xmax = len+15000), xmin = -14999, ymin = -0.01, ymax = 1.01, r = unit(0.1, 'npc'), color = "black", fill = NA, size = 0.3) +
  facet_wrap(~ CHROM, ncol = 18, nrow = 1, scales = "fixed") +
  scale_x_continuous(expand = c(0.07,0.07)) +
  scale_y_discrete(limits=rev, expand = c(0.01,0.01)) +
  scale_color_manual(breaks = c("hap3","hap4"),
                     values = c("#F45587","#FFA900"), name = "Color key", drop = TRUE, limits = c("hap3","hap4")) +
  scale_fill_manual(breaks = c("hap25"),
                    values = c("white"), name = "Color key", drop = TRUE, guide = "none") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(-0.6,1.6)) + 
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_blank(),
        strip.text.x.top = element_blank(),
        strip.background.y = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(c(0.1), units = "in"),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.spacing = unit(0, units = "lines"),
        legend.position = "none",
        legend.background = element_blank(),
        legend.justification = "left",
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1, fill = NA, linetype = 0)))

hist <- ggplot(data = hap25_forplot) +
  geom_rect(data = subset(hap25_forplot, V4 == "hap3"), aes(xmin = V2, xmax = V3, ymax = V5), fill = "#F45587", ymin = 0) +
  geom_rect(data = subset(hap25_forplot, V4 == "hap4"), aes(xmin = V2, xmax = V3, ymax = V5), fill = "#FFA900", ymin = 0) +
  facet_wrap(~ V1, drop = TRUE, scales = "fixed", ncol = 18, nrow = 1) +
  scale_y_continuous(breaks = c(0,1000,2000), position = "left") +
  scale_x_continuous(expand = c(0.07, 0.07)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background.y = element_blank(),
        strip.background.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.01,0.01,0.01,0.01), units = "lines"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey30", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "grey70", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0, units = "lines"))

recombinant_plot2 <- ggpubr::ggarrange(hist, plot, nrow = 2, ncol = 1, heights = c(1,1), align = "v")

ggsave("Supplemental_figure9f.tiff", recombinant_plot2, device = "tiff", width = 6.6, height = 1.6, units = "in", dpi = 600)
