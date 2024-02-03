library(ggplot2)

data <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/GitHub/Pca_pangenome/Figure4/90TX52_haplotypes_mash.txt", header = TRUE, sep = "\t")

data$haplotype <- as.factor(data$haplotype)

tx_plot <- ggplot(data = data) +
  #plot points
  geom_point(data = data, aes(x = ident, y = contain, color = interest), size = 1.5) +
  #make margins around plot 0
  scale_x_continuous(expand = c(0,0), breaks = c(99.98, 100)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(NA,100), breaks = c(99.65, 100)) +
  facet_wrap(~ haplotype, scales = "free_x") +
  #set plot labs
  labs(x = "% k-mer identity", y = "% shared k-mers") +
  scale_color_manual(breaks = c("N","1","2"), values = c("white","#EC3F3F","#3432DE")) +
  coord_cartesian(xlim = c(99.98, 100.001), ylim = c(99.65, 100.01)) +
  #theme stuff. make plot margins 0 to get closer fit between scatterplot and next plots
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 5, color = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1, color = "black"),
  ) 


ggsave("C:/Users/HEN294/OneDrive - CSIRO/Documents/GitHub/Pca_pangenome/Figure4/tx_plot.tiff", tx_plot, device = "tiff", width = 3, height = 2, units = "in")


data2 <- read.delim("C:/Users/HEN294/OneDrive - CSIRO/Documents/GitHub/Pca_pangenome/Figure4/20QLD86_haplotypes_mash.txt", header = TRUE, sep = "\t")

data2$haplotype <- factor(data2$haplotype, levels = c("hap5","hap6","hap23","hap24"))
sub <- subset(data2, interest != "3")

aus_plot <- ggplot(data = sub) +
  #plot points
  geom_point(data = sub, aes(x = ident, y = contain, color = interest), size = 1.6) +
  #make margins around plot 0
  scale_x_continuous(expand = c(0,0), breaks = c(99.98, 100)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(NA,100), breaks = c(99.65, 100)) +
  facet_wrap(~ haplotype, scales = "fixed", ncol = 2, nrow = 2) +
  #set plot labs
  labs(x = "% k-mer identity", y = "% shared k-mers") +
  scale_color_manual(breaks = c("N","1","2","4"), values = c("white","black","grey50", "#D253FF")) +
  coord_cartesian(xlim = c(99.98, 100.001), ylim = c(99.65, 100.001)) +
  #theme stuff. make plot margins 0 to get closer fit between scatterplot and next plots
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 5, color = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1, color = "black"),
  ) 

ggsave("C:/Users/HEN294/OneDrive - CSIRO/Documents/GitHub/Pca_pangenome/Figure4/aus_plot.tiff", aus_plot, device = "tiff", width = 3, height = 3.7, units = "in")
