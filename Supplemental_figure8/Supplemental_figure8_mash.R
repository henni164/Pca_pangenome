library(ggplot2)

data <- read.delim("90TX52_haplotypes_mash.txt", header = TRUE, sep = "\t")

data$haplotype <- as.factor(data$haplotype)

tx_plot <- ggplot(data = data) +
  #plot points
  geom_point(data = data, aes(x = contain, y = ident, color = interest), size = 1.5) +
  #make margins around plot 0
  #scale_x_continuous(expand = c(0,0), breaks = c(99.98, 100)) +
  #scale_y_continuous(expand = c(0.01,0.01), limits = c(NA,100), breaks = c(99.65, 100)) +
  facet_wrap(~ haplotype, scales = "free_x") +
  #set plot labs
  labs(x = "% shared k-mers", y = "% k-mer identity") +
  scale_color_manual(breaks = c("N","1","2"), values = c("grey90","#EC3F3F","#3432DE")) +
  #coord_cartesian(xlim = c(99, 100.001), ylim = c(99, 100.01)) +
  #theme stuff. make plot margins 0 to get closer fit between scatterplot and next plots
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8, color = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "lines"),
        axis.text.x = element_text(size = 8, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
  ) 


ggsave("Figure_S8a.tiff", tx_plot, device = "tiff", width = 7, height = 2, units = "in")