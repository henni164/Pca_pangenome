library(ggplot2)
library(ggpubr)

chr9_lengths <- read.delim("chr9_lengths_STE.txt", header = TRUE)

chr9_lengths$allele <- factor(chr9_lengths$allele, levels = c("STE3.2.2", "STE3.2.3"))
chr9_lengths <- chr9_lengths[,2:3]

chr9_length_fig <- ggboxplot(chr9_lengths, x = "allele", y = "Chr9") +
  stat_compare_means(method = "wilcox.test", method.args = list("two.sided"), comparisons = list(c("STE3.2.2", "STE3.2.3")), size = 3) +
  scale_y_continuous(labels = function(y)y/1000000) +
  labs(x = NULL, y = "Length of chromosome 9 (Mbp)") +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_line(color = "grey90"),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10, face = "italic"),
        axis.ticks.x = element_blank())

ggsave("Figure_S5c_chr9lengths.tiff", chr9_length_fig, device = "tiff", width=2.25, height=4, units = "in")

