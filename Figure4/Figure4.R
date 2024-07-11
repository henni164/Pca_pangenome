library(ggplot2)
library(maps)
library(dplyr)

options(scipen = 999)


## Figure 4A

data2 <- read.delim("20QLD86_haplotypes_mash.txt", header = TRUE, sep = "\t")

data2$haplotype <- factor(data2$haplotype, levels = c("hap5","hap6","hap23","hap24"))

aus_plot <- ggplot(data = data2) +
  #plot points
  geom_point(data = data2, aes(x = contain, y = ident, color = interest), size = 1.6) +
  geom_point(data = subset(data2, interest == "1"), aes(x = contain, y = ident, color = interest), size = 1.6) +
  geom_point(data = subset(data2, interest == "2"), aes(x = contain, y = ident, color = interest), size = 1.6) +
  geom_point(data = subset(data2, interest == "3"), aes(x = contain, y = ident, color = interest), size = 1.6) +
  geom_point(data = subset(data2, interest == "4"), aes(x = contain, y = ident, color = interest), size = 1.6) ++
  facet_wrap(~ haplotype, scales = "fixed", ncol = 4, nrow = 1) +
  #set plot labs
  labs(x = "% shared k-mers", y = "% k-mer identity") +
  scale_color_manual(breaks = c("N","1","2","3","4"), values = c("grey90","black","grey30", "#0D3E5F", "#D253FF")) +
  #theme stuff. make plot margins 0 to get closer fit between scatterplot and next plots
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 8, color = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0,0.1,0,0), units = "lines"),
        axis.text.x = element_text(size = 8, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
  ) 

ggsave("Figure_4A.tiff", aus_plot, device = "tiff", width = 7.5, height = 2, units = "in")


# Figure 4F

dat <- map_data("world") %>%
  select(lon = long, lat, group, region)

mapplot <- ggplot(data= dat, aes(x = lon, y = lat, group = group)) +
  geom_polygon(data = dat, fill = "grey90", color = "grey70") +
  geom_polygon(data = dat[dat$region == "Taiwan",], fill = "grey10", color = "black", alpha = 0.5) +
  geom_polygon(data = dat[dat$region == "Australia",], fill = "grey50", color = "black", alpha = 0.5) +
  theme_void() +
  coord_cartesian(xlim = c(90,180), ylim = c(-50,40))

ggsave("Figure4F_map.tiff", mapplot, device = "tiff", width = 6, height = 6, units = "in", dpi = 600)
