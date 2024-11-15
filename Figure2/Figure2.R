library(ggplot2)
library(ozmaps)
library(scatterpie)
library(dplyr)
library(patchwork)


## Figure 2b
piedat <- read.delim("Figure2_map_piecharts_wide.txt", sep = "\t")
piedat$state <- as.factor(piedat$state)

map <- ozmap("states")


isolate_metadata <- read.csv("Figure2_aus_sequenced_metatadata.csv")
isolate_metadata$Lineage <- as.factor(isolate_metadata$Lineage)

ACT <- abs_lga %>% dplyr::filter(grepl("ACT", NAME))

ACT_box_window <- ggplot(data = ACT) + 
  geom_sf(data = ACT, fill = "white",linewidth = 0.7, color = "grey70") +
  geom_point(data = isolate_metadata[isolate_metadata$State == "ACT",], aes(x = long, y = lat, color = Lineage), size = 1, alpha = 0.6) +
  geom_rect(aes(xmin = 149.45, xmax = 148.7, ymin = -35.95, ymax = -35.1), fill = NA, color = "black", linewidth = 0.2) +
  geom_text(aes(x = 149.27, y = -35.6, label = "ACT\nn = 9"), size = 2) +
  scale_color_manual(values = c("black", "grey30", "red")) +
  theme_void() +
  theme(legend.position = "none")


#Plotting just the map layer, the geom_scatterpie pie charts get warped using coord_sf, but we need coord_sf for this map method...
map_plot <- ggplot(data = map) + 
  geom_sf(data = map, fill = "white", linewidth = 0.7, color = "grey70") +
  #geom_scatterpie(aes(x=lat, y=long, group = state, r = log10(total)*2), 
   #               data = piedat, cols = colnames(piedat[,c(4:(ncol(piedat)-1))]), legend_name = "Lineage", color = NA) +
  geom_point(data = isolate_metadata, aes(x = long, y = lat, color = Lineage), size = 1, alpha = 0.6) +
  geom_text(data = piedat, aes(x = (lat + (log10(total)*2) + 1.6), y = (long + 0.5), label = paste(state,"\nn = ",total, sep ="")), size = 2) +
  scale_color_manual(values = c("black", "grey30", "#23A6FE", "#0D3E5F","#95D4FF", "#558035", "#8EC763", "#B6EA90", "#FDD9A0", "#FF9C00","#895501", "#C75A52","#FE9F98","#D297F1", "#A406F8", "#F939F3","#A50076","red"),labels = c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12","L13","L14","L15","L16","L17","L18")) + 
  scale_fill_manual(values = c("black", "grey30", "#23A6FE", "#0D3E5F", "#95D4FF", "#558035", "#8EC763", "#B6EA90", "#FDD9A0", "#FF9C00", "#895501", "#C75A52", "#FE9F98","#D297F1", "#A406F8","#F939F3","#A50076","red"), guide = "none") +
  coord_sf(xlim = c(112,157), ylim = c(-45,-10)) +
  theme_void() +
  theme(legend.position = c(0.1,0.82),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.width = unit(0.1, units = "in"),
        legend.key.height = unit(0.1, units = "in"),
        legend.direction = "vertical") +
  guides(color = guide_legend(nrow = 9)) +
  geom_segment(aes(x = 148.9, xend = 151.46, y = -35.9, yend = -41.45), color = 'black', linewidth = 0.4) +
  geom_segment(aes(x = 149, xend = 151.5, y = -35.1, yend = -34.4), color = 'black', linewidth = 0.4) +
  inset_element(ACT_box_window, left = 0.8, bottom = 0.13, right = 1, top = 0.33)

# Plotting just the pie charts with coord_fixed so they are perfect circles
map_plot2 <- ggplot(data = map) + 
  geom_scatterpie(aes(x=lat, y=long, group = state, r = log10(total)*2), 
                  data = piedat, cols = colnames(piedat[,c(4:(ncol(piedat)-1))]), legend_name = "Lineage", color = NA) +
  #geom_text(data = piedat, aes(x = (lat + (log10(total)*2) + 1.6), y = (long + 0.5), label = paste(state,"\nn = ",total, sep ="")), size = 2) +
  scale_color_manual(values = c("black", "grey30", "#23A6FE", "#0D3E5F","#95D4FF", "#558035", "#8EC763", "#B6EA90", "#FDD9A0", "#FF9C00","#895501", "#C75A52","#FE9F98","#D297F1", "#A406F8", "#F939F3","#A50076","red"),labels = c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12","L13","L14","L15","L16","L17","L18")) + 
  scale_fill_manual(values = c("black", "grey30", "#23A6FE", "#0D3E5F", "#95D4FF", "#558035", "#8EC763", "#B6EA90", "#FDD9A0", "#FF9C00", "#895501", "#C75A52", "#FE9F98","#D297F1", "#A406F8","#F939F3","#A50076","red"), guide = "none") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = c(0.1,0.82),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.width = unit(0.1, units = "in"),
        legend.key.height = unit(0.1, units = "in"),
        legend.direction = "vertical") +
  guides(color = guide_legend(nrow = 9))


# Save the plots separately and combine in powerpoint.
ggsave(filename = "Figure2b_part1.tiff", plot = map_plot, device = "tiff", width = 6, height = 5, dpi = 600, units = "in")
ggsave(filename = "Figure2b_part2.tiff", plot = map_plot2, device = "tiff", width = 6, height = 5, dpi = 600, units = "in")
