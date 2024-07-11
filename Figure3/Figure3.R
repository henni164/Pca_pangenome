library(ggplot2)
library(sf)
library(giscoR)

# Figure 3A maps

## Oceania map
crs_string <- "+proj=ortho +lon_0=130 +lat_0=0"

ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho

# one of the visible ones red (don't really matter which one :)
world$fill_color <- ifelse(world$ISO3_CODE == "AUS", "interesting", "dull")

# now the action!
aus_globe <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "white", color = "black") + # background first
  geom_sf(aes(fill = fill_color), lwd = .2, color = "black") + # now land over the oceans
  scale_fill_manual(values = c("interesting" = "#FFE043",
                               "dull" = "grey90"),
                    guide = "none") +
  theme_void()

ggsave("aus_globe.tiff", aus_globe, device = "tiff", width = 1, height = 1, units = "in", dpi = 600)

## North America map
crs_string <- "+proj=ortho +lon_0=-90 +lat_0=16"

ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho

# one of the visible ones red (don't really matter which one :)
world$fill_color <- ifelse(world$ISO3_CODE == "USA", "interesting", "dull")

# now the action!
us_globe <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "white", color = "black") + # background first
  geom_sf(aes(fill = fill_color), lwd = .2, color = "black") + # now land over the oceans
  scale_fill_manual(values = c("interesting" = "#9087FF",
                               "dull" = "grey90"),
                    guide = "none") +
  theme_void()
us_globe

ggsave("us_globe.tiff", us_globe, device = "tiff", width = 1, height = 1, units = "in", dpi = 600)


# Figure 3B

hap1_containment <- read.delim("Fig3B_hap1_containment.txt", header = FALSE)

colnames(hap1_containment) <- c("Isolate","Lineage","hap1_ident","hap1_shared", "hap2_ident","hap2_shared")

hap1_containment_plot <- ggplot(hap1_containment) + 
  geom_point(aes(x = hap1_shared, y = hap1_ident), color = "black", alpha = 0.7, size = 1.6) + 
  geom_point(data = hap1_containment[hap1_containment$Isolate=="Pca203",], aes(x = hap1_shared, y = hap1_ident), 
             color = "red", alpha = 1, size = 1.6) +
  geom_hline(yintercept = 99.98, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 99.65, color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(y = expression(paste("% ",italic("k"),"-mer identity", sep = "")), 
       x = expression(paste("% shared ",italic("k"),"-mers", sep = ""))) +
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
        axis.line = element_line(color = "black")
  ) 

ggsave("Figure3B.tiff", hap1_containment_plot, device = "tiff", width = 3, height = 3, units = "in")

# Figure 3C

clones_containment <- read.delim("clones_containment.txt", header = FALSE)

clones_containment_plot <- ggplot(clones_containment) + 
  geom_point(aes(x = V3, y = V2), 
             color = "black", alpha = 0.3, size = 1.6) +
  geom_point(data = clones_containment[clones_containment$V1=="Pca203",], aes(x = V3, y = V2), 
             color = "red", size = 1.6) +
  geom_hline(yintercept = 99.98, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 99.65, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(99.5,100.001)) +
  scale_y_continuous(limits = c(99.97,100.001)) +
  labs(y = expression(paste("% ",italic("k"),"-mer identity", sep = "")), 
       x = expression(paste("% shared ",italic("k"),"-mers", sep = ""))) +
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
        axis.line = element_line(color = "black")
  ) 


ggsave("Figure3C.tiff", clones_containment_plot, device = "tiff", width = 3, height = 3, units = "in")
