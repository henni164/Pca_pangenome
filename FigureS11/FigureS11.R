library(dplyr)
library(ggplot2)
library(stringr)
library(micropan)
library(tidyr)
library(reshape2)
library(ggpubr)

# FigureS11A

#orthogroup_counts_all <- read.delim("Orthogroups_combined_haplotypes.txt", header = TRUE)
#ortho_t <- t(orthogroup_counts_all[,2:ncol(orthogroup_counts_all)])
#colnames(ortho_t) <- paste0("Cluster_",seq(1:ncol(ortho_t)))
#rownames(ortho_t) <- paste0("GID",seq(1:nrow(ortho_t)))
#alpha_estimate <- heaps(ortho_t, n.perm = 1000)
#rarefaction_plot <- rarefaction(ortho_t, n.perm = 1000)
#no_zero <- subset(rarefaction_plot, !(Genome %in% c("0")))
#small_ortho <- no_zero %>% gather(key = "Permutation", value = "Clusters", -Genome)
#write.table(small_ortho, "orthogroups_parsed.txt", quote = FALSE, sep = "\t")

alpha_estimate_ortho <- 1.182893
ortho_small <- read.delim("orthogroups_parsed.txt", header = TRUE, row.names = 1)
fitheaps_ortho <- lm(log10(Clusters) ~ log10(Genome), data = ortho_small)
fitlog_ortho <- lm(Clusters ~ log10(Genome), data = ortho_small)

orthogroup_pangenome_plot <- ggplot(data = ortho_small) +
  geom_point(aes(x = Genome, y = Clusters, group = Permutation)) +
  geom_function(fun = function(x) 10^(summary(fitheaps_ortho)$coef[2]*log10(x) + summary(fitheaps_ortho)$coef[1]), color = "red") +
  geom_function(fun = function(x) summary(fitlog_ortho)$coef[2]*log10(x) + summary(fitlog_ortho)$coef[1], color = "blue") +
  annotate("text", x = 12.5, y = 30000, label = paste0("\u03b1", " = ", round(alpha_estimate_ortho, digits = 2)), size = 3) +
  scale_x_continuous(limits = c(0,(max(ortho_small$Genome) + 1)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,32000), expand = c(0,0)) +
  labs(x = NULL, y = "n orthogroups") +
  theme(plot.margin = margin(0,0.3,0.3,0, unit = "lines"),
        #axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, margin = margin(0,0,0.3,0, unit = "lines")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", hjust = 1, margin = margin(0,0,0,0.3, unit = "lines")),
        axis.title = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.5))


# FigureS11B

#all_blocks <- read.delim("all_haploblocks_complete.txt", header = FALSE)
#casted <- dcast(all_blocks, V1 ~ V2)
#casted[is.na(casted)] <- 0
#rownames(casted) <- casted$V1
#hapblocks <- t(casted[,2:ncol(casted)])
#colnames(hapblocks) <- paste0("Cluster_",seq(1:ncol(hapblocks)))
#rownames(hapblocks) <- paste0("GID",seq(1:nrow(hapblocks)))
#alpha_estimate <- heaps(hapblocks, n.perm = 1000)
#rarefaction_plot <- rarefaction(hapblocks, n.perm = 1000)
#no_zero <- subset(rarefaction_plot, !(Genome %in% c("0")))
#small_blocks <- no_zero %>% gather(key = "Permutation", value = "Clusters", -Genome)
#write.table(small_blocks, "hapblocks_parsed.txt", quote = FALSE, sep = "\t")

alpha_estimate_haplo <- 0.4611783
haplo_small <- read.delim("hapblocks_parsed.txt", header = TRUE, row.names = 1)
fitheaps_haplo <- lm(log10(Clusters) ~ log10(Genome), data = haplo_small)
fitlog_haplo <- lm(Clusters ~ log10(Genome), data = haplo_small)

haplotype_blocks_saturation_plot <- ggplot(data = haplo_small) +
  geom_point(aes(x = Genome, y = Clusters, group = Permutation)) +
  geom_function(fun = function(x) 10^(summary(fitheaps_haplo)$coef[2]*log10(x) + summary(fitheaps_haplo)$coef[1]), color = "red") +
  geom_function(fun = function(x) summary(fitlog_haplo)$coef[2]*log10(x) + summary(fitlog_haplo)$coef[1], color = "blue") +
  annotate("text", x = 12.5, y = 27000, label = paste0("\u03b1", " = ", round(alpha_estimate_haplo, digits = 2)), size = 3) +
  scale_x_continuous(limits = c(0,(max(haplo_small$Genome) + 1)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,33000), expand = c(0,0)) +
  labs(x = NULL, y = "n 100kb haplotype blocks") +
  theme(plot.margin = margin(0,0.3,0.3,0, unit = "lines"),
        #axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, margin = margin(0,0,0.3,0, unit = "lines")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", hjust = 1, margin = margin(0,0,0,0.3, unit = "lines")),
        axis.title = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.5))


#Figure S11C

#run separately on cluster because iterations take forever to run

#all_snps <- read.delim("hap1_snp_matrix_fixed.txt", header = TRUE)
#hapsnps <- t(all_snps[,2:ncol(all_snps)])
#colnames(hapsnps) <- paste0("Cluster_",seq(1:ncol(hapsnps)))
#rownames(hapsnps) <- paste0("GID",seq(1:nrow(hapsnps)))
#alpha_estimate <- heaps(hapsnps, n.perm = 1000)
#alpha_estimate
#rarefaction_plot <- rarefaction(hapsnps, n.perm = 1000)
#no_zero <- subset(rarefaction_plot, !(Genome %in% c("0")))
#small <- no_zero %>% gather(key = "Permutation", value = "Clusters", -Genome)
#write.table(small, file = "SNP_parsed.txt", sep = "\t", quote = FALSE)

alpha_estimate_snp <- 0 #number

snps_small <- read.delim("SNP_parsed.txt", header = TRUE)
fitheaps_snp <- lm(log10(Clusters) ~ log10(Genome), data = snps_small)
fitlog_snp <- lm(Clusters ~ log10(Genome), data = snps_small)

haplotype_snps_saturation_plot <- ggplot(data = snps_small) +
  geom_point(aes(x = Genome, y = Clusters, group = Permutation)) +
  geom_function(fun = function(x) 10^(summary(fitheaps_snp)$coef[2]*log10(x) + summary(fitheaps_snp)$coef[1]), color = "red") +
  geom_function(fun = function(x) summary(fitlog_snp)$coef[2]*log10(x) + summary(fitlog_snp)$coef[1], color = "blue") +
  #annotate("text", x = 12.5, y = 1000000, label = paste0("\u03b1", " = ", round(alpha_estimate_snp, digits = 2)), size = 3) +
  scale_x_continuous(limits = c(0,(max(ortho_small$Genome) + 1)), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1250000), expand = c(0,0)) +
  labs(x = "N haplotypes", y = "n SNPs") +
  theme(plot.margin = margin(0,0.3,0,0, unit = "lines"),
        axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, margin = margin(0,0,0.3,0, unit = "lines")),
        axis.text.y = element_text(size = 8, color = "black", hjust = 1, margin = margin(0,0,0,0.3, unit = "lines")),
        axis.title = element_text(size = 8, color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.5))


#Combine and save

combined_plot <- ggarrange(orthogroup_pangenome_plot, haplotype_blocks_saturation_plot, haplotype_snps_saturation_plot, ncol = 1, nrow = 3, align = "hv")

ggsave("FigureS11.tiff", combined_plot, device = "tiff", width = 5, height = 6, units = "in")


