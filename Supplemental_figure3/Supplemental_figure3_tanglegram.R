library(circlize)
library(dendextend)
library(ape)
library(phytools)
library(ggplot2)

scoring_data <- read.csv("Supplemental_figure3_pathotypes.csv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

key <- setNames(c(0, 0, 0, 0, 0, 0, 0, 2, 2, 2), c("0", "0;", ";", ";C", "1;", "1", "2", "3", "3+", "4"))

scoring_data_numeric <- scoring_data

scoring_data_numeric[,2:ncol(scoring_data_numeric)] <- key[unlist(scoring_data_numeric[,2:ncol(scoring_data_numeric)])]

scoring_matrix <- as.matrix(scoring_data_numeric[,c(2:ncol(scoring_data_numeric)), drop = FALSE])

row_names <- scoring_data$"Differential"

rownames(scoring_matrix) = row_names

pathotype_clust <- as.dendrogram(hclust(dist(scoring_matrix), method = "complete"))

phylo_tree <- read.tree(file = "Figure_S3_pruned_tree.nwk")
phylo_tree$node.label <- NULL
chrono <- chronos(phylo_tree, 1)
plotTree.singletons(chrono)
chrono_collapse <- collapse.singles(chrono)
plot.phylo(chrono, show.node.label = TRUE, use.edge.length = TRUE)
is.rooted(chrono_collapse)
is.ultrametric(chrono_collapse)
is.binary(chrono_collapse)
test <- as.hclust(chrono_collapse)
phylo_dendro <- as.dendrogram(test)
idea <- remove_branches_edgePar(phylo_dendro)

repeat_tests <- dendextend::untangle_step_rotate_1side(pathotype_clust, phylo_dendro, direction = "backward", leaves_matching_method = "labels")

dendextend::tanglegram(repeat_tests, edge.lwd = 1, main_left = "Pathotype clustering", main_right = "Genetic phylogeny",
                       lwd = 2, columns_width = c(3,1,3), margin_inner = 4.7, axes = FALSE, cex_main = 1.4, lab.cex = 1, margin_top = 1, margin_bottom = 1, margin_outer = 1, rank_branches = TRUE)

dev.copy(pdf, "Supplemental_figure3_tanglegram.pdf", width = 8, height = 9)
dev.off()

cor_bakers_gamma(repeat_tests)
