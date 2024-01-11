library(ape)

stetree <- read.tree("STE3.2_tree_pangenome1_final.nwk")

rooted <- root(stetree, "Pgt_STE3.2.1")

write.tree(rooted, file = "STE3.2_tree_pangenome1_final_rooted.nwk")
plot.phylo(rooted)



bWtree <- read.tree("bWtree_unrooted.nwk")

rooted2 <- root(bWtree, "Pgt_bW1-HD1")

write.tree(rooted2, file = "bW_tree_pgt_rooted.nwk")
plot.phylo(rooted2)


bEtree <- read.tree("bEtree_unrooted.nwk")

rooted3 <- root(bEtree, "Pgt_bE1-HD2")

write.tree(rooted3, file = "bE_tree_pgt_rooted.nwk")
plot.phylo(rooted3)
