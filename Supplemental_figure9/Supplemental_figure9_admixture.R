library(ggplot2)
library(pophelper)

Kvals <- c(2,3,4,5,6,7,8,9,10,11)
CV <- c(0.20988,0.16099,0.14692,0.13321,0.11682,0.10669,0.09827,0.11117,0.10077,0.12824)

CV_error <- data.frame(Kvals, CV)
ggplot(CV_error) +
  geom_line(aes(x = Kvals, y = CV))

pop2 <- list.files("path/", pattern = "*.Q", full.names = TRUE)

pop2_dat <- readQ(pop2)

pop2_dat <- lapply(pop2_dat, function(x) { rownames(x) <- c(2,1,116,34,36,113,117,3,4,61,63,52,68,66,51,67,17,32,114,112,35,18,28,118,33,19,71,65,64,115,46,12,10,24,104,11,81,91,31,108,30,22,14,15,6,69,70,88,16,87,9,29,27,94,79,105,90,57,75,85,62,86,45,58,59,84,41,49,47,38,43,44,48,54,53,56,103,98,42,80,73,50,72,60,37,74,76,77,82,102,55,83,101,111,97,106,96,95,20,21,93,107,100,110,109,99,5,7,92,25,26,8,13,40,39,23,89,78); x})
pop2_dat <- lapply(pop2_dat, function(x) {x[order(as.numeric(row.names(x))),]})
pop2_dat <- sortQ(pop2_dat, by = "k")
colnames(pop2_dat[[2]]) <- c("Cluster1","Cluster3","Cluster2")
colnames(pop2_dat[[3]]) <- c("Cluster4", "Cluster2","Cluster3", "Cluster1")
colnames(pop2_dat[[4]]) <- c("Cluster2","Cluster3","Cluster1","Cluster5","Cluster4")
colnames(pop2_dat[[5]]) <- c("Cluster3","Cluster4","Cluster5","Cluster1","Cluster6","Cluster2")
colnames(pop2_dat[[6]]) <- c("Cluster1","Cluster5","Cluster6","Cluster2","Cluster3","Cluster7","Cluster4")
colnames(pop2_dat[[7]]) <- c("Cluster5","Cluster4","Cluster1","Cluster3","Cluster2","Cluster6","Cluster8","Cluster7")
colnames(pop2_dat[[8]]) <- c("Cluster8","Cluster2","Cluster3","Cluster9","Cluster6","Cluster5","Cluster1","Cluster7","Cluster4")
colnames(pop2_dat[[9]]) <- c("Cluster9","Cluster2","Cluster4","Cluster5","Cluster8","Cluster3","Cluster7","Cluster6","Cluster1","Cluster10")
colnames(pop2_dat[[10]]) <- c("Cluster8","Cluster3","Cluster5","Cluster1","Cluster7","Cluster11","Cluster4","Cluster2","Cluster6","Cluster9","Cluster10")
plotQ(pop2_dat, imgoutput = "join", clustercol = c("black","#FF6E6E","#A9171D","#FFB06E","#7BD8DD","#C848BA","#A6E179","hotpink","#7A571F","#D692FC","#9792FC"), sharedindlab = TRUE, exportplot = FALSE, showsubtitle = FALSE, returnplot = TRUE, showlegend = FALSE, legendtextsize = 10)
