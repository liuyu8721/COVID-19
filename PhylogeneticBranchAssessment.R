rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/all.meta.RData")

###################################################################################################
### phylogenetic tree
###################################################################################################

tree <- read.tree("./data/Worldwide/Worldwide_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("223"))
# pdf("worldwide_tree_node.pdf", width = 7, height = 14)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
# dev.off()

library(RColorBrewer)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
tree2 <- groupClade(tree, c("10", "221", "188", "203", "190", "103", "205", 
                            "106", "333", "1185", "334", "1896", "3333", "5067"))
ggtree(tree2, aes(color=group)) +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c("grey",
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10],
                              color.list[4], #G:333
                              color.list[2], #GK: 1185
                              color.list[5], #GH: 334
                              color.list[3], #GV: 1896
                              color.list[6], #GR: 3333
                              color.list[1], #GRY: 5067
                              rep("grey99", 10))) +
  theme(legend.position = "none")

###################################################################################################
### Mark the detected samples 
###################################################################################################

all.meta$Virus.name <- sub(pattern = " ", replacement = "_", x = all.meta$Virus.name)
tree.detail <- subset(all.meta, Virus.name %in% tree$tip.label)

#---Subgroup
subgroup03_P1 <- read.delim2("./data/Worldwide/ValidationData/subgroup03_P1.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03_P1 <- subset(tree.detail, Accession.ID %in% subgroup03_P1$gisaid_epi_isl)
subgroup05_B12 <- read.delim2("./data/Worldwide/ValidationData/subgroup05_B12.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05_B12 <- subset(tree.detail, Accession.ID %in% subgroup05_B12$gisaid_epi_isl)
subgroup06_B11284 <- read.delim2("./data/Worldwide/ValidationData/subgroup06_B11284.txt", 
                              header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup06_B11284 <- subset(tree.detail, Accession.ID %in% subgroup06_B11284$gisaid_epi_isl)
subgroup07_D2 <- read.delim2("./data/Worldwide/ValidationData/subgroup07_D2.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07_D2 <- subset(tree.detail, Accession.ID %in% subgroup07_D2$gisaid_epi_isl)
subgroup08_GKplus <- read.delim2("./data/Worldwide/ValidationData/subgroup08_GKplus.txt", 
                                 header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08_GKplus <- subset(tree.detail, Accession.ID %in% subgroup08_GKplus$gisaid_epi_isl)
subgroup10_GRYplus <- read.delim2("./data/Worldwide/ValidationData/subgroup10_GRYplus.txt", 
                                  header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup10_GRYplus <- subset(tree.detail, Accession.ID %in% subgroup10_GRYplus$gisaid_epi_isl)
subgroup14_UC <- read.delim2("./data/Worldwide/ValidationData/subgroup14_UC.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14_UC <- subset(tree.detail, Accession.ID %in% subgroup14_UC$gisaid_epi_isl)

# pdf("worldwide_tree_highlight.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c("grey",
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10],
                              color.list[4], #G:333
                              color.list[2], #GK: 1185
                              color.list[5], #GH: 334
                              color.list[3], #GV: 1896
                              color.list[6], #GR: 3333
                              color.list[1], #GRY: 5067
                              rep("grey99", 10)))  +
  
  geom_tippoint(aes(subset=((label %in% subgroup03_P1$Virus.name)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[1], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup05_B12$Virus.name)==TRUE)),
                size = 4, fill = "black", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup06_B11284$Virus.name)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[2], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup07_D2$Virus.name)==TRUE)),
                size = 4, fill = "blue", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08_GKplus$Virus.name)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[3], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup10_GRYplus$Virus.name)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[4], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup14_UC$Virus.name)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[5], color= 'grey80', shape = 21) +
  
  theme(legend.position = "none")
# dev.off()
