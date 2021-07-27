rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")
library("treeio")

load("Worldwide_PhyTree.RData")

#########################################################################

tree <- read.tree("Worldwide_PhyTree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("1040"))
plot(tree, show.node.label = T, cex = 0.05, show.tip.label = F, edge.color = "grey90")

ggtree(tree, aes(color=1), mrsd = "2021-02-27", as.Date = TRUE, size = 0.2) + theme_tree2() +
  geom_tippoint(size = 4, fill = 'grey96', color = 'gray60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "17")$strain)==TRUE)), 
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[8], color= 'grey60',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "03")$strain)==TRUE)), 
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[1], color= 'grey60',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "05")$strain)==TRUE)), 
                size = 4, fill = "black", color = 'grey90',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "06")$strain)==TRUE)), 
                size = 4, fill = "yellow", color = 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "07")$strain)==TRUE)), 
                size = 4, fill = "blue", color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "09")$strain)==TRUE)), 
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[3], color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "08")$strain)==TRUE)), 
                size = 4, fill = adjustcolor(RColorBrewer::brewer.pal(9, "Pastel1")[2], 1), color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "11")$strain)==TRUE)), 
                size = 4, fill = adjustcolor(RColorBrewer::brewer.pal(9, "Pastel1")[4], 1), color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "13")$strain)==TRUE)), 
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[5], color= 'grey60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group == "15")$strain)==TRUE)), 
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[7], color= 'grey60', shape = 21) +
  theme(legend.position = "none")
