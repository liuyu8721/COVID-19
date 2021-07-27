rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")
library("treeio")
library(RColorBrewer)

setwd("D:/00SARS-Cov-2/ADS")
color.list <- rev(c(RColorBrewer::brewer.pal(8, "Set2"),  RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(9, "Pastel1")))
load("India_PhyTree.RData")

#########################################################################

tree <- read.tree("India_PhyTree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("90"))
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")

ggtree(tree, aes(color=1), mrsd = NULL, as.Date = FALSE, size = 0.2) + theme_tree2() +
  geom_tippoint(size = 4, fill = 'grey96', color = 'gray60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="03")$strain)==TRUE)), 
                size = 4, fill = color.list[3], color= 'white',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="02")$strain)==TRUE)), 
                size = 4, fill = color.list[2], color= 'grey88',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="05")$strain)==TRUE)), 
                size = 4, fill = color.list[6], color = 'grey90',shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="06")$strain)==TRUE)), 
                size = 4, fill = color.list[7], color = 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="08")$strain)==TRUE)), 
                size = 4, fill = color.list[10], color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="11")$strain)==TRUE)), 
                size = 4, fill = color.list[14], color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="09")$strain)==TRUE)), 
                size = 4, fill = color.list[12], color= 'grey60', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="10")$strain)==TRUE)), 
                size = 4, fill = color.list[13], color= 'black', shape=21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="12")$strain)==TRUE)), 
                size = 4, fill = color.list[16], color= 'grey60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="13")$strain)==TRUE)), 
                size = 4, fill = color.list[17], color= 'grey60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="14")$strain)==TRUE)), 
                size = 4, fill = color.list[19], color= 'grey60', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="16")$strain)==TRUE)), 
                size = 4, fill = color.list[21], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subset(bg.clade, group=="17")$strain)==TRUE)), 
                size = 4, fill = color.list[23], color= 'white', shape = 21) +
  theme(legend.position = "none")

