rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")
# library("treeio")
# suppressMessages(library(tidyverse))
Sys.setlocale("LC_TIME", "English")

setwd("D:/00SARS-Cov-2/ADS")
load("./RData/all.meta.update.RData")

###################################################################################################

tree <- read.tree("D:\\00SARS-Cov-2\\ADS\\20211114\\bg_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("203")) #203
pdf("bg_tree_node.pdf", width = 7, height = 14)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
dev.off()

all.meta.update$Virus.name <- sub(pattern = " ", replacement = "_", x = all.meta.update$Virus.name)
tree.detail <- subset(all.meta.update, Virus.name %in% tree$tip.label)

# l <- subset(tree.detail, Clade=="L", select = Virus.name)
# s <- subset(tree.detail, Clade=="S", select = Virus.name)
# v <- subset(tree.detail, Clade=="V", select = Virus.name)
# o <- subset(tree.detail, Clade=="O", select = Virus.name)
# g <- subset(tree.detail, Clade %in% c("G", "GR", "GV", "GH", "GK", "GRY"), select = Virus.name)
# gh <- subset(tree.detail, Clade=="GH", select = Virus.name)
# gr <- subset(tree.detail, Clade=="GR", select = Virus.name)
# gry <- subset(tree.detail, Clade=="GRY", select = Virus.name)
# gv <- subset(tree.detail, Clade=="GV", select = Virus.name)
# gk <- subset(tree.detail, Clade=="GK", select = Virus.name)
# b117 <- subset(tree.detail, Lineage=="B.1.1.7", select = Virus.name)
# b16171 <- subset(tree.detail, Lineage=="B.1.617.1", select = Virus.name)

load("D:/00SARS-Cov-2/ADS/RData/all.muta.update.list.RData")
tree.muta <- subset(all.muta.update.list, sample %in% tree.detail$Accession.ID)
core.muta <- read.csv("Worldwide_CodeMuta.csv", header = T)
for (t in c("S")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
s <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("V")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
v <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("G")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
g <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("GH")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
gh <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("GV")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
gv <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("GR")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
gr <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("GRY")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
gry <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)
for (t in c("GK")) {
  core <- subset(core.muta, RelatedVariant==t & CoreMutation=="YES")
  sg <- tree.muta
  for (m in core$NucleotideChange) {
    sg <- subset(sg, grepl(m, muta.list))
  }
}
gk <- subset(tree.detail, Accession.ID %in% sg$sample, select = Virus.name)

library(RColorBrewer)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
tree2 <- groupClade(tree, c("10", "205", "240", "196", "224", "48", #S
                            "111", #V
                            "333", #G
                            "2394",#GH
                            "367", #GV
                            "2825",#GR
                            "3896",#GRY
                            "1882" #GK
                            ))
pdf("bg_tree_skeleton.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  geom_tippoint(size = 4, fill = 'grey96', color = 'gray60', shape = 21) +
  # geom_nodelab(size = 1, color = "black") +
  scale_color_manual(values=c("grey",
                              color.list[9],  #S
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10],#V:111
                              color.list[4], #G:333
                              color.list[5], #GH:2394
                              color.list[3], #GV:367
                              color.list[6], #GR:2825
                              color.list[1], #GRY:3896
                              color.list[2], #GK:1882
                              rep("grey99", 10))) +
  theme(legend.position = "none")
dev.off()

pdf("bg_tree.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  geom_tippoint(size = 4, fill = 'grey96', color = 'gray60', shape = 21) +
  # geom_nodelab2(size = 1) +
  scale_color_manual(values=c("grey",
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10],
                              color.list[4], #G:333
                              color.list[5], #GH:2394
                              color.list[3], #GV:367
                              color.list[6], #GR:2825
                              color.list[1], #GRY:3896
                              color.list[2], #GK:1882
                              rep("grey99", 10))) +
  
  geom_tippoint(aes(subset=((label %in% s$Virus.name)==TRUE)),
                size = 4, fill = color.list[9], color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% v$Virus.name)==TRUE)),
                size = 4, fill = color.list[10] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% g$Virus.name)==TRUE)),
                size = 4, fill = color.list[4] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% gh$Virus.name)==TRUE)),
                size = 4, fill = color.list[5], color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% gv$Virus.name)==TRUE)),
                size = 4, fill = color.list[3] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% gr$Virus.name)==TRUE)),
                size = 4, fill = color.list[6] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% gry$Virus.name)==TRUE)),
                size = 4, fill = color.list[1] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% gk$Virus.name)==TRUE)),
                size = 4, fill = color.list[2] , color= color.list[2], shape = 21) +
  # geom_tippoint(aes(subset=((label %in% b16171$Virus.name)==TRUE)),
  #               size = 8, fill = "blue" , color= color.list[2], shape = 21) +
  
  theme(legend.position = "none")
dev.off()

table(subset(tree.detail, Virus.name %in% gh$Virus.name)$Clade)
subset(tree.detail, Virus.name %in% gh$Virus.name & Clade=="G")
