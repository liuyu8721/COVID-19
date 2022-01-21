rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")
library("treeio")
suppressMessages(library(tidyverse))
Sys.setlocale("LC_TIME", "English")

setwd("D:/00SARS-Cov-2/ADS")

library(RColorBrewer)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))

###################################################################################################

tree <- read.tree("D:\\00SARS-Cov-2\\ADS\\India_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("130")) #130
pdf("India_tree_node.pdf", width = 25, height = 50)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
dev.off()

###################################################################################################

fw.detail <- read.delim2("D:/00SARS-Cov-2/ADS/India/metadata.tsv", header = T)

#---Background
load("D:/00SARS-Cov-2/ADS/RData/all.meta.update.RData")
all.meta.update$Virus.name <- sub(" ", "_", x = all.meta.update$Virus.name)
tree.detail <- subset(all.meta.update, Virus.name %in% tree$tip.label)
l <- subset(tree.detail, Clade=="L", select = Virus.name)
s <- subset(tree.detail, Clade=="S", select = Virus.name)
o <- subset(tree.detail, Clade=="O", select = Virus.name)
v <- subset(tree.detail, Clade=="V", select = Virus.name)
g <- subset(tree.detail, Clade %in% c("G", "GH", "GV", "GR", "GRY", "GK"), select = Virus.name)
gh <- subset(tree.detail, Clade=="GH", select = Virus.name)
gv <- subset(tree.detail, Clade=="GV", select = Virus.name)
gr <- subset(tree.detail, Clade %in% c("GR", "GRY"), select = Virus.name)
gry <- subset(tree.detail, Clade=="GRY", select = Virus.name)
gk <- subset(tree.detail, Clade=="GK", select = Virus.name)

pdf("India_Skeleton.pdf", width = 7, height = 14)
ggtree(tree, mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  geom_nodelab() +
  geom_tippoint(aes(subset=((label %in% l$Virus.name)==TRUE)),
                size = 4, fill = "blue", color= 'grey88', shape = 21)  +
  geom_tippoint(aes(subset=((label %in% s$Virus.name)==TRUE)),
                size = 4, fill = color.list[9], color= 'grey88', shape = 21)  +
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
                size = 4, fill = color.list[2] , color= 'grey88', shape = 21)
dev.off()

tree2 <- groupClade(tree, c("130", "84", "271", "590", "5099", "1848", "1153", "6287", "8530"))
pdf("India_Skeleton_ColorEdge.pdf", width = 25, height = 50)
ggtree(tree2, aes(color=group)) +
  # geom_nodelab() +
  scale_color_manual(values=c(color.list[9],  #S: 130
                              "grey",         #L
                              color.list[10], #V:  271
                              color.list[4],  #G:  590
                              color.list[5],  #GH: 5099
                              color.list[3],  #GV: 1848
                              color.list[2],  #GK: 1153
                              color.list[6],  #GR: 6287
                              color.list[1],  #GRY:8530
                              rep("grey99", 10))) +
  geom_tippoint(aes(subset=((label %in% l$Virus.name)==TRUE)),
                size = 4, fill = "blue", color= 'grey88', shape = 21)  +
  geom_tippoint(aes(subset=((label %in% s$Virus.name)==TRUE)),
                size = 4, fill = color.list[9], color= 'grey88', shape = 21)  +
  geom_tippoint(aes(subset=((label %in% v$Virus.name)==TRUE)),
                size = 4, fill = color.list[10] , color= 'grey88', shape = 21) +
  geom_tippoint(aes(subset=((label %in% o$Virus.name)==TRUE)),
                size = 4, fill = color.list[7] , color= 'grey88', shape = 21) +
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
                size = 4, fill = color.list[2] , color= 'grey88', shape = 21)
dev.off()

###################################################################################################

#---Subgroup
subgroup01 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_01.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup01 <- subset(fw.detail, gisaid_epi_isl %in% subgroup01$gisaid_epi_isl)
subgroup02 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_02.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup02 <- subset(fw.detail, gisaid_epi_isl %in% subgroup02$gisaid_epi_isl)
subgroup03 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_03.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03 <- subset(fw.detail, gisaid_epi_isl %in% subgroup03$gisaid_epi_isl)
subgroup04 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_04.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup04 <- subset(fw.detail, gisaid_epi_isl %in% subgroup04$gisaid_epi_isl)
subgroup05 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_05.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05 <- subset(fw.detail, gisaid_epi_isl %in% subgroup05$gisaid_epi_isl)
treex <- ape::extract.clade(tree, node = c("5503"))
subgroup06 <- subset(fw.detail, strain %in% treex$tip.label)
# subgroup06 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_06.txt",
#                           header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
# subgroup06 <- subset(fw.detail, gisaid_epi_isl %in% subgroup06$gisaid_epi_isl)
subgroup07 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_07.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07 <- subset(fw.detail, gisaid_epi_isl %in% subgroup07$gisaid_epi_isl)
subgroup08 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_08.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08 <- subset(fw.detail, gisaid_epi_isl %in% subgroup08$gisaid_epi_isl)
subgroup09 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_09.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup09 <- subset(fw.detail, gisaid_epi_isl %in% subgroup09$gisaid_epi_isl)
subgroup10 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_10.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup10 <- subset(fw.detail, gisaid_epi_isl %in% subgroup10$gisaid_epi_isl)
subgroup11 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_11.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup11 <- subset(fw.detail, gisaid_epi_isl %in% subgroup11$gisaid_epi_isl)
subgroup12 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_12.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup12 <- subset(fw.detail, gisaid_epi_isl %in% subgroup12$gisaid_epi_isl)
subgroup13 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_13.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup13 <- subset(fw.detail, gisaid_epi_isl %in% subgroup13$gisaid_epi_isl)
subgroup14 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_14.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14 <- subset(fw.detail, gisaid_epi_isl %in% subgroup14$gisaid_epi_isl)
subgroup15 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_15.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup15 <- subset(fw.detail, gisaid_epi_isl %in% subgroup15$gisaid_epi_isl)
subgroup16 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_16.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup16 <- subset(fw.detail, gisaid_epi_isl %in% subgroup16$gisaid_epi_isl)
subgroup17 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_17.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup17 <- subset(fw.detail, gisaid_epi_isl %in% subgroup17$gisaid_epi_isl)
subgroup18 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_18.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup18 <- subset(fw.detail, gisaid_epi_isl %in% subgroup18$gisaid_epi_isl)
subgroup19 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_19.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup19 <- subset(fw.detail, gisaid_epi_isl %in% subgroup19$gisaid_epi_isl)
subgroup20 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_20.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup20 <- subset(fw.detail, gisaid_epi_isl %in% subgroup20$gisaid_epi_isl)
subgroup21 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_21.txt",
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup21 <- subset(fw.detail, gisaid_epi_isl %in% subgroup21$gisaid_epi_isl)
subgroup23 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_23.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup23 <- subset(fw.detail, gisaid_epi_isl %in% subgroup23$gisaid_epi_isl)
subgroup24 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_24.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup24 <- subset(fw.detail, gisaid_epi_isl %in% subgroup24$gisaid_epi_isl) %>%
  filter(gisaid_epi_isl %in% subset(all.meta.update, Lineage=="B.1.1.326")$Accession.ID)
subgroup25 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_25.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup25 <- subset(fw.detail, gisaid_epi_isl %in% subgroup25$gisaid_epi_isl)
subgroup26 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_26.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup26 <- subset(fw.detail, gisaid_epi_isl %in% subgroup26$gisaid_epi_isl)
subgroup27 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_27.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup27 <- subset(fw.detail, gisaid_epi_isl %in% subgroup27$gisaid_epi_isl)
subgroup28 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_28.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup28 <- subset(fw.detail, gisaid_epi_isl %in% subgroup28$gisaid_epi_isl)
subgroup29 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_29.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup29 <- subset(fw.detail, gisaid_epi_isl %in% subgroup29$gisaid_epi_isl)
subgroup30 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_30.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup30 <- subset(fw.detail, gisaid_epi_isl %in% subgroup30$gisaid_epi_isl)
subgroup32 <- read.delim2("D:/00SARS-Cov-2/ADS/India/subgroup/subgroup_32.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup32 <- subset(fw.detail, gisaid_epi_isl %in% subgroup32$gisaid_epi_isl)

pdf("India_tree.pdf", width = 7, height = 14)
ggtree(tree2, aes(color = group),mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c(color.list[9], #S: 221
                              "grey", #L
                              color.list[10], #V:  112
                              color.list[4],  #G:  589
                              color.list[5],  #GH: 3472
                              color.list[3],  #GV: 1848
                              color.list[2],  #GK: 1153
                              color.list[6],  #GR: 6287
                              color.list[1],  #GRY:8530
                              rep("grey99", 10))) +
  
  # geom_tippoint(aes(subset=((label %in% subgroup11$strain)==TRUE)),
  #               size = 4, fill = "#E78AC3", color= "white", shape = 21) +       #G/B.1
  # geom_tippoint(aes(subset=((label %in% subgroup17$strain)==TRUE)),
  #               size = 4, fill = "#BEAED4", color= "white", shape = 21) +       #GR/B.1.1
  
  geom_tippoint(aes(subset=((label %in% subgroup04$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[6], color= "white", shape = 21) +   #B.1.672
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[7], color= "white", shape = 21) +   #B.1.672.2
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[7], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup23$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[3], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup25$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup32$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[12], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup12$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup02$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[2], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup26$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[6], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup27$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Accent")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup09$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[3], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup13$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[6], color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup19$strain)==TRUE)),
                size = 4, fill = "#FFFF99", color= "grey88", shape = 21) +   #B.1.672.1
  geom_tippoint(aes(subset=((label %in% subgroup03$strain)==TRUE)),
                size = 4, fill = "#A6D854", color= "grey88", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup01$strain)==TRUE)),
                size = 4, fill = "#66C2A5", color= "white", shape = 21) +   #B.1.1.7
  geom_tippoint(aes(subset=((label %in% subgroup16$strain)==TRUE)),
                size = 4, fill = "#FDC086", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup08$strain)==TRUE)),
                size = 4, fill = "#1F78B4", color= "white", shape = 21) +   #B.1.36
  geom_tippoint(aes(subset=((label %in% subgroup07$strain)==TRUE)),
                size = 4, fill = "#A6CEE3", color= "white", shape = 21) +   #B.1.36.29

  geom_tippoint(aes(subset=((label %in% subgroup29$strain)==TRUE)),
                size = 4, fill = "#D9D9D9", color= "grey66", shape = 21) +   #B.1.1.306
  geom_tippoint(aes(subset=((label %in% subgroup20$strain)==TRUE)),
                size = 4, fill = "#B15928", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup18$strain)==TRUE)),
                size = 4, fill = "#6A3D9A", color= "grey66", shape = 21) +   #B.1.1.326
  geom_tippoint(aes(subset=((label %in% subgroup24$strain)==TRUE)),
                size = 4, fill = "#FB8072", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup10$strain)==TRUE)),
                size = 4, fill = "#33A02C", color= "grey66", shape = 21) +      #B.6
  geom_tippoint(aes(subset=((label %in% subgroup30$strain)==TRUE)),
                size = 4, fill = "#BC80BD", color= "white", shape = 21) +       #B.6.6

  geom_tippoint(aes(subset=((label %in% subgroup28$strain)==TRUE)),
                size = 4, fill = "#FCCDE5", color= "white", shape = 21) +         #B.1.1.216
  geom_tippoint(aes(subset=((label %in% subgroup15$strain)==TRUE)),
                size = 4, fill = "#FF7F00", color= "white", shape = 21) +         #B.1.1.46
  geom_tippoint(aes(subset=((label %in% subgroup14$strain)==TRUE)),
                size = 4, fill = "#FDBF6F", color= "white", shape = 21) +         #B.4
  
  geom_tippoint(aes(subset=((label %in% subgroup21$strain)==TRUE)),
                size = 4, fill = "#8DD3C7", color= "white", shape = 21) +       #B.1.210
  geom_tippoint(aes(subset=((label %in% subgroup06$strain)==TRUE)),
                size = 4, fill = "#B3B3B3", color= "black", shape = 21) +      #B.1.113

  theme(legend.position = "none")
dev.off()

###################################################################################################
### B.617
###################################################################################################

tree3 <- ape::extract.clade(tree, node = c("2658"))
tree3 <- groupClade(tree3, c("2658", "2676", "3826", "3827", "4419", "4637"))
pdf("India_tree_Delta.pdf", width = 7, height = 7)
ggtree(tree3, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c(color.list[2],  #GK: 2658
                              adjustcolor("#FC8D62", 0.6),
                              adjustcolor("#FB9A99", 0.6),
                              adjustcolor("#386CB0", 0.6),
                              adjustcolor("#B2DF8A", 0.6),
                              adjustcolor("#E31A1C", 0.6),
                              rep("grey99", 10))) +
  # geom_tippoint(aes(subset=((label %in% subgroup04$strain)==TRUE)),
  #               size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[6], color= "white", shape = 21) +   #B.1.672
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[7], color= "white", shape = 21) +   #B.1.672.2
  geom_tippoint(aes(subset=((label %in% subgroup23$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[3], color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup25$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup32$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[12], color= "white", shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup12$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup02$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Set2")[2], color= "white", shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup26$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(12, "Set3")[6], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup27$strain)==TRUE)),
              size = 4, fill = RColorBrewer::brewer.pal(8, "Accent")[5], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup09$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[3], color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup13$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(8, "Paired")[6], color= "white", shape = 21) +
  
  theme(legend.position = "none")
dev.off()

###################################################################################################

T1130 <- read_tsv("F:/Genomes/GISAID/NotHighCoverage/BASE/metadata.tsv")
names(T1130) <- c("Virus.name",
               "Type",
               "Accession.ID",
               "Collection.date",
               "Location",
               "Additional.location.information",
               "Sequence.length",
               "Host",
               "Patient.age",
               "Gender",
               "Clade",
               "Pango.lineage",                  
               "Pangolin.version",
               "Variant",
               "AA.Substitutions",
               "Submission.date",
               "Is.reference",
               "Is.complete",
               "Is.high.coverage",
               "Is.low.coverage",
               "N-Content",
               "GC-Content")
T1130$Pangolin.version[1]

setwd("D:/00SARS-Cov-2/ADS/India")
load("./India.meta.update.RData")
T1130.India <- subset(T1130, Accession.ID %in% India.meta.update$Accession.ID)
input <- subset(T1130.India, Accession.ID %in% subgroup32$gisaid_epi_isl)
lineage <- data.frame(table(input$Pango.lineage))
lineage <- lineage[order(lineage$Freq, decreasing = T),]
lineage

