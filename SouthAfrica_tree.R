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

tree <- read.tree("D:\\00SARS-Cov-2\\ADS\\SouthAfrica_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("152")) #152
pdf("SouthAfrica_tree_node.pdf", width = 25, height = 50)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
dev.off()

###################################################################################################

fw.detail <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/metadata.tsv", header = T)

#---Background
load("D:/00SARS-Cov-2/ADS/RData/all.meta.update.RData")
all.meta.update$Virus.name <- sub(" ", "_", x = all.meta.update$Virus.name)
tree.detail <- subset(all.meta.update, Virus.name %in% tree$tip.label)
l <- subset(tree.detail, Clade=="L", select = Virus.name)
s <- subset(tree.detail, Clade=="S", select = Virus.name)
v <- subset(tree.detail, Clade=="V", select = Virus.name)
o <- subset(tree.detail, Clade=="O", select = Virus.name)
g <- subset(tree.detail, Clade %in% c("G", "GH", "GV", "GR", "GRY", "GK"), select = Virus.name)
gh <- subset(tree.detail, Clade=="GH", select = Virus.name)
gv <- subset(tree.detail, Clade=="GV", select = Virus.name)
gr <- subset(tree.detail, Clade %in% c("GR", "GRY"), select = Virus.name)
gry <- subset(tree.detail, Clade=="GRY", select = Virus.name)
gk <- subset(tree.detail, Clade=="GK", select = Virus.name)

pdf("SouthAfrica_Skeleton.pdf", width = 25, height = 50)
ggtree(tree, mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  geom_nodelab() +
  geom_tippoint(aes(subset=((label %in% l$Virus.name)==TRUE)),
                size = 4, fill = "blue", color= 'grey88', shape = 21)  +
  geom_tippoint(aes(subset=((label %in% s$Virus.name)==TRUE)),
                size = 4, fill = color.list[9], color= 'grey88', shape = 21)  +
  geom_tippoint(aes(subset=((label %in% o$Virus.name)==TRUE)),
                size = 4, fill = "green", color= 'grey88', shape = 21)  +
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

tree2 <- groupClade(tree, c("152", "33", "34", "371", "2256", "917", "3448", "5142", "6992"))
pdf("SouthAfrica_Skeleton_ColorEdge.pdf", width = 25, height = 50)
ggtree(tree2, aes(color=group)) +
  geom_nodelab() +
  scale_color_manual(values=c(color.list[9], #S: 231
                              "grey", #L
                              color.list[10], #V:  34
                              color.list[4],  #G:  371
                              color.list[5],  #GH: 2256
                              color.list[3],  #GV: 917
                              color.list[2],  #GK: 3448
                              color.list[6],  #GR: 5142
                              color.list[1],  #GRY:6992
                              rep("grey99", 10))) 
dev.off()

###################################################################################################

#---Subgroup
subgroup01 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_01.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup01 <- subset(fw.detail, gisaid_epi_isl %in% subgroup01$gisaid_epi_isl)
subgroup02 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_02.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup02 <- subset(fw.detail, gisaid_epi_isl %in% subgroup02$gisaid_epi_isl)
subgroup03 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_03.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03 <- subset(fw.detail, gisaid_epi_isl %in% subgroup03$gisaid_epi_isl)
subgroup04 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_04.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup04 <- subset(fw.detail, gisaid_epi_isl %in% subgroup04$gisaid_epi_isl)
subgroup05 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_05.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05 <- subset(fw.detail, gisaid_epi_isl %in% subgroup05$gisaid_epi_isl)
subgroup06 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_06.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup06 <- subset(fw.detail, gisaid_epi_isl %in% subgroup06$gisaid_epi_isl)
subgroup07 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_07.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07 <- subset(fw.detail, gisaid_epi_isl %in% subgroup07$gisaid_epi_isl)
subgroup08 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_08.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08 <- subset(fw.detail, gisaid_epi_isl %in% subgroup08$gisaid_epi_isl)
subgroup09 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_09.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup09 <- subset(fw.detail, gisaid_epi_isl %in% subgroup09$gisaid_epi_isl)
subgroup10 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_10.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup10 <- subset(fw.detail, gisaid_epi_isl %in% subgroup10$gisaid_epi_isl) %>%
  filter(!(strain %in% "hCoV-19/South_Africa/VIDA-KRISP-K012173/2020"))
subgroup11 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_11.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup11 <- subset(fw.detail, gisaid_epi_isl %in% subgroup11$gisaid_epi_isl)
subgroup12 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_12.txt", #20211214/subgroup_15.txt
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup12 <- subset(fw.detail, gisaid_epi_isl %in% subgroup12$gisaid_epi_isl)
subgroup13 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_13.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup13 <- subset(fw.detail, gisaid_epi_isl %in% subgroup13$gisaid_epi_isl)
subgroup14 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_14.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14 <- subset(fw.detail, gisaid_epi_isl %in% subgroup14$gisaid_epi_isl) %>%
  filter(!(strain %in% c("hCoV-19/South_Africa/KRISP-K002686/2020",
                         "hCoV-19/South_Africa/VIDA-KRISP-K012406/2020",
                         "hCoV-19/South_Africa/VIDA-KRISP-K012551/2020")))
subgroup15 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_15.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup15 <- subset(fw.detail, gisaid_epi_isl %in% subgroup15$gisaid_epi_isl) %>%
  filter(!(strain %in% "hCoV-19/South_Africa/VIDA-KRISP-K012173/2020"))
subgroup16 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_16.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup16 <- subset(fw.detail, gisaid_epi_isl %in% subgroup16$gisaid_epi_isl)
subgroup17 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_17.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup17 <- subset(fw.detail, gisaid_epi_isl %in% subgroup17$gisaid_epi_isl)
subgroup18 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_18.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup18 <- subset(fw.detail, gisaid_epi_isl %in% subgroup18$gisaid_epi_isl)
subgroup19 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_19.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup19 <- subset(fw.detail, gisaid_epi_isl %in% subgroup19$gisaid_epi_isl)
subgroup20 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_20.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup20 <- subset(fw.detail, gisaid_epi_isl %in% subgroup20$gisaid_epi_isl)
subgroup23 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_23.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup23 <- subset(fw.detail, gisaid_epi_isl %in% subgroup23$gisaid_epi_isl)
subgroup24 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_24.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup24 <- subset(fw.detail, gisaid_epi_isl %in% subgroup24$gisaid_epi_isl)
subgroup25 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_25.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup25 <- subset(fw.detail, gisaid_epi_isl %in% subgroup25$gisaid_epi_isl) %>%
  filter(!(strain %in% c("hCoV-19/South_Africa/CERI-KRISP-K032191/2021")))
subgroup26 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_26.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup26 <- subset(fw.detail, gisaid_epi_isl %in% subgroup26$gisaid_epi_isl)
subgroup27 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup_27.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup27 <- subset(fw.detail, gisaid_epi_isl %in% subgroup27$gisaid_epi_isl)

pdf("SouthAfrica_tree.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c(color.list[9], #S: 231
                              "grey", #L
                              color.list[10], #V:  34
                              color.list[4],  #G:  371
                              color.list[5],  #GH: 2256
                              color.list[3],  #GV: 917
                              color.list[2],  #GK: 3448
                              color.list[6],  #GR: 5142
                              color.list[1],  #GRY:6992
                              rep("grey99", 10))) +
  
  # geom_tippoint(aes(subset=((label %in% subgroup12$strain)==TRUE)),
  #               size = 4, fill = "#E78AC3", color= "white", shape = 21) +       #G/B.1
  # geom_tippoint(aes(subset=((label %in% subgroup16$strain)==TRUE)),
  #               size = 4, fill = "#BEAED4", color= "white", shape = 21) +       #GR/B.1.1
  
  geom_tippoint(aes(subset=((label %in% subgroup01$strain)==TRUE)),             #Omicron
                size = 4, fill = "red", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = "green", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup03$strain)==TRUE)),             #Alpha
                size = 4, fill = "#66C2A5", color= "grey88", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup15$strain)==TRUE)),
                size = 4, fill = "lightyellow", color= "black", shape = 21) +   #C.1
  geom_tippoint(aes(subset=((label %in% subgroup10$strain)==TRUE)),
                size = 4, fill = "tan", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08$strain)==TRUE)),
                size = 4, fill = "greenyellow", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup17$strain)==TRUE)),
                size = 4, fill = "darkturquoise", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup04$strain)==TRUE)),
                size = 4, fill = "#7FC97F", color= "grey88", shape = 21) +       #B.1.351
  geom_tippoint(aes(subset=((label %in% subgroup09$strain)==TRUE)),
                size = 4, fill = "black", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup07$strain)==TRUE)),
                size = 4, fill = "purple", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup11$strain)==TRUE)),
                size = 4, fill = "#FA8072CC", color= "black", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup24$strain)==TRUE)),
                size = 4, fill = "skyblue", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup25$strain)==TRUE)),
                size = 4, fill = "saddlebrown", color= "white", shape = 21) +   #B.1.672.2
  geom_tippoint(aes(subset=((label %in% subgroup02$strain)==TRUE)),
                size = 4, fill = "#FC8D62", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup26$strain)==TRUE)),
                size = 4, fill = "darkorange", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup18$strain)==TRUE)),
                size = 4, fill = "lightcyan", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup13$strain)==TRUE)),
                size = 4, fill = "grey60", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup06$strain)==TRUE)),
                size = 4, fill = "#8DA0CB", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup23$strain)==TRUE)),
                size = 4, fill = "#666666", color= "white", shape = 21) +

  geom_tippoint(aes(subset=((label %in% subgroup14$strain)==TRUE)),
                size = 4, fill = "lightgreen", color= "white", shape = 21) +    #B.1.1.448
  geom_tippoint(aes(subset=((label %in% subgroup20$strain)==TRUE)),
                size = 4, fill = "cyan", color= "white", shape = 21) +          #B.1.140
  geom_tippoint(aes(subset=((label %in% subgroup19$strain)==TRUE)),
                size = 4, fill = "darkgrey", color= "white", shape = 21) +      #B.1.1.459
  geom_tippoint(aes(subset=((label %in% subgroup27$strain)==TRUE)),
                size = 4, fill = "steelblue", color= "white", shape = 21) +     #B.1.8

  theme(legend.position = "none")
dev.off()

# temp <- subset(tree.detail, Virus.name %in% subgroup15$strain)
# table(temp$Clade)
# table(temp$Lineage)
# subset(temp, !grepl("B.1.160|B.1.36", Lineage))
# 
# test <- subset(tree.detail, Virus.name %in% temp$strain)
# test <- subset(fw.detail, gisaid_epi_isl %in% test$Accession.ID)
# 
# tree$tip.label
# temp <- subset(subgroup15, strain %in% tree$tip.label)


###################################################################################################
### B.1.351
###################################################################################################

tree3 <- ape::extract.clade(tree, node = c("2686"))
tree3 <- groupClade(tree3, c("2686"))
pdf("SouthAfrica_tree_Beta.pdf", width = 7, height = 5)
ggtree(tree3, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = 'gray80', shape = 21) +
  scale_color_manual(values=color.list[5]) +
  geom_tippoint(aes(subset=((label %in% subgroup04$strain)==TRUE)),
                size = 4, fill = "#7FC97F", color= "grey88", shape = 21) +       #B.1.351
  geom_tippoint(aes(subset=((label %in% subgroup09$strain)==TRUE)),
                size = 4, fill = "black", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup07$strain)==TRUE)),
                size = 4, fill = "purple", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup11$strain)==TRUE)),
                size = 4, fill = "#FA8072CC", color= "black", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup24$strain)==TRUE)),
                size = 4, fill = "skyblue", color= "white", shape = 21) +
  theme(legend.position = "none")
dev.off()

tree4 <- ape::extract.clade(tree, node = c("6146"))
tree4 <- groupClade(tree4, c("6146"))
pdf("SouthAfrica_tree_omicron.pdf", width = 7, height = 1)
ggtree(tree4, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=color.list[6]) +
  geom_tippoint(aes(subset=((label %in% subgroup01$strain)==TRUE)),             #Omicron
                size = 5, fill = "red", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = "green", color= "white", shape = 21) +
  theme(legend.position = "none")
dev.off()

tree5 <- ape::extract.clade(tree, node = c("3448"))
tree5 <- groupClade(tree5, c("3448"))
pdf("SouthAfrica_tree_delta.pdf", width = 7, height = 5)
ggtree(tree5, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=color.list[2]) +
  geom_tippoint(aes(subset=((label %in% subgroup25$strain)==TRUE)),
                size = 4, fill = "saddlebrown", color= "white", shape = 21) +   #B.1.672.2
  geom_tippoint(aes(subset=((label %in% subgroup02$strain)==TRUE)),
                size = 4, fill = "#FC8D62", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup26$strain)==TRUE)),
                size = 4, fill = "darkorange", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup18$strain)==TRUE)),
                size = 4, fill = "lightcyan", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup13$strain)==TRUE)),
                size = 4, fill = "grey60", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup06$strain)==TRUE)),
                size = 4, fill = "#8DA0CB", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup23$strain)==TRUE)),
                size = 4, fill = "#666666", color= "white", shape = 21) +
  theme(legend.position = "none")
dev.off()

tree6 <- ape::extract.clade(tree, node = c("6196"))
tree6 <- groupClade(tree6, c("6196"))
pdf("SouthAfrica_tree_C1.pdf", width = 7, height = 3)
ggtree(tree6, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=color.list[6]) +
  geom_tippoint(aes(subset=((label %in% subgroup15$strain)==TRUE)),
                size = 4, fill = "lightyellow", color= "black", shape = 21) +   #C.1
  geom_tippoint(aes(subset=((label %in% subgroup10$strain)==TRUE)),
                size = 4, fill = "tan", color= "white", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08$strain)==TRUE)),
                size = 4, fill = "greenyellow", color= "grey88", shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup17$strain)==TRUE)),
                size = 4, fill = "darkturquoise", color= "white", shape = 21) +
  
  theme(legend.position = "none")
dev.off()

table(c1$Lineage)
