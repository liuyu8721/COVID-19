rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/all.meta.RData")

###################################################################################################
### phylogenetic tree
###################################################################################################

tree <- read.tree("./data/India/India_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("557"))
# pdf("India_tree_node.pdf", width = 7, height = 14)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
# dev.off()

library(RColorBrewer)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
tree2 <- groupClade(tree, c("471", "447", "454", "573", "119",
                            "373", "697", "2793", "1523", "4291", "6375", "8291"))
ggtree(tree2, aes(color=group)) +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c("grey",
                              color.list[9], #S
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10], #V
                              color.list[4],  #G:697
                              color.list[5],  #GH:2793
                              color.list[3],  #GV:1523
                              color.list[2],  #GK:4291
                              color.list[6],  #GR:3590
                              color.list[1],  #GRY:8291
                              rep("grey99", 10)))  +
  theme(legend.position = "none")

###################################################################################################
### Mark the detected samples 
###################################################################################################

all.meta$Virus.name <- sub(pattern = " ", replacement = "_", x = all.meta$Virus.name)
tree.detail <- subset(all.meta, Virus.name %in% tree$tip.label)

# Subgroup
subgroup01 <- read.delim2("./data/India/ValidationData/subgroup01.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup01 <- subset(tree.detail, Accession.ID %in% subgroup01$gisaid_epi_isl)
subgroup02 <- read.delim2("./data/India/ValidationData/subgroup02.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup02 <- subset(tree.detail, Accession.ID %in% subgroup02$gisaid_epi_isl)
subgroup03 <- read.delim2("./data/India/ValidationData/subgroup03.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03 <- subset(tree.detail, Accession.ID %in% subgroup03$gisaid_epi_isl)
subgroup04 <- read.delim2("./data/India/ValidationData/subgroup04.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup04 <- subset(tree.detail, Accession.ID %in% subgroup04$gisaid_epi_isl)
subgroup05 <- read.delim2("./data/India/ValidationData/subgroup05.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05 <- subset(tree.detail, Accession.ID %in% subgroup05$gisaid_epi_isl)
subgroup06 <- read.delim2("./data/India/ValidationData/subgroup06.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup06 <- subset(tree.detail, Accession.ID %in% subgroup06$gisaid_epi_isl)
subgroup07 <- read.delim2("./data/India/ValidationData/subgroup07.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07 <- subset(tree.detail, Accession.ID %in% subgroup07$gisaid_epi_isl)
subgroup08 <- read.delim2("./data/India/ValidationData/subgroup08.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08 <- subset(tree.detail, Accession.ID %in% subgroup08$gisaid_epi_isl)
subgroup09 <- read.delim2("./data/India/ValidationData/subgroup09.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup09 <- subset(tree.detail, Accession.ID %in% subgroup09$gisaid_epi_isl)
subgroup10 <- read.delim2("./data/India/ValidationData/subgroup10.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup10 <- subset(tree.detail, Accession.ID %in% subgroup10$gisaid_epi_isl)
subgroup11 <- read.delim2("./data/India/ValidationData/subgroup11.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup11 <- subset(tree.detail, Accession.ID %in% subgroup11$gisaid_epi_isl)
subgroup12 <- read.delim2("./data/India/ValidationData/subgroup12.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup12 <- subset(tree.detail, Accession.ID %in% subgroup12$gisaid_epi_isl)
subgroup13 <- read.delim2("./data/India/ValidationData/subgroup13.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup13 <- subset(tree.detail, Accession.ID %in% subgroup13$gisaid_epi_isl)
subgroup14 <- read.delim2("./data/India/ValidationData/subgroup14.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14 <- subset(tree.detail, Accession.ID %in% subgroup14$gisaid_epi_isl)
subgroup16 <- read.delim2("./data/India/ValidationData/subgroup16.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup16 <- subset(tree.detail, Accession.ID %in% subgroup16$gisaid_epi_isl)
subgroup17 <- read.delim2("./data/India/ValidationData/subgroup17.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup17 <- subset(tree.detail, Accession.ID %in% subgroup17$gisaid_epi_isl)
subgroup18 <- read.delim2("./data/India/ValidationData/subgroup18.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup18 <- subset(tree.detail, Accession.ID %in% subgroup18$gisaid_epi_isl)
subgroup19 <- read.delim2("./data/India/ValidationData/subgroup19.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup19 <- subset(tree.detail, Accession.ID %in% subgroup19$gisaid_epi_isl)
subgroup20 <- read.delim2("./data/India/ValidationData/subgroup20.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup20 <- subset(tree.detail, Accession.ID %in% subgroup20$gisaid_epi_isl)
subgroup21 <- read.delim2("./data/India/ValidationData/subgroup21.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup21 <- subset(tree.detail, Accession.ID %in% subgroup21$gisaid_epi_isl)
subgroup22 <- read.delim2("./data/India/ValidationData/subgroup22.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup22 <- subset(tree.detail, Accession.ID %in% subgroup22$gisaid_epi_isl)
subgroup23 <- read.delim2("./data/India/ValidationData/subgroup23.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup23 <- subset(tree.detail, Accession.ID %in% subgroup23$gisaid_epi_isl)
subgroup24 <- read.delim2("./data/India/ValidationData/subgroup24.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup24 <- subset(tree.detail, Accession.ID %in% subgroup24$gisaid_epi_isl)
subgroup25 <- read.delim2("./data/India/ValidationData/subgroup25.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup25 <- subset(tree.detail, Accession.ID %in% subgroup25$gisaid_epi_isl)
subgroup26 <- read.delim2("./data/India/ValidationData/subgroup26.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup26 <- subset(tree.detail, Accession.ID %in% subgroup26$gisaid_epi_isl)
subgroup27 <- read.delim2("./data/India/ValidationData/subgroup27.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup27 <- subset(tree.detail, Accession.ID %in% subgroup27$gisaid_epi_isl)
subgroup28 <- read.delim2("./data/India/ValidationData/subgroup28.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup28 <- subset(tree.detail, Accession.ID %in% subgroup28$gisaid_epi_isl)
subgroup29 <- read.delim2("./data/India/ValidationData/subgroup29.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup29 <- subset(tree.detail, Accession.ID %in% subgroup29$gisaid_epi_isl)
subgroup31 <- read.delim2("./data/India/ValidationData/subgroup31.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup31 <- subset(tree.detail, Accession.ID %in% subgroup31$gisaid_epi_isl)

# pdf("India_tree_highlight.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = TRUE, size = 0.2) + 
  geom_tippoint(size = 1, fill = "white", color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  # geom_nodelab(size = 1) +
  scale_color_manual(values=c("grey",
                              color.list[9], #S
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10], #V
                              color.list[4],  #G:697
                              color.list[5],  #GH:2793
                              color.list[3],  #GV:1523
                              color.list[2],  #GK:2897
                              color.list[6],  #GR:3590
                              color.list[1],  #GRY:8291
                              rep("grey99", 10))) +
  
  geom_tippoint(aes(subset=((label %in% subgroup01$Virus.name)==TRUE)),             #B.1.1.7
                size = 3, fill = "#66C2A5", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup28$Virus.name)==TRUE)),             #B.1.1. 
                size = 3, fill = "skyblue", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup06$Virus.name)==TRUE)),             #B.1.1.526
                size = 3, fill = "red", color= 'red', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup24$Virus.name)==TRUE)),             #B.1.1.306:1438
                size = 3, fill = "darkgrey", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup18$Virus.name)==TRUE)),             #B.1.1.306:278
                size = 3, fill = "lightgreen", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup19$Virus.name)==TRUE)),             #B.1.1.216
                size = 3, fill = "lightyellow", color= 'grey66', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup21$Virus.name)==TRUE)),             #B.1.326
                size = 3, fill = "darkred", color= 'white', shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup06$Virus.name)==TRUE)),             #B.1.1.526
                size = 2, fill = "red", color= 'white', shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup05$Virus.name)==TRUE)),             #B.1.617:9340
                size = 2, fill = "green", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup25$Virus.name)==TRUE)),             #B.1.617:2120
                size = 2, fill = "orange", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup02$Virus.name)==TRUE)),             #B.1.617.1
                size = 3, fill = "blue", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup12$Virus.name)==TRUE)),             #B.1.617.1
                size = 3, fill = "tan", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup04$Virus.name)==TRUE)),             #B.1.617.2:4848
                size = 2, fill = "#FC8D62", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup27$Virus.name)==TRUE)),             #B.1.617.2:2764
                size = 2, fill = "purple", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup20$Virus.name)==TRUE)),             #B.1.617.2:4375
                size = 2, fill = "royalblue", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup09$Virus.name)==TRUE)),             #B.1.617.2:2435
                size = 3, fill = "magenta", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup22$Virus.name)==TRUE)),             #B.1.617.2:524
                size = 2, fill = adjustcolor("darkgreen", 0.1), color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup13$Virus.name)==TRUE)),             #B.1.617.2:935
                size = 3, fill = adjustcolor("salmon", alpha.f = 0.6), color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup14$Virus.name)==TRUE)),             #B.1.617.2: 602
                size = 3, fill = "cyan", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup03$Virus.name)==TRUE)),             #B.1.617.2:2209
                size = 3, fill = "brown", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup16$Virus.name)==TRUE)),             #B.1.617.2:1042
                size = 3, fill = "#E78AC3", color= 'white', shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup07$Virus.name)==TRUE)),             #B.1.113 
                size = 2, fill = "black", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup10$Virus.name)==TRUE)),             #B.1.36
                size = 3, fill = "#7FC97F", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08$Virus.name)==TRUE)),             #B.1.36.29
                size = 3, fill = "pink", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup23$Virus.name)==TRUE)),             #B.1.36
                size = 3, fill = "darkturquoise", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup26$Virus.name)==TRUE)),             #B.1.210
                size = 3, fill = "darkorange", color= 'white', shape = 21) +
  
  geom_tippoint(aes(subset=((label %in% subgroup11$Virus.name)==TRUE)),             #B.6
                size = 3, fill = "greenyellow", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup29$Virus.name)==TRUE)),             #B.6.6
                size = 3, fill = "saddlebrown", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup17$Virus.name)==TRUE)),             #B.4
                size = 3, fill = "grey60", color= 'white', shape = 21) +
 
  geom_tippoint(aes(subset=((label %in% subgroup31$Virus.name)==TRUE)),             #A
                size = 3, fill = "#386CB0", color= 'white', shape = 21) +
  
  theme(legend.position = "none")
# dev.off()
