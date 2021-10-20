rm(list=ls())
library("ape")
library("ggplot2")
library("ggtree")

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/all.meta.RData")

###################################################################################################
### phylogenetic tree
###################################################################################################

tree <- read.tree("./data/SouthAfrica/SouthAfrica_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("199"))
# pdf("SouthAfrica_tree_node.pdf", width = 7, height = 14)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
# dev.off()

library(RColorBrewer)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
tree2 <- groupClade(tree, c("186", "197", "200", "223", "191", "37", 
                            "101", "333", "2063", "334", "3590", "5147", "2897"))
ggtree(tree2, aes(color=group)) +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c("grey",
                              color.list[9], #S
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10], #V
                              color.list[4], #G:333
                              color.list[5], #GH:2063
                              color.list[3], #GV:334
                              color.list[6], #GR:3590
                              color.list[1], #GRY:5147
                              color.list[2], #GK:2897
                              rep("grey99", 10)))  +
  theme(legend.position = "none")

###################################################################################################
### Mark the detected samples 
###################################################################################################

all.meta$Virus.name <- sub(pattern = " ", replacement = "_", x = all.meta$Virus.name)
tree.detail <- subset(all.meta, Virus.name %in% tree$tip.label)

# Subgroup
subgroup01 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup01_AY4.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup01 <- subset(tree.detail, Accession.ID %in% subgroup01$gisaid_epi_isl)
subgroup02 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup02_Beta.txt", 
                               header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup02 <- subset(tree.detail, Accession.ID %in% subgroup02$gisaid_epi_isl)
subgroup03 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup03_Alpha.txt", 
                                header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03 <- subset(tree.detail, Accession.ID %in% subgroup03$gisaid_epi_isl)
subgroup04 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup04_B1.1.34.txt", 
                                  header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup04 <- subset(tree.detail, Accession.ID %in% subgroup04$gisaid_epi_isl)
subgroup05 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup05_B1.351+R246I.txt", 
                                       header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05 <- subset(tree.detail, Accession.ID %in% subgroup05$gisaid_epi_isl)
subgroup06 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup06_C1+Y508Y.txt", 
                                   header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup06 <- subset(tree.detail, Accession.ID %in% subgroup06$gisaid_epi_isl)
subgroup07 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup07_B1.1.54.txt", 
                                  header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07 <- subset(tree.detail, Accession.ID %in% subgroup07$gisaid_epi_isl)
subgroup08 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup08_AY11.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08 <- subset(tree.detail, Accession.ID %in% subgroup08$gisaid_epi_isl)
subgroup11 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup11_C1+T723T.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup11 <- subset(tree.detail, Accession.ID %in% subgroup11$gisaid_epi_isl)
subgroup11 <- subset(subgroup11, !(gisaid_epi_isl %in% c("EPI_ISL_2494616", "EPI_ISL_2726855")))
subgroup12 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup12_B1.351+H145Y.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup12 <- subset(tree.detail, Accession.ID %in% subgroup12$gisaid_epi_isl)
subgroup14 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup14_B1.617.2+A446V.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14 <- subset(tree.detail, Accession.ID %in% subgroup14$gisaid_epi_isl)
subgroup17 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup17_B1.8.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup17 <- subset(tree.detail, Accession.ID %in% subgroup17$gisaid_epi_isl)
subgroup18 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup18_B1.1.7+V1770F.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup18 <- subset(tree.detail, Accession.ID %in% subgroup18$gisaid_epi_isl)
subgroup19 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup19_B1.140.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup19 <- subset(tree.detail, Accession.ID %in% subgroup19$gisaid_epi_isl)
subgroup20 <- read.delim2("./data/SouthAfrica/ValidationData/subgroup20_C1+Q384H.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup20 <- subset(tree.detail, Accession.ID %in% subgroup20$gisaid_epi_isl)
subgroup20 <- subset(subgroup20, !(gisaid_epi_isl %in% c("EPI_ISL_2494616", "EPI_ISL_2726855")))
subgroup21 <- read.delim2("D:/00SARS-Cov-2/ADS/SouthAfrica/subgroup/subgroup21_B1.1.273.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup21 <- subset(tree.detail, Accession.ID %in% subgroup21$gisaid_epi_isl)

# pdf("SouthAfrica_tree_highlight.pdf", width = 7, height = 14)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c("grey",
                              color.list[9], #S
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[9],
                              color.list[10], #V
                              color.list[4], #G:333
                              color.list[5], #GH:2063
                              color.list[3], #GV:334
                              color.list[6], #GR:3590
                              color.list[1], #GRY:5147
                              color.list[2], #GK:2897
                              rep("grey99", 10))) +
  
  geom_tippoint(aes(subset=((label %in% subgroup01$Virus.name)==TRUE)),
                size = 4, fill = color.list[2], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup02$Virus.name)==TRUE)),
                size = 4, fill = "#7FC97F", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup03$Virus.name)==TRUE)),
                size = 4, fill = "#66C2A5", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup04$Virus.name)==TRUE)),
                size = 4, fill = "yellow", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup07$Virus.name)==TRUE)),
                size = 4, fill = "black", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup05$Virus.name)==TRUE)),
                size = 4, fill = "green", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup11$Virus.name)==TRUE)),
                size = 4, fill = "greenyellow", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup12$Virus.name)==TRUE)),
                size = 4, fill = "tan", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup14$Virus.name)==TRUE)),
                size = 4, fill = "cyan", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup17$Virus.name)==TRUE)),
                size = 4, fill = "grey60", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup18$Virus.name)==TRUE)),
                size = 4, fill = "lightgreen", color= 'white', shape = 21) +  
  geom_tippoint(aes(subset=((label %in% subgroup19$Virus.name)==TRUE)),
                size = 4, fill = "lightyellow", color= 'black', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup20$Virus.name)==TRUE)),
                size = 4, fill = "royalblue", color= 'black', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup21$Virus.name)==TRUE)),
                size = 4, fill = "darkred", color= 'black', shape = 21) +  
  
  geom_tippoint(aes(subset=((label %in% subgroup06$Virus.name)==TRUE)),
                size = 4, fill = "red", color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08$Virus.name)==TRUE)),
                size = 4, fill = "#8DA0CB", color= 'white', shape = 21) +
  
  theme(legend.position = "none")
# dev.off()
