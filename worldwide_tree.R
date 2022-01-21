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

tree <- read.tree("D:\\00SARS-Cov-2\\ADS\\worldwide_tree.nwk")
tree$node.label <- as.numeric(substr(tree$node.label, 6, 12))

tree <- ape::extract.clade(tree, node = c("38")) #38
pdf("worldwide_tree_node.pdf", width = 25, height = 50)
plot(tree, show.node.label = T, cex = 0.1, show.tip.label = F, edge.color = "grey90")
dev.off()

###################################################################################################

fw.detail <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/metadata.tsv", header = T)

#---Background
load("D:/00SARS-Cov-2/ADS/RData/all.meta.update.RData")
all.meta.update$Virus.name <- sub(" ", "_", x = all.meta.update$Virus.name)
tree.detail <- subset(all.meta.update, Virus.name %in% tree$tip.label)
l <- subset(tree.detail, Clade=="L", select = Virus.name)
s <- subset(tree.detail, Clade=="S", select = Virus.name)
v <- subset(tree.detail, Clade=="V", select = Virus.name)
g <- subset(tree.detail, Clade %in% c("G", "GH", "GV", "GR", "GRY", "GK"), select = Virus.name)
gh <- subset(tree.detail, Clade=="GH", select = Virus.name)
gv <- subset(tree.detail, Clade=="GV", select = Virus.name)
gr <- subset(tree.detail, Clade %in% c("GR", "GRY"), select = Virus.name)
gry <- subset(tree.detail, Clade=="GRY", select = Virus.name)
gk <- subset(tree.detail, Clade=="GK", select = Virus.name)

pdf("Worldwide_Skeleton.pdf", width = 25, height = 50)
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
                size = 4, fill = color.list[2] , color= 'grey88', shape = 21) +
  theme_minimal()
dev.off()

tree2 <- groupClade(tree, c("38", "33", "78", "49", "29", "140", "801", "2864", "1340", "3876", "5592"))
pdf("worldwide_Skeleton_ColorEdge.pdf", width = 25, height = 50)
ggtree(tree2, aes(color=group)) +
  # geom_nodelab() +
  scale_color_manual(values=c(color.list[9], #S: 38
                              "grey", #L
                              color.list[9], #S: 78
                              color.list[10], #V: 49
                              "grey", #L
                              color.list[4],  #G: 139
                              color.list[2],  #GK: 801
                              color.list[5],  #GH: 2864
                              color.list[3],  #GV: 1340
                              color.list[6],  #GR: 3876
                              color.list[1],  #GRY: 5592
                              rep("grey99", 10))) 
dev.off()

###################################################################################################

#---Subgroup
subgroup03 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_03.txt", 
                             header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup03 <- subset(fw.detail, gisaid_epi_isl %in% subgroup03$gisaid_epi_isl)
subgroup05 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_05.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup05 <- subset(fw.detail, gisaid_epi_isl %in% subgroup05$gisaid_epi_isl)
subgroup06 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_06.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup06 <- subset(fw.detail, gisaid_epi_isl %in% subgroup06$gisaid_epi_isl)
subgroup07 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_07.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup07 <- subset(fw.detail, gisaid_epi_isl %in% subgroup07$gisaid_epi_isl)
subgroup08 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_08.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup08 <- subset(fw.detail, gisaid_epi_isl %in% subgroup08$gisaid_epi_isl)
subgroup09 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_09.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup09 <- subset(fw.detail, gisaid_epi_isl %in% subgroup09$gisaid_epi_isl)
subgroup12 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_12.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup12 <- subset(fw.detail, gisaid_epi_isl %in% subgroup12$gisaid_epi_isl)
subgroup14 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_14.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup14 <- subset(fw.detail, gisaid_epi_isl %in% subgroup14$gisaid_epi_isl)
treex <- ape::extract.clade(tree, node = c("6841"))
subgroup15 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_15.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup15 <- subset(fw.detail, gisaid_epi_isl %in% subgroup15$gisaid_epi_isl) %>%
  filter(strain %in% treex$tip.label)
# subgroup15 <- subset(fw.detail, gisaid_epi_isl %in% setdiff(subgroup15$gisaid_epi_isl, 
#                                                             c("EPI_ISL_2250065", "EPI_ISL_2789209", "EPI_ISL_1408159", "EPI_ISL_2658767", "EPI_ISL_1999855", "EPI_ISL_2529573",
#                                                               "EPI_ISL_1268731", "EPI_ISL_2652067", "EPI_ISL_2159757", "EPI_ISL_2492229", "EPI_ISL_2267394", "EPI_ISL_1999838",
#                                                               "EPI_ISL_2503153", "EPI_ISL_6050396", "EPI_ISL_2455462", "EPI_ISL_1995372", "EPI_ISL_2116735", "EPI_ISL_1989024",
#                                                               "EPI_ISL_2558196", "EPI_ISL_2259172")))
subgroup16 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_16.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup16 <- subset(fw.detail, gisaid_epi_isl %in% setdiff(subgroup16$gisaid_epi_isl, "EPI_ISL_590899"))
subgroup19 <- read.delim2("D:/00SARS-Cov-2/ADS/worldwide/subgroup/subgroup_19.txt", 
                          header = F, sep = ";", stringsAsFactors = F, col.names = c("gisaid_epi_isl", "colMonth"))
subgroup19 <- subset(fw.detail, gisaid_epi_isl %in% setdiff(subgroup19$gisaid_epi_isl, "EPI_ISL_590899")) %>%
  filter(!(strain %in% c("hCoV-19/Italy/CAM-TIGEM-IZSM-COLLI-15514/2021",
                         "hCoV-19/Switzerland/BL-ETHZ-521740/2021")))

pdf("worldwide_tree.pdf", width = 4, height = 12)
ggtree(tree2, aes(color=group), mrsd = NULL, as.Date = FALSE, size = 0.2) + 
  # geom_nodelab() +
  geom_tippoint(size = 4, fill = adjustcolor('white', alpha.f = 1), color = adjustcolor('gray80', alpha.f = 1), shape = 21) +
  scale_color_manual(values=c(color.list[9], #S: 38
                              "grey", #L
                              color.list[9], #S: 78
                              color.list[10], #V: 49
                              "grey", #L
                              color.list[4],  #G: 139
                              color.list[2],  #GK: 801
                              color.list[5],  #GH: 2864
                              color.list[3],  #GV: 1340
                              color.list[6],  #GR: 3876
                              color.list[1],  #GRY: 5592
                              rep("grey99", 10))) +
  geom_tippoint(aes(subset=((label %in% subgroup03$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[1], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup16$strain)==TRUE)),
                size = 4, fill = "blue", color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup05$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[2], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup06$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[3], color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup07$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[4], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup08$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[5], color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup09$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[6], color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup12$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[7], color= 'white', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup14$strain)==TRUE)),
                size = 4, fill = RColorBrewer::brewer.pal(9, "Pastel1")[8], color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup15$strain)==TRUE)),
                size = 4, fill = "black", color= 'grey80', shape = 21) +
  geom_tippoint(aes(subset=((label %in% subgroup19$strain)==TRUE)),
                size = 4, fill = "#D95F02", color= 'grey80', shape = 21) +
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
