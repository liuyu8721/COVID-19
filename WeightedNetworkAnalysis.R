rm(list=ls())
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

##################################################################################################
### Weighted network analysis
##################################################################################################

# Load the WGCNA package
suppressMessages(library(WGCNA))
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary. Any error here may be ignored but you may want to update WGCNA if you see one. Caution: skip this line if you run RStudio or other third-party R environments.
enableWGCNAThreads()

load("D:/00SARS-Cov-2/ADS/RData/hm.wna.RData")
input <- hm.wna
data <- scale(t(input))

#------------------------------------------------------------------------------------------------
# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
# # Plot the results:
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#------------------------------------------------------------------------------------------------

# calculate the adjacencies, using the soft thresholding power
adjacency = adjacency(data, power = 17, type="signed")
# pheatmap(adjacency, show_rownames = F)

# Call the hierarchical clustering function
mutaTree = hclust(as.dist(1 - adjacency), method = "average");

# set the minimum module size:
minModuleSize = 1;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = mutaTree, distM = 1 - adjacency, cutHeight = 0.94,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
color.list <- rev(c(RColorBrewer::brewer.pal(3, "Set2"),  RColorBrewer::brewer.pal(6, "Accent")))
cluster.row.seq <- list()
dynamicColors <- rep(NA, ncol(data))
for (i in 0:length(table(dynamicMods))) {
  
  group = colnames(data)[which(dynamicMods==i)]
  group.info <- word(word(group, 2, 2, fixed(":")), 1, 1, fixed(","))
  # print(list(i, group))
  
  if(length(grep("A23063T", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[9]
    print(list(i, "GRY", group))
  }
  if(length(grep("C22227T", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[8] #GV
    print(list(i, "GV", group))
  }
  if(length(grep("C21855T", group.info))) { 
    dynamicColors[dynamicMods==i] <- "yellow" #S98F
    print(list(i, "S98F", group))
  }
  if(length(grep("C10319T", group.info))) {
    dynamicColors[dynamicMods==i] <- "black" #A.E
    print(list(i, "S98F", group))
  }
  if(length(grep("A1163T", group.info))) {
    dynamicColors[dynamicMods==i] <- "blue"
    print(list(i, "S477N", group))
  }
  if(length(grep("A23403G", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[7]
    print(list(i, "G", group))
  }
  if(length(grep("G26144T", group.info))) { 
    dynamicColors[dynamicMods==i] <- color.list[1] # V
    print(list(i, "V", group))
  } 
  if(length(grep("T28144C", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[2]
    print(list(i, "S", group))
  }
  if(length(grep("G25563T", group.info))) { 
    dynamicColors[dynamicMods==i] <- color.list[6]; # GH
    print(list(i, "GH", group))
  }
  if(length(grep("GGG28881AAC", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[5] # GR
    print(list(i, "GR", group))
  }
  
  cluster.row.seq[[i+1]] <- which(dynamicMods==i)
}
dynamicColors[dynamicMods==2] <- RColorBrewer::brewer.pal(9, "Pastel1")[1]
dynamicColors[dynamicMods==7] <- RColorBrewer::brewer.pal(9, "Pastel1")[2]
dynamicColors[dynamicMods==8] <- RColorBrewer::brewer.pal(9, "Pastel1")[3]
dynamicColors[dynamicMods==10] <- RColorBrewer::brewer.pal(9, "Pastel1")[4]
dynamicColors[dynamicMods==12] <- RColorBrewer::brewer.pal(9, "Pastel1")[5]
dynamicColors[dynamicMods==13] <- RColorBrewer::brewer.pal(9, "Pastel1")[6]
dynamicColors[dynamicMods==14] <- RColorBrewer::brewer.pal(9, "Pastel1")[7]
dynamicColors[dynamicMods==16] <- RColorBrewer::brewer.pal(9, "Pastel1")[8]

# Plot the dendrogram and colors underneath (Supplemental Figure 3)
plotDendroAndColors(mutaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = colnames(data), cex.dendroLabels = 0.6,
                    hang = 0.03,
                    addGuide = TRUE, guideAll = TRUE, guideHang = 0.15,
                    main = "")

##################################################################################################
### Visualization using igraph
##################################################################################################

# Load the igraph package
suppressMessages(library(igraph))
TOM = TOMsimilarity(adjacency, TOMType="unsigned")
adj <- TOM
# Change the cutoff for visualization
adj[adj > 0.09] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj, mode = "undirected")
network <- simplify(network)  # removes self-loops
# Plot Figure 3A or Supplemental Figure 4
plot(network, 
     layout = layout_nicely,
     vertex.color = dynamicColors,
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0)
# TOM[dynamicMods==5,dynamicMods==9]

# Plot Figure 3B
loc <- c(which(dynamicMods==0), which(dynamicMods==1), which(dynamicMods==9), which(dynamicMods==17))
plot(induced_subgraph(network, loc), 
     vertex.label = NA,
     vertex.color = dynamicColors[sort(loc)],
     impl = "create_from_scratch")
