rm(list=ls())
library(tidyverse)
# suppressMessages()
# library(ggplot2)
# library(pheatmap)
# library(NbClust)
# suppressMessages(library(reshape2))
# library(RColorBrewer)
# Sys.setlocale("LC_TIME", "English")

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/Worldwide/Worldwide_hm.RData")

##################################################################################################
### Weighted gene co-expression network analysis
##################################################################################################

suppressMessages(library(WGCNA))
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

input <- hm.data
data <- scale(t(input))

#------------------------------------------------------------------------------------------------
# library(factoextra)
# # Choose a set of soft-thresholding powers
# powers = c(1:50)
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
adjacency = adjacency(data, power = 22, type="signed")
adjacency[1:5, 1:5]

# Call the hierarchical clustering function
mutaTree = hclust(as.dist(1 - adjacency), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(6,6)
plot(mutaTree, labels = colnames(data),
     xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     cex = 0.5);

# set the minimum module size:
minModuleSize = 1;
# Module identification using dynamic tree cut:
mutaTree$height[order(mutaTree$height)]
dynamicMods = cutreeDynamic(dendro = mutaTree, distM = 1 - adjacency, cutHeight = 0.9,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
# Cluster adjustment
dynamicMods[grep("C14805T", colnames(data))] <- as.numeric(dynamicMods[grep("G26144T", colnames(data))])
dynamicMods[dynamicMods==12] <- 10
table(dynamicMods)

# Convert numeric lables into colors
# dynamicColors <- numbers2colors(dynamicMods)
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
# pie(rep(1, 10), col = color.list)
cluster.row.seq <- list()
dynamicColors <- rep(NA, ncol(data))
for (i in 0:length(table(dynamicMods))) {
  
  group = colnames(data)[which(dynamicMods==i)]
  group.info <- word(word(group, 2, 2, fixed(":")), 1, 1, fixed(","))
  # print(list(i, group))
  
  if(length(grep("C22995A", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[2]
    print(list(i, "GK", group))
  }
  if(length(grep("A23063T", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[1]
    print(list(i, "GRY", group))
  }
  if(length(grep("C22227T", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[3]
    print(list(i, "GV", group))
  }
  if(length(grep("A23403G", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[4]
    print(list(i, "G", group))
  }
  if(length(grep("G26144T", group.info))) { 
    dynamicColors[dynamicMods==i] <- color.list[10] 
    print(list(i, "V", group))
  } 
  if(length(grep("T28144C", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[9]
    print(list(i, "S", group))
  }
  if(length(grep("G25563T", group.info))) { 
    dynamicColors[dynamicMods==i] <- color.list[5]
    print(list(i, "GH", group))
  }
  if(length(grep("GGG28881AAC", group.info))) {
    dynamicColors[dynamicMods==i] <- color.list[6]
    print(list(i, "GR", group))
  }
  if(length(grep("C10319T", group.info))) {
    dynamicColors[dynamicMods==i] <- "black"
    print(list(i, "L89F", group))
  }
  if(length(grep("A1163T", group.info))) {
    dynamicColors[dynamicMods==i] <- "blue"
    print(list(i, "S477N", group))
  }

  cluster.row.seq[[i+1]] <- which(dynamicMods==i)
}
dynamicColors[is.na(dynamicColors)] <- "white"
dynamicColors[dynamicMods==3] <- RColorBrewer::brewer.pal(9, "Pastel1")[1]
dynamicColors[dynamicMods==6] <- RColorBrewer::brewer.pal(9, "Pastel1")[2]
dynamicColors[dynamicMods==8] <- RColorBrewer::brewer.pal(9, "Pastel1")[3]
dynamicColors[dynamicMods==10] <- RColorBrewer::brewer.pal(9, "Pastel1")[4]
dynamicColors[dynamicMods==14] <- RColorBrewer::brewer.pal(9, "Pastel1")[5]
dynamicColors[dynamicMods==16] <- RColorBrewer::brewer.pal(9, "Pastel1")[6]

# Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
# pdf("SF4.pdf", width = 14, height = 7)
plotDendroAndColors(mutaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = colnames(data), cex.dendroLabels = 0.6,
                    hang = 0.03,
                    addGuide = TRUE, guideAll = TRUE, guideHang = 0.15,
                    main = "")
# dev.off()

write.csv(data.frame(group = dynamicMods,
                     snps = colnames(data),
                     refpos = as.numeric(str_extract_all(word(word(colnames(data), 2, 2, fixed(":")), 1, 1, fixed(",")), "[0-9]+"))),
          file = "group.csv", quote = F, row.names = F)

##################################################################################################
### topological overlap measure (TOM)
##################################################################################################

suppressMessages(library(igraph))
TOM = TOMsimilarity(adjacency, TOMType="unsigned")
adj <- TOM
adj[adj > 0.001] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj, mode = "undirected")
network <- simplify(network)  # removes self-loops
dev.off()
# pdf("test.pdf")
plot(network, 
     layout = layout_nicely,
     vertex.color = dynamicColors,
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0)

adj <- TOM
adj[adj > 0.00003] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj, mode = "undirected")
network <- simplify(network)  # removes self-loops
dev.off()
loc <- c(which(dynamicMods==0), #GR 
         which(dynamicMods==4), #GV
         which(dynamicMods==1), #GK
         which(dynamicMods==2), #GRY
         which(dynamicMods==9), #G
         which(dynamicMods==15) #GH
         )
plot(induced_subgraph(network, loc), 
     layout=layout_nicely, 
     vertex.label = NA,
     vertex.color = dynamicColors[sort(loc)],
     impl = "create_from_scratch")

##################################################################################################
### Module identification / take example
##################################################################################################

load("./data/all.muta_list.RData")
load("./data/all.meta.RData")

table(dynamicMods)
group.row.id = which(dynamicMods %in% c(1))
if (length(group.row.id) > 1) {
  group.adj <- adjacency[group.row.id, group.row.id]
  pheatmap(group.adj, cluster_rows = F, cluster_cols = F)
  diag(group.adj) <- 0
  degree.cutoff <- 0.9
  group.adj.TF <- (group.adj>=degree.cutoff)
  group.adj.high.id <- names(which(rowSums(group.adj.TF)>=1))
  group.adj.high <- group.adj[group.adj.high.id, group.adj.high.id]
  # group.adj.high
  diag(group.adj.high) <- 1
  pheatmap(group.adj.high, cluster_rows = F, cluster_cols = F)
} else {
  group.adj.high.id <- colnames(data)[group.row.id]
}
group.core.info <- word(word(group.adj.high.id, 2, 2, fixed(":")), 1, 1, fixed(","))
group.core.info

group.sample <- c()
group.nsample <- c()
for (m in group.core.info) {
  if (m==group.core.info[1]) {
    group.sample <- all.muta_list %>%
      filter(grepl(m, muta.list)) %>%
      pull(sample)
  } else {
    group.sample <- intersect(group.sample,
                              all.muta_list %>%
                                filter(grepl(m, muta.list)) %>%
                                pull(sample))
  }
  group.nsample <- c(group.nsample, length(group.sample))
}
plot(group.nsample)
group.nsample
group.nsample / nrow(all.muta_list)

group.meta <- subset(all.meta, Accession.ID %in% group.sample)
table(group.meta$Clade)
table(group.meta$region)
table(group.meta$country)
group.lineage <- data.frame(table(subset(all.meta, Accession.ID %in% group.sample)$Lineage))
group.lineage <- group.lineage[order(group.lineage$Freq, decreasing = T),]
head(group.lineage, 10)

### Downsampling
input <- subset(all.meta, Accession.ID %in% group.sample)
set.seed(2021)
subgroup.downsample <- input[sample(1:nrow(input), 200),]
write.table(subgroup.downsample[,c("Accession.ID", "colMonth")], "subgroup03.txt", quote = F, row.names = F, col.names = F, sep = ";")

subgroup.meta <- subgroup.downsample %>%
  mutate(strain = Virus.name,
         virus = "betacoronavirus", 
         gisaid_epi_isl = Accession.ID,
         genbank_accession = "?",
         date = as.Date(paste(collect.year, collect.month, collect.day, sep = "-")),
         location = word(Location, 4, 4, fixed("/"))) %>%
  select(strain, virus, gisaid_epi_isl, genbank_accession, date, region, country, division, location)
write.table(subgroup.meta, "subgroup03.meta.txt", quote = FALSE, sep = '\t', row.names = F)
