###################################################################################################
### We took India data as an example to show how to perform Weighted Network modeling of FTMs
###################################################################################################

rm(list=ls())
library(tidyverse)
library(stringi)
library(ggplot2)
library(RColorBrewer)
Sys.setlocale("LC_TIME", "English")

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/India/hm.data.cpa.RData")
load("./data/India/calendar.RData")

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

input <- hm.data.cpa
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

#------------------------------------------------------------------------------------------------
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
#      labels=powers,cex=cex1,col="red")
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#------------------------------------------------------------------------------------------------

# calculate the adjacencies, using the soft thresholding power
adjacency = adjacency(data, power = 24, type="signed")
adjacency[1:5, 1:5]
pheatmap(adjacency, show_rownames = F)

# Call the hierarchical clustering function
mutaTree = hclust(as.dist(1 - adjacency), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(6,6)
plot(mutaTree, labels = colnames(data),
     xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     cex = 0.5)

# set the minimum module size:
minModuleSize = 1;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = mutaTree, distM = 1 - adjacency, 
                            # cutHeight = 0.999,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

### Module reassignment for the lanscape but no influence of module-based variant identification
# dynamicMods[grep("C23604A", colnames(data))] <- 16
# dynamicMods[grep("C15240T|C21846T|C2453T|T8620.|ATT21993.|G23012A|G23012C|T15096C|C26681T|G21987A|C8683T|G23593T", colnames(data))] <- 0

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
dynamicColors[dynamicMods==0] <- "white"
color.list <- c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))
dynamicColors[dynamicMods==1] <- RColorBrewer::brewer.pal(8, "Set2")[1]         #GRY
dynamicColors[dynamicMods==2] <- RColorBrewer::brewer.pal(8, "Set2")[2]         #GK
dynamicColors[dynamicMods==3] <- RColorBrewer::brewer.pal(8, "Set2")[5]
dynamicColors[dynamicMods==4] <- RColorBrewer::brewer.pal(8, "Set2")[6]
dynamicColors[dynamicMods==5] <- RColorBrewer::brewer.pal(8, "Set2")[7]
dynamicColors[dynamicMods==6] <- RColorBrewer::brewer.pal(8, "Set2")[8]
dynamicColors[dynamicMods==7] <- RColorBrewer::brewer.pal(12, "Paired")[1]
dynamicColors[dynamicMods==8] <- RColorBrewer::brewer.pal(12, "Paired")[2]
dynamicColors[dynamicMods==9] <- RColorBrewer::brewer.pal(12, "Paired")[3]
dynamicColors[dynamicMods==10] <- RColorBrewer::brewer.pal(12, "Paired")[4]
dynamicColors[dynamicMods==11] <- RColorBrewer::brewer.pal(8, "Set2")[4]        #G
dynamicColors[dynamicMods==12] <- RColorBrewer::brewer.pal(12, "Paired")[5]
dynamicColors[dynamicMods==13] <- RColorBrewer::brewer.pal(12, "Paired")[6]
dynamicColors[dynamicMods==14] <- RColorBrewer::brewer.pal(12, "Paired")[7]
dynamicColors[dynamicMods==15] <- RColorBrewer::brewer.pal(12, "Paired")[8]
dynamicColors[dynamicMods==16] <- RColorBrewer::brewer.pal(8, "Accent")[3]
dynamicColors[dynamicMods==17] <- RColorBrewer::brewer.pal(8, "Accent")[2]      #GR
dynamicColors[dynamicMods==18] <- RColorBrewer::brewer.pal(12, "Paired")[10]
dynamicColors[dynamicMods==19] <- RColorBrewer::brewer.pal(8, "Accent")[4]
dynamicColors[dynamicMods==20] <- RColorBrewer::brewer.pal(12, "Paired")[12]
dynamicColors[dynamicMods==21] <- RColorBrewer::brewer.pal(12, "Set3")[1]
dynamicColors[dynamicMods==22] <- RColorBrewer::brewer.pal(12, "Set3")[2]
dynamicColors[dynamicMods==23] <- RColorBrewer::brewer.pal(12, "Set3")[3]
dynamicColors[dynamicMods==24] <- RColorBrewer::brewer.pal(12, "Set3")[4]
dynamicColors[dynamicMods==25] <- RColorBrewer::brewer.pal(12, "Set3")[5]
dynamicColors[dynamicMods==26] <- RColorBrewer::brewer.pal(12, "Set3")[6]
dynamicColors[dynamicMods==27] <- RColorBrewer::brewer.pal(8, "Accent")[5]
dynamicColors[dynamicMods==28] <- RColorBrewer::brewer.pal(12, "Set3")[8]
dynamicColors[dynamicMods==29] <- RColorBrewer::brewer.pal(12, "Set3")[9]
dynamicColors[dynamicMods==30] <- RColorBrewer::brewer.pal(12, "Set3")[10]
dynamicColors[dynamicMods==31] <- RColorBrewer::brewer.pal(12, "Set3")[11]
dynamicColors[dynamicMods==32] <- RColorBrewer::brewer.pal(12, "Set3")[12]

# Plot the dendrogram and colors underneath
# pdf("India_WGCNA.pdf", width = 14, height = 5)
plotDendroAndColors(mutaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = colnames(data), cex.dendroLabels = 0.4,
                    hang = 0.02,
                    addGuide = TRUE, guideAll = TRUE, guideHang = 0.15,
                    main = "")
# dev.off()

write.csv(data.frame(group = dynamicMods,
                     snps = colnames(data),
                     color = dynamicColors,
                     refpos = as.numeric(str_extract_all(word(word(colnames(data), 2, 2, fixed(":")), 1, 1, fixed(",")), "[0-9]+"))),
          file = "./data/India/group.csv", quote = F, row.names = F)

##################################################################################################
### Module identification
##################################################################################################

load("./data/India/India.muta.list.RData")
load("./data/India/India.meta.RData")

all.muta_list <- India.muta.list
rm(India.muta.list)
all.meta <- India.meta
rm(India.meta)

nsample <- calendar[,c("colWeek", "n")]
nsample[nsample$n==0, "n"] <- 0.1

module.monitor <- data.frame(colWeek = calendar[,"colWeek"])
for (group.no in setdiff(as.numeric(names(table(dynamicMods))), 0)) {
  
  ### Extract module feature
  group.row.id = which(dynamicMods %in% group.no)
  if (length(group.row.id) > 1) {
    group.adj <- adjacency[group.row.id, group.row.id]
    # pheatmap(group.adj, cluster_rows = F, cluster_cols = F)
    diag(group.adj) <- 0
    degree.cutoff <- 0.7
    group.adj.TF <- (group.adj>=degree.cutoff)
    group.adj.high.id <- names(which(rowSums(group.adj.TF)>=1))
    group.adj.high <- group.adj[group.adj.high.id, group.adj.high.id]
    diag(group.adj.high) <- 1
    # pheatmap(group.adj.high, cluster_rows = F, cluster_cols = F)
  } 
  if (length(group.adj.high.id)) {
    group.core.info <- word(word(group.adj.high.id, 2, 2, fixed(":")), 1, 1, fixed(","))
    # print(group.core.info)
    
    ### Extract genomes with module feature
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
    
    ### Generate module prevalence
    group.meta <- subset(all.meta, Accession.ID %in% group.sample)
    group.prev <- data.frame(table(group.meta$colWeek))
    names(group.prev) <- c("colWeek", paste0("module", group.no))
    group.prev <- merge(calendar[, c("colWeek", "n")], group.prev, all.x = T)
    group.prev[is.na(group.prev)] <- 0
    group.prev[,3] <- group.prev[,3] / nsample$n
    
    module.monitor <- merge(module.monitor, group.prev[,c(1, 3)])
  }
}

mat <- t(module.monitor[1:99, -1]) * 100
ts.label <- rep("", nrow(calendar))
ts.label[seq(1, 99, 2)] <- format(calendar$Start[seq(1, 99, 2)], "%b %d, %y")
pheatmap(mat, 
         cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(colors = c("white","red"))(100),
         labels_col = as.character(ts.label),
         border_color = "grey88",
         angle_col = "45")
