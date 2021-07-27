rm(list=ls())
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
Sys.setlocale("LC_TIME", "English")

###################################################################################################
### FTM filtreation
###################################################################################################

load("./RData/hm.RData")
hm.TF <- (hm.data>0.1)
hm.TF <- apply(hm.TF, 1, sum)
hm.data <- hm.data[hm.TF>0,]
dim(hm.data)

hm.eucl <- dist(hm.data, method = "euclidean")
hc <- hclust(hm.eucl, method = "ward.D")

### Selecting an appropriate cluster number for dendrogram cutting.
library(factoextra)
set.seed(1234)
fviz_nbclust(hm.data, FUNcluster = hcut, method = "gap_stat", diss = hm.eucl, k.max = 10, nboot = 100) +
  geom_vline(xintercept = c(5, 8), linetype = 2) +
  theme(axis.text.x = element_text(angle = 90))

plot(hc, labels = FALSE)
re <- rect.hclust(hc, k = 5)
re.hm <- data.frame()
for (i in 1:length(re)) {
        re.hm <- rbind(re.hm,
                       data.frame(d = i,
                                  n = length(re[[i]])))
}
re.hm <- re.hm[order(re.hm$n, decreasing = T),]
re.hm

re.hm.rest <- re.hm[-1,]
r.lab <- c()
for (i in 1:nrow(re.hm.rest)) {
        r.lab <- c(r.lab, re[[re.hm.rest$d[i]]])
}

hm.wna <- hm.data[names(r.lab),]
pheatmap(t(hm.wna[,60:1]), cluster_rows = F)

