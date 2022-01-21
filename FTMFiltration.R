###################################################################################################
### We took India data as an example to show how to perform FTM filtration
###################################################################################################

rm(list=ls())
library(tidyverse)

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/India/India.p.data.RData")
load("./data/India/India.hm.data.RData")
load("./data/India/calendar.RData")

###################################################################################################
### FTM filtreation
###################################################################################################

hm.data <- India.hm.data

### Filtering according to the relative abundance
hm.TFr <- (hm.data>10)
dim(hm.TFr)

### Filtering according to the absolute amount of case
hm.count <- pivot_wider(data = p.data,
                        id_cols = anno,
                        names_from = colWeek,
                        names_prefix = "Col.",
                        values_from = Freq)
hm.count2 <- as.matrix(hm.count[,-1])
rownames(hm.count2) <- hm.count$anno
colnames(hm.count2) <- word(colnames(hm.count)[-1], 2, 2, fixed("."))
head(hm.count2)
hm.count <- hm.count2
rm(hm.count2)
hm.count <- hm.count[rownames(hm.data),]
hm.TFc <- (hm.count>10)
dim(hm.TFc)
hm.TFc <- hm.TFc[row.names(hm.TFr),]

hm.filter <- hm.TFr + hm.TFc
hm.filter <- (hm.filter == 2)
hm.TF <- apply(hm.filter, 1, sum)
hm.TF <- hm.TF[hm.TF>0]
hm.data.cpa <- hm.data[names(hm.TF),]
dim(hm.data.cpa)
save(hm.data.cpa, file = "./data/India/hm.data.cpa.RData")

library(pheatmap)
ts.label <- rep("", nrow(calendar))
ts.label[seq(1, 99, 2)] <- format(calendar$Start[seq(1, 99, 2)], "%b %d, %y")
# pdf("India_EpidemicMutations.pdf", height = 7, width = 14)
pheatmap(t(hm.data.cpa[,99:1]), cluster_rows = F, legend = F,
         labels_row = ts.label[99:1], fontsize = 4)
# dev.off()
