###################################################################################################
### We took India data as an example to show how to perform FTM filtration
###################################################################################################

rm(list=ls())
library(tidyverse)

setwd("D:/00SARS-Cov-2/GibHub")
load("./data/calendar.RData")

###################################################################################################
### Summary of FTMs
###################################################################################################

p.data <- data.frame()
seg.list <- c("5'UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9",
              "NSP10", "NSP12a", "NSP12b", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E",
              "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3'UTR")
for (seg in seg.list) {
  x <- read.delim2(paste0("FTM_India_", seg, ".txt"), header = T, sep = ";")
  x$segment2 <- seg
  p.data <- rbind(p.data, x)
}

p.data <- p.data %>%
  filter(!(mutation %in% c("CAT28253.", "CTA28093.", "TAAATT27134.", "TCAGTT27618.", "TGG26465.") & varclass=="deletion")) %>%
  mutate(anno = paste0(segment2, ":", mutation, ",", ifelse(segment2 %in% c("5'UTR", "3'URT", "Intergenic"), "Extragenic", variant)))
p.data <- p.data %>%
  mutate(col.week = factor(colWeek, labels = format(calendar$Start, "%b %d, %y")),
         ratio = as.numeric(ratio),
         segment2 = factor(segment2, levels = seg.list))

###################################################################################################
### FTM filtreation
###################################################################################################

hm.data <- pivot_wider(data = p.data,
                       id_cols = anno,
                       names_from = colWeek,
                       names_prefix = "Col.",
                       values_from = ratio)
hm.data <- as.data.frame(hm.data)
for (c in 2:ncol(hm.data)) {
  hm.data[,c]  <- as.numeric(hm.data[,c])
}
hm.data2 <- as.matrix(hm.data[,-1])
rownames(hm.data2) <- hm.data$anno
colnames(hm.data2) <- word(colnames(hm.data)[-1], 2, 2, fixed("."))
head(hm.data2)
hm.data <- hm.data2
rm(hm.data2)

# Filtering according to the relative abundance
hm.TF <- (hm.data>10)
hm.TF <- apply(hm.TF, 1, sum)
hm.data <- hm.data[hm.TF>0,]
dim(hm.data)

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

hm.TF <- (hm.count>10)
hm.TF <- apply(hm.TF, 1, sum)
hm.count <- hm.count[hm.TF>0,]
hm.data.cpa <- hm.data[rownames(hm.count),]
dim(hm.count)

library(pheatmap)
# pdf("India_EpidemicMutations.pdf", height = 7, width = 14)
pheatmap(t(hm.data.cpa[,81:1]), cluster_rows = F, legend = F,fontsize = 6)
# dev.off()
