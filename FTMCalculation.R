###################################################################################################
### We took India data as an example to show how to calculate frequency trajectories of mutations
###################################################################################################

rm(list=ls())
library(tidyverse)
Sys.setlocale("LC_TIME", "English")

setwd("D:/00SARS-Cov-2/GibHub")

###################################################################################################
### Data praperation
###################################################################################################

# load("all.meta.RData")
# load("all.muta.RData")
# load("all.muta.list.RData")
# load("calendar.RData")
# 
# India.meta <- all.meta %>%
#   filter(colWeek %in% calendar$colWeek & country=="India")
# India.muta <- subset(all.muta, sample %in% India.meta$Accession.ID)
# India.muta.list <- subset(all.muta.list, sample %in% India.meta$Accession.ID)
# 
# dim(India.meta)
# dim(India.muta)
# length(unique(India.muta$sample))
# 
# nsample <- data.frame(table(India.meta$colWeek))
# names(nsample) <- c("colWeek", "n")
# calendar <- merge(calendar[, c("colWeek", "Start", "End")], nsample, all.x = T)
# calendar[is.na(calendar)] <- 0
# 
# save(India.meta, file = "./data/India/India.meta.RData")
# save(India.muta, file = "./data/India/India.muta.RData")
# save(India.muta.list, file = "./data/India/India.muta.list.RData")
# save(calendar, file = "./data/India/calendar.RData")

###################################################################################################
### Sample size
###################################################################################################

load("./data/India/India.meta.RData")
load("./data/India/India.muta.RData")
load("./data/India/calendar.RData")

dim(India.meta)
dim(India.muta)
length(unique(India.muta$sample))

### Rename 
all.meta <- India.meta
all.muta <- India.muta

### Sample size by sampling week
nsample <- data.frame(table(all.meta$colWeek))
names(nsample) <- c("colWeek", "n")
nsample <- merge(calendar[,c("colWeek", "Start")], nsample, by = "colWeek", all.x = T, all.y = F) %>%
  mutate(n = ifelse(is.na(n), 0.1, n))
nsample

ts.label <- rep("", nrow(calendar))
ts.label[seq(1, 99, 2)] <- format(calendar$Start[seq(1, 99, 2)], "%b %d, %y")

# pdf("India_SampleDistribution.pdf", width = 14, height = 2)
ggplot(nsample, aes(x = colWeek, y = n, group = 1)) +
  geom_line(size = 1.5, col = "#386CB0") +  
  geom_point(pch = 16, size = 2, col = "#FC8D62") +
  geom_point(pch = 1, size = 2, col = "#386CB0") +
  scale_x_discrete(labels = ts.label) +
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 500)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
  labs(x = "", y = "") +
  # guides(colour = guide_legend(nrow = 1)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = "grey88"),
        legend.position = c(0.45, 0.9),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(fill = NA),
        # axis.title = element_text(size = 14),
        # axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# dev.off()


###################################################################################################
### Epidemic strains
###################################################################################################

wc <- as.data.frame.matrix(xtabs(~ colWeek + factor(Clade), data = all.meta))
wc <- wc %>%
  mutate(colWeek = row.names(wc))
wc <- merge(nsample, wc, all.x = T, all.y = F)
wc[is.na(wc)] <- 0

library(reshape2)
wc <- melt(wc, id.vars = c("colWeek", "Start", "n")) %>%
  mutate(prop = value / n * 100,
         Clade = factor(variable, levels = rev(c("V", "S", "O", "L", "GR", "GH", "G", "GV", "GK", "GRY"))))

# pdf("India_EpidemicStrains.pdf", width = 10, height = 6.8)
ggplot(data = wc, aes(x = colWeek, fill = Clade)) +
  geom_bar(aes(weight = prop),position="stack") +
  scale_fill_discrete(type = c(RColorBrewer::brewer.pal(4, "Set2"),  RColorBrewer::brewer.pal(6, "Accent"))) +
  scale_x_discrete(labels = ts.label) +
  # scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = c("0", "25", "50", "75", "100")) +
  guides(fill = guide_legend(nrow = 1, reverse = T)) +
  labs(x = "Sampling week", y = "Proportion of epidemic strain (%)", fill = "Clade") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major.x = element_line(colour = "grey88"),
        legend.position = "top")
# dev.off()


###################################################################################################
### Frequency Trajectories of Mutations by segment
###################################################################################################

load("./data/India/muta.anno.RData")
all.muta <- merge(India.muta, 
                         muta.anno[, c("mutation", "refpos", "protein", "variant", "varclass", "annotation")], 
                         all = F)
all.muta <- all.muta %>%
  rename(segment2 = protein) %>%
  mutate(mutation =  paste0(refvar, refpos, qvar))
temp <- subset(all.meta, Accession.ID %in% unique(all.muta$sample), select = c("Accession.ID", "colWeek"))
all.muta <- merge(all.muta, temp, 
                  by.x = "sample", by.y = "Accession.ID", all = F)
rm(temp)

muta.tab <- table(all.muta$segment2)
muta.tab

for (seg in names(muta.tab)) {
  S <- subset(all.muta, segment2==seg)
  S.pos <- as.data.frame(table(S$refpos))
  S.pos <- S.pos[order(S.pos$Freq, decreasing = T),]
  head(S.pos, 20)

  muta.detail <- data.frame()
  y.S <- S
  for (p in S.pos$Var1) {
    x <- subset(y.S, refpos==p)

    muta.c <- as.data.frame(xtabs(~variant+varclass, data = x))
    muta.c <- muta.c[order(muta.c$Freq, decreasing = T),]
    muta.c$percent <- with(muta.c, Freq / sum(Freq) * 100)
    muta.c <- muta.c %>%
      filter(percent > 0)
    for (k in 1:nrow(muta.c)) {
      y <- subset(x, variant==muta.c[k,"variant"] & varclass==muta.c[k,"varclass"])
      muta <- data.frame(table(with(y, paste0(refvar, refpos, qvar))))
      muta <- as.character(muta[order(muta$Freq, decreasing = T),1][1])

      muta.y <- data.frame(xtabs(~ colWeek, y))
      if (nrow(muta.y)>0) {
        # names(muta.y)[1] <- "colWeek"
        muta.y <- merge(nsample, muta.y, all.x = T, all.y = F) %>%
          mutate(Freq = ifelse(is.na(Freq), 0, Freq))
        muta.y$ratio <- with(muta.y, Freq / n * 100)

        muta.detail <- rbind(muta.detail,
                             data.frame(segment = y$segment[1],
                                        annotation = y$annotation[1],
                                        status = muta.c[k, "Freq"] / sum(muta.c$Freq) * 100,
                                        position = p,
                                        mutation = muta,
                                        variant = muta.c[k,"variant"],
                                        varclass = muta.c[k,"varclass"],
                                        muta.y))
      }
    }
    y.S <- subset(y.S, refpos!=p)
  }

  mf.ts <- ggplot(data = muta.detail, aes(x = colWeek, y = ratio, group = mutation)) +
    geom_line(aes(color = mutation)) +
    # scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(x = "", y = "Mutational frequency") +
    ggtitle(paste(seg, "mutations over time")) +
    theme(legend.position = "none",
          axis.text = element_text(size = 9, hjust = 1),
          axis.text.x = element_text(angle = 30),
          axis.title = element_text(size = 14),
          panel.grid = element_line(colour = "grey88"),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))

  pdf(paste0("FTM_India_", seg, ".pdf"), width = 14, height = 7)
  print(mf.ts)
  dev.off()
  write.table(muta.detail, paste0("./data/India/FTM/FTM_India_", seg, ".txt"), row.names = F, quote = F, sep = ";")
}

###################################################################################################
### Summary of FTMs
###################################################################################################

p.data <- data.frame()
seg.list <- c("5'UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9",
              "NSP10", "NSP12a", "NSP12b", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E",
              "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3'UTR")
for (seg in seg.list) {
  x <- read.delim2(paste0("./data/India/FTM/FTM_India_", seg, ".txt"), header = T, sep = ";")
  x$segment2 <- seg
  p.data <- rbind(p.data, x)
}

p.data <- subset(p.data, !(segment2=="3'UTR" & position <= 29674))
p.data <- p.data %>%
  mutate(anno = paste0(segment2, ":", mutation, ",", ifelse(segment2 %in% c("5'UTR", "3'UTR"), "Extragenic", variant))) %>%
  filter(mutation!="A28273.")
save(p.data, file = "./data/India/India.p.data.RData")

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
India.hm.data <- hm.data
save(India.hm.data, file = "./data/India/India.hm.data.RData")
