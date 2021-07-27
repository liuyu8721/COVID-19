rm(list=ls())
library(tidyverse)
library(stringi)
library(ggplot2)

setwd("D:/00SARS-Cov-2/ADS")
all.count.India <- read.delim2("all_count_India.txt", header = T, sep = ";")
all.muta.India <- read.delim2("all_muta_India.txt", header = T, sep = ";")

###################################################################################################

dummy <- data.frame(Var1 = c(paste("2020", stri_reverse(substr(stri_reverse(paste0("0", 1:52)), 1, 2)), sep = "-"),
                             paste("2021", stri_reverse(substr(stri_reverse(paste0("0", 1:52)), 1, 2)), sep = "-")[1:8]))
nsample <- as.data.frame(table(all.count.India$col.week))
nsample <- merge(dummy, nsample, all.x = T, sort = T)
nsample[is.na(nsample)] <- 0.1
names(nsample) <- c("col.week", "n")
nsample
sum(nsample$n)

###################################################################################################
### Frequency Trajectories of Mutations by segment
###################################################################################################
muta.tab <- table(all.muta.India$segment2)
muta.tab

## local mutation
for (seg in names(muta.tab)[29]) {
  S <- subset(all.muta.India, segment2==seg)
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
      
      muta.y <- data.frame(xtabs(~ col.week, y))
      if (nrow(muta.y)>0) {
        # names(muta.y)[1] <- "col.week"
        muta.y <- merge(nsample, muta.y, all.x = T, all.y = T) %>%
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
  
  mf.ts <- ggplot(data = muta.detail, aes(x = col.week, y = ratio, group = mutation)) +
    geom_line(aes(color = mutation)) +
    # scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(x = "", y = "Mutational frequency") +
    ggtitle(paste(seg, "mutations over time")) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10, hjust = 1),
          axis.text.x = element_text(angle = 30),
          axis.title = element_text(size = 14),
          panel.grid = element_line(colour = "grey88"),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))

  print(mf.ts)
  write.table(muta.detail, paste0("FTM_India_", seg, ".txt"), row.names = F, quote = F, sep = ";")
}
