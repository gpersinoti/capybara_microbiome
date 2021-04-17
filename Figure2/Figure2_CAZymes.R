#' Update Capybara paper using the annotation from CAZy database
#' Date: July, 22nd, 2020

library("DescTools")
library("ggplot2")
library("cowplot")
library("reshape2")
library("plyr")
library("data.table")

setwd("/data/Projects/Capibara/github/Figure2/")


cazyFamiliesInterest <- read.delim(file="CAZyFamiliesCapybara.csv")
head(cazyFamiliesInterest)

dataMG <- read.delim(file="MG.TPM.matrix", stringsAsFactors = F)
names(dataMG)[1] <- "GeneID"
dataMG$Cecal <- rowMeans(dataMG[2:4])
dataMG$Rectal <- rowMeans(dataMG[5:7])
dataMG <- dataMG[,c(1,8:9)]
dataMG$omics <- "MG"
head(dataMG)

#' Prepare MT TPM data
dataMT <- read.delim(file="MT.TPM.matrix",stringsAsFactors = F)
names(dataMT)[1] <- "GeneID"
dataMT$Cecal <- rowMeans(dataMT[2:4])
dataMT$Rectal <- rowMeans(dataMT[5:7])
dataMT <- dataMT[,c(1,8:9)]
dataMT$omics <- "MT"
head(dataMT)


dataMG_MT <- rbind(dataMG,dataMT)

head(dataMG_MT)

dataCAZy <- merge(data2,dataMG_MT, all.x=T, by="GeneID")
head(dataCAZy)


#'Figure by Group of activity
head(dataCAZy)
dataCAZy <- merge(dataCAZy, cazyFamiliesInterest,by="CAZy")
dataCAZy$CAZy= factor(dataCAZy$CAZy, levels=unique(gtools::mixedsort(as.character(dataCAZy$CAZy))))

ord <- c("plant cell wall glycans", "sucrose", "starch & glycogen", "microbial glycans", "animal glycans",
         "others")
dataCAZy$Activity <- factor(dataCAZy$Activity, levels = ord)

head(dataCAZy)
dataCAZy <- reshape::melt(dataCAZy)

p <- ggplot(dataCAZy[which(dataCAZy$CAZy %in% x$CAZy),], aes(x=CAZy, fill=variable))+
  geom_bar(aes(weight=value), position="dodge") + 
  facet_grid(Activity ~ omics, scales = "free", space = "free", labeller = label_wrap_gen(width=15)) + 
  xlab("CAZy families") + 
  ylab("Abundance") +
  coord_flip() +  scale_fill_brewer(palette = "Set1") +
  background_grid(major = "xy", minor = "xy") + 
  theme_cowplot(12) + background_grid(major = "x", minor = "") +
  theme(legend.title=element_blank())
p

save_plot("Figure2.pdf",p,base_width = 8, base_height = 9)
