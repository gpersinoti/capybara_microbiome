#' ---
#' title: "Create the Taxonomy Figure Capybara Microbiome"
#' author: "Gabriela F. Persinoti"
#' date: "March 3rd, 2020"
#' output: html_document
#' ---

library(cowplot)
library(plyr)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(genefilter)
library(ggplot2)
library(corrplot)

theme_set(theme_bw())

setwd("/data/Projects/Capibara/github/Figure1//")

#' Taxons that will be kept
taxon <- c("Actinobacteria","Ascomycota", "Bacteroidetes","Euryarchaeota","Fibrobacteres","Firmicutes",
           "Fusobacteria", "Proteobacteria","Spirochaetes","Tenericutes","unclassified")

colors <- c("#E41A1C", "#459D71","#566B9B","#6F8272","#FFF32E","#B45B76",
            "#FFA000","#EDBB99","#B6E3BF","#F3F2BA","#BEBADA","#D5DBDB")

#'MG data merging the results from Kaiju at phylum level
filesMG <- list.files(path="MG/",pattern = ".phylum.summary")
filesMG
data <- read.delim(file=paste0("MG/",filesMG[1]),comment.char = "-",stringsAsFactors = F)

for (i in 2: length(filesMG)){
  data2 <- read.delim(file=paste0("MG/",filesMG[i]),comment.char = "-")
  data <- rbind(data, data2)
}


#' Add variables
data$animal <- gsub("_1_S[0-9]_L00[12].out","",data$file)
data$animal <- gsub("_HKTFVBCX2_L00[12].out","",data$animal)
data$part <- gsub("_[1234]","",data$animal)
head(data)

#' Mean between the two lanes for each animal
data_animals <- ddply(data, c("taxon_name","animal"), summarize, percent=mean(percent))
data_animals[!(data_animals$taxon_name %in% taxon),'taxon_name'] <- "Others"

dataMG <- data_animals


#'########################################################################
#'
#'MT data merging the results from Kaiju at phylum level
filesMT <- list.files(path="MT/",pattern = ".phylum.summary")
data <- read.delim(file=paste0("MT/",filesMT[1]),comment.char = "-",stringsAsFactors = F)

for (i in 2: length(filesMT)){
  data2 <- read.delim(file=paste0("MT/",filesMT[i]),comment.char = "-")
  data <- rbind(data, data2)
}

#' Add variables
data$animal <- gsub("_1_L00[12].nonrRNA.R1.fastq.gz.out","",data$file)
data$part <- gsub("_[1234]","",data$animal)


#' Mean between the two lanes for each animal
data_animalsMT <- ddply(data, c("taxon_name","animal"), summarize, percent=mean(percent))
#data_animalsMT[data_animalsMT$percent < 0.5,'taxon_name'] <- "Others"
data_animalsMT[!(data_animalsMT$taxon_name %in% taxon),'taxon_name'] <- "Others"
data_animalsMT[data_animalsMT$taxon_name =='cannot be assigned to a (non','taxon_name'] <- "Others"
dataMT <- data_animalsMT

#'#######################################################################################################

data <- read.delim(file="16S_Phylum_Relative_abundance.txt", sep=",",stringsAsFactors = F)
head(data)
names(data)[1] <- "taxon_name"
data <- melt(data)
head(data)

data$variable <- gsub("X","",data$variable)
data$taxon_name <- gsub("p:","",data$taxon_name)
data[!(data$taxon_name %in% taxon),'taxon_name'] <- "Others"

names(data) <- c('taxon_name','animal','percent')
data16S <- data


#'#######################################################################################################

# Taxonomy based on 16S reads recovered from MG using Reago. 
#' 16S_MG

#'16S MG reads merging the results from KaijuMG at phylum level
files16SMG <- list.files(path="16S_MG/",pattern = ".phylum.summary")
files16SMG
data <- read.delim(file=paste0("16S_MG/",files16SMG[1]),comment.char = "-",stringsAsFactors = F)

for (i in 2: length(files16SMG)){
  data2 <- read.delim(file=paste0("16S_MG/",files16SMG[i]),comment.char = "-")
  data <- rbind(data, data2)
}

#' Add variables
data$animal <- gsub("_1_S[0-9]_L00[12].out","",data$file)
data$animal <- gsub("kaijuMG/","",data$animal)
data$animal <- gsub("_HKTFVBCX2_L00[12].out","",data$animal)

#' Mean between the two lanes for each animal
data_animals <- ddply(data, c("taxon_name","animal"), summarize, percent=mean(percent))
data_animals[!(data_animals$taxon_name %in% taxon),'taxon_name'] <- "Others"
data_animals[data_animals$taxon_name =='belong to a (non','taxon_name'] <- "Others"

data16SMG <- data_animals


#'###########################################
#' Merge the graphs
#' 
head(dataMG)
head(dataMT)
head(data16S)
head(data16SMG)

dataMG$animal = gsub("AC","Cecum",dataMG$animal)
dataMG$animal = gsub("R","Recto",dataMG$animal)

data16S$animal = gsub("AC","Cecum",data16S$animal)
data16S$animal = gsub("R","Recto",data16S$animal)

data16S$animal = gsub("2Cecum","Cecum_2",data16S$animal)
data16S$animal = gsub("3Cecum","Cecum_3",data16S$animal)
data16S$animal = gsub("4Cecum","Cecum_4",data16S$animal)
data16S$animal = gsub("2Recto","Recto_2",data16S$animal)
data16S$animal = gsub("3Recto","Recto_3",data16S$animal)
data16S$animal = gsub("4Recto","Recto_4",data16S$animal)

dataMT$animal = gsub("Ceco","Cecum",dataMT$animal)
dataMT$animal = gsub("Reto","Recto",dataMT$animal)

data16SMG$animal = gsub("AC","Cecum",data16SMG$animal)
data16SMG$animal = gsub("R","Recto",data16SMG$animal)


data16S$omics="16S"
dataMT$omics = "MT"
data16SMG$omics="16S_MG"
dataMG$omics = "MG"

#' The main plot
dataPlot <- rbind(dataMG,dataMT,data16S, data16SMG)
dataPlot$taxon_name = gsub("p:","",dataPlot$taxon_name)


dataPlot$part <- gsub("_[2,3,4]","",dataPlot$animal)
data_parts <- ddply(dataPlot, c("taxon_name","omics","part"), summarize, abundance=mean(percent), sd=sd(percent,na.rm = T))


p.A <- ggplot(data_parts, aes(x = omics, y = abundance, fill = taxon_name)) +
  geom_bar(stat = 'identity',alpha=0.8, position = "fill") + scale_fill_manual(values=colors) +
  xlab("") + ylab("Relative abundance (%)") + ggtitle("") +
  facet_grid(. ~ part) + background_grid(major = "xy", minor = "xy") +
  theme(legend.title=element_blank(), axis.text = element_text(colour = "black")) 
p.A


################################################################
dataRatio <- rbind(dataMG,dataMT)
dataRatio <- dcast(dataRatio, taxon_name ~ omics + animal , sum, value.var="percent")
head(dataRatio)

dataRatio$c1 = (dataRatio$MT_Cecum_2+1)/(dataRatio$MG_Cecum_2+1)
dataRatio$c2 = (dataRatio$MT_Cecum_3+1)/(dataRatio$MG_Cecum_3+1)
dataRatio$c3 = (dataRatio$MT_Cecum_4+1)/(dataRatio$MG_Cecum_4+1)

dataRatio$r1 = (dataRatio$MT_Recto_2+1)/(dataRatio$MG_Recto_2+1)
dataRatio$r2 = (dataRatio$MT_Recto_3+1)/(dataRatio$MG_Recto_3+1)
dataRatio$r3 = (dataRatio$MT_Recto_4+1)/(dataRatio$MG_Recto_4+1)

head(dataRatio)
dataRatio2 <- dataRatio[,c(1,14:19)]
dataRatio2

#'plot the dots uppon the bar plot
dataRatio3 <- melt(dataRatio2)
head(dataRatio3)
dataRatio3$condition <- gsub("[1-3]","",dataRatio3$variable)
dataRatio3$condition <- gsub("c","Cecum",dataRatio3$condition)
dataRatio3$condition <- gsub("r","Recto",dataRatio3$condition)
dataRatio3 <- dataRatio3[which(dataRatio3$taxon_name != "unclassified" & dataRatio3$taxon_name != "Others"),]


head(dataRatio3)
library(ggsignif)
p.B <- ggplot(dataRatio3, aes(x=taxon_name, y= value, fill=condition))+
  geom_bar(stat = "summary", fun="mean", position=position_dodge()) +
  geom_point(size=1,alpha=0.7, position=position_dodge(width=0.9)) + 
  geom_errorbar(stat = 'summary', position=position_dodge(width=0.9), width=.2) +
  coord_flip() +
  ylab("Ratio of average relative abundance (MT/MG)") + xlab("") +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous(trans='log2', limits = c(-1,4) ) + 
  theme_classic() +
  ggtitle("") + background_grid(major = "xy", minor = "xy") +
  theme(legend.title=element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x= 6.5 , y = 0.18, label = "*", size=10) +
  annotate("text", x= 2.5 , y = 0.35, label = "*", size=10) +
  annotate("text", x= 7.5 , y = 0.27, label = "*", size=10)
p.B

f1 <- plot_grid(p.A, p.B, labels = c("a", "b"), ncol=1)

save_plot("Figure1.pdf", f1,
          base_width = 6, 
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 0.5
)

