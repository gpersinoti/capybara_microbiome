library(tidyr)
library(readr)

setwd("/data/Projects/Capibara/github/Figure3/")
##################### Heatmap #################

#Set cazymes of interest
heat_list <- c("GH1","GH2","GH3","GH4","GH5","GH8","GH9",
               "GH10","GH13","GH16","GH26","GH28","GH38",
               "GH43","GH45","GH51","GH78","GH92","GH97","CE1","CE2","CE4","CE6",
               "CE7", "CE8","CE9","CE11","CE12","CE15")

as.data.frame(heat_list) -> heat_list


for(i in 1:79){
  #Load MAG annotation
  MAG <- read_delim(paste0("MAGs/MAG", i,".txt"),
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
  
  #Separate modules in different rows
  separate_rows(MAG, 2, sep = "-", convert = FALSE) -> MAG
  #Remove subfamilies information
  stringr::str_replace_all(MAG$X2, pattern = "\\_[:digit:]*$","") -> MAG$X2
  table(MAG$X2)
  boolean <- c()
  for(l in 1:nrow(heat_list)){
    #Getting the number for each family if it's present in the annotation
    ifelse(heat_list$heat_list[l] %in% names(table(MAG$X2)), table(MAG$X2)[heat_list$heat_list[l]],0) -> boolean[l]
  }
  #Integrate annotation and give the proper name
  heat_list[,i+1] <- boolean
  colnames(heat_list)[i+1] <- paste0("MAG", i)
  rm(boolean)
}

heat_list$heat_list -> row.names(heat_list)
heat_list$heat_list <- NULL

breaks=seq(from=-1,to=27,by=1)
heat_list[is.na(heat_list)] <- 0
#Remove zero labels from figure
t(heat_list) -> non_zero
non_zero[non_zero == 0] <- ""

#Load phylum annotation
read.table("phylum.tsv") -> phylum
#Add colors for phylum and enzyme class
ann_colors = list(
  Phylum= c(Actinobacteriota = "#E41A1C",
            Bacteroidota= "#566B9B",
            Euryarchaeota ="#6F8272",
            Fibrobacterota= "#FFF32E",
            Firmicutes= "#B45B76",
            Fusobacteriota= "#FFA000",
            Planctomycetota= "#DAC2AE",
            Proteobacteria= "#B6E3BF",
            Spirochaetota= "#F3F2BA",
            Synergistota= "#AEC5D9"
  ),
  Class=c(CE= "#CC4D66",GH= "#89C6AF"))

#Add colors for the class of enzyme
cazy_colors <-c("GH1","GH2","GH3","GH4","GH5","GH8","GH9",
                "GH10","GH13","GH16","GH18","GH26","GH28","GH36","GH38",
                "GH43","GH45","GH51","GH78","GH92","GH97","CE1","CE2","CE4","CE6",
                "CE7", "CE8","CE9","CE11","CE12","CE15")

as.data.frame(cazy_colors) -> cazy_colors
row.names(cazy_colors) <- cazy_colors$cazy_colors
stringr::str_replace(cazy_colors$cazy_colors, pattern = "[:digit:]*$", replacement = "") -> cazy_colors$cazy_colors
colnames(cazy_colors) <- 'Class'

library(RColorBrewer)
my_palette <- colorRampPalette(c("white", "blue", "darkblue"))(n = 27)


head(heat_list)
head(t(heat_list))

library(pheatmap)
figA <- pheatmap(t(heat_list), color= my_palette,
                 cellwidth = 9, cellheight = 8,
                 fontsize_row = 8, breaks=seq(from=-1,to=27,by=1),
                 annotation_colors =  ann_colors, annotation_row = phylum, annotation_col = cazy_colors,
                 cluster_rows = F,cluster_cols = F, display_numbers = non_zero,
                 border_color = F, number_color = ifelse(t(heat_list)> 8,"white","grey30"))
pdf(file="Figure3A.pdf", height = 12)
figA
dev.off()



#Load Supplementary table 6 with PULs and CCs
read.table("new_s6.tsv",header = T,sep = "\t") -> s6_table
s6_table[s6_table$category =="PUL",1] -> Pul_count
s6_table[s6_table$category =="CC",1] -> CC_count

df <- data.frame(PUL=numeric(),
                 CC=numeric())

PUL_boolean <- c()
CC_boolean <- c()
#Count for each MAG the number of times they appear
for(i in 1:79){
  ifelse(paste0("MAG", i) %in% names(table(Pul_count)), table(Pul_count)[paste0("MAG", i)],0) -> PUL_boolean[i]
  ifelse(paste0("MAG", i) %in% names(table(CC_count)), table(CC_count)[paste0("MAG", i)],0) -> CC_boolean[i]
}


cbind(as.data.frame(PUL_boolean),as.data.frame(CC_boolean)) -> pul_heat
colnames(pul_heat) <- c("PULs","CCs")

pul_heat -> matrix_values
matrix_values[matrix_values == 0] <- ""

stringr::str_replace_all(row.names(pul_heat),"^","MAG") -> row.names(pul_heat)

my_palette <- colorRampPalette(c("white", "red", "darkred"))(n = 30)
breaks=seq(from=-1,to=30,by=1)

figB <- pheatmap(pul_heat, color= my_palette,
                 cellwidth = 9, cellheight = 8,
                 fontsize_row = 8, #breaks=seq(from=-1,to=27,by=1),
                 cluster_rows = F,cluster_cols = F, display_numbers = matrix_values,
                 border_color = F, number_color = ifelse(pul_heat > 9,"white","grey30") )
#Figures A and B were manually joined  latter


pdf(file="Figure3.B.pdf", height = 12)
figB
dev.off()
