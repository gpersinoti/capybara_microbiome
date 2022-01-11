#' ---
#' title: "Create the Metabolic reconstruction of Capybara Microbiome"
#' author: "Gabriela F. Persinoti"
#' date: "March 3rd, 2020"
#' output: html_document
#' ---

library(RColorBrewer)
library(phylogram)
library(ape)
library(readr)
library(dplyr)
library(stringr)
library(pheatmap)

setwd("/data/Projects/Capibara/github/Figure4/")
read.table("MAG_to_ID.txt", header = F, stringsAsFactors = F) -> samples

read.table("bin_number.txt", sep = "\t", header = F, stringsAsFactors = F) -> bins
colnames(bins) <- c("V3","V1")

inner_join(samples,bins) -> samples_mag
row.names(samples_mag) <- samples_mag$V1

file.exists(paste("top5/",samples_mag$V2,
                  "/kegg_mapper.tsv", sep = "")) -> exists
samples_mag[exists,] -> samples_mag


#Loading bins phylogeny
read.tree("bins_capybara_new.tree.nwk") -> tree
#using Capybara_CONCOCT.68 as root
read.tree(file="bins_capybara_new.tree_rooted.nwk") -> tree
#root(tree,outgroup = "Recto_4.CONCOCT.0",resolve.root = T,) -> rooted
phylogram::as.dendrogram(tree) -> dend
as.hclust(dend) -> hclust
hclust$labels[hclust$order] -> order_samples
order_samples

#Loading compounds list
all_KEGG <- read.table("compounds4heatmap_artigo.csv", sep ="\t",
                       header = F, stringsAsFactors = F)
all_KEGG$V3 <- NULL
rownames(all_KEGG) <- all_KEGG$V2
all_KEGG$V1 -> full_names

list_boolean <- c()
#Identifying the compounds in the MAGs
for( i in 1:nrow(samples_mag)){
  amon_res <- read_delim(paste("top5/",samples_mag$V2[i],
                               "/kegg_mapper.tsv", sep = ""),
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  for (l in 1:nrow(all_KEGG)){
    ifelse(rownames(all_KEGG)[l] %in% amon_res$...1,1,0) -> list_boolean[[l]]
  }
  
  all_KEGG[,i] <- list_boolean 
  colnames(all_KEGG)[i] <- samples_mag$V1[i]
}
rownames(all_KEGG) <- full_names

#Generating the colours

ann_colors = list(
  Phylum= c(Actinobacteriota = "#E41A1C", Bacteroidetes= "#566B9B", Euryarchaeota ="#6F8272", Fibrobacterota= "#FFF32E",
            Firmicutes= "#B45B76", Fusobacteriota= "#FFA000", Planctomycetota= "#DAC2AE",Proteobacteria= "#B6E3BF", 
            Spirochaetota= "#F3F2BA"))


read.table("mags_phylum.txt", row.names = 1,skip = 1, stringsAsFactors = F) -> phylum

colnames(phylum) <- "Phylum"

#Color according to the phylum
all_KEGG -> KEGG_colors

for(i in 1:ncol(KEGG_colors)){
  if (phylum[colnames(KEGG_colors)[i],] == "Actinobacteriota"){
    ifelse(KEGG_colors[,i] == 1,1,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Firmicutes"){
    ifelse(KEGG_colors[,i] == 1,2,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Synergistota"){
    ifelse(KEGG_colors[,i] == 1,3,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Fusobacteriota"){
    ifelse(KEGG_colors[,i] == 1,4,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Bacteroidota"){
    ifelse(KEGG_colors[,i] == 1,5,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Fibrobacterota"){
    ifelse(KEGG_colors[,i] == 1,6,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Planctomycetota"){
    ifelse(KEGG_colors[,i] == 1,7,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Spirochaetota"){
    ifelse(KEGG_colors[,i] == 1,8,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Proteobacteria"){
    ifelse(KEGG_colors[,i] == 1,9,0) -> KEGG_colors[,i]
  }
  if (phylum[colnames(KEGG_colors)[i],] == "Euryarchaeota"){
    ifelse(KEGG_colors[,i] == 1,10,0) -> KEGG_colors[,i]
  }
}

row_ann <- read.table("compounds4heatmap_artigo.csv", sep ="\t",
                      header = F, stringsAsFactors = F)

rownames(row_ann) <- row_ann$V1
row_ann$V2 <- NULL
row_ann$V1 <- NULL
colnames(row_ann) <- "Metabolites"

breaksList = seq(0, 10, by = 1)
ann_colors = list(
  Phylum= c(Actinobacteriota = "#E41A1C",Bacteroidota= "#566B9B",Euryarchaeota ="#6F8272", Fibrobacterota= "#FFF32E",Firmicutes= "#B45B76",
            Fusobacteriota= "#FFA000", Planctomycetota= "#EDBB99",Proteobacteria= "#B6E3BF",
            Spirochaetota= "#F3F2BA",Synergistota= "#AEC5D9"))

pheatmap(KEGG_colors[,order_samples]  ,cluster_rows = F,
         cluster_cols = hclust, legend = F, border_color = "NA",
         color = c("#000000","#E41A1C", "#B45B76", "#AEC5D9", "#FFA000","#566B9B",
                   "#FFF32E","#EDBB99","#F3F2BA","#B6E3BF", "#6F8272"),
         fontsize_row =11,fontsize_col = 5 ,cellwidth = 5,
         gaps_row = 17, annotation_colors =  ann_colors,
         show_colnames =  T, labels_col =  samples_mag[order_samples,3],
         width = 7.5, height = 8, 
         filename = "Figure4.pdf"
         )


