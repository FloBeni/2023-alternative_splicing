options(stringsAsFactors = F, scipen = 999)
library(cowplot)
library(ggplot2)
library(magick)
library(caper)
library(imager)
library(ape)
library(stringr)
library(png)
library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))

resolution = 3
path_require = "images_library/"
path_figure = "pannels/"
path_pannel = "figure/"

pathData = "~/data/Projet-SplicedVariants/"

lm_eqn <- function(m=lm(Y ~ X,data)){
  if (summary(m)$coefficients[2,4] < 0.01){
    paste(" =", format(summary(m)$r.squared, digits = 2) , "   p =",formatC(summary(m)$coefficients[2,4], format = "e", digits = 1))}
  else {paste(" =", format(summary(m)$r.squared, digits = 2) , "   p =",format(summary(m)$coefficients[2,4], digits = 2))}
}

std <- function(x) sd(x)/sqrt(length(x))

vectorColor=c(Diptera="red",Coleoptera="#0a4413",
              Hymenoptera="#ba8e18",Blattodea="darkblue",
              Hemiptera="#00BFFF",Lepidoptera="#419f51",other="grey",
              Testudines="blue",Mammalia="#66281A",Crocodylia="#5c5c5c"
              ,Aves="#000000",Teleostei="#579392",Anura="#39ABAD",Chondrichthyes="#8b27bc",
              Testudines="blue","Monodelphis domestica"="#66281A",Mammalia="#66281A","Rattus norvegicus"="#66281A",
              "Mus musculus"="#66281A","Macaca mulatta"="#66281A","Oryctolagus cuniculus"="#66281A","Homo sapiens"="#66281A"
              ,Crocodylia="#5c5c5c","Gallus gallus"="#000000",Teleostei="#579392",Anura="#39ABAD",Chondrichthyes="#1A6566")

vectorColor = vectorColor[c("Hymenoptera","Diptera","Lepidoptera","Coleoptera","Blattodea","Hemiptera",
                            "Mammalia","Crocodylia","Aves","Chondrichthyes")]

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))

data_1 = read.delim("data/Data1_supp.tab",comment.char = "#")
data_1 = data_1[data_1$species %in% arbrePhylo$tip.label,]
rownames(data_1) = data_1$species
data_1$clade = factor(data_1$clade,levels=names(vectorColor))

