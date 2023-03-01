options(stringsAsFactors = F, scipen = 999)
library(readxl)
library(ggplot2)
library(caper)
library(imager)
library(stringr)
library(cowplot)
library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))


std <- function(x) sd(x)/sqrt(length(x))

pathData = "~/data/Projet-SplicedVariants/"

resolution = 3
path_require = "images_library/"
path_figure = "pannels/"
path_pannel = "figure_refPCI/"


vectorColor=c(Diptera="red",Coleoptera="#0a4413",
              Hymenoptera="#ba8e18",Blattodea="darkblue",
              Hemiptera="#00BFFF",Lepidoptera="#419f51",other="grey",
              Testudines="blue",Mammalia="#66281A",Crocodylia="#5c5c5c"
              ,Aves="#000000",Teleostei="#579392",Anura="#39ABAD",Chondrichthyes="#8b27bc",
              Testudines="blue","Monodelphis domestica"="#66281A",Mammalia="#66281A","Rattus norvegicus"="#66281A",
              "Mus musculus"="#66281A","Macaca mulatta"="#66281A","Oryctolagus cuniculus"="#66281A","Homo sapiens"="#66281A"
              ,Crocodylia="#5c5c5c","Gallus gallus"="#000000",Teleostei="#579392",Anura="#39ABAD",Chondrichthyes="#8b27bc")



vectorColor = vectorColor[c("Hymenoptera","Diptera","Lepidoptera","Coleoptera","Blattodea","Hemiptera",
                            "Mammalia","Crocodylia","Aves","Chondrichthyes")]


lm_eqn <- function(m=lm(Y ~ X,data)){
  if (summary(m)$coefficients[2,4] < 0.01){
    paste(" =", format(summary(m)$r.squared, digits = 2) , "   p =",formatC(summary(m)$coefficients[2,4], format = "e", digits = 1))}
  else {paste(" =", format(summary(m)$r.squared, digits = 2) , "   p =",format(summary(m)$coefficients[2,4], digits = 2))}
}

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))

data_1 = read.delim("data/Data1_supp.tab",comment.char = "#")
data_1 = data_1[data_1$species %in% arbrePhylo$tip.label,]
rownames(data_1) = data_1$species


data_8 = read.delim(paste("data/Data8_supp.tab",sep=""),comment.char = "#")
data_8$clade = data_1[data_8$species,]$clade
data_8$longevity = data_1[data_8$species,]$longevity
data_8$body_size = data_1[data_8$species,]$body_size
data_8$dNdS = data_1[data_8$species,]$dNdS_200k
data_8$ratio =  data_8$prop_fp_sv_abundant * 100
