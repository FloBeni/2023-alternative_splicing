# Generate Data 7

options(stringsAsFactors = F, scipen = 999)
library(ape)

pathData = "/home/fbenitiere/data/Projet-SplicedVariants/"

species = "Homo_sapiens"
df = read.table(file = paste(pathData,"Analyses/",species,"/sequencing_depth.tab",sep="") ,header = T, sep="\t")
df$species = species
data_7 = df[df$replicate %in% c(1:30) ,]

species = "Drosophila_melanogaster"
df=read.table(file = paste(pathData,"Analyses/",species,"/sequencing_depth.tab",sep="") ,header = T, sep="\t")
df$species = species
df = df[df$replicate %in% c(1:30) ,]
data_7 = rbind(data_7,df)

write.table(data_7 , paste("data/Data7_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
