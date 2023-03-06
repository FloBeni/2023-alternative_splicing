# Generate Data 3

options(stringsAsFactors = F, scipen = 999)
library(ape)

std <- function(x) sd(x)/sqrt(length(x))

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))


data_3 = data.frame()
for (species in arbrePhylo$tip.label){print(species)
  fpkm_cov = read.delim(paste("data/per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  
  SPintron_original = read.delim(paste("data/per_species/",species,"_by_intron_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  
  SPintron_original = SPintron_original[SPintron_original$gene_id  %in% fpkm_cov$gene_id ,]
  SPintron_original = SPintron_original[SPintron_original$into_cds == "True",]
  SPintron_original = SPintron_original[SPintron_original$n1 != 0,]
  
  Xaxis = 1 - SPintron_original$splice_variant_rate
  quantile = seq(0, 1,5/100)
  intervalle = cut(Xaxis, quantile,include.lowest = T,include.higher=T)
  ras = table(intervalle)/sum(table(intervalle)) * 100
  
  Xaxis=1-SPintron_original$nonsplice_variant_rate
  quantile = seq(0, 1,5/100)
  intervalle = cut(Xaxis , quantile , include.lowest = T,include.higher=T)
  rans = table(intervalle)/sum(table(intervalle)) * 100
  
  
  data_3 = rbind(data_3,data.frame(
    species,
    rate = seq(5/100/2, 1,(5/100)),
    ras,
    rans
  ))
}
write.table(data_3 , paste("data/Data3_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)




