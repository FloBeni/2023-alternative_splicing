# Generate Data 2

options(stringsAsFactors = F, scipen = 999)
library(ape)


std <- function(x) sd(x)/sqrt(length(x))
pathData="/home/fbenitiere/data/Projet-SplicedVariants/"
# pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))
sp_studied = arbrePhylo$tip.label

std <- function(x) sd(x)/sqrt(length(x))
data_2 = data.frame()
for (species in sp_studied) {# recupere les especes a analyser
  print(species)
  
  minor_introns = read.delim(paste(pathData,"per_species/",species,"_by_intron_analysis.tab.gz",sep=""),  sep="\t")
  fpkm_cov = read.delim(paste(pathData,"per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t")
  
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  rownames(fpkm_cov) = fpkm_cov$gene_id
  
  # minor_introns = read.delim(file=paste(pathData,"Analyses/",species,"/by_minor_intron.tab",sep=""))
  minor_introns = minor_introns[minor_introns$intron_class == "minor" & minor_introns$into_cds == "True" & minor_introns$gene_id %in% fpkm_cov$gene_id,]
  
  # minor_introns = minor_introns[minor_introns$which_shared_site != "both",]
  minor_introns = minor_introns[minor_introns$criptic_intron == "False",]
  minor_introns = minor_introns[minor_introns$distance_from_major < 30 ,]
  
  xaxis = minor_introns[,"mira"]
  proportion = 5/100
  quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
  intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  X=tapply(xaxis, intervalle, mean)
  XerrorBar=tapply(xaxis, intervalle, std)
 
  
  Y = tapply(minor_introns$frame_shift, intervalle, function(x) {
    return(sum(x=="0")/length(x) )
  })
  
  data_2 = rbind(data_2,data.frame(species,
                                   average_mira=X,
                                   framepreserving_proportion=Y,
                                   XerrorBar=XerrorBar,
                                   table(intervalle)
                                   
  ))
}


write.table(data_2,paste("data/Data2_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

