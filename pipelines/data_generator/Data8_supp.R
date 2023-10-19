# Generate Data 8
options(stringsAsFactors = F, scipen = 999)
library(ape)

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))
sp_studied = arbrePhylo$tip.label

std <- function(x) sd(x)/sqrt(length(x))
data_8 = data.frame()
for (species in sp_studied){
  print(species)
  fpkm_cov = read.delim(paste("data/per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  rownames(fpkm_cov) = fpkm_cov$gene_id
  
  minor_introns = read.delim(paste("data/per_species/",species,"_by_intron_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  minor_introns = minor_introns[minor_introns$intron_class == "minor" & minor_introns$into_cds == "True" & minor_introns$gene_id %in% fpkm_cov$gene_id,]
  minor_introns = minor_introns[minor_introns$criptic_intron == "False" & !is.na(minor_introns$criptic_intron),]
  minor_introns = minor_introns[( minor_introns$distance_from_major < 30 ) ,]
  print(table(minor_introns$which_shared_site))
  
  fpkm_cov = fpkm_cov[fpkm_cov$busco_metazoa ,]
  intron_busco = minor_introns[minor_introns$gene_id %in% fpkm_cov$gene_id,]
  
  prop_fp_sv_abundant = sum(minor_introns$mira > 0.05 & minor_introns$frame_shift == 0) / sum(minor_introns$mira > 0.05)
  prop_fp_sv_rare = sum(minor_introns$mira <= 0.05 & minor_introns$frame_shift == 0) / sum(minor_introns$mira <= 0.05)
  prop_fp_sv_all = sum( minor_introns$frame_shift == 0) / nrow(minor_introns)
  
  prop_fp_sv_abundant_busco = sum(intron_busco$mira > 0.05 & intron_busco$frame_shift == 0) / sum(intron_busco$mira > 0.05)
  prop_fp_sv_rare_busco = sum(intron_busco$mira <= 0.05 & intron_busco$frame_shift == 0) / sum(intron_busco$mira <= 0.05)
  prop_fp_sv_all_busco = sum( intron_busco$frame_shift == 0) / nrow(intron_busco)
  
  
  minor_introns = minor_introns[minor_introns$mira > 0.05,]
  
  xaxis = minor_introns[,"mira"]
  proportion = 20/100
  quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
  intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  X = tapply(xaxis, intervalle, mean)
  XerrorBar = tapply(xaxis, intervalle, std)
  
  
  Y = tapply(minor_introns$frame_shift, intervalle, function(x) {
    return(sum(x == "0")/sum(x != "05") * 100)
  })
  
  value = lm(Y ~ X)$coefficients[1] + lm(Y ~ X)$coefficients[2] * 0.1
  
  
  data_8 = rbind(data_8,data.frame(
    species,
    reg_linear = value,
    prop_fp_sv_abundant,
    prop_fp_sv_abundant_busco,
    prop_fp_sv_rare,
    prop_fp_sv_rare_busco,
    prop_fp_sv_all,
    prop_fp_sv_all_busco
  ))
} 

write.table(data_8 , paste("data/Data8_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)




