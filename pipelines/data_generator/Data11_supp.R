# Generate Data 11
options(stringsAsFactors = F, scipen = 999)
library(ape)

arbrePhylo = read.tree(paste("data/phylogenetic_tree.nwk",sep=""))
sp_studied = arbrePhylo$tip.label

data_11 = data.frame()
for (species in sp_studied){print(species)

  by_intron = read.delim(file=paste("data/per_species/",species,"_by_intron_analysis.tab.gz",sep=""), header=T , sep="\t",comment.char = "#")
  by_intron = by_intron[!is.na(by_intron$phase),]

  fpkm_cov = read.delim(paste("data/per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]

  by_intron = by_intron[by_intron$gene_id %in% fpkm_cov$gene_id,]

  ratio = sum(!grepl(",",by_intron$phase)) / nrow(by_intron)
  by_intron = by_intron[ !grepl(",",by_intron$phase) ,]
  print(table( by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id))
  major_introns = by_intron[by_intron$intron_class == "major" & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]

  fpkm_cov = fpkm_cov[fpkm_cov$busco_metazoa ,]

  major_introns_busco = major_introns[major_introns$gene_id %in% fpkm_cov$gene_id,]

  data_11 = rbind(data_11,data.frame(species,
                                 ratio,
                                 major_introns = table(major_introns$phase) / nrow(major_introns),
                                 major_introns_svr = tapply(major_introns$splice_variant_rate,major_introns$phase,mean) * 100,
                                 major_introns_busco = table(major_introns_busco$phase) / nrow(major_introns_busco),
                                 major_introns_busco_svr = tapply(major_introns_busco$splice_variant_rate,major_introns_busco$phase,mean) * 100
  ))
}

write.table(data_11,paste("data/Data11_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)