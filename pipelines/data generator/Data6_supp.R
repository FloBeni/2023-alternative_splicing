# Generate Data 6

options(stringsAsFactors = F, scipen = 999)
std <- function(x) sd(x)/sqrt(length(x))

pathData = "/home/XXXXX/data//Projet-SplicedVariants/"

n1_n2_treshold = 100

data_6 = data.frame()
for (species in c("Homo_sapiens","Drosophila_melanogaster")){
  fpkm_cov = read.delim(paste(pathData,"/Analyses/",species,"/by_gene_analysis.tab",sep="") , header=T , sep="\t",comment.char = "#")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  rownames(fpkm_cov) = fpkm_cov$gene_id
  all_major = read.delim(paste(pathData,"/Analyses/",species,"/by_intron_major_overlap.tab",sep=""))
  all_major = all_major[all_major$gene_id %in% fpkm_cov$gene_id,]
  all_major = all_major[all_major$intron_class == "major" & all_major$into_cds == "True",]
  all_major$fpkm = fpkm_cov[all_major$gene_id,]$weighted_fpkm
  all_major$n1_n2_treshold =  (all_major$n1 +all_major$n2_spl3+all_major$n2_spl5) >= n1_n2_treshold
  all_major = all_major[all_major$n1_n2_treshold,]
  
  intron = all_major
  xaxis = unlist(intron$fpkm)
  proportion = .05
  
  quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
  intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  X = tapply(xaxis, intervalle, median)
  XerrorBar = tapply(xaxis, intervalle, std)
  
  yaxis = unlist(intron$splice_variant_rate)
  
  Y = tapply(yaxis, intervalle, mean)
  YerrorBar = tapply(yaxis, intervalle, std)
  
  data_6 = rbind(data_6,data.frame(median_gene_expression=X ,average_as=Y*100 ,median_gene_expression_errorBar=XerrorBar ,
                                     average_as_errorBar=YerrorBar*100,effectif=table(intervalle),intron="All",
                                     group=paste("All protein-codingenterN=",nrow(intron))
                                     ,species ))
  
  
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                          !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
  
  intron = all_major[all_major$gene_id %in% busco_tab$gene_id,]
  
  xaxis = unlist(intron$fpkm)
  proportion = .1
  
  quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
  intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  X = tapply(xaxis, intervalle, median)
  XerrorBar = tapply(xaxis, intervalle, std)
  
  yaxis = unlist(intron$splice_variant_rate)
  
  Y = tapply(yaxis, intervalle, mean)
  YerrorBar = tapply(yaxis, intervalle, std)
  
  data_6 = rbind(data_6,data.frame(median_gene_expression=X ,average_as=Y*100 ,median_gene_expression_errorBar=XerrorBar ,
                                     average_as_errorBar=YerrorBar*100,effectif=table(intervalle),intron="All",group=paste("BUSCOenterN=",nrow(intron))
                                     ,species ))
  
  
  intron = all_major[  all_major$have_abundant_sv == "False"   ,]
  
  xaxis = unlist(intron$fpkm)
  proportion = .05
  
  quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
  intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
  
  X = tapply(xaxis, intervalle, median)
  XerrorBar = tapply(xaxis, intervalle, std)
  
  yaxis = unlist(intron$splice_variant_rate)
  
  Y = tapply(yaxis, intervalle, mean)
  YerrorBar = tapply(yaxis, intervalle, std)
  
  data_6 = rbind(data_6,data.frame(median_gene_expression=X ,average_as=Y*100 ,median_gene_expression_errorBar=XerrorBar ,
                                     average_as_errorBar=YerrorBar*100,effectif=table(intervalle),intron="Rare_SV",group=paste("All protein-codingenterN=",nrow(intron))
                                     ,species ))
  
  
}


write.table(data_6 , paste("data/Data6_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
