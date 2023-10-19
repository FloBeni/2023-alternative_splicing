# Generate Data 4
library(stringr)
# pathData = "/home/XXXXX/data/Projet-SplicedVariants/"
options(stringsAsFactors = F, scipen = 999)


data_10 = read.delim(paste("data/Data10_supp.tab",sep=""),comment.char = "#")
data_10$tissue = sapply(data_10$LibraryName,function(x) str_split(x,"\\.")[[1]][3])
data_10$BioProjet = paste(data_10$tissue,data_10$BioProject)

data_4 = data.frame()
for(species in c("Gallus_gallus","Homo_sapiens","Macaca_mulatta","Monodelphis_domestica","Mus_musculus","Oryctolagus_cuniculus","Rattus_norvegicus")){
  print(species)
  intron_all = read.delim(paste("data/per_species/",species,"_by_intron_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  intron_all$id = paste(intron_all$gene_id,intron_all$splice3,intron_all$splice5,intron_all$seqname,intron_all$strand,sep=";")
  rownames(intron_all) = intron_all$id
  
  fpkm_cov_all = read.delim(paste("data/per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  
  sra_table = data_10[data_10$species == species,]
  print(table(sra_table$BioProjet))
  
  for (Bioproject in unique(sra_table$BioProjet)){
    list_rnaseq = paste(sra_table[sra_table$BioProjet == Bioproject,]$Run ,collapse = ";")
    print(Bioproject)
    fpkm_cov = read.delim(paste(pathData,"Analyses/",species,"/bioproject_analysis/by_gene_",Bioproject,".tab",sep=""))
    fpkm_cov = fpkm_cov[fpkm_cov_all$busco_metazoa ,]
    CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) 
    
    
    by_intron = read.delim(paste(pathData,"Analyses/",species,"/bioproject_analysis/by_intron_",Bioproject,".tab",sep=""))
    
    by_intron = by_intron[by_intron$intron_class =="major" & by_intron$gene_id %in% fpkm_cov$gene_id ,]
    by_intron$id = paste(by_intron$gene_id,by_intron$splice3,by_intron$splice5,by_intron$seqname,by_intron$strand,sep=";")
    by_intron$into_cds = intron_all[by_intron$id,]$into_cds
    
    by_intron = by_intron[by_intron$into_cds == "True",]
    
    # Collecte les données
    data_4 = rbind(data_4,data.frame(
      species,
      SVR = mean(by_intron$splice_variant_rate)*100,
      Bioproject,
      couverture=CoverageBuscoExon,
      nbIntron=nrow(by_intron),
      list_rnaseq
    ))
  }
}

data_4$organs = sapply(data_4$Bioproject,function(x) str_split(x," ")[[1]][1])



tissue_list = c( "_ovar","_head","_test")
for(species in c( "Dendroctonus_ponderosae" )){print(species)
  sra_table = data_10[data_10$species == species ,]
  
  intron_all = read.delim(paste("data/per_species/",species,"_by_intron_analysis.tab.gz",sep=""),  sep="\t")
  intron_all$id = paste(intron_all$gene_id,intron_all$splice3,intron_all$splice5,intron_all$seqname,intron_all$strand,sep=";")
  rownames(intron_all) = intron_all$id
  
  fpkm_cov_all = read.delim(paste("data/per_species/",species,"_by_gene_analysis.tab.gz",sep=""),  sep="\t",comment.char = "#")
  
  for (tissue in tissue_list){print(tissue)
    samples = sra_table[grepl(tissue,sra_table$LibraryName),]$Run
    list_rnaseq = paste(samples ,collapse = ";")
    print(samples)
    by_intron = read.delim(paste(pathData,"Analyses/" , species , "/by_intron_db.tab.gz" , sep=""), header=T , sep="\t",comment.char = "#")
    by_intron = by_intron[colnames(by_intron)[grepl(paste(c("gene_id","seqname","strand","splice5","splice3",samples),collapse = "|"),colnames(by_intron))]]
    by_intron$sum_n1 = rowSums(by_intron[grepl("n1",colnames(by_intron))])
    by_intron$sum_n2 = rowSums(by_intron[grepl("n2",colnames(by_intron))])
    by_intron$sum_n3 = rowSums(by_intron[grepl("n3",colnames(by_intron))])
    
    fpkm_cov = read.delim(paste(pathData,"Analyses/",species,"/by_gene_db.tab.gz",sep=""), header=T , sep="\t",comment.char = "#")
    print(table(fpkm_cov$gene_id == fpkm_cov_all$gene_id))
    fpkm_cov = fpkm_cov[fpkm_cov_all$busco_metazoa ,]
    fpkm_cov = fpkm_cov[colnames(fpkm_cov)[grepl(paste(c("gene_id",samples),collapse = "|"),colnames(fpkm_cov))]]
    fpkm_cov$sum_exon_coverage = rowSums(fpkm_cov[grepl("exon_coverage_",colnames(fpkm_cov))])
    
    CoverageBuscoExon = round(median(tapply(fpkm_cov$sum_exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) 
    
    by_intron = by_intron[by_intron$gene_id %in% fpkm_cov$gene_id ,]
    
    by_intron$splice_variant_rate = by_intron$sum_n2 / (by_intron$sum_n1 + by_intron$sum_n2)
    by_intron[(by_intron$sum_n1 + by_intron$sum_n2) < 10,]$splice_variant_rate = NA
    by_intron$nonsplice_variant_rate = by_intron$sum_n3 / (by_intron$sum_n3 + 2*by_intron$sum_n1)
    by_intron[(by_intron$sum_n3/2 + by_intron$sum_n1) < 10,]$nonsplice_variant_rate = NA
    
    by_intron$intron_class = "Unclassified"
    by_intron[sapply((by_intron$sum_n1 > 0) & (by_intron$nonsplice_variant_rate < 0.5) &
                       (by_intron$splice_variant_rate < 0.5),isTRUE),]$intron_class = "major"
    
    by_intron[sapply((by_intron$sum_n1 > 0) & (by_intron$nonsplice_variant_rate >= 0.5) |
                       (by_intron$splice_variant_rate >= 0.5),isTRUE), ]$intron_class = "minor"
    table(by_intron$intron_class)
    
    by_intron = by_intron[by_intron$intron_class =="major"  ,]
    by_intron$id = paste(by_intron$gene_id,by_intron$splice3,by_intron$splice5,by_intron$seqname,by_intron$strand,sep=";")
    by_intron$into_cds = intron_all[by_intron$id,]$into_cds
    
    by_intron = by_intron[by_intron$into_cds == "True",]
    
    data_4 = rbind(data_4 , data.frame(
      species ,
      SVR = mean(by_intron$splice_variant_rate)*100 ,
      Bioproject = tissue ,
      couverture = CoverageBuscoExon ,
      nbIntron = nrow(by_intron) ,
      organs = tissue ,
      list_rnaseq
    ))
  }
}
data_4[data_4$organs == "_ovar",]$organs = "Ovary"
data_4[data_4$organs == "_head",]$organs = "Head"
data_4[data_4$organs == "_test",]$organs = "Testis"


write.table(data_4,paste("data/Data4_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
