# Generate Data 4

library(stringr)
pathData="~/data/Projet-SplicedVariants/"
options(stringsAsFactors = F, scipen = 999)


df = data.frame()
for(species in c("Gallus_gallus","Homo_sapiens","Macaca_mulatta","Monodelphis_domestica","Mus_musculus","Oryctolagus_cuniculus","Rattus_norvegicus")){print(species)
  # species="Homo_sapiens"
  intron_all = read.delim(paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  intron_all$id = paste(intron_all$gene_id,intron_all$splice3,intron_all$splice5,intron_all$seqname,intron_all$strand,sep=";")
  rownames(intron_all) = intron_all$id
  
  listBioproject=read.delim(paste(pathData,"RNAseq_table/",species,"/list_Acc.tab",sep=""))
  for (Bioproject in unique(listBioproject$BioProjet)){
    print(Bioproject)
    
    busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
    
    fpkm_cov = read.delim(paste(pathData,"Analyses/",species,"/bioproject_analysis/by_gene_",Bioproject,".tab",sep=""))
    fpkm_cov = fpkm_cov[fpkm_cov$gene_id %in% busco_tab$gene_id ,]
    CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco
    
    
    by_intron = read.delim(paste(pathData,"Analyses/",species,"/bioproject_analysis/by_intron_",Bioproject,".tab",sep=""))
    
    busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                            !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
    
    by_intron = by_intron[by_intron$intron_class =="major" & by_intron$gene_id %in% busco_tab$gene_id ,]
    by_intron$id = paste(by_intron$gene_id,by_intron$splice3,by_intron$splice5,by_intron$seqname,by_intron$strand,sep=";")
    by_intron$into_cds = intron_all[by_intron$id,]$into_cds
    
    by_intron = by_intron[by_intron$into_cds =="True",]
    
    # Collecte les données
    df=rbind(df,data.frame(
      species,
      SVR=mean(by_intron$splice_variant_rate)*100,
      Bioproject,
      couverture=CoverageBuscoExon,
      nbIntron=nrow(by_intron)
    ))
  }
}

df$organs = sapply(df$Bioproject,function(x) str_split(x," ")[[1]][1])
df$Bioproject = sapply(df$Bioproject,function(x) str_split(x," ")[[1]][2])

#selection des bioprojets interessants, garde pour les mammifères que le forebrain
write.table(df,paste("data/Data4_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
