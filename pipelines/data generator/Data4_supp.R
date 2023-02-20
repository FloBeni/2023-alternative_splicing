# Generate Data 4

library(stringr)
pathData="~/data/Projet-SplicedVariants/"
options(stringsAsFactors = F, scipen = 999)


df = data.frame()
for(species in c("Gallus_gallus","Homo_sapiens","Macaca_mulatta","Monodelphis_domestica","Mus_musculus","Oryctolagus_cuniculus","Rattus_norvegicus")){print(species)
  # species="Gallus_gallus"
  intron_all = read.delim(paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  intron_all$id = paste(intron_all$gene_id,intron_all$splice3,intron_all$splice5,intron_all$seqname,intron_all$strand,sep=";")
  rownames(intron_all) = intron_all$id
  
  sra_table = read.delim(paste(pathData,"Annotations/",species,"/SRAruninfo.tab",sep=""))
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
    df = rbind(df,data.frame(
      species,
      SVR=mean(by_intron$splice_variant_rate)*100,
      Bioproject,
      couverture=CoverageBuscoExon,
      nbIntron=nrow(by_intron)
    ))
  }
}

df$organs = sapply(df$Bioproject,function(x) str_split(x," ")[[1]][1])

#selection des bioprojets interessants, garde pour les mammifères que le forebrain


df = read.delim("data/Data4_supp.tab",comment.char = "#")
df = df[df$species != "Callorhinchus_milii",]




tissue_list = c( "_ovar","_head","_test")
for(species in c("Dendroctonus_ponderosae")){print(species)
  sra_table = read.delim(paste(pathData,"Annotations/",species,"/SRAruninfo.tab",sep=""))
  
  intron_all = read.delim(paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  intron_all$id = paste(intron_all$gene_id,intron_all$splice3,intron_all$splice5,intron_all$seqname,intron_all$strand,sep=";")
  rownames(intron_all) = intron_all$id
  
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  
  for (tissue in tissue_list){print(tissue)
    samples = sra_table[grepl(tissue,sra_table$LibraryName),]$Run
    print(samples)
    by_intron = read.delim(paste(pathData,"Analyses/",species,"/by_intron_db.tab.gz",sep=""), header=T , sep="\t",comment.char = "#")
    by_intron = by_intron[colnames(by_intron)[grepl(paste(c("gene_id","seqname","strand","splice5","splice3",samples),collapse = "|"),colnames(by_intron))]]
    by_intron$sum_n1 = rowSums(by_intron[grepl("n1",colnames(by_intron))])
    by_intron$sum_n2 = rowSums(by_intron[grepl("n2",colnames(by_intron))])
    by_intron$sum_n3 = rowSums(by_intron[grepl("n3",colnames(by_intron))])
    
    fpkm_cov = read.delim(paste(pathData,"Analyses/",species,"/by_gene_db.tab.gz",sep=""), header=T , sep="\t",comment.char = "#")
    fpkm_cov = fpkm_cov[fpkm_cov$gene_id %in% busco_tab$gene_id ,]
    fpkm_cov = fpkm_cov[colnames(fpkm_cov)[grepl(paste(c("gene_id",samples),collapse = "|"),colnames(fpkm_cov))]]
    fpkm_cov$sum_exon_coverage = rowSums(fpkm_cov[grepl("exon_coverage_",colnames(fpkm_cov))])
    
    CoverageBuscoExon = round(median(tapply(fpkm_cov$sum_exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco
    
    
    busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                            !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
    
    by_intron = by_intron[by_intron$gene_id %in% busco_tab$gene_id ,]
    
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
    
    by_intron = by_intron[by_intron$intron_class =="major" & by_intron$gene_id %in% busco_tab$gene_id ,]
    by_intron$id = paste(by_intron$gene_id,by_intron$splice3,by_intron$splice5,by_intron$seqname,by_intron$strand,sep=";")
    by_intron$into_cds = intron_all[by_intron$id,]$into_cds
    
    by_intron = by_intron[by_intron$into_cds =="True",]
    
    # Collecte les données
    df = rbind(df , data.frame(
      species ,
      SVR = mean(by_intron$splice_variant_rate)*100 ,
      Bioproject = tissue ,
      couverture = CoverageBuscoExon ,
      nbIntron = nrow(by_intron) ,
      organs = tissue
    ))
  }
}
df[df$organs=="_ovar",]$organs = "Ovary"
df[df$organs=="_head",]$organs = "Head"
df[df$organs=="_test",]$organs = "Testis"




write.table(df,paste("data/Data4_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
