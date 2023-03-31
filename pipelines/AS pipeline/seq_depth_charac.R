options(stringsAsFactors = F, scipen = 999)

args = (commandArgs(TRUE))
species = args[1]
pathData = args[2]
by_intron_cds_path = args[3]
IntronCoord_path = args[4]
busco_path = args[5]
path_output = args[6]
path_RNAseq_table = args[7]
# species="Homo_sapiens"
# pathData="~/data/Projet-SplicedVariants/"
# # pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"
# by_intron_cds_path = "/home/fbenitiere/data/Projet-SplicedVariants/Analyses/Homo_sapiens/by_intron_cds.tab"
# IntronCoord_path = "/home/fbenitiere/data/Projet-SplicedVariants/Annotations/Homo_sapiens/formatted_data/IntronCoords.tab"
# path_RNAseq_table = paste("/home/fbenitiere/data/Projet-SplicedVariants/","RNAseq_table/Homo_sapiens/list_Acc.tab",sep="")
# busco_path = "/home/fbenitiere/data/Projet-SplicedVariants/Annotations/Homo_sapiens/busco_analysis/busco_to_gene_id_metazoa"


RNAseq_table <- read.delim(path_RNAseq_table) # Lecture du fichier RNA-seq
busco_to_gene = read.delim(busco_path)


# SELECTION DES INTRONS AVEC INTRONLIBRARY QUI CHANGAIENT
IntronCoords = read.delim(IntronCoord_path)

library(dplyr)
library(tidyr)
IntronCoords = IntronCoords %>%
  mutate(Genes = strsplit(as.character(Genes), ",")) %>%
  unnest(Genes) %>%
  filter(Genes != "")

IntronCoords$id = apply(IntronCoords,1,function(x){
  if (as.numeric(x["Strand"]) == 1){
    paste(as.numeric(x["Start"]) ,as.numeric(x["End"]) ,x["Chr"],x["Genes"],as.numeric(x["Strand"]),sep=";")
  } else if (as.numeric(x["Strand"]) == -1){
    paste(as.numeric(x["End"]) ,as.numeric(x["Start"]) ,x["Chr"],x["Genes"],as.numeric(x["Strand"]),sep=";")
  }
})


by_intron_cds = read.table(by_intron_cds_path, header=T , sep="\t")
by_intron_cds$id = paste(by_intron_cds$splice5 , by_intron_cds$splice3 , by_intron_cds$seqname , by_intron_cds$gene_id , by_intron_cds$strand,sep=";")
rownames(by_intron_cds) = by_intron_cds$id


busco_unique = busco_to_gene[!(duplicated(busco_to_gene$busco_id,fromLast = FALSE) | duplicated(busco_to_gene$busco_id,fromLast = TRUE)) &
                               !(duplicated(busco_to_gene$gene_id,fromLast = FALSE) | duplicated(busco_to_gene$gene_id,fromLast = TRUE)) ,]

by_intron_cds$busco_id = by_intron_cds$gene_id %in% busco_unique$gene_id




by_gene = read.delim(paste(pathData , "Analyses-RNAseq/",species,"/",RNAseq_table[1,1],"/coverage_fpkm_by_genes.tab",sep=""))
by_gene$busco_id = by_gene$gene_id %in% busco_to_gene$gene_id


summaryTable = read.delim(paste(pathData , "Analyses-RNAseq/",species,"/",RNAseq_table[1,1],"/n1_n2_n3_by_intron.tab",sep=""))
summaryTable$id = paste(summaryTable$splice5,summaryTable$splice3,summaryTable$seqname,summaryTable$gene_id,summaryTable$strand,sep=";")
rownames(summaryTable) = summaryTable$id

summaryTable$intoIntronLibrary = summaryTable$id %in% by_intron_cds$id
summaryTable$intoIntronCoord = summaryTable$id %in% IntronCoords$id


summaryTable$into_cds = by_intron_cds[summaryTable$id,"into_cds"]
summaryTable$busco_id = by_intron_cds[summaryTable$id,"busco_id"]

library = summaryTable



intron_SELECTION = library$into_cds == "True" & library$intoIntronLibrary & library$busco_id
# intron_SELECTION = library$into_cds == "True" & library$intoIntronLibrary 
intron_SELECTION_annotated = library$intoIntronCoord
intron_SELECTION_annotated = intron_SELECTION_annotated[ intron_SELECTION ]


nb_introns = sum(intron_SELECTION)

n1_tab = data.frame(matrix(ncol = 0, nrow = nb_introns))
n2_tab = data.frame(matrix(ncol = 0, nrow = nb_introns))
n3_tab = data.frame(matrix(ncol = 0, nrow = nb_introns))

coverage = data.frame(matrix(ncol = 0, nrow = sum(by_gene$busco_id)))


for ( RNAsample in RNAseq_table$SRA_accession_ID ){
  print(RNAsample)
  
  by_gene = read.delim(paste(pathData , "Analyses-RNAseq/",species,"/",RNAsample,"/coverage_fpkm_by_genes.tab",sep=""))
  by_gene = by_gene[by_gene$gene_id %in% busco_to_gene$gene_id,]
  coverage[RNAsample] = by_gene$exon_coverage
  
  
  summaryTable = read.delim(paste(pathData,"Analyses-RNAseq/",species,"/",RNAsample,"/n1_n2_n3_by_intron.tab",sep=""))
  summaryTable = summaryTable[intron_SELECTION,]
  
  n1_tab[RNAsample] = summaryTable$n1
  n2_tab[RNAsample] = summaryTable$n2_spl5 + summaryTable$n2_spl3
  n3_tab[RNAsample] = summaryTable$n3_spl5 + summaryTable$n3_spl3
  
}

if (length( RNAseq_table$SRA_accession_ID ) >= 20 ){ nb = 20 } else { nb = length( RNAseq_table$SRA_accession_ID )}
df = data.frame()
list_major_intron = list()

for ( i in 1:nb ){ # nombre de RNA-seqs compilés
  no_replicate = 30
  if ( i == 1 ){
    no_replicate = nrow( RNAseq_table )
  }
  
  for ( j in 1:no_replicate ){ # nombre de rééchantillonnages
    if ( i == 1 ){
      intervalle = RNAseq_table$SRA_accession_ID[j]
      
      print(intervalle)
      n1_value = n1_tab[,intervalle]
      n2_value = n2_tab[,intervalle]
      n3_value = n3_tab[,intervalle]
      cov = coverage[,intervalle]
    } else {
      intervalle = sample( RNAseq_table$SRA_accession_ID ,i )
      
      print(intervalle)
      n1_value = apply(n1_tab[,intervalle],1,function(x) sum(x,na.rm=T))
      n2_value = apply(n2_tab[,intervalle],1,function(x) sum(x,na.rm=T))
      n3_value = apply(n3_tab[,intervalle],1,function(x) sum(x,na.rm=T))
      cov = apply(coverage[,intervalle],1,function(x) sum(x,na.rm=T))
    }
    
    SVR = n2_value / (n1_value + n2_value)
    SVR[(n1_value + n2_value) < 10] = NA
    
    IR = n3_value / (2 * n1_value + n3_value)
    IR[ ( n1_value + n3_value/2) < 10] = NA
    
    major =  sum( n1_value != 0 & SVR < 0.5 & IR < 0.5  ,na.rm = T)  
    
    minor = sum(n1_value != 0 & (SVR >= 0.5 | IR >= 0.5) ,na.rm = T)
    
    unclassified = sum( n1_value != 0 ) - (minor + major) 
    
    annot_N1N2_sup10 =  sum( n1_value != 0 & (n1_value + n2_value) >= 10 & intron_SELECTION_annotated )

    
    
    
    
    list_major_intron[[paste(i,j)]] = which(!is.na(SVR))
    
    da = data.frame(
      annot_N1N2_sup10,
      annotated_intron = sum(intron_SELECTION_annotated),
      sequencing_depth = median(cov,na.rm = T),
      average_svr_busco = mean( SVR[n1_value != 0 & SVR < 0.5 & IR < 0.5],na.rm=T ),
#       average_svr_busco = mean( SVR[n1_value != 0 & SVR < 0.5 & IR < 0.5 & intron_SELECTION_annotated],na.rm=T ),
      no_major_intron = major,
      no_compiled_rna_seq = i,
      replicate = j,
      echantillon = "all introns",
      selection = 0
    )
    da$samples = list(intervalle)
    df = rbind(df,da)
  }
}















# if ( length(unlist(list_major_intron)) != 0 ){
#   
#   vector_count = data.frame(intron_position = unique(unlist(list_major_intron)))
#   rownames(vector_count) = vector_count$intron_position
#   vector_count$rep = 0
#   
#   lapply(list_major_intron,function(x) {
#     vector_count[as.character(x),"rep"] <<- vector_count[as.character(x),"rep"] + 1
#     return(NA)
#   }
#   )
#   
#   
#   selected = vector_count[vector_count$rep > nrow(df) * .9,]$intron_position
#   length(selected)
#   
#   
#   previous_df = df
#   
#   for (i in 1:nb){ # nombre de RNA-seqs compilés
#     no_replicate = 10
#     if (i == 1){no_replicate = length( RNAseq_table$SRA_accession_ID )}
#     
#     for (j in 1:no_replicate){ # nombre de rééchantillonnages
#       if (i == 1){
#         intervalle = previous_df[previous_df$no_compiled_rna_seq == i & previous_df$replicate == j,]$samples[[1]]
#         
#         print(intervalle)
#         n1_value = n1_tab[,intervalle]
#         n2_value = n2_tab[,intervalle]
#         n3_value = n3_tab[,intervalle]
#         cov = coverage[,intervalle]
#       }    else {
#         intervalle=previous_df[previous_df$no_compiled_rna_seq == i & previous_df$replicate == j,]$samples[[1]]
#         
#         print(intervalle)
#         n1_value = apply(n1_tab[,intervalle],1,function(x) sum(x,na.rm=T))
#         n2_value = apply(n2_tab[,intervalle],1,function(x) sum(x,na.rm=T))
#         n3_value = apply(n3_tab[,intervalle],1,function(x) sum(x,na.rm=T))
#         cov = apply(coverage[,intervalle],1,function(x) sum(x,na.rm=T))
#       }
#       
#       svr = n2_value / (n1_value+n2_value)
#       ir = n3_value / (2*n1_value+n3_value)
#       n2_value[svr >= 0.5 | ir >= 0.5 | n1_value < 10]=NA
#       n1_value[svr >= 0.5 | ir >= 0.5 | n1_value < 10]=NA
#       
#       
#       svr = n2_value / (n1_value + n2_value)
#       print(table(is.na(svr)))
#       
#       svr = svr[selected]
#       
#       da = data.frame(
#         N1_sup10 = NA,
#         annot_N1N2_sup10 = NA,
#         minor_N1_sup10 = NA,
#         annotated_intron = NA,
#         annot_N1_sup10 = NA,
#         annot_minor_N1_sup10 = NA,
#         sequencing_depth = median(cov,na.rm = T),
#         average_svr_busco = mean(svr,na.rm=T),
#         no_major_intron = length(svr[!is.na(svr)]),
#         no_compiled_rna_seq = i,
#         replicate = j,
#         echantillon = "introns shared across datasets",
#         selection = length(selected)
#       )
#       
#       da$samples = list(intervalle)
#       
#       df = rbind(df,da)
#     }
#   }
# }
df$samples = unlist(lapply(df$samples,function(x) paste(x,collapse=",")))


write.table(df, file=path_output,row.names=F, col.names=T, sep="\t", quote=F)
