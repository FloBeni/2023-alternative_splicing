# Generate Data 1

options(stringsAsFactors = F, scipen = 999)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)

get_rsquared_slope = function(prop.quantile = 0.1,Xaxis,Yaxis){
  quantile = unique(quantile(Xaxis, probs = seq(0, 1,prop.quantile),na.rm=T))
  intervalle = cut(Xaxis, quantile,include.lowest = T,include.higher=T)
  X = tapply(Xaxis, intervalle, median)
  if ( !any(is.na(X)) ){
    Y = tapply(Yaxis, intervalle, mean)
    X = log10(X)
    pearson_method = cor.test(X, Y,method="pearson")
    cor = pearson_method$estimate
    pval_cor = pearson_method$p.value
    return( c(cor,pval_cor) )
  } else {
    return(
      c(NA,NA)
    )
  }
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


pathData="/home/fbenitiere/data//Projet-SplicedVariants/"
pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_species.xls",sep=""))


sp_studied = c()
for (species in names(mysheets)){
  if (!is.na(mysheets[[species]]$Group_study[1]) & (mysheets[[species]]$Group_study[1] == "53_sp" | mysheets[[species]]$Group_study[1] == "69_sp")){
    sp_studied = append(sp_studied,species)
  }
}

data_dNdS_subset_200k_v2 = read.delim(paste(pathData,"DnDs/Metazoa_species_filtered_v53.2/subset_200_ksites_GC3_v2/data_calculation.tab",sep=""))

get_CM_dNdS<-function(D) {
  # Compute the cumulated number of substitutions over all genes
  cum_KS=sum(D$num_dS)
  cum_KN=sum(D$num_dN)
  cum_OS=sum(D$den_dS/D$branch_length)
  cum_ON=sum(D$den_dS/D$branch_length)
  # Compute cumulated DN, DS
  cum_dS = cum_KS/cum_OS
  cum_dN = cum_KN/cum_ON
  # Compute cumulated dN/dS
  cum_dNdS=cum_dN/cum_dS
  
  return(cum_dNdS)
}


all_data = data.frame()
species = "Drosophila_melanogaster"
for (species in sp_studied ){
  print(species)
  con <- file(paste(pathData,"/Annotations/",species,"/data_source/annotation.gff",sep=""),"r")
  first_line <- readLines(con,n=5)
  close(con)
  print(first_line[5])
  genome_assembly = first_line[5]
  genome_assembly = str_replace(genome_assembly,"#!genome-build-accession NCBI_Assembly:","")
  
  
  clade = mysheets[[species]]$Clade[1]
  body_size = mysheets[[species]]$`Length (cm)`[1]
  longevity = mysheets[[species]]$`Longevity (days)`[1]
  bioproj = read.delim(paste(pathData,"RNAseq_table/",species,"/list_Acc.tab",sep="" ) )
  nb_rnaseq = nrow(bioproj)
  list_rnaseq = paste(bioproj$SRA_accession_ID ,collapse = ";")
  
  
  ## INTRON
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  
  fpkm_cov = read.delim(paste(pathData,"/Analyses/",species,"/by_gene_analysis.tab",sep="") , header=T , sep="\t",comment.char = "#")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  rownames(fpkm_cov) = fpkm_cov$gene_id
  
  
  minor_introns = read.delim(file=paste(pathData,"Analyses/",species,"/by_minor_intron.tab",sep=""))
  minor_introns = minor_introns[minor_introns$gene_id %in% fpkm_cov$gene_id,]
  print(table(minor_introns$into_cds))
  minor_introns = minor_introns[minor_introns$into_cds == "True",]
  # minor_introns = minor_introns[minor_introns$criptic_intron == "False",]
  
  print(table(minor_introns$criptic_intron))
  
  all_major = read.delim(paste(pathData,"/Analyses/",species,"/by_intron_major_overlap.tab",sep=""))
  all_major = all_major[all_major$intron_class == "major" & all_major$into_cds == "True" & all_major$gene_id %in% fpkm_cov$gene_id,]
  all_major$fpkm = fpkm_cov[all_major$gene_id,]$weighted_fpkm
  
  
  fpkm_cov = fpkm_cov[fpkm_cov$gene_id %in% busco_tab$gene_id ,]
  CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco
  
  
  
  busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                          !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
  
  busco_major = all_major[all_major$gene_id %in% busco_tab$gene_id,]
  minor_introns_busco = minor_introns[minor_introns$gene_id %in% busco_tab$gene_id,]
  
  
  prop_major_nt_sup100 =  sum(all_major$have_abundant_sv == "False" & (all_major$n1 + all_major$n2_spl3 + all_major$n2_spl5) >= 100) / nrow(all_major)
  
  intron_svr_fpkm = all_major[all_major$have_abundant_sv == "False" & (all_major$n1 + all_major$n2_spl3 + all_major$n2_spl5) >= 100,]
  if (nrow(intron_svr_fpkm) != 0){
    rlowSVR_slope_low_as = get_rsquared_slope(prop.quantile = 0.05, Xaxis=intron_svr_fpkm$fpkm, Yaxis= unlist(intron_svr_fpkm$splice_variant_rate))
  } else { rlowSVR_slope_low_as=NA}
  intron_svr_fpkm = all_major[(all_major$n1 + all_major$n2_spl3 + all_major$n2_spl5) >= 100,]
  if (nrow(intron_svr_fpkm) != 0){
    rlowSVR_slope = get_rsquared_slope(prop.quantile = 0.05, Xaxis=intron_svr_fpkm$fpkm, Yaxis= unlist(intron_svr_fpkm$splice_variant_rate))
  } else { rlowSVR_slope=NA}
  
  gene_n1 =  tapply(all_major$n1,all_major$gene_id,sum)
  gene_n2 =   tapply(all_major$n2_spl3,all_major$gene_id,sum) +  tapply(all_major$n2_spl5,all_major$gene_id,sum)
  gene_as_proteincoding_all = 1 - (1 - (gene_n2 / (gene_n1 + gene_n2)))^(tapply(all_major$n1,all_major$gene_id,length))
  
  gene_n1 =  tapply(all_major[all_major$have_abundant_sv == "False",]$n1,all_major[all_major$have_abundant_sv == "False",]$gene_id,sum)
  gene_n2 =   tapply(all_major[all_major$have_abundant_sv == "False",]$n2_spl3,all_major[all_major$have_abundant_sv == "False",]$gene_id,sum) + tapply(all_major[all_major$have_abundant_sv == "False",]$n2_spl5,all_major[all_major$have_abundant_sv == "False",]$gene_id,sum)
  gene_as_proteincoding_lowas = 1 - (1 - ( gene_n2 / ( gene_n1 + gene_n2 )))^(tapply(all_major$n1,all_major$gene_id,length)[names(gene_n2)])
  
  
  gene_n1 =  tapply(busco_major$n1,busco_major$gene_id,sum)
  gene_n2 =   tapply(busco_major$n2_spl3,busco_major$gene_id,sum) +  tapply(busco_major$n2_spl5,busco_major$gene_id,sum)
  gene_as_busco_all = 1 - (1 - (gene_n2 / (gene_n1 + gene_n2)))^(tapply(busco_major$n1,busco_major$gene_id,length))
  
  gene_n1 =  tapply(busco_major[busco_major$have_abundant_sv == "False",]$n1,busco_major[busco_major$have_abundant_sv == "False",]$gene_id,sum)
  gene_n2 =   tapply(busco_major[busco_major$have_abundant_sv == "False",]$n2_spl3,busco_major[busco_major$have_abundant_sv == "False",]$gene_id,sum) + tapply(busco_major[busco_major$have_abundant_sv == "False",]$n2_spl5,busco_major[busco_major$have_abundant_sv == "False",]$gene_id,sum)
  gene_as_busco_lowas = 1 - (1 - ( gene_n2 / ( gene_n1 + gene_n2 )))^(tapply(busco_major$n1,busco_major$gene_id,length)[names(gene_n2)])
  
  
  
  all_intron = read.delim(paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  all_intron = all_intron[all_intron$into_cds == "True"  ,]
  
  IntronLibrary = read.delim(paste(pathData,"/Analyses/",species,"/IntronLibrary_inclusive.txt",sep=""))
  IntronLibrary = IntronLibrary %>%
    mutate(Gene = strsplit(as.character(Gene), ",")) %>%
    unnest(Gene) %>%
    filter(Gene != "")
  IntronLibrary$id = paste(IntronLibrary$Chr,IntronLibrary$Gene,IntronLibrary$Splice5,IntronLibrary$Splice3,IntronLibrary$Strand,sep=";")
  IntronLibrary$Annotation = grepl("Annotation",IntronLibrary$Source)
  rownames(IntronLibrary) = IntronLibrary$id
  
  all_intron$id = paste(all_intron$seqname,all_intron$gene_id,all_intron$splice5,all_intron$splice3,all_intron$strand,sep=";")
  
  all_intron$SpliceSignal3 = IntronLibrary[all_intron$id,]$SpliceSignal3
  all_intron$SpliceSignal5 = IntronLibrary[all_intron$id,]$SpliceSignal5
  all_intron$Annotation = IntronLibrary[all_intron$id,]$Annotation
  all_intron$splicesite = paste(all_intron$SpliceSignal5,all_intron$SpliceSignal3)
  
  all_intron_busco = all_intron[all_intron$gene_id %in% busco_tab$gene_id  ,]
  
  
  
  all_data = rbind(all_data, data.frame(
    species ,
    genome_assembly,
    clade,
    body_size,
    longevity,
    nb_rnaseq,
    list_rnaseq,
    dNdS_200k = get_CM_dNdS( data_dNdS_subset_200k_v2[data_dNdS_subset_200k_v2$species == species,]),
    
    nb_busco = nrow(busco_tab),
    CoverageBuscoExon,
    
    prop_analyzable_busco = sum(all_intron_busco$Annotation & (all_intron_busco$n1 + all_intron_busco$n2_spl3 + all_intron_busco$n2_spl5) >= 10 ) / sum(all_intron_busco$Annotation),
    prop_analyzable_proteincoding = sum(all_intron$Annotation & (all_intron$n1 + all_intron$n2_spl3 + all_intron$n2_spl5) >= 10 ) / sum(all_intron$Annotation,na.rm = T),
    
    splsite_gtag_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("GT AG")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_gtag_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("GT AG")) / sum(all_intron_busco$intron_class == "major"),
    splsite_gcag_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("GC AG")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_gcag_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("GC AG")) / sum(all_intron_busco$intron_class == "major"),
    splsite_atac_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("AT AC")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_atac_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("AT AC")) / sum(all_intron_busco$intron_class == "major"),
    
    prop_major_sv_busco = sum((busco_major$n2_spl3 + busco_major$n2_spl5) > 0) /  nrow(busco_major),
    prop_major_sv_proteincoding = sum((all_major$n2_spl3 + all_major$n2_spl5) > 0) /  nrow(all_major),
    
    mean_gene_busco_as_lowas = mean(gene_as_busco_lowas),
    mean_gene_busco_as = mean(gene_as_busco_all),
    mean_gene_proteincoding_as_lowas = mean(gene_as_proteincoding_lowas),
    mean_gene_proteincoding_as = mean(gene_as_proteincoding_all),
    
    
    mean_as_proteincoding = mean( all_major$splice_variant_rate ) * 100 ,
    mean_as_proteincoding_high_as = mean(all_major[all_major$have_abundant_sv == "True"  ,]$splice_variant_rate)*100,
    mean_as_proteincoding_low_as = mean(all_major[all_major$have_abundant_sv == "False",]$splice_variant_rate)*100,
    
    mean_as_busco = mean( busco_major$splice_variant_rate ) * 100 ,
    mean_as_busco_high_as = mean(busco_major[busco_major$have_abundant_sv == "True"  ,]$splice_variant_rate)*100,
    mean_as_busco_low_as = mean(busco_major[busco_major$have_abundant_sv == "False",]$splice_variant_rate)*100,
    
    prop_rare_sv_protein_coding = sum( minor_introns$mira <= 0.05 ) / nrow(minor_introns) * 100,
    prop_rare_sv_busco = sum( minor_introns_busco$mira <= 0.05 ) / nrow(minor_introns_busco) * 100,
    
    mean_intron_per_gene_busco = mean(table(busco_major$gene_id)),
    mean_intron_per_gene_proteincoding = mean(table(all_major$gene_id)),
    
    prop_major_nt_sup100,
    cor_as_fpkm_all_as = rlowSVR_slope[1],
    pval_cor_as_fpkm_all_as = rlowSVR_slope[2],
    cor_as_fpkm_low_as = rlowSVR_slope_low_as[1],
    pval_cor_as_fpkm_low_as = rlowSVR_slope_low_as[2]
  ))
}

write.table(all_data,paste("data/Data1_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)



