
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

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_v53.xls",sep=""))


sp_studied = c()
for (species in names(mysheets)){
  if (!is.na(mysheets[[species]]$Group_study[1]) & (mysheets[[species]]$Group_study[1] == "53_sp" | mysheets[[species]]$Group_study[1] == "69_sp")){
    sp_studied = append(sp_studied,species)
  }
}

data_dNdS_subset_200k_v2 = read.delim(paste(pathData , "DnDs/Metazoa_species_filtered_v53.2/subset_200_ksites_GC3_v2/data_calculation.tab",sep=""))

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


data_1 = data.frame()
for (species in sp_studied){
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
  
  
  by_intron = read.delim(paste(pathData,"per_species/",species,".tab.gz",sep=""),  sep="\t")
  fpkm_cov = read.delim(paste(pathData,"per_species/",species,"fpkm.tab.gz",sep=""),  sep="\t")
  
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  
  major_introns = by_intron[by_intron$intron_class == "major" & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]
  minor_introns = by_intron[by_intron$intron_class == "minor" & !is.na(by_intron$mira) & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]
  
  fpkm_cov = fpkm_cov[fpkm_cov$busco_metazoa ,]
  CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco
  
  major_introns_busco = major_introns[major_introns$gene_id %in% fpkm_cov$gene_id,]
  minor_introns_busco = minor_introns[minor_introns$gene_id %in% fpkm_cov$gene_id,]
  all_intron = by_intron[ by_intron$into_cds == "True" ,]
  all_intron_busco = by_intron[by_intron$gene_id %in% fpkm_cov$gene_id & by_intron$into_cds == "True"  ,]
  
  
  
  prop_major_nt_sup100 =  sum(major_introns$have_abundant_sv == "False" & (major_introns$n1 + major_introns$n2_spl3 + major_introns$n2_spl5) >= 100) / nrow(major_introns)
  
  intron_svr_fpkm = major_introns[major_introns$have_abundant_sv == "False" & (major_introns$n1 + major_introns$n2_spl3 + major_introns$n2_spl5) >= 100,]
  if (nrow(intron_svr_fpkm) != 0){
    rlowSVR_slope_low_as = get_rsquared_slope(prop.quantile = 0.05, Xaxis=intron_svr_fpkm$fpkm, Yaxis= unlist(intron_svr_fpkm$splice_variant_rate))
  } else { rlowSVR_slope_low_as=NA}
  intron_svr_fpkm = major_introns[(major_introns$n1 + major_introns$n2_spl3 + major_introns$n2_spl5) >= 100,]
  if (nrow(intron_svr_fpkm) != 0){
    rlowSVR_slope = get_rsquared_slope(prop.quantile = 0.05, Xaxis=intron_svr_fpkm$fpkm, Yaxis= unlist(intron_svr_fpkm$splice_variant_rate))
  } else { rlowSVR_slope=NA}
  
  gene_n1 =  tapply(major_introns$n1,major_introns$gene_id,sum)
  gene_n2 =   tapply(major_introns$n2_spl3,major_introns$gene_id,sum) +  tapply(major_introns$n2_spl5,major_introns$gene_id,sum)
  gene_as_proteincoding_all = 1 - (1 - (gene_n2 / (gene_n1 + gene_n2)))^(tapply(major_introns$n1,major_introns$gene_id,length))
  
  gene_n1 =  tapply(major_introns[major_introns$have_abundant_sv == "False",]$n1,major_introns[major_introns$have_abundant_sv == "False",]$gene_id,sum)
  gene_n2 =   tapply(major_introns[major_introns$have_abundant_sv == "False",]$n2_spl3,major_introns[major_introns$have_abundant_sv == "False",]$gene_id,sum) + tapply(major_introns[major_introns$have_abundant_sv == "False",]$n2_spl5,major_introns[major_introns$have_abundant_sv == "False",]$gene_id,sum)
  gene_as_proteincoding_lowas = 1 - (1 - ( gene_n2 / ( gene_n1 + gene_n2 )))^(tapply(major_introns$n1,major_introns$gene_id,length)[names(gene_n2)])
  
  
  gene_n1 =  tapply(major_introns_busco$n1,major_introns_busco$gene_id,sum)
  gene_n2 =   tapply(major_introns_busco$n2_spl3,major_introns_busco$gene_id,sum) +  tapply(major_introns_busco$n2_spl5,major_introns_busco$gene_id,sum)
  gene_as_busco_all = 1 - (1 - (gene_n2 / (gene_n1 + gene_n2)))^(tapply(major_introns_busco$n1,major_introns_busco$gene_id,length))
  
  gene_n1 =  tapply(major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$n1,major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$gene_id,sum)
  gene_n2 =   tapply(major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$n2_spl3,major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$gene_id,sum) + tapply(major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$n2_spl5,major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$gene_id,sum)
  gene_as_busco_lowas = 1 - (1 - ( gene_n2 / ( gene_n1 + gene_n2 )))^(tapply(major_introns_busco$n1,major_introns_busco$gene_id,length)[names(gene_n2)])
  
  
  
  data_1 = rbind(data_1, data.frame(
    species ,
    genome_assembly,
    clade,
    body_size,
    longevity,
    nb_rnaseq,
    list_rnaseq,
    dNdS_200k = get_CM_dNdS( data_dNdS_subset_200k_v2[data_dNdS_subset_200k_v2$species == species,]),
    
    nb_busco = sum(fpkm_cov$busco_metazoa),
    CoverageBuscoExon,
    
    annotated_intron_proteincoding = sum(all_intron$Annotation),
    analyzable_intron_proteincoding = sum(all_intron$Annotation & (all_intron$n1 + all_intron$n2_spl3 + all_intron$n2_spl5) >= 10 ),
    annotated_intron_busco = sum(all_intron_busco$Annotation),
    analyzable_intron_busco = sum(all_intron_busco$Annotation & (all_intron_busco$n1 + all_intron_busco$n2_spl3 + all_intron_busco$n2_spl5) >= 10 ),
    
    splsite_gtag_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("GT AG")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_gtag_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("GT AG")) / sum(all_intron_busco$intron_class == "major"),
    splsite_gcag_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("GC AG")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_gcag_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("GC AG")) / sum(all_intron_busco$intron_class == "major"),
    splsite_atac_minor_busco = sum(all_intron_busco[all_intron_busco$intron_class == "minor",]$splicesite %in% c("AT AC")) / sum(all_intron_busco$intron_class == "minor"),
    splsite_atac_major_busco = sum(all_intron_busco[all_intron_busco$intron_class == "major",]$splicesite %in% c("AT AC")) / sum(all_intron_busco$intron_class == "major"),
    
    prop_major_sv_busco = sum((major_introns_busco$n2_spl3 + major_introns_busco$n2_spl5) > 0) /  nrow(major_introns_busco),
    prop_major_sv_proteincoding = sum((major_introns$n2_spl3 + major_introns$n2_spl5) > 0) /  nrow(major_introns),
    
    mean_gene_busco_as_lowas = mean(gene_as_busco_lowas),
    mean_gene_busco_as = mean(gene_as_busco_all),
    mean_gene_proteincoding_as_lowas = mean(gene_as_proteincoding_lowas),
    mean_gene_proteincoding_as = mean(gene_as_proteincoding_all),
    
    
    mean_as_proteincoding = mean( major_introns$splice_variant_rate ) * 100 ,
    mean_as_proteincoding_high_as = mean(major_introns[major_introns$have_abundant_sv == "True"  ,]$splice_variant_rate)*100,
    mean_as_proteincoding_low_as = mean(major_introns[major_introns$have_abundant_sv == "False",]$splice_variant_rate)*100,
    
    mean_as_busco = mean( major_introns_busco$splice_variant_rate ) * 100 ,
    mean_as_busco_high_as = mean(major_introns_busco[major_introns_busco$have_abundant_sv == "True"  ,]$splice_variant_rate)*100,
    mean_as_busco_low_as = mean(major_introns_busco[major_introns_busco$have_abundant_sv == "False",]$splice_variant_rate)*100,
    
    prop_rare_sv_protein_coding = sum( minor_introns$mira <= 0.05 ) / nrow(minor_introns) * 100,
    prop_rare_sv_busco = sum( minor_introns_busco$mira <= 0.05 ) / nrow(minor_introns_busco) * 100,
    
    mean_intron_per_gene_busco = mean(table(major_introns_busco$gene_id)),
    mean_intron_per_gene_proteincoding = mean(table(major_introns$gene_id)),
    
    prop_major_nt_sup100,
    cor_as_fpkm_all_as = rlowSVR_slope[1],
    pval_cor_as_fpkm_all_as = rlowSVR_slope[2],
    cor_as_fpkm_low_as = rlowSVR_slope_low_as[1],
    pval_cor_as_fpkm_low_as = rlowSVR_slope_low_as[2]
  ))
}

# write.table(data_1 , paste("data/Data1_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

