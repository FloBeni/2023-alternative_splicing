# Generate Data 1

options(stringsAsFactors = F, scipen = 999)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)

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
# pathData="/beegfs/data/XXXXX/Projet-SplicedVariants/"

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_species.xls",sep=""))


sp_studied = c()
for (species in names(mysheets)){
  if (!is.na(mysheets[[species]]$Group_study[1]) & (mysheets[[species]]$Group_study[1] == "53_sp" | mysheets[[species]]$Group_study[1] == "69_sp")){
    sp_studied = append(sp_studied,species)
  }
}


species = "Drosophila_melanogaster"
for (species in sp_studied ){
  print(species)
  
  ## INTRON
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  
  fpkm_cov = read.delim(paste(pathData,"/Analyses/",species,"/by_gene_analysis.tab",sep="") , header=T , sep="\t",comment.char = "#")
  fpkm_cov$busco_metazoa = fpkm_cov$gene_id %in% busco_tab$gene_id
  rownames(fpkm_cov) = fpkm_cov$gene_id
  
  by_intron = read.delim(file=paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  by_intron$id = paste(by_intron$seqname,by_intron$gene_id,by_intron$splice5,by_intron$splice3,by_intron$strand,sep=";")
  by_intron$fpkm = fpkm_cov[by_intron$gene_id,]$weighted_fpkm
  
  IntronLibrary = read.delim(paste(pathData,"/Analyses/",species,"/IntronLibrary_inclusive.txt",sep=""))
  IntronLibrary = IntronLibrary %>%
    mutate(Gene = strsplit(as.character(Gene), ",")) %>%
    unnest(Gene) %>%
    filter(Gene != "")
  IntronLibrary$id = paste(IntronLibrary$Chr,IntronLibrary$Gene,IntronLibrary$Splice5,IntronLibrary$Splice3,IntronLibrary$Strand,sep=";")
  IntronLibrary$Annotation = grepl("Annotation",IntronLibrary$Source)
  rownames(IntronLibrary) = IntronLibrary$id
  by_intron$Annotation = IntronLibrary[by_intron$id,]$Annotation
  by_intron$splicesite = paste(IntronLibrary[by_intron$id,]$SpliceSignal5,IntronLibrary[by_intron$id,]$SpliceSignal3)
  
  
  minor_introns = read.delim(file=paste(pathData,"Analyses/",species,"/by_minor_intron.tab",sep=""), header=T , sep="\t",comment.char = "#")
  minor_introns$id = paste(minor_introns$seqname,minor_introns$gene_id,minor_introns$splice5,minor_introns$splice3,minor_introns$strand,sep=";")
  rownames(minor_introns) = minor_introns$id
  by_intron$mira = minor_introns[by_intron$id,]$mira
  
  major_introns = read.delim(paste(pathData,"/Analyses/",species,"/by_intron_major_overlap.tab",sep=""))
  major_introns$id = paste(major_introns$seqname,major_introns$gene_id,major_introns$splice5,major_introns$splice3,major_introns$strand,sep=";")
  rownames(major_introns) = major_introns$id
  by_intron$have_abundant_sv = major_introns[by_intron$id,]$have_abundant_sv
  
}




### Analyses
fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]

all_major = by_intron[by_intron$intron_class == "major" & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]
minor_introns = by_intron[by_intron$intron_class == "minor" & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]


fpkm_cov = fpkm_cov[fpkm_cov$busco_metazoa ,]
CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco


busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                        !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]

busco_major = all_major[all_major$gene_id %in% busco_tab$gene_id,]
minor_introns_busco = minor_introns[minor_introns$gene_id %in% busco_tab$gene_id,]
all_intron_busco = by_intron[by_intron$gene_id %in% busco_tab$gene_id & by_intron$into_cds == "True"  ,]




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
  prop_analyzable_proteincoding = sum(by_intron$Annotation & (by_intron$n1 + by_intron$n2_spl3 + by_intron$n2_spl5) >= 10 ) / sum(by_intron$Annotation),
  
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
