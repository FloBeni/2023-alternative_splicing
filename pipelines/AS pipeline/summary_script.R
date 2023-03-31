# Create a summary for a species
options(stringsAsFactors = F, scipen = 999)

library(rgbif)
library(ape)
library(ggpubr)
library(ggrepel)
library(readxl)
library(seqinr)
library(stringr)
library(dplyr)
library(tidyr)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # fonction de lecture de fichier excel
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

add_charac <- function(data_summary,label,description,value){
  data_summary = rbind(data_summary,data.frame(
    label = label,
    description = description,
    value=value
  ))
  return(data_summary)
}



get_rsquared_slope = function(prop.quantile = 0.2,Xaxis,Yaxis){
  print(Xaxis)
  print(Yaxis)
  if (length(Xaxis) != 0){ # Pas de vecteur vide
    quantile = unique(quantile(Xaxis, probs = seq(0, 1,prop.quantile),na.rm=T))
    if (length(quantile) > 2){ # au moins un intervalle
      intervalle = cut(Xaxis, quantile,include.lowest = T,include.higher=T)
      X = tapply(Xaxis, intervalle,  function(x) mean(x, na.rm=T))
      Y = tapply(Yaxis, intervalle, function(x) mean(x, na.rm=T))
      print(Y)
      print(X)
      formule = lm( formula = Y ~ X)
      R.signed = summary(formule)$coefficients["X","Estimate"] / abs(summary(formule)$coefficients["X","Estimate"])*round(summary(formule)$adj.r.squared, 3)
      slope = summary(formule)$coefficients["X","Estimate"]

      return(c(R.signed,slope))
    } else {
      return(c(NA,NA))
    }} else {
      return(c(NA,NA))
    }
}

args = (commandArgs(TRUE))
pathData = args[1]
species = args[2]
excel_input = args[3]
gff_path = args[4]
gc_table_path = args[5]
library_path = args[6]
by_intron_cds_path = args[7]
by_gene_analysis_path = args[8]
by_intron_overlap_path = args[9]
prot_path = args[10]
sumv3_path = args[11]


mysheets <- read_excel_allsheets(excel_input)


if (species == "Musca_domestica") {
  name_sp = "Mdim"
} else {
  name_sp = paste(substr(strsplit(species,"_")[[1]][1],1,1),
                  substr(strsplit(species,"_")[[1]][2],1,3),sep="")
}


data_summary = data.frame()
data_summary = add_charac(data_summary,'species',"",species)
data_summary = add_charac(data_summary,'short_name',"",name_sp)

con <- file(gff_path,"r")
first_line <- readLines(con,n=10)
close(con)
print(first_line[5])
genome_assembly = first_line[5]
genome_assembly = str_replace(genome_assembly,"#!genome-build-accession NCBI_Assembly:","")
data_summary = add_charac(data_summary,'genome_assembly',"",genome_assembly)

data_summary = add_charac(data_summary,'clade;qual',"",mysheets[[species]]$Clade[1])
data_summary = add_charac(data_summary,'socialite;qual',"",mysheets[[species]]$`Socialite`[1])
data_summary = add_charac(data_summary,'longevity;quant',"",mysheets[[species]]$`Longevity (days)`[1])
data_summary = add_charac(data_summary,'length;quant',"",mysheets[[species]]$`Length (cm)`[1])

key = name_backbone(name=str_replace(species,"_"," "),rank="species")$usageKey
RGBIF_count = occ_search(key,limit=0)$meta$count

data_summary = add_charac(data_summary,'RGBIF_observation;quant',"",RGBIF_count)


genome_character = read.table(gc_table_path,header=T)
genome_size_var =  sum(genome_character[genome_character$genome_character %in% c("A","T","G","C"), ]$Freq) / 1000000
GC_content_proportion_var = round( sum(genome_character[genome_character$genome_character %in% c("G","C"), ]$Freq) /
                                     sum(genome_character[genome_character$genome_character %in% c("A","T","G","C"), ]$Freq) ,2)

data_summary = add_charac(data_summary,'genome_size;quant',"Mb",genome_size_var)
data_summary = add_charac(data_summary,'GC_content_proportion;quant',"",GC_content_proportion_var)
data_summary = add_charac(data_summary,'No_prot_annot;quant',"",length(read.fasta(prot_path)))

print("la1")
#### Data_set
IntronLibrary = read.delim(library_path)
IntronLibrary = IntronLibrary %>%
  mutate(Gene = strsplit(as.character(Gene), ",")) %>%
  unnest(Gene) %>%
  filter(Gene != "")

IntronLibrary = IntronLibrary[grepl("Annotation",IntronLibrary$Source),]
IntronLibrary$id = paste(IntronLibrary$Chr,IntronLibrary$Strand,IntronLibrary$Splice5,IntronLibrary$Splice3,IntronLibrary$Gene,sep=";")

print("la1")
by_intron_overlap = read.delim( by_intron_overlap_path, header=T , sep="\t",comment.char = "#")
by_intron_overlap$id = paste(by_intron_overlap$seqname,by_intron_overlap$strand,by_intron_overlap$splice5,by_intron_overlap$splice3,by_intron_overlap$gene_id,sep=";")
rownames(by_intron_overlap) = by_intron_overlap$id

print("la2")
by_intron =  read.delim( by_intron_cds_path, header=T , sep="\t",comment.char = "#")
by_intron$id = paste(by_intron$seqname,by_intron$strand,by_intron$splice5,by_intron$splice3,by_intron$gene_id,sep=";")
by_intron$isAnnotated = by_intron$id %in% IntronLibrary$id
by_intron$length = abs(by_intron$splice3 - by_intron$splice5)
by_intron$overlap = by_intron_overlap[by_intron$id,"overlap"]

print("la3")
by_gene =  read.delim(by_gene_analysis_path , header=T , sep="\t",comment.char = "#")
by_gene = by_gene[by_gene$type == "gene" & grepl("gene_biotype=protein_coding" , by_gene$attributes),]
rownames(by_gene) = by_gene$gene_id
by_gene = by_gene[by_gene$type == "gene",] # FILTRE PSEUDOGENE
by_intron = by_intron[by_intron$gene_id %in% by_gene$gene_id,] # FILTRE PSEUDOGENE
by_intron$median_fpkm = by_gene[by_intron$gene_id,]$median_fpkm

print("la4")
for (busco_group in c("metazoa","embryophyta","eukaryota","None")){ print(busco_group)
  can_analyse = T
  if ( busco_group != "None" ){
    if (file.exists(paste(pathData , "Annotations/",species,"/busco_analysis/busco_to_gene_id_",busco_group,sep=""))){
      busco_to_gene = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_",busco_group,sep=""))

      busco_to_gene = busco_to_gene[!(duplicated(busco_to_gene$busco_id,fromLast = FALSE) | duplicated(busco_to_gene$busco_id,fromLast = TRUE)) &
                          !(duplicated(busco_to_gene$gene_id,fromLast = FALSE) | duplicated(busco_to_gene$gene_id,fromLast = TRUE)) ,]

      rownames(busco_to_gene) = busco_to_gene$gene_id
      by_gene$busco_id = by_gene$gene_id %in% busco_to_gene$gene_id
      by_intron$busco_id = by_intron$gene_id %in% busco_to_gene$gene_id

      by_gene_selected = by_gene[by_gene$busco_id,]
      by_intron_selected = by_intron[by_intron$busco_id & by_intron$intron_class == "major" & by_intron$into_cds == "True",]
      print(nrow(by_intron_selected))

    } else { can_analyse = F }
  } else {
    by_gene_selected = by_gene
    by_intron_selected = by_intron[ by_intron$intron_class == "major" & by_intron$into_cds == "True",]
  }

  if (can_analyse){
    data_summary = add_charac(data_summary,paste("No_genes;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_gene_selected))
    data_summary = add_charac(data_summary,paste("No_introns;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_intron_selected))
    data_summary = add_charac(data_summary,paste("median_length_introns;buscodataset_",busco_group,";quant",sep=""),"",median(by_intron_selected$length))
    data_summary = add_charac(data_summary,paste("No_multi_exonic_genes;buscodataset_",busco_group,";quant",sep=""),"",length(unique(by_intron_selected$gene_id)))
    data_summary = add_charac(data_summary,paste("No_intron_per_gene;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_intron_selected) / length(unique(by_intron_selected$gene_id)))
    data_summary = add_charac(data_summary,paste("prop_intron_annotated;buscodataset_",busco_group,";quant",sep=""),"",sum(by_intron_selected$isAnnotated) / nrow(by_intron_selected))
    data_summary = add_charac(data_summary,paste("median_coverage_exon;buscodataset_",busco_group,";quant",sep=""),"",median(by_gene_selected$exon_coverage,na.rm = T))
    data_summary = add_charac(data_summary,paste("median_gene_fpkm;buscodataset_",busco_group,";quant",sep=""),"",median(by_gene_selected$median_fpkm,na.rm = T))
    data_summary = add_charac(data_summary,paste("prop_overlap;buscodataset_",busco_group,";quant",sep=""),"",sum(by_intron_selected$overlap != 0 ) /nrow(by_intron_selected) )

    for (svr_class in c("all" , "high" , "low")){
      by_intron_selected_svr = by_intron_selected
      if (svr_class == "all"){print("all")} else if ( svr_class == "high" ){
        print("high")
        by_intron_selected_svr = by_intron_selected[  by_intron_selected$splice_variant_rate > 0.05 ,]
        print(nrow(by_intron_selected_svr))
      } else if (svr_class == "low") {
        print("low")
        by_intron_selected_svr = by_intron_selected[  by_intron_selected$splice_variant_rate <= 0.05 ,]
        print(nrow(by_intron_selected_svr))
      }

      data_summary = add_charac(data_summary,paste("average_svr;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",mean(by_intron_selected_svr$splice_variant_rate))
      data_summary = add_charac(data_summary,paste("average_nsvr;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",mean(by_intron_selected_svr$nonsplice_variant_rate))

      Rsvrfpkm = get_rsquared_slope(prop.quantile=0.2,Xaxis=log10(by_intron_selected_svr[ by_intron_selected_svr$median_fpkm != 0 & !is.na(by_intron_selected_svr$median_fpkm),]$median_fpkm),
                                    Yaxis=by_intron_selected_svr[ by_intron_selected_svr$median_fpkm != 0 & !is.na(by_intron_selected_svr$median_fpkm),]$splice_variant_rate)

      Rnsvrfpkm = get_rsquared_slope(prop.quantile=0.2,Xaxis=log10(by_intron_selected_svr[ by_intron_selected_svr$median_fpkm != 0 & !is.na(by_intron_selected_svr$median_fpkm),]$median_fpkm),
                                     Yaxis=by_intron_selected_svr[ by_intron_selected_svr$median_fpkm != 0 & !is.na(by_intron_selected_svr$median_fpkm),]$nonsplice_variant_rate)

      data_summary = add_charac(data_summary,paste("rsquared_svr_fpkm;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",Rsvrfpkm[1]  )
      data_summary = add_charac(data_summary,paste("slope_svr_fpkm;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",Rsvrfpkm[2]  )
      data_summary = add_charac(data_summary,paste("rsquared_nsvr_fpkm;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",Rnsvrfpkm[1]  )
      data_summary = add_charac(data_summary,paste("slope_nsvr_fpkm;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",Rnsvrfpkm[2]  )
    }
  }
}

write.table(data_summary, file = sumv3_path,row.names=F, col.names=T, sep="\t", quote=F)
