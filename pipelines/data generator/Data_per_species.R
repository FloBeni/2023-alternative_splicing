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
pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_species.xls",sep=""))


sp_studied = c()
for (species in names(mysheets)){
  if (!is.na(mysheets[[species]]$Group_study[1]) & (mysheets[[species]]$Group_study[1] == "53_sp" | mysheets[[species]]$Group_study[1] == "69_sp")){
    sp_studied = append(sp_studied,species)
  }
}


# species = "Drosophila_melanogaster"
for (species in rev(sp_studied) ){
  print(species)
  ## INTRON
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  
  busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                          !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
  
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
  
  
  write.table(by_intron,paste(pathData,species,".tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(fpkm_cov,paste(pathData,species,"fpkm.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}
