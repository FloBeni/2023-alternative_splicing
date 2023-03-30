# Generate per_species tables

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


pathData="/home/XXXXX/data/Projet-SplicedVariants/"
pathData="/beegfs/data/XXXXX/Projet-SplicedVariants/"

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_69species.xls",sep=""))
sp_studied = names(mysheets)


# species = "Drosophila_melanogaster"
for (species in sp_studied ){
  print(species)
  ## INTRON
  busco_tab = read.delim(paste(pathData,"Annotations/",species,"/busco_analysis/busco_to_gene_id_metazoa",sep="" ) )
  
  busco_tab = busco_tab[!(duplicated(busco_tab$busco_id,fromLast = FALSE) | duplicated(busco_tab$busco_id,fromLast = TRUE)) &
                          !(duplicated(busco_tab$gene_id,fromLast = FALSE) | duplicated(busco_tab$gene_id,fromLast = TRUE)) ,]
  
  fpkm_cov = read.delim(paste(pathData,"Analyses/",species,"/by_gene_analysis.tab",sep="") , header=T , sep="\t",comment.char = "#")
  fpkm_cov$busco_metazoa = fpkm_cov$gene_id %in% busco_tab$gene_id
  rownames(fpkm_cov) = fpkm_cov$gene_id
  
  by_intron = read.delim(file=paste(pathData,"Analyses/",species,"/by_intron_cds.tab",sep=""), header=T , sep="\t",comment.char = "#")
  by_intron$id = paste(by_intron$seqname,by_intron$gene_id,by_intron$splice5,by_intron$splice3,by_intron$strand,sep=";")
  by_intron$fpkm = fpkm_cov[by_intron$gene_id,]$weighted_fpkm
  
  IntronLibrary = read.delim(paste(pathData,"Analyses/",species,"/IntronLibrary_inclusive.txt",sep=""))
  IntronLibrary = IntronLibrary %>%
    mutate(Gene = strsplit(as.character(Gene), ",")) %>%
    unnest(Gene) %>%
    filter(Gene != "")
  IntronLibrary$id = paste(IntronLibrary$Chr,IntronLibrary$Gene,IntronLibrary$Splice5,IntronLibrary$Splice3,IntronLibrary$Strand,sep=";")
  IntronLibrary$Annotation = grepl("Annotation",IntronLibrary$Source)
  rownames(IntronLibrary) = IntronLibrary$id
  # by_intron$Annotation = IntronLibrary[by_intron$id,]$Annotation
  by_intron$splicesite = paste(IntronLibrary[by_intron$id,]$SpliceSignal5,IntronLibrary[by_intron$id,]$SpliceSignal3)
  
  annotation_gtf = read.delim(paste(pathData,"Annotations/",species,"/data_source/annotation.gtf",sep=""),header=F)
  annotation_gtf$exon_splice3 = -100
  annotation_gtf[annotation_gtf$V7 == "+",]$exon_splice3 = annotation_gtf[annotation_gtf$V7 == "+",]$V5
  annotation_gtf[annotation_gtf$V7 == "-",]$exon_splice3 = annotation_gtf[annotation_gtf$V7 == "-",]$V4
  annotation_gtf$transcrit = sapply( annotation_gtf$V9 , function(x) str_split(x,";")[[1]][1])
  annotation_gtf$transcrit = str_replace_all(annotation_gtf$transcrit,"transcript_id ","")
  rownames(annotation_gtf) = paste(annotation_gtf$transcrit,annotation_gtf$V3,annotation_gtf$exon_splice3,sep=":")
  
  IntronCoord_original = read.delim(paste(pathData,"Annotations/",species,"/formatted_data/IntronCoords.tab",sep=""))
  IntronCoord_original$Splice5 = NA
  IntronCoord_original$Splice3 = NA
  IntronCoord_original[IntronCoord_original$Strand == 1,]$Splice5 = IntronCoord_original[IntronCoord_original$Strand == 1,]$Start
  IntronCoord_original[IntronCoord_original$Strand == 1,]$Splice3 = IntronCoord_original[IntronCoord_original$Strand == 1,]$End
  IntronCoord_original[IntronCoord_original$Strand == -1,]$Splice3 = IntronCoord_original[IntronCoord_original$Strand == -1,]$Start
  IntronCoord_original[IntronCoord_original$Strand == -1,]$Splice5 = IntronCoord_original[IntronCoord_original$Strand == -1,]$End
  
  IntronCoord = IntronCoord_original %>%
    mutate(Transcripts = strsplit(as.character(Transcripts), ",")) %>%
    unnest(Transcripts) %>%
    filter(Transcripts != "")
  IntronCoord$exon_splice3 = -100
  IntronCoord[IntronCoord$Strand == 1,]$exon_splice3 = IntronCoord[IntronCoord$Strand == 1,]$Start-1
  IntronCoord[IntronCoord$Strand == -1,]$exon_splice3 = IntronCoord[IntronCoord$Strand == -1,]$End+1
  IntronCoord$id = paste(IntronCoord$Transcripts,IntronCoord$exon_splice3,sep=":")
  IntronCoord$phase = annotation_gtf[IntronCoord$id,]$V8
  
  
  df = IntronCoord[grepl(":CDS",IntronCoord$Transcripts),]
  df$id = paste(df$Chr,df$Genes,df$Splice5,df$Splice3,df$Strand,sep=";")
  
  by_intron$phase = tapply(df$phase,df$id,function(x) paste(unique(x),collapse = ","))[by_intron$id]
  
  
  IntronCoord = IntronCoord_original %>%
    mutate(Genes = strsplit(as.character(Genes), ",")) %>%
    unnest(Genes) %>%
    filter(Genes != "")
  IntronCoord$id = paste(IntronCoord$Chr,IntronCoord$Genes,IntronCoord$Splice5,IntronCoord$Splice3,IntronCoord$Strand,sep=";")
  by_intron$Annotation = by_intron$id %in% IntronCoord$id 
  
  minor_introns = read.delim(file=paste(pathData,"Analyses/",species,"/by_minor_intron.tab",sep=""), header=T , sep="\t",comment.char = "#")
  minor_introns$id = paste(minor_introns$seqname,minor_introns$gene_id,minor_introns$splice5,minor_introns$splice3,minor_introns$strand,sep=";")
  rownames(minor_introns) = minor_introns$id
  by_intron$mira = minor_introns[by_intron$id,]$mira
  by_intron$criptic_intron = minor_introns[by_intron$id,]$criptic_intron
  by_intron$distance_from_major = minor_introns[by_intron$id,]$distance_from_major
  by_intron$frame_shift = minor_introns[by_intron$id,]$frame_shift
  
  major_introns = read.delim(paste(pathData,"/Analyses/",species,"/by_intron_major_overlap.tab",sep=""))
  major_introns$id = paste(major_introns$seqname,major_introns$gene_id,major_introns$splice5,major_introns$splice3,major_introns$strand,sep=";")
  rownames(major_introns) = major_introns$id
  by_intron$have_abundant_sv = major_introns[by_intron$id,]$have_abundant_sv
  
  
  write.table(by_intron,paste(pathData,"per_species/",species,"_by_intron_analysis.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(fpkm_cov,paste(pathData,"per_species/",species,"_by_gene_analysis.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}
