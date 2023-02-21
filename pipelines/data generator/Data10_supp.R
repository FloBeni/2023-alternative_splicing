# Generate Data 8

options(stringsAsFactors = F, scipen = 999)
library(readxl)

pathData="/home/fbenitiere/data/Projet-SplicedVariants/"
# pathData="/beegfs/data/fbenitiere/Projet-SplicedVariants/"

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

mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/metazoa_species.xls",sep=""))

sp_studied = c()
for (species in names(mysheets)){
  if (!is.na(mysheets[[species]]$Group_study[1]) & (mysheets[[species]]$Group_study[1] == "53_sp" | mysheets[[species]]$Group_study[1] == "69_sp")){
    sp_studied = append(sp_studied,species)
  }
}

data_10 = data.frame()
for (species in sp_studied ){
  table_sra = read.delim(paste(pathData,"Annotations/",species,"/SRAruninfo.tab",sep=""))
  table_sra$species = species
  data_10 = rbind(data_10 , table_sra)
}


write.table(data_10,paste("data/Data10_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
