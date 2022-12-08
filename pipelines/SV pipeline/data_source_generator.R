# Generate
# Generate RNAseq_list files for each species and the Species.txt file which is used to follow the advancement of the pipeline
# .libPaths()
# Library / Parameters / Arguments
options(stringsAsFactors = F, scipen = 999)
# .libPaths(c( "/beegfs/data/soft/R-3.5.2/lib/R/library" , .libPaths() ))
# .libPaths()
library(readxl)
args <- commandArgs(TRUE)
pathData = args[1]
pathHome = args[2]
Excel_dataset = args[3]

# This function permit to read an excel file
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

# Reading of the excel file
mysheets <- read_excel_allsheets(Excel_dataset)

# Extract SRA and species names
SpeciesTXT = c()
for (Species in names(mysheets)) {
  dir.create(paste(pathData, Species , sep = ''))
    dt = data.frame(SRA_accession_ID = mysheets[[Species]]$`SRA Accession ID`,
               BioProjet = mysheets[[Species]]$`BioProject`)
    dt = dt[!is.na(dt$SRA_accession_ID),]

  write.table(dt

    ,paste(pathData, Species , "/list_Acc.tab", sep = ''),row.names=F, col.names=T, sep="\t", quote=F)
  #if (mysheets[[Species]]$`Analyse (O/N)`[1] == 'O') {
  #  SpeciesTXT = append(SpeciesTXT, mysheets[[Species]]$Species[1])
  #}
}

# Save the Species.txt file setting the step to 0
#SpeciesTXT = gsub(" ", "_", SpeciesTXT)
#SpeciesTXT = paste(SpeciesTXT, '_0', sep = '')
#write(SpeciesTXT, paste(pathData,'/Species.txt',sep=''))
