options(stringsAsFactors = F, scipen = 999)

args = (commandArgs(TRUE))
SRAruninfo_path = args[1]
SRAruninfo_output = args[2]
path_RNAseq_table = args[3]


RNAseq_table <- read.delim(path_RNAseq_table)# Lecture du fichier RNA-seq


runinfo = read.delim(SRAruninfo_path)
if (ncol(runinfo) == 1){
    runinfo <- read.csv(SRAruninfo_path)
}


# runinfo = runinfo[runinfo$LibrarySource == "TRANSCRIPTOMIC",]

runinfo = runinfo[runinfo$Run %in% RNAseq_table$SRA_accession_ID,]

write.table(runinfo,SRAruninfo_output, row.names=F, col.names=T, sep="\t", quote=F)

