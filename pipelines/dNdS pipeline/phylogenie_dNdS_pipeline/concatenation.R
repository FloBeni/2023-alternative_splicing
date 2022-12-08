
options(stringsAsFactors = F, scipen = 999)

args = (commandArgs(TRUE))
pathPhylo = args[1]
library(seqinr)


all.files=system(paste("ls ",pathPhylo, "PRANK_CDS/", sep=""), intern=T)


passage=0
for(file in all.files){
  if (passage == 0){
    protein=read.fasta(paste(pathPhylo, "PRANK_CDS/", file, sep=""), seqtype="DNA")
    ids.protein=names(protein)

    prot=data.frame(protein=substr(names(protein),1,4))
    prot$sequence=protein
    rownames(prot)=substr(names(protein),1,4)
    passage=1
    for(id in ids.protein){
      prot[substr(id,1,4),'sequence'][[1]] = list(prot[substr(id,1,4),'sequence'][[1]][prot[substr(id,1,4),'sequence'][[1]] != " "])
    }
  } else {
    protein=read.fasta(paste(pathPhylo, "PRANK_CDS/", file, sep=""), seqtype="DNA")
    ids.protein=names(protein)
    for(id in ids.protein){
      prot[substr(id,1,4),'sequence'][[1]] = list(append(prot[substr(id,1,4),'sequence'][[1]],protein[[id]]))
      prot[substr(id,1,4),'sequence'][[1]] = list(prot[substr(id,1,4),'sequence'][[1]][prot[substr(id,1,4),'sequence'][[1]] != " "])
    }
  }
}

write.fasta(prot$sequence, prot$protein, file.out= paste(pathPhylo,"/prank_concatenat/concatenatCDS.aln", sep = ""))

