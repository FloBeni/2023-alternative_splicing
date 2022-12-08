# Ce script forme un fichier d'alignement pour chaque gene busco.

args = (commandArgs(TRUE))
pathPhylo = args[1]
pathData = args[2]
pathScripts = args[3]
buscoSample = args[4]
fileName = args[5]

options(stringsAsFactors = F, scipen = 999)
library(readxl)
library(seqinr)

speciesStudied = read.delim(paste(pathScripts,"list_sp",sep=""),header = F)
speciesStudied = speciesStudied$V1

print(speciesStudied)
# Collecte des IDbusco des arthropodes
IDBusco = read.delim(paste(pathData, "Busco/", buscoSample, "/info/ogs.id.info", sep = ""))# Collecte des IDbusco
IDBusco = IDBusco[duplicated(IDBusco[, 2]) == F, 2]

allSequence = data.frame(matrix(ncol = 6, nrow = 0))
colnames(allSequence) = c('name', 'sequenceAAS', 'sequenceCDS', 'busco_id', 'protein', "species")

information_table = data.frame()


for (sp in speciesStudied) {
  print(sp)
  pathAnalyses = paste(pathData, "Analyses/", sp, sep = "")
  pathAnnotations = paste(pathData, "Annotations/", sp, sep = "")
  buscoRef = read.delim(paste(pathAnnotations,"/busco_analysis/", fileName, sep = "/"))
  buscoRef = buscoRef[!(duplicated(buscoRef$busco_id,fromLast = FALSE) | duplicated(buscoRef$busco_id,fromLast = TRUE)) &
                          !(duplicated(buscoRef$gene_id,fromLast = FALSE) | duplicated(buscoRef$gene_id,fromLast = TRUE)) ,]
  rownames(buscoRef) = buscoRef[,"sequence"]

  AAsequence = read.fasta(paste(pathAnnotations, "/data_source/protein.faa", sep = ""))
  aas = sapply(AAsequence, function(x) attr(x, 'name')) #recupere le code de la proteine
  aas = data.frame(protein = aas)
  aas$sequence = AAsequence #recupere la sequence de la proteine


  CDSsequence = read.fasta(paste(pathAnnotations, "/data_source/cds_from_genomic.fna", sep = ""))
  cds = sapply(CDSsequence, function(x) attr(x, 'Annot'))
  cds = data.frame(protein = sub(".*protein_id=*(.*?) *].*", "\\1", cds)) #recupere le code de la proteine
  cds$sequence = CDSsequence #recupere la sequence de la proteine
  cds$longueur = lapply(cds$sequence, function(x) length(x)) #recupere la longueur de la cds
  cds = cds[order(unlist(cds$longueur), decreasing = T),]
  cds = cds[duplicated(cds[, 1]) == F,] #elimine les proteines ayant plus d'une cds en gardant la plus grande
  rownames(cds) = cds$protein

  aas$busco_id=buscoRef[aas$protein,"busco_id"]
  cds$busco_id=buscoRef[cds$protein,"busco_id"]

  rownames(buscoRef) = buscoRef$busco_id

  information_table=rbind(information_table,data.frame(
    species=sp,
    No_AAS_total = nrow(aas),
    No_CDS_total = nrow(cds),
    No_AAS_Busco = nrow(aas[!is.na(aas$busco_id),]),
    No_CDS_Busco = nrow(cds[!is.na(cds$busco_id),])
  ))


  print(information_table[information_table$species==sp,])

  for (busco_id in IDBusco) { #pour chaque gene busco recupere la cds et aas de sa proteine la plus longue
    if (busco_id %in% buscoRef$busco_id) {
      allSequence[paste(sp, "_buscoid:", busco_id, "_longestProt:", buscoRef[busco_id, "sequence"], sep = ""),] =
        list(paste(sp, "_buscoid:", busco_id, "_longestProt:", buscoRef[busco_id, "sequence"], sep = ""), aas[buscoRef[busco_id, "sequence"], 2], cds[buscoRef[busco_id, "sequence"], 2], busco_id, buscoRef[busco_id, "sequence"], sp)
    }
  }
}

allSequence$sequenceCDS = lapply(allSequence$sequenceCDS, function(x) if (translate(x[(length(x) - 2) : length(x)]) == '*')
{toupper(x[1 : (length(x) - 3)])} else {toupper(x[1 : (length(x))])}) #supprime le dernier codon s'il est STOP=*



############################################################################


for (i in rownames(allSequence)) { #verifie que la CDS correspond a la AAS
  this.translation = translate(allSequence[i, 'sequenceCDS'][[1]])
  if (length(allSequence[i, 'sequenceCDS'][[1]])/3 != length(toupper(allSequence[i, 'sequenceAAS'][[1]])) |
    ! all(this.translation == toupper(allSequence[i, 'sequenceAAS'][[1]])))
  {allSequence = allSequence[- which(rownames(allSequence) == i),]
  }
}

# Write per species the number of genes for which we got the cds and aas, corresponding
information_table$No_AAS_CDS_corresponding = table(allSequence$species)[information_table$species]



rownames(information_table) = information_table$species
information_table$nb_gene_85 = 0


buscoAppearance = table(allSequence$busco_id)
length( which( buscoAppearance == length(speciesStudied) ) )
for (busco_id in IDBusco) {
  if (!is.na(buscoAppearance[busco_id]) ) {
    write.fasta(allSequence[allSequence$busco_id == busco_id, 'sequenceCDS'],
                allSequence[allSequence$busco_id == busco_id, 'name'],
                file.out = paste(pathPhylo, "CDS/", busco_id, ".fa", sep = ""))

    write.fasta(allSequence[allSequence$busco_id == busco_id, 'sequenceAAS'],
                allSequence[allSequence$busco_id == busco_id, 'name'],
                file.out = paste(pathPhylo, "AAS/", busco_id, ".fa", sep = ""))

    if (buscoAppearance[busco_id] >= floor(.85 * nrow(information_table) )){
      information_table[sapply(allSequence[allSequence$busco_id == busco_id, 'name'],function(x) strsplit(x,"_busco")[[1]][1]),"nb_gene_85"] =
        information_table[sapply(allSequence[allSequence$busco_id == busco_id, 'name'],function(x) strsplit(x,"_busco")[[1]][1]),"nb_gene_85"] + 1
    }
  }
}

write.table(information_table,paste(pathPhylo,"/gene_No_aas_cds",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
