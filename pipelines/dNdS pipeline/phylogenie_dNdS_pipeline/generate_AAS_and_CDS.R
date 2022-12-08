# Ce script forme un fichier d'alignement pour chaque gene busco.

args = (commandArgs(TRUE))
pathPhylo = args[1]
pathData = args[2]
pathHome = args[3]
buscoSample = args[4]
fileName = args[5]

options(stringsAsFactors = F, scipen = 999)
library(readxl)
library(seqinr)
read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- readxl::excel_sheets(filename)
    x <-
    lapply(sheets, function(X)
    readxl::read_excel(filename, sheet = X))
    if (! tibble)
    x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

# Extraire les espèces
# Lecture du fichier excel
mysheets <- read_excel_allsheets(paste(pathData,"Fichiers-data/69_sp.xls", sep = ""))# Lecture du fichier excel
speciesStudied = c()# Extraction des espèces
for (Species in names(mysheets)) {
    if (mysheets[[Species]]$`dN/dS (O/N)`[1] == 'O') {
        speciesStudied = append(speciesStudied, Species)
    }
}
speciesStudied = names(mysheets)
print(speciesStudied)
# Collecte des IDbusco des arthropodes
IDBusco = read.delim(paste(pathData, "/Busco/", buscoSample, "/info/ogs.id.info", sep = ""))# Collecte des IDbusco
IDBusco = IDBusco[duplicated(IDBusco[, 2]) == F, 2]

allSequence = data.frame(matrix(ncol = 6, nrow = 0))
colnames(allSequence) = c('name', 'sequenceAAS', 'sequenceCDS', 'buscoID', 'protein', "species")

information_table=data.frame()


for (sp in speciesStudied) {
    print(sp)
    pathAnalyses = paste(pathData, "/Analyses/Species_", sp, sep = "")
    pathAnnotations = paste(pathData, "/Annotations/Species_", sp, sep = "")
    buscoRef = read.delim(paste(pathAnnotations, "/", fileName, sep = ""))
    # rownames(buscoRef) = buscoRef$Protein
    rownames(buscoRef) = buscoRef$NA..1

    AAsequence = read.fasta(paste(pathAnnotations, "/proteins.faa", sep = ""))
    aas = sapply(AAsequence, function(x) attr(x, 'name')) #recupere le code de la proteine
    aas = data.frame(protein = aas)
    aas$sequence = AAsequence #recupere la sequence de la proteine


    CDSsequence = read.fasta(paste(pathAnnotations, "/CDS.txt", sep = ""))
    cds = sapply(CDSsequence, function(x) attr(x, 'Annot'))
    cds = data.frame(protein = sub(".*protein_id=*(.*?) *].*", "\\1", cds)) #recupere le code de la proteine
    cds$sequence = CDSsequence #recupere la sequence de la proteine
    cds$longueur = lapply(cds$sequence, function(x) length(x)) #recupere la longueur de la cds
    cds = cds[order(unlist(cds$longueur), decreasing = T),]
    cds = cds[duplicated(cds[, 1]) == F,]#elimine les proteines ayant plus d'une cds en gardant la plus grande
    rownames(cds) = cds$protein

    aas$BuscoID=buscoRef[aas$protein,"buscoID"]
    cds$BuscoID=buscoRef[cds$protein,"buscoID"]

    rownames(buscoRef) = buscoRef$buscoID

    information_table=rbind(information_table,data.frame(
    species=sp,
    No_AAS_total = nrow(aas),
    No_CDS_total = nrow(cds),
    No_AAS_Busco = nrow(aas[!is.na(aas$BuscoID),]),
    No_CDS_Busco = nrow(cds[!is.na(cds$BuscoID),])
    ))


    print(information_table[information_table$species==sp,])

    for (buscoID in IDBusco) { #pour chaque gene busco recupere la cds et aas de sa proteine la plus longue
        if (buscoID %in% buscoRef$buscoID) {
            allSequence[paste(sp, "_buscoid:", buscoID, "_longestProt:", buscoRef[buscoID, 3], sep = ""),] =
            list(paste(sp, "_buscoid:", buscoID, "_longestProt:", buscoRef[buscoID, 3], sep = ""), aas[buscoRef[buscoID, 3], 2], cds[buscoRef[buscoID, 3], 2], buscoID, buscoRef[buscoID, 3], sp)}
    }
}

allSequence$sequenceCDS = lapply(allSequence$sequenceCDS, function(x) if (translate(x[(length(x) - 2) : length(x)]) == '*')
{toupper(x[1 : (length(x) - 3)])}else {toupper(x)}) #supprime le dernier codon s'il est STOP=*



for (i in rownames(allSequence)) { #verifie que la CDS correspond a la AAS
    this.translation = translate(allSequence[i, 'sequenceCDS'][[1]])
    if (length(allSequence[i, 'sequenceCDS'][[1]])/3 != length(toupper(allSequence[i, 'sequenceAAS'][[1]])) |
    ! all(this.translation == toupper(allSequence[i, 'sequenceAAS'][[1]])))
    {allSequence = allSequence[- which(rownames(allSequence) == i),]
    }
}

# Write per species the number of genes for which we got the cds and aas, corresponding
information_table$No_AAS_CDS_corresponding=table(allSequence$species)[information_table$species]

write.table(information_table,paste(pathPhylo,"/gene_No_aas_cds",sep=""), row.names=F, col.names=T, sep="\t", quote=F)


buscoAppearance = table(allSequence$buscoID)
length(which(buscoAppearance == length(speciesStudied)))
for (buscoid in IDBusco) {#pour chaque gene busco partagé par toutes les especes ecrit un fichier fasta contenant les CDS et AAS
    if (! is.na(buscoAppearance[buscoid]) & buscoAppearance[buscoid] == length(speciesStudied)) {
        write.fasta(allSequence[allSequence$buscoID == buscoid, 'sequenceCDS'],
        allSequence[allSequence$buscoID == buscoid, 'name'],
        file.out = paste(pathPhylo, "/CDS/", buscoid, ".aln", sep = ""))

        write.fasta(allSequence[allSequence$buscoID == buscoid, 'sequenceAAS'],
        allSequence[allSequence$buscoID == buscoid, 'name'],
        file.out = paste(pathPhylo, "/AAS/", buscoid, ".aln", sep = ""))
    }
}
