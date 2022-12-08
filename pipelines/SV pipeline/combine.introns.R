################################################################################

options(stringsAsFactors = F, scipen = 999)
# .libPaths(c( "/beegfs/data/soft/R-3.5.2/lib/R/library" , .libPaths() ))
# .libPaths()
library('stringr')
library('seqinr')
args = (commandArgs(TRUE))
print(args)
species = args[1]
pathData = args[2]
pathIntronCoords = args[3]
path_genome = args[4]
output_intronlibrary = args[5]
pathRNASEQ = args[6]
print(species)


pathAnalyses=paste(pathData, "Analyses-RNAseq/", species, sep="")

################################################################################

## first read annotated junctions

annot=read.table(pathIntronCoords, h=T, stringsAsFactors=F, sep="\t", quote="\"")

if(length(grep(",", annot$Chr))!=0){
    stop("Weird: there are commas in chromosome name, this will cause problems")
}

annot$Strand[which(annot$Strand=="+")]="1"
annot$Strand[which(annot$Strand=="-")]="-1"

## remove NA introns if any

annot=annot[which(!is.na(annot$Genes)),]


annot$ID=paste(annot$Chr, annot$Start, annot$End, annot$Strand, sep=",")

all.introns=c(annot$ID)
all.sources=rep("Annotation", length(annot$ID))
all.genes=annot$Genes

###############################################################################

## then read RNAseq splice junctions
samples=read.delim(pathRNASEQ)

for(sample in samples[,1]){
    if(file.exists(paste(pathAnalyses, "/", sample,"/junctions_gene_ids_inclusive.txt", sep=""))){
        print(sample)

        junctions.all=read.table(paste(pathAnalyses, "/", sample,"/junctions_gene_ids_inclusive.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
        junctions.all$ID=paste(junctions.all$Chr, junctions.all$Start, junctions.all$End, junctions.all$Strand, sep=",")

        if(length(grep(",", junctions.all$Chromosome)) != 0){
            stop("Weird: there are commas in chromosome name, this will cause problems")
        }

        all.introns=c(all.introns, junctions.all$ID)
        ## all.sources=c(all.sources, rep("RNASeq", length(junctions.all$ID)))

        all.sources=c(all.sources, rep(paste("RNASeq", sample, sep="_"), length(junctions.all$ID)))
        all.genes=c(all.genes, junctions.all$Gene)
    }
}

###############################################################################

unique.introns=levels(as.factor(all.introns))
unique.sources=tapply(all.sources, as.factor(all.introns), function(x) paste(unique(x), collapse=","))
unique.genes=tapply(all.genes, as.factor(all.introns), function(x) {y=unique(unlist(lapply(x, function(z) unlist(strsplit(x, split=","))))); return(paste(y, collapse=","))})
print(unique.introns)
unique.info=lapply(unique.introns, function(x) unlist(strsplit(x,split=",")))

chr=unlist(lapply(unique.info, function(x) x[1]))
start=unlist(lapply(unique.info, function(x) x[2]))
end=unlist(lapply(unique.info, function(x) x[3]))
strand=unlist(lapply(unique.info, function(x) x[4]))

###############################################################################

results=data.frame("Chr"=chr, "Start"=start, "End"=end, "Strand"=strand, "Source"=unique.sources, "Gene"=unique.genes, stringsAsFactors=F)

results$Splice5 = NA
results$Splice3 = NA

results[results$Strand == 1,]$Splice5 = results[results$Strand == 1,]$Start
results[results$Strand == -1,]$Splice5 = results[results$Strand == -1,]$End
results[results$Strand == 1,]$Splice3 = results[results$Strand == 1,]$End
results[results$Strand == -1,]$Splice3 = results[results$Strand == -1,]$Start

genome = read.fasta(path_genome)

SpliceSignal = apply(results,1,function(x){
      if (!is.na(x["Strand"])){
        if (as.numeric(x["Strand"]) == 1){
          return( toupper(c(genome[[x["Chr"]]][as.numeric(x["Splice5"]) : (as.numeric(x["Splice5"])+1) ] ,
                            genome[[x["Chr"]]][ (as.numeric(x["Splice3"])-1) : as.numeric(x["Splice3"])    ] )))
        } else if (as.numeric(x["Strand"]) == -1){
          return( chartr("ATGC","TACG",toupper(c(genome[[x["Chr"]]][ as.numeric(x["Splice5"]) : (as.numeric(x["Splice5"])-1) ] ,
                                                 genome[[x["Chr"]]][  (as.numeric(x["Splice3"])+1) : as.numeric(x["Splice3"])  ] ))))
        }} else{
          return( c(NA,NA,NA,NA))
        } })

results$SpliceSignal5 = paste(SpliceSignal[1,],SpliceSignal[2,],sep="")
results$SpliceSignal3 = paste(SpliceSignal[3,],SpliceSignal[4,],sep="")

composition = apply(results,1,function(x){
  if (!is.na(x["Strand"])){
    if (as.numeric(x["Strand"]) == 1){
      freq = table( toupper(genome[[x["Chr"]]][(as.numeric(x["Splice5"])+2) : (as.numeric(x["Splice3"])-2)]))
      return( freq[c("A","T","C","G")])
    } else if (as.numeric(x["Strand"]) == -1){
      freq = table( chartr("ATGC","TACG",toupper(c(genome[[x["Chr"]]][ (as.numeric(x["Splice5"])-2) : (as.numeric(x["Splice3"])+2)]))))
      return( freq[c("A","T","C","G")])
        }} else{
      return( c(NA,NA,NA,NA) )
    }
    }
    )

results$nb_A = composition[1,]
results$nb_T = composition[2,]
results$nb_C = composition[3,]
results$nb_G = composition[4,]


results$length = abs(as.numeric(results$Splice3) - as.numeric(results$Splice5))
results$nb_nucl = results$nb_G + results$nb_T + results$nb_A + results$nb_C

print( table( results$nb_nucl == results$length - 3 ) )

write.table(results, file=output_intronlibrary, row.names=F, col.names=T, sep="\t", quote=F)

