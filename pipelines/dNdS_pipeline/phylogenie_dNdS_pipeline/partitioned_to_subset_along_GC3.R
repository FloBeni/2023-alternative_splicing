# Créer des concatenats de plus de 200000 sites
options(stringsAsFactors = F, scipen = 999)
library(seqinr)

args = (commandArgs(TRUE))
prank_path = args[1]
output_path = args[2]
data_sequence = args[3]
nb_sites = as.numeric(args[4])

data_sequence = read.table(data_sequence,header=T)
data_sequence = data_sequence[data_sequence$species_filtered & data_sequence$busco_filtered,]
list_species = unique(data_sequence$species)

gene_dt = data.frame(
  busco_gene = names(tapply(data_sequence$GC3_ratio,data_sequence$busco_gene,mean)),
  mean_gc3 = tapply(data_sequence$GC3_ratio,data_sequence$busco_gene,mean),
  length_prank = tapply(data_sequence$length_prank,data_sequence$busco_gene,mean),
  med_AA_length = tapply(data_sequence$AA_length,data_sequence$busco_gene,median)
  )

gene_dt$name_align = paste(gene_dt$busco_gene,".fa-prank.aln",sep="")
rownames(gene_dt) = gene_dt$name_align


k = nb_sites / mean(gene_dt$length_prank)
nb_group = floor(nrow(gene_dt)/k)

group = rep(1:nb_group,500)[1:nrow(gene_dt)]
gene_dt = gene_dt[order(gene_dt$mean_gc3),]
gene_dt$concatenate = group[order(group)]

all.files = paste(gene_dt$busco_gene,".fa-prank.aln",sep="")

dir.create(paste(output_path,sep=""))

for (replicate.id in unique(gene_dt$concatenate)){
  dir.create(paste(output_path, replicate.id,sep=""))
  output.file = paste(output_path, replicate.id,"/concatenatCDS.aln", sep = "")
  info.file = paste(output_path, replicate.id,"/info", sep = "")



  prot = data.frame(protein = list_species)
  rownames(prot) = prot$protein
  list.gene = c()
  size.gene = c()
  GC3.gene = c()
  list.sp.by.gene = c()

#   while (length(all.files) != 0 & length(prot[sapply(rownames(prot)[1],function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]]) <= 200000 ) {
  for (file in gene_dt[gene_dt$concatenate == replicate.id,]$name_align){
    print(length(prot[sapply(rownames(prot)[1],function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]]))
#     file = all.files[1]
    list.gene=append(list.gene,substr(file,1,11))
#     all.files = all.files[-which(all.files == file)]
    protein = read.fasta(paste(prank_path , file, sep = ""), seqtype = "DNA")
    ids.protein = names(protein)

    list.sp.by.gene = append(list.sp.by.gene,length(ids.protein))
    size.gene = append(size.gene, length(protein[[ids.protein[1]]]))
    GC3.gene = append(GC3.gene,gene_dt[file,"mean_gc3"] )

    for (name in prot$protein) {
      if (name %in% sapply(ids.protein,function(x) strsplit(x,"_busco")[[1]][1])) {
        id = ids.protein[which(name == sapply(ids.protein,function(x) strsplit(x,"_busco")[[1]][1]))]
        prot[sapply(id,function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]] = list(append(prot[sapply(id,function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]], protein[[id]]))
        prot[sapply(id,function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]] = list(prot[sapply(id,function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]][prot[sapply(id,function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]] != " "])
      } else {
        id_length = ids.protein[1]
        prot[name, 'sequence'][[1]] = list(append(prot[name, 'sequence'][[1]], rep('-', length(protein[[id_length]]))))
      }
    }
  }

  print(list.gene)
  print(size.gene)
  print(list.sp.by.gene)

  write.fasta(prot$sequence, prot$protein, file.out = output.file)

  info.table=data.frame(gene.id=list.gene,
                        gene.size.nucl=size.gene,
                        average.GC3 = GC3.gene,
                        nb.sp=list.sp.by.gene)

  info.table=rbind(info.table,data.frame(
    gene.id=paste("No genes:",length(list.gene),
                  "  total size:",print(length(prot[sapply(rownames(prot)[1],function(x) strsplit(x,"_busco")[[1]][1]), 'sequence'][[1]])),"pb"),
    gene.size.nucl="",
    average.GC3="",
    nb.sp=""))

  write.table(info.table,info.file, row.names=F, col.names=T, sep="\t", quote=F)
}
