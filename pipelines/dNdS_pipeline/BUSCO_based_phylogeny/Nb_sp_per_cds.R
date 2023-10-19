

library(seqinr)
library(ggplot2)


args = (commandArgs(TRUE))
pathCDS = args[1]
pathCDS_85percent = args[2]
species_list = args[3]
path_Folder = args[4]

gene_quantity = data.frame()

list_species = read.table(species_list,header=T)$species

for (file in list.files(pathCDS)){
  print(file)
  r = read.fasta(paste(pathCDS, file, sep = ""))

  gene_quantity=rbind(gene_quantity,data.frame(
    gene=file,
    nb_species=length(r)
  ))
}

total_genes_number = nrow(gene_quantity)

p = ggplot(gene_quantity,aes(x = nb_species)) + geom_histogram(bins=100)+
  xlab("Number of species per genes") + ylab("Number of genes") + theme_bw()+
  geom_vline(xintercept = length(list_species), linetype="dashed",color = "red", size=0.5) +
  geom_vline(xintercept = floor(.85 * length(list_species)), linetype="dashed",color = "blue", size=0.5) +
  ggtitle(paste("Total number of genes :",total_genes_number))

jpeg(paste(path_Folder,"distribution_genes.jpg",sep=""), width = 5000/2, height = 3000/2,res=500/3)
print(p)
dev.off()


gene_quantity = gene_quantity[gene_quantity$nb_species >= floor(.85 * length(list_species) ),]


gene_per_species = data.frame(species=list_species)
rownames(gene_per_species) = gene_per_species$species
gene_per_species$nb_gene = 0

#### Same set gene then Metazoa_v9 for Nematoda_v10
gene_list_met = list.files("/beegfs/XXXXXX/PRANK_CDS_at_least_85percent/")
gene_list_met = gene_list_met[gene_list_met %in% gene_quantity$gene]
####
print(gene_list_met)
#for(gene in gene_quantity$gene){print(gene)
for(gene in gene_list_met){print(gene)

  file.copy(paste(pathCDS,gene,sep="")
            , paste(pathCDS_85percent,gene,sep=""))
  r = read.fasta(paste(pathCDS_85percent,gene,sep=""))
  gene_per_species[sapply(names(r),function(x) strsplit(x,"_busco")[[1]][1]),"nb_gene"] = gene_per_species[sapply(names(r),function(x) strsplit(x,"_busco")[[1]][1]),"nb_gene"] + 1

} #pour charger en local les fichiers nécessaires

p = ggplot(gene_per_species,aes(x = nb_gene)) + geom_histogram(bins=100,col="white")+
  xlab("Number of genes per species") + ylab("Number of species") + theme_bw()+
  geom_vline(xintercept = floor(.8 * length(list.files(pathCDS_85percent))), linetype="dashed",color = "blue", size=0.5) +
  ggtitle(paste("Total number of genes :",total_genes_number))

jpeg(paste(path_Folder,"distribution_species.jpg",sep=""), width = 5000/2, height = 3000/2,res=500/3)
print(p)
dev.off()
