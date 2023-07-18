################################################################################

options(stringsAsFactors = F, scipen = 999)
library('stringr')
library('seqinr')
library(dplyr)
library(tidyr)

args = (commandArgs(TRUE))
codon_usage_path = args[1]
by_gene_path = args[2]
intron_content_path = args[3]
intronic_UC_path = args[4]


by_gene = read.delim(by_gene_path , header=T , sep="\t",comment.char = "#")
rownames(by_gene) = by_gene$gene_id


codon_usage = read.delim(codon_usage_path)
codon_usage$weighted_fpkm = by_gene[codon_usage$gene_id,"weighted_fpkm"]
codon_usage$mean_fpkm = by_gene[codon_usage$gene_id,"mean_fpkm"]
codon_usage$median_fpkm = by_gene[codon_usage$gene_id,"median_fpkm"]


intron_content = read.delim(intron_content_path)


intron_content = intron_content %>%
  mutate(protein_id = strsplit(as.character(protein_id), ";")) %>%
  unnest(protein_id) %>%
  filter(protein_id != "")


intron_content_db = aggregate(intron_content[,8:75], by=list(protein_id=intron_content$protein_id), FUN=sum)

codon_usage = merge(x=codon_usage,y=intron_content_db,by.x="protein_id",by.y = "protein_id",all.x=T,all.y=F,suffixes=c("","_intronic"))

write.table(codon_usage, file=intronic_UC_path, row.names=F, col.names=T, sep="\t", quote=F)



# args = (commandArgs(TRUE))
# codon_usage_path = args[1]
# by_gene_path =  args[2]
# major_path =  args[3]
# intronic_UC_path = args[4]
# 
# by_gene = read.delim(by_gene_path)
# rownames(by_gene) = by_gene$gene_id
# 
# 
# codon_usage = read.delim(codon_usage_path)
# codon_usage$weighted_fpkm = by_gene[codon_usage$gene_id,"weighted_fpkm"]
# codon_usage$mean_fpkm = by_gene[codon_usage$gene_id,"mean_fpkm"]
# codon_usage$median_fpkm = by_gene[codon_usage$gene_id,"median_fpkm"]
# 
# major = read.delim(major_path)
# major = major[major$into_cds == "True",]
# 
# major = aggregate(major[,18:85], by=list(gene_id=major$gene_id), FUN=sum)
# 
# codon_usage = merge(x=codon_usage,y=major,by.x="gene_id",by.y = "gene_id",all.x=T,all.y=F,suffixes=c("","_intronic"))
# 
# write.table(codon_usage, file=intronic_UC_path, row.names=F, col.names=T, sep="\t", quote=F)

