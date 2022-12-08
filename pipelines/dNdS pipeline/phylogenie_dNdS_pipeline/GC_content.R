# Créer des concatenats de plus de 200000 sites
.libPaths()
options(stringsAsFactors = F, scipen = 999)
library(ape)
library(seqinr)
library(stringr)

std <- function(x) sd(x)/sqrt(length(x))

args = (commandArgs(TRUE))
cds_dir = args[1]
prank_dir = args[2]
tree_input = args[3]
output_path = args[4]

data_sequence = data.frame()

list_busco <- str_replace(list.files(paste(cds_dir,sep="")),".fa","")

# GC3 CONTENT
for ( busco_gene in list_busco ){ print(busco_gene)# through fasta file
  fasta_file <- read.fasta(paste(cds_dir,busco_gene,".fa",sep=""))
  prank_file <- read.fasta(paste(prank_dir,"/",busco_gene,".fa-prank.aln",sep=""))
  length_prank = length(prank_file[[1]])

  for (seq_id in names(fasta_file) ){ # through sequences
    species = str_split(seq_id,"_buscoid")[[1]][1]
    sequence = fasta_file[[seq_id]]
    seq_pos3 = sequence[seq(3,length(sequence),3)] # Position 3 de la sequence
    GC3 = sum(seq_pos3 %in% c("c","g")) # Count GC

    data_sequence = rbind(data_sequence,data.frame(
      species ,
      busco_gene ,
      AA_length = length(seq_pos3) ,
      length_prank,
      GC3_count = GC3,
      GC3_ratio = GC3 / (length(seq_pos3))
    ))
  }
}

# Busco filtered
list_busco <- list.files(paste(prank_dir,"_at_least_85percent/",sep=""), full.names = F)
list_busco = str_replace(list_busco,".fa-prank.aln","")
data_sequence$busco_filtered = data_sequence$busco_gene %in% list_busco

# Species filtered
tree = read.tree(paste(tree_input,sep=""))
data_sequence$species_filtered = data_sequence$species %in% tree$tip.label


write.table(data_sequence,paste(output_path,sep=""), row.names=F, col.names=T, sep="\t", quote=F)
