source("figure generator/library_path.R")

library(stringr)
library(dplyr)
library(tidyr)

pathData="/home/fbenitiere/data//Projet-SplicedVariants/"

dg = read.table("/home/fbenitiere/2023-alternative_splicing/data/Data1_supp.tab",header=T)
rownames(dg) = dg$species
subdt = data.frame()
for (species in rev(dg$species)){print(species)
  
  by_intron = read.delim(file=paste(pathData,species,".tab",sep=""), header=T , sep="\t",comment.char = "#")
  by_intron = by_intron[!is.na(by_intron$phase),]
  
  fpkm_cov = read.delim(paste(pathData,"per_species/",species,"fpkm.tab.gz",sep=""),  sep="\t")
  fpkm_cov = fpkm_cov[fpkm_cov$type == "gene" & grepl("gene_biotype=protein_coding" , fpkm_cov$attributes),]
  
  by_intron = by_intron[by_intron$gene_id %in% fpkm_cov$gene_id,]
  
  ratio = sum(!grepl(",",by_intron$phase)) / nrow(by_intron)
  by_intron = by_intron[ !grepl(",",by_intron$phase) ,]
  print(table( by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id))
  major_introns = by_intron[by_intron$intron_class == "major" & by_intron$into_cds == "True" & by_intron$gene_id %in% fpkm_cov$gene_id,]
  
  fpkm_cov = fpkm_cov[fpkm_cov$busco_metazoa ,]
  CoverageBuscoExon = round(median(tapply(fpkm_cov$exon_coverage,fpkm_cov$gene_id,mean),na.rm=T)) # detection de la couverture médiane des gènes Busco
  
  major_introns_busco = major_introns[major_introns$gene_id %in% fpkm_cov$gene_id,]
  
  table(major_introns$phase) / (nrow(major_introns)/ 3)
  
  major_introns_have_abundant = major_introns[major_introns$have_abundant_sv == "True",]
  major_introns_not_abundant = major_introns[major_introns$have_abundant_sv == "False",]
  # if (length(unique(major_introns_have_abundant$phase)) == 3){
  subdt = rbind(subdt,data.frame(species,
                                 ratio,
                                 major_introns = table(major_introns$phase) / nrow(major_introns),
                                 major_introns_svr = tapply(major_introns$splice_variant_rate,major_introns$phase,mean) * 100,
                                 major_introns_busco = table(major_introns_busco$phase) / nrow(major_introns_busco),
                                 major_introns_busco_svr = tapply(major_introns_busco$splice_variant_rate,major_introns_busco$phase,mean) * 100
                                 # major_introns_have_abundant = table(major_introns_have_abundant$phase) / nrow(major_introns_have_abundant),
                                 # major_introns_have_abundant_svr = tapply(major_introns_have_abundant$splice_variant_rate,major_introns_have_abundant$phase,mean) * 100,
                                 # major_introns_not_abundant = table(major_introns_not_abundant$phase) / nrow(major_introns_not_abundant),
                                 # major_introns_not_abundant_svr = tapply(major_introns_not_abundant$splice_variant_rate,major_introns_not_abundant$phase,mean) * 100
  ))
  # }
}

library(ggplot2)
source("figure_supp generator/library_path.R")
subdt = subdt[subdt$species %in% data_1$species,]

p1 = ggplot(subdt,aes(y=major_introns.Freq * 100,x=(major_introns.Var1) ))  + geom_boxplot(fill=set_color[2]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) +  scale_y_continuous(breaks=seq(0,100,20), labels=paste(seq(0,100,20),"%"),limits = c(0,70)) + 
  ggtitle("(all protein-coding genes)") + theme(legend.position = "none") + 
  ylab("Fraction of introns") + xlab("Phase") 
p1

jpeg(paste(path_figure,"p1_diff_phase.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p1)
dev.off()

p2 = ggplot(subdt,aes(y=major_introns_busco.Freq*100,x=(major_introns_busco.Var1) ))  + geom_boxplot(fill=set_color[1]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(BUSCO genes)") + theme(legend.position = "none") + ylab("Fraction of introns")+ xlab("Phase") + 
  scale_y_continuous(breaks=seq(0,100,20), labels=paste(seq(0,100,20),"%"),limits = c(0,70)) 
p2

jpeg(paste(path_figure,"p2_diff_phase_busco.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p2)
dev.off()


dt = data.frame(
  values = subdt[grepl("0",subdt$major_introns.Var1),]$major_introns_svr - subdt[grepl("1",subdt$major_introns.Var1),]$major_introns_svr,
  group = "Phase 0 - phase 1"
)
dt = rbind(dt,data.frame(
  values = subdt[grepl("0",subdt$major_introns.Var1),]$major_introns_svr - subdt[grepl("2",subdt$major_introns.Var1),]$major_introns_svr,
  group = "Phase 0 - phase 2"
))
dt = rbind(dt,data.frame(
  values = subdt[grepl("2",subdt$major_introns.Var1),]$major_introns_svr - subdt[grepl("1",subdt$major_introns.Var1),]$major_introns_svr,
  group = "Phase 2 - phase 1"
))

p3 = ggplot(dt,aes(y=values ,x=group ))  + geom_boxplot(fill=set_color[2]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(all protein-coding genes)") + theme(legend.position = "none") + ylab("Average AS rate ")+ xlab("") + 
  scale_y_continuous(breaks=seq(-1,1,.2), labels=paste(seq(-1,1,.2),"%")) +
labs(y=expression(paste("Average AS rate ",italic("per")," intron differences")))
p3

jpeg(paste(path_figure,"p3_as_diff_phase.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p3)
dev.off()


dt = data.frame(
  values = subdt[grepl("0",subdt$major_introns_busco.Var1),]$major_introns_busco_svr - subdt[grepl("1",subdt$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "Phase 0 - phase 1"
)
dt = rbind(dt,data.frame(
  values = subdt[grepl("0",subdt$major_introns_busco.Var1),]$major_introns_busco_svr - subdt[grepl("2",subdt$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "Phase 0 - phase 2"
))
dt = rbind(dt,data.frame(
  values = subdt[grepl("2",subdt$major_introns_busco.Var1),]$major_introns_busco_svr - subdt[grepl("1",subdt$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "Phase 2 - phase 1"
))

p4 = ggplot(dt,aes(y=values ,x=group ))  + geom_boxplot(fill=set_color[1]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(BUSCO genes)") + theme(legend.position = "none") + ylab("Average AS rate ")+ xlab("") + 
  scale_y_continuous(breaks=seq(-1,1,.2), labels=paste(seq(-1,1,.2),"%")) +
labs(y=expression(paste("Average AS rate ",italic("per")," intron differences")))
p4

jpeg(paste(path_figure,"p4_as_diff_phase_busco.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p4)
dev.off()

############## Supplementary Figure 1

imgA = load.image(paste(path_figure,"p1_diff_phase.jpg",sep=""))
imgB = load.image(paste(path_figure,"p2_diff_phase_busco.jpg",sep=""))
imgC = load.image(paste(path_figure,"p3_as_diff_phase.jpg",sep=""))
imgD = load.image(paste(path_figure,"p4_as_diff_phase_busco.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure1x_supp.pdf",sep=""), width=5.4, height=4)
  
  m=matrix(c(1,3,2,4), nrow=2)
  
  
  m
  layout(m)
  
  par(mar=c(1, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=0,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=0,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgC, axes = F)
  mtext("C", side=2,at=50,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgD, axes = F)
  mtext("D", side=2,at=60,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  dev.off()
}









