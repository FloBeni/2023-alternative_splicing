source("figure/figure_refPCI_generator/library_path.R")

data_1$mean_gene_proteincoding_as = data_1$mean_gene_proteincoding_as * 100
data_1$mean_gene_proteincoding_as_v2 = data_1$mean_gene_proteincoding_as_v2 * 100
data_1$mean_gene_busco_as = data_1$mean_gene_busco_as * 100
data_1$mean_gene_busco_as_v2 = data_1$mean_gene_busco_as_v2 * 100

############## PCI referee Pannel 1 A
ylabel="mean_gene_busco_as"
xlabel="mean_gene_busco_as_v2"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],
                                                     xlabel=data_1[,xlabel],
                                                     ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1A = ggplot(data_1, aes_string(x=xlabel,y=ylabel,fill="clade"))  + theme_bw() + 
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  labs(x=expression(paste("Average AS rate ",italic("per")," gene version 2")))+
  geom_point(pch=21,alpha=0.7, size=7) + ggtitle("BUSCO genes") + 
  scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%")) +
  scale_x_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%")) +
  scale_fill_manual("Clades",values=vectorColor) +labs(fill="Clades") +
  scale_color_manual("Clades",values=vectorColor) + 
  theme(
    axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
      " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                 lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  ) + theme(legend.position = "none")

p1A


jpeg(paste(path_figure,"supp_p1.2A.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p1A)
dev.off()



############## PCI referee Pannel 1 B

ylabel="mean_gene_proteincoding_as"
xlabel="mean_gene_proteincoding_as_v2"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],
                                                     xlabel=data_1[,xlabel],
                                                     ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1B = ggplot(data_1, aes_string(x=xlabel,y=ylabel,fill="clade"))  + theme_bw() + 
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  labs(x=expression(paste("Average AS rate ",italic("per")," gene version 2")))+
  geom_point(pch=21,alpha=0.7, size=7) + 
  scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%")) +
  scale_x_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%")) +
  scale_fill_manual("Clades",values=vectorColor) +labs(fill="Clades") +
  scale_color_manual("Clades",values=vectorColor) +  ggtitle("all protein-coding genes")+ 
  theme(
    axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
      " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                 lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  ) 

p1B


jpeg(paste(path_figure,"supp_p1.2B.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p1B)
dev.off()



############## Supplementary Figure 7

imgA = load.image(paste(path_figure,"supp_p1.2A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p1.2B.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure1.2_refPCI.pdf",sep=""), width=6.75, height=2.75)
  
  m=matrix(rep(NA,10*1), nrow=1)
  
  m[1,]=c(rep(1,5),rep(2,5))
  
  m
  layout(m)
  
  
  par(mar=c(0.5, 3, 1.2, 3.5))
  plot(imgA, axes=F)
  mtext("A",at=20,adj=0.5, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0, 0, .7, 0))
  plot(imgB, axes=F)
  mtext("B", side=2,adj=0.5,at=30, line=1, font=2, cex=1.7,las=2)
  
  dev.off()
}
