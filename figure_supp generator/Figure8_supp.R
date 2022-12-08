source("figure_supp generator/library_path.R")



############## Supplementary Pannel 8 A
shorebird <- comparative.data(arbrePhylo, data.frame(species=data_8[,"species"],
                                                     xlabel=data_8[,"body_size"],
                                                     ylabel=data_8[,"ratio"]), species, vcv=TRUE)


p4 = ggplot(data_8, aes(x=body_size,y=ratio,fill=clade))  + theme_bw() + 
  xlab("Body length (cm, log scale)") + ylab("Proportion of frame-preserving SVs")+
  geom_point(pch=21,alpha=0.7, size=7) + 
  scale_y_continuous(breaks=seq(0,100,10), labels=paste(seq(0,100,10),"%")) +
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000))+
  scale_fill_manual("Clades",values=vectorColor) +labs(fill="Clades") +
  scale_color_manual("Clades",values=vectorColor) + ggtitle("Abundant SVs (all protein-coding genes)")+
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + theme(legend.position = "none")

p4


jpeg(paste(path_figure,"p27_ratio_fp_busco_body.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p4)
dev.off()



############## Supplementary Pannel 8 B
shorebird <- comparative.data(arbrePhylo, data.frame(species=data_8[,"species"],
                                                     xlabel=data_8[,"dNdS"],
                                                     ylabel=data_8[,"ratio"]), species, vcv=TRUE)


p4 = ggplot(data_8, aes(x=dNdS,y=ratio,fill=clade))  + theme_bw() + 
  xlab("dN/dS")  +  ylab("Proportion of frame-preserving SVs")+
  geom_point(pch=21,alpha=0.7, size=7) + 
  scale_y_continuous(breaks=seq(0,100,10), labels=paste(seq(0,100,10),"%")) +
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12)) +
  scale_fill_manual("Clades",values=vectorColor) +labs(fill="Clades") +
  scale_color_manual("Clades",values=vectorColor) + ggtitle("Abundant SVs (all protein-coding genes)")+
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 

p4

jpeg(paste(path_figure,"p28_ratio_fp_busco_dNdS.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p4)
dev.off()



############## Supplementary Figure  8

imgA = load.image(paste(path_figure,"p27_ratio_fp_busco_body.jpg",sep=""))
imgB = load.image(paste(path_figure,"p28_ratio_fp_busco_dNdS.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure8_supp.pdf",sep=""), width=6.75, height=2.75)
  
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
