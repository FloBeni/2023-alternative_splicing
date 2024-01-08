source("figure/figure_supp_generator/library_path.R")


############## Supplementary Pannel Figure 11 A
ylabel="mean_as_busco"
xlabel="longevity"

data_vertebrate = data_1[data_1$clade %in% c("Mammalia","Crocodylia","Aves","Chondrichthyes" ),]

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_vertebrate[,"species"],xlabel=data_vertebrate[,xlabel],ylabel=data_vertebrate[,ylabel]), species, vcv=TRUE)

p2A = ggplot(  data_vertebrate,aes(data_vertebrate[,xlabel],data_vertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  geom_abline(lwd=1,slope = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[2], intercept = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[1])+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(1000,51000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
                              " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                                           lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")

p2A



jpeg(paste(path_figure,"refPCI_p2A.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p2A)
dev.off()


############## Supplementary Pannel Figure 11 B
ylabel="mean_as_busco"
xlabel="body_size"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_vertebrate[,"species"],xlabel=data_vertebrate[,xlabel],ylabel=data_vertebrate[,ylabel]), species, vcv=TRUE)

p2B = ggplot(  data_vertebrate,aes(data_vertebrate[,xlabel],data_vertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(5,1000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm, log scale)") +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
                              " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                                             lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")

p2B


jpeg(paste(path_figure,"refPCI_p2B.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p2B)
dev.off()


############## Supplementary Pannel Figure 11 C
ylabel="mean_as_busco"
xlabel="dNdS_200k"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_vertebrate[,"species"],xlabel=data_vertebrate[,xlabel],ylabel=data_vertebrate[,ylabel]), species, vcv=TRUE)

p2C = ggplot(  data_vertebrate,aes(data_vertebrate[,xlabel],data_vertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  geom_abline(lwd=1,slope = coef(lm((ylabel)~(xlabel),data=shorebird$data))[2], intercept = coef(lm((ylabel)~(xlabel),data=shorebird$data))[1])+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("dN/dS")  +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
                              " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird)),
                                                             lm_eq=lm_eqn(lm((ylabel)~(xlabel),data=shorebird$data))))
  ) 

p2C


jpeg(paste(path_figure,"refPCI_p2C.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p2C)
dev.off()



############## Supplementary Pannel Figure 11 D
ylabel = "mean_as_busco"
xlabel="longevity"

data_invertebrate = data_1[!data_1$clade %in% c("Mammalia","Crocodylia","Aves","Chondrichthyes" ),]
# data_invertebrate = data_1[!data_1$clade %in% c("Mammalia","Crocodylia","Aves","Chondrichthyes" ) & data_1$species != "Megachile_rotundata",]

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_invertebrate[,"species"],xlabel=data_invertebrate[,xlabel],ylabel=data_invertebrate[,ylabel]), species, vcv=TRUE)

p2D = ggplot(  data_invertebrate,aes(data_invertebrate[,xlabel],data_invertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  geom_abline(lwd=1,slope = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[2], intercept = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[1])+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(7,50000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
                              " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                                             lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")
p2D


jpeg(paste(path_figure,"refPCI_p2D.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p2D)
dev.off()


############## Supplementary Pannel Figure 11 E
ylabel = "mean_as_busco"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_invertebrate[,"species"],xlabel=data_invertebrate[,xlabel],ylabel=data_invertebrate[,ylabel]), species, vcv=TRUE)

p2E = ggplot(  data_invertebrate,aes(data_invertebrate[,xlabel],data_invertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  geom_abline(lwd=1,slope = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[2], intercept = coef(lm((ylabel)~log10(xlabel),data=shorebird$data))[1])+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,10)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm, log scale)") +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste(
                              " LM "," R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird)),
                                                             lm_eq=lm_eqn(lm((ylabel)~log10(xlabel),data=shorebird$data))))
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")

p2E


jpeg(paste(path_figure,"refPCI_p2E.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p2E)
dev.off()



############## Supplementary Pannel Figure 11 F
ylabel = "mean_as_busco"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_invertebrate[,"species"],xlabel=data_invertebrate[,xlabel],ylabel=data_invertebrate[,ylabel]), species, vcv=TRUE)

p2F = ggplot(  data_invertebrate,aes(data_invertebrate[,xlabel],data_invertebrate[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  geom_abline(lwd=1,slope = coef(lm((ylabel)~(xlabel),data=shorebird$data))[2], intercept = coef(lm((ylabel)~(xlabel),data=shorebird$data))[1])+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major-isoform introns (BUSCO genes)")+ 
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("dN/dS")  +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"  
  )+
  labs(
    caption = substitute(paste(
                               " LM ","R"^2,lm_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird)),
                                                               lm_eq=lm_eqn(lm((ylabel)~(xlabel),data=shorebird$data))))
  )

p2F
jpeg(paste(path_figure,"refPCI_p2F.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p2F)
dev.off()



############## Supplementary Figure 11

imgA = load.image(paste(path_figure,"refPCI_p2A.jpg",sep=""))
imgB = load.image(paste(path_figure,"refPCI_p2B.jpg",sep=""))
imgC = load.image(paste(path_figure,"refPCI_p2C.jpg",sep=""))
imgD = load.image(paste(path_figure,"refPCI_p2D.jpg",sep=""))
imgE = load.image(paste(path_figure,"refPCI_p2E.jpg",sep=""))
imgF = load.image(paste(path_figure,"refPCI_p2F.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure11_supp.pdf",sep=""), width=6.75*3/2, height=2.75*2)
  
  m=matrix(rep(NA,15*2), nrow=2)
  
  m[1,]=c(rep(1,5),rep(2,5),rep(3,5))
  m[2,]=c(rep(4,5),rep(5,5),rep(6,5))
  m
  layout(m)
  
  par(mar=c(0, 4.6, 0.5, 2))
  plot(imgA, axes=F)
  mtext("A",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 2.6, 0.5, 4))
  plot(imgB, axes=F)
  mtext("B",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 0, 0.5, 0))
  plot(imgC, axes=F)
  mtext("C", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 4.6, 0.5, 2))
  plot(imgD, axes=F)
  mtext("D",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 2.6, 0.5, 4))
  plot(imgE, axes=F)
  mtext("E",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 0, 0.5, 0))
  plot(imgF, axes=F)
  mtext("F", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  
  dev.off()
}
