source("figure_supp generator/library_path.R")



############## Supplementary Pannel 7 A
ylabel="mean_as_proteincoding"
xlabel="longevity"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (all protein-coding genes)")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,51000)) + theme_bw() +
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + theme(legend.position = "none")

p1



jpeg(paste(path_figure,"p11_proteincoding_as_longevity.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()


############## Supplementary Pannel 7 B
ylabel="mean_as_proteincoding"
xlabel="body_size"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (all protein-coding genes)")+ 
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000)) + theme_bw() +
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + theme(legend.position = "none")

p1



jpeg(paste(path_figure,"p10_proteincoding_as_body.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()


############## Supplementary Pannel 7 C
ylabel="mean_as_proteincoding"
xlabel="dNdS_200k"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (all protein-coding genes)")+ 
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 

p1



jpeg(paste(path_figure,"p9_proteincoding_as_dNdS.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()



############## Supplementary Pannel 7 D
ylabel = "mean_as_proteincoding_low_as"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("Low-AS major introns (all protein-coding genes)")+
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,51000)) +
  scale_y_continuous(breaks=seq(0,10,0.2), labels=paste("",seq(0,10,0.2),"%")) +
  theme(
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")
p6


jpeg(paste(path_figure,"p19_proteincoding_lowas_as_longevity.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()


############## Supplementary Pannel 7 E
ylabel = "mean_as_proteincoding_low_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("Low-AS major introns (all protein-coding genes)")+
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm, log scale)") +
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000))+
  scale_y_continuous(breaks=seq(0,10,0.2), labels=paste("",seq(0,10,0.2),"%")) +
  theme(
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")
p6


jpeg(paste(path_figure,"p21_proteincoding_lowas_as_body.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()



############## Supplementary Pannel 7 F
ylabel = "mean_as_proteincoding_low_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("Low-AS major introns (all protein-coding genes)")+
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("dN/dS") +
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12))+
  scale_y_continuous(breaks=seq(0,10,0.2), labels=paste("",seq(0,10,0.2),"%")) +
  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=21 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 
p6


jpeg(paste(path_figure,"p20_proteincoding_lowas_as_dNdS.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()


############## Supplementary Pannel 7 G
ylabel = "mean_as_proteincoding_high_as"
xlabel = "longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("High-AS major introns (all protein-coding genes)") +
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,51000)) +
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%")) +
  theme(
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")

p



jpeg(paste(path_figure,"p22_proteincoding_highas_as_longevity.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p)
dev.off()



############## Supplementary Pannel 7 H
ylabel = "mean_as_proteincoding_high_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("High-AS major introns (all protein-coding genes)") +
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm, log scale)")  +
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000))+
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%")) +
  theme(
    
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")

p



jpeg(paste(path_figure,"p24_proteincoding_highas_as_body.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p)
dev.off()


############## Supplementary Pannel 7 I

ylabel = "mean_as_proteincoding_high_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("High-AS major introns (all protein-coding genes)") +
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("dN/dS")  +
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12))+
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%")) +
  theme(
    
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
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  )

p


jpeg(paste(path_figure,"p23_proteincoding_highas_as_dNdS.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p)
dev.off()



############## Supplementary Figure 7

imgA = load.image(paste(path_figure,"p11_proteincoding_as_longevity.jpg",sep=""))
imgB = load.image(paste(path_figure,"p10_proteincoding_as_body.jpg",sep=""))
imgC = load.image(paste(path_figure,"p9_proteincoding_as_dNdS.jpg",sep=""))
imgD = load.image(paste(path_figure,"p19_proteincoding_lowas_as_longevity.jpg",sep=""))
imgE = load.image(paste(path_figure,"p21_proteincoding_lowas_as_body.jpg",sep=""))
imgF = load.image(paste(path_figure,"p20_proteincoding_lowas_as_dNdS.jpg",sep=""))
imgG = load.image(paste(path_figure,"p22_proteincoding_highas_as_longevity.jpg",sep=""))
imgH = load.image(paste(path_figure,"p24_proteincoding_highas_as_body.jpg",sep=""))
imgI = load.image(paste(path_figure,"p23_proteincoding_highas_as_dNdS.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure7_supp.pdf",sep=""), width=6.75*3/2, height=2.75*3)
  
  m=matrix(rep(NA,15*3), nrow=3)
  
  m[1,]=c(rep(1,5),rep(2,5),rep(3,5))
  m[2,]=c(rep(4,5),rep(5,5),rep(6,5))
  m[3,]=c(rep(7,5),rep(8,5),rep(9,5))
  
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
  par(mar=c(0, 4.6, 0.5, 2))
  plot(imgG, axes=F)
  mtext("G",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 2.6, 0.5, 4))
  plot(imgH, axes=F)
  mtext("H",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0, 0, 0.5, 0))
  plot(imgI, axes=F)
  mtext("I", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  
  dev.off()
}
