source("figure_supp generator/library_path.R")


############## Supplementary Pannel 4 A
ylabel="mean_as_busco"
xlabel="body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(.5,4)) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm log scale)") +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.72, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")

p1



jpeg(paste(path_figure,"p7_busco_as_body.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()


############## Supplementary Pannel 4 B
ylabel="mean_as_busco"
xlabel="dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(.5,4)) +
   labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("dN/dS")  +  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 

p1



jpeg(paste(path_figure,"p8_busco_as_dNdS.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()



############## Supplementary Pannel 4 C
ylabel = "mean_as_busco_low_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Low-AS major introns (BUSCO genes)")+
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm log scale)") +
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000))+
  scale_y_continuous(breaks=seq(0,10,0.2), labels=paste("",seq(0,10,0.2),"%")) +
  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .72, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none")
p6


jpeg(paste(path_figure,"p13_busco_lowas_as_body.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()


############## Supplementary Pannel 4 D
ylabel = "mean_as_busco_low_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Low-AS major introns (BUSCO genes)")+
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
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 
p6


jpeg(paste(path_figure,"p12_busco_lowas_as_dNdS.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()






############## Supplementary Pannel 4 E
ylabel = "mean_as_busco_high_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("High-AS major introns (BUSCO genes)") +
  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm log scale)")  +
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000))+
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%")) +
  theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
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
  ) + theme(legend.position = "none")

p



jpeg(paste(path_figure,"p15_busco_highas_as_body.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p)
dev.off()


############## Supplementary Pannel 4 F
ylabel = "mean_as_busco_high_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("High-AS major introns (BUSCO genes)") +
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
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  ) 

p


jpeg(paste(path_figure,"p14_busco_highas_as_dNdS.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p)
dev.off()



############## Supplementary Figure 4

imgA = load.image(paste(path_figure,"p7_busco_as_body.jpg",sep=""))
imgB = load.image(paste(path_figure,"p8_busco_as_dNdS.jpg",sep=""))
imgC = load.image(paste(path_figure,"p13_busco_lowas_as_body.jpg",sep=""))
imgD = load.image(paste(path_figure,"p12_busco_lowas_as_dNdS.jpg",sep=""))
imgE = load.image(paste(path_figure,"p15_busco_highas_as_body.jpg",sep=""))
imgF = load.image(paste(path_figure,"p14_busco_highas_as_dNdS.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure4_supp.pdf",sep=""), width=6.75, height=2.75*3)
  
  m=matrix(rep(NA,10*3), nrow=3)
  
  m[1,]=c(rep(1,5),rep(2,5))
  m[2,]=c(rep(3,5),rep(4,5))
  m[3,]=c(rep(5,5),rep(6,5))
  
  m
  layout(m)
  
  
  par(mar=c(0.5, 3, 1.2, 3.5))
  plot(imgA, axes=F)
  mtext("A",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0, 0, .7, 0))
  plot(imgB, axes=F)
  mtext("B", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0.5, 3, 1.2, 3.5))
  plot(imgC, axes=F)
  mtext("C",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)

  par(mar=c(0, 0, .7, 0))
  plot(imgD, axes=F)
  mtext("D", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  par(mar=c(0.5, 3, 1.2, 3.5))
  plot(imgE, axes=F)
  mtext("E",at=20,adj=0, side=2, line=1, font=2, cex=1.7,las=2)

  par(mar=c(0, 0, .7, 0))
  plot(imgF, axes=F)
  mtext("F", side=2,adj=0,at=30, line=1, font=2, cex=1.7,las=2)
  dev.off()
}
