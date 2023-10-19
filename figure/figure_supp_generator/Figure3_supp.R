source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 3 A
ylabel="mean_as_busco"
xlabel="body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3A = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.01,0.1,0.5,1,5,10,100,1000),labels=c(0.01,0.1,0.5,1,5,10,100,1000),limits = c(0.01,1000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(.5,4)) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Body length (cm, log scale)") +  theme(
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
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")

p3A



jpeg(paste(path_figure,"supp_p3A.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p3A)
dev.off()


############## Supplementary Pannel 3 B
ylabel="mean_as_busco"
xlabel="dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3B = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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

p3B



jpeg(paste(path_figure,"supp_p3B.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p3B)
dev.off()



############## Supplementary Pannel 3 C
ylabel = "mean_as_busco_low_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3C = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Low-AS major introns (BUSCO genes)")+
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
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .72, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")
p3C


jpeg(paste(path_figure,"supp_p3C.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p3C)
dev.off()


############## Supplementary Pannel 3 D
ylabel = "mean_as_busco_low_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3D = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
p3D


jpeg(paste(path_figure,"supp_p3D.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p3D)
dev.off()






############## Supplementary Pannel 3 E
ylabel = "mean_as_busco_high_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3E = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("High-AS major introns (BUSCO genes)") +
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
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")

p3E



jpeg(paste(path_figure,"supp_p3E.jpg",sep=""), width = 6000/resolution, height = 5500/resolution,res=700/resolution)
print(p3E)
dev.off()


############## Supplementary Pannel 3 F
ylabel = "mean_as_busco_high_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3F = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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

p3F


jpeg(paste(path_figure,"supp_p3F.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p3F)
dev.off()



############## Supplementary Figure 3

imgA = load.image(paste(path_figure,"supp_p3A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p3B.jpg",sep=""))
imgC = load.image(paste(path_figure,"supp_p3C.jpg",sep=""))
imgD = load.image(paste(path_figure,"supp_p3D.jpg",sep=""))
imgE = load.image(paste(path_figure,"supp_p3E.jpg",sep=""))
imgF = load.image(paste(path_figure,"supp_p3F.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure3_supp.pdf",sep=""), width=6.75, height=2.75*3)
  
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
