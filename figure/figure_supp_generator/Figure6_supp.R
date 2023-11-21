source("figure/figure_supp_generator/library_path.R")


############## Supplementary Pannel 6 A
ylabel="mean_as_proteincoding"
xlabel="longevity"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6A = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")

p6A



jpeg(paste(path_figure,"supp_p6A.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6A)
dev.off()


############## Supplementary Pannel 6 B
ylabel="mean_as_proteincoding"
xlabel="body_size"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6B = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  )  + theme(legend.position = "none")+ annotation_logticks(sides="b")

p6B



jpeg(paste(path_figure,"supp_p6B.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6B)
dev.off()


############## Supplementary Pannel 6 C
ylabel="mean_as_proteincoding"
xlabel="dNdS_200k"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6C = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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

p6C



jpeg(paste(path_figure,"supp_p6C.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p6C)
dev.off()



############## Supplementary Pannel 6 D
ylabel = "mean_as_proteincoding_low_as"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6D = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")
p6D


jpeg(paste(path_figure,"supp_p6D.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6D)
dev.off()


############## Supplementary Pannel 6 E
ylabel = "mean_as_proteincoding_low_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6E = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  ) + theme(legend.position = "none")+ annotation_logticks(sides="b")
p6E


jpeg(paste(path_figure,"supp_p6E.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6E)
dev.off()



############## Supplementary Pannel 6 F
ylabel = "mean_as_proteincoding_low_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6F = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
p6F


jpeg(paste(path_figure,"supp_p6F.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p6F)
dev.off()


############## Supplementary Pannel 6 G
ylabel = "mean_as_proteincoding_high_as"
xlabel = "longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6G = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  ) + theme(legend.position = "none")+ annotation_logticks(sides="b")

p6G



jpeg(paste(path_figure,"supp_p6G.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6G)
dev.off()



############## Supplementary Pannel 6 H
ylabel = "mean_as_proteincoding_high_as"
xlabel = "body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6H = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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
  ) + theme(legend.position = "none")+ annotation_logticks(sides="b")

p6H



jpeg(paste(path_figure,"supp_p6H.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6H)
dev.off()


############## Supplementary Pannel 6 I

ylabel = "mean_as_proteincoding_high_as"
xlabel = "dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6I = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
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

p6I


jpeg(paste(path_figure,"supp_p6I.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p6I)
dev.off()



############## Supplementary Figure 6

imgA = load.image(paste(path_figure,"supp_p6A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p6B.jpg",sep=""))
imgC = load.image(paste(path_figure,"supp_p6C.jpg",sep=""))
imgD = load.image(paste(path_figure,"supp_p6D.jpg",sep=""))
imgE = load.image(paste(path_figure,"supp_p6E.jpg",sep=""))
imgF = load.image(paste(path_figure,"supp_p6F.jpg",sep=""))
imgG = load.image(paste(path_figure,"supp_p6G.jpg",sep=""))
imgH = load.image(paste(path_figure,"supp_p6H.jpg",sep=""))
imgI = load.image(paste(path_figure,"supp_p6I.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure6_supp.pdf",sep=""), width=6.75*3/2, height=2.75*3)
  
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
