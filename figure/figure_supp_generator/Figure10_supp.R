source("figure/figure_supp_generator/library_path.R")


############## Supplementary Pannel A
data_1$error = data_1$mean_gene_busco_as*100

ylabel="error"
xlabel="longevity"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8A = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~log10(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~log10(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("BUSCO genes")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(7,51000)) + theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("Longevity (days, log scale)") + theme(
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
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"),limits=c(1,25)) + theme(legend.position = "none")+ annotation_logticks(sides="b")

p8A

jpeg(paste(path_figure,"supp_p8A.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p8A)
dev.off()




############## Supplementary Pannel B
data_1$error = data_1$mean_gene_busco_as*100

ylabel="error"
xlabel="body_size"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8B = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~log10(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~log10(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("BUSCO genes")+ 
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=c(0.01,0.1,1,10,100,1000),limits = c(0.01,1000))+  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("Body length (cm, log scale)") + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"),limits=c(1,25)) + theme(legend.position = "none")+ annotation_logticks(sides="b")

p8B

jpeg(paste(path_figure,"supp_p8B.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p8B)
dev.off()




############## Supplementary Pannel C
data_1$error = data_1$mean_gene_busco_as*100
ylabel="error"
xlabel="dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8C = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("BUSCO genes")+ 
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12))+  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("dN/dS") + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"),limits=c(1,25)) 

p8C

jpeg(paste(path_figure,"supp_p8C.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p8C)
dev.off()





############## Supplementary Pannel D
data_1$error = data_1$mean_gene_proteincoding_as*100

ylabel="error"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8D = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~log10(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~log10(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("all protein-coding genes")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(7,51000))+  theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("Longevity (days, log scale)") + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"))+ theme(legend.position = "none")+ annotation_logticks(sides="b")
p8D
jpeg(paste(path_figure,"supp_p8D.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p8D)
dev.off()





############## Supplementary Pannel E
data_1$error = data_1$mean_gene_proteincoding_as*100

ylabel="error"
xlabel="body_size"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8E = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~log10(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~log10(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("all protein-coding genes")+ 
  scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=c(0.01,0.1,1,10,100,1000),limits = c(0.01,1000))+ theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("Body length (cm, log scale)") + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .7, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"))+ theme(legend.position = "none")+ annotation_logticks(sides="b")
p8E
jpeg(paste(path_figure,"supp_p8E.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p8E)
dev.off()



############## Supplementary Pannel F
data_1$error = data_1$mean_gene_proteincoding_as*100

ylabel="error"
xlabel="dNdS_200k"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p8F = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~(xlabel) , shorebird))[1])+geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("all protein-coding genes")+ 
  scale_x_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12))+theme_bw() +
  labs(y=expression(paste("Average AS rate ",italic("per")," gene")))+
  xlab("dN/dS") + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=24 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~(xlabel),shorebird))))
  )  + scale_y_continuous(breaks=seq(0,50,5), labels=paste(seq(0,50,5),"%"))
p8F
jpeg(paste(path_figure,"supp_p8F.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p8F)
dev.off()











############## Supplementary Figure 10

imgA = load.image(paste(path_figure,"supp_p8A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p8B.jpg",sep=""))
imgC = load.image(paste(path_figure,"supp_p8C.jpg",sep=""))
imgD = load.image(paste(path_figure,"supp_p8D.jpg",sep=""))
imgE = load.image(paste(path_figure,"supp_p8E.jpg",sep=""))
imgF = load.image(paste(path_figure,"supp_p8F.jpg",sep=""))


{
  pdf(file= paste(path_pannel,"Figure10_supp.pdf",sep=""), width=6.75*3/2, height=2.75*2)
  
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
