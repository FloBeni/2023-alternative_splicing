source("figure/figure_main generator/library_path.R")

lm_eqn <- function(m=cor.test(X, Y,method="pearson")){
  cor = m$estimate
  pval_cor = m$p.value
  paste(" =", format(cor, digits = 3) , "   p =",formatC(pval_cor,format = "e", digits = 2))
}



############## Pannel 6 A
xmin = min(data_1$cor_as_fpkm_low_as[which(data_1$pval_cor_as_fpkm_low_as >= 0.05)])

p6A = ggplot(  data_1,aes(x=cor_as_fpkm_low_as) ) + geom_histogram(col="black",fill="#A6CEE3",bins = 25,  boundary = 0) +
  geom_text(label="p < 0.05",x=xmin-.17,y=10,cex=10,family="serif")+
  geom_text(label="p \u2265 0.05",x=-0.15,y=10,cex=10,family="serif")+
  geom_vline(xintercept = xmin,linetype="dashed",size=1)+  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.2, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + ylab("Frequency") + xlab("Pearson coefficient correlation (R)")+ xlim(-1,1)+ ylim(0,25) + ggtitle("Low-AS major introns (all protein-coding genes)")
p6A


jpeg(paste(path_figure,"p6A.jpg",sep="") , width = 8500/resolution, height = 5000/resolution,res=700/resolution)
print(p6A)
dev.off()


############## Pannel 6 B
xmin = min(data_1$cor_as_fpkm_all_as[which(data_1$pval_cor_as_fpkm_all_as>=0.05)])
xmax = max(data_1$cor_as_fpkm_all_as[which(data_1$pval_cor_as_fpkm_all_as>=0.05)])


p6B = ggplot(  data_1,aes(x=cor_as_fpkm_all_as) ) + geom_histogram(col="black",fill="#A6CEE3",bins = 25,  boundary = 0) +
  geom_text(label="p < 0.05",x=xmin-.17,y=10,cex=10,family="serif")+
  geom_text(label="p < 0.05",x=xmax+.17,y=10,cex=10,family="serif")+
  geom_text(label="p \u2265 0.05",x=-0.04,y=10,cex=10,family="serif")+
  geom_vline(xintercept = xmin,linetype="dashed",size=1)+ 
  geom_vline(xintercept = xmax,linetype="dashed",size=1)+ theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.2, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + ylab("Frequency") + xlab("Pearson coefficient correlation (R)")+ xlim(-1,1)+ ylim(0,25)+
  ggtitle("Major introns (all protein-coding genes)")
p6B


jpeg(paste(path_figure,"p6B.jpg",sep="") , width = 8500/resolution, height = 5000/resolution,res=700/resolution)
print(p6B)
dev.off()



data_6 = read.delim(paste("data/Data6_supp.tab",sep=""),comment.char = "#")
data_6$group = str_replace(data_6$group,"enter","\n")


############## Pannel 6 C
data_sp = data_6[data_6$species=="Homo_sapiens" & data_6$intron == "Rare_SV",]

p6C = ggplot(data_sp,aes(x=median_gene_expression+10^-10,y=average_as,fill=group))  +theme_bw() +
  geom_errorbar(aes(ymin=average_as-average_as_errorBar, ymax=average_as+average_as_errorBar),size=0.5,width = 0) +
  geom_point(pch=21,size=5) +
  theme_bw()+ scale_fill_manual("Gene set",values=set_color[c(2,2,4,8)])+
  theme( 
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=22, family="serif"),
    plot.caption = element_text(hjust = .65, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + xlab("Gene expression level (FPKM, log scale)") + labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10(
    breaks=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    limits=c(0.1,100)
  ) + ggtitle("Low-AS major introns (all protein-coding genes)") +  labs(
    caption = substitute(paste("Pearson correlation:"," R",pgls_eq), list(pgls_eq=lm_eqn(
      m=cor.test(log10(data_sp$median_gene_expression), data_sp$average_as,method="pearson")))))+
  scale_y_continuous(breaks=seq(0,100,.25), labels=paste(seq(0,100,.25),"%")) + theme(legend.position = "none")
print(p6C)
p6C = ggdraw() + draw_plot(p6C, 0, 0, 1, 1)+
  draw_image(paste(path_require,"human.png",sep=""),.85,.7,.1,.2)

jpeg(paste(path_figure,"p6C.jpg",sep="") , width = 8500/resolution, height = 5000/resolution,res=700/resolution)
print(p6C)
dev.off()



############## Pannel 6 D
data_sp = data_6[data_6$species=="Drosophila_melanogaster" & data_6$intron == "Rare_SV",]

p6D = ggplot(data_sp,aes(x=median_gene_expression+10^-10,y=average_as,fill=group))  +theme_bw() +
  geom_errorbar(aes(ymin=average_as-average_as_errorBar, ymax=average_as+average_as_errorBar),size=0.5,width = 0) +
  geom_point(pch=21,size=5) +   theme_bw()+ scale_fill_manual("Gene set",values=set_color[c(2)])+
  theme( 
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=22, family="serif"),
    plot.caption = element_text(hjust = .63, face= "italic", size=23),
    plot.caption.position =  "plot"
    
  ) + ylab("") + xlab("Gene expression level (FPKM, log scale)")  + 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
                labels=c(0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000)
                ,limits=c(1,500)
  ) + ggtitle("Low-AS major introns (all protein-coding genes)")+ ylab("Average AS rate per intron") +labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+  labs(
    caption = substitute(paste("Pearson correlation:"," R",pgls_eq), list(pgls_eq=lm_eqn(
      m=cor.test(log10(data_sp$median_gene_expression), data_sp$average_as,method="pearson")))))+
  scale_y_continuous(breaks=seq(0,100,.1), labels=paste(seq(0,100,0.1),"%")) + theme(legend.position = "none")
print(p6D)

p6D = ggdraw() + draw_plot(p6D, 0, 0, 1, 1)+
  draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.85,.7,.07,.2)
p6D


jpeg(paste(path_figure,"p6D.jpg",sep="") , width = 8500/resolution, height = 5000/resolution,res=700/resolution)
print(p6D)
dev.off()




#### Figure 6

imgA = load.image(paste(path_figure,"p6A.jpg",sep=""))
imgB = load.image(paste(path_figure,"p6B.jpg",sep=""))
imgC = load.image(paste(path_figure,"p6C.jpg",sep="") )
imgD = load.image(paste(path_figure,"p6D.jpg",sep="") )


{
  pdf(file=paste(path_pannel,"Figure6.pdf",sep=""), width=6.75, height=4)
  
  m=matrix(rep(NA,2*2), nrow=2)
  
  for(i in 1:2){
    m[i,]=c(rep(1,1),rep(2,1))
  }
  
  for(i in 2){
    m[i,]=c(rep(3,1),rep(4,1))
  }
  m
  layout(m)
  
  par(mar=c(0, 2, 0, 1))
  plot(imgA, axes=F)
  mtext("A",at=-20,adj=0.6, side=2, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 2, 0, 1))
  plot(imgB, axes=F)
  mtext("B",at=-20,adj=0.6, side=2, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgC, axes=F)
  mtext("C",at=-20,adj=-1, side=2, line=1, font=2, cex=1,las=2)
  
  par(mar=c(0, 1, 0, 1))
  plot(imgD, axes=F)
  mtext("D",at=-20,adj=-1, side=2, line=1, font=2, cex=1,las=2)
  
  dev.off()
}
