source("figure/figure_supp_generator/library_path.R")

lm_eqn <- function(m=cor.test(X, Y,method="pearson")){
  cor = m$estimate
  pval_cor = m$p.value
  paste(" =", format(cor, digits = 3) , "   p =",formatC(pval_cor,format = "e", digits = 2))
}


data_6 = read.delim(paste("data/Data6_supp.tab",sep=""),comment.char = "#")
data_6$group = str_replace(data_6$group,"enter","\n")
data_6$group = sapply(data_6$group,function(x){
  split = str_split(x,"\n")[[1]]
  paste(split[1]," genes\n",split[2]," introns",sep="")
} )


############## Supplementary Pannel 5 A
data_sp = data_6[data_6$species=="Homo_sapiens" & data_6$intron == "All"  ,]


p5A = ggplot(data_sp,aes(x=median_gene_expression+10^-10,y=average_as,fill=group))  +
  geom_errorbar(aes(ymin=average_as-average_as_errorBar, ymax=average_as+average_as_errorBar),size=0.5,width = 0) +
  geom_point(pch=21,size=5) +   theme_bw()+ scale_fill_manual("Introns set",values=set_color[c(2,1)])+
  theme( 
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=23, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=22, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .41, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + xlab("Gene expression level (FPKM, log scale)") + 
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10(
    breaks=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    limits=c(0.1,100)
  ) +ggtitle("Major introns") +  labs(
    caption = substitute(paste("Pearson correlation:"," R",pgls_eq), list(pgls_eq=lm_eqn(
      m=cor.test(log10(data_sp[grepl("All protein-coding",data_sp$group),]$median_gene_expression),
                 data_sp[grepl("All protein-coding",data_sp$group),]$average_as_errorBar,method="pearson")))))+
  scale_y_continuous(breaks=seq(0,100,1), labels=paste(seq(0,100,1),"%")) + annotation_logticks(sides="b")
p5A

p5A = ggdraw() + draw_plot(p5A, 0, 0, 1, 1)+
  draw_image(paste(path_require,"human.png",sep=""),.65,.7,.1,.2)

jpeg(paste(path_figure,"supp_p5A.jpg",sep="") , width = 10000/resolution, height = 5000/resolution,res=700/resolution)
print(p5A)
dev.off()



############## Supplementary Pannel 5 B
data_sp = data_6[data_6$species=="Drosophila_melanogaster" & data_6$intron == "All"  ,]

p5B = ggplot(data_sp,aes(x=median_gene_expression+10^-10,y=average_as,fill=group))  +theme_bw() +
  geom_errorbar(aes(ymin=average_as-average_as_errorBar, ymax=average_as+average_as_errorBar),size=0.5,width = 0) +
  geom_point(pch=21,size=5) +   theme_bw()+ scale_fill_manual("Introns set",values=set_color[c(2,1)])+
  theme( 
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=23, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=22, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = .41, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + xlab("Gene expression level (FPKM, log scale)") + 
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10(
    breaks=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    labels=c(0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000,10000,50000),
    limits=c(1,500)
  ) +ggtitle("Major introns") +  labs(
    caption = substitute(paste("Pearson correlation:"," R",pgls_eq), list(pgls_eq=lm_eqn(
      m=cor.test(log10(data_sp[grepl("All protein-coding",data_sp$group),]$median_gene_expression),
                 data_sp[grepl("All protein-coding",data_sp$group),]$average_as,method="pearson")))))+
  scale_y_continuous(breaks=seq(0,100,1), labels=paste(seq(0,100,1),"%")) + annotation_logticks(sides="b")
p5B


p5B = ggdraw() + draw_plot(p5B, 0, 0, 1, 1)+
  draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.7,.7,.05,.2)
p5B


jpeg(paste(path_figure,"supp_p5B.jpg",sep="") , width = 10000/resolution, height = 5000/resolution,res=700/resolution)
print(p5B)
dev.off()



############## Supplementary Figure 5

imgA = load.image(paste(path_figure,"supp_p5A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p5B.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure5_supp.pdf",sep=""), width=4, height=4.5)
  
  m=matrix(c(1,2), nrow=2)
  
  
  m
  layout(m)
  
  par(mar=c(1, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=0,adj=-2, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=0,adj=-2, line=1, font=2, cex=1.1,las=2)
  dev.off()
}


