source("figure/figure_supp_generator/library_path.R")


############## Supplementary Pannel 2 A
ylabel="prop_major_sv_busco"
xlabel="CoverageBuscoExon"

p2A = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel]*100, fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  theme_bw() + ylab(expression(paste("Proportion with ", N[a], " > 0")))+
  xlab("Median read coverage on BUSCO genes\n(reads/bp, log scale)" ) +
  scale_y_continuous(breaks=seq(0,100,20), labels=paste(seq(0,100,20),"%"),limits=c(0,100)) + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=22 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("LM model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(lm((data_1[,ylabel])~log10(data_1[,xlabel])))))
  ) + scale_x_log10() + annotation_logticks(sides="b")
p2A

jpeg(paste(path_figure,"supp_p2A.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p2A)
dev.off()


############## Supplementary Figure 2
imgA = load.image(paste(path_figure,"supp_p2A.jpg",sep=""))
{
  pdf(file= paste(path_pannel,"Figure2_supp.pdf",sep=""), width=7.75*1/1, height=2.75)
  
  m=matrix(rep(1,15*2), nrow=2)
  
  m
  layout(m)
  
  
  par(mar=c(0, 0.5, 0.5, 0.5))
  plot(imgA, axes=F)
  
  dev.off()
}
