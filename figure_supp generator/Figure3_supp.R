source("figure_supp generator/library_path.R")


############## Supplementary Pannel 3 A
ylabel="prop_major_sv_busco"
xlabel="CoverageBuscoExon"

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel]*100, fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  theme_bw() + ylab("Proportion with N2 > 0")+
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
  ) + scale_x_log10()
p1

jpeg(paste(path_figure,"p18_prop_major_sv_busco_coverage.jpg",sep=""), width = 8200/resolution, height = 5500/resolution,res=700/resolution)
print(p1)
dev.off()


############## Supplementary Figure 3
imgA = load.image(paste(path_figure,"p18_prop_major_sv_busco_coverage.jpg",sep=""))
{
  pdf(file= paste(path_pannel,"Figure3_supp.pdf",sep=""), width=6.75*1/2, height=2.75)
  
  m=matrix(rep(1,15*2), nrow=2)
  
  m
  layout(m)
  
  
  par(mar=c(0, 0.5, 0.5, 0.5))
  plot(imgA, axes=F)
  # mtext("A",at=-100,adj=-2, side=2, line=1, font=2, cex=1.2,las=2)
  
  dev.off()
}
