source("figure/figure_main_generator/library_path.R")


############## Pannel 3 A
ylabel="mean_as_busco"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p3A = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ 
  geom_abline(lwd=1,slope = coef(pgls((ylabel)~log10(xlabel) , shorebird))[2], intercept = coef(pgls((ylabel)~log10(xlabel) , shorebird))[1])+
  geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(7,50000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(.5,4)) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + annotation_logticks(sides="b")
p3A



jpeg(paste(path_figure,"p3A.jpg",sep=""), width = 8500/resolution, height = 6000/resolution,res=700/resolution)
print(p3A)
dev.off()




############## Pannel 3 B

point_color = c("Brain"="#5c89b2","Cerebellum"="#7bc7eb","Heart"="#8c1b1f","Liver"="#54833b","Kidney"="#c19e40",
                "Testis"="#e16527","Ovary"="#ac3495","Head"="#A6CEE3")

point_shape=c("Gallus gallus"=5,
              "Homo sapiens"=0,
              "Macaca mulatta"=8,
              "Monodelphis domestica"=2,
              "Mus musculus"=1,
              "Oryctolagus cuniculus"=4,
              "Rattus norvegicus"=6
              # "Dendroctonus ponderosae"=9
)

data_4 = read.delim("data/Data4_supp.tab",comment.char = "#")
data_4$species_name = str_replace_all(data_4$species,"_"," ")
data_4$longevity = data_1[data_4$species,]$longevity

data_4 = data_4[data_4$organs %in% names(point_color),]
data_4 = data_4[data_4$species_name %in% names(point_shape),]

life_span_order = tapply(data_4$longevity,data_4$species_name,mean)
life_span_order = life_span_order[order(life_span_order,decreasing = T)]

data_4$species_name = factor( data_4$species_name, levels = names(life_span_order))



p3B = ggplot(font.label = c(50, "plain"),font.legend= c(20, "plain"),font.x= c(20, "plain"),font.y= c(20, "plain"),
             data_4, aes(x=longevity, y=SVR,group=species_name,shape=species_name)) + ggtitle("Major introns (BUSCO genes)")+ 
  geom_abline(lwd=1,slope = coef(lm((data_4$SVR)~log10(data_4$longevity)))[2], intercept = coef(lm((data_4$SVR)~log10(data_4$longevity)))[1])+
  geom_point(size=5,alpha=1,aes(color=organs),stroke=1.5)+ scale_color_manual("Organs",values=point_color) +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5),labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(0.5,4)) +
  scale_shape_manual("Species",values=point_shape) + labs(fill="Species") +
  xlab("Longevity (days, log scale)") + theme_bw() + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif"),
    plot.caption = element_text(hjust = 0.365, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10( breaks = c(0.05,0.1,0.5,1,5,10,100,1000,10000,50000), limits=c(7,50000)) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 6),order=2, 
                              label.theme = element_text(color="black", size=26, family="serif",vjust = 1.1,margin = margin(t = 5))),
         shape = guide_legend(override.aes = list(stroke=1.1), 
                              label.theme = element_text(color="black", size=26,face="italic", family="serif",vjust = 1.5,margin = margin(t = 5)))
  ) +
  labs(
    caption = substitute(paste("LM model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(lm((data_4$SVR)~log10(data_4$longevity)))))
  ) + annotation_logticks(sides="b")
p3B

jpeg(paste(path_figure , "p3B.jpg",sep=""), width = 8000/resolution, height = 5100/resolution,res=600/resolution)
print(p3B)
dev.off()


#### Figure 3

imgA = load.image(paste(path_figure,"p3A.jpg",sep=""))
imgB = load.image(paste(path_figure,"p3B.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure3.pdf",sep=""), width=6.75, height=2)
  
  m=matrix(rep(1,10*7), nrow=7)
  
  for(i in 1:7){
    m[i,]=c(rep(1,5),rep(2,5))
  }
  m
  layout(m)
  
  par(mar=c(0, 3, 0, 1))
  plot(imgA, axes = F)
  mtext("A",at=80,adj=1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 0, 0, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=70,adj=1, line=1, font=2, cex=1.2,las=2)
  dev.off()
}
