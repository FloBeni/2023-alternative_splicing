source("figure/figure_supp generator/library_path.R")


############## Supplementary Pannel 9

point_color = c("Brain"="#5c89b2","Cerebellum"="#7bc7eb","Heart"="#8c1b1f","Liver"="#54833b","Kidney"="#c19e40",
                "Testis"="#e16527","Ovary"="#ac3495","Head"="#A6CEE3")

point_shape=c("Gallus gallus"=5,
              "Homo sapiens"=0,
              "Macaca mulatta"=8,
              "Monodelphis domestica"=2,
              "Mus musculus"=1,
              "Oryctolagus cuniculus"=4,
              "Rattus norvegicus"=6,
              "Dendroctonus ponderosae"=9
)

data_4 = read.delim("data/Data4_supp.tab",comment.char = "#")
data_4$species_name = str_replace_all(data_4$species,"_"," ")
data_4$longevity = data_1[data_4$species,]$longevity

data_4 = data_4[data_4$organs %in% names(point_color),]
data_4 = data_4[data_4$species_name %in% names(point_shape),]

life_span_order = tapply(data_4$longevity,data_4$species_name,mean)
life_span_order = life_span_order[order(life_span_order,decreasing = T)]

data_4$species_name = factor( data_4$species_name, levels = names(life_span_order))



p9A = ggplot(font.label = c(50, "plain"),font.legend= c(20, "plain"),font.x= c(20, "plain"),font.y= c(20, "plain"),
            data_4, aes(x=longevity, y=SVR,group=species_name,shape=species_name)) + ggtitle("Major introns (BUSCO genes)")+ 
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
    plot.caption = element_text(hjust = 0.337, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10( breaks = c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,50000)) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 6),order=2, 
                              label.theme = element_text(color="black", size=26, family="serif",vjust = 1.1,margin = margin(t = 5))),
         shape = guide_legend(override.aes = list(stroke=1.1), 
                              label.theme = element_text(color="black", size=26,face="italic", family="serif",vjust = 1.5,margin = margin(t = 5)))
  ) +
  labs(
    caption = substitute(paste("LM model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(lm((data_4$SVR)~log10(data_4$longevity)))))
  ) + annotation_logticks(sides="b")
p9A

jpeg(paste(path_figure , "supp_p9A.jpg",sep=""), width = 8000/resolution, height = 5100/resolution,res=600/resolution)
print(p9A)
dev.off()


############## Supplementary Figure 9
imgA = load.image(paste(path_figure,"supp_p9A.jpg",sep=""))
{
  pdf(file= paste(path_pannel,"Figure9_supp.pdf",sep=""), width=7.75*1/2, height=2.75)
  
  m=matrix(rep(1,15*2), nrow=2)
  
  m
  layout(m)
  
  
  par(mar=c(0, 0.5, 0.5, 0.5))
  plot(imgA, axes=F)
  
  dev.off()
}
