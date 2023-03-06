source("figure/figure_main generator/library_path.R")

get_rescaled_density<-function(dist_err_raw, fact, plot=TRUE) {
  tmp=dist_err_raw
  tmp$SVRraw=tmp$SVRraw*fact
  tmp=tmp[which(tmp$SVRraw<0.5),]
  tmp$density=tmp$density/(sum(tmp$density))
  
  if(plot==TRUE) {
    main=sprintf("Factor=%.2f", fact)
    plot_dist(tmp, main)
  }
  return(tmp)
}

plot_dist<-function(dist_err_raw, main="", col="black", xlim=c(0,0.5)) {
  test=sum(dist_err_raw$density)
  if(test != 1) warning("Error: sum density != 1")
  #sub=sprintf("Check sum density = %.2f; N=%d", test, nrow(dist_err_raw))
  sub=NULL
  plot(dist_err_raw$SVRraw, dist_err_raw$density, type = "l", xlab="SVR", ylab="Density", sub=sub, main=main, cex.lab=1.5, cex.axis=1.2, xlim=xlim, col=col, lwd=2)
  
  dist_err_raw$SVRweighted=dist_err_raw$SVRraw*dist_err_raw$density
  svr_mean=sum(dist_err_raw$SVRweighted)
  svr_mode=dist_err_raw$SVRraw[which(dist_err_raw$density==max(dist_err_raw$density))]
  frq_highSVR=sum(dist_err_raw$density[which(dist_err_raw$SVRraw>0.05)])
  abline(v=0.05, lty=3, lwd=2)
  
  lsvr_mean=sprintf("Mean SVR = %.3f", svr_mean)
  lsvr_mode=sprintf("Mode SVR = %.3f", svr_mode)
  frq_highSVR=sprintf("Frq. high SVR = %.1f%%", 100*frq_highSVR)
  #legend("topright", legend = c(lsvr_mode, lsvr_mean, frq_highSVR))
  legend("topright", legend = c(lsvr_mean, frq_highSVR))
}

# draw_svr(dist_svr=d_err, nb_introns=nb_introns_tot-nb_introns_fn_SVR, seq_depth=seq_depth)
draw_svr<-function(dist_svr, nb_introns, seq_depth=100) {
  n=nrow(dist_svr)
  NbSVR=NULL
  for(ii in 1:n) {
    ni=round(dist_svr$density[ii]*nb_introns,0)
    if(ni > 0) {
      x=rbinom(n=ni, size = seq_depth, prob = dist_svr$SVRraw[ii])/seq_depth
      NbSVR=c(NbSVR, x)
    }
  }
  ni=length(NbSVR)
  if(ni<nb_introns) {
    for(ii in (ni+1):nb_introns) {
      NbSVR=c(NbSVR, 0)
    }
  }
  return(NbSVR[(1:nb_introns)])
}

myhist<-function(x, main="", breaks=(c(0:101)/100)) {
  n=length(x)
  main=sprintf("%s N=%d introns", main, n)
  hist(x, breaks = breaks, main=main, xlab="SVR", xlim=c(0,0.5))
  svr_median=median(x)
  svr_mean=mean(x)
  frq_highSVR=length(x[which(x>0.05)])/n
  lsvr_mean=sprintf("Mean SVR = %.3f", svr_mean)
  lsvr_mode=sprintf("Median SVR = %.3f", svr_median)
  frq_highSVR=sprintf("Frq. high SVR = %.1f%%", 100*frq_highSVR)
  legend("topright", legend = c(lsvr_mode, lsvr_mean, frq_highSVR))
  abline(v=svr_median, lty=3)
  abline(v=0.05, lty=3, col="red")
}

combine_hist<-function(SVRerr, SVRfn, main="", ymax=0) {
  ner=length(SVRerr)
  nfn=length(SVRfn)
  SVR_tot=c(SVRerr, SVRfn)
  SVR_high=SVR_tot[which(SVR_tot>0.05)]
  SVR_low=SVR_tot[which(SVR_tot<=0.05)]
  n=ner+nfn
  
  svr_mean_er=mean(SVRerr)
  svr_mean_fn=mean(SVRfn)
  svr_mean=mean(SVR_tot)
  svr_mean_high=mean(SVR_high)
  svr_mean_low=mean(SVR_low)
  
  n_highSVR=length(SVR_tot[which(SVR_tot>0.05)])
  n_highSVR_err=length(SVRerr[which(SVRerr>0.05)])
  n_highSVR_fn=length(SVRfn[which(SVRfn>0.05)])
  
  n_SVR_err_detected=length(SVRerr[which(SVRerr>0)])
  n_SVR_fn_detected=length(SVRfn[which(SVRfn>0)])
  frq_SVR_detected = (n_SVR_err_detected + n_SVR_fn_detected)/n
  frq_SVR_err_detected = n_SVR_err_detected /ner
  frq_SVR_fn_detected = n_SVR_fn_detected /nfn
  
  prop_SVRfn_among_detected = (n_SVR_fn_detected)/(n_SVR_err_detected + n_SVR_fn_detected)
  prop_SVRfn_highSVR = (n_highSVR_fn)/(n_highSVR_err + n_highSVR_fn)
  
  
  n_lowSVR=length(SVR_tot[which(SVR_tot<=0.05)])
  n_lowSVR_err=length(SVRerr[which(SVRerr<=0.05)])
  n_lowSVR_fn=length(SVRfn[which(SVRfn<=0.05)])
  
  frq_highSVR=n_highSVR/n
  frq_highSVR_err=n_highSVR_err/ner
  
  f_modulo3_err=0.33
  f_modulo3_fn=1
  f_modulo3_highSVR=(n_highSVR_fn*f_modulo3_fn + n_highSVR_err*f_modulo3_err)/n_highSVR
  f_modulo3_lowSVR=(n_lowSVR_fn*f_modulo3_fn + n_lowSVR_err*f_modulo3_err)/n_lowSVR
  
  
  
  breaks=(c(0:101)/100)
  h_err=hist(SVRerr, breaks = breaks, plot = FALSE)
  h_fn=hist(SVRfn, breaks = breaks, plot = FALSE)
  h_tot=hist(SVR_tot, breaks = breaks, plot = FALSE)
  
  cols=c("orange", "green", "red")
  
  
  ym=max(h_tot$counts/n)
  plot(h_tot$mids, h_tot$counts/n, type = "l", xlab="SVR", ylab="Density", xlim=c(0,0.5),  col=cols[1], lwd=3, main=main, cex.lab=1.5, cex.axis=1.2)
  points(h_fn$mids, h_fn$counts/n, type = "l", col=cols[2], lty=3, lwd=2)
  points(h_err$mids, h_err$counts/n, type = "l", col=cols[3], lty=3, lwd=2)
  abline(v=0.05, lty=3)
  leg1=sprintf("Mean SVR=%.1f%%", svr_mean*100)
  leg2=sprintf("Prop. high SVR=%.1f%%", frq_highSVR*100)
  #legend("topright", title = "All AS", legend = c(leg1, leg2))
  legend(x=0.2,y=ym, title = "All introns", legend = c(leg1, leg2))
  
  
  leg1=sprintf("Mean SVR=%.1f%%", svr_mean_low*100)
  leg2=sprintf("Prop. modulo3=%.0f%%", f_modulo3_lowSVR*100)
  #legend("right", title = "Low SVR", legend = c(leg1, leg2))
  legend(x=0.2,y=ym*2/3, title = "Low SVR", legend = c(leg1, leg2))
  
  leg1=sprintf("Mean SVR=%.1f%%", svr_mean_high*100)
  leg2=sprintf("Prop. modulo3=%.0f%%", f_modulo3_highSVR*100)
  #legend("bottomright", title = "High SVR", legend = c(leg1, leg2), bty="n")
  legend(x=0.2,y=ym/3, title = "High SVR", legend = c(leg1, leg2))
  
  text(x=0.06, y=ym*0.9, labels = "High\nSVR", pos = 4)
  text(x=0.045, y=ym*0.9, labels = "Low\nSVR", pos = 2)
  
  if(ymax>0) {
    main=sprintf("%s (zoom Y-axis)", main)
    plot(h_tot$mids, h_tot$counts/n, type = "l", xlab="SVR", ylab="Density", xlim=c(0,0.5), ylim=c(0,ymax), col=cols[1], lwd=2, main=main, cex.lab=1.5, cex.axis=1.2)
    points(h_fn$mids, h_fn$counts/n, type = "l", col=cols[2], lty=3, lwd=2)
    points(h_err$mids, h_err$counts/n, type = "l", col=cols[3], lty=3, lwd=2)
    abline(v=0.05, lty=3)
    legend("topright", legend = c("All variants", "Functional var.", "Splicing errors"), col=cols, lwd=c(3,2,2), lty=c(1,3,3))
    
    leg1=sprintf("All variants: %.1f%%", svr_mean*100)
    leg2=sprintf("Functional variants: %.1f%%", svr_mean_fn*100)
    leg3=sprintf("Splicing errors: %.1f%%", svr_mean_er*100)
    #legend("bottomright", title = "High SVR", legend = c(leg1, leg2), bty="n")
    legend("right", title = "Mean SVR:", legend = c(leg1, leg2, leg3))
    
  }
  
  
  return(list(svr_mean=svr_mean, svr_mean_low=svr_mean_low, 
              f_modulo3_lowSVR=f_modulo3_lowSVR, svr_mean_high=svr_mean_high, 
              f_modulo3_highSVR=f_modulo3_highSVR, frq_SVR_detected=frq_SVR_detected, 
              prop_SVRfn_among_detected=prop_SVRfn_among_detected, prop_SVRfn_highSVR=prop_SVRfn_highSVR, 
              frq_SVR_err_detected=frq_SVR_err_detected, frq_SVR_fn_detected=frq_SVR_fn_detected, svr_mean_er=svr_mean_er))
  
}






















nb_introns_tot= 400000

# Proportion of introns with functional splice variants
f_fn_svr = 0.05

# Mean (and sd) SVR for functional splice variants
mean_fn_SVR = 0.25
sd_fn_SVR = 0.05

# Parameters of the gamma distribution for erroneous splice variants
gamma_shape = 0.8
gamma_scale = 0.1


# RNAseq sequencing depth (affects the sensitivity of detection of splice variants)
seq_depth = 300






nb_introns_fn_SVR = nb_introns_tot * f_fn_svr

# Number of bins to discretize the density distribution
nbin=1000

# Draw nsim introns with functional variants 
nsim = 1e7
SVRfn = rnorm(n = nsim, mean = mean_fn_SVR, sd = sd_fn_SVR)

# Exclude cases with with SVR<0 (possible with a normal distribution) or > 0.5
SVRfn=SVRfn[which(SVRfn>=0 & SVRfn<0.5)]
# Check that we exclude no more than 10% of the initial dataset
if(length(SVRfn)<(0.9*nsim)) {
  mess=sprintf("WARNING: parameters not consistent with SVR in [0,0.5[")
  warning(mess)
  return(NULL)
}


d_fn_raw=hist(SVRfn, breaks = c(0:nbin)/nbin, plot = F)

dist_fn_raw=data.frame(SVRraw=d_fn_raw$mids)
dist_fn_raw$density=d_fn_raw$density/(sum(d_fn_raw$density))

plot_dist(dist_fn_raw, main="Functional variants", col="green")














# Shape of the distribution of error rates
vmax=20

y = rgamma(n = 200000, shape = gamma_shape, scale = gamma_scale)
if(max(y) > vmax) {
  warning("Error: increase VMAX")
}
d_er_raw=hist(y, breaks = c(0:nbin)*vmax/nbin, plot = F)

dist_err_raw=data.frame(SVRraw=d_er_raw$mids/vmax/2)
dist_err_raw$density=d_er_raw$density/(sum(d_er_raw$density))

plot_dist(dist_err_raw, main="Erroneous variants (very high Ne)", col="red")











NeS=c("High Ne","Medium Ne","Low Ne")



par(mfrow=c(1,2))
# Combine functional and erroneous AS
Ne_vector = c(3,6,10)
Ne_vector = c( 1,3,6)
tab_res=data.frame(factor=Ne_vector)
tab_res$svr_mean=NA
tab_res$svr_mean_low=NA
tab_res$svr_mean_high=NA
tab_res$f_modulo3_lowSVR=NA
tab_res$f_modulo3_highSVR=NA
tab_res$frq_SVR_detected=NA
tab_res$prop_SVRfn_among_detected=NA
tab_res$prop_SVRfn_highSVR=NA
tab_res$frq_SVR_err_detected=NA
tab_res$frq_SVR_fn_detected=NA
tab_res$svr_mean_er=NA


for(ii in 1:nrow(tab_res)) {
  d_err=get_rescaled_density(dist_err_raw, fact=tab_res$factor[ii],plot=FALSE)
  SVRerr=draw_svr(dist_svr=d_err, nb_introns=nb_introns_tot-nb_introns_fn_SVR, seq_depth=seq_depth)
  SVRfn=draw_svr(dist_svr=dist_fn_raw, nb_introns=nb_introns_fn_SVR, seq_depth=seq_depth)
  SVR_tot=c(SVRerr, SVRfn)
  breaks=(c(0:101)/100)
  h_err=hist(SVRerr, breaks = breaks, plot = FALSE)
  h_fn=hist(SVRfn, breaks = breaks, plot = FALSE)
  h_tot=hist(SVR_tot, breaks = breaks, plot = FALSE)
  n = length(SVR_tot)
  
  
  
  res=combine_hist(SVRerr, SVRfn, main=sprintf("Longevity = %.0f", tab_res$factor[ii]), ymax=0.04)
  tab_res$svr_mean[ii]=res$svr_mean
  tab_res$svr_mean_low[ii]=res$svr_mean_low
  tab_res$svr_mean_high[ii]=res$svr_mean_high
  tab_res$f_modulo3_lowSVR[ii]=res$f_modulo3_lowSVR
  tab_res$f_modulo3_highSVR[ii]=res$f_modulo3_highSVR
  tab_res$frq_SVR_detected[ii]=res$frq_SVR_detected
  tab_res$prop_SVRfn_among_detected[ii]=res$prop_SVRfn_among_detected
  tab_res$prop_SVRfn_highSVR[ii]=res$prop_SVRfn_highSVR
  tab_res$frq_SVR_err_detected[ii]=res$frq_SVR_err_detected
  tab_res$frq_SVR_fn_detected[ii]=res$frq_SVR_fn_detected
  tab_res$svr_mean_er[ii]=res$svr_mean_er
  
  
  svr_mean_er=mean(SVRerr)
  svr_mean_fn=mean(SVRfn)
  
  vectorColor=c("Functional variants"="#0a4413","da"="red","non_functional_Gne"="red","Splicing errors"="red","non_functional_Mne"="#ba8e18")
  vectorColorFill=c("All variants"="orange","Functional variants"="#0a4413","da"="red","non_functional_Gne"="red","Splicing errors"="red","non_functional_Mne"="#ba8e18")
  
  dc =  data.frame(density =  h_tot$counts/n,svr  =h_tot$mids, group = rep("All variants"))
  dc =rbind(dc, data.frame(density = h_fn$counts/n,svr =h_fn$mids, group = rep("Functional variants")))
  dc =rbind(dc,data.frame(density = h_err$counts/n,svr =h_err$mids, group = rep("Splicing errors")))
  
  titre = "Mean AS rate:"
  text=paste("All variants:",round(tab_res$svr_mean[ii]*100,1),
             "%\nFunctional variants:",round(svr_mean_fn*100,1),
             "%\nSplicing errors:",round(svr_mean_er*100,1),"%")
  
  
  text_all = paste("All variants:",round(tab_res$svr_mean[ii]*100,1),"%")
  text_fn = paste( "Functional variants:",round(svr_mean_fn*100,1),"%")
  text_err = paste("Splicing errors:",round(svr_mean_er*100,1),"%")
  
  dc$density=dc$density*100
  
  {
    p1 = ggplot(dc,aes(x=svr*100,y=density,col=group)) + geom_line(data = dc[ dc$group=="All variants",],aes(x=svr*100,y=density,col=group),lwd=1,linetype = "dashed")  + theme_bw() +
      ylab("Percentage of introns") + 
      coord_cartesian(ylim=c(0,80),xlim=c(0,50))+
      scale_x_continuous(breaks=seq(0,50,10), labels=paste(seq(0,50,10),"%")) +
      scale_color_manual(values=vectorColorFill) +
      geom_text(label="5%",x=8,y=78,color="black", size=8, family="serif")+
      theme(
        axis.title.x = element_text(color="black", size=31,family="serif"),
        axis.title.y = element_text(color="black", size=31, family="serif"),
        axis.text.y =  element_text(color="black", size=23, family="serif"),
        axis.text.x =  element_text(color=NA, size=NA, family="serif"),
        title =  element_text(color="black", size=31, family="serif"),
        text =  element_text(color="black", size=31, family="serif"),
        legend.title  =  element_text(color=NA, size=26, family="serif"),
        legend.text =  element_text(color="black", size=26, family="serif"),
        # legend.box.background = element_rect(colour = "black")
      ) + xlab("")  +
      scale_y_continuous(breaks=c(0,10,seq(25,80,25)), labels=paste(c(0,10,seq(25,80,25)),"%")) +
      
      scale_fill_manual(values=vectorColor)+geom_area(aes(fill = group,col=NA, group = group),alpha=0.3, position = 'identity')+ guides(group = guide_legend(order = 1)) +
      geom_vline(xintercept=5, linetype="dashed", color = "black", size=1,alpha=0.7)+theme(legend.position = "None") 
    # geom_hline(yintercept=1, linetype="dashed", color = "black", size=1,alpha=0.7)+theme(legend.position = "None")
    
    
    p1
    
    if (ii == 3){
      p1 = p1 +  
        annotate(family="serif",x=28, y=30+14,cex=8.5, label=titre,geom="text")+
        # annotate(family="serif",x=30, y=17+14,cex=7, label=text,geom="text") +
        annotate(family="serif",col="black",x=32-3, y=24+14,cex=7, label=text_all,geom="text") +
        geom_segment(data=dc[1,],x=32-15,xend=32-20,y=24+14,yend=24+14,col="orange",lwd=1.5,linetype = "dashed",alpha=.7) +
        annotate(family="serif",col="black",x=36-3, y=20+14,cex=7, label=text_fn,geom="text") +
        geom_point(data=dc[1,],size=7,pch=21,fill="#0a4413",col="black",x=36-20,y=20+14,alpha=0.7) +
        annotate(family="serif",col="black",x=32.5-1.5, y=16+14,cex=7, label=text_err,geom="text") +
        geom_point(data=dc[1,],size=7,pch=21,fill="red",col="black",x=36-20,y=16+14,alpha=0.7) +
        annotate(family="serif",x=40, y=75,cex=10, label=NeS[ii],geom="text") + labs(x=expression(paste("AS rate ",italic("per")," intron")))
    }else {
      p1 = p1 +  
        annotate(family="serif",x=28, y=30+14,cex=8.5, label=titre,geom="text")+
        # annotate(family="serif",x=30, y=17+14,cex=7, label=text,geom="text") +
        annotate(family="serif",col="black",x=32-3, y=24+14,cex=7, label=text_all,geom="text") +
        geom_segment(data=dc[1,],x=32-15,xend=32-20,y=24+14,yend=24+14,col="orange",lwd=1.5,linetype = "dashed",alpha=.7) +
        annotate(family="serif",col="black",x=36-3, y=20+14,cex=7, label=text_fn,geom="text") +
        geom_point(data=dc[1,],size=7,pch=21,fill="#0a4413",col="black",x=36-20,y=20+14,alpha=0.7) +
        annotate(family="serif",col="black",x=34-3, y=16+14,cex=7, label=text_err,geom="text") +
        geom_point(data=dc[1,],size=7,pch=21,fill="red",col="black",x=36-20,y=16+14,alpha=0.7) +
        annotate(family="serif",x=40, y=75,cex=10, label=NeS[ii],geom="text") + labs(x=expression(paste("AS rate ",italic("per")," intron")))
    }
    
    p1
    
    
    p2  = ggplot(dc,aes(x=svr*100,y=density,col=group)) + geom_line(data = dc[dc$group=="All variants",],aes(x=svr*100,y=density,col=group),lwd=1.5,linetype = "dashed")  + theme_bw() +
      ylab("") + coord_cartesian(ylim=c(0,1),xlim=c(0,50))+  scale_x_continuous(breaks=seq(0,50,10), labels=paste(seq(0,50,10),"%")) +
      scale_color_manual(values=vectorColorFill) +scale_y_continuous(breaks=seq(0,100,0.5), labels=paste(seq(0,100,0.5),"%")) +
      theme(
        axis.title.x = element_text(color="black", size=31,family="serif"),
        axis.title.y = element_text(color="black", size=31, family="serif"),
        axis.text.y =  element_text(color="black", size=23, family="serif"),
        axis.text.x =  element_text(color="black", size=23, family="serif"),
        title =  element_text(color="black", size=31, family="serif"),
        text =  element_text(color="black", size=31, family="serif"),
        legend.title  =  element_text(color=NA, size=26, family="serif"),
        legend.text =  element_text(color="black", size=26, family="serif")
      ) + xlab("") +
      
      
      scale_fill_manual(values=vectorColor)+geom_area(aes(fill = group,col=NA, group = group),alpha=0.3, position = 'identity') + guides(group = guide_legend(order = 1)) +
      geom_vline(xintercept=5, linetype="dashed", color = "black", size=1,alpha=0.7)+theme(legend.position = "None")
    
    p2 = p2 + xlab("AS rate per intron") + labs(x=expression(paste("AS rate ",italic("per")," intron")))+ theme(
      axis.title.x = element_text(color="black", size=31,family="serif"),
      axis.title.y = element_text(color="black", size=31, family="serif"),
      axis.text.y =  element_text(color="black", size=23, family="serif"),
      axis.text.x =  element_text(color="black", size=23, family="serif"),
      title =  element_text(color="black", size=31, family="serif"),
      text =  element_text(color="black", size=31, family="serif"),
      legend.title  =  element_text(color=NA, size=26, family="serif"),
      legend.text =  element_text(color="black", size=26, family="serif")
    )
    
    
    p2
    
    p = ggdraw() + draw_plot(p1, 0.014, .3, 0.955, 0.7) + draw_plot(p2, 0, 0, 0.97, .24) 
    p
    
    
    
    jpeg(paste(path_figure,tab_res$factor[ii],"Pnes.jpg",sep=""), width = 5100/resolution, height = 7000/resolution,res=700/resolution)
    print(p)
    dev.off()
    
    
    ############## Pannel 7 AB
    if (ii > 1){
      p2 = p2 +  theme(
        axis.title.x = element_text(color="black", size=31,family="serif"),
        axis.title.y = element_text(color=NA, size=NA, family="serif"),
        axis.text.y =  element_text(color=NA, size=NA, family="serif"),
        axis.text.x =  element_text(color="black", size=23, family="serif"),
        title =  element_text(color="black", size=31, family="serif"),
        text =  element_text(color="black", size=31, family="serif"),
        legend.title  =  element_text(color=NA, size=26, family="serif"),
        legend.text =  element_text(color="black", size=26, family="serif")
      )
      
      p1 = p1 +  theme(
        axis.title.x = element_text(color="black", size=31,family="serif"),
        axis.title.y = element_text(color=NA, size=NA, family="serif"),
        axis.text.y =  element_text(color=NA, size=NA, family="serif"),
        axis.text.x =  element_text(color=NA, size=NA, family="serif"),
        title =  element_text(color="black", size=31, family="serif"),
        text =  element_text(color="black", size=31, family="serif"),
        legend.title  =  element_text(color=NA, size=26, family="serif"),
        legend.text =  element_text(color="black", size=26, family="serif")
      )+ labs(x=expression(paste("AS rate ",italic("per")," intron")))
      
      
      p = ggdraw() + draw_plot(p1, 0.0075, .3, 0.962, .7) + draw_plot(p2, 0, 0, 0.97, .24) 
      p
      
      
      
      jpeg(paste(path_figure,tab_res$factor[ii],"Pnes.jpg",sep=""), width = 4700/resolution, height = 7000/resolution,res=700/resolution)
      print(p)
      dev.off()
    }
  }
}




############## Pannel 7 AB

p = ggdraw()+ draw_image(paste(path_figure,Ne_vector[3],"Pnes.jpg",sep=""),0.3,0,1,1)  + 
  draw_image(paste(path_figure,Ne_vector[2],"Pnes.jpg",sep=""),0,0,1,1) + 
  draw_image(paste(path_figure,Ne_vector[1],"Pnes.jpg",sep=""),-0.32,0,1,1)
p
jpeg(paste(path_figure,"F1.jpg",sep=""), width = 5300/1, height = 2500/1,res=350)
print(p)
dev.off()




############## Pannel 7 E
data = data.frame(selection = tab_res$factor, svr = tab_res$svr_mean_low*100)

p3  = ggplot(data,aes(x=selection,y=svr)) +geom_point(shape=21,fill="grey",size=9,alpha=1, stroke = 1.5)+
  theme_bw() +  scale_y_continuous(breaks=seq(0,1.5,0.25), labels=paste(seq(0,1.5,0.25),"%"))+ 
  scale_x_continuous(breaks=Ne_vector, labels=NeS)+ coord_cartesian(ylim=c(0,1.5),xlim=c(1,7))+ 
  ylab("Average AS rate per intron")+ labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+ggtitle("Low-AS introns") +
  xlab("")+  theme(
    axis.title.x = element_text(color="black", size=40,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=45, family="serif"),
    axis.text.y =  element_text(color="black", size=35, family="serif"),
    axis.text.x =  element_text(color="black", size=40, family="serif",vjust=-1),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.title  =  element_text(color=NA, size=26, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) 


p3

jpeg(paste(path_figure,"Rare splice variants.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p3)
dev.off()



############## Pannel 7 F
data = data.frame(selection = tab_res$factor, svr = tab_res$svr_mean_high*100)

p4  = ggplot(data,aes(x=selection,y=svr)) +geom_point(shape=21,fill="grey",size=9,alpha=1, stroke = 1.5)+
  theme_bw() +  scale_y_continuous(breaks=seq(10,25,2.5), labels=paste(seq(10,25,2.5),"%"))+ 
  scale_x_continuous(breaks=Ne_vector, labels=NeS)+ coord_cartesian(ylim=c(15,25),xlim=c(1,7))+ 
  ylab("Average AS rate per intron")+ labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+ggtitle("High-AS introns") +
  xlab("")+  theme(
    axis.title.x = element_text(color="black", size=40,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=45, family="serif"),
    axis.text.y =  element_text(color="black", size=40, family="serif"),
    axis.text.x =  element_text(color="black", size=40, family="serif",vjust=-1),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.title  =  element_text(color=NA, size=26, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) 


p4

jpeg(paste(path_figure,"Common splice variants.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p4)
dev.off()



############## Pannel 7 C
data = data.frame(selection = tab_res$factor, svr = tab_res$svr_mean*100)

p6  = ggplot(data,aes(x=selection,y=svr)) +geom_point(shape=21,fill="grey",size=9,alpha=1, stroke = 1.5)+
  theme_bw() +  scale_y_continuous(breaks=seq(1,100,1), labels=paste(seq(1,100,1),"%"))+ 
  scale_x_continuous(breaks=Ne_vector, labels=NeS)+ coord_cartesian(ylim=c(1,3),xlim=c(1,7))+ 
  ylab("Average AS rate per intron")+ labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("")+  theme(
    axis.title.x = element_text(color="black", size=40,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=45, family="serif"),
    axis.text.y =  element_text(color="black", size=40, family="serif"),
    axis.text.x =  element_text(color="black", size=40, family="serif",vjust=-1),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.title  =  element_text(color=NA, size=26, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) 


p6

jpeg(paste(path_figure,"mean_SVR.jpg",sep=""), width = 8000/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()



############## Pannel 7 D
data = data.frame(selection = tab_res$factor, svr = tab_res$prop_SVRfn_highSVR*100)

p5  = ggplot(data,aes(x=selection,y=svr)) +geom_point(shape=21,fill="grey",size=9,alpha=1, stroke = 1.5)+
  theme_bw() +  scale_y_continuous(breaks=seq(10,100,10), labels=paste(seq(10,100,10),"%"))+ 
  scale_x_continuous(breaks=Ne_vector, labels=NeS)+ coord_cartesian(ylim=c(55,100),xlim=c(1,7))+ 
  ylab("Proportion of\nfunctional variants")+ ggtitle("High-AS introns") +
  xlab("")+  theme(
    axis.title.x = element_text(color="black", size=40,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 15, b = 0, l = 0), size=45, family="serif"),
    axis.text.y =  element_text(color="black", size=40, family="serif"),
    axis.text.x =  element_text(color="black", size=40, family="serif",vjust=-1),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.title  =  element_text(color=NA, size=26, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) 


p5

jpeg(paste(path_figure,"mod3_High_SVR.jpg",sep=""), width = 8500/resolution, height = 5500/resolution,res=700/resolution)
print(p5)
dev.off()





#### Figure 7

imgAB = load.image(paste(path_figure,"F1.jpg",sep=""))
imgC = load.image(paste(path_figure,"mean_SVR.jpg",sep=""))
imgD = load.image(paste(path_figure,"mod3_High_SVR.jpg",sep=""))
imgE = load.image(paste(path_figure,"Rare splice variants.jpg",sep=""))
imgF = load.image(paste(path_figure,"Common splice variants.jpg",sep=""))


{
  res=.8
  pdf(file= paste(path_pannel,"Figure7.pdf",sep=""), width=6.5/res, height=7/res)
  
  m=matrix(rep(NA,6*20), nrow=20)
  
  for(i in 1:10){
    m[i,]=c(rep(1,6))
  }
  for(i in 11:17){
    m[i,]=c(rep(2,3),rep(3,3))
  }
  for(i in 16:20){
    m[i,]=c(rep(4,3),rep(5,3))
  }
  m
  layout(m)
  
  
  axesShow = F
  
  
  par(mar=c(0, 0, 0, 0))
  plot(imgAB, axes=axesShow)
  mtext("A",at=0,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)
  mtext("B",at=1850,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgC, axes=axesShow)
  mtext("C",at=50,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)

  par(mar=c(0, 0, 1, 0))
  plot(imgD, axes=axesShow)
  mtext("D",at=50,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)

  plot(imgE, axes=axesShow)
  mtext("E",at=50,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)

  par(mar=c(0, 0, 1, 0))
  plot(imgF, axes=axesShow)
  mtext("F",at=50,adj=-1, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}
