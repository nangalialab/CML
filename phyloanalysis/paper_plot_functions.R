
## Plots a 
do_tree_paper_plot=function(PD,b.add.lines=FALSE, cex.annot=0,b.extend.colours=TRUE,ff=0.15,prioritise.group=NULL,
  b.add.bs=FALSE,min.clone.count=-1,vspace.reserve=NA,force.count=NA,expand.top.tree="none",ymax=-1,...){
  if(expand.top.tree=="none"){
    tree=do_tree_paper_plot2(PD,b.add.lines = b.add.lines,cex.annot=cex.annot,b.extend.colours = b.extend.colours,ff=ff,prioritise.group=prioritise.group,b.add.bs=FALSE,
                             min.clone.count = min.clone.count,vspace.reserve = vspace.reserve,force.count = force.count,ymax=ymax,...)
  }else{
    if(expand.top.tree=="ultra"){
      # alpha=0.1
      a=0.75
      H=ymax#max(nodeHeights(PD$pdx$tree_ml))
      MH=10*floor(H/10)
      sf=12#alpha*(H-a)/(a*(1-alpha))
      PD$pdx=magnify_top_of_tree_pdx(PD$pdx,a=a,scale.factor = sf)
      ymax=ymax+a*(sf-1)
      tree=do_tree_paper_plot2(PD,b.add.lines = b.add.lines,cex.annot=cex.annot,b.extend.colours = b.extend.colours,ff=ff,prioritise.group=prioritise.group,b.add.bs=FALSE,
                               min.clone.count = min.clone.count,vspace.reserve = vspace.reserve,force.count = force.count,b.add.scale=FALSE,ymax=ymax,...)
      add_scale(tree,a =a,scale.factor = sf,
                breaks1 = c(0),breaks2=seq(0,MH,10)+a,
                labels1 = c("Zygote"),labels2=c("Birth",seq(10,MH,10))
      )
      
    }else{
      H=ymax#max(nodeHeights(PD$pdx$tree_ml))
      PD$pdx=set_color_by_age(PD$pdx)
      MH=100*floor(H/100)
      a=50
      sf=3
      PD$pdx=magnify_top_of_tree_pdx(PD$pdx,a=a,scale.factor = sf)
      ymax=ymax+a*(sf-1)
      tree=do_tree_paper_plot2(PD,b.add.lines = b.add.lines,cex.annot=cex.annot,b.extend.colours = b.extend.colours,ff=ff,prioritise.group=prioritise.group,b.add.bs=FALSE,
                               min.clone.count = min.clone.count,vspace.reserve = vspace.reserve,force.count = force.count,b.add.scale=FALSE,ymax=ymax,...)
      add_scale(tree,a = a,scale.factor = sf,breaks1 = c(0,50),
                breaks2=seq(200,MH,200),labels1 = c("Zygote","Birth"))
      
    }
  }
  tree
}

do_tree_paper_plot2=function(PD,b.add.lines=FALSE, cex.annot=0,b.extend.colours=TRUE,ff=0.15,
                             prioritise.group=NULL,
                             b.add.bs=FALSE,min.clone.count=-1,vspace.reserve=NA,force.count=NA,...){
  lm=ff/(1-ff)
  #browser()
  ds=get_driver_scheme()
  dff=get_all_tree_drivers(PD$pdx,genes=GENES,cv = CV)
  ## Some label transformations.
  dff=dff %>% mutate(label=ifelse(label=="t(9;22)(q34;q11)","BCR::ABL1",label))
  dff=dff %>% mutate(label=ifelse(label=="Tri8_AlleleA","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(label=="Tri8_AlleleB","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(label=="Tri8_AlleleC","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(label=="Tri8_AlleleD","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(label=="Trisomy8A","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(label=="Trisomy8B","Trisomy 8",label))
  dff=dff %>% mutate(label=ifelse(grepl("loY",label),"Loss of Y",label))
  
  dff$group=gsub(":.+","",dff$label)
  dff$group=gsub("_[A-Z]$","",dff$group)
  dff=dff %>% left_join(ds,by=c("group"))
  dff=dff[order(dff$group,dff$cv),]
  print(dff)
  if(b.add.bs){
    PD=add_bs(PD,get_nexus(PD))
  }
  tree=PD$pdx$tree_ml
  ID=PD$patient
  if(dim(dff)[1]==0){
    if(is.na(vspace.reserve)){
      vspace.reserve=ifelse(dim(dff)[1]>=10,0.8,0.3)
    }
    if(is.na(force.count)){
      force.count=ifelse(dim(dff)[1]>=10,-1,8)
    }
    pdx=set_color_by_age(PD$pdx)
    tree=plot_basic_tree(pdx,ID,
                         style="vertical",
                         genes=GENES,
                         cex.label=0,
                         cex.annot.label = 1,b.show.just.group=TRUE,
                         cv = CV ,left.margin.prop = lm,cex.terminal.dots = 0.5,seed=123576,legpos = NULL,vspace.reserve=vspace.reserve,
                         b.extend.colours = b.extend.colours,...
    );
    cord=tree$coords[match(1:length(tree$tip.label),tree$edge[,2]),]
    points(cord$b1,tree$top-cord$a1,col=tree$tip.color,cex=1,pch=19)
    legend("topleft",legend = sprintf("%3.1f",pdx$age.df$age),col=pdx$age.df$color,pch=19,title = "Age at Sample",bty = "n")
  }else{
    dff$nsample=sapply(dff$node,function(node) length(get_samples_in_clade(node,PD$pdx$tree_ml)))
    dff=dff[order(-dff$nsample),]
    nodes=dff$node
    labels=dff$label
    hm=matrix("white",ncol=length(tree$tip.label),nrow=length(nodes))
    colnames(hm)=tree$tip.label
    rownames(hm)=labels;
    for(i in 1:length(nodes)){
      tips=get_samples_in_clade(nodes[i],tree);hm[i,]=adjustcolor(dff$colour[i],alpha.f = 0.15)
      hm[i,match(tips,tree$tip.label)]=dff$colour[i]
    }
    print(hm)
    if(min.clone.count>0){
      cols=c("black","black")
      lineages=get_lineages(tree,min.clone.count)
      hsc=rep("white",length(tree$tip.label))
      filler=hsc
      if(dim(lineages)[1]>0){
        inf=lapply(lineages$node,function(node) get_samples_in_clade(node,tree))
        k=1
        for(x in inf){
          hsc[match(x,tree$tip.label)]=cols[(k %% length(cols))+1]
          k=k+1
        }
      }
      hsc=rbind(hsc,filler)
      colnames(hsc)=colnames(hm)
      hm=rbind(hsc,hm)
      rownames(hm)[1:2]=c("Clonal Fraction","")
      txtcols=c("black","black",dff$colour)
    }else{
      txtcols=dff$colour
    } 
    pdx=set_color_by_age(PD$pdx)
    if(is.na(vspace.reserve)){
      vspace.reserve=ifelse(dim(dff)[1]>=10,0.8,0.3)
    }
    if(is.na(force.count)){
      force.count=ifelse(dim(dff)[1]>=10,-1,8)
    }
    #a=50
    #sf=50
    tree=plot_basic_tree(pdx,ID,#plot_basic_tree(magnify_top_of_tree_pdx(pdx,a=50,scale.factor = 5),ID,b.add.scale=FALSE,
                         style="vertical",
                         genes=GENES,
                         cex.label=0,
                         cex.annot.label = cex.annot,b.show.just.group=TRUE,
                         cv = CV ,left.margin.prop = lm,cex.terminal.dots = 0.5,
                         seed=123576,legpos = NULL,
                         vspace.reserve=vspace.reserve,
                         b.extend.colours = b.extend.colours,
                         prioritise.group = prioritise.group,lwd=1.2,...
    );
    cord=tree$coords[match(1:length(tree$tip.label),tree$edge[,2]),]
    points(cord$b1,tree$top-cord$a1,col=tree$tip.color,cex=0.5,pch=19)
    legend("topleft",legend = sprintf("%2.0f",pdx$age.df$age),col=pdx$age.df$color,pch=19,title = "Age at Sample",bty = "n")
    ##browser()
    tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7,border=NA,txtcols=txtcols,force.count = force.count,b.add.lines = b.add.lines)#,font=1)
    #tree=add_heatmap(tree,heatmap=hm,cex.label = 0.7,border=NA,txtcols=txtcols,force.count = ifelse(dim(dff)[1]>10,-1,10),b.add.lines = b.add.lines,font=1)
    #lm=3*ff/(1-3*ff)
    if(b.add.bs){
      bs_labels2(tree,80)
    }
  }
  text(PD$INTERNAL_ID,x=length(tree$tip.label)/2,y=tree$top,pos = 3,cex=1.5,xpd=TRUE)
  tree
}



plot_all_trees_scaled_v2=function(PDD,bw,mar=c(0,2,0,2),m=c(1200,1400,1800),d=c(5,10,15),expand.top.tree="none",figure.path=sprintf("../export/figure_2%s.pdf",VERSION)){## box width fraction of total image height
  ## hard coding here ###
  pdat=data.frame(donor=names(PDD),n=sapply(PDD,function(PD) length(PD$pdx$tree_ml$tip.label)+1))
  pdat$row=c(1,1,1,2,2,2,3,3,3)
  rdat=data.frame(row=1:3,twidth=c(1,1,0.8))
  
  MAXCOUNT= pdat %>% group_by(row) %>% summarise(N=sum(n)) %>% (function(x) max(x$N)) ## = 202..
  u=1/(MAXCOUNT*1.20) ## set the unit width of each tip or colony
  pconfig=pdat %>%  
    left_join(pdat %>% group_by(row) %>% summarise(N=sum(n),c=n()) %>% mutate(W=N*u) %>% left_join(rdat) %>% mutate(lmargin=twidth-W)) %>% 
    mutate(w=twidth*(n+0.05*N)/(N+c*0.05*N),ff=w/(n*u)-1)#ff=lmargin/(c*n*u))
  print(pconfig)
  ###
  ypm=(1-bw*sum(d))/sum(m)
  
  ## relative heights of rows
  h=m*ypm+bw*d
  h=h/sum(h)
  vsr=bw*d/(m*ypm)
  #cat("vsr",vsr,"\n")

  vec=do.call("c",lapply(1:3, function(thisrow){
    pp=pconfig %>% filter(row==thisrow)
    w=pp$w
    if(pp$twidth[1]<1){
      w=c(w,1-sum(w))
    }
    cnt=round(100*w)
    cnt[length(cnt)]=100-sum(cnt[-length(cnt)])
    cnt
  }))
  
  vec2=do.call("c",lapply(1:length(vec),function(i) rep(i,vec[i])))
  ppp=matrix(vec2,nrow=3,byrow = TRUE)
  
  pdf(figure.path,w=12,h=8,pointsize = 8)
  layout(ppp,width=rep(1/100,100),heights = h)
  ff=pconfig$ff/(1+pconfig$ff) ## 
  zz=lapply(1:3,function(i) do_tree_paper_plot(PDD[[i]],ff=ff[i],force.count =d[1],ymax=m[1],vspace.reserve = vsr[1],mar=mar,expand.top.tree=expand.top.tree))
  
  zz=lapply(4:6,function(i) do_tree_paper_plot(PDD[[i]],ff=ff[i],force.count =d[2],ymax=m[2],vspace.reserve = vsr[2],mar=mar,expand.top.tree=expand.top.tree))
  zz=lapply(7:9,function(i) do_tree_paper_plot(PDD[[i]],ff=ff[i],force.count =d[3],ymax=m[3],vspace.reserve = vsr[3],mar=mar,expand.top.tree=expand.top.tree))
  add_legend(cex=1.5,bty="n")
  ## Add legend in final box
  dev.off()
}

add_legend=function(...){
  plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
  cdf=get_driver_scheme() %>% mutate(group=ifelse(number==5,"Copy Number Change",group)) %>% 
    mutate(group=ifelse(number==2,"BCR::ABL1",group)) %>% filter(number != 11) %>% 
    unique() %>% left_join(data.frame(group=c("BCR::ABL1","Copy Number Change","DNMT3A"),priority2=c(3,2,1))) %>% (function(x) x[order(x$priority2,decreasing = TRUE,na.last = TRUE),])
  
  #cdf=readRDS("driverscheme.RDS")
  legend("center",cdf$group,col=cdf$colour,pch=15,...)
}
