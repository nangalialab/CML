### Part Thre
## Generate signature related signature for paper.
## Generate figures for supplementary note 3. 

source("load_and_annotate_tree.R")
source("revised_sigs.R")
source("cml_sigs.R")
source("PHYLOFIT.R")
require(cloneRate)
require(broom)
require(patchwork)
require(ggpubr)
VERSION=readLines("../export/VERSION.txt")
phylores=readRDS("../cache/phylores.RDS")
phyloresNULL=readRDS("../cache/phyloresNULL.RDS")
PDD=readRDS("../cache/PDDxs.RDS")


add_bracket=function(gp,DFF,field1,field2,index=1){
  
  ##. Need to do change to support ggplot + function..
  dat=DFF %>% group_by(signature) %>% summarise(ub=max(ub)+index*(1-max(ub))/4)
  
  gp+geom_segment(data=dat,aes(x=ub,xend=ub),y=field1,yend=field2,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field1,yend=field1,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field2,yend=field2,col="black")
  
}
META_ANALYSIS=function(df){ranm=rma(df$mean,sei=(df$ub-df$lb)/(2*1.96));sfit=summary(ranm);data.frame(mean=sfit$beta[1],lb=sfit$ci.lb,ub=sfit$ci.ub,se=sfit$se)}

create_forest_plot=function(DF,XMAX=1.1){
  
  DF=DF %>% left_join(data.frame(category=c("Mutant","MutantI","Top100","Trunk","WildType"),
                                 newcat=c("CML","Early CML","Early life","Pre-CML Lineage","Wild Type"))) %>%
    mutate(category=newcat)
  
  DF$category=factor(DF$category,levels=rev(c("Wild Type","Pre-CML Lineage","Early life","Early CML","CML")))
  gp=ggplot(DF,aes(y=category,col=category))+geom_errorbar(aes(xmin=lb,xmax=ub),width=0.5)+
    geom_point(aes(x=mean))+facet_grid(rows=vars(signature))+theme_bw()+ggtitle("Cohort Level Signatures vs Branch Category")+
    ylab("")+xlab("Proportion")+scale_x_continuous(limits=c(0,XMAX),breaks=seq(0,1,0.1))
  
  index=1;
  dat=DF %>% group_by(signature) %>% summarise(ub=max(ub)+index*(XMAX-max(ub))/4)
  field1="Wild Type";field2="Pre-CML Lineage"
  D1=DF %>% filter(category==field1) %>% mutate(value1=mean,se1=se) %>% dplyr::select(signature,value1,se1)
  D1=D1 %>% left_join(DF %>% filter(category==field2) )
  D1=D1 %>% mutate(p=pnorm(abs(value1-mean)/sqrt(se1**2+se**2),lower.tail = FALSE))
  
  
  g5=gp+geom_segment(data=dat,aes(x=ub,xend=ub),y=field1,yend=field2,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field1,yend=field1,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field2,yend=field2,col="black")+
    geom_text(data=D1 %>% dplyr::select(-ub) %>% inner_join(dat),aes(x=ub,y=field1,label=sprintf("P=%3.2g",p)),nudge_x = 0.01,nudge_y=-0.5,hjust="left",col="black",size=2.5)
  
  if(!("CML" %in% DF$category)){
    return(g5)
  }
  index=3;dat=DF %>% group_by(signature) %>% summarise(ub=max(ub)+index*(XMAX-max(ub))/4)
  field1="Wild Type";field2="CML"
  D1=DF %>% filter(category==field1) %>% mutate(value1=mean,se1=se) %>% dplyr::select(signature,value1,se1)
  D1=D1 %>% left_join(DF %>% filter(category==field2) )
  D1=D1 %>% mutate(p=pnorm(abs(value1-mean)/sqrt(se1**2+se**2),lower.tail = FALSE))
  g6=g5+geom_segment(data=dat,aes(x=ub,xend=ub),y=field1,yend=field2,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field1,yend=field1,col="black")+
    geom_segment(data=dat,aes(x=ub-0.02,xend=ub),y=field2,yend=field2,col="black")+
    geom_text(data=D1 %>% dplyr::select(-ub) %>% inner_join(dat),aes(x=ub,y=field1,label=sprintf("P=%3.2g",p)),nudge_x = 0.01,nudge_y=-0.5,col="black",hjust="left",size=2.5)
  
  
  g6
}

compare_phylores=function(p1,p2,label1,label2){
  donors=names(p1)
  do.call("rbind",lapply(donors,function(donor){
    s1=summary(p1[[donor]]$res)$summary
    s2=summary(p2[[donor]]$res)$summary
    data.frame(S_A=s1["S",'50%'],S_lb_A=s1["S",'2.5%'],S_ub_A=s1["S",'97.5%'],
               S_B=s2["S",'50%'],S_lb_B=s2["S",'2.5%'],S_ub_B=s2["S",'97.5%'],donor=donor
    )
  }
  )
  )
}

plot_SBS18_rcsm=function(PDD){
  pcawg=get_extended_pcawg()
  zz2=lapply(PDD,get_clade_sigs,prior=pcawg[,c("SBS1","SBSblood")],extracats=c("MutantI","Top50","Top100"))#readRDS("../cache/sigs2.RDS")#
  zz3=lapply(PDD,get_clade_sigs,prior=pcawg[,c("SBS1","SBSblood","SBS18")],extracats=c("MutantI","Top50","Top100"))#readRDS("../cache/sigs3.RDS")#
  bb2=do.call("rbind",lapply(zz2, function(PD) PD$sigstats))
  bb3=do.call("rbind",lapply(zz3, function(PD) PD$sigstats))
  bba=rbind(bb3 %>% mutate(palette="1,Blood & 18"),bb2 %>% mutate(palette="1 & Blood"))
  BB=bba %>% filter(signature=="SBS1") %>% dplyr::select(category,donor,rec_csm,palette,total) %>% pivot_wider(names_from=palette,values_from=rec_csm)
  BB %>% filter(category!="Top50") %>% ggplot(aes(x=total,xend=total,y=`1 & Blood`,yend=`1,Blood & 18`,col=category,label=donor))+geom_segment(arrow = arrow(length = unit(0.2, "cm"),type="closed"))+theme_bw()+ylim(c(0.9,1))+ylab("Cosine Similarity")+scale_x_log10(limits=c(100,1e5),breaks=c(100,500,1000,5000,10000,50000,100000),labels=c(100,500,1000,5000,10000,"","100000"))+xlab("Total SNVs in Attribution")+ggtitle("Reconstruction Cosine Similarity vs SNV count:\nChange from Sig=SBS1,SBSblood to Sig=SBS1,SBSblood,SBS18")
  
}


CAT_ORDER=c("Mutant","MutantI","Top100","Trunk","WildType")
if(TRUE){
  dfcomb=read.table(sprintf("../export/per_donor_summary_v%s.txt",VERSION),head=TRUE,stringsAsFactors = FALSE,comment.char = "")
  dfcomb=dfcomb %>% mutate(donor=Donor)
  #CAT_ORDER=c("Mutant","MutantI","Top100","Trunk","WildType")
  #S2=ctsummary %>% mutate(signature=ctx) %>% dplyr::select(donor,category,signature,mean,lb,ub)
  S2=do.call("rbind",lapply(PDD,function(PD) get_ctx_counts(PD) )) %>% mutate(signature="C>T at CpG",category=type,mean=prop) %>% dplyr::select(donor,category,signature,mean,lb,ub,total)
  sbssummary=do.call("rbind",lapply( PDD,function(PD) PD$pdx$dat$info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`) ))
  
  S1=sbssummary %>% dplyr::select(donor,category,signature,mean,lb,ub,total)
  df=rbind(S1,S2) %>% mutate(type=factor(category,levels=rev(CAT_ORDER)))
  
  df2=df %>% filter(category %in% c("Mutant","MutantI")) %>% left_join(dfcomb %>% mutate(donor=Donor),by="donor") 
  donor.exclude=c("PD57333")
  gml=ggplot(df2,aes(x=latency_median,y=mean))+geom_point(aes(col=Donor))+geom_smooth(method="lm")+
    scale_color_manual(values=get_colors(dfcomb %>% filter(!(Donor %in% donor.exclude )),field = "Donor"))+theme_bw()+
    facet_grid(cols=vars(category),rows=vars(signature),scales="free_y")+ylab("Proportion")+xlab("Latency(Years)")
  tmp=df2 %>% as.data.frame() %>% group_by(signature,category) %>% summarise(tidy(lm(mean~latency_median)) %>% filter(term=="latency_median") %>% mutate(xmax=max(mean),xmin=min(mean)))
  #tmp=df2 %>% as.data.frame() %>% group_by(signature,category) %>% summarise(tidy(lm(Logit(mean)~latency_median)) %>% filter(term=="latency_median") %>% mutate(xmax=max(mean),xmin=min(mean)))
  
  gml=gml+geom_text(data=tmp %>% filter(p.value<0.05),aes(y=xmin+1.2*(xmax-xmin),x=4,label=sprintf("P=%3.2g",p.value)),hjust="left")
  pdf(sprintf("../figures/latency_vs_signature_v%s.pdf",VERSION),w=6,h=8)
  print(gml)
  dev.off()
  
  donor.by.S=dfcomb[order(dfcomb$S_ltt_median),] %>% pull("donor")
  g1=ggplot(df %>% filter(!(donor %in% donor.exclude)) %>% mutate(donor=factor(donor,levels=donor.by.S)),aes(y=donor,col=donor))+
    geom_errorbar(aes(xmin=lb,xmax=ub),width=0.5)+geom_point(aes(x=mean))+
    facet_grid(rows=vars(type),cols=vars(signature))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
    theme_bw()+ggtitle("Signature vs Donor within Branch Category")+ylab("")+xlab("Proportion")+scale_color_manual(values=get_colors(df %>% mutate(Donor=donor) %>% filter(!(donor %in% donor.exclude )),"Donor"))
  
  g2=ggplot(df %>% filter(!(donor %in% donor.exclude)),
            aes(y=type,col=type))+geom_errorbar(aes(xmin=lb,xmax=ub),width=0.5)+geom_point(aes(x=mean))+
    facet_grid(rows=vars(donor),cols=vars(signature))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+theme_bw()+ggtitle("Signature vs Category within Donor")+ylab("")+xlab("Proportion")
  
  g1all=ggplot(df %>% mutate(donor=factor(donor,levels=donor.by.S)),aes(y=donor,col=donor))+
    geom_errorbar(aes(xmin=lb,xmax=ub),width=0.5)+geom_point(aes(x=mean))+
    facet_grid(rows=vars(type),cols=vars(signature))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
    theme_bw()+ggtitle("Signature vs Donor within Branch Category")+ylab("")+xlab("Proportion")+scale_color_manual(values=get_colors(df %>% mutate(Donor=donor) %>% filter(!(donor %in% donor.exclude )),"Donor"))
  
  g2all=ggplot(df ,
               aes(y=type,col=type))+geom_errorbar(aes(xmin=lb,xmax=ub),width=0.5)+geom_point(aes(x=mean))+
    facet_grid(rows=vars(donor),cols=vars(signature))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+theme_bw()+ggtitle("Signature vs Branch Category within Donor")+ylab("")+xlab("Proportion")
  
  DF=df %>% group_by(signature,category) %>% summarise(META_ANALYSIS(.data)) #%>% left_join(dfcomb %>% mutate(donor=Donor),by="donor") 
  donor.exclude=c("PD57333","PD57332")
  gfa=create_forest_plot(DF)+
    ggtitle("Cohort Level Signatures vs Branch Category (All Donors)")
  gf=create_forest_plot(df %>% filter(!(donor %in% donor.exclude)) %>% group_by(signature,category) %>% summarise(META_ANALYSIS(.data)),XMAX = 1.13)+
    ggtitle("Cohort Level Signatures vs Branch Category (Exc. PD57333 and PD57332)")
  
  
  pdf(sprintf("../figures/figure_3e_signature_vs_category_cohort_v%s.pdf",VERSION),w=8,h=6)
  #print(gfa)
  print(gf)
  dev.off()
  if(TRUE){
    pdf(sprintf("../figures/signature_vs_donor_within_category_v%s.pdf",VERSION),w=10,h=10)
    print(g1)
    dev.off()
    
    pdf(sprintf("../figures/signature_vs_category_within_donor_v%s.pdf",VERSION),w=10,h=6)
    print(g2)
    dev.off()
    pdf(sprintf("../figures/signature_vs_donor_within_category_ALL_v%s.pdf",VERSION),w=10,h=12)
    print(g1all)
    dev.off()
    pdf(sprintf("../figures/signature_vs_category_within_donor_ALL_v%s.pdf",VERSION),w=10,h=8)
    print(g2all)
    dev.off()
    
  }
  pdf(sprintf("../figures/extended_figure_3a_reconstruction_csm_vs_sbspalette_v%s.pdf",VERSION),w=10,h=8)
  g1=plot_SBS18_rcsm(PDD)
  print(g1)
  dev.off()
  # plot vs upper bound
  pp=readRDS("../cache/phylores_smax10.RDS")
  dd=compare_phylores(phylores,pp)
  g1=ggplot(dd,aes(x=S_A,xmin=S_lb_A,xmax=S_ub_A,y=S_B,ymin=S_ub_B,ymax=S_lb_B,col=donor))+geom_point() +
    geom_errorbar(width=0.1,size=1)+geom_errorbarh(height=0.1,size=1)+scale_x_log10(labels = scales::percent,limits=c(0.1,1e11))+
    scale_y_log10(labels = scales::percent,limits=c(0.1,1e11))+geom_abline(slope=1,intercept=0)+theme_bw()+
    scale_color_manual(values=get_colors(dd %>% mutate(Donor=donor),"Donor"))+ggtitle("Growth Estimates: Max s=30 vs Max s=10")+
    ylab("S(Phylofit):Max s=10")+xlab("S(Phylofit):Max s=30")+theme(legend.title = element_blank())
  pdf(sprintf("../figures/extended_figure_5a_phylofit10_vs_phylofit30_v%s.pdf",VERSION),w=10,h=8)
  print(g1)
  dev.off()
  sel2=read.table(sprintf("../export/growth_rates_ALL_%s.txt",VERSION),head=TRUE,sep="\t")
  tmp=sel2 %>% filter(grepl(" BCR::ABL1$",LABEL))%>% mutate(Donor=donor)
  
  g1=ggplot(tmp ,aes(y=exp(estimate_cloneRateBirthDeathMCMC)-1,ymin=exp(lowerBound_cloneRateBirthDeathMCMC)-1,ymax=exp(upperBound_cloneRateBirthDeathMCMC)-1,x=exp(estimate_phylofit)-1,xmin=exp(lowerBound_phylofit)-1,xmax=exp(upperBound_phylofit)-1,col=Donor))+geom_point() +
    geom_errorbar(width=0.1,size=1)+geom_errorbarh(height=0.1,size=1)+scale_x_log10(labels = scales::percent,limits=c(0.1,1e11))+scale_y_log10(labels = scales::percent,limits=c(0.1,1e11))+
    geom_abline(slope=1,intercept=0)+theme_bw()+scale_color_manual(values=get_colors(tmp %>% mutate(Donor=donor),"Donor"))+xlab("S(Phylofit)")+ylab("S(BirthDeathMCMC)")
  
  pdf(sprintf("../figures/extended_figure_5b_phylofit_vs_johnsonBirthDeathMCMC_v%s.pdf",VERSION),w=10,h=8)
  print(g1)
  dev.off()
  
  sel=rbind(
    do.call("rbind",lapply(names(phyloresNULL),function(id) get_all_results(phyloresNULL[[id]]) %>% mutate(model="null",donor=id))),
    do.call("rbind",lapply(names(phylores),function(id) get_all_results(phylores[[id]]) %>% mutate(model="alt",donor=id)))
  )
  g1=sel %>% mutate(label=sprintf("%s_%s",method,model),type=ifelse(grepl("Moltime",method),"Molecular Time Based","Time Based"))  %>% ggplot(aes(y=label,xmin=lowerBound,xmax=upperBound,x=estimate,col=model,shape=type))+geom_errorbarh()+geom_point()+facet_wrap(~donor,scale="free_x")+theme_bw()+theme(legend.title=element_blank())+xlab("s: Instantaneous growth rate (per year)")+ylab("")+scale_color_manual(labels = c("Clade Specific Rate", "Global Rate"), values = c(alt="gray70",null="gray10"))
  pdf(sprintf("../figures/figure_S2_%s.pdf",VERSION),w=10,h=10)
  print(g1)
  dev.off()
}
## More supplementary figures
qqq=readRDS("../data/supp_note_2_benchmark_results.RDS")
## cap the difference - will only affect cloneRate
qqq=qqq %>% mutate(cap=ifelse(s_true>3,30,10)) %>% mutate(diff=ifelse(log(median+1)>cap,cap-s_true,log(median+1)-s_true)) %>% mutate(s_true=sprintf("s=%5.1f (S=%d%%)",s_true,S_true*100))
#qqq=qqq %>% mutate(diff=log(median+1)-s_true)
XLAB="n(Sampled Mutant Clade Size)"
grms=qqq  %>% group_by(s_true,method,N) %>% 
  summarise(rms=sqrt(mean(diff**2))) %>% filter(N>3) %>% 
  ggplot(aes(x=N,y=rms,col=method))+scale_y_log10()+geom_line()+xlab(XLAB)+
  geom_point()+facet_wrap(~s_true,scales="free_y")+theme_bw()+scale_x_continuous(breaks=unique(qqq$N))+
  ylab("RMS Error")+ggtitle("Root Mean Square Error for inference of s")+
  annotation_logticks(sides = "l")
gadiff=qqq %>% group_by(s_true,method,N) %>% 
  summarise(ad=median(abs(diff))) %>% filter(N>4) %>% 
  ggplot(aes(x=N,y=ad,col=method))+scale_y_log10()+geom_line()+xlab(XLAB)+
  geom_point()+facet_wrap(~s_true,scales="free_y")+theme_bw()+scale_x_continuous(breaks=unique(qqq$N))+
  ylab("Median Absolute Error")+ggtitle("Root Mean Square Error for inference of s")+
  annotation_logticks(sides = "l")
gmean=qqq %>% group_by(s_true,method,N) %>% 
  summarise(me=mean(diff)) %>% filter(N>4) %>% 
  ggplot(aes(x=N,y=me,col=method))+geom_line()+xlab(XLAB)+
  geom_point()+facet_wrap(~s_true,scales="free_y")+theme_bw()+scale_x_continuous(breaks=unique(qqq$N))+
  ylab("Mean Error")+ggtitle("Mean Error for inference of s")+
  annotation_logticks(sides = "l")+geom_hline(yintercept = 0)
pdf(sprintf("../figures/figure_S4_benchmark_rms_error_%s.pdf",VERSION),w=10,h=8)
print(grms)
dev.off()
pdf(sprintf("../figures/figure_S5_benchmark_mean_error_%s.pdf",VERSION),w=10,h=8)
print(gmean)
dev.off()
pdf(sprintf("../figures/figure_S5_benchmark_median_abs_error_%s.pdf",VERSION),w=10,h=8)
print(gadiff)
dev.off()
png(sprintf("../figures/figure_S4_benchmark_rms_error_%s.png",VERSION),w=480*1.5,h=480*1.2)
print(grms+theme_bw(base_size = 14))
dev.off()
png(sprintf("../figures/figure_S5_benchmark_mean_error_%s.png",VERSION),w=480*1.5,h=480*1.2)
print(gmean+theme_bw(base_size = 14))
dev.off()
png(sprintf("../figures/figure_S5_benchmark_median_abs_error_%s.png",VERSION),w=480*1.5,h=480*1.2)
print(gadiff+theme_bw(base_size = 14))
dev.off()
## rebuttal figures.
EXTRASEL=readRDS("../cache/extrasel.RDS")
rr=do.call("rbind",lapply(names(EXTRASEL),function(y) get_all_results(EXTRASEL[[y]]) %>% mutate(label=y))) %>% dplyr::select(-warningMessage)
rr=rr %>% mutate(LABEL=sprintf("%s:%s",label,label2)
)
rr=rr %>% mutate(S_ltt_median=exp(estimate)-1,   S_ltt_mean=NA,  S_ltt_lb=exp(lowerBound)-1,S_ltt_ub=exp(upperBound)-1)
R1=do.call("rbind",lapply(names(phylores),function(x) get_all_results(phylores[[x]]) %>% mutate(donor=x)))
R1=R1 %>% mutate(S_ltt_median=exp(estimate)-1,   S_ltt_mean=NA,  S_ltt_lb=exp(lowerBound)-1,S_ltt_ub=exp(upperBound)-1)
R1=R1 %>% mutate(label=sprintf("%s: BCR::ABL1",donor)) %>% mutate(label2=label)
tmp=rbind(R1 %>% dplyr::select(-warningMessage) %>% mutate(LABEL=label),rr %>% mutate(donor=LABEL) %>% mutate(LABEL=gsub("_[0-9]+","",LABEL)))
tmp=tmp %>% mutate(LABEL=gsub("_[0-9]+","",LABEL))

zz=read.table(sprintf("../export/growth_rates_ALL_%s.txt",VERSION),head=TRUE,sep="\t")
##
g0=ggplot(tmp %>% filter(!grepl("Moltime",method)) %>% mutate(method=gsub("cloneRate","",method)) %>% mutate(method=gsub("ML","maxLikelihood",method)) %>%
            filter(grepl(": BCR",LABEL)) %>% mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","maxLikelihood"))),
          aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method))+
  geom_point()+geom_errorbar()+facet_wrap(~LABEL,nrow=4,scales="free_y")+theme_bw()+
  theme(axis.text.x=element_blank())+ylab("s")+xlab("")

plts=lapply(list(list(methodA="phylofit",methodB="cloneRateBirthDeathMCMC"),
                 list(methodA="phylofit",methodB="cloneRateML"),list(methodA="cloneRateML",methodB="cloneRateBirthDeathMCMC")),function(x){ methodA=x$methodA;methodB=x$methodB;
zz %>% filter(grepl(": BCR::ABL1",LABEL)) %>%
dplyr::rename(estimate_A=paste0("estimate_",methodA),lowerBound_A=paste0("lowerBound_",methodA), upperBound_A=paste0("upperBound_",methodA),
estimate_B=paste0("estimate_",methodB),lowerBound_B=paste0("lowerBound_",methodB), upperBound_B=paste0("upperBound_",methodB)) %>% ggplot(aes(x=estimate_A,xmin=lowerBound_A,xmax=upperBound_A,y=estimate_B,ymin=lowerBound_B,ymax=upperBound_B,col=donor))+geom_point() +
geom_errorbar(width=0.1,size=1)+geom_errorbarh(height=0.1,size=1)+
geom_abline(slope=1,intercept=0)+theme_bw()+scale_color_manual(values=get_colors(tmp %>% mutate(Donor=donor),"Donor"))+
                   xlab(sprintf("s(%s)",gsub("ML","maxLikelihood",gsub("cloneRate","",methodA))))+
                   ylab(sprintf("s(%s)",gsub("ML","maxLikelihood",gsub("cloneRate","",methodB))))})
leg=get_legend(plts[[1]])
g2=ggarrange(g0+ theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=0.5)),
             plts[[1]],plts[[2]],plts[[3]],ncol=2,nrow=2,legend.grob = leg,
             legend = "right",
             labels=letters[1:4])

pdf(sprintf("../figures/figure_S2_%s.pdf",VERSION),w=10,h=10)
#plt=annotate_figure(g2,top="Comparison of methods for estimating Growth Rates of all BCR::ABL Expansions (s)")
print(g2)
dev.off()

png(sprintf("../figures/figure_S2_%s.png",VERSION),w=480*1.5,h=480*1.5)
#plt=annotate_figure(g2,top="Comparison of methods for estimating Growth Rates of all BCR::ABL Expansions (s)")
print(g2)
dev.off()
if(FALSE){
relabel=data.frame(LABEL=c("PD57335: BCR::ABL1","PD57335_REMOVE_LATE:BCR::ABL:Early Time Point: Remove 3 late coal.","PD57335yr:BCR::ABL:Time point 2"),
                   newlabel=c("BCR::ABL1: Early Time Point","BCR::ABL1:Early Time Point:\nRemove 3 late coalescences","BCR::ABL1:Later Time Point"),order=c(1,2,3))
relabel$newlabel=factor(relabel$newlabel,levels=relabel$newlabel[relabel$order])
gch57335=ggplot(tmp %>% filter(grepl("PD57335",LABEL)) %>% left_join(relabel) %>%
                  mutate(method=gsub("ML","maxLikelihood",gsub("cloneRate","",method))) %>% 
                  filter(method!="Moltime") %>% mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","maxLikelihood"))),
                aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method,col=method))+geom_point()+geom_errorbar()+facet_wrap(~newlabel,ncol=3)+theme_bw()+
  theme(axis.text.x=element_blank(),legend.position = "bottom")+ylab("s")+xlab("")
png("../figures/rebuttal_PD57335_2timepoints.png",w=480,h=480/1.5)
gch57335+scale_y_continuous(breaks=2:10)
dev.off()
gch51633=ggplot(tmp %>% filter(grepl("PD51633",LABEL)) %>% 
                  mutate(method=gsub("cloneRate","",method)) %>% 
                  filter(method!="Moltime") %>% 
                  mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","ML"))),
                aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method,col=method))+
  geom_point()+geom_errorbar()+facet_wrap(~LABEL,ncol=3)+
  theme_bw()+theme(axis.text.x=element_blank())+ylab("s")+xlab("")
pdf("../figures/rebuttal_PD51633_methods.pdf",w=6,h=4)
print(gch51633)
dev.off()
pdf("../figures/rebuttal_PD51633_pairs.pdf",w=6,h=6)
pairs(phylores$PD51633$res,pars=c("s","tm"))
dev.off()
## redo 
if(FALSE){
  param=PARAM %>% filter(donor=="PD51633")
  tree=PDD$PD51633$fit$poisson_tree$altmodel$ultratree
  r2=fit_clade(tree,
               node=param$node,
               nmutcolony =-1,
               nwtcolony =-1,
               maxt=60,
               stan.control = list(adapt_delta=param$adapt_delta,max_treedepth=14),
               maxSYear = param$maxS,
               niter = 10000)
}
gch57334=ggplot(tmp %>% filter(grepl("PD57334",LABEL)) %>% 
                  mutate(method=gsub("cloneRate","",method)) %>% 
                  filter(method!="Moltime") %>% 
                  mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","ML"))),
                aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method,col=method))+
  geom_point()+geom_errorbar()+facet_wrap(~LABEL,ncol=3)+
  theme_bw()+theme(axis.text.x=element_blank())+ylab("s")+xlab("")
pdf("../figures/rebuttal_PD57334_methods.pdf",w=6,h=4)
print(gch57334)
dev.off()

gch57332=ggplot(tmp %>% filter(grepl("PD57332",LABEL)) %>% 
                  mutate(method=gsub("cloneRate","",method)) %>% 
                  filter(method!="Moltime") %>% 
                  mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","ML"))),
                aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method,col=method))+
  geom_point()+geom_errorbar()+facet_wrap(~LABEL,ncol=3)+
  theme_bw()+theme(axis.text.x=element_blank(),legend.position = "bottom",legend.title = element_blank())+ylab("s")+xlab("")

pdf("../figures/rebuttal_PD57332_methods_clades.pdf",w=8,h=4)
print(gch57332)
dev.off()

gch51635=ggplot(tmp %>% filter(grepl("PD51635",LABEL) & grepl("BCR",LABEL)) %>% 
                  mutate(method=gsub("cloneRate","",method)) %>% 
                  filter(method!="Moltime") %>% 
                  mutate(method=factor(method,levels=c("phylofit","BirthDeathMCMC","ML"))),
                aes(ymin=lowerBound,y=estimate,ymax=upperBound,x=method,col=method))+
  geom_point()+geom_errorbar()+facet_wrap(~LABEL,ncol=3)+
  theme_bw()+theme(axis.text.x=element_blank(),legend.position = "bottom",legend.title = element_blank())+ylab("s")+xlab("")

pdf("../figures/rebuttal_PD51635_subclades.pdf",w=6,h=4)
print(gch51635)
dev.off()
}

