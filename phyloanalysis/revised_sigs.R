require("MutationalPatterns")
require("hdp")
require("dplyr")
require("tibble")
require("stringr")
if(!exists("SPECIES_BUILD")){
  stop("Please create SPECIES (Hsapiens, Mmusculus ) SPECIES_BUILD and set to e.g. hg19,mm10")
}
if(!exists("ref_genome") ){
  ref_genome= sprintf("BSgenome.%s.UCSC.%s",SPECIES,SPECIES_BUILD)
  library(ref_genome, character.only = TRUE)
}

if(SPECIES_BUILD=="hg38"){
  csvcosmic="https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh38.txt"
}else{
  if(SPECIES_BUILD=="hg19"){
    csvcosmic="https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.2_SBS_GRCh37.txt"
  }else{
    stop(sprintf("unsupported SPECIES_BUILD=%s",SPECIES_BUILD))
  }
}

get_pcawg60=function(csvfile=csvcosmic){
  res=read.table(csvfile,head=TRUE,stringsAsFactors = FALSE,sep="\t")  %>% column_to_rownames("Type")
  res[MutationalPatterns:::TRIPLETS_96,]
}

##Get mutation count matrix from VCFs
get_mut_matrix_from_vcfs=function(vcf_files,sample_names){
  vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
  mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
}

extract_details_from_hdpres=function(res,priorsigs,mutcount,b.is.denovo=TRUE,threshold=0.9){
  
  N=dim(mutcount)[2]
  ##check integrity of mutcount
  dps=dp(final_hdpState(chains(res)[[1]]))
  chkcount=tail(sapply(dps, function(x) x@numdata),N)
  idx=which(colSums(mutcount)!=chkcount)
  if(length(idx)>0){
    browser()
    stop("matrix mismatch!")
  }
  comp_distn= comp_categ_distn(res)
  dp_distn=comp_dp_distn(res)
  sigs=comp_distn$mean
  cdc=comp_distn$cred.int
  sigs_lb95=t(sapply(1:length(cdc),function(x) cdc[[x]][1,]))
  sigs_ub95=t(sapply(1:length(cdc),function(x) cdc[[x]][2,]))
  contr=tail(dp_distn$mean,N)
  cdp=dp_distn$cred.int
  contr_lb95=tail(t(sapply(1:length(cdp),function(x) cdp[[x]][1,])),N)
  contr_ub95=tail(t(sapply(1:length(cdp),function(x) cdp[[x]][2,])),N)
  #Only keep signatures that have non-zero lower bound in at least one sample...
  idx=which(apply(contr_lb95,2,max)>0)
  if(length(idx)==0){
    stop("No significant sigs found!")
  }
  if(length(idx)==1){
    sigs_lb95=matrix(sigs_lb95[idx,],nrow=1)
    sigs_ub95=matrix(sigs_ub95[idx,],nrow=1)
    sigs=matrix(sigs[idx,],nrow=1)
    contr_lb95=matrix(contr_lb95[,idx],ncol=1)
    contr_ub95=matrix(contr_ub95[,idx],ncol=1)
    contr=matrix(contr[,idx],ncol=1)
  }else{
    sigs_lb95=sigs_lb95[idx,]
    sigs_ub95=sigs_ub95[idx,]
    sigs=sigs[idx,]
    contr_lb95=contr_lb95[,idx]
    contr_ub95=contr_ub95[,idx]
    contr=contr[,idx]
  }
  rownames(contr)=colnames(mutcount)
  if(b.is.denovo){
    ##Try ma
    csm=cos_sim_matrix(t(sigs),priorsigs)
    labels=apply(csm,1,function(x){
      idx=which.max(x)
      if(x[idx]>threshold){
        colnames(csm)[idx]
      }else{
        NA
      }})
    idx=which(!is.na(labels))
    if(length(idx)>0){
      rownames(sigs)[idx]=labels[idx]
    }
    idx2=which(is.na(labels))
    if(length(idx2)>0){
      rownames(sigs)[idx2]=sprintf("N%s",1:length(idx2))
    }
    csmmax=apply(csm,1,function(x){max(x)})
    names(csmmax)=rownames(sigs)
  }else{
    csmmax=NULL
  }
  colnames(sigs)=rownames(mutcount)
  colnames(contr)=rownames(sigs)
  rownames(contr_lb95)=rownames(contr)
  rownames(contr_ub95)=rownames(contr)
  colnames(contr_lb95)=colnames(contr)
  colnames(contr_ub95)=colnames(contr)
  cnt=colSums(mutcount)
  sigcount=t(sigs) %*% diag(as.vector(cnt %*% contr))
  activities=diag(cnt) %*% contr
  rownames(activities)=rownames(contr)
  reconstruction=contr %*% sigs
  sample_csm=sapply(1:dim(reconstruction)[1],function(i) cos_sim(reconstruction[i,],mutcount[,i]))
  list(signatures=t(sigs),
       sigs_lb95=t(sigs_lb95),
       sigs_ub95=t(sigs_ub95),
       contr=contr,
       contr_lb95=contr_lb95,
       contr_ub95=contr_ub95,
       reconstruction=reconstruction,
       csm_max=csmmax,
       sample_csm=sample_csm,
       mutcount=mutcount,
       sigcount=sigcount,
       activities=activities
  )
}

##### HDP functions...

##Denove
runHDP_denovo=function(mut_count,burnin=10000,spacing=500,n_per_chain=250,seed=-1,mtype=96){
  mut_count=t(mut_count)
  n_samp=dim(mut_count)[1]
  n_cat=dim(mut_count)[2]
  hdp_mut <- hdp_init(ppindex=c(0,1,rep(2,n_samp)),
                      cpindex=c(1,2,rep(3,n_samp)),
                      hh=rep(1,mtype),
                      alphaa=rep(1,3),
                      alphab=rep(1,3))
  hdp_mut = hdp_setdata(hdp_mut,
                        dpindex=3:(n_samp+2),
                        mut_count)
  if(seed<0){
    seed_offset=round(as.numeric(Sys.time()))
  }else{
    seed_offset=seed
  }
  ##Now run the sampling chains
  chlist=mclapply(1:4,function(i){
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=seed_offset+i*2000)
    hdp_posterior(hdp_activated,
                  burnin=burnin,
                  n=n_per_chain,
                  space=spacing,
                  cpiter=3,
                  seed=seed_offset+i*1e5)
  },mc.cores = 4)
  mut_multi=hdp_multi_chain(chlist)
  hdp_extract_components(mut_multi)
}




runHDP_projection=function(mut_count,prior_sigs=get_cosmic(),prior_strength=10000,seed=-1,spacing=250,n_per_chain=250,burnin=10000){
  mut_count=t(mut_count)
  n_samp=dim(mut_count)[1]
  n_cat=dim(mut_count)[2]
  nps=ncol(prior_sigs)
  if(seed<0){
    seed_offset=round(as.numeric(Sys.time()))
  }else{
    seed_offset=seed
  }
  ##browser()
  mpn_prior=hdp_prior_init(prior_distn = prior_sigs,
                           prior_pseudoc = rep(prior_strength, nps),
                           hh=rep(1, 96),
                           alphaa=c(1, 1),
                           alphab=c(1, 1))
  
  
  mpn_prior=hdp_addconparam(mpn_prior,
                            alphaa = c(1,1),
                            alphab = c(1,1))
  
  mpn_prior=hdp_adddp(mpn_prior,
                      numdp = dim(mut_count)[1]+1,
                      ppindex = c(1, rep(1+nps+1, dim(mut_count)[1])),
                      cpindex = c(3, rep(4, dim(mut_count)[1])))
  
  mpn_prior=hdp_setdata(mpn_prior,
                        dpindex = (1+nps+1)+1:dim(mut_count)[1],
                        mut_count)
  
  chlist =mclapply(1:4,function(i){
    cat("#")
    mpn_activated <- dp_activate(mpn_prior,
                                 dpindex = (1+nps+1)+0:dim(mut_count)[1],
                                 initcc = nps+5,
                                 seed = seed_offset+i*10000)
    cat("processing posterior\n")
    hdp_posterior(mpn_activated,
                  burnin = burnin,  ##5000
                  n = n_per_chain,
                  space = spacing,
                  cpiter = 3,
                  seed = seed_offset+i*1e6)
  },mc.cores = 4)
  cat("finished processing posterior\n")
  mpn_multi=hdp_multi_chain(chlist)
  hdp_extract_components(mpn_multi)
  
}
##Denove
runHDP_denovo_structured=function(mut_count,id,burnin=10000,spacing=500,n_per_chain=250,seed=-1){
  mut_count=t(mut_count)
  n_samp=dim(mut_count)[1]
  n_cat=dim(mut_count)[2]
  df=as.data.frame(table(id))
  #browser()
  dfr=df %>% filter(Freq>1)
  n_ind=length(df$id)
  ##ppindices associated with groups of counts associated with one individual
  ppex=rep(seq(3,2+length(dfr$id)),dfr$Freq)
  ccex=rep(seq(3,2+length(dfr$id)),dfr$Freq)+1
  hdp_mut <- hdp_init(ppindex=c(0,1,rep(2,n_ind),ppex),
                      cpindex=c(1,2,rep(3,n_ind),ccex),
                      hh=rep(1,96),
                      alphaa=rep(1,3+length(dfr$id)),
                      alphab=rep(1,3+length(dfr$id)))
  hdp_mut = hdp_setdata(hdp_mut,
                        dpindex=(3+length(dfr$id)):(3+length(dfr$id)+n_samp-1),
                        mut_count)
  if(seed<0){
    seed_offset=round(as.numeric(Sys.time()))
  }else{
    seed_offset=seed
  }
  ##Now run the sampling chains
  chlist=mclapply(1:4,function(i){
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=seed_offset+i*2000)
    hdp_posterior(hdp_activated,
                  burnin=burnin,
                  n=n_per_chain,
                  space=spacing,
                  cpiter=3,
                  seed=seed_offset+i*1e5)
  },mc.cores = 4)
  mut_multi=hdp_multi_chain(chlist)
  hdp_extract_components(mut_multi)
}


runHDP_denovo_structured2=function(mut_count,id,burnin=10000,spacing=500,n_per_chain=250,seed=-1){
  mut_count=t(mut_count)
  n_samp=dim(mut_count)[1]
  n_cat=dim(mut_count)[2]
  df=as.data.frame(table(id))
  #browser()
  dfr=df %>% filter(Freq>1)
  n_ind=length(df$id)
  ##ppindices associated with groups of counts associated with one individual
  ppex=rep(seq(3,2+length(dfr$id)),dfr$Freq)
  ccex=rep(seq(3,2+length(dfr$id)),dfr$Freq)+1
  hdp_mut <- hdp_init(ppindex=c(0,1,rep(2,n_ind),ppex),
                      cpindex=c(1,2,rep(3,n_ind),ccex),
                      hh=rep(1,96),
                      alphaa=c(rep(1,3),rep(100,length(dfr$id))), ## Wack up the concentration prior hopefully to better pool the sample..
                      #alphaa=rep(1,3+length(dfr$id)),
                      alphab=rep(1,3+length(dfr$id)))
  hdp_mut = hdp_setdata(hdp_mut,
                        dpindex=(3+length(dfr$id)):(3+length(dfr$id)+n_samp-1),
                        mut_count)
  if(seed<0){
    seed_offset=round(as.numeric(Sys.time()))
  }else{
    seed_offset=seed
  }
  ##Now run the sampling chains
  chlist=mclapply(1:4,function(i){
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=seed_offset+i*2000)
    hdp_posterior(hdp_activated,
                  burnin=burnin,
                  n=n_per_chain,
                  space=spacing,
                  cpiter=3,
                  seed=seed_offset+i*1e5)
  },mc.cores = 4)
  mut_multi=hdp_multi_chain(chlist)
  hdp_extract_components(mut_multi)
}

get_signatures_color_scheme=function(signames){
  cols=RColorBrewer::brewer.pal(9,"Set1")
  cols[9]="lightgrey"
  scheme=data.frame(
    signames=c("SBS1","SBS5","SBS32","SBS23","SBS2","SBS8","SBS9","SBS40","SBS19","SBS_5_40","SBS_19_23","SBS_5_9_40"),
    col=c(cols,cols[2],cols[9],cols[2]),stringsAsFactors = FALSE
  )
  idx=match(signames,scheme$signames)
  #browser()
  if(length(which(is.na(idx))>0)){
    cat("Can't match all signatures to preferred colour scheme. Using General one instead...\n")
    return(hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0,direction = 1)(length(signames)))
  }else{
    return(scheme$col[idx])
  }
}



expanded_plot=function(contribution,title,extracol=NULL,css=NULL,lb=NULL,ub=NULL,cex.label=1,renorm=TRUE){
  k=dim(contribution)[1]
  N=dim(contribution)[2]
  if(renorm){
    qm=sweep(contribution,2,colSums(contribution),"/")
    if(!is.null(lb)){
      lb=sweep(lb,2,colSums(contribution),"/")
      ub=sweep(ub,2,colSums(contribution),"/")
    }
  }else{
    qm=contribution
  }
  cols=get_signatures_color_scheme(rownames(contribution))
  
  
  par(mar=c(5, 4, 4, 2) + 0.1+c(2,4,0,0))
  plot(NULL,xlim=c(0,N),ylim=c(0,k),xaxt="n",yaxt="n",
       main=title,xlab="",ylab="",cex.main=2)
  for(i in 1:k){
    for(j in 1:N){
      rect(xleft=j-0.5,xright=j+0.5,ybottom=i-1,ytop=i+qm[i,j]-1,col=cols[i],border=NA)
      if(!is.null(lb)){
        rect(col="black",xleft=j-0.02,xright=j+0.02,ybottom = i+lb[i,j]-1,ytop=i+ub[i,j]-1,border=NA)
      }
    }
  }
  axis(2,at=seq(0.5,k),labels=rownames(qm),las=2)
  if(!is.null(extracol)){
    axis(1,at=seq(1,N),labels=FALSE,las=2)
    mtext(colnames(qm),side=1,line=1,at=seq(1,N),col=extracol,las=2,cex.axis=cex.label)
  }else{
    axis(1,at=seq(1,N),labels=colnames(qm),las=2,cex.axis=cex.label)
  }
  arrows(y0=par("usr")[3],y1=par("usr")[4],x0=0.5:(N+1),length=0,lwd=0.5)
  arrows(y0=0:k,y1=0:k,x1=par("usr")[1],x0=N+0.5,length=0,lwd=0.5)
}

addSpacedAxisVert=function(at,labels,text.size=1,min.space.factor=0.4,xfac=0.05,x0=NULL){
  #browser()
  xl=par("usr")[1:2]
  yl=par("usr")[3:4]
  par(xpd=NA)
  if(is.null(x0)){
    x0=xl[1]
  }
  xm=x0-(xl[2]-x0)*xfac
  
  n=length(at)
  
  
  #d.to=seq(xl[1]+((xl[2]-xl[1])/(2*(n))),xl[2],(xl[2]-xl[1])/(n))
  tick.at=at
  
  d.to=applyMinimumSpacing(tick.at,min.space.factor*(yl[2]-yl[1])/length(at))
  arrows(y0=tick.at,y1=d.to,x1=xm,x0=xl[1],length=0)
  #TODO!
  
  text(y=d.to,x=xm,labels=labels,srt=0,cex=text.size,pos=2)
}
applyMinimumSpacing=function(test,min.space,test.min=min(test)-min.space,test.max=max(test)+2*min.space){
  prev=0
  n=length(test)
  for(i in 1:n){
    if((test[i]-prev)>min.space){
      prev=test[i]
    }else{
      test[i]=prev+min.space
      prev=test[i]
    }
  }
  if(max(test)>test.max){
    ##TODO deal with overextension on the RHS of plot
  }
  test
}


add_signature_decomp=function(pdx,##<< PDX object or list including details matrix
                              tree,##<< enhanced phylo returned from plot_tree
                              node,##<< Tree node - maps to pdx$details$node
                              control,##<< Control parameters.. Requires bfield
                              ...)
{
  ##tree,node,res,ptarget){
 # browser()
  info=get_edge_info(pdx,tree,node)
  i=which(tree$edge[,2]==node)
  
  decomp=pdx$sigbynode[[i]]$contr
  csm=pdx$sigbynode[[i]]$csm
  control=add_defaults(control,defaults=list( maxlen=tree$top*0.2),
                       mandatory_fields=c("col.scheme"))
  
  
  y0=info$yb
  y1=info$yt
  
  x=info$x
  mh=y1-y0
  maxlen=control$maxlen
  #minlen=mh;##100
  maxlen=min(maxlen,mh)
  
  
  #if(mh<minlen){
  #  return(NULL)
  #}
  #frac=minlen/mh
  ##dd=0.5*(1-frac) ##0.01
  start=y1-0.01*mh
  start=y1-(mh-maxlen)*0.5
  start=y0+maxlen##ifelse(maxlen<(mh-0.1*maxlen),1.1*maxlen,maxlen)
  
  unit=0.98*mh##frac
  unit=maxlen
  #browser()
  if(is.na(decomp[1])){
    rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*1,col="grey",border=NA)
    return(NULL)
  }
  sigs=rownames(decomp)
  col.scheme=control$col.scheme
  cols=col.scheme$COL[match(sigs,col.scheme$SIG)]
  for(i in 1:length(decomp)){
    if(decomp[i]>0){
      rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*decomp[i],col=cols[i],border=NA)
      start=start-unit*decomp[i]
    }
  }
  if(control$b.include.csm.text){
    text(x-0.25,y0+maxlen-0.5*unit,sprintf("%3.2f",csm))
  }
  NULL
}

add_sig_to_pdx=function(PD,prior,mutprob=NULL,minmut=50){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  ##data.frame(SIG=sigs,count=NA)
  if(is.null(mutprob)){
    ##fit using mutational patterns...
    for(i in 1:length(PD$pdx$tree_ml$edge.length)){
      cat("#")
      node=PD$pdx$tree_ml$edge[i,2]
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc) & sum(mutc)>minmut){
        qf=quickfit(as.vector(mutc),prior)
        sigbynode[[i]]=list(node=node,contr=qf$contribution,csm=qf$csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
    }
    cat("\n")
    PD$pdx$dat$sigbynode=sigbynode
    PD$pdx$dat$siglist=sigs
    return(PD)
  }else{
    stop("Please check this is applicable to current setup!")
    sigs=colnames(mutprob)[-c(1,2)]
    prior=prior[,sigs]
    ##Expect mutprob to have form sample <PATIENT>.<node index> (wildtype =0, others are presented in same order as PD$nodes$node)
    descdf=data.frame(id=sprintf("%s.%s",PD$patient,0:length(PD$nodes$driver3[PD$nodes$status>=0])),desc=c("WT",PD$nodes$driver3[PD$nodes$status>=0]))
    mutprob=mutprob %>% filter(sample %in% descdf$id)
    ##Save in the form 
    tmp=markup_tree(PD$pdx$tree_ml,PD$nodes$node[PD$nodes$status>=0])
    tmp$sample=sprintf("%s.%s",PD$patient,tmp$type)
    sigbynode=list()
    for(i in 1:length(tmp$sample)){
      cat("#")
      idx=which(mutprob$sample==tmp$sample[i])
      node=PD$pdx$tree_ml$edge[i,2]
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% as.matrix(mutprob[idx,sigs]))
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      
    }
    cat("\n")
    PD$pdx$dat$sigbynode=sigbynode
    PD$pdx$dat$siglist=sigs
    return(PD)
  }
  
}

get_ctx_config=function(){
  df=data.frame(mut_type=c("C>A","C>G","C>T other","C>T at CpG","T>A","T>C","T>G","INDEL"),
                mut_type7=c("Other","Other","Other","C>T at CpG","Other","Other","Other","INDEL"),
              name=c("CA","CG","CGC","CGO","TA","TC","TG","INDEL"),
              colors=c(MutationalPatterns:::COLORS7,"orange"))
  df=rbind(df %>% filter(mut_type=="C>T at CpG"),df %>% filter(mut_type!="C>T at CpG"))
}
add_ctx_to_pdx=function(PD,b.all=FALSE,b.include.indel=TRUE){
  ## Just use the signatures template to add in context 
  PD$pdx$dat$details=add_trin(PD$pdx$dat$details)
  PD$pdx$dat$details=PD$pdx$dat$details %>% mutate(mut_type=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",
                                                                   ifelse(TYPE=="INDEL",TYPE,
                                                                          ifelse(grepl("\\[C>T\\]",trin),"C>T other",substr(trin,3,5)))
                                                                   )
                                                   )
  sigbynode=list()
  df=get_ctx_config() 
  if(!b.include.indel){
    df=df %>% filter(mut_type!="INDEL")
  }
  sigs=df$label
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  ##data.frame(SIG=sigs,count=NA)
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
      cat("#")
      NODE=PD$pdx$tree_ml$edge[i,2]
      df2=df %>% left_join(PD$pdx$dat$details %>% filter(node==NODE) %>% group_by(mut_type) %>% summarise(N=n()),by="mut_type") %>%
        mutate(N=ifelse(is.na(N),0,N))
      N=sum(df2$N)
      df2$contr=df2$N/sum(df2$N)
      #browser()
      contr=matrix(df2$contr,ncol=1)
      rownames(contr)=df2$mut_type
      sigbynode[[i]]=list(node=NODE,contr=contr,csm=1,count=N)
  }
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  PD
}

get_cohort_count_matrix_with_short_branches_grouped=function(PDD,min.count=50,early.branch.threshold=Inf,mc.cores=8){
  do.call("cbind",mclapply(PDD,function(PD) 
    get_pd_count_matrix_with_short_branches_grouped(
      PD,
      min.count=min.count,
      early.branch.threshold=early.branch.threshold),mc.cores = mc.cores)
    )
}

get_pd_count_matrix_with_short_branches_grouped=function(PD,min.count=50,early.branch.threshold=Inf){
  nodes=PD$pdx$tree_ml$edge[,2]
  mats=lapply(nodes,function(x) get_mut_matrix_for_node(PD$pdx,x))
  nmut=sapply(mats,function(x) sum(x))
  ## Identify which nodes
  nh=nodeHeights(PD$pdx$tree_ml)
  idx.young=which(nmut<min.count & nmut>0 & nh[,2]<early.branch.threshold)
  idx.old=which(nmut<min.count & nmut>0 & nh[,2]>=early.branch.threshold)
  if(length(idx.young)>0){
    young=rowSums(do.call("cbind",mats[idx.young]))
    mats[[paste0(PD$patient,".youngSB")]]=young
  }
  if(length(idx.old)>0){
    old=rowSums(do.call("cbind",mats[idx.old]))
    mats[[paste0(PD$patient,".oldSB")]]=old
  }
  idx.drop=c(idx.young,idx.old)
  if(length(idx.drop)>0){
    mats=mats[-idx.drop]
  }
  out=do.call("cbind",mats)
  out
}


get_cohort_count_matrix=function(PDD){
  do.call("cbind",lapply(PDD,get_pd_count_matrix))
}

get_pd_count_matrix=function(PD){
  nodes=PD$pdx$tree_ml$edge[,2]
  out=do.call("cbind",lapply(nodes,function(x) get_mut_matrix_for_node(PD$pdx,x)))
  out
}

get_mut_matrix_for_node=function(pdx,node){
  idx=get_idx_for_node(pdx$dat,node = node)
  if(length(idx)==0){
    return(NULL)
  }
  out=with(pdx$dat$details[idx,] %>% filter(TYPE=="SNV"),get_mut_matrix(Chrom,Pos,Ref,Alt))
  colnames(out)=sprintf("%s.%s",pdx$meta$prefix,node)
  out
}

quickfit=function(sig,priorsigs,method="MP"){
  if(method=="MP"){
    fitfunc=MutationalPatterns:::fit_to_signatures
  }else{
    fitfunc=fixed_weight_em
  }
  ##browser()
  decomp=fitfunc(matrix(sig,ncol=1),priorsigs)
  decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
  decomp$csm=cos_sim(decomp$reconstructed[,1],sig)
  decomp
}


extract_details_from_hdpx=function(hdpxres,mc,priorsigs,threshold=0.9){
  rownames(hdpxres$signature)=rownames(mc)
  hdpxres$exposureCounts=hdpxres$exposureProbs %*% diag(colSums(mc))
  colnames(hdpxres$exposureCounts)=colnames(hdpxres$exposureProbs)
  rownames(hdpxres$signature)=rownames(mc)
  tmp=relabel_sigs(hdpxres$signature,priorsigs,threshold)
  hdpxres$signature=t(tmp$sig)
  rownames(hdpxres$exposureCounts)=colnames(hdpxres$signature)
  #browser()
  reconstruction= hdpxres$signature %*% hdpxres$exposureCounts
  sample_csm=sapply(1:dim(reconstruction)[2],function(i) cos_sim(reconstruction[,i],mc[,i]))
  names(sample_csm)=colnames(hdpxres$exposureCounts)
  hdpxres$reconstruction=reconstruction
  hdpxres$sample_csm=sample_csm
  hdpxres$priorsigs=tmp
  hdpxres
}


relabel_sigs=function(sigs,priorsigs,threshold){
  ##Try ma
  sigs=t(sigs)
  csm=cos_sim_matrix(t(sigs),priorsigs)
  labels=apply(csm,1,function(x){
    idx=which.max(x)
    if(x[idx]>threshold){
      colnames(csm)[idx]
    }else{
      NA
    }})
  idx=which(!is.na(labels))
  if(length(idx)>0){
    rownames(sigs)[idx]=labels[idx]
  }
  idx2=which(is.na(labels))
  if(length(idx2)>0){
    rownames(sigs)[idx2]=sprintf("N%s",1:length(idx2))
  }
  csmmax=apply(csm,1,function(x){max(x)})
  names(csmmax)=rownames(sigs)
  colnames(sigs)=rownames(priorsigs)
  list(sig=sigs,
       csmmax=csmmax,
       csmmax_name=colnames(priorsigs)[apply(csm,1,function(x){which.max(x)})]
  )
}

get_extended_pcawg=function(){
  lym=read.table("../data/lympho_sbs_sigs.txt",head=T,stringsAsFactors = FALSE)
  pcawg=get_pcawg60()
  pcawg=cbind(pcawg,SBSblood=lym[,"SBSblood"]) %>% as.matrix()
  pcawg
}

get_sig_color_scheme_ORIG=function(){
  ## Save this to a file
  if(file.exists("../data/signature_scheme.txt")){
    scheme=read.table("../data/signature_scheme.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t",comment.char = "")
    if(any(is.na(match(c("COL","SIG"),names(scheme))))){
      stop("Invalid signature scheme - need COL and SIG columns")
    }
  }else{
    cols=RColorBrewer::brewer.pal(8,"Set1")
    scheme=data.frame(SIG=c("SBS1","SBSblood","SBS18"),COL=cols[c(1:2,8)])
  }
  scheme
                    
}

plot_cml_sig_tree=function(PD,prior=get_extended_pcawg()[,c("SBS1","SBSblood")],b.ultra=FALSE,txtitle=PD$patient){
  PD2=add_cml_sig_to_pdx(PD,prior=prior)
  colorscheme=get_sig_color_scheme()
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = 0.3)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  browser()
  tree=plot_signature_annotated_tree(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  ##
  inf=get_sig_info(PD)
  idx=which(rownames(inf$counts)=="C>T at CpG")
  ct_at_cpg=inf$counts[idx,]/colSums(inf$counts)
  bounds=sapply(1:3,function(i) binom.test(inf$counts[idx,i],colSums(inf$counts)[i])$conf.int)
  colnames(bounds)=colnames(inf$counts)
  xm=abs(par("usr")[1])
  ym=par("usr")[c(3,4)]
  ymm=mean(ym)
  h=diff(ym)*0.5
  w=0.2
  gap=0.1
  col.err="black"
  col.cpg="darkorange"
  col.other="white"
  f=ct_at_cpg[1,"WildType"]
  ##browser()
  segments(x0 =-xm+0.09*xm,y0=ymm-(h/2),y1=ymm+(h/2),lwd=1)
  for(x in c(0,0.25,0.5,0.75,1)){
    segments(x0=-xm+0.05*xm,x1=-xm+0.09*xm,y0=ymm-(h/2)+(h*x),lwd=1)
    text(sprintf("%3.2f",x),x=-xm+0.05*xm,y=ymm-(h/2)+(h*x),pos=2,cex=0.8,xpd=NA,offset=0.1)
  }
  
  rect(xleft=-xm+0.1*xm,xright=-xm+(0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
  rect(xleft=-xm+0.1*xm,xright=-xm+(0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
  segments(x0 = -xm+(0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"WildType"],y1 = ymm-(h/2)+h*bounds[2,"WildType"],lwd=2,lend=2,col=col.err)
  text("WT",x=-xm+(0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
  
  f=ct_at_cpg[1,"Trunk"]
  rect(xleft=-xm+(0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
  rect(xleft=-xm+(0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
  text("Trunk",x=-xm+(0.1+w+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
  segments(x0 = -xm+(0.1+w+0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"Trunk"],y1 = ymm-(h/2)+h*bounds[2,"Trunk"],lwd=2,lend=2,col=col.err)
  text("C>T @ CpG",x=-xm+(0.1+w+0.1+0.5*w)*xm,y=ymm+(h/2),pos=3,cex=1.2)
  
  f=ct_at_cpg[1,"Mutant"]
  rect(xleft=-xm+(0.1+w+0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
  rect(xleft=-xm+(0.1+w+0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w+0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
  text("BCR::ABL1",x=-xm+(0.1+w+0.1+w+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
  segments(x0 = -xm+(0.1+w+0.1+w+0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"Mutant"],y1 = ymm-(h/2)+h*bounds[2,"Mutant"],lwd=2,lend=2,col=col.err)
  tree
}


plot_cml_sig_tree_posterior=function(PD,prior,b.ultra=FALSE,txtitle=PD$patient,b.add.sigs=TRUE){
  #browser()
  info=PD$pdx$dat$info
  if(is.null(info)){
    info=get_clade_sigs(PD,prior,extracats = c())
  }
  if(is.null(PD$pdx$dat$sigbynode)){
    PD2=add_cml_sig_to_pdx2(PD,prior=prior,info = info)
  }else{
    PD2=PD
  }
  colorscheme=get_sig_color_scheme()
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  lmp=ifelse(b.add.sigs,0.25,0.35)
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = lmp)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  tree=plot_signature_annotated_tree(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  ## Convert data to form ## For 
  if(!b.add.sigs){
    agg1=info$ctx %>% group_by(category) %>% summarise(total=sum(count))
    df=info$ctx %>% filter(ctx=="C>T at CpG") %>% left_join(agg1)
    bounds=t(sapply(1:dim(df)[1],function(i) binom.test(df$count[i],df$total[i])$conf.int))
    colnames(bounds)=c("lb","ub")
    df=cbind(df %>% mutate(value=count/total),as.data.frame(bounds))
    df=df %>% left_join(data.frame(category=c("Mutant","MutantI","Trunk","WildType"),label=c("CML","Early CML","Pre-CML","WT")))#left_join(data.frame(category=c("Mutant","WildType","Trunk","MutantI"),label=c("Mutant","WT","Trunk","MutantI")))
    df=df %>% filter(category %in% c("Mutant","MutantI","Trunk","WildType"))
    
    df$category=factor(df$category,levels=c("Mutant","MutantI","Trunk","WildType"))
    df=df %>% filter(!is.na(category))
    df=df[order(df$category),]
    v=par("usr")
    ncat=length(df$category)
    ## approx normalis the width for number of categories
    v=c( v[1]+abs(v[1])*0.1, v[1]+abs(v[1])*(0.1+0.9*((ncat+0.5)/(4+0.5))), v[3]+0.15*(v[4]-v[3]), v[3]+0.65*(v[4]-v[3]) )
    #subplot(do_cpg_barplot(df),x=v[1:2],y=v[3:4])
    #browser()
    #add_bar_chart(df %>% mutate(type=category),col2 =c("darkturquoise","white"),txttitle="C>T @ CpG")
    add_bar_chart2(df %>% mutate(type=category),col2 =c("darkturquoise","white"),txttitle="C>T @ CpG")
  }else{
    #df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean) %>% filter(signature=="SBS1") %>% 
    #  left_join(data.frame(category=c("Mutant","WildType","Trunk","MutantI"),label=c("BCR::ABL","WT","Trunk","MTI")))
    #browser()
    #add_bar_chart(df %>% mutate(type=category),colorscheme$COL)
    df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean)  
    df=df %>% filter(category %in% c("Mutant","MutantI","Trunk","WildType"))
    df$category=factor(df$category,levels=c("Mutant","MutantI","Trunk","WildType"))
    
    df=df %>% left_join(data.frame(category=c("Mutant","MutantI","Trunk","WildType"),label=c("CML","Early CML","Pre-CML","WT")))
    df=df[order(df$category),]
    df$category=as.character(df$category)
    #browser()
    DF=sapply(c("mean","lb","ub"),function(field) df %>% 
                #dplyr::select(c("category","signature",field)) %>% mutate(category=ifelse(category=="WildType","WT",category)) %>%
                #pivot_wider(names_from = "category",values_from=field) %>% 
                dplyr::select(c("label","signature",field)) %>%
                pivot_wider(names_from = "label",values_from=field) %>% 
                column_to_rownames("signature") %>% as.matrix(),USE.NAMES = TRUE,simplify = FALSE)
    v=par("usr")
    ncat=dim(DF$mean)[2]
    ## approx normalis the width for number of categories
    v=c( v[1]+abs(v[1])*0.1, v[1]+abs(v[1])*(0.1+0.9*((ncat+0.5)/(4+0.5))), v[3]+0.15*(v[4]-v[3]), v[3]+0.65*(v[4]-v[3]) )
    
    #browser()
    subplot(do_barplot(DF,colorscheme),x=v[1:2],y=v[3:4])
    
  }
  tree
}


do_barplot=function(DF,colorscheme){
  ###browser()
  bp=barplot(DF$mean,col=data.frame(SIG=rownames(DF$mean)) %>% inner_join(colorscheme) %>% pull(COL) ,beside = TRUE,ylim=c(0,1),xaxt="n")
  for(i in 1:dim(DF$lb)[1]){
    segments(x0=bp[i,],y0=as.numeric(DF$lb[i,]),y1 = as.numeric(DF$ub[i,]),lwd=2)
  }
  text(x=apply(bp,2,max),y=-0.05,colnames(DF$mean),srt=60,pos=2,xpd=TRUE)
}

do_cpg_barplot=function(df){
  bp=barplot(df$value,ylim=c(0,1),names.arg = df$label,col="darkorange")
  for(i in 1:dim(bp)[1]){
    segments(x0=bp[i,],y0=df$lb[i],y1 = df$ub[i],lwd=2)
  }
}






plot_cml_ctx_tree=function(PD,b.ultra=FALSE,txtitle=PD$patient,axis.ticks=seq(0,1,0.1),b.include.indel=FALSE){
  
  PD2=add_ctx_to_pdx(PD,b.include.indel = b.include.indel )
  cols=RColorBrewer::brewer.pal(3,"Set1")
  colorscheme=get_ctx_config() %>% mutate(COL=colors,SIG=mut_type)#data.frame(COL=cols[1:2],SIG=c("SBS1","SBSblood"))
  if(!b.include.indel){
    colorscheme=colorscheme %>% filter(mut_type != "INDEL")
  }
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = 0.3)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  tree=plot_signature_annotated_tree(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  info=PD2$pdx$dat$info
  agg1=info$ctx %>% group_by(category) %>% summarise(total=sum(count))
  df=info$ctx %>% filter(ctx=="C>T at CpG") %>% left_join(agg1)
  df=df %>% filter(!is.na(total))
  bounds=t(sapply(1:dim(df)[1],function(i) binom.test(df$count[i],df$total[i])$conf.int))
  colnames(bounds)=c("lb","ub")
  df=cbind(df %>% mutate(value=count/total),as.data.frame(bounds))
  df=df %>% left_join(data.frame(category=c("Mutant","WildType","Trunk","MutantI"),label=c("Mutant","WT","Trunk","MutantI")))
  df=df %>% filter(category %in% c("Mutant","MutantI","Trunk","WildType"))
  df$category=factor(df$category,levels=c("Mutant","MutantI","Trunk","WildType"))
  df=df[order(df$category),]
  #browser()
  add_bar_chart(df %>% mutate(type=category),col2 =c("darkorange","white"),txttitle="C>T @ CpG")
  
  if(FALSE){
    inf=get_sig_info(PD)
    idx=which(rownames(inf$counts)=="C>T at CpG")
    ct_at_cpg=inf$counts[idx,]/colSums(inf$counts)
    bounds=sapply(1:3,function(i) binom.test(inf$counts[idx,i],colSums(inf$counts)[i])$conf.int)
    colnames(bounds)=colnames(inf$counts)
    xm=abs(par("usr")[1])
    ym=par("usr")[c(3,4)]
    ymm=mean(ym)
    h=diff(ym)*0.5
    w=0.2
    gap=0.1
    col.err="black"
    col.cpg="darkorange"
    col.other="white"
    f=ct_at_cpg[1,"WildType"]
    ##browser()
    segments(x0 =-xm+0.09*xm,y0=ymm-(h/2),y1=ymm+(h/2),lwd=1)
    for(x in axis.ticks){
      segments(x0=-xm+0.05*xm,x1=-xm+0.09*xm,y0=ymm-(h/2)+(h*x),lwd=1)
      text(sprintf("%3.2f",x),x=-xm+0.05*xm,y=ymm-(h/2)+(h*x),pos=2,cex=0.8,xpd=NA,offset=0.1)
    }
    
    rect(xleft=-xm+0.1*xm,xright=-xm+(0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
    rect(xleft=-xm+0.1*xm,xright=-xm+(0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
    segments(x0 = -xm+(0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"WildType"],y1 = ymm-(h/2)+h*bounds[2,"WildType"],lwd=2,lend=2,col=col.err)
    text("WT",x=-xm+(0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    
    f=ct_at_cpg[1,"Trunk"]
    rect(xleft=-xm+(0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
    rect(xleft=-xm+(0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
    text("Trunk",x=-xm+(0.1+w+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    segments(x0 = -xm+(0.1+w+0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"Trunk"],y1 = ymm-(h/2)+h*bounds[2,"Trunk"],lwd=2,lend=2,col=col.err)
    text("C>T @ CpG",x=-xm+(0.1+w+0.1+0.5*w)*xm,y=ymm+(h/2),pos=3,cex=1.2)
    
    f=ct_at_cpg[1,"Mutant"]
    rect(xleft=-xm+(0.1+w+0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*f),col=col.cpg)
    rect(xleft=-xm+(0.1+w+0.1+w+0.1)*xm,xright=-xm+(0.1+w+0.1+w+0.1+w)*xm,ybottom =ymm-(h/2)+(h*f),ytop=ymm+(h/2),col=col.other)
    text("BCR::ABL1",x=-xm+(0.1+w+0.1+w+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    segments(x0 = -xm+(0.1+w+0.1+w+0.1+0.5*w)*xm,y0=ymm-(h/2)+h*bounds[1,"Mutant"],y1 = ymm-(h/2)+h*bounds[2,"Mutant"],lwd=2,lend=2,col=col.err)
  }
  tree
}



plot_signature_annotated_tree=function(tree,pdx,sig.col.scheme,maxlen=NULL,b.include.csm.text=TRUE){
  if(is.null(pdx$sigbynode)){
    stop("Need to build sigbynode first: see add_sig_to_pdx")
  }
  ##need df of SIG,col
  control=list(col.scheme=sig.col.scheme)
  if(!is.null(maxlen)){
    control$maxlen=maxlen
    control$b.include.csm.text=b.include.csm.text
  }
  res=add_annotation(pdx,
                     tree,
                     add_signature_decomp,control=control)
  leg=legend("topleft",sig.col.scheme$SIG,col = sig.col.scheme$COL,pch=15,pt.cex=2)$rect
  tree
}


add_bar_chart_ORIG=function(df,col2,col.err="black",txttitle=""){
  ## df has the form value,lb,ub,type,label
  xm=abs(par("usr")[1])
  ym=par("usr")[c(3,4)]
  ymm=mean(ym)
  h=diff(ym)*0.5
  w=0.15
  gap=0.1
  inc=0.05
  #col.err="black"
  col.mut=col2[1]
  col.other=col2[2]
  segments(x0 =-xm+0.09*xm,y0=ymm-(h/2),y1=ymm+(h/2),lwd=1)
  for(x in c(0,0.25,0.5,0.75,1)){
    segments(x0=-xm+0.05*xm,x1=-xm+0.09*xm,y0=ymm-(h/2)+(h*x),lwd=1)
    text(sprintf("%3.2f",x),x=-xm+0.05*xm,y=ymm-(h/2)+(h*x),pos=2,cex=0.8,xpd=NA,offset=0.1)
  }
  xd=0
  for(TYPE in unique(df$type)){
    tmp=df %>% filter(type==TYPE) 
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*tmp$value),col=col.mut)
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom =ymm-(h/2)+(h*tmp$value),ytop=ymm+(h/2),col=col.other)
    segments(x0 = -xm+(0.1+xd+0.5*w)*xm,y0=ymm-(h/2)+h*tmp$lb,
             y1 = ymm-(h/2)+h*tmp$ub,lwd=2,lend=2,col=col.err)
    text(tmp$label,x=-xm+(xd+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    xd=xd+inc+w
  }
  text(txttitle,x=-xm+(0.1*w)*xm,y=ymm+1.1*(h/2),pos=4,cex=1.2)
}

convert_tree_to_hdp=function(PD,node){
  
}

get_private_vaf_vs_signature_DF=function(PD,signame="SBS18"){
  
  details=add_trin(PD$pdx$dat$details %>% mutate(IDX=1:length(Chrom)))
  node=PD$nodes %>% filter(driver=="BCR::ABL1") %>% pull(node)
  if(length(node)==0){
    return(NULL)
  }
  nodes=intersect(get_all_node_children(node,PD$pdx$tree_ml),1:length(PD$pdx$tree_ml$tip.label))
  d1=details %>% filter(TYPE=="SNV" & is_localx_excluded==0 & node %in% nodes)
  d1$P_SIG=PD$pdx$dat$info$MP$Mutant[match(d1$trin,rownames(PD$pdx$dat$info$MP$Mutant)),signame]
  d1=d1 %>% filter(node<=length(PD$pdx$tree_ml$tip.label))
  d1$sample=PD$pdx$tree_ml$tip.label[d1$node]
  vaf=PD$pdx$dat$mtr/PD$pdx$dat$dep
  VAF=vaf %>% as.data.frame() %>% mutate(IDX=1:dim(vaf)[1]) %>% pivot_longer(cols=colnames(vaf),names_to = "sample",values_to = "VAF")
  D1=d1 %>% left_join(VAF)
  D1 %>% mutate(SIG=signame,donor=PD$patient) %>% dplyr::select(Chrom,Pos,P_SIG,node,SIG,VAF,donor)
}


get_sig_info_counts=function(PD){
  
  vals=sapply(
    c("Mutant","WildType","Trunk","MutantI"),
    function(categ) {
      details=get_details_for_category(PD,categ)
      if(dim(details)[1]==0){
        return(NULL)
      }else{
        mcm=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
        mto=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
        return(list(mcm=mcm,mto=mto))
      }
    },USE.NAMES = TRUE,simplify = FALSE
  )
  vals=vals[!sapply(vals,is.null)]
  if(length(vals)>1){
    mcc=do.call("cbind",lapply(vals,function(x) x$mcm)) %>% (function(x){colnames(x)=names(vals);x})
    alt=do.call("cbind",lapply(vals,function(x) t(x$mto))) %>% (function(x){colnames(x)=names(vals);x})
  }else{
    mcc=vals[[1]]$mcm %>% (function(x){colnames(x)=names(vals);x})
    alt=vals[[1]]$mto %>% (function(x){colnames(x)=names(vals);x})
  }
  mcc2=mcc %>% as.data.frame() %>% tibble:::rownames_to_column(var="trin")
  mcc2=mcc2 %>% mutate(ctx=substr(trin,3,5)) %>% 
    mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
    mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx))
  mcc3= mcc2 %>% pivot_longer(cols=names(vals),names_to = "group",values_to = "count")
  counts=mcc3 %>% group_by(group,ctx) %>% summarise(N=n(),count=sum(count)) %>% pivot_wider(names_from = "group",values_from = "count") %>% 
    tibble::column_to_rownames(var="ctx") %>% dplyr::select(names(vals))
  ##print(p1)
  list(counts=counts,alt=alt)
}


get_clade_sigs=function(PD,prior,sigmethod="SIGFIT",corecats=c("Mutant","WildType","Trunk"),extracats=c("MutantI","Top50","Top100"),min.muts=100){
  stats=NULL
  CTX=NULL
  mp=list()
  for(categ in c(corecats,extracats)){
    cat("Processing",categ,"..\n")
    details=get_details_for_category(PD,categ)
    if(dim(details)[1]<min.muts){
      next
    }
    mutc=with(details,get_mut_matrix(Chrom,Pos,Ref,Alt))
    ctx=as.data.frame(mutc) %>% mutate(count=gr,trin=rownames(.)) %>% mutate(ctx=substr(trin,3,5)) %>% 
      mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
      mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx)) %>% group_by(ctx) %>% summarise(N=n(),count=sum(count)) %>%
      mutate(donor=PD$patient,category=categ)
    CTX=rbind(CTX,ctx)
    sig=matrix(as.vector(mutc),ncol=1)
    if(sigmethod=="SIGFIT"){
      sf=fit_signatures(t(sig),t(prior),iter=4000,refresh=0,control=list(adapt_delta=0.95))
      exposures=retrieve_pars(sf,par="exposures")
      decomp=list(contribution=t(as.matrix(exposures$mean)))
      colnames(decomp$contribution)=categ
      # retrieve distribution csm.
      zz=rstan::extract(sf$result)
      exposures$sd=apply(zz$exposures[,1,],2,sd)
      names(exposures$sd)=names(exposures$mean)
      tmp=zz$exposures[,1,] %*% t(prior)
      CSM=sapply(1:8000,function(i) cos_sim(tmp[i,],as.vector(sig)))
      bounds=quantile(CSM,prob=c(0.025,0.975))
      decomp$csm=mean(CSM) ##cos_sim(as.vector(retrieve_pars(sf,par="reconstructions")$mean),as.vector(sig))
      decomp$csm_lb=bounds[1]
      decomp$csm_ub=bounds[2]
      ## Following is a bit ugly - but emulates the structure emitted by MP method below.
      qq=rbind(do.call("rbind",exposures) %>% (function(x){rownames(x)=c("mean","2.5%","97.5%","sd");x}),
               total=rep(sum(mutc),dim(exposures$mean)[2]))
    }else{
      ## Mutational patterns fit
      decomp=MutationalPatterns:::fit_to_signatures(sig,prior)
      decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
      decomp$csm=cos_sim(decomp$reconstructed[,1],as.vector(mutc))
      colnames(decomp$contribution)=categ
      ## Next do bootstrapping to get CIs.
      bs=fit_to_signatures_bootstrapped(sig,prior,verbose = FALSE,method = "regular") %>% (function(probs) sweep(probs,MARGIN=1,rowSums(probs),"/"))
      qq=sapply(colnames(bs),function(x) {qq=quantile(bs[,x],c(0.025,0.5,0.975));c(qq,mean=mean(bs[,x],),sd=sd(bs[,x]),total=sum(mutc))})
      
    }
    stats=rbind(stats,as.data.frame(qq) %>% mutate(field=rownames(.),category=categ,donor=PD$patient,rec_csm=decomp$csm,csm_lb=decomp$csm_lb,csm_ub=decomp$csm_ub) )
    MP=as.matrix(get_mutation_probability(prior,t(decomp$contribution))[,-(1:2)])
    mp[[categ]]=MP
    cat("\n")
  }
  list(corecats=corecats,extracats=extracats,sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
         pivot_wider(names_from="field",values_from = "proportion"),
       ctx=CTX,MP=mp)
  
}


add_bar_chart2=function(df,col2,col.err="black",txttitle=""){
  ## df has the form value,lb,ub,type,label
  xm=abs(par("usr")[1])
  ym=par("usr")[c(3,4)]
  ymm=mean(ym)
  h=diff(ym)*0.5
  w=0.15
  gap=0.1
  inc=0.05
  #col.err="black"
  col.mut=col2[1]
  col.other=col2[2]
  segments(x0 =-xm+0.09*xm,y0=ymm-(h/2),y1=ymm+(h/2),lwd=1)
  for(x in seq(0,1,0.1)){
    segments(x0=-xm+0.05*xm,x1=-xm+0.09*xm,y0=ymm-(h/2)+(h*x),lwd=1)
    text(sprintf("%3.2f",x),x=-xm+0.05*xm,y=ymm-(h/2)+(h*x),pos=2,cex=0.8,xpd=NA,offset=0.1)
  }
  xd=0
  for(TYPE in unique(df$type)){
    tmp=df %>% filter(type==TYPE) 
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*tmp$value),col=col.mut,border=NA)
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom =ymm-(h/2)+(h*tmp$value),ytop=ymm+(h/2),col=col.other,border=NA)
    segments(x0 = -xm+(0.1+xd+0.5*w)*xm,y0=ymm-(h/2)+h*tmp$lb,
             y1 = ymm-(h/2)+h*tmp$ub,lwd=2,lend=2,col=col.err)
    text(tmp$label,x=-xm+(xd+0.1+w)*xm,y=ymm-(h/2),cex=1,srt=60,pos=2,xpd=TRUE)
    xd=xd+inc+w
  }
  text(txttitle,x=-xm+(0.1*w)*xm,y=ymm+1.1*(h/2),pos=4,cex=1.2)
}

