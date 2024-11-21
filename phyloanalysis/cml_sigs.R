require(RVAideMemoire)
require(MutationalPatterns)
require(sigfit)

get_cml_details=function(PD,tree.cat){
  ## 
  tnode=PD$nodes %>% filter(driver=="BCR::ABL1") %>% pull(node)
  if(length(tnode)==0){
    if(tree.cat=="WildType"){
      return(PD$pdx$dat$details %>% filter(TYPE=="SNV" & is_localx_excluded==0))
    }else{
      cat("no mutant node defined - returning empty matrix\n")
    }
    return(PD$pdx$dat$details %>% head(0))
  }
  if(!(tree.cat %in% c("Mutant","WildType","Trunk","MutantI"))){
    stop("invalid tree category")
  }
  if(tree.cat=="Trunk"){
    PD$pdx$dat$details %>% filter(node==tnode & TYPE=="SNV" & is_localx_excluded==0)
  }else if(tree.cat=="WildType"){
    PD$pdx$dat$details %>% filter((!(node %in% c(tnode,get_all_node_children(tnode,PD$pdx$tree_ml))) & TYPE=="SNV" & is_localx_excluded==0))
  }else if(tree.cat=="Mutant"){
    PD$pdx$dat$details %>% filter((node %in% get_all_node_children(tnode,PD$pdx$tree_ml) )& TYPE=="SNV" & is_localx_excluded==0 )
  }else if(tree.cat=="MutantI"){
    PD$pdx$dat$details %>% filter((node %in% setdiff(get_all_node_children(tnode,PD$pdx$tree_ml),1:length(PD$pdx$tree_ml$tip.label) ) ) & TYPE=="SNV" & is_localx_excluded==0)
  }
}




get_sig_info=function(PD,node=PD$nodes %>% filter(driver=="BCR::ABL1") %>% (function(x)(x$node))){
  tnode=node
  details=PD$pdx$dat$details %>% filter((node %in% get_all_node_children(tnode,PD$pdx$tree_ml) )& TYPE=="SNV")
  mcm=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcmg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  
  details=PD$pdx$dat$details %>% filter((node %in% setdiff(get_all_node_children(tnode,PD$pdx$tree_ml),1:length(PD$pdx$tree_ml$tip.label) ) ) & TYPE=="SNV")
  mci=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcmi=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  
  
  details=PD$pdx$dat$details %>% filter((!(node %in% c(tnode,get_all_node_children(tnode,PD$pdx$tree_ml))) & TYPE=="SNV"))
  mcw=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcwg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  details=PD$pdx$dat$details %>% filter(node==tnode & TYPE=="SNV")
  mct=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mctg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  mcc=cbind(mcm,mcw,mct,mci)
  colnames(mcc)=c("Mutant","WildType","Trunk","MutantI")
  #plot_96_profile(mcc)
  p1=MutationalPatterns::plot_spectrum(mcwg,error_bars = "none",CT=TRUE,legend = FALSE)+ggtitle("WildType")+ylim(c(0,0.6))
  p2=MutationalPatterns::plot_spectrum(mctg,error_bars = "none",CT=TRUE,legend = FALSE)+ggtitle("Trunk")+ylim(c(0,0.6))
  p3=MutationalPatterns::plot_spectrum(mcmg,error_bars = "none",CT=TRUE,legend=TRUE)+ggtitle("Mutant")+ylim(c(0,0.6))
  p4=MutationalPatterns::plot_spectrum(mcmi,error_bars = "none",CT=TRUE,legend=TRUE)+ggtitle("MutantI")+ylim(c(0,0.6))
  
  #zz=grid.arrange(p1,p2,p3,ncol=3,widths=c(1.75,1.75,3),top=textGrob(PD$patient))
  mcc2=mcc %>% as.data.frame() %>% tibble:::rownames_to_column(var="trin")
  mcc2=mcc2 %>% mutate(ctx=substr(trin,3,5)) %>% 
    mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
    mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx))
  mcc3= mcc2 %>% pivot_longer(cols=c("Mutant","WildType","Trunk","MutantI"),names_to = "group",values_to = "count")
  counts=mcc3 %>% group_by(group,ctx) %>% summarise(N=n(),count=sum(count)) %>% pivot_wider(names_from = "group",values_from = "count") %>% 
    tibble::column_to_rownames(var="ctx") %>% dplyr::select(Mutant,Trunk,WildType,MutantI)
  ## tests - include all contexts..
  c1=chisq.test(t(as.matrix(counts)))
  g1=G.test(t(as.matrix(counts)))
  gp=pairwise.G.test(t(as.matrix(counts)))
  print(c1)
  print(g1)
  ##print(p1)
  list(c1=c1,g1=g1,gp=gp,plots=list(wt=p1,trunk=p2,mutant=p3,mutanti=p4),mcc=mcc,counts=counts)
}

get_sig_info2=function(PD,node=PD$nodes %>% filter(driver=="BCR::ABL1") %>% (function(x)(x$node)),extra.node=-1){
  tnode=node
  details=PD$pdx$dat$details %>% filter((node %in% get_all_node_children(tnode,PD$pdx$tree_ml) )& TYPE=="SNV")
  mcm=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcmg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  ##
  ##if(extra.node>0){
  browser()
  details=PD$pdx$dat$details %>% filter((node %in% get_all_node_children(extra.node,PD$pdx$tree_ml) )& TYPE=="SNV")
  
  mcmx=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcmgx=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  details=PD$pdx$dat$details %>% filter((node %in% setdiff(get_all_node_children(tnode,PD$pdx$tree_ml),get_all_node_children(extra.node,PD$pdx$tree_ml) )) & TYPE=="SNV")
  mcmy=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcmgy=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  
  
  details=PD$pdx$dat$details %>% filter((!(node %in% c(tnode,get_all_node_children(tnode,PD$pdx$tree_ml))) & TYPE=="SNV"))
  
  
  
  
  
  mcw=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mcwg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  details=PD$pdx$dat$details %>% filter(node==tnode & TYPE=="SNV")
  mct=get_mut_matrix(details$Chrom,details$Pos,details$Ref,details$Alt)
  mctg=mut_type_occurrences(get_grange_df(details),ref_genome = ref_genome)
  
  
  
  
  mcc=cbind(mcm,mcw,mct,mcmx,mcmy)
  colnames(mcc)=c("Mutant","WildType","Trunk","MutantDiagnosis","MutantRemission")
  #plot_96_profile(mcc)
  p1=MutationalPatterns::plot_spectrum(mcwg,error_bars = "none",CT=TRUE,legend = FALSE)+ggtitle("WildType")+ylim(c(0,0.6))
  p2=MutationalPatterns::plot_spectrum(mctg,error_bars = "none",CT=TRUE,legend = FALSE)+ggtitle("Trunk")+ylim(c(0,0.6))
  p3=MutationalPatterns::plot_spectrum(mcmg,error_bars = "none",CT=TRUE,legend=TRUE)+ggtitle("Mutant")+ylim(c(0,0.6))
  p4=MutationalPatterns::plot_spectrum(mcmgx,error_bars = "none",CT=TRUE,legend=TRUE)+ggtitle("Mutant: Diagnosis Timepoint")+ylim(c(0,0.6))
  p5=MutationalPatterns::plot_spectrum(mcmgy,error_bars = "none",CT=TRUE,legend=TRUE)+ggtitle("Mutant: Remission Timepoint")+ylim(c(0,0.6))
  
  #zz=grid.arrange(p1,p2,p3,ncol=3,widths=c(1.75,1.75,3),top=textGrob(PD$patient))
  mcc2=mcc %>% as.data.frame() %>% tibble:::rownames_to_column(var="trin")
  mcc2=mcc2 %>% mutate(ctx=substr(trin,3,5)) %>% 
    mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
    mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx))
  mcc3= mcc2 %>% pivot_longer(cols=c("Mutant","WildType","Trunk","MutantDiagnosis","MutantRemission"),names_to = "group",values_to = "count")
  counts=mcc3 %>% group_by(group,ctx) %>% summarise(N=n(),count=sum(count)) %>% pivot_wider(names_from = "group",values_from = "count") %>% tibble::column_to_rownames(var="ctx") %>% dplyr::select(Mutant,Trunk,WildType,MutantDiagnosis,MutantRemission)
  ## tests - include all contexts..
  c1=chisq.test(t(as.matrix(counts)))
  g1=G.test(t(as.matrix(counts)))
  gp=pairwise.G.test(t(as.matrix(counts)))
  print(c1)
  print(g1)
  ##print(p1)
  list(c1=c1,g1=g1,gp=gp,plots=list(wt=p1,trunk=p2,mutant=p3,mutantx=p4,mutanty=p5),mcc=mcc,counts=counts)
}
add_sig_to_pdx_cml=function(PD,prior,mutprob=NULL,minmut=50){
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
    ##We have Trunk,Mutant, WT
    sigs=colnames(mutprob)[-c(1:2)]
    prior=prior[,sigs]
    tree=PD$pdx$tree_ml
    tree$color=rep("black",length(tree$edge.length))
    ## child Node of transition assumed to happen at the end of this branch.
    node=PD$nodes %>% filter(driver=="BCR::ABL1") %>% (function(x)(x$node))
    mt.samples=get_samples_in_clade(node,tree)
    #wt.samples=setdiff(tree$tip.label,c(mt.samples,"zeros"))
    #wt.nodes=match(wt.samples,tree$tip.label)
    mt.nodes=get_all_node_children(node,tree)
    wt.nodes=setdiff(tree$edge[,2],c(mt.nodes,node))
    trunk.nodes=node
    wt.idx=match(wt.nodes,tree$edge[,2])
    mt.idx=match(mt.nodes,tree$edge[,2])
    trunk.idx=match(trunk.nodes,tree$edge[,2])
    sigbynode=list()
    tmp=list(WT=list(label="WildType",idx=wt.idx,color="black"),
             TRUNK=list(label="Trunk",idx=trunk.idx,color="blue"),
             MT=list(label="Mutant",idx=mt.idx,color="red"))
    for(x in tmp){
      tree$color[x$idx]=x$color
      pidx=which(mutprob$sample==x$label)
      cat("running",x$label,"\n")
      for(i in x$idx){
        cat("#")
        node=PD$pdx$tree_ml$edge[i,2]
        mutc=get_mut_matrix_for_node(PD$pdx,node)
        if(!is.null(mutc)){
          
          contr=t(t(mutc) %*% as.matrix(mutprob[pidx,sigs]))
          contr=contr/sum(contr)
          csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
          sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
        }else{
          sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
        }
      }
    }
    cat("\n")
    PD$pdx$dat$sigbynode=sigbynode
    PD$pdx$dat$siglist=sigs
    PD$test_tree=tree
    return(PD)
  }
  
}


do_sigs=function(PD,sigs){
  siginf=get_sig_info(PD)
  mcmcfit = sigfit::fit_signatures(
    counts = t(siginf$mcc),
    signatures = t(sigs),
    iter = 10000,warmup = 5000,chains = 1,seed = 1234567
    )
  
  zzz=retrieve_pars(mcmcfit,par="activities",hpd_prob=0.95)
  mutprob=get_mutation_probability(sigs,contr = as.matrix(zzz$mean))
  PD=add_sig_to_pdx_cml(PD,prior = sigs,mutprob = mutprob)
  list(PD=PD,siginf=siginf,mcmcfit=mcmcfit)
}


get_ctx_counts=function(PD){
  if(!is.null(PD$pdx$dat$info$ctx)){
    inf=list(counts=PD$pdx$dat$info$ctx %>% dplyr::select(ctx,count,category) %>% 
               pivot_wider(names_from = category,values_from=count) %>% column_to_rownames("ctx"))
  }else{
    inf=get_sig_info_counts(PD) 
  }
  idx=which(rownames(inf$counts)=="C>T at CpG")
  ct_at_cpg=inf$counts[idx,]/colSums(inf$counts)
  bounds=sapply(1:ncol(inf$counts),function(i) binom.test(inf$counts[idx,i],colSums(inf$counts)[i])$conf.int)
  colnames(bounds)=colnames(inf$counts)
  bounds=as.data.frame(t(bounds) %>% (function(x){colnames(x)=c("lb","ub");x})) %>% mutate(donor=PD$patient,type=rownames(.))
  bounds$count=as.numeric(inf$counts["C>T at CpG",bounds$type])
  bounds$total=colSums(inf$counts[,bounds$type,drop=FALSE])
  bounds$prop=bounds$count/bounds$total
  bounds
}


get_clade_sigs_ORIG=function(PD,prior){
  stats=NULL
  CTX=NULL
  mp=list()
  #browser()
  for(categ in c("WildType","Trunk","Mutant","MutantI")){
    cat("Processing",categ,"..\n")
    details=get_cml_details(PD,categ)
    browser()
    if(dim(details)[1]==0){
      next
    }
    mutc=with(details,get_mut_matrix(Chrom,Pos,Ref,Alt))
    ctx=as.data.frame(mutc) %>% mutate(count=gr,trin=rownames(.)) %>% mutate(ctx=substr(trin,3,5)) %>% 
      mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
      mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx)) %>% group_by(ctx) %>% summarise(N=n(),count=sum(count)) %>%
      mutate(donor=PD$patient,category=categ)
    CTX=rbind(CTX,ctx)
    ## Do the 
    #qf=quickfit(as.vector(mutc),prior)
    sig=matrix(as.vector(mutc),ncol=1)
    ## Mutational patterns fit
    decomp=MutationalPatterns:::fit_to_signatures(sig,prior)
    decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
    decomp$csm=cos_sim(decomp$reconstructed[,1],as.vector(mutc))
    colnames(decomp$contribution)=categ
    ## Next do bootstrapping to get CIs.
    bs=fit_to_signatures_bootstrapped(sig,prior,verbose = FALSE,method="regular") %>% (function(probs) sweep(probs,MARGIN=1,rowSums(probs),"/"))
    
    #quantile(bs[,1]/rowSums(bs))
    qq=sapply(colnames(bs),function(x) {qq=quantile(bs[,x],c(0.025,0.5,0.975));c(qq,mean=mean(bs[,x],),sd=sd(bs[,x]),total=sum(mutc))})
    stats=rbind(stats,as.data.frame(qq) %>% mutate(field=rownames(.),category=categ,donor=PD$patient,rec_csm=decomp$csm) )
    MP=as.matrix(get_mutation_probability(prior,t(decomp$contribution))[,-(1:2)])
    mp[[categ]]=MP
    cat("\n")
  }
  list(sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
    pivot_wider(names_from="field",values_from = "proportion"),
    ctx=CTX,MP=mp)
  
}



## Adds grouped category level signature to the tree 
add_cml_sig_to_pdx=function(PD,prior){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  stats=NULL
  for(categ in c("WildType","Trunk","Mutant")){
    cat("Processing",categ,"..\n")
    details=get_cml_details(PD,categ)
    mutc=with(details,get_mut_matrix(Chrom,Pos,Ref,Alt))
  
    ## Do the 
    #qf=quickfit(as.vector(mutc),prior)
    sig=matrix(as.vector(mutc),ncol=1)
    ## Mutational patterns fit
    decomp=MutationalPatterns:::fit_to_signatures(sig,prior)
    decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
    decomp$csm=cos_sim(decomp$reconstructed[,1],as.vector(mutc))
    colnames(decomp$contribution)=categ
    ## Next do bootstrapping to get CIs.
    bs=fit_to_signatures_bootstrapped(sig,prior,verbose = FALSE) %>% (function(probs) sweep(probs,MARGIN=1,rowSums(probs),"/"))
    
    #quantile(bs[,1]/rowSums(bs))
    qq=sapply(colnames(bs),function(x) {qq=quantile(bs[,x],c(0.025,0.5,0.975));c(qq,mean=mean(bs[,x],),sd=sd(bs[,x]))})
    #browser()
    stats=rbind(stats,as.data.frame(qq) %>% mutate(field=rownames(.),category=categ,donor=PD$patient) )
    MP=as.matrix(get_mutation_probability(prior,t(decomp$contribution))[,-(1:2)])
    for(node in unique(details$node)){
      i=which(PD$pdx$tree_ml$edge[,2]==node)
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% MP)
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      cat("#")
    }
    cat("\n")
  }
  ## Fill in any empty nodes
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
    if(i>length(sigbynode) || is.null(sigbynode[[i]])){
      sigbynode[[i]]=list(node=PD$pdx$tree_ml$edge[i,2],contr=empty,csm=NA,count=0)
    }
  }
  #browser()
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  PD$pdx$dat$sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
    pivot_wider(names_from="field",values_from = "proportion")
  return(PD)
}


## Adds grouped category level signature to the tree 
add_cml_sig_to_pdx2=function(PD,prior,info=get_clade_sigs(PD,prior)){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  #stats=NULL
  for(categ in c("WildType","Trunk","Mutant")){
    cat("Processing",categ,"..\n")
    details=get_cml_details(PD,categ)
    MP=info$MP[[categ]]
    for(node in unique(details$node)){
      i=which(PD$pdx$tree_ml$edge[,2]==node)
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% MP)
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      cat("#")
    }
    cat("\n")
  }
  ## Fill in any empty nodes
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
    if(i>length(sigbynode) || is.null(sigbynode[[i]])){
      sigbynode[[i]]=list(node=PD$pdx$tree_ml$edge[i,2],contr=empty,csm=NA,count=0)
    }
  }
  #browser()
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  #PD$pdx$dat$sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
  #  pivot_wider(names_from="field",values_from = "proportion")
  return(PD)
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
      sf=fit_signatures(t(sig),t(prior),iter=4000,refresh=0,control=list(adapt_delta=0.995))
      exposures=retrieve_pars(sf,par="exposures")
      decomp=list(contribution=t(as.matrix(exposures$mean)))
      colnames(decomp$contribution)=categ
      # retrieve distribution csm.
      zz=rstan::extract(sf$result)
      tmp=zz$exposures[,1,] %*% t(prior)
      CSM=sapply(1:8000,function(i) cos_sim(tmp[i,],as.vector(sig)))
      bounds=quantile(CSM,prob=c(0.025,0.975))
      decomp$csm=mean(CSM) ##cos_sim(as.vector(retrieve_pars(sf,par="reconstructions")$mean),as.vector(sig))
      decomp$csm_lb=bounds[1]
      decomp$csm_ub=bounds[2]
      ## Following is a bit ugly - but emulates the structure emitted by MP method below.
      qq=rbind(do.call("rbind",exposures) %>% (function(x){rownames(x)=c("mean","2.5%","97.5%");x}),
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

plot_sig_tree_posterior=function(PD,prior,b.ultra=FALSE,txtitle=PD$patient,b.add.sigs=TRUE,focal.signature="SBS1"){
  info=PD$pdx$dat$info
  if(is.null(info)){
    stop('prior to running please do: PD$pdx$dat$info=get_clade_sigs(PD,prior=pcawg[, c("SBS1","SBSblood","SBS18")])')
    ###info=get_clade_sigs_mouse(PD,prior)
  }
  if(is.null(PD$pdx$dat$sigbynode)){
    PD2=add_v2_sig_to_pdx(PD,prior=prior,info = info)
  }else{
    PD2=PD
  }
  colorscheme=get_sig_color_scheme() %>% filter(SIG %in% colnames(prior))
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = 0.3)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  tree=plot_signature_annotated_tree_v2(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  ## Convert data to form ## For 
  if(!b.add.sigs){
    agg1=info$ctx %>% group_by(category) %>% summarise(total=sum(count))
    df=info$ctx %>% filter(ctx=="C>T at CpG") %>% left_join(agg1)
    bounds=t(sapply(1:dim(df)[1],function(i) binom.test(df$count[i],df$total[i])$conf.int))
    colnames(bounds)=c("lb","ub")
    df=cbind(df %>% mutate(value=count/total),as.data.frame(bounds))
    df=df %>% left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
    add_bar_chart(df %>% mutate(type=category),col2 =c("darkorange","white"),txttitle="C>T @ CpG")
  }else{
    df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean) %>% filter(signature==focal.signature) %>% mutate(label=category)
    #  left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
    
    add_bar_chart(df %>% mutate(type=category),c(colorscheme %>% filter(SIG==focal.signature) %>% pull(COL),"white"))
  }
  tree
}

## Adds grouped category level signature to the tree 
add_v2_sig_to_pdx=function(PD,prior,info){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  for(categ in info$corecats){
    cat("Processing",categ,"..\n")
    details=get_details_for_category(PD,categ)
    MP=info$MP[[categ]]
    for(node in unique(details$node)){
      i=which(PD$pdx$tree_ml$edge[,2]==node)
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% MP)
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      cat("#")
    }
    cat("\n")
  }
  ## Fill in any empty nodes
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
    if(i>length(sigbynode) || is.null(sigbynode[[i]])){
      sigbynode[[i]]=list(node=PD$pdx$tree_ml$edge[i,2],contr=empty,csm=NA,count=0)
    }
  }
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  return(PD)
}



get_details_for_category=function(PD,categ,TYPE="SNV"){
  ## 
  tnode=PD$nodes %>% filter(driver=="BCR::ABL1") %>% pull(node)
  if(categ %in% c("Mutant","Trunk","MutantI") && length(tnode)==0){
    cat("no mutant node defined - returning empty matrix\n")
    return(PD$pdx$dat$details %>% head(0))
  }
  ntip=length(PD$pdx$tree_ml$tip.label)
  details=PD$pdx$dat$details %>% filter(is_localx_excluded==0)
  if(TYPE=="SNV"){
    details=details %>% filter(TYPE=="SNV")
  }
  if(categ=="All"){
    details %>% filter(TYPE=="SNV")
  }else if(categ=="Trunk"){
    details %>% filter(node==tnode)
  }else if(categ=="WildType"){
    details %>% filter((!(node %in% c(tnode,get_all_node_children(tnode,PD$pdx$tree_ml))) ))
  }else if(categ=="Mutant"){
    details %>% filter((node %in% get_all_node_children(tnode,PD$pdx$tree_ml) ) )
  }else if(categ=="MutantI"){
    details %>% filter((node %in% setdiff(get_all_node_children(tnode,PD$pdx$tree_ml),1:length(PD$pdx$tree_ml$tip.label) ) ) )
  }else if(categ=="Shared"){
    details %>% filter(node > ntip )
  }else if(categ=="Private"){
    details %>% filter(node <= ntip )
  }else if(grepl("^Top",categ)){
    threshold=as.numeric(gsub("Top","",categ))
    nh=nodeHeights(PD$pdx$tree_ml)
    idx=which(nh[,2]<=threshold)
    nodes=unique(PD$pdx$tree_ml$edge[idx,2])
    details %>% filter(node %in% nodes )
  }else{
    stop("invalid tree category")
  }
}




get_sig_color_scheme=function(){
  ## Save this to a file
  if(file.exists("../data/signature_scheme.txt")){
    scheme=read.table("../data/signature_scheme.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t",comment.char = "")
    if(any(is.na(match(c("COL","SIG"),names(scheme))))){
      stop("Invalid signature scheme - need COL and SIG columns")
    }
  }else{
    cols=RColorBrewer::brewer.pal(8,"Set1")
    scheme=data.frame(SIG=c("SBS1","SBSblood","SBS18"),COL=c("#008080","#DEB887","darkred"))
  }
  scheme
  
}

add_bar_chart=function(df,col2,col.err="black",txttitle=""){
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
    text(tmp$label,x=-xm+(xd+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    xd=xd+inc+w
  }
  text(txttitle,x=-xm+(0.1*w)*xm,y=ymm+1.1*(h/2),pos=4,cex=1.2)
}

plot_signature_annotated_tree_v2=function(tree,pdx,sig.col.scheme,maxlen=NULL,b.include.csm.text=TRUE){
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

get_mutation_probability=function(signatures,### Each column is a signature, each row a context
                                  contr ## Each column is a signature, row as sample
){
  signames=colnames(signatures)
  
  ##Check the mutation context ordering
  idx=match(MutationalPatterns:::TRIPLETS_96,rownames(signatures))
  if(length(which(is.na(idx)))>0){
    stop("ERROR:Unexpected rownames")
  }
  signatures=signatures[idx,]
  mutcontext=rownames(signatures)
  contr=sweep(contr,MARGIN = 1,rowSums(contr),"/")## Make sure contr normalised
  signatures==sweep(signatures,MARGIN = 2,colSums(signatures),"/")
  #samples=rownames(contr)
  ###browser()
  do.call("rbind",
          lapply(rownames(contr),function(x) {
            
            probs=signatures %*% diag(contr[x,])
            probs=sweep(probs,MARGIN=1,rowSums(probs),"/")
            colnames(probs)=colnames(signatures)
            cbind(data.frame(
              sample=rep(x,dim(signatures)[1]),
              mutcontext=mutcontext),
              as.data.frame(probs))})
  )
  
}

get_sig_burdens=function(PD,priors=pcawg[,c("SBS1","SBSblood")],ultramodel="nullmodel",corecats=c("Shared","Private")){
  if(is.null(PD$pdx$dat$sigbynode)){
    PD$pdx$dat$info=get_clade_sigs(PD,prior=priors,extracats = c(),corecats = corecats)
    PD=add_v2_sig_to_pdx(PD,prior=priors,info = PD$pdx$dat$info)
  }
  ## Get burden as tree tip height (tree is adjusted and corresponds to SNV burden)
  nh=nodeHeights(PD$pdx$tree_ml)
  tmp=data.frame(tip.label=PD$pdx$tree_ml$tip.label,
                 nsub_adj=nh[match(1:length(PD$pdx$tree_ml$tip.label),PD$pdx$tree_ml$edge[,2]),2]
  )
  nht=nodeHeights(PD$fit$poisson_tree[[ultramodel]]$ultratree)
  traj=data.frame(start=nht[,1],mid=0.5*(nht[,1]+nht[,2]),end=nht[,2],nsub_adj_start=nh[,1],nsub_adj_mid=0.5*(nh[,1]+nh[,2]),nsub_adj_end=nh[,2],
                  parent=PD$pdx$tree_ml$edge[,1],child=PD$pdx$tree_ml$edge[,2])
  for(SIG in colnames(priors)){
    frac=sapply(1:length(PD$pdx$tree_ml$edge.length),function(i) if(i<=length(PD$pdx$dat$sigbynode)){PD$pdx$dat$sigbynode[[i]]$contr[SIG,]}else{0})
    frac=ifelse(is.na(frac),0,frac)
    tree=PD$pdx$tree_ml
    tree$edge.length=tree$edge.length*frac
    nh=nodeHeights(tree)
    ht=nh[match(1:length(tree$tip.label),tree$edge[,2]),2]
    #tmp2=data.frame(tip.label=tree$tip.label)
    tmp[[paste0("nsub_adj_",SIG)]]=ht
    traj[[paste0(SIG,"_start")]]=nh[,1]
    traj[[paste0(SIG,"_mid")]]=0.5*(nh[,1]+nh[,2])
    traj[[paste0(SIG,"_end")]]=nh[,2]
  }
  
  df=PD$pdx$agedf %>% filter(tip.label!="zeros")  %>% left_join(tmp) %>%
    left_join(PD$pdx$cfg %>% mutate(tip.label=SHORT_LABEL,colony=LABEL) %>% dplyr::select(tip.label,colony)) %>% mutate(donor=PD$patient)
  list(traj=traj,df=df)
}


### DNDS 
do_category_dnds=function(PDD,categ,gene_list=NULL){
  muts=do.call("rbind",lapply(PDD,function(PD){
    muts=get_details_for_category(PD,categ) %>% 
      dplyr::rename(chr=Chrom,pos=Pos,ref=Ref,alt=Alt) %>% 
      mutate(sampleID=sprintf("%s:%s",PD$patient,node)) %>% 
      dplyr::select(chr,pos,ref,alt,sampleID)
  }))
  if(SPECIES_BUILD=="hg38"){
    load("/nfs/casm/team273jn/nw14/projects/pipe2.0/cml_hg38v2/data/covariates_hg19_hg38_epigenome_pcawg.rda")
    load("/nfs/casm/team273jn/nw14/projects/pipe2.0/cml_hg38v2/data/RefCDS_human_GRCh38_GencodeV18_recommended.rda")
    gr_genes<<-gr_genes
    MUT=muts %>% (function(x) x[,c("sampleID","chr","pos","ref","alt")])
    list(DNDS=dndscv(MUT,refdb = RefCDS,cv=covs,gene_list=gene_list),muts=MUT)
  }else{
    dndscv(muts %>% (function(x) x[,c("sampleID","chr","pos","ref","alt")]),gene_list=gene_list)
  }
}
### Function specifying parameters of phylofit.
get_sel_params=function(PDD){
  patients=read.table("../data/patients.txt",head=TRUE,sep="\t") %>% mutate(donor=PDID)
  ## Use the same params as used for phylofit
  params=data.frame(donor=c("PD57334","PD56961","PD57335","PD51634","PD51633","PD51632","PD51635","PD57332"),
                    node=c(75 ,244,71 ,131 ,98 ,89  ,185 ,38),
                    maxS=rep(exp(30)-1,8),
                    adapt_delta=c(0.99,0.999,0.99,0.99,0.999,0.99,0.99,0.99995),
                    max_treedepth=c(14,14,14,14,16,14,14,16)
  )
  params=params %>% left_join(patients %>% dplyr::select(donor,DIAGNOSIS_GAP))
  params$maxt=sapply(PDD[params$donor],function(x) min(x$pdx$agedf$age_at_sample_pcy[x$pdx$agedf$age_at_sample_exact>1],na.rm=TRUE))
  params$maxt=params$maxt-params$DIAGNOSIS_GAP
  params$model="altmodel"
  params
}

##stop("done for now")
get_all_results=function(out){
  PR=get_phylofit_summary(out)
  msg=sprintf("ndiv=%d:rhatmax=%4.3g",PR$ndivt,PR$rhatmax)
  s=log(PR$S+1)
  res=data.frame(lowerBound=s["2.5%"],estimate=s["50%"],upperBound=s["97.5%"],warningMessage=msg)
  fields=names(res)
  if(!is.null(out$cloneRateSharedMuts)){
    rbind(res %>% mutate(method="phylofit"),
          out$cloneRateML[fields[1:3]] %>% mutate(warningMessage="<NA>",method="cloneRateML"),
          out$cloneRateBdMCMC[fields] %>% mutate(method="cloneRateBirthDeathMCMC"),
          out$cloneRateSharedMuts[fields[1:3]] %>% mutate(warningMessage="<NA>",method="cloneRateMoltime")
    ) %>% mutate(ntip= out$cloneRateML$n,label2=out$label)
  }else{
    rbind(res %>% mutate(method="phylofit"),
          out$cloneRateML[fields[1:3]] %>% mutate(warningMessage="<NA>",method="cloneRateML"),
          out$cloneRateBdMCMC[fields] %>% mutate(method="cloneRateBirthDeathMCMC")
    ) %>% mutate(ntip= out$cloneRateML$n,label2=out$label)
  }
}

get_stored_benchmark_results=function(S,maxS=1e8,
                                      stub="/lustre/scratch127/casm/team273jn/nw14/CML_SELECTION_BENCHMARK_REVISION",
                                      nn=c(4,11,32,58,75,91),
                                      DOFILTER=FALSE,
                                      P=0.9999){
  #browser()
  zz=mclapply(nn,
              function(i) do.call("c",
                                  lapply(list.files(sprintf(
                                    "%s/maxs_%d__divperyear_1.00__ageacquistion_10__N_100000__nmut_%d_P_%5.4f/S_%d",stub,round(log(maxS+1)),
                                    i,P,S),full.names = TRUE),function(x){
                                      y=readRDS(x)
                                      if(length(y)>0){ ### or can list expected number of sims. length(y)==10
                                        y=lapply(y,function(w){w$file=x;w})
                                      }else{
                                        NULL
                                        }
                                      })
                                  )
              ,mc.cores = 4)
  idx=which(sapply(zz,length)>0)
  zz=zz[idx]
  names(zz)=paste0("n",nn)[idx]
  test=NULL
  # browser()
  for(FIELD in c("S","LN","tm")){
    test=rbind(test,do.call("rbind",lapply(nn[idx],function(n) as.data.frame(
      do.call("rbind",lapply(zz[[paste0("n",n)]],
                             function(x) 
                             {if(!DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5 && x$no_acf$rhatmax<1.05 )){x$no_acf[[FIELD]]}else{NULL}}))[,c(1,3,5)]) %>% 
        (function(y){names(y)=c("lb","median","ub")
        y %>% mutate(field=FIELD,i=1:dim(y)[1])}) %>% mutate(N=n))
    )
    
    )
  }
  test=test %>% mutate(method="phylofit")
  ##browser()
  test=rbind(test ,
             do.call("rbind",lapply(nn[idx],function(n)
               do.call("rbind",lapply(zz[[paste0("n",n)]],function(x) if(!DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5  )){exp(x$clonerateres[,1:3])-1}else{NULL} ))[,c(1,2,3)] %>% 
                 as.data.frame() %>%(function(y){names(y)=c("lb","median","ub");y %>% mutate(field="S",i=1:dim(y)[1])}) %>% 
                 mutate(N=n,method="cloneRate")
             )
             )
  )
  
  test=rbind(test,
             do.call("rbind",lapply(nn[idx],function(n)
               do.call("rbind",lapply(zz[[paste0("n",n)]],function(x) if(DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5 )){exp(x$bdres[,1:3])-1}else{NULL} ))[,c(1,2,3)] %>% 
                 as.data.frame() %>%(function(y){names(y)=c("lb","median","ub");y %>% mutate(field="S",i=1:dim(y)[1])}) %>% 
                 mutate(N=n,method="birthDeathMCMC")
             )
             )
  )
  test %>% filter(field=="S") %>% mutate(S_true=S,s_true=log(1+S),diff=log(median+1)-s_true,status=ifelse(lb>S_true | ub<S_true,"FAIL","PASS"))
}

get_stored_benchmark_results_ORIG=function(S,maxS=1e8,stub="/lustre/scratch127/casm/team273jn/nw14/CML_REVISION_SELECTION_BENCHMARK_NO_OFFSET_V3",nn=c(4,11,32,58,75,91),DOFILTER=FALSE){
  zz=mclapply(nn,
              function(i) do.call("c",
                                  lapply(list.files(sprintf(
                                    "%s/maxs_%d__divperyear_1.00__ageacquistion_10__N_100000__nmut_%d/S_%d",stub,round(log(maxS+1)),
                                    i,S),full.names = TRUE),function(x){y=readRDS(x);if(length(y)==20){y=lapply(y,function(w){w$file=x;w})}else{NULL}}))
              ,mc.cores = 4)
  idx=which(sapply(zz,length)>0)
  zz=zz[idx]
  names(zz)=paste0("n",nn)[idx]
  test=NULL
  # browser()
  for(FIELD in c("S","LN","tm")){
    test=rbind(test,do.call("rbind",lapply(nn[idx],function(n) as.data.frame(
      do.call("rbind",lapply(zz[[paste0("n",n)]],
                             function(x) 
                             {if(!DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5 && x$no_acf$rhatmax<1.05 )){x$no_acf[[FIELD]]}else{NULL}}))[,c(1,3,5)]) %>% 
        (function(y){names(y)=c("lb","median","ub")
        y %>% mutate(field=FIELD,i=1:dim(y)[1])}) %>% mutate(N=n))
    )
    
    )
  }
  test=test %>% mutate(method="phylofit")
  ##browser()
  test=rbind(test ,
             do.call("rbind",lapply(nn[idx],function(n)
               do.call("rbind",lapply(zz[[paste0("n",n)]],function(x) if(!DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5  )){exp(x$clonerateres[,1:3])-1}else{NULL} ))[,c(1,2,3)] %>% 
                 as.data.frame() %>%(function(y){names(y)=c("lb","median","ub");y %>% mutate(field="S",i=1:dim(y)[1])}) %>% 
                 mutate(N=n,method="cloneRate")
             )
             )
  )
  
  test=rbind(test,
             do.call("rbind",lapply(nn[idx],function(n)
               do.call("rbind",lapply(zz[[paste0("n",n)]],function(x) if(DOFILTER || (x$no_acf$rhatmax<1.05 && x$no_acf$ndivt<5 )){exp(x$bdres[,1:3])-1}else{NULL} ))[,c(1,2,3)] %>% 
                 as.data.frame() %>%(function(y){names(y)=c("lb","median","ub");y %>% mutate(field="S",i=1:dim(y)[1])}) %>% 
                 mutate(N=n,method="birthDeathMCMC")
             )
             )
  )
  test %>% filter(field=="S") %>% mutate(S_true=S,s_true=log(1+S),diff=log(median+1)-s_true,status=ifelse(lb>S_true | ub<S_true,"FAIL","PASS"))
}