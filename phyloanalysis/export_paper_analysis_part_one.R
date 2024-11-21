##  Script generating derived tree based results.
source("load_and_annotate_tree.R")
source("PHYLOFIT.R")
source("revised_sigs.R")
source("cml_sigs.R")
library("sigfit")
library(RVAideMemoire)
library(DescTools)
library(ggpubr)
library(cloneRate)
VERSION="REVISION_V2"
writeLines(VERSION,"../export/VERSION.txt")
pdd_cache="../cache/PDD.RDS"
## We set a separate seed for each stage of the analysis to facilitate 
## rerunning from saved interim results

# Section 1
## Read in trees with assigned mutations
## Convert molecular time trees to time-based trees with clade-specific mutation rates
set.seed(1234567)
## Split PDD file into git friendly file size.
if(file.exists("../cache/PDD_A.RDS") && !file.exists(pdd_cache)){
  pdd=c(readRDS("../cache/PDD_A.RDS"),readRDS("../cache/PDD_B.RDS"))
  saveRDS(pdd,pdd_cache)
  rm(pdd)
}
if(file.exists(pdd_cache)){
  cat("loading stored PDD data\n")
  PDD=readRDS(pdd_cache)
}else{
  PDD=get_PDD()  ### Relies on data with germline variants so this execution path not supported.
  cat("removing germline data...")
  PDD=lapply(PDD,function(x){x[["inf"]]=NULL;x})
  names(PDD)=sapply(PDD,function(PD) PD$patient)
  cat("caching PDD to ",pdd_cache,"\n")
  saveRDS(PDD,file=pdd_cache)
}
if(file.exists("../cache/PDDx.RDS")){
  PDD=readRDS("../cache/PDDx.RDS")
}else{
  PDDo=PDD
  PATIENTS=setdiff(PATIENTS,c("PD57332","PD57333"))
  PDD=PDD[PATIENTS]
  
  PDD=lapply(PDD,function(PD) {PD$nodes =PD$nodes %>% mutate(status=ifelse(driver=="BCR::ABL1",1,-1));PD})
  PDD=lapply(PDD,function(x) wraptreefit(x,niter=10000,xcross = 0.99,b.fit.null = FALSE,stan_control=list(adapt_delta=0.99)))
  PDD=lapply(PDD,function(x) wraptreefit(x,niter=10000,xcross = 0.99,b.fit.null =TRUE,stan_control=list(adapt_delta=0.99)))
  #browser()
  ## Special treatment for PD57322
  PD=PDDo$PD57332 ##get_PD("PD57332")
  PD$nodes =PD$nodes %>% mutate(status=ifelse(driver=="BCR::ABL1",1,-1))
  idx=which(PD$pdx$agedf$tip.label=="zeros")
  PD2=PD
  PD2$pdx$agedf$age_at_sample_pcy[idx]=1000 ### We are making a long outgroup
  idx2=match(idx,PD2$pdx$tree_ml$edge[,2])
  PD2$pdx$tree_ml$el.snv.local.filtered[idx2]=round(PD2$trunk_params$mean.WT.rate*1000+PD2$trunk_params$additional.muts.at.birth)
  PD2=wraptreefit(PD2,niter=10000,xcross = 0.99,b.fit.null = FALSE,stan_control=list(adapt_delta=0.99))
  PD2=wraptreefit(PD2,niter=10000,xcross = 0.99,b.fit.null =TRUE,stan_control=list(adapt_delta=0.99))
  
  PD2$fit$poisson_tree$altmodel$ultratree$edge.length[idx2]=0.00001
  PD2$fit$poisson_tree$nullmodel$ultratree$edge.length[idx2]=0.00001
  PD$fit=PD2$fit
  PDD$PD57332=PD
  ## PD57333 is all wild type - so just fit the null model
  PD=PDDo$PD57333
  PD=wraptreefit(PD,niter=10000,xcross = 0.99,b.fit.null =TRUE,stan_control=list(adapt_delta=0.99))
  PD$fit$poisson_tree$altmodel=PD$fit$poisson_tree$nullmodel
  PDD$PD57333=PD
  #cat("saving PDD\n")
  #saveRDS(PDD,"../cache/PDD_step1.RDS")
  ages=sapply(PDD,function(x) min(x$pdx$agedf$age_at_sample_exact[x$pdx$agedf$age_at_sample_exact>1],na.rm=TRUE))
  PDD=PDD[order(ages)]
  #PDD=lapply(PDD,function(PD){PD$long_id=sprintf("%s:%s",PD$patient,patients$SHORTNAME[match(PD$patient,patients$PDID)]);PD})
  #PDD2=lapply(PDD,function(x){x$inf=NULL;x$fit$poisson_tree$altmodel$fullres=NULL;x$fit$poisson_tree$nullmodel$fullres=NULL;x})
  saveRDS(PDD,"../cache/PDDx.RDS")
}
# Section 2
## Estimate growth rates all available methods 




set.seed(1234567)
smax=30
maxSyear=exp(smax)-1
PARAM=get_sel_params(PDD) %>% mutate(smax=smax,Smax=exp(smax)-1)
if(file.exists("../cache/phylores.RDS")){
  phylores=readRDS("../cache/phylores.RDS")
}else{

  donors=PARAM %>% filter(donor!="PD57335") %>% pull("donor")
  phylores=sapply(donors,function(thisdonor){
    param=PARAM %>% filter(donor==thisdonor)
    PD=PDD[[param$donor]]
    context=paste0(PD$patient,"_selection_altmodel_")
    tree=PD$fit$poisson_tree[[param$model]]$ultratree
    ### Now run cloneRate, birthDeathMCMC and birthDeathMRCA
    
    hashme=list(tree=tree,param=param)
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=param$node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = log(param$maxS+1),nCores = 4,chainLength = 4000)
    
    muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
    rate=PD$fit$poisson_tree[[param$model]]$lambda$median[2]
    cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
    
    res=get_cached_result(context,hashme)
    if(!is.null(res)){
      cat(param$donor,":Using cached result..\n")
      return(res)
    }else{
      out=fit_clade(tree,
                    node=param$node,
                    nmutcolony =-1,
                    nwtcolony =-1,
                    maxt=param$maxt,
                    stan.control = list(adapt_delta=param$adapt_delta,max_treedepth=14),
                    maxSYear = param$maxS,
                    niter = 10000)
      out$cloneRateML=cloneRateML
      out$cloneRateBdMCMC=cloneRateBdMCMC
      out$cloneRateSharedMuts=cloneRateShared
      
      
      cache_result(context,hashme,out)
      return(out)
    }
  },USE.NAMES = TRUE,simplify = FALSE
  )
  ## Fix PD57335 because of multiple timepoints
  param=PARAM %>% filter(donor=="PD57335")
  PD=PDD$PD57335
  tree=PD$fit$poisson_tree[[param$model]]$ultratree
  tree=drop.tip(tree,PD$pdx$agedf %>% filter(age_at_sample_pcy>47)  %>% pull("tip.label"))
  test=fit_clade(tree,node=71,nmutcolony =-1,nwtcolony=-1,maxt=46.72827,nchain = 4,stan.control = list(adapt_delta=0.999,max_treedepth=14),maxSYear = maxSyear,niter = 10000)
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=71,tree = tree))
  muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
  rate=PD$fit$poisson_tree[[param$model]]$lambda$median[2]
  cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
  test$cloneRateML=cloneRateML
  test$cloneRateBdMCMC=cloneRateBdMCMC
  test$cloneRateSharedMuts=cloneRateShared
  phylores[["PD57335"]]=test
  ### Include the extra time point
  
  
  saveRDS(phylores,"../cache/phylores.RDS")
}

## Estimate additional growth rates of interest, e.g. for CH.
set.seed(1234567)
EXTRASEL=list()
if(!file.exists("../cache/extrasel.RDS")){
  smax=30
  maxSyear=exp(smax)-1
  tree=PDD$PD57335$fit$poisson_tree$altmodel$ultratree
  tree=drop.tip(tree,PDD$PD57335$pdx$agedf %>% filter(age_at_sample_pcy<47)  %>% pull("tip.label"))
  test=fit_clade(tree,node=40,nchain=4,
                 nmutcolony =-1,
                 nwtcolony=-1,
                 maxt=46.72827,
                 stan.control = list(adapt_delta=0.999,max_treedepth=14),
                 maxSYear = maxSyear,
                 niter = 10000)
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=40,tree = tree))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 4000)
  test$cloneRateML=cloneRateML
  test$cloneRateBdMCMC=cloneRateBdMCMC
  test$label="BCR::ABL:Time point 2"
  EXTRASEL[["PD57335_49yr"]]=test
  smax=5
  maxSyear=exp(smax)-1
  ### PD51635
  ### Age 70..
  ## 127: TET2:p.P1889L
  ## 163: DNMT3A: o1310T
  ## 214: DNMT3A:p.D845fs*8
  labels=list("127"="TET2:p.P1889L",
              "163"="DNMT3A: o1310T",
              "214"= "DNMT3A:p.D845fs*8",
              "162"="DNMT3A: o1310T",
              "82"="DNMT3A:p.F868L",
              "111"="No Driver",
              "39"="RUNX1:BCR::ABL1",
              "192"="BCR::ABL1:Subclone",
              "8"="RUNX1(Singleton):BCR::ABL1",
              "68"="BCR::ABL1:Early Timepoint.Remove 3 late coalescences"
              )
  tree=PDD$PD51635$fit$poisson_tree$altmodel$ultratree
  nodes=c(127,163,214)
  for(node in nodes){ ## 127,163,214
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
    out=fit_clade(tree,
                  node=node,nchain = 4,
                  nmutcolony =-1,
                  nwtcolony =-1,
                  maxt=log(1e6)/0.1,
                  stan.control = list(adapt_delta=0.999,max_treedepth=14),
                  maxSYear = maxSyear,
                  niter = 10000)
    out$cloneRateML=cloneRateML
    out$cloneRateBdMCMC=cloneRateBdMCMC
    out$label=labels[[as.character(node)]]
    EXTRASEL[[sprintf("PD51635_%d",node)]]=out
  }
  node=192
  smax=5
  maxSyear=exp(smax)-1
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
  out=fit_clade(tree,
                node=node,nchain = 4,
                nmutcolony =-1,
                nwtcolony =-1,
                maxt=PARAM %>% filter(donor=="PD51635") %>% pull("maxt"),
                stan.control = list(adapt_delta=0.999,max_treedepth=14),
                maxSYear = maxSyear,
                niter = 10000)
  out$cloneRateML=cloneRateML
  out$cloneRateBdMCMC=cloneRateBdMCMC
  out$label=labels[[as.character(node)]]
  EXTRASEL[[sprintf("PD51635_%d",node)]]=out
  ## Remove late coalescence
  tree=drop.tip(PDD$PD51635$fit$poisson_tree$altmodel$ultratree,"lo0066")
  nodes=162
  for(node in nodes){ 
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
    out=fit_clade(tree,
                  node=node,nchain = 4,
                  nmutcolony =-1,
                  nwtcolony =-1,
                  maxt=log(1e6)/0.1,
                  stan.control = list(adapt_delta=0.999,max_treedepth=14),
                  maxSYear = maxSyear,
                  niter = 10000)
    out$cloneRateML=cloneRateML
    out$cloneRateBdMCMC=cloneRateBdMCMC
    out$label=labels[[as.character(node)]]
    EXTRASEL[[sprintf("PD51635_%d_REMOVE_LATE",node)]]=out
  }
  ### PD57333
  smax=5
  maxSyear=exp(smax)-1
  ## 82: DNMT3A:p.F868L
  ## 111: <NONE>
  tree=PDD$PD57333$fit$poisson_tree$altmodel$ultratree
  for(node in c(82,111)){
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
    out=fit_clade(tree,
                  node=node,nchain=4,
                  nmutcolony =-1,
                  nwtcolony =-1,
                  maxt=log(1e6)/0.1,
                  stan.control = list(adapt_delta=0.999,max_treedepth=14),
                  maxSYear = maxSyear,
                  niter = 10000)
    out$cloneRateML=cloneRateML
    out$cloneRateBdMCMC=cloneRateBdMCMC
    out$label=labels[[as.character(node)]]
    EXTRASEL[[sprintf("PD57733_%d",node)]]=out
  }
  smax=30
  maxSyear=exp(smax)-1
  tree=PDD$PD57332$fit$poisson_tree$altmodel$ultratree
  for(node in c(39)){
    MAXT=PARAM %>% filter(donor=="PD57332") %>% pull("maxt")
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
    out=fit_clade(tree,
                  node=node,nchain=4,
                  nmutcolony =-1,
                  nwtcolony =-1,
                  maxt=MAXT,
                  stan.control = list(adapt_delta=0.999,max_treedepth=14),
                  maxSYear = maxSyear,
                  niter = 10000)
    out$cloneRateML=cloneRateML
    out$cloneRateBdMCMC=cloneRateBdMCMC
    out$label=labels[[as.character(node)]]
    EXTRASEL[[sprintf("PD57332_%d",node)]]=out
  }
  ### Some additional tests for Merge RUNX1 and non-RUNX1 representing RUNX1 with just
  ## one clone
  set.seed(12345667)
  keep.tips.1=c("lo0003","lo0007","lo0017","lo0016","lo0023","zeros")
  smax=30
  maxSyear=exp(smax)-1
  tree=keep.tip(PDD$PD57332$fit$poisson_tree$altmodel$ultratree,keep.tips.1)
  for(node in c(8)){
    MAXT=PARAM %>% filter(donor=="PD57332") %>% pull("maxt")
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
    out=fit_clade(tree,
                  node=node,nchain=4,
                  nmutcolony =-1,
                  nwtcolony =-1,
                  maxt=MAXT,
                  stan.control = list(adapt_delta=0.999,max_treedepth=14),
                  maxSYear = maxSyear,
                  niter = 10000)
    out$cloneRateML=cloneRateML
    out$cloneRateBdMCMC=cloneRateBdMCMC
    out$label=labels[[as.character(node)]]
    EXTRASEL[[sprintf("PD57332_%d",node)]]=out
  }
  #### Remove late coalescence experiment.
  tree=PDD$PD57335$fit$poisson_tree$altmodel$ultratree
  tree=keep.tip(tree,PDD$PD57335$pdx$agedf %>% filter(age_at_sample_pcy<47)  %>% pull("tip.label"))
  tree=drop.tip(tree,c("a64","a75","a07"))
  test=fit_clade(tree,node=68,nchain=4,
                 nmutcolony =-1,
                 nwtcolony=-1,
                 maxt=46.72827,
                 stan.control = list(adapt_delta=0.999,max_treedepth=14),
                 maxSYear = maxSyear,
                 niter = 10000)
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=68,tree = tree))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 4000)
  test$cloneRateML=cloneRateML
  test$cloneRateBdMCMC=cloneRateBdMCMC
  test$label="BCR::ABL:Early Time Point: Remove 3 late coal."
  EXTRASEL[["PD57335_REMOVE_LATE"]]=test
  saveRDS(EXTRASEL,"../cache/extrasel.RDS")
}else{
  EXTRASEL=readRDS("../cache/extrasel.RDS")
}
## Estimate growth rates under the assumption that the wild-type and mutant mutation rates are the same.
set.seed(1234567)
smax=30
maxSyear=exp(smax)-1
PARAM=get_sel_params(PDD) %>% mutate(model="nullmodel",smax=smax,Smax=exp(smax)-1)
if(file.exists("../cache/phyloresNULL.RDS")){
  phylores=readRDS("../cache/phyloresNULL.RDS")
}else{
  donors=PARAM %>% filter(donor!="PD57335") %>% pull("donor")
  phylores=sapply(donors,function(thisdonor){
    param=PARAM %>% filter(donor==thisdonor)
    PD=PDD[[param$donor]]
    context=paste0(PD$patient,"_selection_nullmodel_")
    tree=PD$fit$poisson_tree[[param$model]]$ultratree
    ### Now run cloneRate, birthDeathMCMC and birthDeathMRCA
    
    hashme=list(tree=tree,param=param)
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=param$node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = log(param$maxS+1),nCores = 4,chainLength = 4000)
    muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
    rate=PD$fit$poisson_tree[[param$model]]$lambda$median[1]
    cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
    res=get_cached_result(context,hashme)
    if(!is.null(res)){
      cat(param$donor,":Using cached result..\n")
      return(res)
    }else{
      out=fit_clade(tree,
                    node=param$node,
                    nmutcolony =-1,
                    nwtcolony =-1,
                    maxt=param$maxt,
                    stan.control = list(adapt_delta=param$adapt_delta,max_treedepth=14),
                    maxSYear = param$maxS,
                    niter = 10000)
      out$cloneRateML=cloneRateML
      out$cloneRateBdMCMC=cloneRateBdMCMC
      out$cloneRateSharedMuts=cloneRateShared
      cache_result(context,hashme,out)
      return(out)
    }
  },USE.NAMES = TRUE,simplify = FALSE
  )
  ## Fix PD57335 because of multiple timepoints
  param=PARAM %>% filter(donor=="PD57335")
  PD=PDD$PD57335
  tree=PD$fit$poisson_tree[[param$model]]$ultratree
  tree=drop.tip(tree,PD$pdx$agedf %>% filter(age_at_sample_pcy>47)  %>% pull("tip.label"))
  test=fit_clade(tree,node=71,nmutcolony =-1,nwtcolony=-1,maxt=46.72827,nchain = 4,stan.control = list(adapt_delta=0.999,max_treedepth=14),maxSYear = maxSyear,niter = 10000)
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=71,tree = tree))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
  
  muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
  rate=PD$fit$poisson_tree[[param$model]]$lambda$median[1]
  cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
  
  test$cloneRateML=cloneRateML
  test$cloneRateBdMCMC=cloneRateBdMCMC
  test$cloneRateSharedMuts=cloneRateShared
  phylores[["PD57335"]]=test
  ### Include the extra time point
  saveRDS(phylores,"../cache/phyloresNULL.RDS")
}
## Finally estimate growth rate with max prior s=10
set.seed(1234567)
smax=5
maxSyear=exp(smax)-1
PARAM=get_sel_params(PDD) %>% mutate(smax=smax,maxS=exp(smax)-1)
if(file.exists("../cache/phylores_smax10.RDS")){
  phylores=readRDS("../cache/phylores_smax10.RDS")
}else{
  
  donors=PARAM %>% filter(donor!="PD57335") %>% pull("donor")
  phylores=sapply(donors,function(thisdonor){
    param=PARAM %>% filter(donor==thisdonor)
    PD=PDD[[param$donor]]
    context=paste0(PD$patient,"_selection_altmodel_")
    tree=PD$fit$poisson_tree[[param$model]]$ultratree
    ### Now run cloneRate, birthDeathMCMC and birthDeathMRCA
    
    hashme=list(tree=tree,param=param)
    justmutcladetree=keep.tip(tree,get_samples_in_clade(node=param$node,tree = tree))
    cloneRateML=maxLikelihood(justmutcladetree)
    cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = log(param$maxS+1),nCores = 4,chainLength = 4000)
    
    muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
    rate=PD$fit$poisson_tree[[param$model]]$lambda$median[2]
    cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
    
    res=get_cached_result(context,hashme)
    if(!is.null(res)){
      cat(param$donor,":Using cached result..\n")
      return(res)
    }else{
      out=fit_clade(tree,
                    node=param$node,
                    nmutcolony =-1,
                    nwtcolony =-1,
                    maxt=param$maxt,
                    stan.control = list(adapt_delta=param$adapt_delta,max_treedepth=14),
                    maxSYear = param$maxS,
                    niter = 10000)
      out$cloneRateML=cloneRateML
      out$cloneRateBdMCMC=cloneRateBdMCMC
      out$cloneRateSharedMuts=cloneRateShared
      
      
      cache_result(context,hashme,out)
      return(out)
    }
  },USE.NAMES = TRUE,simplify = FALSE
  )
  ## Fix PD57335 because of multiple timepoints
  param=PARAM %>% filter(donor=="PD57335")
  PD=PDD$PD57335
  tree=PD$fit$poisson_tree[[param$model]]$ultratree
  tree=drop.tip(tree,PD$pdx$agedf %>% filter(age_at_sample_pcy>47)  %>% pull("tip.label"))
  test=fit_clade(tree,node=71,nmutcolony =-1,nwtcolony=-1,maxt=46.72827,nchain = 4,stan.control = list(adapt_delta=0.999,max_treedepth=14),maxSYear = maxSyear,niter = 10000)
  justmutcladetree=keep.tip(tree,get_samples_in_clade(node=71,tree = tree))
  muttree=keep.tip(PD$pdx$tree_ml,justmutcladetree$tip.label)
  rate=PD$fit$poisson_tree[[param$model]]$lambda$median[2]
  cloneRateShared=suppressWarnings(sharedMuts(muttree,nu=rate))
  cloneRateML=maxLikelihood(justmutcladetree)
  cloneRateBdMCMC=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 10000)
  test$cloneRateML=cloneRateML
  test$cloneRateBdMCMC=cloneRateBdMCMC
  test$cloneRateSharedMuts=cloneRateShared
  phylores[["PD57335"]]=test
  ### Include the extra time point
  saveRDS(phylores,"../cache/phylores_smax10.RDS")
}
phylores=readRDS("../cache/phylores.RDS")




if(!file.exists("../cache/PDDxs.RDS")){
  pcawg=get_extended_pcawg()
  PDD=lapply(PDD,function(PD){
    PD$pdx$dat$info=get_clade_sigs(PD,prior=pcawg[, c("SBS1","SBSblood","SBS18")],extracats = c("MutantI","Top100"))
    PD
  }
  )
  PDD=lapply(PDD,function(PD){
    add_v2_sig_to_pdx(PD,
                      prior=pcawg[, c("SBS1","SBSblood","SBS18")],
                      info=PD$pdx$dat$info
    )
  }
  )
  saveRDS(PDD,"../cache/PDDxs.RDS")
}else{
  PDD=readRDS("../cache/PDDxs.RDS")
}
## Wrangle data into a fat per patient data.frame.
source("lineages.R")
PATIENTS=names(PDD)
#browser()
MT_PATIENTS=setdiff(PATIENTS,"PD57333")
## The following interrogate the rtreefit results to get the timing plus CI of the BCR-ABL branch
bcatimings=do.call("rbind",lapply(PDD[MT_PATIENTS],function(PD) {
  lineages=get_lineages2(PD$pdx$tree_ml,PD$nodes$node[which(PD$nodes$driver=="BCR::ABL1")])
  lineages=lineages %>% mutate(start_mutcount=start,end_mutcount=end, n_samples_in_clade=N) %>% dplyr::select(node,start_mutcount,end_mutcount, n_samples_in_clade) 
  tmp=do_tree_branch_length_ci_etc(PD$pdx$tree_ml,PD$fit$poisson_tree$altmodel,lineages)$lineages
  tmp$patient=PD$patient
  tmp})
)

## The diagnosis timing is the age at first sample time point for this cohort.
## but there are a couple of exceptions (DIAGNOSIS_GAP) encoded in patients.txt
## Estimate aternce based on timing of mrca of bcr-abl colonies
pinf=read.table("../data/patients.txt",header = T,stringsAsFactors = FALSE)
## Add in age at first sample
pinf=pinf %>% mutate(age_at_first_sample_pcy=sapply(PDD[pinf$PDID],function(PD) min(PD$pdx$agedf$age_at_sample_pcy[PD$pdx$agedf$age_at_sample_pcy>0.1])))
bcatimings=pinf %>% left_join(bcatimings,by=c("PDID"="patient")) %>% mutate(Donor=PDID) %>% mutate(t_diagnosis=ifelse(is.na(t_diagnosis),age_at_first_sample_pcy,t_diagnosis))
bcatimings=bcatimings %>% mutate(t_first_sample=t_diagnosis) %>% mutate(t_diagnosis=t_diagnosis-DIAGNOSIS_GAP)
bcatimings=bcatimings %>% mutate(latency_median=t_diagnosis-t_upper_median,latency_lb=t_diagnosis-t_upper_ub95,latency_ub=t_diagnosis-t_upper_lb95)

## Extract Wildtype (WT) bcr-abl (bcrabl( specific estimates of mutation rate (calculated assuming rate operates after mcra of BCRABL colonies)
PDD$PD57333$fit$poisson_tree$altmodel$lambda=lapply(PDD$PD57333$fit$poisson_tree$altmodel$lambda,function(x){x[2]=NA;names(x)=c("lambda[1]","lambda[2]");x})
lmean=sapply(PDD,function(PD) PD$fit$poisson_tree$altmodel$lambda$mean)
lb=sapply(PDD,function(PD) PD$fit$poisson_tree$altmodel$lambda$lb)
ub=sapply(PDD,function(PD) PD$fit$poisson_tree$altmodel$lambda$ub)


rownames(ub)=c("WT_ub","BCRABL_ub")
rownames(lb)=c("WT_lb","BCRABL_lb")
rownames(lmean)=c("WT_mean","BCRABL_mean")
rates=as.data.frame(t(rbind(lmean,ub,lb))) %>% rownames_to_column("Donor")
### Add the rates to bcatimings ##
bcadat=bcatimings %>% left_join(rates)

ratediff=sapply(PDD[MT_PATIENTS],
                function(PD){
                  tt=rstan::extract(PD$fit$poisson_tree$altmodel$fullres,pars=c("lambda"));
                  c(mean=mean(tt$lambda[,2]-tt$lambda[,1]),quantile(tt$lambda[,2]-tt$lambda[,1],c(0.025,0.5,0.975)))})
ratediff=as.data.frame(t(as.data.frame(ratediff))) %>% rownames_to_column(var="Donor") %>% (function(x){colnames(x)[-1]=sprintf("lambda_diff_%s",c("mean","lb","median","ub"));x})
bcadat=bcadat %>% left_join(ratediff)
bcadat$age_at_diagnosis=round(bcadat$t_diagnosis-0.72827)
bcadat=bcadat %>% mutate(label=sprintf("%s",age_at_diagnosis)) %>% mutate(label2=sprintf("%s(%s)",Donor,age_at_diagnosis))

patients=read.table("../data/patients.txt",head=TRUE,sep="\t")
patients$color="" ##gg_color_hue(dim(patients)[1])
patients$Donor=patients$PDID
patients=patients %>% inner_join(bcadat %>% dplyr::select(Donor,label,label2))
patients=patients[order(match(patients$PDID,PATIENTS)),]

dfsel=do.call("rbind",
              lapply(names(phylores),
                     function(x) as.data.frame(summary(phylores[[x]]$res)$summary) %>% rownames_to_column(var="field") %>% mutate(Donor=x)
              )
)
dfsel=dfsel %>% filter(field=="S") %>% left_join(patients)
colnames(dfsel)[c(2,5,9)]=c("S_ltt","S_ltt_lb","S_ltt_ub")
dfsel=dfsel %>% dplyr::rename("S_ltt_mean"="S_ltt") %>% dplyr::rename("S_ltt_median"="50%") %>% mutate(S_ltt=S_ltt_median) 

zzz=lapply(PDD,function(PD) list(counts=PD$pdx$dat$info$ctx %>% dplyr::select(ctx,count,category) %>% 
                                   pivot_wider(names_from = category,values_from=count) %>% column_to_rownames("ctx")))

tmp=do.call("rbind",lapply(names(zzz),function(donor) zzz[[donor]]$counts %>% as.matrix() %>% (function(x){data.frame(n_C_T=x["C>T at CpG",],n_total=colSums(x),Donor=donor) %>% rownames_to_column("category")})))
cats=unique(tmp$category)
tmp=tmp %>% pivot_wider(names_from = category,values_from = c("n_C_T","n_total"))
for(CAT in cats){tmp=tmp %>% dplyr::rename(!!paste0("n_",CAT,"_total"):=paste0("n_total_",CAT));tmp[[paste0("C_T_",CAT)]]=tmp[[paste0("n_C_T_",CAT)]]/tmp[[paste0("n_",CAT,"_total")]]}
c_to_t_at_cpg=tmp


dfcomb=bcadat %>% left_join(c_to_t_at_cpg)
dfcomb=dfcomb %>% left_join(dfsel %>% dplyr::select(-setdiff(colnames(patients),"Donor")))
## Make sure everything is ordered correctly.
dfcomb=dfcomb %>% mutate(Patient=factor(Donor,levels=patients$PDID)) %>% 
  mutate(Donor=factor(Patient,levels=patients$PDID)) %>% mutate(label=factor(label,levels=patients$label)) %>%
  mutate(label2=factor(label2,levels=patients$label2))
dfcomb=dfcomb %>% left_join(patients %>% dplyr::select(Donor,color))

fields1=c("Donor","DIAGNOSIS_GAP","age_at_diagnosis","node","start_mutcount","end_mutcount","n_samples_in_clade","t_lower_lb95","t_lower_ub95","t_lower_median","t_upper_lb95", "t_upper_median","t_upper_ub95")
fields2=c("t_first_sample","t_diagnosis","latency_median","latency_lb","latency_ub","S_ltt_median","S_ltt_mean","S_ltt_lb","S_ltt_ub","WT_mean","BCRABL_mean","WT_ub","BCRABL_ub","WT_lb", "BCRABL_lb","lambda_diff_mean","lambda_diff_median","lambda_diff_lb","lambda_diff_ub","C_T_Mutant","C_T_Trunk","C_T_WildType","n_C_T_Mutant","n_C_T_Trunk","n_C_T_WildType","n_Mutant_total","n_Trunk_total","n_WildType_total")

## Add in s

write.table(dfcomb[,c(fields1,fields2)] %>% 
              mutate(s_median=log(1+S_ltt_median),s_lb=log(1+S_ltt_lb),s_ub=log(1+S_ltt_ub)),
            file =sprintf("../export/per_donor_summary_v%s.txt",VERSION) ,sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

dfs=do.call("rbind",
            lapply(PDD,function(PD) 
              PD$pdx$agedf %>% filter(tip.label != "zeros") %>% left_join(PD$pdx$cfg,c("tip.label"="SHORT_LABEL")))
)
write.table(dfs %>% mutate(sample=LABEL) %>% dplyr::select(patient,LABEL,tip.label,age_at_sample_exact,meanvaf,rho.bb,meandepth,sensitivity,driver3),
            file =sprintf("../export/retained_sample_summary_v%s.txt" ,VERSION),sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE
)
## Get burden as tree tip height (tree is adjusted and corresponds to SNV burden)
dfox=get_burden(PDD,chk.ctcpg = FALSE)
write.table(dfox,file =sprintf("../export/per_colony_burden_etc_v%s.txt",VERSION) ,sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

ch=do.call("c",lapply(get_lineages(PDD$PD51635$pdx$tree_ml,100)$node,function(node){get_samples_in_clade(node,PDD$PD51635$pdx$tree_ml)}))
tmp=PDD$PD51635$pdx$agedf %>% left_join(data.frame(tip.label=ch,status="CH")) %>% mutate(status=ifelse(is.na(status),"NOCH",status)) %>% dplyr::select(tip.label,driver3,status) %>% inner_join(PDD$PD51635$pdx$cfg %>% dplyr::select(LABEL,SHORT_LABEL),by=c("tip.label"="SHORT_LABEL"))
write.table(tmp,file = sprintf("../export/PD51635.CH.v%s.txt",VERSION),row.names = FALSE,quote=FALSE,sep="\t",col.names=TRUE)

### Generate tables from EXTRASEL and phylores
rr=do.call("rbind",lapply(names(EXTRASEL),function(y) get_all_results(EXTRASEL[[y]]) %>% mutate(label=y)
                          )
  ) %>% 
  dplyr::select(-warningMessage)
rr=rr %>% mutate(LABEL=sprintf("%s:%s",label,label2))
rr=rr %>% mutate(S_ltt_median=exp(estimate)-1,   S_ltt_mean=NA,  S_ltt_lb=exp(lowerBound)-1,S_ltt_ub=exp(upperBound)-1)
R1=do.call("rbind",lapply(names(phylores),function(x) get_all_results(phylores[[x]]) %>% mutate(donor=x)))
R1=R1 %>% mutate(S_ltt_median=exp(estimate)-1,   S_ltt_mean=NA,  S_ltt_lb=exp(lowerBound)-1,S_ltt_ub=exp(upperBound)-1)
R1=R1 %>% mutate(label=sprintf("%s: BCR::ABL1",donor)) %>% mutate(label2=label)
tmp=rbind(R1 %>% dplyr::select(-warningMessage) %>% mutate(LABEL=label),rr %>% mutate(donor=LABEL) %>% mutate(LABEL=gsub("_[0-9]+","",LABEL)))
tmp=tmp %>% mutate(LABEL=gsub("_[0-9]+","",LABEL))
tmp=tmp %>% mutate(donor=substr(donor,1,7))
#write.table(tmp,sprintf("../export/growth_rates_clonal_haem_etc_ALL_METHODS_%s.txt",VERSION),
#            col.names = TRUE,row.names=FALSE,sep="\t")
zz=tmp %>% dplyr::select(lowerBound,estimate,upperBound,method,ntip,donor,LABEL) %>% 
  mutate(s_text=sprintf("%3.1f(%3.1f-%3.1f)",
                        estimate,lowerBound,upperBound),
         S_text=sprintf("%3.1f(%3.1f-%3.1f)",
                        exp(estimate)-1,exp(lowerBound)-1,exp(upperBound)-1),
         method=factor(method,levels=c("phylofit","cloneRateBirthDeathMCMC","cloneRateML","cloneRateMoltime"))) %>% 
  pivot_wider(values_from = c("estimate","lowerBound","upperBound","s_text","S_text"),names_from = method)
phylo.summary=zz %>% dplyr::select(donor,LABEL,ntip,estimate_phylofit,lowerBound_phylofit,upperBound_phylofit) %>% 
  mutate(s=estimate_phylofit,s_lb=lowerBound_phylofit,s_ub=upperBound_phylofit) %>% 
  mutate(doubling_time=log(2)/s,S=exp(s)-1,S_lb=exp(s_lb)-1,S_ub=exp(s_ub)-1) %>% dplyr::select(donor,LABEL,ntip,S,S_lb,S_ub,s,s_lb,s_ub,doubling_time)
write.table(phylo.summary,sprintf("../export/growth_rates_PHYLOFIT_%s.txt",VERSION),col.names = TRUE,row.names=FALSE,sep="\t")
write.table(zz,sprintf("../export/growth_rates_ALL_%s.txt",VERSION),col.names = TRUE,row.names=FALSE,sep="\t")
