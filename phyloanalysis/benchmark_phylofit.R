require("rsimpop")
require("parallel")
require("rstan")
require("cloneRate")
#source("birthDeathStanModels.R")
PHYLOFIT_NITER=10000
source("PHYLOFIT.R")
run_benchmark_sim=function(S=0.4,N=5e4,divrate_per_year=1,nyears_driver_acquisition = 10,nyears = 33,minprop=0.05){
  run_selection_sim(0.1,
                    final_division_rate = divrate_per_year/365,
                    target_pop_size = N,
                    nyears_driver_acquisition = nyears_driver_acquisition,
                    nyears = nyears,fitness=log(1+S)/(divrate_per_year),
                    minprop = minprop)
}

get_subtree=function(selsim,ncolony,nmut_colony=-1){
  nmutcolony=-1
  kk=0
  if(nmut_colony<0){
    while(nmutcolony<2 && kk<100){
      st=get_subsampled_tree(selsim,ncolony)
      idx=which(st$events$driverid==1)
      kk=kk+1
      if(length(idx)!=1){
        next()
      }
      node=st$events$node[idx]
      nmutcolony=length(get_samples_in_clade(node,st))
      cat("node=",node,"\n")
    }
  }else{
    st=get_subsampled_tree_fixed_prop(selsim,ncolony,nmut_colony)
    idx=which(st$events$driverid==1)
    if(length(idx)!=1){
      stop("No mutant clade")
      #next()
    }
    node=st$events$node[idx]
    nmutcolony=length(get_samples_in_clade(node,st))
  }
  #if(nmutcolony<2){
  #  return(NULL)
  #}
  if(length(selsim$cfg$info$population)>2){
    acf=selsim$cfg$info$population[3]/(selsim$cfg$info$population[2]+selsim$cfg$info$population[3])
  }else{
    acf=1
  }
  #st=get_subsampled_tree(selsim,ncolony)
  idx=which(st$events$driverid==1)
  if(length(idx)!=1){
    stop("No mutant clade")
  }
  node=st$events$node[idx]
  ts=st$events$ts[idx]
  idx=which(st$edge[,2]==node)
  nmutcolony=length(get_samples_in_clade(node,st))
  nwtcolony=length(st$tip.label)-1-nmutcolony ##subtract 1 for the outgroup
  T=max(st$timestamp)/365
  
  #Nhsc=selsim$cfg$compartment$popsize[2]
  #if(b.use.perfect.priors){
  #  LNmin=log(0.99*Nhsc)
  #  LNmax=log(1.01*Nhsc)
  #}
  
  st=get_elapsed_time_tree(st)
  st$edge.length=st$edge.length/365
  list(st=st,nmutcolony=nmutcolony,nwtcolony=nwtcolony,node=node,acf=acf)
}

get_subsampled_tree_fixed_prop=function(selsim,N=100,nmut=10){
  mutants=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$driverid>0 & selsim$state==1)
  wt=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$driverid==0 & selsim$state==1)
  outgroup=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$state==0)
  tips=selsim$edge[c(sample(mutants,nmut,replace = FALSE),sample(wt,N-nmut,replace=FALSE),outgroup),2]
  get_subsampled_tree(selsim,-1,tips)
}

bench=function(N,S,sid,dpy=1,P=0.9999,OUTDIR="../benchmark_results/",maxSYear=10,M=5,nyd=10,mutperyear=20,nmut=20,bsavefull=FALSE){
  initSimPop(seed = as.integer(sid),bForce = TRUE)
  
  #nyd=5
  s=log(1+S)
  smax=log(1+maxSYear)
  #Following has expected ACF=0.1
  #set nyear such that 99% of the time the pure birthdeath would have a clone size=N.  
  ##. N ~ exp(st)*Exp(rate=s/(1+s))
  ## Solve cumulative distribution=p: 1-exp(rate*N*exp(-st))=1-p
  gett=function(s,N,p=0.99) -(1/s)*log(-log(p)/(N*(s/(1+s))))
  nyear=nyd+gett(s,N,p=P)
  top.label=sprintf("maxs_%d__divperyear_%3.2f__ageacquistion_%d__N_%3.0f__nmut_%d_P_%5.4f",round(log(maxSYear+1)),dpy,nyd,N,nmut,P)
  outdir=sprintf("%s/%s/S_%d",OUTDIR,top.label,S)
  dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
  label=sprintf("sim%d",sid)
  outfile=sprintf("%s/%s.RDS",outdir,label)
  cat("Saving results to",outfile,"\n")
  PHYLOFIT_NITER=10000
  N1=4
  N2=6
  z1=lapply(1:M,function(i){
    selsim=run_benchmark_sim(S=S,divrate_per_year = dpy,N = N,nyears_driver_acquisition = nyd,nyear=nyear,minprop = 0.5)
    
    st=get_subtree(selsim,nmut,nmut)
    justmutcladetree=drop.tip(st$st,st$st$tip.label[st$st$edge[st$st$state==0,2]])
    clonerateres=maxLikelihood(justmutcladetree)
    bdres=birthDeathMCMC(justmutcladetree,maxGrowthRate = smax,nCores = 4,chainLength = 4000)
    out=lapply(list(
      # with_acf_tightN_tightS=fit_clade(st$st,st$node,st$nmutcolony,st$nwtcolony,nchain=3,maxt=NA,minLN=N1,maxLN=N2,stan.control=list(adapt_delta=0.99,max_treedepth = 15),maxSYear = 1,mutperyear = 1000),
      no_acf=fit_clade(st$st,
                       node=st$node,
                       nmutcolony =-1,
                       nwtcolony =-1,
                       maxt=nyear,
                       minLN=N1,
                       maxLN=N2,
                       stan.control = list(adapt_delta=0.999,max_treedepth = 15),
                       maxSYear = maxSYear,
                       niter = 10000)    
      ),get_phylofit_summary,b.full=bsavefull,extra=list(acf=st$acf)
    )
    out$clonerateres=clonerateres
    out$bdres=bdres
    out$dt=nyear-nyd-max(nodeHeights(justmutcladetree))
    out$traj=data.frame(ts=selsim$timestamp,pop=selsim$pop.size,mutcellcount=selsim$totaldrivercount)
    out$trees=list(s=s,S=S,justmutcladetree=justmutcladetree,phylofit_tree=st$st[c("edge","tip.label","Nnode","edge.length")] %>% (function(x){class(x)="phylo";x}))
    out
  })
  saveRDS(z1,outfile)
  z1
}
