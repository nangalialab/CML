require("rstan")
### Inference stuff
PHYLOFIT_NITER=10000
phylofit.stan="functions{

real glogistic_lpdf(real[] t, real s, real tm, real T,real N,real offset){
int n;
real ll;
n=size(t);

ll=0.0;
for(k in 2:n){
  ll=ll+log(1+exp(s*(t[k-1]-offset-T+tm)));
  ll=ll-(1/(s*N))*choose(k,2)*(exp(-s*(T-tm-t[k]+offset)+s*(t[k-1]-t[k]))-exp(-s*(T-tm-t[k]+offset)));
  ll=ll+choose(k,2)*(t[k]-t[k-1])/N+log(choose(k,2))-log(N);
}
//print(\"LL=\",ll,\"s=\",s,\";offset=\",offset,\";N=\",N,\";tm=\",tm,\";T=\",T);
return(ll);
}

}

data{
int N;      //num of coalescences
real t[N];  //timing of coalescences
real T;
real maxT;
real minT;
real maxLN;
real minLN;
real mutPerYear;
real smax;
int nmutcolony; //Could be number of mutant colonies or reads.
int nwtcolony; //Could be number of wild type colonies or reads.
}

parameters {
real<lower=0.0001,upper=smax> s; //instantaneous growth rate-> exp(s)-1 per Year. log(2)  0.6931472 
real <lower=minT,upper=maxT> tm; //midpoint
real<lower=minLN,upper=maxLN> LN; //log pop size
}

model {
  s ~ uniform(0.001,smax);
  LN ~ uniform(minLN,maxLN);
  tm ~ uniform(minT,maxT);
  t ~ glogistic(s,tm,T,pow(10,LN),0);//offset);
  if(nmutcolony>0)
    nmutcolony ~ binomial(nmutcolony+nwtcolony,1/(1+exp(-s*(T-tm))));
}
generated quantities {
real S;
S=exp(s)-1;
}
"
modelfile="../data/PHYLOFIT_STAN_MODEL.RDS"
if(file.exists(modelfile)){
  PHYLOFIT_STAN_MODEL=readRDS(modelfile)
}else{
  PHYLOFIT_STAN_MODEL=stan_model(model_code=phylofit.stan,model_name = "phylofit")
  saveRDS(PHYLOFIT_STAN_MODEL,modelfile)
}

get_phylologistic_dat=function(ultratree,node,maxt=NA,mutperyear=20){
  nc=get_all_node_children(node,ultratree)
  idx=match(nc,ultratree$edge[,2])
  nh=nodeHeights(ultratree)
  ages=nh[match(match(get_samples_in_clade(node,ultratree),ultratree$tip.label),ultratree$edge[,2]),2]
  T=max(ages)
  clade=keep.tip(ultratree,get_samples_in_clade(node=node,tree = ultratree))
  ##
  if(!is.ultrametric(clade)){
    stop("mutant clade must be ultrametric..")
  }
  if(!is.binary(clade)){
    clade=multi2di(clade);
  }
  tc=c(sort(branching.times(clade),decreasing = TRUE),0)
  if(is.na(maxt)){
    maxt=2*T
  }
  list(N=length(tc),t=tc,T=T,maxT=maxt,mutPerYear=mutperyear)
}

##Fits the clade with the specified ancestral branch (node)
fit_clade=function(ultratree,node,nmutcolony,nwtcolony,nchain=4,maxt=NA,minLN=4,maxLN=6,mutperyear=20,maxSYear=1,niter=PHYLOFIT_NITER,stan.control=list(adapt_delta=0.95),...){
  dat=get_phylologistic_dat(ultratree,node,maxt = maxt,mutperyear = mutperyear )
  dat$nmutcolony=nmutcolony
  dat$nwtcolony=nwtcolony
  dat$smax=log(1+maxSYear)
  dat$smin=0.001
  if(nmutcolony>=0){
    tvaf=nmutcolony/(nmutcolony+nwtcolony)
    ##Get reasonable range for mid-point prior (tm).  
    #Firstly get range of likely true ACF based on colony sampling.
    ci=binom.test(nmutcolony,nmutcolony+nwtcolony,conf.level = 0.999)$conf.int
    #Assuming a maximum instantaneous growth rate (s) of 2 and a minimum of 0.05.   
    #Then solve ACF=1/(1+exp(-s(t-tm)))
    #Giving tm=(1/s)*log((1/ACF)-1)+t
    #The following considers the 4 extremal possibilities of minACF,maxACF x minS,maxS 
    #and uses the range as the bounds for the uniform prior of tm. 
    vlog1=log((1/ci[2])-1)
    vlog2=log((1/ci[1])-1)
    ismin=1/0.05
    ismax=1/2
    tmrange=c(ismin*vlog1,ismax*vlog1,ismin*vlog2,ismin*vlog2)
    mint=max(min(tmrange)+dat$T,0)
    maxt=max(max(tmrange)+dat$T,10)
    if(mint>maxt){
      stop("Inconsistency in prior range for tm!")
    }
    if(dat$T-dat$t[1]>mint){
      mint=mint
    }
    dat$maxT=maxt
    dat$minT=mint
    cat("maxt=",maxt,"mint=",mint,"\n")
  }else{
    dat$minT=dat$T-dat$t[1]
  }
  dat$minLN=minLN
  dat$maxLN=maxLN
  ## When initialisied too high the markov chain tends to get stuck. Initialise towards the start range.
  init.s=lapply(1:nchain,function(i) list(s=runif(1,
                                                  min=dat$smin,
                                                  max=min(dat$smin+0.2*(dat$smax-dat$smin),8)
  )
  )
  )
  if(FALSE){
  cat("initialising s as follows:")
  print(init.s)
  }
  
  stanr=sampling(PHYLOFIT_STAN_MODEL,
                 data=dat,iter = niter,control=stan.control,chains = nchain,cores = nchain,init=init.s,...)
  list(posterior=rstan::extract(stanr),
       res=stanr,
       ultratree=ultratree,
       dat=dat
  )
}

get_phylofit_summary=function(zzz,b.full=FALSE,extra=NULL){
  if(b.full){
    res=zzz$res
  }else{
    res=NULL
  }
  list(
    S=quantile(zzz$posterior$S,prob=c(0.025,0.25,0.5,0.75,0.975)),
    Smean=mean(zzz$posterior$S),
    LN=quantile(zzz$posterior$LN,prob=c(0.025,0.25,0.5,0.75,0.975)),
    LNmean=mean(zzz$posterior$LN),
    tm=quantile(zzz$posterior$tm,prob=c(0.025,0.25,0.5,0.75,0.975)),
    dat=zzz$dat,
    ndivt=get_num_divergent(zzz$res),
    lowbfmichains=get_low_bfmi_chains(zzz$res),
    rhatmax=max(summary(zzz$res)$summary[,"Rhat"]),
    res=res,
    extra=extra
  )
}