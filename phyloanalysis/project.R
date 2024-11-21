##
SPECIES="Hsapiens" # Hsapiens  Mmusculus
SPECIES_BUILD="hg38" # hg19 hg38 mm10
patients=read.table("../data/patients.txt",head=T)
PATIENTS=patients$PDID
CV=c("missense","nonsense","ess_splice","frameshift","inframe","loh","start_lost","cna","stop_lost")
GENES=readLines("GENES.txt")
MUTCOUNTBIN=sprintf("../cache/mutcountbin.%s.RDS",SPECIES_BUILD)

get_driver_scheme2=function(){
  driver.scheme=read.table("driver_scheme_simple.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  n=max(driver.scheme$number)
  pallete=c(RColorBrewer::brewer.pal(9,"Set1")[-6],RColorBrewer::brewer.pal(8,"Dark2"))
  driver.scheme$colour=pallete[driver.scheme$number]
  driver.scheme
}

##Import age of sample data
##Add all relevant meta
##  Could also add columns to CFG here for individul specific info

add_agedf=function(PD){
  tree=PD$pdx$tree_ml
  agefile="../data/sample_ages.txt"
  if(!file.exists(agefile)){
    cat("can't find",agefile,"\n")
    cat("need to creat tab delimited file ",agefile," with column headings sample and age_at_sample.  Note that sample column should contain sample long names")
    stop("unable to find age file ")
  }
  sample=PD$pdx$cfg$LABEL[match(PD$pdx$tree_ml$tip.label,PD$pdx$cfg$SHORT_LABEL)]
  if(length(which(is.na(sample)==1))){
    sample[is.na(sample)]="zeros"
  }else{
    stop("unexpected lookup error!")
  }
  
  agedf=data.frame(tip.label=tree$tip.label,sample=sample)
  af=read.table(agefile,stringsAsFactors=FALSE,sep="\t",head=TRUE)
  agedf=agedf %>% left_join(af) %>% dplyr::rename(age_at_sample_exact=age_at_sample)
  agedf$age_at_sample_pcy=agedf$age_at_sample_exact+(MEAN_AGE_AT_DELIVERY_DAYS)/365.25
  agedf$age_at_sample_pcy[which(agedf$sample=="zeros")]=1e-6
  PD$pdx$agedf=agedf
  PD
}

### CML specific nodes.


do_donor_specific_filtering_and_annotation=function(PD){
  if(PD$patient=="PD51634"){
    idx=length(PD$pdx$meta$CNA)+1
    PD$pdx$meta$CNA[[idx]]=list(LABEL="BCR::ABL1",chr=9,start=0,end=1,ploidy=2,samples=get_samples_in_clade(131,PD$pdx$tree_ml))
    if(length(PD$pdx$meta$CNA[[idx]]$samples)!=75){
      browser()
      stop("node was ascertained on a different version of the tree!")
    }
  }
  if(PD$patient=="PD51635"){
    idx=length(PD$pdx$meta$CNA)+1
    PD$pdx$meta$CNA[[idx]]=list(LABEL="BCR::ABL1",chr=9,start=0,end=1,ploidy=2,samples=get_samples_in_clade(183,PD$pdx$tree_ml))
    if(length(PD$pdx$meta$CNA[[idx]]$samples)!=32){
      stop("node was ascertained on a different version of the tree!")
    }
  }
  if(PD$patient=="PD57334"){
    ##browser()
    PD$pdx$meta$CNA[[1]]$LABEL="BCR::ABL1"
  }
  if(PD$patient=="PD57332"){
    PD$pdx$meta$CNA[[1]]$LABEL="BCR::ABL1"
  }
  if(PD$patient=="PD56961"){
    idx=length(PD$pdx$meta$CNA)+1
    PD$pdx$meta$CNA[[idx]]=list(LABEL="BCR::ABL1",chr=9,start=0,end=1,ploidy=2,samples=get_samples_in_clade(244,PD$pdx$tree_ml))
    if(length(PD$pdx$meta$CNA[[idx]]$samples)!=58){
      browser()
      stop("node was ascertained on a different version of the tree!")
    }
  }
  if(PD$patient=="PD57335"){
    idx=length(PD$pdx$meta$CNA)+1
    PD$pdx$meta$CNA[[idx]]=list(LABEL="BCR::ABL1",chr=9,start=0,end=1,ploidy=2,samples=get_samples_in_clade(109,PD$pdx$tree_ml))
    if(length(PD$pdx$meta$CNA[[idx]]$samples)!=91){
      browser()
      stop("node was ascertained on a different version of the tree!")
    }
  }
  if(PD$patient=="PD51633"){
    for(sample in c("lo0120","lo0127")){
      PD=fixContaminatedSample(PD,sample)
    }
  }
  return(PD)
}

fix_trunk_PD57332=function(PD){
  if(PD$patient!="PD57332"){
    ## Keep germline variants as stub
    return(PD)
  }
  #  ##PD=get_pd("PD57332",b.remove.germline = TRUE)
  #tree=plot_basic_tree(PD$pdx,PD$patient,genes = GENES,cv=CV,cex.annot.label = 1,cex.label=0)
  tree=PD$pdx$tree_ml
  nh=nodeHeights(tree)
  #
  mean.BCRABL.clade.burden=mean(nh[which(tree$edge[,2]<length(tree$tip.label)),2])
  #bcadat %>% filter(Donor!="PD57335") %>% (function(x) mean(x$BCRABL_mean))
  # Mean BCRABL rate   
  mean.BCRABL.rate=26.10834# See above
  mean.WT.rate=17.78061 # See above
  additional.muts.at.birth=35 ## Additional mutates above that expected for 9 months accumulaton.
  
  ##warning("CHECK cohort mean WT and BCRABLE rates are correct")
  BCRABL.duration=mean.BCRABL.clade.burden/mean.BCRABL.rate
  
  burden=round(additional.muts.at.birth+mean.WT.rate*(PD$pdx$agedf$age_at_sample_pcy[2]-BCRABL.duration))
  idx=get_germline_idx(PD$pdx)
  for( field in c("edge.length","el.snv","el.snv.local.filtered")){
    PD$pdx$tree_ml[[field]][idx]=burden
  }
  sens.fields= names(PD$pdx$tree_ml)
  sens.fields=sens.fields[grepl("sens",sens.fields)]
  for( field in sens.fields){
    PD$pdx$tree_ml[[field]][idx]=0.999
  }
  PD$trunk_params=list(mean.BCRABL.rate=mean.BCRABL.rate,
                       mean.WT.rate=mean.WT.rate,
                       additional.muts.at.birth=additional.muts.at.birth)
  PD
}

fixContaminatedSample=function(PD,sample){
  ## get parents and estimate VAF based on shared branches.
  ## Then assume that assigned variants are due to 2 clones VAF and 0.5-VAF.
  tree=PD$pdx$tree_ml
  tip.node=match(sample,tree$tip.label)
  parents=get_parents(tip.node,tree$edge)
  parents=setdiff(parents,tip.node)
  idx=with(PD$pdx$dat$details,which(node %in% parents & TYPE %in% "SNV" & is_localx_excluded==0))
  depth=PD$pdx$dat$dep[idx,sample]
  mtr=PD$pdx$dat$mtr[idx,sample]
  vaf=sum(mtr)/sum(depth)
  cat("estimated VAF for",sample,"is",vaf,"\n")
  #browser()
  ## Now identify those that have been assigned to the private branch and estimate the proportion (or number "n") that should
  ## be retained and then assign the n most probable mutations - the others are then removed.  
  idx=with(PD$pdx$dat$details,which(node == tip.node & TYPE %in% "SNV" & is_localx_excluded==0))
  private.branch=data.frame(mtr=PD$pdx$dat$mtr[idx,sample],dep=PD$pdx$dat$dep[idx,sample],IDX=idx)
  private.branch$p2=sapply(1:dim(private.branch)[1],function(i) dbinom(private.branch$mtr[i],private.branch$dep[i],0.5-vaf))
  private.branch$p1=sapply(1:dim(private.branch)[1],function(i) dbinom(private.branch$mtr[i],private.branch$dep[i],vaf))
  private.branch=private.branch %>% mutate(p=p1/(p1+p2))
  N=length(idx)
  #browser()
  p=sum(private.branch$p)
  n=floor(p)
  cat("Retaining a proportion",
      p/N,
      "(n=",
      n,
      ") of private",
      sample,
      "variants\n"
  )
  ## Check original SNV 
  remove.private.branch=head(private.branch[order(private.branch$p),],N-n)
  ## Update el.
  orig.snv.el=N #with(PD$pdx$dat$details[private.branch$IDX,],length(which(TYPE=="SNV" & is_localx_excluded==0)))
  eidx=match(tip.node,tree$edge[,2])
  if(!is.null(tree$per.branch.sensitivity.reg)){
    stop("Run before add_adjustment_models")
  }
  cat("Original private SNV EL=",tree$el.snv.local.filtered[eidx], " chk=",orig.snv.el)
  idx=remove.private.branch$IDX
  if(length(idx)>0){
    cat(PD$patient,":Removing",length(idx)," variants from",sample,"..\n")
    PD$pdx$summary=PD$pdx$summary[-idx,]
    N=dim(PD$pdx$dat$details)[1]
    for(x in names(PD$pdx$dat)){
      if(dim(PD$pdx$dat[[x]])[1]!=N){
        stop("Unexpected dimension pdx$dat")
      }
      PD$pdx$dat[[x]]=PD$pdx$dat[[x]][-idx,]
    }
  }
  idx=with(PD$pdx$dat$details,which(node == tip.node & TYPE %in% "SNV" & is_localx_excluded==0))
  cat(sample,":New SNV EL=",length(idx)," for \n")
  PD$pdx$tree_ml$el.snv.local.filtered[eidx]=length(idx)
  PD
}

get_PD=function(patient,b.force.reload=FALSE,...){
  cfile=sprintf("%s/%s.RDS",CACHE,patient)
  cat("looking for",cfile,"\n")
  if(file.exists(cfile) && b.force.reload==FALSE){
    cat("reading cached data")
    return(readRDS(cfile))
  }else{
    PD=get_pd(patient,...)
    if(patient=="PD57333"){
      PD$gmiss.correction=get_gmiss_correction(PD)
      PD$localx.correction2=PD$localx.correction2*PD$gmiss.correction
    }else{
      ## Minimal correction for the others
    }
    
    PD=add_adjustment_models(PD)
    
    PD=fix_trunk_PD57332(PD)
    cat("caching data")
    PD$nodes=PD$nodes %>% dplyr::select(node,driver,status,driver2,driver3,child_count)
    saveRDS(PD,cfile)
    
    return(PD)
  }
}

get_gmiss_correction=function(PD,filter.ex.field="filter_max_gmiss"){
  fields=colnames(PD$inf$snv$details)
  filters_ex_gmiss=setdiff(fields[grepl("^filter_",fields)],filter.ex.field)
  ##cat("Filters ex gmiss are",paste(filters_ex_gmiss,collapse=","),"\n")
  df=PD$inf$snv$details
  df$FILTER_EX_GMISS=ifelse(rowSums(df[,filters_ex_gmiss],na.rm=TRUE)>0,1,0)
 # browser()
  sum(df$FILTER_EX_GMISS==0)/sum(df$FILTER_EX_GMISS==0 & df[[filter.ex.field]]==0)
}
