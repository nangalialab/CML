
get_lineages=function(tree,threshold,min.count=2){
  nh=nodeHeights(tree)
  idx=which(nh[,1]<threshold & nh[,2]>=threshold)
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>1)
}

get_lineages2=function(tree,nodes){
  nh=nodeHeights(tree)
  idx=match(nodes,tree$edge[,2])
  if(any(is.na(idx))){
    browser()
    stop("bad node list")
  }
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>1)
}
get_lineages3=function(tree,threshold,min.count=2){
  nh=nodeHeights(tree)
  idx=which(nh[,1]<threshold & nh[,2]>=threshold)
  nodes=tree$edge[idx,2]
  out=data.frame(node=nodes,start=nh[idx,1],end=nh[idx,2])
  if(length(nodes)>0){
    out$N=sapply(nodes,function(node) length(get_samples_in_clade(node,tree)))
  }else{
    out$N=integer()
  }
  out %>% filter(N>=min.count)
}

do_tree_branch_length_ci_etc=function(mtree,treefitres,lineages){
  #lineages=get_lineages(mtree,threshold=threshold,min.count=min.count )
  tree=get_treefit_ci(treefitres)
  #browser()
  nh=nodeHeights(tree)
  idx=match(lineages$node,tree$edge[,2])
  ##Wild type = non-expanded 
  mt.clones=do.call("c",
                    lapply(lineages$node,
                           function(node) get_samples_in_clade(node,tree
                           )
                    )
  )
  wt.clones=setdiff(tree$tip.label,c("zeros",mt.clones))
  nh[,2]=round(nh[,2],6)
  ts=unique(nh[match(1:length(tree$tip.label),tree$edge[,2]),2],digits = 6)
  ts=ts[which(ts>1e-3)]
  cat("Found",length(ts),"timepoints.. just using the first time point - assumed at diagnosis")
  ts=min(ts)
  #browser()
  tmp=data.frame(samples=wt.clones,age=nh[match(match(wt.clones,tree$tip.label),tree$edge[,2]),2])
  nwtcolony=sapply(ts,function(t) length(which(tmp$age==t)))
  
  lineages$lower=tree$upper_lb95[idx]
  lineages$upper=tree$upper_ub95[idx]
  for(field in c("lower_lb95","lower_ub95","lower_median","upper_lb95","upper_median","upper_ub95")){
    lineages[[sprintf("t_%s",field)]]=tree[[field]][idx]
  }
  lineages$t_diagnosis=ts
  return(list(lineages=lineages))
}

get_treefit_ci=function(treefit){
  tree=treefit$ultratree
  tmp = tree
  post.dist=rstan::extract(treefit$fullres)
  nodetdist = matrix(apply(post.dist$ta, 1, function(x) {
    tmp$edge.length = x
    nodeHeights(tmp)[,1]
  }), nrow = length(tmp$edge.length))
  parent_bounds = apply(nodetdist, 1, function(x) quantile(x, prob = c(0.025, 
                                                              0.5, 0.975), na.rm = TRUE))
  nodetdist = matrix(apply(post.dist$ta, 1, function(x) {
    tmp$edge.length = x
    nodeHeights(tmp)[,2]
  }), nrow = length(tmp$edge.length))
  child_bounds = apply(nodetdist, 1, function(x) quantile(x, prob = c(0.025, 
                                                                       0.5, 0.975), na.rm = TRUE))
  
  tree$lower_lb95=parent_bounds[1,]
  tree$lower_median=parent_bounds[2,]
  tree$lower_ub95=parent_bounds[3,]
  tree$upper_lb95=child_bounds[1,]
  tree$upper_median=child_bounds[2,]
  tree$upper_ub95=child_bounds[3,]
  tree
}




mlog2=function(splus1){
  log(2)/log(splus1)
}

