dAA=function(tL,l=5){
  require(Biostrings)
  require(parallel)
  require(future)
  require(future.apply)
  
  a=AA_STANDARD
  aad=array(0, dim=c(20,20,4), dimnames=list(a,a,c(1:4)))
  ij=t(combn(l,2))
  ij=cbind(ij,ij[,2]-ij[,1])
  ij=ij[ij[,3]<=l,]
  L=names(tL)
  L=strsplit(L,split="")
  #cl=makeCluster(12, timeout = 60)
  plan(multisession, workers=12)
  taa=rowSums(future_vapply(seq_along(L), function(i){
      p=L[[i]]
      n=tL[i]
      ij12=aggregate(rep(1,nrow(ij)), by=list(p[ij[,1]],p[ij[,2]],ij[,3]), "sum")
      aad[as.matrix(ij12[,1:3])]=n*ij12$x
      return(aad)
  }, aad), dims=3)
  closeAllConnections()
  return(taa)
}