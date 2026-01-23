diamap=function(tL){
  require(Biostrings)
  require(igraph)
  require(parallel)
  require(future.apply)
  require(pbapply)
  require(reshape2)
  require(matrixStats)
  require(Matrix)
  
  L=names(tL)
  aa=AA_STANDARD; n=unique(nchar(L))
  if (length(n)>1) stop("The sequences should be of equal length!")
  x=c(sapply(aa,\(a1) sapply(aa,\(a2) paste(a1,a2,sep = ""))))
  x=paste(rep(x,each=4),rep(0:(n-2),length(x)),sep="")
  dpp=rep(0,length(x)); names(dpp)=x
  
  # ncores=20
  # plan(multisession, workers=ncores)
  dienc=future_sapply(L,displit, future.chunk.size=length(L)%/%20)

  diencsp=apply(dienc,1,\(x) {
    x=x[x>0]
    x=x*tL[names(x)]
    return(x)
  })
 
  AM=future_sapply(diencsp,\(s1){
       sapply(diencsp,\(s2){
         il=intersect(names(s1), names(s2))
         if (length(il)>0) return(sum(s1[il]*s2[il])) else return(0)
       })
  }, future.chunk.size=length(diencsp) %/% ncores)  #
  #plan(sequential)
  
  diag(AM)=0
  AM=AM/sum(AM)

  return(AM)
}

displit=function(p,dpp_=dpp){

  ij=combn(nchar(p),2)
  dij=colDiffs(ij)-1
  sp=unlist(strsplit(p,split=""))
  dsp=sapply(seq_along(dij),\(i) {
      paste(paste(sp[ij[,i]], collapse=""),dij[i],sep="")
  })
  dummy=sapply(dsp, \(dp) {
      dpp_[dp]<<-dpp_[dp]+1
  })

  return(dpp_)
}
