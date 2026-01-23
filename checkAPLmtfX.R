checkAPLmtfX=function(m){
  require(Biostrings)
  aa=AA_STANDARD
  l=length(m)
  f=c("A", "V", "I", "L", "M", "F", "W", "C", "P", "G")
  z=c("Y", "T", "S", "H", "K", "R", "E", "D", "Q", "N")
  ff="F"
  ai=c("f","f","z","z","ff") #c("f","f","f","f","z","z","ff","aa","f")
  lm0=length(ai)
  m0=matrix(0,20,lm0)
  rownames(m0)=aa
  for (i in 1:lm0) {
    m0[aa %in% get(ai[i]),i]=1/length(get(ai[i]))
  }
  m0=m0+0.05/10; m0=apply(m0,2,\(cl) cl/sum(cl))
  m=consensusMatrix(AAStringSet(m), as.prob = T)[aa,]
  n=ncol(m)
  m=m+0.05/l; m=apply(m,2,\(cl) cl/sum(cl))
  mx=t(m) %*% m0
  i=(2-n):lm0
  max(sapply(i,\(ii){
    prod(sapply(1:n,\(jj){
      if((ii+jj-1)>0 & (ii+jj-1)<(lm0+1)) mx[jj,ii+jj-1] else 0.05
    }))
  }))
}