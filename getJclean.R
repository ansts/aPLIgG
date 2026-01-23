getJclean=function(path0, path, f="vdjserver.tsv", batch=5000000){
  require(parallel)
  require(pbapply)
  require(future.apply)
  require(Biostrings)
  
  aa=AA_STANDARD
  
  n=1;X=c()    
  fn=paste(path0, path, f, sep="")
  repeat {
    print(n)
    cx=read.delim(fn,nrows = batch, skip = (n-1)*batch, header=(n==1))
    if (n==1) clnm=colnames(cx) else colnames(cx)=clnm
    x=cx[,c("junction_aa","repertoire_id")]

      x=x[nchar(x[,1])>6,]
      x[,1]=substring(x[,1],2)
      i=grep("[^[ARNDCQEGHILKMFPSTWYV]",x[,1])
      x=x[-i,]
      rownames(x)=NULL
      y=x[,2];names(y)=x[,1];x=y;rm(y)

    X=c(X,x)  
    print(length(X))
    pres=nrow(cx)
    rm(cx)
    gc()
    n=n+1
    if (pres<batch) break
  }
  rep=unique(X)
  X=unlist(lapply(rep, function(rx) unique(names(X)[X==rx])))
  return(X)
}