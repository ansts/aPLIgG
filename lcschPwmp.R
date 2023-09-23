lcschPwmp <- function(L, M=pwm) {
  require(parallel)
  require(pbapply)
  require(gtools)
  require(matrixStats)
  
  cl=makeCluster(10)
  clusterExport(cl, c("L", "aa", "M"),envir = environment())
  clusterEvalQ(cl,{
      require(matrixStats)
      require(gtools)
    })
  Res1=pbsapply(L[1:(length(L)-1)], function(x) { #
    xs=as.numeric(unlist(strsplit(x, split="")))
    n=max(xs)
    xsu=sort(unique(xs))
    xi=sapply(xsu,function(y) {
      j=which(xs==y)
      if (length(j)>1) {
        p=rowProds(M[,j])
        names(p)=rownames(M)
      } else p=M[,j]
      return(p)
    })
    if (n==1) return(sum(xi))
    colnames(xi)=unique(xs)
    ij=permutations(20,n)
    res=sum(apply(ij,1,function(prm){
        prod(xi[cbind(prm,1:n)])
    }))
    rm(ij)
    gc()
    return(res)
  }, cl=cl)
  stopCluster(cl)
  gc()
  print("Start the 7 cycle.............")
  x=L[length(L)]
    xs=as.numeric(unlist(strsplit(x, split="")))
    n=max(xs)
    xsu=sort(unique(xs))
    xi=sapply(xsu,function(y) {
      j=which(xs==y)
      if (length(j)>1) {
        p=rowProds(M[,j])
        names(p)=rownames(M)
      } else p=M[,j]
      return(p)
    })
    print("xi ready")
    colnames(xi)=unique(xs)
    #N=choose(20,n)*factorial(n)
    ij=permutations(20,n)
    print("permutaion list calculated")
    
    nsl=80
    msl=15
    J=splitIndices(nrow(ij), nsl)
    print("indices split")
    Res2=0
    for (j in 1:nsl){
      print(j)
      X=ij[J[[j]],, drop=F]
      J1=splitIndices(nrow(X), msl)
      cl=makeCluster(15)
      clusterExport(cl, c("n","xi","X","J1"),envir = environment())
      clusterEvalQ(cl, require(matrixStats))
      rs=sum(pbsapply(1:length(J1), function(cj){
        sum(rowProds(t(apply(X[J1[[cj]],, drop=F],1,function(l) xi[cbind(l,1:n)]))))
        },cl=cl))
      rm(X,J1)
      closeAllConnections()
      gc()
      Res2=Res2+rs
    }

    rm(ij)
    closeAllConnections()
  
  return(c(Res1,Res2))
}
