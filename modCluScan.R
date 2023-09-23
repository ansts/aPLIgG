modCluScan <- function(G, ii, niter=-1, N=20, ncore=4) {
  require(igraph)
  require(parallel)
  require(pbapply)
  
  vG=names(V(G))
  nlim=length(vG)*0.01
  cl=makeCluster(ncore)
  clusterExport(cl, c("G","vG"), envir = environment())
  clusterEvalQ(cl, require(igraph))
  clusterEvalQ(cl, require(stringi))
  modsc=pbsapply(ii, function(resol){
    print(resol)
    X=sapply(1:3,function(i){
      x=cluster_leiden(G, objective_function = "modularity", resolution_parameter = resol, n_iterations = niter)
      x=x$membership
      list(modularity(G,as.numeric(as.factor(x))), x)
    })
    Xm=unlist(X[1,])
    Xmm=X[2,]
    m=median(Xm)
    i=which(Xm==m)
    n=length(table(Xmm[[i]])[table(Xmm[[i]])>1])
    BSmod=c()
    for (ij in 1:N) {
      labv=rep("_",vcount(G))
      names(labv)=vG
      seed=sample(names(V(G)),n)
      seedn=paste("m",1:n,sep="")
      labv[seed]=seedn
      
      repeat{
        #proct=proc.time()
        #print(length(seed)/vcount(G))
        av=adjacent_vertices(G,V(G)[seed])
        v=sample(seedn)
        for (i in v){
          ss=names(labv)[labv==i]
          jj=which(names(labv) %in% unique(unlist(stri_extract_all(names(unlist(av[ss])), regex="(?<=\\.)\\w+"))) & labv=="_")
          labv[jj]=i
        }
        #print(table(labv))
        seed=V(G)[vG %in% names(labv)[labv!="_"] & !(vG %in% names(seed))] 
        if (sum(labv=="_")<nlim) break
        #print(proc.time()-proct)
      }
      for (l in names(labv[labv=="_"])){
        tx=table(labv[names(unlist(ego(G,1,V(G)[names(labv[labv=="_"])])))])
        tx=tx[names(tx)!="_"]
        labv[l]=names(tx)[which.max(tx)]
      }
      mo=modularity(G,as.numeric(as.factor(labv[vG])))
      BSmod=c(BSmod,mo)
    }
    res=m/mean(BSmod)
    rm(m,BSmod,av,labv,seed,seedn,x,n)
    gc()
    return(res)
  }, cl=cl)
  closeAllConnections()
  return(cbind(ii,modsc))
}