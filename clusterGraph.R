# Takes a sequence graph and clusters it to about 2000 clusters 
# and contracts it to a graph of the clusters which is returned. 
# The number of clusters is controlled by the resolution parameter
# which has to be established in advance (default=32).
# The edges are summed and the Group attribute is returned as a 
# table of the values found in the clusters. The edge weight is 
# the log(modularity) of the subgraph formed of the pair of 
# clusters of the initial graph connected in the graph of 
# clusters. Outputs a graphml file with the luster graph
# and returns the cluster graph as an igraph object. Size is
# the size of the clusters in the initial graph.

clusterGraph <- function(G, resol=32, fnm=NULL, nodes=14) {
  require(igraph)
  require(parallel)
  require(pbapply)
  require(numbers)
  
  dGm=degree(G)
  dGmi=cut(dGm, 50, labels = F)
  dGmh=aggregate(dGm, by=list(dGmi), "length")
  
  pdf(file=paste("DegreeDistr",fnm,".pdf",sep=""),width=5, height=5)
  plot(dGmh, log="y", main=paste("Degree Distribution of Gctr",fnm,sep=""))
  dev.off()
  
  if (is.null(fnm)) fnm=round(runif(1,1,1e4),0)
  print("Louvain clustering.........")
  Glouv=cluster_louvain(G, resolution = resol)
  n=Glouv$names
  Glouv=Glouv$memberships[2,]
  names(Glouv)=n
  save(Glouv,file=paste("Glouv",fnm,sep="")) 
  Grpsall=table(V(G)$Group)  
  gnm=names(Grpsall)
  x0=rep(0,length(gnm))
  names(x0)=gnm
  print("Contraction.........")
  Gctr=contract(G,mapping=Glouv,vertex.attr.comb=list(name="concat", Group=function(x) {
    x3=table(x)
    x0[names(x3)]=x3
    return(x0)
  }))
  Gctr= simplify(Gctr, edge.attr.comb = "sum")
  clnm=paste("C",seq_along(V(Gctr)), sep="")
  Gpeps=vertex_attr(Gctr)$name
  names(Gpeps)=clnm

  Gctr=set_vertex_attr(Gctr, "name", value=clnm)
  Gctr=set_vertex_attr(Gctr, "size", value=lengths(Gpeps))
  Gctr=set_vertex_attr(Gctr, "edges", value=sapply(Gpeps,function(x){
    g=induced_subgraph(G,x)
    length(E(g))
  }))
  rm(G)
  gc()
  print("Calculating distances.........")
  
  n=length(E(Gctr))
  ei=seq_along(E(Gctr))
  nmax=max(sapply(14:1, function(i) GCD(n,i)))
  if (nmax==1) nmax=14
  ni=nmax *10
  ij=cut(ei,ni,labels=F)
  
  cl=makeCluster(nmax)
  clusterExport(cl, c("Gctr","ij"), envir = environment())
  clusterEvalQ(cl, require(igraph))
  
  x=pbsapply(1:ni, function(j){ # using a log(1/estimate) of the modularity as a dissimilarity measure
    sapply(seq_along(E(Gctr))[ij==j], function(i){
      es=ends(Gctr,i, names=F)   # indices of incident vertices
      ebtw=E(Gctr)$weight[i]     # number of between edges (weight of edge after contraction)
      s=vertex_attr(Gctr)$size[es[1,]]               # sizes of the two respective clusters in G
      allein=sum(c(vertex_attr(Gctr)$edges[es[1,]])) # No of edges inside both clusters
      d=(s[1]*(s[1]-1)/2+s[2]*(s[2]-1)/2)
      if (d <= 0) d=1
      pein=allein/d              # graph density inside clusters
      if (pein<=0) pein=0.01
      pebtw=ebtw/prod(s)         # edge density between clusters 
      log10(pein/pebtw)          # log of the ratio of within/between edge density
    })
  },cl=cl)
  
  stopCluster(cl)
  
  w=x[[1]]
  for (i in 2:ni){
    w=c(w,x[[i]])
  }
  
  x=w[order(w)]
  x=x[length(x)]-x[length(x)-1]
  w=w-min(w)+x/2
  
  Gctr=set.edge.attribute(Gctr, name="weight", value=w)
  save(Gctr,file=paste("Gctr",fnm,sep=""))
  
  print("Prune large distances down to the smallest connected graph...")
  x=mst(Gctr)
  thr=max(E(x)$weight)
  Gctrsm=delete.edges(Gctr,E(Gctr)[E(Gctr)$weight>thr])
  Gctrsm=set.vertex.attribute(Gctrsm,"size",value=3*log(vertex_attr(Gctrsm)$size+0.5,2))
  save(Gctrsm,file=paste("Gctrsm",fnm,sep=""))
 
  write_graph(Gctrsm,format = "graphml", file=paste("Gctrsm",fnm, ".graphml",sep=""))
 
  return(Gctr)
}