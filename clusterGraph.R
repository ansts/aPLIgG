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

clusterGraph <- function(G, resol=30, fnm=NULL, nodes=14) {
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
  print("Leiden clustering.........")
  Glouv=cluster_leiden(G, objective_function = "modularity", resolution_parameter = resol, n_iterations = 5)
  n=Glouv$names
  #d=nrow(Glouv$memberships)
  #if (d>1) d=d-1
  Glouv=Glouv$membership
  names(Glouv)=n
  print(paste(c(max(Glouv)," clusters and " ,sum(table(Glouv)==1), " singlets"), sep="", collapse=""))
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
  vertex_attr(Gctr)$name=unlist(vertex_attr(Gctr)$name)
  Gctr=set_vertex_attr(Gctr, "size", value=lengths(Gpeps))
  Gctr=set_vertex_attr(Gctr, "edges", value=sapply(Gpeps,function(x){
    g=induced_subgraph(G,x)
    length(E(g))
  }))
  rm(G)
  gc()
  print("Calculating distances.........")

  n=length(V(Gctr))
  A=matrix(0,n,n)
  rownames(A)=names(V(Gctr))
  colnames(A)=names(V(Gctr))
  eG=apply(ends(Gctr, E(Gctr), E(Gctr)),2,unlist) # matrix of edge ends
  w=edge.attributes(Gctr)$weight                  # number of edges between clusters
  A[eG]=w                                         # adjacency mx of clusters
  A=t(A)
  A[eG]=w
  Ein=vertex_attr(Gctr)$edges                     # number of edges in each cluster
  S=vertex_attr(Gctr)$size                        # number of vertices in each cluster
  ex=(Ein+1*(S==1))/(S*(S-1)/2+1*(S==1))
  ex=matrix(rep(ex,n),n,n)
  m=(ex + t(ex))/2                                # mean density
  w=(S%*%t(S))*m/(A+1*(A==0))                     # expected no of edges between/ actual
  rownames(w)=names(V(Gctr))                      # w is a dissimilarity measure
  colnames(w)=names(V(Gctr))
  w=w[eG]

  w=log10(w+min(w[w>0])*(w==0)) 
  x=w[order(w)]
  x=x[length(x)]-x[length(x)-1]
  w=w-min(w)+x/2
  
  Gctr=set.edge.attribute(Gctr, name="weight", value=w)
  save(Gctr,file=paste("Gctr",fnm,sep=""))
  save(Gpeps,file="Gctrpeps")
  
  print("Prune large distances ...")
  thrng=range(edge.attributes(Gctr)$weight)
  n=length(V(Gctr))
  thrscan=t(sapply(seq(thrng[1],thrng[2],length.out=200), function(th){
    Gctrsm=delete.edges(Gctr,E(Gctr)[E(Gctr)$weight>th])
    cmp=components(Gctrsm)
    c(th,cmp$no,sum(cmp$csize<5), graph.density(Gctrsm))
  }))
  # x=mst(Gctr)
  # thr=max(E(x)$weight)
  thr=min(thrscan[thrscan[,3]<round(0.05*n),1])
  Gctrsm=delete.edges(Gctr,E(Gctr)[E(Gctr)$weight>thr])
  #Gctrsm=set.edge.attribute(Gctrsm,"weight", value=x)
  #Gctrsm=set.vertex.attribute(Gctrsm,"size",value=3*log(vertex_attr(Gctrsm)$size+0.5,2))
  save(Gctrsm,file=paste("Gctrsm",fnm,sep=""))
 
  write_graph(Gctrsm,format = "graphml", file=paste("Gctrsm",fnm, ".graphml",sep=""))
 
  return(Gctrsm)
}