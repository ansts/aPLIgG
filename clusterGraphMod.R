# Takes Am sequence graph and clusters it to about 2000 clusters 
# and contracts it to Am graph of the clusters which is returned. 
# The number of clusters is controlled by the resolution parameter
# which has to be established in advance (default=32).
# The edges are summed and the Group attribute is returned as Am 
# table of the values found in the clusters. The edge weight is 
# the log(modularity) of the subgraph formed of the pair of 
# clusters of the initial graph connected in the graph of 
# clusters. Outputs a graphml file with the cluster graph
# and returns the cluster graph as an igraph object. Size is
# the size of the clusters in the initial graph.
# Needs the matrix "gmatwtmx" calculated in iGrawate.

clusterGraphMod <- function(G, resol=1, fnm=NULL, nodes=14) {
  require(igraph)
  require(parallel)
  #require(pbapply)
  require(future.apply)
  require(numbers)
  source("repptrnalize.R")
  vrepptrnalize=Vectorize(repptrnalize)
  
  dGm=degree(G)
  dGmi=cut(dGm, 50, labels = F)
  dGmh=aggregate(dGm, by=list(dGmi), "length")
  
  pdf(file=paste("DegreeDistr",fnm,".pdf",sep=""),width=5, height=5)
  plot(dGmh, log="y", main=paste("Degree Distribution of Gctr",fnm,sep=""))
  dev.off()
  
  load("gmatwtmx")
  
  if (is.null(fnm)) fnm=round(runif(1,1,1e4),0)
  print("Leiden clustering.........")
  Gleiden=cluster_leiden(G, objective_function = "modularity", resolution_parameter = resol, n_iterations = 15)
  n=Gleiden$names
  #d=nrow(Gleiden$memberships)
  #if (d>1) d=d-1
  Gleiden=Gleiden$membership
  names(Gleiden)=n
  print(paste(c(max(Gleiden)," clusters and " ,sum(table(Gleiden)==1), " singlets"), sep="", collapse=""))
  save(Gleiden,file=paste("Gleiden",fnm,sep="")) 
  Grpsall=table(V(G)$Group)  
  gnm=names(Grpsall)
  x0=rep(0,length(gnm))
  names(x0)=gnm
  print("Contraction.........")
  Gctr=contract(G,mapping=Gleiden,vertex.attr.comb=list(name="concat", Group=function(x) {
    x3=table(x)
    x0[names(x3)]=x3
    return(x0)
  }))
  Gctr= simplify(Gctr, edge.attr.comb = list(weight="sum", LCS="ignore"))
  clnm=paste("C",seq_along(V(Gctr)), sep="")
  Gpeps=vertex_attr(Gctr)$name
  names(Gpeps)=clnm
  Gctr=delete_vertex_attr(Gctr, "name")
  Gctr=set_vertex_attr(Gctr, "name", value=clnm)
  Gctr=set_vertex_attr(Gctr, "size", value=lengths(Gpeps)*1.25+5)
  Gctr=set_vertex_attr(Gctr, "edges", value=sapply(Gpeps,function(x){
    g=induced_subgraph(G,x)
    sum(edge_attr(g)$weight)
  }))

  print("Calculating distances.........")

  vrp=lapply(Gpeps, function(l) vrepptrnalize(l))
  plan(multisession, workers=8)
  mtw=future_sapply(vrp, function(l1){                # maximal theoretical sums of weights 
      sapply(vrp, function(l2){
        mx=1/gmatwtmx[l1,l2]
        if (length(l1)==length(l2)){
          if (all(l1==l2)) res=sum(mx[lower.tri(mx)]+diag(mx)) else res=sum(mx)
        } else res=sum(mx)
        return(res)
      })
  })
  closeAllConnections()
  print("Maximal stength calculated.")
  
  n=length(V(Gctr))
  Am=matrix(0,n,n)
  rownames(Am)=names(V(Gctr))
  colnames(Am)=names(V(Gctr))
  eG=apply(ends(Gctr, E(Gctr), E(Gctr)),2,unlist)  # matrix of edge ends
  wG=edge.attributes(Gctr)$weight                  # sums of weights of edges between clusters
  Am[eG]=wG                                        # adjacency mx of clusters
  Am=t(Am)
  Am[eG]=wG
  wGin=array(rep(vertex_attr(Gctr)$edges, n), dim=c(n,n))
  wGin=wGin+t(wGin)                                # sums of weights inside all pairs of clusters 
  diag(wGin)=0
  mtwin=array(rep(diag(mtw), n), dim=c(n,n))       # maximal theoretical sums of weights inside clusters 
  mtwin=mtwin+t(mtwin)
  diag(mtwin)=1
  densin=wGin/mtwin
  densbtw=Am/mtw                                    # density - the diagonal is inside each cluster and the rest are for pairs of clusters
  wfinal=densin/densbtw                             # ratio of densities between and within
  rownames(wfinal)=names(V(Gctr))                       # wG is Am dissimilarity measure
  colnames(wfinal)=names(V(Gctr))
  wfinal=wfinal[eG]
  
  wfinal=log10(wfinal+min(wfinal[wfinal>0])*(wfinal==0)) 
  x=wfinal[order(wfinal)]
  x=x[length(x)]-x[length(x)-1]
  wfinal=wfinal-min(wfinal)+x/2
  
  Gctr=set.edge.attribute(Gctr, name="weight", value=wfinal)
  save(Gctr,file=paste("Gctr",fnm,sep=""))
  save(Gpeps,file="Gctrpeps")
  
  print("Prune large distances ...")
  thrng=range(edge.attributes(Gctr)$weight)
  if (ecount(Gctr)>200) seqthr=seq(thrng[1],thrng[2],length.out=200) else seqthr=sort(edge.attributes(Gctr)$weight)
  thrscan=t(sapply(seqthr, function(th){
    Gctrsm=delete.edges(Gctr,E(Gctr)[E(Gctr)$weight>th])
    cmp=components(Gctrsm)
    c(th,cmp$no,sum(cmp$csize<5), graph.density(Gctrsm))
  }))
 
  thr=min(thrscan[thrscan[,3]<(components(Gctr)$no*1.1),1])
  Gctrsm=delete.edges(Gctr,E(Gctr)[E(Gctr)$weight>thr])
  Gctrsm=set.edge.attribute(Gctrsm,"weight", value=1/edge_attr(Gctrsm)$weight) # Back to weight as Am similarity measure

  save(Gctrsm,file=paste("Gctrsm",fnm,sep=""))
 
  write_graph(Gctrsm,format = "graphml", file=paste("Gctrsm",fnm, ".graphml",sep=""))
 
  return(Gctrsm)
}
