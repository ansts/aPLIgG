grEmbed=function(G, att=NULL,attsc=NULL,a_e_ratio=0.67, fl="all"){
  require(igraph)
  require(MASS)
  require(uwot)
  
  if (is.null(att)) att=length(vertex_attr(G))
  if (is.null(attsc)) attsc=att
  G=set_edge_attr(G, name="weight", value=1/edge_attr(G)$weight)
  # cmp=components(G)$csize
  # i=which.max(cmp)
  # G=induced.subgraph(G, vids=V(G)[cmp==i])
  # ajm=distances(G, to=V(G))
  ajm=as.matrix(as_adjacency_matrix(G, attr = "weight"))
  #ajm[is.infinite(ajm)]=max(ajm[!is.infinite(ajm)]*3)
  ajm[ajm==0]=max(ajm)*1000
  diag(ajm)=0
  x=summary(prcomp(ajm))$importance[2,]
  x=sapply(seq_along(x), function(i) sum(x[1:i]))
  k=min(which(x>0.95))
  emb=cmdscale(ajm, k=k)
  colnames(emb)=paste("Emb", seq_along(emb[1,]), sep="")
  
  x=as.matrix(apply(as.data.frame(vertex_attr(G)[att]), 2, as.numeric))
  x[is.na(x)]=0
  x[,att %in% attsc]=scale(x[,att %in% attsc])
  if (fl=="all") emb=cbind(emb,a_e_ratio*x)
  if (fl=="par") emb=x
  
  xy=umap(emb, n_neighbors = round(sqrt(vcount(G))))    
  xy=500*xy/diff(range(xy))
  G=set.vertex.attribute(G, "x", value=xy[,1])
  G=set.vertex.attribute(G, "y", value=xy[,2])
  G=set.vertex.attribute(G, "Label", value=vertex.attributes(G)$name)
  G=set.vertex.attribute(G, "size", value=scale(vertex.attributes(G)$size, center=F))
  return(G)
}