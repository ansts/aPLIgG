grEmbed=function(G, att=NULL,a_e_ratio=0.67, k=20){
  require(igraph)
  require(MASS)
  require(uwot)
  
  if (is.null(att)) att=length(vertex_attr(G))
  # cmp=components(G)$csize
  # i=which.max(cmp)
  # G=induced.subgraph(G, vids=V(G)[cmp==i])
  ajm=distances(G, to=V(G))
  ajm[is.infinite(ajm)]=max(ajm[!is.infinite(ajm)]*10)
  emb=isoMDS(ajm, k=k, tol=1e-4)$points
  colnames(emb)=paste("Emb", seq_along(emb[1,]), sep="")
  x=as.matrix(apply(as.data.frame(vertex_attr(G)[att]), 2, as.numeric))
  emb=cbind(emb,a_e_ratio*scale(x))
  xy=umap(emb, n_neighbors = min(nrow(emb),50))
  xy=500*xy/diff(range(xy))
  G=set.vertex.attribute(G, "x", value=xy[,1])
  G=set.vertex.attribute(G, "y", value=xy[,2])
  G=set.vertex.attribute(G, "Label", value=vertex.attributes(G)$name)
  return(G)
}