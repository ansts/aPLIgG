grEmbed=function(G, att=NULL,attsc=NULL,a_e_ratio=0.67, k=20, fl="all"){
  require(igraph)
  require(MASS)
  require(uwot)
  
  if (is.null(att)) att=length(vertex_attr(G))
  # cmp=components(G)$csize
  # i=which.max(cmp)
  # G=induced.subgraph(G, vids=V(G)[cmp==i])
  ajm=distances(G, to=V(G))
  ajm[is.infinite(ajm)]=max(ajm[!is.infinite(ajm)]*10)
  emb=cmdscale(ajm, k=k)
  colnames(emb)=paste("Emb", seq_along(emb[1,]), sep="")
  x=as.matrix(apply(as.data.frame(vertex_attr(G)[att]), 2, as.numeric))
  x[is.na(x)]=0
  x[,att %in% attsc]=scale(x[,att %in% attsc])
  if (fl=="all") emb=cbind(emb,a_e_ratio)
  xy=umap(emb, n_neighbors = min(nrow(emb),50))
  xy=500*xy/diff(range(xy))
  G=set.vertex.attribute(G, "x", value=xy[,1])
  G=set.vertex.attribute(G, "y", value=xy[,2])
  G=set.vertex.attribute(G, "Label", value=vertex.attributes(G)$name)
  return(G)
}