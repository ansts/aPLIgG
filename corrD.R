corrD=function(G){
  require(igraph)
  require(parallel)
  require(pbapply)
  require(Rfast)
  
  cmp=components(G)$membership
  i=which.max(table(cmp))
  G=induced.subgraph(G, vids=V(G)[cmp==i])
  N=length(V(G))
  mx=distances(G, to=V(G), weights = NA)
  emb=cmdscale(mx,k=N-3, eig=T, list.=T)
  emb=emb$points[,1:sum(emb$eig>0)]
  
  D=sapply(seq_along(emb[,1]), function(i){
    #print(i)
    d=dista(t(emb[i,]), emb, trans=F)[-i]
    s=cut(d,100, labels=F)
    ts=table(s)
    ds=sapply(as.numeric(names(ts)), function(si) max(d[s==si]))
    crv=cbind(log(ds),log(cumsum(ts)/(sum(ts))))
    xlo=diff(sapply(seq_along(crv[,1])[-(1:2)], function(i){
      res=summary(lm(crv[1:i,2]~crv[1:i,1]))
      res$r.squared
    }))
    xlo=min(which(xlo>0))
    xup=sapply(seq_along(crv[,1])[(xlo+4):length(crv[,1])], function(i){
      res=summary(lm(crv[xlo:i,2]~crv[xlo:i,1]))
      res$r.squared
    })
      xup=xlo+3+which.max(xup)
    #plot(crv, main=i)
    #abline(lm(crv[xlo:xup,2]~crv[xlo:xup,1]))
    summary(lm(crv[xlo:xup,2]~crv[xlo:xup,1]))$coefficients[2,1]
  })
  return(median(D))
}

