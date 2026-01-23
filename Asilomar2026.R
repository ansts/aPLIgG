require(ape)
require(purrr)
require(MixGHD)
require(Biostrings)
require(pwalign)
require(parallel)
require(pbapply)
require(future.apply)
require(stringdist)
require(ggseqlogo)
require(reticulate)
require(chisq.posthoc.test)
require(FactoMineR)
require(factoextra)
require(scales)
require(pheatmap)
require(Matrix)
require(uwot)
require(matrixStats)
require(qualV)
require(rgl)
require(latexpdf)
require(tools)
require(MASS)
require(genieclust)
require(eulerr)
require(vioplot)
require(dunn.test)
library(ggplot2)
library(rstatix)
library(ggpubr)
require(corrplot)
require(stringi)
require(Rfast)
require(dbscan)
require(RandPro)
require(mixtools)
require(dendextend)
require(data.tree)
require(igraph)
require(reshape2)
require(LaplacesDemon)
require(gplm)
require(msa)
require(alluvial)
require(fastcluster)
require(data.tree)
require(dendextend)
require(dynamicTreeCut)
require(proxyC)
require(irlba)
require(bio3d)
require(alluvial)
require(vcd)
require(vcdExtra)

options(future.globals.maxSize=10*1024^3)
arpopt=list(maxiter=500000, tol=1e-8)

cpl=colorRampPalette(c("#000040FF","#0050AA9F","#10AA109F","#99FF5080","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
cplcb=colorRampPalette(c("#000000","#0072B2", "#56B4E9",  "#009E73","#F0E442", "#E69F00",  "#D55E00", "#CC79A7"), alpha=T)




# Direct lcs grouping ---------------------------------------------------------
gi=list()
flcm=rownames(fl[fl$cM,]); flcm=flcm[flcm %in% vgp]; gi[[1]]=induced_subgraph(Gpos, V(Gpos)[flcm])
flam=rownames(fl[fl$aM,]); flam=flam[flam %in% vgp]; gi[[2]]=induced_subgraph(Gpos, V(Gpos)[flam])
flcg=rownames(fl[fl$cG,]); flcg=flcg[flcg %in% vgp]; gi[[3]]=induced_subgraph(Gpos, V(Gpos)[flcg])
flag=rownames(fl[fl$aG,]); flag=flag[flag %in% vgp]; gi[[4]]=induced_subgraph(Gpos, V(Gpos)[flag])


grbylcs1=lapply(gi, function(g){
  ee=ends(g, E(g))
  ee=cbind(ee, edge_attr(g)$LCS)
  eeag=aggregate(ee[,1:2], by=list(ee[,3]), c)
  eeagpp=pbapply(eeag,1,function(l) {
    unlist(c(l[[2]],l[[3]]))
  })
  names(eeagpp)=eeag[,1]
  l=lengths(eeagpp)
  l5=l[nchar(names(l))==5]
  l6=l[nchar(names(l))==6]
  l65=unlist(lapply(names(l6),\(li) {
    x=sapply(1:6,\(i) {
      if (i>1 & i<6) paste0(substr(li,1,(i-1)),substr(li,(i+1),6),collapse="") else {
        if(i==1) return(substr(li,2,6)) else return(substr(li,1,5))
      }
    })
    lj=rep(l6[li],length(x)); names(lj)=x
    return(lj)
  }))
  l5=c(l5,l65);l5=aggregate(l5,by=list(names(l5)), sum)
  x=l5$x;names(x)=l5$Group.1
  return(x)
})


eeagpp1=lapply(gi, function(g){
  ee=ends(g, E(g))
  ee=cbind(ee, edge_attr(g)$LCS)
  eeag=aggregate(ee[,1:2], by=list(ee[,3]), c)
  eeagpp=pbapply(eeag,1,function(l) {
    unlist(c(l[[2]],l[[3]]))
  })
  names(eeagpp)=eeag[,1]
  return(eeagpp)
})
names(eeagpp1)=c("CM","AM","CG","AG")

names(grbylcs1)=c("CM","AM","CG","AG")

allcs=unique(unlist(rownames(grbylcs1)))
x=matrix(0,length(allcs), 4); rownames(x)=allcs; colnames(x)=names(grbylcs1);allcs=x
x=sapply(names(grbylcs1), function(n){
  y=allcs
  y[names(grbylcs1[[n]]),n]=grbylcs1[[n]]
  return(y[,n])
})
grbylcs1=x; rm(x)

tt=as.table(colSums(grbylcs1))
plan("multisession", workers=14)
grbylcs1_chsq=future_apply(grbylcs1,1, function(l){
  x=chisq.posthoc.test(as.table(cbind(l,tt-l)), method="BH", simulate.p.value=T)
  return(x[,3])
}, future.seed = T)
plan("sequential")

grbylcs1_chsq=t(grbylcs1_chsq)
grbylcs1_chsq=list(z=grbylcs1_chsq[,(1:(ncol(grbylcs1_chsq)/2))*2-1],
                   p=grbylcs1_chsq[,(1:(ncol(grbylcs1_chsq)/2))*2])
grbylcs1_chsq[[2]]=apply(grbylcs1_chsq[[2]],2,p.adjust, method="BH")
lvp=min(abs(grbylcs1_chsq[[1]])[grbylcs1_chsq[[2]]<0.05])
grbylcs1_sig=grbylcs1_chsq[[1]][pbapply(grbylcs1_chsq[[2]],1,function(l) any(l<0.05)),]
grbylcs1_sig[abs(grbylcs1_sig)<lvp]=0

colnames(grbylcs1_sig)=c("CM","AM","CG","AG")

X=grbylcs1_sig; X[X>0]=1; X[X<0]=-1; #X[X[,1]==-1,1]=0
lcsig1_prfls=aggregate(rownames(X), by=as.data.frame(X), c)
clcsig1=lcsig1_prfls$x; names(clcsig1)=paste("C",seq_along(clcsig1), sep="")
rownames(lcsig1_prfls)=names(clcsig1); clcsigsz1=lengths(clcsig1)
rm(X);gc()

plan("multisession", workers=4)
mnId1=future_lapply(eeagpp1, function(L) sapply(L, function(l) pslg$inverse(mean(pslg$transform(Idmapmepi[l])))))
plan("sequential")

clcsig1Idmn=t(sapply(clcsig1, function(l) sapply(l, function(il) {
  if (sum(l %in% names(il))>0) {
    x=il[l]
    pslg$inverse(mean(pslg$transform(x[!is.na(x)])))
  } else return(0)
})))

#### "################"--------------------------------------------------------
### new Id check --------------------------------------------------------------
clcsig1m=melt(clcsig1);x=clcsig1m$L1; names(x)=clcsig1m$value; clcsig1m=x; rm(x)
x=qgrams(names(clcsig1m)[nchar(names(clcsig1m))==5],nJ3,q=5)
x6=qgrams(names(clcsig1m)[nchar(names(clcsig1m))==6],nJ3,q=6)
x=x[,colProds(x)>0];x6=x6[,colProds(x6)>0]
clcsigId5_6=c(x[2,], x6[2,])
x=names(clcsig1m)[!(names(clcsig1m) %in% names(clcsigId5_6))]
y=rep(0,length(x)); names(y)=x
clcsigId5_6=c(clcsigId5_6,y)
clcsig1mId=clcsigId5_6[names(clcsig1m)]

clcsig1mIdbyCl=aggregate(clcsig1mId, by=list(clcsig1m), mean)
clcsIdunn=dunn.test(clcsig1mId,g=clcsig1m,method="bh")
x=clcsIdunn$comparisons;y=clcsIdunn$P.adjusted;z=clcsIdunn$Z
clcsIdunnsig=data.frame(Comparison=x[y<0.05],Z=z[y<0.05])
clcsIdunnsig=clcsIdunnsig[order(clcsIdunnsig$Z),]
xy=t(sapply(clcsIdunnsig[,1], \(dif) unlist(strsplit(dif, split=" - "))))
clcsIdunnsig=cbind(xy,clcsIdunnsig[,2])
lcsig1_prfl=cbind(lcsig1_prfls,lengths(lcsig1_prfls$x))
nclcs=nrow(lcsig1_prfl)
X=lcsig1_prfl[,1:4];X[X==-1]="Down";X[X==1]="Up";X[X==0]="NS"
Y=clcsig1mIdbyCl$x; names(Y)=clcsig1mIdbyCl$Group.1;Y=Y[rownames(X)]
alluvial(X, freq=log(lcsig1_prfl[,6]), col=cpl(nclcs))
alluvial(X, freq=log(lcsig1_prfl[,6]), col=cpl(nclcs)[rank(Y)])


ij=combn(6,5)
eeagpp1n5=lapply(eeagpp1,\(l1){
  x=names(l1)
  i6=nchar(x)==6
  n5=pblapply(x[i6],\(x6) {
    x6=unlist(strsplit(x6, split=""))
    apply(ij,2,\(j) paste(x6[j], collapse=""))
  })
  x[i6]=n5
  return(x)
})

plan(multisession, workers=14)
clcspp=future_sapply(lcsig1_prfl$x,\(l){
  as.character(unique(unlist(sapply(names(eeagpp1),\(nl) {
    j=unlist(sapply(l,\(li) grep(li,eeagpp1n5[[nl]])))
    eeagpp1[[nl]][j]
  }))))
})
plan(sequential)

lcsig1_prfl=cbind(lcsig1_prfl,NPep=lengths(clcspp))

qs=quantile(lcsig1_prfl[,6], c(0.33,0.67,1))
j=cut(lcsig1_prfl[,6],c(0,qs[1:2],28000), labels=F)
co=rank(Y)
alluvial(X[j==1,], freq=lcsig1_prfl[j==1,6], col=cplcb(nclcs)[co[j==1]],alpha=0.67)
alluvial(X[j==2,], freq=lcsig1_prfl[j==2,6], col=cplcb(nclcs)[co[j==2]],alpha=0.67)
alluvial(X[j==3,], freq=lcsig1_prfl[j==3,6], col=cplcb(nclcs)[co[j==3]],alpha=0.67)
z=1:100;barplot(rep(1,100),col=cplcb(100)[1:100], border=cplcb(100)[1:100], space=0)


proct=proc.time()
allumotfs=lapply(seq(clcspp),\(i){
  l=clcspp[[i]]
  print(c(i,length(l)))
  print(proc.time()-proct)
  if (length(l)>20) g=induced_subgraph(Gposu, V(Gposu)[l[l %in% names(V(Gposu))]]) else return(l)
  gcmp=components(g)
  lapply(seq(gcmp$no), \(i) {
    nm=names(gcmp$membership)[gcmp$membership==i]
    if (gcmp$csize[i]>1000) {
      gi=induced_subgraph(g, V(g)[nm])
      clstnhldn(gi)
    } else return(nm)
  })
})

x=lapply(allumotfs,\(l){
  lapply(l, \(l){
    if (is.numeric(l)){
      x=melt(l)
      aggregate(rownames(x), by=list(x[,1]), c)[,2]
    } else return(l)
  })  
})

x=lapply(seq(x),\(i){
  list_flatten(x[[i]])
})

allumotfs_flat=x; rm(x,xi)

pdf(file="allumtfs.pdf", width=25, height=25)
for (x in allumotfs_flat){
  xi=ggseqlogo(x)
  print(xi)
}
dev.off()

# diamap-ping -----------------------------------------------------------------

ij=combn(7,5)
mataa=t(pbsapply(matpep,\(p) unlist(strsplit(p, split=""))))
plan(multisession, workers=20)
matlcs5=future_apply(ij,2,\(j) list(apply(mataa[,j],1,\(l) paste(l,sep="",collapse=""))))
plan(sequential)
matlcs5=unlist(matlcs5, recursive = F)
matlcs5=do.call(cbind,matlcs5)
matlcs5=apply(matlcs5,1,unique)
matlcs5=unlist(matlcs5)
matlcs5t=table(matlcs5)
lmat5=length(matlcs5)
rm(mataa)

aa=AA_STANDARD; n=unique(nchar(L))
if (length(n)>1) stop("The sequences should be of equal length!")
x=c(sapply(aa,\(a1) sapply(aa,\(a2) paste(a1,a2,sep = ""))))
x=paste(rep(x,each=4),rep(0:(n-2),length(x)),sep="")
dpp=rep(0,length(x)); names(dpp)=x
AM0=array(0,dim=c(1600,1600), dimnames=list(names(dpp),names(dpp)))
nallcs=length(clcsig1m)

ncores=20
plan(multisession, workers=ncores)
matlcsAMs=pbsapply(1:100,\(i){
  tX=sample(matlcs5t,nallcs)
  AM=diamap(tX)
  AM0[rownames(AM),colnames(AM)]=AM
  gc()
  return(AM0)
}, simplify="array")
plan(sequential)

matlcsAM50p=pbapply(matlcsAMs,1:2,quantile, prob=0.5)
matlcsAM90p=pbapply(matlcsAMs,1:2,quantile, prob=0.9)
matlcsAM50p[matlcsAM50p==0]=min(matlcsAMs[matlcsAMs>0])/2
matlcsAM90p[matlcsAM90p==0]=min(matlcsAMs[matlcsAMs>0])/2

clcsig1mt=rowSums(grbylcs1)[names(clcsig1m)]
allcs5AM=diamap(clcsig1mt)

pslg2=pseudo_log_trans(base=2)
allcs5AMr=allcs5AM/matlcsAM90p
diag(allcs5AMr)=0
allcs5AMr[allcs5AMr<1.5]=0
allcs5AMr=allcs5AMr[rowSums(allcs5AMr)>0,colSums(allcs5AMr)>0]
gdi0=graph_from_adjacency_matrix(allcs5AMr,mode="undirected", weighted=T)
gdi0=simplify(gdi0)

plan(multisession, workers=20)
X=future_sapply(names(clcsig1m),displit, future.chunk.size=length(L) %/% 20)
plan(sequential)
X=melt(X)
X=X[X$value>0,]
allcs5pp=data.frame(Feature=as.character(X[,1]),Lcs5=as.character(X[,2]))

gdi=graph_from_adjacency_matrix(allcs5AMr,mode="undirected", weighted=T)
gdi=simplify(gdi)
gdism=induced_subgraph(gdi, V(gdi)[components(gdi)$membership==which.max(components(gdi)$csize)])
diag(matlcsAM50p)=0
gdimat=graph_from_adjacency_matrix(matlcsAM50p,mode="undirected", weighted=T)
gdimat=simplify(gdimat)
diag(allcs5AM)=0
gdineat=graph_from_adjacency_matrix(allcs5AM,mode="undirected", weighted=T)
gdineat=simplify(gdineat)
gdineat=delete_edges(gdineat, E(gdineat)[edge_attr(gdineat)$weight<quantile(edge_attr(gdineat)$weight, 0.9885)])

comsm=components(gdi)$membership %in% which(components(gdi)$csize<max(components(gdi)$csize))
aggregate(names(V(gdi))[comsm], by=list(components(gdi)$membership[comsm]),c)
clcompgdi=components(gdi)$membership[comsm]
clcompgdi=paste("co",clcompgdi, sep="_")
names(clcompgdi)=names(V(gdi))[comsm]

## Clique-based -----------------------------------------------------------
# mxqgdism=max_cliques(gdism, min=3)
# mxqgdism=lapply(mxqgdism,names)
clq3gdism=lapply(cliques(gdism,3,3), names)

plan(multisession, workers=10)
clq3_2pep=future_lapply(clq3gdism,\(qi){
  ppi=allcs5pp[allcs5pp[,1] %in% qi,2]
  tppi=table(ppi)
  good=names(tppi[tppi==3])  
  n=which(names(clcsig1m) %in% good)
  return(n)
}, future.chunk.size=length(clq3gdism) %/% 10, future.seed = T)
plan(sequential)

j=lengths(clq3_2pep)>1
clq3_2pep=clq3_2pep[j]
clq3gdsh=clq3gdism[j]

plan(multisession, workers=14)
clq3AdjM=future_sapply(clq3_2pep,\(l1){
  sapply(clq3_2pep,\(l2){
    length(intersect(l1,l2))/min(length(l1), length(l2))
  })
}, future.chunk.size=length(clq3_2pep) %/% 14)
plan(sequential)

diag(clq3AdjM)=0
dimnames(clq3AdjM)=list(seq_along(clq3_2pep),seq_along(clq3_2pep))
Mx=clq3AdjM;Mx[Mx<1]=0
gclq3=graph_from_adjacency_matrix(Mx, mode="undirected")
ajlgclq3=as_adj_list(gclq3)
jstay=pbsapply(names(ajlgclq3),\(ni){
  li=names(ajlgclq3[[ni]])
  if (length(li)>0) {
    ni=as.numeric(ni); li=as.numeric(li)
    xni=length(clq3_2pep[[ni]])
    xli=lengths(clq3_2pep[li])
    return(xni>max(xli))
  } else return(TRUE)
})
clq3_5OK=lapply(clq3_2pep[as.numeric(names(V(gclq3))[jstay])],\(ni) names(clcsig1m)[ni])


plan(multisession, workers=14)
clq3mtfs=future_lapply(clq3_5OK,\(ppi){
  rivtab=aggregate(names(clcsig1m[ppi]),by=list(clcsig1m[ppi]),list)
  clcspp=lengths(apply(rivtab,1,\(li) {
    as.character(unique(unlist(sapply(names(eeagpp1),\(nl) {
      j=unlist(sapply(li$x,\(lj) grep(lj,eeagpp1n5[[nl]])))
      eeagpp1[[nl]][j]
    }))))
  }, simplify = F))
  names(clcspp)=rivtab$Group.1
  tb=as.table(cbind(clcspp,lcsig1_prfl[rivtab$Group.1, "NPep"]-clcspp))
  if (nrow(tb)>1) {
    chsqt=chisq.posthoc.test(tb, simulate.p.value=T)
    i=(1:length(clcspp))*2; rivi=chsqt[i,][chsqt[i,3]<0.01&chsqt[i-1,3]>0,1]
  } else rivi=rownames(tb)
  aln=apply(as.matrix(msa(AAStringSet(ppi), method="Muscle", gapOpening=1.4)),1,paste,collapse="")
  return(list(rivi, aln))
}, future.chunk.size=length(clq3_5OK) %/% 14, future.seed = T)
plan(sequential)

i=sapply(clq3mtfs,\(l) (length(l[[1]])>0)&(length(l[[2]])>2))
clq3_touse=clq3gdism[i]
clq3mtfs_ok=clq3mtfs[i]

Ajmtfs=pbsapply(clq3mtfs_ok,\(m1){
  m1=m1[[2]]
  sapply(clq3mtfs_ok,\(m2){
    m2=m2[[2]]
    length(intersect(m1,m2))/min(lengths(list(m1,m2)))
  })
})
diag(Ajmtfs)=0
dimnames(Ajmtfs)=list(seq_along(clq3mtfs_ok),seq_along(clq3mtfs_ok))
gmtfs=graph_from_adjacency_matrix(as.matrix(Ajmtfs),  mode="undirected", weighted=T)
gmtfs05=delete_edges(gmtfs, E(gmtfs)[E(gmtfs)$weight<=0.5])

pdf(file="gmtfs05.pdf", width=15, height=15)
plot(gmtfs05, vertex.size=1, vertex.label=NA)
dev.off()

cmpmtfs05=components(gmtfs05)
jlc=which(cmpmtfs05$membership==which.max(cmpmtfs05$csize))
ggseqlogo(unique(unlist(sapply(jlc,\(i) clq3mtfs_ok[[i]][[2]]))))
table(unlist(sapply(jlc,\(i) clq3mtfs_ok[[i]][[1]])))
jlc=which(cmpmtfs05$membership==which(cmpmtfs05$csize==28))
ggseqlogo(unique(unlist(sapply(jlc,\(i) clq3mtfs_ok[[i]][[2]]))))
table(unlist(sapply(jlc,\(i) clq3mtfs_ok[[i]][[1]])))

gmtfs05lc=induced_subgraph(gmtfs05, V(gmtfs05)[jlc])
clgmtfs05lc=stbLeiden(gmtfs05lc)
j1s=cmpmtfs05$membership %in% which(cmpmtfs05$csize==1)
j1s=names(V(gmtfs05))[j1s]
jcmp=names(V(gmtfs05))[cmpmtfs05$membership %in%  which(cmpmtfs05$csize>1 & cmpmtfs05$csize<500)]
jcmp=aggregate(jcmp, 
               by=list(cmpmtfs05$membership[cmpmtfs05$membership %in%  which(cmpmtfs05$csize>1 & cmpmtfs05$csize<500)]), c)$x
jlcm=names(V(gmtfs05lc))
jlcm=aggregate(jlcm, 
               by=list(clgmtfs05lc), c)$x
jmtfnl=c(j1s,jcmp,jlcm)

plan(multisession, workers=14)
Motfnl=future_lapply(jmtfnl,\(l){
  L=clq3mtfs_ok[as.numeric(l)]
  if (length(l)>1) {
    # riv=unique(unlist(lapply(L,\(l1) l1[[1]])))
    mbr=unique(unlist(lapply(L,\(l1) l1[[2]])))
    mbr=gsub("-","",mbr)
    
    rivtab=aggregate(names(clcsig1m[mbr]),by=list(clcsig1m[mbr]),list)
    clcspp=lengths(apply(rivtab,1,\(li) {
      as.character(unique(unlist(sapply(names(eeagpp1),\(nl) {
        j=unlist(sapply(li$x,\(lj) grep(lj,eeagpp1n5[[nl]])))
        eeagpp1[[nl]][j]
      }))))
    }, simplify = F))
    names(clcspp)=rivtab$Group.1
    tb=as.table(cbind(clcspp,lcsig1_prfl[rivtab$Group.1, "NPep"]-clcspp))
    if (nrow(tb)>1) {
      chsqt=chisq.posthoc.test(tb, simulate.p.value=T)
      i=(1:length(clcspp))*2; rivi=chsqt[i,][chsqt[i,3]<0.01&chsqt[i-1,3]>0,1]
    } else rivi=rownames(tb)
    
    aln=as.character(apply(as.matrix(msa(AAStringSet(mbr), method="Muscle", gapOpening=1.4)),1,paste,collapse=""))
    L=list(rivi,aln)
  }
  return(unlist(L, recursive=F))
})
plan(sequential)

j=sapply(Motfnl,\(l){
  length(l[[1]])>0
})

Motfnl=Motfnl[j]

xy=t(sapply(Motfnl, \(mi) {
  c(Ncats=length(mi[[1]]), Nseqs=length(mi[[2]]))
}))

plot(xy, log="xy", pch=16, col=rgb(0,0,0,0.2))

barplot(table(xy[,"Ncats"]), ylab="N Categories")

catMotfnl=lapply(Motfnl,\(l){
  l[[1]]
})

x=catMotfnl[lengths(catMotfnl)>1]
z=matrix(0,55,55, dimnames=list(rownames(lcsig1_prfl),rownames(lcsig1_prfl)))
catcrossr=pbsapply(x,\(l){
  ij=combn(length(l),2)
  apply(ij,2,\(ii) {
    z[l[ii[1]],l[ii[2]]]<<-z[l[ii[1]],l[ii[2]]]+1
    z[l[ii[2]],l[ii[1]]]<<-z[l[ii[2]],l[ii[1]]]+1
  })
})


Y=clcsig1mIdbyCl$x; names(Y)=clcsig1mIdbyCl$Group.1;Y=Y[rownames(X)]
co=rank(Y);names(co)=names(Y)
z=max(z)-z
xyz=cmdscale(z);xyz=xyz+runif(nrow(xyz)*2,-diff(range(xyz))/20,diff(range(xyz))/20)
plot(xyz, cex=0, xlab="D1", ylab="D2", main="Expression Category Crossreactivty Map")
pointLabel(xyz, labels=rownames(z), cex=rank(lcsig1_prfl[,6])/40+0.5, 
           col=cplcb(nclcs)[co])

CMorAM=rownames(lcsig1_prfl[(lcsig1_prfl[,1])==1 & (lcsig1_prfl[,2])<1 ,])
crosclmtfs=list(c("C52","C42","C50"),
                c("C41","C55","C49","C54"),
                c("C12","C36","C35"),
                c("C3","C4","C5"),
                c("C1","C24","C21","C37","C29"))
rN=rank(lcsig1_prfl[,6]); names(rN)=rownames(lcsig1_prfl)
for (i in seq_along(crosclmtfs)){
  plot(xyz[crosclmtfs[[i]],], cex=0, xlab="D1", ylab="D2", 
       xlim=range(xyz),ylim=range(xyz))
  pointLabel(xyz[crosclmtfs[[i]],], labels=crosclmtfs[[i]], 
             cex=rN[crosclmtfs[[i]]]/40+0.5, 
             col=cplcb(nclcs)[co[crosclmtfs[[i]]]])
}

X=lcsig1_prfl[,1:4];X[X==-1]="Down";X[X==1]="Up";X[X==0]="NS"

for (i in seq_along(crosclmtfs)){
  col=rep("#FFFFFF", 55);names(col)=rownames(X)
  col[crosclmtfs[[i]]]=cplcb(nclcs)[co[crosclmtfs[[i]]]]
  trns=rep(0,55); names(trns)=rownames(X)
  trns[crosclmtfs[[i]]]=0.7
  alluvial(X[,1:4], 
           freq=(lcsig1_prfl[,7])/10,
           col=col, alpha = trns, border=col,
           layer = !(rownames(X) %in% crosclmtfs[[i]]))
}

commtf=lapply(crosclmtfs,\(li) {
  j=lapply(li,\(lj) grep(lj,Motfnl))
  tj=table(unlist(j))
  ntj=names(tj)[tj==length(li)]
  prop=length(ntj)/lengths(j)
  ls=sapply(Motfnl[as.numeric(ntj)],\(mi) length(mi[[2]]))
  ntj[order(ls, decreasing = T)]
  return(list(ntj[1:min(length(ntj),5)], prop))
})

pdf(file="comtflogos.pdf", width=7, height=7)
for (l in commtf){
  l=l[[1]]
  l=lapply(Motfnl[as.numeric(l)],\(l1) l1[[2]])
  print(ggseqlogo(l))
}
dev.off()

unimtft=table(unlist(catMotfnl[lengths(catMotfnl)==1]))
barplot(sort(as.array(unimtft),decreasing = T), las=2, ylab="N specific motifs")

rn=rownames(lcsig1_prfl);unmn=names(unimtft);rnunmn=rn[rn %in% unmn]
XY=cbind((lcsig1_prfl[rnunmn,7]),(c(unimtft)[rnunmn]))
plot(XY,cex=0, log="xy", xlab="N peptides in the category", ylab="N specific motifs")
pointLabel(XY, labels=rnunmn)

unimtf=Motfnl[lengths(sapply(Motfnl,\(x) x[[1]]))==1]
lunm=lapply(unimtf,\(l) l[[2]])
unimtf=aggregate(seq_along(lunm),by=list(sapply(unimtf,\(l) l[[1]])), list)
x=lapply(unimtf$x,\(l) lunm[l])
names(x)=unimtf$Group.1
unimtf=x

pdf(file="unimotifs,pdf", width=15, height=15)
for (nm in names(unimtf)){
  print(ggseqlogo(unimtf[[nm]])+ggtitle(nm))
}
dev.off()

# Gapped dimer feature maps ------------------------------------------------

pslg3=pseudo_log_trans(sigma=1e-9, base=10)
w=edge_attr(gdineat)$weight;w=pslg3$transform(w)
gdineat=set_edge_attr(gdineat, name="weight", value=w)
w=edge_attr(gdi0)$weight;w=pslg3$transform(w)
gdi0=set_edge_attr(gdi0, name="weight", value=w)
w=edge_attr(gdimat)$weight;w=pslg3$transform(w)
gdimat=set_edge_attr(gdimat, name="weight", value=w)

GDI=union(gdineat,gdimat,gdi0)
w=sapply(1:3, \(i){
  w=paste("weight",i,sep="_")
  w0=edge_attr(GDI)[[w]]
  GDI<<-delete_edge_attr(GDI, name=w)
  w0[is.na(w0)]=0
  return(w0)
})

w=rowMaxs(w, value=T)
GDI=set_edge_attr(GDI, name="weight", value=w)
GDI=delete_edges(GDI, E(GDI)[edge_attr(GDI)$weight<7.5])
GDI=simplify(GDI)
GDI=induced_subgraph(GDI, V(GDI)[components(GDI)$membership==1])
LGDI=embed_laplacian_matrix(GDI, no=800, which="sa", type="I-DAD")
uLGDI=umap(LGDI$X[,2:46], verbose = T)
cuLGDI=uLGDI/(sqrt(rowSums(uLGDI^2)))
pdf(file="cGDI.pdf", width=20, height=20)
plot(GDI, layout=cuLGDI, vertex.size=1, vertex.label=NA, edge.width=0.05)
dev.off()

clgdi0=stbLeiden(gdi0, thrmn = 0.9, thrmx = 0.95)
clgdineat=stbLeiden(gdineat)
clgdimat=stbLeiden(gdimat, thrmn = 0.9, thrmx = 0.95)

clmx_gdineatmat=as.matrix(table(clgdineat,clgdimat))
r1=rowSums(clmx_gdineatmat)>2;c1=colSums(clmx_gdineatmat)>2
pheatmap(clmx_gdineatmat[r1,c1], col=cplcb(11), clustering_method = "ward.D2")

cplcb1=colorRampPalette(c("#000000","#0000FF","#0055FF","#0072B2", "#56B4E9",  "#009E73","#F0E442", "#E69F00",  "#D55E00", "#CC79A7", bias=2), alpha=T)

clmx_gdi0mat=as.matrix(table(clgdi0,clgdimat))
r1=rowSums(clmx_gdi0mat)>2;c1=colSums(clmx_gdi0mat)>2
pheatmap(clmx_gdi0mat[r1,c1], col=cplcb1(59),  clustering_method = "ward.D2")

gdi0large=graph_from_adjacency_matrix(allcs5AM/matlcsAM90p,mode="undirected", weighted=T)
gdi0large=simplify(gdi0large)
x=components(gdi0large)$membership
gdi0large=induced_subgraph(gdi0large, V(gdi0large)[x==1])
gdi0large=delete_edges(gdi0large,E(gdi0large)[edge_attr(gdi0large)$weight<3.5e-3])

Lgdism=embed_laplacian_matrix(gdism, no=200, which="sa", type="I-DAD")
D=dim_select(Lgdism$D)
uLgdism=umap(Lgdism$X[,2:D], verbose = T, min_dist = 0.01)
uLgdism5=umap(Lgdism$X[,2:D], verbose = T, min_dist = 0.01, n_components = 5, n_neighbors=25)
rownames(uLgdism5)=names(V(gdism))
dx=sqrt(rowSums(uLgdism^2))*cbind(runif(nrow(uLgdism),-0.025,0.025),runif(nrow(uLgdism),-0.025,0.025))

cLeigismre=as.numeric(cLeigism[-grep("co_", cLeigism)])
names(cLeigismre)=names(cLeigism[-grep("co_", cLeigism)])
tcLeigismre=table(cLeigismre)
clcol=cLeigismre*(cLeigismre %in% as.numeric(names(tcLeigismre))[tcLeigismre>3])+1
Mx=(uLgdism+dx)
rownames(Mx)=names(V(gdism))
Mx=Mx[names(clcol),]
pdf(file="gdi0L.pdf",gdi0largepdf(file="gdi0L.pdf", width=20, height=20)
    plot(Mx, cex=0.75, pch=16, col=clcol)
dev.off()
    
cluLgdism=hclust(dist(uLgdism5),method="ward.D2")
cLuLmmb=cutreeDynamic(cluLgdism)
names(cLuLmmb)=names(V(gdism))

cLei2m=melt(cLei2); x=cLei2m$L1; names(x)=cLei2m$value;cLei2m=x
crossclgdism=table(cLei2m, cLuLmmb[names(cLei2m)])
colx=cplcb1(max(crossclgdism))
pdf(file="crossclusters_gdism.pdf", width=6, height=16)
pheatmap(crossclgdism, col=c("#FFFFFF",colx), clustering_method = "ward.D2")
dev.off()

# Final rivers -----------------------------------------------------------
cplcb=colorRampPalette(c("#0000FF","#0072B2", "#56B4E9",  "#009E73","#F0E442", "#E69F00",  "#D55E00", "#CC79A7"), alpha=T,bias=2)
X=lcsig1_prfl[,1:4];X[X==-1]="Down";X[X==1]="Up";X[X==0]="NS"
Y=clcsig1mIdbyCl$x; names(Y)=clcsig1mIdbyCl$Group.1;Y=Y[rownames(X)]
co=Y-min(Y)+1
j=(X[,2]=="NS"&X[,4]=="NS")|rownames(X) %in% c("C12","C35") 
Xj=X[!j,]; Ns=sqrt(lcsig1_prfl$NPep)[!j]; coj=co[!j]
ij=order(Ns, decreasing = T)

Xjij=Xj[ij,]; Nsij=Ns[ij]
ij1=(Xjij[,1] =="Up" & Xjij[,2] %in% c("NS","Down"))|(Xjij[,1] =="NS" & Xjij[,2] =="Down")
col=cplcb(max(Y)-min(Y)+1)[co];col=col[!j][ij]

alluvial(Xjij, 
         freq=Nsij+25,
         col=col,alpha=0.67,hide=!ij1)

ij1=(Xjij[,1] =="Down" & Xjij[,2] %in% c("NS","Up"))|(Xjij[,1] =="NS" & Xjij[,2] =="Up")
col=cplcb(max(Y)-min(Y)+1)[coj[ij]]

alluvial(Xjij, 
         freq=Nsij,
         col=col,alpha=0.67, hide=!ij1)

ij1=(Xjij[,4] =="Up" & Xjij[,3] %in% c("NS","Down"))|(Xjij[,4] =="NS" & Xjij[,3] =="Down")
col=cplcb(max(Y)-min(Y)+1)[coj[ij]]

alluvial(Xjij, 
         freq=Nsij+15,
         col=col,alpha=0.67, hide=!ij1)


idj=quantile(coj, c(0.33, 0.67))
ij1=coj[ij]>idj[1]&coj[ij]<idj[2]
col=cplcb(max(Y)-min(Y)+1)[coj[ij]]

alluvial(Xjij, 
         freq=Nsij+15,
         col=col,alpha=0.67, hide=!ij1)


colorbar=barplot(array(rep((max(Y)-min(Y)+1)/length(co),length(co)),dim=c(length(co),1)),
                 col=cplcb(length(co)), 
                 border=NA, ylab="Frequency of Idiotopes")

# Check motifs against known epitopes -----------------------------------------

knownEpi=read.csv("APLepis.csv")

mtfXcr=pbsapply(Motfnl,\(mi) checkAPLmtfX(mi[[2]]))

plan(multisession, workers=14)
bckgmtfX=future_sapply(names(matlcs5t),\(mi) {
  checkAPLmtfX(mi)
}, future.chunk.size=length(matlcs5t) %/% 14)
plan(sequential)     

n=length(mtfXcr)

maxbcgmtfX=sapply(1:1000, \(i) max(sample(bckgmtfX,n)))
h=hist(maxbcgmtfX, breaks=30)
j=which(mtfXcr>max(maxbcgmtfX))

mtfX=Motfnl[j]
mtfXchk=sapply(mtfX, \(mi){
  j=sapply(mi[[2]],\(mj) checkAPLmtfX_sh(mj))
  if (any(j)) return(c(mi[1],mi[[2]][j]))
})

# 40 5mers in 15 motifs in 12 categories - 7/15 motifs in C31  

Xlogo=ggseqlogo(unlist(sapply(mtfXchk, \(x) x[-1])))
print(Xlogo)

allpep=unique(unlist(eeagpp1))
epi=knownEpi$Sequence[-5]
x=qgrams(allpep,epi, q=5)
x=x[,colProds(x)>0]
epi5mer=colnames(x)

x=qgrams(unique(unlist(lcsig1_prfl$x)),epi, q=5)
x=x[,colProds(x)>0]
epi5lcs=colnames(x)


mtfs_per_cat=table(unlist(catMotfnl))
tmtfxC=table(unlist(sapply(mtfXchk,\(x) x[[1]])))
X=lcsig1_prfl[,1:4];X[X==-1]="Down";X[X==1]="Up";X[X==0]="NS"
Y=clcsig1mIdbyCl$x; names(Y)=clcsig1mIdbyCl$Group.1;Y=Y[rownames(X)]
co=Y-min(Y)+1
j=rownames(X) %in% names(tmtfxC)
Xj=X[j,]; Ns=sqrt(lcsig1_prfl$NPep)

col=cplcb(max(Y)-min(Y)+1)[co]

alluvial(X, 
         freq=Ns+25,
         col=col,alpha=0.67, hide=(!j))

text(x=1.35,y=0.875,labels="C31 (7/199)", col="#851000")
text(x=1.33,y=0.835,labels="C8 (2/221)", col="#851000")
text(x=1.33,y=0.74,labels="C6 (1/149)", col="#851000")
text(x=1.33,y=0.69,labels="C33 (1/54)", col="#851000")
text(x=1.33,y=0.5,labels="C16 (1/211)", col="#851000")
text(x=1.33,y=0.44,labels="C7 (2/123)", col="#851000")
text(x=1.33,y=0.415,labels="C23 (1/86)", col="#851000")
text(x=1.33,y=0.35,labels="C36 (1/113)", col="#851000")
text(x=1.33,y=0.315,labels="C13 (1/102)", col="#851000")
text(x=1.33,y=0.22,labels="C42 (1/318)", col="#851000")
text(x=3.6,y=0.77,labels="C54 (1/360)", col="#851000")
text(x=3,y=0.7,labels="C12 (1/104)", col="#851000")

# River Stat ------------------------------------------------------------------

clcsIdunnsig1=t(apply(clcsIdunnsig,1,\(l){
  x=as.numeric(l[3])
  if (x<0) {
    y=l[2]
    l[2]=l[1]
    l[1]=y
    x=-x
  }
  return(c(l[1],l[2], Z=x))
}))
rownames(clcsIdunnsig1)=NULL
clcsIdunnsig1=as.data.frame(clcsIdunnsig1)
clcsIdunnsig1=clcsIdunnsig1[order(clcsIdunnsig1[,1],clcsIdunnsig1[,2]),]

X=grbylcs1_sig; X[X>0]=1; X[X<0]=-1
X=cbind(X,clcsig1mId[rownames(X)])
X[,5]=cut(X[,5], c(-1,0.5,10.5,1e6), labels = F)
X=apply(X, 2, factor)
colnames(X)[5]="Id"
Xa=aggregate(rownames(X), by=as.data.frame(X), length)
Xt=xtabs(x~CM+AM+CG+AG+Id, data=Xa)
model_indep=loglm(~ CM+AM+CG+AG+Id, data = Xt)

stresid=residuals(model_indep, type = "pearson")
table_allways=xtabs(x~CM+AM+CG+AG+Id, data = Xa)
observed=as.vector(table_allways); n_total=sum(observed)
expected=as.vector(fitted(model_indep))
margins=sapply(1:5,\(i) apply(table_allways, i, sum) / n_total)
dims=dim(table_allways)
indices=expand.grid(1:dims[1], 1:dims[2], 1:dims[3], 1:dims[4], 1:dims[5])
var_adjustment=rep(1, length(observed))
for(i in 1:nrow(indices)) {
  var_adjustment[i]=(1 - margins[indices[i, 1],1]) *
    (1 - margins[indices[i, 2],2]) *
    (1 - margins[indices[i, 3],3]) *
    (1 - margins[indices[i, 4],4]) *
    (1 - margins[indices[i, 5],5])
}

stresid=stresid/sqrt(expected * var_adjustment)
presid=array(p.adjust(1-pnorm(abs(stresid))), dim=dim(stresid), dimnames=dimnames(Xt))
j=which(presid<0.05, arr.ind = T)
sigresid=cbind(j,stresid[j])
sigresid[order(sigresid[,6], decreasing=T),]

# 
wilcox.test(clcsig1mId[rownames(X)][X[,"CM"]=="1"],clcsig1mId[rownames(X)][X[,"CM"]=="-1"])
wilcox.test(clcsig1mId[rownames(X)][X[,"CM"]=="1"],clcsig1mId[rownames(X)][X[,"CM"]=="0"])
boxplot(rank(clcsig1mId[rownames(X)])~X[,"CM"], notch=T, xlab="CM", ylab="Rank of N Id Hits")
wilcox.test(clcsig1mId[rownames(X)][X[,"AM"]=="1"],clcsig1mId[rownames(X)][X[,"AM"]=="-1"])
wilcox.test(clcsig1mId[rownames(X)][X[,"AM"]=="-1"],clcsig1mId[rownames(X)][X[,"AM"]=="0"])
boxplot(rank(clcsig1mId[rownames(X)])~X[,"AM"], notch=T, xlab="AM", ylab="Rank of N Id Hits")
wilcox.test(clcsig1mId[rownames(X)][X[,"AG"]=="1"],clcsig1mId[rownames(X)][X[,"AG"]=="-1"])
boxplot(rank(clcsig1mId[rownames(X)])~X[,"AG"], notch=T, xlab="AG", ylab="Rank of N Id Hits")
wilcox.test(clcsig1mId[rownames(X)][X[,"CG"]=="1"],clcsig1mId[rownames(X)][X[,"CG"]=="-1"])
wilcox.test(clcsig1mId[rownames(X)][X[,"CG"]=="1"],clcsig1mId[rownames(X)][X[,"CG"]=="0"])
boxplot(rank(clcsig1mId[rownames(X)])~X[,"CG"], notch=T, xlab="CG", ylab="Rank of N Id Hits")

