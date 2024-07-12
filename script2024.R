# packages ---------------------------------------------------------------------

require(parallel)
require(igraph)
require(reshape2)
require(stringdist)
require(matrixStats)
require(chisq.posthoc.test)
require(heatmap3)
require(Biostrings)
require(uwot)
require(corrplot)
require(stringi)
require(qualV)
require(eulerr)
require(ggplot2)
require(ggseqlogo)
require(ggalluvial)
require(MASS)
require(multcomp)
require(smallstuff)
require(pbapply)
require(infotheo)
require(entropy)
require(rsetse)
require(R.utils)
require(vioplot)
require(matlib)
require(splines)
require(future.apply)
require(dunn.test)
require(rgl)
require(ggseqlogo)
require(PST)


load("mixnew")
load("flnew")

# Isotype Graphs --------------------------------------------------------------

## IgM ------------------------------------------------------------------------

X=mix[fl$cM|fl$aM]
x=adjL(X)
AL=lapply(list.files(pattern="ALm_ "),function(f) read.delim(file=f, stringsAsFactors = F, header=F))
AL=unlist(AL, recursive = F)
AL=unlist(AL, recursive = F)
AL=lapply(AL, strsplit, split=" ")
AL=unlist(AL, recursive = F)
nms=sapply(AL, function(l) l[1])
names(AL)=nms
AL=lapply(AL,function(l) l[-1])
save(AL,file="ALM")

load("ALM")
EL=melt(AL)
Gm=graph_from_edgelist(as.matrix(EL), directed = F)
Gm=simplify(Gm)
dGm=degree(Gm)

options(future.globals.maxSize = 1.0 * 1e9)
plan("multisession", workers=16)
x=t(future_sapply(1:35, function(i){    
  z=cluster_leiden(Gm, objective_function = "modularity", resolution_parameter = i)$quality
  return(c(i,z))
}))
closeAllCosigections()

clusGm=cluster_leiden(Gm, objective_function = "modularity", resolution_parameter = 50)

## IgG ------------------------------------------------------------------------

X=mix[fl$cG|fl$aG]
x=adjL(X)
AL=lapply(list.files(pattern="ALm_ "),function(fl) read.delim(file=fl, stringsAsFactors = F, header=F))
AL=unlist(AL, recursive = F)
AL=unlist(AL, recursive = F)
AL=lapply(AL, strsplit, split=" ")
AL=unlist(AL, recursive = F)
AL=unlist(AL, recursive = F)
nms=sapply(AL, function(l) l[1])
names(AL)=nms
AL=lapply(AL,function(l) l[-1])
save(AL,file="ALG")

load("ALG")
EL=melt(AL)
Gg=graph_from_edgelist(as.matrix(EL), directed = F)
Gg=simplify(Gg)
dGg=degree(Gg)
save(Gg, file="Gg")
rm(AL,EL)
gc()

plan("multisession", workers=16)
x=t(future_sapply(1:35, function(i){    
  z=cluster_leiden(Gg, objective_function = "modularity", resolution_parameter = i)$quality
  return(c(i,z))
}))
closeAllCosigections()
clusGG=cluster_leiden(Gg, objective_function = "modularity", resolution_parameter = 50)

# Distribution of disease specific mimotopes ---------------------------------
vm=names(V(Gm)); vg=names(V(Gg))
tm=table(Cluster=clusGm$membership,Controls=fl[vm,2],APLS=fl[vm,3])
tm=cbind(A=tm[,1,2],AC=tm[,2,2],C=tm[,2,1])
chtm=chisq.posthoc.test(tm, simulate.p.value=T)
n=nrow(chtm)/2; jp=2*(1:n); jz=jp-1
chtm=chtm[jz,][rowAnys(chtm[jp,]<0.05),c(1,3:5)]; colnames(chtm)[1]="Cluster"
chtm=chtm[order(chtm[,2]),]

tg=table(Cluster=clusGG$membership,Controls=fl[vg,4],APLS=fl[vg,5])
tg=cbind(A=tg[,1,2],AC=tg[,2,2],C=tg[,2,1])
chtg=chisq.posthoc.test(tg, simulate.p.value=T)
n=nrow(chtg)/2; jp=2*(1:n); jz=jp-1
chtg=chtg[jz,][rowAnys(chtg[jp,]<0.05),c(1,3:5)]; colnames(chtg)[1]="Cluster"
chtg=chtg[order(chtg[,2]),]

sigcl=lapply(list(IgM=chtm, IgG=chtg), function(l){
  list(C=as.numeric(l[l[,2]<(-3) | l[,4]>3,1]),A=as.numeric(l[l[,4]<(-3) | l[,2]>3,1]))
})

msigcl=melt(sigcl)
pepsig=apply(msigcl, 1, function(l){
  if (l[3] == "IgM") v=vm[clusGm$membership==as.numeric(l[1])] 
          else v=vg[clusGG$membership==as.numeric(l[1])]
  return(v)
})

X=aggregate(lengths(pepsig), by=list(msigcl[,2],msigcl[,3]), "sum")
X$x/rev(colSums(fl[,-1]))

aggregate(seq_along(msigcl[,1]), by=list(msigcl[,2],msigcl[,3]), "length")


# Cluster comparison ---------------------------------------------------------

allpepsig=unique(unlist(pepsig))

clmx=pbapply(msigcl,1,function(l1){
      if (l1[3]=="IgM") pp1=vm[clusGm$membership==as.numeric(l1[1])] else pp1=vg[clusGG$membership==as.numeric(l1[1])] 
      apply(msigcl,1,function(l2){
        if (l2[3]=="IgM") pp2=vm[clusGm$membership==as.numeric(l2[1])] else pp2=vg[clusGG$membership==as.numeric(l2[1])] 
        length(intersect(pp1,pp2))
      })
})
diag(clmx)=0

vsigtbl=sapply(1:4, function(i) {
  l=unlist(sigcl, recursive = F)[[i]]
  pp=unlist(sapply(l, function(cl){
    if (i %in% 1:2) x=vm[clusGm$membership %in% l] else x=vg[clusGG$membership %in% l]
    return(x)
  }))
  allpepsig %in% pp
})
colnames(vsigtbl)=names(unlist(sigcl, recursive = F))
rownames(vsigtbl)=allpepsig

vsig=euler(vsigtbl, shape="ellipse")
plot(vsig, quantities=T)

intrscts=rownames(vsigtbl)[rowSums(vsigtbl)>1]



# Merged graph ---------------------------------------------------------------

G=graph.union(Gm,Gg)
G=simplify(G)
VG=names(V(G))
mnD=mean(degree(G))
dG=degree(G)
save(G, file="G")
AjL=as_adj_list(G)
save(AjL,file="AjL")
rm(G)


tb=melt(table(fl[,-1]))
nms=apply(tb[,1:4],1, function(l) paste(colnames(tb[,-5])[l], collapse = ""))
tb=tb[,5]; names(tb)=nms
tb[1]=0;#tb=tb/sum(tb)
N=as.numeric(nrow(fl[rowSums(fl[,-1])>0,-1])) #=sum(tb)
clszM=sapply(sigcl[[1]], function(L) sapply(unlist(L), function(cl) sum(clusGm$membership==cl)))
clszG=sapply(sigcl[[2]], function(L) sapply(unlist(L), function(cl) sum(clusGG$membership==cl)))
# 
# ORDM=lapply(clszM, function(L) array(L%x%tb, dim=c(length(tb),length(L))))
# ORDG=lapply(clszG, function(L) array(L%x%tb, dim=c(length(tb),length(L))))

ddM=pblapply(1:2, function(i){
  L=sigcl[[1]][[i]]
  sapply(L, function(cl) {
    pp=vm[clusGm$membership==cl]
    n=length(pp)
    c(mean(log10(dG[pp])),n)
  }) 
})

ddG=pblapply(1:2, function(i){
  L=sigcl[[2]][[i]]
  sapply(L, function(cl) {
    pp=vg[clusGG$membership==cl]
    n=length(pp)
    c(mean(log10(dG[pp])),n)
  }) 
})

dd=unlist(lapply(list(ddM,ddG),function(L) sapply(L,function(l) l[1,])),recursive = F)
names(dd)=c("IgM C","IgM A","IgG C","IgG A")
ns=unlist(lapply(list(ddM,ddG),function(L) sapply(L,function(l) l[2,])),recursive = F)

vioplot(dd, h=0.075, ylab="Mean Log Degree ", names=names(dd))
x=melt(dd)[,1]; cc=melt(dd)[,2]
dusig.test(x,cc)
vioplot(ns, h=700, ylab="Cluster Size ")

x=melt(ns)[,1]; cc=melt(ns)[,2]
dunn.test(x,cc)
x=melt(dd)[,1]; cc=melt(dd)[,2]
dunn.test(x,cc)

tb23=melt(table(fl[,2:3]))
nms23=tb23[,-3]
tb23=tb23[,3]
nms23=apply(as.matrix(nms23),1, function(l) paste(colnames(nms23)[l], collapse = ""))
names(tb23)=nms23
tb23[1]=0
Nm=sum(tb23)

tb45=melt(table(fl[,4:5]))
nms45=tb45[,-3]
tb45=tb45[,3]
nms45=apply(as.matrix(nms45),1, function(l) paste(colnames(nms45)[l], collapse = ""))
names(tb45)=nms45
tb45[1]=0
Ng=sum(tb45)

sigM=pblapply(1:2, function(i){
  L=sigcl[[1]][[i]]
  X=sapply(L, function(cl) {
    pp=vm[clusGm$membership==cl]
    rt=(Nm-length(pp))/length(pp)
    tx=melt(table(fl[pp,2:3]))[,3]
    names(tx)=c("","c","a","ca")
    tx[-1]=rt*tx[-1]/(tb23[-1]-tx[-1])
    return(tx[-1])
  }) 
  colnames(X)=L
  return(X)
})

sigG=pblapply(1:2, function(i){
  L=sigcl[[2]][[i]]
  X=sapply(L, function(cl) {
    pp=vg[clusGG$membership==cl]
    rt=(Ng-length(pp))/length(pp)
    tx=melt(table(fl[pp,4:5]))[,3]
    names(tx)=c("","c","a","ca")
    tx[-1]=rt*tx[-1]/(tb45[-1]-tx[-1])
  }) 
  colnames(X)=L
  return(X)
})
names(sigM)=paste("M",c("C","A"),sep="")
names(sigG)=paste("G",c("C","A"),sep="")

X=unlist(list(sigM,sigG), recursive = F)
for(i in 1:4){
  colnames(X[[i]])=paste(names(X)[i],colnames((X[[i]])), sep="_")
}
Y=X[[1]]
for (i in 2:4) Y=cbind(Y,X[[i]])
sigMG=t(Y)

pdf(file="hmsigMG.pdf", width = 50, height=50)
hmsigMG=heatmap3(sigMG, method="ward.D2", cexCol=5, cexRow = 2,
                 scale="row", margins = c(5,10))
dev.off()

# vioplot(t(sigM[[1]][order(rowMedians(sigM[[1]])),]), h=0.001, las=2, main="IgM - C")
# vioplot(t(sigG[[1]][order(rowMedians(sigG[[1]])),]), h=0.003, las=2, main="IgG - C", ylim=c(0,0.04))
# vioplot(t(sigM[[2]][order(rowMedians(sigM[[2]])),]), h=0.001, las=2, main="IgM - A")
# vioplot(t(sigG[[2]][order(rowMedians(sigG[[2]])),]), h=0.003, las=2, main="IgG - A", ylim=c(0,0.04))

NNM=pblapply(1:2, function(i){
  L=sigcl[[1]][[i]]
  X=sapply(L, function(cl) {
    pp=vm[clusGm$membership==cl]
    n=length(pp)
    y=unlist(sapply(ego(G,1,V(G)[pp]),names))
    y=setdiff(y,pp)
    tx0=melt(table(fl[y,c(2:3)]))[,3]
    tx=melt(table(fl[y,c(4:5)]))[,3]
    tx=tx/tx0
    names(tx)=c("","c","a","ca")
    tx=tx[-1]*Nm/Ng
    return(tx)
  }) 
  colnames(X)=L
  return(X)
})

NNG=pblapply(1:2, function(i){
  L=sigcl[[2]][[i]]
  X=sapply(L, function(cl) {
    pp=vg[clusGG$membership==cl]
    n=length(pp)
    y=unlist(sapply(ego(G,1,V(G)[pp]),names))
    tx0=melt(table(fl[y,c(4:5)]))[,3]
    tx=melt(table(fl[y,c(2:3)]))[,3]
    tx=tx/tx0
    names(tx)=c("","c","a","ca")
    tx=tx[-1]*Nm/Ng
    return(tx)
  })
  colnames(X)=L
  return(X)
})

names(NNM)=paste("M",c("C","A"),sep="")
names(NNG)=paste("G",c("C","A"),sep="")

X=unlist(list(NNM,NNG), recursive = F)
for(i in 1:4){
  colnames(X[[i]])=paste(names(X)[i],colnames((X[[i]])), sep="_")
}
Y=X[[1]]
for (i in 2:4) Y=cbind(Y,X[[i]])
NNMG=t(Y)

pdf(file="hmNNMG.pdf", width = 50, height=50)
hmsigMG=heatmap3(NNMG, method="ward.D2", cexCol=5, cexRow = 2,
                 scale="row", margins = c(5,10))
dev.off()


x1=sigMG
rownames(x1)=paste(rownames(sigMG),"S", sep="_")
x2=NNMG
rownames(x2)=paste(rownames(x2),"NN", sep="_")
X=rbind(x1,x2)
pdf(file="hmsigNNMG.pdf", width = 75, height=75)
hmsigNNMG=heatmap3(X, method="ward.D2", cexCol=5, cexRow = 2,
scale="row", margins = c(5,10))
dev.off()

x1=sigMG
colnames(x1)=paste(colnames(x1),"S", sep="_")
x2=NNMG
colnames(x2)=paste(colnames(x2),"NN", sep="_")
X=cbind(x1,x2)
pdf(file="hmsigNNMGcol.pdf", width = 75, height=75)
hmsigNNMG=heatmap3(X, method="ward.D2", cexCol=5, cexRow = 2,
                   scale="column", margins = c(15,10))
dev.off()

x=as.matrix(fl[,2:3])
x=x[names(dG),]
#dGm=dG[rowAnys(x)]
mdGm=10^(mean(log10(dGm)))

x=as.matrix(fl[,4:5])
x=x[names(dG),]
#dGg=dG[rowAnys(x)]
mdGg=10^(mean(log10(dGg)))

nd=names(dG)
Dsbygroup=aggregate(dG, 
                    by=list(cM=fl[nd,2],aM=fl[nd,3],cG=fl[nd,4],aG=fl[nd,5]), 
                    function(x) 10^(mean(log10(x))))
nmsbg=apply(Dsbygroup[,-5],1, function(l){
  paste(colnames(Dsbygroup)[-5][l],collapse="_")
})
Dbg=Dsbygroup$x; names(Dbg)=nmsbg
par(mar=c(7.5, 4.1, 4.1, 2.1))
barplot(sort(Dbg), las=2, ylab="Mean Degree")

nDbg=names(Dbg)
rownames(Dsbygroup)=nDbg

fld=fl[names(V(G)),-1]
modutypes=apply(Dsbygroup[,-5],1,function(l){
  l=matrix(rep(l,nrow(fld)), nrow=nrow(fld), ncol=4, byrow = T)
  j=1*rowAlls(as.matrix(fld)==l)+1
  modularity(G,j)
})
barplot(sort(modutypes, decreasing = T), las=2, ylab="Modularity")

fldiso=rowSums(fld[,1:2])-rowSums(fld[,3:4])
assortiso=assortativity(G,fldiso)
modiso=modularity(G,fldiso+3)
fldia=rowSums(fld[,c(1,3)])-rowSums(fld[,c(2,4)])
assordia=assortativity(G,fldia)
moddia=modularity(G,fldia+3)

# merged graph clusters ------------------------------------------------------

length(vm)
length(vg)
length(intersect(vm,vg))
length(vm)/1e8*length(vg) #expected overlap
length(vm)/1e9*length(vg)
length(vm)/(length(intersect(vm,vg))/length(vg)) # 4e6 - the whole IgM igome?

# the resolution parameter is selected to yield distribution of cluster sizes
# similar to the degree distribution, so the clusters will have sizes
# in the same range as the sizes of the neighborhoods

clus=cluster_leiden(G, objective_function = "modularity", resolution=70)
tg=melt(table(Cluster=clus$membership,cM=fl[VG,2],aM=fl[VG,3],cG=fl[VG,4],aG=fl[VG,5]))
tg=tg[rowAnys(as.matrix(tg[,2:5])),]
x=apply(tg[,2:5],1,function(l) paste(colnames(tg[,2:5])[l],collapse='_'))
tg=data.frame(tg[,1],x,tg[,6]); colnames(tg)=c("Cluster","Fraction","N")
tg=acast(data=tg,Cluster~Fraction)
j=nchar(colnames(tg))
tg=cbind(tg[,j<8],rowSums(tg[,j>5])); colnames(tg)[11]="Mixed"

chtg=chisq.posthoc.test(tg, simulate.p.value=T)
n=nrow(chtg)/2; jp=2*(1:n); jz=jp-1
chtg=chtg[jz,][rowAnys(chtg[jp,]<0.05),c(1,3:13)]; colnames(chtg)[1]="Cluster"
chtg=chtg[order(chtg[,2]),]

pdf(file="hmclus.pdf", width = 50, height=50)
hmclus=heatmap3(tg, method="ward.D2", cexCol=5, cexRow = 2,
                 scale="row", margins = c(20,10))
dev.off()
pdf(file="hmclusRNK.pdf", width = 50, height=50)
hmclus=heatmap3(t(apply(tg,1,rank)), method="ward.D2", cexCol=5, cexRow = 2,
                scale="row", margins = c(20,10))
dev.off()

vmg=intersect(vm,vg)
Gint=induced_subgraph(G,V(G)[vmg])
components(Gint)
Gint=induced_subgraph(Gint, V(Gint)[components(Gint)$membership==1])
vin=names(V(Gint))

clint=cluster_leiden(Gint, objective_function = "modularity", resolution=18)
tig=melt(table(Cluster=clint$membership,cM=fl[vin,2],aM=fl[vin,3],cG=fl[vin,4],aG=fl[vin,5]))
tig=tig[rowAnys(as.matrix(tig[,2:5])) & tig[,6]>0,]
x=apply(tig[,2:5],1,function(l) paste(colnames(tig[,2:5])[l],collapse='_'))
tig=data.frame(tig[,1],x,tig[,6]); colnames(tig)=c("Cluster","Fraction","N")
tig=acast(data=tig,Cluster~Fraction)
tig[is.na(tig)]=0
j=nchar(colnames(tig))
tig=cbind(tig[,j<8],rowSums(tig[,j>5])); colnames(tig)[5]="Mixed"

chtig=chisq.posthoc.test(tig, simulate.p.value=T)
n=nrow(chtig)/2; jp=2*(1:n); jz=jp-1
chtig=chtig[jz,][rowAnys(chtig[jp,]<0.05),c(1,3:7)]; colnames(chtig)[1]="Cluster"
chtig=chtig[order(chtig[,2]),]



# NN ------------------------------------------------------------------------
tb=melt(table(fl[,2:3]))
nms=tb[,-3]
tb=tb[,3]
nms=apply(as.matrix(nms),1, function(l) paste(colnames(nms)[l], collapse = ""))
names(tb)=nms
tb0=rep(0,)


# too slow
options(future.globals.maxSize= 10000000000)
plan("multisession", workers=14)
proct=proc.time()
NNall=t(pbsapply(VG[1:1000], function(p){
  nn=VG[AjL[[p]]]
  tb=melt(table(fl[nn,]))
  nms=tb[,-3]
  tb=tb[,3]
  nms=apply(as.matrix(nms),1, function(l) paste(colnames(nms)[l], collapse = ""))
  names(tb)=nms
  return(tb[-1])
}))
print(proc.time()-proct)
closeAllConnections()

# SigCl Graph ------------------------------------------------------------

#Gsig0=induced_subgraph(G, V(G)[allpepsig])
x=adjL(allpepsig, fln = "ALsig_")
AL=lapply(list.files(pattern="ALsig_"),function(f) read.delim(file=f, stringsAsFactors = F, header=F))
AL=unlist(AL, recursive = F)
AL=unlist(AL, recursive = F)
AL=lapply(AL, strsplit, split=" ")
AL=unlist(AL, recursive = F)
nms=sapply(AL, function(l) l[1])
names(AL)=nms
AL=lapply(AL,function(l) l[-1])
Gsig=adjL2Glite(AL)

save(Gsig, file="Gsig")
load("Gsig")

#Ls=embed_laplacian_matrix(Gsig,no=1000, which = "sa") # too much
dGsig=degree(Gsig)
median(log10(dGsig))

clG=cluster_leiden(Gsig, objective_function = "modularity", resolution_parameter = 8)
median(log10(table(clG$membership)))
hist(log10(dGsig), xlim=c(0,3))
hist(log10(table(clG$membership)))

Gsig=set_vertex_attr(Gsig, name="cluster", value=clG$membership)
typ=pbapply(fl[names(V(Gsig)),-1],1, function(l){
  rownames(Dsbygroup)[rowAlls(as.matrix(Dsbygroup[,-5])==matrix(rep(l,15), nrow=nrow(Dsbygroup), byrow=T))]
})
Gsig=set_vertex_attr(Gsig, name="type", value=typ)

distr=table(clG$membership,typ)
chsqdstr=chisq.posthoc.test(distr, simulate.p.value=T)
n=nrow(chsqdstr)/2; jp=2*(1:n); jz=jp-1
chsqdstr=chsqdstr[jz,][rowAnys(chsqdstr[jp,]<0.05),]; colnames(chsqdstr)[1]="Cluster"
rownames(chsqdstr)=chsqdstr[,1]
chsqdstr=chsqdstr[,-c(1:2)]

write.graph(Gsig, format="graphml", file="Gsig.graphml")
badcli=as.numeric(names(table(clG$membership)))[table(clG$membership)<20]

Gclsi=pbsapply((2:clG$nb_clusters), function(i1){
        sapply(1:(i1-1), function(i2){
          if ((i1 %in% badcli) | (i2 %in% badcli)) return(c(NA,NA,NA))
          g=induced_subgraph(Gsig, V(Gsig)[clG$membership %in% c(i1,i2)])
          c(i1,i2,modularity(g, vertex_attr(g)$cluster))
    })
})
Gclsi=t(array(unlist(Gclsi, recursive=F), dim=c(3,clG$nb_clusters*(clG$nb_clusters-1)/2)))
Gclsi=Gclsi[!(rowAnys(is.na(Gclsi))),]
X=apply(Gclsi[,1:2],2,as.character)
Gclsig=graph_from_edgelist(as.matrix(X), directed=F)         
Gclsig=set_edge_attr(Gclsig, name="weight", value=Gclsi[,3])
x=mst(Gclsig);mw=max(edge_attr(x)$weight)
Gclsigsm=delete_edges(Gclsig, E(Gclsig)[edge_attr(Gclsig)$weight>mw]) #mw
sz=sapply(names(V(Gclsigsm)), function(i) sum(clG$membership==i))
Gclsigsm=set_vertex_attr(Gclsigsm, name="size", value=sqrt(sz))
ttb0=table(vertex_attr(Gsig)$type)/vcount(Gsig)
typetb=aggregate(vertex_attr(Gsig)$type, by=list(vertex_attr(Gsig)$cluster), 
                 function(x) {
                   t1=table(x)
                   t0=rep(0,length(ttb0));names(t0)=names(ttb0)
                   t0[names(t1)]=t1
                   t0=t0/length(x)
                   t0/ttb0
                 })
typetb=typetb[as.numeric(names(V(Gclsigsm))),]
for (i in 1:15){
 Gclsigsm=set_vertex_attr(Gclsigsm, name=colnames(typetb$x)[i], value=typetb$x[,i]) 
}
Gclsigsm=set_vertex_attr(Gclsigsm, name="Comb_pure", 
              value=vertex_attr(Gclsigsm)$cM+vertex_attr(Gclsigsm)$aM+vertex_attr(Gclsigsm)$cG+vertex_attr(Gclsigsm)$aG) 
j=c(5:7,9,11:17)
X=as.data.frame(vertex_attr(Gclsigsm)[j])
Gclsigsm=set_vertex_attr(Gclsigsm, name="Comb_comb", 
                         value=rowSums(X))

X=typetb$x; rownames(X)=typetb$Group.1
pdf(file="hmtypetb.pdf", width = 15, height=15)
heatmap3(X,  method="ward.D2", margins = c(10,10), scale = "column", cexRow = 0.75)
dev.off()

Lgclsm=embed_laplacian_matrix(Gclsigsm, no=87, type="I-DAD")

plot(Lgclsm$D[75:87])
plot3d(Lgclsm$X[,79:87])
plotG3D(Gclsigsm, coff=5, att="cM", grclr = T)
snapshot3d(filename="Glsigsm_cM.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="aM", grclr = T)
snapshot3d(filename="Glsigsm_aM.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="cG", grclr = T)
snapshot3d(filename="Glsigsm_cG.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="aG", grclr = T)
snapshot3d(filename="Glsigsm_aG.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="Comb_pure")
snapshot3d(filename="Glsigsm_pure.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="Comb_comb")
snapshot3d(filename="Glsigsm_comb.png")

ec=eigen_centrality(Gclsigsm)
Gclsigsm=set_vertex_attr(Gclsigsm, name="EigenCEntr", value=(ec$vector)^4) 
Gclsigsm=set_vertex_attr(Gclsigsm, name="Degree", value=degree(Gclsigsm)) 

plotG3D(Gclsigsm, coff=5, att="EigenCEntr", grclr = T)
snapshot3d(filename="Glsigsm_EiC.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="Degree", grclr = T)
snapshot3d(filename="Glsigsm_Degree.png")
clear3d()

write.graph(Gclsigsm, format="graphml", file="Gclsigsm.graphml")

# degrgee hub ----------------------------------------------------------------

centrHubs=names(degree(Gclsigsm))[degree(Gclsigsm)>60]
centrHubseq=lapply(centrHubs, function(h) names(V(Gsig))[clG$membership==as.numeric(h)])
for (i in seq_along(centrHubs)) print(ggseqlogo(centrHubseq[[i]]))
centrhubspp=unique(unlist(centrHubseq))
vsigHubs=vsigtbl[centrhubspp,]
plot(euler(vsigHubs), quantities=T)

Gclsigsm=set_vertex_attr(Gclsigsm, name="25", 
                         value=1+1*(degree(Gclsigsm)==69))
Gclsigsm=set_vertex_attr(Gclsigsm, name="83", 
                         value=1+1*(degree(Gclsigsm)==74))
clear3d()
plotG3D(Gclsigsm, coff=5, att="25")
snapshot3d(filename="Glsigsm_25.png")
clear3d()
plotG3D(Gclsigsm, coff=5, att="83")

vsigHubs1=vsigtbl[centrHubseq[[1]],]
plot(euler(vsigHubs1), quantities=T)

vsigHubs2=vsigtbl[centrHubseq[[2]],]
plot(euler(vsigHubs2), quantities=T)

AjLm=as_adj_list(Gm)
AjLm=lapply(AjLm, function(l) vm[l])
AjLg=as_adj_list(Gg)
AjLg=lapply(AjLg, function(l) vg[l])

LCSm=pblapply(centrHubseq, function(L){
  X=AjLm[L]
  lapply(seq_along(X),function(i){
    l=X[[i]]; p=names(X)[i]
    s1=unlist(strsplit(p, split=""))
    sapply(l, function(pp){
      s2=unlist(strsplit(pp, split=""))
      paste(LCS(s1,s2)$LCS, collapse="")
    })
  })
})

LCSm=lapply(LCSm, function(L) sort(table(unlist(L)),decreasing = T))


LCSg=pblapply(centrHubseq, function(L){
  X=AjLg[L]
  lapply(seq_along(X),function(i){
    l=X[[i]]; p=names(X)[i]
    s1=unlist(strsplit(p, split=""))
    sapply(l, function(pp){
      s2=unlist(strsplit(pp, split=""))
      paste(LCS(s1,s2)$LCS, collapse="")
    })
  })
})

LCSg=lapply(LCSg, function(L) sort(table(unlist(L)),decreasing = T))

X1=intersect(names(LCSm[[1]]),names(LCSg[[1]]))
X1=cbind(LCSm[[1]][X1],LCSg[[1]][X1])
j=rowProds(X1)
X1=X1[order(j, decreasing = T),];X1=X1[j>10,]
X2=intersect(names(LCSm[[2]]),names(LCSg[[2]]))
X2=cbind(LCSm[[2]][X2],LCSg[[2]][X2])
j=rowProds(X2)
X2=X2[order(j, decreasing = T),];X2=X2[j>10,]

x2=rownames(X2)
ggseqlogo(x2)
x2aa=t(sapply(x2, function(p) unlist(strsplit(p, split=""))))
x2aa=t(apply(x2aa,1,function(l) rev(l)))
colnames(x2aa)=paste("P",1:5, sep="")
z=seqdef(x2aa)
zpst=pstree(z); #zpst.p1=prune(zpst, gain = "G1", C = 1.3, delete = FALSE)
zres2=generate(zpst, 7, 50000, method="pmax")
zmine2=pmine(zpst,zres2)
plot(predict(zpst,zmine2))
zprob2=apply(as.matrix(zmine2),1,paste, collapse="")[predict(zpst,zmine2)>5e-3]
zprob2=sapply(zprob2,function(p) paste(rev(unlist(strsplit(p,split=""))), collapse=""))

#1 - NA

#2 - 
# Query  1   HFGLG  5     2X
#            HFGLG
# Sbjct  7   HFGLG  11
# Query  1   HFGLGQH  7 
#            H+GLG H
# Sbjct  13  HYGLGDH  19
#3 - 
# Query  1   HENVHE  6
#            HENVH+
# Sbjct  11  HENVHD  16

#4 -

x1=rownames(X1)
ggseqlogo(x1)
x1aa=t(sapply(x1, function(p) unlist(strsplit(p, split=""))))
x1aa=t(apply(x1aa,1,function(l) rev(l)))
colnames(x1aa)=paste("P",1:5, sep="")
z=seqdef(x1aa)
zpst=pstree(z); #zpst.p1=prune(zpst, gain = "G1", C = 1.3, delete = FALSE)
zres1=generate(zpst, 7, 50000, method="pmax")
zmine1=pmine(zpst,zres1)
plot(predict(zpst,zmine1))
zprob1=apply(as.matrix(zmine1),1,paste, collapse="")[predict(zpst,zmine1)>5e-3]
zprob1=sapply(zprob1,function(p) paste(rev(unlist(strsplit(p,split=""))), collapse=""))
                                               
# Query  1   KTITSST  7
#            KTIT ST
# Sbjct  6   KTITGST  12


# EHub ----------------------------------------------------------------------

EcentrHub=names(ec$vector)[ec$vector>0.95]
EcentrHubseq=lapply(EcentrHub, function(n) names(V(Gsig))[clG$membership==as.numeric(n)])
for (L in EcentrHubseq) print(ggseqlogo(L))
Ecentrhubspp=lapply(EcentrHubseq, function(n) unique(unlist(n)))
vsigEHub=lapply(EcentrHubseq, function(n) vsigtbl[n,])
plot(euler(vsigEHub[[1]], quantities=T))

X=AjLm[EcentrHubseq[[1]]]
LCSem=lapply(seq_along(X),function(i){
    l=X[[i]]; p=names(X)[i]
    s1=unlist(strsplit(p, split=""))
    sapply(l, function(pp){
      s2=unlist(strsplit(pp, split=""))
      paste(LCS(s1,s2)$LCS, collapse="")
    })
  })

LCSem=sort(table(unlist(LCSem)),decreasing = T)

X=AjLg[EcentrHubseq[[1]]]
LCSeg=lapply(seq_along(X),function(i){
    l=X[[i]]; p=names(X)[i]
    s1=unlist(strsplit(p, split=""))
    sapply(l, function(pp){
      s2=unlist(strsplit(pp, split=""))
      paste(LCS(s1,s2)$LCS, collapse="")
    })
  })

LCSeg=sort(table(unlist(LCSeg)),decreasing = T)

X1=intersect(names(LCSem),names(LCSeg))
X1=cbind(LCSem[X1],LCSeg[X1])

x3=rownames(X1)[nchar(rownames(X1))==5]
x3aa=t(sapply(x3, function(p) unlist(strsplit(p, split=""))))
x35aa=t(apply(x3aa,1,function(l) rev(l)))
colnames(x35aa)=paste("P",1:5, sep="")
z=seqdef(x35aa)
zpst=pstree(z); #zpst.p1=prune(zpst, gain = "G1", C = 1.3, delete = FALSE)
zres3=generate(zpst, 7, 50000, method="pmax")
zmine3=pmine(zpst,zres3)
plot(predict(zpst,zmine3))
zprob3=apply(as.matrix(zmine3),1,paste, collapse="")[predict(zpst,zmine3)>1e-3]
zprob3=sapply(zprob3,function(p) paste(rev(unlist(strsplit(p,split=""))), collapse=""))

ggseqlogo(zprob3)
#1 - 
#2 - 
# Query  1   GSNSMN  6
#            GSNSM+
# Sbjct  13  GSNSMD  18
#3 -
#4 - 
#5 -
#6 -
#7 -
# Query  1  MSNSM  5
#           MSNSM
# Sbjct  5  MSNSM  9
#8 -
# Query  1   ISNSMN  6
#            ISN MN
# Sbjct  6   ISNGMN  11


