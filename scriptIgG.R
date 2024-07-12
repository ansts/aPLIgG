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

options(future.globals.maxSize= 10^10)

source("clusterGraph.R")
source("corrD.R")
source("grEmbed.R")
source("adjL.R")
source("adjL2Glite.R")
# Data Acquisition and Cleaning ----------------------------------------------

cpl3=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))
cpl1=colorRampPalette(c("#0000FFFF","#00FF00FF","#FFFF00FF","#FF0000FF"), alpha=T)
aa=AA_ALPHABET[1:20]
load("mixIgM")
load("flIgM")
mixIgM=mix
flIgM=fl
rm(mix, fl)
aa=AA_STANDARD

cntr_all=lapply(list.files(pattern="*.cntr.txt"),function(fl) read.table(file=fl, header = FALSE,  quote=NULL, stringsAsFactors = F))
aPL_all=lapply(list.files(pattern="*T.txt"),function(fl) read.table(file=fl, header = FALSE,  quote=NULL, stringsAsFactors = F))

aPL=aPL_all[[1]]
for (i in 2:length(aPL_all)){
  aPL=rbind(aPL,aPL_all[[i]])
}
xa=aggregate(aPL[,2], by=list(as.character(aPL[,1])), FUN=sum)
xa1=xa[xa[,2]==1,1]
xa2=xa[xa[,2]==2,1]
xa3up=xa[xa[,2]>2&xa[,2]<101,1]

cntr=cntr_all[[1]]
for (i in 2:length(cntr_all)){
  cntr=rbind(cntr,cntr_all[[i]])
}

xc=aggregate(cntr[,2], by=list(as.character(cntr[,1])), FUN=sum)
xc1=xc[xc[,2]==1,1]
xc2=xc[xc[,2]==2,1]
xc3up=xc[xc[,2]>2&xc[,2]<101,1]

# cross-multiplicates

xac12=intersect(xa1,xc2)
xac21=intersect(xa2,xc1)
xac123=intersect(union(xa1,xa2),xc3up)
xac321=intersect(union(xc1,xc2),xa3up)
xam123=intersect(union(xa1,xa2),mixIgM)
xcm123=intersect(union(xc1,xc2),mixIgM)

aPL=unique(c(xa3up,xac12,xac21,xac123,xam123))
cntr=unique(c(xc3up,xac12,xac21,xac321,xcm123))

taPL=table(xa[,2])
xi=log10(as.double(names(taPL)))
yi=log10(taPL)
plot(xi,yi)

tcntr=table(xc[,2])
xi=log10(as.double(names(tcntr)))
yi=log10(tcntr)
plot(xi,yi)

ix=grep("X",cntr)
write(cntr, file="cntr.txt")
write(aPL, file="aPL.txt")

aPLTUPS=read.table(file="aPLTUP.txt", header = T, sep="\t",quote = NULL, stringsAsFactors = F)
excpt=unique(aPLTUPS[,3])
aPLbad=aPLTUPS[!(aPLTUPS[,3] %in% excpt[c(4,5,8)]),1]
aPL=setdiff(aPL,aPLbad)
cntrTUPS=read.table(file="cntrTUP.txt", header = T, sep="\t",quote = NULL, stringsAsFactors = F)
excpt=unique(cntrTUPS[,3])
cntrbad=cntrTUPS[!(cntrTUPS[,3] %in% excpt[c(3,5)]),1]
cntr=setdiff(cntr,cntrbad)
write(cntr, file="cntr.txt")
write(aPL, file="aPL.txt")

X=unique(c(aPL,cntr))

# from here on mix will denote the mixture of 
# the IgM mix database and the IgG results

mix=union(mixIgM,X)
mxcM=mixIgM[flIgM %in% c(1,3,5,7)]
mxaM=mixIgM[flIgM %in% c(2,3,6,7)]
mxMM=mixIgM[flIgM %in% c(4,5,6,7)]
fl=data.frame(Pub=mix %in% mxMM, cM=mix %in% mxcM,aM=mix %in% mxaM, 
                                  cG=mix %in% cntr, aG=mix %in% aPL)
rownames(fl)=mix
save(mix, file="mixnew")
save(fl, file="flnew")
rm(list=ls())

# Venn -------------------------------------------------------------------


load("mixnew")
load("flnew")

plot(euler(fl, shape="ellipse"), quantitites=T)

tbl=reshape2::melt(table(fl))
szs=colSums(fl)
tblprop=apply(tbl,1,function(l) {
  y=unlist(l[1:(length(l)-1)])
  x=szs[as.logical(y)]
  x=rep(l[[length(l)]],length(x))/x
  if(any(y)) return(list(x))  else return(0)
})
sink(file="venn_table.txt")
print(tblprop)
sink()
write.csv(tbl,file="venn_table_raw.csv")

x=fl[,-1]
plot(venn(x), quantities=T)
colSums(fl)

# The Graph ---------------------------------------------------------------
## Construction -----

load("mixnew")
load("flnew")
source("adjL.R")
source("adjL2Glite.R")
source("treeCrawLouv.R")

save(mix,file="varAdjL")

X=mix[fl$cM|fl$aM|fl$cG|fl$aG]
x=adjL(X)
AL=lapply(list.files(pattern="AL_ "),function(fl) read.delim(file=fl, stringsAsFactors = F, header=F))
AL=unlist(AL, recursive = F)
AL=unlist(AL, recursive = F)
AL=lapply(AL, strsplit, split=" ")
AL=unlist(AL, recursive = F)
nms=sapply(AL, function(l) l[1])
names(AL)=nms
AL=lapply(AL,function(l) l[-1])
save(AL,file="AL")

load("AL")
nms=names(AL)
nmn=seq_along(nms)
names(nmn)=nms

cl=makeCluster(14)
clusterExport(cl, "nmn")
Al=pblapply(AL,function(l) nmn[l], cl=cl)
stopCluster(cl)
names(Al)=NULL
Al=lapply(Al, function(l) {
  names(l)=NULL
  return(l)
})
G=graph_from_adj_list(adjlist=Al, mode="all")
G=set_vertex_attr(G, "name", value=nms)
x=apply(fl,1, function(l) paste(1*l,collapse=""))
G=set_vertex_attr(G, "Group", value=x[names(V(G))])
E(G)$weight=1
save(G, file="G")
write.graph(G, format = "graphml", file="G.graphml")

loners=names(V(G))[components(G)$membership>1]
fl[loners,]

# the disconnected sequences are mimotopes of APL associated antibodies 
# 5/8 from the IgG sample and 2/8 from the IgM sample. 1 is from the control IgM
# which is found in the public repertoire too as is one of the IgG's.
# Note the KCCVFQV sequence similar to a motif of 2G12 mimotopes.

load("G")
G=induced_subgraph(G, components(G)$membership==1)

es=ends(G,E(G))

ednd1=strsplit(es[,1],split="")
ednd1=t(sapply(ednd1, unlist))
ednd2=strsplit(es[,2],split="")
ednd2=t(sapply(ednd2, unlist))
endp=cbind(ednd1,ednd2)
rm(es,ednd1,ednd2)
ij=combn(7,5)

j=cut(1:nrow(endp),30,labels=F)
edp=sapply(unique(j), function(i) {
  M=endp[j==i,]
  cl=makeCluster(14)
  clusterExport(cl,c("ij","M"))
  X=pbapply(M,1,function(l){
    p1=l[1:7]
    p2=l[8:14]
    l1=nrow(unique(t(apply(ij,2,function(ji) (p1[ji])))))
    l2=nrow(unique(t(apply(ij,2,function(ji) (p2[ji])))))
    l1*l2
  },cl=cl)
  closeAllConnections()
  gc()
  return(X)
})


cl=makeCluster(16)
clusterEvalQ(cl, require(qualV))
lcsG=pbapply(es,1,function(eg) {
  eg=strsplit(eg,split="")
  res=paste(LCS(eg[[1]],eg[[2]])$LCS, collapse="")
  return(res)
}, cl=cl)
stopCluster(cl)

G=set_edge_attr(G, name="LCS",value=lcsG)

dG=degree(G)
ccG=transitivity(G, "localundirected")
ecG=eigen_centrality(G)

save(G,file="Ggiant")

## Cluster -----

load("Ggiant")
G70=clusterGraph(G, resol=70)
load("Gctrpeps")
G70peps=Gpeps
save(G70peps, file="G70peps")
save(G70, file="G70")
G=NULL
rm(G)

## Map Idiotypes -----------------------------------------------------------

load("IgJT")
v=names(dG)
# As qgrams - only complete 7-mers
qq7I=t(qgrams(v,IgJtrim,q=7))
qq7I=qq7I[v,2]
save(qq7I,file="qq7I")

## Add attributes ------

G70gr=t(sapply(vertex_attr(G70)$Group, unlist))
rownames(G70gr)=vertex_attr(G70)$name

n=names(dG)
A=cbind(fls[n,],bkgNN[n],qq7I[n])
colnames(A)[6:7]=c("Bkg","Id")
cl=makeCluster(14)
clusterExport(cl,"A")
AttrG=t(pbsapply(G70peps, function(x) colSums(A[x,]), cl=cl))
stopCluster(cl)
AttrG[,1:6]=scale(sweep(AttrG[,1:6],1,lengths(G70peps),"/"))

vertex_attr(G70)=cbind(as.data.frame(vertex_attr(G70)[2:4]),as.data.frame(AttrG))

## Stats ----

X=G70gr[,1:15]+G70gr[,16:30]
X[,7]=rowSums(X[,c(7,11,13:15)])
X=X[,-c(11,13:15)]

chsqX=chisq.posthoc.test(X, simulate.p.value = T)
j2=(1:(nrow(chsqX)/2))*2
j1=j2-1
G_good=(chsqX[j2,3:13]<0.05)*sign(chsqX[j1,3:13])
rownames(G_good)=chsqX$Dimension[j2]

for (n in colnames(G_good)) G70=set.vertex.attribute(G70,name=n,value=as.data.frame(G_good)[,n])

## Logos of Clusters -----------------------------------------------------------

logosG70=lapply(G70peps,ggseqlogo)
save(logosG70,file="LOGOSG70") 

for (i in 1:100) print(logosG70[[i]])

# Background frequencies ------------

load("matpep")
matpep=sample(matpep, length(V(Gg)))
x=adjL(matpep)

# AL=lapply(list.files(pattern="AL_ "),function(fl) read.delim(file=fl, stringsAsFactors = F, header=F))
# AL=unlist(AL, recursive = F)
# AL=unlist(AL, recursive = F)
# AL=lapply(AL, strsplit, split=" ")
# AL=unlist(AL, recursive = F)
nms=names(x)
AL=lapply(x,function(l) l[-1])
save(AL,file="AdjLmtckh")
nmn=seq_along(nms)
names(nmn)=nms

cl=makeCluster(14)
clusterExport(cl, "nmn")
Al=pblapply(AL,function(l) {
  l=nmn[l]
  names(l)=NULL
  return(l)
}, cl=cl)
stopCluster(cl)
names(Al)=NULL

Gmtchk=graph_from_adj_list(adjlist=Al, mode="all")
Gmtchk=set_vertex_attr(Gmtchk, "name", value=nms)
E(Gmtchk)$weight=1
Gbg32sm=clusterGraph(Gmtchk)

ij=cut(seq_along(vG),1600,labels = F)

cl=makeCluster(16)
clusterExport(cl,c("vG","matpep","ij"), envir = environment())
clusterEvalQ(cl, require(stringdist))
bkgNN=pbsapply(1:1600, function(i) {
  sapply(vG[ij==i], function(p){
    x=stringdist(p, matpep, method = "lcs", nthread = 1 )
    sum(x<5)
  })
},cl=cl)
stopCluster(cl)

# Isotype graphs -----

load("Ggiantw")
load("flnew")

fls=fl[names(V(G)),]

## IgM ----
vG=names(V(G))
vG=rownames(fls)[rownames(fls) %in% vG]
flss=fls[vG,]
gm=flss$cM|flss$aM
Gmw=induced.subgraph(G,gm)
save(Gmw, file="Gmw")
write.graph(Gmw, format = "graphml", file = "Gmw.graphml")

Gm39w=clusterGraphMod(Gmw, resol=39)
load("Gctrpeps")
Gm39peps=Gpeps
save(Gm39peps, file="Gm39peps")
rm(Gmw)

grps=names(vertex_attr(Gm39w)$Group[[1]])
cim=as.numeric(substr(grps,2,3))==10
aim=as.numeric(substr(grps,2,3))==1
caim=as.numeric(substr(grps,2,3))==11
d=cim+aim+2*caim

Gm39CA=t(sapply(vertex_attr(Gm39w)$Group, function(x){
  x=as.data.frame(x)
  cntr=sum(x[cim,])
  apl=sum(x[aim,])
  all=sum(x[caim,])
  return(c(cntr, apl, all))
}))
colnames(Gm39CA)=c("Control","APLS","All")
rownames(Gm39CA)=names(V(Gm39w))

y=Gm39CA
chsqmCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=(1:(nrow(chsqmCA)/2))*2
j1=j2-1
Gmw_good=(chsqmCA[j2,3:5]<0.05)*sign(chsqmCA[j1,3:5])
rownames(Gmw_good)=chsqmCA$Dimension[j2]

jc=chsqmCA[j1,1][rowSums(abs(chsqmCA[j1,3:5])<0.25)==3]
negm=Gm39peps[jc]

Gm39P=sapply(vertex_attr(Gm39w)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})


## IgG ----
vG=names(V(G))
vG=rownames(fls)[rownames(fls) %in% vG]
flss=fls[vG,]
gg=flss$cG|flss$aG
Ggw=induced.subgraph(G,gg)
save(Ggw, file="Ggw")

Gg40w=clusterGraphMod(Ggw, resol=40)
load("Gctrpeps")
Gg40peps=Gpeps
save(Gg40peps, file="Gg45peps")
rm(Ggw)

grps=names(vertex_attr(Gg40w)$Group[[1]])
cig=as.numeric(substr(grps,4,5))==10
aig=as.numeric(substr(grps,4,5))==1
caig=as.numeric(substr(grps,4,5))==11
d=cig+aig+2*caig

Gg40CA=t(sapply(vertex_attr(Gg40w)$Group, function(x){
  x=as.data.frame(x)
  cntr=sum(x[cig,])
  apl=sum(x[aig,])
  all=sum(x[caig,])
  return(c(cntr, apl, all))
}))
colnames(Gg40CA)=c("Control","APLS","All")
rownames(Gg40CA)=names(V(Gg40w))

y=Gg40CA
chsqgCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=(1:(nrow(chsqgCA)/2))*2
j1=j2-1
Ggw_good=(chsqgCA[j2,3:5]<0.05)*sign(chsqgCA[j1,3:5])
rownames(Ggw_good)=chsqgCA$Dimension[j2]

jc=chsqgCA[j1,1][rowSums(abs(chsqgCA[j1,3:5])<0.25)==3]
negg=Gg40peps[jc]

Gg40P=sapply(vertex_attr(Gg40w)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})

load("qq7I")
qq7Im=sapply(Gm39peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Im=log10(qq7Im+1e-5)
qq7Ig=sapply(Gg40peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Ig=log10(qq7Ig+1e-5)

Gm39w=set.vertex.attribute(Gm39w,name="C",value=Gm39CA[,1]/vertex_attr(Gm39w)$size)
Gm39w=set.vertex.attribute(Gm39w,name="C_sig",value=Gmw_good[,1])
Gm39w=set.vertex.attribute(Gm39w,name="A",value=Gm39CA[,2]/vertex_attr(Gm39w)$size)
Gm39w=set.vertex.attribute(Gm39w,name="A_sig",value=Gmw_good[,2])
Gm39w=set.vertex.attribute(Gm39w,name="all_sig",value=Gmw_good[,3])
Gm39w=set.vertex.attribute(Gm39w,name="Pub",value=Gm39P/vertex_attr(Gm39w)$size)
Gm39w=set.vertex.attribute(Gm39w,name="Id",value=qq7Im)

Gg40w=set.vertex.attribute(Gg40w,name="C",value=Gg40CA[,1]/vertex_attr(Gg40w)$size)
Gg40w=set.vertex.attribute(Gg40w,name="C_sig",value=Ggw_good[,1])
Gg40w=set.vertex.attribute(Gg40w,name="A",value=Gg40CA[,2]/vertex_attr(Gg40w)$size)
Gg40w=set.vertex.attribute(Gg40w,name="A_sig",value=Ggw_good[,2])
Gg40w=set.vertex.attribute(Gg40w,name="all_sig",value=Ggw_good[,3])
Gg40w=set.vertex.attribute(Gg40w,name="Pub",value=Gg40P/vertex_attr(Gg40w)$size)
Gg40w=set.vertex.attribute(Gg40w,name="Id",value=qq7Ig)
save(Gm39w, file="Gm39w")
save(Gg40w, file="Gg40w")

write_graph(Gm39w,format = "graphml", file="Gm39w.graphml")
write_graph(Gg40w,format = "graphml", file="Gg40w.graphml")

# Fix corrD, anyway - not really useful
# d=round(corrD(Gm39w))
# x=grEmbed(Gm39w, att=c(5,7,10,11),  k=d)
# write_graph(x,format = "graphml", file="Gm39we.graphml")
# d=round(corrD(Gg40w))
# x=grEmbed(Gg40w, att=c(5,7,10,11), k=d)
# write_graph(x,format = "graphml", file="Gg40we.graphml")

# Fusion graph -----
load("Ggiantw")
vGg=names(V(G))
allposcl=list(Gm39peps[rowSums(abs(Gmw_good[,1:2]))>0],
              Gg40peps[rowSums(abs(Ggw_good[,1:2]))>0])
allpos=unique(unlist(allposcl))

ubiqp=names(which(rowSums(fl)==5))
crosubiq=unlist(lapply(ego(G, 1, ubiqp), names))
allpos=setdiff(allpos, crosubiq)

Gpos=induced.subgraph(G,V(G)[allpos])
Gpos=induced.subgraph(Gpos,components(Gpos)$membership==1)
ecGp=eigen_centrality(Gpos)
dGp=degree(Gpos)
pepGpos=names(V(Gpos))
rm(G)
closeAllConnections()


origcl=pbsapply(pepGpos, function(p){
  i=grep(p, allposcl)
  sapply(i,function(j) grep(p,allposcl[[j]]))
})

save(Gpos, file="Gpos")

x=modCluScan(Gpos,seq(1,60,0.5), ncore=15)
plot(x, ty="b", xlab="Resolution", ylab="Normalized modularity")
plot(smooth.spline(x,spar=0.35), ty="l", xlab="Resolution", ylab="Normalized modularity")

Gposm1=clusterGraph(Gpos, resol=27)
load("Gctrpeps")
Gposmpeps=Gpeps
save(Gposm, file="Gposm")
save(Gposmpeps,file="Gposmpeps")
ecGpsm=sapply(Gposmpeps, function(p){
  g=induced_subgraph(Gpos,V(Gpos)[names(V(Gpos)) %in% p])
  names(V(g))[which.max(eigen_centrality(g)$vector)]
})


## Gpos cluster signif ----
ns=lengths(Gposmpeps)
grpsGposm=t(sapply(Gposmpeps, function(pp){
  colSums(fl[pp,])
}))
justgr=grpsGposm[,2:5]
chsqgrpsGpsm=chisq.posthoc.test(justgr, simulate.p.value=T)
chsqgrpsGpsm_V=chsqgrpsGpsm[(1:(nrow(chsqgrpsGpsm)/2))*2-1,]
rownames(chsqgrpsGpsm_V)=chsqgrpsGpsm_V[,1]
chsqgrpsGpsm_V=chsqgrpsGpsm_V[,-(1:2)]
chsqgrpsGpsm_p=chsqgrpsGpsm[(1:(nrow(chsqgrpsGpsm)/2))*2,]
rownames(chsqgrpsGpsm_p)=chsqgrpsGpsm_p[,1]
chsqgrpsGpsm_p=chsqgrpsGpsm_p[,-(1:2)]

justPub=data.frame(Pub=grpsGposm[,1], Not=rowSums(grpsGposm[,2:5])-grpsGposm[,1])
chsqPubGpsm=chisq.posthoc.test(justPub, simulate.p.value=T)
chsqPubGpsm_V=chsqPubGpsm[(1:(nrow(chsqPubGpsm)/2))*2-1,]
rownames(chsqPubGpsm_V)=chsqPubGpsm_V[,1]
chsqPubGpsm_V=chsqPubGpsm_V[,-(1:2)]
chsqPubGpsm_p=chsqPubGpsm[(1:(nrow(chsqPubGpsm)/2))*2,]
rownames(chsqPubGpsm_p)=chsqPubGpsm_p[,1]
chsqPubGpsm_p=chsqPubGpsm_p[,-(1:2)]

Gposm=set.vertex.attribute(Gposm, name=colnames(chsqPubGpsm_V)[1], value=chsqPubGpsm_V[,1])
for (i in 1:4){
  Gposm=set.vertex.attribute(Gposm, name=colnames(chsqgrpsGpsm_V)[i], value=chsqgrpsGpsm_V[,i])
}

M=cbind(chsqgrpsGpsm_p<0.05&chsqgrpsGpsm_V>0,chsqgrpsGpsm_p<0.05&chsqgrpsGpsm_V<0)
colnames(M)=paste(colnames(M),c(rep("Up",4),rep("Down",4)), sep="_")

vennGpsm=euler(M, shape="ellipse", control=list(extraopt_threshold=0.0001))
plot(vennGpsm, legend=T, quantiities=T, adjust_labels=T)

l=seq_along(chsqgrpsGpsm_p[,1])
irnd=sapply(1:4,function(j) l %in% sample(l,sum(chsqgrpsGpsm_p[,j])))
irnd=cbind(irnd,sapply(1:4,function(j) l %in% unique(sample(which(!(l %in% l[irnd[,j]])),sum(chsqgrpsGpsm_p[,j]), replace=T))))
colnames(irnd)=colnames(M)
vennRnd=euler(irnd, shape="ellipse", control=list(extraopt_threshold=0.0001))
vennRnd
plot(vennRnd, legend=F, quantiities=T, adjust_labels=T)

corrplot(cor(chsqgrpsGpsm_V), method="color")

#Gposmemb=grEmbed(Gposm,att=4:9)
XY=embed_laplacian_matrix(Gposm,no=50)$X
Gposmemb=Gposm
Gposmemb=set_vertex_attr(Gposmemb, name="X", value = XY[,1])
Gposmemb=set_vertex_attr(Gposmemb, name="Y", value = XY[,2])

write.graph(Gposmemb,file="Gposmemb.graphml", format="graphml")

Vsig=sapply(1:4,function(j){
    x=chsqgrpsGpsm_V[,j]
    y=chsqgrpsGpsm_p[,j]
    range(x[y>=0.05])
})
colnames(Vsig)=colnames(chsqgrpsGpsm_p)
table(as.data.frame(chsqgrpsGpsm_V[rowSums(chsqgrpsGpsm_p<0.05)==4,]>0))
table(rowSums(chsqgrpsGpsm_p<0.05))
clgood=rownames(chsqgrpsGpsm_p)[!(rowProds((chsqgrpsGpsm_p<0.05)[,c(1,3)]==(chsqgrpsGpsm_p<0.05)[,c(2,4)])==1)]
Gposmppfr=Gposmpeps[clgood]
ecGpsm_good=ecGpsm[clgood]

## Entropy (not calculated last) -----

vGgaa=strsplit(vGg, split="")
En=pbsapply(vGgaa, function(p){
  x=table(p)
  entropy(x, method="ML")
})

names(En)=vGg
vioplot(En~as.factor((names(En) %in% unlist(Gposmpeps))))

vGmx=melt(Gm39peps)
vGmxaa=strsplit(vGmx$value, split="")
Enm=pbsapply(vGmxaa, function(p){
  x=table(p)
  round(entropy(x, method="ML"),2)
})
names(Enm)=vGmx$value
Gm_ggr=(apply(Gmw_good, 1,paste, collapse=""))
Gm_ggr[Gm_ggr=="0-11"]=NA
Gm_ggr[Gm_ggr=="100"]="1-10"
Gm_ggr=as.factor(Gm_ggr[!is.na(Gm_ggr)])
x=Gm39peps[names(Gm39peps) %in% names(Gm_ggr)]
x=melt(x)
x=cbind(x,Gm_ggr[x$L1])
y=table(Gm_ggr)
yalf=rep(1,nrow(x))/y[(x$`Gm_ggr[x$L1]`)]
yalf=yalf/max(yalf)
names(yalf)=NULL
vioplot(Enm[x$value]~x$`Gm_ggr[x$L1]`, col="lightgrey")
stripchart((Enm[x$value]+runif(nrow(x),-0.033,0.033))~x$`Gm_ggr[x$L1]`, 
           method="jitter",vertical=T,add=T, pch=16, col=rgb(0,0.4,1,1), 
           cex=0.01)
distrat=sapply(names(y), function(n) {
  x=table(Enm[x$value[x$`Gm_ggr[x$L1]`==n]])
  x=x/sum(x)
  x0=table(Enm)
  x0=x0/sum(x0)
  x1=rep(0,length(x0))
  names(x1)=names(x0)
  x1[names(x)]=x
  return(x1/x0)
})

cats=c("Under in Control","Over in APLS and Under in Control","Under in APLS","NS","Over in Ubiquitous","Over in Control and Under in APLS")
cats=paste("IgM - ", cats, sep="")
names(cats)=colnames(distrat)

sapply(colnames(distrat),function(j) {
  co=distrat[,j]
  barplot(log10(co+0.07),ylim=range(log10(distrat+0.07)), main=cats[j], xlab="Entropy in nats", ylab="Log fold change relative to the whole set")
})

# and for IgG

vGgx=melt(Gg40peps)
vGgxaa=strsplit(vGgx$value, split="")
Eng=pbsapply(vGgxaa, function(p){
  x=table(p)
  round(entropy(x, method="ML"),2)
})
names(Eng)=vGgx$value
Gg_ggr=(apply(Ggw_good, 1,paste, collapse=""))
Gg_ggr[Gg_ggr=="-100"]="-110"
Gg_ggr[Gg_ggr=="100"]="1-10"
Gg_ggr=as.factor(Gg_ggr[!is.na(Gg_ggr)])
x=Gg40peps[names(Gg40peps) %in% names(Gg_ggr)]
x=melt(x)
x=cbind(x,Gg_ggr[x$L1])
y=table(Gg_ggr)
yalf=rep(1,nrow(x))/y[(x$`Gg_ggr[x$L1]`)]
yalf=yalf/max(yalf)
names(yalf)=NULL
vioplot(Eng[x$value]~x$`Gg_ggr[x$L1]`, col="lightgrey")
stripchart((Eng[x$value]+runif(nrow(x),-0.033,0.033))~x$`Gg_ggr[x$L1]`, 
           method="jitter",vertical=T,add=T, pch=16, col=rgb(0,0.4,1,1), 
           cex=0.01)
distratg=sapply(names(y), function(n) {
  x=table(Eng[x$value[x$`Gg_ggr[x$L1]`==n]])
  x=x/sum(x)
  x0=table(Eng)
  x0=x0/sum(x0)
  x1=rep(0,length(x0))
  names(x1)=names(x0)
  x1[names(x)]=x
  return(x1/x0)
})
catsg=c("Over in APLS and Under in Control","Under in APLS","NS","Over in Ubiquitous","Over in APLS","Over in Control and Under in APLS")
catsg=paste("IgG - ", catsg, sep="")
names(catsg)=colnames(distratg)

sapply(colnames(distratg),function(j) {
  co=distratg[,j]
  barplot(log10(co+0.01),ylim=range(log10(distratg+0.01)), main=catsg[j], xlab="Entropy in nats", ylab="Log fold change relative to the whole set")
})

## Line graph fusion ----
X=induced.subgraph(Gpos, V(Gpos)[unlist(Gposmppfr)])
lGpos=make_line_graph(X)
lGpos=set_vertex_attr(lGpos, name="lcs", value=unlist(edge_attr(Gpos)$LCS))
lGpos=set_vertex_attr(lGpos, name="N", value=1)
lGpos=set_vertex_attr(lGpos, name="W", value=edge_attr(Gpos)$weight)
es=ends(lGpos, E(lGpos))
wt=vertex_attr(lGpos)$W
names(wt)=names(V(lGpos))
wt=rowMeans(array(wt[es], dim=dim(es)))
lGpos=set_edge_attr(lGpos, name="weight", value=wt)
x=as.numeric(as.factor(vertex_attr(lGpos)$lcs))
lGpos=contract(lGpos,x, vertex.attr.comb = list(N ="sum", lcs="first", W="sum"))
lGpos=simplify(lGpos, edge.attr.comb = "sum")
lGpos=set_vertex_attr(lGpos, name="name", value=vertex_attr(lGpos)$lcs)
lGpos=delete_vertex_attr(lGpos,"lcs")
save(lGpos,file="lGpos")

dlG=degree(lGpos)
vlG=names(V(lGpos))
eclG=eigen_centrality(lGpos)
lGpos=set_vertex_attr(lGpos, name="size", value=log10(vertex_attr(lGpos)$N)*10+10)

lGpsm=clusterGraph_l(lGpos, resol=6)
load("Gctrpeps_l")
lGpsmpeps=Gpeps
x=edge_attr(lGpsm)$weight
x=log10(x)
x=(x-min(x)+0.1)/diff(range(x))
lGpsm=set_edge_attr(lGpsm, name="weight",value=x)
lGpsm=set_vertex_attr(lGpsm, name="sizeR",value=sapply(vertex_attr(lGpsm)$N, function(x) sqrt(sum(x))))
save(lGpsm, file="lGpsm")
write.graph(lGpsm, format = "graphml", file="lGpsm.graphml")

for (j in seq_along(lGpsmpeps)){
  p=lGpsmpeps[[j]]
  seq2fasta(S=p,fnm=paste("lGpsmpeps/pep",j,collapse = "_"))
}

## complexity comparison ----
plan("multisession", workers=10)
complx=future_sapply(vlG,function(p) mean(table(unlist(strsplit(p,split="")))))
complx_deg=cbind(complx,dlG)
complx_N=cbind(complx,N=vertex_attr(lGpos)$N)
boxplot(data=complx_deg,log10(dlG)~complx,notch=T)
boxplot(data=complx_N,log10(N)~complx,notch=T)

complx=future_sapply(unlist(lGpsmpeps),function(p) mean(table(unlist(strsplit(p,split="")))))
complx_deg=cbind(complx,dlGt=degree(lGpsm))
complx_N=cbind(complx,N=vertex_attr(lGpsm)$N)
boxplot(data=complx_deg,log10(dlGt)~complx,notch=T)
boxplot(data=complx_N,log10(N)~complx,notch=T)

# There is some negative correlation between the mean degree or the 
# mean N in the lGpos and the entropy of the lcs-s for the length 5 
# but it is lost in lGpsm.

## Eigencentrality based selection ----

nx=vertex_attr(lGpos)$N
plan("multisession", workers=16)
names(nx)=names(V(lGpos))
ecchng=future_lapply(lGpsmpeps, function(l) {
  if (length(l)>1) {
    eclst=eigen_centrality(induced.subgraph(lGpos,l))$vector
    x0=names(which.max(eclst))
  } else {
    x0=l
  }
  x=names(ego(induced.subgraph(lGpos,l),1,x0)[[1]])
  xn=sort(nx[x], decreasing = T)
  thr=sapply(seq_along(xn), function(i) sum(xn[1:i])/sum(xn))
  return(list(x0,xn[thr<0.33]))
})
closeAllConnections()
lcs_g=sapply(ecchng, function(l) names(l[[2]]))
lcs1_g=sapply(ecchng, function(l) l[[1]])

dx=distances(lGpos,lcs1_g,lcs1_g, weights=NA)
range(dx[dx>0])
diag(dx)=100
table(rowMins(dx))
ecmaxntwk=components(induced_subgraph(lGpos,lcs1_g))$membership[lcs1_g]

logos_lGpsm5=lapply(lGpsmpeps,function(s) {
  ggseqlogo(s[nchar(s)==5])
})
save(logos_lGpsm5,file="logos_lGpsm5") 

pdf(file="logos_lGpsm.pdf", width=5,height = 5)
  for (i in seq_along(lGpsmpeps)) plot(logos_lGpsm5[[i]])
dev.off()

logos_lGpsm6=lapply(lGpsmpeps,function(s) {
  if (any(nchar(s)==6)) ggseqlogo(s[nchar(s)==6]) else NULL
})
save(logos_lGpsm6,file="logos_lGpsm6") 

pdf(file="logos_lGpsm6.pdf", width=5,height = 5)
for (i in seq_along(lGpsmpeps)) plot(logos_lGpsm6[[i]])
dev.off()

## connectivity study ----

ownpeps=future_lapply(Gposmpeps, function(l) {
  ex=ego(Gpos,1,l)
  L=length(l)
  prp=sapply(ex, function(ei){
    length(setdiff(names(ei)[-1],l))/L
  })
  return(prp)
})

y=table(sapply(ownpeps,function(x) sum(x==0)))

plot(lengths(Gposmpeps)+0.5, sapply(ownpeps,function(x) sum(x==0))+0.5, log="xy")
pepins=unlist(lapply(seq_along(ownpeps),function(i) Gposmpeps[[i]][ownpeps[[i]]==0]))

boxplot(log10(ecGp$vector[pepins]+0.5),log10(ecGp$vector[!(names(V(Gpos)) %in% pepins)]+0.5), notch=T)
ecGp$vector[pepins]
f=ecdf(ecGp$vector)
hist(log10(f((ecGp$vector[pepins]))))
10^(mean(log10(f((ecGp$vector[pepins])))))
10^(mean(log10(f((ecGp$vector)))))
hist(log10(degree(Gpos)[pepins]), xlim=c(0,3), col=rgb(1,0,0,0.3))
hist(log10(degree(Gpos)))

# Only ?20% of the clusters have some strictly inside sequences, but one has 429. To
# have non-zero probability of an inside peptide, the size of the cluster should
# be more than 40 and to have 100% probability for at least one - more than
# 1000. The inside peptides also have on the average about 10 times lower
# degree. 
# Those parameters were derived from a graph without edge weight
# adjustment to compensate the dependence of the degree on entropy. After that adjustment
# the new Gpos and its derivatives yield positive clusters 24% of which have some
# strictly inside sequences, the maximum is 79. To
# have a non-zero probability of an inside peptide, the size of the cluster should
# be more than 100 and to have 100% probability for at least one - more than
# 1010.

ownpepsflat=unlist(lapply(seq_along(ownpeps), function(i) {
  p=ownpeps[[i]]
  names(p)=Gposmpeps[[i]]
  return(p)
}))

plot(ecGp$vector[names(ownpepsflat)]+min(ecGp$vector[ecGp$vector>0])/2, 
     ownpepsflat+min(ownpepsflat[ownpepsflat>0])/2, 
     log="xy", pch=16, cex=0.5) #, col=rgb(0,0,0,0.15)

# Peptides chooser ------

# only for lGpos
es=ends(Gpos, E(Gpos))
es=cbind(apply(es,2,as.character),LCS=as.character(edge_attr(Gpos)$LCS))
es1=es[,-2]
es2=es[,-1]
es=rbind(es1,es2)
x=aggregate(data=es, es[,1]~es[,2], "list")
y=sapply(x[,2], table)
names(y)=x[,1]
lcs2pep=y
rm(x,y,es,es1,es2)
save(lcs2pep,file="lcs2pep")

plan("multisession", workers=14)
lgrps=t(future_sapply(lGpsmpeps, function(pp){
  x=lcs2pep[pp]
  rowMeans(sapply(x, function(y){
    colSums(fls[names(y),]*(y/sum(y)))
  }))
}))
closeAllConnections()

nletlvg=sapply(vlG, function(p){
  length(table(unlist(strsplit(p,split=""))))/nchar(p)
})

j5=nletlvg %in% ((1:5)/5)
j6=nletlvg %in% ((1:6)/6)
x5=nletlvg[j5]
X5=bs(x5, knots = 3)
Y5=log10(vertex_attr(lGpos)$N[j5])
XY5lm=lm(Y5~X5)
prdXY5=Y5-predict(XY5lm, as.data.frame(X5))
boxplot(prdXY5~nletlvg[j5], notch=T)
x6=nletlvg[j6]
X6=bs(x6, knots = 3)
Y6=log10(vertex_attr(lGpos)$N[j6])
XY6lm=lm(Y6~X6)
prdXY6=Y6-predict(XY6lm, as.data.frame(X6))
boxplot(prdXY6~nletlvg[j6], notch=T)
Nwt=rep(0, length(vlG))
Nwt[j5]=prdXY5
Nwt[j6]=prdXY6

boxplot(Nwt~nletlvg, notch=T)

# recalculate this -----
plan("multisession", workers=14)
lgrpsbest=future_sapply(lGpsmpeps, function(pp){
  j=vlG %in% pp
  x=max(Nwt[j])
  vlG[j&Nwt==x][1]
})
closeAllConnections()

lgrpsbestpp=sapply(lgrpsbest, function(l) names(lcs2pep[[l]])[which.max(lcs2pep[[l]])])
intersect(lgrpsbestpp,ecGpsm_good)
lgrpppercl=sapply(Gposmpeps, function(p) sum(p %in% lgrpsbestpp))
lgrpppercl_g=sapply(Gposmppfr, function(p) sum(p %in% lgrpsbestpp))
jb=!(names(lgrpppercl) %in% names(lgrpppercl_g))
table(lgrpppercl)
table(lgrpppercl_g)
table(lgrpppercl[jb])
hist(log10(lgrpppercl_g/lengths(Gposmppfr)), xlim=c(-4,-1),main="", xlab="")
par(new=T)
hist(log10(lgrpppercl[jb]/lengths(Gposmpeps)[jb]), xlim=c(-4,-1), col=rgb(1,0,0,0.3), main="", xlab="Log10(Proportion of LCS based)")
legend("topleft", fill=c("grey", rgb(1,0,0,0.3)), legend=c("Significant Clusters", "Non-significant Clusters"), bty="n")

Xp=as.numeric(table(sapply(ecGpsm_good, function(x) length(table(strsplit(x, split=""))))))
Xl=as.numeric(table(sapply(lgrpsbestpp, function(x) length(table(strsplit(x, split=""))))))

cols0=lgrps
cols1=(apply(cols0, 2, function(x) (x-min(x))/diff(range(x))))

for (i in 1:5){
  lGpsm=set_vertex_attr(lGpsm, name=colnames(cols0)[i], value=cols1[,i])
}

x=grEmbed(lGpsm,att=9:13, fl="emb")
nms=unlist(vertex_attr(x)$name)
x=delete_vertex_attr(x, "name")
x=set_vertex_attr(x,name="name",value=nms)
write.graph(x, format = "graphml", file="lGpsm_emb_only.graphml")

np=unique(unlist(c(negm,negg)))
pepneg=unique(unlist(sapply(c(negm,negg), function(pp) sample(pp,round(0.0179*length(pp))))))
pepposlg=unique(lgrpsbestpp)
pepposG=unique(ecGpsm_good)
peptotest=unique(c(pepposG,pepposlg,pepneg))
peptotest=data.frame(peptotest,code=(peptotest %in% pepposG)+2*(peptotest %in% pepposlg)+4*(peptotest %in% pepneg))
peptotest[,2]=cut(peptotest[,2],c(0,1.5,3,5), labels=c("G","lG","N"))
write.csv(peptotest, file="peptotest.csv")

plan("multisession", workers=16)
ppb=peptotest[peptotest[,2]!="N",1]
ppttlab=t(future_sapply(ppb,function(p){
  pm=unlist(Gmw_good[grep(p,Gm39peps),1:2])
  if (sum(lengths(pm))==0) pm=c(Control=0,APLS=0)
  names(pm)=paste(c("C","A"),"IgM",sep="_")
  pg=unlist(Ggw_good[grep(p,Gg40peps),1:2])
  if (sum(lengths(pg))==0) pg=c(Control=0,APLS=0)
  names(pg)=paste(c("C","A"),"IgG",sep="_")
  data.frame(t(pm),t(pg))
}))
closeAllConnections()

x=matrix(unlist(ppttlab),nrow=length(ppb))
colnames(x)=colnames(ppttlab)
peptypes=aggregate(ppb, by=as.data.frame(x), "list")
peptypes=peptypes[order(lengths(peptypes$x), decreasing = T),]
xn=lengths(peptypes$x)
peptypes=cbind(peptypes,xn)

#pptpshrt=read.csv("peptypesshorter.csv")

MCu=(ppttlab[,1]>0)
MCd=(ppttlab[,1]<0)
MAu=(ppttlab[,2]>0)
MAd=(ppttlab[,2]<0)
GCu=(ppttlab[,3]>0)
GCd=(ppttlab[,3]<0)
GAu=(ppttlab[,4]>0)
GAd=(ppttlab[,4]<0)
venn=data.frame(MCu,MCd,MAu,MAd,GCu,GCd,GAu,GAd)
venn=euler(venn, shape="ellipse")
plot(venn, quantities=T)
venn=data.frame(MAd,MAu,GAd,GAu)
venn=euler(venn, shape="ellipse")
plot(venn, quantities=T)

write.csv(peptotest, file="peptotest.csv")
seq2fasta(peptotest, fnm="pptt.fasta")



## Distribution of the selected peptides in the giant Graph -----

Gp=Gpos
Gp=set.edge.attribute(Gp, "weight", value=1/edge_attr(Gp)$weight)
pp=peptotest[peptotest$code!="N",1]
dpepch=distances(Gp, V(Gp)[unlist(pp)])
hist(dpepch, breaks=1000)
dppch=distances(Gp, V(Gp)[pp], V(Gp)[pp])
hist(dppch, breaks=1000, xlim=range(dpepch))
tdpep=table(cut(dpepch[lower.tri(dpepch)], c(0,100,800,2000,4000), labels=F))/length(dpepch[lower.tri(dpepch)])
tdpp=table(cut(dppch[lower.tri(dppch)], c(0,100,800,2000,4000), labels=F))/length(dppch[lower.tri(dppch)])
plan("multisession", workers=16)
tpep=t(future_sapply(1:32, function(i){
  print(i)
  x=sample(names(V(Gp)), 558)
  d=distances(Gp, V(Gp)[x], V(Gp)[x])
  table(cut(d[lower.tri(d)], c(0,100,800,2000,4000), labels=F))/length(d[lower.tri(d)])
}))
closeAllConnections()

boxplot(tpep, log="y")
par(new=T)
plot(tdpp, ylim=range(tpep), ty="b", col=2, log="y")

cbind(colRanges(tpep), tdpp)

Gptt=induced.subgraph(Gpos, V(Gpos)[pp])
components(Gptt)
tdgptt=degree(Gptt)
plot(table(tdgptt))
ttdgptt=table(tdgptt)
plot((as.numeric(names(ttdgptt))), log10(ttdgptt))
ppsrc=(pp %in% lgrpsbestpp) +2*(pp %in% ecGpsm_good)
names(ppsrc)=pp
Gptt=set_vertex_attr(Gptt, name="Source", value=ppsrc[names(V(Gptt))])
write.graph(Gptt,format="graphml", file="Gptt.graphml")
modularity(Gptt,ppsrc[names(V(Gptt))])
assortativity(Gptt,ppsrc[names(V(Gptt))])
table(ppsrc, components(Gptt)$membership)
x=rbind(c(184,79),table(ppsrc)-c(184,79))
chisq.posthoc.test(as.data.frame(x))
