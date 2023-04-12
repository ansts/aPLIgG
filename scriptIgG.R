require(parallel)
require(igraph)
require(reshape2)
require(stringdist)
require(matrixStats)
require(chisq.posthoc.test)
require(heatmap3)
require(uwot)
require(corrplot)
require(stringi)
require(qualV)
require(eulerr)
require(ggplot2)
require(ggseqlogo)
require(MASS)
require(multcomp)
require(smallstuff)
require(pbapply)
require(infotheo)
require(rsetse)
require(R.utils)
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

euler(fl, shape="ellipse")
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

load("Ggiant")
load("flnew")

fls=fl[names(V(G)),]

## IgM ----
gm=fls$cM|fls$aM
Gm=induced.subgraph(G,gm)
save(Gm, file="Gm")

Gm38=clusterGraph(Gm, resol=38)
load("Gctrpeps")
Gm38peps=Gpeps
save(Gm38peps, file="Gm38peps")
rm(Gm)

grps=names(vertex_attr(Gm38)$Group[[1]])
cim=as.numeric(substr(grps,2,3))==10
aim=as.numeric(substr(grps,2,3))==1
caim=as.numeric(substr(grps,2,3))==11
d=cim+aim+2*caim

Gm38CA=t(sapply(vertex_attr(Gm38)$Group, function(x){
  x=as.data.frame(x)
  cntr=sum(x[cim,])
  apl=sum(x[aim,])
  all=sum(x[caim,])
  return(c(cntr, apl, all))
}))
colnames(Gm38CA)=c("Control","APLS","All").
rownames(Gm38CA)=names(V(Gm38))

y=Gm38CA
chsqmCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=(1:(nrow(chsqmCA)/2))*2
j1=j2-1
Gm_good=(chsqmCA[j2,3:5]<0.05)*sign(chsqmCA[j1,3:5])
rownames(Gm_good)=chsqmCA$Dimension[j2]

Gm38P=sapply(vertex_attr(Gm38)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})


## IgG ----

gg=fls$cG|fls$aG
Gg=induced.subgraph(G,gg)
save(Gg, file="Gg")

Gg45=clusterGraph(Gg, resol=45)
load("Gctrpeps")
Gg45peps=Gpeps
save(Gg45peps, file="Gg45peps")
rm(Gg)

grps=names(vertex_attr(Gg45)$Group[[1]])
cig=as.numeric(substr(grps,4,5))==10
aig=as.numeric(substr(grps,4,5))==1
caig=as.numeric(substr(grps,4,5))==11
d=cig+aig+2*caig

Gg45CA=t(sapply(vertex_attr(Gg45)$Group, function(x){
  x=as.data.frame(x)
  cntr=sum(x[cig,])
  apl=sum(x[aig,])
  all=sum(x[caig,])
  return(c(cntr, apl, all))
}))
colnames(Gg45CA)=c("Control","APLS","All")
rownames(Gg45CA)=names(V(Gg45))

y=Gg45CA
chsqgCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=(1:(nrow(chsqgCA)/2))*2
j1=j2-1
Gg_good=(chsqgCA[j2,3:5]<0.05)*sign(chsqgCA[j1,3:5])
rownames(Gg_good)=chsqgCA$Dimension[j2]

Gg45P=sapply(vertex_attr(Gg45)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})

load("qq7I")
qq7Im=sapply(Gm38peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Im=log10(qq7Im+1e-5)
qq7Ig=sapply(Gg45peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Ig=log10(qq7Ig+1e-5)

Gm38=set.vertex.attribute(Gm38,name="C",value=Gm38CA[,1]/vertex_attr(Gm38)$size)
Gm38=set.vertex.attribute(Gm38,name="C_sig",value=Gm_good[,1])
Gm38=set.vertex.attribute(Gm38,name="A",value=Gm38CA[,2]/vertex_attr(Gm38)$size)
Gm38=set.vertex.attribute(Gm38,name="A_sig",value==Gm_good[,2])
Gm38=set.vertex.attribute(Gm38,name="all_sig",value=Gm_good[,3])
Gm38=set.vertex.attribute(Gm38,name="Pub",value=Gm38P/vertex_attr(Gm38)$size)
Gm38=set.vertex.attribute(Gm38,name="Id",value=qq7Im)

Gg45=set.vertex.attribute(Gg45,name="C",value=Gg45CA[,1]/vertex_attr(Gg45)$size)
Gg45=set.vertex.attribute(Gg45,name="C_sig",value=Gg_good[,1])
Gg45=set.vertex.attribute(Gg45,name="A",value=Gg45CA[,2]/vertex_attr(Gg45)$size)
Gg45=set.vertex.attribute(Gg45,name="A_sig",value=Gg_good[,2])
Gg45=set.vertex.attribute(Gg45,name="all_sig",value=Gg_good[,3])
Gg45=set.vertex.attribute(Gg45,name="Pub",value=Gg45P/vertex_attr(Gg45)$size)
Gg45=set.vertex.attribute(Gg45,name="Id",value=qq7Ig)
save(Gm38, file="Gm38")
save(Gg45, file="Gg45")

write_graph(Gm38,format = "graphml", file="Gm38.graphml")
write_graph(Gg45,format = "graphml", file="Gg45.graphml")

d=round(corrD(Gm38))
x=grEmbed(Gm38, att=c(5,7,10,11),  k=d)
write_graph(x,format = "graphml", file="Gm38e.graphml")
d=round(corrD(Gg45))
x=grEmbed(Gg45, att=c(5,7,10,11), k=d)
write_graph(x,format = "graphml", file="Gg45e.graphml")
d=round(corrD(G70))
x=grEmbed(G70, att=c(4:10), k=d)
write_graph(x,format = "graphml", file="G70e.graphml")

# Fusion graph -----
load("Ggiant")
allposcl=list(G70peps[rowSums(abs(G_good))>0],
              Gm38peps[rowSums(abs(Gm_good[,2:3]))>0],
              Gg45peps[rowSums(abs(Gg_good[,2:3]))>0])
allpos=unique(unlist(allposcl))
Gpos=induced.subgraph(G,V(G)[allpos])
Gpos=induced.subgraph(Gpos,components(Gpos)$membership==1)
ecGp=eigen_centrality(Gpos)
rm(G)
gc()


## Line graph fusion ----
lGpos=make_line_graph(Gpos)
lGpos=set_vertex_attr(lGpos, name="name", value=edge_attr(Gpos)$LCS)
lGpos=set_vertex_attr(lGpos, name="N", value=1)
lGpos=set_edge_attr(lGpos, name="weight", value=1)
lGpos=contract(lGpos,as.numeric(as.factor(vertex_attr(lGpos)$name)), vertex.attr.comb = list(N ="sum", name="first"))
lGpos=simplify(lGpos, edge.attr.comb = "sum")

dlG=degree(lGpos)
vlG=names(V(lGpos))
eclG=eigen_centrality(lGpos)
lGpos=set_vertex_attr(lGpos, name="size", value=log10(vertex_attr(lGpos)$N))

lGpsm=clusterGraph_l(lGpos, resol=20)
load("Gctrpeps_l")
lGpsmpeps=Gpeps
x=edge_attr(lGpsm)$weight
x=log10(x)
x=(x-min(x)+0.1)/diff(range(x))
lGpsm=set_edge_attr(lGpsm, name="weight",value=x)
lGpsm=set_vertex_attr(lGpsm, name="sizeR",value=sapply(vertex_attr(lGpsm)$N, function(x) sqrt(sum(x))))
write.graph(lGpsm, format = "graphml", file="lGpsm.graphml")

lGptrim=induced_subgraph(lGpos, V(lGpos)[vertex_attr(lGpos)$N>3])
write.graph(lGptrim, format = "graphml", file="lGptrim.graphml")

lGptrimsm=clusterGraph_l(lGptrim, resol=20)
load("Gctrpeps_l")
lGptrimsmpeps=Gpeps
write.graph(lGptrimsm,format = "graphml", file="lGptrimsm.graphml")

logoslGptrim=lapply(lGptrimsmpeps,function(s) {
  x=s[nchar(s)==6]
  s=s[nchar(s)==5]
  x=c(sapply(x,function(i) c(substr(i,1,5),substr(i,2,6))))
  s=c(s,x)
  ggseqlogo(s)
})
save(logoslGptrim,file="logoslGptrim") 

pdf(file="logoslGptrim.pdf", width=5,height = 5)
for (i in seq_along(logoslGptrim)) plot(logoslGptrim[[i]])
dev.off()

complx=future_sapply(vlG,function(p) mean(table(unlist(strsplit(p,split="")))))
complx_deg=cbind(complx,dlG)
complx_N=cbind(complx,N=vertex_attr(lGpos)$N)
boxplot(data=complx_deg,log10(dlG)~complx,notch=T)
boxplot(data=complx_N,log10(N)~complx,notch=T)

complx=future_sapply(unlist(lGptrimsmpeps),function(p) mean(table(unlist(strsplit(p,split="")))))
complx_deg=cbind(complx,dlGt=degree(lGptrim))
complx_N=cbind(complx,N=vertex_attr(lGptrim)$N)
boxplot(data=complx_deg,log10(dlGt)~complx,notch=T)
boxplot(data=complx_N,log10(N)~complx,notch=T)

# There is some negative correlation between the mean degree or the 
# mean N in the lGpos and the entropy of the lcs-s for the length 5 
# but it is lost in lGptrim.

## Eigencentrality based selection ----

# nx=vertex_attr(lGpos)$N  
# names(nx)=names(V(lGpos))
# y=eclG$vector
# y[is.na(y)|y==0]=min(y[!is.na(y)&y>0])/2

ecchng=future_lapply(lGpsmpeps, function(l) {
  #eclst=closeness(induced.subgraph(lGpos,l))
  # mn=min(eclst[!is.na(eclst)&eclst>0])
  # eclst[is.na(eclst)|eclst==0]=mn/2
  #x0=names(which.max(eclst*y[names(eclst)]))
  if(length(l)>1) {
    eclst=eigen_centrality(induced.subgraph(lGpos,l))$vector
    x0=names(which.max(eclst))
  }
  else {
    x0=l
  }
  x=names(ego(induced.subgraph(lGpos,l),1,x0)[[1]])
  xn=nx[x]
  xn=sort(xn, decreasing = T)
  i=seq_along(xn)
  #thr=length(xn)
  thr=sapply(1:length(xn),function(i) sum(xn[1:i])/sum(xn))
  thr=length(thr[thr<=0.33])+1
  #clr=round(log10(nx[names(eclst$vector)]+0.5))+1
  #plot(log10(eclst$vector), log10(eclG$vector[names(eclst$vector)]), cex=0)
  #text(log10(eclst$vector), log10(eclG$vector[names(eclst$vector)]), cex=0.2*clr,labels = names(eclst$vector), col=clr)
  return(list(x0,xn[1:thr]))
})
lcs_g=sapply(ecchng, function(l) names(l[[2]]))

x=unlist(sapply(lcs_g, function(y) y[1]))
dx=distances(lGpos, x,x, weights=NA)
range(dx[dx>0])
diag(dx)=100
table(rowMins(dx))
ecmaxntwk=components(induced_subgraph(lGpos,x))$membership[x]


logos_lcs_g=lapply(lcs_g,function(s) {
  ggseqlogo(s)
})
save(logos_lcs_g,file="logos_lcs_g") 

pdf(file="logos_lcs_g.pdf", width=5,height = 5)
for (i in seq_along(logos_lcs_g)) plot(logos_lcs_g[[i]])
dev.off()


# Peptides chooser ------

load("Ggiant")
es=ends(G, E(G))
lcs=edge_attr(G)$LCS
G=NULL
rm(G)
gc()
es=cbind(es,lcs)
es1=es[,-2]
es2=es[,-1]
es=rbind(es1,es2)
rm(es1,es2)
es=aggregate(es[,2], by=list(es[,1]), "list")
x=sapply(es[,2], table)
names(x)=es[,1]
pep2lcs=x

# only for lGpos
es=es[es[,2] %in% names(V(lGpos)) & es[,1] %in% names(V(Gpos)),]
x=aggregate(es[,1], by=list(es[,2]), "list")
y=sapply(x[,2], table)
names(y)=x[,1]
lcs2pep=y
rm(x,y,es)

plan("cluster", workers=makeCluster(16))
pepbe=future_sapply(lcs_g, function(l){
  z=melt(lcs2pep[l])
  z=aggregate(z$value, by=list(z$Var1), "sum")
  x=z$x
  names(x)=z$Group.1
  z=x[order(x, decreasing = T)] # peptides having lcs from the ith cluster 
  return(z)   # ordered by the # of edges with their NN named with these lcs-s 
})
closeAllConnections()

for (j in seq_along(pepbe)){
  p=names(pepbe[[j]])
  seq2fasta(S=p,fnm=paste("pepsel/pepsel",j,collapse = "_"))
}

pepbingood=future_sapply(1:80, function(j){
    x=t(sapply(seq_along(pepbe), function(i){
      l=pepbe[[i]]
      print(i)
      if (length(l)<j) return(rep(0,8))
      p=names(l)[j] 
      pall=c(G_good[grep(p,G70peps),][,c(8,4,2,1)])
      pm=unlist(c(Gm_good[grep(p,Gm38peps),1:2],0,0))
      pg=unlist(c(0,0,Gg_good[grep(p,Gg45peps),1:2]))
      m=rbind(pall,pm,pg)
      mup=colSums(m==1)
      mdown=colSums(m==-1)
      return(c(mup,-mdown))
    }))
    return(colSums(x))
})

x=sapply(1:ncol(pepbingood),function(j) pepbingood[,j]/sum(lengths(pepbe)>(j-1)))

boxplot(t(x), notch=T, ylim=c(-1,1))

xm=rowMeans(x)
xm[5:8]=-xm[5:8]
xsd=rowSds(x)
plot(xm[1:4], xm[5:8], cex=0, xlim=c(-0.01,0.8),ylim=c(-0.01,0.8),xlab="Proportion overexressed", ylab="Proportion underexpressed")
par(new=T)
for (i in 1:4) {
  plot(c(xm[i]-xsd[i],xm[i]+xsd[i]), c(xm[i+4],xm[i+4]), ty="l",xlim=c(-0.01,0.8),ylim=c(-0.01,0.8), xlab="",ylab="")
  par(new=T)
}
for (i in 1:4) {
  plot(c(xm[i],xm[i]), c(xm[i+4]-xsd[i],xm[i+4]+xsd[i]), ty="l",xlim=c(-0.01,0.8),ylim=c(-0.01,0.8), xlab="",ylab="")
  par(new=T)
}
text(xm[1:4]+0.05, xm[5:8]+0.05,labels=c("IgM_cntr","IgM_APS","IgG_cntr","IgG_APS"))

peptotest=unlist(lapply(pepbe, function(p) {
  n=length(p) %/% 100
  if(n==0) n=1 
  names(p)[1:n]
  #names(p[1])
}))
ppttlab=t(future_sapply(peptotest,function(p){
  pall=unlist(G_good[grep(p,G70peps),])
  pm=unlist(Gm_good[grep(p,Gm38peps),1:2])
  if (sum(lengths(pm))==0) pm=c(Control=0,APLS=0)
  names(pm)=paste(c("C","A"),"IgM",sep="_")
  pg=unlist(Gg_good[grep(p,Gg45peps),1:2])
  if (sum(lengths(pg))==0) pg=c(Control=0,APLS=0)
  names(pg)=paste(c("C","A"),"IgG",sep="_")
  data.frame(t(pall),t(pm),t(pg))
}))
x=hclust(dist(ppttlab), method = "ward.D2")
plot(x)
x$order
x=matrix(unlist(ppttlab),nrow=length(peptotest))
colnames(x)=colnames(ppttlab)
peptypes=aggregate(names(peptotest), by=as.data.frame(x), "list")
peptypes=peptypes[order(lengths(peptypes$x), decreasing = T),]
xn=lengths(peptypes$x)
peptypes=cbind(peptypes,xn)

#pptpshrt=read.csv("peptypesshorter.csv")

MCu=((ppttlab[,8]>0)+(ppttlab[,12]>0))>0
MCd=((ppttlab[,8]<0)+(ppttlab[,12]<0))>0
MAu=((ppttlab[,4]>0)+(ppttlab[,13]>0))>0
MAd=((ppttlab[,4]<0)+(ppttlab[,13]<0))>0
GCu=((ppttlab[,2]>0)+(ppttlab[,14]>0))>0
GCd=((ppttlab[,2]<0)+(ppttlab[,14]<0))>0
GAu=((ppttlab[,1]>0)+(ppttlab[,15]>0))>0
GAd=((ppttlab[,1]<0)+(ppttlab[,15]<0))>0
venn=data.frame(MCu,MCd,MAu,MAd,GCu,GCd,GAu,GAd)
venn=euler(venn, shape="ellipse")
plot(venn, quantities=T)
venn=data.frame(MAd,MAu,GAd,GAu)
venn=euler(venn, shape="ellipse")
plot(venn, quantities=T)

write.csv(peptotest, file="peptotest.csv")
seq2fasta(peptotest, fnm="pptt.fasta")

j=dlG %in% 1:51
i=cut(log10(vertex_attr(lGpos)$N[j]), 25, labels = F)
boxplot(log10(eclG$vector[j])~i, notch=T)
aggregate(vertex_attr(lGpos)$N[j], by=list(i), "range")

j=dlG>100
i=cut(log10(vertex_attr(lGpos)$N[j]), 50, labels = F)
boxplot(log10(eclG$vector[j])~i, notch=T)
aggregate(vertex_attr(lGpos)$N[j], by=list(i), "range")

# Relation between ec in Gpos and in lGpos
topecGp=ecGp$vector[ecGp$vector>0.2]
n=names(topecGp)
l=pep2lcs[n]
plot(lengths(l), ty="l")

topec_lcl=future_lapply(n,function(p){
  sort(table(unlist(apply(sapply(l[p],names),1,grep,lGpsmpeps))),decreasing=T)
})

# crosscheck Gpos and lGpos clustering

Gposm=clusterGraph(Gpos, resol=25)
load("Gctrpeps")
Gposmpeps=Gpeps
x=Gposmpeps
x=melt(x)
y=x[,2]
names(y)=x[,1]
modularity(Gpos,as.factor(y[names(V(Gpos))]))

ppperlcl=lapply(lGpsmpeps, function(l) table(unlist(lapply(lcs2pep[l],names))))
tppperl=table(unlist(lapply(ppperlcl, names)))
length(tppperl)
table(tppperl)
lclrepp=melt(ppperlcl)
x=aggregate(lclrepp[,2:3],by=list(lclrepp$Var1),"c")
x[[2]]=sapply(seq_along(x[[3]]),function(i) x[[3]][[i]][which.max(x[[2]][[i]])])            
x=x[,-3]            
y=x[,2]            
names(y)=x[,1]          
x=y
modularity(Gpos,as.factor(x[names(V(Gpos))]))

tpppcl=melt(x)
tpppclpeps=aggregate(rownames(tpppcl), by=list(tpppcl$value), "list")
x=tpppclpeps$x
names(x)=tpppclpeps$Group.1
tpppclpeps=x
mxlcl=sapply(Gposmpeps, function(p1){
      sapply(tpppclpeps, function(p2){ 
        length(intersect(p1,p2))
      })
})

require(RcppHungarian)

sln=HungarianSolver(log10(3097-mxlcl))
o=sln$pairs
ow=mxlcl[o]/rowMins(cbind(lengths(tpppclpeps)[o[,1]],lengths(Gposmpeps)[o[,2]]))
image(log10(mxlcl+1)[sln$pairs[order(ow,decreasing=T) ,1], sln$pairs[order(ow,decreasing=T),2]])

# ec of Gposm

ppeccGpos=sapply(Gposmpeps, function(l) {
  eclst=eigen_centrality(induced.subgraph(Gpos,l))
  names(which.max(eclst$vector*ecGp$vector[names(eclst$vector)]))
})

ppgood=intersect(ppeccGpos,peptotest)           

ppecc_=setdiff(ppeccGpos,peptotest)
ppptt_=setdiff(peptotest,ppeccGpos)

dpppp=distances(Gpos,v=ppecc_,to=ppptt_)
ppplus=ppecc_[rowMins(dpppp)==2]

ppNL=c(ppgood,ppptt_,ppplus)
dx=distances(Gpos,ppNL,ppNL)
x=dx
diag(x)=10
rm=rowMins(x)

# Connectivity in lGpos and the lGpsm cluster boundaries
y=eclG$vector
y[is.na(y)|y==0]=1e-20

ecmxmncon=t(future_sapply(lGpsmpeps, function(l) {
  L=length(l)
  eclst=eigen_centrality(induced.subgraph(lGpos,l))$vector
  eclst[is.na(eclst)|eclst==0]=1e-20
  ec=eclst*y[names(eclst)]
  ec2=max(ec)
  ec1=min(ec)
  ecmax=names(which.max(ec))
  ecmin=names(which.min(ec))
  nn2=length(intersect(l,names(ego(induced.subgraph(lGpos,l),1,ecmax)[[1]])))
  nn1=length(intersect(l,names(ego(induced.subgraph(lGpos,l),1,ecmin)[[1]])))
  return(c(L=L, ecmin=ec1, ecmax=ec2, NNmin=nn1, NNmax=nn2))
}))
plot(as.data.frame(log10(ecmxmncon)), cex=0.7, pch=16, col=rgb(0,0,0,0.3))
for (i in 1:ncol(ecmxmncon)) lGpsm=set_vertex_attr(lGpsm,name=colnames(ecmxmncon)[i], value=ecmxmncon[,i])

lGpsm=set_vertex_attr(lGpsm, name="ClofCl", value=ecmaxntwk)
write.graph(lGpsm, format = "graphml", file="lGpsm.graphml")
