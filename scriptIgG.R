require(parallel)
require(pbapply)
require(igraph)
require(reshape2)
require(stringdist)
require(matrixStats)
require(chisq.posthoc.test)
require(heatmap3)
require(uwot)
require(corrplot)
require(stringi)

require(ggseqlogo)
require(infotheo)
require(rsetse)
require(R.utils)
require(eulerr)

# Data Acquisition and Cleaning ----------------------------------------------

cpl3=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))
cpl1=colorRampPalette(c("#0000FF80","#00FF0080","#FFFF0080","#FF000080"))
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

J=2
cl=makeCluster(J)
#clusterExport(cl, "G")
clusterEvalQ(cl, require(igraph))
GlouvTable=pblapply(1:(6*J),function(i){
  load("G")
  n=(i-1)*5+1
  cluster_louvain(G, resolution = n)
},cl=cl)
stopCluster(cl)
save(GlouvTable,file="GlouvTable") #exists as GlouvTable (1:12) and GlouvTable1 (13:24)
GlouvT=sapply(1:24, function(i){
  ij=((i-1) %/% 2)+1
  ji=(i %% 2)+1
  GlouvTable[[ij]]$memberships[ji,]
})
cnms=c((0:5)*5+1,(0:5)*5+31 )
cnms=paste(rep(cnms, each=2),rep(1:2,12), sep="_")
colnames(GlouvT)=cnms
save(GlouvT, file="GlouvT_final")

crosT=lapply(1:11,function(i){
  print(i)
  table(GlouvT[,(i*2)],GlouvT[,(i*2)+2])
})

mutI=sapply(1:11,function(i){
  print(i)
  mutinformation(GlouvT[,(i*2)],GlouvT[,(i*2)+2])
})

load("GlouvT_final")
Gl24=GlouvT[,24]
n=max(Gl24)
load("G")
names(Gl24)=names(V(G))
proct=proc.time()
G24=contract(G,mapping=Gl24,"ignore")
G24=simplify(G24, edge.attr.comb = "sum")
print(proc.time()-proct)
save(G24, file="G24")

Groupsall=table(V(G)$Group)  
gnm=names(Groupsall)
x0=rep(0,length(gnm))
names(x0)=gnm
propG24=aggregate(Gl24, by=list(Gl24), function(x){
  i=Gl24==unique(x)
  x1=length(x)
  g=induced_subgraph(G,i)
  x2=length(E(g))
  x3=table(V(G)$Group[i])
  x0[names(x3)]=x3
  c(Size=x1,Edges=x2,x0)
})
propG24=propG24[,-1]
save(propG24,file="propG24")
chisq.posthoc.test(propG24[,-c(1,2)])

x=propG24[,-c(1,2)]
xnm=substr(colnames(x),2,5)
xrn=seq_along(x[,1])
y0=rep(0,length(xnm)/2)
names(y0)=xnm[1:(length(xnm)/2)]
x=t(apply(x,1,function(l){
  y=aggregate(l, by=list(xnm), "sum")
  y0[y[,1]]=y[,2]
  return(y0)
}))
rownames(x)=xrn
x=as.table(x)
x=chisq.posthoc.test(as.table(x), simulate.p.value=T)
px=x[(1:(nrow(x)/2))*2,3:17]
rownames(px)=seq_along(rownames(px))
ex=x[(1:(nrow(x)/2))*2-1,3:17]
rownames(ex)=seq_along(rownames(ex))
pxgood=apply(px, 1,function(l) any(l<0.05))
sigcl=cbind(Eff=ex[pxgood,],P=px[pxgood,])
sigclnm=rownames(sigcl)

x=propG24[,-c(1,2)]
x=t(apply(x,1,function(l){
  y1=sum(l[1:15])
  y2=sum(l[16:30])
  return(c(y1,y2))
}))
x=as.table(x)
pubchsq=chisq.posthoc.test(x, simulate.p.value=T)
ppub=pubchsq[(1:(nrow(pubchsq)/2))*2,4]
names(ppub)=seq_along(ppub)
epub=pubchsq[(1:(nrow(pubchsq)/2))*2-1,4]
names(epub)=seq_along(epub)
pubsig=epub[which(ppub<0.05)]
pubsig=intersect(pubsig, as.numeric(sigclnm))
hilit=which(px[pxgood,]<0.05, arr.ind=T)
hilit=data.frame(hilit,rep(1, nrow(hilit)))
pdf(file="hmsigcl.pdf", width=20,height = 20)
hmsigcl=heatmap3(sigcl[,1:15], method = "ward.D2", col=cpl1(1000), cexRow = 0.5, highlightCell = hilit)
dev.off()

NBcl=c(350,248,246,1809,351,91,1766,1027,915,481,387,372,318,477,486,1886,400,
       1808,497,1794,485,219,483,521,355,480,1799,1798,2212,1028,1797,1813,
       1566,348,523,1796,1561,479,460,1800,365,220,1810,1707,364,361,50,599,
       1801,819,3,10,1404,450,636,1339,1349,532,1579,2064,187,1344,491,817,2293,
       310,527,1555,1852,296,1460,1729,1728,881,530,1290,2346,590,188,508,426,
       329,288,290,591,1834,723,284,2087,178,1348,223,271,1730,1440,1531,582,
       392,1842,515,1750,1772,2221,1308,467,101,1340,1546,249,1582,1713,
       1841,275,1450,339,1732,1469,1319,1457,593,577,445,571,572)
ip=apply(px[pxgood,],1,function(l) colnames(px)[l<0.05])
sigcl1=lapply(seq_along(ip), function(i) {
  x=ex[pxgood,][i, ip[[i]]]
  n=ip[[i]]
  names(x)=n
  return(x)
}) 
t0=c(0,0)
names(t0)=c(-1,1)
syngcl=aggregate(unlist(sigcl1), by=list(names(unlist(sigcl1))), function(x) {
  t1=table(sign(x))
  t0[names(t1)]=t1
  return(t0)
})

attrG=array(0,dim=c(nrow(propG24),30))
colnames(attrG)=colnames(ex)
for (rn in rownames(sigcl)){
  print(rn)
  j=sigcl[rn,31:60]<0.05
  attrG[as.numeric(rn),]=as.numeric(sign(sigcl[,1:30][rn,]*j))
}
save(attrG,file="attrG")

for (j in 1:ncol(attrG)){
  G24=set.vertex.attribute(G24, name=colnames(attrG)[j], value=attrG[,j])
}

ij=cut(seq_along(E(G24)),140,labels=F)

cl=makeCluster(14)
clusterExport(cl, c("G24","propG24","ij"))
clusterEvalQ(cl, require(igraph))

x=pbsapply(1:140, function(j){ # using a log(1/estimate) of the modularity as a dissimilarity measure
  sapply(seq_along(E(G24))[ij==j], function(i){
    es=ends(G24,i)            # indices of incident vertices
    ebtw=E(G24)$weight[i]     # number of between edges (weight of edge after contraction)
    s=propG24[es[1,],1]       # sizes of the two respective clusters in G
    allein=sum(c(propG24[es[1,],2])) # No of edges inside both clusters
    pein=allein/(s[1]*(s[1]-1)/2+s[2]*(s[2]-1)/2) # graph density inside clusters
    pebtw=ebtw/prod(s)        # edge density between clusters 
    log10(pein/pebtw)        #  log of the ratio of within/between edge density
  })
},cl=cl)
  
stopCluster(cl)

w=x[[1]]
for (i in 2:140){
  w=c(w,x[[i]])
}
x=w[order(w)]
x=x[length(x)]-x[length(x)-1]
w=w-min(w)+x/2

G24=set.edge.attribute(G24, name="weight", value=w)
save(G24,file="G24")

propG24[components(G24)$membership>1,]
loners=names(V(G))[components(G)$membership>1]
fl[loners,]
# the disconnected sequences are mimotopes of APL associated antibodies 
# 5/8 from the IgG sample and 2/8 from the IgM sample. 1 is from the control IgM
# which is found in the public repertoire too as is one of the IgG's.
# Note the KCCVFQV sequence similar to a motif of 2G12 mimotopes.

x=mst(G24)
thr=range(E(x)$weight)
cl=makeCluster(14)
clusterExport(cl,c("G24"))
clusterEvalQ(cl, require(igraph))
edgewscan=t(pbsapply(seq(0,2.5,by=0.001), function(th){
  X=delete.edges(G24,E(G24)[E(G24)$weight>th])
  X=components(X)
  x=length(E(delete_edges(G24, E(G24)[E(G24)$weight>th])))
  mxCs=max(X$csize)
  c(Thr=th,ENo=x, Ncomp=X$no, Giant=mxCs)
}, cl=cl))
stopCluster(cl)
plot(as.data.frame(edgewscan), pch=16,cex=0.1)
plot(scale(edgewscan[,2])*scale(edgewscan[,3]), ty="l")
j=(scale(edgewscan[,2])*scale(edgewscan[,3]))
plot(j)
j=which.max(j)
edgewscan[j,]
# select the threshold optimizing between edge No. and components No
G24sm=delete.edges(G24,E(G24)[E(G24)$weight>1.214])
G24sm=set.vertex.attribute(G24sm,name = "Size", value=propG24[,1])
louvG24sm=cluster_louvain(G24sm)
lvG24sm=louvG24sm$memberships[1,]
G24sm=set.vertex.attribute(G24sm,name = "Louv", value=lvG24sm)
write.graph(G24, format="graphml", file = "G24.graphml")
write.graph(G24sm, format="graphml", file = "G24sm.graphml")
save(G24sm, file="G24sm")
save(edgewscan,file="edgewscan")


gs# Smooth Property Mapping -------------------------------------------------

load("propG24")
x=as.data.frame(propG24)
xG=(x$`00001`+x$`00010`+x$`00011`+x$`00101`/2+x$`00110`/2+x$`00111`/2+
           x$`01001`/2+x$`01010`/2+x$`01011`/2+x$`01101`/2+x$`01110`/2+
           x$`01111`/2+x$`10001`+x$`10010`+x$`10011`+x$`10101`/2+x$`10110`/2+
           x$`10111`/2+x$`11001`/2+x$`11010`/2+x$`11011`/2+x$`11101`/2+
           x$`11110`/2+x$`11111`/2)
xM=(x$`00100`+x$`01000`+x$`01100`+x$`00101`/2+x$`00110`/2+x$`00111`/2+
      x$`01001`/2+x$`01010`/2+x$`01011`/2+x$`01101`/2+x$`01110`/2+
      x$`01111`/2+x$`10100`+x$`11000`+x$`11100`+x$`10101`/2+x$`10110`/2+
      x$`10111`/2+x$`11001`/2+x$`11010`/2+x$`11011`/2+x$`11101`/2+
      x$`11110`/2+x$`11111`/2)
xD=(x$`00001`+x$`00100`+x$`00011`/2+x$`00101`+x$`00110`/2+2*x$`00111`/3+
      x$`01100`/2+ x$`01001`/2+x$`01011`/3+2*x$`01101`/3+x$`01110`/3+
      x$`01111`/2+x$`10001`+x$`10011`/2+x$`10100`+x$`10101`+x$`11100`/2+x$`10110`/2+
      2*x$`10111`/3+x$`11001`/2+x$`11011`/3+2*x$`11101`/3+
      x$`11110`/3+x$`11111`/2)
xC=(x$`00010`+x$`01000`+x$`00011`/2+x$`01010`+x$`00110`/2+x$`00111`/3+
      x$`01100`/2+ x$`01001`/2+2*x$`01011`/3+x$`01101`/3+2*x$`01110`/3+
      x$`01111`/2+x$`10010`+x$`10011`/2+x$`11000`+x$`10101`+x$`11100`/2+x$`10110`/2+
      x$`10111`/3+x$`11001`/2+2*x$`11011`/3+x$`11101`/3+
      2*x$`11110`/3+x$`11111`/2)
xP=(x$`10001`+x$`10010`+x$`10011`+x$`11000`+x$`10100`++x$`10101`+
    x$`11100`+x$`10110`+x$`10111`+x$`11001`+x$`11010`++x$`11011`+x$`11101`+
      x$`11110`+x$`11111`)

propG24GMCD=cbind(xG,xM,xC,xD,xP)
z=sweep(propG24GMCD, 1, rowSums(propG24GMCD), "/")
z=t(apply(z,2,function(zx) {
  attributes(zx)=NULL
  return(zx)
}))

xGC=(x$`00010`+  x$`00011`+  x$`00110`/2+x$`00111`/2+x$`01010`/2+x$`01011`/2+
     x$`01110`/2+x$`01111`/2+x$`10010`  +x$`10011`  +x$`10110`/2+
     x$`10111`/2+x$`11010`/2+x$`11011`/2+x$`11110`/2+x$`11111`/2)
xGA=(x$`00001`+  x$`00011`+  x$`00101`/2+x$`00111`/2+x$`01001`/2+x$`01011`/2+
     x$`01101`/2+x$`01111`/2+x$`10001`+  x$`10011`+  x$`10101`/2+
     x$`10111`/2+x$`11001`/2+x$`11011`/2+x$`11101`/2+x$`11111`/2)
xMC=(x$`01000`+  x$`01100`+  x$`01001`/2+x$`01010`/2+x$`01011`/2+x$`01101`/2+
     x$`01110`/2+x$`01111`/2+x$`11000`  +x$`11100`+  x$`11001`/2+
     x$`11010`/2+x$`11011`/2+x$`11101`/2+x$`11110`/2+x$`11111`/2)
xMA=(x$`00100`+  x$`01100`+  x$`00101`/2+x$`00110`/2+x$`00111`/2+x$`01101`/2+
     x$`01110`/2+x$`01111`/2+x$`10100`  +x$`11100`+  x$`10101`/2+
     x$`10110`/2+x$`10111`/2+x$`11101`/2+x$`11110`/2+x$`11111`/2)
propG24GcMd=cbind(xGC,xGA,xMC,xMA,xP)
z=sweep(propG24GcMd, 1, x$Size, "/")

colnames(z)=colnames(propG24GcMd)
uprGcMd=umap(z,n_neighbors = 50, verbose = T )

sz=2*(propG24[,1]-min(propG24[,1]))/diff(range(propG24[,1]))+0.5
plot(uprGcMd, pch=16, cex=sz, col=rgb(recol(z[,1]),recol(z[,2]),recol(z[,5]),0.75), main="IgG Contr vs APL / Public")
plot(uprGcMd, pch=16, cex=sz, col=rgb(recol(z[,4]),recol(z[,3]),recol(z[,5]),0.75), main="IgM Contr vs APL / Public")
plot(uprGcMd, pch=16, cex=sz, col=rgb(recol(z[,4]),recol(z[,3]),recol(z[,1]),recol(z[,2])), main="IgM vs IgG / Contr vs APL")
plot(uprGcMd, pch=16, cex=sz, col=rgb(0,0,0,recol(xP)))

minz=min(z[z>0])/2
z1=apply(z,2,function(x){
  x=log10(x+minz)
  (x-min(x))/(diff(range(x))+0.02)
})


G24sm=set.vertex.attribute(G24sm,name = "IgM", value=z1[,2])
G24sm=set.vertex.attribute(G24sm,name = "Cntr", value=z1[,3])
G24sm=set.vertex.attribute(G24sm,name = "APL", value=z1[,4])
G24sm=set.vertex.attribute(G24sm,name = "Pub", value=z1[,5])
write.graph(G24sm, format="graphml", file = "G24sm.graphml")

# Logos of Gl24 -----------------------------------------------------------

seqs24=aggregate(names(V(G)), by=list(Gl24), "list")
x=seqs24$Group.1
seqs24=seqs24$x
names(seqs24)=x
logosG24=lapply(seqs24,ggseqlogo)
save(logosG24,file="LOGOSG24") 

xlv=V(G24sm)$Louv
x=sapply(seq_along(xlv), function(i) {
  cbind(Seq=unlist(seqs24[[i]]), Lvcl=xlv[i])
})
x0=x[[1]]
for (i in 2:length(x)) x0=rbind(x0,x[[i]]) 
seqs24lv=aggregate(x0[,1], by=list(x0[,2]), list)
x=seqs24lv$Group.1
seqs24lv=seqs24lv$x
names(seqs24lv)=x
logosG24lv=lapply(seqs24lv,ggseqlogo)
save(logosG24lv,file="LOGOSG24lv")  

print(logosG24lv[lengths(seqs24lv)>10000])

# Map Idiotypes -----------------------------------------------------------

load("IgJT")
vG=names(V(G))

# As qgrams - only complete 7-mers
qq7I=t(qgrams(vG,IgJtrim,q=7))
qq7I=qq7I[vG,2]
save(qq7I,file="qq7I")

x=aggregate(qq7I, by=list(Gl24), "mean")
q7_24=x$x
G24sm=set.vertex.attribute(G24sm,name = "Mean_q7I", value=q7_24)
write.graph(G24sm, format="graphml", file = "G24sm.graphml")


# Using the same LCS conditions as for the rest of the graph
L=seq_along(vG)
q7=colnames(qgrams(IgJtrim,q=7))
ij=cut(L,1400, labels=F)
cl=makeCluster(14)
clusterExport(cl, c("q7","vG","ij"))
clusterEvalQ(cl, require(stringdist))
LcsI5=pbsapply(1:1400, function(i){
  j=ij==i
  sapply(vG[j], function(p){
    x=stringdist(p,q7,method="lcs",nthread = 1)
    sum(x<5)
  })
}, cl=cl)
stopCluster(cl)
save(LcsI5,file="LcsI5")

x=LcsI5[[1]]
for (i in 2:length(LcsI5)) x=c(x,LcsI5[[i]])
LcsI5=x
LcsI5_24=aggregate(LcsI5, by=list(Gl24), function(l) mean(log10(l+0.5)))
LcsI5_24=LcsI5_24$x
hist(LcsI5_24)

save(LcsI5_24,file="LcsI5_24")

G24sm=set.vertex.attribute(G24sm,name = "LcsI5_24", value=LcsI5_24)
write.graph(G24sm, format="graphml", file = "G24sm.graphml")


# All features embedding ---------------------------------------------------

x=mst(G24)
thr=max(E(x)$weight)
G24smcon=delete.edges(G24,E(G24)[E(G24)$weight>thr])


for (n in vertex_attr_names(G24smcon)) G24smcon=delete_vertex_attr(G24smcon, n)
G24smcon=set.vertex.attribute(G24smcon,name = "IgG-C", value=z1[,1])
G24smcon=set.vertex.attribute(G24smcon,name = "IgG-A", value=z1[,2])
G24smcon=set.vertex.attribute(G24smcon,name = "IgM-C", value=z1[,3])
G24smcon=set.vertex.attribute(G24smcon,name = "IgM-A", value=z1[,4])
#G24smcon=set.vertex.attribute(G24smcon,name = "IgMvG", value=log10((z1[,3]+z1[,4]+0.2)/(z1[,1]+z1[,2]+0.2)))
#G24smcon=set.vertex.attribute(G24smcon,name = "CvA", value=log10((z1[,1]+z1[,3]+0.2)/(z1[,2]+z1[,4]+0.2)))
G24smcon=set.vertex.attribute(G24smcon,name = "Pub", value=z1[,5])
G24smcon=set.vertex.attribute(G24smcon,name = "Size", value=propG24[,1])
G24smcon=set.vertex.attribute(G24smcon,name = "Louv", value=lvG24sm)
G24smcon=set.vertex.attribute(G24smcon,name = "Mean_q7I", value=recol(q7_24))

G24smcon=remove_small_components(G24smcon)
AjM24sm=distances(G24smcon, to=V(G24smcon))
selfembG24=cmdscale(AjM24sm,k=58)
x=as.matrix(as.data.frame(vertex_attr(G24smcon)))
x[,6]=recol(x[,6])
x=x[,-7]
#corx=cor(selfembG24,x)

#uattr=umap(x, n_neighbors = 50, verbose = T)
#uself=umap(selfembG24,n_neighbors = 50, verbose = T)

mx=cbind(x,selfembG24/4)
#mx=scale(mx)
uall=umap(mx,n_neighbors = 50, verbose = T)

G24smcon=set.vertex.attribute(G24smcon, name="x",value=uall[,1])
G24smcon=set.vertex.attribute(G24smcon, name="y",value=uall[,2])

plot(G24smcon, vertex.label=NULL, vertex.size=0.2)
save(G24smcon, file="G24smcon")

write.graph(G24smcon, format="graphml", file="G24smcon.graphml")

# Separate IgG and IgM Graphs ---------------------------------------------

load("G")
grps=unique(vertex_attr(G)$Group)
gig=as.numeric(substr(grps,4,5))>0
gim=as.numeric(substr(grps,2,3))>0
ig=vertex_attr(G)$Group %in% grps[gig]
im=vertex_attr(G)$Group %in% grps[gim]

Gg=induced.subgraph(G,ig)
Gm=induced.subgraph(G,im)
save(Gg, file="Ggw")
save(Gm, file="Gmw")

load("Ggw")
Gg=set_edge_attr(Gg,"weight",value=1)
Gg30sm=clusterGraph(Gg, resol=30)
rm(Gg)
load("Gmw")
Gm=set_edge_attr(Gm,"weight",value=1)
Gm30sm=clusterGraph(Gm)
rm(Gm)
#Gg1s m=clusterGraph(Gg, resol = 1)

### Background frequencies ------------

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
theorNN=(choose(7,5)*19^2+6*19)*length(matpep)/(20^7)

load("Gmw")
bkgNNm=bkgNN[names(V(Gm))]
rm(Gm)
load("Ggw")
bkgNNg=bkgNN[names(V(Gg))]
rm(Gg)

## Cross graph similarities -----
cl=makeCluster(14)
clusterExport(cl, "Gm30peps")
Ggxm=pbsapply(Gg30peps,function(p1){
      sapply(Gm30peps, function(p2){
        length(intersect(p1,p2))
      })
}, cl=cl)
stopCluster(cl)

Gginm=apply(Ggxm,1,function(l){
  max(l)/(sum(l)+0.1)
})
Gm30sm=set.vertex.attribute(Gm30sm,name="IgGovrlp",value=Gginm)
Gming=apply(Ggxm,2,function(l){
  max(l)/(sum(l)+0.1)
})
Gg30sm=set.vertex.attribute(Gg30sm,name="IgMovrlp",value=Gming)

Gbging=sapply(Gg30peps, function(p){
  median(bkgNNg[p])
})
Gbginm=sapply(Gm30peps, function(p){
  median(bkgNNm[p])
})
Gm30sm=set.vertex.attribute(Gm30sm,name="Bkgovrlp",value=Gbginm)
Gg30sm=set.vertex.attribute(Gg30sm,name="Bkgovrlp",value=Gbging)
write_graph(Gm30sm,format = "graphml", file="Gm30sm.graphml")
write_graph(Gg30sm,format = "graphml", file="Gg30sm.graphml")

qq7Ig=sapply(Gg30peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Ig=log10(qq7Ig+1e-5)
qq7Im=sapply(Gm30peps, function(l){
  sum(qq7I[l]>0)/length(l)
})
qq7Im=log10(qq7Im+1e-5)


grps=names(vertex_attr(Gg30sm)$Group[[1]])
gic=as.numeric(substr(grps,4,4))>0
gia=as.numeric(substr(grps,5,5))>0
d=gic+gia

Gg60CA=t(sapply(vertex_attr(Gg30sm)$Group, function(x){
  x=as.data.frame(x)/d
  cntr=sum(x[gic,])
  apls=sum(x[gia,])
  return(c(cntr, apls))
}))
colnames(Gg60CA)=c("Control","APLS")
y=Gg60CA+0.1
chsqgCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=unique(2*(seq_along(chsqgCA[,1])[-1] %/% 2))
j1=j2-1
ii=which(chsqgCA[j2,4]<0.05)
Gg_good=rep(0,length(chsqgCA[,4]))
Gg_good[2*ii-1]=sign(chsqgCA[2*ii-1,4])
Gg_good=Gg_good[j1]

Gg60P=sapply(vertex_attr(Gg30sm)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})

grps=names(vertex_attr(Gm30sm)$Group[[1]])
mic=as.numeric(substr(grps,2,2))>0
mia=as.numeric(substr(grps,3,3))>0
d=mic+mia

Gm60CA=t(sapply(vertex_attr(Gm30sm)$Group, function(x){
  x=as.data.frame(x)/d
  cntr=sum(x[mic,])
  apls=sum(x[mia,])
  return(c(cntr, apls))
}))
colnames(Gm60CA)=c("Control","APLS")
y=Gm60CA+0.1
chsqmCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=unique(2*(seq_along(chsqmCA[,1])[-1] %/% 2))
j1=j2-1
ii=which(chsqmCA[j2,4]<0.05)
Gm_good=rep(0,length(chsqmCA[,4]))
Gm_good[2*ii-1]=sign(chsqmCA[2*ii-1,4])
Gm_good=Gm_good[j1]

Gm60P=sapply(vertex_attr(Gm30sm)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})

Gm30sm=set.vertex.attribute(Gm30sm,name="A/C",value=(Gm60CA[,2]+0.5)/(Gm60CA[,1]+0.5))
Gm30sm=set.vertex.attribute(Gm30sm,name="A/C_sig",value=Gm_good)
Gm30sm=set.vertex.attribute(Gm30sm,name="Pub",value=Gm60P)
Gm30sm=set.vertex.attribute(Gm30sm,name="Id",value=qq7Im)
Gg30sm=set.vertex.attribute(Gg30sm,name="A/C",value=(Gg60CA[,2]+0.5)/(Gg60CA[,1]+0.5))
Gg30sm=set.vertex.attribute(Gg30sm,name="A/C_sig",value=Gg_good)
Gg30sm=set.vertex.attribute(Gg30sm,name="Pub",value=Gg60P)
Gg30sm=set.vertex.attribute(Gg30sm,name="Id",value=qq7Ig)

write_graph(Gm30sm,format = "graphml", file="Gm30sm.graphml")
write_graph(Gg30sm,format = "graphml", file="Gg30sm.graphml")

x=cbind(rowAnys(as.matrix(fl[,2:3])),rowAnys(as.matrix(fl[,4:5])))
table(rowSums(x))

GMreun=unlist(c(Gm30peps[Gm_good!=0],Gg30peps[Gg_good!=0]))
table(table(GMreun))
nm=length(unlist(Gm30peps))
ng=length(unlist(Gg30peps))
ij=sample(ng,nm)
mc=unlist(Gm30peps[Gm_good==-1])
ma=unlist(Gm30peps[Gm_good==1])
gc=unlist(Gg30peps[ij][Gg_good[ij]==-1])
ga=unlist(Gg30peps[ij][Gg_good[ij]==1])
mgca=data.frame(IgM_contr=GMreun %in% mc,IgM_APLS=GMreun %in% ma,
                IgG_contr=GMreun %in% gc,IgG_APLS=GMreun %in% ga)
mgcavenn=euler(mgca, shape="ellipse")
plot(mgcavenn, quantities = T)

# The DEx clusters in IgG and IgM are anticorrelated - 1.4% peptides overlap
# vs 5.8% for the whole mixed library.
### Dex graph fused -------
GMreun=unique(GMreun)
Ggmdex=induced.subgraph(G, GMreun)
Ggmdexsm=clusterGraph(Ggmdex, resol=30)

load("GGMdexpepnames9047")
Ggmpeps=Gpeps
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="size",value = sqrt(vertex_attr(Ggmdexsm)$size))

grps=names(vertex_attr(Ggmdexsm)$Group[[1]])
mic=as.numeric(substr(grps,2,2))>0
mia=as.numeric(substr(grps,3,3))>0
gic=as.numeric(substr(grps,4,4))>0
gia=as.numeric(substr(grps,5,5))>0
dm=mic+mia
dg=gic+gia
dc=mic+gic
da=mia+gia
GgmdexsmCA=t(sapply(vertex_attr(Ggmdexsm)$Group, function(x){
  xm=as.data.frame(x)/dm
  cntrm=sum(x[mic])
  aplsm=sum(x[mia])
  xg=as.data.frame(x)/dg
  cntrg=sum(x[gic])
  aplsg=sum(x[gia])
  return(c(cntrm, aplsm,cntrg, aplsg))
}))
colnames(GgmdexsmCA)=c("Control_IgM","APLS_IgM","Control_IgG","APLS_IgG")
GgmdexsmP=sapply(vertex_attr(Ggmdexsm)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})
qq7Igmdex=sapply(Ggmpeps, function(l){
  sum(qq7I[l]>0)/length(l)
})
Gbgdex=sapply(Ggmpeps, function(p){
  median(bkgNN[p])
})

Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="Control_IgM", value=GgmdexsmCA[,1])
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="APLS_IgM", value=GgmdexsmCA[,2])
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="Control_IgG", value=GgmdexsmCA[,3])
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="APLS_IgG", value=GgmdexsmCA[,4])
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="Pub", value=GgmdexsmP)
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="Id", value=qq7Igmdex)
Ggmdexsm=set.vertex.attribute(Ggmdexsm, name="Bkg", value=Gbgdex)

x=grEmbed(Ggmdexsm)

write_graph(Ggmdexsm, format = "graphml", file="Ggmdexsm.graphml")

### Dex graph combined ----
dexmpeps=Gm30peps[Gm_good!=0]
dexmpeps=melt(dexmpeps)
x=dexmpeps[,2]
names(x)=dexmpeps[,1]
dexmpeps=x
dexgpeps=Gg30peps[Gg_good!=0]
dexgpeps=melt(dexgpeps)
x=dexgpeps[,2]
names(x)=dexgpeps[,1]
dexgpeps=x
clnm=Gm_good[Gm_good!=0]
names(clnm)=names(Gm30peps[Gm_good!=0])
clnm[clnm==-1]="C"
clnm[clnm==1]="A"
clng=Gg_good[Gg_good!=0]
names(clng)=names(Gg30peps[Gg_good!=0])
clng[clng==-1]="C"
clng[clng==1]="A"
y=rep(0,length(GMreun))
names(y)=GMreun
y1=y
y1[names(dexgpeps)]=dexgpeps
y2=y
y2[names(dexmpeps)]=dexmpeps
x=aggregate(GMreun, by=list(y1,y2), "list")
nms=paste("Cr",seq_along(x[,1]), sep="")
reunpeps=x[,3]  
names(reunpeps)=nms
reungroup=apply(x[,1:2],1,function(l){
  paste("G",clng[l[1]],"M",clnm[l[2]], sep="")
})

GGMdex=induced.subgraph(G, GMreun)
GGMdex=set.edge.attribute(GGMdex,name="weight", value=1)

v=melt(x[,3])
x0=v[,2]
names(x0)=v[,1]
v=x0
GGMdexsm=contract(GGMdex,mapping=v,vertex.attr.comb="ignore")
GGMdexsm=simplify(GGMdexsm, edge.attr.comb = "sum")
GMdexsmpeps=x[,3]
GGMdexsm=set_vertex_attr(GGMdexsm, "edges", value=sapply(GMdexsmpeps,function(x){
  g=induced_subgraph(G,x)
  length(E(g))
}))
grdex=vertex.attributes(GGMdex)$Group
names(grdex)=names(V(GGMdex))
pub=sapply(GMdexsmpeps, function(p){
  sum(substr(grdex[p],1,1)==1)/length(p)
})
id=sapply(GMdexsmpeps, function(p){
  sum(qq7I[p])/length(p)
})

GGMdexsm=set.vertex.attribute(GGMdexsm, name="size", value=lengths(GMdexsmpeps))
GGMdexsm=set.vertex.attribute(GGMdexsm, name="group", value=reungroup)
GGMdexsm=set.vertex.attribute(GGMdexsm, name="pub", value=pub)
GGMdexsm=set.vertex.attribute(GGMdexsm, name="id", value=id)

n=length(V(GGMdexsm))
A=matrix(0,n,n)
rownames(A)=names(V(GGMdexsm))
colnames(A)=names(V(GGMdexsm))
eG=apply(ends(GGMdexsm, E(GGMdexsm), E(GGMdexsm)),2,unlist) # matrix of edge ends
w=edge.attributes(GGMdexsm)$weight                  # number of edges between clusters
A[eG]=w                                         # adjacency mx of clusters
A=t(A)
A[eG]=w
Ein=vertex_attr(GGMdexsm)$edges                     # number of edges in each cluster
S=vertex_attr(GGMdexsm)$size                        # number of vertices in each cluster
ex=(Ein)/(S*(S-1)/2+1*(S==1))
ex=matrix(rep(ex,n),n,n)
m=(ex + t(ex))/2                                # mean density
w=(S%*%t(S))*m/(A+1*(A==0))                     # expected no of edges between/ actual
rownames(w)=names(V(GGMdexsm))                      # w is a dissimilarity measure
colnames(w)=names(V(GGMdexsm))
w=w[eG]

w=log10(w+min(w[w>0])*(w==0)) 
x=w[order(w)]
x=x[length(x)]-x[length(x)-1]
w=w-min(w)+x/2

GGMdexsm=set.edge.attribute(GGMdexsm, name="weight", value=w)

x=mst(GGMdexsm)
thr=max(E(x)$weight)
GGMdexsm=delete.edges(GGMdexsm,E(GGMdexsm)[E(GGMdexsm)$weight>thr])
x=edge_attr(GGMdexsm)$weight               # convert weights from dissimilarity to similarity
x=scale(1/x, center=F)
GGMdexsm=set.edge.attribute(GGMdexsm,"weight", value=x)

write_graph(GGMdexsm, format = "graphml", file="GGM_dexsm.graphml")
# Modularity scan =====

load("Gg")
cl=makeCluster(3)
clusterExport(cl, "Gg")
clusterEvalQ(cl, require(igraph))
mod=pbsapply(18+((1:9)*4), function(i){
  x=cluster_louvain(Gg, resolution = 60)
  c(i,max(x$modularity),length(table(x$membership)))
},cl=cl)
stopCluster(cl)

