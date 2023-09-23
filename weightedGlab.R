require(parallel)
require(pbapply)
require(igraph)
require(matrixStats)
require(reshape2)
require(stringdist)
require(heatmap3)
require(Biostrings)
require(uwot)
require(corrplot)
require(stringi)
require(qualV)
require(chisq.posthoc.test)

load("Ggiantw")
#G=induced_subgraph(G, components(Gw)$membership==1)



## IgM ----
load("Ggiantw")
vG=names(V(G))
vG=rownames(fls)[rownames(fls) %in% vG]
flss=fls[vG,]
gm=flss$cM|flss$aM
Gmw=induced.subgraph(G,gm)
save(Gmw, file="Gmw")
rm(G)
vGm=names(V(Gmw))
closeAllConnections()
# modularity scan

cl=makeCluster(4)
clusterExport(cl, c("G","vG"))
clusterEvalQ(cl, require(igraph))
modsc=pbsapply(c(5,5.5,6,6.5,7,7.5,8,9), function(resol){
  print(resol)
    x=cluster_leiden(G, objective_function = "modularity", resolution_parameter = resol, n_iterations = 3)
      x=x$membership
      m=modularity(G,as.numeric(as.factor(x)))
      n=length(table(x)[table(x)>1])
      BSmod=sapply(1:10, function(j){
        labv=rep("_",vcount(G))
        names(labv)=vG
        seed=sample(names(V(G)),n)
        seedn=paste("m",1:n,sep="")
        labv[seed]=seedn
        proct=proc.time()
        repeat{
          print(length(seed)/vcount(G))
          av=adjacent_vertices(G,V(G)[seed])
          for (i in sample(seedn)){
            ss=names(labv)[labv==i]
            jj=which((names(labv) %in% unique(unlist(sapply(av[ss],names), use.names = F))) & labv=="_")
            labv[jj]=i
          }
          print(table(labv))
          seed=V(G)[vG %in% names(labv)[labv!="_"] & !(vG %in% names(seed))] 
          if (sum(labv=="_")<100) break
        }
        print(proc.time()-proct)
        modularity(G,as.numeric(as.factor(labv[vG])))
    })
    res=m/median(BSmod)
    rm(m,BSmod,av,labv,seed,seedn,x,n)
    gc()
    return(res)
}, cl=cl)
closeAllConnections()


Gm39w=clusterGraph(Gmw, resol=39)
load("Gctrpeps")
Gm39wpeps=Gpeps
save(Gm39wpeps, file="Gm39wpeps")
rm(Gmw)

grps=names(vertex_attr(Gm39w)$Group[[1]])
cim=as.numeric(substr(grps,2,3))==10
aim=as.numeric(substr(grps,2,3))==1
caim=as.numeric(substr(grps,2,3))==11
d=cim+aim+2*caim

Gm39wCA=t(sapply(vertex_attr(Gm38w)$Group, function(x){
  x=as.data.frame(x)
  cntr=sum(x[cim,])
  apl=sum(x[aim,])
  all=sum(x[caim,])
  return(c(cntr, apl, all))
}))
colnames(Gm38wCA)=c("Control","APLS","All")
rownames(Gm38wCA)=names(V(Gm38w))

y=Gm38wCA
chsqmCA=chisq.posthoc.test(as.table(y), simulate.p.value=T)
j2=(1:(nrow(chsqmCA)/2))*2
j1=j2-1
Gmw_good=(chsqmCA[j2,3:5]<0.05)*sign(chsqmCA[j1,3:5])
rownames(Gmw_good)=chsqmCA$Dimension[j2]

Gm38wP=sapply(vertex_attr(Gm38w)$Group, function(x){
  x=as.data.frame(x)
  n=rownames(x)
  j=as.numeric(substr(n,1,1))>0
  return(sum(x[j,]))
})


## IgG ----
load("Ggiantw")
vG=names(V(G))
vG=rownames(fls)[rownames(fls) %in% vG]
flss=fls[vG,]
gg=flss$cG|flss$aG
Ggw=induced.subgraph(G,V(G)[gg])
save(Ggw, file="Ggw")
rm(G)
vGg=names(V(Ggw))
closeAllConnections()


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

# Gmat -----

load("Gmat")

vG=names(V(Gmat))
x=sample(vG, 4.5e5)
G=induced.subgraph(Gmat, names(V(Gmat)) %in% x)
vG=names(V(G))
rm(Gmat)
closeAllConnections()


# comparissons between weighted and unweighted clustering -----

p38=melt(Gm38peps)
p39=melt(Gm39peps)
p3839=array(0,dim = dim(p38))
rownames(p3839)=sort(p38$value)
x=p38$L1
names(x)=p38$value
p38=x
x=p39$L1
names(x)=p39$value
p39=x
p3839[,1]=p38[rownames(p3839)]
p3839[,2]=p39[rownames(p3839)]
colnames(p3839)=c("N", "W")
x=table(as.data.frame(p3839))
x=as.data.frame(x)

x=acast(data=x,W~N, value.var = "Freq")
pdf(file="p38vsp38_corrplot.pdf", width=15, height=15)
  corrplot(log10(x+0.5), is.corr = F, method="color", hclust.method = "ward.D2")
dev.off()