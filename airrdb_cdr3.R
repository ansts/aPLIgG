require(pbapply)
require(parallel)
require(future.apply)
require(Biostrings)
require(stringi)
require(reshape2)
require(protr)
require(bioseq)
require(seqinr)
require(stringdist)
require(pheatmap)
require(furrr)

pth0="~/pth0="~/pth0="~/Documents/IgGenes/"
pth0=""
pth1="/media/anastas/Ext2/IgGenes/"
aa=AA_STANDARD

cpl=colorRampPalette(c("black","blue","green","yellow","red"))

JVergannaive=getJclean(pth0, "Vergani_B_naive/", f="ireceptor-public-archive.tsv")
JDeKoskymem2=getJclean(pth0, "", f="vdjserver.tsv")

save(JDeKoskymem2, file="JDeKoskymem2")
save(JDeKoskymem, file="JDeKoskymem")
save(JRobinsmem, file="JRobinsmem")
save(JMarkDavismem, file="JMarkDavismem")
save(JMarkDavisnaive, file="JMarkDavisnaive")
save(JRobinsnaive, file="JRobinsnaive")
save(JVergannaive, file="JVergannaive")

Jnaive=c(JRobinsnaive,JVergannaive,JMarkDavisnaive)
tJn=table(Jnaive)
save(tJn,file="tJn")
Jnpub=names(tJn)[tJn>1]
Jnpubhi=names(tJn)[tJn>5]
Jnaive=names(tJn)[tJn==1]


Jmem=c(JDeKoskymem2,JRobinsmem,JMarkDavismem)
tJm=table(Jmem)
save(tJm,file="tJm")
Jmpub=names(tJm)[tJm>1]
Jmpubhi=names(tJm)[tJm>5]
Jmem=names(tJm)[tJm==1]

Jmntrsct=intersect(Jnaive, Jmem)
Jnaive=setdiff(Jnaive,Jmntrsct)
Jmem=setdiff(Jmem,Jmntrsct)
Jpubttl=unique(c(Jmpub,Jnpub,Jmntrsct))

save(list=c("Jnaive","Jnpub","Jnpubhi"), file="Jn")
save(list=c("Jmem","Jmpub","Jmpubhi"), file="Jm")
save(Jmntrsct, file="Jmntrsct")
save(Jpubttl, file="Jpubttl")

rm(JDeKoskynaive,JRobinsnaive,JVergannaive,JMarkDavisnaive,JDeKoskymem,JDeKoskymem2,JRobinsmem,JMarkDavismem)

Jall7=t(qgrams(setdiff(Jnaive,Jmntrsct),c(Jnpub,Jmntrsct),Jnpubhi, setdiff(Jmem,Jmntrsct), c(Jmpub,Jmntrsct), Jmpubhi, q=7))
colnames(Jall7)=c("N","NP","NPhi","M","MP","MPhi")
save(Jall7,file="Jall7")

load("Jall7")
Jall7=scale(log10(Jall7+0.5))

J7hist=pbapply(Jall7,2,function(co) cut(co,100,labels = F))
J7nxm=table(as.data.frame(J7hist[,c(1,4)]))
pheatmap(log10(J7nxm+0.5), col=cpl(100), cluster_rows = F, cluster_cols = F)

J7npubxmpub=table(as.data.frame(J7hist[,c(2,5)]))
pheatmap(log10(J7npubxmpub+0.5), col=cpl(100), cluster_rows = F, cluster_cols = F)

J7npubhxmpubh=table(as.data.frame(J7hist[,c(3,6)]))
pheatmap(log10(J7npubhxmpubh+0.5), col=cpl(100), cluster_rows = F, cluster_cols = F)

J7nxpubT=table(as.data.frame(J7hist[,c(1,7)]))
pheatmap(log10(J7nxpubT+0.5), col=cpl(100), cluster_rows = F, cluster_cols = F)

J7mxpubT=table(as.data.frame(J7hist[,c(1,7)]))
pheatmap(log10(J7mxpubT+0.5), col=cpl(100), cluster_rows = F, cluster_cols = F)

Fr4="WGQGTLVTVSS"
Fr4_7=colnames(qgrams(Fr4,q=7))

# get VDJ calls from the data =-----------------------------------------

Pthsi=rbind(c("DeKosky_B_mem2/", "vdjserver.tsv"),
            c("DeKosky_B_mem/", "vdjserver.tsv"),
            c("Robins_B_mem/", "vdjserver.tsv"),
            c("MarkDavis_B_mem/", "vdjserver.tsv"),
            c("MarkDavis_B_naive/", "vdjserver.tsv"),
            c("Robins_B_naive/", "vdjserver.tsv"),
            c("Vergani_B_naive/", "ireceptor-public-archive.tsv"))


load("Jcl2_2")

L=unlist(Jcl2_2)

LOR=paste(L, collapse="|")

proct=proc.time()
vdjcl2_2n=grep(LOR,Jnaive, value = T)
print(proc.time()-proct)

proct=proc.time()
vdjcl2_2m=grep(LOR,Jmem, value = T)
print(proc.time()-proct)

proct=proc.time()
vdjcl2_2p=grep(LOR,Jpubttl, value = T)
print(proc.time()-proct)

L=c(vdjcl2_2n,vdjcl2_2m,vdjcl2_2p)

vdjgen2_2n=c()
for (i in 5:7) {
    P=Pthsi[i,]
    print(c(P[1],P[2]))
    vdjgen2_2n=c(vdjgen2_2n, list(getJgene(L=vdjcl2_2n, path0=pth1, path=P[1], f=P[2])))
}
save(vdjgen2_2n, file="vdjgen2_2n")

vdjgen2_2m=c()
for (i in 1:4) {
  P=Pthsi[i,]
  print(c(P[1],P[2]))
  vdjgen2_2m=c(vdjgen2_2m, list(getJgene(L=vdjcl2_2m, path0=pth1, path=P[1], f=P[2])))
}
save(vdjgen2_2m, file="vdjgen2_2m")

vdjgen2_2p=c()
for (i in 1:7) {
  P=Pthsi[i,]
  print(c(P[1],P[2]))
  vdjgen2_2p=c(vdjgen2_2p, list(getJgene(L=vdjcl2_2p, path0=pth1, path=P[1], f=P[2])))
}
save(vdjgen2_2p, file="vdjgen2_2p")

X=c()

Mall=c(vdjgen2_2n,vdjgen2_2m,vdjgen2_2p)

for (Mi in Mall){
    X=rbind(X,Mi)
}

vdjcl2_2=aggregate(X[,1], by=list(X[,2],X[,3],X[,4]), function(x) array(c(length(x),list(x)), dim=c(1,2)))
d2_2=table(unlist(strsplit(vdjcl2_2$Group.2, split=c("\\,| |/OR"))))
j2_2=table(unlist(strsplit(vdjcl2_2$Group.3, split=c("\\,| |/OR"))))

# IGHJ4 - YFDYWGQGTLVTVSS
# IGHJ6 - YYYYYGMDVWGQGTTVTVSS
# IGHJ5 - NWFDSWGQGTLVTVSS

# IGHD6-13 - GYSSSWY
# IGHD6-6 - EYSSSS
# IGHD2-2 - GYCSSTSCYA