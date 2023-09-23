# Packages -----
require(igraph)
require(reshape2)
require(stringdist)
require(matrixStats)
require(Biostrings)
require(stringi)
require(qualV)
require(R.utils)
require(parallel)
require(pbapply)
require(future.apply)
require(combinat)
require(entropy)
require(gmp)
require(DEoptimR)
require(gtools)

load("lcslab_vars")
load("matpep")
aa=AA_STANDARD
x="ACDEFGA"
ij=t(combn(7,5))
q=consensusMatrix(AAStringSet(pepGpos), as.prob = T)
xs=unlist(strsplit(x,split=""))
lcss=apply(ij,1,function(ji) paste(xs[ji], collapse=""))
table(lcss)
cmplxlcs=matrix(0,nrow=7,ncol=21)
colnames(cmplxlcs)=1:21
rownames(cmplxlcs)=1:7

XY=pbsapply(1:1e6,function(i){
  x=sapply(1:7, function(i) sample(aa, 1, replace=T, prob=q[,i]))
  cmpl=length(table(x))
  lcss=apply(ij,1,function(ji) paste(x[ji], collapse=""))
  lcss=length(unique(lcss))
  cmplxlcs[cmpl,lcss] <<- cmplxlcs[cmpl,lcss]+1
  return(NULL)
})
closeAllConnections()

# Entropies ----

L2P=unlist(lcs2pep)

L=strsplit(lcs_lG,split="")
P=strsplit(pepGpos,split="")
cl=makeCluster(15)
clusterEvalQ(cl, require(entropy))
el=pbsapply(L,function(l) entropy(table(l), method="ML", unit="log2", verbose=F), cl=cl)
closeAllConnections()
cl=makeCluster(15)
clusterEvalQ(cl, require(entropy))
ep=pbsapply(P,function(l) entropy(table(l), method="ML", unit="log2", verbose=F), cl=cl)
closeAllConnections()
names(el)=lcs_lG
names(ep)=pepGpos
el=round(el,2)
ep=round(ep,2)

# N of combinations of all AA with repetitions of every kind for 7-mer peptides

thrcxN=c(choose(20,7),
      20*choose(19,5),
      choose(20,2)*choose(18,3),
      choose(20,3)*17,
      20*choose(19,4),
      choose(20,2)*18,
      20*choose(19,3),
      20*choose(19,2),
      20*19,
      choose(20,2)*2*choose(18,2),
      choose(20,3)*3,
      choose(20,2)*2*18,
      choose(20,2)*2,
      choose(20,2)*2,
      20) 

n0=c(7,6,5,4,5,3,4,3,2,4,3,3,2,2,1) # Number of elements per string
reps=list(0,2,c(2,2),c(2,2,2),3,c(3,3),4,5,6,c(2,3),c(2,2,3),c(2,4),c(3,4),c(2,5),7)

# Number of permutations with repetitions of different kinds
ns=sapply(1:15, function(i){
  factorial(7)/(prod(sapply(reps[[i]],factorial)))
})

# Permutations in 7-mers with repetitions by combinations out of 20 letters give all possible permutations
sum(ns*thrcxN)

(ns*thrcxN)/sum(ns*thrcxN)
round((ns*thrcxN)/sum(ns*thrcxN)*100,6)

# n for dmultinom                                # rep       unique     ttl
thrcmplx=t(array(c(c(rep(1,7),rep(0,13)),         # 0           7         7
              c(rep(1,5),2,rep(0,14)),            # 2           5         6
              c(rep(1,3),2,2,rep(0,15)),          # 2,2         3         5
              c(1,rep(2,3),rep(0,16)),            # 2,2,2       1         4
              c(rep(1,4),3,rep(0,15)),            # 3           4         5  
              c(1,3,3,rep(0,17)),                 # 3,3         1         3
              c(rep(1,3),4,rep(0,16)),            # 4           3         4
              c(1,1,5,rep(0,17)),                 # 5           2         3  
              c(1,6,rep(0,18)),                   # 6           1         2
              c(1,1,2,3,rep(0,16)),               # 2,3         2         4
              c(2,2,3,rep(0,17)),                 # 2,2,3       0         3
              c(1,2,4,rep(0,17)),                 # 2,4         1         3
              c(3,4,rep(0,18)),                   # 3,4         0         2 
              c(2,5,rep(0,18)),                   # 2,5         0         2
              c(7,rep(0,19))), dim=c(20,15)))     # 7           0         1 

# likelihood of each type of repetition with uniform distribution of AA frequencies
X=apply(thrcmplx,1,function(n) dmultinom(n,prob=rep(0.05,20)))

round((X*thrcxN),8)==round((ns*thrcxN)/sum(ns*thrcxN),8)
round(X,8)==round(ns/sum(ns*thrcxN),8)

# LCS distributions

p_nlbye=pbsapply(sort(unique(ep), decreasing = T), function(en){
  p=t(sapply(names(ep[ep==en]), function(pp) unlist(strsplit(pp, split=""))))
  ij=combn(7,5)
  l=apply(p,1,function(pp) length(unique(apply(array(pp[ij], dim=c(5,21)),2,paste,collapse=""))))
  table(l)/sum(table(l))
})

names(p_nlbye)=sort(unique(ep), decreasing = T)
wmnlbye=sapply(p_nlbye,function(x) sum(x*as.numeric(names(x))))

# N of combinations of the repetitions ----

tchoose=c(1,
          choose(7,2),
          (choose(7,2)*choose(5,2))/2,
          (choose(7,2)*choose(5,2)*choose(3,2))/3,
          choose(7,3),
          (choose(7,3)*choose(4,2)),
          (choose(7,3)*choose(4,3))/2,
          2*choose(7,4),
          (choose(7,4)*choose(3,2)),
          2*choose(7,5),
          7,1)


vgrep=Vectorize("grep",vectorize.args = "pattern")

# N of LCS by type of scheme ----

p=rndpp
px=sapply(strsplit(p, split=""), unlist)
cl=makeCluster(15)
clusterExport(cl,c("vgrep"))
clusterEvalQ(cl, require(reshape2))
sch=pbapply(px,2,function(p) {
  s=vgrep(unique(p),p)
  names(s)=NULL
  if (is.array(s)) return("1111111")
  if (length(s)<7) {
    ij=melt(s)
    ij=ij[order(ij$value),2]
    return(paste(ij, collapse = ""))
  } else {return(paste(s, collapse = ""))}
},cl=cl)
closeAllConnections()

sch0=sort(unique(sch))

ij=combn(7,5)
lcsch0=sapply(sch0, function(s) length(unique(apply(ij,2,function(j) paste(unlist(strsplit(s,split=""))[j], collapse="")))))
tlcsch0=aggregate(names(lcsch0), by=list(lcsch0), "list")
lcschn=names(lcsch0)

# lcschn are all 877 possible arrangements of peptide seqs with residues labeled 
# from left to right with numbers 1:7 and with repetitions labeled with the same 
# number - the number of the first occurrence of the repeated residue. 
# The homoheptamer is omitted as too rare and trivial. Lcsch are all unique 
# subsequences of each lcschn recoded to unique patterns of repetition rather than 
# the successive numbers of residues in the peptide. x - the matrix of the 5-letter
# subsequences of each lcschn. tx - their tables indicating how many repetitions 
# exist of each label (position related numbers) with z - individual tables. 
# The names of z (nm) are the unique positions found in individual lcschn and,
# respectively, in their subsequences. 
# t0 - where the recoded x is kept. The code is 0 - unique letters, B,D - 
# possible pairs, T - triplicate, Q- quadriplicate, and P - all 5 letters 
# are the same.

# Calculating tables of P ------

X1=sample(lcschn,30)
X2=sample(lcschn,30)
ij5=combn(7,5)
ij6=combn(7,6)
Lcsch=pblapply(lcschn, function(n){
  n=unlist(strsplit(n, split=""))

  cmbix=apply(ij5,2,function(i) paste(i, collapse=""))
  al5=apply(ij5,2,function(i) (n[i]))
  aln5=apply(al5,2,paste, collapse="")
  colnames(al5)=aln5
  x5=t(unique(t(al5)))
  nl=aggregate(1:21, by=list(aln5), "length")
  g=nl$x
  names(g)=nl$Group.1
  x5=x5[,names(g)]
  ni=aggregate(cmbix, by=list(aln5), "list")
  v5=ni$x
  names(v5)=ni$Group.1
  if (is.null(dim(x5))) tx5=table(x5) else tx5=lapply(seq_along(x5[1,]), function(i) table(x5[,i]))
  
  tx5=sapply(seq_along(tx5), function(j) {
    t0=rep(0,5)
    names(t0)=1:5
    z=tx5[[j]]
    nm=as.numeric(names(z))
    if (length(z)==1) t0=rep("P",5)
    if (length(z)==2) {   
      if (max(z)==4) {                     
        t0[x5[,j]==nm[z==4]]="Q"            
      }                                    
      if (max(z)==3) {
        t0[x5[,j]==nm[z==3]]="T"
        t0[x5[,j]==nm[z==2]]="B"
      }
    }
    if (length(z)==3) {
      if (max(z)==3) {
        t0[x5[,j]==nm[z==3]]="T"
      }
      if (max(z)==2) {
        nj=as.character(unique(x5[,j])[unique(x5[,j]) %in% nm[z==2]])
        for (i in 1:2){
          t0[x5[,j]==nj[i]]=c("B","D")[i]
        }
      }
    }
    if (length(z)==4) {
      t0[x5[,j]==nm[z==2]]="B"
    }
    return(paste(t0,collapse=""))
  })
  v5=lapply(v5, function(x){
    if (length(x)==1) return(x)
    x=sapply(x,function(l) unlist(strsplit(l, split="")))
    print(x)
    tx=table(c(x))
    tx=names(tx)[tx>1]
    jbad=which(colnames(x) %in% 
            apply(
              combn(length(tx),5),2,function(cmb) 
              paste(tx[cmb], collapse="")))
    }
    print(jbad)
    return(x[-jbad])
  })
  y5=cbind(tx5,g)
  y5=aggregate(cbind(rownames(y5),y5[,2],v5[rownames(y5)]), by=list(y5[,1]), "list")
  y5=data.frame(y5, N=sapply(y5$V2,function(z) sum(as.numeric(z))), n=sapply(y5$V1,function(z) length(z)))
  colnames(y5)[1:4]=c("Reptrn0","Reptrn1","n_1per0","indices")
  
  return(y5)
})
names(Lcsch)=lcschn

thprob=thProb(X1,X2)
ij0=which(thprob==0, arr.ind=T)
ij0=cbind(X2[ij0[,1]],X1[ij0[,2]])
ijp=which(thprob!=0, arr.ind=T)

jj=combn(877,2)
jj=apply(jj,2,sort)
jj=t(unique(t(jj)))
pij=pbapply(jj, 2, function(i) {
  thProb(lcschn[i[1]],lcschn[i[2]])
})
jj=rbind(jj,pij)
jj=jj[,order(pij)]
# jj=rbind(jj,pij[order(pij)])
ij=combn(7,5)

#clrs=(probBs6==0)*1
plot(probBs,thprob, log="xy", pch=16, cex=.5, col=rgb(0,0,0,0.3))
lmprobs=lm(log10(c(thprob+5e-6))~log10(c(probBs+5e-6)))
summary(lmprobs)
lmprobs$residuals
x=melt(thprob)
#x=x[x[,3]>0,]
x=cbind(x,lmprobs$residuals)
x=x[order(x[,4]),]
y=sapply(seq_along(x[,1]),function(i) {
  x1=Lcsch[[as.character(x[i,1])]]
  x2=Lcsch[[as.character(x[i,2])]]
  z=intersect(x1[,1],x2[,1])
  z1=x1[x1[,1] %in% z,2]
  z2=x2[x2[,1] %in% z,2]
  z3=length(intersect(z1,z2))
  unlist(c(as.character(x[i,1:2]),x[i,3],z,z1,"-___",z2,z3))
})

# Calculate complete probability table ----

jj=combn(877,2)
jj=cbind(jj,rbind(1:877,1:877))
cl=makeCluster(16)
clusterExport(cl, c("lcschn", "thProb1", "Lcsch"))
clusterEvalQ(cl, require(reshape2))
pij1=pbapply(jj, 2, function(i) {
  thProb1(lcschn[i[1]],lcschn[i[2]])
}, cl=cl)
closeAllConnections()
pij1=data.frame(S1=lcschn[jj[1,]],S2=lcschn[jj[2,]],P=pij1, stringsAsFactors = F)

# Matochko graph -----

load("matpep")
Gmat=simplify(adjL2Glite(adjL(matpep)))
save(Gmat, file="Gmat")
load("Gmat")
vrepptrnalize=Vectorize(repptrnalize,"s")

lcschnf=sapply(lcschn, function(x) {
  n=max(as.numeric(unlist(strsplit(x, split=""))))
  (20^-7)*factorial(20)/factorial(20-n)
})


lcschent=sapply(lcschn, function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))

# Small graph ----

lcschlt=sapply(lcschn, function(s) length(table(unlist(strsplit(s,split="")))))

# cl=makeCluster(15)
# clusterExport(cl, "aa")
# rndscflat=c(pbsapply(lcschn, function(s) {
#   x=as.numeric(unlist(strsplit(s,split="")))
#   u=length(unique(x))
#   sapply(1:100,function(j) {
#     paste(sample(aa,u)[x],collapse="")
#   })
# }, cl=cl))
# closeAllConnections()
# 
# cl=makeCluster(6)
# clusterExport(cl, c("lcschlt","lcschn", "aa"))
# rndlcflat=unlist(pbsapply(2:7, function(i) {
#   sc=lcschn[lcschlt==i]
#   n=length(sc)
#   sapply(sc, function(sci){
#     sapply(1:(16000 %/% n),function(j) {
#       s=sample(aa,i)
#       x=as.numeric(unlist(strsplit(sci, split="")))
#       paste(s[x],collapse="")
#     })
#   })
# }, cl=cl))
# closeAllConnections()
# 

x=sapply(1:7,function(z){
  y=rep(aa,each=10000)
  sample(y)
})
x=apply(x,1,function(l) paste(l,collapse="",sep=""))
ss=x

#ss=sample(names(V(Gmat)),250000)
#ss=sapply(1:250000, function(i) paste(sample(aa,7,replace = T),collapse="", sep=""))
#ss=sample(rndscflat)
 # core=sample(V(Gmat),1)
 # egg=names(ego(Gmat,2,core)[[1]])
 # g=induced.subgraph(Gmat,egg)
# 
# al=adjL(ss)
# g=adjL2Glite(al)
# g=induced_subgraph(g, V(g)[components(g)$membership==1])


load("Gmat")
g=Gmat
Gmat=NULL
rm(Gmat)
gc()

vg=names(V(g))
ltrcnt=pbsapply(vg,function(x) length(table(unlist(strsplit(x, split="")))))
dg=degree(g)
vrp=vrepptrnalize(vg)
esg=ends(g, E(g))

x=rep(0,length(lcschn))
names(x)=lcschn
tmpn=table(vrp)
x[names(tmpn)]=tmpn
tmpn=x

p1=tmpn%x%tmpn
p1=array(p1,dim=c(length(tmpn),length(tmpn)), dimnames = list(names(tmpn),names(tmpn)))
diag(p1)=tmpn*(tmpn-1)/2

esp=cbind(vrp[esg[,1]],vrp[esg[,2]])
rm(esg)

pwm=consensusMatrix(AAStringSet(vg))/length(vg)

lcschnfmat=lcschPwmp(L=lcschn)
save(lcschnfmat,file="lcschnfmat")
plot(lcschnf, lcschnfmat, xlim=c(1e-8,1), ylim=c(1e-8,1), log="xy")
lines(c(1e-8,1),c(1e-8,1))

# x=apply(pwm+3e-4,2,function(l) sum(l*log2(l)))
# x=x/20
# y=0.05*log2(0.05)
# c0=2^y
# exp(lambertWn(log(c0)))
# xx=exp(lambertWn(log(2^x)))
# plike=mean(1/xx)

cl=makeCluster(16)
clusterExport(cl, c("lcschn", "thprobEK", "Lcsch7"), envir = environment())
Bymg=pbsapply(lcschn, function(sc1) {
        sapply(lcschn, function(sc2) {
          thprobEK(sc1,sc2)
        })
}, cl=cl)
closeAllConnections()

Bymg=(Bymg+t(Bymg))/2

lcschlt_n=factorial(20)/factorial(20-lcschlt)
p1_=lcschlt_n%x%lcschlt_n
p1_=array(p1_,dim=dim(Bymg), dimnames=dimnames(Bymg))
diag(p1_)=lcschlt_n*(lcschlt_n-1)/2
lcschlt_r=(p1*20^7*(20^7-1))/(2*p1_*ecount(g))

gentr=pbsapply(names(V(g)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
pentr=pbsapply(lcschn,function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))

#save(esp, file="esp")
el=esp
cl=makeCluster(16)
el=t(pbapply(el,1,sort,cl=cl))
closeAllConnections()

vagg=aggregate(el[,1], by=list(el[,1],el[,2]), "length")
ij=as.matrix(vagg[,1:2])
y=as.numeric(vagg$x)
vaggm=array(0, dim=dim(Bymg), dimnames = dimnames(Bymg))
vaggm[ij]=y
vaggm=vaggm+t(vaggm)
diag(vaggm)=diag(vaggm)/2
rv=rowSums(vaggm)
vaggmm=pbsapply(rv, function(n1){
          sapply(rv, function(n2){
            n1+n2
          })
})
vaggmm=vaggmm-vaggm

nsqpp=pbsapply(seq_along(tmpn), function(i1){
      sapply(seq_along(tmpn), function(i2){
        tmpn[i1]+tmpn[i2]
  })
})
diag(nsqpp)=diag(nsqpp)/2
dimnames(nsqpp)=dimnames(vaggmm)

image(log10(vaggm[order(pentr),order(pentr)]/p1[order(pentr),order(pentr)]), col=cpl1(2000))
image((vaggm[order(pentr),order(pentr)]/p1[order(pentr),order(pentr)])>0)
corrplot(vaggm>0, method = "color", is.corr = F, order="AOE",hclust.method = "ward.D2")

# empirical
# pseudocounts have to be added to vaggmm and nsqpp
nes=ecount(g)
p4=aggregate((dg), by=list(vrp), "mean")
x=p4$x
names(x)=p4$Group.1
p4=x
p4=pbsapply(names(p4), function(n1){
      sapply(names(p4), function(n2){
          if (n1!=n2) p4[n1]*p4[n2] else p4[n1]*(p4[n1]-1)/2
      })
})
rownames(p4)=colnames(p4)
x=array(0,dim=dim(Bymg), dimnames=dimnames(Bymg))
x[rownames(p4),colnames(p4)]=p4
gmatwtmx=(x)*nsqpp/(2*vaggmm)

wt=1/gmatwtmx[el]

# theoretical ?
p4t=Bymg*p1
p4t=rowSums(p4t)/(tmpn)
p4t=pbsapply(names(p4t), function(n1){
  sapply(names(p4t), function(n2){
    if (n1!=n2) p4t[n1]*p4t[n2] else p4t[n1]*(p4t[n1]-1)/2
  })
})
rownames(p4t)=colnames(p4t)
#wt=(p4t[el]*nsqpp[el])/(2*vaggmm[el])

g=set_edge_attr(g, "weight", value=wt)

boxplot(log10(strength(g))~round(gentr,2), notch=T, xlab="Entropy [bits]", ylab="Log Strength")
boxplot(log10(dg)~round(gentr,2), notch=T, xlab="Entropy [bits]", ylab="Log Strength")
(log10(strength(g))~ltrcnt, notch=T)
boxplot(log10(dg)~ltrcnt, notch=T)
boxplot(log10(strength(g))~vrp[order(gentr)], notch=T)
boxplot(log10(dg)~vrp[order(gentr)], notch=T)

bigX=aggregate(Yg, by=list(esp[,1],esp[,2]), "mean")
bigX[,1]=lcschent[bigX[,1]]
bigX[,2]=lcschent[bigX[,2]]
#bigX=aggregate(bigX[,3], by=list(bigX[,1], bigX[,2]), "c")

clgw=cluster_leiden(g, objective_function = "modularity", resolution_parameter = 10)
clg=cluster_leiden(g, objective_function = "modularity", weights=rep(1, ecount(g)), resolution_parameter = 10)
clcegv=pbsapply(1:clg$nb_clusters, function(i){
  gg=induced_subgraph(g, V(g)[clg$membership==i])
  entropy(table(unlist(strsplit(names(which.max(eigen_centrality(gg)$vector)), split=""))), method="ML")
})
clcegvw=pbsapply(1:clgw$nb_clusters, function(i){
  gg=induced_subgraph(g, V(g)[clgw$membership==i])
  entropy(table(unlist(strsplit(names(which.max(eigen_centrality(gg)$vector)), split=""))), method="ML")
})

hist(clcegv, xlim=c(0,2))
hist(clcegvw, xlim=c(0,2), main="Weights")
mean(clcegv)
mean(clcegvw)
mean(gentr)

ecg=eigen_centrality(g)


# Igome graph -----

load("G")
dG=degree(G)
ventr=pbsapply(names(V(G)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
boxplot((dG)~round(ventr[names(dG)],2), notch=T, varwidth=T)

clG=cluster_leiden(G, objective_function = "modularity", resolution_parameter = 10)

cntrv=sapply(1:146, function(i){
  g=induced_subgraph(G, V(G)[clG$membership==i])
  V(g)[which.max(eigen_centrality(g)$vector)]
})

ventrc=pbsapply(names(cntrv),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
hist(ventr,xlim=c(0,2))
par(new=T)
hist(ventrc,xlim=c(0,2), col=rgb(1,0,0,0.3))
par(new=F)

esG=ends(G, E(G))
vG=names(V(G))
vrp=vrepptrnalize(vG)
x=rep(0,length(lcschn))
names(x)=lcschn
tmpn=table(vrp)
x[names(tmpn)]=tmpn
tmpn=x
esp=cbind(vrp[esG[,1]],vrp[esG[,2]])
rm(esG)
closeAllConnections()
el=esp
cl=makeCluster(16)
el=t(pbapply(el,1,sort,cl=cl))
rm(esp)
closeAllConnections()
load("gmatwtmx")
wt=1/gmatwtmx[el]
gentr=pbsapply(names(V(G)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))


G=set_edge_attr(G, "weight", value=wt)
ij=sapply(unique(gentr), function(n) sample(which(gentr==n),500, replace = T))
boxplot(log10(strength(G)[ij])~round(gentr[ij],2), xlab="Entropy [bits]", ylab="Log Strength")
boxplot(log10(dG[ij])~round(gentr[ij],2), xlab="Entropy [bits]", ylab="Log Strength")

X=lcschnf*vcount(G)
table(Expect=X<1,Actual=tmpn==0)
table(X[X<1],tmpn[X<1])

# define pseudocounts for 3 schemes which are 0 in Gmat and not 0 in G (one is GGGGGGG)