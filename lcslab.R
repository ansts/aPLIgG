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
tabentr=array(0, dim=c(length(unique(el)),length(unique(ep))))
rownames(tabentr)=sort(unique(el))
colnames(tabentr)=sort(unique(ep))

tbntr=pbsapply(seq_along(L2P), function(i){
  if ((i %% 10000)==0) print(list(i/length(L2P),tabentr))
  n=L2P[i]
  p=names(L2P)[i]
  lp=strsplit(p,split="\\.")[[1]]
  tabentr[as.character(el[lp[1]]),as.character(ep[lp[2]])]<<-tabentr[as.character(el[lp[1]]),as.character(ep[lp[2]])]+n
  return(NULL)
})

ptrentr_p=lapply(sort(unique(ep)), function(en) c(en,head(names(ep[ep==en]))))
x=c("XXXXXX_","XXXXXYY","XXXXYYY","XXXXX__","XXXXYY_","XXXYYY_","XXXYYZZ","XXXX___","XXYYY__","XXYYZZ_","XXX____","XXYY___","XX_____","_______")
ptrentr_p=lapply(seq_along(ptrentr_p),function(i) c(ptrentr_p[[i]][1],x[i]))
ptrentr_p=as.data.frame(t(as.data.frame(ptrentr_p)))
x=as.numeric(ptrentr_p[,1])
names(x)=ptrentr_p[,2]
ptrentr_p=x
ptrentr_l=lapply(sort(unique(el)), function(en) c(en,head(names(el[el==en]))))
x=c("XXXXX","XXXXXY","XXXXY","XXXXYY","XXXYY","XXXYYY","XXXX__","XXX__","XXXYY_","XXYY_","XXYYZZ","XXX___","XXYY__","XX___","XX____","_____","______")
ptrentr_l=lapply(seq_along(ptrentr_l),function(i) c(ptrentr_l[[i]][1],x[i]))
ptrentr_l=as.data.frame(t(as.data.frame(ptrentr_l)))
x=as.numeric(ptrentr_l[,1])
names(x)=ptrentr_l[,2]
ptrentr_l=x



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

# Bootstrap of entropies of peptides and edge lcs

lrnglcsBS=c()
for (z in 1:20){
  proct=proc.time()
  
  # cl=makeCluster(14)
  # clusterExport(cl,c("aa","q"))
  # rndpp=pbsapply(1:100000, function(i){
  #   paste(sapply(1:7, function(i) sample(aa, 1, replace=T, prob=q[,i])), collapse="")
  # },cl=cl)
  # closeAllConnections()
  
  # rndpp=sample(matpep,1e5)
  
  rndpp=rndpp[-duplicated(rndpp)]
  AL=adjL(rndpp)
  rndppG=adjL2Glite(AL)
  rndppG=simplify(rndppG)
  vrng=names(V(rndppG))
  drng=degree(rndppG)
  ecrng=eigen_centrality(rndppG)
  erng=ends(rndppG, E(rndppG))
  cmpxrng=sapply(vrng, function(p) length(table(unlist(strsplit(p, split="")))))
  # boxplot(log10(drng)~cmpxrng, notch=T)
  # boxplot(ecrng$vector~cmpxrng, notch=T)
  # boxplot(log10(ecrng$vector/(drng^2))~cmpxrng, notch=T)
  # 
  cl=makeCluster(14)
  clusterEvalQ(cl, require(qualV))
  lcsrng=pbapply(erng,1,function(l) paste(LCS(unlist(strsplit(l[1], split="")),unlist(strsplit(l[2], split="")))$LCS, collapse=""), cl=cl)
  closeAllConnections()
  
  #lcsrng5=lcsrng[nchar(lcsrng)==5]
  cl=makeCluster(14)
  cmpxlcsrng=pbsapply(lcsrng, function(p) length(table(unlist(strsplit(p, split="")))), cl=cl)
  closeAllConnections()
  
  lrnG=make_line_graph(rndppG)
  lrnG=set_vertex_attr(lrnG, name="lcs", value=lcsrng)
  lrnG=set_vertex_attr(lrnG, name="N", value=1)
  lrnG=set_edge_attr(lrnG, name="weight", value=1)
  x=as.numeric(as.factor(vertex_attr(lrnG)$lcs))
  lrnG=contract(lrnG,x, vertex.attr.comb = list(N ="sum", lcs="first"))
  lrnG=simplify(lrnG, edge.attr.comb = "sum")
  lrnG=set_vertex_attr(lrnG, name="name", value=vertex_attr(lrnG)$lcs)
  lrnG=delete_vertex_attr(lrnG,"lcs")
  vlrng=vertex_attr(lrnG)$N
  names(vlrng)=vertex_attr(lrnG)$name
  gc()
  lrngcl=cluster_leiden(lrnG, objective_function = "modularity", resolution_parameter = 2)
  
  cl=makeCluster(14)
  clusterExport(cl, c("lrnG","lrngcl"), envir = environment())
  clusterEvalQ(cl, require(igraph))
  eclrngcl=pbsapply(1:29, function(i){
    g=induced.subgraph(lrnG,V(lrnG)[lrngcl$membership==i])
    n=names(V(g))[which.max(eigen_centrality(g)$vector)]
    names(ego(g,1,n)[[1]])
  }, cl=cl)
  closeAllConnections()
  eclrngcl=unique(unlist(eclrngcl))
  cmpxeclrngcl=sapply(eclrngcl, function(p) length(table(unlist(strsplit(p, split="")))))
  tl=table(cmpxlcsrng)/sum(table(cmpxlcsrng))
  tecl=table(cmpxeclrngcl)/sum(table(cmpxeclrngcl))
  # plot(as.numeric(names(tl)),tl, xlim=c(1,7), ylim=c(0,1), ty="b")
  # par(new=T)
  # plot(as.numeric(names(tecl)),tecl, xlim=c(1,7), ylim=c(0,1), col=rgb(1,0,0,0.3), ty="b")
  # par(new=F)
  t0=rep(0,5)
  names(t0)=1:5
  t0i=t0
  t0i[names(tl)]=tl
  tl=t0i
  t0i=t0
  t0i[names(tecl)]=tecl
  tecl=t0i
  print(proc.time()-proct)
  lrnglcsBS=cbind(lrnglcsBS,c(tl,tecl))
}
save(rndpp, rndppG, lcsrng, cmpxlcsrng,cmpxrng,lrnG, file="lcslab_vars")

pwmlcs=consensusMatrix(AAStringSet(lcsrng[nchar(lcsrng)==6]), as.prob = T)
cmpxvrng=sapply(vrng, function(p) length(table(unlist(strsplit(p, split="")))))
vrng5mer=t(qgrams(vrng,q=5))
vrng5mer=unlist(sapply(rownames(vrng5mer),function(n) rep(n,vrng5mer[n,])))
cmpxvrng5mer=sapply(vrng5mer, function(p) length(table(unlist(strsplit(p, split="")))))
tv5m=table(cmpxvrng5mer)/sum(table(cmpxvrng5mer))
tl5=table(cmpxlcsrng[nchar(names(cmpxlcsrng))==5])/sum(table(cmpxlcsrng[nchar(names(cmpxlcsrng))==5]))
plot(as.numeric(names(tv5m)),tv5m, ty="b", xlim=c(1,5.5))
par(new=T)
plot(as.numeric(names(tl5)),tl5,col=2, ty="b", xlim=c(1,5.5))
par(new=F)



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

ALmatp=adjL(matpep)
Gmat=adjL2Glite(ALmatp)
vGmat=names(V(Gmat))
es=ends(Gmat,E(Gmat))
schms=c()
proct=proc.time()
jj=sample(seq_along(es[,1]))
i=1
repeat {
  p=es[jj[i],]
  px=sapply(strsplit(p, split=""), unlist)
  sch=apply(px,2,function(p) {
    s=vgrep(unique(p),p)
    names(s)=NULL
    if (length(s)<7) {
      ij=melt(s)
      ij=ij[order(ij$value),2]
      return(ij)
    } else {return(s)}
  })
  s=paste(apply(sch,2,function(l) paste(l, collapse = "")), collapse=".")
  if (length(grep(s,schms))==0) {
    schms=c(schms,s)
    prob=mean(pbsapply(1:10,function(j){
      psch=apply(sch,2,function(ij){
        n=length(unique(ij))
        sapply(1:10000,function(i){
          q=(sample(aa,n))
          paste(q[ij],collapse="")
        })
      })
      sum(stringdistmatrix(psch[,1],psch[,2], method="lcs",)<5)/1e8
    }))
  }
  print(list(i,length(schms),proc.time()-proct))
  i=i+1
}



# lcschn are all 876 possible arrangements of peptide seqs with residues labeled 
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


# just a random graph of sequences -----

cl=makeCluster(18)
clusterExport(cl,"aa")
rndpp=t(pbsapply(1:500000, function(i) {
  x=sample(aa,7,replace = T)
  paste(x,collapse="")
}, cl=cl))
closeAllConnections()

Grnd=simplify(adjL2Glite(adjL(rndpp)))
singl=names(V(Grnd)[components(Grnd)$membership>1])
dGrnd=degree(Grnd)
ventr=pbsapply(names(V(Grnd)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
boxplot(log10(dGrnd)~round(ventr,2), notch=T, varwidth=T)


clGrnd=cluster_leiden(Grnd, objective_function = "modularity", resolution_parameter = 10)

cntrv=sapply(1:391, function(i){
  g=induced_subgraph(Grnd, V(Grnd)[clGrnd$membership==i])
  V(g)[which.max(eigen_centrality(g)$vector)]
})

ventrc=pbsapply(names(cntrv),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
hist(ventr,xlim=c(0,2))
par(new=T)
hist(ventrc,xlim=c(0,2), col=rgb(1,0,0,0.3))
par(new=F)

Ghom=set_edge_attr(Ghom,"weight", value=1)
boxplot(log10(strength(Ghom))~round(ventr,2), notch=T)
boxplot(log10(degree(Ghom)+0.5)~round(ventr,2), notch=T)

# Homogeneous random sequence graph ? ----

cl=makeCluster(18)
clusterExport(cl,"aa")
rndpp0=t(pbsapply(1:10000000, function(i) {
  x=sample(aa,7,replace = T)
  tx=paste(sort(table(x)), collapse="")
  x=paste(x,collapse="")
  c(x,tx)
}, cl=cl))
closeAllConnections()

ptx0=sort(table(rndpp0[,2])/1e7, decreasing = T)
tx0=names(ptx0)
rm(rndpp0)
gc()

cl=makeCluster(15)
clusterExport(cl, c("tx0", "aa"))
rndppflat=pbsapply(1:500000, function(i) {
  sc=as.numeric(unlist(strsplit(sample(tx0,1),split="")))
  x=sample(aa,length(sc))
  paste(sample(unlist(mapply(rep,x,each=sc))),collapse="")
}, cl=cl)
closeAllConnections()
rndppflat=c(rndppflat,"GGGGGGG")

Ghom=simplify(adjL2Glite(adjL(rndppflat)))
singl=names(V(Ghom)[components(Ghom)$membership>1])
dGhom=degree(Ghom)
ventr=pbsapply(names(V(Ghom)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
boxplot(log10(dGhom)~round(ventr,2), notch=T, varwidth=T)
plot(table(ventr))


clGhom=cluster_leiden(Ghom, objective_function = "modularity", resolution_parameter = 10)

cntrv=sapply(1:208, function(i){
  g=induced_subgraph(Ghom, V(Ghom)[clGhom$membership==i])
  V(g)[which.max(eigen_centrality(g)$vector)]
})

ventrc=pbsapply(names(cntrv),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
plot(table((ventr)))
par(new=T)
plot(table(ventrc),xlim=c(0,2), col=rgb(1,0,0,0.3))
par(new=F)

save(Ghom, file="Ghom")
save(Grnd, file="Grnd")

# Matochko aa distr random -----
load("matpep")
mmatp=t(sapply(matpep, function(p) unlist(strsplit(p,split=""))))
rndmat=apply(apply(mmatp,2,sample),1,paste, collapse="")
Grmat=simplify(adjL2Glite(adjL(rndmat)))
save(Grmat, file="Grmat")

load("Grmat")
singl=names(V(Grmat)[components(Grmat)$membership>1])
Grmat=induced_subgraph(Grmat, V(Grmat)[!(names(V(Grmat)) %in% singls]))
dGrmat=degree(Grmat)
ventrr=pbsapply(names(V(Grmat)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
boxplot(log10(dGrmat)~round(ventrr,2), notch=T)
plot(table(ventrr))
x=aggregate(dGmat,by=list(ventr), median)
y=aggregate(dGrmat,by=list(ventrr), median)
ventrboth=cbind(x[,2],y[,2])
plot(ventrboth, cex=0)
text(ventrboth, labels = round(x$Group.1,3))
ppx=names(ventr[round(ventr,3)==0.598])
ppy=names(ventrr[round(ventrr,3)==0.598])

# Matochko graph -----

load("matpep")
Gmat=simplify(adjL2Glite(adjL(matpep)))
save(Gmat, file="Gmat")
load("Gmat")
singl=names(V(Gmat)[components(Gmat)$membership>1])
dGmat=degree(Gmat)
ventr=pbsapply(names(V(Gmat)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
boxplot(log10(dGmat)~round(ventr,2), notch=T)
plot(table(ventr))
ventrj=cut(ventr,c(0,unique(ventr)+1e-4), labels=F)
names(ventrj)=names(ventr)

clGmat=cluster_leiden(Gmat, objective_function = "modularity", resolution_parameter = 10)

cntrv=sapply(1:146, function(i){
  g=induced_subgraph(Gmat, V(Gmat)[clGmat$membership==i])
  V(g)[which.max(eigen_centrality(g)$vector)]
})

ventrc=pbsapply(names(cntrv),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
hist(ventr, xlim=c(0,2), col=rgb(0,1,0,0.2))
par(new=T)
hist(ventrc,xlim=c(0,2), col=rgb(1,0,0,0.3))
par(new=F)



# More homogenous graph -----  

load("Gmat")
Gmat=induced.subgraph(Gmat,V(Gmat)[components(Gmat)$membership==1])
v=names(V(Gmat))
dGmat=degree(Gmat)
vrepptrnalize=Vectorize(repptrnalize,"s")
vrptrn=vrepptrnalize(names(V(Gmat)))
Gmat=set_vertex_attr(Gmat, "repptrn", value=vrptrn)
es=ends(Gmat, E(Gmat))
esmnd=rowMeans(apply(es,2,function(j) dGmat[j]))
es=apply(es,2,function(n) vertex_attr(Gmat, "repptrn", V(Gmat)[n]))
eshist=aggregate(es, by=as.data.frame(es), "length")

x=as.numeric(factorize(length(es)))
es=array(es, dim=c(5653,41161,2))

lcschnf=sapply(lcschn, function(x) {
  n=max(as.numeric(unlist(strsplit(x, split=""))))
  (20^-7)*factorial(20)/factorial(20-n)
})

lcschent=sapply(lcschn, function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))

x=rep(0,length(lcschn))
names(x)=lcschn
tmatreptrn=table(vrptrn)
x[names(tmatreptrn)]=tmatreptrn
tmatreptrn=x

p1=tmatreptrn/sum(tmatreptrn)
x=p1[pij1[,1]]*p1[pij1[,2]]
pij1psc=pij1
pij1psc[,3]=x

pij1m=dcast(data=pij1, S1~S2)
rownames(pij1m)=pij1m[,1]
pij1m=pij1m[,-1]
pij1m=as.matrix(pij1m)
pij1m[lower.tri(pij1m)]=t(pij1m)[lower.tri(t(pij1m))]

pij1pscm=dcast(data=pij1psc, S1~S2)
rownames(pij1pscm)=pij1pscm[,1]
pij1pscm=pij1pscm[,-1]
pij1pscm=as.matrix(pij1pscm)
pij1pscm[lower.tri(pij1pscm)]=t(pij1pscm)[lower.tri(t(pij1pscm))]
Bym=pij1m*pij1pscm
Bym=Bym/(sum(Bym[lower.tri(Bym)])+sum(diag(Bym)))

save(pij1pscm,file="pij1pscm")



ij=combn(7,5)
cl=makeCluster(10)
clusterExport(cl, c("ij","pwm"))
clusterEvalQ(cl,require(qualV))
esfold=pbapply(es,1,function(x){
  x=strsplit(x,split="")
  z=LCS(x[[1]],x[[2]])$LCS
  p=sum(apply(ij,2,function(j){10^sum(pwm[cbind(sapply(z, function(a) which(rownames(pwm)==a)),j)])}))
  p*p/4.306640625e-11
}, cl=cl)
closeAllConnections()
esprb=array(esprb, dim=c(5653,41161,2))
save(esprb, file="esprb")

AL=as_adj_list(Gmat)
AL=pblapply(AL, names)
n=names(AL)
AL=pbmapply(c,n,AL)
cl=makeCluster(10)
clusterExport(cl, "vrepptrn")
AL=pblapply(AL,function(l) vrepptrn(names(l)), cl=cl)
closeAllConnections()
names(AL)=vrepptrn(names(AL))
save(AL,file="AL_gmat")

cl=makeCluster(14)
clusterExport(cl,"Bym")
wGmat=pbapply(es, 2, function(j) {
  Bym[j] 
}, cl=cl)
closeAllConnections()
save(wGmat,file="wGmat")
load("wGmat")

wGmat=wGmat*(10^esprb)
#wGmat1=wGmat1/sum(wGmat1)
cl=makeCluster(14)
clusterExport(cl,"pij1m")
wpij=pbapply(es, 2, function(j) {
  pij1m[j] 
}, cl=cl)
closeAllConnections()

rm(es)
gc()  

wGmat=c(wGmat)  
wpij=c(wpij)


# USE es just out of the ends function
es=cbind(es,wGmat)
colnames(es)[1:2]=c("P1","P2")
es=as.data.frame(es)
es$wGmat=as.numeric(es$wGmat)
es$wGmat=log10(es$wGmat)
es1=es[,c(1,3)]
es2=es[,2:3]
es=rbind(as.matrix(es1),as.matrix(es2))
es=as.data.frame(es)
es[,2]=as.numeric(es[,2])
save(es, file="esw_long")

load("esw_long")
Bgx=aggregate(es[,2], by=list(es[,1]), "mean")
x=Bgx$x
names(x)=Bgx$Group.1
Bgx=x
j=intersect(names(Bgx),names(dGmat))


vrj=vrptrn[j]
lmBxjD=lm(log10(dGmat[j])~Bgx[j])
lmBxjDe=lm(log10(dGmat[j][vrj=="1112345" ])~Bgx[j][vrj=="1112345"])
summary(lmBxjDe)

BxDcf=lmBxjD$coefficients
BxDcf=c(10^BxDcf[1],BxDcf[2])

plot(Bgx, ventr[names(Bgx)])
plot(Bgx[j], dGmat[j], pch=16, cex=0.1, col=rgb(0,0,0,0.7), log="y", xlab="P", ylab="D")
Y=log10(dGmat)
x=cut(Bgx, c(-16,-12,-11,-9.5,-8,-6), labels=F)
X=x[names(dGmat)]
boxplot(Y~X, notch=T, varwidth=T)


tvrj=table(vrj)
vrji=vrj[vrj %in% names(tvrj[tvrj>3])]
vrji=vrji[names(vrji) %in% j]
rptrnvrji=sapply(unique(vrji), function(l) entropy(table(unlist(strsplit(l,split=""))), method="ML"))
colrptr=as.numeric(as.factor(rptrnvrji))
vrjprms=t(pbsapply(unique(vrji), function(rp){
  lmBxjDe=lm(log10(dGmat[names(vrji)[vrji==rp]])~Bgx[names(vrji)[vrji==rp]])
  n=sum(vrji==rp)
  # print(rp)
  # print(n)
  # print(summary(lmBxjDe))
  #plot(Bgx[j][vrj==rp], dGmat[j][vrj==rp], log="y", pch=16, cex=0.3, col=rgb(0,0,0,0.7), xlim=range(Bgx), ylim=range(dGmat[j]), main=summary(lmBxjDe)$coefficients[2,1])
  #abline(lmBxjDe)
  return(c(n,summary(lmBxjDe)$coefficients[,1]))
}))

plot(vrjprms[,1:2], log="x", col=cpl1(14)[colrptr], pch=16, cex=colrptr/8,xlim=range(vrjprms[,1]), ylim=range(vrjprms[,2]))
plot(Bgx[j], dGmat[j], pch=16, cex=0.1, col=rgb(0,0,0,0.7), log="y", xlab="P", ylab="D")
plot(Bijx[j], dGmat[j], pch=16, cex=0.1, col=rgb(0,0,0,0.7), log="y", xlab="P", ylab="D")

zz=names(rptrnvrji)[colrptr==4]

for (n in zz){
  plot(Bgx[names(vrji)[vrji==n]], dGmat[names(vrji)[vrji==n]], pch=16, cex=0.7, log="y", xlab="P", ylab="D", xlim=range(Bgx), ylim=range(dGmat), main=n)
  lmBxjDe=lm(log10(dGmat[names(vrji)[vrji==n]])~Bgx[names(vrji)[vrji==n]])
  abline(lmBxjDe)
  legend("bottomright",legend = summary(lmBxjDe)$coefficients[,1])
}


uvntr=sort(unique(ventr))
vrprms=t(pbsapply(uvntr, function(p){
  #if (table(vrj)[p]<2) return(c(0,0,0))
  ii=ventr[j]==p
  lmj=summary(lm(log10(dGmat[j][ii])~Bgx[j][ii]))
  return(c(length(dGmat[j][ii]),lmj$coefficients[,1]))
}))
closeAllConnections()

gcjj=vrj %in% names(table(vrj)[table(vrj)==20])
plot(Bgx[j][ventr[j]==unique(ventr)[8]], dGmat[j][ventr[j]==unique(ventr)[8]], log="y", pch=16, cex=0.3, col=rgb(0,0,0,0.7), xlim=range(Bgx), ylim=range(dGmat[j]), main=i)


W=Bgx[v]
W=(10^-W)/1e10

Gmat=set_edge_attr(Gmat, "weight", value=W)
save(Gmat,file="Gmatw")

jj=c(-11,seq(min(log10(wGmat)), max(log10(wGmat)),length.out=10000),-1)
jj=cut(log10(wGmat), jj, labels=F)
esmndagg=aggregate(cbind(esmd=log10(esmnd),wGm=log10(wGmat)), by=list(jj), "mean")
plot(esmndagg[,c(3,2)], pch=16, cex=1, col=rgb(0,0,0,0.7))
Y=lm(data=esmndagg, esmd~wGm)
abline(Y)

boxplot(log10(strength(Gmat))~round(ventr,2), notch=T)
boxplot(log10(dGmat)~round(ventr,2), notch=T)

# all v mx ----

Bymrptrn=pbsapply(rownames(Bym),function(p) Bym[p,vrptrn])
save(Bymrptrn,file="Bymrptrn")
Pijrptrn=pbsapply(rownames(pij1m),function(p) pij1m[p,vrptrn])
save(Pijrptrn,file="Pijrptrn")

load("AL_gmat")
load("Gmat")
load("es")
load("Bymrptrn")
esptrn=pbapply(es,2,function(p) vrptrn[p])
rm(es)
gc()

BymrptrnS=colSums(Bymrptrn)
PijrptrnS=colSums(Pijrptrn)
BymrptrnS=array(BymrptrnS%x%BymrptrnS, dim=c(length(BymrptrnS),length(BymrptrnS)), dimnames=dimnames(Bym))
#PijrptrnS=array(PijrptrnS%x%PijrptrnS, dim=c(length(PijrptrnS),length(PijrptrnS)), dimnames=dimnames(pij1m0))
#BymrptrnS=pbapply(Bymrptrn, 2, function(s1){
#   apply(Bymrptrn,2,function(s2){
#     sum(s1)*sum(s2)
#   })
# })

#vfield=pij1m[esptrn]/(PijrptrnS[esptrn[,1]]+PijrptrnS[esptrn[,2]])

cl=makeCluster(14)
clusterExport(cl, c("Bym"))
adje=pbsapply(AL, function(l){
  p=l[1]
  pp=l[-1]
  #p1=sum(Bym[p,vrptrn])
  p2=sum(Bym[p,pp])
  return(p2)
},cl=cl)
closeAllConnections()

BymrptrnSv=BymrptrnS[vrptrn]
PijrptrnSv=PijrptrnS[vrptrn]

names(BymrptrnSv)=NULL

vfield=adje[,1]/(BymrptrnSv)
names(vfield)=v

Y=(2/(vfield[es[,1]]+vfield[es[,2]]))

Gmat=set_edge_attr(Gmat, "weight", value=Y)
boxplot(log10(strength(Gmat))~round(ventr,2), notch=T)

lapply(sapply(sort(unique(ventr)), function(x) names(ventr)[which(ventr==x)[1]]), function(y) table(unlist(strsplit(y,split=""))))


clsGmat=cluster_leiden(Gmat, objective_function = "modularity")
aggregate(ventr, by=list(clsGmat$membership), "mean")
j=which(table(clsGmat$membership)>100)

for (i in j) {
  g=induced_subgraph(Gmat, v[clsGmat$membership==i])
  ecw=eigen_centrality(g)
  ec0=eigen_centrality(g, weights = rep(1,ecount(g)))
  print(c(ventr[names(which.max(ecw$vector))],ventr[names(which.max(ec0$vector))]))
}

pij1noz=pij1[pij1[,3]>0,]

pjsmpl1=sample(lcschn[1:100],500, replace=T)
pjsmpl2=sample(lcschn[787:877],500, replace=T)
pjjs=rbind(pjsmpl1,pjsmpl2)
pjjs=t(unique(t(pjjs)))
lsc0=array(seq_along(lcschn), dimnames=list(lcschn))
pjjs=t(apply(pjjs,2,function(l) {
  l=l[order(lsc0[l])]
  unlist(pij1[pij1[,1]==l[1]&pij1[,2]==l[2],])
}))
pjjs=data.frame(S1=pjjs[,1],S2=pjjs[,2],P=as.numeric(pjjs[,3]))
pjjs=pjjs[pjjs$P>0,]

cl=makeCluster(10)
clusterExport(cl, c("aa"))
clusterEvalQ(cl, require(stringdist))
clusterEvalQ(cl, require(stringi))
probBs=pbapply(pjjs,1,function(ii){
    s10=as.numeric(unlist(strsplit(ii[[1]], split="")))
    s20=as.numeric(unlist(strsplit(ii[[2]], split="")))
    pr=as.numeric(ii[[3]])
    n1=length(unique(s10))
    n2=length(unique(s20))
    N=round(50/pr)
    if (N>5e7) N=5e7
    p1=sapply(1:N,function(z1) {
      z1=sample(aa,n1)
      stri_join(z1[s10], collapse="")
      })
    p2=sapply(1:N,function(z1) {
      z1=sample(aa,n2)
      stri_join(z1[s20], collapse="")
    })
    m=stringdist(p1,p2, nthread=2, method="lcs")
    p1=NULL;p2=NULL;rm(p1,p2);gc()
    s=sum(m<5)
    m=NULL;rm(m);gc()
    
    return(s/N)
}, cl=cl)

closeAllConnections()
plot(probBs,as.numeric(pjjs[,3]), log="xy", pch=16, cex=.5, col=rgb(0,0,0,0.3))

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
# "
#ss=sample(names(V(Gmat)),250000)
#ss=sapply(1:250000, function(i) paste(sample(aa,7,replace = T),collapse="", sep=""))
#ss=sample(rndscflat)
core=sample(V(Gmat),1)
egg=names(ego(Gmat,2,core)[[1]])
g=induced.subgraph(Gmat,egg)

# al=adjL(ss)
# g=adjL2Glite(al)
# g=induced_subgraph(g, V(g)[components(g)$membership==1])
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

pwm=consensusMatrix(AAStringSet(vg))/length(vg)
x=apply(pwm+3e-4,2,function(l) sum(l*log2(l)))
x=x/20
y=0.05*log2(0.05)
c0=2^y
exp(lambertWn(log(c0)))
xx=exp(lambertWn(log(2^x)))
plike=mean(1/xx)

cl=makeCluster(16)
clusterExport(cl, c("lcschn", "thProb1", "Lcsch","plike"), envir = environment())
Bymg=pbsapply(lcschn, function(sc1) {
        sapply(lcschn, function(sc2) {
          thProb1(sc1,sc2, A=20)
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

vagg=aggregate(esp[,1], by=list(esp[,1],esp[,2]), "length")

p4=aggregate(log10(dg), by=list(vrp), "mean")
x=p4$x
names(x)=p4$Group.1
p4=10^x
vaggF=cbind(D1=p4[vagg$Group.1],D2=p4[vagg$Group.2], N=vagg$x)
plot()
p4=sqrt(p4[esp[,1]]*p4[esp[,2]])

p4t=Bymg*p1
p4t=rowSums(p4t)/(tmpn+0.5)
p4t=sqrt(p4t%x%p4t)
p4t=array(p4t, dim = dim(Bymg), dimnames=dimnames(Bymg))

#p4t=p4t/(2*ecount(g))
#p4t=log10((p4t[esp[,1]]*p4t[esp[,2]])/ecount(g))

p3=Bymg*p1_
#p3=rowSums(p3)/lcschlt_n
p3=p3*vcount(g)/(plike^7)

p5=sqrt(p3[esp[,1]]*p3[esp[,2]])

Yg=1/p4  #1/(p3[esp[,1:2]])          /(p1[esp[,1:2]]*)   (1/p2[esp[,1:2]])  espf[esp[,1:2]]/(sum(espf))  1/(60*ecount(g)*MX[esp[,1:2]]
#Yg[is.infinite(Yg)]=0

g=set_edge_attr(g, "weight", value=Yg)

boxplot(log10(strength(g))~round(gentr,2), notch=T)
boxplot(log10(dg)~round(gentr,2), notch=T)
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

# Small Graph Edgewise 

n=length(vg)
  #w=1/(12*p1)+1/(ecount(g)*pij1m)
  w=(pij1m+p2)/(pij1m*p1)
  w[is.infinite(w)]=0
  #w=w/sum(w)

w=w[esptrn]  

g=set_edge_attr(g, "weight", value=w)
boxplot(log10(strength(g))~round(gentr,2), notch=T)
boxplot(log10(dg)~round(gentr,2), notch=T)

gn0=gentr
gn0[gn0<0.9]=0
optx=JDEoptim(1,100,thiggraw)

# Igome graph -----

load("G")
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



# Scrapbook -----

Bymgpn=pbsapply(rownames(Bymg),function(p) Bymg[p,vrp])
Pijprn=pbsapply(rownames(pij1m),function(p) pij1m[p,vrp])
esptrn=pbapply(esg,2,function(p) vrp[p])
BymgpnS=colSums(Bymgpn)
PijprnS=colSums(Pijprn)

ALg=as_adj_list(g)
ALg=lapply(ALg, names)
n=names(ALg)
ALg=pbmapply(c,n,ALg)
cl=makeCluster(10)
clusterExport(cl, "vrepptrn")
ALg=pblapply(ALg,function(l) vrepptrn(l), cl=cl)
closeAllConnections()
names(ALg)=vrepptrn(names(ALg))

cl=makeCluster(14)
clusterExport(cl, c("Bymg"))
adge=pbsapply(ALg, function(l){
  p=l[1]
  pp=l[-1]
  #p1=sum(Bym[p,vrptrn])
  p2=sum(Bymg[p,pp])
  return(p2)
},cl=cl)
closeAllConnections()
gfield=adge[1]/(BymgpnS)


