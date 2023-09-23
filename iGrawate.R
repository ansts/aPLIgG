# Calculates sequence pattern of repetition specific edge weights based on a
# particular background (control) IgOme graph

iGrawate=function(g,lcschn=NULL){
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
  require(gtools)
  
  aa=AA_STANDARD
  rm(Gmat)
# convert peptide seqs in rep. patterns
  vrepptrnalize=Vectorize(repptrnalize,"s") 
  
# All patterns of repetitions can be calculated by -----
  print("Rep. patterns calculation........")
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
  print(sum(tchoose))
  vgrep=Vectorize("grep",vectorize.args = "pattern")
  
# lcschn - all patterns of repetitions ----
# lcschlt - number of letters of the patterns
# lcschnf - theoretical distribution of the frequencies in uniformly distrbuted
# aminoacids by positions
  
  if (is.null(lcschn)) {
      sch0=c()
      print("Exctracting all possible patterns of repetition by simulation......")
      repeat {
        px=t(sapply(1:1000000, function(i) sample(aa,7,replace = T)))
        cl=makeCluster(15)
        clusterExport(cl,c("vgrep"))
        clusterEvalQ(cl, require(reshape2))
        sch=pbapply(px,1,function(p) {
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
        sch0=c(sch0,sch)
        sch0=(unique(sch0))
        l=length(sch0)
        print(l)
        
        if (l==877) break
      }
      lcschn=sch0
  }
  lcschlt=sapply(lcschn, function(s) length(table(unlist(strsplit(s,split="")))))
  lcschn=lcschn[order(lcschlt)]
  lcschnf=sapply(lcschn, function(x) {
    n=max(as.numeric(unlist(strsplit(x, split=""))))
    (20^-7)*factorial(20)/factorial(20-n)
  })
  
  # Basic graph parameters, reppatternalizing edges -----------
  print("Basic graph parameters extraction........")
  vg=names(V(g))
  ltrcnt=pbsapply(vg,function(x) length(table(unlist(strsplit(x, split="")))))
  dg=degree(g)
  vrp=vrepptrnalize(vg)
  esg=ends(g, E(g))
  esp=cbind(vrp[esg[,1]],vrp[esg[,2]])
  rm(esg)
  el=esp
  cl=makeCluster(16)
  el=t(pbapply(el,1,sort,cl=cl))
  rm(esp)
  closeAllConnections()
  gentr=pbsapply(names(V(g)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))

  # PSPM of sequences in the graph ----
  pwm=consensusMatrix(AAStringSet(vg))/length(vg)
  
  # Distribution of rep. patterns based on the PSPM
  print("Calculation of rep. pattern distribution as a function of PSPM.......")
  lcschnfmat=lcschPwmp(L=lcschn)
  save(lcschnfmat,file="lcschnfmat")
  
  plot(lcschnf, lcschnfmat, xlim=c(1e-8,1), ylim=c(1e-8,1), log="xy", 
       xlab="Probability based on uniform distribution of aa", 
       ylab="Probability based on the PSPM")
  lines(c(1e-8,1),c(1e-8,1))
  
  print("Weight calculation......")
  
  # theoretical number of vertices of each rep. pattern
  thvc=lcschnfmat*vcount(g)/100  

  # Distribution of patterns in the vertices
  x=rep(0,length(lcschn))
  names(x)=lcschn
  tvrp=table(vrp)
  x[names(tvrp)]=tvrp
  tvrp=x
  
  # with pseudocounts
  
  tvrpc=vcount(g)*(tvrp+thvc)/sum(tvrp+thvc)
  
  # Matrix of sums of vertices by pattern - the vertices in all subgraphs
  # connected (egonetwork) to all the vertices of each pair of repetitions
  nsqpp=pbsapply(seq_along(tvrpc), function(i1){
    sapply(seq_along(tvrpc), function(i2){
      tvrpc[i1]+tvrpc[i2]
    })
  })
  diag(nsqpp)=diag(nsqpp)/2
  dimnames(nsqpp)=dimnames(eaggmm) 
  
  # Empirical distribution of number of edges by rep.pattern pairs ----
  eagg=aggregate(el[,1], by=list(el[,1],el[,2]), "length")
  ij=as.matrix(eagg[,1:2])
  y=as.numeric(eagg$x)
  eaggm=array(0, dim=c(length(lcschn),length(lcschn)), dimnames = list(lcschn,lcschn))
  eaggm[ij]=y
  eaggm=eaggm+t(eaggm)
  diag(eaggm)=diag(eaggm)/2
  
  # with pseudocounts
  eaggmpc=sum(eaggm)*(eaggm+min(eaggm[eaggm>0])/100)/sum(eaggm+min(eaggm[eaggm>0])/100)
                                                                           
  thdg=rowSums(eaggmpc)/tvrpc
  
  # Theoretic number of edges between rep. patterns based on the graph degrees
  # p4 - potential number of edges within each pair of rep. pattern
  p4=aggregate((dg), by=list(vrp), "mean")
  x=array(0,dim=length(lcschn), dimnames = list(lcschn))
  x[p4$Group.1]=p4$x
  p4=x
  
  #with pseudocaounts
  p4pc=sum(p4)*(p4+thdg/100)/sum(p4+thdg/100)
    
  p4pc=pbsapply(names(p4pc), function(n1){
    sapply(names(p4pc), function(n2){
      if (n1!=n2) p4pc[n1]*p4pc[n2] else p4pc[n1]*(p4pc[n1]-1)/2
    })
  })
  rownames(p4pc)=colnames(p4pc)
  # Expanding the matrix to include the 0 freq rep patterns
  x=array(0,dim=dim(eaggm), dimnames=dimnames(eaggm))
  x[rownames(p4pc),colnames(p4pc)]=p4pc
  p4pc=x
  rm(x) 

  # Matrix of sums of edges by pattern - the edges in all subgraphs connected
  # (egonetwork) to all the vertices of each pair of repetitions
  rv=rowSums(eaggmpc)
  eaggmm=pbsapply(rv, function(n1){
    sapply(rv, function(n2){
      n1+n2
    })
  })
  eaggmm=eaggmm-eaggm
  
  # Calculation of the weights based on the characteristic edge densities -----
  gmatwtmx=(p4pc)*nsqpp/(2*eaggmm)
  save(gmatwtmx, file="gmatwtmx")
  
  # x/(2*eaggmm) is the probability of pair of peptides to have an edge between
  # them (based on the empirical distribution from this graph) as a function of
  # their rep. patterns. Multiplied by the number of vertices in this pair of
  # rep. patterns nsqpp, it gives the expected degree of a vertex in this
  # subgrpah (slice). Its reciprocal value is the weight which normalizes the
  # variance of this rep. pattern defined characteristic degree.
  
  wt=1/gmatwtmx[el]
  
  g=set_edge_attr(g, "weight", value=wt)
  
  boxplot(log10(strength(g))~round(gentr,2), notch=T, 
          xlab="Entropy [bits]", ylab="Log Strength", main = "Corrected")
  boxplot(log10(dg)~round(gentr,2), notch=T, 
          xlab="Entropy [bits]", ylab="Log Strength", main = "Raw")
  rm(g,el,wt)
  closeAllConnections()
  
  load("G")
  dG=degree(G)
  ventr=pbsapply(names(V(G)),function(x) entropy(table(unlist(strsplit(x, split=""))), method="ML"))
  esG=ends(G, E(G))
  vG=names(V(G))
  vrp=vrepptrnalize(vG)
  esp=cbind(vrp[esG[,1]],vrp[esG[,2]])
  rm(esg)
  el=esp
  cl=makeCluster(16)
  el=t(pbapply(el,1,sort,cl=cl))
  rm(esp)
  closeAllConnections()
  wt=1/gmatwtmx[el]
  Gw=set_edge_attr(G, "weight", value=wt)
  save("Gw", file="Gw")
  return(NULL)
}  