thProb1 = function(X1,X2,co=NULL, L=Lcsch, A=20) {
  require(reshape2)
    p1=X1
    n1=L[[p1]]
    n1_=n1$n
    n10=n1$N
    l1=n1$Reptrn0
    p2=X2
    n2=L[[p2]]
    n2_=n2$n
    n20=n2$N
    l2=n2$Reptrn0
    l=intersect(l1,l2)
    if (length(l)==0) return(0)
    i=which(n1$Reptrn0 %in% l)
    j=which(n2$Reptrn0 %in% l)
    # Reshaping the data
    ll1=data.frame(cbind(Rep0=rep(n1[i,1],lengths(n1[i,2])),Rep1=unlist(n1[i,2])),Source=rep(1,sum(lengths(n1[i,2]))))
    ll2=data.frame(cbind(Rep0=rep(n2[j,1],lengths(n2[j,2])),Rep1=unlist(n2[j,2])),Source=rep(2,sum(lengths(n2[j,2]))))
    ovr=rbind(ll1,ll2)
    
    # No of unique letters
    al=sapply(unique(as.character(ovr$Rep1)),function(x) length(unique(unlist(strsplit(x, split=""))))) 
    
    n_=unique(as.character(ovr$Rep0))
    x=lapply(n_, function(m) {
      gr=expand.grid(ovr$Rep1[ovr$Source==1&ovr$Rep0==m],ovr$Rep1[ovr$Source==2&ovr$Rep0==m])
      n=apply(gr,1,function(l) al[l[1]])
      cbind(gr,n)
    })
    names(x)=n_
    x=suppressMessages({
      melt(x)
    })
    x=x[,-3]
    colnames(x)=c("Sch1","Sch2","al","n_")
    x=x[order(x$al),]
    
    if (nrow(x)>1){
      # Some combination of letters (usually few) generate more than one of the
      # observed common patterns because the 1234567 (Rep1) patterns show that
      # same letters are expected at the same place in both peptides while
      # generating more than one pattern of the Rep0 type. Each such combination
      # of letters has to be represented finally by only one Rep0 pattern and it
      # has to be one of a maximal number of different letter. To this end we
      # construct a matrix of correspondences by position in Rep1 patterns
      # identify and to eliminate redundant Rep0 patterns.
      
        p12=unique(apply(x,1,function(y){
            apply(sapply(strsplit(y[1:2],split=""),unlist),1,paste,collapse ="")
        }))
        p12u=unique(c(p12))
        p12=t(sapply(p12u,function(x){
          apply(p12,2,function(y)( x %in% y)*1) 
        }))
        j=c()                         # Sweep the matrix with each column from large to small no of letters
        for (c0 in ncol(p12):2) {                                                              
          c1=p12[,c0]
          mx=p12[,1:(c0-1), drop=F]
          j12=colSums(mx-mx*c1)==0
          if (any(j12)==T) j=c(j,colnames(p12)[c0]) # Mark the columns including a smaller pattern (of fewer letters)
        }
        j=!(colnames(p12) %in% j)
        p12=p12[,j, drop=F]
        ij=colnames(p12)  # ij are the patterns after filtering out the more complex patterns covered by simpler ones
    } else {
      ij=1
    }
    x=x[ij,]
    n1_=table(unique(x[,c(1,4)])[,2])
    nn=names(n1_)
    n2_=table(unique(x[,c(2,4)])[,2])
    n2_=n2_[nn]
    l  =unique(x[,3:4])
    l=array(l$al,dimnames=list(l$n_))[nn]
    pp=sapply(l,function(li) {
      1/(factorial(A)/factorial(A-li))
    })
    sum(pp*n1_*n2_)                                       
}


