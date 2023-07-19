thprobEK = function(X1,X2,L=Lcsch7,A=20) {
  
  
    p1=X1
    n1=L[[p1]]
    n1_=n1$n
    l1=n1$Group.1
    
      p2=X2
      n2=L[[p2]]
      n2_=n2$n
      l2=n2$Group.1
      
      l=intersect(l1,l2)
      if (length(l)==0) return(0)
      
      pp=sapply(l,function(li) {
        n=unlist(strsplit(li,split=""))
        n0=sum(n==0)
        n1=sum(unique(n)!=0)
        n=n0+n1
        factorial(A-n)/factorial(A)
      })
      # order patterns according to their probabilities
      l = l[order(pp,decreasing = TRUE)]
      i=which(l1 %in% l)[order(pp,decreasing = TRUE)]
      j=which(l2 %in% l)[order(pp,decreasing = TRUE)]
      nprod = (n1_[i]*n2_[j])
      pp = pp[order(pp,decreasing = TRUE)]
      
      allcomb1 = c()
      allcomb2 = c()
      
      # start from the pattern with highest probability(and lowest diversity),
      # in order to exclude probabilities which are a subset of the bigger ones:
      for (ii in 1:length(l)){
        
        if(length(l) == 1){
          break
        } else if(l[1] == "PPPPP"){
          nprod[2:length(nprod)] = 0
          break
        }
        else{
          # write all pairs of potentially corresponding aas between X1 and X2:
          temp = sapply(n1$V3[[i[ii]]], function(x){
            sapply(n2$V3[[j[ii]]], function(y){ 
              x = unlist(strsplit(x,split=""))
              y = paste0(unlist(strsplit(y,split=""))[order(x)], collapse = "")
              x = paste0(x[order(x)], collapse = "")
              c(x,y)
            })
          })
          # separate variable for each string, temp1[i] corresponds to temp2[i]
          temp1 = sapply(1:(length(temp)/2), function(z) temp[2*z-1]) 
          temp2 = sapply(1:(length(temp)/2), function(z) temp[2*z])
          
          if(ii>1) {
            # search which correspondences from the previous patterns occur in the new one:
            m = unlist(sapply(1:length(allcomb1),function(k){
              g1 = regexpr(allcomb1[k],temp1)
              g2 = regexpr(allcomb2[k],temp2)
              m = which(g1 != -1 & g1==g2 & attr(g1,"match.length") == attr(g2,"match.length"))
              return(m)
            }))
            
            
            if(length(m)>0)
            { # remove the duplicated correspondences:
              temp1 = temp1[-m]
              temp2 = temp2[-m]
              nprod[ii] = length(temp1)
            }
          }
          # convert to regular expressions, in order to search in the next patterns:
          temp1 = sapply(temp1, function(x) paste0(unlist(strsplit(x,"")),collapse = "\\d*"))
          temp2 = sapply(temp2, function(x) paste0(unlist(strsplit(x,"")),collapse = "\\d*"))
          allcomb1 = c(allcomb1,temp1)
          allcomb2 = c(allcomb2,temp2)
        }
      }
      sum(pp*nprod)                   
}


