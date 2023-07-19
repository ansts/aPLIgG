repptrnalize=function(s){
  x=unlist(strsplit(s,split=""))
  y=unique(x)
  j=seq_along(y)
  names(j)=y
  paste(j[x], collapse="")
}