checkAPLmtfX_sh=function(s){
  require(Biostrings)
  aa=AA_STANDARD
  l=length(m)
  f=c("A", "V", "I", "L", "M", "F", "W", "C", "P", "G")
  z=c("Y", "T", "S", "H", "K", "R", "E", "D", "Q", "N")
  ff="F"
  ai=c("f","f","z","z","ff") #c("f","f","f","f","z","z","ff","aa","f")
  r0=c("f","z","ff")
  x=unlist(strsplit(s, split=""))
  all(sapply(1:5,\(i){
    x[i] %in% get(ai[i])
  }))
}