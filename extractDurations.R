# extracts (preceding) durations of the m largest exceedances
Tell <- function(TT,idxJ,m){
  m=ceiling(m)
  thinT <- sort(TT[idxJ[1:m]])
  diff(c(0,thinT))
}
