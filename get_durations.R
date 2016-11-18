# extracts (preceding) durations of the m largest exceedances
get_durations <- function(TT,idxJ,m){
  m=ceiling(m)
  thinT <- sort(TT[idxJ[1:m]])
  diff(thinT)
}
