#functions used in QQPlot App

get_durations <- function(TT,idxJ,m){
  m=ceiling(m)
  thinT <- sort(TT[idxJ[1:m]])
  diff(thinT)
}

threshold <- function(a){
  Vipaka[Vipaka$mag >= a, ]
}

EULER.C = 0.57721566490153286;

ml.est = function (T) {
  log.T = log(T)
  m = mean(log.T)
  s.2 = var(log.T)
  tail = pi/sqrt(3*(s.2 + pi^2/6))
  scale = exp(m + EULER.C)
  
  return(list(tail = tail, scale = scale))
}