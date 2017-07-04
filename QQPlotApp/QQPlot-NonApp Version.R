require(parsedate)
require(POT)
require(plyr)
require(MittagLeffleR)
source("get_durations.R")
Vipaka <- read.csv("Data/Vipaka.csv")
Vipaka <- Vipaka[-(1:27), ]

parse_iso_8601(Vipaka$time) -> Vipaka$time
n=dim(Vipaka)[1]
plot(Vipaka$time, Vipaka$mag, type='h', ylim=c(0,9))

threshold <- function(a){
  Vipaka[Vipaka$mag >= a, ]
}

ml.est = function (T) {
  log.T = log(T)
  m = mean(log.T)
  s.2 = var(log.T)
  tail = pi/sqrt(3*(s.2 + pi^2/6))
  scale = exp(m + EULER.C)
  
  return(list(tail = tail, scale = scale))
}

##TT Times, JJ magnitudes
TT=as.vector(Vipaka$time)
JJ=Vipaka$mag

# set lowest threshold at this many values
k_max = n
idxJ <- order(JJ, decreasing = T)[1:k_max]
# set highest threshold at this many values
k_min = 5



#set K here
k=750
x=get_durations(TT,idxJ,k)
estimates=ml.est(x)

y=qml(ppoints(k-1),estimates$tail, estimates$scale)
x=sort(x)

plot(x,y)
