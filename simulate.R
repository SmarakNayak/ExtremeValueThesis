require(stabledist)
n = 1000
beta = 0.9
TT = rstable(n,beta, 1, pm=1)
require(fExtremes)
JJ = rgev(n, xi = 0.3, mu = 0, beta = 1)

# This is what the original data could look like
par(mfrow=c(2,2))
plot(TT,JJ,type='h', xlim = c(0,max(TT)), ylim=c(0,max(JJ)))

# Now threshold the data
# Take the m largest observations
m=500
ii = order(JJ, decreasing = T)[1:m]
plot(TT[ii],JJ[ii], type='h',xlim = c(0,max(TT)), ylim=c(0,max(JJ)))

## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs
amountOfObservations = 1
hillEstimate=matrix(0,nrow=m)
while(amountOfObservations<m){
  orderedObservations=order(TT,decreasing = T)[1:amountOfObservations]
  hillEstimate[amountOfObservations]=(1/amountOfObservations)*sum(log(TT[orderedObservations]))-log(TT[orderedObservations][amountOfObservations])
  amountOfObservations=amountOfObservations+1
}
par(mfrow=c(2,2))
plot(hillEstimate,type="l")
grid()

##or use inbuilt

hillPlot(TT,start=10)
gevFit(JJ,block=10,type = "pwm")