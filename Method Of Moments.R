require(stabledist)
n = 100000
beta = 0.7
sumWaitingTimes = cumsum(rstable(n,beta, 1, pm=1))
require(fExtremes)
jumpSizes = rgev(n, xi = 0.3, mu = 0, beta = 1)

# This is what the original data could look like
par(mfrow=c(1,2))
plot(sumWaitingTimes,jumpSizes,type='h', xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumpSizes)))

# Now threshold the data
# Take the m largest observations
m = 50
indexesOfLargestJumps = order(jumpSizes, decreasing = T)[1:m]
plot(sumWaitingTimes[indexesOfLargestJumps],jumpSizes[indexesOfLargestJumps], type='h',xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumpSizes)))

a=jumpSizes[indexesOfLargestJumps[m]]

##Cheating here to calculate actual F(a) rather than estimating it
Fa=pgev(a,xi = 0.3, mu = 0, beta = 1)
#Calculating actual Ta's
cumulativeDurations=c(0,sort(sumWaitingTimes[indexesOfLargestJumps]))
Ta=diff(cumulativeDurations)
lnTa=log(Ta)

## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs

ec=0.57721566490153286

# #Generate strechable Mittag-Leffler (GML (a,b, lam))
# dat=( rexp(n, rate=lam)^(1/a)  )*stbl(a,n) #MITTAG-LEFFLER data
# dat2= dat^(a/g) #strechable MITTAG-LEFFLER data

####################
#Estimation part
 
mu1=mean(lnTa)  
mu2=var(lnTa)  
mu3=mean( ( lnTa-mu1 )^3 )   
c= ( (mu3^2)^(1/3) )/mu2

require(VGAM)
ah=sqrt( (c*pi^2)/( 3*( (2*zeta(3))^(2/3) + (c*pi^2)/6 ) ) ) #nu  estimate
ah #stability parameter
gh=sqrt(  ( (ah^2)/mu2 )*(pi^2)*( 1 /(3*ah^2) -1/6 )  ) #gamma estimate
gh # gamma shape parameter
lamh= exp( -(mu1*gh  + ec*ah) )  #intensity lambda
lamh #exponential parameter/gamma scale parameter

#or

betaEst=sqrt(2*trigamma(1)/(mu2+trigamma(1)))
betaEst
# f1=function(x,y,beta) {-log(Fa)*dexp(-log(Fa)*exp((x-y)*beta),1)*beta*exp((x-y)*beta)*dstable(exp(y),beta,1,pm=1)*exp(y)}
# f2=function(x,y) {f1(x,y,beta = beta)}
# likelihood=vector("numeric",length = m)
# totalLikelihood=vector("numeric",length = 41)
# betaVec=seq(0.1,0.9,0.02)
# 
# for (j in c(1:41))
# {
#   for (i in c(1:m))
#   {
#     betaHat=betaVec[j]
#     likelihood[i]=as.numeric(integrate(f=function(y){f1(lnTa[i],y,betaHat)},lower=-10,upper=10)[1])
#   }
#   totalLikelihood[j]=prod(likelihood)
# }
# par(mfrow=c(1,1))
# plot(betaVec,totalLikelihood,main=bquote("MLE estimate for data generated with "*beta == .(beta)),xlab=expression(hat(beta)))
