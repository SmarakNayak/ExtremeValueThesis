require(stabledist)
require(fExtremes)
n=1000
betaGen = 0.4
sumWaitingTimes = cumsum(rstable(n,betaGen, 1, pm=1))

jumps = rgev(n, xi = 0.3, mu = 0, beta = 1)

# Now threshold the data
# Take the m largest observations
m = 50
indexesOfLargestJumps = order(jumps, decreasing = T)[1:m]
#plot(sumWaitingTimes[indexesOfLargestJumps],jumps[indexesOfLargestJumps], type='h',xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumps)))

a=jumps[indexesOfLargestJumps[m]]

##Cheating here to calculate actual F(a) rather than estimating it
Fa=pgev(a,xi = 0.3, mu = 0, beta = 1)
#Calculating actual Ta's
cumulativeDurations=c(0,sort(sumWaitingTimes[indexesOfLargestJumps]))
Ta=diff(cumulativeDurations)
lnTa=log(Ta)

## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs

##g(x,y;beta)
gxy=function(x,y,beta) {-log(Fa)*dexp(-log(Fa)*exp((x-y)*beta),1)*beta*exp((x-y)*beta)}
##g(x,y;beta) evaluated at x=0
g0y=function(y,beta) {gxy(0,y,beta)}
##f(y)
fy=function(y,beta) {dstable(exp(y),beta,1,pm=1)*exp(y)}

##bounds of integration, dy is the spacing between each quadrature point (1 must be a multiple of dy)
dy=0.2
intLB=-20
intUB=40
integrationMeshForY=seq(intLB,intUB,dy)
noPointsInMesh = length(integrationMeshForY)

##we need to round the empirical log exceedances to the nearest dy
meshedLnTa=round(lnTa/dy)*dy
#Since y is bounded by intLB and intUB, -y is bounded by -intUB and -intLB.
#Thus x-y is bounded by -intUB+min(meshedLnTa) and -intLB+max(meshedLnTa)
extendedIntegrationMeshForXminusY=seq(-intUB+min(meshedLnTa),-intLB+max(meshedLnTa),dy)

##without for loop
beta=0.3
allPossibleGvalues=g0y(extendedIntegrationMeshForXminusY,beta)
xmin=min(meshedLnTa)
G_Entries=allPossibleGvalues
GRowcreator=function(x)
{
  startIndex=(x-xmin)/dy+1
  endIndex=(x-xmin)/dy+noPointsInMesh 
  GRow=G_Entries[startIndex:endIndex]
  return(GRow)
}

GRowcreator(xmin)

G=lapply(meshedLnTa, GRowcreator)

##with for loop
l=41
totalLikelihood=vector("numeric",length = l)
betaHatVec=seq(0.1,0.9,length.out = l)

for (j in c(1:l))
{
  beta=betaHatVec[j]
  F=fy(integrationMeshForY,beta)
  allPossibleGvalues=g0y(extendedIntegrationMeshForXminusY,beta)
  
  xmin=min(meshedLnTa)
  xm=max(meshedLnTa)

  G=matrix(,m,noPointsInMesh)
  
  
  
  itr=1
  for (x in meshedLnTa) 
  {
    startIndex=(x-xmin)/dy+1
    endIndex=(x-xmin)/dy+noPointsInMesh 
    G[itr,]=allPossibleGvalues[startIndex:endIndex]
    itr=itr+1
  }
  likelihoodVector=G%*%F
  logLikelihoodVector=log(likelihoodVector)
  totalLikelihood[j]=sum(logLikelihoodVector)
  #or 
  #totLikelihood=apply(G%*%F, 1, +)
}
plot(betaHatVec,totalLikelihood,main=bquote("MLE estimate for data generated with "*beta == .(betaGen)),xlab=expression(hat(beta)),ylab="Likelihood")
