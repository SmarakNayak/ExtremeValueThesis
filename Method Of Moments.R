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

EULER.C = 0.57721566490153286;

# T = data, 1-a= 1- alpha =confidence level/coefficient
ml.par.est = function (T, a) {
  log.T = log(T)
  m = mean(log.T)
  s.2 = var(log.T)
  nu = pi/sqrt(3*(s.2 + pi^2/6))
  mu = exp(-nu*(m + EULER.C))
  n=length(T)
  se.nu=sqrt( (nu^2)*(32-20*nu^2-nu^4)/(40*n)  )
  zcv=qnorm(1-a/2,0,1)  
  l.nu= nu -zcv*se.nu
  u.nu =  nu + zcv*se.nu
  
  se.mu = sqrt( ((mu^2)/(120*(pi^2)*n) )*( 20*(pi^4)*(2-nu^2)-3*(pi^2)*(nu^4+20*nu^2-32)*(log(mu))^2 - 720*(nu^3)*log(mu)*1.20206 )  ) #MU
  l.mu= mu -zcv*se.mu
  u.mu =  mu + zcv*se.mu
  return(list(nu = nu, mu = mu, CInu=c(l.nu, u.nu), CImu=c(l.mu, u.mu)) )      
}

ml.par.est(Ta,0.05)