#GEV ESTIMATE
#install.packages("fExtremes")
require(fExtremes)
n=1000
jumps = rgev(n, xi = 0.3, mu = 0, beta = 1)

# ##BlockSizing - not necessary, can just change the block size parameter in gevFit
# b=2
# jumpMatrix=matrix(jumps,nrow=b,ncol=floor(n/b))
# jumps=apply(jumpMatrix,2,max)


# Now threshold the data
# Take the m largest observations
m = 50
indexesOfLargestJumps = order(jumps, decreasing = T)[1:m]

a=jumps[indexesOfLargestJumps[m]]

model=gevFit(jumps,block=1)
xi=model@fit$par.ests[1]
mu=model@fit$par.ests[2]
beta=model@fit$par.ests[3]

##Works out to be the proportion m/n
Fa=pgev(a,xi,mu,beta)