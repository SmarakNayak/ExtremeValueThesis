require(stabledist)
n=1000
betaVec=seq(0.1,0.9,0.1)
y=seq(-10,20,1)
par(mfrow=c(4,3))
for (beta in betaVec)
{
  for (x in c(0:10))
  {
    FirstDens=dexp(exp((x-y)*beta),1)*beta*exp((x-y)*beta)
    SecondDens=dstable(exp(y),beta,1,pm=1)*exp(y)
    fxy_integrand=FirstDens*SecondDens
    plot(y,fxy_integrand,main=sprintf("x=%d,b=%1.1f",x,beta))
  }
}
beta=0.4
par(mfrow=c(1,1))
x1=c(-10:10)
y1=seq(-20,10,0.1)
f=function(x,y) {dexp(exp((x-y)*beta),1)*beta*exp((x-y)*beta)*dstable(exp(y),beta,1,pm=1)*exp(y)}
z<-outer(x1,y1,f)
par(mfrow=c(1,1))
persp(x1,y1,z,theta = 30, phi = 30, expand = 0.5, col = "lightblue",ticktype = "detailed")