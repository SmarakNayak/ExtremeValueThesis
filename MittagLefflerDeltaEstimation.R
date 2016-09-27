# From Dexter Cahoy, http://www2.latech.edu/~dcahoy/
#######################
#Random number generation and  estimation for the Mittag-Leffler distribution  
#and fractional Poisson process
#SOURCE: (1) Parameter estimation for fractional Poisson processes (with V Uchaikin and W Woyczynski). 
#Journal of Statistical Planning and Inference, 140(11), 3106-3120, Nov 2010.
#(2)Estimation of Mittag-Leffler parameters. Communications in Statistics-Simulation and Computation, 42(2), 303-315, 2013
#Updated last 02/03/2016.
#Email me at dcahoy@latech.edu if you have questions, suggestions, etc.
########################

EULER.C = 0.57721566490153286;

# n - number of data needed
# nu, mu - the parameters of Mittag-Leffler dist.
mittag.leffler = function(n, nu, mu){
  U.1 = runif(n)
  U.2 = runif(n)
  U.3 = runif(n)
  T = ( (abs(log(U.1))/mu)^(1/nu) )*( ( sin(nu*pi*U.2)*sin((1-nu)*pi*U.2)^(1/nu-1) )/ (  ( sin(pi*U.2)^(1/nu) )*abs(log(U.3))^(1/nu-1))  )
  return(T)
}


# T = data, 1-a= 1- alpha =confidence level/coefficient, F.l = cumulative prob of threshold level
ml.par.est.delta = function (T, a, F.l) {
  log.T = log(T)
  m = mean(log.T)
  s.2 = var(log.T)
  nu = pi/sqrt(3*(s.2 + pi^2/6))
  delta = exp(m + EULER.C)
  n=length(T)
  
  se.nu=sqrt( (nu^2)*(32-20*nu^2-nu^4)/(40*n)  )
  zcv=qnorm(1-a/2,0,1)  
  l.nu= nu -zcv*se.nu
  u.nu =  nu + zcv*se.nu
  
  se.delta = sqrt(((pi^2*delta^2)/(6*n))*((2/nu^2) - 1)) #delta
  l.delta= delta -zcv*se.delta
  u.delta =  delta + zcv*se.delta
  
  b=delta*(-log(F.l))^(1/nu)
  var.b=(b^2*s.2 + (b^2/nu^2)*log(-log(F.l)) * (6*sqrt(2)*pi/(6*s.2+pi^2)^(3/2) + se.nu^2/(nu^2)*log(-log(F.l))))/n
  se.b=sqrt(var.b)
  l.b=b-zcv*se.b
  u.b=b+zcv*se.b
  
  
  return(list(nu = nu, delta = delta, b=b, CInu=c(l.nu, u.nu), CIdelta=c(l.delta, u.delta), CIb=c(l.b,u.b)) )      
}

#Generate Mittag-Leffler distributed data 
#or inter-event or  inter-jump times of a fractional Poisson process
# dat=mittag.leffler(n=500,nu=0.2, mu=2)

#Point and interval estimates of  nu and mu. 
# ml.par.est(dat, a=0.05)

