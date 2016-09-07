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

#Generate Mittag-Leffler distributed data 
#or inter-event or  inter-jump times of a fractional Poisson process
# dat=mittag.leffler(n=500,nu=0.2, mu=2)

#Point and interval estimates of  nu and mu. 
# ml.par.est(dat, a=0.05)

  
