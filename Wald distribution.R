#################### Working with WALD distribution ################################

# Author: Bram Timmers
# Use: testing WALD distributional forms
# dependencies: already sourced by dependencies.r

###########################################################################

source("dependencies.R")
# Single Accumulation process with random drift variability - adopted from DMC software-package
# Setting sv to 1 as scaling parameter

rPNRWald <- function(n,a,v,sv,t0)
  # random function for single acumulator
{
  drifts <- matrix(rtnorm(n = n , mean = v, sd = sv, lower = 0), 
                   nrow = length(v))
  t0+statmod::rinvgauss(n,mean=a/drifts,shape=a^2)        
}

# i.e. sampling 100 reaction times from SW distribution:
RT <- rPNRWald(1000,2,2,0,1)


#############################################################################

# simulation design  (example)

#############################################################################

# name HMM states
statesNames <- c("strategy 1","strategy 2")

#N
N <- 1000

# Scenario 1 - independent states 
RT <- vector(length=N)

#initiate a markov chain 
HMM.1 <- new("markovchain", transitionMatrix=matrix(c(0.9,0.1,0.1,0.9), nrow=2, byrow=TRUE,
                                                    dimnames = list(statesNames,statesNames)))

#sample a random sequence
HMM.1.seq <- markovchainSequence( N, HMM.1, t0 = sample(HMM.1@states,1))

#make numeric
HMM.1.seq <- as.numeric(factor(HMM.1.seq))

N.state.1 <- sum(HMM.1.seq==1)
N.state.2 <- sum(HMM.1.seq==2)

#initial parameter values of accumulators function(n,a,v,sv,t0)
#p.vector.1 <- c(a=2,v=1,sv=0,t0=1) 
#p.vector.2 <- c(a=2,v=0.1,sv=0,t0=3)

# simulate response times 
state.1 <- rPNRWald(N.state.1,a=2,v=1,sv=0,t0=1)
state.2 <- rPNRWald(N.state.2,a=2,v=1,sv=0,t0=10)

RT[HMM.1.seq==1] <- state.1
RT[HMM.1.seq==2] <- state.2

data.test <- cbind(RT,HMM.1.seq)
data.test <- data.frame(data.test)
names(data.test) <- c("RT","states")
rt <- data.test$RT
############################################################################

source("optimalization obj functions.R")
source("DepmixS4 - Shifted_Wald.R")
source("gamlsss stuff.R") # contains Swald functions
# Optimization checks
lb <- rep(0.1,3)
ub <- c(20,20,min(state.1))
res.1 <- optim(c(1,1,1),obj.seq,y=state.1,w=rep(1,length(state.1)))
res.1b <- optim(c(10,20,1),obj.own,y=state.1)
res.1c <- optim(c(1,20,1),obj,y=state.1)
res.1
# correct <- a=2,v=1,sv=0,t0=1
res.2 <- optim(c(1,1,1),obj.seq,y=state.2,w=rep(1,length(state.2)))
res.2
###############################################################################

# with gamlss

gamlss(RT~1,weights=rep(1,length(data.test$RT)),data=as.data.frame(data.test$RT),family=Swald(),
       control=gamlss.control(c.crit=1e-5,n.cyc=100,trace=FALSE),
       mu.start=2,
       sigma.start=1,
       nu.start=2)

# EM - algorithm - doesn't work now
rt <- data.test$RT
#EM algorithm 
rModels <- list(
  list(
    Swald(rt,pstart=c(2,1,1))),
  list(
    Swald(rt,pstart=c(2,1,1))
  )
)

#rt<-as.data.frame(data.test[,1])
transition <- list()
trstart <- c(0.7,0.3,0.3,0.7)
transition[[1]] <- transInit(~1,nstates=2,data=data.test,pstart=trstart[1:2])#,family=multinomial("identity"))
transition[[2]] <- transInit(~1,nstates=2,data=data.test,pstart=trstart[3:4])#,family=multinomial("identity"))#,family=multinomial("identity"))

instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,data=data.frame(rep(1,1)),family=multinomial("identity"))
  
mod2 <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=N,homogeneous=FALSE)
fm3 <- fit(mod2,emc=em.control(rand=F))

lb <- c(0,0,rep(0,4),-Inf,-Inf,-Inf,-Inf,-Inf,-Inf)
ub <- c(1,1,rep(1,4),Inf,Inf,Inf,Inf,Inf,Inf)

# uses optimatization on the basis of smoother functions
fm3.r <- fit(mod2,method='rsolnp',solnpcntrl=list(rho = 2, outer.iter = 600, inner.iter = 1200, 
                                                  delta = 1e-10, tol = 1e-10,trace=T),conrows=diag(12),
             conrows.lower=lb,
             conrows.upper=ub)

single.fit <- solnp(c(2,1,1),obj.seq,y=state.2,w=1,LB=c(0.1,0.1,0.1),UB=c(100,100,100))
#################################################################################################


##########################################################################################################

# extra stuff

#reparameterization function
#repara <- function(par){
#  new.par <- vector(length=length(par))
# new.par[1] <- sqrt(par[2])/par[1]
#  new.par[2] <- sqrt(par[2])
# new.par[3] <- par[3]
#  return(new.par)
# }

#check simulations - correct
#fit <- fitdist(data=a,distr="SWald",start=list(mu=exp(object@parameters$mu),lambda=exp(object@parameters$lambda),threshold=object@parameters$threshold),method="mle")

dSWald <- function(x,mu,lambda,threshold){
  res <- statmod::dinvgauss(x-threshold, mean=mu,shape=lambda)
  res
}

a <- fitdist(state.1, distr="SWald",method="mle",
             start=list(mu=1,
                        lambda=1,
                        threshold=1)) 
?fitdistr


#density function
dPNRWald <- function(t,a,v,sv,posdrift=TRUE)
  # density for single accumulator
{
  # Convert to Desmond and Yang's notation
  d <- v/a # mean of delta
  l <- a^2 # lambda assuming diffusive variance is 1
  v <- sv^2 # variance of delta
  sqrt(l/(2*pi*t^3*(v*t+1)))*exp(-(l*(d*t-1)^2)/(2*t*(v*t+1)))*
    pnorm((v+d)*sqrt(l/(t*v^2+v)))/pnorm(d*sqrt(l/v)) # normalize
}

t <- seq(0,140,0.1)
hist(RT,freq=F,ylim=c(0,0.3))
lines(dPNRWald(t,a=2,v=0.1,sv=0,posdrift=TRUE))

# from github
library(seqmodels)

# with dinvgauss - note that it uses slightly different parameterization: kappa refers to threshold
# xi to drift rate and tau is the shift in response times. 
# optional sigma for the within-trial variability of the drift rate
?dinvgauss
dinvgauss(t, kappa=2, xi, tau = 1, sigma = 0,
          ln = FALSE)

