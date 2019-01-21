#################### Working with markovchain ################################

# Author: Bram Timmers
# Use: testing/how to simulate
# dependencies: already sourced by dependencies.r

###########################################################################

setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Thesis project")
source("dependencies.R")
source("DepmixS4 - lnorm distribution.R")
# source("DepmixS4 - Shifted_lognormal.R") - contains GAMLSS specifications
source("DepmixS4- Shifted_lognormal_fitdist trials.R")

# Initiate a single acccumulator LNR model
model.single <- model.dmc(type="lnr",constants=c(st0=0),
                          p.map=list(meanlog="1",sdlog="1",t0="1",st0="1"),
                          match.map=list(M=list(resp1="RESP1")),
                          factors=list(S=c("resp1")),
                          responses=c("RESP1"))

# name HMM states
statesNames <- c("strategy 1","strategy 2")

#N
N <- c(10000,1000,10000)

# Scenario 1 - independent states 
RT <- vector(length=N[2])

#initiate a markov chain 
HMM.1 <- new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.7,0.3), nrow=2, byrow=TRUE,
                                                    dimnames = list(statesNames,statesNames)))

#sample a random sequence
HMM.1.seq <- markovchainSequence( N[2], HMM.1, t0 = sample(HMM.1@states,1))

#make numeric
HMM.1.seq <- as.numeric(factor(HMM.1.seq))

N.state.1 <- sum(HMM.1.seq==1)
N.state.2 <- sum(HMM.1.seq==2)

#initial parameter values
p.vector.1 <- c(meanlog=-2,sdlog=2,t0=0) 
p.vector.2 <- c(meanlog=1,sdlog=2,t0=0)

# simulate response times 
state.1 <- data.model.dmc(simulate.dmc(p.vector.1,model.single,n=N.state.1),model.single)[,3]
state.2 <- data.model.dmc(simulate.dmc(p.vector.2,model.single,n=N.state.2),model.single)[,3]

RT[HMM.1.seq==1] <- state.1
RT[HMM.1.seq==2] <- state.2

data.test <- cbind(RT,HMM.1.seq)
data.test <- data.frame(data.test)
names(data.test) <- c("RT","states")

##### fit dist check #############################

# lnorm
fitdistr(data.test[(as.integer(factor(data.test[,2]))==1),1],dens=dlnorm,lower=c(-Inf,0),start=list(meanlog=0,sdlog=1))
fitdistr(data.test[(as.integer(factor(data.test[,2]))==2),1],dens=dlnorm,lower=c(-Inf,0),start=list(meanlog=0,sdlog=1))
fitdist(data.test[(as.integer(factor(data.test[,2]))==1),1],distr="lnorm")
fitdist(data.test[(as.integer(factor(data.test[,2]))==2),1],distr="lnorm")


#norm
fitdistr(log(data.test[(as.integer(factor(data.test[,2]))==1),1]),dens=dnorm,lower=c(-Inf,0),start=list(mean=0,sd=1))
fitdistr(log(data.test[(as.integer(factor(data.test[,2]))==2),1]),dens=dnorm,lower=c(-Inf,0),start=list(mean=0,sd=1))
fitdist(log(data.test[(as.integer(factor(data.test[,2]))==1),1]),distr="norm")
fitdist(log(data.test[(as.integer(factor(data.test[,2]))==2),1]),distr="norm")

#shifted lnorm - works when bounding the boundary with upper limits
fitdistr(data.test[(as.integer(factor(data.test[,2]))==1),1],dens=dSLOGNO,lower=c(-Inf,0.1,0),upper=c(Inf,Inf,min(state.1)),start=list(meanlog=0,sdlog=1,shift=0),method="L-BFGS-B",control=list(maxit=500,pgtol=1e-20,trace=TRUE))
fitdistr(data.test[(as.integer(factor(data.test[,2]))==2),1],dens=dSLOGNO,lower=c(-Inf,0.1,-1),upper=c(Inf,Inf,min(state.2)),start=list(meanlog=0,sdlog=1,shift=0),method="L-BFGS-B",control=list(maxit=500,pgtol=1e-20,trace=TRUE))


# with GAMLSSML - needs better strating values - doesn't work
#fit <- gamlssML(data.test[(as.integer(factor(data.test[,2]))==1),1]~1,weights=NULL,family=SLOGNO(),
 #               control=gamlss.control(c.crit=1e-10,n.cyc=10000,trace=FALSE),
  #              mu.start=0,
   #             sigma.start=1,
    #            nu.start=0)

# ps optim 
source("obj.R")

# limits
lb <- c(-10,0,-1e-1)
ub <- c(10,10,min(state.1))

#controls
max.its <- 500
nparticles <- 100
#ps optim - lower bound includes zero
res<-psoptim(rep(NA,3),obj,lower=lb,upper=ub,rt=state.1,w=rep(1,length(state.1)),control=list(trace=T,hybrid="improved",hybrid.control=list(pgtol=1e-8)))
        
#DEoptim - lower bound excludes zero
res <-DEoptim(obj, lower=lb, upper=ub, control=DEoptim.control(
  itermax=max.its, trace=T, NP=nparticles), rt=state.1)
       
        



############# depmix fit check #####################
RT <- data.test$RT

rModels2 <- list(
  list(
    lnorm(RT,pstart=c(-1,1))
  ),
  list(
    lnorm(RT,pstart=c(1,1))
  )
)

transition <- list()
transition[[1]] <- transInit(~1,nstates=2,data=data.test)
transition[[2]] <- transInit(~1,nstates=2,data=data.test)

instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart, data=data.frame(rep(1,1)))

mod2 <- makeDepmix(response=rModels2,transition=transition,prior=inMod,ntimes=N[2], 
                   homogeneous=FALSE)

fm3.r <- fit(mod2,emc=em.control(rand=FALSE))

# for plotting purposes
hist(data.test[,1],breaks=100,freq=FALSE,ylim=c(0,4),xlim=c(0,4))
lines(density(data.test[data.test[,2]==1,1]),lwd=2)
lines(density(data.test[data.test[,2]==2,1]),lwd=2,col="red")

#for plotting purposes
par.state1 <- fm3.r@response[[1]][[1]]@parameters
par.state2 <- fm3.r@response[[2]][[1]]@parameters

x <- seq(0,4,0.001)
lines(x,dlnorm(x,par.state1[[1]],exp(par.state1[[2]])),col="blue",lwd=2)
lines(x,dlnorm(x,par.state2[[1]],exp(par.state2[[2]])),col="green",lwd=2)

start <- c(0.5,1.1,0.1)
names(start) <- c("meanlog","sdlog","shift") 
start <- as.list(start)
fitdist(data.test[(as.integer(factor(data.test[,2]))==1),1],"SLOGNO",start=start)
fitdist(data.test[(as.integer(factor(data.test[,2]))==2),1],"shifted_lnorm",start=start,lower=c(-Inf,0,0),upper=rep(Inf,3))

fitdist(data.test[(as.integer(factor(data.test[,2]))==1),1],distr="lnorm")
fitdist(data.test[(as.integer(factor(data.test[,2]))==2),1],distr="lnorm")

