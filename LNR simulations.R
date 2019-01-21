############## Modeling with LNR ###############################################

# Author: Bram Timmers -- mostly from tutorial DMC
# use: experimenting with LNR, creating functions for simulations 1,2,3
# Dependencies: DMC -> install_packages.R, dmc/dmc.R, load_model ("LNR","lnr.R")

################################################################################

setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Thesis project")
source("dependencies.R")
source("DepmixS4 - lnorm distribution.R")
source("DepmixS4 - Shifted_lognormal.R")

# start with a two-response design 
factors=list(S=c("resp1","resp2"))
responses=c("RESP1","RESP2")
match.map=list(M=list(resp1="RESP1",resp2="RESP2"))

# NB: The LNR model has parameters meanlog, sdlog (see ?plnorm) and t0. 
#     The LNR model is not implemented in the rtdists package but 
#     "rtdists_extras.R" provides the necessary functions (rlnr,   
#     a random function, and n1PDF.lnr). 

model <- model.dmc(type="lnr",constants=c(st0=0),
                   p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(resp1="RESP1",resp2="RESP2")),
                   factors=list(S=c("resp1","resp2")),
                   responses=c("RESP1","RESP2"))

model.extended.t0 <- model.dmc(type="lnr",constants=c(st0=0),
                              p.map=list(meanlog="M",sdlog="M",t0="M",st0="1"),
                              match.map=list(M=list(resp1="RESP1",resp2="RESP2")),
                              factors=list(S=c("resp1","resp2")),
                              responses=c("RESP1","RESP2"))


# Parameter vector names are: ( see attr(,"p.vector") ) --> extension by t0 only leads to one accumulator winning more.
# [1] "meanlog.true"  "meanlog.false" "sdlog.true"    "sdlog.false"   "t0"           

# 1) In this example accuracy turns out to be around 75% - this is the same ex. as in tutorial of DMC
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=0.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e3),model)

plot.score.dmc(data.model)

## extended

p.vector.ex  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0.true=2,t0.false=.2)
data.model.ex <- data.model.dmc(simulate.dmc(p.vector.ex,model.extended.t0,n=1e3),model)

plot.score.dmc(data.model.ex)


#################### Working with depmixS4 ################################

# Author: Bram Timmers
# Use: testing/how to use depmixS4
# dependencies: already sourced by dependencies.r

###########################################################################

# data.model <- data.model.ex #change to see extended model as simulated above

# the normal fitdistr on true/false response.
fitdist(data.model[as.integer(factor(data.model[,2]))==as.integer(factor(data.model[,1])),3],"lnorm")
fitdist(data.model[as.integer(factor(data.model[,2]))!=as.integer(factor(data.model[,1])),3],"lnorm")

# add error category 
data.model$error <- ifelse(as.integer(factor(data.model[,2]))==as.integer(factor(data.model[,1])),1,0)

#transform from lognormal to normal
data.model$RTNorm <- log(data.model$RT)

#simple plot of transformed RTs
hist(data.model$RTNorm,breaks=20,freq=FALSE,ylim=c(0,2))
lines(density(data.model[as.integer(factor(data.model[,2]))==as.integer(factor(data.model[,1])),5]),lwd=2)
lines(density(data.model[as.integer(factor(data.model[,2]))!=as.integer(factor(data.model[,1])),5]),lwd=2,col="red")

# note that RTs are merged in data.model[,3], now we specify our depmixS4 element # seems to work for normal distribution
mod.test <- depmix(list(RTNorm~1,error~1), data = data.model, nstates = 2,family = list(gaussian(),multinomial("identity")))

# only RTs
mod2.test <- depmix(RTNorm~1, data = data.model, nstates = 2, family=gaussian(),instart = c(0.99, 0.01))

# fit
fm.mod1 <- fit(mod.test)
fm.mod2 <- fit(mod2.test)
summary(fm.mod1)
summary(fm.mod2)

#for plotting purposes
par.state1 <- fm.mod1@response[[1]][[1]]@parameters
par.state2 <- fm.mod1@response[[2]][[1]]@parameters

x <- seq(-2,5,0.001)
lines(x,dnorm(x,par.state1[[1]],par.state1[[2]]),col="blue",lwd=2)
lines(x,dnorm(x,par.state2[[1]],par.state2[[2]]),col="green",lwd=2)

# We need to work on some functions to include the lognormal distribution in the depmixS4 package

# This works and is correct! :-)

RT <- data.model$RT

# this gives weird transition matrix # something goes wrong with 1st model. 

rModels <- list(
  list(
    lnorm(RT,pstart=c(-1,1)),
    GLMresponse(formula=error~1, data=data.model, 
                family=multinomial("identity"), pstart=c(0.1,0.9))
  ),
  list(
    lnorm(RT,pstart=c(0,1)),
    GLMresponse(formula=error~1, data=data.model, 
                family=multinomial("identity"), pstart=c(0.9,0.1))
  )
)

rModels2 <- list(
  list(
    lnorm(RT,pstart=c(-0.5,1))
  ),
  list(
    lnorm(RT,pstart=c(0.5,1))
  )
)

transition <- list()
transition[[1]] <- transInit(~1,nstates=2,data=data.model)
transition[[2]] <- transInit(~1,nstates=2,data=data.model)


instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,family=multinomial("identity"), data=data.frame(rep(1,1)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=2000, 
                  homogeneous=FALSE)

fm3 <- fit(mod,emc=em.control(rand=FALSE))

inMod <- transInit(~1,ns=2,ps=instart, data=data.frame(rep(1,1)))

mod2 <- makeDepmix(response=rModels2,transition=transition,prior=inMod,ntimes=2000, 
                  homogeneous=FALSE)

fm3.r <- fit(mod2,emc=em.control(rand=FALSE))

# for plotting purposes
hist(data.model$RT,breaks=40,freq=FALSE,ylim=c(0,4),xlim=c(0,4))
lines(density(data.model[as.integer(factor(data.model[,2]))==as.integer(factor(data.model[,1])),3]),lwd=2)
lines(density(data.model[as.integer(factor(data.model[,2]))!=as.integer(factor(data.model[,1])),3]),lwd=2,col="red")

#for plotting purposes
par.state1 <- fm3@response[[1]][[1]]@parameters
par.state2 <- fm3@response[[2]][[1]]@parameters

x <- seq(0,4,0.001)
lines(x,dlnorm(x,par.state1[[1]],exp(par.state1[[2]])),col="blue",lwd=2)
lines(x,dlnorm(x,par.state2[[1]],exp(par.state2[[2]])),col="green",lwd=2)

