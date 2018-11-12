############## Modeling with LNR ###############################################

# Author: Bram Timmers
# use: experimenting with LNR, creating functions for simulations 1,2,3
# Dependencies: DMC -> install_packages.R, dmc/dmc.R, load_model ("LNR","lnr.R")

################################################################################

setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Thesis project")
source("dependencies.R")

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

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "meanlog.true"  "meanlog.false" "sdlog.true"    "sdlog.false"   "t0"           

# 1) In this example accuracy turns out to be around 75% - this is the same ex. as in tutorial of DMC
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
plot.score.dmc (data.model)

#################### Working with depmixS4 ################################

# Author: Bram Timmers
# Use: testing/how to use depmixS4
# dependencies: already sourced by dependencies.r

###########################################################################

# the normal fitdistr on true/false response.
fitdist(data.model[as.integer(factor(data.model[,2]))==as.integer(factor(data.model[,1])),3],"lnorm")
fitdist(data.model[as.integer(factor(data.model[,2]))!=as.integer(factor(data.model[,1])),3],"lnorm")

# note that RTs are merged in data.model[,3], now we specify our depmixS4 element # seems to work for lognormal distribution
mod.test <- depmix(list(RT~1,R~1), data = data.model, nstates = 2,family = list(gaussian(),multinomial()))

# only RTs
mod2.test <- depmix(RT~1, data = data.model, nstates = 2, family=gaussian(),instart = c(0.99, 0.01))

# fit
fm.mod1 <- fit(mod.test)
fm.mod2 <- fit(mod2.test)
summary(fm.mod1)
summary(fm.mod2)

# We need to work on some functions to include the lognormal distribution in the depmixS4 package
# Work on this in weekend

source("DepmixS4 - lnorm distribution.R")

RT <- data.model$RT

rModels <- list(
  list(
    lnorm(RT,pstart=c(-0.5,.1)),
    GLMresponse(formula=R~1, data=data.model, 
                family=multinomial("identity"), pstart=c(0.5,0.5))
  ),
  list(
    lnorm(RT,pstart=c(-1,.1)),
    GLMresponse(formula=R~1, data=data.model, 
                family=multinomial("identity"), pstart=c(.1,.9))
  )
)

transition <- list()
transition[[1]] <- transInit(~1,nstates=2,data=data.model)
transition[[2]] <- transInit(~1,nstates=2,data=data.model)


instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,family=multinomial("identity"), data=data.frame(rep(1,1)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=20000, 
                  homogeneous=FALSE)

fm3 <- fit(mod,emc=em.control(rand=FALSE))
summary(fm3)

## End(Not run)




