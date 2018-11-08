############## Modeling with LNR ###############################################

# Author: Bram Timmers
# use: experimenting with LNR, creating functions for simulations 1,2,3
# Dependencies: DMC -> install_packages.R, dmc/dmc.R, load_model ("LNR","lnr.R")

################################################################################

# start with a two-response design 
factors=list(S=c("resp1","resp2"))
responses=c("RESP1","RESP2")
match.map=list(M=list(left="RESP1",right="RESP2"))

# NB: The LNR model has parameters meanlog, sdlog (see ?plnorm) and t0. 
#     The LNR model is not implemented in the rtdists package but 
#     "rtdists_extras.R" provides the necessary functions (rlnr,   
#     a random function, and n1PDF.lnr). 

model <- model.dmc(type="lnr",constants=c(st0=0),
                   p.map=list(meanlog="M",sdlog="M",t0="1",st0="1"),
                   match.map=list(M=list(left="RESP1",right="RESP2")),
                   factors=list(S=c("resp1","resp2")),
                   responses=c("RESP1","RESP2"))

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "meanlog.true"  "meanlog.false" "sdlog.true"    "sdlog.false"   "t0"           

# 1) In this example accuracy turns out to be around 75%
p.vector  <- c(meanlog.true=-1,meanlog.false=0,
               sdlog.true=1,sdlog.false=1,t0=.2)
data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e4),model)
