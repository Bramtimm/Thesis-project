#################### Working with markovchain ################################

# Author: Bram Timmers
# Use: testing/how to simulate - Function
# dependencies: already sourced by dependencies.r

###########################################################################

source("dependencies.R")

model.single <- model.dmc(type="lnr",constants=c(st0=0),
                          p.map=list(meanlog="1",sdlog="1",t0="1",st0="1"),
                          match.map=list(M=list(resp1="RESP1")),
                          factors=list(S=c("resp1")),
                          responses=c("RESP1"))

# create function with arguments: p.vector.1,p.vector.2,N

simulate.HMM.LNR <- function(p.vector.1,p.vector.2,N,trans,NoS){
  
#checks for function arguments
  if(missing(p.vector.1) | missing(p.vector.2)){
    stop("model parameters are not specified for the accumulators;
         please specify parameters for both p.vector.1 and p.vector.2")
  }
  if(missing(N)){
    stop("please specify the number of trials,N, to be simulated")
  }
  if(length(p.vector.1) !=3 | length(p.vector.2)!=3){
    stop("please specify only mu0, sd0, and t0/ter for every accumulator")
  }
  if(missing(trans) | identical(dim(trans),rep(2,2))){
    stop("please make sure that the dimensions of the transition matrix are equal to all
         transition probabilities between states")
  } 
  if(missing(NoS)){
    stop("please specify the number of sequences of chains to be simulated")
  }
  

#set default parameters
statesNames <- c("strategy 1","strategy 2")
dimnames(trans) <- list(statesNames,statesNames)
N <- N
NoS <- NoS
chains <- list()

#initiate Markov Chain
HMM <- new("markovchain", transitionMatrix = trans)

#simulate NoS chains
chains <- lapply(1:NoS,markovchainSequence,n=N,markovchain=HMM,t0=sample(HMM@states,1))

#binarize
chains <- lapply(chains,function(chains)as.numeric(as.factor(chains)))

chains
}
