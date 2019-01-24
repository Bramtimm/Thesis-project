####################################################################
## Project: RM Thesis - Final Project
## Script purpose: Simulates HMM's according to given transition matrix
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

simulate.HMM <- function(trans,N,NoS){
  statesNames <- c("strategy 1","strategy 2")
  dimnames(trans) <- list(statesNames,statesNames)
  chains <- matrix(NA,nrow=NoS,ncol=N)
  
  #initiate Markov Chain
  HMM <- new("markovchain", transitionMatrix = trans)
  
  #simulate NoS chains
  chains <- lapply(1:NoS,markovchainSequence,n=N,markovchain=HMM,t0=sample(HMM@states,1),include.t0=FALSE)
  
  #binarize
  chains <- lapply(chains,function(chains)as.numeric(as.factor(chains)))
  chains <- matrix(unlist(chains),ncol=N,byrow=T)
  chains
}