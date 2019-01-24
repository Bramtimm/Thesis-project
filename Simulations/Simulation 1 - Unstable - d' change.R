####################################################################
## Project: RM Thesis - Final Project
## Script purpose: Simulation 1 - (unstable condition - varying d')
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

# make sure that dependencies are loaded
source('dependencies.R')
source('functions.R')

# set transition matrix, N and Number of Sequences
trans <- matrix(c(0.3,0.7,0.7,0.3),ncol=2,nrow=2,byrow=T)
N <- 10000
NoS <- 100

# simulates HMM chains according to trans,N,NoS
chains <- simulate.HMM(trans,N,NoS)


