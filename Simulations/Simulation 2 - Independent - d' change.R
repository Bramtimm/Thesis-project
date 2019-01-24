####################################################################
## Project: RM Thesis - Final Project
## Script purpose: Simulation 2 - (independent state - varying d')
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

# make sure that dependencies are loaded
source('dependencies.R')
source('functions.R')

# set transition matrix, N and Number of Sequences
trans <- matrix(rep(0.5,4),ncol=2,nrow=2,byrow=T)
N <- 10000
NoS <- 100

# simulates HMM chains according to trans,N,NoS
chains <- simulate.HMM(trans,N,NoS)
