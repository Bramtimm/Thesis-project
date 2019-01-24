####################################################################
## Project: RM Thesis - Final Project
## Script purpose: Simulation 4 - (unstable state - varying N)
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

# make sure that dependencies are loaded
source('dependencies.R')
source('functions.R')

# set transition matrix, N and Number of Sequences
trans <- matrix(c(0.3,0.7,0.7,0.3),ncol=2,nrow=2,byrow=T)
N <- seq(50,10000,50)
NoS <- 200

# simulates HMMs with different N and saves in list