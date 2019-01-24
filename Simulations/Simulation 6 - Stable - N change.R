####################################################################
## Project: RM Thesis - Final Project
## Script purpose: Simulation 6 - (stable states - varying N)
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

# make sure that dependencies are loaded
source('dependencies.R')
source('functions.R')

# set transition matrix, N and Number of Sequences
trans <- matrix(c(0.9,0.1,0.1,0.9),ncol=2,nrow=2,byrow=T)
N <- seq(50,10000,50)
NoS <- 200

# simulates HMMs with different N and saves in list
