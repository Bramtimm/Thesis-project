############## Dependecies of Thesis Projects ##################################

# Author: Bram Timmers
# use: Dependencies of Thesis Project
# Dependencies: none, are dependencies self

################################################################################

setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Dependencies/DMC/DMC-MBN18")

# dependencies of DMC: Dynamics of Choice software
source("packages/install_packages.R")

# clear directory and par
rm(list=ls()); par(mfrow = c(1,1))

# source functions needed for LNR
source("dmc/dmc.R")

load_model ("LNR","lnr.R")
