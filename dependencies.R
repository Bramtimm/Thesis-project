####################################################################
## Project: RM Thesis - Final Project
## Script purpose: General software dependencies and packages for scripts
## Date: 24-1-2019
## Author: Bram Timmers
####################################################################

setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Dependencies/DMC/DMC-MBN18")

# dependencies of DMC: Dynamics of Choice software
source("packages/install_packages.R")

# clear directory and par
rm(list=ls()); par(mfrow = c(1,1))

# source functions needed for LNR
source("dmc/dmc.R")
load_model ("LNR","lnr.R")

# function for loading and checking packages
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# install packages - note that seqmodels needs to be imported directly from github
packages <- c('depmixS4','fitdistrplus','gamlss','gamlss.dist','markovchain','brms','seqmodels')
sapply(packages,pkgTest)

# reset working directory
setwd("C:/Users/Bram/Documents/RM Thesis/Analyses/Thesis project")