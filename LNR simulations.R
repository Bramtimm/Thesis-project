############## Modeling with LNR ###############################################

# Author: Bram Timmers
# use: experimenting with LNR, creating functions for simulations 1,2,3
# Dependencies: DMC -> install_packages.R, dmc/dmc.R, load_model ("LNR","lnr.R")

################################################################################

# start with a two-response design 
factors=list(S=c("resp1","resp2"))
responses=c("RESP1","RESP2")
match.map=list(M=list(left="RESP1",right="RESP2"))
