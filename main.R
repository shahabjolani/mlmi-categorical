pack <- lapply(c("mvtnorm", "mice", "MASS", "lme4", "msm", "broom.mixed", "micemd", "tidyverse"), require, character.only = TRUE)

source("functions.R")
source("mlmi-binomial.R")
source("mlmi-poisson.R")

# this function runs the specified simulation scenario with 2 replications 
res.sim <- sim(rep = 2, k = 20, model = "bin", outcome = "bin", sys = .1, omega = "moderate")
# this functino extracts the corresponding results 
extract(res.sim)
