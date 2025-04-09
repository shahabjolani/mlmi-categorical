pack <- lapply(c("mvtnorm", "mice", "MASS", "lme4", "msm", "broom.mixed", "micemd", "tidyverse"), require, character.only = TRUE)

res.sim <- sim(rep = 2, k = 20, model = "bin", outcome = "bin", sys = .1, omega = "moderate")
extract(res.sim)
