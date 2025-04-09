datagen <- function(k = 30, omega = "weak", outcome = "bin", sigma = 1, ...){
  n <- 3000/k # of subjects per study
  
  # model specification
  # (nu_1k, nu_2k) ~ N_2(0, lambda) where lambda_11 = lambda_22 = .2 and lambda_12 = 0.005 (corr = 0.025) 
  # x1 ~ binary with logit{Pr(x1 = 1)} = gamma_01 + nu_1k, gamma_01 = 0.4
  # x2 ~ count with log{E(x2)} = gamma_02 + nu_2k, gamma_02 = 0
  # y ~ cont or binary, y = a0+ a1x1 + a2x2 + u0 + u1x1 + e with a = (-1.85, 1.05, -0.04) and sigma = 1 (for cont)
  # (u0, u1) ~ N_2(0, omega) where omega = (0.10, -0.024, 0.25) for weak (corr = -0.15) and
  #                                omega = (0.75, -0.3, 0.85) for moderate (corr = -0.375) and 
  #                                omega = (2.75, -1.82, 2.85) for strong (corr = -0.65) 
  
  lambda <- matrix(c(0.2, 0.005, 0.005, 0.2), 2, 2) 
  gamma <- c(0.4, 0.0)
  alpha <- c(-1.85, 1.05, -0.04)
  omega <- switch(omega, weak = matrix(c(0.10, -0.024, -0.024, 0.25), 2, 2), 
                  moderate = matrix(c(0.75, -0.30, -0.30, 0.85), 2, 2),
                  strong = matrix(c(2.75, -1.82, -1.82, 2.85), 2, 2))
  
  # generating data
  data <- NULL
  for (i in 1:k){
    tempy <- 0 # to prevent studies with zero event
    nu <- rmvnorm(1, mean = rep(0, length(diag(lambda))), sigma = lambda)
    link1 <- gamma[1] + nu[1]
    link2 <- gamma[2] + nu[2]
    # logit link/poisson link
    x1 <- rbinom(n, 1, exp(link1)/(1+exp(link1)))
    x2 <- rpois(n, exp(link2))
    
    tomega <- rmvnorm(1, mean = rep(0, length(diag(omega))), sigma = omega)
    link <- cbind(1, x1, x2)%*%alpha + cbind(1, x1)%*%t(tomega)
    if (outcome == "bin"){
      while (sum(tempy) == 0) tempy <- rbinom(n, 1, 1/(1 + exp(-link))) #logit link
    }else{tempy <- link + rnorm(n, sd = sigma)}
    data <- rbind(data, cbind(i, tempy, x1, x2))
  }
  colnames(data) <- c("st", "y", "x1", "x2")
  setup <- list(alpha = alpha, omega = c(diag(omega), omega[1,2], sigma))
  data <- list (setup = setup, data = data)
  return(data)
}

misgen <- function(data, sys = .30, spo = .15, mech = "mcar",...){
  delta <- matrix(0, nrow = 4, ncol = 2)
  if (mech == "mar"){
    delta[2,1] <- 1 
    delta[2,2] <- 0.5
  }
  delta[1,1] <- log(spo/(1 - spo)) - delta[2,1]*mean(data[,"y"]) - delta[3,1]*mean(data[,"x1"]) - delta[4,1]*mean(data[,"x2"])
  delta[1,2] <- log(spo/(1 - spo)) - delta[2,2]*mean(data[,"y"]) - delta[3,2]*mean(data[,"x1"]) - delta[4,2]*mean(data[,"x2"])
  logit <- as.matrix(cbind(1, data[,-1]))%*%delta
  p <- exp(logit)/(1 + exp(logit))
  data[(runif(nrow(p)) > p[,1]) == FALSE, "x1"] <- NA
  data[(runif(nrow(p)) > p[,2]) == FALSE, "x2"] <- NA
  sys1 <- sample(unique(data[,"st"]), round(length(unique(data[,"st"]))*sys))
  sys2 <- sample(unique(data[,"st"]), round(length(unique(data[,"st"]))*sys)) 
  data[data[,"st"] %in% unique(data[,"st"])[sys1], "x1"] <- NA
  data[data[,"st"] %in% unique(data[,"st"])[sys2], "x2"] <- NA
  return(data)
}

ana.cont <- function(...){
  out <- datagen(...)
  cdata <- as.data.frame(out$data)
  # REF - no missing data
  fit.ref <- lmer(y ~ x1 + x2 + (1 + x1|st), data = cdata)
  res.ref <- list(fix = fixef(fit.ref), sefix = sqrt(diag(vcov(fit.ref))), ran = as.data.frame(VarCorr(fit.ref))[,"vcov"])
  # creating missing data
  mdata <- data <- as.data.frame(misgen(cdata, ...))
  # CC, here 'data' is complete cases
  data <- data[apply(!is.na(data), 1, sum) == 4, ]
  fit.cc <- lmer(y ~ x1 + x2 + (1 + x1|st), data = data)
  res.cc <- list(fix = fixef(fit.cc), sefix = sqrt(diag(vcov(fit.cc))), ran = as.data.frame(VarCorr(fit.cc))[,"vcov"])
  # SMI - stratified MI with dummies
  data <- mdata
  data[,"st"] <- as.factor(data[,"st"])
  data[,"x1"] <- as.factor(data[,"x1"])
  imp1 <- mice(data, print = F)
  suppressMessages(fit1 <- with(imp1, lmer(y ~ x1 + x2 + (1 + x1|st))))
  fit.smi <- pool(fit1)
  ran.smi <- fit1$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.smi <- list(fix = fit.smi$pooled$estimate, sefix = sqrt(fit.smi$pooled$t), ran = ran.smi)
  # MLMI
  data <- mdata
  pred <- make.predictorMatrix(data) 
  pred[pred[,"st"] != 0,"st"] <- -2
  pred[pred[,"y"] != 0,"y"] <- 2
  meth <- c("", "", "2l.glmer.binomial", "2l.glmer.poisson")
  suppressMessages(imp <- mice(data, print = F, pred=pred, method = meth, maxit = 10))
  suppressMessages(fit <- with(imp, lmer(y ~ x1 + x2 + (1 + x1|st))))
  fit.mlmi <- pool(fit)
  ran.mlmi <- fit$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.mlmi <- list(fix = fit.mlmi$pooled$estimate, sefix = sqrt(fit.mlmi$pooled$t), ran = ran.mlmi)
  # 2stage
  data <- mdata
  data[,"x1"] <- as.factor(data[,"x1"])
  pred <- make.predictorMatrix(data) 
  pred[pred[,"st"] != 0,"st"] <- -2
  pred[pred[,"y"] != 0,"y"] <- 2
  meth <- c("", "", "2l.2stage.bin", "2l.2stage.pois")
  suppressMessages(imp <- mice(data, print = F, pred=pred, method = meth, maxit = 10))
  suppressMessages(fit <- with(imp, lmer(y ~ x1 + x2 + (1 + x1|st))))
  fit.stage <- pool(fit)
  ran.stage <- fit$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.stage <- list(fix = fit.stage$pooled$estimate, sefix = sqrt(fit.stage$pooled$t), ran = ran.stage)
  # all results
  res <- list(ref = res.ref, cca = res.cc, smi = res.smi, mlmi = res.mlmi, twostage = res.stage)
  output <- list(cdata = cdata, mdata = mdata, setup = out$setup, results = res)
}

ana.bin <- function(...){
  out <- datagen(...)
  cdata <- as.data.frame(out$data)
  # REF - no missing data
  fit.ref <- glmer(y ~ x1 + x2 + (1 + x1|st), family = binomial(link = "logit"), data = cdata)
  res.ref <- list(fix = fixef(fit.ref), sefix = sqrt(diag(vcov(fit.ref))), ran = as.data.frame(VarCorr(fit.ref))[,"vcov"])
  # creating missing data
  mdata <- data <- as.data.frame(misgen(cdata, ...))
  # CC, here 'data' is complete cases
  data <- data[apply(!is.na(data), 1, sum) == 4, ]
  fit.cc <- glmer(y ~ x1 + x2 + (1 + x1|st), family = binomial(link = "logit"), data = data)
  res.cc <- list(fix = fixef(fit.cc), sefix = sqrt(diag(vcov(fit.cc))), ran = as.data.frame(VarCorr(fit.cc))[,"vcov"])
  # SMI - stratified MI with dummies 
  data <- mdata
  data[,"st"] <- as.factor(data[,"st"])
  data[,"x1"] <- as.factor(data[,"x1"])
  imp1 <- mice(data, print = F)
  suppressMessages(fit1 <- with(imp1, glmer(y ~ x1 + x2 + (1 + x1|st), family = binomial(link = "logit"))))
  fit.smi <- pool(fit1)
  ran.smi <- fit1$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.smi <- list(fix = fit.smi$pooled$estimate, sefix = sqrt(fit.smi$pooled$t), ran = ran.smi)
  # MLMI
  data <- mdata
  pred <- make.predictorMatrix(data) 
  pred[pred[,"st"] != 0,"st"] <- -2
  pred[pred[,"y"] != 0,"y"] <- 2
  meth <- c("", "", "2l.glmer.binomial", "2l.glmer.poisson")
  suppressMessages(imp <- mice(data, print = F, pred=pred, method = meth, maxit = 10))
  suppressMessages(fit <- with(imp, glmer(y ~ x1 + x2 + (1 + x1|st), family = binomial(link = "logit"))))
  fit.mlmi <- pool(fit)
  ran.mlmi <- fit$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.mlmi <- list(fix = fit.mlmi$pooled$estimate, sefix = sqrt(fit.mlmi$pooled$t), ran = ran.mlmi)
  # 2stage
  data <- mdata
  data[,"x1"] <- as.factor(data[,"x1"])
  pred <- make.predictorMatrix(data) 
  pred[pred[,"st"] != 0,"st"] <- -2
  pred[pred[,"y"] != 0,"y"] <- 2
  meth <- c("", "", "2l.2stage.bin", "2l.2stage.pois")
  suppressMessages(imp2 <- mice(data, print = F, pred=pred, method = meth, maxit = 10))
  suppressMessages(fit2 <- with(imp2, glmer(y ~ x1 + x2 + (1 + x1|st), family = binomial(link = "logit"))))
  fit.stage <- pool(fit2)
  ran.stage <- fit2$analyses %>% map(VarCorr) %>% as.data.frame %>% select(starts_with("vcov")) %>% apply(., 1, mean)
  res.stage <- list(fix = fit.stage$pooled$estimate, sefix = sqrt(fit.stage$pooled$t), ran = ran.stage)
  # all results
  res <- list(ref = res.ref, cca = res.cc, smi = res.smi, mlmi = res.mlmi, twostage = res.stage)
  output <- list(cdata = cdata, mdata = mdata, setup = out$setup, results = res)
}

extract <- function(tst, alpha = .05){
  true <- tst[[1]]$setup$alpha
  trueran <- sqrt(tst[[1]]$setup$omega[-3])
  zed <- qnorm(1 - alpha/2)
  t <- length(tst)
  methods <- names(tst[[t]][["results"]])
  
  # summarizing the simulation results
  fixef <- sefixef <- sdfixef <- cr <- rmsefixef <- NULL 
  ranef <- sqranef <- rmseranef <- sqrmseranef <- NULL
  for (j in methods){
    # FIXED
    # estimate: fixed-effect coefficients
    fixef <- rbind(fixef, apply(sapply(1:t, function(i) tst[[i]][["results"]][[j]][["fix"]]), 1, mean))
    # model SE: sqrt of mean of estimated standard error (SE)^2
    sefixef <- rbind(sefixef, sqrt(apply(sapply(1:t, function(i) tst[[i]][["results"]][[j]][["sefix"]]^2), 1, mean)))
    # emp SE: standard deviation of the fixed-effect coefficients (SD)
    sdfixef <- rbind(sdfixef, apply(sapply(1:t, function(i) tst[[i]][["results"]][[j]][["fix"]]), 1, sd))
    # coverage rate %
    cr <- rbind(cr, apply(sapply(1:t, function(i) ifelse(
      true > tst[[i]][["results"]][[j]][["fix"]] - zed*tst[[i]][["results"]][[j]][["sefix"]] & 
        true < tst[[i]][["results"]][[j]][["fix"]] + zed*tst[[i]][["results"]][[j]][["sefix"]], 1, 0)), 1, mean))
    # rmse
    rmsefixef <- rbind(rmsefixef, sqrt(apply((sapply(1:t, function(i) tst[[i]][["results"]][[j]][["fix"]]) - true)^2, 1, mean)))
    
    
    # RANDOM
    if (length(tst[[1]][["results"]][["ref"]][["ran"]]) == 3) trueran <- trueran[1:2]
    # random-effects coefficients (between study heterogeneity) - sqrt(diag(omega)) + sigam
    sqranef <- rbind(sqranef, apply(sapply(1:t, function(i) sqrt(tst[[i]][["results"]][[j]][["ran"]][-3])), 1, mean))
    # rmse of random-effects coefficients (between study heterogeneity) omega + sigam
    sqrmseranef <- rbind(sqrmseranef, sqrt(apply((sapply(1:t, function(i) sqrt(tst[[i]][["results"]][[j]][["ran"]][-3])) - trueran)^2, 1, mean)))
  }
  
  # bias and relative bias
  bfixef <- sweep(fixef, 2, true)
  branef <- sweep(sqranef, 2, trueran)

  # fixed
  beta <- rbind(
    cbind(fixef[,1], bfixef[,1], sefixef[,1], sdfixef[,1], round(100*cr[,1],3), rmsefixef[,1]),
    cbind(fixef[,2], bfixef[,2], sefixef[,2], sdfixef[,2], round(100*cr[,2],3), rmsefixef[,2]),
    cbind(fixef[,3], bfixef[,3], sefixef[,3], sdfixef[,3], round(100*cr[,3],3), rmsefixef[,3])
  )
  rownames(beta) <- rep(methods, 3)
  colnames(beta) <- c("Estimate", "Bias", "Model SE", "Emp SE", "CR", "RMSE")
  
  # random
  omega <- rbind(
    cbind(sqranef[,1], branef[,1], sqrmseranef[,1]),
    cbind(sqranef[,2], branef[,2], sqrmseranef[,2])
  )
  rownames(omega) <- rep(methods, 2)
  colnames(omega) <- c("Estimate", "Bias", "RMSE")
  list(beta = beta, omega = omega)
}

sim <- function(rep = 1, model = "bin", ...){
  output <- vector("list", rep)
  if (model == "bin"){
    for (i in 1:rep) output[[i]] <- ana.bin(...)  
  }
  if (model == "cont"){
    for (i in 1:rep) output[[i]] <- ana.cont(...)  
  }
  return(output)
}
