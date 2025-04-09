mice.impute.2l.glmer.binomial <- function (y, ry, x, type, wy = NULL, intercept = TRUE, ...) 
{
  #install.on.demand("lme4", ...)
  if (is.null(wy)) 
    wy <- !ry
  if (intercept) {
    x <- cbind(1, as.matrix(x))
    type <- c(2, type)
    names(type)[1] <- colnames(x)[1] <- "(Intercept)"
  }
  clust <- names(type[type == -2])
  rande <- names(type[type == 2])
  fixe <- names(type[type > 0])
  lev <- unique(x[, clust])
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]
  fr <- ifelse(length(rande) > 1, paste("+ ( 1 +", paste(rande[-1L], 
                                                         collapse = "+")), "+ ( 1 ")
  randmodel <- paste("yobs ~ ", paste(fixe[-1L], collapse = "+"), 
                     fr, "|", clust, ")")
  
  #step1: fitting glmer to the observed data
  suppressWarnings(fit <- try(lme4::glmer(formula(randmodel), data = data.frame(yobs, xobs), 
                                          nAGQ = 0, family = binomial(link = "logit")),
  silent = TRUE))
  if (!is.null(attr(fit, "class"))) {
    if (attr(fit, "class") == "try-error") {
      warning("glmer does not run. Simplify imputation model")
      return(y[wy])
    }
  }

  # calculation of bound M for the accept-reject sampling
  probM <- predict(fit, type = "response")
  probM <- probM*yobs + (1 - probM)*(1-yobs) #likelihood of each individual (for a binary variable)

  #step2: extracting and drawing parameters
  beta <- lme4::fixef(fit)
  RX <- lme4::getME(fit, "RX")
  covmat <- chol2inv(RX)
  rv <- t(chol(covmat))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  
  rancoef <- as.matrix(lme4::ranef(fit)[[1]])
  lambda <- t(rancoef) %*% rancoef
  df.psi <- nrow(rancoef)
  temp.psi.star <- stats::rWishart(1, df.psi, diag(nrow(lambda)))[, , 1]
  temp <- MASS::ginv(lambda)
  ev <- eigen(temp)
  if (sum(ev$values > 0) == length(ev$values)) {
    deco <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
    psi.star <- MASS::ginv(deco %*% temp.psi.star %*% t(deco))
  } else {
    try(temp.svd <- svd(lambda))
    if (class(temp.svd) != "try-error") {
      deco <- temp.svd$u %*% diag(sqrt(temp.svd$d), nrow = length(temp.svd$d))
      psi.star <- MASS::ginv(deco %*% temp.psi.star %*% 
                               t(deco))
    } else {
      psi.star <- temp
      warning("psi fixed to estimate")
    }
  }
  
  #step3: imputation of b and yobs 
  for (jj in lev){
    myi <- matrix(0, nrow = nrow(psi.star), ncol = 1)
    vyi <- psi.star
    vyi <- vyi - upper.tri(vyi) * vyi + t(lower.tri(vyi) * vyi)
    deco1 <- eigen(vyi)
    if (sum(deco1$values > 0) == length(deco1$values)) {
      A <- deco1$vectors %*% sqrt(diag(deco1$values, nrow = length(deco1$values)))
      bi.star <- myi + A %*% rnorm(length(myi))
    } else {
      try(deco1 <- svd(vyi))
      if (class(deco1) != "try-error") {
        A <- deco1$u %*% sqrt(diag(deco1$d, nrow = length(deco1$d)))
        bi.star <- myi + A %*% rnorm(length(myi))
      }
      else {
        bi.star <- myi
        warning("b_", jj, " fixed to estimate")
      }
    }
    if (jj %in% unique(xobs[, clust])){
      reject <- accept <- 0
      while(accept == 0 & reject < 250){
        w <- runif(1)
        f <- as.vector(as.matrix(X[!wy & x[, clust] == jj, , drop = FALSE]) %*% beta.star + 
                         as.matrix(Z[!wy & x[, clust] == jj, , drop = FALSE]) %*% as.matrix(bi.star))
        f <- exp(f)/(1+exp(f))
        f <- f*y[!wy & x[, clust] == jj] + (1-f)*(1- y[!wy & x[, clust] == jj]) # above formula is fine too
        M <- probM[x[!wy,clust] == jj]
        ratio <- prod(f/M)
        if (is.na(ratio)) ratio <- 0 
        if(w > ratio){
          reject <- reject + 1
          if(!all(bi.star == myi)) bi.star <- myi + A %*% rnorm(length(myi))
        }else{accept <- 1} 
      }
      if (reject == 250) warning("b_", jj, " fixed to unconditional draw")
    }
    logit <- as.vector(as.matrix(X[wy & x[, clust] == jj, , drop = FALSE]) %*% beta.star + 
                         as.matrix(Z[wy & x[, clust] == jj, , drop = FALSE]) %*% as.matrix(bi.star))
    prob <- exp(logit)/(1+exp(logit))
    y[wy & x[, clust] == jj] <- rbinom(length(prob), 1, prob)
  }
  return(y[wy])
}