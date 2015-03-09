require(VGAM)


# # # # #
# Common 
# # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # #

get.O <- function(y){
    O <- function(Z.i, T.i){
        return(which(y$assign==Z.i & y$treat==T.i))
    }
    return(O)
}

norm.ratio <- function(Wi, bet1, sig1, bet2, sig2, Xi){
    return( sig2/sig1 * exp( -(Wi - sum(Xi * bet1))^2/(2*sig1^2) +
                            (Wi - sum(Xi * bet2))^2/(2*sig2^2)))
}

update.one.posterior.strata.proba <- function(stratum, strata, pis.t.idx, Wi, Xi, bet.list, sig.list){
  denom <- 1
  for(alt.stratum in strata[strata!=stratum]){

    pis.ratio <- pis.t.idx[[alt.stratum]] / pis.t.idx[[stratum]]

    denom <- denom + pis.ratio * norm.ratio(Wi,
                                            bet.list[[alt.stratum]], sig.list[[alt.stratum]],
                                            bet.list[[stratum]], sig.list[[stratum]],
                                            Xi)
    
  }
  return(1/denom)
}

estimate.strata.proportions <- function(pis, wgt){
  # pis = list( "cEE"= vector, ..., "nNN"= vector), where each vector has length nrow(X)
  # wgt = the design weights (a vector of length nrow(X) that can be extracted from the X matrix)
  res = list()
  total.wgt <- sum(wgt)
  for(stratname in names(pis)){
    res[[stratname]] <- sum(pis[[stratname]] * wgt) / total.wgt
  }
  return(res)
}




initialize.params <- function(W.obs, y, X){

    # getting the indices
    
    # AT
    idx.all <- which(y$assign==1 & y$treat==1)
    idx.AT.1 <- sample(idx.all, size=floor(length(idx.all)/2))
    idx.1.1 <- setdiff(idx.all, idx.AT.1)

    # C
    idx.all <- which(y$assign==0 & y$treat==0)
    idx.C.0 <- sample(idx.all, size=floor(length(idx.all)/2))
    idx.0.0 <- setdiff(idx.all, idx.C.0)
    idx.C.1 <- idx.1.1 

    # NT
    idx.all <- which(y$assign==1 & y$treat == 0)
    idx.NT.1 <- sample(idx.all, size=floor(length(idx.all)/2))
    idx.1.0 <- setdiff(idx.all, idx.NT.1)
    idx.NT.0 <- idx.0.0

    # D
    idx.all <- which(y$assign==0 & y$treat == 1)
    idx.D.0 <- sample(idx.all, size=floor(length(idx.all)/2))
    idx.0.1 <- setdiff(idx.all, idx.D.0)
    idx.D.1 <- idx.1.0

    # AT (end)
    idx.AT.0 <- idx.0.1

    # setting up the beta params
    lm.AT.1 <- lm(W.obs[idx.AT.1] ~ X[idx.AT.1,] - 1)
    bet.AT.1 <- as.vector(lm.AT.1$coefficients)
    sig.AT.1 <- as.vector(summary(lm.AT.1)$sigma)

    lm.AT.0 <- lm(W.obs[idx.AT.0] ~ X[idx.AT.0,] - 1)
    bet.AT.0 <- as.vector(lm.AT.0$coefficients)
    sig.AT.0 <- as.vector(summary(lm.AT.0)$sigma)

    lm.C.1 <- lm(W.obs[idx.C.1] ~ X[idx.C.1,] - 1)
    bet.C.1 <- as.vector(lm.C.1$coefficients)
    sig.C.1 <- as.vector(summary(lm.C.1)$sigma)

    lm.C.0 <- lm(W.obs[idx.C.0] ~ X[idx.C.0,] - 1)
    bet.C.0 <- as.vector(lm.C.0$coefficients)
    sig.C.0 <- as.vector(summary(lm.C.0)$sigma)
    
    lm.NT.1 <- lm(W.obs[idx.NT.1] ~ X[idx.NT.1,] - 1)
    bet.NT.1 <- as.vector(lm.NT.1$coefficients)
    sig.NT.1 <- as.vector(summary(lm.NT.1)$sigma)

    lm.NT.0 <- lm(W.obs[idx.NT.0] ~ X[idx.NT.0,] - 1)
    bet.NT.0 <- as.vector(lm.NT.0$coefficients)
    sig.NT.0 <- as.vector(summary(lm.NT.0)$sigma)

    lm.D.1 <- lm(W.obs[idx.D.1] ~ X[idx.D.1,] - 1)
    bet.D.1 <- as.vector(lm.D.1$coefficients)
    sig.D.1 <- as.vector(summary(lm.D.1)$sigma)

    lm.D.0 <- lm(W.obs[idx.D.0] ~ X[idx.D.0,] - 1)
    bet.D.0 <- as.vector(lm.D.0$coefficients)
    sig.D.0 <- as.vector(summary(lm.D.0)$sigma)

    theta.init <- NULL
    theta.init$beta <- list(AT=list(bet.AT.0,bet.AT.1),
                            NT=list(bet.NT.0,bet.NT.1),
                            C=list(bet.C.0,bet.C.1),
                            D=list(bet.D.0,bet.D.1))
    theta.init$sigma <- list(AT=list(sig.AT.0,sig.AT.1),
                             NT=list(sig.NT.0,sig.NT.1),
                             C=list(sig.C.0,sig.C.1),
                             D=list(sig.D.0,sig.D.1))
    theta.init$alpha <- list(AT=rep(1,5),
                             NT=rep(1,5),
                             C=rep(1,5),
                             D=rep(1,5))
    return(theta.init)
}



# # # # #
# Full
# # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # #

# E-step

# checked
update.posterior.strata.proba.full <- function(theta, X, W.obs, O.fn){
    strata.probas <- list("AT" = rep(0,nrow(X)),
                          "NT" = rep(0,nrow(X)),
                          "C"  = rep(0,nrow(X)),
                          "D"  = rep(0,nrow(X)))

    # O(1,1)  =>  G \in { AT, C }
    indices <- O.fn(1,1)
    strata <- c("AT", "C")

    for(index in indices){
        pis.t.idx <- list(AT = theta[["pi"]][["AT"]][index],
                          C  = theta[["pi"]][["C"]][index])
        bet.list <- list(AT = theta[["beta"]][["AT"]][[1 + 1]],
                         C = theta[["beta"]][["C"]][[1 + 1]])
        sig.list <- list(AT = theta[["sigma"]][["AT"]][[1 + 1]],
                         C = theta[["sigma"]][["C"]][[1 + 1]])

        strata.probas[["AT"]][index] <- 1 / (1 + pis.t.idx[["C"]] / pis.t.idx[["AT"]] *
                                             exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["C"]]),  sig.list[["C"]], log=T)
                                                 - dnorm(W.obs[index], sum(X[index,] * bet.list[["AT"]]),  sig.list[["AT"]], log=T)))

        strata.probas[["C"]][index] <- 1 / (1 + pis.t.idx[["AT"]] / pis.t.idx[["C"]] *
                                            exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["AT"]]),  sig.list[["AT"]], log=T)
                                                - dnorm(W.obs[index], sum(X[index,] * bet.list[["C"]]),  sig.list[["C"]], log=T)))
    }

    
    # O(1,0)  =>   G \in { NT, D }
    indices <- O.fn(1,0)
    strata <- c("NT", "D")

    for(index in indices){
        pis.t.idx <- list(NT = theta[["pi"]][["NT"]][index],
                          D  = theta[["pi"]][["D"]][index])
        bet.list <- list(NT = theta[["beta"]][["NT"]][[1 + 1]],
                         D = theta[["beta"]][["D"]][[1 + 1]])
        sig.list <- list(NT = theta[["sigma"]][["NT"]][[1 + 1]],
                         D = theta[["sigma"]][["D"]][[1 + 1]])

        strata.probas[["NT"]][index] <- 1 / (1 + pis.t.idx[["D"]] / pis.t.idx[["NT"]] *
                                             exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["D"]]),  sig.list[["D"]], log=T)
                                                 - dnorm(W.obs[index], sum(X[index,] * bet.list[["NT"]]),  sig.list[["NT"]], log=T)))

        strata.probas[["D"]][index] <- 1 / (1 + pis.t.idx[["NT"]] / pis.t.idx[["D"]] *
                                            exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["NT"]]),  sig.list[["NT"]], log=T)
                                                - dnorm(W.obs[index], sum(X[index,] * bet.list[["D"]]),  sig.list[["D"]], log=T)))
    }
    
    # O(0,1)  =>   G \in { AT, D }
    indices <- O.fn(0,1)
    strata <- c("AT", "D")

    for(index in indices){
        pis.t.idx <- list(AT = theta[["pi"]][["AT"]][index],
                          D  = theta[["pi"]][["D"]][index])
        bet.list <- list(AT = theta[["beta"]][["AT"]][[0 + 1]],
                         D = theta[["beta"]][["D"]][[0 + 1]])
        sig.list <- list(AT = theta[["sigma"]][["AT"]][[0 + 1]],
                         D = theta[["sigma"]][["D"]][[0 + 1]])

        strata.probas[["AT"]][index] <- 1 / (1 + pis.t.idx[["D"]] / pis.t.idx[["AT"]] *
                                             exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["D"]]),  sig.list[["D"]], log=T)
                                                 - dnorm(W.obs[index], sum(X[index,] * bet.list[["AT"]]),  sig.list[["AT"]], log=T)))
        
        strata.probas[["D"]][index] <- 1 / (1 + pis.t.idx[["AT"]] / pis.t.idx[["D"]] *
                                            exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["AT"]]),  sig.list[["AT"]], log=T)
                                                - dnorm(W.obs[index], sum(X[index,] * bet.list[["D"]]),  sig.list[["D"]], log=T)))

        
    }
    
    # O(0,0)  =>   G \in { NT, C }
    indices <- O.fn(0,0)
    strata <- c("NT", "C")

    for(index in indices){
        pis.t.idx <- list(NT = theta[["pi"]][["NT"]][index],
                          C  = theta[["pi"]][["C"]][index])
        bet.list <- list(NT = theta[["beta"]][["NT"]][[0 + 1]],
                         C = theta[["beta"]][["C"]][[0 + 1]])
        sig.list <- list(NT = theta[["sigma"]][["NT"]][[0 + 1]],
                         C = theta[["sigma"]][["C"]][[0 + 1]])

        strata.probas[["NT"]][index] <- 1 / (1 + pis.t.idx[["C"]] / pis.t.idx[["NT"]] *
                                             exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["C"]]),  sig.list[["C"]], log=T)
                                                 - dnorm(W.obs[index], sum(X[index,] * bet.list[["NT"]]),  sig.list[["NT"]], log=T)))
        
        strata.probas[["C"]][index] <- 1 / (1 + pis.t.idx[["NT"]] / pis.t.idx[["C"]] *
                                            exp(dnorm(W.obs[index], sum(X[index,] * bet.list[["NT"]]),  sig.list[["NT"]], log=T)
                                                - dnorm(W.obs[index], sum(X[index,] * bet.list[["C"]]),  sig.list[["C"]], log=T)))        }

    return(strata.probas)
}

# M step

fitfn.multinom <- function(alph.arg, X, strata.probas){
        alpha.AT <- c(alph.arg[seq(1,5)])
        alpha.NT <- c(alph.arg[seq(6,10)])
        alpha.C <- c(alph.arg[seq(11,15)])
        alpha.D <- c(0, 0, 0, 0, 0)
                    
        res <- 0
        eta.AT <- exp(as.vector(X %*% alpha.AT))
        eta.NT <- exp(as.vector(X %*% alpha.NT))
        eta.C <- exp(as.vector(X %*% alpha.C))
        eta.D <- exp(as.vector(X %*% alpha.D))

        denom <- eta.AT + eta.NT + eta.C + eta.D

        p.AT <- eta.AT/denom
        p.NT <- eta.NT/denom
        p.C <- eta.C/denom
        p.D <- eta.D/denom


        res <- sum(strata.probas[["AT"]] * log(p.AT) +
                   strata.probas[["NT"]] * log(p.NT) +
                   strata.probas[["C"]] * log(p.C) +
                   strata.probas[["D"]] * log(p.D))

        return(res)
    }


# checked (works with BFGS)
optimize.multinom.model.optim <- function(X, strata.probas, method="Nelder-Mead"){

    fitfn <- function(alph.arg){
        alpha.AT <- c(alph.arg[seq(1,5)])
        alpha.NT <- c(alph.arg[seq(6,10)])
        alpha.C <- c(alph.arg[seq(11,15)])
        alpha.D <- c(0, 0, 0, 0, 0)
                    
        res <- 0
        eta.AT <- exp(as.vector(X %*% alpha.AT))
        eta.NT <- exp(as.vector(X %*% alpha.NT))
        eta.C <- exp(as.vector(X %*% alpha.C))
        eta.D <- exp(as.vector(X %*% alpha.D))

        denom <- eta.AT + eta.NT + eta.C + eta.D

        p.AT <- eta.AT/denom
        p.NT <- eta.NT/denom
        p.C <- eta.C/denom
        p.D <- eta.D/denom


        res <- sum(strata.probas[["AT"]] * log(p.AT) +
                   strata.probas[["NT"]] * log(p.NT) +
                   strata.probas[["C"]] * log(p.C) +
                   strata.probas[["D"]] * log(p.D))

        return(res)
    }

    res <- optim(par=rep(0,15), fn=fitfn, method=method, control=list(fnscale=-1, maxit=10000))
    return(list(optim.res = res,
                alphas =
                list(AT = res$par[seq(1,5)],
                     NT = res$par[seq(6,10)],
                     C  = res$par[seq(11,15)],
                     D  = rep(0, 5))))
}


fitfn.normal <- function(pars, W.obs.is, X.is, strata.probas.is){
    bet.arg <- pars[seq(5)]
    sig <- pars[6]
    dens.ls <- (dnorm(W.obs.is, as.vector(X.is %*% bet.arg), sig, log=T) *
                strata.probas.is)        
    return(sum(dens.ls))        
}

# *!* gives same result for beta, but not for sigma!!
optimize.normal.model.optim <- function(W.obs.is, X.is, strata.probas.is, method="Nelder-Mead"){
    fitfn <- function(pars){
        bet.arg <- pars[seq(5)]
        sig <- pars[6]
        dens.ls <- (dnorm(W.obs.is, as.vector(X.is %*% bet.arg), sig, log=T) *
                    strata.probas.is)        
        return(sum(dens.ls))        
    }

    res <- optim(par=c(rep(0,5),1), fn=fitfn, method=method, control=list(fnscale=-1, maxit=10000))
    return(list(optim.res=res,
                beta=res$par[seq(5)],
                sigma=res$par[6]))
}

m.step.optim <- function(assign, W.obs, theta, X, O.fn, strata.probas, method="Nelder-Mead"){

    # optimizing the multinomial model
    opt.res <- optimize.multinom.model.optim(X, strata.probas, method)
    theta[["alpha"]][["AT"]] <- opt.res$alphas[["AT"]]
    theta[["alpha"]][["NT"]] <- opt.res$alphas[["NT"]]
    theta[["alpha"]][["C"]] <- opt.res$alphas[["C"]]
    theta[["alpha"]][["D"]] <- opt.res$alphas[["D"]]

    denom <- (exp( X %*% theta[["alpha"]][["AT"]]) +
              exp( X %*% theta[["alpha"]][["NT"]]) +
              exp( X %*% theta[["alpha"]][["C"]]) +
              exp( X %*% theta[["alpha"]][["D"]]))

    theta[["pi"]][["AT"]] <- as.vector(exp( X %*% theta[["alpha"]][["AT"]]) / denom)
    theta[["pi"]][["NT"]] <- as.vector(exp( X %*% theta[["alpha"]][["NT"]]) / denom)
    theta[["pi"]][["C"]] <- as.vector(exp( X %*% theta[["alpha"]][["C"]]) / denom)
    theta[["pi"]][["D"]] <- as.vector(exp( X %*% theta[["alpha"]][["D"]]) / denom)

    
    ## Taking care of the normal model
    ### AT
    indices <- O.fn(1,1)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["AT"]][indices],method)
    theta[["beta"]][["AT"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["AT"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,1)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["AT"]][indices],method)
    theta[["beta"]][["AT"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["AT"]][[0+1]] <- opt.res$sigma

    ### C
    indices <- O.fn(1,1)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["C"]][indices],method)
    theta[["beta"]][["C"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["C"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,0)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["C"]][indices],method)
    theta[["beta"]][["C"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["C"]][[0+1]] <- opt.res$sigma

    ### NT
    indices <- O.fn(1,0)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["NT"]][indices],method)
    theta[["beta"]][["NT"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["NT"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,0)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["NT"]][indices],method)
    theta[["beta"]][["NT"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["NT"]][[0+1]] <- opt.res$sigma

    ### D
    indices <- O.fn(1,0)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["D"]][indices],method)
    theta[["beta"]][["D"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["D"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,1)
    opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,],
                                           strata.probas[["D"]][indices],method)
    theta[["beta"]][["D"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["D"]][[0+1]] <- opt.res$sigma

    return(theta)
}

em.algo.optim.stopping <- function(W.obs, theta.init, assign, X, O.fn, eps=0.001, maxiter=10){
    theta <- theta.init
    res.strata.probas <- NULL
    strata.probas <- NULL
    obs.lik.list <- NULL

    flag=FALSE
    i = 0
    while(!flag){
        i <- i+1
        cat("Iteration:", i, "\n")

        cat("   E-step:\n")
        strata.probas <- update.posterior.strata.proba.full(theta, X, W.obs, O.fn)
        strata.probas <- add.small.noise(strata.probas)

        cat("   M-step:\n")
        theta <- m.step.optim(assign, W.obs, theta, X, O.fn, strata.probas, "BFGS")

        #summary
        est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,1])))
        res.strata.probas <- rbind(res.strata.probas, est.strata.prop)
        obs.lik.i <- obs.lik.full(X, W.obs, theta, O.fn)
        obs.lik.list <- c(obs.lik.list, obs.lik.i)
        
        print(est.strata.prop)
        cat("  LOWER BOUND LIK:\n")
        print(obs.lik.i)

        cat("  DELTA LIK:\n")
        if(i>1){
            print(obs.lik.i - obs.lik.list[i-1])
        }
        
        # stopping
        if(i > maxiter){
            flag=TRUE
        }
        if(i>1){
            if(obs.lik.i - obs.lik.list[i-1] < eps){
                flag=TRUE
            }
        }
    }
    return(list(theta, strata.probas, res.strata.probas, obs.lik.list))
}



em.algo.optim <- function(W.obs, theta.init, assign, X, O.fn, niter=10){
    theta <- theta.init
    res.strata.probas <- NULL
    strata.probas <- NULL
    obs.lik.list <- NULL

    for(i in seq(1, niter)){
        cat("Iteration:", i, "\n")

        cat("   E-step:\n")
        strata.probas <- update.posterior.strata.proba.full(theta, X, W.obs, O.fn)
        strata.probas <- add.small.noise(strata.probas)

        cat("   M-step:\n")
        theta <- m.step.optim(assign, W.obs, theta, X, O.fn, strata.probas, "BFGS")

        #summary
        est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,1])))
        res.strata.probas <- rbind(res.strata.probas, est.strata.prop)
        obs.lik.i <- obs.lik.full(X, W.obs, theta, O.fn)
        obs.lik.list <- c(obs.lik.list, obs.lik.i)
        
        print(est.strata.prop)
        cat("  LOWER BOUND LIK:\n")
        print(obs.lik.i)
    }
    return(list(theta, strata.probas, res.strata.probas, obs.lik.list))
}


# checked 
obs.lik.full <- function(X, W, theta, O.fn){
    res <- 0

    # O(1,1)  =>  G \in { AT, C }
    indices <- O.fn(1,1)
    for(index in indices){
        bet.AT.1 <- theta[["beta"]]$AT[[1+1]]
        sig.AT.1 <- theta[["sigma"]]$AT[[1+1]]
        pi.AT <- theta[["pi"]]$AT[index]
        bet.C.1 <- theta[["beta"]]$C[[1+1]]
        sig.C.1 <- theta[["sigma"]]$C[[1+1]]
        pi.C <- theta[["pi"]]$C[index]

        Wi <- W[index]
        Xi <- X[index,]
        
        res <- res + log(pi.AT * dnorm(Wi, sum(Xi*bet.AT.1), sig.AT.1) +
                         pi.C * dnorm(Wi, sum(Xi*bet.C.1), sig.C.1))
    }

    # O(1,0)  =>   G \in { NT, D }
    indices <- O.fn(1,0)
    for(index in indices){
        bet.NT.1 <- theta[["beta"]]$NT[[1+1]]
        sig.NT.1 <- theta[["sigma"]]$NT[[1+1]]
        pi.NT <- theta[["pi"]]$NT[index]
        bet.D.1 <- theta[["beta"]]$D[[1+1]]
        sig.D.1 <- theta[["sigma"]]$D[[1+1]]
        pi.D <- theta[["pi"]]$D[index]

        Wi <- W[index]
        Xi <- X[index,]
        
        res <- res + log(pi.NT * dnorm(Wi, sum(Xi*bet.NT.1), sig.NT.1) +
                         pi.D * dnorm(Wi, sum(Xi*bet.D.1), sig.D.1))
    }

    # O(0,1)  =>   G \in { AT, D }
    indices <- O.fn(0,1)
    for(index in indices){
        bet.AT.0 <- theta[["beta"]]$AT[[0+1]]
        sig.AT.0 <- theta[["sigma"]]$AT[[0+1]]
        pi.AT <- theta[["pi"]]$AT[index]
        bet.D.0 <- theta[["beta"]]$D[[0+1]]
        sig.D.0 <- theta[["sigma"]]$D[[0+1]]
        pi.D <- theta[["pi"]]$D[index]

        Wi <- W[index]
        Xi <- X[index,]
        
        res <- res + log(pi.AT * dnorm(Wi, sum(Xi*bet.AT.0), sig.AT.0) +
                         pi.D * dnorm(Wi, sum(Xi*bet.D.0), sig.D.0))
    }

    
    # O(0,0)  =>   G \in { NT, C }
    indices <- O.fn(0,0)
    for(index in indices){
        bet.NT.0 <- theta[["beta"]]$NT[[0+1]]
        sig.NT.0 <- theta[["sigma"]]$NT[[0+1]]
        pi.NT <- theta[["pi"]]$NT[index]
        bet.C.0 <- theta[["beta"]]$C[[0+1]]
        sig.C.0 <- theta[["sigma"]]$C[[0+1]]
        pi.C <- theta[["pi"]]$C[index]

        Wi <- W[index]
        Xi <- X[index,]
        
        res <- res + log(pi.NT * dnorm(Wi, sum(Xi*bet.NT.0), sig.NT.0) +
                         pi.C * dnorm(Wi, sum(Xi*bet.C.0), sig.C.0))
    }

    return(res)
}



# # # # #
# Reduced
# # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # #

# utilities
PS <- function(Z.i, T.i){
    O.vals <- paste(as.character(c(Z.i, T.i)), collapse=',')
    res <- switch(O.vals,
                  "1,1"=c("AT", "C"),
                  "1,0"=c("NT", "D"),
                  "0,1"=c("AT", "D"),
                  "0,0"=c("NT", "C"))
    return(res)
}

# returns a mapping between strata, and which O(.,.) they can come from
S2O <- function(stratum){
    res <- switch(stratum,
                  "AT" = list(c(1,1), c(0,1)),
                  "NT" = list(c(1,0), c(0,0)),
                  "C"  = list(c(1,1), c(0,0)),
                  "D"  = list(c(1,0), c(0,1)))
    return(res)
}


# for each stratum, returns whether the model is the same for
# treatment and control (returns 1) or different (returns 2). This
# is useful for setting null hypotheses

MOD <- function(stratum){
    res <- switch(stratum,
                  "AT"=2,
                  "NT"=2,
                  "C"=2,
                  "D"=2)
    return(res)
}


# E-step

extract.current.bet.sig <- function(theta, strata, assign, MOD.fn){
    bet.list <- list()
    sig.list <- list()
    for(stratum in strata){
        if(MOD.fn(stratum) == 1){
            bet.list[[stratum]] <- theta[["beta"]][[stratum]]
            sig.list[[stratum]] <- theta[["sigma"]][[stratum]]
        } else {
            bet.list[[stratum]] <- theta[["beta"]][[stratum]][[assign+1]]
            sig.list[[stratum]] <- theta[["sigma"]][[stratum]][[assign+1]]
        }
    }
    return(list(bet = bet.list, sig = sig.list))
}


update.posterior.strata.proba <- function(theta, X, W.obs, O.fn, CP.fn, PS.fn, MOD.fn){
    init.vals <- rep(0, nrow(X))
    strata.probas <- list("AT" = init.vals,
                          "NT" = init.vals,
                          "C"  = init.vals,
                          "D"  = init.vals)
    O.list <- list(c(1,1), c(1,0), c(0,1), c(0,0))

    for(O.val in O.list){
        # get the indices, and potential strata of units corresponding to that O val
        indices <- do.call(O.fn, as.list(O.val))
        strata <- do.call(PS.fn, as.list(O.val))

        for(idx in indices){
            pis.t.idx <- as.list(sapply(strata, function(stratum) theta[["pi"]][[stratum]][idx])) # get current pi values for idx
            names(pis.t.idx) <- strata
            assign = O.val[1]
            bet.sig <- extract.current.bet.sig(theta, strata, assign, MOD.fn) # get current value for beta and sigma
            bet.list <- bet.sig$bet ; sig.list <- bet.sig$sig
            for(stratum in strata){
                # update the posterior strata proba, at iteration t
                strata.probas[[stratum]][idx] <- update.one.posterior.strata.proba(stratum, strata, pis.t.idx, W.obs[idx], X[idx,],
                                                                                   bet.list, sig.list)
            }            
        }               
    }
    return(strata.probas)
}


# M step

setup.glm <- function(x, N, K){
  y.mod <- matrix(rep(diag(K), N), nrow=K*N, byrow=T)
  x.mod <- matrix(rep(x, each=K), nrow=K*N)
  return(list(ymod=y.mod, xmod=x.mod))
}

add.small.noise <- function(strata.probas, noise = 0.0000000001){
  # adds small noise to the strata probabilities to prevent them from
  # being exactly zero. Zero probabilities seem (I don't why yet) to throw
  # off the multinomial optimization step (vglm raises an error)

  for(i in seq(length(strata.probas[[1]]))){
    for(stratname in names(strata.probas)){
      if(strata.probas[[stratname]][i] == 0){
        strata.probas[[stratname]][i] <- noise 
      }
    }
  }
  return(strata.probas)
}

optimize.multinom.model <- function(X, strata.probas){
  # G and X are not subsetted. We update the entire thing in one go!
  # what we need is to select the right weights in strata.probas.
  W <- as.vector(t(do.call('cbind', strata.probas)))
  K <- length(strata.probas)
  print("setting up glm")
  mod.input <- setup.glm(X, nrow(X), K)
  print("starting regression")

  #if(!("coefstart" %in% ls(envir=.GlobalEnv))){
  #    coefstart <- NULL
  #}
  res <- vglm(mod.input$ymod~mod.input$xmod - 1,
              family=multinomial,
              weights=W,
              trace=TRUE, epsilon=0.0001)

  coefs <- res@coefficients
  
  #assign("coefstart", coef(res), envir=.GlobalEnv) # ugly as hell, but just to try something
  
  alphas <- list()
  idx <- 0
  for(stratum in names(strata.probas)[-K]){
    idx <- idx + 1
    slice <- rep(0, K-1)
    slice[idx] <- 1
    alphas[[stratum]] <- coefs[as.logical(rep(slice, ncol(X)))]
  }
  alphas[[names(strata.probas)[K]]] <- rep(0,ncol(X)) 
  
  return(alphas) # alphas is a list of lists. Beware the order.
}


optimize.normal.model <- function(W.obs.is, X.is, strata.probas.is){
  # W.obs.is and X.is have already been subseted, as well as the
  # weights, strata.probas.is. Note that W.obs.is denotes the observed outcome
  res <- lm(W.obs.is ~ X.is - 1, weights = strata.probas.is)
  res.beta <- res$coefficients
  res.sigma <- summary(res)$sigma
  return(list("beta" = res.beta, "sigma" = res.sigma))
}


m.step <- function(assign, W.obs, theta, X, PS.fn, MOD.fn, O.fn, S2O.fn, strata.probas){
    O.list <- list(c(1,1), c(1,0), c(0,1), c(0,0))

    # we first optimize the multinomial model and update the alphas
    opt.res <- optimize.multinom.model(X, strata.probas) 
    denom <- 0
    for(stratum in c("AT", "C", "NT", "D")){
        theta[["alpha"]][[stratum]] <- opt.res[[stratum]]
        denom <- denom + exp(X %*% opt.res[[stratum]])
    }
    for(stratum in c("AT", "C", "NT", "D")){
        theta[["pi"]][[stratum]] <- as.vector(exp(X %*% opt.res[[stratum]]) / denom)
    }

    # we then take care of normal model
    for(stratum in c("AT", "C", "NT", "D")){
        indices <- NULL
        for(oval in S2O.fn(stratum)){
            indices <- c(indices, do.call(O.fn, as.list(oval)))
        }
        # dealing with control first
        indices.0 <- intersect(indices, which(assign==0))
        opt.res <- optimize.normal.model(W.obs[indices.0], X[indices.0,], strata.probas[[stratum]][indices.0])
        theta[["beta"]][[stratum]][[0+1]] <- opt.res$beta
        theta[["sigma"]][[stratum]][[0+1]] <- opt.res$sigma

        # dealing with treatment
        indices.1 <- intersect(indices, which(assign==1))
        opt.res <- optimize.normal.model(W.obs[indices.1], X[indices.1,], strata.probas[[stratum]][indices.1])
        theta[["beta"]][[stratum]][[1+1]] <- opt.res$beta
        theta[["sigma"]][[stratum]][[1+1]] <- opt.res$sigma
    }
    return(theta)
}


m.step.mixte <- function(assign, W.obs, theta, X, O.fn, strata.probas, method="Nelder-Mead"){

    # optimizing the multinomial model
    opt.res <- optimize.multinom.model.optim(X, strata.probas, method)
    theta[["alpha"]][["AT"]] <- opt.res$alphas[["AT"]]
    theta[["alpha"]][["NT"]] <- opt.res$alphas[["NT"]]
    theta[["alpha"]][["C"]] <- opt.res$alphas[["C"]]
    theta[["alpha"]][["D"]] <- opt.res$alphas[["D"]]

    denom <- (exp( X %*% theta[["alpha"]][["AT"]]) +
              exp( X %*% theta[["alpha"]][["NT"]]) +
              exp( X %*% theta[["alpha"]][["C"]]) +
              exp( X %*% theta[["alpha"]][["D"]]))

    theta[["pi"]][["AT"]] <- as.vector(exp( X %*% theta[["alpha"]][["AT"]]) / denom)
    theta[["pi"]][["NT"]] <- as.vector(exp( X %*% theta[["alpha"]][["NT"]]) / denom)
    theta[["pi"]][["C"]] <- as.vector(exp( X %*% theta[["alpha"]][["C"]]) / denom)
    theta[["pi"]][["D"]] <- as.vector(exp( X %*% theta[["alpha"]][["D"]]) / denom)

    
    ## Taking care of the normal model
    ### AT
    indices <- O.fn(1,1)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["AT"]][indices])
    theta[["beta"]][["AT"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["AT"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,1)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["AT"]][indices])
    theta[["beta"]][["AT"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["AT"]][[0+1]] <- opt.res$sigma

    ### C
    indices <- O.fn(1,1)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["C"]][indices])
    theta[["beta"]][["C"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["C"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,0)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["C"]][indices])
    theta[["beta"]][["C"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["C"]][[0+1]] <- opt.res$sigma

    ### NT
    indices <- O.fn(1,0)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["NT"]][indices])
    theta[["beta"]][["NT"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["NT"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,0)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                           strata.probas[["NT"]][indices])
    theta[["beta"]][["NT"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["NT"]][[0+1]] <- opt.res$sigma

    ### D
    indices <- O.fn(1,0)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["D"]][indices])
    theta[["beta"]][["D"]][[1+1]] <- opt.res$beta
    theta[["sigma"]][["D"]][[1+1]] <- opt.res$sigma

    indices <- O.fn(0,1)
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,],
                                     strata.probas[["D"]][indices])
    theta[["beta"]][["D"]][[0+1]] <- opt.res$beta
    theta[["sigma"]][["D"]][[0+1]] <- opt.res$sigma

    return(theta)
}


em.algo.for <- function(W.obs, theta.init, assign, X, ps, mod, O.fn, cp, S2O.fn, niter=10){
  theta <- theta.init
  res.strata.probas <- NULL
  strata.probas <- NULL #strata.probas.init
  lower.bound.lik.list <- NULL
  
  for(i in seq(1, niter)){
      cat("Iteration:", i, "\n")
      
      cat("   E-step:\n")
      strata.probas <- update.posterior.strata.proba(theta, X, W.obs, O.fn, cp, ps, mod)
      strata.probas <- add.small.noise(strata.probas)
      
      cat("   M-step:\n")
      theta <- m.step(assign, W.obs, theta, X, ps, mod, O.fn, S2O.fn, strata.probas)
      
      est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,1]))) # all obs equally weighted
      print(est.strata.prop)
      res.strata.probas <- rbind(res.strata.probas, est.strata.prop)
      
      #lower.bound.i <- lower.bound.lik(X, W.obs, theta, strata.probas, ps, O.fn)
      lower.bound.i <- obs.lik.full(X, W.obs, theta, O.fn)#obs.lik(X, W.obs, theta, ps, O.fn)
      lower.bound.lik.list <- c(lower.bound.lik.list, lower.bound.i)
      cat("  LOWER BOUND LIK:\n")
      print(lower.bound.i)
  }
  return(list(theta, strata.probas, res.strata.probas, lower.bound.lik.list))
}



# likelihood
    
obs.lik.contribution <- function(index, assign, strata, X, W, theta){
    res <- 0
    Xi <- X[index,]
    for(stratum in strata){
        bet <- theta[["beta"]][[stratum]][[1+assign]]
        sig <- theta[["sigma"]][[stratum]][[1+assign]]
        res <- res + theta[["pi"]][[stratum]][index] * dnorm(W[index], sum(X[index,] * bet), sig)
    }
    return(log(res))
}

obs.lik <- function(X, W, theta, PS.fn, O.fn){
    O.list <- list(c(1,1), c(0,1), c(1,0), c(0,0))
    res <- 0
    for(O.val in O.list){
        indices <- do.call(O.fn, as.list(O.val))
        strata <- do.call(PS.fn, as.list(O.val))
        #print(strata)
        for(index in indices){
            Zi <- O.val[1]
            res <- res + obs.lik.contribution(index, Zi, strata, X, W, theta)
        }
    }
    return(res)
}


# # # # # #
# Direct optimization
# # # # # #

#direct.optim <- function(X, W, theta.init, O.fn){
#    fitfn <- function(pars){
#        bet <- pars[seq(1, 40)]
#        sig <- pars[seq(41, 48)]
#        alpha <- pars[seq(49, 63)]
#    }
#
#}
