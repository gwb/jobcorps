

# # # # # # #
# Utilities #
# # # # # # #

get.O <- function(y){
    O <- function(Z.i, T.i){
        return(which(y$assign==Z.i & y$treat==T.i))
    }
    return(O)
}


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
# E-step #
# # # # #

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



# # # # # 
# M-step #
# # # # #

optimize.multinom.model.optim <- function(X, strata.probas, method="Nelder-Mead"){
    K <- ncol(X)
    fitfn <- function(alph.arg){
        alpha.AT <- c(alph.arg[seq(1,K)])
        alpha.NT <- c(alph.arg[seq(K+1,2*K)])
        alpha.C <- c(alph.arg[seq(2*K+1,3*K)])
        alpha.D <- rep(0, K)
                    
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

    res <- optim(par=rep(0,3*K), fn=fitfn, method=method, control=list(fnscale=-1, maxit=10000))
    return(list(optim.res = res,
                alphas =
                list(AT = res$par[seq(1,K)],
                     NT = res$par[seq(K+1,2*K)],
                     C  = res$par[seq(2*K+1,3*K)],
                     D  = rep(0, K))))
}

optimize.normal.model.optim <- function(W.obs.is, X.is, strata.probas.is, method="Nelder-Mead"){
    K <- ncol(X.is)
    fitfn <- function(pars){
        bet.arg <- pars[seq(K)]
        sig <- pars[K+1]
        dens.ls <- (dnorm(W.obs.is, as.vector(X.is %*% bet.arg), sig, log=T) *
                    strata.probas.is)        
        return(sum(dens.ls))        
    }

    res <- optim(par=c(rep(0,K),1), fn=fitfn, method=method, control=list(fnscale=-1, maxit=10000))
    return(list(optim.res=res,
                beta=res$par[seq(K)],
                sigma=res$par[K+1]))
}


m.step.optim.old <- function(assign, W.obs, theta, X, PS.fn, MOD.fn, O.fn, S2O.fn, strata.probas, method="BFGS"){
    O.list <- list(c(1,1), c(1,0), c(0,1), c(0,0))

    # we first optimize the multinomial model and update the alphas
    opt.res <- optimize.multinom.model.optim(X, strata.probas, method=method)$alphas 
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
        opt.res <- optimize.normal.model.optim(W.obs[indices.0], X[indices.0,], strata.probas[[stratum]][indices.0], method)
        theta[["beta"]][[stratum]][[0+1]] <- opt.res$beta
        theta[["sigma"]][[stratum]][[0+1]] <- opt.res$sigma

        # dealing with treatment
        indices.1 <- intersect(indices, which(assign==1))
        opt.res <- optimize.normal.model.optim(W.obs[indices.1], X[indices.1,], strata.probas[[stratum]][indices.1], method)
        theta[["beta"]][[stratum]][[1+1]] <- opt.res$beta
        theta[["sigma"]][[stratum]][[1+1]] <- opt.res$sigma
    }
    return(theta)
}


m.step.optim <- function(assign, W.obs, theta, X, PS.fn, MOD.fn, O.fn, S2O.fn, strata.probas, method="BFGS"){
    O.list <- list(c(1,1), c(1,0), c(0,1), c(0,0))

    # we first optimize the multinomial model and update the alphas
    opt.res <- optimize.multinom.model.optim(X, strata.probas, method=method)$alphas 
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
        if(MOD.fn(stratum) == 2){
            # dealing with control first
            indices.0 <- intersect(indices, which(assign==0))
            opt.res <- optimize.normal.model.optim(W.obs[indices.0], X[indices.0,], strata.probas[[stratum]][indices.0], method)
            theta[["beta"]][[stratum]][[0+1]] <- opt.res$beta
            theta[["sigma"]][[stratum]][[0+1]] <- opt.res$sigma
            
            # dealing with treatment
            indices.1 <- intersect(indices, which(assign==1))
            opt.res <- optimize.normal.model.optim(W.obs[indices.1], X[indices.1,], strata.probas[[stratum]][indices.1], method)
            theta[["beta"]][[stratum]][[1+1]] <- opt.res$beta
            theta[["sigma"]][[stratum]][[1+1]] <- opt.res$sigma
        } else {
            opt.res <- optimize.normal.model.optim(W.obs[indices], X[indices,], strata.probas[[stratum]][indices], method)
            theta[["beta"]][[stratum]] <- opt.res$beta
            theta[["sigma"]][[stratum]] <- opt.res$sigma            
        }
    }
    return(theta)
}

# EM

em.algo.optim.stopping.multistart <- function(W.obs, y, X, O.fn, CP.fn, PS.fn, S2O.fn, MOD.fn, n.restart=5, eps=0.001, maxiter=10, method="BFGS"){
    res <- vector('list', n.restart)
    for(i in seq(n.restart)){
        N <- length(W.obs)
        theta.init <- initialize.params(W.obs, y, X)
        #names(strata.probas.list) <- c("AT", "NT", "C", "D")
        theta.init[["pi"]] <- list(AT=rep(0.25, N),
                                   NT=rep(0.25, N),
                                   C=rep(0.25,N),
                                   D=rep(0.25,N))
        em.res <- em.algo.optim.stopping(W.obs, theta.init, y$assign, X, O.fn, CP.fn, PS.fn, S2O.fn, MOD.fn, maxiter=maxiter, eps=eps, method=method)
        res[[i]] <- em.res
    }
    mls <- do.call('c',lapply(res, function(x) tail(x[[4]],1)))
    return(res[[which.max(mls)]])
}

em.algo.optim.stopping <- function(W.obs, theta.init, assign, X, O.fn, CP.fn, PS.fn, S2O.fn, MOD.fn, eps=0.001, maxiter=10,
                                   method="BFGS"){
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
        strata.probas <- update.posterior.strata.proba(theta, X, W.obs, O.fn, CP.fn, PS.fn, MOD.fn)
        #strata.probas <- add.small.noise(strata.probas)

        cat("   M-step:\n")
        theta <- m.step.optim(assign, W.obs, theta, X, PS.fn, MOD.fn, O.fn, S2O.fn, strata.probas, method)

        #summary
        est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, rep(1, nrow(X)))))
        res.strata.probas <- rbind(res.strata.probas, est.strata.prop)
        obs.lik.i <- obs.lik(X, W.obs, theta, PS.fn, O.fn, MOD.fn)
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


# single covariate optim

lik.optim <- function(W, PS.fn, O.fn, method="BFGS"){
    X <- matrix(1, nrow=length(W), ncol=1)
    fitfn <- function(pars){
        beta <- list(AT = list(pars[1], pars[2]),
                     NT = list(pars[3], pars[4]),
                     C = list(pars[5], pars[6]),
                     E = list(pars[7], pars[8]))
        alpha <- list(AT=pars[9],
                      NT=pars[10],
                      C=pars[11],
                      D=pars[12])
        sigma <- list(AT=list(exp(pars[13]), exp(pars[14])),
                      NT=list(exp(pars[15]), exp(pars[16])),
                      C=list(exp(pars[17]), exp(pars[18])),
                      D=list(exp(pars[19]), exp(pars[20])))
        theta <- list(alpha = alpha,
                      beta=beta,
                      sigma=sigma)
        alph <- do.call('cbind', alpha)
        lin.comp <- X %*% alph
        strata.probas <- exp(lin.comp) / rowSums(exp(lin.comp))

        theta$pi <- lapply(seq_len(ncol(strata.probas)), function(i) strata.probas[,i])
        names(theta$pi) <- c("AT", "NT", "C", "D")
        return(obs.lik(X, W, theta, PS.fn, O.fn))
    }
    return(optim(par=runif(20), fitfn, method=method, control=list(fnscale=-1, trace=3)))
}


# likelihood
obs.lik.contribution.old <- function(index, assign, strata, X, W, theta){
    res <- 0
    Xi <- X[index,]
    for(stratum in strata){
        bet <- theta[["beta"]][[stratum]][[1+assign]]
        sig <- theta[["sigma"]][[stratum]][[1+assign]]
        res <- res + theta[["pi"]][[stratum]][index] * dnorm(W[index], sum(X[index,] * bet), sig)
    }
    return(log(res))
}

obs.lik.contribution <- function(index, assign, strata, X, W, theta, MOD.fn=NULL){
    if(is.null(MOD.fn)){
        MOD.fn <- function(x) return(2)
    }
    res <- 0
    Xi <- X[index,]
    for(stratum in strata){
        if(MOD.fn(stratum)==2){
            bet <- theta[["beta"]][[stratum]][[1+assign]]
            sig <- theta[["sigma"]][[stratum]][[1+assign]]
        }else{
            bet <- theta[["beta"]][[stratum]]
            sig <- theta[["sigma"]][[stratum]]
        }
        res <- res + theta[["pi"]][[stratum]][index] * dnorm(W[index], sum(X[index,] * bet), sig)
    }
    return(log(res))
}



obs.lik <- function(X, W, theta, PS.fn, O.fn, MOD.fn){
    O.list <- list(c(1,1), c(0,1), c(1,0), c(0,0))
    res <- 0
    for(O.val in O.list){
        indices <- do.call(O.fn, as.list(O.val))
        strata <- do.call(PS.fn, as.list(O.val))
        #print(strata)
        for(index in indices){
            Zi <- O.val[1]
            res <- res + obs.lik.contribution(index, Zi, strata, X, W, theta, MOD.fn)
        }
    }
    return(res)
}

