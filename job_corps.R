require(VGAM)

# Data structures

# theta - nested list of the following form:
# theta:
#       alpha:
#             cEE: vector
#             cEN: vector
#             cNE: vector
#             cNN: vector
#             nEE: vector
#             nNN: vector
#       beta:
#             cEE:
#                 1: vector
#                 2: vector
#             cEN: vector
#             cNE: vector
#             nEE: vector
#       sigma:
#             cEE:
#                 1: float
#                 2: float
#             cEN: float
#             cNE: float
#             nEE: float
#       pi:
#             cEE: vector (of length N=13987)
#             cEN: vector (of length N=13987)
#             cNE: vector (of length N=13987)
#             cNN: vector (of length N=13987)
#             nEE: vector (of length N=13987)
#             nNN: vector (of length N=13987)


# Globals

# returns the mapping betwen values of "O" and principal strata
PS <- function(Z.i, D1.i, S.i){
  
  O.vals <- paste(as.character(c(Z.i,D1.i,S.i)), collapse=",")

  res <- switch(O.vals,
                "1,1,1"=c("cEE","cEN"),
                "1,1,0"=c("cNN","cNE"),
                "1,0,1"=c("nEE"),
                "1,0,0"=c("nNN"),
                "1,NA,1"=c("nEE","cEE","cEN"),
                "1,NA,0"=c("nNN","cNN","cNE"),
                "0,NA,1"=c("cEE","cNE","nEE"),
                "0,NA,0"=c("cEN","cNN","nNN"),
                "1,1,NA"=c("cEE","cEN","cNE","cNN"),
                "1,0,NA"=c("nEE","nNN"),
                "0,NA,NA"=c("cEE","cEN","cNE","cNN","nEE","nNN")
                )
  return(res)
}

# returns a mapping between strata, and which O(.,.,.) they can come from
S2O <- function(stratum){
  res <- switch(stratum,
                "cEE" = list(c(1,1,1), c(1,NA,1), c(0,NA,1), c(1,1,NA), c(0,NA,NA)),
                "cEN" = list(c(1,1,1), c(1,NA,1), c(0,NA,0), c(1,1,NA), c(0,NA,NA)),
                "cNE" = list(c(1,1,0), c(1,NA,0), c(0,NA,1), c(1,1,NA), c(0,NA,NA)),
                "cNN" = list(c(1,1,0), c(1,NA,0), c(0,NA,0), c(1,1,NA), c(0,NA,NA)),
                "nEE" = list(c(1,0,1), c(1,NA,1), c(1,0,NA), c(0,NA,1)),
                "nNN" = list(c(1,0,0), c(1,NA,0), c(0,NA,0), c(1,0,NA), c(0,NA,NA)))
  return(res)
}

# returns a mapping like above, but returns only the O(.,.,.) corresponding to normal optimization
S2O.normal <- function(stratum){
  res <- switch(stratum,
                "cEE1" = list(c(1,1,1), c(1,NA,1)),
                "cEE0" = list(c(0,NA,1)),
                "cEN" = list(c(1,1,1), c(1,NA,1)),
                "cNE" = list(c(0,NA,1)),
                "cNN" = list(),
                "nEE" = list(c(1,0,1), c(1,NA,1), c(0,NA,1)),
                "nNN" = list())
  return(res)
}

# for each stratum, returns whether the model is the same
# for treatment and control (returns 1) or different (returns 2)

MOD <- function(stratum){
  res <- switch(stratum,
                "cEE"=2,
                "cEN"=1,
                "cNE"=1,
                "nEE"=1)
  return(res)
}


# for each possible o value, wether the proba of belonging to specific latent strata
# should be updated using regular formula (including normal mixtures) or short formula,
# in the E step.
CP <- function(Z.i, D1.i, S.i){

  O.vals <- paste(as.character(c(Z.i,D1.i,S.i)), collapse=",")

  res <- switch(O.vals,
                "1,1,1"="regular",
                "1,1,0"="short",
                "1,0,1"="short",
                "1,0,0"="short",
                "1,NA,1"="regular",
                "1,NA,0"="short",
                "0,NA,1"="regular",
                "0,NA,0"="short",
                "1,1,NA"="short",
                "1,0,NA"="short"
                #"0,NA,NA"=c("cEE","cEN","cNE","cNN","nEE","nNN")
                )
  return(res)
}


# Returns a closure that indicates which units belong to a certain O category
get.O <- function(y, post=52){
  O <- function(Z.i, D1.i, S.i){

    # if Z.i = 0, then D1.i = NA corresponds to treat = 0, so we make
    # a slight modification
    if(Z.i == 0){
      if(is.na(D1.i)){
        D1.i <- 0
      }
    }
    
    assign.p <- y$assign == Z.i

    if(is.na(D1.i)){ treat.p <- is.na(y$treat) }
    else { treat.p <- y$treat == D1.i} 

    if(is.na(S.i)){
      if(post==52){ work.p <- is.na(y$work52) }
      else if(post==130){ work.p <- is.na(y$work130)}
      else if(post==208){ work.p <- is.na(y$work208)}
    } else {
      if(post==52){ work.p <- y$work52 == S.i }
      else if(post==130){ work.p <- y$work130 == S.i }
      else if(post==208){ work.p <- y$work208 == S.i}
    }
    
    raw.res <- seq(nrow(y))[assign.p & treat.p & work.p]
    return(raw.res[!is.na(raw.res)])
  }
  return(O)
}


# likelihood helper functions

lower.bound.contribution.normal.all <- function(indices, stratum, X, W, theta, strata.probas){
    res <- 0
    for(i in indices){
        Xi <- X[i,]
        bet <- theta[["beta"]][[stratum]]
        sig <- theta[["sigma"]][[stratum]]
        res <- res + strata.probas[[stratum]][i] * (log(theta[["pi"]][[stratum]][i]) +
                                                    dnorm(log(W[i]), sum(Xi * bet), sig, log=T))   
    }
    return(res)
}

lower.bound.contribution.normal.cEE <- function(indices, treat, X, W, theta, strata.probas){
    res <- 0
    for(i in indices){
        #if(i == 4138){
        #    browser()
        #}
        Xi <- X[i,]
        bet <- theta[["beta"]][["cEE"]][[1+treat]]
        sig <- theta[["sigma"]][["cEE"]][[1+treat]]
        res <- res + strata.probas[["cEE"]][i] * (log(theta[["pi"]][["cEE"]][i]) +
                                                    dnorm(log(W[i]), sum(Xi * bet), sig, log=T))
        #print(paste(i, "-", res))
        #if(is.infinite(res)){
        #    browser()
        #}
    }
    return(res)
}

lower.bound.contribution.prob <- function(indices, stratum, theta, strata.probas){
    res <- 0
    for(i in indices){
        res <- res + strata.probas[[stratum]][i] * log(theta[["pi"]][[stratum]][i])
    }
    return(res)
}

lower.bound.lik <- function(X, W, theta, strata.probas, cp, ps, O.fn){
      O.list <- list(c(1,1,1),
                     c(1,1,0),
                     c(1,0,1),
                     c(1,0,0),
                     c(1,NA,1),
                     c(1,NA,0),
                     c(0,NA,1),
                     c(0,NA,0),
                     c(1,1,NA),
                     c(1,0,NA))
      res <- 0
      for(O.val in O.list){
          indices <- do.call(O.fn, as.list(O.val)) # we get the indices of the units corresponding to that group
          strata <- do.call(ps, as.list(O.val))
          for(stratum in strata){
              if(do.call(cp, as.list(O.val)) == "regular"){
                  Zi <- O.val[1]
                  if(stratum=="cEE"){
                      res.tmp <- lower.bound.contribution.normal.cEE(indices, Zi, X, W, theta, strata.probas)
                  } else {
                      res.tmp <- lower.bound.contribution.normal.all(indices, stratum, X, W, theta, strata.probas)
                  }                  
              } else {
                  res.tmp <- lower.bound.contribution.prob(indices, stratum, theta, strata.probas)
              }
              res <- res + res.tmp
          }
      }
      return(res)
}

log.lik.contribution <- function(indices, strata, Z, G, W.obs, mod.fn, X, theta, S){
  log.lik <- 0

  # TODO: rewrite the Z[idx] part, since this should be the same for all the indices
  #       (since indices = O(zi, di(1), si) )
  
  for(idx in indices){
    for(stratum in strata){

      # extract pis parameters
      pis <- theta[["pi"]][[stratum]]

      # Adds contribution of (stratum x idx) to log likelihood
      if(G[idx] == stratum){
        if(!is.na(S) && S == 1){
          # if S = 1 (the unit is salaried) then we look at the wage contribution

          # extract parameters this has to be done only if S=1 (bet do not exist if there is no salary)
          if(mod.fn(stratum) == 1){
            bet <- theta[["beta"]][[stratum]]
            sig <- theta[["sigma"]][[stratum]]
          }else{
            bet <- theta[["beta"]][[stratum]][[Z[idx]+1]]
            sig <- theta[["sigma"]][[stratum]][[Z[idx]+1]]
          }

          # adds contribution
          log.lik <- log.lik + log(pis[idx]) + dnorm(log(W.obs[idx]), mean= sum(X[idx,] * bet), sd=sig, log=T)
          
        } else {
          # if S = 0 then it is simply a pi contribution
          log.lik <- log.lik + log(pis[idx])
        }
      }
    }
    
  }
  return(log.lik)
}


complete.data.loglik <- function(theta, W.obs, G, Z, X, ps, mod.fn, O.fn){
  # parameters:
  #    - mod: a function returning an integer (1 or 2) indicating wether the wage model
  #           is the same in control and treatment (1) or different (2). If it is different,
  #           then we have {sigma,beta}->stratum->0or1, else, we have {sigma,beta}->stratum
  #      mod stands for *model*
  #    - ps: a function taking Zi, D1i, Si as input, and returning the list of possible principle strata
  #          to which it belongs (mixture).
  #    - O.fn: takes Z, D1, S as input and returns a vector of all indices such that (Zi,D1i,Si) == (Z,D1,S)
  #    - G: a vector linking a specific index to its stratum
  #    - Wobs: observed wage
  #    - X: covariates
  #    - Z: indicator for assignment to control (Z=0) or treatment (Z=1)
  
  N <- ncol(Z)
  O.list <- list(c(1,1,1),
                 c(1,1,0),
                 c(1,0,1),
                 c(1,0,0),
                 c(1,NA,1),
                 c(1,NA,0),
                 c(0,NA,1),
                 c(0,NA,0),
                 c(1,1,NA),
                 c(1,0,NA)) # corresponds to (Z,D,S)
  log.lik <- 0
  for(O.val in O.list){
    print(O.val)
    log.lik <- log.lik + log.lik.contribution(O.fn(O.val[[1]], O.val[[2]], O.val[[3]]),
                                              ps(O.val[[1]], O.val[[2]], O.val[[3]]),
                                              Z, G, W.obs, mod.fn, X, theta, O.val[[3]])
  }
  return(log.lik)
}

# # # # # # # # # #
# helpers for E-step
# # # # # # # # # #

update.one.regular.strata.proba.old <- function(stratum, strata, pis.t.idx, Wi, Xi, bet.list, sig.list){
  numerator <- pis.t.idx[[stratum]] * dnorm(log(Wi), mean= sum(Xi * bet.list[[stratum]]), sd=sig.list[[stratum]])
  denom <- 0
  for(alt.stratum in strata){
    denom <- denom + pis.t.idx[[alt.stratum]] * dnorm(log(Wi),
                                                      mean= sum(Xi * bet.list[[alt.stratum]]),
                                                      sd=sig.list[[alt.stratum]])

  }

  return(numerator/denom)
}

norm.ratio <- function(Wi, bet1, sig1, bet2, sig2, Xi){
  return( sig2/sig1 * exp( -(log(Wi) - sum(Xi * bet1))^2/(2*sig1^2) +
                          (log(Wi) - sum(Xi * bet2))^2/(2*sig2^2)))
}

update.one.regular.strata.proba <- function(stratum, strata, pis.t.idx, Wi, Xi, bet.list, sig.list){
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

update.one.short.strata.proba <- function(stratum, strata, pis.t.idx){
  numerator <- pis.t.idx[[stratum]]
  denom <- sum(sapply(strata, function(alt.stratum) pis.t.idx[[alt.stratum]]))
  return(numerator/denom)
}

# Computes the conditionals, e.g P(G_i = "cEE" | theta)
update.strata.proba <- function(theta, X, W.obs, O.fn, cp, ps, mod){
  init.vals <- rep(0, nrow(X))
  strata.probas <- list("cEE"=init.vals,
                        "cEN"=init.vals,
                        "cNE"=init.vals,
                        "cNN"=init.vals,
                        "nEE"=init.vals,
                        "nNN"=init.vals)
  O.list <- list(c(1,1,1),
                 c(1,1,0),
                 c(1,0,1),
                 c(1,0,0),
                 c(1,NA,1),
                 c(1,NA,0),
                 c(0,NA,1),
                 c(0,NA,0),
                 c(1,1,NA),
                 c(1,0,NA))
  
  for(O.val in O.list){
    indices <- do.call(O.fn, as.list(O.val)) # we get the indices of the units corresponding to that group
    strata <- do.call(ps, as.list(O.val)) 
    #print(O.val)
    #print(length(indices))
    for(idx in indices){

      # we first get the old pis for idx, for every stratum valid for O.val
      pis.t.idx <- as.list(sapply(strata, function(stratum) theta[["pi"]][[stratum]][idx]))
      names(pis.t.idx) <- strata

      # there are two types of conditionals:
      # -> "regular" involving pis_t and N(log(W) | XB, sig)
      # -> "short" involving only the pis_t
      # we treat them differently here
      if(do.call(cp, as.list(O.val))== "regular"){

        # we compute all the betas and sigmas for this idx, for all strata
        bet.list <- list()
        sig.list <- list()
        for(stratum in strata){
          # there are two cases:
          # -> either the distribution of wages are the same params in control and treat (mod ==1)
          # -> or they the distributions are different (mod==2)
          # we treat the two differently
          if(mod(stratum) == 1){
            bet.list[[stratum]] <- theta[["beta"]][[stratum]]
            sig.list[[stratum]] <- theta[["sigma"]][[stratum]]
          } else {
            bet.list[[stratum]] <- theta[["beta"]][[stratum]][[O.val[1]+1]]
            sig.list[[stratum]] <- theta[["sigma"]][[stratum]][[O.val[1]+1]]
          }
        }

        # we now compute the new pis for this idx, and for each strata (regular case)
        for(stratum in strata){
          #print("update regular")
          #  if(idx == 1){browser()}
          strata.probas[[stratum]][idx] <- update.one.regular.strata.proba(stratum,
                                                                           strata,
                                                                           pis.t.idx,
                                                                           W.obs[idx],
                                                                           X[idx,],
                                                                           bet.list,
                                                                           sig.list)
        }
        
      } else {
        
        # we compute the new pis (short case)
        for(stratum in strata){
          #print("update short")
          strata.probas[[stratum]][idx] <- update.one.short.strata.proba(stratum, strata, pis.t.idx)
        }
        
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

compare.ovalues <- function(oval1, oval2){
  for(i in seq(3)){
    if(! is.na(oval1)[i] & ! is.na(oval2)[i] ){
      if(oval1[i] != oval2[i]){
        return(FALSE)
      }
    }
    
    if( is.na(oval1)[i] & ! is.na(oval2)[i] ){
      return(FALSE)
    }
    
    if( !is.na(oval1)[i] & is.na(oval2)[i] ){
      return(FALSE)
    }
  }
  return(TRUE)
}

check.optimize.normal <- function(oval){

  return(compare.ovalues(oval, c(1,1,1)) |
         compare.ovalues(oval, c(1,0,1)) |
         compare.ovalues(oval, c(1,NA,1))|
         compare.ovalues(oval, c(0,NA,1)))

}

optimize.normal.model <- function(W.obs.is, X.is, strata.probas.is){
  # W.obs.is and X.is have already been subseted, as well as the
  # weights, strata.probas.is. Note that W.obs.is denotes the observed wages.
  res <- lm(log(W.obs.is) ~ X.is - 1, weights = strata.probas.is)
  res.beta <- res$coefficients
  res.sigma <- summary(res)$sigma
  return(list("beta" = res.beta, "sigma" = res.sigma))
}

optimize.multinom.model <- function(X, strata.probas){
  # G and X are not subsetted. We update the entire thing in one go!
  # what we need is to select the right weights in strata.probas.
  W <- as.vector(t(do.call('cbind', strata.probas)))
  K <- length(strata.probas)
  print("setting up glm")
  mod.input <- setup.glm(X, nrow(X), K)
  print("starting regression")
  res <- vglm(mod.input$ymod~mod.input$xmod - 1,
              family=multinomial,
              weights=W,
              trace=TRUE, epsilon=0.001)
  
  coefs <- res@coefficients
  
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

m.step <- function(W.obs, theta, Z, X, ps, mod, O.fn, strata.probas){
  N <- ncol(Z)
  O.list <- list(c(1,1,1),
                 c(1,1,0),
                 c(1,0,1),
                 c(1,0,0),
                 c(1,NA,1),
                 c(1,NA,0),
                 c(0,NA,1),
                 c(0,NA,0),
                 c(1,1,NA),
                 c(1,0,NA))

  # we first optimize the multinomial, and update the \alphas
  opt.res <- optimize.multinom.model(X, strata.probas)
 
  denom <- 0
  for(stratum in c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")){
    theta[["alpha"]][[stratum]] <- opt.res[[stratum]]
    denom <- denom + exp(X %*% opt.res[[stratum]])
  }
  for(stratum in c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")){
    theta[["pi"]][[stratum]] <- as.vector(exp(X %*% opt.res[[stratum]]) / denom)
  }


  for(stratum in c("cEN", "cNE", "nEE")){
    print(paste("processing: ", stratum))
    indices <- c()

    for(oval in S2O.normal(stratum)){
      indices <- c(indices, do.call(O.fn, as.list(oval)))
    }

    # we run the weighted OLS routine
    opt.res <- optimize.normal.model(W.obs[indices], X[indices,], strata.probas[[stratum]][indices])

    # => update just B_{g} and sig_{g}
    theta[["beta"]][[stratum]] <- opt.res$beta
    theta[["sigma"]][[stratum]] <- opt.res$sigma
  }

  # deals separately with cEE

  # cEE1
  stratum="cEE"

  print(paste("processing: ", stratum))
  
  indices <- c(O.fn(1,1,1), O.fn(1,NA,1))
  opt.res <- optimize.normal.model(W.obs[indices], X[indices,], strata.probas[[stratum]][indices])
  theta[["beta"]][[stratum]][[1+1]] <- opt.res$beta
  theta[["sigma"]][[stratum]][[1+1]] <- opt.res$sigma

  # cEE0
  indices <- O.fn(0,NA,1)
  opt.res <- optimize.normal.model(W.obs[indices], X[indices,], strata.probas[[stratum]][indices])
  theta[["beta"]][[stratum]][[0+1]] <- opt.res$beta
  theta[["sigma"]][[stratum]][[0+1]] <- opt.res$sigma

  return(theta) # returns the updated theta
}


m.step.old <- function(W.obs, theta, Z, X, ps, mod, O.fn, strata.probas){
  N <- ncol(Z)
  O.list <- list(c(1,1,1),
                 c(1,1,0),
                 c(1,0,1),
                 c(1,0,0),
                 c(1,NA,1),
                 c(1,NA,0),
                 c(0,NA,1),
                 c(0,NA,0),
                 c(1,1,NA),
                 c(1,0,NA))

  # we first optimize the multinomial, and update the \alphas
  opt.res <- optimize.multinom.model(X, strata.probas)
 
  denom <- 0
  for(stratum in c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")){
    theta[["alpha"]][[stratum]] <- opt.res[[stratum]]
    denom <- denom + exp(X %*% opt.res[[stratum]])
  }
  for(stratum in c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")){
    theta[["pi"]][[stratum]] <- as.vector(exp(X %*% opt.res[[stratum]]) / denom)
  }

  
  for(O.val in O.list){
    #print(O.val)
    indices <- do.call(O.fn, as.list(O.val))
    strata <- do.call(ps, as.list(O.val))

    if(check.optimize.normal(O.val)){
      for(stratum in strata){

      # we run the weighted OLS routine
        opt.res <- optimize.normal.model(W.obs[indices], X[indices,], strata.probas[[stratum]][indices])
      
      # then depending on the mod, we update \betas and \sigmas
        if(mod(stratum) == 1){
        # => update just B_{g} and sig_{g}
          theta[["beta"]][[stratum]] <- opt.res$beta
          theta[["sigma"]][[stratum]] <- opt.res$sigma
        
        } else {
        # => update B_{g,0}, sig_{g,0} or B_{g,1}, sig_{g,1}
        # depending on the first element of O(zi,d1.i,si)
          theta[["beta"]][[stratum]][[O.val[1]+1]] <- opt.res$beta
          theta[["sigma"]][[stratum]][[O.val[1]+1]] <- opt.res$sigma
        }
      }
    }
  }
  
  return(theta) # returns the updated theta
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

em.algo <- function(W.obs, theta.init, Z, X, ps, mod, O.fn, cp, strata.probas.init, niter=10){
  theta <- theta.init
  res.strata.probas <- NULL
  strata.probas <- strata.probas.init
  #res.list <- vector(mode='list', niter)
  for(i in seq(1, niter)){
    cat("Iteration:", i, "\n")
    
    cat("   E-step:\n")
    strata.probas <- update.strata.proba(theta, X, W.obs, O.fn, cp, ps, mod)
    strata.probas <- add.small.noise(strata.probas)
    
    cat("   M-step:\n")
    theta <- m.step(W.obs, theta, Z, X, ps, mod, O.fn, strata.probas)

    #res.list[[i]] <- list(theta, strata.probas)
    #print(as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,"wgt"]))))
    est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,"wgt"])))
    print(est.strata.prop)
    res.strata.probas <- rbind(res.strata.probas, est.strata.prop)
  }
  return(list(theta, strata.probas, res.strata.probas))
  #return(res.list)
}


# Post-processing utilities

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


estimate.causal.effect <- function(pis, wgt){
  strata.proportions <- estimate.strata.proportions(pis, wgt)
  pi.c <- with(strata.proportions, cEE + cEN + cNE + cNN)
  zs <-  with(strata.proportions, (cEN - cNE))
  return( list("ZD" = pi.c,
               "ZS" = zs,
               "DS" = zs / pi.c))
}


estimate.all.strata.proportions <- function(em.res, wgt){
  res <- sapply(seq(1, length(em.res)),
                function(i) as.vector(do.call('cbind', estimate.strata.proportions(em.res[[i]][[1]]$pi, wgt))))
  return(res)
}



