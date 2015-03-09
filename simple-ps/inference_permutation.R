
#source("estimation-enhanced.R")
source("estimation-automated.R")

get.CACE.estimate <- function(theta, X){
    numerator <- 0
    denom <- 0
    for(i in seq(nrow(X))){
        num.tmp.1 <- theta$pi$C[i] * sum(X[i,] * theta$beta$C[[2]])
        num.tmp.0 <- theta$pi$C[i] * sum(X[i,] * theta$beta$C[[1]])
        numerator <- numerator + num.tmp.1 - num.tmp.0
        denom <- denom + theta$pi$C[i]
    }
    return(numerator / denom)
}


reobserve.data <- function(treat.assignment, strata.assignment){
    assign.y <- treat.assignment
    treat.y <- NULL
    for(i in seq(N)){
        if(strata.assignment[i] == "AT"){
            treat.y <- c(treat.y, 1)
        }
        if(strata.assignment[i] == "NT"){
            treat.y <- c(treat.y,0)
        }
        if(strata.assignment[i] == "C"){
            treat.y <- c(treat.y, assign.y[i])
        }
        if(strata.assignment[i] == "D"){
            treat.y <- c(treat.y, 1-assign.y[i])
        }
    }
    y <- data.frame(assign=assign.y, treat=treat.y)
    return(y)
}

.do.permutation.add <- function(W, X, strata.probas.post, treat.assignment, tau=0){

    # impute potential outcomes under null
    W.0 <- W
    W.0[treat.assignment == 1] <- W[treat.assignment == 1] - tau

    W.1 <- W
    W.1[treat.assignment==0] <- W[treat.assignment == 0] + tau

    # impute strata
    strata.probas.post <- do.call('cbind', strata.probas.post)
    strata.assignment <- apply(strata.probas.post, 1,
                               function(x) sample(c("AT", "NT", "C", "D"), 1, prob=x))

    # rerandomize
    n.treat.assignment <- sample(treat.assignment, length(treat.assignment))

    # reobserve
    n.W <- rep(NA, length(W))
    n.W[n.treat.assignment == 0] <- W.0[n.treat.assignment == 0]
    n.W[n.treat.assignment == 1] <- W.1[n.treat.assignment == 1]

    Y <- reobserve.data(n.treat.assignment, strata.assignment)
    n.O.52 <- get.O(Y)
    return(list(W=n.W, Y=Y, n.O.52=n.O.52))
}

do.permutation.add <- function(W, X, strata.probas.post, treat.assignment, maxiter=10, tau=0, eps=0.01){
    # performs the permutation
    perm <- .do.permutation.add(W, X, strata.probas.post, treat.assignment, tau)
    
    # initialize parameters for EM based on new permutation
    theta.init <- initialize.params(perm$W, perm$Y, X)
    theta.init[["pi"]] <- list(AT=rep(0.25, N),
                           NT=rep(0.25, N),
                           C=rep(0.25,N),
                           D=rep(0.25,N))
    
    res.em <- em.algo.optim.stopping(perm$W, theta.init, perm$Y$assign, as.matrix(X), perm$n.O.52, maxiter=maxiter,eps=eps)
    res <- get.CACE.estimate(res.em[[1]], X)
    return(res)
}

# different initialization scheme
do.permutation.add.alt <- function(W, X, strata.probas.post, treat.assignment, theta.init, maxiter=10, tau=0, eps=0.01){
    # performs the permutation
    perm <- .do.permutation.add(W, X, strata.probas.post, treat.assignment, tau)
    
    res.em <- em.algo.optim.stopping(perm$W, theta.init, perm$Y$assign, as.matrix(X), perm$n.O.52, maxiter=maxiter,eps=eps)
    res <- get.CACE.estimate(res.em[[1]], X)
    return(res)
}

# with multiple restarts!
do.permutation.add.robust <- function(W, y, X, strata.probas.post, CP.fn, PS.fn, S2O.fn, MOD.fn, tau=0, n.restart=5, eps=0.001, maxiter=10, method="BFGS"){
    
    perm <- .do.permutation.add(W, X, strata.probas.post, y$assign, tau)
    n.O <- get.O(perm$Y)
    res.em <- em.algo.optim.stopping.multistart(perm$W, perm$Y, X, n.O, CP.fn, PS.fn, S2O.fn, MOD.fn, n.restart=n.restart,maxiter=maxiter, eps=eps, method=method)
    res <- get.CACE.estimate(res.em[[1]],X)
    return(res)
}
