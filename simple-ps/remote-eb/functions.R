
source("../generate.R")

draw.new.data <- function(theta.hat, X, assign.y){
    strata.probas.hat <- get.strata.proba(X, theta.hat)
    strata.assignment.hat <- apply(strata.probas.hat, 1, function(x) sample(c("AT", "NT", "C", "D"), 1, prob=x))
    
    table(strata.assignment.hat)
    
    #assign.y <- sample(c(0,1), N, replace=T)
    treat.y <- NULL
    for(i in seq(N)){
        if(strata.assignment.hat[i] == "AT"){
            treat.y <- c(treat.y, 1)
        }
        if(strata.assignment.hat[i] == "NT"){
            treat.y <- c(treat.y, 0)
        }
        if(strata.assignment.hat[i] == "C"){
            treat.y <- c(treat.y, assign.y[i])
        }
        if(strata.assignment.hat[i] == "D"){
            treat.y <- c(treat.y, 1-assign.y[i])
        }
    }
    
    y <- list(assign=assign.y,
              treat=treat.y)
    
    W.y <- sapply(seq(N), function(i) gen.wages(X[i,], strata.assignment.hat[i], assignment=assign.y[i], theta.hat))
    return(list(y, W.y))
}

impute.mis.W <- function(X, assign, strata.probas.post, theta.hat){
    Gs <- apply(strata.probas.post, 1, function(x) sample(c("AT", "NT", "C", "D"), 1, prob=x))
    W.mis <- sapply(seq(length(Gs)),
                    function(i) gen.wages(X[i,], Gs[i], 1-y$assign[i], theta.hat))
    return(list(Gs, W.mis))
}

get.tau.compliers <- function(W.complete, Gs){
    idx.C <- which(Gs=="C")
    return(mean(W.complete$treat[idx.C]) - mean(W.complete$control[idx.C]))
}

.eb.post.estimate <- function(W.y, assign, strata.probas.post, theta.hat){
    impute.res <- impute.mis.W(X, assign, strata.probas.post, theta.hat)
    Gs <- impute.res[[1]]
    W.mis <- impute.res[[2]]
    W.complete <- list(treat=ifelse(assign==1, W.y, W.mis),
                       control=ifelse(assign==0, W.y, W.mis))
    return(get.tau.compliers(W.complete, Gs))
}

eb.post.estimate <- function(W.y, assign, strata.probas.post, theta.hat, nsimul=100){
    return(sapply(seq(nsimul), function(i) .eb.post.estimate(W.y, assign, strata.probas.post, theta.hat)))
}



draw.tau.post.distr.old <- function(Y, strata.probas.post, assign.y, niter=100){
    tau.post.list <- vector('numeric', length=niter)
    for(i in seq(niter)){
        strata.assignment.post <- apply(strata.probas.post, 1, function(x) sample(c("AT", "NT", "C", "D"), 1, prob=x))
        Y.1 <- Y[which(strata.assignment.post == "C" & assign.y == 1)]
        Y.0 <- Y[which(strata.assignment.post == "C" & assign.y == 0)]
        tau.post.list[i] <- mean(Y.1) - mean(Y.0)
    }
    return(tau.post.list)
}
