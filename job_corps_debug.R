
#set.seed(3341)

set.seed(123)


# The synthetic data we use will have same variance (to simplify)

K <- 10 # number of covariates including intercept and weights
#N <- 1000
N <- 6000 # total number of units

X <- matrix(0, nrow=N, ncol=K)
X[,1] <- 1
X[,2] <- runif(N)
X[,3] <- rnorm(N)
X[,4] <- rnorm(N)
X[,5] <- rnorm(N)
X[,6] <- rnorm(N)
X[,7] <- rnorm(N)
X[,8] <- rnorm(N)
X[,9] <- rnorm(N)
X[,10] <- rnorm(N)

alpha <- list(cEE=runif(K, 0, 5),
              cEN=runif(K, 0, 5),
              cNE=runif(K, 0, 5),
              cNN=runif(K, 0, 5),
              nEE=runif(K, 0, 5),
              nNN=runif(K, 0, 5))

strong.idx <- sample(seq(3, K), 6)
for(i in seq(6)){
    alpha[[i]][strong.idx[i]] <- runif(1,8,10)
}

beta <- list(cEE=list(
                 runif(K, -3, 3),
                 runif(K, -3, 3)),
             cEN=runif(K, -3, 3),
             cNE=runif(K, -3, 3),
             nEE=runif(K, -3, 3))

sigma <- list(cEE=list(1,1),
              cEN=1,
              cNE=1,
              nEE=1)

theta <- list(alpha=alpha,
              beta=beta,
              sigma=sigma)

get.strata.proba <- function(X, theta){
    alph <- do.call('cbind', theta$alpha)
    lin.comp <- X %*% alph
    probas <- exp(lin.comp) / rowSums(exp(lin.comp))
    return(probas)
}


check.stratum.treat <- function(stratum, treat){
    if(stratum == "cEN" && treat == 0){
        #stop("cEN and treat == 0 are not compatible for wages")
        return(FALSE)
    }
    if(stratum == "cNE" && treat == 1){
        #stop("cNE and treat == 1 are not compatible for wages")
        return(FALSE)
    }
    return(TRUE)
}



# generates synthetic log wages from covariates, and parameters
gen.log.wages <- function(Xi, stratum, treat, theta){
    if(!check.stratum.treat(stratum, treat)){
        return(NA)
    }

    if(stratum == "cNN" || stratum == "nNN"){return(NA)}
    
    if(stratum == "cEE"){
        #bet <- theta$beta[["cEE"]][[2-treat]]
        #sig <- theta$sigma[["cEE"]][[2-treat]]
        bet <- theta$beta[["cEE"]][[1+treat]]
        sig <- theta$sigma[["cEE"]][[1+treat]]
    } else {
        bet <- theta$beta[[stratum]]
        sig <- theta$sigma[[stratum]]
    }
    return(rnorm(1, sum(Xi * bet), sig))
}


strata.probas <- get.strata.proba(X, theta)
strata.assignment <- apply(strata.probas, 1, function(x) sample(c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN"), 1, prob=x))

table(strata.assignment)


assign.y <- sample(c(0,1), N, replace=T)

#treat.y <- sapply(assign.y, function(x) ifelse(x==1, sample(c(0,1), 1, prob=c(0.2, 0.8)), 0))

treat.y <- NULL
for(i in seq(N)){
    if(assign.y[i] == 1 && strsplit(strata.assignment[i],"")[[1]][1] == "c"){
        treat.y <- c(treat.y, 1)
    } else {
        treat.y <- c(treat.y, 0)
    }
}

work.y <- sapply(seq(N), function(i) ifelse(strsplit(strata.assignment[i],"")[[1]][3-treat.y[i]] == "E", 1, 0))

    

y <- list(assign = assign.y,
          treat = treat.y,
          work52 = work.y)



log.wages <- sapply(seq(N), function(i) gen.log.wages(X[i,], strata.assignment[i], treat.y[i], theta))

wages <- exp(log.wages)


