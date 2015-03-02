
set.seed(123)

N <- 5000
X <- matrix(0, nrow=N, ncol=5)

X[,1] <- 1
X[,2] <- rnorm(N)
X[,3] <- rnorm(N)
X[,4] <- rnorm(N)
X[,5] <- rnorm(N)

alpha <- list(AT = c(0.5, 1, 2, 3, 4),
              NT = c(0.5, 2, 3, 4, 1),
              C = c(0.5, 3, 4, 1, 2),
              D = c(0.5, 4, 1, 2, 3))

beta <- list(AT = list(
                 runif(5, -3, 3),
                 runif(5, -3, 3)),
             NT = list(
                 runif(5, -3, 3),
                 runif(5, -3, 3)),
             C = list(
                 runif(5, -3, 3),
                 runif(5, -3, 3)),
             D = list(
                 runif(5, -3, 3),
                 runif(5, -3, 3)))

sigma <- list(AT = list(0.1,0.1),
              NT = list(0.1,0.1),
              C = list(0.1,0.1),
              D = list(0.1,0.1))

theta <- list(alpha=alpha,
              beta=beta,
              sigma=sigma)


get.strata.proba <- function(X, theta){
    alph <- do.call('cbind', theta$alpha)
    lin.comp <- X %*% alph
    probas <- exp(lin.comp) / rowSums(exp(lin.comp))
    return(probas)
}


gen.wages <- function(Xi, stratum, assignment, theta){
    bet <- theta$beta[[stratum]][[1+assignment]]
    sig <- theta$sigma[[stratum]][[1+assignment]]
    return(rnorm(1, sum(Xi * bet), sig))
}


#

strata.probas <- get.strata.proba(X, theta)
strata.assignment <- apply(strata.probas, 1, function(x) sample(c("AT", "NT", "C", "D"), 1, prob=x))

table(strata.assignment)

assign.y <- sample(c(0,1), N, replace=T)
treat.y <- NULL
for(i in seq(N)){
    if(strata.assignment[i] == "AT"){
        treat.y <- c(treat.y, 1)
    }
    if(strata.assignment[i] == "NT"){
        treat.y <- c(treat.y, 0)
    }
    if(strata.assignment[i] == "C"){
        treat.y <- c(treat.y, assign.y[i])
    }
    if(strata.assignment[i] == "D"){
        treat.y <- c(treat.y, 1-assign.y[i])
    }
}

y <- list(assign=assign.y,
          treat=treat.y)

W.y <- sapply(seq(N), function(i) gen.wages(X[i,], strata.assignment[i], assignment=assign.y[i], theta))
