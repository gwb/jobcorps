source("job_corps_debug.R")
source("job_corps.R")



# first simul -- initializing at true value
Y <- data.frame(assign=y$assign,
                treat=y$treat,
                work52=y$work52)
colnames(X) <- c("X0", "wgt", paste("X", seq(ncol(X)-2), sep=""))
O.52 <- get.O(Y)
strata.probas.full <- lapply(seq_len(ncol(strata.probas)), function(i) strata.probas[,i])
names(strata.probas.full) <- c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")
theta[["pi"]] <- strata.probas.full

bli <- em.algo(wages, theta, Y$assign, as.matrix(X), PS, MOD, O.52, CP, strata.probas.full, niter=30)

est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,"wgt"])))

# second simul -- initializing at modified value

theta.0 <- list(
    alpha = list(cEE=runif(K, 0, 5),
        cEN=runif(K, 0, 5),
        cNE=runif(K, 0, 5),
        cNN=runif(K, 0, 5),
        nEE=runif(K, 0, 5),
        nNN=runif(K, 0, 5)),
    beta = list(cEE=list(
                    runif(K, -3, 3),
                    runif(K, -3, 3)),
        cEN=runif(K, -3, 3),
        cNE=runif(K, -3, 3),
        nEE=runif(K, -3, 3)),
    sigma = list(cEE=list(1.2,1.5),
        cEN=2,
        cNE=0.5,
        nEE=1.3))

strata.probas.0 <- get.strata.proba(X, theta.0)
strata.probas.full.0 <- lapply(seq_len(ncol(strata.probas.0)), function(i) strata.probas.0[,i])
names(strata.probas.full.0) <- c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")
theta.0[["pi"]] <- strata.probas.full.0

bli.2 <- em.algo(wages, theta.0, Y$assign, as.matrix(X), PS, MOD, O.52, CP, strata.probas.full.0, niter=30)




# second simul

strata.assignment <- apply(strata.probas, 1, function(x) sample(c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN"), 1, prob=x))
assign.y <- sample(c(0,1), 1000, replace=T)
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


Y <- data.frame(assign=y$assign,
                treat=y$treat,
                work52=y$work52)

colnames(X) <- c("X0", "wgt", paste("X", seq(ncol(X)-2), sep=""))

O.52 <- get.O(Y)


strata.probas.full <- lapply(seq_len(ncol(strata.probas)), function(i) strata.probas[,i])
names(strata.probas.full) <- c("cEE", "cEN", "cNE", "cNN", "nEE", "nNN")

theta[["pi"]] <- strata.probas.full

est.strata.prop <- as.vector(do.call('cbind', estimate.strata.proportions(theta$pi, X[,"wgt"])))

bli <- em.algo(wages, theta, Y$assign, as.matrix(X), PS, MOD, O.52, CP, strata.probas.full, niter=5)


# debug stuff

W <- as.vector(t(do.call("cbind", strata.probas.full)))
M <- length(strata.probas.full)

mod.input <- setup.glm(X, nrow(X), M)




# more debug

foo <- update.strata.proba(theta, X, wages, O.52, CP, PS, MOD)
