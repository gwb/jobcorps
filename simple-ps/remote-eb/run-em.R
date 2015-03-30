source("../generate.R")
source("../estimation-automated.R")

source("functions.R")

O.fn <- get.O(y)

em.res <- em.algo.optim.stopping.multistart(W.y, y, X, O.fn, CP, PS, S2O, MOD, maxiter=30)

strata.probas.post <- do.call("cbind", em.res[[2]])
theta.hat <- em.res[[1]]


save(theta.hat, strata.probas.post, y, W.y, X, file="hat.Rdata")
