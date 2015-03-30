
source("../generate.R")
source("../estimation-automated.R")
source("functions.R")

args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])


# loads bootstrapped data
load(paste("data-replicates/bs-sample_",idx,".Rdata", sep=""))

# finds mle of model for bootstrapped data
# note: the y and W.y are loaded from the bootstratpped data file!
O.fn <- get.O(y)
em.res <- em.algo.optim.stopping.multistart(W.y, y, X, O.fn, CP, PS, S2O, MOD, maxiter=30, n.restart=3)
strata.probas.post <- do.call("cbind", em.res[[2]])
theta.hat <- em.res[[1]]

# computes tau posterior distriubtion conditional on theta = \hat{theta}
tau.hat.distr <- eb.post.estimate(W.y, y$assign, strata.probas.post, theta.hat, nsimul=100)


save(tau.hat.distr, em.res, W.y, y, file = paste("results-tau/res_",idx,".Rdata", sep=""))





