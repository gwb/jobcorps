
source("generate.R")
source("estimation-enhanced.R")

#
# Comprehensive tests checking whether the "shortcut" implementation is correct
#


O.fn <- get.O(y)
strata.probas.list <- lapply(seq_len(ncol(strata.probas)), function(i) strata.probas[,i])
names(strata.probas.list) <- c("AT", "NT", "C", "D")
theta[["pi"]] <- strata.probas.list


# Testing likelihoods implementations

obs.lik.full(X, W.y, theta, O.fn)
obs.lik(X, W.y, theta, PS, O.fn)


# testing E-step

theta.post.full <- update.posterior.strata.proba.full(theta, X, W.y, O.fn)
theta.post <- update.posterior.strata.proba(theta, X, W.y, O.fn, CP, PS, MOD)

all(abs(do.call('cbind', theta.post.full) - do.call('cbind', theta.post)) < 0.00000000001)


# testing multinomial optimization

## Default optim is Nelder-Mead. We also try BFGS
theta.post.noise <- add.small.noise(theta.post)
alphas.opt <- optimize.multinom.model(X, theta.post.noise)
alphas.opt.optim <- optimize.multinom.model.optim(X, theta.post.noise)
alphas.opt.optim.nonoise <- optimize.multinom.model.optim(X, theta.post)
alphas.opt.optim.bfgs <- optimize.multinom.model.optim(X, theta.post.noise, "BFGS")

# With BFGS we get exactly the same results as with the multinomial logistic regression (except much faster)
alphas.opt.par <- as.vector(c(alphas.opt$AT, alphas.opt$NT, alphas.opt$C))
alphas.opt.optim.par.bfgs <- alphas.opt.optim.bfgs[[1]]$par
alphas.opt.optim.par <- alphas.opt.optim[[1]]$par
all(abs(alphas.opt.par - alphas.opt.optim.par.bfgs) < 0.001)

# BFGS and multinom perform better than Nelder-Mead
fitfn.multinom(alphas.opt.par, X, theta.post.noise)
fitfn.multinom(alphas.opt.optim.par, X, theta.post.noise)
fitfn.multinom(alphas.opt.optim.par.bfgs, X, theta.post.noise)

# Timing comparison
system.time(alphas.opt <- optimize.multinom.model(X, theta.post.noise))   # 1.7 sec
system.time(alphas.opt.optim.bfgs <- optimize.multinom.model.optim(X, theta.post.noise, "BFGS")) # 0.7 sec


# testing normal optimization

## the methods find the same beta, but not the same sigma.
tmp.indices <- O.fn(0,1)
res.norm <- optimize.normal.model(W.y[tmp.indices], X[tmp.indices,],
                                  theta.post[["AT"]][tmp.indices])
res.norm.optim <- optimize.normal.model.optim(W.y[tmp.indices], X[tmp.indices,],
                                              theta.post[["AT"]][tmp.indices],
                                              method="BFGS")
res.norm.optim$beta
res.norm.optim$sigma
res.norm

## the optim method gives higher likelihood. see warning for the lm method
res.norm.pars <- c(res.norm$beta, res.norm$sigma)
res.norm.optim.pars <- c(res.norm.optim$beta, res.norm.optim$sigma)
fitfn.normal(res.norm.pars, W.y[tmp.indices], X[tmp.indices,],
             theta.post[["AT"]][tmp.indices])
fitfn.normal(res.norm.optim.pars, W.y[tmp.indices], X[tmp.indices,],
             theta.post[["AT"]][tmp.indices])



tmp.indices <- O.fn(1,1)
res.norm <- optimize.normal.model(W.y[tmp.indices], X[tmp.indices,],
                                  theta.post[["AT"]][tmp.indices])
res.norm.optim <- optimize.normal.model.optim(W.y[tmp.indices], X[tmp.indices,],
                                              theta.post[["AT"]][tmp.indices],
                                              method="BFGS")
res.norm.optim$beta
res.norm.optim$sigma
res.norm


# testing m-step

## the mixte one checks that if I use the linear regression in the long form, then I get the same
## as the previous implementation, which is the case.
res.mstep.optim <- m.step.optim(y$assign, W.y, theta, X, O.fn, theta.post.noise, method="BFGS")
res.mstep <- m.step(y$assign, W.y, theta, X, PS, MOD, O.fn, S2O, theta.post.noise)
res.mstep.mixte <- m.step.mixte(y$assign, W.y, theta, X, O.fn, theta.post.noise, method="BFGS")


# testing full algo

## first notes: the observed likelihood for the optim versions is much higher (almost x 3),
## is pretty much monotonous (to 3 digit precision) while the non optim version has a sharp
## drop in likelihood
res.em.optim <- em.algo.optim(W.y, theta, y$assign, X, O.fn, niter=2)
res.em <- em.algo.for(W.y, theta, y$assign, X, PS, MOD, O.fn, CP, S2O, niter=2)

## As hinted in all the above, the problem might come from the *lm* estimate of the variance which
## is pretty bad for some reason. Here, we evaluate the final likelihood of the non optim version,
## but we set the sigma values to the sigma values obtained from the optim version. Hurray, we recover
## pretty much the same likelihood as that of the optim version!!

theta.res.em <- res.em[[1]]
theta.res.em.tmp <- res.em[[1]]
theta.res.em.tmp$sigma <- res.em.optim[[1]]$sigma

obs.lik.full(X, W.y, theta.res.em, O.fn)
obs.lik.full(X, W.y, theta.res.em.tmp, O.fn)


# Now we test the optim algo when initializing at the wrong place

## It works! Converges to the truth. Note however that it is fairly sensitive to the initialization. So
## I really need to think about good initialization schemes. Since the function initialize.params is
## partially random, I could try multiple random initializations.
N <- length(W.y)
theta.init <- initialize.params(W.y, y, X)
names(strata.probas.list) <- c("AT", "NT", "C", "D")
theta.init[["pi"]] <- list(AT=rep(0.25, N),
                           NT=rep(0.25, N),
                           C=rep(0.25,N),
                           D=rep(0.25,N))

n.res.em.optim <- em.algo.optim(W.y, theta.init, y$assign, X, O.fn, niter=10)
