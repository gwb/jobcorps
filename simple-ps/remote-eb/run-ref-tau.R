source("../generate.R")
source("functions.R")

load("hat.Rdata")

tau.hat.distr <- eb.post.estimate(W.y, y$assign, strata.probas.post, theta.hat, nsimul=100)

save(tau.hat.distr, file="tau-ref-distr.Rdata")
