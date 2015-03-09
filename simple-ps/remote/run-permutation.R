#source("../generate.R", chdir=T)
source("../estimation-automated.R", chdir=T)
source("../inference_permutation.R", chdir=T)


args <- commandArgs(trailingOnly = TRUE)
tau <- as.numeric(args[2])

load("simdata.Rdata")
load("theta.Rdata")


res.permutation <- do.permutation.add.robust(W.y, y, X, post.strata.proba, CP, PS, S2O, MOD, tau=tau, maxiter=30,n.restart=5,eps=0.005)

save(res.permutation, file=paste("results/res.",args[1],".Rdata", sep=""))
