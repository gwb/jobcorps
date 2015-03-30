
source("functions.R")
load("hat.Rdata")

set.seed(123)

for(i in seq(100)){
    print(i)
    res.i <- draw.new.data(theta.hat, X, y$assign)
    y <- res.i[[1]]
    W.y <- res.i[[2]]
    save(y, W.y, file=paste("data-replicates/bs-sample_",i,".Rdata",sep=""))
}

