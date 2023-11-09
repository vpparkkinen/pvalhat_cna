# check that the p-val test does what it should, i.e. every p-val
# equally probable,
# , i.e. multiple tests produce p-val <0.05 about 5% of time, 
# when null of pure noise data is really true

if(!require(rstudioapi)){
  setwd(system2("pwd", stdout = TRUE)) # R must run in a shell
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}


library(cnasimtools)
library(doParallel)
source("pval_hat_cna.R")

n.cores <- 5
dcluster <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl = dcluster)


ndat <- bs_dat_create(Nsets = 1, type = "cs")[[1]]
# WARNING: will take time to run
N <- 200 # increase this to get a better estimate of the p-val distr.
models <- replicate(N, randomCsf(ndat))

fits <- mclapply(models, function(x) condition(x, ndat))
fits <- mclapply(fits, function(x) attributes(x)$info[,c("consistency", "coverage")]) 

pvals <- foreach(i = 1:N) %dopar% {
  pval_hat_single(models[[i]],
                  ndat,
                  obs_con = fits[[i]]$consistency,
                  obs_cov = fits[[i]]$coverage,
                  nulltype = "iid",
                  bs_samples = 1000)
}


parallel::stopCluster(dcluster)

pvals_con <- unlist(lapply(pvals, '[', 1))
pvals_cov <- unlist(lapply(pvals, '[', 2))

ks.test(pvals_con, runif(1e6))
ks.test(pvals_cov, runif(1e6))

# length(which(pvals_con >= 0.99))
# length(which(pvals_cov <= 0.05))

hist(pvals_con, breaks = 40, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
hist(pvals_cov, breaks = 20, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
