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

#library(cnaOpt)
library(cnasimtools)
library(doParallel)
source("pval_hat_cna.R")

n.cores <- 5
dcluster <- makeCluster(n.cores, type = "FORK")
registerDoParallel(cl = dcluster)


N <- 1000 # increase this to get a better estimate of the p-val distr.
vn <- 7
#ndat <- bs_dat_create(Nsets = 1, type = "cs")[[1]]
ndat <- bs_dat_create(Nsets = N, varnum = vn ,type = "cs")
# WARNING: will take time to run
#models <- replicate(N, randomCsf(ndat))
models <- replicate(N, randomCsf(vn))

#fits <- mclapply(models, function(x) condition(x, ndat))
fits <- mcmapply(condition, models, ndat, SIMPLIFY = FALSE)

fits <- mclapply(fits, function(x) attributes(x)$info[,c("consistency", "coverage")]) 
fits <- lapply(fits, 
              function(f) {f[is.na(f)] <- 0L ; return(f)})


# pvals <- foreach(i = 1:N) %dopar% {
#   pval_hat_single(models[[i]],
#                   ndat,
#                   obs_con = fits[[i]]$consistency,
#                   obs_cov = fits[[i]]$coverage,
#                   nulltype = "iid",
#                   bs_samples = 100)
# }


pvals <- foreach(i = 1:N) %dopar% {
  pval_hat_single(models[[i]],
                  ndat[[i]],
                  obs_con = fits[[i]]$consistency,
                  obs_cov = fits[[i]]$coverage,
                  nulltype = "iid",
                  bs_samples = 10000)
}


parallel::stopCluster(dcluster)

pvals_con <- unlist(lapply(pvals, '[', 1))
pvals_cov <- unlist(lapply(pvals, '[', 2))

pv_con_no <- na.omit(pvals_con)
pv_cov_no <- na.omit(pvals_cov)

# any(duplicated(pv_con_no))
# which(duplicated(pv_cov_no))
 pv_con_no <- pv_con_no + runif(length(pv_con_no), min = -0.001, max = 0.001)
#  
 pv_cov_no <- pv_cov_no + runif(length(pv_cov_no), min = -0.001, max = 0.001)

# kolmogorov-smirnov, is the distr of p-values different from uniform
ks.test(pv_con_no, "punif")
ks.test(pv_cov_no, "punif")

# chi-squared. Why are the results so different? Which one, if either, is appropriate?
chisq.test(pv_con_no)
chisq.test(pv_cov_no, p = rep(1/length(pv_cov_no), length(pv_cov_no)))

# length(which(pvals_con >= 0.99))
# length(which(pvals_cov <= 0.05))

hist(pvals_con, breaks = 40, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
hist(pvals_cov, breaks = 20, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
