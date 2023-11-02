# check that the p-val test does what it should, i.e. every p-val
# equally probable,
# , i.e. multiple tests produce p-val <0.05 about 5% of time, 
# when null of pure noise data is really true

source(pvsl_hat_cna.R)

ndat <- bs_dat_create(Nsets = 1, type = "cs")[[1]]
# WARNING: will take time to run
N <- 1000 # increase this to get a better estimate of the p-val distr.
models <- replicate(N, randomCsf(ndat))

fits <- lapply(models, function(x) condition(x, ndat))
fits <- lapply(fits, function(x) attributes(x)$info[,c("consistency", "coverage")]) 

pvals <- vector("list", length(models))
for(i in seq_along(pvals)){
  pvals[[i]] <- pval_hat_single(models[[i]],
                                ndat,
                                obs_con = fits[[i]]$consistency,
                                obs_cov = fits[[i]]$coverage,
                                nulltype = "iid")
}
#pvals

pvals_con <- unlist(lapply(pvals, '[', 1))
pvals_cov <- unlist(lapply(pvals, '[', 2))

length(which(pvals_con <= 0.05))
length(which(pvals_cov <= 0.05))

hist(pvals_con, breaks = 20, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
hist(pvals_cov, breaks = 20, xaxt = 'n')
axis(side = 1, at=seq(0,1,0.05))
