# install cnasimtools from vpparkkinen/cnasimtools
# and sort out your working directory before running

library(cnasimtools)
source(pval_hat_cna.R)


# generate cs data sets with 12,5% noise
ds <- replicate(200, noisyDat(6, noisefraction = .125), simplify = FALSE)

# grab dgs's
targets <- lapply(ds, function(x) attributes(x)$target)

# check fit of dgs to each data set
m_fits <- mapply(condition, targets, ds, SIMPLIFY = FALSE)
m_fits <- lapply(m_fits, 
                        function(x) attributes(x)$info[, c("consistency", "coverage")])


# calculate p-vals for dgs model for each data set
res <- vector("list", length(m_fits))
for(i in seq_along(m_fits)){
  res[[i]] <- pval_hat_single(model = targets[[i]],
                              obs_con = m_fits[[i]][,1],
                              obs_cov = m_fits[[i]][,2],
                              dat = ds[[i]],
                              nulltype = "perm.outcome")  
}
# what proportion of dgs's fail to reject null in either
reu <- unlist(lapply(res, function(x) any(unlist(x) > 0.05)))
length(which(reu)) / length(reu)
