# install cnasimtools from vpparkkinen/cnasimtools
# and sort out your working directory before running

library(cnasimtools)
source(pval_hat_cna.R)


# generate cs data sets with 12,5% noise
ds <- replicate(200, noisyDat(6, noisefraction = .125), simplify = FALSE)



# grab dgs's
targets <- lapply(ds, function(x) attributes(x)$target)

outs <- lapply(targets, getoutcomes)

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
wreu <- which(reu)
length(wreu) / length(reu)

# now check if there are discoverable submodels
# of the non-significant targets that would be significant
nonsig_targets <- targets[wreu]
nst_fits <- m_fits[wreu]
nst_outs <- outs[wreu]
nst_dat <- ds[wreu]

# generate potential submodels by 
pot_subs <- vector("list", length(nonsig_targets))
for(i in seq_along(pot_subs)){
 pot_subs[[i]] <- cna(nst_dat[[i]], 
                 con = (nst_fits[[i]]$consistency - 0.4),
                 cov = (nst_fits[[i]]$consistency - 0.4),
                 outcome = nst_outs[[i]]) 
}

pot_subs <- lapply(pot_subs, csf)
pot_s_m <- lapply(pot_subs, function(x) x$condition)
any_sub <- mapply(is.submodel, pot_s_m, nonsig_targets)
w_has_sub <- lapply(any_sub, which)
idx <- unlist(lapply(w_has_sub, any))
w_has_sub <- w_has_sub[idx]
subs_w_info <- mapply(function(x, y) x[y,], 
                      pot_subs[idx], 
                      w_has_sub, 
                      SIMPLIFY = FALSE)
hassub_target <- nonsig_targets[idx]
hassub_dat <- nst_dat[idx]

ch_all_pvals <- function(x,ds){
  out <- vector("list", nrow(x))
  for(i in seq_along(out)){
    out[[i]] <- pval_hat_single(x$condition[i],
                                obs_con = x$consistency[i],
                                obs_cov = x$coverage[i],
                                dat = ds,
                                nulltype = "perm.outcome",
                                bs_samples = 1000)
  }
  return(out)
}

ch_sub_pvals <- vector("list", length(hassub_dat))
for(i in seq_along(hassub_dat)){
  ch_sub_pvals[[i]] <- ch_all_pvals(subs_w_info[[i]], ds = hassub_dat[[i]])
}

