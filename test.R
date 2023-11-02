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
                 cov = (nst_fits[[i]]$coverage - 0.4),
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

# p-values for the submodels. But now we are doing multiple tests, what do
# these p-values mean?
ch_sub_pvals <- vector("list", length(hassub_dat))
for(i in seq_along(hassub_dat)){
  ch_sub_pvals[[i]] <- ch_all_pvals(subs_w_info[[i]], ds = hassub_dat[[i]])
}

uch_sub_pvals <- unlist(ch_sub_pvals, recursive = FALSE)
sub_signif_idx <- unlist(lapply(uch_sub_pvals, function(x) all(x < 0.05)))
length(which(sub_signif_idx)) / length(sub_signif_idx)
# generate overfitted models  
pot_overfit <- vector("list", length(nonsig_targets))
for(i in seq_along(pot_overfit)){
  pot_overfit[[i]] <- csf(cna(nst_dat[[i]], 
                       con = ifelse(nst_fits[[i]]$consistency < 1,
                                     (nst_fits[[i]]$consistency + 
                                       (1-nst_fits[[i]]$consistency)/2),1),
                       cov = ifelse(nst_fits[[i]]$coverage < 1,
                                     (nst_fits[[i]]$coverage + 
                                       (1-nst_fits[[i]]$coverage)/2),1),
                       outcome = nst_outs[[i]])) 
}

pot_overfit_head <- lapply(pot_overfit, head)

nonempty_idx <- unlist(lapply(pot_overfit_head, function(x) nrow(x)>0))

pot_of_nonempty <- pot_overfit_head[nonempty_idx]
pot_of_nonempty_targets <- nonsig_targets[nonempty_idx]

pot_minus_submodels <- vector("list", length(pot_of_nonempty))
for(i in seq_along(pot_minus_submodels)){
  pot_minus_submodels[[i]] <- pot_of_nonempty[[i]][
    sapply(pot_of_nonempty[[i]]$condition, 
             function(x) !is.submodel(x, pot_of_nonempty_targets[[i]])),
  ]
}


nst_dat_nonempty <- nst_dat[nonempty_idx]

of_pvals <- vector("list", length(pot_minus_submodels))
for(i in seq_along(pot_minus_submodels)){
  of_pvals[[i]] <- ch_all_pvals(pot_minus_submodels[[i]], 
                                ds = nst_dat_nonempty[[i]])
}

uof_pvals <- unlist(of_pvals, recursive = FALSE)
signif_idx <- unlist(lapply(uof_pvals, function(x) all(x < 0.05)))
length(which(signif_idx)) / length(signif_idx)
