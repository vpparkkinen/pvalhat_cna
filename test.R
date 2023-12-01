# install cnasimtools from vpparkkinen/cnasimtools
# and sort out your working directory before running
if(!require(rstudioapi)){
  setwd(system2("pwd", stdout = TRUE)) # R must run in a shell
} else {
  path <- rstudioapi::getActiveDocumentContext()$path
  Encoding(path) <- "UTF-8"
  setwd(dirname(path))
}

library(cnasimtools)
# library(parallel)
# library(foreach)
library(doParallel)
source("pval_hat_cna.R")

n.cores <- detectCores() - 5
dcluster <- makeCluster(n.cores, type = "FORK")
doParallel::registerDoParallel(cl = dcluster)
# generate cs data sets with noise
ds <- replicate(100, noisyDat(6, noisefraction = .2), simplify = FALSE)



# grab dgs's
targets <- mclapply(ds, function(x) attributes(x)$target)

outs <- mclapply(targets, getoutcomes)

# check fit of dgs to each data set
m_fits <- mcmapply(condition, targets, ds, SIMPLIFY = FALSE)
m_fits <- mclapply(m_fits, 
                        function(x) attributes(x)$info[, c("consistency", "coverage")])

res <- foreach(i = 1:length(targets)) %dopar% {
  pval_hat_single(model = targets[[i]],
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

pot_subs <- foreach(i = seq_along(nonsig_targets)) %dopar% {
  cna(nst_dat[[i]], 
      con = (nst_fits[[i]]$consistency * 0.9),
      cov = (nst_fits[[i]]$coverage * 0.9),
      outcome = nst_outs[[i]])
}



pot_subs <- mclapply(pot_subs, csf)
pot_s_m <- mclapply(pot_subs, function(x) x$condition)
any_sub <- mcmapply(is.submodel, pot_s_m, nonsig_targets, SIMPLIFY = FALSE)
w_has_sub <- mclapply(any_sub, which)
idx <- unlist(mclapply(any_sub, any))
w_has_sub <- w_has_sub[idx]
subs_w_info <- mcmapply(function(x, y) x[y,], 
                      pot_subs[idx], 
                      w_has_sub, 
                      SIMPLIFY = FALSE)

#randomly pick one submodel of the target for each data set
sub_picks <- mclapply(subs_w_info, 
                      function(x) x[sample(1:nrow(x), 1),])

hassub_target <- nonsig_targets[idx]
hassub_dat <- nst_dat[idx]

ch_sub_pvals <- foreach(i = seq_along(sub_picks)) %dopar% {
  pval_hat_single(sub_picks[[i]][1,]$condition,
                  obs_con = sub_picks[[i]][1,]$consistency,
                  obs_cov = sub_picks[[i]][1,]$coverage,
                  dat = hassub_dat[[i]],
                  nulltype = "perm.outcome",
                  bs_samples = 1000)
}


# p-values for the submodels. But now we are doing multiple tests, what do
# these p-values mean?
# ch_sub_pvals <- vector("list", length(hassub_dat))
# for(i in seq_along(hassub_dat)){
#   ch_sub_pvals[[i]] <- ch_all_pvals(subs_w_info[[i]], ds = hassub_dat[[i]])
# }

#uch_sub_pvals <- unlist(ch_sub_pvals, recursive = FALSE)
sub_signif_idx <- unlist(lapply(ch_sub_pvals, function(x) all(x < 0.05)))
length(which(sub_signif_idx)) / length(wreu)



# generate incorrect models by analyzing the non-sig data sets
# with higher fit-threshold than the target's fit

pot_overfit <- foreach(i = seq_along(nonsig_targets)) %dopar% {
  csf(cna(nst_dat[[i]], 
          con = ifelse(nst_fits[[i]]$consistency < 1,
                       (nst_fits[[i]]$consistency + 
                          (1-nst_fits[[i]]$consistency)/2),1),
          cov = ifelse(nst_fits[[i]]$coverage < 1,
                       (nst_fits[[i]]$coverage + 
                          (1-nst_fits[[i]]$coverage)/2),1),
          outcome = nst_outs[[i]]))
}




nonempty_idx <- unlist(mclapply(pot_overfit, function(x) nrow(x)>0))

pot_overfit_nonempty <- pot_overfit[nonempty_idx]
nonempty_nst_dat <- nst_dat[nonempty_idx]


pot_of_nonempty_targets <- nonsig_targets[nonempty_idx]

# make sure to collect only incorrect models
pot_minus_submodels <- foreach(i = seq_along(pot_overfit_nonempty)) %dopar% {
  pot_overfit_nonempty[[i]][
    sapply(pot_overfit_nonempty[[i]]$condition, 
           function(x) !is.submodel(x, pot_of_nonempty_targets[[i]])),
  ]
}

pot_minus_submodels_nonempty_idx <- unlist(mclapply(pot_minus_submodels, 
                                             function(x) nrow(x) > 0))

pot_of_nonempty_targets <- pot_of_nonempty_targets[pot_minus_submodels_nonempty_idx]

pot_minus_submodels_nempty <- pot_minus_submodels[pot_minus_submodels_nonempty_idx]

# sample one model
pot_overfit_picks <- mclapply(pot_minus_submodels_nempty, 
                              function(x) x[sample(1:nrow(x), 1),])


nst_dat_nonempty <- nst_dat[nonempty_idx][pot_minus_submodels_nonempty_idx]

of_pvals <- foreach(i = seq_along(pot_overfit_picks)) %dopar% {
  pval_hat_single(pot_overfit_picks[[i]][1,]$condition,
                  obs_con = pot_overfit_picks[[i]][1,]$consistency,
                  obs_cov = pot_overfit_picks[[i]][1,]$coverage,
                  dat = nst_dat_nonempty[[i]],
                  nulltype = "perm.outcome",
                  bs_samples = 1000)
}


#uof_pvals <- unlist(of_pvals, recursive = FALSE)
signif_idx_of <- unlist(lapply(of_pvals, function(x) all(x < 0.05)))
length(which(signif_idx_of)) / length(wreu)

parallel::stopCluster(cl = dcluster)
