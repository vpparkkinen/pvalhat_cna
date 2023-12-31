library(cna)

# create a list of fs or cs data sets w/ iid uniform/binom variables
bs_dat_create <- function(Nsets = 1e3, 
                          size = 30, 
                          varnum = 7,
                          type = c("cs", "fs"),
                          varnames = LETTERS[1:varnum]){
  type = match.arg(type)
  if (type == "fs"){
    c <- quote(runif(size, min = 0, max = 1))
  } 
  if (type == "cs"){
    c <- quote(rbinom(n = size, size = 1, prob = 0.5))
  }
  dsets <- vector("list", Nsets)
  for(i in 1:Nsets){
    dsets[[i]] <- data.frame(setNames(
      replicate(varnum, eval(c), simplify = FALSE), varnames))
  }
  return(dsets)
}

# permutate_vec <- function(x){
#   out <- sample(x, length(x), replace = FALSE)
#   return(out)
# }

permutate_cols <- function(dat, colnames = NA){
  dat <- ct2df(dat)
  l <- nrow(dat)
  what_perm <- if(all(is.na(colnames))){names(dat)} else {colnames}
  for(i in what_perm){
    dat[,i] <- sample(dat[,i], l, replace = FALSE)
  }
  return(dat)
}


# full_perm_dist <- function(dat, outcome){
#   dat <- ct2df(dat)
#   oi <- which(names(dat) == outcome)
#   dmo <- dat[,-oi]
#   ouc <- dat[, oi]
#   ol <- length(ouc)
#   perms <- gtools::permutations(ol, ol, 1:ol)
# }


getoutcomes <- function(x){
  x <- noblanks(x)
  asfs <- unlist(extract_asf(x))
  out <- unlist(lapply(asfs, rhs))
  return(out)
}

# create "bootstrap" data sets. 
# Currently only works for fs or cs.
# x is the actual analyzed data, used to detect type and sample size.
sim_null <- function(x,
                     model,
                     outcomes = NA,
                     Nsets = 1e3, 
                     nulltype){
  #ntype <- match.arg(nulltype)
  if(nulltype == "perm.outcome"){outcomes <- getoutcomes(model)}
  facs <- names(x)
  N <- nrow(x)
  type <- attributes(configTable(x))$type
  if(nulltype=="iid"){
  out <- bs_dat_create(Nsets = Nsets, 
                       size = N,
                       #size = N*40,
                       varnames = facs,
                       varnum = length(x),
                       type = type)
  } else {
    if(any(is.na(outcomes))){outcomes <- names(x)}
    out <- replicate(Nsets, permutate_cols(x, colnames = outcomes), simplify = FALSE)
  }
  return(out)
}


# simulate empirical distribution of con/cov for a model
# under null. x is the analyzed data set.

ccov_dist_single_model <- function(dat, 
                                   model,
                                   nulltype = c("perm.outcome", "iid", "perm.all"),
                                   Nsets = 1e3){
  dats <- sim_null(x = dat, 
                   model = model, 
                   Nsets = Nsets, 
                   nulltype = match.arg(nulltype))
  ccov <- lapply(dats, function(x) condition(model, x))
  ccov <- lapply(ccov, 
                 function(x) attributes(x)$info[,c("consistency", "coverage")])
  out <- do.call(rbind, ccov)
  return(out)
}


# simulate empirical p-val for consistency/coverage of a single model
# dat is the analyzed data set. obs_stat is observed actual con or cov.
# bs_samples is number of bootstrap samples.


pval_hat_single <- function(model, 
                            dat,
                            #stat_type = c("consistency", "coverage"),
                            obs_con = attributes(condition(model, dat))$info[,"consistency"],
                            obs_cov = attributes(condition(model, dat))$info[,"coverage"],
                            nulltype = c("perm.outcome", "iid", "perm.all"),
                            bs_samples = 1e3){
  
  nulltype <- match.arg(nulltype)
  simdist <- ccov_dist_single_model(dat, 
                                    model = model, 
                                    Nsets = bs_samples,
                                    nulltype = nulltype)
  concov <- c(obs_con, obs_cov)
  #stat <- match.arg(stat_type)
  out <- vector("numeric", 2)
  ccn <- c("consistency", "coverage")
  for (i in 1:2){
    if(is.na(concov[i])){out[i] <- NA} else {
      what_dist <- simdist[,which(names(simdist) == ccn[i])]
      n <- length(what_dist)
      r <- length(what_dist[what_dist >= concov[i]])
      #r <- length(what_dist[what_dist > concov[i]]) 
      #pval <- (r + 1) / (n + 1) # think this over!
      #out[i] <- r / n
      out[i] <- (r+1) / (n+1)
    }

  }
  
  names(out) <- c("p-val_con", "p-val_cov")
  #out <- pval
  attr(out, "model") <- model
  attr(out, "observed con and cov") <- concov
  return(out)
}

msc_pval_hat <- function(model, dat){
  mc <- frscore:::decompose_model(model)
  mscs <- lapply(mc$lhss, function(x) gsub("\\+", "", x))
  mscs <- lapply(mscs, function(x) unlist(strsplit(x, "")))
  for(i in seq_along(mscs)){
    mscs[[i]] <- sapply(mscs[[i]], function(x) paste0(x, "->", mc$rhss[i]), 
                        USE.NAMES = FALSE)
  }
  mscs <- unlist(mscs)
  ccl <- lapply(mscs, function(x) 
    attributes(condition(x, dat))$info[,c("consistency", "coverage")])
  #names(ccl) <- mscs
  ccl <- do.call(rbind, ccl)
  out <- mapply(pval_hat_single, 
                model = mscs, 
                obs_con = ccl$consistency,
                obs_cov = ccl$coverage,
                MoreArgs = list(dat = dat, nulltype = "perm.outcome"),
                SIMPLIFY = FALSE)
  return(out)
}


#EXAMPLE

# get some models
# jsre <- csf(cna(d.jobsecurity, con = .6, cov = .8, outcome = "JSR"))
#   
# 
# ## test one of them
# pval_hat_single("C + R <-> JSR", d.jobsecurity, obs_con = 0.722, obs_cov = 0.968,
#                 nulltype = "perm.outcome")
# 
# pval_hat_single("C + R <-> JSR", d.jobsecurity, obs_con = 0.722, obs_cov = 0.968,
#                 nulltype = "perm.all")
# 
# pval_hat_single("C + R <-> JSR", d.jobsecurity, obs_con = 0.722, obs_cov = 0.968,
#                 nulltype = "iid")

  