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




# create "bootstrap" data sets where everything is independent
# for the inference test. Currently only works for fs or cs.
# x is the actual analyzed data, used to detect type and sample size.
sim_null <- function(x, Nsets = 1e3){
  facs <- names(x)
  N <- nrow(x)
  type <- attributes(configTable(x))$type
  # if(type=="fs"){
  out <- bs_dat_create(Nsets = Nsets, 
                       size = N,
                       varnames = facs,
                       varnum = length(x),
                       type = type)
  #}
  return(out)
}


# simulate empirical distribution of con/cov for a model
# under null. x is the analyzed data set.

ccov_dist_single_model <- function(dat, model, Nsets = 1e3){
  dats <- sim_null(x = dat)
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
                            stat_type = c("consistency", "coverage"),
                            obs_stat,
                            bs_samples = 1e3){
  simdist <- ccov_dist_single_model(dat, model = model, Nsets = bs_samples)
  stat <- match.arg(stat_type)
  what_dist <- simdist[,which(names(simdist)==stat)]
  n <- length(what_dist)
  r <- length(what_dist[what_dist >= obs_stat])
  #pval <- (r + 1) / (n + 1)
  pval <- r / n
  names(pval) <- "p-value"
  out <- pval
  attr(out, "model") <- model
  attr(out, "observed fit") <- obs_stat
  return(out)
}

#EXAMPLE

# get some models
jsre <- csf(cna(d.jobsecurity, con = .6, cov = .8, outcome = "JSR"))
  

# try one of them for consistency pval
pval_hat_single("C + R <-> JSR", d.jobsecurity, obs_stat = 0.722)
  