findpairs <- function(x) {
  pairs <- list()
  for (w in 1:floor(sqrt(x))) {
    if (x %% w == 0) {
      y <- x / w
      pairs <- append(pairs, list(c(w, y)))
    }
  }
  return(pairs)
}

setPhenoGCA2 = function (pop,
                         testers,
                         use = "pheno",
                         h2 = NULL,
                         inbred = FALSE,
                         nenvs_tpe,
                         nreps,
                         nenvs_met,
                         nCycles,
                         sampling = "random")
{
  if (sampling == "random")
    pick = ceiling(runif(1, 1, nCycles))
  else
    pick = sampling
  stopifnot("You picked the wrong number when set 'sampling'" = pick <= nCycles)
  simParam = get("SP", envir = .GlobalEnv)
  if (any(duplicated(pop@id))) {
    stop("This function does not work with duplicate IDs")
  }
  stopifnot(class(pop) == "Pop", class(testers) == "Pop")
  use = tolower(use)
  tmp = hybridCross(
    females = pop,
    males = testers,
    crossPlan = "testcross",
    returnHybridPop = inbred,
    simParam = simParam
  )
  if (use == "pheno") {
    gv_tpe <- unstr_asr_output(
      pop = tmp,
      ntraits = 1,
      nenvs = nenvs_tpe,
      nreps = nreps
    )
    envs_met <- samples_met[[pick]][1:nenvs_met]
    gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met, ])
    
    s2e = var(gv_met$gv.Trait1) / h2 - var(gv_met$gv.Trait1)
    if (nenvYT1 > 1)
      s2ej = diag(skew_diag_mat(
        n = nenvs_met,
        mean.var = s2e,
        shape = 3,
        scale = 1.2
      ))
    else
      s2ej = s2e
    
    pp = findpairs(nlevels(gv_met$id) * nreps)
    
    if (pp[[length(pp)]][1] == pp[[length(pp)]][2])
      ncol = nrow = pp[[length(pp)]][1]
    else{
      ncol = ifelse(pp[[length(pp)]][1] %% nreps == 0, pp[[length(pp)]][1], pp[[length(pp)]][2])
      nrow = pp[[length(pp)]][which(pp[[length(pp)]] != ncol)]
    }
    error_met <- field_trial_error(
      nenvs = nenvs_met,
      nblocks = nreps,
      varR = s2ej,
      spatial.model = "AR1",
      ncols = ncol,
      nrows = nrow
    )
    pheno_met <- make_phenotypes(gv.df = gv_met,
                                 error.df = error_met,
                                 randomise = TRUE)
    mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
    mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
    
    tmp@pheno[, 1] <- mean_pheno_met
    
    cat("h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
    
    y = tmp@pheno
  }
  else if (use == "gv") {
    y = tmp@gv
  }
  else {
    stop(paste0("use=", use, " is not a valid option"))
  }
  if (ncol(y) == 0) {
    stop(paste("No values for", use))
  }
  female = factor(tmp@mother, levels = unique(tmp@mother))
  if (nlevels(female) == 1) {
    GCAf = matrix(colMeans(y), nrow = 1)
  }
  else {
    if (testers@nInd == 1) {
      GCAf = y
    }
    else {
      tmp = aggregate(y[, 1] ~ female, FUN = mean)
      GCAf = unname(as.matrix(tmp[, -1, drop = F]))
    }
  }
  pop@pheno = GCAf
  # pop@fixEff = rep(as.integer(fixEff), pop@nInd)
  return(pop)
}

setPheno2 = function (pop,
                      h2 = NULL,
                      inbred = FALSE,
                      nenvs_tpe,
                      nreps,
                      nenvs_met,
                      nCycles,
                      sampling = "random")
{
  if (sampling == "random")
    pick = ceiling(runif(1, 1, nCycles))
  else
    pick = sampling
  stopifnot("You picked the wrong number when set 'sampling'" = pick <= nCycles)
  simParam = get("SP", envir = .GlobalEnv)
  
  gv_tpe <- unstr_asr_output(
    pop = pop,
    ntraits = 1,
    nenvs = nenvs_tpe,
    nreps = nreps
  )
  envs_met <- samples_met[[pick]][1:nenvs_met]
  gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met, ])
  
  s2e = var(gv_met$gv.Trait1) / h2 - var(gv_met$gv.Trait1)
  if (nenvs_met > 1)
    s2ej = diag(skew_diag_mat(
      n = nenvs_met,
      mean.var = s2e,
      shape = 3,
      scale = 1.2
    ))
  else
    s2ej = s2e
  
  pp = findpairs(nlevels(gv_met$id) * nreps)
  if (pp[[length(pp)]][1] == pp[[length(pp)]][2])
    ncol = nrow = pp[[length(pp)]][1]
  else{
    ncol = ifelse(pp[[length(pp)]][1] %% nreps == 0, pp[[length(pp)]][1], pp[[length(pp)]][2])
    nrow = pp[[length(pp)]][which(pp[[length(pp)]] != ncol)]
  }
  
  error_met <- field_trial_error(
    nenvs = nenvs_met,
    nblocks = nreps,
    varR = s2ej,
    spatial.model = "AR1",
    ncols = ncol,
    nrows = nrow
  )
  pheno_met <- make_phenotypes(gv.df = gv_met,
                               error.df = error_met,
                               randomise = TRUE)
  mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
  mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
  
  pop@pheno[, 1] = mean_pheno_met
  
  cat("h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
  
  return(pop)
}
