# Script to generate three-way crosses ------------------------------------

simu_met = function(malepop,
                    femalepop,
                    nenv,
                    h2,
                    mumu,
                    sdmu,
                    nrept,
                    selfamf1,
                    seed,
                    rank = 7, 
                    skew = 0){
  require(FieldSimR)
  require(AlphaSimR)
  require(dplyr)
  require(tidyr)
  
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
  
  output = list()

  ## Generate F1s
  malepop@id = paste0("Ma", malepop@id)
  femalepop@id = paste0("Fe", femalepop@id)
  
  ### Cross plan and generating the F1s
  planf1 = as.matrix(expand.grid(femalepop@id, malepop@id))
  f1s = makeCross2(males = malepop, females = femalepop, crossPlan = planf1)
  
  gv_tpe <- unstr_asr_output(pop = f1s, ntraits = 1, nenvs = nEnvTPE, nreps = nrept)
  envs_met <- samples_met[[13]][1:nenv]
  gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
  
  s2e = var(gv_met$gv.Trait1) / h2 - var(gv_met$gv.Trait1)
  if (nenv > 1){
    s2ej = diag(skew_diag_mat(
      n = nenv,
      mean.var = s2e,
      shape = 3,
      scale = 1.2
    ))
  } else s2ej = s2e
  pp = findpairs(nlevels(gv_met$id) * nrept)
  if (pp[[length(pp)]][1] == pp[[length(pp)]][2]){
    ncol = nrow = pp[[length(pp)]][1]
    } else{
    ncol = ifelse(pp[[length(pp)]][1] %% nrept == 0, pp[[length(pp)]][1], pp[[length(pp)]][2])
    nrow = pp[[length(pp)]][which(pp[[length(pp)]] != ncol)]
  }
  
  error_met <- field_trial_error(
    ntraits = 1,
    nenvs = nenv,
    nblocks = nrept,
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
  f1s@pheno[,1] = mean_pheno_met
  selef1s = selectFam(
    pop = f1s,
    nFam = selfamf1,
    # trait = mean,
    use = 'pheno'
  )
  selef1s@id = paste(selef1s@mother, selef1s@father, sep = '_')
  
  message("===> Simulated F1s")
  
  ## Generate three-way hybrids
  plan = expand.grid(selef1s@id, c(malepop@id, femalepop@id))
  plan$control = apply(as.matrix(plan), 1, function(x) {
    temp = do.call(c, strsplit(x[1], split = "_"))
    x[2] %in% temp
  })
  
  plan = as.matrix(plan[which(!plan$control), -3])
  plan = plan[sample(nrow(plan), 100),]
  hybtr = makeCross2(
    females = selef1s,
    males = c(malepop, femalepop),
    crossPlan = plan
  )
  
  message("===> Simulated three-way hybrids")
  
  ## Get the true genetic values (g + ge)
  set.seed(seed+30)
  gv_df = unstr_asr_output(pop = hybtr, ntraits = 1, nenvs = nenv, nreps = nrept)
  gv_df = gv_df |> dplyr::left_join(data.frame(
    'id' = hybtr@id,
    'mum' = hybtr@mother,
    'dad' = hybtr@father
  ), by = 'id') |> dplyr::mutate(id = paste(mum, dad, sep = '//')) |>
    dplyr::reframe(
      gv.Trait1 = mean(gv.Trait1),
      .by = c(env, id, rep)
    ) |> dplyr::mutate(parents = id) |> 
    tidyr::separate(parents, into = c("mum", 'dad'), sep = '//')
  
  numfam = length(unique(gv_df$id)) * nrept
  pp = findpairs(numfam)
  if (pp[[length(pp)]][1] == pp[[length(pp)]][2]){
    ncol = nrow = pp[[length(pp)]][1]
  } else{
    ncol = ifelse(pp[[length(pp)]][1] %% nrept == 0, pp[[length(pp)]][1], pp[[length(pp)]][2])
    nrow = pp[[length(pp)]][which(pp[[length(pp)]] != ncol)]
  }

  ## Simulate plot errors
  s2e = mean(diag(De))/h2 - mean(diag(De))
  s2ej = diag(skew_diag_mat(n = nenv, mean.var = s2e, shape = 3, scale = 1.2))
  h2j = c(diag(De)/(diag(De) + s2e))
  
  set.seed(seed+40)
  error_df = field_trial_error(
    nenvs = nenv, 
    nblocks = nrept, 
    spatial.model = "AR1",
    varR = s2e, 
    ncols = ncol,
    nrows = nrow,
    return.effects = TRUE,
    col.cor = 0, 
    row.cor = 0
  )
  
  ## Simulate the phenotypes
  met_df = make_phenotypes(
    gv.df = gv_df[,c('env', 'id', 'rep', 'gv.Trait1')],
    error.df = error_df$error.df,
    randomise = TRUE,
    return.effects = FALSE
  )
  
  met_df = left_join(met_df, gv_df)

  output$simudf = met_df
  
  # output$param = data.frame(
  #   corpred = cor(with(met_df, tapply(y.Trait1, id, mean)), with(gv_df, tapply(gv.Trait1, id, mean))),
  #   cortrue = sqrt((Ge_vars[1,2] + Ge_vars[2,2]/nenv)/(Ge_vars[1,2] + Ge_vars[2,2]/nenv + mean(s2e)/nenv/nrept))
  # )
  # 
  # output$ge_met_true <- with(gv_df, tapply(gv.Trait1, list(id, env), mean))
  # output$g_met_true <- with(gv_df, tapply(gv.Trait1, id, mean))
  
  id = unique(met_df$id)
  shyb = gsub("//.*", "", id)
  output$simuped = rbind(
    data.frame(
      id = unique(gsub(".*//", "", id)),
      mum = 0,
      dad = 0
    ),
    data.frame(id = unique(shyb), temp = unique(shyb)) |> separate(temp, into = c("mum", 'dad'), sep = "_"),
    data.frame(
      id = id,
      mum = shyb,
      dad = gsub(".*//", "", id)
    )
  )
  
  message("===> Simulated multi-environment trials \n \n @@@@@@@@@@@@@@@ \n")
  
  return(output)
  
#   g_tpe_true <- with(gv_tpe, tapply(gv.Trait1, id, mean))
#   g_met_true <- with(gv_met, tapply(gv.Trait1, id, mean))
#   ge_tpe_true <- with(gv_tpe, tapply(gv.Trait1, list(id, env), mean))
#   ge_met_true <- with(gv_met, tapply(gv.Trait1, list(id, env), mean))
#   (Ge_vars <- measure_variances(Ge))
#   sigm2e
#   
#   # Main effect accuracy in TPE
#     sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/nlevels(gv_df$env) + s2e/nlevels(gv_df$env)/nlevels(gv_df$rep)))
#   # Main effect accuracy in MET
#   
#     sqrt((Ge_vars[1,2] + Ge_vars[2,2]/nlevels(gv_df$env))/(Ge_vars[1,2] + Ge_vars[2,2]/nlevels(gv_df$env) + s2e/nlevels(gv_df$env)/nlevels(gv_df$rep)))
#   # MET-TPE alignment
#   cor(g_tpe_true, g_met_true) # observed
#  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/nlevels(gv_df$env)))
#   # Accuracy for GE effects in MET
# sqrt((Ge_vars[1,2] + Ge_vars[2,2])/(Ge_vars[1,2] + Ge_vars[2,2] + s2e/nlevels(gv_df$rep)))
  
}

### low heritability ===========
LH = simu_met(
  malepop = MaleParents,
  femalepop = FemaleParents,
  nenv = nenvYT5,
  h2 = 0.2,
  mumu = MeanG,
  sdmu = sqrt(initVarg),
  nrept = 3,
  selfamf1 = 40,
  seed = 847,
  rank = 7,
  skew = 0
)

### High heritability ===========
HH = simu_met(
  malepop = MaleParents,
  femalepop = FemaleParents,
  nenv = nenvYT5,
  h2 = 0.8,
  mumu = MeanG,
  sdmu = sqrt(initVarg),
  nrept = 3,
  selfamf1 = 60,
  seed = 847,
  rank = 7,
  skew = 0
)

