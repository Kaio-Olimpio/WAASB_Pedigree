# Store variables

# Set up dataframe to store simulation variables
output = data.frame(year = rep(1:nCycles, times = nStages),
                    rep  = numeric(nCycles*nStages),
                    scenario    = Scenario,
                    GEI         = GEI,
                    stage       = rep(c("Parents","F1","DH/YT1","YT2","YT3","YT4","YT5"),
                                      each = nCycles),
                    mean_gv_tpe_male = numeric(nCycles*nStages),
                    mean_gv_met_male = numeric(nCycles*nStages),
                    var_gv_tpe_male = numeric(nCycles*nStages),
                    var_gv_met_male = numeric(nCycles*nStages),
                    mean_gv_tpe_female = numeric(nCycles*nStages),
                    mean_gv_met_female = numeric(nCycles*nStages),
                    var_gv_tpe_female = numeric(nCycles*nStages),
                    var_gv_met_female = numeric(nCycles*nStages),
                    acc_tpe = numeric(nCycles*nStages),
                    acc_met = numeric(nCycles*nStages),
                    met_tpe = numeric(nCycles*nStages))

# Function to automate saving output
save_output <- function(df, malepop, femalepop, stage, year, nenvs = NULL, gv_tpe_male, gv_tpe_female, 
                        GCAf = NULL, GCAm = NULL) {
  # Note: if nenvs not specified, then MET values returned as NA
  # Note: if pheno_met not specified, then phenotypic values returned as NA
  require(data.table)

  # Genotype main effects
  # gv_means <- with(gv_tpe, tapply(gv.Trait1, list(id,env), mean)) # slow
  setDT(gv_tpe_male)
  setDT(gv_tpe_female)
  gv_means_male <- gv_tpe_male[, .(mean_gv = mean(gv.Trait1)), by = .(env, id)]
  gv_means_male <- matrix(gv_means_male$mean_gv, nrow = malepop@nInd, ncol = uniqueN(gv_tpe_male$env))
  mean_gv_tpe_male <- rowMeans(gv_means_male)
  mean_gv_met_male <- if (is.null(nenvs)) {
    NA
  } else if (nenvs == 1) {
    gv_means_male[, samples_met[[year]][1]]
  } else{
    rowMeans(gv_means_male[, samples_met[[year]][1:nenvs]])
  }
  gv_means_female <- gv_tpe_female[, .(mean_gv = mean(gv.Trait1)), by = .(env, id)]
  gv_means_female <- matrix(gv_means_female$mean_gv, nrow = femalepop@nInd, ncol = uniqueN(gv_tpe_female$env))
  mean_gv_tpe_female <- rowMeans(gv_means_female)
  mean_gv_met_female <- if (is.null(nenvs)) {
    NA
  }else if (nenvs == 1) {
    gv_means_female[, samples_met[[year]][1]]
  }else{
    rowMeans(gv_means_female[, samples_met[[year]][1:nenvs]])
  }

  # phenotypic main effects when available
  # if (is.null(pheno_met)) {
  #   mean_pheno_met <- NA
  # } else {
  #   pheno_means <- with(pheno_met, tapply(y.Trait1, list(id, env), mean))
  #   mean_pheno_met <- if (is.null(nenvs)) NA
  #                     else if (nenvs == 1) pheno_means[,1]
  #                     else rowMeans(pheno_means[,1:nenvs])
  # }

  index <- which(df$year == year & df$stage == stage)
  df[index, c("mean_gv_tpe_male", "mean_gv_met_male")] <- c(mean(mean_gv_tpe_male), mean(mean_gv_met_male))
  df[index, c("var_gv_tpe_male", "var_gv_met_male")]   <- c(var(mean_gv_tpe_male), var(mean_gv_met_male))
  df[index, c("mean_gv_tpe_female", "mean_gv_met_female")] <- c(mean(mean_gv_tpe_female), mean(mean_gv_met_female))
  df[index, c("var_gv_tpe_female", "var_gv_met_female")]   <- c(var(mean_gv_tpe_female), var(mean_gv_met_female))
  
  df[index, c("acc_tpe", "acc_met", "met_tpe")] <-
    c(if (is.null(GCAf) | is.null(GCAm)) NA else cor(c(mean_gv_met_male,mean_gv_met_female), c(GCAm, GCAf)),
      if (is.null(GCAf) | is.null(GCAm)) NA else cor(c(mean_gv_tpe_male, mean_gv_tpe_female), c(GCAm, GCAf)),
      if (is.null(GCAf) | is.null(GCAm))     NA else cor(c(mean_gv_met_male,mean_gv_met_female), c(GCAm, GCAf)))

  return(df)
}
