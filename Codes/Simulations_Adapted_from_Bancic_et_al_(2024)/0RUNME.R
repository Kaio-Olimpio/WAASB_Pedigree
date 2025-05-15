# ======================================================================
#
# Script: Simulate phenotypic line breeding program with moderate GEI
#
# Script authors: J. Bancic and  D.J. Tolhurst (Adapted by Saulo Chaves)
#
# ======================================================================
#
# Description:
# This script demonstates how to simulate a phenotypic line breeding program
# with GEI using FieldSimR functionality. AlphaSimR is used to generate an
# additive genetic trait with population structure and FieldSimR is used to
# generate correlated plot errors based on a randomised complete block design.
# Each year, the simulation tracks the breeding progress by summarizing measures
# of genetic mean (average of genotype main effects in MET and TPE), variance
# (variance of genotype main effects in MET and TPE), accuracy of main genotype
# effects in MET and TPE, and MET-TPE alignment.
#

# Load necessary libraries
rm(list = ls())  # Clean up the environment
library(AlphaSimR)
library(FieldSimR)
library(data.table)

# Scenario name
Scenario <- "Pheno"
GEI <- "High"

# Import simulation parameters
source("Codes/Jon_S4_BreedingComparison/GlobalParameters.R")

# Import simulated target population of environments (TPE)
# Ge <- readRDS("Ge_ModerateGEI.rds") # Load presimulated Ge for moderate GEI from the github
# Ce <- cov2cor(Ge)
# De <- diag(diag(Ge))

# Run 20 years of breeding across multiple replicates
# High GEI: basecor = 0.2
basecor = 0.2
HGEI = list()
for (Rep in 1:nReps) {
  cat("\nWorking on replicate:", Rep, "\n")

  # Simulate founders
  source("Codes/Jon_S4_BreedingComparison/CreateFounders.R")

  # Import dataframe for storing variables
  source("Codes/Jon_S4_BreedingComparison/StoreVariables.R")
  output$rep <- Rep

  # Sample MET envrionments from TPE for each simulation year
  # Note: Environments within each sample are sampled at random
  # and without ordering
  covs_tpe = input_asr$corA
  samples_met <- sample_met(nenvs = nEnvTPE,
                            nsamples = nCycles,
                            sample.size = nenvYT5,
                            replace = TRUE,
                            cov.mat = input_asr$corA)$sample

  # Fill breeding pipeline
  source("Codes/Jon_S4_BreedingComparison/FillPipeline.R")

  # Run simulation
  for (year in 1:nCycles) {
    time <- timestamp(prefix = "", suffix = "", quiet = TRUE)
    cat("Working on year: ", year, "    (", time, ")\n", sep = "")

    # Select new parents
    source("Codes/Jon_S4_BreedingComparison/UpdateParents.R")

    # Advance the breeding year
    source("Codes/Jon_S4_BreedingComparison/AdvanceYear.R")
  }

  # Save results for this replicate
  cat("Saving results \n")
  file_name <- paste0("Output_", Scenario, GEI, ".csv")
  write.table(output, file_name, sep = ",",
    col.names = !file.exists(file_name), row.names = FALSE, append = TRUE)
  
  source("Codes/Jon_S4_BreedingComparison/Make_3hybs.R")
  
  HGEI[[Rep]] = list(HHHG = HH, LHHG = LH)
}

# Low GEI: basecor = 0.8
Scenario <- "Pheno"
GEI <- "Low"
basecor = 0.8
LGEI = list()
for (Rep in 1:nReps) {
  cat("\nWorking on replicate:", Rep, "\n")
  T
  # Simulate founders
  source("Codes/Jon_S4_BreedingComparison/CreateFounders.R")
  
  # Import dataframe for storing variables
  source("Codes/Jon_S4_BreedingComparison/StoreVariables.R")
  output$rep <- Rep
  
  # Sample MET envrionments from TPE for each simulation year
  # Note: Environments within each sample are sampled at random
  # and without ordering
  covs_tpe = input_asr$corA
  samples_met <- sample_met(nenvs = nEnvTPE,
                            nsamples = nCycles,
                            sample.size = nenvYT5,
                            replace = TRUE,
                            cov.mat = input_asr$corA)$sample
  
  # Fill breeding pipeline
  source("Codes/Jon_S4_BreedingComparison/FillPipeline.R")
  
  # Run simulation
  for (year in 1:nCycles) {
    time <- timestamp(prefix = "", suffix = "", quiet = TRUE)
    cat("Working on year: ", year, "    (", time, ")\n", sep = "")
    
    # Select new parents
    source("Codes/Jon_S4_BreedingComparison/UpdateParents.R")
    
    # Advance the breeding year
    source("Codes/Jon_S4_BreedingComparison/AdvanceYear.R")
  }
  
  # Save results for this replicate
  cat("Saving results \n")
  file_name <- paste0("Output_", Scenario, GEI, ".csv")
  write.table(output, file_name, sep = ",",
              col.names = !file.exists(file_name), row.names = FALSE, append = TRUE)
  
  source("Codes/Jon_S4_BreedingComparison/Make_3hybs.R")
  
  LGEI[[Rep]] = list(HHLG = HH, LHLG = LH)
}

save(HGEI, LGEI, file = "saves/simu_data.RDA")


# end of script

