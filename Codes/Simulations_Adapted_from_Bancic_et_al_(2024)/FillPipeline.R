# Populate breeding stages with unique genotypes
source("Codes/Jon_S4_BreedingComparison/Aux_fun.R")


for(year in 1:nStages){
  cat("FillPipeline year:", year, "of", nStages,"\n")

  # Stage 1
  MaleF1   = randCross(MaleParents, nCrosses)
  FemaleF1 = randCross(FemaleParents, nCrosses)
  
  # Stage 2
  if(year < 6){
    # p = P[6-cohort]
    
    MaleDH   = makeDH(MaleF1, nDH)
    FemaleDH = makeDH(FemaleF1, nDH)
    
    MaleYT1 = setPhenoGCA2(
      pop = MaleDH,
      testers = FemaleTester1,
      use = "pheno",
      h2 = h2YT1,
      inbred = TRUE,
      nenvs_tpe = nEnvTPE,
      nenvs_met = nenvYT1,
      nreps = nreplYT1,
      nCycles = nCycles,
      sampling = 4
    )
    
    FemaleYT1 = setPhenoGCA2(
      pop = FemaleDH,
      testers = MaleTester1,
      use = "pheno",
      h2 = h2YT1,
      inbred = TRUE,
      nenvs_tpe = nEnvTPE,
      nenvs_met = nenvYT1,
      nreps = nreplYT1,
      nCycles = nCycles,
      sampling = 4
    )
    
  }
  if(year<5){
    #Stage 3
    MaleYT2   = selectInd(MaleYT1, nInbred2, use = "pheno")
    FemaleYT2 = selectInd(FemaleYT1, nInbred2, use = "pheno")
    
    MaleYT2 = setPhenoGCA2(
      pop = MaleYT2,
      testers = FemaleTester2,
      use = "pheno",
      h2 = h2YT2,
      inbred = TRUE,
      nenvs_tpe = nEnvTPE,
      nenvs_met = nenvYT2,
      nreps = nreplYT2,
      nCycles = nCycles,
      sampling = 5
    )
    
    FemaleYT2 = setPhenoGCA2(
      pop = FemaleYT2,
      testers = MaleTester2,
      use = "pheno",
      h2 = h2YT2,
      inbred = TRUE,
      nenvs_tpe = nEnvTPE,
      nenvs_met = nenvYT2,
      nreps = nreplYT2,
      nCycles = nCycles,
      sampling = 5
    )
    
  }
  if(year<4){
    #Stage 4
    
    MaleInbredYT3   = selectInd(MaleYT2, nInbred3)
    FemaleInbredYT3 = selectInd(FemaleYT2, nInbred3)
    
    MaleHybridYT3   = hybridCross(MaleInbredYT3, FemaleElite)
    FemaleHybridYT3 = hybridCross(FemaleInbredYT3, MaleElite)
    
    MaleHybridYT3 = setPheno2(
      pop = MaleHybridYT3,
      h2 = h2YT3,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT3,
      nenvs_met = nenvYT3,
      nCycles = nCycles,
      sampling = 6
    )
    
    FemaleHybridYT3 = setPheno2(
      pop = FemaleHybridYT3,
      h2 = h2YT3,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT3,
      nenvs_met = nenvYT3,
      nCycles = nCycles,
      sampling = 6
    )
  }
  if(year<3){
    #Stage 5
    
    MaleHybridYT4   = selectInd(MaleHybridYT3, nYT4)
    FemaleHybridYT4 = selectInd(FemaleHybridYT3, nYT4)
    
    MaleHybridYT4 = setPheno2(
      pop = MaleHybridYT4,
      h2 = h2YT4,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT4,
      nenvs_met = nenvYT4,
      nCycles = nCycles,
      sampling = 7
    )
    
    FemaleHybridYT4 = setPheno2(
      pop = FemaleHybridYT4,
      h2 = h2YT4,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT4,
      nenvs_met = nenvYT4,
      nCycles = nCycles,
      sampling = 7
    )
    
    MaleInbredYT4 = MaleInbredYT3[MaleInbredYT3@id%in%MaleHybridYT4@mother]
    FemaleInbredYT4 = FemaleInbredYT3[FemaleInbredYT3@id%in%FemaleHybridYT4@mother]
    
  }
  if(year<2){
    #Stage 6
    
    MaleHybridYT5   = selectInd(MaleHybridYT4, nYT5)
    FemaleHybridYT5 = selectInd(FemaleHybridYT4, nYT5)
    
    MaleHybridYT5 = setPheno2(
      pop = MaleHybridYT5,
      h2 = h2YT5,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT5,
      nenvs_met = nenvYT5,
      nCycles = nCycles,
      sampling = 8
    )
    
    FemaleHybridYT5 = setPheno2(
      pop = FemaleHybridYT5,
      h2 = h2YT5,
      inbred = FALSE,
      nenvs_tpe = nEnvTPE,
      nreps = nreplYT5,
      nenvs_met = nenvYT5,
      nCycles = nCycles,
      sampling = 8
    )
    
    MaleInbredYT5 = MaleInbredYT4[MaleInbredYT4@id %in% MaleHybridYT5@mother]
    FemaleInbredYT5 = FemaleInbredYT4[FemaleInbredYT4@id %in% FemaleHybridYT5@mother]
  }
  # if(year<1){
  #   #Year 7
  #   Release = selectInd(pop = EYT, nInd = 1)
  #   cat(" Release", "\n")
  # }
}
