#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
cat(" Advance year \n")

#-- Year 7

#Release hybrid
# Release <- selectInd(pop = EYT, nInd = 1)



#-- Year 6
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

gv_tpe_male <- unstr_asr_output(
  pop = MaleHybridYT5,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT5
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleHybridYT5,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT5
)

output <- save_output(
  df = output,
  malepop = MaleHybridYT5,
  femalepop = FemaleHybridYT5,
  year = year,
  stage = "YT5",
  nenvs = nenvYT5,
  gv_tpe_male = gv_tpe_male,
  gv_tpe_female = gv_tpe_female,
  GCAf = FemaleHybridYT5@pheno[,1],
  GCAm = MaleHybridYT5@pheno[,1]
)

#-- Year 5
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

gv_tpe_male <- unstr_asr_output(
  pop = MaleHybridYT4,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT4
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleHybridYT4,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT4
)

output <- save_output(
  df = output,
  malepop = MaleHybridYT4,
  femalepop = FemaleHybridYT4,
  year = year,
  stage = "YT4",
  nenvs = nenvYT4,
  gv_tpe_male = gv_tpe_male,
  gv_tpe_female = gv_tpe_female,
  GCAf = FemaleHybridYT4@pheno[,1],
  GCAm = MaleHybridYT4@pheno[,1]
)

#-- Year 4
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

gv_tpe_male <- unstr_asr_output(
  pop = MaleHybridYT3,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT3
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleHybridYT3,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT3
)

output <- save_output(
  df = output,
  malepop = MaleHybridYT3,
  femalepop = FemaleHybridYT3,
  year = year,
  stage = "YT3",
  nenvs = nenvYT3,
  gv_tpe_male = gv_tpe_male,
  gv_tpe_female = gv_tpe_female,
  GCAf = FemaleHybridYT3@pheno[,1],
  GCAm = MaleHybridYT3@pheno[,1]
)

#-- Year 3
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

gv_tpe_male <- unstr_asr_output(
  pop = MaleYT2,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT2
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleYT2,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT2
)

output <- save_output(
  df = output,
  malepop = MaleYT2,
  femalepop = FemaleYT2,
  year = year,
  stage = "YT2",
  nenvs = nenvYT2,
  gv_tpe_male = gv_tpe_male,
  gv_tpe_female = gv_tpe_female,
  GCAf = FemaleYT2@pheno[,1],
  GCAm = MaleYT2@pheno[,1]
)

#-- Year 2
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

gv_tpe_male <- unstr_asr_output(
  pop = MaleYT1,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT2
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleYT1,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT2
)

output <- save_output(
  df = output,
  malepop = MaleYT1,
  femalepop = FemaleYT1,
  year = year,
  stage = "DH/YT1",
  nenvs = nenvYT2,
  gv_tpe_male = gv_tpe_male,
  gv_tpe_female = gv_tpe_female,
  GCAf = FemaleYT1@pheno[,1],
  GCAm = MaleYT1@pheno[,1]
)

#-- Year 1
MaleF1   = randCross(MaleParents, nCrosses)
FemaleF1 = randCross(FemaleParents, nCrosses)

#------------------------
# rm(envs_met, gv_tpe, gv_met, error_met, pheno_met, mean_pheno_met)