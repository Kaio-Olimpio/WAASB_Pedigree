# Script removes 10 oldest parents each year and replaces them 
# with most recent genotypes from the EYT stage to form a new crossing block
cat(" Update parents \n")

# Update with new 10 parents from the YT4 stage
newmale = selectInd(MaleInbredYT4,10)
newfemale = selectInd(FemaleInbredYT4,10)

tmp = matrix(NA, nrow = length(newmale@pheno), ncol = nEnvTPE)
tmp[,1] = newmale@pheno
newmale@pheno = tmp

tmp = matrix(NA, nrow = length(newfemale@pheno), ncol = nEnvTPE)
tmp[,1] = newfemale@pheno
newfemale@pheno = tmp

MaleParents   = c(MaleParents[11:nParents], newmale)
FemaleParents = c(FemaleParents[11:nParents], newfemale)

# Report parameters
gv_tpe_male <- unstr_asr_output(
  pop = MaleParents,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT4
)
gv_tpe_female <- unstr_asr_output(
  pop = FemaleParents,
  ntraits = 1,
  nenvs = nEnvTPE,
  nreps = nreplYT4
)

output <- save_output(df = output,
                      malepop = MaleParents,
                      femalepop = FemaleParents,
                      stage="Parents",
                      year = year,
                      gv_tpe_male = gv_tpe_male,
                      gv_tpe_female = gv_tpe_female)
