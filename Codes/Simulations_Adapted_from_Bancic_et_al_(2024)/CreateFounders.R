##### Hybrid breeding ---------------

# Create founders
cat(" Simulating founder parents \n")

# Generate initial haplotypes
founderPop <- runMacs(nInd     = nParents * 2, 
                      nChr     = 10, 
                      segSites = nSnp + nQtl,
                      inbred   = TRUE, 
                      split    = nGenSplit,
                      species  = "MAIZE")
SP = SimParam$new(founderPop)

# Add SNP chip
SP$restrSegSites(nQtl,nSnp)
if (nSnp > 0) {
  SP$addSnpChip(nSnp)
}

# Create an input object with FieldSimR for AlphaSimR
# to obtain desired genotype slopes 
# set.seed(seed)
De = skew_diag_mat(n = nEnvTPE, shape = 3, scale = 1.2, mean.var = initVarg)

## Correlations between environments (High GEI scenario)
Ce = struc_cor_mat(
  n = nEnvTPE,
  base.cor = basecor,
  rank = rank,
  skew = skew,
  pos.def = TRUE
)

## GE matrix
Ge = sqrt(De) %*% Ce %*% sqrt(De)

input_asr <- unstr_asr_input(nenvs  = nEnvTPE, 
                             mean   = MeanG,
                             var    = diag(De), 
                             corA   = Ce, 
                             ntraits = 1,
                             varDD = VarDD,
                             meanDD = MeanDD,
                             corDD = diag(1, nrow = nEnvTPE))

# Add additive trait
SP$addTraitAD(
  nQtlPerChr = nQtl,
  mean       = input_asr$mean,
  var        = input_asr$var,
  corA       = input_asr$corA,
  meanDD = input_asr$meanDD,
  varDD = input_asr$varDD,
  corDD = input_asr$corDD,
  useVarA = FALSE
)

# Split heterotic pools to form initial parents
FemaleParents = newPop(founderPop[1:nParents])
MaleParents   = newPop(founderPop[(nParents+1):(nParents*2)])

# Set hybrid parents for later yield trials
MaleElite   = selectInd(MaleParents, nElite, use = "gv")
FemaleElite = selectInd(FemaleParents, nElite, use = "gv")

# Reverse order to keep best parent in longer
MaleElite = MaleElite[nElite:1]
FemaleElite = FemaleElite[nElite:1]

# Set initial testers for YT1 and YT2
# Requires nTesters to be smaller than nElite
MaleTester1   = MaleElite[1:nTester1]
FemaleTester1 = FemaleElite[1:nTester1]
MaleTester2   = MaleElite[1:nTester2]
FemaleTester2 = FemaleElite[1:nTester2]
rm(founderPop)
