### Simulating for hybrid breeding
#=======================================================================
##-- Simulation parameters
#=======================================================================

nReps    = 50   # Number of simulation replicates
nCycles = 20  # Number of simulation cycles
nBurnin = 20

# Trait parameters
MeanG     = 8.3   # Overall trait mean
initVarg  = 0.1905995  # Initial genetic variance
nSnp      = 1700     # Number of SNPs per chromosome
nQtl      = 300 # Number of QTLs per chromosome
k         = 7   # Number of environmental covariate terms
MeanDD    = 0.92   # mean
VarDD     = 0.3    # variance
nGenSplit = 0
nrept     = 3

# TPE parameters
# basecor = 0.2
rank = 7
skew = 0

#=======================================================================
##-- Stage-specific parameters
#=======================================================================

nEnvTPE = 1000 # Target population of environments
nStages = 7

# Parents
nParents  = 20  # Number of parents (and founders)
nCrosses  = 100 # Number of crosses per year
famMax   = 15  # The maximum number of DH lines per cross
nDH      = 50  # DH lines produced per cross

# Number of trials per stage
nenvYT1 = 1 
nenvYT2 = 3    
nenvYT3 = 10    
nenvYT4 = 20    
nenvYT5 = 50  

# Number of replicates per trial
nreplYT1 = 1    
nreplYT2 = 2    
nreplYT3 = 3    
nreplYT4 = 5   
nreplYT5 = 15

# heritability per trial
h2YT1 = 0.06
h2YT2 = 0.11
h2YT3 = 0.20
h2YT4 = 0.34
h2YT5 = 0.80

#  ----Selection on GCA ----
# Number of inbreds per heterotic pool per stage
nInbred1 = nCrosses*nDH #Do not change
nInbred2 = 400
nInbred3 = 40

# Number of testers per heterotic pool per stage
# Values must be smaller than nElite
nTester1 = 2
nTester2 = 4

# Yield trial entries
nYT1 = nInbred1*nTester1 #Do not change
nYT2 = nInbred2*nTester2 #Do not change

#  ---- Selection on SCA ----

# Elite parents per heterotic pool
nElite = 10

# Elite YT size
nYT3 = nInbred3*nElite #Do not change
nYT4 = 20
nYT5 = 7

