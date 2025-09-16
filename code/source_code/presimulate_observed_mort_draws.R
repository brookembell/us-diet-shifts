# Title: Presimulate observed mortality draws
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-14-2025

############################################################################

### Documentation for R script ###

# presimulate.observed.mort.draws_cluster.r creates simulation draws for 
# morality/incidence counts of the subgroups of interest based on estimates 
# for mean and corresponding standard errors. I call it "pre"-simulate because 
# I'm creating the draws before doing the simulating of PIFs. It's useful to 
# create these sims at the beginning (as opposed to redoing them when we 
# calculate PIFs for each outcome) so we can use a consistent set of sims for all 
# health outcomes.  

############################################################################

# First, read in input file with CMD/Cancer incidence/death data, and do change 
# some some processing (change some variable names and subset to variables on 
# interest). Note that, if you decide to change the strata names or end up 
# adding/subtracing strata, you'll have to make corresponding changes to the code 
# (specfically, where you define vars.to.keep).

raw.file <- read_csv(paste0(file_location, "cvd_cancer_merged_", version.date, ".csv"))

# remove BMImed and SBPmed for now
index <- grep(pattern = "_medBMI|_medSBP", x = raw.file$diseases)
raw.file <- raw.file[-index,]

raw.file$Sex <- factor(raw.file$Sex)
raw.file$se <- raw.file$crude_se # set se to crude_se

# rename disease
raw.file <- raw.file %>% rename("disease" = "diseases")

vars.to.keep <- c("disease", "subgroup", "Sex", "Age_label", "Race_label", 
                  "count", "se")

mort <- raw.file[, vars.to.keep]

# Loop through each row and generate mortality/incidence sims based on mean and 
# standard error for each subgroup. 

observed.mort.draws <- matrix(data = NA, nrow = dim(mort)[1], ncol = nsim1)

for(i in 1:dim(mort)[1]) {
  
  temp <- rnorm(n = nsim1, mean = mort$count[i], sd = mort$se[i])
  temp[temp < 0] <- 0
  observed.mort.draws[i,] <- temp
  
}

observed.mort.draws <- cbind(mort, observed.mort.draws)

# Duplicate the sims to use for estimating effects from mediated pathways, as 
# well as combined effects (as in, combining direct and mediated pahways). 
# Combine the duplicated sims into one file.  

medBMI.observed.mort.draws <- observed.mort.draws
medBMI.observed.mort.draws$disease <- paste(observed.mort.draws$disease, 
                                            "_medBMI", sep = "")

medSBP.observed.mort.draws <- observed.mort.draws
medSBP.observed.mort.draws$disease <- paste(observed.mort.draws$disease, 
                                            "_medSBP", sep = "")

total.observed.mort.draws <- observed.mort.draws
total.observed.mort.draws$disease <- paste(observed.mort.draws$disease, 
                                           "_total", sep = "")

observed.mort.draws <- rbind(observed.mort.draws, 
                             medBMI.observed.mort.draws, 
                             medSBP.observed.mort.draws, 
                             total.observed.mort.draws)

# fix variable name
observed.mort.draws <- observed.mort.draws %>% 
  rename("X1" = `1`,
         "X2" = `2`,
         "X3" = `3`)

# Save sims as a csv file. Don't be confused by the file name, it's just another 
# artifact from the days the code was used for the cancer CRA project. 
# "observed.mort.draws" contains simulation draws for morality/incidence counts 
# of the subgroups for CMD/Cancer, by subgroup of interest.

write_csv(x = observed.mort.draws, 
          file = paste0(file_location, "observed.cancer.mortality.draws.csv"))




