# Title: Presimulate observed disease burden (cancer incidence; CMD mortality) draws
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

disease <- raw.file[, vars.to.keep]

# Loop through each row and generate mortality/incidence sims based on mean and 
# standard error for each subgroup. 

observed.disease.draws_mat <- matrix(data = NA, nrow = dim(disease)[1], ncol = nsim1)

for(i in 1:dim(disease)[1]) {
  
  temp <- rnorm(n = nsim1, mean = disease$count[i], sd = disease$se[i])
  temp[temp < 0] <- 0
  observed.disease.draws_mat[i,] <- temp
  
}

print(observed.disease.draws_mat)

# change to data frame
observed.disease.draws_df <- as.data.frame(observed.disease.draws_mat)

print(observed.disease.draws_df)

# change variable names to start with "X"
observed.disease.draws_df1 <- observed.disease.draws_df %>% 
  rename_with(~gsub("V", "X", .x))

# bind
observed.disease.draws_df2 <- cbind(disease, observed.disease.draws_df1)

print(observed.disease.draws_df2)

# Duplicate the sims to use for estimating effects from mediated pathways, as 
# well as combined effects (as in, combining direct and mediated pahways). 
# Combine the duplicated sims into one file.  

medBMI.observed.disease.draws <- observed.disease.draws_df2
medBMI.observed.disease.draws$disease <- paste(observed.disease.draws_df2$disease, 
                                            "_medBMI", sep = "")

medSBP.observed.disease.draws <- observed.disease.draws_df2
medSBP.observed.disease.draws$disease <- paste(observed.disease.draws_df2$disease, 
                                            "_medSBP", sep = "")

total.observed.disease.draws <- observed.disease.draws_df2
total.observed.disease.draws$disease <- paste(observed.disease.draws_df2$disease, 
                                           "_total", sep = "")

# final dataset
observed.disease.draws <- rbind(observed.disease.draws_df2, 
                             medBMI.observed.disease.draws, 
                             medSBP.observed.disease.draws, 
                             total.observed.disease.draws)

# Save sims as a csv file. Don't be confused by the file name, it's just another 
# artifact from the days the code was used for the cancer CRA project. 
# "observed.disease.draws" contains simulation draws for morality/incidence counts 
# of the subgroups for CMD/Cancer, by subgroup of interest.

write_csv(x = observed.disease.draws, 
          file = paste0(file_location, "observed.disease.burden.draws.csv"))




