# Title: Presimulate population draws
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-14-2025

############################################################################

### Documentation for R script ###
  
# presimulate_pop_draws.r creates simulation draws for population 
# counts of the subgroups of interest based on estimates for mean and corresponding 
# standard errors. We call it "pre"-simulate because we're creating the draws before 
# doing the simulating of PIFs and impacts (the `4-pillar` results). It's useful to 
# get these sims at the beginning (as opposed to redoing them when we calculate 
# impacts for each outcome) so we can use a consistent set of sims we use for all 
# impacts and outcomes.   

############################################################################

# Create variables to match the specific naming and coding of other input files. 
# Note that "pop" is the dataframe with data on population counts. That data has 
# already been read in in	"run_models.rmd". 

# If you make changes to the age/sex/race groups (or 
# add/subtract strata), then you will need to change this part of the code to match.

pop_edit <- pop

pop_edit$agecat <- as.numeric(as.factor(pop_edit$Age))
pop_edit$female <- as.numeric(pop_edit$Sex)
pop_edit$female[as.numeric(pop_edit$Sex) == 2] <- 0
pop_edit$race <- pop_edit$Race

# add vars used for environment models
pop_edit <- pop_edit %>% mutate(subgroup_id = subgroup,
                        age_gp = case_when(agecat == 1 ~ '20-34',
                                           agecat == 2 ~ '35-44',
                                           agecat == 3 ~ '45-54',
                                           agecat == 4 ~ '55-64',
                                           agecat == 5 ~ '65-74',
                                           agecat == 6 ~ '75+'),
                        #race_gp = case_when(Race_label == 1 ~ 'NHW',
                        #                   Race_label == 2 ~ 'NHB',
                        #                   Race_label == 3 ~ 'HIS',
                        #                   Race_label == 4 ~ 'OTH'),
                        sex_gp = case_when(female == 0 ~ 'Male',
                                           female == 1 ~ 'Female')
)

# note to self: this line of code is needed for calculate_change.r to run
pop_edit$race_gp <- pop_edit$Race_label

pop_edit$age <- 0
pop_edit$age[pop_edit$agecat == 1] <- 25
pop_edit$age[pop_edit$agecat == 2] <- 35
pop_edit$age[pop_edit$agecat == 3] <- 45
pop_edit$age[pop_edit$agecat == 4] <- 55
pop_edit$age[pop_edit$agecat == 5] <- 65
pop_edit$age[pop_edit$agecat == 6] <- 75

pop_edit$race <- 0
pop_edit$race[pop_edit$Race_label == "NHW"] <- 1
pop_edit$race[pop_edit$Race_label == "NHB"] <- 2
pop_edit$race[pop_edit$Race_label == "HIS"] <- 3
pop_edit$race[pop_edit$Race_label == "OTH"] <- 4

# Loop through each row and generate population sims based on mean and standard 
# error of population for each subgroup. 

# create empty matrix
observed.pop.draws_mat <- matrix(data = NA, nrow = dim(pop_edit)[1], ncol = nsim1)

for(i in 1:dim(pop_edit)[1]) {
  
  temp <- rnorm(n = nsim1, mean = pop_edit$population[i], sd = pop_edit$population.se[i])
  temp[temp < 0] <- 0
  observed.pop.draws_mat[i,] <- temp
  
}

print(observed.pop.draws_mat)

# change to data frame
observed.pop.draws_df <- as.data.frame(observed.pop.draws_mat)

print(observed.pop.draws_df)

# change variable names to start with "X"
observed.pop.draws_df1 <- observed.pop.draws_df %>% 
  rename_with(~gsub("V", "X", .x))

# final dataset
observed.pop.draws <- cbind(pop_edit, observed.pop.draws_df1)

print(observed.pop.draws)

# Save sims as a csv file
write_csv(x = observed.pop.draws,
          file = paste0(file_location, "observed.pop.draws.csv"))
