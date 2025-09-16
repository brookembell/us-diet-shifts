# Title: Presimulate population draws
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-14-2025

############################################################################

### Documentation for R script ###
  
# presimulate.pop.draws_cluster.r, which creates simulation draws for population 
# counts of the subgroups of interest based on estimates for mean and corresponding 
# standard errors. I call it "pre"-simulate because I'm creating the draws before 
# doing the simulating of PIFs and impacts (the 4 pillar results). It's useful to 
# get these sims at the beginning (as opposed to redoing them when we calculate 
# impacts for each outcome) so we can use a consistent set of sims we use for all 
# impacts and outcomes.   

############################################################################

# Create variables to match the specific naming and coding of other input files. 
# Note that "pop" is the dataframe with data on population counts. That data has 
# already been read in in	"LASTING_cluster_w_master_input.r". Also, don't be 
# confused by the "pop" being essentially renamed to "mort". This is actually just 
# legacy from previous versions of the code which created mortality draws instead 
# of population draws. I consider this code to be very "under the hood" so 
# hopefully you don't have to navigate this part of the code very often. 
# That being said, if you make changes to the age/sex/race groups (or 
# add/subtract strata), then you will need to change this part of the code to match.

mort <- pop

mort$agecat <- as.numeric(as.factor(mort$Age))
mort$female <- as.numeric(mort$Sex)
mort$female[as.numeric(mort$Sex) == 2] <- 0
mort$race <- mort$Race

# add vars used for price/environment models
mort <- mort %>% mutate(subgroup_id = subgroup,
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

mort$race_gp <- mort$Race_label

mort$age <- 0
mort$age[mort$agecat == 1] <- 25
mort$age[mort$agecat == 2] <- 35
mort$age[mort$agecat == 3] <- 45
mort$age[mort$agecat == 4] <- 55
mort$age[mort$agecat == 5] <- 65
mort$age[mort$agecat == 6] <- 75

mort$race <- 0
mort$race[mort$Race_label == "NHW"] <- 1
mort$race[mort$Race_label == "NHB"] <- 2
mort$race[mort$Race_label == "HIS"] <- 3
mort$race[mort$Race_label == "OTH"] <- 4

# Loop through each row and generate population sims based on mean and standard 
# error of population for each subgroup. 

observed.pop.draws <- matrix(data = NA, nrow = dim(mort)[1], ncol = nsim1)

print(observed.pop.draws)

for(i in 1:dim(mort)[1]) {
  
  temp <- rnorm(n = nsim1, mean = mort$population[i], sd = mort$population.se[i])
  temp[temp < 0] <- 0
  observed.pop.draws[i,] <- temp
  
}

print(observed.pop.draws)

observed.pop.draws <- cbind(mort, observed.pop.draws)

print(observed.pop.draws)

# fix variable name
observed.pop.draws <- observed.pop.draws %>% 
  rename("X1" = `1`,
         "X2" = `2`,
         "X3" = `3`)

# Save sims as a csv file.

write_csv(x = observed.pop.draws,
          file = paste0(file_location, "observed.pop.draws.csv"))



