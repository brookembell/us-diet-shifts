# Title: Combine Food-At-Home (FAH) and Food-Away-From-Home (FAFH) output
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-14-2025

############################################################################

### Documentation for R script ###

# combine.FAH.FAFH_cluster.r is used to combine food at home (FAH) and food away 
# from home (FAFH) results. This script is called in "run_models.r" to combine 
# these two results which are generated separately in LASTING model. To better 
# understand the context, in which this script is being used, it may be useful to 
# review the section in "run_models.r" that calls this script.

############################################################################

# Load necessary libraries and read in FAH and FAFH output files (all simulated 
# draws) for the diet pattern and outcome type of interest (defined in 
# "LASTING_cluster_w_master_input.r" where this script is called). Give each file 
# indicator variable "X" so we can differentiate between FAH and FAFH results 
# after merging.

library(data.table)

FAH.output <- read_csv(paste0("output/envecosoc/", 
                              diet, "_diet_Gro/", 
                              output.type.sims, ".csv"))


FAFH.output <- read_csv(paste0("output/envecosoc/", 
                               diet, "_diet_Oth/", 
                               output.type.sims, ".csv"))

FAH.output$X <- "FAH"
FAFH.output$X <- "FAFH"

# Merge FAH and FAFH sims files, sum FAH and FAFH sims, and save results.

combined.output <- rbind(FAH.output, FAFH.output)
cols.to.sum <- paste("X", 1:n.sims, sep = "")
sims.data.table <- as.data.table(combined.output)

keep <- c("Foodgroup", "sex_gp", "race_gp","age_gp", "subgroup_id", "population",
          "Intake_unit", "outcome", "outcome_unit")

all.output <- sims.data.table[,lapply(.SD, sum), by = keep, .SDcols = cols.to.sum]

# first create directory if it doesn't already exist
ifelse(!dir.exists(file.path(paste0("output/envecosoc/", diet, "_diet_both"))),
       dir.create(file.path(paste0("output/envecosoc/", diet, "_diet_both"))),
       "Directory Exists")

# export
write_csv(x = all.output, 
          file = paste0("output/envecosoc/", diet, 
                      "_diet_both/", output.type.sims, 
                      ".csv"))

# Read in pop sims, use it to calculate per capita output for summed sims 
# (and then save output).

pop.sims <- read_csv(paste0(file_location, "observed.pop.draws.csv"))

percapita.all <- 
  get.percapita.sims(impact.sims = as.data.frame(all.output), 
                     population.sims = pop.sims, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))


# first create directory if it doesn't already exist
ifelse(!dir.exists(file.path(paste0("output/envecosoc/", diet, 
                                    "_diet_both/per_capita"))),
       dir.create(file.path(paste0("output/envecosoc/", diet, 
                                   "_diet_both/per_capita"))),
       "Directory Exists")

fwrite(x = all.output, 
       file = paste0("output/envecosoc/", diet, 
                   "_diet_both/per_capita/", 
                   output.type.sims, 
                   ".percapita.csv"))

# Calculate and save summary stats (mean, median, sd, 2.5th and 9.5th percentiles) from sims. 

combined.summary <- get.summary.stats.simple(impact.sims = as.data.frame(all.output), 
                                             vars = paste("X", 1:n.sims, sep = ""))

write_csv(x = combined.summary, 
          file = paste0("output/envecosoc/", diet, 
                      "_diet_both/", output.type.summary, ".csv"))

combined.summary.percapita <- 
  get.summary.stats.simple(impact.sims = as.data.frame(percapita.all), 
                           vars = paste("X", 1:n.sims, sep = ""))

write_csv(x = combined.summary.percapita, 
          file = paste0("output/envecosoc/", diet, 
                      "_diet_both/per_capita/", output.type.summary, 
                      ".percapita.csv"))

# Next, we want to get these results for all strata combos. First, we define the 
# vector "strata" which has strings of all the strata you are using for your 
# analysis. Note that, if the names of your strata change or you are using your 
# own strata, you need to change this vector to match your analysis. 

strata <- c("age_gp", "sex_gp", "race_gp", "Foodgroup")

# The following code gets all the strata combos of interest in the a list object 
# named "strata.combos", based on strata vector defined above. Here, "i" is a 
# possible number of a strata that a strata-combination is made out of. In our 
# case, we have 4 possibilities. 1-strata "combos" like age, sex, etc.., 2-strata 
# combos like age/sex, age/race, etc.., 3-strata combos like age/sex/race etc... 
# and one four strata combo: age/sex/race/Foodgroup. So, You loop through from 
# i =1 to i=4, and at each iteration, you get all possible strata combinations 
# with i strata using combn function. Then, you add  to strata.combos.list. 
# Finally, you want add a "null" element to your list so you when you loop through 
# strata.combos, you also create results that sum over all subgroups.

strata.combos <- list()

for(i in 1:length(strata)){
  
  x <- combn(strata, i)
  combos.for.i <- split(x, rep(1:ncol(x), each = nrow(x)))
  strata.combos <- c(strata.combos, combos.for.i)
  
}

# add null element to list to use for overall numbers
strata.combos["null"] <- list(NULL) 

# Use summary.stats.by.strata function defined in "cost_change_functions.r" to get 
# summary stats for all strata combinations, which takes the strata.combos list 
# you just created and the FAH+FAFH combined sims as input. summ.stats.by.strata.cost 
# will be a list of 4 lists, each corresponding to an output type: 1 = summary stats, 
# 2 = sims, 3 = per capita summary stats, 4 = per capita sims. Each list is a list 
# with an element for each strata combo. 

output_location <- paste0("output/envecosoc/", diet, "_diet_both/")

summ.stats.by.strata.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = all.output, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

# Next, loop through the strata combos and save the sim output, rename variables 
# names for the summary output and include them in summary.list and 
# summary.list.percapita created in "LASTING_cluster_w_masterinput.r". 
# Finally, Some fine-tuning on names of the list elements (adding underscores 
# between strata, as well as adding "all" label to last element of lists which 
# corresponds the "null" strata combo summing over all subgroups).

# first create directory if it doesn't already exist
ifelse(!dir.exists(file.path(paste0("output/envecosoc/", 
                                    diet, "_diet_both/By_Subgroup"))),
       dir.create(file.path(paste0("output/envecosoc/", 
                                   diet, "_diet_both/By_Subgroup"))),
       "Directory Exists")

ifelse(!dir.exists(file.path(paste0("output/envecosoc/", 
                                    diet, "_diet_both/By_Subgroup/full_sims"))),
       dir.create(file.path(paste0("output/envecosoc/", 
                                   diet, "_diet_both/By_Subgroup/full_sims"))),
       "Directory Exists")

ifelse(!dir.exists(file.path(paste0("output/envecosoc/", 
                                    diet, "_diet_both/per_capita/By_SubGroup"))),
       dir.create(file.path(paste0("output/envecosoc/", 
                                   diet, "_diet_both/per_capita/By_SubGroup"))),
       "Directory Exists")

ifelse(!dir.exists(file.path(paste0("output/envecosoc/", 
                                    diet, "_diet_both/per_capita/By_SubGroup/full_sims"))),
       dir.create(file.path(paste0("output/envecosoc/", 
                                   diet, "_diet_both/per_capita/By_SubGroup/full_sims"))),
       "Directory Exists")

for(j in 1:length(strata.combos)){
  
  fwrite(x = summ.stats.by.strata.envecosoc[[2]][[j]], 
         file = paste0(output_location, "By_SubGroup/full_sims/", 
                     output.type.sims.bystrata,".sims.output_by_", 
                     paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.envecosoc[[4]][[j]], 
         file = paste0(output_location, "per_capita/By_SubGroup/full_sims/", 
                     output.type.sims.bystrata, ".sims.percapita.output_by_", 
                     paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  names(
    summ.stats.by.strata.envecosoc[[1]][[j]])[which(names(
      summ.stats.by.strata.envecosoc[[1]][[j]]) %in% c("lower_bound (2.5th percentile)", 
                                                  "median", 
                                                  "upper_bound  (97.5th percentile)", 
                                                  "mean", 
                                                  "SD"))] <- 
    c(LB.name, median.name, UB.name, mean.name, SD.name)
  
  names(
    summ.stats.by.strata.envecosoc[[3]][[j]])[which(names(
      summ.stats.by.strata.envecosoc[[3]][[j]]) %in% c("lower_bound (2.5th percentile)", 
                                                  "median", 
                                                  "upper_bound  (97.5th percentile)", 
                                                  "mean", 
                                                  "SD"))] <- 
    c(LB.name, median.name, UB.name, mean.name, SD.name)
  
  summary.list[[k]][[j]] <- summ.stats.by.strata.envecosoc[[1]][[j]]
  
  summary.list.percapita[[k]][[j]] <- summ.stats.by.strata.envecosoc[[3]][[j]]
  
}

names(summary.list[[k]]) <- lapply(strata.combos, paste, collapse = "_")
names(summary.list[[k]])[length(summary.list[[k]])] <- "all"
names(summary.list.percapita[[k]]) <- lapply(strata.combos, paste, collapse = "_")
names(summary.list.percapita[[k]])[length(summary.list[[k]])] <- "all"

# Note that we're not saving the summary stats output just yet. There is more 
# processing done in "run_models.rmd" after this script is called.

