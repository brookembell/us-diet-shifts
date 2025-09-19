# Title: Calculate change
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-14-2025

############################################################################

### Documentation for R script ###

# Calculate_change.r  is called by "run_models.rmd" and is used to calculate 
# the environmental, economic, and social ("env-eco-soc") impacts for a given 
# dietary type and food source (food at home vs food away from home). It will 
# calculate population-level and per capita results, and output the full sims, as 
# well as summary stats, into separate output files.

############################################################################

# First, read in functions used to calculate impacts.

source(paste0(code_location, "functions/change_functions.r"))

# Define output location.

output_location <- paste0("output/envecosoc/", envecosoc_output_loc.vec[k], "/")

print(output_location)

# Read in relevant inputs (population draws, env-eco-soc impact factors of dietary 
# factors) and do some minor data processing.

pop.sims <- read_csv(paste0(file_location, "observed.pop.draws.csv"))

envecosoc.inputs <- envecosoc_inputs.vec[[k]]

merged <- merge(envecosoc.inputs, pop.sims, by = "subgroup_id")
merged$Sex <- merged$Sex.y
merged$population <- merged$population.y

# Create empty list objects to store outputs. Each element of the lists corresponds to 
# a specific subgroup. There are 31 lists in total, corresponding to some combination 
# of (1) what we are looking at (impact (as in, change in impact), substitution impact, 
# combined impact (change in impact + substitution impact), population, delta 
# (aka, change in intake), current intake impact, counterfactual intake impact)*, 
# and (2) level of disaggregation (total, consumed, unconsumed, inedible, wasted). 
# Hopefully, what each list corresponds to is clear from name.

# See "change_functions.r" for more detail on what these things mean

list.of.impact.sims.total <- list()
list.of.substitution.impact.sims.total <- list()
list.of.combined.impact.sims.total <- list()
list.of.population.sims <- list()
list.of.delta.sims.total <- list()

list.of.impact.sims.consumed <- list()
list.of.substitution.impact.sims.consumed <- list()
list.of.combined.impact.sims.consumed <- list()
list.of.delta.sims.consumed <- list()

list.of.impact.sims.unconsumed <- list()
list.of.substitution.impact.sims.unconsumed <- list()
list.of.combined.impact.sims.unconsumed <- list()
list.of.delta.sims.unconsumed <- list()

list.of.current.intake.impact.sims.total <- list()
list.of.CF.intake.impact.sims.total <- list()
list.of.current.intake.impact.sims.consumed <- list()
list.of.CF.intake.impact.sims.consumed <- list()
list.of.current.intake.impact.sims.unconsumed <- list()
list.of.CF.intake.impact.sims.unconsumed <- list()

list.of.impact.sims.inedible <- list()
list.of.substitution.impact.sims.inedible <- list()
list.of.combined.impact.sims.inedible <- list()
list.of.delta.sims.inedible <- list()

list.of.impact.sims.wasted <- list()
list.of.substitution.impact.sims.wasted <- list()
list.of.combined.impact.sims.wasted <- list()
list.of.delta.sims.wasted <- list()

list.of.current.intake.impact.sims.inedible <- list()
list.of.CF.intake.impact.sims.inedible <- list()
list.of.current.intake.impact.sims.wasted <- list()
list.of.CF.intake.impact.sims.wasted <- list()

# Rename some variables.

n.sims <- nsim1
impact.factor.names <- envecosoc.outcomes.vec

# Start for-loop. We are going to be traversing the "merged" file row-by-row and 
# calculating impacts. Each row should correspond to a subgroup and dietary factor.

# i=1

for(i in 1:dim(merged)[1]){ # for each row in the input file (corresponding to fcid for a given subgroup)
  
  # Extract label variables and put them into their own data frames 
  # (to be merged later with the results).
  
  print(i)
  
  vars.to.keep.for.df.a <- c("Foodgroup", "sex_gp", "race_gp", "age_gp",
                             "subgroup_id", "population", 
                             "Mean_Intake", "SE_Intake", "CF_Mean_Intake", 
                             "CF_SE_Intake", "Intake_unit", 
                             paste0(all.outcomes.vec, "_unit"),
                             paste0(all.outcomes.vec, "_to_intake_conversion"))
  
  df.a <- merged[i, vars.to.keep.for.df.a] 
  
  # Here, we are just creating a data frame that repeats df.a for however many 
  # outcomes we are looking at (six, in our case). 
  # We are going to have results for each pair of dietary factor and outcome.  
  
  df.aa <- df.a[rep(seq_len(nrow(df.a)), each = length(impact.factor.names)),]
  
  # For the sake of readability of code, we extract various relevant inputs 
  # from merged[i,] and give them their own object with own name.
  
  current.mean <- merged[i, "Mean_Intake"]
  current.se <- merged[i, "SE_Intake"]
  counterfactual.mean <- merged[i, "CF_Mean_Intake"]
  counterfactual.se <- merged[i, "CF_SE_Intake"]
  impact.factor.means <- merged[i, envecosoc.outcomes.vec.mn]
  impact.factor.ses <- merged[i, envecosoc.outcomes.vec.se]
  
  RRunit <- merged %>% 
    select(c(paste0(impact.factor.names, "_to_intake_conversion"))) %>% 
    slice(1)
  
  pop.sims.names <- paste("X", c(1:n.sims), sep = "")
  
  population.sims <- matrix(rep(x = as.numeric(merged[i, pop.sims.names]), 
                                times = length(impact.factor.names)),
                            nrow = length(impact.factor.names))
  
  print(population.sims)
  
  substitution.impact.factor.means <- merged[i, envecosoc.outcomes.substitution.vec.mn]
  substitution.impact.factor.ses <- merged[i, envecosoc.outcomes.substitution.vec.se]
  substitution_unit <- merged[i, "substitution_unit"] 
  
  current_inedible_p <- merged[i, envecosoc.current.inedible.p.mn]
  current_inedible_p_se <- merged[i, envecosoc.current.inedible.p.se]
  current_foodwaste_p <- merged[i, envecosoc.current.foodwaste.p.mn]
  current_foodwaste_p_se <- merged[i, envecosoc.current.foodwaste.p.se]
  counterfactual_inedible_p <- merged[i, envecosoc.counterfactual.inedible.p.mn]
  counterfactual_inedible_p_se <- merged[i, envecosoc.counterfactual.inedible.p.se]
  counterfactual_foodwaste_p <- merged[i, envecosoc.counterfactual.foodwaste.p.mn]
  counterfactual_foodwaste_p_se <- merged[i, envecosoc.counterfactual.foodwaste.p.se]
  
  # We take the above inputs and feed it into a function that simulates the impact 
  # (simulate.impact). The resulting "b" is a list of 31 matrices with simulation 
  # results that correspond to the 31 lists we made earlier.
  
  b <- simulate.impact(current.mean = current.mean, 
                       current.se = current.se, 
                       counterfactual.mean = counterfactual.mean,
                       counterfactual.se = counterfactual.se, 
                       impact.factor.names = impact.factor.names, 
                       impact.factor.means = as.numeric(impact.factor.means), 
                       impact.factor.ses = as.numeric(impact.factor.ses),
                       RRunit = RRunit, 
                       substitution.impact.factor.means = as.numeric(substitution.impact.factor.means),
                       substitution.impact.factor.ses = as.numeric(substitution.impact.factor.ses), 
                       substitution_unit = substitution_unit, 
                       current_inedible_p = as.numeric(current_inedible_p),
                       current_inedible_p_se = as.numeric(current_inedible_p_se),
                       current_foodwaste_p = as.numeric(current_foodwaste_p),
                       current_foodwaste_p_se = as.numeric(current_foodwaste_p_se),
                       counterfactual_inedible_p = as.numeric(counterfactual_inedible_p),
                       counterfactual_inedible_p_se = as.numeric(counterfactual_inedible_p_se),
                       counterfactual_foodwaste_p = as.numeric(counterfactual_foodwaste_p),
                       counterfactual_foodwaste_p_se = as.numeric(counterfactual_foodwaste_p_se),
                       n.sims = n.sims,
                       population.sims = population.sims)
  
  # Extract results from "b" (which itself is a list of 31 matrices) and 
  # create data frames with proper labels.
  
  df.b1 <- data.frame(outcome = rownames(b[[1]]), b[[1]])
  df.b2 <- data.frame(outcome = rownames(b[[2]]), b[[2]])
  df.b3 <- data.frame(outcome = rownames(b[[3]]), b[[3]])
  df.b4 <- data.frame(outcome = rownames(b[[4]]), b[[4]])
  df.b5 <- data.frame(outcome = rownames(b[[5]]), b[[5]])
  df.b6 <- data.frame(outcome = rownames(b[[6]]), b[[6]])
  df.b7 <- data.frame(outcome = rownames(b[[7]]), b[[7]])
  df.b8 <- data.frame(outcome = rownames(b[[8]]), b[[8]])
  df.b9 <- data.frame(outcome = rownames(b[[9]]), b[[9]])
  df.b10 <- data.frame(outcome = rownames(b[[10]]), b[[10]])
  df.b11 <- data.frame(outcome = rownames(b[[11]]), b[[11]])
  df.b12 <- data.frame(outcome = rownames(b[[12]]), b[[12]])
  df.b13 <- data.frame(outcome = rownames(b[[13]]), b[[13]])
  df.b14 <- data.frame(outcome = rownames(b[[14]]), b[[14]])
  df.b15 <- data.frame(outcome = rownames(b[[15]]), b[[15]])
  df.b16 <- data.frame(outcome = rownames(b[[16]]), b[[16]])
  df.b17 <- data.frame(outcome = rownames(b[[17]]), b[[17]])
  df.b18 <- data.frame(outcome = rownames(b[[18]]), b[[18]])
  df.b19 <- data.frame(outcome = rownames(b[[19]]), b[[19]])
  df.b20 <- data.frame(outcome = rownames(b[[20]]), b[[20]])
  df.b21 <- data.frame(outcome = rownames(b[[21]]), b[[21]])
  df.b22 <- data.frame(outcome = rownames(b[[22]]), b[[22]])
  df.b23 <- data.frame(outcome = rownames(b[[23]]), b[[23]])
  df.b24 <- data.frame(outcome = rownames(b[[24]]), b[[24]])
  df.b25 <- data.frame(outcome = rownames(b[[25]]), b[[25]])
  df.b26 <- data.frame(outcome = rownames(b[[26]]), b[[26]])
  df.b27 <- data.frame(outcome = rownames(b[[27]]), b[[27]])
  df.b28 <- data.frame(outcome = rownames(b[[28]]), b[[28]])
  df.b29 <- data.frame(outcome = rownames(b[[29]]), b[[29]])
  df.b30 <- data.frame(outcome = rownames(b[[30]]), b[[30]])
  df.b31 <- data.frame(outcome = rownames(b[[31]]), b[[31]])
  
  # Bind df.aa to the df.bs. So now, the df.cs have the relevant variables 
  # for interpreting/labeling results (foodgroup, age/sex/race, etc..., all 
  # the variables listed in vars.to.keep.for.df.a)
  
  df.c1 <- cbind(df.aa, df.b1)
  df.c2 <- cbind(df.aa, df.b2)
  df.c3 <- cbind(df.aa, df.b3)
  df.c4 <- cbind(df.aa, df.b4)
  df.c5 <- cbind(df.aa, df.b5)
  df.c6 <- cbind(df.aa, df.b6)
  df.c7 <- cbind(df.aa, df.b7)
  df.c8 <- cbind(df.aa, df.b8)
  df.c9 <- cbind(df.aa, df.b9)
  df.c10 <- cbind(df.aa, df.b10)
  df.c11 <- cbind(df.aa, df.b11)
  df.c12 <- cbind(df.aa, df.b12)
  df.c13 <- cbind(df.aa, df.b13)
  df.c14 <- cbind(df.aa, df.b14)
  df.c15 <- cbind(df.aa, df.b15)
  df.c16 <- cbind(df.aa, df.b16)
  df.c17 <- cbind(df.aa, df.b17)
  df.c18 <- cbind(df.aa, df.b18)
  df.c19 <- cbind(df.aa, df.b19)
  df.c20 <- cbind(df.aa, df.b20)
  df.c21 <- cbind(df.aa, df.b21)
  df.c22 <- cbind(df.aa, df.b22)
  df.c23 <- cbind(df.aa, df.b23)
  df.c24 <- cbind(df.aa, df.b24)
  df.c25 <- cbind(df.aa, df.b25)
  df.c26 <- cbind(df.aa, df.b26)
  df.c27 <- cbind(df.aa, df.b27)
  df.c28 <- cbind(df.aa, df.b28)
  df.c29 <- cbind(df.aa, df.b29)
  df.c30 <- cbind(df.aa, df.b30)
  df.c31 <- cbind(df.aa, df.b31)
  
  # Create "impact_unit_strings" string vector, which is just the "outcome" 
  # variable with "_unit" appended to it. Then, create "outcome_unit" variable 
  # for all df.cs. Then, for each df.c, loop through each row and populate 
  # "outcome_unit" with the the corresponding unit for that outcome for that row 
  # (which is extracted from the row j and column impact_unit_strings[j] for the df.c).    
  
  impact_unit_strings <- paste(df.b1$outcome, "_unit", sep = "")
  
  df.c1$outcome_unit <- NA
  df.c2$outcome_unit <- NA
  df.c3$outcome_unit <- NA
  df.c4$outcome_unit <- NA
  df.c5$outcome_unit <- NA
  df.c6$outcome_unit <- NA
  df.c7$outcome_unit <- NA
  df.c8$outcome_unit <- NA
  df.c9$outcome_unit <- NA
  df.c10$outcome_unit <- NA
  df.c11$outcome_unit <- NA
  df.c12$outcome_unit <- NA
  df.c13$outcome_unit <- NA
  df.c14$outcome_unit <- NA
  df.c15$outcome_unit <- NA
  df.c16$outcome_unit <- NA
  df.c17$outcome_unit <- NA
  df.c18$outcome_unit <- NA
  df.c19$outcome_unit <- NA
  df.c20$outcome_unit <- NA
  df.c21$outcome_unit <- NA
  df.c22$outcome_unit <- NA
  df.c23$outcome_unit <- NA
  df.c24$outcome_unit <- NA
  df.c25$outcome_unit <- NA
  df.c26$outcome_unit <- NA
  df.c27$outcome_unit <- NA
  df.c28$outcome_unit <- NA
  df.c29$outcome_unit <- NA
  df.c30$outcome_unit <- NA
  df.c31$outcome_unit <- NA
  
  for(j in 1:dim(df.c1)[1]){
    
    df.c1$outcome_unit[j] <- df.c1[j,impact_unit_strings[j]]
    df.c2$outcome_unit[j] <- df.c2[j,impact_unit_strings[j]]
    df.c3$outcome_unit[j] <- df.c3[j,impact_unit_strings[j]]
    df.c4$outcome_unit[j] <- df.c4[j,impact_unit_strings[j]]
    df.c5$outcome_unit[j] <- df.c5[j,impact_unit_strings[j]]
    df.c6$outcome_unit[j] <- df.c6[j,impact_unit_strings[j]]
    df.c7$outcome_unit[j] <- df.c7[j,impact_unit_strings[j]]
    df.c8$outcome_unit[j] <- df.c8[j,impact_unit_strings[j]]
    df.c9$outcome_unit[j] <- df.c9[j,impact_unit_strings[j]]
    df.c10$outcome_unit[j] <- df.c10[j,impact_unit_strings[j]]
    df.c11$outcome_unit[j] <- df.c11[j,impact_unit_strings[j]]
    df.c12$outcome_unit[j] <- df.c12[j,impact_unit_strings[j]]
    df.c13$outcome_unit[j] <- df.c13[j,impact_unit_strings[j]]
    df.c14$outcome_unit[j] <- df.c14[j,impact_unit_strings[j]]
    df.c15$outcome_unit[j] <- df.c15[j,impact_unit_strings[j]]
    df.c16$outcome_unit[j] <- df.c16[j,impact_unit_strings[j]]
    df.c17$outcome_unit[j] <- df.c17[j,impact_unit_strings[j]]
    df.c18$outcome_unit[j] <- df.c18[j,impact_unit_strings[j]]
    df.c19$outcome_unit[j] <- df.c19[j,impact_unit_strings[j]]  
    df.c20$outcome_unit[j] <- df.c20[j,impact_unit_strings[j]]
    df.c21$outcome_unit[j] <- df.c21[j,impact_unit_strings[j]]
    df.c22$outcome_unit[j] <- df.c22[j,impact_unit_strings[j]]
    df.c23$outcome_unit[j] <- df.c23[j,impact_unit_strings[j]]
    df.c24$outcome_unit[j] <- df.c24[j,impact_unit_strings[j]]
    df.c25$outcome_unit[j] <- df.c25[j,impact_unit_strings[j]]
    df.c26$outcome_unit[j] <- df.c26[j,impact_unit_strings[j]]
    df.c27$outcome_unit[j] <- df.c27[j,impact_unit_strings[j]]
    df.c28$outcome_unit[j] <- df.c28[j,impact_unit_strings[j]]
    df.c29$outcome_unit[j] <- df.c29[j,impact_unit_strings[j]]
    df.c30$outcome_unit[j] <- df.c30[j,impact_unit_strings[j]]
    df.c31$outcome_unit[j] <- df.c31[j,impact_unit_strings[j]]
    
    }
  
  # Move the new "outcome_unit" column to right after "outcome", and get rid of 
  # the impact factor specific columns for units.
  
  df.c1 <- df.c1 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c2 <- df.c2 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c3 <- df.c3 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c4 <- df.c4 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c5 <- df.c5 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c6 <- df.c6 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c7 <- df.c7 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c8 <- df.c8 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c9 <- df.c9 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c10 <- df.c10 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c11 <- df.c11 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c12 <- df.c12 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c13 <- df.c13 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c14 <- df.c14 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c15 <- df.c15 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c16 <- df.c16 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c17 <- df.c17 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c18 <- df.c18 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c19 <- df.c19 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c20 <- df.c20 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c21 <- df.c21 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c22 <- df.c22 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c23 <- df.c23 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c24 <- df.c24 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c25 <- df.c25 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c26 <- df.c26 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c27 <- df.c27 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c28 <- df.c28 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c29 <- df.c29 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c30 <- df.c30 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  df.c31 <- df.c31 %>% relocate(outcome_unit, .after = outcome) %>% 
    select(-impact_unit_strings)
  
  # We assign each data frame we made to be the i-th element of their 
  # corresponding lists. Recall, earlier in the code, we made 31 different lists 
  # to store output (one for each kind of output type). We are now populating 
  # that list, the i-th element corresponding to the appropriate df.c. (e.g the 
  # current version of df.c1 for this iteration within the for loop is being 
  # stored as the i-th element of list.of.impact.sims.total, and so on).
  
  # total
  list.of.impact.sims.total[[i]] <- df.c1
  list.of.substitution.impact.sims.total[[i]] <- df.c2
  list.of.combined.impact.sims.total[[i]] <- df.c3
  list.of.population.sims[[i]] <- df.c4
  list.of.delta.sims.total[[i]] <- df.c5
  
  # consumed
  list.of.impact.sims.consumed[[i]] <- df.c6
  list.of.substitution.impact.sims.consumed[[i]] <- df.c7
  list.of.combined.impact.sims.consumed[[i]] <- df.c8
  list.of.delta.sims.consumed[[i]] <- df.c9
  
  # unconsumed
  list.of.impact.sims.unconsumed[[i]] <- df.c10
  list.of.substitution.impact.sims.unconsumed[[i]] <- df.c11
  list.of.combined.impact.sims.unconsumed[[i]] <- df.c12
  list.of.delta.sims.unconsumed[[i]] <- df.c13
  
  list.of.current.intake.impact.sims.total[[i]] <- df.c14
  list.of.CF.intake.impact.sims.total[[i]] <- df.c15
  list.of.current.intake.impact.sims.consumed[[i]] <- df.c16
  list.of.CF.intake.impact.sims.consumed[[i]] <- df.c17
  list.of.current.intake.impact.sims.unconsumed[[i]] <- df.c18
  list.of.CF.intake.impact.sims.unconsumed[[i]] <- df.c19
  
  # inedible
  list.of.impact.sims.inedible[[i]] <- df.c20
  list.of.substitution.impact.sims.inedible[[i]] <- df.c21
  list.of.combined.impact.sims.inedible[[i]] <- df.c22
  list.of.delta.sims.inedible[[i]] <- df.c23
  
  # wasted
  list.of.impact.sims.wasted[[i]] <- df.c24
  list.of.substitution.impact.sims.wasted[[i]] <- df.c25
  list.of.combined.impact.sims.wasted[[i]] <- df.c26
  list.of.delta.sims.wasted[[i]] <- df.c27
  
  list.of.current.intake.impact.sims.inedible[[i]] <- df.c28
  list.of.CF.intake.impact.sims.inedible[[i]] <- df.c29
  list.of.current.intake.impact.sims.wasted[[i]] <- df.c30
  list.of.CF.intake.impact.sims.wasted[[i]] <- df.c31
  
} # end of loop

# For each list, row-bind the elements of the list (corresponding to a subgroup 
# and dietary factor) together into a data frame. Each of these data frames has 
# results for all dietary factors and subgroups, for one outcome (specified in 
# name of data frame). NAs are also replaced with 0. 

# Recall there are 31 lists,each corresponding to some combination of (1) what we 
# are looking at (impact (as in, impact of change), substitution impact, combined 
# impact, population, delta (aka, change in intake), current intake impact, 
# counterfactual intake impact), and (2) level of disaggregation (total, consumed, 
# unconsumed, inedible, wasted).

impact.sims.total.allgroups <- 
  Reduce(rbind,list.of.impact.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

substitution.impact.sims.total.allgroups <- 
  Reduce(rbind,list.of.substitution.impact.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

combined.impact.sims.total.allgroups <- 
  Reduce(rbind,list.of.combined.impact.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

population.sims.allgroups <- 
  Reduce(rbind,list.of.population.sims) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

delta.sims.total.allgroups <- 
  Reduce(rbind,list.of.delta.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

impact.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.impact.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

substitution.impact.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.substitution.impact.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

combined.impact.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.combined.impact.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

delta.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.delta.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

impact.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.impact.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

substitution.impact.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.substitution.impact.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

combined.impact.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.combined.impact.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

delta.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.delta.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

current.intake.impact.sims.total.allgroups <- 
  Reduce(rbind,list.of.current.intake.impact.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

CF.intake.impact.sims.total.allgroups <- 
  Reduce(rbind,list.of.CF.intake.impact.sims.total) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

current.intake.impact.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.current.intake.impact.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

CF.intake.impact.sims.consumed.allgroups <- 
  Reduce(rbind,list.of.CF.intake.impact.sims.consumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

current.intake.impact.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.current.intake.impact.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

CF.intake.impact.sims.unconsumed.allgroups <- 
  Reduce(rbind,list.of.CF.intake.impact.sims.unconsumed) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

impact.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.impact.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

substitution.impact.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.substitution.impact.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

combined.impact.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.combined.impact.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

delta.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.delta.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

impact.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.impact.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

substitution.impact.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.substitution.impact.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

combined.impact.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.combined.impact.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

delta.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.delta.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

current.intake.impact.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.current.intake.impact.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

CF.intake.impact.sims.inedible.allgroups <- 
  Reduce(rbind,list.of.CF.intake.impact.sims.inedible) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

current.intake.impact.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.current.intake.impact.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

CF.intake.impact.sims.wasted.allgroups <- 
  Reduce(rbind,list.of.CF.intake.impact.sims.wasted) %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

# Save dataframes into csv files.

# first create directory if it doesn't already exist
ifelse(!dir.exists(file.path(output_location)),
       dir.create(file.path(output_location)),
       "Directory Exists")

# export as CSVs

fwrite(x = impact.sims.total.allgroups, 
       file = paste0(output_location, "envecosoc.sims.granular.csv"))

fwrite(x = substitution.impact.sims.total.allgroups, 
       file = paste0(output_location, "substitution.envecosoc.sims.granular.csv"))

fwrite(x = combined.impact.sims.total.allgroups, 
       file = paste0(output_location, "combined.envecosoc.sims.granular.csv"))

fwrite(x = delta.sims.total.allgroups, 
       file = paste0(output_location, "delta.envecosoc.sims.granular.csv"))

fwrite(x = impact.sims.consumed.allgroups, 
       file = paste0(output_location, "envecosoc.sims.granular.consumed.csv"))

fwrite(x = substitution.impact.sims.consumed.allgroups, 
       file = paste0(output_location, "substitution.envecosoc.sims.granular.consumed.csv"))

fwrite(x = combined.impact.sims.consumed.allgroups, 
       file = paste0(output_location, "combined.envecosoc.sims.granular.consumed.csv"))

fwrite(x = impact.sims.unconsumed.allgroups, 
       file = paste0(output_location, "envecosoc.sims.granular.unconsumed.csv"))

fwrite(x = substitution.impact.sims.unconsumed.allgroups, 
       file = paste0(output_location, "substitution.envecosoc.sims.granular.unconsumed.csv"))

fwrite(x = combined.impact.sims.unconsumed.allgroups, 
       file = paste0(output_location, "combined.envecosoc.sims.granular.unconsumed.csv"))

fwrite(x = current.intake.impact.sims.total.allgroups, 
       file = paste0(output_location, "current.intake.envecosoc.sims.granular.csv"))

fwrite(x = CF.intake.impact.sims.total.allgroups, 
       file = paste0(output_location, "CF.envecosoc.intake.sims.granular.csv"))

fwrite(x = current.intake.impact.sims.consumed.allgroups, 
       file = paste0(output_location, "current.intake.envecosoc.sims.granular.consumed.csv"))

fwrite(x = CF.intake.impact.sims.consumed.allgroups, 
       file = paste0(output_location, "CF.intake.envecosoc.sims.granular.consumed.csv"))

fwrite(x = current.intake.impact.sims.unconsumed.allgroups, 
       file = paste0(output_location, "current.intake.envecosoc.sims.granular.unconsumed.csv"))

fwrite(x = CF.intake.impact.sims.unconsumed.allgroups, 
       file = paste0(output_location, "CF.intake.envecosoc.sims.granular.unconsumed.csv"))

fwrite(x = impact.sims.inedible.allgroups, 
       file = paste0(output_location, "envecosoc.sims.granular.inedible.csv"))

fwrite(x = substitution.impact.sims.inedible.allgroups, 
       file = paste0(output_location, "substitution.envecosoc.sims.granular.inedible.csv"))

fwrite(x = combined.impact.sims.inedible.allgroups, 
       file = paste0(output_location, "combined.envecosoc.sims.granular.inedible.csv"))

fwrite(x = impact.sims.wasted.allgroups, 
       file = paste0(output_location, "envecosoc.sims.granular.wasted.csv"))

fwrite(x = substitution.impact.sims.wasted.allgroups, 
       file = paste0(output_location, "substitution.envecosoc.sims.granular.wasted.csv"))

fwrite(x = combined.impact.sims.wasted.allgroups, 
       file = paste0(output_location, "combined.envecosoc.sims.granular.wasted.csv"))

fwrite(x = current.intake.impact.sims.inedible.allgroups, 
       file = paste0(output_location, "current.intake.envecosoc.sims.granular.inedible.csv"))

fwrite(x = CF.intake.impact.sims.inedible.allgroups, 
       file = paste0(output_location, "CF.intake.envecosoc.sims.granular.inedible.csv"))

fwrite(x = current.intake.impact.sims.wasted.allgroups, 
       file = paste0(output_location, "current.intake.envecosoc.sims.granular.wasted.csv"))

fwrite(x = CF.intake.impact.sims.wasted.allgroups, 
       file = paste0(output_location, "CF.intake.envecosoc.sims.granular.wasted.csv"))

# All of what we just did pertains to population level impact. Next. we get 
# per capita sims for all outputs of interest. 

# Recall there are 31 outcome types of interest, each corresponding to some 
# combination of (1) what we are looking at (impact (as in, impact of change), 
# substitution impact, combined impact, population, delta (aka, change in intake), 
# current intake impact, counterfactual intake impact), and (2) level of 
# disaggregation (total, consumed, unconsumed, inedible, wasted).

# PER CAPITA 

impact.sims.percapita.total.allgroups <- 
  get.percapita.sims(impact.sims = impact.sims.total.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

substitution.impact.sims.percapita.total.allgroups <- 
  get.percapita.sims(impact.sims = substitution.impact.sims.total.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

combined.impact.percapita.total.allgroups <- 
  get.percapita.sims(impact.sims = combined.impact.sims.total.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

delta.percapita.allgroups <- 
  get.percapita.sims(impact.sims = delta.sims.total.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

impact.sims.percapita.consumed.allgroups <- 
  get.percapita.sims(impact.sims = impact.sims.consumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

substitution.impact.sims.percapita.consumed.allgroups <- 
  get.percapita.sims(impact.sims = substitution.impact.sims.consumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

combined.impact.percapita.consumed.allgroups <- 
  get.percapita.sims(impact.sims = combined.impact.sims.consumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

impact.sims.percapita.unconsumed.allgroups <- 
  get.percapita.sims(impact.sims = impact.sims.unconsumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

substitution.impact.sims.percapita.unconsumed.allgroups <- 
  get.percapita.sims(impact.sims = substitution.impact.sims.unconsumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

combined.impact.percapita.unconsumed.allgroups <- 
  get.percapita.sims(impact.sims = combined.impact.sims.unconsumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

current.intake.impact.sims.percapita.total.allgroups <- 
  get.percapita.sims(impact.sims = current.intake.impact.sims.total.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

CF.intake.impact.sims.percapita.total.allgroups <- 
  get.percapita.sims(impact.sims = CF.intake.impact.sims.total.allgroups, 
                     population.sims = population.sims.allgroups,  
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

current.intake.impact.sims.percapita.consumed.allgroups <- 
  get.percapita.sims(impact.sims = current.intake.impact.sims.consumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

CF.intake.impact.sims.percapita.consumed.allgroups <- 
  get.percapita.sims(impact.sims = CF.intake.impact.sims.consumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

current.intake.impact.sims.percapita.unconsumed.allgroups <- 
  get.percapita.sims(impact.sims = current.intake.impact.sims.unconsumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

CF.intake.impact.sims.percapita.unconsumed.allgroups <- 
  get.percapita.sims(impact.sims = CF.intake.impact.sims.unconsumed.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

impact.sims.percapita.inedible.allgroups <- 
  get.percapita.sims(impact.sims = impact.sims.inedible.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

substitution.impact.sims.percapita.inedible.allgroups <- 
  get.percapita.sims(impact.sims = substitution.impact.sims.inedible.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

combined.impact.percapita.inedible.allgroups <- 
  get.percapita.sims(impact.sims = combined.impact.sims.inedible.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

impact.sims.percapita.wasted.allgroups <- 
  get.percapita.sims(impact.sims = impact.sims.wasted.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

substitution.impact.sims.percapita.wasted.allgroups <- 
  get.percapita.sims(impact.sims = substitution.impact.sims.wasted.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

combined.impact.percapita.wasted.allgroups <- 
  get.percapita.sims(impact.sims = combined.impact.sims.wasted.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

current.intake.impact.sims.percapita.inedible.allgroups <- 
  get.percapita.sims(impact.sims = current.intake.impact.sims.inedible.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

CF.intake.impact.sims.percapita.inedible.allgroups <- 
  get.percapita.sims(impact.sims = CF.intake.impact.sims.inedible.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

current.intake.impact.sims.percapita.wasted.allgroups <- 
  get.percapita.sims(impact.sims = current.intake.impact.sims.wasted.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

CF.intake.impact.sims.percapita.wasted.allgroups <- 
  get.percapita.sims(impact.sims = CF.intake.impact.sims.wasted.allgroups, 
                     population.sims = population.sims.allgroups, 
                     pop.sims.names = paste("X", c(1:n.sims), sep = ""), 
                     new.pop.sims.names = paste("pop", c(1:n.sims), sep = ""), 
                     sims.names = paste("X", c(1:n.sims), sep = ""))

# Save them into csv files.

output_capita_location <- paste0(output_location, "per_capita/")

# first create directory if it doesn't already exist
ifelse(!dir.exists(file.path(output_capita_location)),
       dir.create(file.path(output_capita_location)),
       "Directory Exists")

fwrite(x = impact.sims.percapita.total.allgroups, 
       file = paste0(output_capita_location, 
                   "envecosoc.sims.percapita.granular.csv"))

fwrite(x = substitution.impact.sims.percapita.total.allgroups, 
       file = paste0(output_capita_location, 
                   "substitution.envecosoc.sims.percapita.granular.csv"))

fwrite(x = combined.impact.percapita.total.allgroups, 
       file = paste0(output_capita_location, 
                   "combined.envecosoc.sims.percapita.granular.csv"))

fwrite(x = delta.percapita.allgroups, 
       file = paste0(output_capita_location, 
                   "delta.envecosoc.sims.percapita.granular.csv"))

fwrite(x = impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "envecosoc.sims.percapita.consumed.granular.csv"))

fwrite(x = substitution.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "substitution.envecosoc.sims.percapita.consumed.granular.csv"))

fwrite(x = combined.impact.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "combined.envecosoc.sims.percapita.consumed.granular.csv"))

fwrite(x = impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "envecosoc.sims.percapita.unconsumed.granular.csv"))

fwrite(x = substitution.impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "substitution.envecosoc.sims.percapita.unconsumed.granular.csv"))

fwrite(x = combined.impact.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "combined.envecosoc.sims.percapita.unconsumed.granular.csv"))

fwrite(x = current.intake.impact.sims.percapita.total.allgroups, 
       file = paste0(output_capita_location, 
                   "current.intake.envecosoc.sims.percapita.granular.csv"))

fwrite(x = CF.intake.impact.sims.percapita.total.allgroups, 
       file = paste0(output_capita_location, 
                   "CF.envecosoc.intake.sims.percapita.granular.csv"))

fwrite(x = current.intake.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "current.intake.envecosoc.sims.percapita.granular.consumed.csv"))

fwrite(x = CF.intake.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "CF.intake.envecosoc.sims.percapita.granular.consumed.csv"))

fwrite(x = current.intake.impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "current.intake.envecosoc.sims.percapita.granular.unconsumed.csv"))

fwrite(x = CF.intake.impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "CF.intake.envecosoc.sims.percapita.granular.unconsumed.csv"))

fwrite(x = impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "envecosoc.sims.percapita.inedible.granular.csv"))

fwrite(x = substitution.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "substitution.envecosoc.sims.percapita.inedible.granular.csv"))

fwrite(x = combined.impact.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "combined.envecosoc.sims.percapita.inedible.granular.csv"))

fwrite(x = impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "envecosoc.sims.percapita.wasted.granular.csv"))

fwrite(x = substitution.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "substitution.envecosoc.sims.percapita.wasted.granular.csv"))

fwrite(x = combined.impact.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "combined.envecosoc.sims.percapita.wasted.granular.csv"))

fwrite(x = current.intake.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "current.intake.envecosoc.sims.percapita.granular.inedible.csv"))

fwrite(x = CF.intake.impact.sims.percapita.consumed.allgroups, 
       file = paste0(output_capita_location, 
                   "CF.intake.envecosoc.sims.percapita.granular.inedible.csv"))

fwrite(x = current.intake.impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "current.intake.envecosoc.sims.percapita.granular.wasted.csv"))

fwrite(x = CF.intake.impact.sims.percapita.unconsumed.allgroups, 
       file = paste0(output_capita_location, 
                   "CF.intake.envecosoc.sims.percapita.granular.wasted.csv"))

# Get summary stats from simulations (mean, median, standard deviation, 
# uncertainty intervals, i.e., the things we report in the manuscript) for 
# all outcome of interests and save (for both total and per capita impact).

# get summary stats 

all.impact.sims.summary <- 
  get.summary.stats(
    impact.sims = impact.sims.total.allgroups, 
    substitution.impact.sims = substitution.impact.sims.total.allgroups, 
    combined.impact.sims = combined.impact.sims.total.allgroups, 
    delta.sims = delta.sims.total.allgroups, 
    impact.sims.consumed = impact.sims.consumed.allgroups, 
    substitution.impact.sims.consumed = substitution.impact.sims.consumed.allgroups, 
    combined.impact.sims.consumed = combined.impact.sims.consumed.allgroups,
    impact.sims.unconsumed = impact.sims.unconsumed.allgroups, 
    substitution.impact.sims.unconsumed = substitution.impact.sims.unconsumed.allgroups, 
    combined.impact.sims.unconsumed = combined.impact.sims.unconsumed.allgroups, 
    current.intake.impact.sims = current.intake.impact.sims.total.allgroups, 
    current.intake.impact.sims.consumed = current.intake.impact.sims.consumed.allgroups, 
    current.intake.impact.sims.unconsumed = current.intake.impact.sims.unconsumed.allgroups, 
    CF.intake.impact.sims = CF.intake.impact.sims.total.allgroups, 
    CF.intake.impact.sims.consumed = CF.intake.impact.sims.consumed.allgroups, 
    CF.intake.impact.sims.unconsumed = CF.intake.impact.sims.unconsumed.allgroups,
    impact.sims.inedible = impact.sims.inedible.allgroups, 
    substitution.impact.sims.inedible = substitution.impact.sims.inedible.allgroups, 
    combined.impact.sims.inedible = combined.impact.sims.inedible.allgroups,
    impact.sims.wasted = impact.sims.wasted.allgroups, 
    substitution.impact.sims.wasted = substitution.impact.sims.wasted.allgroups, 
    combined.impact.sims.wasted = combined.impact.sims.wasted.allgroups, 
    current.intake.impact.sims.inedible = current.intake.impact.sims.inedible.allgroups, 
    current.intake.impact.sims.wasted = current.intake.impact.sims.wasted.allgroups, 
    CF.intake.impact.sims.inedible = CF.intake.impact.sims.inedible.allgroups, 
    CF.intake.impact.sims.wasted = CF.intake.impact.sims.wasted.allgroups,
    vars = paste("X", 1:n.sims, sep = ""))

write_csv(x = all.impact.sims.summary, 
          file = paste0(output_location, "envecosoc.summary.granular.csv"))

# now use the function for per capita results
all.impact.sims.percapita.summary <- 
  get.summary.stats(
    impact.sims = impact.sims.percapita.total.allgroups, 
    substitution.impact.sims = substitution.impact.sims.percapita.total.allgroups, 
    combined.impact.sims = combined.impact.percapita.total.allgroups,
    delta.sims = delta.percapita.allgroups, 
    impact.sims.consumed = impact.sims.percapita.consumed.allgroups, 
    substitution.impact.sims.consumed = substitution.impact.sims.percapita.unconsumed.allgroups, 
    combined.impact.sims.consumed = combined.impact.percapita.unconsumed.allgroups,
    impact.sims.unconsumed = impact.sims.percapita.unconsumed.allgroups, 
    substitution.impact.sims.unconsumed = substitution.impact.sims.percapita.unconsumed.allgroups, 
    combined.impact.sims.unconsumed = combined.impact.percapita.unconsumed.allgroups, 
    current.intake.impact.sims = current.intake.impact.sims.percapita.total.allgroups, 
    current.intake.impact.sims.consumed = current.intake.impact.sims.percapita.consumed.allgroups, 
    current.intake.impact.sims.unconsumed = current.intake.impact.sims.percapita.unconsumed.allgroups, 
    CF.intake.impact.sims = CF.intake.impact.sims.percapita.total.allgroups,
    CF.intake.impact.sims.consumed = CF.intake.impact.sims.percapita.consumed.allgroups, 
    CF.intake.impact.sims.unconsumed = CF.intake.impact.sims.percapita.unconsumed.allgroups,
    impact.sims.inedible = impact.sims.percapita.inedible.allgroups, 
    substitution.impact.sims.inedible = substitution.impact.sims.percapita.inedible.allgroups, 
    combined.impact.sims.inedible = combined.impact.percapita.inedible.allgroups,
    impact.sims.wasted = impact.sims.percapita.wasted.allgroups, 
    substitution.impact.sims.wasted = substitution.impact.sims.percapita.wasted.allgroups, 
    combined.impact.sims.wasted = combined.impact.percapita.wasted.allgroups,
    current.intake.impact.sims.inedible = current.intake.impact.sims.percapita.inedible.allgroups, 
    current.intake.impact.sims.wasted = current.intake.impact.sims.percapita.wasted.allgroups, 
    CF.intake.impact.sims.inedible = CF.intake.impact.sims.percapita.inedible.allgroups, 
    CF.intake.impact.sims.wasted = CF.intake.impact.sims.percapita.wasted.allgroups,
    vars = paste("X", 1:n.sims, sep = ""))

write_csv(x = all.impact.sims.percapita.summary, 
          file = paste0(output_capita_location, 
                      "envecosoc.percapita.summary.granular.csv"))

# Next, we want to get summary stat for all possible strata combinations (not 
# just the most granular, age/sex/race/foodgroup).

# First we create a list of all possible strata combinations. Note, the string 
# vector "strata" contains all the strata being used for the analysis. If you want 
# to use different strata, you need to change the line defining "strata" accordingly. 

# get all possible combinations of strata
strata <- c("age_gp", "sex_gp", "race_gp", "Foodgroup")

strata.combos <- list()

for(i in 1:length(strata)) {
  
  x <- combn(strata, i)
  combos.for.i <- split(x, rep(1:ncol(x), each = nrow(x)))
  strata.combos <- c(strata.combos, combos.for.i)
  
}

# add null element to list to use for overall numbers
strata.combos["null"] <- list(NULL) 

# Get summary stats by all strata combinations for all outcomes of interest. Note 
# that the objects being created here are lists of lists. The second level are the 
# different strata combinations. The first level denotes outcome type: 1 = total 
# summary stats, 2 = total sims, 3 = per capita summary stats, 4 = per capita sims. 
# That means, for example, that summ.stats.by.strata.cost is a list, and each 
# element within that list is also a list. Specifically, it is a list of 4 lists 
# (each one denoting outcome type), and each of those 4 lists is itself a list of 
# 48 data frames, one for each subgroup. 

summ.stats.by.strata.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = impact.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.substitution.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = substitution.impact.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.combined.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = combined.impact.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.delta_forenvecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = delta.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.combined.consumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = combined.impact.sims.consumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.combined.unconsumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = combined.impact.sims.unconsumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.current.intake.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = current.intake.impact.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.CF.intake.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = CF.intake.impact.sims.total.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.current.intake.consumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = current.intake.impact.sims.consumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.CF.intake.consumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = CF.intake.impact.sims.consumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.current.intake.unconsumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = current.intake.impact.sims.unconsumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.CF.intake.unconsumed.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = CF.intake.impact.sims.unconsumed.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.combined.inedible.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = combined.impact.sims.inedible.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.combined.wasted.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = combined.impact.sims.wasted.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.current.intake.inedible.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = current.intake.impact.sims.inedible.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.CF.intake.inedible.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = CF.intake.impact.sims.inedible.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.current.intake.wasted.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = current.intake.impact.sims.wasted.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

summ.stats.by.strata.CF.intake.wasted.envecosoc <- 
  summary.stats.by.strata(impact.sims.allgroups = CF.intake.impact.sims.wasted.allgroups, 
                          strata.combos.list = strata.combos, 
                          pop.sims = pop.sims)

# Next, we're going to go through a for-loop and output all outcomes into their 
# own csv files for each strata combination (In my opinion, this is a better way 
# to go about it than manually saving each file. Less prone to error, and robust 
# to future changes). 

output_subgroup_location <- paste0(output_location, "By_SubGroup/")

# create 'By_SubGroup' output directory
ifelse(!dir.exists(file.path(output_subgroup_location)),
       dir.create(file.path(output_subgroup_location)),
       "Directory Exists")

# create 'By_SubGroup/full_sims' output directory
ifelse(!dir.exists(file.path(paste0(output_subgroup_location, "full_sims"))),
       dir.create(file.path(paste0(output_subgroup_location, "full_sims"))),
       "Directory Exists")

# create 'per_capita/By_SubGroup' output directory
ifelse(!dir.exists(file.path(paste0(output_capita_location, "By_SubGroup"))),
       dir.create(file.path(paste0(output_capita_location, "By_SubGroup"))),
       "Directory Exists")

# create 'per_capita/By_SubGroup/full_sims' output directory
ifelse(!dir.exists(file.path(paste0(output_capita_location, "By_SubGroup/full_sims"))),
       dir.create(file.path(paste0(output_capita_location, "By_SubGroup/full_sims"))),
       "Directory Exists")

# First, start the for-loop.

for(j in 1:length(strata.combos)) {
  
  # First, save full sim outputs for total and per capita
  
  fwrite(x = summ.stats.by.strata.envecosoc[[2]][[j]], 
         file = paste0(output_subgroup_location, 
                     "full_sims/full.envecosoc.sims.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.substitution.envecosoc[[2]][[j]], 
         file = paste0(output_subgroup_location, 
                     "full_sims/full.envecosoc.sims.substitution.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.combined.envecosoc[[2]][[j]], 
         file = paste0(output_subgroup_location, 
                     "full_sims/full.envecosoc.sims.combined.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.envecosoc[[4]][[j]], 
         file = paste0(output_capita_location, 
                     "By_SubGroup/full_sims/full.envecosoc.sims.percapita.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.substitution.envecosoc[[4]][[j]], 
         file = paste0(output_capita_location, 
                     "By_SubGroup/full_sims/full.envecosoc.sims.percapita.substitution.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  fwrite(x = summ.stats.by.strata.combined.envecosoc[[4]][[j]], 
         file = paste0(output_capita_location, 
                     "By_SubGroup/full_sims/full.envecosoc.sims.percapita.combined.output_by_", 
                    paste(strata.combos[[j]], collapse = '_'), ".csv"))
  
  # Next, create outcome names to be used as variable names, and put relevant 
  # outcome data frames into two lists (for total and per capita). Within the 
  # k-loop, the summary stat variable names for the dataframes listed in 
  # outcome.data" are being changed so that they all consistently use the 
  # naming convention outlined in "outcome.names". 
  
  to.rename <- c("lower_bound (2.5th percentile)", 
                 "median", 
                 "upper_bound  (97.5th percentile)", 
                 "mean", 
                 "SD")
  
  outcome.names <- c("impact", 
                     "substitution_impact", 
                     "combined_impact", 
                     "delta", 
                     "combined_consumed_impact", 
                     "combined_unconsumed_impact", 
                     "current_intake_impact", 
                     "CF_intake_impact", 
                     "current_intake_consumed_impact", 
                     "CF_intake_consumed_impact", 
                     "current_intake_unconsumed_impact", 
                     "CF_intake_unconsumed_impact",
                     
                     # inedible, wasted
                     "combined_inedible_impact", 
                     "combined_wasted_impact", 
                     "current_intake_inedible_impact", 
                     "CF_intake_inedible_impact", 
                     "current_intake_wasted_impact", 
                     "CF_intake_wasted_impact")
  
  outcome.data <- list(summ.stats.by.strata.envecosoc[[1]][[j]],
                       summ.stats.by.strata.substitution.envecosoc[[1]][[j]],
                       summ.stats.by.strata.combined.envecosoc[[1]][[j]],
                       summ.stats.by.strata.delta_forenvecosoc[[1]][[j]],
                       summ.stats.by.strata.combined.consumed.envecosoc[[1]][[j]],
                       summ.stats.by.strata.combined.unconsumed.envecosoc[[1]][[j]],
                       summ.stats.by.strata.current.intake.envecosoc[[1]][[j]],
                       summ.stats.by.strata.CF.intake.envecosoc[[1]][[j]],
                       summ.stats.by.strata.current.intake.consumed.envecosoc[[1]][[j]],
                       summ.stats.by.strata.CF.intake.consumed.envecosoc[[1]][[j]],
                       summ.stats.by.strata.current.intake.unconsumed.envecosoc[[1]][[j]],
                       summ.stats.by.strata.CF.intake.unconsumed.envecosoc[[1]][[j]],
                       
                       # inedible, wasted
                       summ.stats.by.strata.combined.inedible.envecosoc[[1]][[j]],
                       summ.stats.by.strata.combined.wasted.envecosoc[[1]][[j]],
                       summ.stats.by.strata.current.intake.inedible.envecosoc[[1]][[j]],
                       summ.stats.by.strata.CF.intake.inedible.envecosoc[[1]][[j]],
                       summ.stats.by.strata.current.intake.wasted.envecosoc[[1]][[j]],
                       summ.stats.by.strata.CF.intake.wasted.envecosoc[[1]][[j]]
                       )
  
  outcome.data.per.capita <- list(
    summ.stats.by.strata.envecosoc[[3]][[j]],
    summ.stats.by.strata.substitution.envecosoc[[3]][[j]],
    summ.stats.by.strata.combined.envecosoc[[3]][[j]],
    summ.stats.by.strata.delta_forenvecosoc[[3]][[j]],
    summ.stats.by.strata.combined.consumed.envecosoc[[3]][[j]],
    summ.stats.by.strata.combined.unconsumed.envecosoc[[3]][[j]],
    summ.stats.by.strata.current.intake.envecosoc[[3]][[j]],
    summ.stats.by.strata.CF.intake.envecosoc[[3]][[j]],
    summ.stats.by.strata.current.intake.consumed.envecosoc[[3]][[j]],
    summ.stats.by.strata.CF.intake.consumed.envecosoc[[3]][[j]],
    summ.stats.by.strata.current.intake.unconsumed.envecosoc[[3]][[j]],
    summ.stats.by.strata.CF.intake.unconsumed.envecosoc[[3]][[j]],
                                  
    # inedible, wasted
    summ.stats.by.strata.combined.inedible.envecosoc[[3]][[j]],
    summ.stats.by.strata.combined.wasted.envecosoc[[3]][[j]],
    summ.stats.by.strata.current.intake.inedible.envecosoc[[3]][[j]],
    summ.stats.by.strata.CF.intake.inedible.envecosoc[[3]][[j]],
    summ.stats.by.strata.current.intake.wasted.envecosoc[[3]][[j]],
    summ.stats.by.strata.CF.intake.wasted.envecosoc[[3]][[j]]
    )
  
  for(k in 1:length(outcome.names)){
    
    names(
      outcome.data[[k]])[which(names(
        outcome.data[[k]]) %in% c("lower_bound (2.5th percentile)", 
                                  "median", 
                                  "upper_bound  (97.5th percentile)", 
                                  "mean", 
                                  "SD"))] <- 
      paste(outcome.names[k], c("lower_bound (2.5th percentile)", 
                                "median", 
                                "upper_bound  (97.5th percentile)", 
                                "mean", 
                                "SD"), sep = "_")
    
    names(
      outcome.data.per.capita[[k]])[which(names(
        outcome.data.per.capita[[k]]) %in% c("lower_bound (2.5th percentile)", 
                                             "median", 
                                             "upper_bound  (97.5th percentile)", 
                                             "mean", 
                                             "SD"))] <- 
      paste(outcome.names[k], c("lower_bound (2.5th percentile)", 
                                "median", 
                                "upper_bound  (97.5th percentile)", 
                                "mean", 
                                "SD"), sep = "_")
    
    }
  
  # Then, for both lists, merge their contents into a data frame.
  
  list.names.intersect <- Reduce(intersect, lapply(outcome.data, colnames))
  
  all.summ.stats.by.strata <- 
    Reduce(function(...) merge(..., by = list.names.intersect), outcome.data)
  
  all.summ.stats.by.strata.per.capita <- 
    Reduce(function(...) merge(..., by = list.names.intersect), outcome.data.per.capita)
  
  # Save the two outputs created (one for total and one for per capita) as csvs. 
  # Each of these data frames will have summary stats for all outputs of 
  # interest for subgroup j.
  
  write_csv(x = all.summ.stats.by.strata, 
            file = paste0(output_subgroup_location, 
                        "summary.output_by_", 
                        paste(strata.combos[[j]], collapse = '_'), 
                        ".envecosoc.csv"))
  
  write_csv(x = all.summ.stats.by.strata.per.capita, 
            file = paste0(output_capita_location, 
                        "By_SubGroup/summary.output.percapita_by_", 
                        paste(strata.combos[[j]], collapse = '_'), 
                        ".envecosoc.csv"))
  
} # end of loop

