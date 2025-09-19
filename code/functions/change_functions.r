# Title: Change functions
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-15-2025

############################################################################

### Documentation for R script ###

# change_functions.r contains various functions used in calculating and 
# summarizing the impacts for the non-health (env-eco-soc) outcomes. 
# All of the functions in this script are used in calculate_change.r.

############################################################################

# The first function here, simulate.impact, is the core of the entire operation. 
# It runs monte carlo simulations for env-eco-soc impact calculations and 
# outputs the resulting distributions. The inputs are going to be specific to the 
# population and dietary factor of interest. In "calculate_change.r", 
# this function is being called for each pair of subgroup/dietary factor. 

# The meaning of the 21 inputs (stated below) should be mostly self explanatory, 
# but 2 things may need explanation.
# 
# (1) "substitution.impact.factor.means", "substitution.impact.factor.ses", 
# "substitution_unit" refer to substitution impact. While we don't use this 
# right now (all these inputs are zero), we can calculate substitution effects 
# as well. Say, if we change the diet for just one diet factor, we need replace 
# the calories from other foods. So, we need to calculate the impact of the 
# increase/decrease that comes from the substitution (probably using the impact 
# per unit of the average diet as substitution impact factor). Since we keep energy 
# constant for the five alternate scenarios, we don't actually use this for the 
# 4 pillars paper. But if you're not keeping energy constant when going from 
# current to counterfactual intake (i.e., just changing intake for one dietary 
# factor), you'd want to use substitution effects (and by "not use", I mean 
# setting all substitution impact factors to zero, not commenting out code 
# related to substitution factors).
# 
# (2) "current_foodwaste_p", "current_foodwaste_p_se, "counterfactual_foodwaste_p", 
# "counterfactual_foodwaste_p_se": Don't forget that food waste proportion is 
# conditional on edible. It's the proportion of the food wasted out of the 
# edible portion. It is NOT food waste / total food including inedible portion.

simulate.impact <- function(current.mean, current.se, 
                          counterfactual.mean, counterfactual.se = 0,
                          impact.factor.names, impact.factor.means,
                          impact.factor.ses, RRunit, 
                          substitution.impact.factor.means,
                          substitution.impact.factor.ses, substitution_unit, 
                          current_inedible_p, current_inedible_p_se,
                          current_foodwaste_p, current_foodwaste_p_se,
                          counterfactual_inedible_p, counterfactual_inedible_p_se,
                          counterfactual_foodwaste_p, counterfactual_foodwaste_p_se,
                          n.sims,
                          # population.sims must be a vector, so make sure to 
                          # convert to vector (from single row datafame) if necessary
                          population.sims) {
  
  # First simulate intake means for current and counterfactual scenarios. We 
  # are creating a matrix of sims because we are going to use them to do 
  # calculations for all outcomes from three pillars (with each row corresponding 
  # to an outcome of interest). 
  
  current.mean.sims <- rnorm(n = n.sims, 
                           mean = current.mean, 
                           sd = current.se)
  
  current.mean.sims.matrix <- t(matrix(rep(current.mean.sims, 
                                           times = length(impact.factor.names)),
                                       ncol = length(impact.factor.names)))
  
  counterfactual.mean.sims <- rnorm(n = n.sims, 
                                  mean = counterfactual.mean, 
                                  sd = counterfactual.se)
  
  counterfactual.mean.sims.matrix <- t(matrix(rep(counterfactual.mean.sims, 
                                                  times = length(impact.factor.names)),
                                              ncol = length(impact.factor.names)))
  
  # Calculate correction factors (from inedible and waste proportion) to go 
  # from current intake (our input) to edible amount and total produced 
  # (edible + inedible), which is what we need to use to calculate impact 
  # on these three pillars. 
  
  current.correction <- 1 / ((1 - rnorm(length(impact.factor.names), 
                                    mean = current_inedible_p, 
                                    sd = current_inedible_p_se)) * 
                             (1 - rnorm(length(impact.factor.names), 
                                      mean = current_foodwaste_p, 
                                      sd = current_foodwaste_p_se)))
  
  counterfactual.correction <- 1 / ((1 - rnorm(length(impact.factor.names), 
                                           mean = counterfactual_inedible_p, 
                                           sd = counterfactual_inedible_p_se)) * 
                                    (1 - rnorm(length(impact.factor.names), 
                                             mean = counterfactual_foodwaste_p, 
                                             sd = counterfactual_foodwaste_p_se)))
  
  # this is a correction factor to get total edible. 
  # "current.correction" above is for getting total produced 
  # (edible (consumed and wasted) + inedible)
  
  current.correction.for.waste <- 1 / (1 - rnorm(length(impact.factor.names), 
                                             mean = current_foodwaste_p, 
                                             sd = current_foodwaste_p_se))
  
  counterfactual.correction.for.waste <- 1 / (1 - rnorm(length(impact.factor.names), 
                                                    mean = counterfactual_foodwaste_p, 
                                                    sd = counterfactual_foodwaste_p_se))
  
  # Apply corrections, and based on that, calculate unconsumed, edible, 
  # inedible and wasted amounts. 
  
  current.mean.sims.matrix.after.correction <- 
    current.mean.sims.matrix * current.correction
  
  counterfactual.mean.sims.matrix.after.correction <- 
    counterfactual.mean.sims.matrix * counterfactual.correction
  
  current.mean.sims.matrix.unconsumed <- 
    current.mean.sims.matrix.after.correction - current.mean.sims.matrix
  
  counterfactual.mean.sims.matrix.unconsumed <- 
    counterfactual.mean.sims.matrix.after.correction - 
    counterfactual.mean.sims.matrix
  
  current.mean.sims.matrix.edible <- 
    current.mean.sims.matrix * current.correction.for.waste
  
  counterfactual.mean.sims.matrix.edible <- 
    counterfactual.mean.sims.matrix * counterfactual.correction.for.waste
  
  current.mean.sims.matrix.inedible <- 
    current.mean.sims.matrix.after.correction - current.mean.sims.matrix.edible
  
  counterfactual.mean.sims.matrix.inedible <- 
    counterfactual.mean.sims.matrix.after.correction - 
    counterfactual.mean.sims.matrix.edible
  
  current.mean.sims.matrix.wasted <- 
    current.mean.sims.matrix.edible - current.mean.sims.matrix
  
  counterfactual.mean.sims.matrix.wasted <- 
    counterfactual.mean.sims.matrix.edible - counterfactual.mean.sims.matrix
  
  # Calculate deltas (difference between current and counterfactual)
  
  delta.sims.matrix <- 
    counterfactual.mean.sims.matrix.after.correction - 
    current.mean.sims.matrix.after.correction
  
  delta.sims.matrix.consumed <- 
    counterfactual.mean.sims.matrix - 
    current.mean.sims.matrix
  
  delta.sims.matrix.unconsumed <- 
    counterfactual.mean.sims.matrix.unconsumed - 
    current.mean.sims.matrix.unconsumed
  
  delta.sims.matrix.total <- 
    delta.sims.matrix.consumed + 
    delta.sims.matrix.unconsumed
  
  delta.sims.matrix.inedible <- counterfactual.mean.sims.matrix.inedible - 
    current.mean.sims.matrix.inedible
  
  delta.sims.matrix.wasted <- counterfactual.mean.sims.matrix.wasted - 
    current.mean.sims.matrix.wasted
  
  # Simulate impact factors.
  
  impact.factor.sims <- 
    matrix(rnorm(n = n.sims * length(impact.factor.names), 
                 mean = as.numeric(impact.factor.means), 
                 sd = as.numeric(impact.factor.ses)), 
           nrow = length(impact.factor.names))
  
  substitution.impact.factor.sims <- 
    matrix(rnorm(n = n.sims * length(impact.factor.names), 
                 mean = as.numeric(substitution.impact.factor.means), 
                 sd = as.numeric(substitution.impact.factor.ses)), 
           nrow = length(impact.factor.names))
  
  # Calculate impact and assign row names to the output matrices (for 
  # consumed, unconsumed, total, inedible, and wasted)
  
  impact.sims.consumed <- 
    impact.factor.sims * 
    delta.sims.matrix.consumed * 
    population.sims / 
    as.numeric(RRunit)
  
  impact.sims.unconsumed <- 
    impact.factor.sims * 
    delta.sims.matrix.unconsumed * 
    population.sims / 
    as.numeric(RRunit)
  
  impact.sims.total <- 
    impact.sims.consumed + 
    impact.sims.unconsumed
  
  impact.sims.inedible <- 
    impact.factor.sims * 
    delta.sims.matrix.inedible * 
    population.sims / 
    as.numeric(RRunit)
  
  impact.sims.wasted <- 
    impact.factor.sims * 
    delta.sims.matrix.wasted * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(impact.sims.consumed) <- impact.factor.names
  
  rownames(impact.sims.unconsumed) <- impact.factor.names
  
  rownames(impact.sims.total) <- impact.factor.names
  
  rownames(impact.sims.inedible) <- impact.factor.names
  
  rownames(impact.sims.wasted) <- impact.factor.names
  
  # Calculate substitution impact and assign row names to the output matrices
  
  substitution.sims.consumed <- 
    substitution.impact.factor.sims * 
    -delta.sims.matrix.consumed * 
    population.sims / 
    substitution_unit
  
  substitution.sims.unconsumed <- 
    substitution.impact.factor.sims * 
    -delta.sims.matrix.unconsumed * 
    population.sims / 
    substitution_unit
  
  substitution.sims.total <- 
    substitution.sims.consumed + 
    substitution.sims.unconsumed
  
  substitution.sims.inedible <- 
    substitution.impact.factor.sims * 
    -delta.sims.matrix.inedible * 
    population.sims / 
    substitution_unit
  
  substitution.sims.wasted <- 
    substitution.impact.factor.sims * 
    -delta.sims.matrix.wasted * 
    population.sims / 
    substitution_unit
  
  rownames(substitution.sims.consumed) <- impact.factor.names
  
  rownames(substitution.sims.unconsumed) <- impact.factor.names
  
  rownames(substitution.sims.total) <- impact.factor.names
  
  rownames(substitution.sims.inedible) <- impact.factor.names
  
  rownames(substitution.sims.wasted) <- impact.factor.names
  
  # Combine impact and substitution impact to get total impact
  
  combined.impact.sims.consumed <- 
    impact.sims.consumed + substitution.sims.consumed
  
  combined.impact.sims.unconsumed <- 
    impact.sims.unconsumed + substitution.sims.unconsumed
  
  combined.impact.sims.total <- 
    combined.impact.sims.consumed + combined.impact.sims.unconsumed
  
  combined.impact.sims.inedible <- 
    impact.sims.consumed + substitution.sims.inedible
  
  combined.impact.sims.wasted <- 
    impact.sims.unconsumed + substitution.sims.wasted
  
  # Calculate impact of current and counterfactual intake on the outcomes 
  # (as opposed to the the impact of the shift from current to counterfactual, 
  # which is what we've been doing up to this point).
  
  combined.impact.sims.consumed <- 
    impact.sims.consumed + substitution.sims.consumed
  
  combined.impact.sims.unconsumed <- 
    impact.sims.unconsumed + substitution.sims.unconsumed
  
  combined.impact.sims.total <- 
    combined.impact.sims.consumed + combined.impact.sims.unconsumed
  
  combined.impact.sims.inedible <- 
    impact.sims.consumed + substitution.sims.inedible
  
  combined.impact.sims.wasted <- 
    impact.sims.unconsumed + substitution.sims.wasted
  
  # calculate total impact, not just difference between current and counterfactual
  
  current.impact.sims <- 
    impact.factor.sims * 
    current.mean.sims.matrix.after.correction * 
    population.sims / 
    as.numeric(RRunit)
  
  counterfactual.impact.sims <- 
    impact.factor.sims * 
    counterfactual.mean.sims.matrix.after.correction * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(current.impact.sims) <- impact.factor.names
  
  rownames(counterfactual.impact.sims) <- impact.factor.names
  
  # consumed
  current.impact.sims.consumed <- 
    impact.factor.sims * 
    current.mean.sims.matrix * 
    population.sims / 
    as.numeric(RRunit)
  
  counterfactual.impact.sims.consumed <- 
    impact.factor.sims * 
    counterfactual.mean.sims.matrix * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(current.impact.sims.consumed) <- impact.factor.names
  
  rownames(counterfactual.impact.sims.consumed) <- impact.factor.names
  
  # unconsumed
  current.impact.sims.unconsumed <- 
    impact.factor.sims * 
    current.mean.sims.matrix.unconsumed * 
    population.sims / 
    as.numeric(RRunit)
  
  counterfactual.impact.sims.unconsumed <- 
    impact.factor.sims * 
    counterfactual.mean.sims.matrix.unconsumed * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(current.impact.sims.unconsumed) <- impact.factor.names
  
  rownames(counterfactual.impact.sims.unconsumed) <- impact.factor.names
  
  # inedible
  current.impact.sims.inedible <- 
    impact.factor.sims * 
    current.mean.sims.matrix.inedible * 
    population.sims / 
    as.numeric(RRunit)
  
  counterfactual.impact.sims.inedible <- 
    impact.factor.sims * 
    counterfactual.mean.sims.matrix.inedible * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(current.impact.sims.inedible) <- impact.factor.names
  
  rownames(counterfactual.impact.sims.inedible) <- impact.factor.names
  
  # wasted
  current.impact.sims.wasted <- 
    impact.factor.sims * 
    current.mean.sims.matrix.wasted * 
    population.sims / 
    as.numeric(RRunit)
  
  counterfactual.impact.sims.wasted <- 
    impact.factor.sims * 
    counterfactual.mean.sims.matrix.wasted * 
    population.sims / 
    as.numeric(RRunit)
  
  rownames(current.impact.sims.wasted) <- impact.factor.names
  
  rownames(counterfactual.impact.sims.wasted) <- impact.factor.names
  
  # Calculating total deltas at population levels by multiplying our mean 
  # deltas with population counts
  
  population.sims.output <- as.matrix(population.sims)
  
  rownames(population.sims.output) <- impact.factor.names 
  
  delta.sims.output.consumed <- 
    as.matrix(delta.sims.matrix.consumed * population.sims)
  
  rownames(delta.sims.output.consumed) <- impact.factor.names
  
  delta.sims.output.unconsumed <- 
    as.matrix(delta.sims.matrix.unconsumed * population.sims)
  
  rownames(delta.sims.output.unconsumed) <- impact.factor.names
  
  delta.sims.output.total <- 
    as.matrix(delta.sims.matrix.total * population.sims)
  
  rownames(delta.sims.output.total) <- impact.factor.names
  
  delta.sims.output.inedible <- 
    as.matrix(delta.sims.matrix.inedible * population.sims)
  
  rownames(delta.sims.output.inedible) <- impact.factor.names
  
  delta.sims.output.wasted <- 
    as.matrix(delta.sims.matrix.wasted * population.sims)
  
  rownames(delta.sims.output.wasted) <- impact.factor.names
  
  # Take all output created and put it into a list of 31 matrices, each 
  # corresponding to simulations for a specific output of interest for the 
  # particular dietary factor - counterfactual diet pair, and return that list. 
  # Let's spell out what all the outcomes are here:
  # 
  # 1. impact.sims.total = impact of change from current to counterfactual for 
  # dietary factor of interest
  # 
  # 2. substitution.sims.total = impact from substitution of calories needed 
  # to maintain energy going from current to counterfactual for dietary factor 
  # of interest
  # 
  # 3. combined.impact.sims.total = combined impact of dietary factor shift 
  # and corresponding substitution impact (impact.sim.total + substitution.sims.total)
  # 
  # 4. population.sims.total = population count
  # 
  # 5. delta.sims.total = difference in dietary factor amount produced for 
  # current vs counterfactual
  # 
  # 6. impact.sims.consumed = impact of change from current to counterfactual 
  # for dietary factor of interest, but just for the consumed part of the dietary 
  # factor
  # 
  # 7. substitution.sims.consumed = impact from substitution of calories needed 
  # to maintain energy going from current to counterfactual for dietary factor of 
  # interest, but just for the consumed part of the dietary factor
  # 
  # 8. combined.impact.sims.consumed = combined impact of dietary factor shift 
  # and corresponding substitution impact (impact.sim.consumed + substitution.sims.consumed), 
  # but just for the consumed part of the dietary factor
  # 
  # 9. delta.sims.consumed = difference in dietary factor amount produced for 
  # current vs counterfactual, but just for the consumed part of the dietary factor
  # 
  # 10. impact.sims.unconsumed = impact of change from current to counterfactual 
  # for dietary factor of interest, but just for the unconsumed part of the dietary factor
  # 
  # 11. substitution.sims.consumed = impact from substitution of calories needed 
  # to maintain energy going from current to counterfactual for dietary factor of 
  # interest, but just for the consumed part of the dietary factor
  # 
  # 12. combined.impact.sims.consumed = combined impact of dietary factor shift 
  # and corresponding substitution impact (impact.sim.unconsumed + substitution.sims.unconsumed), 
  # but just for the unconsumed part of the dietary factor
  # 
  # 13. delta.sims.unconsumed = difference in dietary factor amount produced 
  # for current vs counterfactual, but just for the unconsumed part of the dietary factor
  # 
  # 14. current.intake.impact.total = impact of current intake for dietary 
  # factor of interest
  # 
  # 15. CF.intake.impact.total = impact of counterfactual intake for dietary 
  # factor of interest
  # 
  # 16. current.intake.impact.consumed = impact of current intake for dietary 
  # factor of interest, but just for the consumed part
  # 
  # 17. CF.intake.impact.consumed = impact of counterfactual intake for dietary 
  # factor of interest, but just for the consumed part
  # 
  # 18. current.intake.impact.unconsumed = impact of current intake for dietary 
  # factor of interest, but just for the unconsumed part
  # 
  # 19. CF.intake.impact.unconsumed = impact of counterfactual intake for dietary 
  # factor of interest, but just for the unconsumed part
  # 
  # 20. impact.sims.inedible = impact of change from current to counterfactual for 
  # dietary factor of interest, but just for the inedible part of the dietary factor
  # 
  # 21. substitution.sims.inedible = impact from substitution of calories needed 
  # to maintain energy going from current to counterfactual for dietary factor of 
  # interest, but just for the inedible part of the dietary factor
  # 
  # 22. combined.impact.sims.inedible = combined impact of dietary factor shift 
  # and corresponding substitution impact (impact.sim.inedible + substitution.sims.inedible), 
  # but just for the inedible part of the dietary factor
  # 
  # 23. delta.sims.inedible = difference in dietary factor amount produced for 
  # current vs counterfactual, but just for the inedible part of the dietary factor
  # 
  # 24. impact.sims.wasted = impact of change from current to counterfactual for 
  # dietary factor of interest, but just for the wasted part of the dietary factor
  # 
  # 25. substitution.sims.wasted = impact from substitution of calories needed to 
  # maintain energy going from current to counterfactual for dietary factor of 
  # interest, but just for the wasted part of the dietary factor
  # 
  # 26. combined.impact.sims.wasted = combined impact of dietary factor shift and 
  # corresponding substitution impact (impact.sim.wasted + substitution.sims.wasted), 
  # but just for the wated part of the dietary factor
  # 
  # 27. delta.sims.wasted = difference in dietary factor amount produced for 
  # current vs counterfactual, but just for the inedible part of the dietary factor
  # 
  # 28. current.intake.impact.inedible = impact of current intake for dietary 
  # factor of interest, but just for the inedible part
  # 
  # 29. CF.intake.impact.inedible = impact of counterfactual intake for dietary 
  # factor of interest, but just for the inedible part
  # 
  # 30. current.intake.impact.wasted = impact of current intake for dietary 
  # factor of interest, but just for the wasted part
  # 
  # 31. CF.intake.impact.wasted = impact of counterfactual intake for dietary 
  # factor of interest, but just for the wasted part
  
  impact.list <- list("impact.sims.total" = impact.sims.total,
                    "substitution.sims.total" = substitution.sims.total,
                    "combined.impact.sims.total" = combined.impact.sims.total,
                    "population.sims.total" = population.sims.output,
                    "delta.sims.total" = delta.sims.output.total,
                    
                    "impact.sims.consumed" = impact.sims.consumed,
                    "substitution.sims.consumed" = substitution.sims.consumed,
                    "combined.impact.sims.consumed" = combined.impact.sims.consumed,
                    "delta.sims.consumed" = delta.sims.output.consumed,
                    
                    "impact.sims.unconsumed" = impact.sims.unconsumed,
                    "substitution.sims.unconsumed" = substitution.sims.unconsumed,
                    "combined.impact.sims.unconsumed" = combined.impact.sims.unconsumed,
                    "delta.sims.unconsumed" = delta.sims.output.unconsumed,
                    
                    "current.intake.impact.total" = current.impact.sims,
                    "CF.intake.impact.total" = counterfactual.impact.sims,
                    "current.intake.impact.consumed" = current.impact.sims.consumed,
                    "CF.intake.impact.consumed" = counterfactual.impact.sims.consumed,
                    "current.intake.impact.unconsumed" = current.impact.sims.unconsumed,
                    "CF.intake.impact.unconsumed" = counterfactual.impact.sims.unconsumed,
                    
                    "impact.sims.inedible" = impact.sims.inedible,
                    "substitution.sims.inedible" = substitution.sims.inedible,
                    "combined.impact.sims.inedible" = combined.impact.sims.inedible,
                    "delta.sims.inedible" = delta.sims.output.inedible,
                    
                    "impact.sims.wasted" = impact.sims.wasted,
                    "substitution.sims.wasted" = substitution.sims.wasted,
                    "combined.impact.sims.wasted" = combined.impact.sims.wasted,
                    "delta.sims.wasted" = delta.sims.output.wasted,
                    
                    "current.intake.impact.inedible" = current.impact.sims.inedible,
                    "CF.intake.impact.inedible" = counterfactual.impact.sims.inedible,
                    "current.intake.impact.wasted" = current.impact.sims.wasted,
                    "CF.intake.impact.wasted" = counterfactual.impact.sims.wasted
                    
                    )
  
  return(impact.list)
  
}

# Next is a function that calculates per capita sims by dividing impact sims with population count. 

get.percapita.sims <- function(impact.sims, population.sims, 
                               pop.sims.names, new.pop.sims.names, 
                               sims.names) {
  
  # replace X's with "pop"s
  names(population.sims)[names(population.sims) %in% pop.sims.names] <- 
    new.pop.sims.names 
  
  impact.sims.plus.pop <- merge(impact.sims, population.sims)
  
  impact.sims.percapita <- impact.sims.plus.pop
  
  impact.sims.percapita[, sims.names] <- 
    impact.sims.percapita[, sims.names] / 
    impact.sims.percapita[, new.pop.sims.names]
  
  impact.sims.percapita <- 
    impact.sims.percapita[, !(names(impact.sims.percapita) %in% new.pop.sims.names)]
  
  return(impact.sims.percapita)
  
}

# Next, a function that calculates summary stats from simulation outputs, 
# taking in simulation results for a specific diet factor as an input. 
# Returns a list of 22 data tables, each corresponding to a different outcome 
# (deltas, impact, substitution impact, combined impact, etc...). 
# The 22 data tables are spelled out below:
#   
# 1. delta.sims.summary = summary stats for change in amount of food produced 
# for dietary factor of interest
# 
# 2. impact.sims.summary = summary stats for impact of change from current to 
# counterfactual for dietary factor of interest
# 
# 3. substitution.impact.sims.summary = summary stats for impact from 
# substitution of calories needed to maintain energy going from current to 
# counterfactual for dietary factor of interest
# 
# 4. combined.impact.sims.summary = summary stats for combined impact of dietary 
# factor shift and corresponding substitution impact
# 
# 5. impact.sims.consumed.summary = summary stats for impact of change from 
# current to counterfactual for dietary factor of interest, but just for the 
# consumed part of the dietary factor
# 
# 6. substitution.impact.sims.consumed.summary = summary stats for combined 
# impact of dietary factor shift and corresponding substitution impact, but just 
# for the consumed part of the dietary factor 
# 
# 7. combined.impact.sims.consumed.summary = summary stats for combined impact 
# of dietary factor shift and corresponding substitution impact, but just for 
# the consumed part of the dietary factor
# 
# 8. impact.sims.unconsumed.summary = summary stats for impact of change from 
# current to counterfactual for dietary factor of interest, but just for the 
# unconsumed part of the dietary factor
# 
# 9. substitution.impact.sims.unconsumed.summary = summary stats for impact 
# from substitution of calories needed to maintain energy going from current to 
# counterfactual for dietary factor of interest, but just for the unconsumed 
# part of the dietary factor
# 
# 10. combined.impact.sims.unconsumed.summary = summary stats for combined 
# impact of dietary factor shift and corresponding substitution impact, but 
# just for the unconsumed part of the dietary factor
# 
# 11. current.intake.impact.sims.summary = summary stats for impact of current 
# intake for dietary factor of interest 
# 
# 12. current.intake.impact.consumed.sims.summary = summary stats for impact 
# of current intake for dietary factor of interest, but just for the consumed 
# part of the dietary factor
# 
# 13. current.intake.impact.unconsumed.sims.summary = summary stats for impact 
# of current intake for dietary factor of interest, but just for the unconsumed 
# part of the dietary factor
# 
# 14. CF.intake.impact.sims.summary = summary stats for impact of 
# counterfactual intake for dietary factor of interest 
# 
# 15. CF.intake.impact.consumed.sims.summary = summary stats for impact of 
# counterfactual intake for dietary factor of interest, but just for the 
# consumed part of the dietary factor
# 
# 16. CF.intake.impact.unconsumed.sims.summary = summary stats for impact of 
# counterfactual intake for dietary factor of interest, but just for the 
# unconsumed part of the dietary factor
# 
# 17. impact.sims.inedible.summary = summary stats for impact of change from 
# current to counterfactual for dietary factor of interest, but just for the 
# inedible part of the dietary factor
# 
# 18. substitution.impact.sims.inedible.summary = summary stats for impact 
# from substitution of calories needed to maintain energy going from current 
# to counterfactual for dietary factor of interest, but just for the inedible 
# part of the dietary factor
# 
# 19. combined.impact.sims.inedible.summary = summary stats for combined impact 
# of dietary factor shift and corresponding substitution impact, but just for the 
# inedible part of the dietary factor
# 
# 20. impact.sims.wasted.summary = summary stats for impact of change from 
# current to counterfactual for dietary factor of interest, but just for the 
# wasted part of the dietary factor
# 
# 21. substitution.impact.sims.wasted.summary = summary stats for impact from 
# substitution of calories needed to maintain energy going from current to 
# counterfactual for dietary factor of interest, but just for the wasted part 
# of the dietary factor
# 
# 22. combined.impact.sims.wasted.summary = summary stats for combined impact 
# of dietary factor shift and corresponding substitution impact, but just for 
# the wasted part of the dietary factor

get.summary.stats <- 
  function(impact.sims, substitution.impact.sims, combined.impact.sims, 
           delta.sims, impact.sims.consumed, substitution.impact.sims.consumed, 
           combined.impact.sims.consumed,impact.sims.unconsumed, 
           substitution.impact.sims.unconsumed, combined.impact.sims.unconsumed, 
           current.intake.impact.sims, current.intake.impact.sims.consumed, 
           current.intake.impact.sims.unconsumed, CF.intake.impact.sims, 
           CF.intake.impact.sims.consumed, CF.intake.impact.sims.unconsumed,
           impact.sims.inedible, substitution.impact.sims.inedible, 
           combined.impact.sims.inedible, impact.sims.wasted, 
           substitution.impact.sims.wasted, combined.impact.sims.wasted, 
           current.intake.impact.sims.inedible, current.intake.impact.sims.wasted, 
           CF.intake.impact.sims.inedible, CF.intake.impact.sims.wasted,
           vars, convert.to.kcal.units = "convert.to.kcal.units") {
    
  impact.sims.dt <- as.data.table(impact.sims)
  
  impact.sims.summary  <-  
    impact.sims.dt[, ":="(mean_impact = 
                            rowMeans(.SD, na.rm = TRUE),
                          
                          sd_impact = 
                            apply(.SD, 1, sd),
                          
                          median_impact = 
                            apply(.SD, 1, median),
                          
                          lowerci_95_impact = 
                            apply(.SD, 1, function(x) 
                              quantile(x, .025, na.rm = TRUE)),
                          
                          upperci_95_impact = 
                            apply(.SD, 1, function(x) 
                              quantile(x, .975, na.rm = TRUE))),
                   
                   .SDcols = vars
                   
  ]
  
  impact.sims.summary <- impact.sims.summary[, (vars) := NULL]
  
  substitution.impact.sims.dt <- as.data.table(substitution.impact.sims)
  
  substitution.impact.sims.summary <- 
    substitution.impact.sims.dt[, ":="(mean_substitution_impact = 
                                         rowMeans(.SD, na.rm = TRUE),
                                       
                                       sd_substitution_impact = 
                                         apply(.SD, 1, sd),
                                       
                                       median_substitution_impact = 
                                         apply(.SD, 1, median),
                                       
                                       lowerci_95_substitution_impact = 
                                         apply(.SD, 1, function(x) 
                                           quantile(x, .025, na.rm = TRUE)),
                                       
                                       upperci_95_substitution_impact = 
                                         apply(.SD, 1, function(x) 
                                           quantile(x, .975, na.rm = TRUE))),
                                       
                                .SDcols = vars
                                
  ]
  
  substitution.impact.sims.summary <- 
    substitution.impact.sims.summary[, (vars) := NULL]
  
  combined.impact.sims.dt <- 
    as.data.table(combined.impact.sims)
  
  combined.impact.sims.summary <- 
    combined.impact.sims.dt[, ":="(mean_combined_impact = 
                                     rowMeans(.SD, na.rm = TRUE),
                                                                 
                                   sd_combined_impact = 
                                     apply(.SD, 1, sd),
                                                                 
                                   median_combined_impact = 
                                     apply(.SD, 1, median),
                                                                 
                                   lowerci_95_combined_impact = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .025, na.rm = TRUE)),
                                                                 
                                   upperci_95_combined_impact = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .975, na.rm = TRUE))),
                                                          
                            .SDcols = vars
                            
  ]
  
  combined.impact.sims.summary <- 
    combined.impact.sims.summary[, (vars) := NULL]
  
  delta.sims.dt <- 
    as.data.table(delta.sims)
  
  delta.sims.summary <- 
    delta.sims.dt[, ":="(mean_delta = 
                           rowMeans(.SD, na.rm = TRUE),
                                             
                         sd_delta = apply(.SD, 1, sd),
                                             
                         median_delta = apply(.SD, 1, median),
                                             
                         lowerci_95_delta = 
                           apply(.SD, 1, function(x) 
                             quantile(x, .025, na.rm = TRUE)),
                                             
                         upperci_95_delta = 
                           apply(.SD, 1, function(x) 
                             quantile(x, .975, na.rm = TRUE))),
                                      
                  .SDcols = vars
  ]
  
  delta.sims.summary <- delta.sims.summary[, (vars) := NULL]
  
  # consumed
  impact.sims.consumed.dt <- as.data.table(impact.sims.consumed)
  
  impact.sims.consumed.summary <- 
    impact.sims.consumed.dt[, ":="(mean_impact_consumed = 
                                     rowMeans(.SD, na.rm = TRUE),
                                                                 
                                   sd_impact_consumed = 
                                     apply(.SD, 1, sd),
                                                                 
                                   median_impact_consumed = 
                                     apply(.SD, 1, median),
                                                                 
                                   lowerci_95_impact_consumed = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .025, na.rm = TRUE)),
                                                                 
                                   upperci_95_impact_consumed = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .975, na.rm = TRUE))),
                                                          
                            .SDcols = vars
                            
  ]
  
  impact.sims.consumed.summary <- 
    impact.sims.consumed.summary[, (vars) := NULL]
  
  substitution.impact.sims.consumed.dt <- 
    as.data.table(substitution.impact.sims.consumed)
  
  substitution.impact.sims.consumed.summary <- 
    substitution.impact.sims.consumed.dt[, ":="(mean_substitution_impact_consumed = 
                                                  rowMeans(.SD, na.rm = TRUE),
                                                                                           
                                                sd_substitution_impact_consumed = 
                                                  apply(.SD, 1, sd),
                                                                                          
                                                median_substitution_impact_consumed = 
                                                  apply(.SD, 1, median),
                                                                                           
                                                lowerci_95_substitution_impact_consumed = 
                                                  apply(.SD, 1, function(x) 
                                                    quantile(x, .025, na.rm = TRUE)),
                                                                                           
                                                upperci_95_substitution_impact_consumed = 
                                                  apply(.SD, 1, function(x) 
                                                    quantile(x, .975, na.rm = TRUE))),
                                                                                    
                                         .SDcols = vars
                                         
  ]
  
  substitution.impact.sims.consumed.summary <- 
    substitution.impact.sims.consumed.summary[, (vars) := NULL]
  
  combined.impact.sims.consumed.dt <- 
    as.data.table(combined.impact.sims.consumed)
  
  combined.impact.sims.consumed.summary <- 
    combined.impact.sims.consumed.dt[, ":="(mean_combined_impact_consumed = 
                                              rowMeans(.SD, na.rm = TRUE),
                                                                                   
                                            sd_combined_impact_consumed = 
                                              apply(.SD, 1, sd),
                                                                                   
                                            median_combined_impact_consumed = 
                                              apply(.SD, 1, median),
                                                                                   
                                            lowerci_95_combined_impact_consumed = 
                                              apply(.SD, 1, function(x) 
                                                quantile(x, .025, na.rm = TRUE)),
                                                                                   
                                            upperci_95_combined_impact_consumed = 
                                              apply(.SD, 1, function(x) 
                                                quantile(x, .975, na.rm = TRUE))),
                                                                            
                                     .SDcols = vars
                                     
  ]
  
  combined.impact.sims.consumed.summary <- 
    combined.impact.sims.consumed.summary[, (vars) := NULL]
  
  # unconsumed
  impact.sims.unconsumed.dt <- as.data.table(impact.sims.unconsumed)
  
  impact.sims.unconsumed.summary <- 
    impact.sims.unconsumed.dt[, ":="(mean_impact_unconsumed = 
                                       rowMeans(.SD, na.rm = TRUE),
                                                                     
                                     sd_impact_unconsumed = 
                                       apply(.SD, 1, sd),
                                                                     
                                     median_impact_unconsumed = 
                                       apply(.SD, 1, median),
                                                                     
                                     lowerci_95_impact_unconsumed = 
                                       apply(.SD, 1, function(x) 
                                         quantile(x, .025, na.rm = TRUE)),
                                                                     
                                     upperci_95_impact_unconsumed = 
                                       apply(.SD, 1, function(x) 
                                         quantile(x, .975, na.rm = TRUE))),
                                                              
                              .SDcols = vars
                              
  ]
  
  impact.sims.unconsumed.summary <- 
    impact.sims.unconsumed.summary[, (vars) := NULL]
  
  substitution.impact.sims.unconsumed.dt <- 
    as.data.table(substitution.impact.sims.unconsumed)
  
  substitution.impact.sims.unconsumed.summary <- 
    substitution.impact.sims.unconsumed.dt[, ":="(mean_substitution_impact_unconsumed = 
                                                    rowMeans(.SD, na.rm = TRUE),
                                                                                               
                                                  sd_substitution_impact_unconsumed = 
                                                    apply(.SD, 1, sd),
                                                                                               
                                                  median_substitution_impact_unconsumed = 
                                                    apply(.SD, 1, median),
                                                                                               
                                                  lowerci_95_substitution_impact_unconsumed = 
                                                    apply(.SD, 1, function(x) 
                                                      quantile(x, .025, na.rm = TRUE)),
                                                                                               
                                                  upperci_95_substitution_impact_unconsumed = 
                                                    apply(.SD, 1, function(x) 
                                                      quantile(x, .975, na.rm = TRUE))),
                                                                                        
                                           .SDcols = vars
                                           
  ]
  
  substitution.impact.sims.unconsumed.summary <- 
    substitution.impact.sims.unconsumed.summary[, (vars) := NULL]
  
  combined.impact.sims.unconsumed.dt <- 
    as.data.table(combined.impact.sims.unconsumed)
  
  combined.impact.sims.unconsumed.summary <- 
    combined.impact.sims.unconsumed.dt[, ":="(mean_combined_impact_unconsumed = 
                                                rowMeans(.SD, na.rm = TRUE),
                                                                                       
                                              sd_combined_impact_unconsumed = 
                                                apply(.SD, 1, sd),
                                                                                       
                                              median_combined_impact_unconsumed = 
                                                apply(.SD, 1, median),
                                                                                       
                                              lowerci_95_combined_impact_unconsumed = 
                                                apply(.SD, 1, function(x) 
                                                  quantile(x, .025, na.rm = TRUE)),
                                                                                       
                                              upperci_95_combined_impact_unconsumed = 
                                                apply(.SD, 1, function(x) 
                                                  quantile(x, .975, na.rm = TRUE))),
                                                                                
                                       .SDcols = vars
                                       
  ]
  
  combined.impact.sims.unconsumed.summary <- 
    combined.impact.sims.unconsumed.summary[, (vars) := NULL]
  
  # current intake
  current.intake.impact.sims.dt <- as.data.table(current.intake.impact.sims)
  
  current.intake.impact.sims.summary <- 
    current.intake.impact.sims.dt[, ":="(mean_impact_current_intake_ = 
                                           rowMeans(.SD, na.rm = TRUE),
                                                                             
                                         sd_impact_current_intake = 
                                           apply(.SD, 1, sd),
                                                                             
                                         median_impact_current_intake = 
                                           apply(.SD, 1, median),
                                                                             
                                         lowerci_95_impact_current_intake = 
                                           apply(.SD, 1, function(x) 
                                             quantile(x, .025, na.rm = TRUE)),
                                                                             
                                         upperci_95_impact_current_intake = 
                                           apply(.SD, 1, function(x) 
                                             quantile(x, .975, na.rm = TRUE))),
                                                                      
                                  .SDcols = vars
                                  
  ]
  
  current.intake.impact.sims.summary <- 
    current.intake.impact.sims.summary[, (vars) := NULL]
  
  current.intake.impact.consumed.sims.dt <- 
    as.data.table(current.intake.impact.sims.consumed)
  
  current.intake.impact.consumed.sims.summary <- 
    current.intake.impact.consumed.sims.dt[, ":="(mean_impact_current_intake_consumed = 
                                                    rowMeans(.SD, na.rm = TRUE),
                                                                                               
                                                  sd_impact_current_intake_consumed = 
                                                    apply(.SD, 1, sd),
                                                                                               
                                                  median_impact_current_intake_consumed = 
                                                    apply(.SD, 1, median),
                                                                                               
                                                  lowerci_95_impact_current_intake_consumed = 
                                                    apply(.SD, 1, function(x) 
                                                      quantile(x, .025, na.rm = TRUE)),
                                                                                               
                                                  upperci_95_impact_current_intake_consumed = 
                                                    apply(.SD, 1, function(x) 
                                                      quantile(x, .975, na.rm = TRUE))),
                                                                                        
                                           .SDcols = vars
                                           
  ]
  
  current.intake.impact.consumed.sims.summary <- 
    current.intake.impact.consumed.sims.summary[, (vars) := NULL]
  
  current.intake.impact.unconsumed.sims.dt <- 
    as.data.table(current.intake.impact.sims.unconsumed)
  
  current.intake.impact.unconsumed.sims.summary <- 
    current.intake.impact.unconsumed.sims.dt[, ":="(mean_impact_current_intake_unconsumed = 
                                                      rowMeans(.SD, na.rm = TRUE),
                                                                                                   
                                                    sd_impact_current_intake_unconsumed = 
                                                      apply(.SD, 1, sd),
                                                                                                   
                                                    median_impact_current_intake_unconsumed = 
                                                      apply(.SD, 1, median),
                                                                                                   
                                                    lowerci_95_impact_current_intake_unconsumed = 
                                                      apply(.SD, 1, function(x) 
                                                        quantile(x, .025, na.rm = TRUE)),
                                                                                                   
                                                    upperci_95_impact_current_intake_unconsumed = 
                                                      apply(.SD, 1, function(x) 
                                                        quantile(x, .975, na.rm = TRUE))),
                                                                                            
                                             .SDcols = vars
                                             
  ]
  
  current.intake.impact.unconsumed.sims.summary <- 
    current.intake.impact.unconsumed.sims.summary[, (vars) := NULL]
  
  # Counterfactual intake
  CF.intake.impact.sims.dt <- as.data.table(CF.intake.impact.sims)
  
  CF.intake.impact.sims.summary <- 
    CF.intake.impact.sims.dt[, ":="(mean_impact_CF_intake = 
                                      rowMeans(.SD, na.rm = TRUE),
                                                                   
                                    sd_impact_CF_intake = 
                                      apply(.SD, 1, sd),
                                                                   
                                    median_impact_CF_intake = 
                                      apply(.SD, 1, median),
                                                                   
                                    lowerci_95_impact_CF_intake = 
                                      apply(.SD, 1, function(x) 
                                        quantile(x, .025, na.rm = TRUE)),
                                                                   
                                    upperci_95_impact_CF_intake = 
                                      apply(.SD, 1, function(x) 
                                        quantile(x, .975, na.rm = TRUE))),
                                                            
                             .SDcols = vars
                             
  ]
  
  CF.intake.impact.sims.summary <- 
    CF.intake.impact.sims.summary[, (vars) := NULL]
  
  CF.intake.impact.consumed.sims.dt <- 
    as.data.table(CF.intake.impact.sims.consumed)
  
  CF.intake.impact.consumed.sims.summary <- 
    CF.intake.impact.consumed.sims.dt[, ":="(mean_impact_CF_intake_consumed = 
                                               rowMeans(.SD, na.rm = TRUE),
                                                                                     
                                             sd_impact_current_CF_consumed = 
                                               apply(.SD, 1, sd),
                                                                                     
                                             median_impact_CF_intake_consumed = 
                                               apply(.SD, 1, median),
                                                                                     
                                             lowerci_95_impact_CF_intake_consumed = 
                                               apply(.SD, 1, function(x) 
                                                 quantile(x, .025, na.rm = TRUE)),
                                                                                     
                                             upperci_95_impact_CF_intake_consumed = 
                                               apply(.SD, 1, function(x) 
                                                 quantile(x, .975, na.rm = TRUE))),
                                                                              
                                      .SDcols = vars
                                      
  ]
  
  CF.intake.impact.consumed.sims.summary <- 
    CF.intake.impact.consumed.sims.summary[, (vars) := NULL]
  
  CF.intake.impact.unconsumed.sims.dt <- 
    as.data.table(CF.intake.impact.sims.unconsumed)
  
  CF.intake.impact.unconsumed.sims.summary <- 
    CF.intake.impact.unconsumed.sims.dt[, ":="(mean_impact_CF_intake_unconsumed = 
                                                 rowMeans(.SD, na.rm = TRUE),
                                                                                         
                                               sd_impact_CF_intake_unconsumed = 
                                                 apply(.SD, 1, sd),
                                                                                         
                                               median_impact_CF_intake_unconsumed = 
                                                 apply(.SD, 1, median),
                                                                                         
                                               lowerci_95_impact_CF_intake_unconsumed = 
                                                 apply(.SD, 1, function(x) 
                                                   quantile(x, .025, na.rm = TRUE)),
                                                                                         
                                               upperci_95_impact_CF_intake_unconsumed = 
                                                 apply(.SD, 1, function(x) 
                                                   quantile(x, .975, na.rm = TRUE))),
                                                                                  
                                        .SDcols = vars
                                        
  ]
  
  CF.intake.impact.unconsumed.sims.summary <- 
    CF.intake.impact.unconsumed.sims.summary[, (vars) := NULL]
  
  # inedible
  impact.sims.inedible.dt <- as.data.table(impact.sims.inedible)
  
  impact.sims.inedible.summary <- 
    impact.sims.inedible.dt[, ":="(mean_impact_inedible = 
                                     rowMeans(.SD, na.rm = TRUE),
                                                                 
                                   sd_impact_inedible = 
                                     apply(.SD, 1, sd),
                                                                 
                                   median_impact_inedible = 
                                     apply(.SD, 1, median),
                                                                 
                                   lowerci_95_impact_inedible = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .025, na.rm = TRUE)),
                                                                 
                                   upperci_95_impact_inedible = 
                                     apply(.SD, 1, function(x) 
                                       quantile(x, .975, na.rm = TRUE))),
                                                          
                            .SDcols = vars
                            
  ]
  
  impact.sims.inedible.summary <- 
    impact.sims.inedible.summary[, (vars) := NULL]
  
  substitution.impact.sims.inedible.dt <- 
    as.data.table(substitution.impact.sims.inedible)
  
  substitution.impact.sims.inedible.summary <- 
    substitution.impact.sims.inedible.dt[, ":="(mean_substitution_impact_inedible = 
                                                  rowMeans(.SD, na.rm = TRUE),
                                                                                           
                                                sd_substitution_impact_inedible = 
                                                  apply(.SD, 1, sd),
                                                                                           
                                                median_substitution_impact_inedible = 
                                                  apply(.SD, 1, median),
                                                                                           
                                                lowerci_95_substitution_impact_inedible = 
                                                  apply(.SD, 1, function(x) 
                                                    quantile(x, .025, na.rm = TRUE)),
                                                                                           
                                                upperci_95_substitution_impact_inedible = 
                                                  apply(.SD, 1, function(x) 
                                                    quantile(x, .975, na.rm = TRUE))),
                                                                                    
                                         .SDcols = vars
                                         
  ]
  
  substitution.impact.sims.inedible.summary <- 
    substitution.impact.sims.inedible.summary[, (vars) := NULL]
  
  combined.impact.sims.inedible.dt <- 
    as.data.table(combined.impact.sims.inedible)
  
  combined.impact.sims.inedible.summary <- 
    combined.impact.sims.inedible.dt[, ":="(mean_combined_impact_inedible = 
                                              rowMeans(.SD, na.rm = TRUE),
                                                                                   
                                            sd_combined_impact_inedible = 
                                              apply(.SD, 1, sd),
                                                                                   
                                            median_combined_impact_inedible = 
                                              apply(.SD, 1, median),
                                                                                   
                                            lowerci_95_combined_impact_inedible = 
                                              apply(.SD, 1, function(x) 
                                                quantile(x, .025, na.rm = TRUE)),
                                                                                   
                                            upperci_95_combined_impact_inedible = 
                                              apply(.SD, 1, function(x) 
                                                quantile(x, .975, na.rm = TRUE))),
                                                                            
                                     .SDcols = vars
                                     
  ]
  
  combined.impact.sims.inedible.summary <- 
    combined.impact.sims.inedible.summary[, (vars) := NULL]
  
  # wasted
  impact.sims.wasted.dt <- as.data.table(impact.sims.wasted)
  
  impact.sims.wasted.summary <- 
    impact.sims.wasted.dt[, ":="(mean_impact_wasted = 
                                   rowMeans(.SD, na.rm = TRUE),
                                                             
                                 sd_impact_wasted = 
                                   apply(.SD, 1, sd),
                                                             
                                 median_impact_wasted = 
                                   apply(.SD, 1, median),
                                                             
                                 lowerci_95_impact_wasted = 
                                   apply(.SD, 1, function(x) 
                                     quantile(x, .025, na.rm = TRUE)),
                                                             
                                 upperci_95_impact_wasted = 
                                   apply(.SD, 1, function(x) 
                                     quantile(x, .975, na.rm = TRUE))),
                                                      
                          .SDcols = vars
                          
  ]
  
  impact.sims.wasted.summary <- 
    impact.sims.wasted.summary[, (vars) := NULL]
  
  substitution.impact.sims.wasted.dt <- 
    as.data.table(substitution.impact.sims.wasted)
  
  substitution.impact.sims.wasted.summary <- 
    substitution.impact.sims.wasted.dt[, ":="(mean_substitution_impact_wasted = 
                                                rowMeans(.SD, na.rm = TRUE),
                                                                                       
                                              sd_substitution_impact_wasted = 
                                                apply(.SD, 1, sd),
                                                                                       
                                              median_substitution_impact_wasted = 
                                                apply(.SD, 1, median),
                                                                                       
                                              lowerci_95_substitution_impact_wasted = 
                                                apply(.SD, 1, function(x) 
                                                  quantile(x, .025, na.rm = TRUE)),
                                                                                       
                                              upperci_95_substitution_impact_wasted = 
                                                apply(.SD, 1, function(x) 
                                                  quantile(x, .975, na.rm = TRUE))),
                                                                                
                                       .SDcols = vars
                                       
  ]
  
  substitution.impact.sims.wasted.summary <- 
    substitution.impact.sims.wasted.summary[, (vars) := NULL]
  
  combined.impact.sims.wasted.dt <- 
    as.data.table(combined.impact.sims.wasted)
  
  combined.impact.sims.wasted.summary <- 
    combined.impact.sims.wasted.dt[, ":="(mean_combined_impact_wasted = 
                                            rowMeans(.SD, na.rm = TRUE),
                                                                               
                                          sd_combined_impact_wasted = 
                                            apply(.SD, 1, sd),
                                                                               
                                          median_combined_impact_wasted = 
                                            apply(.SD, 1, median),
                                                                               
                                          lowerci_95_combined_impact_wasted = 
                                            apply(.SD, 1, function(x) 
                                              quantile(x, .025, na.rm = TRUE)),
                                                                               
                                          upperci_95_combined_impact_wasted = 
                                            apply(.SD, 1, function(x) 
                                              quantile(x, .975, na.rm = TRUE))),
                                                                        
                                   .SDcols = vars
                                   
  ]
  
  combined.impact.sims.wasted.summary <- 
    combined.impact.sims.wasted.summary[, (vars) := NULL]
  
  all.impact.sims.summary <- 
    Reduce(function(...) merge(..., all = TRUE), 
           list(delta.sims.summary, 
                impact.sims.summary, 
                substitution.impact.sims.summary, 
                combined.impact.sims.summary,
                impact.sims.consumed.summary, 
                substitution.impact.sims.consumed.summary, 
                combined.impact.sims.consumed.summary,
                impact.sims.unconsumed.summary, 
                substitution.impact.sims.unconsumed.summary, 
                combined.impact.sims.unconsumed.summary,
                current.intake.impact.sims.summary, 
                current.intake.impact.consumed.sims.summary, 
                current.intake.impact.unconsumed.sims.summary,
                CF.intake.impact.sims.summary, 
                CF.intake.impact.consumed.sims.summary, 
                CF.intake.impact.unconsumed.sims.summary,
                impact.sims.inedible.summary, 
                substitution.impact.sims.inedible.summary, 
                combined.impact.sims.inedible.summary,
                impact.sims.wasted.summary, 
                substitution.impact.sims.wasted.summary, 
                combined.impact.sims.wasted.summary))
  
  return(all.impact.sims.summary)
  
}

# Simple version of above function to use on just one set of outputs 
# (as opposed to 22). Inputs are:
# 
# 1. impact.sims: the data frame with output sims of interest, and 
# 
# 2. vars: the variable names for the sims (e.g, V1, V2, .... V1000).  

get.summary.stats.simple <- function(impact.sims, vars) {
  
  impact.sims.dt <- as.data.table(impact.sims)
  
  impact.sims.summary <- 
    impact.sims.dt[, ":="(mean_impact = rowMeans(.SD, na.rm = TRUE),
                                               
                          sd_impact = apply(.SD, 1, sd),
                                               
                          median_impact = apply(.SD, 1, median),
                                               
                          lowerci_95_impact = apply(.SD, 1, function(x) 
                            quantile(x, .025, na.rm = TRUE)),
                                               
                          upperci_95_impact = apply(.SD, 1, function(x) 
                            quantile(x, .975, na.rm = TRUE))),
                                        
                   .SDcols = vars
                   
  ]
  
  impact.sims.summary <- impact.sims.summary[, (vars) := NULL]
  
  return(impact.sims.summary)
  
}

# And finally, a function that sums the output by all strata combinations 
# of interest and outputs summary stats. Note that input "strata.combos.list" 
# is a list of the strata combinations of interest, with each element being a 
# string vector denoting the strata for that combination. "impact.sims.allgroups" 
# is data frame with output sims that you want to summarize.

summary.stats.by.strata <- function(impact.sims.allgroups, strata.combos.list, pop.sims) {
  
  # Convert sims to data tables, create some string vectors for labeling 
  # simulated values, and create lists for storing results.
  
  cols.to.sum <- paste("X", 1:n.sims, sep = "")
  pop.cols.to.sum <- paste("pop", 1:n.sims, sep = "")
  
  sims.data.table <- as.data.table(impact.sims.allgroups)
  pop.sims.data.table <- as.data.table(pop.sims)
  
  strata.combos.sims <- list()
  strata.combos.sims.long <- list()
  strata.combos.summary.stats.long <- list()
  strata.combos.summary.stats.wide <- list()
  strata.combos.summary.stats.misc <- list()
  
  pop.strata.combos.sims <- list() 
  
  percapita.sims <- list()
  percapita.sims.long <- list()
  percapita.summary.stats.long <- list()
  percapita.summary.stats.wide <- list()
  percapita.summary.stats.misc <- list()
  
  # Start a for-loop to traverse through each strata combination of interest.
  
  for(j in 1:length(strata.combos.list)) {
    
    # Sum by combinations of strata in strata.combos.list[[j]]
    
    labels <- c(strata.combos[[j]], "outcome", "outcome_unit")
    
    strata.combos.sims[[j]] <- 
      sims.data.table[, lapply(.SD, sum), 
                      by = c(strata.combos.list[[j]], "outcome", "outcome_unit"), 
                      .SDcols = cols.to.sum]
    
    # Sum population sims by combinations of strata in strata.combos.list[[j]]. 
    # Note that we can't sum population counts by "Foodgroup" or "broad_foodgroup" 
    # since we would be double (or triple, quadruple) counting (these are not 
    # population attributes) if we did that, so we remove those from the strata 
    # combination when summing population counts. If the strata combination of 
    # interest is just the foodgroups, then we sum over all subgroups.
    
    pop.strata <- 
      strata.combos[[j]][!(strata.combos.list[[j]] %in% c("Foodgroup", 
                                                          "broad_foodgroup"))]
    
    if(is.null(pop.strata)) {
      
      pop.strata.combos.sims[[j]] <- 
        pop.sims.data.table[, lapply(.SD, sum), 
                            .SDcols = cols.to.sum]
    
      } else {
      
        pop.strata.combos.sims[[j]] <- 
          pop.sims.data.table[, lapply(.SD, sum), 
                              by = c(pop.strata), 
                              .SDcols = cols.to.sum]
      
      }
    
    names(pop.strata.combos.sims[[j]]) <- 
      gsub(x = names(pop.strata.combos.sims[[j]]), pattern = "X", replacement = "pop")
    
    # merge the summed inputs and pop counts into one file
    
    if(length(intersect(names(strata.combos.sims[[j]]), names(pop.strata.combos.sims[[j]]))) == 0) {
      
      merged <- cbind(strata.combos.sims[[j]], pop.strata.combos.sims[[j]])
      
    }
    
    if(length(intersect(names(strata.combos.sims[[j]]), names(pop.strata.combos.sims[[j]]))) != 0) {
      
      merged <- merge(strata.combos.sims[[j]], pop.strata.combos.sims[[j]])
      
    }
    
    # Calculate per capita results of sims summed by strata
    
    merged.names <- merged[,..labels]
    percapita.sims[[j]] <- merged[,..cols.to.sum] / merged[,..pop.cols.to.sum]
    percapita.sims[[j]] <- cbind(merged.names, percapita.sims[[j]])
    
    # Reshape sims into long format, calculate summary stats, and label 
    # outcomes appropriately, then reshape back to wide format so that each 
    # summary stat has its own column.
    
    strata.combos.sims.long[[j]] <- 
      melt(strata.combos.sims[[j]], 
           id.vars = c(strata.combos[[j]], "outcome", "outcome_unit"), 
           measure.vars = cols.to.sum)
    
    strata.combos.summary.stats.long[[j]] <- 
      strata.combos.sims.long[[j]][, .(quantile(value, c(.025, .5, .975))), 
                                   by = c(strata.combos[[j]], "outcome", "outcome_unit")]
    
    strata.combos.summary.stats.misc[[j]] <- 
      strata.combos.sims.long[[j]][, .(mean(value), sd(value)), 
                                   by = c(strata.combos[[j]], "outcome")]
    
    names(strata.combos.summary.stats.misc[[j]])[(dim(strata.combos.summary.stats.misc[[j]])[2] - 1):
          dim(strata.combos.summary.stats.misc[[j]])[2]] <- 
      c("mean", "SD")
    
    strata.combos.summary.stats.long[[j]] <- 
      cbind(strata.combos.summary.stats.long[[j]], 
            as.factor(rep(c("lower_bound (2.5th percentile)", "median", 
                            "upper_bound  (97.5th percentile)"), 
                          dim(strata.combos.summary.stats.long[[j]])[1] / 3)))
    
    strata.combos.summary.stats.wide[[j]] <- 
      dcast(strata.combos.summary.stats.long[[j]], 
            paste(paste(strata.combos[[j]], collapse = ' + '), 
                  " + outcome + outcome_unit", " ~ V2", sep = ""),
            value.var = "V1")    

    strata.combos.summary.stats.wide[[j]] <- 
      merge(strata.combos.summary.stats.wide[[j]], 
            strata.combos.summary.stats.misc[[j]], 
            by = intersect(names(strata.combos.summary.stats.wide[[j]]), 
                         names(strata.combos.summary.stats.misc[[j]])))
    
    # Repeat for per capita sims.

    percapita.sims.long[[j]] <- 
      melt(percapita.sims[[j]], 
           id.vars = c(strata.combos.list[[j]], "outcome", "outcome_unit"), 
           measure.vars = cols.to.sum)
    
    percapita.summary.stats.long[[j]] <- 
      percapita.sims.long[[j]][, .(quantile(value, c(.025, .5, .975))), 
                               by = c(strata.combos.list[[j]], "outcome", "outcome_unit")]
    
    percapita.summary.stats.misc[[j]] <- 
      percapita.sims.long[[j]][, .(mean(value), sd(value)), 
                               by = c(strata.combos.list[[j]], "outcome")]
    
    names(percapita.summary.stats.misc[[j]])[(dim(strata.combos.summary.stats.misc[[j]])[2] - 1):
                                               dim(strata.combos.summary.stats.misc[[j]])[2]] <- 
      c("mean", "SD")
    
    percapita.summary.stats.long[[j]] <- 
      cbind(percapita.summary.stats.long[[j]], 
            as.factor(rep(c("lower_bound (2.5th percentile)", "median", 
                            "upper_bound  (97.5th percentile)"), 
                          dim(percapita.summary.stats.long[[j]])[1] / 3)))
    
    percapita.summary.stats.wide[[j]] <- 
      dcast(percapita.summary.stats.long[[j]], 
            paste(paste(strata.combos.list[[j]], collapse = ' + '), 
                  " + outcome + outcome_unit", " ~ V2", sep = ""),
            value.var = "V1")    
    
    percapita.summary.stats.wide[[j]] <- 
      merge(percapita.summary.stats.wide[[j]], 
            percapita.summary.stats.misc[[j]], 
            by = intersect(names(percapita.summary.stats.wide[[j]]), 
                         names(percapita.summary.stats.misc[[j]])))

  } # end loop
  
  # After looping through all strata combos of interest, we return the result 
  # as a 4 element list. 1 = summary output, 2 = full sims, 
  # 3 = per capita summary output, 4 = per capita full sims. 
  # Each of these 4 elements are themselves lists with with element 
  # corresponding to a strata combination of interest.
  
  return(list("summary.output" = strata.combos.summary.stats.wide, 
              "full.sims.output" = strata.combos.sims,
              "summary.output.percapita" = percapita.summary.stats.wide,
              "full.sims.output.percapita" = percapita.sims)
         )
  
}


