# Title: Part 3 - Run PIF analysis
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-15-2025

############################################################################

### Documentation for R script ###

# You'll read in this script in "Complete.PIF.Analysis.partTwo.cluster.r". It 
# defines a bunch of function used for summarizing the output (and a couple 
# obsolete functions that have seen been replaced. Since this code is just 
# defining functions, there is no harm in not commenting them out). 

############################################################################

# This function takes in alldisease (simulation of attributable mortality/incidence) 
# and total disease (simulation of total mortality/incidence) to calculate 
# "reverse-engineered" PIFs, which is simply the attributable mortality/incidence 
# divided by total mortality/incidence. What's the point when we calculate the PIFs 
# directly in part one? We want to be able to report PIFs for broader subgroups 
# (not just by age/sex/race) and this is the way to do that. First step in the 
# function is two merge totaldisease and alldisease (to ensure everything matches up 
# by subgroup), next is to identify column indices for the totaldisease and alldisease 
# simulations, and then finally use those indices to take alldisease sims and divide 
# by total disease sims. Note that this function only works when your subgroups 
# are age/sex/race specific. Output is sims for reverse-engineered PAFs for 
# all subgroups.

# function of calculating RE PAFs sum stats
# for now, only works for age/sex/race
get.RE.PAFs.draws <- function(alldisease, totaldisease) {
  
  merged <- merge(alldisease, totaldisease)
  
  alldisease.starting.point <- which(names(merged) == "V1")
  alldisease.ending.point <- which(names(merged) == "riskfactor") - 1
  totaldisease.starting.point <- which(names(merged) == "X1")
  totaldisease.ending.point <- dim(merged)[2]
  
  RE.PAFs.info <- 
    merged[,which(names(merged) %in% c("age", "female", "race", "outcome", 
                                       "mean..food.", "se..food", "sd..food.", 
                                       "riskfactor", "count", "se"))]
  
  RE.PAFs <- cbind(merged[,which(names(merged) %in% c("age", "female", "race", 
                                                      "outcome", "mean..food.", 
                                                      "se..food", "sd..food.", 
                                                      "riskfactor", "count", 
                                                      "se"))], 
                 merged[,alldisease.starting.point:alldisease.ending.point] / 
                   merged[,totaldisease.starting.point:totaldisease.ending.point])
  
  RE.PAFs <- RE.PAFs[order(RE.PAFs$riskfactor, RE.PAFs$outcome, 
                           RE.PAFs$female, RE.PAFs$age, RE.PAFs$race),]
  
  return(RE.PAFs)
  
}

# Literally the same thing as above, but with the added step of getting summary 
# stats from the simulated RE PIFs.

# for now, only works for age/sex/race
get.RE.PAFs <- function(alldisease, totaldisease) {
  
  merged <- merge(alldisease, totaldisease)
  
  alldisease.starting.point <- which(names(merged) == "V1")
  alldisease.ending.point <- which(names(merged) == "riskfactor") - 1
  totaldisease.starting.point <- which(names(merged) == "X1")
  totaldisease.ending.point <- dim(merged)[2]
  
  RE.PAFs.info <- 
    merged[,which(names(merged) %in% c("age", "female", "race", "outcome", 
                                       "mean..food.", "se..food", "sd..food.", 
                                       "riskfactor", "count", "se"))]
  
  RE.PAFs <- cbind(merged[,which(names(merged) %in% c("age", "female", "race", 
                                                      "outcome", "mean..food.", 
                                                      "se..food", "sd..food.", 
                                                      "riskfactor", "count", 
                                                      "se"))], 
                 merged[,alldisease.starting.point:alldisease.ending.point] / 
                   merged[,totaldisease.starting.point:totaldisease.ending.point])
  
  RE.PAFs <- RE.PAFs[order(RE.PAFs$riskfactor, RE.PAFs$outcome, 
                           RE.PAFs$female, RE.PAFs$age, RE.PAFs$race),]
  
  RE.PAFs.sum.stats <- sum.stats.maker.RE.PAFs(RE.PAFs, n.covar)
  
  return(RE.PAFs.sum.stats)
  
}

# Next is a function we'll use often. It'll read in a simulation file (like 
# attributable mortality sims) and output summary stats (mean, standard 
# deviation, median, 2.5th and 97.5th percentile).

# function to get summary stats
sum.stats.maker <- function(alldisease) {
  
  alldisease.dt <- as.data.table(alldisease)
  
  vars <- paste("V", 1:n.sims, sep = "")
  
  disease.summary <- 
    alldisease.dt[, ":="(means = rowMeans(.SD, na.rm = TRUE),
                      sd.devs = apply(.SD, 1, sd),
                      medians = apply(.SD, 1, median),
                      LB = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                      UB = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))
  ),
  
  .SDcols = vars
  
  ]
  
  disease.summary <- disease.summary[, (vars) := NULL]
  
  return(disease.summary)
  
}

# Next, a similar function that is used for Reverse engineered PAF draws. In 
# addition for getting the five summary stats mentioned above, it also outputs 
# median, 2.5th percentile, and 97.5th percentile as a percent (as opposed to as 
# a proportion) for convenience. 

sum.stats.maker.RE.PAFs <- function(RE.PAFs) {
  
  RE.PAFs.dt <- as.data.table(RE.PAFs)
  vars <- paste("V", 1:n.sims, sep = "")
  
  RE.PAFs.summary <- 
    RE.PAFs.dt[, ":="(RE.PAF.means = rowMeans(.SD, na.rm = TRUE),
                      RE.PAF.sd.devs = apply(.SD, 1, sd),
                      RE.PAF.medians = apply(.SD, 1, median),
                      RE.PAF.LB = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                      RE.PAF.UB = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T)),
                      RE.PAF.medians.as.percent = apply(.SD, 1, median) * 100,
                      RE.PAF.LB.as.percent = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)) * 100,
                      RE.PAF.UB.as.percent = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T)) * 100
  ),
  
  .SDcols = vars
  
  ]
  
  RE.PAFs.summary <- RE.PAFs.summary[, (vars) := NULL]
  
  return(RE.PAFs.summary)
  
}

# Next is Sum.by.strata. Very important function that produces output by strata 
# of interest (e.g, by age instead of by age/sex/race, etc...). Inputs are 
# simulation output for attributable disease (alldisease), total disease (total disease), 
# total population (pop), and finally, a string vector of the strata of interest 
# (covar). This function is a bit more steps than function above so I'll break it 
# down chunk by chunk, by the general idea is simple: Step 1 is to sum over 
# strata of interest, Step 2 is calculate summary stats for those strata of interest.

# debug
# alldisease = alldisease
# totaldisease = totaldisease
# pop = pop.draws
# covar = strata.combos[[i]]

Sum.by.strata <- function(alldisease, totaldisease, pop, covar) {
  
  # First, convert input data frames to data tables and get string of variable 
  # names for the simulations for each file. (note, we also rename our string 
  # vector of strata to "strata")
  
  alldisease.dt <- as.data.table(alldisease)
  totaldisease.dt <- as.data.table(totaldisease)
  pop.dt <- as.data.table(pop)
  
  strata <- covar
  
  cols.to.sum.alldisease <- paste("V", 1:n.sims, sep = "")
  cols.to.sum.totaldisease <- paste("X", 1:n.sims, sep = "")
  cols.to.sum.pop <- paste("Y", 1:n.sims, sep = "")
  
  # Following code sums the simulations over strata of interest. 
  # 
  # ".SDcols: are your variables of interest (in our case, the simulated 
  # values: V1, V2, etc...). 
  # 
  # "lapply(.SD, sum)" is the function you're applying: in our case summing 
  # values in each column. Whatever you specify in "by = " is the grouping of 
  # interest. Let's say it's age, as an example. For attributable mortality/incidence, 
  # you're summing by age and risk factor, because you don't want to sum results 
  # for different foods together. For total mortality, you only have data by 
  # age/sex/race and outcome/disease_type so whatever you listed in "covar", 
  # you are limited to those five groupings (that's what 
  # c(strata[strata %in% c("age", "female", "race", "outcome", "disease_type")] 
  # does). Total population numbers aren't divisible by outcome/disease type, so 
  # you're strata are limited to some combination of age/sex/race.
  
  strata.sims.alldisease <- 
    alldisease.dt[,lapply(.SD, sum), 
               by = c(strata, "riskfactor"), 
               .SDcols = cols.to.sum.alldisease]
  
  strata.sims.totaldisease <- 
    totaldisease.dt[,lapply(.SD, sum), 
                 by = c(strata[strata %in% c("age", "female", "race", "outcome", "disease_type", "disease_type1", "disease_type2", "disease_type3")]), 
                 .SDcols = cols.to.sum.totaldisease]
  
  strata.sims.pop <- pop.dt[,lapply(.SD, sum), 
                            by = c(strata[strata %in% c("age", "female", "race")]), 
                            .SDcols = cols.to.sum.pop]
  
  # Next, we merge three output data tables together. How to merge them 
  # correctly depends on the what strata you are interested, hence the four 
  # if statements. Everything in the if statement is doing the same thing 
  # (merging the three outputs together), but they need to be implemented 
  # differently depending on whether one or more of the data tables are one 
  # dimensional or not (merge function basically doesn't work if one of the 
  # things to be merged is one dimensional). So, the first if statement, where 
  # none of the three are one dimentionsal, then we can use the merge function 
  # without issue. In the other three cases, where either strata.sims.totaldisease 
  # is one dimensional, we need to use cbind instead o of merge to bind the 
  # files together. So if strata.sims.pop is one dimensional 
  # (if length(strata[strata %in% c("age", "female", "race")]) == 0  is true) 
  # and/or if strata.sims.alldisease is one dimensional 
  # (if dim(strata.sims.totaldisease)[1]>1), 
  # then we need to use cbind to combine those to the rest of the files. 
  
    if(length(strata[strata %in% c("age", "female", "race")]) > 0 & 
       dim(strata.sims.totaldisease)[1] > 1) {
      
      merged <- 
        merge(merge(strata.sims.alldisease, strata.sims.totaldisease), 
              strata.sims.pop, 
              by = strata[strata %in% c("age", "female", "race")])
      
    }
  
    if(length(strata[strata %in% c("age", "female", "race")]) > 0 & 
       dim(strata.sims.totaldisease)[1] == 1) {
      
      merge.a <- cbind(strata.sims.alldisease, strata.sims.totaldisease)
      
      merged <- merge(merge.a, 
                      strata.sims.pop, 
                      by = strata[strata %in% c("age", "female", "race")])
      
    }
  
    if(length(strata[strata %in% c("age", "female", "race")]) == 0 & 
       dim(strata.sims.totaldisease)[1] > 1) {
      
      merged <- cbind(merge(strata.sims.alldisease, strata.sims.totaldisease), 
                      strata.sims.pop)
      
    }
  
    if(length(strata[strata %in% c("age", "female", "race")]) == 0 & 
       dim(strata.sims.totaldisease)[1] == 1){
      
      merged <- cbind(strata.sims.alldisease, 
                      strata.sims.totaldisease, 
                      strata.sims.pop)
      
    }
  
  # Next is identifying the index for starting and ending of simulation values 
  # for each of the three simulated outputs.
  
  alldisease.starting.point <- which(names(merged) == "V1")
  alldisease.ending.point <- which(names(merged) == paste("V", nsim1, sep = ""))
  totaldisease.starting.point <- which(names(merged) == "X1")
  totaldisease.ending.point <- which(names(merged) == paste("X", nsim1, sep = ""))
  pop.starting.point <- which(names(merged) == "Y1")
  pop.ending.point <- which(names(merged) == paste("Y", nsim1, sep = ""))
  
  # Along the same lines as described previously (in part two documentation): 
  # calculate "reverse-engineered" PIFs and "standardized disease", 
  # and get summary stats.
  
  RE.PAFs <- cbind(merged[, c(strata, "riskfactor"), with = FALSE], 
                 merged[, alldisease.starting.point:alldisease.ending.point] / 
                   merged[, totaldisease.starting.point:totaldisease.ending.point])
  
  setorderv(RE.PAFs, cols=c("riskfactor", strata))
  
  sum.stats.attr.disease.PAF <- sum.stats.maker(strata.sims.alldisease)
  
  RE.PAFs.sum.stats <- sum.stats.maker.RE.PAFs(RE.PAFs)
  
  standardized.disease <- cbind(merged[, c(strata, "riskfactor"), with = FALSE], 
                           merged[, alldisease.starting.point:alldisease.ending.point] / 
                             merged[,pop.starting.point:pop.ending.point] * 100000)
  
  setorderv(standardized.disease, cols = c("riskfactor", strata))
  
  standardized.disease.sum.stats <- sum.stats.maker(standardized.disease)
  
  names(standardized.disease.sum.stats)[(length(names(standardized.disease.sum.stats)) - 4):
                                       length(names(standardized.disease.sum.stats))] <- 
    paste("disease.per.100k.", 
          names(standardized.disease.sum.stats)[(length(names(standardized.disease.sum.stats)) - 4):
                                               length(names(standardized.disease.sum.stats))], 
          sep = "")
  
  # Merge all the summary output we created (attributable disease, RE.PIF, 
  # standardized disease) into one super output file which we will call 
  # "strata.summary.stats". 
  
  #merged.results <- merge(sum.stats.attr.disease.PAF, RE.PAFs.sum.stats)
  merged.results <- reduce(list(sum.stats.attr.disease.PAF, 
                                RE.PAFs.sum.stats, 
                                standardized.disease.sum.stats), 
                           merge)
  
  strata.summary.stats <- merged.results
  
  # Finally, return all output of interest as a list of five elements: 
  # 1. summary stats for all three by specified strata of interest, 
  # 2. simulated values for attributable mortality/incidence by specified strata of interest, 
  # 3. simulated values total mortality/incidence by specified strata of interest, 
  # 4. "Reverse-engineered" PIFs by specified strata of interest, 
  # 5. standardized attributable mortality by specified strata of interest,
  
  return(list(strata.summary.stats, 
              strata.sims.alldisease, 
              strata.sims.totaldisease, 
              RE.PAFs, 
              standardized.disease))
  
}

# function to use to get summary stats by strata for joint PIFs (that is, 
# for combination of shifts of multiple dietary factors). It is now identical 
# to the Sum.by.strata function. (It used to be different, before sum.by.strata 
# was made to be more flexible). Easier to change the function than change 
# every instance of Sum.by.strata.joint being used in part three.

Sum.by.strata.joint <- function(alldisease, totaldisease, pop, covar) {
  
  alldisease.dt <- as.data.table(alldisease)
  totaldisease.dt <- as.data.table(totaldisease)
  pop.dt <- as.data.table(pop)
  
  strata <- covar
  
  cols.to.sum.alldisease <- paste("V", 1:n.sims, sep = "")
  cols.to.sum.totaldisease <- paste("X", 1:n.sims, sep = "")
  cols.to.sum.pop <- paste("Y", 1:n.sims, sep = "")
  
  strata.sims.alldisease <- alldisease.dt[, lapply(.SD, sum), 
                                    by = c(strata), 
                                    .SDcols = cols.to.sum.alldisease]
  
  strata.sims.totaldisease <- 
    totaldisease.dt[,lapply(.SD, sum), 
                 by = c(strata[strata %in% c("age", "female", "race", 
                                           "outcome", "disease_type",
                                           "disease_type1", "disease_type2", 
                                           "disease_type3")]), 
                 .SDcols = cols.to.sum.totaldisease]
  
  strata.sims.pop <- 
    pop.dt[,lapply(.SD, sum), 
           by = c(strata[strata %in% c("age", "female", "race")]), 
           .SDcols = cols.to.sum.pop]
  
  if(length(strata[strata %in% c("age", "female", "race")]) > 0 & 
     dim(strata.sims.totaldisease)[1] > 1) {
    
    merged <- merge(merge(strata.sims.alldisease, strata.sims.totaldisease), 
                    strata.sims.pop, 
                    by = strata[strata %in% c("age", "female", "race")])
    
  }
  
  if(length(strata[strata %in% c("age", "female", "race")]) > 0 & 
     dim(strata.sims.totaldisease)[1] == 1) {
    
    merge.a <- cbind(strata.sims.alldisease, strata.sims.totaldisease)
    
    merged <- merge(merge.a, 
                    strata.sims.pop, 
                    by = strata[strata %in% c("age", "female", "race")])
    
  }
  
  if(length(strata[strata %in% c("age", "female", "race")]) == 0 & 
     dim(strata.sims.totaldisease)[1] > 1) {
    
    merged <- cbind(merge(strata.sims.alldisease, strata.sims.totaldisease), 
                    strata.sims.pop)
    
  }
  
  if(length(strata[strata %in% c("age", "female", "race")]) == 0 & 
     dim(strata.sims.totaldisease)[1] == 1) {
    
    merged <- cbind(strata.sims.alldisease, strata.sims.totaldisease, strata.sims.pop)
    
  }
  
  alldisease.starting.point <- which(names(merged) == "V1")
  alldisease.ending.point <- which(names(merged) == paste("V", nsim1, sep = ""))
  totaldisease.starting.point <- which(names(merged) == "X1")
  totaldisease.ending.point <- which(names(merged) == paste("X", nsim1, sep = ""))
  pop.starting.point <- which(names(merged) == "Y1")
  pop.ending.point <- which(names(merged) == paste("Y", nsim1, sep = ""))
  
  RE.PAFs <- cbind(merged[, c(strata), with = FALSE], 
                 merged[, alldisease.starting.point:alldisease.ending.point] / 
                   merged[, totaldisease.starting.point:totaldisease.ending.point])
  
  
  setorderv(RE.PAFs, cols = c(strata))
  
  sum.stats.attr.disease.PAF <- sum.stats.maker(strata.sims.alldisease)
  
  RE.PAFs.sum.stats <- sum.stats.maker.RE.PAFs(RE.PAFs)
  
  standardized.disease <- cbind(merged[, c(strata), with = FALSE], 
                           merged[, alldisease.starting.point:alldisease.ending.point] / 
                             merged[, pop.starting.point:pop.ending.point] * 100000)
  
  setorderv(standardized.disease, cols = c(strata))
  
  standardized.disease.sum.stats <- sum.stats.maker(standardized.disease)
  
  names(standardized.disease.sum.stats)[(length(names(standardized.disease.sum.stats)) - 4):
                                       length(names(standardized.disease.sum.stats))] <- 
    paste("disease.per.100k.", 
          names(standardized.disease.sum.stats)[(length(names(standardized.disease.sum.stats)) - 4):
                                               length(names(standardized.disease.sum.stats))], 
          sep = "")
  
  merged.results <- 
    reduce(list(sum.stats.attr.disease.PAF, RE.PAFs.sum.stats, standardized.disease.sum.stats), 
           merge, 
           by = strata)
  
  strata.summary.stats <- merged.results
  
  return(list(strata.summary.stats, 
              strata.sims.alldisease, 
              strata.sims.totaldisease, 
              RE.PAFs, 
              standardized.disease))
  
}

