# Title: Part 3 - Run PIF analysis
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-15-2025

############################################################################

### Documentation for R script ###

# This is documentation for "part 3" of the CRA code, which is for calculating 
# the joint PIFs (attributable mort/incidence of multiple concurrent 
# counterfactual exposures).

############################################################################

# Define covar.vec.string (string vector of covariates of interest) so we read 
# in the correct file (defining n.covar is obsolete).

covar.vec.string <- paste(covar.vec, sep = "", collapse = "")

n.covar <- length(covar.vec)

if(identical(covar.vec, c("Overall")))
  n.covar <- 0

# Next, we'll define "calc.joint.PAFs", the function that is actually going to 
# calculate the joint PIFs. It reads in a PIF simulations for a given dietary 
# pattern, and calculates joint PIFs for the defined strata of interest.

calc.joint.PAFs <- function(strata) {
  
  # Read in PAF simulation file of interest. 
  # Note that the function itself only reads in strata, so output.location, 
  # diet_pattern, year.vec.string, and covar.vec.string need to be defined 
  # before calling the function.
  
  PAF.draws <- read_csv(file = paste0(output_location, diet_pattern, 
                                   "/all.PIFs_", year.vec.string, "_",
                                   covar.vec.string, "_", diet_pattern, ".csv"))
  
  # If strata of interest includes one of the 4 ways of classifying disease 
  # (disease_type, disease_type1, disease_type2, disease_type3) then we get rid 
  # of the others from PAF.draws. These four are not mutually exclusive, they 
  # are different ways to the disease. Therefore, it doesn't make sense to 
  # stratify by more than one of these (Code won't work unless you take out the 
  # ones you don't use).
  
  # if contains disease type, then get rid of others
  if("disease_type" %in% strata) {
    
    PAF.draws <- PAF.draws %>% 
      select(-c("disease_type1", "disease_type2", "disease_type3"))
    
  } else if("disease_type1" %in% strata) {
    
    PAF.draws <- PAF.draws %>% 
      select(-c("disease_type", "disease_type2", "disease_type3"))
    
  } else if("disease_type2" %in% strata) {
    
    PAF.draws <- PAF.draws %>% 
      select(-c("disease_type", "disease_type1", "disease_type3"))
    
  } else if("disease_type3" %in% strata) {
    
    PAF.draws <- PAF.draws %>% 
      select(-c("disease_type", "disease_type1", "disease_type2"))
    
  } else {
    
    PAF.draws <- PAF.draws %>% 
      select(-c("disease_type1", "disease_type2", "disease_type3"))
    
  }
  
  # Define number of sims (by counting up number of variables that start with "V")
  n.draws <- length(grep(pattern = "V", x = names(PAF.draws)))
  
  PAF.draws$strata <- interaction(PAF.draws[,c(strata)])
  
  PAF.draws.by.strata.outcome <- split(x = PAF.draws, f = list(PAF.draws$strata))
  
  joint.PAFs.all.draws <- matrix(nrow = length(PAF.draws.by.strata.outcome), 
                                 # adding plus 1 for pathway
                                 ncol = (length(strata)+n.draws)) 
  
  joint.PAFs.all.draws <- as.data.frame(joint.PAFs.all.draws)
  
  names(joint.PAFs.all.draws)[1:(length(strata)+n.draws)] <- 
    c(strata, paste("V", 1:n.draws, sep = ""))  
  
  start <- which(names(PAF.draws.by.strata.outcome[[1]]) == "V1")
  
  end <- start+n.draws-1
  
  for(i in 1:length(PAF.draws.by.strata.outcome)) {
    
    joint.PAFs.all.draws[i,1:(length(strata))] <- 
      PAF.draws.by.strata.outcome[[i]][1, strata]
   
     joint.PAFs.all.draws[i, length(strata)+c(1:n.draws)] <- 
       1-apply(1-PAF.draws.by.strata.outcome[[i]][, start:end], 
               MARGIN = 2, FUN = prod)
    
    }
  
  return(joint.PAFs.all.draws)
  
}

# caclulates joint mort, by multilpying joint paf with total mort
calc.joint.mort <- function(joint.paf.draws, total.mort.draws, strata) {
  
  merged <- merge(joint.paf.draws, total.mort.draws)
  
  Joint.PAF.starting.point <- which(names(merged) == "V1")
  Joint.PAF.ending.point <- which(names(merged) == paste("V", nsim1, sep = ""))
  totalmort.starting.point <- which(names(merged) == "X1")
  totalmort.ending.point <- which(names(merged) == paste("X", nsim1, sep = ""))
  
  joint.mort <- cbind(merged[,c(strata)],
                      merged[,Joint.PAF.starting.point:Joint.PAF.ending.point] * 
                        merged[,totalmort.starting.point:totalmort.ending.point])
  return(joint.mort)
  
}

# create new combinations of strata without disease type overlap
strata0 <- c("age", "female", "race", "pathway", "disease_type")
strata1 <- c("age", "female", "race", "pathway", "disease_type1")
strata2 <- c("age", "female", "race", "pathway", "disease_type2")
strata3 <- c("age", "female", "race", "pathway", "disease_type3")
strata4 <- c("age", "female", "race", "outcome", "pathway")
strata5 <- c("age", "female", "race", "pathway")
strata_list <- list(strata0, strata1, strata2, strata3, strata4, strata5)

strata_list

# Once we know which dietary factors we will keep to calculate joint PIFs, 
# we will subset PAF.draws to keep only those dietary factors

library(magrittr)

joint.PAFs.all.draws <- list()
total.mort <- list()
joint.mort.all.draws <- list()

# create 'joint' directory if it doesn't exist
ifelse(!dir.exists(file.path(paste0(output_location, diet_pattern, "/joint"))),
       dir.create(file.path(paste0(output_location, diet_pattern, "/joint"))),
       "Directory Exists")

for (i in 1:length(strata_list)) {
  
  strata <- strata_list[[i]]
  
  print(strata)
  
  joint.PAFs.all.draws[[i]] <- calc.joint.PAFs(strata = strata)
  
  write.csv(x = joint.PAFs.all.draws[[i]], 
            file = paste(output_location, diet_pattern, 
                       "/joint/joint_PIFs_all_draws_", paste(strata, collapse = "_"), 
                       "_", year.vec.string, "_", diet_pattern, ".csv", sep = ""))
  
  total.mort[[i]] <- read.csv(paste0(output_location, 
                                    diet_pattern, 
                                    "/total_USmortality_draws_by_", 
                                    strata %>% 
                                      # setdiff(., "pathway") %>%
                                      paste(collapse = "_"),
                                    "_", year.vec, "_", 
                                    diet_pattern, ".csv"))
  
  joint.mort.all.draws[[i]] <- 
    calc.joint.mort(joint.paf.draws = joint.PAFs.all.draws[[i]], 
                    total.mort.draws = total.mort[[i]], 
                    strata = strata)
  
  write.csv(x = joint.mort.all.draws[[i]], 
            file = paste0(output_location, diet_pattern, 
                       "/joint/joint_mort_all_draws_", 
                       paste(strata, collapse = "_"), "_", 
                       year.vec.string, "_", diet_pattern, 
                       ".csv"))
  
}

# calculate combos of interest
# need to create new strata.combos.of.interest.index for each disease type category

# first look at all strata combos
print(strata.combos)

# create new list of relevant strata combos
strata.combos.new <- strata.combos[strata.combos.of.interest.index]

print(strata.combos.new) # good

strata.combos.new0 <- keep(strata.combos.new, function(x) "disease_type" %in% x)
strata.combos.new1 <- keep(strata.combos.new, function(x) "disease_type1" %in% x)
strata.combos.new2 <- keep(strata.combos.new, function(x) "disease_type2" %in% x)
strata.combos.new3 <- keep(strata.combos.new, function(x) "disease_type3" %in% x)
strata.combos.new4 <- keep(strata.combos.new, function(x) "outcome" %in% x)

strata.combos.new5 <- discard(strata.combos.new, 
                              function(x) "disease_type" %in% x | 
                                "disease_type1" %in% x | 
                                "disease_type2" %in% x | 
                                "disease_type3" %in% x | 
                                "outcome" %in% x)

# new function
calculate_combos <- function(my.strata.combos, x) {
  
  for(i in 1:length(my.strata.combos)) {
    
    print(i)
    print(my.strata.combos[[i]])
    print("")
    
    sum.stats.byX.joint.attr.mort <- 
      Sum.by.strata.joint(allmort = joint.mort.all.draws[[x]], 
                          totalmort = totalmort, 
                          pop = pop.draws, 
                          covar = my.strata.combos[[i]])
    
    write.csv(x = sum.stats.byX.joint.attr.mort[[1]], 
              file = paste0(output_location, diet_pattern, 
                         "/joint/summarystats_joint_attributable_USmortality_", 
                         "by_", paste(my.strata.combos[[i]], sep = "", collapse = "_"), 
                         "_", year.vec.string, "_", diet_pattern, "_new.csv"), 
              row.names = FALSE)
    
    write.csv(x = sum.stats.byX.joint.attr.mort[[2]], 
              file = paste0(output_location, diet_pattern, 
                         "/joint/joint_attributable_USmortality_draws_", 
                         "by_", paste(my.strata.combos[[i]], sep = "", collapse = "_"), 
                         "_", year.vec.string, "_", diet_pattern, ".csv"), 
              row.names = FALSE)
    
    write.csv(x = sum.stats.byX.joint.attr.mort[[3]], 
              file = paste0(output_location, diet_pattern, 
                         "/joint/RE_joint_PIFs_draws_", 
                         "by_", paste(my.strata.combos[[i]], sep = "", collapse = "_"), 
                         "_", year.vec.string, "_", diet_pattern, ".csv"), 
              row.names = FALSE)
    
    write.csv(x = sum.stats.byX.attr.mort[[4]], 
              file = paste0(output_location, diet_pattern, 
                         "/joint/RE_joint_PIFs_draws_", 
                         "by_", paste(my.strata.combos[[i]], sep = "", collapse = "_"), 
                         "_", year.vec.string, "_", diet_pattern, ".csv"), 
              row.names = FALSE)
    
  }
  
}

# apply
calculate_combos(my.strata.combos = strata.combos.new0, x = 1)
calculate_combos(my.strata.combos = strata.combos.new1, x = 2)
calculate_combos(my.strata.combos = strata.combos.new2, x = 3)
calculate_combos(my.strata.combos = strata.combos.new3, x = 4)
calculate_combos(my.strata.combos = strata.combos.new4, x = 5)
calculate_combos(my.strata.combos = strata.combos.new5, x = 6)

