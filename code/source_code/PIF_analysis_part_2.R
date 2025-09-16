# Title: Part 2 - Run PIF analysis
# Authors: Fred Cudhea & Brooke Bell
# Last Updated: 09-15-2025

############################################################################

### Documentation for R script ###

# This script runs after "part one". After you've generated a distribution of 
# PIFs/attributable mortality/incidence for each subgroup/dietary factor/disease, 
# you run this to generate summary statistics (median, uncertainty intervals) of 
# the output generated in part one.

############################################################################

# First, read in relevant files, including the PIFs/attributable mortality sims 
# file (allmort), observed mortality/incidence sims file (total mort), and 
# observed population sims file (pop.draws).

allmort <- read_csv(file = paste0(output_location, diet_pattern, 
                               "/all.incidence.draws_", year.vec.string, 
                               "_", covar.vec.string, "_", diet_pattern, 
                               ".csv"))


allmort <- allmort[order(allmort$riskfactor, allmort$outcome, 
                         allmort$female, allmort$age, allmort$race),]

totalmort <- read_csv(paste0(file_location, "observed.cancer.mortality.draws.csv"))

pop.draws <- read_csv(paste0(file_location, "observed.pop.draws.csv"))

# For pop.draws, rename sim variables (specifically change "X123", to "Y123", etc...). 
# This will come in handy when we later merge with files that use "X123" naming scheme. 

names(pop.draws) <- gsub("X", "Y", names(pop.draws))

# Recode demographic variables in totalmort (to match with "allmort", the 
# output created in part one).

totalmort <- totalmort[order(totalmort$disease, totalmort$subgroup),]
totalmort$age <- 0
totalmort$age[totalmort$Age_label == "20-34"] <- 25
totalmort$age[totalmort$Age_label == "35-44"] <- 35
totalmort$age[totalmort$Age_label == "45-54"] <- 45
totalmort$age[totalmort$Age_label == "55-64"] <- 55
totalmort$age[totalmort$Age_label == "65-74"] <- 65
totalmort$age[totalmort$Age_label == "75+"] <- 75

totalmort$female[totalmort$Sex == 1] <- 1
totalmort$female[totalmort$Sex == 2] <- 0

totalmort$race[totalmort$Race_label == "NHW"] <- 1
totalmort$race[totalmort$Race_label == "NHB"] <- 2
totalmort$race[totalmort$Race_label == "HIS"] <- 3
totalmort$race[totalmort$Race_label == "OTH"] <- 4

# Recall that, for the purposes of calculating the PIFs, we treat each 
# outcome/pathway as a distinct outcome (direct effect, mediated by BMI, and 
# mediated by SBP). But for purposes of collapsing outcomes (e.g, calculating 
# joint PIFs), it's better to have the outcomes and pathway separated into two 
# different variables. Below is code do that. 

totalmort[!grepl("_", totalmort$disease, fixed = TRUE),]$disease <- 
  paste(totalmort[!grepl("_", totalmort$disease, fixed = TRUE),]$disease, "_direct", sep = "")

totalmort <-totalmort %>% separate(col = disease, into = c("outcome", "pathway"), sep = "_")
totalmort <- totalmort[totalmort$pathway  == "direct",]
totalmort <- subset(totalmort, outcome %in% diseases.vec, select = -c(pathway))

# Recall, we also have three broader outcome categories of interest (denoted as 
# disease type 1, diseaes type 2, disease type 3), and have mapped each outcome 
# to each category (for all three) in run_models.r, and that mapping is stored 
# in disease.table. Code below uses that mapping to create variables for the 
# three disease type categories.

totalmort$disease_type <- 
  disease.type.vec[match(totalmort$outcome[totalmort$outcome %in% diseases.vec], 
                         disease.table)]

totalmort$disease_type1 <- 
  disease.type1.vec[match(totalmort$outcome[totalmort$outcome %in% diseases.vec], 
                          disease.table)]

totalmort$disease_type2 <- 
  disease.type2.vec[match(totalmort$outcome[totalmort$outcome %in% diseases.vec], 
                          disease.table)]

totalmort$disease_type3 <- 
  disease.type3.vec[match(totalmort$outcome[totalmort$outcome %in% diseases.vec], 
                          disease.table)]

# Get rid of unused demo vars (we'll be using the new recoded ones we created earlier).

totalmort <- 
  totalmort[, -which(names(totalmort) %in% c("X", "disease", "Sex", 
                                            "age_gp", "race_gp", "agecat"))]

# Read in PIF sims file.

pif <- read_csv(file = paste0(output_location, diet_pattern, "/all.PIFs_", 
                            year.vec.string, "_", covar.vec.string, "_", 
                            diet_pattern, ".csv"))

# Read in script that defines functions for summing up by demographic variables.

source(paste0(code_location, "functions/sum.by.strata.redone.r"))

# Get summary stats of attributable mortality/incidence for each 
# subgroup/dietary-factor+pathway/outcome. Rename certain variables of 
# resulting variables for clarity.

sum.stats.attr.mort.PIF <- sum.stats.maker(allmort = allmort)

sum.stats.PIF <- sum.stats.maker(allmort = pif)

names(sum.stats.attr.mort.PIF)[(1 + length(covar.vec)):(6 + length(covar.vec))] <- 
  c("mean (food)", "se (food)", "sd (food)", "mean (food, counterfactual)", 
    "se (food, counterfactual)", "sd (food, counterfactual)")

names(sum.stats.PIF)[(1 + length(covar.vec)):(6 + length(covar.vec))] <- 
  c("mean (food)", "se (food)", "sd (food)", "mean (food, counterfactual)", 
    "se (food, counterfactual)", "sd (food, counterfactual)")

# Read in original exposure file, strictly for the purpose of calculating "l", 
# the number of subgroups. Note this is only valid assuming you subgroups are 
# defined by age, sex and race.

# Read CSV file of exposure
mort.dat <- read_csv(file = paste0(file_location, "nhanes1518_agesexrace_merged_", 
                                 version.date ,".csv"))

mort.dat <- mort.dat[order(mort.dat$diet),]

if(identical(covar.vec, c("Age", "Sex", "Race"))) {
  
  l <- length(table(mort.dat$agecat)) *
    length(table(mort.dat$female)) *
    length(table(mort.dat$race))
  
}

# Now we will calculate "reverse engineered PIFs", meaning we are going to 
# back-calculte PIFs from taking the attibutable mortality/incidence and dividing 
# my total number of deaths/incidence. This is how we calculate PIFs for broader 
# subgroups (results by age, results by sex, results by disease type, overall 
# results for the population, etc... as opposed results by each 
# age/sex/race/outcome/etc..). 

# First, we merge "allmort", "totalmort" and "pop.draws" into one file, and 
# identify the column index where part starts in this new dataframe.

merged <- reduce(list(allmort, totalmort, pop.draws), merge)

allmort.starting.point <- which(names(merged) == "V1")
allmort.ending.point <- which(names(merged) == paste("V", nsim1, sep = ""))
totalmort.starting.point <- which(names(merged) == "X1")
totalmort.ending.point <- which(names(merged) == paste("X", nsim1, sep = ""))
pop.starting.point <- which(names(merged) == "Y1")
pop.ending.point <- which(names(merged) == paste("Y", nsim1, sep = ""))

# Create a new matrix with "reverse-engineered" PIFs calculated by dividing 
# "allmort" values (attributable mort) by "totalmort" values (observed mort). 
# Order it (by risk factor, disease type, outcome, sex, age, race). If any NAs 
# (from observed mort values being zero for whatever reason), change to 0. Save. 

RE.PIFs <- cbind(merged[, which(names(merged) %in% c("age", 
                                                    "female", 
                                                    "race", 
                                                    "outcome", 
                                                    "pathway", 
                                                    "disease_type", 
                                                    "disease_type1", 
                                                    "disease_type2", 
                                                    "disease_type3",
                                                    "riskfactor"))], 
                 merged[, allmort.starting.point:allmort.ending.point] / 
                   merged[, totalmort.starting.point:totalmort.ending.point])

RE.PIFs <- RE.PIFs[order(RE.PIFs$riskfactor, 
                         RE.PIFs$disease_type, 
                         RE.PIFs$outcome, 
                         RE.PIFs$female, 
                         RE.PIFs$age, 
                         RE.PIFs$race),]

RE.PIFs[is.na(RE.PIFs)] <- 0

write_csv(x = RE.PIFs, 
          file = paste0(output_location, diet_pattern, "/all.mort.draws_RE.PIFs_", 
                     year.vec.string, "_", covar.vec.string, "_", diet_pattern, 
                     ".csv"))

# Next, we use the same idea to get "standardized mortality", as in, the number 
# of people saved (if at the counterfactual) per 100,000 people. So now, we take 
# attributable mort, divide by total population, then multiply by 100,000.

standardized.mort.info <- merged[,which(names(merged) %in% c("age", 
                                                             "female", 
                                                             "race", 
                                                             "outcome", 
                                                             "pathway", 
                                                             "disease_type", 
                                                             "disease_type1", 
                                                             "disease_type2", 
                                                             "disease_type3",
                                                             "riskfactor"))]

standardized.mort <- 
  cbind(standardized.mort.info, 
        merged[, allmort.starting.point:allmort.ending.point] /
          merged[, pop.starting.point:pop.ending.point] * 100000)

standardized.mort <- 
  standardized.mort[order(standardized.mort$riskfactor, 
                          standardized.mort$outcome, 
                          standardized.mort$disease_type, 
                          standardized.mort$female, 
                          standardized.mort$age, 
                          standardized.mort$race),]

standardized.mort[is.na(standardized.mort)] <- 0

write_csv(x = standardized.mort, 
          file = paste0(output_location, diet_pattern, "/standardized.mort_", 
                     year.vec.string, "_", covar.vec.string, "_", 
                     diet_pattern, ".csv"))

# Now we have a distribution of "reverse-engineered" PIFs and events prevented 
# per 100,000 for all subgroups. Now we take those and get summary stats (and 
# edit var names for standardized.mort.sum.stats file for clarity).

RE.PIFs.sum.stats <- sum.stats.maker.RE.PAFs(RE.PIFs)

standardized.mort.sum.stats <- sum.stats.maker(standardized.mort)

names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats)) - 4):
                                     length(names(standardized.mort.sum.stats))] <- 
  paste0("mort.per.100k.", 
         names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats)) - 4):
                                              length(names(standardized.mort.sum.stats))])

# Merge the summary stats files for attributable mort/incidence, 
# "reverse-engineered" PIFs, and standardized attributable mort/incidence, and 
# save it. Also, save the summary stats file for the actual PIFs 
# (not "reverse-engineered") we made earlier.

sum.stats.attr.mort.merged <- 
  reduce(list(sum.stats.attr.mort.PIF, 
              RE.PIFs.sum.stats, 
              standardized.mort.sum.stats), 
         merge)

write.csv(x = sum.stats.attr.mort.merged, 
          file = paste0(output_location, diet_pattern, 
                     "/summarystats_attributable_USmortality_", 
                     covar.vec.string, "_",year.vec.string, "_", 
                     diet_pattern, ".csv"), row.names = FALSE)

write.csv(x = sum.stats.PIF, 
          file = paste0(output_location, diet_pattern, "/summarystats_PIFs_", 
                     covar.vec.string, "_",year.vec.string, "_", 
                     diet_pattern, ".csv"), row.names = FALSE)

# Ok, now that we've done that. Let's get to making output for the manuscript. 
# That is, summary statistics for broader subgroups (results by age, results by 
# sex, results by disease type, etc... as opposed results by each age/sex/race/outcome/etc..).

# We will create results for all possible combinations that could be of interest 
# here but you will choose what to report based on what your research focus is.

# First we put all strata of interest into a vector called "strata". Then we 
# loop through the length of that vector. For each iteration i of that loop, we 
# extract all combinations of that vector with i element using the combn function, 
# convert the resulting matrix to a list, and append those list elements to 
# "strata.combos" list which will have an element for each possible strata combination.

strata.combos <- list()
strata <- c("age", "female", "race", "outcome", "pathway", 
            "disease_type", "disease_type1", "disease_type2", "disease_type3")

# order should match with covar.vec
for(i in 1:length(strata)) {
  
  x <- combn(strata, i)
  combos.for.i <- split(x, rep(1:ncol(x), each = nrow(x))) 
  strata.combos <- c(strata.combos, combos.for.i)
  
}

# Not all combinations are of interest. We want to, at the very least, have 
# outcome pathway or a disease type as a strata (Results are event specific 
# after all. We aren't reporting all preventable events by age when the types 
# of events are not even consistent by outcome). Also, not all the strata are 
# mutually exclusive. It's redundant to stratify by disease type 1 if you're 
# already stratifying by outcome so we don't need any strata combos that has 
# both outcome and disease type 1, for example. We want to filter out unnecessary 
# strata combos, which is what the code below does (we loop through the list of 
# strata combos and identify the index the strata combos that qualify as of interest).
                                                                                                                         
# An astute reader might note that stratifying just by "pathway" and not 
# outcome/disease type will mix up incidence (from cancer) and mortality 
# (from CMDs). Good point! You got me! (Generally, I recommend not mixing 
# cancer mortality and CMD incidence when reporting results)

strata.combos.of.interest.boolean <- c()

for(i in 1:length(strata.combos)) {
  
  disease_overlap_count <- 
    sum("outcome" %in% strata.combos[[i]], 
        "disease_type" %in% strata.combos[[i]], 
        "disease_type1" %in% strata.combos[[i]], 
        "disease_type2" %in% strata.combos[[i]], 
        "disease_type3" %in% strata.combos[[i]])
  
  strata.combos.of.interest.boolean[i] <- 
    any(strata.combos[[i]] %in% c("outcome", 
                                  "pathway", 
                                  "disease_type", 
                                  "disease_type1", 
                                  "disease_type2", 
                                  "disease_type3")) & 
    (disease_overlap_count < 2)
  
}

strata.combos.of.interest.index <- which(strata.combos.of.interest.boolean)

# check
print(strata.combos.of.interest.index)

# Having identifying the index for each strata combo of interest, we can loop 
# through them and apply the Sum.by.strata command which takes as inputs the 
# attributable mort sims, total mort sims, total pop sims, and the stratas of 
# interest, and uses that to calculate 4 output files  of interest: 
#   
# 1: summary stats of attributable mort/incidence  +  "reverse engineered" PIFs 
#  +  standardized preventable mort/incidence by specified strata, 
# 
# 2: entire simulation of attributable mort/incidence by specified strata, 
# 
# 3: entire simulation of observed deaths by specified strata, 
# 
# 4: entire simulation of "reverse engineered" PIFs by specified strata.

for(i in strata.combos.of.interest.index) {
  
  print(i)
  print(strata.combos[[i]])
  print("")
  
  sum.stats.byX.attr.mort <- 
    Sum.by.strata(allmort = allmort, 
                  totalmort = totalmort, 
                  pop = pop.draws, 
                  covar = strata.combos[[i]])
  
  write.csv(x = sum.stats.byX.attr.mort[[1]], 
            file = paste0(output_location, diet_pattern, 
                       "/summarystats_attributable_USmortality_", 
                       "by_", paste(strata.combos[[i]], sep = "", collapse = "_"), 
                       "_", year.vec.string, "_", diet_pattern,".csv"), 
            row.names = FALSE)
  
  write.csv(x = sum.stats.byX.attr.mort[[2]], 
            file = paste0(output_location, diet_pattern, 
                       "/attributable_USmortality_draws_", 
                       "by_", paste(strata.combos[[i]], sep = "", collapse = "_"), 
                       "_", year.vec.string, "_", diet_pattern, ".csv"), 
            row.names = FALSE)
  
  write.csv(x = sum.stats.byX.attr.mort[[3]], 
            file = paste0(output_location, diet_pattern, 
                       "/total_USmortality_draws_", 
                       "by_", paste(strata.combos[[i]], sep = "", collapse = "_"), 
                       "_", year.vec.string, "_", diet_pattern, ".csv"), 
            row.names = FALSE)
  
  write.csv(x = sum.stats.byX.attr.mort[[4]], 
            file = paste0(output_location, diet_pattern, 
                       "/RE_PIFs_draws_", 
                       "by_", paste(strata.combos[[i]], sep = "", collapse = "_"), 
                       "_", year.vec.string, "_", diet_pattern, ".csv"), 
            row.names = FALSE)
  
}


