#only works on specific allmort file created by my code
##use covar="overall" for all

##function of calculating RE PAFs sum stats
get.RE.PAFs.draws<-function(allmort, totalmort) ##for now, only works for age/sex/race
{
  merged<-merge(allmort, totalmort)
  
  allmort.starting.point<-which(names(merged)=="V1")
  allmort.ending.point<-which(names(merged)=="riskfactor")-1
  totalmort.starting.point<-which(names(merged)=="X1")
  totalmort.ending.point<-dim(merged)[2]
  
  RE.PAFs.info<-merged[,which(names(merged) %in% c("age", "female", "race", "outcome", "mean..food.", "se..food", "sd..food.", "riskfactor", "count", "se"))]
  
  RE.PAFs<-cbind(merged[,which(names(merged) %in% c("age", "female", "race", "outcome", "mean..food.", "se..food", "sd..food.", "riskfactor", "count", "se"))], 
                 merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
  RE.PAFs<-RE.PAFs[order(RE.PAFs$riskfactor, RE.PAFs$outcome, RE.PAFs$female, RE.PAFs$age, RE.PAFs$race),]
  
  return(RE.PAFs)
}
  
get.RE.PAFs<-function(allmort, totalmort) ##for now, only works for age/sex/race
{
  merged<-merge(allmort, totalmort)

  allmort.starting.point<-which(names(merged)=="V1")
  allmort.ending.point<-which(names(merged)=="riskfactor")-1
  totalmort.starting.point<-which(names(merged)=="X1")
  totalmort.ending.point<-dim(merged)[2]

  RE.PAFs.info<-merged[,which(names(merged) %in% c("age", "female", "race", "outcome", "mean..food.", "se..food", "sd..food.", "riskfactor", "count", "se"))]

  RE.PAFs<-cbind(merged[,which(names(merged) %in% c("age", "female", "race", "outcome", "mean..food.", "se..food", "sd..food.", "riskfactor", "count", "se"))], 
               merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
  RE.PAFs<-RE.PAFs[order(RE.PAFs$riskfactor, RE.PAFs$outcome, RE.PAFs$female, RE.PAFs$age, RE.PAFs$race),]

  RE.PAFs.sum.stats<-sum.stats.maker.RE.PAFs(RE.PAFs, n.covar)
  return(RE.PAFs.sum.stats)
}

#function to get summary stats##################################
sum.stats.maker<-function(allmort)
{
  allmort.dt<-as.data.table(allmort)
  vars<-paste("V", 1:n.sims, sep="")
  
  
  mort.summary <- allmort.dt[, ":="(means = rowMeans(.SD, na.rm = TRUE),
                                    sd.devs = apply(.SD, 1, sd),
                                    medians = apply(.SD, 1, median),
                                    LB = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                    UB = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))
  ),
  .SDcols = vars
  ]
  
  mort.summary <- mort.summary[, (vars) := NULL]
  
  return(mort.summary)
}

sum.stats.maker.RE.PAFs<-function(RE.PAFs)
{
  RE.PAFs.dt<-as.data.table(RE.PAFs)
  vars<-paste("V", 1:n.sims, sep="")
  
  RE.PAFs.summary <- RE.PAFs.dt[, ":="(RE.PAF.means = rowMeans(.SD, na.rm = TRUE),
                                       RE.PAF.sd.devs = apply(.SD, 1, sd),
                                       RE.PAF.medians = apply(.SD, 1, median),
                                       RE.PAF.LB = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                       RE.PAF.UB = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T)),
                                       RE.PAF.medians.as.percent = apply(.SD, 1, median)*100,
                                       RE.PAF.LB.as.percent = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T))*100,
                                       RE.PAF.UB.as.percent = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))*100
  ),
  .SDcols = vars
  ]
  
  RE.PAFs.summary <- RE.PAFs.summary[, (vars) := NULL]
  
  return(RE.PAFs.summary)
}


#########################################

# debug
# pop=pop.draws
# covar=strata.combos[[i]]

Sum.by.strata<-function(allmort, totalmort, pop, covar)
{
  allmort.dt<-as.data.table(allmort)
  totalmort.dt<-as.data.table(totalmort)
  pop.dt<-as.data.table(pop)
  
  strata<-covar
  
  cols.to.sum.allmort<-paste("V", 1:n.sims, sep="")
  cols.to.sum.totalmort<-paste("X", 1:n.sims, sep="")
  cols.to.sum.pop<-paste("Y", 1:n.sims, sep="")
  

    strata.sims.allmort<-allmort.dt[,lapply(.SD, sum), by=c(strata, "riskfactor"), .SDcols=cols.to.sum.allmort]
    strata.sims.totalmort<-totalmort.dt[,lapply(.SD, sum), 
                                        by=c(strata[strata %in% c("age", "female", "race", "outcome", "disease_type", "disease_type1", "disease_type2", "disease_type3")]), 
                                        .SDcols=cols.to.sum.totalmort]
    strata.sims.pop<-pop.dt[,lapply(.SD, sum), by=c(strata[strata %in% c("age", "female", "race")]), .SDcols=cols.to.sum.pop]
    
    
    #merged<-reduce(list(strata.sims.allmort, strata.sims.totalmort, strata.sims.pop), merge)
    if(length(strata[strata %in% c("age", "female", "race")])>0 & dim(strata.sims.totalmort)[1]>1)
      merged<-merge(merge(strata.sims.allmort, strata.sims.totalmort), strata.sims.pop, by=strata[strata %in% c("age", "female", "race")])
    if(length(strata[strata %in% c("age", "female", "race")])>0 & dim(strata.sims.totalmort)[1]==1)
    {
      merge.a<-cbind(strata.sims.allmort, strata.sims.totalmort)
      merged<-merge(merge.a, strata.sims.pop, by=strata[strata %in% c("age", "female", "race")])
    }
    if(length(strata[strata %in% c("age", "female", "race")])==0 & dim(strata.sims.totalmort)[1]>1)
    {
      merged<-cbind(merge(strata.sims.allmort, strata.sims.totalmort), strata.sims.pop)
    }
    if(length(strata[strata %in% c("age", "female", "race")])==0 & dim(strata.sims.totalmort)[1]==1)
      merged<-cbind(strata.sims.allmort, strata.sims.totalmort, strata.sims.pop)
    
    allmort.starting.point<-which(names(merged)=="V1")
    allmort.ending.point<-which(names(merged)==paste("V", nsim1, sep=""))
    totalmort.starting.point<-which(names(merged)=="X1")
    totalmort.ending.point<-which(names(merged)==paste("X", nsim1, sep=""))
    pop.starting.point<-which(names(merged)=="Y1")
    pop.ending.point<-which(names(merged)==paste("Y", nsim1, sep=""))
    
    #RE.PAFs.info<-merged[,c(unlist(strata.combos[j]), "outcome"), with=FALSE]
    RE.PAFs<-cbind(merged[,c(strata, "riskfactor"), with=FALSE], 
                   merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
    

    setorderv(RE.PAFs, cols=c("riskfactor", strata))
    
    sum.stats.attr.mort.PAF<-sum.stats.maker(strata.sims.allmort)
    RE.PAFs.sum.stats<-sum.stats.maker.RE.PAFs(RE.PAFs)
    
    standardized.mort<-cbind(merged[,c(strata, "riskfactor"), with=FALSE], 
                             merged[,allmort.starting.point:allmort.ending.point]/merged[,pop.starting.point:pop.ending.point]*100000)
    
    setorderv(standardized.mort, cols=c("riskfactor", strata))
    standardized.mort.sum.stats<-sum.stats.maker(standardized.mort)
    names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))]<-
      paste("mort.per.100k.",names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))], sep="")
    
    #merged.results<-merge(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats)
    merged.results<-reduce(list(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats, standardized.mort.sum.stats), merge)
    
    strata.summary.stats<-merged.results
    
    return(list(strata.summary.stats, strata.sims.allmort, strata.sims.totalmort, RE.PAFs, standardized.mort)) 
    

}

Sum.by.overall<-function(allmort, totalmort, pop)
{
  allmort.dt<-as.data.table(allmort)
  totalmort.dt<-as.data.table(totalmort)
  pop.dt<-as.data.table(pop)
  
  cols.to.sum.allmort<-paste("V", 1:n.sims, sep="")
  cols.to.sum.totalmort<-paste("X", 1:n.sims, sep="")
  cols.to.sum.pop<-paste("Y", 1:n.sims, sep="")
  
  strata.sims.allmort<-allmort.dt[,lapply(.SD, sum), by=c("outcome", "disease_type", "riskfactor"), .SDcols=cols.to.sum.allmort]
  strata.sims.totalmort<-totalmort.dt[,lapply(.SD, sum), by=c("outcome", "disease_type"), .SDcols=cols.to.sum.totalmort]
  strata.sims.pop<-pop.dt[,lapply(.SD, sum), .SDcols=cols.to.sum.pop]
  

  #merged<-merge(strata.sims.allmort, strata.sims.totalmort)
  #merged<-reduce(list(strata.sims.allmort, strata.sims.totalmort, strata.sims.pop), merge)
  #merged<-merge(merge(strata.sims.allmort, strata.sims.totalmort), strata.sims.pop, by="")
  merged<-merge(strata.sims.allmort, strata.sims.totalmort)
  #replicate to get n rows of same numbers
  strata.sims.pop<-map_dfr(seq_len(dim(merged)[1]), ~strata.sims.pop)
  merged<-cbind(merged, strata.sims.pop)
  
  allmort.starting.point<-which(names(merged)=="V1")
  allmort.ending.point<-which(names(merged)==paste("V", nsim1, sep=""))
  totalmort.starting.point<-which(names(merged)=="X1")
  totalmort.ending.point<-which(names(merged)==paste("X", nsim1, sep=""))
  pop.starting.point<-which(names(merged)=="Y1")
  pop.ending.point<-which(names(merged)==paste("Y", nsim1, sep=""))
  
  #RE.PAFs.info<-merged[,c(unlist(strata.combos[j]), "outcome"), with=FALSE]
  RE.PAFs<-cbind(merged[,c("outcome", "riskfactor"), with=FALSE], 
                 merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
  
  setorder(RE.PAFs, riskfactor, outcome)
  sum.stats.attr.mort.PAF<-sum.stats.maker(strata.sims.allmort)
  RE.PAFs.sum.stats<-sum.stats.maker.RE.PAFs(RE.PAFs)
  
  standardized.mort<-cbind(merged[,c("outcome", "riskfactor"), with=FALSE], 
                           merged[,allmort.starting.point:allmort.ending.point]/merged[,pop.starting.point:pop.ending.point]*100000)
  setorder(standardized.mort, riskfactor, outcome)
  
  standardized.mort.sum.stats<-sum.stats.maker(standardized.mort)
  names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))]<-
    paste("mort.per.100k.",names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))], sep="")
  
  #merged.results<-merge(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats)
  merged.results<-reduce(list(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats, standardized.mort.sum.stats), merge)
  
  strata.summary.stats<-merged.results
  
  return(list(strata.summary.stats, strata.sims.allmort, strata.sims.totalmort, RE.PAFs, standardized.mort)) 
  
  
}

# debug
# allmort=joint.mort.all.draws[[4]]
# totalmort=totalmort[[4]]
# pop=pop.draws
# covar=strata.combos[[i]]

# allmort=joint.mort.all.draws[[x]]
# totalmort=total.mort[[x]]
# pop=pop.draws
# covar=my.strata.combos[[i]]

Sum.by.strata.joint<-function(allmort, totalmort, pop, covar)
{
  allmort.dt<-as.data.table(allmort)
  totalmort.dt<-as.data.table(totalmort)
  pop.dt<-as.data.table(pop)
  
  strata<-covar
  
  cols.to.sum.allmort<-paste("V", 1:n.sims, sep="")
  cols.to.sum.totalmort<-paste("X", 1:n.sims, sep="")
  cols.to.sum.pop<-paste("Y", 1:n.sims, sep="")
  
  #strata.sims.allmort<-allmort.dt[,lapply(.SD, sum), by=c(strata, "outcome", "outcome_all", "disease_type"), .SDcols=cols.to.sum.allmort]
  #strata.sims.totalmort<-totalmort.dt[,lapply(.SD, sum), by=c(strata, "outcome", "outcome_all", "disease_type"), .SDcols=cols.to.sum.totalmort]
  
  strata.sims.allmort<-allmort.dt[,lapply(.SD, sum), by=c(strata), .SDcols=cols.to.sum.allmort]
  
  # this was the line causing the problem - fixed now
  # 4-11-24 against my better judgement, i'm changing this back to see if it gets fixed
  strata.sims.totalmort<-totalmort.dt[,lapply(.SD, sum), by=c(strata[strata %in% c("age", "female", "race", "outcome", 
                                                                                   # "pathway", 
                                                                                   "disease_type", 
                                                                                   "disease_type1", "disease_type2", "disease_type3"
                                                                                   )]), 
                                      .SDcols=cols.to.sum.totalmort]

  strata.sims.pop<-pop.dt[,lapply(.SD, sum), by=c(strata[strata %in% c("age", "female", "race")]), .SDcols=cols.to.sum.pop]
  
  
  #merged<-merge(strata.sims.allmort, strata.sims.totalmort)
  if(length(strata[strata %in% c("age", "female", "race")])>0 & dim(strata.sims.totalmort)[1]>1)
    merged<-merge(merge(strata.sims.allmort, strata.sims.totalmort), strata.sims.pop, by=strata[strata %in% c("age", "female", "race")])
  if(length(strata[strata %in% c("age", "female", "race")])>0 & dim(strata.sims.totalmort)[1]==1)
  {
    merge.a<-cbind(strata.sims.allmort, strata.sims.totalmort)
    merged<-merge(merge.a, strata.sims.pop, by=strata[strata %in% c("age", "female", "race")])
  }
  if(length(strata[strata %in% c("age", "female", "race")])==0 & dim(strata.sims.totalmort)[1]>1)
  {
    merged<-cbind(merge(strata.sims.allmort, strata.sims.totalmort), strata.sims.pop)
  }
  if(length(strata[strata %in% c("age", "female", "race")])==0 & dim(strata.sims.totalmort)[1]==1)
    merged<-cbind(strata.sims.allmort, strata.sims.totalmort, strata.sims.pop)
  
  allmort.starting.point<-which(names(merged)=="V1")
  allmort.ending.point<-which(names(merged)==paste("V", nsim1, sep=""))
  totalmort.starting.point<-which(names(merged)=="X1")
  totalmort.ending.point<-which(names(merged)==paste("X", nsim1, sep=""))
  pop.starting.point<-which(names(merged)=="Y1")
  pop.ending.point<-which(names(merged)==paste("Y", nsim1, sep=""))
  
  #RE.PAFs.info<-merged[,c(unlist(strata.combos[j]), "outcome"), with=FALSE]
  RE.PAFs<-cbind(merged[,c(strata), with=FALSE], 
                 merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
  
  
  #setorder(RE.PAFs, outcome)
  setorderv(RE.PAFs, cols=c(strata))
  
  
  sum.stats.attr.mort.PAF<-sum.stats.maker(strata.sims.allmort)
  RE.PAFs.sum.stats<-sum.stats.maker.RE.PAFs(RE.PAFs)
  
  standardized.mort<-cbind(merged[,c(strata), with=FALSE], 
                           merged[,allmort.starting.point:allmort.ending.point]/merged[,pop.starting.point:pop.ending.point]*100000)
  
  setorderv(standardized.mort, cols=c(strata))
  standardized.mort.sum.stats<-sum.stats.maker(standardized.mort)
  names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))]<-
    paste("mort.per.100k.",names(standardized.mort.sum.stats)[(length(names(standardized.mort.sum.stats))-4):length(names(standardized.mort.sum.stats))], sep="")
  
  
  #merged.results<-merge(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats)
  merged.results<-reduce(list(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats, standardized.mort.sum.stats), merge, by=strata)
  
  strata.summary.stats<-merged.results
  
  #return(list(strata.summary.stats, strata.sims.allmort, RE.PAFs)) 
  return(list(strata.summary.stats, strata.sims.allmort, strata.sims.totalmort, RE.PAFs, standardized.mort))
}
  
Sum.by.overall.joint<-function(allmort, totalmort)
{
  allmort.dt<-as.data.table(allmort)
  totalmort.dt<-as.data.table(totalmort)
  
  cols.to.sum.allmort<-paste("V", 1:n.sims, sep="")
  cols.to.sum.totalmort<-paste("X", 1:n.sims, sep="")
  
  
  strata.sims.allmort<-allmort.dt[,lapply(.SD, sum), by=c("outcome", "outcome_all", "disease_type"), .SDcols=cols.to.sum.allmort]
  strata.sims.totalmort<-totalmort.dt[,lapply(.SD, sum), by=c("outcome", "outcome_all", "disease_type"), .SDcols=cols.to.sum.totalmort]
  
  merged<-merge(strata.sims.allmort, strata.sims.totalmort)
  
  allmort.starting.point<-which(names(merged)=="V1")
  allmort.ending.point<-which(names(merged)==paste("V", nsim1, sep=""))
  totalmort.starting.point<-which(names(merged)=="X1")
  totalmort.ending.point<-which(names(merged)==paste("X", nsim1, sep=""))
  
  #RE.PAFs.info<-merged[,c(unlist(strata.combos[j]), "outcome"), with=FALSE]
  RE.PAFs<-cbind(merged[,c("outcome"), with=FALSE], 
                 merged[,allmort.starting.point:allmort.ending.point]/merged[,totalmort.starting.point:totalmort.ending.point])
  
  setorder(RE.PAFs, outcome)
  
  sum.stats.attr.mort.PAF<-sum.stats.maker(strata.sims.allmort)
  RE.PAFs.sum.stats<-sum.stats.maker.RE.PAFs(RE.PAFs)
  
  merged.results<-merge(sum.stats.attr.mort.PAF, RE.PAFs.sum.stats, by="outcome", "outcome_all", "disease_type")
  
  strata.summary.stats<-merged.results
  
  return(list(strata.summary.stats, strata.sims.allmort, RE.PAFs)) 
}

###function to create "sensitivity" version of an output file
# create.sensitivity.output<-function(output_file, prop)
# {
#   output<-output_file
#   output[,c("means", "sd.devs", "medians", "LB", "UB", "mort.per.100k.means", "mort.per.100k.sd.devs", "mort.per.100k.medians", "mort.per.100k.LB", "mort.per.100k.UB")]<-
#     output[,c("means", "sd.devs", "medians", "LB", "UB", "mort.per.100k.means", "mort.per.100k.sd.devs", "mort.per.100k.medians", "mort.per.100k.LB", "mort.per.100k.UB")]*prop
#   return(output)
# }
