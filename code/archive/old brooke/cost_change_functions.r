#cost_change_functions.r

simulate.impact<-function(current.mean, current.se, 
                          counterfactual.mean, counterfactual.se=0,
                          impact.factor.names, impact.factor.means,
                          impact.factor.ses, RRunit, 
                          substitution.impact.factor.means,
                          substitution.impact.factor.ses, substitution_unit, 
                          current_inedible_p, current_inedible_p_se,
                          current_foodwaste_p, current_foodwaste_p_se,
                          counterfactual_inedible_p, counterfactual_inedible_p_se,
                          counterfactual_foodwaste_p, counterfactual_foodwaste_p_se,
                          n.sims,
                          population.sims ##<- must be a vector, so make sure to convert to vector (from single row datafame) if necessary
)
{
  # current.mean.sims<-rnorm(n=n.sims, 
  #                          mean=current.mean, 
  #                          sd=current.se)
  # counterfactual.mean.sims<-rnorm(n=n.sims, 
  #                                 mean=counterfactual.mean, 
  #                                 sd=counterfactual.se)
  
  current.mean.sims<-rnorm(n=n.sims, 
                            mean=current.mean, 
                            sd=current.se)
  current.mean.sims.matrix<-t(matrix(rep(current.mean.sims, times=length(impact.factor.names)),
                           ncol=length(impact.factor.names)))
  counterfactual.mean.sims<-rnorm(n=n.sims, 
                                  mean=counterfactual.mean, 
                                  sd=counterfactual.se)
  counterfactual.mean.sims.matrix<-t(matrix(rep(counterfactual.mean.sims, times=length(impact.factor.names)),
                                     ncol=length(impact.factor.names)))
  
  #current.correction<-1/((1-rnorm(n.sims, mean=current_inedible_p, sd=current_inedible_p_se))*(1-rnorm(n.sims, mean=current_foodwaste_p, sd=current_foodwaste_p_se)))
  #counterfactual.correction<-1/((1-rnorm(n.sims, mean=counterfactual_inedible_p, sd=counterfactual_inedible_p_se))*(1-rnorm(n.sims, mean=counterfactual_foodwaste_p, sd=counterfactual_foodwaste_p_se)))
  
  current.correction<-1/((1-rnorm(length(impact.factor.names), mean=current_inedible_p, sd=current_inedible_p_se))*(1-rnorm(length(impact.factor.names), mean=current_foodwaste_p, sd=current_foodwaste_p_se)))
  counterfactual.correction<-1/((1-rnorm(length(impact.factor.names), mean=counterfactual_inedible_p, sd=counterfactual_inedible_p_se))*(1-rnorm(length(impact.factor.names), mean=counterfactual_foodwaste_p, sd=counterfactual_foodwaste_p_se)))
  
  # UNCOMMENT LATER 
  # ##this is a correction factor to get total edible. "current.correction" above is for getting total produced (edible (consumed and wasted) + inedible)
   current.correction.for.waste<-1/(1-rnorm(length(impact.factor.names), mean=current_foodwaste_p, sd=current_foodwaste_p_se))
   counterfactual.correction.for.waste<-1/(1-rnorm(length(impact.factor.names), mean=counterfactual_foodwaste_p, sd=counterfactual_foodwaste_p_se))
   
  
  #current.mean.sims.after.correction<-current.mean.sims*current.correction
  #counterfactual.mean.sims.after.correction<-counterfactual.mean.sims*counterfactual.correction
  
  current.mean.sims.matrix.after.correction<-current.mean.sims.matrix*current.correction
  counterfactual.mean.sims.matrix.after.correction<-counterfactual.mean.sims.matrix*counterfactual.correction
  
  current.mean.sims.matrix.unconsumed<-current.mean.sims.matrix.after.correction-current.mean.sims.matrix
  counterfactual.mean.sims.matrix.unconsumed<-counterfactual.mean.sims.matrix.after.correction-counterfactual.mean.sims.matrix
  
  ## UNCOMMENT LATER 
   current.mean.sims.matrix.edible<-current.mean.sims.matrix*current.correction.for.waste
   counterfactual.mean.sims.matrix.edible<-counterfactual.mean.sims.matrix*counterfactual.correction.for.waste
  
   current.mean.sims.matrix.inedible <- current.mean.sims.matrix.after.correction - current.mean.sims.matrix.edible
   counterfactual.mean.sims.matrix.inedible <- counterfactual.mean.sims.matrix.after.correction - counterfactual.mean.sims.matrix.edible
  
   current.mean.sims.matrix.wasted <- current.mean.sims.matrix.edible - current.mean.sims.matrix
   counterfactual.mean.sims.matrix.wasted <- counterfactual.mean.sims.matrix.edible - counterfactual.mean.sims.matrix
  
  #delta.sims<-counterfactual.mean.sims-current.mean.sims
  #delta.sims<-counterfactual.mean.sims.after.correction-current.mean.sims.after.correction

  delta.sims.matrix<-counterfactual.mean.sims.matrix.after.correction-current.mean.sims.matrix.after.correction

  delta.sims.matrix.consumed<-counterfactual.mean.sims.matrix-current.mean.sims.matrix
  delta.sims.matrix.unconsumed<-counterfactual.mean.sims.matrix.unconsumed-current.mean.sims.matrix.unconsumed
  delta.sims.matrix.total<-delta.sims.matrix.consumed+delta.sims.matrix.unconsumed
  
  # UNCOMMENT LATER 
   delta.sims.matrix.inedible<-counterfactual.mean.sims.matrix.inedible-current.mean.sims.matrix.inedible
   delta.sims.matrix.wasted<-counterfactual.mean.sims.matrix.wasted-current.mean.sims.matrix.wasted
  
  

  #population.sims<-rnorm(n=n.sims, mean=population, sd=population.se)
  
  impact.factor.sims<-matrix(rnorm(n=n.sims*length(impact.factor.names), 
                                   mean=as.numeric(impact.factor.means), 
                                   sd=as.numeric(impact.factor.ses)), 
                             nrow=length(impact.factor.names))
  
  substitution.impact.factor.sims<-matrix(rnorm(n=n.sims*length(impact.factor.names), 
                                                mean=as.numeric(substitution.impact.factor.means), 
                                                sd=as.numeric(substitution.impact.factor.ses)), 
                                          nrow=length(impact.factor.names))
  
  #impact.sims <- t(t(impact.factor.sims)*delta.sims*t(population.sims))/as.numeric(RRunit)
  
  #impact.sims <- impact.factor.sims*delta.sims.matrix*population.sims/as.numeric(RRunit)
  
  impact.sims.consumed <- impact.factor.sims*delta.sims.matrix.consumed*population.sims/as.numeric(RRunit)
  impact.sims.unconsumed <- impact.factor.sims*delta.sims.matrix.unconsumed*population.sims/as.numeric(RRunit)
  impact.sims.total<-impact.sims.consumed+impact.sims.unconsumed
  
  # UNCOMMENT LATER 
   impact.sims.inedible <- impact.factor.sims*delta.sims.matrix.inedible*population.sims/as.numeric(RRunit)
   impact.sims.wasted <- impact.factor.sims*delta.sims.matrix.wasted*population.sims/as.numeric(RRunit)
  
  
  #rownames(impact.sims)<-impact.factor.names
  rownames(impact.sims.consumed)<-impact.factor.names
  rownames(impact.sims.unconsumed)<-impact.factor.names
  rownames(impact.sims.total)<-impact.factor.names
  
  # UNCOMMENT LATER 
   rownames(impact.sims.inedible)<-impact.factor.names
   rownames(impact.sims.wasted)<-impact.factor.names
  
  #substitution.sims<-t(t(substitution.impact.factor.sims)*(-delta.sims)*t(population.sims))/substitution_unit
  
  #substitution.sims<-substitution.impact.factor.sims*-delta.sims.matrix*population.sims/substitution_unit
  
  substitution.sims.consumed<-substitution.impact.factor.sims*-delta.sims.matrix.consumed*population.sims/substitution_unit
  substitution.sims.unconsumed<-substitution.impact.factor.sims*-delta.sims.matrix.unconsumed*population.sims/substitution_unit
  substitution.sims.total<-substitution.sims.consumed+substitution.sims.unconsumed
  
  # UNCOMMENT LATER 
   substitution.sims.inedible<-substitution.impact.factor.sims*-delta.sims.matrix.inedible*population.sims/substitution_unit
   substitution.sims.wasted<-substitution.impact.factor.sims*-delta.sims.matrix.wasted*population.sims/substitution_unit
  
  #rownames(substitution.sims)<-impact.factor.names
  rownames(substitution.sims.consumed)<-impact.factor.names
  rownames(substitution.sims.unconsumed)<-impact.factor.names
  rownames(substitution.sims.total)<-impact.factor.names
  
  # UNCOMMENT LATER 
   rownames(substitution.sims.inedible)<-impact.factor.names
   rownames(substitution.sims.wasted)<-impact.factor.names

  
  #combined.impact.sims<-impact.sims+substitution.sims
  
  combined.impact.sims.consumed<-impact.sims.consumed+substitution.sims.consumed
  combined.impact.sims.unconsumed<-impact.sims.unconsumed+substitution.sims.unconsumed
  combined.impact.sims.total<-combined.impact.sims.consumed+combined.impact.sims.unconsumed
  
  # UNCOMMENT LATER 
   combined.impact.sims.inedible<-impact.sims.consumed+substitution.sims.inedible
   combined.impact.sims.wasted<-impact.sims.unconsumed+substitution.sims.wasted
  
  ###total impact/cost, not just difference between current and counterfactual
  current.impact.sims<-impact.factor.sims*current.mean.sims.matrix.after.correction*population.sims/as.numeric(RRunit)
  counterfactual.impact.sims<-impact.factor.sims*counterfactual.mean.sims.matrix.after.correction*population.sims/as.numeric(RRunit)
  rownames(current.impact.sims)<-impact.factor.names
  rownames(counterfactual.impact.sims)<-impact.factor.names
  
  ##consumed
  current.impact.sims.consumed<-impact.factor.sims*current.mean.sims.matrix*population.sims/as.numeric(RRunit)
  counterfactual.impact.sims.consumed<-impact.factor.sims*counterfactual.mean.sims.matrix*population.sims/as.numeric(RRunit)
  rownames(current.impact.sims.consumed)<-impact.factor.names
  rownames(counterfactual.impact.sims.consumed)<-impact.factor.names
  
  ##unconsumed
  current.impact.sims.unconsumed<-impact.factor.sims*current.mean.sims.matrix.unconsumed*population.sims/as.numeric(RRunit)
  counterfactual.impact.sims.unconsumed<-impact.factor.sims*counterfactual.mean.sims.matrix.unconsumed*population.sims/as.numeric(RRunit)
  rownames(current.impact.sims.unconsumed)<-impact.factor.names
  rownames(counterfactual.impact.sims.unconsumed)<-impact.factor.names
  
  # UNCOMMENT LATER 
  ## inedible
  current.impact.sims.inedible<-impact.factor.sims*current.mean.sims.matrix.inedible*population.sims/as.numeric(RRunit)
  counterfactual.impact.sims.inedible<-impact.factor.sims*counterfactual.mean.sims.matrix.inedible*population.sims/as.numeric(RRunit)
  rownames(current.impact.sims.inedible)<-impact.factor.names
  rownames(counterfactual.impact.sims.inedible)<-impact.factor.names
  ## wasted
  current.impact.sims.wasted<-impact.factor.sims*current.mean.sims.matrix.wasted*population.sims/as.numeric(RRunit)
  counterfactual.impact.sims.wasted<-impact.factor.sims*counterfactual.mean.sims.matrix.wasted*population.sims/as.numeric(RRunit)
  rownames(current.impact.sims.wasted)<-impact.factor.names
  rownames(counterfactual.impact.sims.wasted)<-impact.factor.names

  # population.sims.output<-as.matrix(population.sims)
  # if(dim(population.sims.output)[1]>1)
  # {
  #   population.sims.output<-population.sims.output[rep(1, times=length(impact.factor.names)),] ##repeat the one row for each impact factor, to match with other sims output
  #   rownames(population.sims.output)<-impact.factor.names 
  # }
  # if(dim(population.sims.output)[1]==1)
  # {
  #   rownames(population.sims.output)<-impact.factor.names 
  # }
  # 
  # 
  # delta.sims.output<-as.matrix(delta.sims.matrix*population.sims)
  # if(dim(delta.sims.output)[1]>1)
  # {
  #   delta.sims.output<-delta.sims.output[rep(1, times=length(impact.factor.names)),] ##repeat the one row for each impact factor, to match with other sims output
  #   rownames(delta.sims.output)<-impact.factor.names
  # }
  # if(dim(delta.sims.output)[1]==1)
  # {
  #   rownames(delta.sims.output)<-impact.factor.names    
  # }
  
  ##Now population.sims.output and detla.sims.output are already matrices, don't need above code to replicate first row. Can just use as is
  ##later, can decide if we want to keep pure deltas, or the modified deltas we use now
  
  population.sims.output<-as.matrix(population.sims)
  rownames(population.sims.output)<-impact.factor.names 
  
  #delta.sims.output<-as.matrix(delta.sims.matrix*population.sims)
  #rownames(delta.sims.output)<-impact.factor.names
  
  delta.sims.output.consumed<-as.matrix(delta.sims.matrix.consumed*population.sims)
  rownames(delta.sims.output.consumed)<-impact.factor.names
  delta.sims.output.unconsumed<-as.matrix(delta.sims.matrix.unconsumed*population.sims)
  rownames(delta.sims.output.unconsumed)<-impact.factor.names
  delta.sims.output.total<-as.matrix(delta.sims.matrix.total*population.sims)
  rownames(delta.sims.output.total)<-impact.factor.names
  
  # UNCOMMENT LATER 
  delta.sims.output.inedible<-as.matrix(delta.sims.matrix.inedible*population.sims)
  rownames(delta.sims.output.inedible)<-impact.factor.names
  delta.sims.output.wasted<-as.matrix(delta.sims.matrix.wasted*population.sims)
  rownames(delta.sims.output.wasted)<-impact.factor.names
  
  # impact.list<-list("impact.sims" = impact.sims,
  #                   "substitution.sims" = substitution.sims,
  #                   "combined.impact.sims" = combined.impact.sims,
  #                   "population.sims" = population.sims.output,
  #                   "delta.sims" = delta.sims.output)
  
  impact.list<-list("impact.sims.total" = impact.sims.total,
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
                    
                    #UNCOMMENT LATER,
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

get.percapita.sims<-function(impact.sims, population.sims, pop.sims.names, new.pop.sims.names, sims.names)
{
  names(population.sims)[names(population.sims) %in% pop.sims.names]<-new.pop.sims.names ##replace X's with "pop"s
  impact.sims.plus.pop<-merge(impact.sims, population.sims)
  
  impact.sims.percapita<-impact.sims.plus.pop
  impact.sims.percapita[,sims.names]<-impact.sims.percapita[,sims.names]/impact.sims.percapita[,new.pop.sims.names]
  impact.sims.percapita<-impact.sims.percapita[,!(names(impact.sims.percapita) %in% new.pop.sims.names)]
  
  return(impact.sims.percapita)
}

#let's make a functino out of this, since we're gonna repeat this whole thing for per capita
get.summary.stats<-function(impact.sims, substitution.impact.sims, combined.impact.sims, delta.sims, 
                            impact.sims.consumed, substitution.impact.sims.consumed, combined.impact.sims.consumed,
                            impact.sims.unconsumed, substitution.impact.sims.unconsumed, combined.impact.sims.unconsumed, 
                            current.intake.impact.sims, current.intake.impact.sims.consumed, current.intake.impact.sims.unconsumed, 
                            CF.intake.impact.sims, CF.intake.impact.sims.consumed, CF.intake.impact.sims.unconsumed,
                            impact.sims.inedible, substitution.impact.sims.inedible, combined.impact.sims.inedible,
                            impact.sims.wasted, substitution.impact.sims.wasted, combined.impact.sims.wasted, 
                            current.intake.impact.sims.inedible, current.intake.impact.sims.wasted,
                            CF.intake.impact.sims.inedible, CF.intake.impact.sims.wasted,
                            vars, convert.to.kcal.units="convert.to.kcal.units")
{
  impact.sims.dt<-as.data.table(impact.sims)
  impact.sims.summary <- impact.sims.dt[, ":="(mean_impact = rowMeans(.SD, na.rm = TRUE),
                                               sd_impact = apply(.SD, 1, sd),
                                               median_impact = apply(.SD, 1, median),
                                               lowerci_95_impact = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                               upperci_95_impact = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                        .SDcols = vars
  ]
  impact.sims.summary <- impact.sims.summary[, (vars) := NULL]
  
  substitution.impact.sims.dt<-as.data.table(substitution.impact.sims)
  substitution.impact.sims.summary <- substitution.impact.sims.dt[, ":="(mean_substitution_impact = rowMeans(.SD, na.rm = TRUE),
                                                                         sd_substitution_impact = apply(.SD, 1, sd),
                                                                         median_substitution_impact = apply(.SD, 1, median),
                                                                         lowerci_95_substitution_impact = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                         upperci_95_substitution_impact = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                  #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                  #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                  .SDcols = vars
  ]
  substitution.impact.sims.summary <- substitution.impact.sims.summary[, (vars) := NULL]
  
  combined.impact.sims.dt<-as.data.table(combined.impact.sims)
  combined.impact.sims.summary <- combined.impact.sims.dt[, ":="(mean_combined_impact = rowMeans(.SD, na.rm = TRUE),
                                                                 sd_combined_impact = apply(.SD, 1, sd),
                                                                 median_combined_impact = apply(.SD, 1, median),
                                                                 lowerci_95_combined_impact = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                 upperci_95_combined_impact = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                          #lowerci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                          #upperci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                          .SDcols = vars
  ]
  combined.impact.sims.summary <- combined.impact.sims.summary[, (vars) := NULL]
  
  delta.sims.dt<-as.data.table(delta.sims)
  delta.sims.summary <- delta.sims.dt[, ":="(mean_delta = rowMeans(.SD, na.rm = TRUE),
                                             sd_delta = apply(.SD, 1, sd),
                                             median_delta = apply(.SD, 1, median),
                                             lowerci_95_delta = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                             upperci_95_delta = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                      #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                      #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                      .SDcols = vars
  ]
  delta.sims.summary <- delta.sims.summary[, (vars) := NULL]
  
  ######_consumed
  impact.sims.consumed.dt<-as.data.table(impact.sims.consumed)
  impact.sims.consumed.summary <- impact.sims.consumed.dt[, ":="(mean_impact_consumed = rowMeans(.SD, na.rm = TRUE),
                                               sd_impact_consumed = apply(.SD, 1, sd),
                                               median_impact_consumed = apply(.SD, 1, median),
                                               lowerci_95_impact_consumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                               upperci_95_impact_consumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                        .SDcols = vars
  ]
  impact.sims.consumed.summary <- impact.sims.consumed.summary[, (vars) := NULL]
  
  substitution.impact.sims.consumed.dt<-as.data.table(substitution.impact.sims.consumed)
  substitution.impact.sims.consumed.summary <- substitution.impact.sims.consumed.dt[, ":="(mean_substitution_impact_consumed = rowMeans(.SD, na.rm = TRUE),
                                                                         sd_substitution_impact_consumed = apply(.SD, 1, sd),
                                                                         median_substitution_impact_consumed = apply(.SD, 1, median),
                                                                         lowerci_95_substitution_impact_consumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                         upperci_95_substitution_impact_consumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                  #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                  #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                  .SDcols = vars
  ]
  substitution.impact.sims.consumed.summary <- substitution.impact.sims.consumed.summary[, (vars) := NULL]
  
  combined.impact.sims.consumed.dt<-as.data.table(combined.impact.sims.consumed)
  combined.impact.sims.consumed.summary <- combined.impact.sims.consumed.dt[, ":="(mean_combined_impact_consumed = rowMeans(.SD, na.rm = TRUE),
                                                                 sd_combined_impact_consumed = apply(.SD, 1, sd),
                                                                 median_combined_impact_consumed = apply(.SD, 1, median),
                                                                 lowerci_95_combined_impact_consumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                 upperci_95_combined_impact_consumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                          #lowerci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                          #upperci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                          .SDcols = vars
  ]
  combined.impact.sims.consumed.summary <- combined.impact.sims.consumed.summary[, (vars) := NULL]
  
  ########unconsumed
  impact.sims.unconsumed.dt<-as.data.table(impact.sims.unconsumed)
  impact.sims.unconsumed.summary <- impact.sims.unconsumed.dt[, ":="(mean_impact_unconsumed = rowMeans(.SD, na.rm = TRUE),
                                                                 sd_impact_unconsumed = apply(.SD, 1, sd),
                                                                 median_impact_unconsumed = apply(.SD, 1, median),
                                                                 lowerci_95_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                 upperci_95_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                          #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                          #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                          .SDcols = vars
  ]
  impact.sims.unconsumed.summary <- impact.sims.unconsumed.summary[, (vars) := NULL]
  
  substitution.impact.sims.unconsumed.dt<-as.data.table(substitution.impact.sims.unconsumed)
  substitution.impact.sims.unconsumed.summary <- substitution.impact.sims.unconsumed.dt[, ":="(mean_substitution_impact_unconsumed = rowMeans(.SD, na.rm = TRUE),
                                                                                           sd_substitution_impact_unconsumed = apply(.SD, 1, sd),
                                                                                           median_substitution_impact_unconsumed = apply(.SD, 1, median),
                                                                                           lowerci_95_substitution_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                           upperci_95_substitution_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                    #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                    #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                    .SDcols = vars
  ]
  substitution.impact.sims.unconsumed.summary <- substitution.impact.sims.unconsumed.summary[, (vars) := NULL]
  
  combined.impact.sims.unconsumed.dt<-as.data.table(combined.impact.sims.unconsumed)
  combined.impact.sims.unconsumed.summary <- combined.impact.sims.unconsumed.dt[, ":="(mean_combined_impact_unconsumed = rowMeans(.SD, na.rm = TRUE),
                                                                                   sd_combined_impact_unconsumed = apply(.SD, 1, sd),
                                                                                   median_combined_impact_unconsumed = apply(.SD, 1, median),
                                                                                   lowerci_95_combined_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                   upperci_95_combined_impact_unconsumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                            #lowerci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                            #upperci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                            .SDcols = vars
  ]
  combined.impact.sims.unconsumed.summary <- combined.impact.sims.unconsumed.summary[, (vars) := NULL]
  
  ##### current.intake
 
  current.intake.impact.sims.dt<-as.data.table(current.intake.impact.sims)
  current.intake.impact.sims.summary <- current.intake.impact.sims.dt[, ":="(mean_impact_current_intake_ = rowMeans(.SD, na.rm = TRUE),
                                                                         sd_impact_current_intake = apply(.SD, 1, sd),
                                                                         median_impact_current_intake = apply(.SD, 1, median),
                                                                         lowerci_95_impact_current_intake = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                         upperci_95_impact_current_intake = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                  #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                  #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                  .SDcols = vars
  ]
  current.intake.impact.sims.summary <- current.intake.impact.sims.summary[, (vars) := NULL]
  
  current.intake.impact.consumed.sims.dt<-as.data.table(current.intake.impact.sims.consumed)
  current.intake.impact.consumed.sims.summary <- current.intake.impact.consumed.sims.dt[, ":="(mean_impact_current_intake_consumed = rowMeans(.SD, na.rm = TRUE),
                                                                                               sd_impact_current_intake_consumed = apply(.SD, 1, sd),
                                                                                               median_impact_current_intake_consumed = apply(.SD, 1, median),
                                                                                               lowerci_95_impact_current_intake_consumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                               upperci_95_impact_current_intake_consumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                        .SDcols = vars
  ]
  current.intake.impact.consumed.sims.summary <- current.intake.impact.consumed.sims.summary[, (vars) := NULL]
  
  
  current.intake.impact.unconsumed.sims.dt<-as.data.table(current.intake.impact.sims.unconsumed)
  current.intake.impact.unconsumed.sims.summary <- current.intake.impact.unconsumed.sims.dt[, ":="(mean_impact_current_intake_unconsumed = rowMeans(.SD, na.rm = TRUE),
                                                                                               sd_impact_current_intake_unconsumed = apply(.SD, 1, sd),
                                                                                               median_impact_current_intake_unconsumed = apply(.SD, 1, median),
                                                                                               lowerci_95_impact_current_intake_unconsumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                               upperci_95_impact_current_intake_unconsumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                        .SDcols = vars
  ]
  current.intake.impact.unconsumed.sims.summary <- current.intake.impact.unconsumed.sims.summary[, (vars) := NULL]
  
  #####CF.intake. _CF_intake
  CF.intake.impact.sims.dt<-as.data.table(CF.intake.impact.sims)
  CF.intake.impact.sims.summary <- CF.intake.impact.sims.dt[, ":="(mean_impact_CF_intake = rowMeans(.SD, na.rm = TRUE),
                                               sd_impact_CF_intake = apply(.SD, 1, sd),
                                               median_impact_CF_intake = apply(.SD, 1, median),
                                               lowerci_95_impact_CF_intake = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                               upperci_95_impact_CF_intake = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                        .SDcols = vars
  ]
  CF.intake.impact.sims.summary <- CF.intake.impact.sims.summary[, (vars) := NULL]
  
  CF.intake.impact.consumed.sims.dt<-as.data.table(CF.intake.impact.sims.consumed)
  CF.intake.impact.consumed.sims.summary <- CF.intake.impact.consumed.sims.dt[, ":="(mean_impact_CF_intake_consumed = rowMeans(.SD, na.rm = TRUE),
                                                                                               sd_impact_current_CF_consumed = apply(.SD, 1, sd),
                                                                                               median_impact_CF_intake_consumed = apply(.SD, 1, median),
                                                                                               lowerci_95_impact_CF_intake_consumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                               upperci_95_impact_CF_intake_consumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                        .SDcols = vars
  ]
  CF.intake.impact.consumed.sims.summary <- CF.intake.impact.consumed.sims.summary[, (vars) := NULL]
  
  
  CF.intake.impact.unconsumed.sims.dt<-as.data.table(CF.intake.impact.sims.unconsumed)
  CF.intake.impact.unconsumed.sims.summary <- CF.intake.impact.unconsumed.sims.dt[, ":="(mean_impact_CF_intake_unconsumed = rowMeans(.SD, na.rm = TRUE),
                                                                                                   sd_impact_CF_intake_unconsumed = apply(.SD, 1, sd),
                                                                                                   median_impact_CF_intake_unconsumed = apply(.SD, 1, median),
                                                                                                   lowerci_95_impact_CF_intake_unconsumed = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                                   upperci_95_impact_CF_intake_unconsumed = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                            #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                            #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                            .SDcols = vars
  ]
  CF.intake.impact.unconsumed.sims.summary <- CF.intake.impact.unconsumed.sims.summary[, (vars) := NULL]
  
  ####
  
  ########inedible
  impact.sims.inedible.dt<-as.data.table(impact.sims.inedible)
  impact.sims.inedible.summary <- impact.sims.inedible.dt[, ":="(mean_impact_inedible = rowMeans(.SD, na.rm = TRUE),
                                                                     sd_impact_inedible = apply(.SD, 1, sd),
                                                                     median_impact_inedible = apply(.SD, 1, median),
                                                                     lowerci_95_impact_inedible = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                     upperci_95_impact_inedible = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                              #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                              #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                              .SDcols = vars
  ]
  impact.sims.inedible.summary <- impact.sims.inedible.summary[, (vars) := NULL]
  
  substitution.impact.sims.inedible.dt<-as.data.table(substitution.impact.sims.inedible)
  substitution.impact.sims.inedible.summary <- substitution.impact.sims.inedible.dt[, ":="(mean_substitution_impact_inedible = rowMeans(.SD, na.rm = TRUE),
                                                                                               sd_substitution_impact_inedible = apply(.SD, 1, sd),
                                                                                               median_substitution_impact_inedible = apply(.SD, 1, median),
                                                                                               lowerci_95_substitution_impact_inedible = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                               upperci_95_substitution_impact_inedible = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                        #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                        #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                        .SDcols = vars
  ]
  substitution.impact.sims.inedible.summary <- substitution.impact.sims.inedible.summary[, (vars) := NULL]
  
  combined.impact.sims.inedible.dt<-as.data.table(combined.impact.sims.inedible)
  combined.impact.sims.inedible.summary <- combined.impact.sims.inedible.dt[, ":="(mean_combined_impact_inedible = rowMeans(.SD, na.rm = TRUE),
                                                                               sd_combined_impact_inedible = apply(.SD, 1, sd),
                                                                               median_combined_impact_inedible = apply(.SD, 1, median),
                                                                               lowerci_95_combined_impact_inedible = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                               upperci_95_combined_impact_inedible = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                        #lowerci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                        #upperci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                        .SDcols = vars
  ]
  combined.impact.sims.inedible.summary <- combined.impact.sims.inedible.summary[, (vars) := NULL]
  
  
  ########wasted
  impact.sims.wasted.dt<-as.data.table(impact.sims.wasted)
  impact.sims.wasted.summary <- impact.sims.wasted.dt[, ":="(mean_impact_wasted = rowMeans(.SD, na.rm = TRUE),
                                                                 sd_impact_wasted = apply(.SD, 1, sd),
                                                                 median_impact_wasted = apply(.SD, 1, median),
                                                                 lowerci_95_impact_wasted = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                 upperci_95_impact_wasted = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                          #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                          #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                          .SDcols = vars
  ]
  impact.sims.wasted.summary <- impact.sims.wasted.summary[, (vars) := NULL]
  
  substitution.impact.sims.wasted.dt<-as.data.table(substitution.impact.sims.wasted)
  substitution.impact.sims.wasted.summary <- substitution.impact.sims.wasted.dt[, ":="(mean_substitution_impact_wasted = rowMeans(.SD, na.rm = TRUE),
                                                                                           sd_substitution_impact_wasted = apply(.SD, 1, sd),
                                                                                           median_substitution_impact_wasted = apply(.SD, 1, median),
                                                                                           lowerci_95_substitution_impact_wasted = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                           upperci_95_substitution_impact_wasted = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                    #lowerci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                    #upperci_90_substitution_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                    .SDcols = vars
  ]
  substitution.impact.sims.wasted.summary <- substitution.impact.sims.wasted.summary[, (vars) := NULL]
  
  combined.impact.sims.wasted.dt<-as.data.table(combined.impact.sims.wasted)
  combined.impact.sims.wasted.summary <- combined.impact.sims.wasted.dt[, ":="(mean_combined_impact_wasted = rowMeans(.SD, na.rm = TRUE),
                                                                                       sd_combined_impact_wasted = apply(.SD, 1, sd),
                                                                                       median_combined_impact_wasted = apply(.SD, 1, median),
                                                                                       lowerci_95_combined_impact_wasted = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                                                                       upperci_95_combined_impact_wasted = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                                                                #lowerci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                                                                #upperci_90_combined_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                                                                .SDcols = vars
  ]
  combined.impact.sims.wasted.summary <- combined.impact.sims.wasted.summary[, (vars) := NULL]
  
  
  all.impact.sims.summary<-Reduce(function(...) merge(..., all = T),list(delta.sims.summary, 
                                                                         impact.sims.summary, substitution.impact.sims.summary, combined.impact.sims.summary,
                                                                         impact.sims.consumed.summary, substitution.impact.sims.consumed.summary, combined.impact.sims.consumed.summary,
                                                                         impact.sims.unconsumed.summary, substitution.impact.sims.unconsumed.summary, combined.impact.sims.unconsumed.summary,
                                                                         current.intake.impact.sims.summary, current.intake.impact.consumed.sims.summary, current.intake.impact.unconsumed.sims.summary,
                                                                         CF.intake.impact.sims.summary, CF.intake.impact.consumed.sims.summary, CF.intake.impact.unconsumed.sims.summary,
                                                                         impact.sims.inedible.summary, substitution.impact.sims.inedible.summary, combined.impact.sims.inedible.summary,
                                                                         impact.sims.wasted.summary, substitution.impact.sims.wasted.summary, combined.impact.sims.wasted.summary))
  
  return(all.impact.sims.summary)
}

###simple version of get summary stats to use on just one set of outputs (as oppposed to 4)
get.summary.stats.simple<-function(impact.sims, vars)
{
  impact.sims.dt<-as.data.table(impact.sims)
  impact.sims.summary <- impact.sims.dt[, ":="(mean_impact = rowMeans(.SD, na.rm = TRUE),
                                               sd_impact = apply(.SD, 1, sd),
                                               median_impact = apply(.SD, 1, median),
                                               lowerci_95_impact = apply(.SD, 1, function(x) quantile(x, .025, na.rm = T)),
                                               upperci_95_impact = apply(.SD, 1, function(x) quantile(x, .975, na.rm = T))),
                                        #lowerci_90_impact = apply(.SD, 1, function(x) quantile(x, .05, na.rm = T)),
                                        #upperci_90_impact = apply(.SD, 1, function(x) quantile(x, .95, na.rm = T))),
                                        .SDcols = vars
  ]
  impact.sims.summary <- impact.sims.summary[, (vars) := NULL]
  
  return(impact.sims.summary)
}

summary.stats.by.strata<-function(impact.sims.allgroups, strata.combos.list, pop.sims)
{
  cols.to.sum<-paste("X", 1:n.sims, sep="")
  pop.cols.to.sum<-paste("pop", 1:n.sims, sep="")
  
  sims.data.table<-as.data.table(impact.sims.allgroups)
  pop.sims.data.table<-as.data.table(pop.sims)
  
  strata.combos.sims<-list()
  strata.combos.sims.long<-list()
  strata.combos.summary.stats.long<-list()
  strata.combos.summary.stats.wide<-list()
  strata.combos.summary.stats.misc<-list()
  
  
  pop.strata.combos.sims<-list() 
  
  percapita.sims<-list()
  percapita.sims.long<-list()
  percapita.summary.stats.long<-list()
  percapita.summary.stats.wide<-list()
  percapita.summary.stats.misc<-list()
  
  
  for(j in 1:length(strata.combos.list))
  {
    labels<-c(strata.combos[[j]], "outcome", "outcome_unit")
    
    strata.combos.sims[[j]]<-sims.data.table[,lapply(.SD, sum), by=c(strata.combos.list[[j]], "outcome", "outcome_unit"), 
                                             .SDcols=cols.to.sum]
    
    pop.strata<-strata.combos[[j]][!(strata.combos.list[[j]] %in% c("Foodgroup", "broad_foodgroup"))]
    
    if(is.null(pop.strata))
    {
      pop.strata.combos.sims[[j]]<-pop.sims.data.table[,lapply(.SD, sum), .SDcols=cols.to.sum]
    }else
    {
      pop.strata.combos.sims[[j]]<-pop.sims.data.table[,lapply(.SD, sum), by=c(pop.strata), .SDcols=cols.to.sum]
    }
    names(pop.strata.combos.sims[[j]])<-gsub(x=names(pop.strata.combos.sims[[j]]), pattern="X", replacement="pop")
    
    if(length(intersect(names(strata.combos.sims[[j]]), names(pop.strata.combos.sims[[j]])))==0)
    {
      merged<-cbind(strata.combos.sims[[j]], pop.strata.combos.sims[[j]])
    }
    if(length(intersect(names(strata.combos.sims[[j]]), names(pop.strata.combos.sims[[j]])))!=0)
    {
      merged<-merge(strata.combos.sims[[j]], pop.strata.combos.sims[[j]])
    }
    
    
    merged.names<-merged[,..labels]
    percapita.sims[[j]]<-merged[,..cols.to.sum]/merged[,..pop.cols.to.sum]
    percapita.sims[[j]]<-cbind(merged.names, percapita.sims[[j]])
    
    strata.combos.sims.long[[j]]<-melt(strata.combos.sims[[j]], id.vars = c(strata.combos[[j]], "outcome", "outcome_unit"), measure.vars=cols.to.sum)
    strata.combos.summary.stats.long[[j]]<-strata.combos.sims.long[[j]][, .(quantile(value, c(.025, .5, .975))), by=c(strata.combos[[j]], "outcome", "outcome_unit")]
    strata.combos.summary.stats.misc[[j]]<-strata.combos.sims.long[[j]][, .(mean(value), sd(value)), by=c(strata.combos[[j]], "outcome")]
    names(strata.combos.summary.stats.misc[[j]])[(dim(strata.combos.summary.stats.misc[[j]])[2]-1):dim(strata.combos.summary.stats.misc[[j]])[2]]<-c("mean", "SD")
    
    strata.combos.summary.stats.long[[j]]<-cbind(strata.combos.summary.stats.long[[j]], 
                                                 as.factor(rep(c("lower_bound (2.5th percentile)", "median", "upper_bound  (97.5th percentile)"), dim(strata.combos.summary.stats.long[[j]])[1]/3)))
    #strata.combos.summary.stats.wide[[j]]<-dcast(strata.combos.summary.stats.long[[j]], "Comorbid_count + Insurance + Race ~ V2", value.var="V1")
    strata.combos.summary.stats.wide[[j]]<-dcast(strata.combos.summary.stats.long[[j]], paste(paste(strata.combos[[j]], collapse=' + '), " + outcome + outcome_unit", " ~ V2", sep=""),
                                                 value.var="V1")    
    #strata.combos.summary.stats.wide[[j]]<-merge(strata.combos.summary.stats.wide[[j]], strata.combos.summary.stats.misc[[j]], by=c("age_gp", "sex_gp", "outcome"))
    strata.combos.summary.stats.wide[[j]]<-merge(strata.combos.summary.stats.wide[[j]], strata.combos.summary.stats.misc[[j]], 
                                            by=intersect(names(strata.combos.summary.stats.wide[[j]]), names(strata.combos.summary.stats.misc[[j]])))
    
    
    ##nperhaps in future, make it so that the function doesn't rely on specific names for intake_units, outcome, outcome_unit 
    
    percapita.sims.long[[j]]<-melt(percapita.sims[[j]], id.vars = c(strata.combos.list[[j]], "outcome", "outcome_unit"), measure.vars=cols.to.sum)
    percapita.summary.stats.long[[j]]<-percapita.sims.long[[j]][, .(quantile(value, c(.025, .5, .975))), by=c(strata.combos.list[[j]], "outcome", "outcome_unit")]
    percapita.summary.stats.misc[[j]]<-percapita.sims.long[[j]][, .(mean(value), sd(value)), by=c(strata.combos.list[[j]], "outcome")]
    names(percapita.summary.stats.misc[[j]])[(dim(strata.combos.summary.stats.misc[[j]])[2]-1):dim(strata.combos.summary.stats.misc[[j]])[2]]<-c("mean", "SD")
    
    percapita.summary.stats.long[[j]]<-cbind(percapita.summary.stats.long[[j]], 
                                             as.factor(rep(c("lower_bound (2.5th percentile)", "median", "upper_bound  (97.5th percentile)"), 
                                                           dim(percapita.summary.stats.long[[j]])[1]/3)))
    #strata.combos.summary.stats.wide[[j]]<-dcast(strata.combos.summary.stats.long[[j]], "Comorbid_count + Insurance + Race ~ V2", value.var="V1")
    percapita.summary.stats.wide[[j]]<-dcast(percapita.summary.stats.long[[j]], paste(paste(strata.combos.list[[j]], collapse=' + '), " + outcome + outcome_unit", " ~ V2", sep=""),
                                             value.var="V1")    
    percapita.summary.stats.wide[[j]]<-merge(percapita.summary.stats.wide[[j]], percapita.summary.stats.misc[[j]], 
                                             by=intersect(names(percapita.summary.stats.wide[[j]]), names(percapita.summary.stats.misc[[j]])))
    
    
  }
  
  return(list("summary.output" = strata.combos.summary.stats.wide, 
              "full.sims.output" = strata.combos.sims,
              "summary.output.percapita"= percapita.summary.stats.wide,
              "full.sims.output.percapita" = percapita.sims)
  )
}