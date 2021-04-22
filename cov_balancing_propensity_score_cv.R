library(tidyverse)
library(CBPS)
library(rstanarm)

# This script documents some of the cross-validation for the Bayesian logistic regression model in the final submission.
# These Bayesian logistic regression are to some extent deliberately underfit to avoid overfitting odd artefacts in 
# the training data (and also, because it is easier to debug what is going on). They  were initially just a testbed 
# for trying out adding random effects, regularization and some prior insights into super-simple by-prediction-task 
# models that seemed to do okay. Then they turned out to perform really well, once we added the post-processing of predictions
#
# The model structure is from previous elastic net regression that followed a similar idea (see stratify_the_models.R)
# of having effectively different models for different prediction tasks (generics, flu vaccines/insulin, new drugs,
# LCM etc.). Those previous experiments more or less determined the model structure seen below.
# The problem had been that disease specific information was hard to account for, so I tried adding
# random therapeutic area effects (newta variable). However, trying to do that with the lme4 package runs into two problems:
#   1. convergence issues / hierarchical standard deviation estimated as 0 (but we know there should be some across-TA variation)
#   2. no ability to regularize easily (like LASSO, elastic net etc.)
# Both of these key issues can be addressed by Bayesian logistic regression, which in fact allows us to specify some of our
# prior beliefs on what direction certain effects should go.

#####################################################
# Generate the data
#####################################################
adat1 <- read_csv(file="/files/feat_bjoern_trial.csv") %>%
  dplyr::select(-row_id, -starts_with("intsponsor"), -starts_with("intpers"), -intidentifiedsites) %>%
  group_by(DrugKey, indicationkey) %>%
  summarize_if(is.numeric, median) %>% # categorical: names(sort(summary(as_factor(c(0,1,1,2,3))), decreasing = T))[1]
  ungroup() %>%
  dplyr::select(-newta)

adat2 <- read_csv(file="/files/feat_bjoern_trial.csv") %>%
  dplyr::select(-row_id, -starts_with("intsponsor"), -starts_with("intpers"), -intidentifiedsites) %>%
  #mutate(newta = ifelse(!is.na(newta), newta, ifelse( strDiseaseType == "[Type 2 Diabetes|NAFLD]", )))
  group_by(DrugKey, indicationkey) %>%
  summarize(number_of_studies=n(),
            newta = as.numeric(names(sort(summary(as_factor(newta)), decreasing = T, na.last = T)[1]))) %>%
  ungroup() 

adat <- left_join(adat1, adat2, by=c("DrugKey", "indicationkey")) %>%
  ungroup() %>%
  mutate(prediction_task = case_when(
      insulin + fluvacc>0 ~ 1,
      is_a_generic>0  ~ 2,
      time_since_first_outcome==0 & prior_approval==0 ~ 3,
      prior_approval == 1 ~ 4,
      time_since_first_outcome>0 ~ 5,
      TRUE ~ 6),
    efficacy_assessed = ifelse(is.na(efficacy_assessed),0 , efficacy_assessed),
    pk_aspects_na = 1L*is.na(pk_aspects),
    pk_aspects = ifelse(is.na(pk_aspects),mean(pk_aspects, na.rm=T) , pk_aspects),
    ta7byorph = ta7*orphan_indication,
    myclass = as_factor(1L*is.na(outcome)),
    dtmeanclu50_t345 = scale( dtmeanclu50 * (prediction_task)>=3),
    dtmeancll50_t345 = scale( dtmeancll50 * (prediction_task)>=3),
    tameanclu50_t345 = scale( tameanclu50* (prediction_task)>=3),
    tameancll50_t345 = scale( tameancll50* (prediction_task)>=3),
    orphan_indication_t345 = orphan_indication *  (prediction_task>=3),
    rel_ph2_size_dis_t345 = rel_ph2_size_dis *  (prediction_task>=3),
    ta7_t345 = ta7 * (prediction_task>=3),
    ta7_orph_t345 = ta7*any_oncology_orphan_designation *  (prediction_task>=3),
    newta_t345 = newta * (prediction_task>=3),
    monoclonal_no_inn_t35 = monoclonal_no_inn * (prediction_task %in% c(3,5)),
    prediction_task = as_factor(prediction_task),
    newta = as_factor(newta),
    mean_trialendscore = scale(mean_trialendscore),
    rel_ph2_size_dis_t345 = scale(rel_ph2_size_dis_t345)    
    ) %>% 
  mutate_if(is.numeric , replace_na, replace = 0) %>%
  left_join( distinct(dplyr::select(read_csv(file="/files/feat_bjoern_trial.csv"),
                                    DrugKey, indicationkey, GenericName, strDiseaseType)), by=c("DrugKey", "indicationkey")) %>%
  ungroup() %>%
  arrange(DrugKey, indicationkey)
  
#####################################################
# Obtain Covariate-balancing propensity score
#####################################################
cbps1 <- CBPS( myclass ~ newta1 + newta2 + newta3 + newta4 + newta5 + newta6 + newta7 + newta8 + newta9 + newta10 +
                 newta11 + newta12 + newta13 + newta14 + newta15 + newta16 + newta17 + newta18 + newta19 + newta20 +
                 + newta21 + newta22 + newta23 +  pk_aspects + time_since_first_outcome + 
                 prediction_task*(mean_trialendscore + spons4) + tameancll50 + tameanclu50 + orphan_indication + any_oncology_orphan_designation + 
                 ta7 + ta7*any_oncology_orphan_designation + any_nononc_orphan_designation +
        + rel_ph2_size_dis + safety + mean_trialendscore + monoclonal_no_inn + pk_aspects + after2007 +
        drugtype1 + drugtype2 + drugtype4 + drugtype5 + drugtype6 + drugtype7 + drugtype8 + drugtype9 + unwilling_to_pay_12k +
           + fluvacc + insulin + rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + tameanclu50_t345 + tameancll50_t345,
      data=filter(adat, phaseendyear<2018))

adat <- adat %>%
  mutate( propensity = 0)

adat$propensity[adat$phaseendyear<2018] = cbps1$fitted.values

# adat %>% 
#   ggplot(aes(x=phaseendyear, y=propensity, col=as_factor(myclass))) + geom_jitter(alpha=0.6)

# # This is what the model looked like in lme4
# glmerfit1 <- glmer(outcome ~ (1 | newta_t345) + prediction_task*(mean_trialendscore + spons4) + fluvacc + insulin +
#          rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7_orph_t345 + orphan_indication_t345,
#       family=binomial, 
#       data=filter(adat, myclass==0), 
#       weights = propensity,
#       control=glmerControl(optCtrl=list(maxfun=50000)))

#####################################################
# Set-up past-future cross-validation
#####################################################

# Create list with what records are in the validation set for each (overlapping) CV-fold
new_cv_splits <- filter(adat, myclass==0) %>%
  dplyr::select(DrugKey, indicationkey) %>%
  left_join(read_csv("/files/new_cv_splits.csv"), by=c("DrugKey", "indicationkey")) %>%
  dplyr::select(foldid, set, DrugKey, indicationkey) %>%
  arrange(foldid, DrugKey, indicationkey)

cv_index = list()
cixl = 22:26
for (fi in 1:length(cixl)){ #1:5){
  cv_index[[fi]] = (1:length(filter(new_cv_splits, foldid==cixl[fi])$set))[filter(new_cv_splits, foldid==cixl[fi])$set=="val"]
}
#table(new_cv_splits$foldid, new_cv_splits$set)

# Probably just want to use folds 22 to 26 (5 fold for quick evaluation?)?

#####################################################
# Function for assessing on a CV split
#####################################################
sglmerfitting <- function( ){
  # prior_location = c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5),
  # prior_scale = rep(1, length(prior_location)),
  # int_location = -2.5,
  # int_sd = 0.5,
  # sd_exp_rate = 1, 
  # val_indices=1:10,
  # yearcutoff = 2007
  # Fit Bayesian random effects logistic regression model with somewhat informative priors
  sglmerfit1 <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                             rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                           family=binomial(link="logit"), 
                           data=filter(adat, myclass==0)[-val_indices, ] %>% 
                             filter(phaseendyear>yearcutoff), 
                           weights = propensity,
                           chains=4, cores=4,
                           iter=2000, # Reduced number of iterations to speed up CV
                           seed=2020,
                           prior_intercept=normal(int_location, int_sd, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                           # (last year of training data has 6.45%, but previous year higher)
                           # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                           prior_aux = exponential(sd_exp_rate, autoscale=FALSE),
                           prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                             # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                             # -2.5 for stuff where I more or less assume deterministic relationship)
                             # Things I have opionions on:
                             # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                             # very bad: no INN requsted for a monoclonal antibody
                             # good: large Phase 2 by standards of the disease area
                             df=4, 
                             location=prior_location,  #-0.5
                             scale=prior_scale))
  
  preds = posterior_predict(object=sglmerfit1, 
                            newdata=filter(adat, myclass==0)[val_indices, ],
                            seed=2020)
  
  filter(adat, myclass==0)[val_indices, ] %>%
    mutate(pred=pmin(0.999, pmax(0.001, colMeans(preds))),
           logloss = - number_of_studies * ( log(pred)*outcome + log(1-pred)*(1-outcome) )) %>%
    summarize( logloss = sum(logloss)/sum(number_of_studies) )
    
}
         
#####################################################
# Run experiments
#####################################################
                  
experiments = list()
experimentid = 0

# >2009 cutoff, weakly informative priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,0,  0,0,0,0,  0,0,0,0)
  prior_scale = rep(2, length(prior_location))
  int_location = 0
  int_sd = 2
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2009
  
  experiments[[experimentid]] = sglmerfitting()
}

# >2007 cutoff, weakly informative priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,0,  0,0,0,0,  0,0,0,0)
  prior_scale = rep(2, length(prior_location))
  int_location = 0
  int_sd = 2
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2007
  
  experiments[[experimentid]] = sglmerfitting()
}

# no cutoff, weakly informative priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,0,  0,0,0,0,  0,0,0,0)
  prior_scale = rep(2, length(prior_location))
  int_location = 0
  int_sd = 2
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 1990
  
  experiments[[experimentid]] = sglmerfitting()
}

# >2009 cut-off, not so opinionated prior
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,-0.5,  -0.5,0.5,-0.5,0,  0,0,0,0)
  prior_scale = rep(1, length(prior_location))
  int_location = -2
  int_sd = 1
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2009
  
  experiments[[experimentid]] = sglmerfitting()
}


# >2009 cut-off, somewhat informative more opinionated priors -0.5 on mab-INN feature
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,-0.5,  -0.5,0.5,-0.5,0,  -0.5,-0.5,-0.5,-0.5)
  prior_scale = rep(1, length(prior_location))
  int_location = -2.5
  int_sd = 0.5
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2009
  
  experiments[[experimentid]] = sglmerfitting()
}


# >2009 cut-off, somewhat informative more opinionated priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5)
  prior_scale = rep(1, length(prior_location))
  int_location = -2.5
  int_sd = 0.5
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2009
  
  experiments[[experimentid]] = sglmerfitting()
}

# >2007 cut-off, somewhat informative more opinionated priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5)
  prior_scale = rep(1, length(prior_location))
  int_location = -2.5
  int_sd = 0.5
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 2007
  
  experiments[[experimentid]] = sglmerfitting()
}


# all data, somewhat informative more opinionated priors
for (fi in 1:length(cv_index)){
  experimentid = experimentid + 1
  # Putting function parameters in global envir because of silly behavior of rstanarm
  prior_location = c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5)
  prior_scale = rep(1, length(prior_location))
  int_location = -2.5
  int_sd = 0.5
  sd_exp_rate = 1
  val_indices=cv_index[[fi]]
  yearcutoff = 1990
  
  experiments[[experimentid]] = sglmerfitting()
}

#####################################################
# Summarize experiments
#####################################################
tibble(scenario=rep(1:(length(experiments)/5), each=5),
       foldid = rep(1:5, length(experiments)/5),
       logloss = as_vector(map(1:length(experiments), function(x) experiments[[x]]$logloss))) %>%
  group_by(scenario) %>%
  summarize(logloss = mean(logloss))

#####################################################
# Conclusions
##################################################################################################################
# These cross-validation scores are pretty close together, but only using >2009 data seems best 
# (that seems like relatively consistent thing across all models in cross-validation & on the leaderboard).
# Almost same ordering seen on public LB, but somewhat informative seemed to outperform weakly informative.
# Different ways of setting priors seem to perform similarly in cross validation, but logically several of
# the informative priors we specified should be meaningful (and do not seem so strong that they should be badly
# miscalibrated, which seems supported by the outcomes on the public LB). Especially, the mab with no INN feature
# should be a very strong (almost deterministic) signal, so went with a strong prior since it does not hurt 
# in cross-validation (and LB seems to be better that way).
##################################################################################################################
# weakly informative >2009:    0.446
# not so opinionated somewhat informative  > 2009: 0.484
# like scenario below, but less strong prior on mab/no-INN effect: 0.483
# opinionated somewhat informative > 2009: 0.456
# opinionated somewhat informative > 2007: 0.484
# opinionated somewhat informative all data:  0.495
