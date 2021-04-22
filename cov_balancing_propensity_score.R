library(tidyverse)
library(CBPS)
library(rstanarm)

# This script trains some Bayesian logistic regression models used in the final submission (see the same file name 
# with "_cv" included for some cross-validation. 
#
# The Bayesian logistic regression is the second best model in our ensemble. xgboost is a good bit better, but not
# by so much that it was worth it to ignore BLR (especially since predictions had a low correlation of 0.5 - pretty
# low for two models trained on the same data) in the final ensemble.
#
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
#
# We kind of ran out of time to do more with this model (and focussed more on xgboost due to its promise):
# - There's probably plenty of effects that should be in the model, candidates: Alzheimer's disease, drugs without INN (other than mab)...
# - Exploring using the horseshoe prior might have been interesting, it is meant to be good at dealing with a
#   high dimensional feature vector
# - Other ways of creating CBPS-types of weights might be attractive. The `cobalt` R package offers a lot of options.`

#################################################################
# Read data, some model specific feature engineering
#################################################################

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
  summarize(newta = as.numeric(names(sort(summary(as_factor(newta)), decreasing = T, na.last = T)[1]))) %>%
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
    
    dtmeanclu50_t345n = scale( dtmeanclu50 * (prediction_task>=3)),
    dtmeancll50_t345n = scale( dtmeancll50 * (prediction_task>=3)),
    tameanclu50_t345n = scale( tameanclu50* (prediction_task>=3)),
    tameancll50_t345n = scale( tameancll50* (prediction_task>=3)),
    
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
  ungroup()

####################################################  
# Obtain Covariate-balancing propensity score
####################################################
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

adat %>% 
  ggplot(aes(x=phaseendyear, y=propensity, col=as_factor(myclass))) + geom_jitter(alpha=0.6)


################################################################
# For background the related lme4 model that had various issues
################################################################

# glmerfit1 <- glmer(outcome ~ (1 | newta_t345) + prediction_task*(mean_trialendscore + spons4) + fluvacc + insulin +
#          rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7_orph_t345 + orphan_indication_t345,
#       family=binomial, 
#       data=filter(adat, myclass==0), 
#       weights = propensity,
#       control=glmerControl(optCtrl=list(maxfun=50000)))

######################################################
# Bayesian logistic regression fitting using rstanarm
######################################################

## Main model used in final sumbission, used more MCMC samples than other fits below (seems to help a little)
# Fit Bayesian random effects logistic regression model with somewhat informative priors
# even limit to 2010 or later
sglmerfit1b <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7 + tameanclu50_t345, 
                          #  + ta7_orph_t345 + orphan_indication_t345#spons4) + # fluvacc + insulin +
                          family=binomial(link="logit"), 
                          data=filter(adat, myclass==0 & phaseendyear>=2010), 
                          weights = propensity,
                          chains=4, cores=4,
                          iter=10000,
                          seed=2020,
                          prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                          # (last year of training data has 6.45%, but previous year higher)
                          # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                          prior_aux = exponential(1, autoscale=FALSE),
                          prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                            # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                            # -2.5 for stuff where I more or less assume deterministic relationship)
                            # Things I have opionions on:
                            # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                            # very bad: no INN requsted for a monoclonal antibody
                            # good: large Phase 2 by standards of the disease area
                            df=4, 
                            location=c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                            scale=1))
write_rds(sglmerfit1b, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1b2.rds", compress = "bz2")

# Fit Bayesian random effects logistic regression model with somewhat informative priors
sglmerfit1 <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                           rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                         family=binomial(link="logit"), 
                         data=filter(adat, myclass==0), 
                         weights = propensity,
                         chains=4, cores=4,
                         seed=2020,
                         prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                                                            # (last year of training data has 6.45%, but previous year higher)
                                                            # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                         prior_aux = exponential(1, autoscale=FALSE),
                         prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                           # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                           # -2.5 for stuff where I more or less assume deterministic relationship)
                           # Things I have opionions on:
                           # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                           # very bad: no INN requsted for a monoclonal antibody
                           # good: large Phase 2 by standards of the disease area
                           df=4, 
                           location=c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                           scale=1))

# # This is how one gets predictions for the testset, probably want to add fixed seed
# preds = posterior_predict(object=sglmerfit1, newdata=filter(adat, myclass==1))
# 
# preds = adat %>%
#   filter(myclass==1) %>%
#   mutate(pred=colMeans(preds))
# 
# preds %>%
#   group_by(phaseendyear) %>%
#   summarize(n=n(),
#             meanpred=mean(pred))

write_rds(sglmerfit1, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1.rds", compress = "bz2")

# Fit Bayesian random effects logistic regression model with somewhat informative priors
# Limit to post 2007 records
sglmerfit1a <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                           rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                         family=binomial(link="logit"), 
                         data=filter(adat, myclass==0 & phaseendyear>2007), 
                         weights = propensity,
                         chains=4, cores=4,
                         seed=2020,
                         prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                         # (last year of training data has 6.45%, but previous year higher)
                         # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                         prior_aux = exponential(1, autoscale=FALSE),
                         prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                           # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                           # -2.5 for stuff where I more or less assume deterministic relationship)
                           # Things I have opionions on:
                           # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                           # very bad: no INN requsted for a monoclonal antibody
                           # good: large Phase 2 by standards of the disease area
                           df=4, 
                           location=c(0,0,0,-0.5,  -0.5,0.5,-2.5,0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                           scale=1))
write_rds(sglmerfit1a, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1a.rds", compress = "bz2")


# Fit Bayesian random effects logistic regression model with somewhat informative priors
# leave out the monoclonal by INN feature
sglmerfit1c <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                           rel_ph2_size_dis_t345 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                         family=binomial(link="logit"), 
                         data=filter(adat, myclass==0), 
                         weights = propensity,
                         chains=4, cores=4,
                         seed=2020,
                         prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                         # (last year of training data has 6.45%, but previous year higher)
                         # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                         prior_aux = exponential(1, autoscale=FALSE),
                         prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                           # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                           # -2.5 for stuff where I more or less assume deterministic relationship)
                           # Things I have opionions on:
                           # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                           # very bad: no INN requsted for a monoclonal antibody
                           # good: large Phase 2 by standards of the disease area
                           df=4, 
                           location=c(0,0,0,-0.5,  -0.5,0.5, 0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                           scale=1))
write_rds(sglmerfit1c, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1c.rds", compress = "bz2")

# Fit Bayesian random effects logistic regression model with somewhat informative prior
# omit INN based feature, train on post-2007
sglmerfit1d <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                            rel_ph2_size_dis_t345 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                          family=binomial(link="logit"), 
                          data=filter(adat, myclass==0 & phaseendyear>2007), 
                          weights = propensity,
                          chains=4, cores=4,
                          seed=2020,
                          prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                          # (last year of training data has 6.45%, but previous year higher)
                          # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                          prior_aux = exponential(1, autoscale=FALSE),
                          prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                            # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                            # -2.5 for stuff where I more or less assume deterministic relationship)
                            # Things I have opionions on:
                            # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                            # very bad: no INN requsted for a monoclonal antibody
                            # good: large Phase 2 by standards of the disease area
                            df=4, 
                            location=c(0,0,0,-0.5,  -0.5,0.5, 0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                            scale=1))
write_rds(sglmerfit1d, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1d.rds", compress = "bz2")

# Fit Bayesian random effects logistic regression model with somewhat informative prior
# omit INN based feature, train on post-2010
sglmerfit1e <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + #spons4) + # fluvacc + insulin +
                            rel_ph2_size_dis_t345 + ta7 + tameanclu50_t345, #  + ta7_orph_t345 + orphan_indication_t345
                          family=binomial(link="logit"), 
                          data=filter(adat, myclass==0 & phaseendyear>=2010), 
                          weights = propensity,
                          chains=4, cores=4,
                          seed=2020,
                          prior_intercept=normal(-2.5, 0.5, autoscale=FALSE), # The mean of the prior corresponds to about 7.5% approvals,
                          # (last year of training data has 6.45%, but previous year higher)
                          # SD of 0.5 leads to + SD being ~ 18% (reasonable somewhat informed prior)
                          prior_aux = exponential(1, autoscale=FALSE),
                          prior = student_t( # Follow prior recommendations of https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations 
                            # for student_t priors, but make them somehwat informative (location=0.5 or -0.5 if I have a good guess, 
                            # -2.5 for stuff where I more or less assume deterministic relationship)
                            # Things I have opionions on:
                            # bad: drug with previous failed indications, poor Phase 2 outcome (e.g. stopped for safety)
                            # very bad: no INN requsted for a monoclonal antibody
                            # good: large Phase 2 by standards of the disease area
                            df=4, 
                            location=c(0,0,0,-0.5,  -0.5,0.5, 0,  -0.5,-0.5,-0.5,-0.5),  #-0.5
                            scale=1))
write_rds(sglmerfit1e, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1e.rds", compress = "bz2")

# # Note explored: Could one even take the predictions for the test data and refit with those?! 
# preds = posterior_predict(object=sglmerfit1, newdata=filter(adat, myclass==1))
# 
# adat3 <- filter(adat, myclass==1) %>%
#   mutate(outcome=colMeans(preds),
#          propensity=0.2) %>%
#   bind_rows(filter(adat, myclass==0))
# 
# sglmerfit2 <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*(mean_trialendscore + spons4) + fluvacc + insulin +
#                            rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7_orph_t345 + orphan_indication_t345,
#                          family=binomial(link="logit"), 
#                          data=adat3, 
#                          weights = propensity,
#                          prior_intercept=normal(-2, 2),
#                          prior_aux = exponential(1),
#                          prior = student_t(4, 0, 2))
# write_rds(sglmerfit2, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit2.rds", compress = "bz2")

##############################################################
# Exploring the posterior and model coefficients further
##############################################################
library(bayesplot)

bayesplot_theme_set(theme_default(base_size = 12,
                                  base_family = "sans"))
color_scheme_set("red")
mcmc_intervals(sglmerfit1b ,
               pars=rev(c("(Intercept)",
                          "prediction_task2",
                          "prediction_task3",
                          "prediction_task4",
                          "prediction_task5",
                          "rel_ph2_size_dis_t345",
                          "monoclonal_no_inn_t35",
                          "mean_trialendscore",
                          "prediction_task2:mean_trialendscore",
                          "prediction_task3:mean_trialendscore",
                          "prediction_task4:mean_trialendscore",
                          "prediction_task5:mean_trialendscore"))) +
  scale_shape_manual(values=c(19)) +
  geom_vline(xintercept=0) +
  scale_y_discrete(
    breaks=rev(c("(Intercept)",
                 "prediction_task2",
                 "prediction_task3",
                 "prediction_task4",
                 "prediction_task5",
                 "rel_ph2_size_dis_t345",
                 "monoclonal_no_inn_t35",
                 "mean_trialendscore",
                 "prediction_task2:mean_trialendscore",
                 "prediction_task3:mean_trialendscore",
                 "prediction_task4:mean_trialendscore",
                 "prediction_task5:mean_trialendscore")),
    labels=rev(c("Log-odds for generics (intercept)",
                 "Flu vaccine/insulin",
                 "Innovative drug (no prior outcomes)",
                 "Innovative with prior approval",
                 "Innovative with prior failed indications",
                 "Relative Phase 2 size for disease",
                 "Monoclonal without INN",
                 "Mean trial outcome score (MTOS)",
                 "MTOS adjustment for flu vacc/ins)",
                 "MTOS adjustment for (no prior outcomes)",
                 "MTOS adjustment for with prior approval",
                 "MTOS adjustment for prior failed ind."
    ))) +
  xlab("Coefficients on logit scale")

#names(sglmerfit1$coefficients)
mcmc_intervals(sglmerfit1b , regex_pars = " newta|ta7") +
  geom_vline(xintercept=0) +
  scale_y_discrete(
    breaks=rev(c("b[(Intercept) newta_t345:0]",
                 "b[(Intercept) newta_t345:1]",  "b[(Intercept) newta_t345:2]",
                 "b[(Intercept) newta_t345:3]",  "b[(Intercept) newta_t345:4]",
                 "b[(Intercept) newta_t345:5]",  "b[(Intercept) newta_t345:6]",
                 "b[(Intercept) newta_t345:7]",  "b[(Intercept) newta_t345:8]",
                 "b[(Intercept) newta_t345:9]",  "b[(Intercept) newta_t345:10]",
                 "ta7",
                 "b[(Intercept) newta_t345:11]", "b[(Intercept) newta_t345:12]",
                 
                 "b[(Intercept) newta_t345:13]", "b[(Intercept) newta_t345:14]",
                 "b[(Intercept) newta_t345:15]", "b[(Intercept) newta_t345:16]",
                 "b[(Intercept) newta_t345:17]", "b[(Intercept) newta_t345:18]",
                 "b[(Intercept) newta_t345:19]", "b[(Intercept) newta_t345:20]",
                 "b[(Intercept) newta_t345:21]", "b[(Intercept) newta_t345:22]",
                 "b[(Intercept) newta_t345:23]")),
    labels=rev(c("Generics/flu vaccines/insulin",
                 "Women's health",
                 "Vaccines",
                 "Urinary",
                 "Transplant",
                 "Respiratory",
                 "Pain",
                 "Other",
                 "Renal",
                 "Hepatic",
                 "Ophthalmology",
                 "Oncology",
                 "Solid tumours (diff. to Oncology)",
                 "Liquid tumours (diff. to Oncology)",
                 
                 "Neuroscience",
                 "Neurodegenerative",
                 "Metabolic",
                 "Viral infections",
                 "(Bacterial) infections",
                 "Hematology",
                 "Gastrointestinal",
                 "Dermatology",
                 "Cardiovascular",
                 "Bone",
                 "Autoimmune"
    ))
    # , labels=rev(c("Log-odds for generics (intercept)",
    #              "Flu vaccine/insulin",
    #              "Innovative drug (no prior outcomes)",
    #              "Innovative with prior approval",
    #              "Innovative with prior failed indications",
    #              "Relative Phase 2 size for disease",
    #              "Monoclonal without INN",
    #              "Oncology",
    #              "Mean trial outcome score (MTOS)",
    #              "MTOS adjustment for flu vacc/ins)",
    #              "MTOS adjustment for (no prior outcomes)",
    #              "MTOS adjustment for with prior approval",
    #              "MTOS adjustment for prior failed ind."))
  ) +
  xlab("Coefficients on logit scale")
