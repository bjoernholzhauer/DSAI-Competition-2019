library(tidyverse)
library(CBPS)
library(rstanarm)

# This script trains some Bayesian logistic regression models used in the final submission (see the same file name 
# with "_cv" included for some cross-validation: compared to the other similar script, some extra covariates are added 
# (further explanations below remain the same)

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
    unwilling_to_pay_12k_t35 = unwilling_to_pay_12k  * (prediction_task %in% c(3,5)),
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
        drugtype1 + drugtype2 + drugtype4 + drugtype5 + drugtype6 + drugtype7 + drugtype8 + drugtype9 + unwilling_to_pay_12k_t35 +
           + fluvacc + insulin + rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + tameanclu50_t345n + tameancll50_t345n,
      data=filter(adat, phaseendyear<2018))

adat <- adat %>%
  mutate( propensity = 0)
adat$propensity[adat$phaseendyear<2018] = cbps1$fitted.values

adat %>% 
  ggplot(aes(x=phaseendyear, y=propensity, col=as_factor(myclass))) + geom_jitter(alpha=0.6)

######################################################
# Bayesian logistic regression fitting using rstanarm
######################################################

## Main model used in final sumbission, used more MCMC samples than other fits below (seems to help a little)
# Fit Bayesian random effects logistic regression model with somewhat informative priors
# even limit to 2010 or later
sglmerfit1b <- stan_glmer(outcome ~ (1 | newta_t345) + prediction_task*mean_trialendscore + rel_ph2_size_dis_t345 + monoclonal_no_inn_t35 + ta7 + 
                            tameanclu50_t345n + unwilling_to_pay_12k_t35, 
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
                            location=c(0,0,0,-0.5,  -0.5,0.5,-2.5,0, 0.5,-0.5, -0.5,-0.5,-0.5,-0.5),  #-0.5
                            scale=1))
write_rds(sglmerfit1b, "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/sglmerfit1f.rds", compress = "bz2")