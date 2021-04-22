##########################################################################
##########################################################################
# prediction R script of Team Insight-Out based on example script by
# Sharada Mohanty.
##########################################################################
##########################################################################

library(tidyverse)
library(glmnet)
library(xgboost)
library(rstanarm)

################################################################################################################
################################################################################################################
## Expected ENVIRONMENT Variables
################################################################################################################
# 
# * AICROWD_TEST_DATA_PATH : Absolute Path to the Test CSV file
# * AICROWD_PREDICTIONS_OUTPUT_PATH  : path where you are supposed to write the output predictions.csv
################################################################################################################
AICROWD_TRAIN_DATA_PATH <- Sys.getenv(c("AICROWD_TRAIN_DATA_PATH"), unset="/shared_data/data/training_data/") #training_data_2015_split_on_outcome.csv
AICROWD_TEST_DATA_PATH <- Sys.getenv(c("AICROWD_TEST_DATA_PATH"), unset="/shared_data/data/test_data_full/testing_phase2_release.csv")
AICROWD_PREDICTIONS_OUTPUT_PATH <- Sys.getenv(c("AICROWD_PREDICTIONS_OUTPUT_PATH"), unset="random_prediction.csv")

print(paste("AICROWD_TEST_DATA_PATH : " , AICROWD_TEST_DATA_PATH))
print(paste("AICROWD_TRAIN_DATA_PATH : " , AICROWD_TRAIN_DATA_PATH))
print(paste("AICROWD_PREDICTIONS_OUTPUT_PATH : " , AICROWD_PREDICTIONS_OUTPUT_PATH))

#############################################################################################
# Read data
#############################################################################################
# Load the Test Data
print("Load the Test Data")
df <- read.csv(AICROWD_TEST_DATA_PATH)

# The AI Crowd evaluation server cannot handle our data wrangling code (apparently a 5GB RAM 
# limit and a higher limit was not possible for cost reasons). Thus, we had to read the wrangled
# data in from csv and merge them in by DrugKey/indicationkey. See 
# https://discourse.aicrowd.com/t/submissions-get-killed-without-any-error-message/2544/6
# https://discourse.aicrowd.com/t/evaluation-error/2521/19
# for details. Checked by email with Nick that will have to do it this way.
# For the programs to create the data wrangled csv files, see the data_wrangling.tar file
# in the repository.

#############################################################################################
# Useful helper functions
#############################################################################################
# log(1+x)
log1p <- function(x) {
    y = 1 + x
    z = y - 1
    return( ifelse( z == 0, x, x * log(y) / z))
}

# log(1-x)
log1m <- function(x){
    return(log1p(-x))
}

# log(1-exp(x))
log1m_exp <- function(x){
    # Note: expm1(x) := exp(x)-1
    return( ifelse(x > -0.693147, log(-expm1(x)), log1m(exp(x)) ) )
}

# log(exp(x)-exp(y))
log_diff_exp <- function(x,y){
    return( x + log1m_exp(y - x) )
}

# Logit function in numerically stable form
logit <- function(x){
    return( -log_diff_exp(-log(x),0) )
}

# log(1+exp(x))
log1p_exp <- function(x){
    return( ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x))) )
}

# Numerically stable inverse logit function
inv_logit <- function(x){
    return(  ifelse(x==Inf, 1, ifelse(x==-Inf, 0, exp(x - log1p_exp(x)))) )
}

# Helper function to re-scale a set of probabilities by adding a constant
# on the logit scale in order to obtain a target mean probability.
revise_probs <- function(probs, target_mean_prob = 0.5){
    add_on_logit <- function(probs, x){
        r = logit(probs) + x # log(probs)-log(1-probs) + x
        inv_logit(r) # exp(r)/(1+exp(r))
    }
    psumdelta2 <- function(x){
        ( mean( add_on_logit(probs, x) ) - target_mean_prob )^2
    }
    optimres = optim(par=0, fn=psumdelta2, method="Brent", lower=-100, upper=100)
    add_value = optimres$par
    add_on_logit(probs, add_value)
}

 
#############################################################################################
# xgboost model follwed by averaging mean, max and min predicted logits
#############################################################################################

# print("Obtain features")
adatx <- read_csv(file="feat_bjoern_trial.csv") %>%
    filter(is.na(outcome)) %>%
    dplyr::select(-outcome, -GenericName, -strDiseaseType, 
                  -predgroup, -casewgt1, -casewgt2, -newta, -sponstype, -logitoffset, 
                  -intidentifiedsites, -starts_with("intpersonid"), -starts_with("intsponsorid"))

xgb_model = read_rds("xgb_trials2_7jan2020.rds")

the_predictions = adatx %>%
    dplyr::select(DrugKey, indicationkey, ta7, phaseendyear) %>%
    mutate( pred = predict(xgb_model, 
               newdata=data.matrix(dplyr::select(adatx, 
                                                 -row_id, -DrugKey, -indicationkey)))) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize(phaseendyear = round(mean(phaseendyear)),
              ta7 = round(mean(ta7)),
              meanpred = mean(pred),
              minlogit = min(logit(pred)),
              meanlogit = mean(logit(pred)),
              maxlogit = max(logit(pred)),
              meanlogit_back = inv_logit(meanlogit)) %>%
    ungroup()

elnet3 <- read_rds("xgboost_meta_avg2_7jan2020.rds")

the_predictions2a <- the_predictions %>%
    mutate(model="xgboost-ridge",
           meta_logit = predict(elnet3, as.matrix(dplyr::select(the_predictions, minlogit, meanlogit, maxlogit))),
           prob_approval = inv_logit(meta_logit)) #(meanpred + meanlogit_back + exp(meta_logit)/(1+exp(meta_logit)))/3 )# meanpred) #= exp(meta_logit)/(1+exp(meta_logit)))

#############################################################################################
# 2nd xgboost model follwed by averaging mean, max and min predicted logits
#############################################################################################

# print("Obtain features")
adatx <- read_csv(file="feat_bjoern_trial.csv") %>%
    filter(is.na(outcome)) %>%
    dplyr::select(-outcome, -GenericName, -strDiseaseType, 
                  -predgroup, -casewgt1, -casewgt2, -newta, -sponstype, -logitoffset, 
                  -intidentifiedsites, -starts_with("intpersonid"), -starts_with("intsponsorid"))

#xgb_model = read_rds("xgb_trials3_10jan2020.rds")
xgb_model = read_rds("xgb_trials3_10jan2020.rds")

the_predictions = adatx %>%
    dplyr::select(DrugKey, indicationkey, ta7, phaseendyear) %>%
    mutate( pred = predict(xgb_model, 
                           newdata=data.matrix(dplyr::select(adatx, 
                                                             -row_id, -DrugKey, -indicationkey)))) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize(phaseendyear = round(mean(phaseendyear)),
              ta7 = round(mean(ta7)),
              meanpred = mean(pred),
              minlogit = min(logit(pred)),
              meanlogit = mean(logit(pred)),
              maxlogit = max(logit(pred)),
              meanlogit_back = inv_logit(meanlogit)) %>%
    ungroup()

elnet3 <- read_rds("xgboost_meta_avg3_10jan2020.rds")

the_predictions2b <- the_predictions %>%
    mutate(model="xgboost-ridge",
           meta_logit = predict(elnet3, as.matrix(dplyr::select(the_predictions, minlogit, meanlogit, maxlogit))),
           prob_approval = inv_logit(meta_logit)) 
    # %>% mutate(prob_approval = ifelse(prob_approval<0.1, prob_approval,
    #                               0.1+(prob_approval-0.1)*0.944)) # This particular xgboost is super aggressive and assigns very high probabilities 
    #                               # (99.8% is clearly absurd) and perhaps (?) should be down-scaled perhaps above 0.1 to at most 0.95? 

    
#############################################################################################
# Code for Bayesian logistic regression
#############################################################################################

# print("Obtain features")
adat1 <- read_csv(file="feat_bjoern_trial.csv") %>%
    dplyr::select(-row_id, -starts_with("intsponsor"), -starts_with("intpers"), -intidentifiedsites) %>%
    group_by(DrugKey, indicationkey) %>%
    summarize_if(is.numeric, median) %>% # categorical: names(sort(summary(as_factor(c(0,1,1,2,3))), decreasing = T))[1]
    ungroup() %>%
    dplyr::select(-newta)

adat2 <- read_csv(file="feat_bjoern_trial.csv") %>%
    dplyr::select(-row_id, -starts_with("intsponsor"), -starts_with("intpers"), -intidentifiedsites) %>%
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
        
        
        time_to_deadline = case_when(
            myclass==0 & outcome==1 ~ 2015-phaseendyear,
            myclass==0 & outcome==0 ~ 2014-phaseendyear,
            myclass==1 ~ 2018.5-phaseendyear,
            TRUE ~ 999),
        ttdemax=time_to_deadline/(time_to_deadline+2),
        
        time_to_deadline=scale(time_to_deadline),
        ttdemax=scale(ttdemax),

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
    mutate_if(is.numeric , replace_na, replace = 0)

# Get model
sglmerfit1 <- read_rds("sglmerfit1b.rds")

# Get predictions
# Could not quickly find an option to get more predictive draws
# from rstanarm than the number MCMC samples.
# Simple fix for the moment: just loop to get 100,000 samples per
# indication - should give a decent enough approx. even for low probabilities
maxdraws = 25 
for (draw in 1:maxdraws){ 
    preds = posterior_predict(object=sglmerfit1,
                              newdata=filter(adat, myclass==1),
                              #draws=4000, #Draws will be 4000 by default (so comment out), if we have more samples, doing more samples may be worthwhile
                              seed=2019+draw)
    if (draw==1) {
        # each column of preds is a MCMC sample predicting 0 or 1,
        # so take mean of columns to get probability
        pa=colMeans(preds)
    } else {
        pa = pa + colMeans(preds)
    }
}
pa = pa/maxdraws

the_predictions3 <- filter(adat, myclass==1) %>%
    dplyr::select(DrugKey, indicationkey, phaseendyear, ta7) %>%
    mutate(model="BLR-CBPS",
           prob_approval = pa) 

#############################################################################################
# Model averaging/stacking
#############################################################################################

# # "Real" model based stacking performed worse on public-LB in logloss by 0.02 than just taking the mean of probabilities
# # Since that was a quite deep nesting of cross-validations, the stacking appeared rather unreliable to us and we
# # relied on the public-LB to decide to discard it.
# elnet_stack1 = read_rds("elnet_stack_int.rds")
# elnet_stack2 = read_rds("elnet_stack_noint.rds")
# 
# the_predictions4a <- the_predictions2 %>%
#     dplyr::select(DrugKey, indicationkey, phaseendyear, ta7, model, prob_approval) %>%
#     bind_rows(the_predictions3) %>%
#     group_by(DrugKey, indicationkey, phaseendyear, ta7) %>%
#     pivot_wider(names_from = model, values_from=prob_approval) %>%
#     mutate(cvpred = logit(`xgboost-ridge`),
#            pred = logit( ifelse(`BLR-CBPS`<=1e-06, 1e-06, ifelse(`BLR-CBPS`>=1-1e-06, 1-1e-06, `BLR-CBPS`)) ),
#            meanpred = (`xgboost-ridge`+`BLR-CBPS`)/2,
#            logit_mean = logit(meanpred)) %>%
#     ungroup()
# 
# the_predictions4 <- the_predictions4a %>%
#     mutate( simple_mean = (`xgboost-ridge` + `BLR-CBPS`)/2,
#             stack_logit1 = predict(elnet_stack1,
#                                   newx = as.matrix( dplyr::select(the_predictions4a, pred, cvpred, logit_mean) )),
#             stack_logit2 = predict(elnet_stack2,
#                                    newx = as.matrix( dplyr::select(the_predictions4a, pred, cvpred, logit_mean) )),
#             prob_approval = (simple_mean + inv_logit(stack_logit1) + inv_logit(stack_logit2))/3,
#             model="Stacked BLR-CBPS + xbgboost-ridge",
#              )

# This is the for just calculating the simple mean:
the_predictions4 <- the_predictions2a %>%
    dplyr::select(DrugKey, indicationkey, phaseendyear, ta7, model, prob_approval) %>%
    bind_rows(dplyr::select(the_predictions2b, DrugKey, indicationkey, phaseendyear, ta7, model, prob_approval)) %>%
    group_by(DrugKey, indicationkey, phaseendyear, ta7) %>%
    summarize(prob_approval=mean(prob_approval)) %>% # Average the xgboosts first (this way they get wgts of 0.25 and do not excessively push predictions)
    ungroup() %>%
    bind_rows(the_predictions3) %>%
    group_by(DrugKey, indicationkey, phaseendyear, ta7) %>%
    summarize(prob_approval=mean(prob_approval)) %>% # Now average in BLR
    # avg. on logit scale inv_logit(mean(logit(prob_approval))) less good on pub-LB
    mutate(model="Model-averaged xbgboost-ridge + genetic-xgboost-ridge + BLR-CBPS") %>%
    ungroup()

#############################################################################################
# Post-processing of predictions
#############################################################################################    

# Step 1:
# Reflects that nothing ending Phase 2 in 2018 can realistically be approved by mid-2019. 
# There might be exceptions for Oncology indications with high unmet need, and log-loss 
# cannot be calculated for probability=0, so set to something low like 0.001
the_predictions5 <- the_predictions4 %>%
    mutate(prob_approval = ifelse( phaseendyear<2018, 
                                   prob_approval, 
                                   ifelse(phaseendyear==2018,
                                          pmin(prob_approval, pmax(0.001, 0.1 * prob_approval)), 
                                          pmin(prob_approval, 0.0005)))) # 0.0005: small gain in cae of 0 approvals with EoPh2 in 2019, okay if 1-3, bad thereafter

# Step 2:
# Reflects that Phase 2 programs without development were declared abandoned after 540 days (=1.5 years),
# so programs with Phase 2 end before 2014 (assuming test set projects are approved 2016 or later)
# must have had a Phase 3 program. Thus, for these projects Phase 3 to approval success rates apply.
# Based on Lo et al. (no leaked data only up to 2015) these are 35.5% (Oncology) and 63.6% (non-Onc).
# However, some of these probabilities are getting pretty high (>90%), which seems implausibly certain.
# So, I decdied it might be better to be cautious and average between the predictions and the re-scaled predictions.
the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==1 ] = 
    (the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==1 ] +
         revise_probs(the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==1 ],
                      target_mean_prob=0.355))/2

the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==0 ] = 
    (the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==0 ] +
         revise_probs(the_predictions5$prob_approval[ the_predictions5$phaseendyear<2014 & the_predictions5$ta7==0 ],
                      target_mean_prob=0.636))/2

# espi. xgboost seems too excitable and predicts very high probabilities for
#    some projects. If these were generics or fluc vaccines, maybe that would be okay,
#    but for new drug developments that seems crazy. Thus, capped predicted probabilities as 0.75.

the_predictions5 <- the_predictions5 %>%
    mutate(prob_approval = ifelse(prob_approval<0.5, prob_approval,
                                  0.5+(prob_approval-0.5)/1.65)) %>%
    mutate(prob_approval = case_when( # Never make an absuredly overconfident prediction, limit probabilities [0.001, 0.999]
        prob_approval>0.999 ~ 0.999,
        prob_approval<0.001 & phaseendyear<=2018 ~ 0.001,
        prob_approval<0.0005 & phaseendyear>2018 ~ 0.0005, # Special exception: for 2019 allow 0.05% instead of 0.1% at lower end.
        TRUE ~ prob_approval))

################################################################################
# Create predictions dataframe
################################################################################

print("Create predictions dataframe")
prediction_df <- df %>%
    dplyr::select(DrugKey, indicationkey, row_id) %>%
    left_join(dplyr::select(the_predictions5, DrugKey, indicationkey, prob_approval), 
              by=c("DrugKey", "indicationkey") ) %>%
    dplyr::select(row_id, prob_approval) %>%
    arrange(row_id) %>%
    data.frame()

# Write final prediction_df to a CSV file at the appropriate location
print("write predictions csv file")
write.csv(prediction_df, file=AICROWD_PREDICTIONS_OUTPUT_PATH, quote=FALSE, row.names=FALSE)
print(c("Prediction File Written to : ", AICROWD_PREDICTIONS_OUTPUT_PATH))
