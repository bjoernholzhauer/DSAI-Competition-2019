library(Matrix)
library(tidyverse)
library(xgboost)
library(DEoptim)

# Script to tune xgboost using a differential evolution algorithm for global optimization of the
# leader-board style loss function using xgboost (taking approximate blending into account)

###########################################################################################################
# Read trial level training data, keep only the variables we hope to be reliable, limit to data after 2007
###########################################################################################################
adat <- read_csv(file="/files/feat_bjoern_trial.csv") %>%
  arrange(row_id) %>%
  filter(!is.na(outcome) & phaseendyear>2007)

# Remove any variables that xgboost cannot or should not use 
# (e.g. character, non-one-hot-encoded categroical like newta, DrugKey, indicationkey), or should not use (like row_id, casewgt1, etc. )
# or that I simply do not trust, because I do not understand them (the intpersonid... and intsponsorid... variables)
traindat = adat %>% 
  dplyr::select(-row_id, -DrugKey, -indicationkey, -predgroup, -outcome, -GenericName, -strDiseaseType, 
                -predgroup, -casewgt1, -casewgt2, -newta, -sponstype, -logitoffset, -intidentifiedsites, -starts_with("intpersonid"), -starts_with("intsponsorid"))
#%>%   mutate(newta = as_factor(newta),sponstype = as_factor(sponstype)) # instead of using factor, use existing multi-class membership indictors

###########################################################################################################
# Create list with what records are in the validation set for each (overlapping) CV-fold
###########################################################################################################
new_cv_splits <- dplyr::select(adat, row_id, DrugKey, indicationkey) %>%
  left_join(read_csv("/files/new_cv_splits.csv"), by=c("DrugKey", "indicationkey")) %>%
  dplyr::select(foldid, set, row_id) %>%
  arrange(foldid, row_id)

cv_index = list()
pick_folds = 22:26
for (fi in 1:5){ 
  cv_index[[fi]] = (1:length(filter(new_cv_splits, foldid==pick_folds[fi])$set))[filter(new_cv_splits, foldid==pick_folds[fi])$set=="val"]
}

# Define function that fits xgboost, does predictions and calculates actual target loss 
call_xgbcv <- function(x){
  params = list(booster = "gbtree", 
                objective = "binary:logistic", 
                eta=0.05, 
                gamma=x[1], 
                max_depth=x[2], 
                min_child_weight=x[3], 
                subsample=x[4], 
                colsample_bytree=x[5])
  
  # Fit xgboost with custom cross-valdiation
  xc1 <- xgb.cv( data = data.matrix(traindat),
                 label=adat$outcome,
                 weight = adat$casewgt2, # Case wgts that down weight 2008 to 2011 proportional to how much the success rate is too high compared to 2012+
                 params = params,
                 nrounds = 10000, 
                 folds=cv_index, 
                 metrics = "logloss",
                 showsd = F, 
                 print_every_n = 100, 
                 early_stopping_rounds = 20, 
                 maximize = F,
                 #verbose=F,
                 nthread = 8,
                 callbacks = list(cb.cv.predict(save_models = TRUE)))
  
  # Get cross-validation prediction results, 
  # aggregate predictions for a project across the trials (min or mean of probabilities),
  # calcuilate logLoss on this basis
  cv_results = tibble(fold=1:length(cv_index)) %>% 
    mutate(res=map(fold, function(x) filter(adat, !is.na(outcome))[cv_index[[x]], c("predgroup", "outcome")] %>% 
                     bind_cols( tibble( pred = predict(xc1$models[[x]], newdata=data.matrix(traindat[cv_index[[x]],])) ) ) )) %>% 
    unnest(res) %>%
    group_by(fold,predgroup) %>%
    summarize(n=n(),
              outcome = max(outcome),
              minpred=min(pred),
              meanpred=mean(pred),
              maxpred=max(pred),
              meanlogit = mean( log(pred)-log(1-pred) ),
              ridge1 = -0.2294281*(log(minpred)-log(1-minpred)) + 0.2751232*meanlogit + 0.8999464*(log(maxpred)-log(1-maxpred)),
              ridge_mean = ( exp(ridge1)/(1+exp(ridge1)) + meanpred )/2 ,
              logloss = -( log(ridge_mean)*outcome + log(1-ridge_mean)*(1-outcome) ) * n ) %>%
    ungroup() %>%
    group_by(fold) %>%
    summarize(logloss = sum(logloss) / sum(n)) %>%
    ungroup()
  
  # Summarize across folds (note that SD would not make so much sense, because of differing fold sizes)
  cv_results2 <- cv_results %>%
    summarize(logloss=mean(logloss))
  # %>%
  #   mutate(fold=0) %>%
  #   bind_rows(cv_results) %>%
  #   bind_rows(cv_results %>%
  #               summarize(loglossmin=mean(ifelse(fold>=22, loglossmin, NA_real_), na.rm=T),
  #                         loglossmean=mean(ifelse(fold>=22, loglossmean, NA_real_), na.rm=T)) %>%
  #               mutate(fold=-1))
  
  return( cv_results2$logloss )
}

#call_xgbcv(c(0.5, 5, 3, 0.5, 0.2))

# for (current_depth in 3:11){
#   print(paste("current depth ", current_depth))
#   experimentid = experimentid + 1
#   params$max_depth = current_depth
#   experiments[[experimentid]] = call_xgbcv()
# }

fnmap_f <- function(x){
  c(x[1], round(x[2]), round(x[3]), x[4], x[5])
}

ptm = proc.time()
set.seed(1234)
resdeopt <- DEoptim(fn=call_xgbcv, 
        lower=c(0, 3, 2, 0.2, 0.1), # gamma, max_depth, min_child_weight, subsample, colsample_bytree
        upper=c(0.6, 10, 10, 1, 0.75), 
        control = DEoptim.control(initialpop = matrix(c(0, 8, 4, 1, 0.1,
                                             0, 8, 4, 1, 0.2,
                                             0.001, 8, 4, 1, 0.25,
                                             0.01, 8, 3, 0.5, 0.15,
                                             0.1, 7, 8, 0.9, 0.12,
                                             0, 5, 6, 0.9, 0.4,
                                             0.001, 6, 3, 0.85, 0.149,
                                             0.1, 5, 2.5, 0.5, 0.19,
                                             0.05, 4, 4, 0.4, 0.3,
                                             0.03, 5, 9, 0.7, 0.3,
                                             0, 10, 4, 1, 0.3,
                                             0.1, 6, 2, 0.9, 0.18,
                                             0.01, 9, 4, 1, 0.12,
                                             0, 7, 4, 1, 0.18,
                                             0.1, 6, 4, 1, 0.26,
                                             0.015, 9, 3, 0.5, 0.14,
                                             0.12, 6, 8, 0.9, 0.13,
                                             0, 5, 5, 0.9, 0.42,
                                             0.002, 4, 3, 0.85, 0.1,
                                             0.12, 6, 2.5, 0.5, 0.05,
                                             0.06, 5, 4, 0.4, 0.25,
                                             0, 4, 9, 0.68, 0.29,
                                             0.01, 11, 4, 0.99, 0.28,
                                             0.101, 5, 2, 0.85, 0.17,
                                             0, 6, 4, 0.75, 0.125), 
                                           nrow = 25, ncol=5, byrow = T),
                                  NP=25, itermax=100, storepopfrom = 1),
        fnMap = fnmap_f)
proc.time() - ptm

# write_rds(resdeopt, 
#          "/files/genetic_opt_xgboost_DEoptim1.rds", 
#          compress = "bz2")

# Results:
# Iteration: 100 bestvalit: 0.228899 bestmemit:    0.024047    6.000000    2.000000    0.950466    0.736181

params = list(booster = "gbtree", 
              objective = "binary:logistic", 
              eta=0.01, 
              gamma=0.024047, 
              max_depth=6, 
              min_child_weight=2, 
              subsample=0.950466, 
              colsample_bytree=0.736181)

###########################################################################################################
# Use full CV scheme now 
###########################################################################################################
new_cv_splits <- dplyr::select(adat, row_id, DrugKey, indicationkey) %>%
  left_join(read_csv("/files/new_cv_splits.csv"), by=c("DrugKey", "indicationkey")) %>%
  dplyr::select(foldid, set, row_id) %>%
  arrange(foldid, row_id)

cv_index = list()
for (fi in sort(unique(new_cv_splits$foldid))){ #1:5){
  cv_index[[fi]] = (1:length(filter(new_cv_splits, foldid==fi)$set))[filter(new_cv_splits, foldid==fi)$set=="val"]
}


##################################################################
# Now switch to a lower learning rate (eta=0.01) to find a 
# good number of rounds to use when refitting with the whole data
##################################################################
set.seed(2020)
fin_train = xgb.cv( data = data.matrix(traindat),
                    label=adat$outcome,
                    weight = adat$casewgt2,
                    params = params,
                    nrounds = 10000, 
                    folds=cv_index, 
                    metrics = "logloss",
                    showsd = T, 
                    print_every_n = 10, 
                    early_stopping_rounds = 50, 
                    maximize = F,
                    nthread = 8,
                    callbacks = list(cb.cv.predict(save_models = TRUE)))

# Best iteration: [3228]	train-logloss:0.004151+0.000384	test-logloss:0.027630+0.054358

# best iteration 3228, add 10% so 3551

##################################################
# Train a final model
##################################################
# Maybe this version trained too long (training error dropping to zero is concerning)
main_model = xgboost( data = data.matrix(traindat),
                      label = adat$outcome,
                      weight = adat$casewgt2,
                      params = params,
                      nrounds = 3551,
                      metrics = "logloss",
                      print_every_n = 10,
                      maximize = F)

# write_rds(main_model, 
#          "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/xgb_trials3_10jan2020.rds", 
#          compress = "bz2")

# 1800 seems to get the train loss below the training loss in CV, perhaps that's a good try
# However, public LB shows no improvement (in fact slight worsening, but it's all within the scope of chance)
main_model2 = xgboost( data = data.matrix(traindat),
                      label = adat$outcome,
                      weight = adat$casewgt2,
                      params = params,
                      nrounds = 1800,
                      metrics = "logloss",
                      print_every_n = 10,
                      maximize = F)
# write_rds(main_model2,
#          "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/xgb_trials4_10jan2020.rds",
#          compress = "bz2")



######################################################
# Now look into how we average the predictions that
# are currently by trial.
######################################################

# Create dataset for optimal averaging of trial predictions based on validation-fold predictions.
# It is very critical to respect that the same project must either 
# always be in training or in validaiton, toherwise we will get overfitting.
foravging = tibble(fold=1:length(cv_index)) %>% 
  mutate(res=map(fold, function(x) filter(adat, !is.na(outcome))[cv_index[[x]], c("predgroup", "outcome")] %>% 
                   bind_cols( tibble( pred = predict(fin_train$models[[x]], newdata=data.matrix(traindat[cv_index[[x]],])) ) ) )) %>% 
  unnest(res) %>%
  group_by(fold,predgroup) %>%
  summarize(n=n(),
            logn=log(n),
            outcome = max(outcome),
            minpred=min(pred),
            minlogit = min( log(pred) - log(1-pred)),
            meanpred=mean(pred),
            meanlogit = mean( log(pred) - log(1-pred)),
            maxpred=max(pred),
            maxlogit = max( log(pred) - log(1-pred)),
            sd_pred = ifelse(n==1, 0, sd(pred)),
            sd_logit = ifelse(n==1, 0, sd( log(pred) - log(1-pred)))) %>%
  ungroup() %>%
  mutate(avgfold = sample(x=1:10, size=n(), replace=T)) %>% # Hand-specifying folds is important, don't want a project in multiple folds
  group_by(predgroup) %>%
  mutate(avgfold = first(avgfold)) %>%
  ungroup()

#write_rds(foravging, "/home/desktop1/Documents/holzhbj1/for_xgboost_averaging_10jan2020.rds", compress = "bz2")  
#foravging <- read_rds("/home/desktop1/Documents/holzhbj1/for_xgboost_averaging_10jan2020.rds")


library(glmnet)
library(doParallel)

# This seems like a decent  model to use for averaging min, max and mean predicted probability in a project across studies.
# I explored some other thing, but it seems woryrying that I might be overfitting noise in the cross-validation,
# if I add too many other things to the model (like untransformed probabilities, number of studies, other project level features),
# and any improvement seemed minimal. Thus, better to keep it simple.
# Relaxed-ridge regression makes sense, because I feel like everything should be in the model (possibly even slightly underfit)
# and should implicitly do something like Platt scaling (which is actually super-convenient, because we wish to optimize logLoss,
# not accuracy or F1-score and xgboost is sometimes not well calibrated otherwise).
cl <- makePSOCKcluster(12)
registerDoParallel(cl)
set.seed(2020)
elnet4 <- cv.glmnet(y=foravging$outcome,
                    x=as.matrix(foravging %>% dplyr::select(minlogit, meanlogit, maxlogit)), 
                    family = "binomial",
                    nfolds=10,
                    foldid=foravging$avgfold, # This is the super important bit (forgot it at first) - keep a whole project always in train or validation fold
                    alpha=0,
                    gamma=seq(0, 1, 0.05),
                    type.measure = "deviance", 
                    relax = T,
                    intercept=F,
                    parallel=TRUE,
                    keep=F)
stopCluster(cl)

write_rds(elnet4,
          "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/xgboost_meta_avg3_10jan2020.rds",
          compress="bz2")

###########################################################################################################################
# Try Platt scaling - abandoned due to non-intuitive behavior, it seems like ridge-regression alone already does a good job
###########################################################################################################################
# 
# cl <- makePSOCKcluster(12)
# registerDoParallel(cl)
# set.seed(2020)
# elnet4 <- cv.glmnet(y=foravging$outcome,
#                     x=as.matrix(foravging %>% dplyr::select(minlogit, meanlogit, maxlogit)), 
#                     family = "binomial",
#                     nfolds=10,
#                     foldid=foravging$avgfold, # This is the super important bit (forgot it at first) - keep a whole project always in train or validation fold
#                     alpha=0,
#                     gamma=seq(0, 1, 0.05),
#                     type.measure = "deviance", 
#                     relax = T,
#                     intercept=F,
#                     parallel=TRUE,
#                     keep=T)
# stopCluster(cl)
# 
# 
# foravging2 <- foravging %>% 
#   mutate(cvpred=elnet4$fit.preval$`g:0.15`[,1],
#          cvpredprob = exp(cvpred)/(1+exp(cvpred))) 
# 
# foravging2 %>%
#   ggplot(aes(x=outcome, y=cvpred)) + 
#   geom_jitter()
# 
# # platt<-glm(y~x,data=fitdata,family=binomial) 
# # calib_pred<-predict(platt,newdata=temp,type="response")
# 
# platt <- glm(outcome ~ cvpredprob,
#              data=foravging2,
#              family=binomial) 
# write_rds(platt,
#           "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/platt1_10jan2020.rds",
#           compress="bz2")
# # library(lme4)
# # platt2 <- glmer(outcome ~ (1|predgroup) + cvpredprob, 
# #                 data=foravging2,
# #                 family=binomial) 
# require(splines)
# platt3 <- glm(outcome ~ bs(cvpred, knots=c(-2,0,1)),
#              data=foravging2 , #%>% group_by(predgroup,outcome) %>% summarize( cvpred=mean(cvpred)) %>% ungroup()
#              family=binomial) 
# inv_logit( predict(platt3, tibble(cvpred=logit(c(0.01,0.2,0.5,0.8,0.999)))) )
# 
# inv_logit(predict(platt, tibble(cvpredprob=c(0.01,0.2,0.5,0.8,0.999) )))