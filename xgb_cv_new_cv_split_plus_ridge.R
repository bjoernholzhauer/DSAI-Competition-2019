library(Matrix)
library(tidyverse)
library(xgboost)

# This R code is training a xgboost model to make predictions per trial and 
# then averaging the predictions on a project level using a relaxed-ridge regression
# fit to projects on tha validation folds.
#
# This is the top performing model in our final ensemble. It brutually uses all features
# except for those that we distrust or do not understand. We failed to get a clear understanding
# of the features in the provded dataset that start with "intpersonid" or "intsponsorid" and worry
# that these increment with time for the same sponsor, which might make them very fragile for
# an extrapolation task.
#
# The cross-validation for xgboost hyperparameter optimization uses a past-predicting-the-future
# 26-fold split. However, validation folds do overlap, but they are according to a couple of 
# different ways of splitting past-vs-future (sometimes with some future projects in the training data).
# One of the key tricks is to calculate the logloss for the cross-validation folds in the exact same
# way as for the leaderboard (but with the simplifying assumption that one just takes the mean of the
# predictions for each trial - one could perhap shave iterated and averaged using the ridge regression
# coefficients).
# 
# Subsequently, ridge regression is used to combine predictions for each trial into an overall prediction
# for a drug-indication-pair. This has the nice side effect of improving calibration (it has much the
# same role as Platt scaling, which is also sometimes used in Kaggle competitions to adjust xgboost
# predictions when logloss is of primary interest).
#
# NOTE 1: The code to save the models is currently commented out to prevent accidentally overwriting them.
# NOTE 2: Script runs 6+ hours unless you have a proper ML platform
# NOTE 3: The model is not particularly hyper-parameter tuned (mostly following a basic recipe that
#   has worked okay in some past Kaggle competitions), because the Aridhia platform was too slow to get it 
#   done quickly (and we were running out of time). 
#   I REALLY WISH WE HAD HAD A GPU AND rapids.ai AVAILABLE! This would have been so much faster & easier
#   on my own machine at home.
#   Yes, LightGBM would have been an obvious alternative, for speeding things up, but we decided to 
#   prioritize other issues in the remaining time.

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
for (fi in sort(unique(new_cv_splits$foldid))){ #1:5){
  cv_index[[fi]] = (1:length(filter(new_cv_splits, foldid==fi)$set))[filter(new_cv_splits, foldid==fi)$set=="val"]
}

# Define function that fits xgboost, does predictions and calculates actual target loss 
# based on taking mean or min of predicted trial-level probabilities (we will average more cleverly, but this is a good proxy)
call_xgbcv <- function(){
  # Fit xgboost with custom cross-valdiation
  xc1 <- xgb.cv( data = data.matrix(traindat),
            label=adat$outcome,
            weight = adat$casewgt2, # Case wgts that down weight 2008 to 2011 proportional to how much the success rate is too high compared to 2012+
            params = params,
            nrounds = 10000, 
            folds=cv_index, 
            metrics = "logloss",
            showsd = T, 
            print_every_n = 100, 
            early_stopping_rounds = 20, 
            maximize = F,
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
              loglossmin = -( log(minpred)*outcome + log(1-minpred)*(1-outcome) ) * n,
              loglossmean = -( log(meanpred)*outcome + log(1-meanpred)*(1-outcome) ) * n ) %>%
    ungroup() %>%
    group_by(fold) %>%
    summarize(loglossmin = sum(loglossmin) / sum(n),
              loglossmean = sum(loglossmean) / sum(n)) %>%
    ungroup()
  
  # Summarize across folds (note that SD would not make so much sense, because of differing fold sizes)
  cv_results %>%
    summarize(loglossmin=mean(loglossmin),
              loglossmean=mean(loglossmean)) %>%
    mutate(fold=0) %>%
    bind_rows(cv_results) %>%
    bind_rows(cv_results %>%
                summarize(loglossmin=mean(ifelse(fold>=22, loglossmin, NA_real_), na.rm=T),
                          loglossmean=mean(ifelse(fold>=22, loglossmean, NA_real_), na.rm=T)) %>%
                mutate(fold=-1))
  
}

######################################## Start actual xgboost training here ######################################

# Strategy: Follow basic approach for creating a sensible, if not totally optimal xgboost model
# by choosing sensible defaults and then optimizing first the parameters that depend less on others: first max_depth
# We do this with a relatively high learning rate (eta=0.1)
experiments = list()
experimentid = 0

params = list(booster = "gbtree", 
              objective = "binary:logistic", 
              eta=0.1, 
              gamma=0, 
              max_depth=10, 
              min_child_weight=4, 
              subsample=1.0, 
              colsample_bytree=0.3)

for (current_depth in 3:11){
  print(paste("current depth ", current_depth))
  experimentid = experimentid + 1
  params$max_depth = current_depth
  experiments[[experimentid]] = call_xgbcv()
  }

res1 <- tibble(current_depth = 3:11, 
               exid = 1:experimentid) %>%
  mutate(res = map(exid, function(x) experiments[[x]])) %>%
  unnest(res)

params$max_depth = filter(res1, fold==0)$current_depth[ which.min(filter(res1, fold==0)$loglossmean) ] #6


#### Now optimize subsample
experiments = list()
experimentid = 0

for (current_subsample in seq(0.1,1,0.05)){
  experimentid = experimentid + 1
  params$subsample = current_subsample
  experiments[[experimentid]] = call_xgbcv()
}

res2 <- tibble(subsample = seq(0.1,1,0.05), 
               exid = 1:experimentid) %>%
  mutate(res = map(exid, function(x) experiments[[x]])) %>%
  unnest(res)
  

params$subsample = filter(res2, fold==0)$subsample[ which.min(filter(res2, fold==0)$loglossmean) ] 

#### Now optimize min_child_weight
experiments = list()
experimentid = 0

for (current_min_child_weight in seq(0.5,10,0.5)){
  experimentid = experimentid + 1
  params$min_child_weight = current_min_child_weight
  experiments[[experimentid]] = call_xgbcv()
}

res3 <- tibble(min_child_weight = seq(0.5,10,0.5),
               exid = 1:experimentid) %>%
  mutate(res = map(exid, function(x) experiments[[x]])) %>%
  unnest(res)

  
params$min_child_weight = filter(res3, fold==0)$min_child_weight[ which.min(filter(res3, fold==0)$loglossmean) ] 
  
#### Now optimize current_colsample_bytree
experiments = list()
experimentid = 0
for (current_colsample_bytree in seq(0.05,1,0.05)){
  experimentid = experimentid + 1
  params$colsample_bytree = current_colsample_bytree
  experiments[[experimentid]] = call_xgbcv()
}

res4 <-  tibble(colsample_bytree = seq(0.05,1,0.05),
                exid = 1:experimentid) %>%
  mutate(res = map(exid, function(x) experiments[[x]])) %>%
  unnest(res)

# res4 <-  tibble(colsample_bytree = seq(0.05,0.95,0.05),
#                 exid = 1:experimentid) %>%
#   mutate(res = map(exid, function(x) experiments[[x]])) %>%
#   unnest(res)


params$colsample_bytree = filter(res4, fold==0)$colsample_bytree[ which.min(filter(res4, fold==0)$loglossmean) ] 


# MAIN CHOICE: Limit to 2007 + casewgts2
# $booster
# [1] "gbtree"
# 
# $objective
# [1] "binary:logistic"
# 
# $eta
# [1] 0.05
# 
# $gamma
# [1] 0
# 
# $max_depth
# [1] 8 # 5 could also be quite good
# 
# $min_child_weight
# [1] 4 # Could also go as high as 8
# 
# $subsample
# [1] 1 # Could also do 0.7 or so
# 
# $colsample_bytree
# [1] 0.1 # Could be as high as 0.2 or 0.25

#  Another attempt (discarded) was all Phase 2 data, no case weights:
# $booster
# [1] "gbtree"
# 
# $objective
# [1] "binary:logistic"
# 
# $eta
# [1] 0.05
# 
# $gamma
# [1] 0
# 
# $max_depth (5 may also be good)
# [1] 6
# 
# $min_child_weight (could however go as high as, say, 6 to 8)
# [1] 3
# 
# $subsample (could also do 0.5, but 0.85 seems okay)
# [1] 0.85
# 
# $colsample_bytree
# [1] 0.15

##################################################################
# Now switch to a lower learning rate (eta=0.05) to find a 
# good number of rounds to use when refitting with the whole data
##################################################################
params$eta = 0.05
fin_train = xgb.cv( data = data.matrix(traindat),
                    label=adat$outcome,
                    weight = adat$casewgt2,
                    params = params,
                    nrounds = 10000, 
                    folds=cv_index, 
                    metrics = "logloss",
                    showsd = T, 
                    print_every_n = 10, 
                    early_stopping_rounds = 20, 
                    maximize = F,
                    callbacks = list(cb.cv.predict(save_models = TRUE)))


# fin_train : best iteration 738, add 10% so 810

# Plots to check that iteration number is not so critical
fin_train$evaluation_log %>%
  filter(iter>400) %>%
  ggplot(aes(x=iter, y=test_logloss_mean)) +
  geom_point() +
  #geom_errorbar() +
  geom_line() +
  scale_y_log10()

fin_train$evaluation_log %>% as_tibble() %>% 
  filter(iter>600) %>%
  ggplot(aes(x=iter, y=train_logloss_mean)) + geom_line(col="blue") +
  geom_line(aes(y=test_logloss_mean), col="red")

##################################################
# Train a final model
##################################################
main_model = xgboost( data = data.matrix(traindat),
                      label = adat$outcome,
                      weight = adat$casewgt2,
                      params = params,
                      nrounds = 810, 
                      metrics = "logloss",
                      print_every_n = 40, 
                      maximize = F)

# write_rds(main_model, 
#           "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/xgb_trials2_7jan2020.rds", 
#           compress = "bz2")

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

#write_rds(foravging, "/home/desktop1/Documents/holzhbj1/for_xgboost_averaging_7jan2020.rds", compress = "bz2")  
#foravging <- read_rds("/home/desktop1/Documents/holzhbj1/for_xgboost_averaging_7jan2020.rds")

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
elnet3 <- cv.glmnet(y=foravging$outcome,
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

# write_rds(elnet3,
#           "/home/desktop1/Documents/holzhbj1/novartis-dsai-challenge-starter-kit/xgboost_meta_avg2_7jan2020.rds",
#           compress="bz2")

# Measure: Binomial Deviance 
# 
# Gamma Lambda Measure       SE Nonzero
# min  0.00 0.0001 0.04185 0.002340       2
# 1se  0.15 0.4233 0.04408 0.002283       2

# cl <- makePSOCKcluster(12)
# registerDoParallel(cl)
# elnet1 <- cv.glmnet(y=foravging$outcome,
#           x=as.matrix(foravging %>% dplyr::select(logn, minpred, minlogit, meanpred, meanlogit, maxpred, maxlogit, sd_pred, sd_logit)), 
#           nfolds=10,
#           foldid=foravging$avgfold,
#           family = "binomial",
#           alpha=0.5,
#           gamma=c(0, 1, 0.1),
#           type.measure = "deviance", 
#           relax = T,
#           parallel=TRUE)
# stopCluster(cl)

# # This is how one would access the validation fold predicted values, if one wanted to do anything further (e.g. stacking).
#foravging %>% mutate(cvpred=elnet3$fit.preval$`g:0.15`[,2])

# # Some further exploration of the optimization of xgboost hyper-parameters
# res5 = tibble(fold=1:length(cv_index)) %>% 
#   mutate(res=map(fold, function(x) filter(adat, !is.na(outcome))[cv_index[[x]], c("predgroup", "outcome")] %>% 
#                    bind_cols( tibble( pred = predict(fin_train$models[[x]], newdata=data.matrix(traindat[cv_index[[x]],])) ) ) )) %>% 
#   unnest(res) %>%
#   group_by(fold,predgroup) %>%
#   summarize(n=n(),
#             outcome = max(outcome),
#             minpred=min(pred),
#             meanpred=mean(pred),
#             loglossmin = -( log(minpred)*outcome + log(1-minpred)*(1-outcome) ) * n,
#             loglossmean = -( log(meanpred)*outcome + log(1-meanpred)*(1-outcome) ) * n ) %>%
#   ungroup() %>%
#   group_by(fold) %>%
#   summarize(loglossmin = sum(loglossmin) / sum(n),
#             loglossmean = sum(loglossmean) / sum(n)) %>%
#   ungroup() 
# 
# res5 <- res5 %>%
#   summarize(loglossmin=mean(loglossmin),
#             loglossmean=mean(loglossmean)) %>%
#   mutate(fold=0) %>%
#   bind_rows(res5) %>%
#   bind_rows(res5 %>%
#               summarize(loglossmin=mean(ifelse(fold>=22, loglossmin, NA_real_), na.rm=T),
#                         loglossmean=mean(ifelse(fold>=22, loglossmean, NA_real_), na.rm=T)) %>%
#               mutate(fold=-1))
# 
# print(params)
#write_csv(res5, "/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29_1.csv")
# #write_rds(list(res1, res2, res3), "/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29_res1_3.rds", compress = "bz2")
# # reslist <- read_rds("/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29_res1_3.rds")
# # res1 <- reslist[[1]]
# # res2 <- reslist[[2]]
# # res3 <- reslist[[3]]
# res1 %>% #filter(fold %in% c(-1,0)) %>%
#   ggplot(aes(x=current_depth, y=loglossmean, col=as_factor(fold))) +
#   geom_line()  + facet_wrap(~as_factor(fold), scale="free") 
# res2 %>% #filter(fold %in% c(-1,0)) %>%
#   ggplot(aes(x=subsample, y=loglossmean, col=as_factor(fold))) +
#   geom_line() + facet_wrap(~as_factor(fold), scale="free")
# res3 %>% #filter(fold %in% c(-1,0)) %>%
#   ggplot(aes(x=min_child_weight, y=loglossmean, col=as_factor(fold))) +
#   geom_line() + facet_wrap(~as_factor(fold), scale="free")
# res4 %>% #filter(fold %in% c(-1,0)) %>%
#   ggplot(aes(x=colsample_bytree, y=loglossmean, col=as_factor(fold))) +
#   geom_line() + facet_wrap(~as_factor(fold), scale="free")


# write_rds(experiments, "/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29_experi.rds", compress = "bz2")
# experiments = read_rds("/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29_experi.rds")
#write_rds(main_model, "/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29.rds", compress = "bz2")
# main_model <- read_rds("/home/desktop1/Documents/holzhbj1/xgb_cv_2019_12_29.rds")

############################################################################################
################################ Model interpretation stuff ################################
############################################################################################

# impmatrix <- xgb.importance(colnames(traindat), model = main_model)
# 
# impmatrix[ (impmatrix$Importance/max(impmatrix$Importance)>0.05), ] %>%
#   xgb.plot.importance(., rel_to_first=T, xlab="Relative importance")
# 
# library(iml)
# 
# mypredict = function(model, newdata){
#   predict(model, newdata=data.matrix(newdata))
# }
# predictor = Predictor$new(model = main_model, data=traindat, y=adat$outcome, predict=mypredict)
# 
# imp = FeatureImp$new(predictor, loss = "logLoss")
# plot(imp)
# 
# ale = FeatureEffect$new(predictor, feature = "rel_ph2_size_ta")
# ale$plot()  + theme_bw(base_size=18) + geom_hline(yintercept=0) + ylab("Accumulated local effects")  + ylab("Relative size (for TA) of Ph2")
# ale = FeatureEffect$new(predictor, feature = "rel_ph2_size_dis")
# ale$plot() + theme_bw(base_size=18) + geom_hline(yintercept=0)  + ylab("Accumulated local effects") + ylab("Relative size (for disease type) of Ph2")
# 
# tree = TreeSurrogate$new(predictor, maxdepth = 2)
# plot(tree)
# #predict(tree, newdata = tibble(rel_ph2_size_ta=1, dtmeanclu50 = 0, dtmeancll50=0))
# 
# 
# 
# # serelaxin: 5300 5301 5302 5303 5304
# shapley = Shapley$new(predictor, x.interest = traindat[5300,])
# shapley$plot()
# 
# shapley$results %>%
#   filter(abs(phi)>1e-02) %>%
#   arrange(desc(phi)) %>%
#   mutate(color = 1*(phi>0)) %>%
#   ggplot(aes(x=phi, y=feature, col=as_factor(color))) +
#   geom_vline(xintercept=0) +
#   geom_point() + theme_bw(base_size=18) + theme(legend.position="none") 
# 
# # 
# # shapley = Shapley$new(predictor, x.interest = X[1,])
# # shapley$plot()


