# Winning solution Novartis Data Science and Artificial Intelligence Competition 2019/2020

This repository contains the code for the final solution of the winning team (Team Insight-Out) of the Novartis Data Science and Artificial Intelligence 2019/2020 competition. This was an company internal data science and artificial intelligence challenge held by Novartis. Teams were challenged to predict the [probability of success](https://arxiv.org/abs/2102.02752) (PoS) for drugs (whether a drug would be approved for a specific indication) based on the data available from [Phase 2 clinical trials](https://www.fda.gov/patients/drug-development-process/step-3-clinical-research). The public leaderboard during the competition and the final private leaderboard at the end of the competition was based on the binary log-loss metric on a set of data from after the period covered by the training data. Our team achieved the top spot in the competition with a public leaderboard log-loss of 0.196	(F1-score 0.740) and final private leaderboard log-loss of 0.202 (F1-score 0.739).

A manuscript describing the competition and our team's solution in more detail is currently under review at a scientific journal.

# Brief overview of solution

Our final solution (shown in the figure below) was based on extensive feature engineering, an ensemble of three models and post-processing of the predictions. The ensemble consisted of two [xgboost](https://xgboost.readthedocs.io/en/latest/) models that predict probabilities of approval per trial (combined into a single drug-indication prediction using [ridge regression](https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf)) and one [Bayesian logistic regression](https://avehtari.github.io/modelselection/diabetes.html) with [covariate balancing propensity](https://cran.r-project.org/web/packages/CBPS/index.html) weights. [Given the time series nature of the data, our validation approach](https://www.fast.ai/2017/11/13/validation-sets/) used several different past-vs.-future (cross-)-validation to combine trial-level xgboost predictions, optimize hyperparameters and make other modeling decisions. For further details see the forthcoming manuscript that is currently under review at a scientific journal.
![image](https://user-images.githubusercontent.com/18594459/115703334-4c7c2f80-a36a-11eb-91be-3ebb194bfd90.png)

# Will this code run on my computer without modification?

This code was intended to run in the competition environment, in which teams developed their solutions (provided by [Aridhia](https://www.aridhia.com/)), and the solution submission system (provided by [AIcrowd](https://www.aicrowd.com/)), respectively. You will have to modify it (especially with respect to file paths) to run in your environment. Additionally, the training data was based on two proprietary pharmaceutical pipeline databases provided by [Informa](https://pharmaintelligence.informa.com/)&copy; ([Pharmaprojects](https://pharmaintelligence.informa.com/products-and-services/data-and-analysis/pharmaprojects) and [Trialtrove](https://pharmaintelligence.informa.com/clinical-trial-data)), and we therefore do not  provide the raw or processed training and inference data in this repository.

# Description of files

Below, we describe the key files with data wrangling, model training, inference code and model files that are included in this repository. 

Each of the `.R` programs described below comes with additional comments within the program that describe what the program does and what assumptions are made. In particular, note that the code assumes the existence of the raw data and the created `.csv` files with features, which cannot be shared due to being based on a proprietary database. The main created dataset used for model training was `feat_bjoern_trial.csv`. While the contents cannot be shared due to the proprietary nature of the data, the features are described in the code, as well as in a manuscript (and its supplmentary materials) that is currently under review at a journal. The contents of the `.csv` file are summarized using the skimr R package in the `skim_feat_bjoern_trial.html` file.

## Files for inference

The following are the key files for obtaining inference for a set of Phase 2 trial results:
* `predict.R`: runs predictions on the evaluation server of the competition
* Trained models used by `predict.R`:
 * `xgb_trials2_7jan2020.rds`: First xgboost model using trial-level level features going into the final ensemble (can be loaded with the `xgboost` `R` package, see `predict.R`)
 * `xgboost_meta_avg2_7jan2020.rds`: Relaxed ridge for combining trial-level predictions from xgboost model #1 into project level predictions (can be loaded with the `glmnet` `R` package, see `predict.R`)
 * `xgb_trials3_10jan2020.rds`: Second xgboost model using trial-level level features going into the final ensemble (can be loaded with `xgboost` `R` package, see `predict.R`)
 * `xgboost_meta_avg3_10jan2020.rds`: Relaxed ridge for combining trial-level predictions from xgboost model #2 into project level predictions (can be loaded with `glmnet` `R` package, see `predict.R`)
 * `sglmerfit1b.rds`: fitted Bayesian logistic regression model using project level features that goes into the final ensemble (can be loaded with the `rstanarm` `R` package, see `predict.R`)

## Files for model training

There are three main model training programs, each of which trained one of the main models in our ensemble for our best public leaderboard submission:
* xgboost model 1: `xgb_cv_new_cv_split_plus_ridge.R`
* xgboost model 2 (automatic hyperparameter tuning using differential evolution): `tune_xgb_more.R`
* Bayesian logistic regression: `cov_balancing_propensity_score.R` (additional attempts around this model in the files with the additional suffic `*_cv` and `*_further_try`)

## Files for data wrangling

* `new_cv_splits.R`: Program to create the cross-validation splits that we used to tune the hyperparameters of xgboost models and to combine trial-level predictions into drug-indication-level predictions.
* `get_orphan.R`: Uses additional data (i.e. in addition to the main merged file) provided by competition organizers in order to obtain information on orphan drug status.
* `features_bjoern.R`: Feature engineering and data wrangling on a drug-indication level
* `feat_bjoern_trial.R`: Feature engineering and data wrangling at the individual trial level (produces `feat_bjoern_trial.csv`)
