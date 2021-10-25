## File containing functions for postprocessing ensemble forecasts via EMOS-GB

# Distribution: (zero-)truncated logistic distribution
# Estimation: MLE

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Prediction ####
# Function for prediction based on the EMOS parameters obtained via boosting #
emos_bst_pred <- function(X, emos_bst_train, pred_vars, n_ens = 20, 
                          scores_pp = TRUE, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #X................Ensemble data for prediction including predictors (and obs.) (n x n_preds (+ 1) data.frame)
  #emos_bst_train...Output of crch function used for prediction
  #pred_vars........Predictors used for EMOS-bst (vector of strings)
  #n_ens............Ensemble size (integer)
  #.................Default: 20 member (COSMO)
  #scores_ens/pp....Should scores of ensemble and Emos (boosting) forecasts, interval lengths, ranks be calculated? (logical)
  #.................Data needs to include variables 'ens_1',...,'ens_(n_ens)', 'obs'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............EMOS forecasts (i.e. location and scale) based on emos_bst_train (n x 2 matrix)
  #......runtime.........Prediction time (numeric)
  #......n_test..........Number of test samples (integer)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank/pit.....Ranks of ensemble forecasts / PIT values of EMOS-bst (n vector)
  #.........crps.........CRPS of ensemble/EMOS-bst forecasts (n vector)
  #.........logs.........Log-Score of ensemble/EMOS-bst forecasts (n vector)
  #.........lgt..........Ensemble range / Length of EMOS-bst prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  library(crch)
  
  # Relevant variables for prediction
  test_vars <- pred_vars
  
  # Input check
  if(scores_ens){
    # Ensemble members
    ens_str <- paste0("ens_", 1:n_ens)
    
    # Check
    if(any(!is.element(ens_str, names(X)))){
      print(paste0("EMOS-bst-pred: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(test_vars, names(X)))){ 
    print("EMOS-bst-pred: Data does not include all of the relevant variables.") }
  
  # Cut data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(X))){
    print("EMOS-bst-pred: Data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(is.element("ens_sd", test_vars)){ if(any(X[["ens_sd"]] < 0)){ 
    print("EMOS-bst-pred: At least one ensemble standard deviation is negative!") }}
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(X)
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Calculate forecasts
  f <- predict(object = emos_bst_train,
               newdata = X,
               type = "parameter")
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Evaluation ####
  # Calculate evaluation measure of EMOS forecasts
  if(scores_pp){ scores_pp <- fn_scores_distr(f = f,
                                              y = X[["obs"]]) }
  
  # Calculate evaluation measure of ensemble forecasts
  scores_ens <- fn_scores_ens(ens = as.matrix(X[,paste0("ens_", 1:n_ens)]),
                              y = X[["obs"]],
                              scores_ens = scores_ens)
  
  #### Output ####
  return(list(f = f, 
              scores_pp = scores_pp, 
              scores_ens = scores_ens,
              n_test = nrow(X),
              runtime = runtime))
}

#### Estimation ####
# Function for estimating EMOS parameters via boosting #
emos_bst_est <- function(train, pred_vars = c("ens_mean", "ens_sd"), 
                         bst_ls = list(), n_ens = 20){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predcitors and obs. (n_train x (n_preds + 1) data.frame)
  #pred_vars.......Predictors used for EMOS boosting (vector of strings)
  #................Default: c("ens_mean", "ens_sd") -> Use only mean and variance
  #bst_ls..........List that may contain the following variables:
  #...nu...........Step size (positive scalar)
  #................Default: 0.05
  #...maxit........Maximum number of iterations (integer)
  #................Default: 1,000
  #...mstop........Stopping criterion ("max", "aic", "bic", "cv")
  #................Default: AIC (Aikake Information criterion)
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......emos_bst_train...EMOS-bst estimation (crch-object)
  #......bst_ls...........Hyperparameters (list)
  #......pred_vars........Predictors (string vector)
  #......n_preds..........Number of predictors used (integer)
  #......n_train..........Number of training samples (integer)
  #......runtime..........Estimation time (numeric)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(crch)
  
  # Relevant variables for training
  train_vars <- c("obs", pred_vars)
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("EMOS-bst-est: Training data does not include relevant variables.") }
  
  # Cut data to relevant variables
  train <- train[,train_vars]
  
  # Input check
  if(any(is.na(train))){
    print("EMOS-bst-est: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(is.element("ens_sd", names(train))){ if(any(train[["ens_sd"]] < 0)){ 
    print("EMOS-bst-est: At least one ensemble standard deviation is negative!") }}
  
  # Threshold for (almost) constant predictors
  t_c <- 1e-4
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(nu = 0.05,
                  maxit = 1000,
                  mstop = "aic")
  
  # Update hyperparameters
  bst_ls <- update_hpar(hpar_ls = hpar_ls,
                        in_ls = bst_ls)
  
  #### Data preparation ####
  # Remove constant predictors
  pred_vars <- rm_const(data = train,
                        cols = pred_vars,
                        t_c = t_c)
  
  #### Estimation ####
  # Use all predictors for both location and scale
  temp_fml <- paste0("obs ~ ", paste0(pred_vars, collapse = " + "), " | ", 
                     paste0(pred_vars, collapse = " + "))
  
  # Take time
  start_tm <- Sys.time()
  
  # Boosting
  est <- crch(formula = temp_fml,
              data = train,
              link.scale = "log",
              dist = "logistic",
              truncated = TRUE,
              left = 0,
              type = "ml",
              control = crch.boost(maxit = bst_ls$maxit,
                                   nu = bst_ls$nu,
                                   mstop = bst_ls$mstop))
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Output ####
  return(list(emos_bst_train = est,
              bst_ls = bst_ls,
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = nrow(train),
              runtime = runtime))
}


