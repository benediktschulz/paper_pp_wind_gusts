## File containing functions for postprocessing ensemble forecasts via QRF

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Prediction ####
# Function for prediction based on QRF #
qrf_pred <- function(X, qrf_train, pred_vars, 
                     q_levels = NULL, n_ens = 20, n_cores = NULL, 
                     scores_pp = TRUE, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #X...............Ensemble data for prediction including predictors (and obs.) (n x n_preds (+ 1) data.frame)
  #qrf_train.......Output of qrf function used for prediction
  #pred_vars.......Predictors used for QRF (vector of strings)
  #q_levels........Quantile levels used for output and evaluation (probability vector)
  #................Default: NULL -> At least 100 member, incl. median and COSMO coverage
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #n_cores.........Number of cores used in predict.ranger (integer)
  #................Default: NULL -> Use one less than available
  #scores_ens/pp...Should scores of ensemble and QRF forecasts, interval lengths, ranks be calculated? (logical)
  #................Data needs to include variables 'ens_1',...,'ens_(n_ens)', 'obs'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............QRF forecasts (i.e. quantiles) based on qrf_train (n x n_q matrix)
  #......runtime.........Prediction time (numeric)
  #......n_test..........Number of test samples (integer)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank.........Ranks of ensemble/QRF forecasts (n vector)
  #.........crps.........CRPS of ensemble/QRF forecasts (n vector)
  #.........logs.........Log-Score of ensemble/QRF forecasts (n vector)
  #.........lgt..........Ensemble range / Length of QRF prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(ranger)
  library(scoringRules)
  
  # Relevant variables for prediction
  test_vars <- pred_vars
  
  # Input check
  if(scores_ens){
    # Ensemble members
    ens_str <- paste0("ens_", 1:n_ens)
    
    # Check
    if(any(!is.element(ens_str, names(X)))){
      print(paste0("QRF-pred: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(test_vars, names(X)))){ 
    print("QRF-pred: Data does not include all of the relevant variables.") }

  # Cut data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(X))){
    print("QRF-pred: Data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(is.element("ens_sd", test_vars)){ if(any(X[["ens_sd"]] < 0)){ 
    print("QRF-pred: At least one ensemble standard deviation is negative!") }}
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores() - 1 }
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/(6*(n_ens + 1)),
                                         to = (6*(n_ens + 1) - 1)/(6*(n_ens + 1)),
                                         by = 1/(6*(n_ens + 1))) }
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(X)
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Calculate quantiles
  q <- predict(object = qrf_train,
               data = X,
               type = "quantiles",
               quantiles = q_levels,
               num.threads = n_cores)$predictions
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Evaluation ####
  # Calculate evaluation measure of QRF forecasts
  scores_pp <- fn_scores_ens(ens = q,
                             y = X[["obs"]],
                             scores_ens = scores_pp)
  
  # Calculate evaluation measure of ensemble forecasts
  scores_ens <- fn_scores_ens(ens = as.matrix(X[,paste0("ens_", 1:n_ens)]),
                              y = X[["obs"]],
                              scores_ens = scores_ens)
  
  # Transform ranks to n_(ens + 1) bins (for multiples of (n_ens + 1) exact)
  if(ncol(q) != n_ens){ scores_pp[["rank"]] <- ceiling(scores_pp[["rank"]]*(n_ens + 1)/(ncol(q) + 1)) }
  
  #### Output ####
  return(list(f = q, 
              scores_pp = scores_pp, 
              scores_ens = scores_ens,
              n_test = nrow(X),
              runtime = runtime))
}

#### Estimation ####
# Function for estimating QRF #
qrf_est <- function(train, pred_vars = c("ens_mean", "ens_sd"), 
                    n_ens = 20, n_cores = NULL, qrf_ls = list()){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #pred_vars.......Predictors used for QRF (vector of strings)
  #................Default: c("ens_mean", "ens_sd") -> Use only mean and variance
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #n_cores.........Number of cores used in ranger (integer)
  #................Default: NULL -> Use one less than available
  #qrf_ls..........List that may contain the following variables:
  #...console......Query, if output should be shown (logical)
  #................Default: FALSE
  #...importance...Importance setting for ranger (string)
  #................Default: "permutation"
  #...n_mtry.......Number of variables considered at each split (1,...,length(pred_vars))
  #................Default: -1 -> NULL (< 10 preds) resp. half of predictors
  #...n_trees......Number of trees (integer)
  #................Default: 1,000
  #...min_node.....Minimal node size (integer)
  #................Default: 10
  #...max_depth....Maximal tree depth (integer)
  #................Default: 0 (default) -> unlimited
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......qrf_train.......Estimated QRF ('ranger'-object)
  #......qrf_ls..........Hyperparameters (list)
  #......pred_vars........Predictors (string vector)
  #......n_preds.........Number of predictors used (integer)
  #......n_train.........Number of training samples (integer)
  #......runtime.........Estimation time (numeric)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(ranger)
  
  # Relevant variables for training
  train_vars <- c("obs", pred_vars)
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("QRF-est: Training data does not include relevant variables.") }
  
  # Cut data to relevant variables
  train <- train[,train_vars]
  
  # Input check
  if(any(is.na(train))){
    print("QRF-est: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(is.element("ens_sd", train_vars)){ if(any(train[["ens_sd"]] < 0)){ 
    print("QRF-est: At least one ensemble standard deviation is negative!") }}
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores() - 1 }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values (for global QRF)
  hpar_ls <- list(console = 0,
                  importance = "permutation",
                  n_trees = 1000,
                  min_node = 10,
                  max_depth = 0,
                  n_mtry = -1) # -2 for ranger-default, -1 for own default
  
  # Update hyperparameters
  qrf_ls <- update_hpar(hpar_ls = hpar_ls,
                        in_ls = qrf_ls)
  
  # Set mtry to ranger-default (-2) or half of predictors (-1)
  if(qrf_ls$n_mtry == -1){ qrf_ls$n_mtry <- floor(length(pred_vars)/2) }
  else if(qrf_ls$n_mtry == -2){ qrf_ls$n_mtry <- NULL }
  
  #### Data preparation ####
  # Remove constant predictors
  pred_vars <- rm_const(data = train,
                        cols = pred_vars,
                        t_c = 0)
  
  #### Estimation ####
  # Get formula from predictors
  qrf_formula <- paste0("obs ~ ", paste0(pred_vars, collapse = " + "))
  
  # Take time
  start_tm <- Sys.time()
  
  # Quantile regression forest
  est <- ranger(formula = qrf_formula,
                data = train,
                num.trees = qrf_ls$n_trees,
                mtry = qrf_ls$n_mtry,
                min.node.size = qrf_ls$min_node,
                max.depth = qrf_ls$max_depth,
                quantreg = TRUE,
                num.threads = n_cores,
                importance = qrf_ls$importance,
                verbose = qrf_ls$console)
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Output ####
  return(list(qrf_train = est,
              qrf_ls = qrf_ls,
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = nrow(train),
              runtime = runtime))
}

