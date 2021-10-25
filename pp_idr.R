## File containing functions for postprocessing ensemble forecasts via IDR

# Comment: Includes subagging

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Estimation and Prediction ####
# Function for estimating IDR #
idr_pp <- function(train, X, pred_vars = NULL, q_levels = NULL, 
                   idr_ls = list(), 
                   n_ens = 20, scores_ens = TRUE, scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data (-"-)
  #pred_vars.......Predictors used for IDR (vector of strings)
  #................Default: NULL -> Use ensemble member
  #q_levels........Quantile levels used for output and evaluation (probability vector)
  #................Default: NULL -> At least 100 member, incl. median and COSMO coverage
  #idr_ls..........List that may contain the following variables:
  #...groups.......Groups of predictors corresponding to a partial order (named vector)
  #................Default: NULL -> All predictors in one group
  #...orders.......Define partial orders of given groups (named vector)
  #................Default: NULL -> All increasing convex order (exchangeable ensemble)
  #...n_sbg........Number of subsamples for subagging (integer)
  #................Default: NULL -> No subagging
  #...n_sample.....Size of subsamples in subagging (integer)
  #................Default: NULL -> floor(n_train/n_sbg)
  #...max_iter.....OSQP parameter: Number of maximum iterations (integer)
  #................Default: 1000L
  #...eps_abs......OSQP parameter: Absolute convergence tolerance
  #................Default: 1e-3
  #...eps_rel......OSQP parameter: Relative convergence tolerance
  #................Default: 1e-3
  #...console......Show output of IDR function (logical)
  #................Default: FALSE (Use FALSE and not 0 !!)
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #scores_ens/pp...Should CRPS and Log-Score of ensemble and QRF forecasts be calculated (logical)
  #................Training data needs to include variables 'ens_1',...,'ens_(n_ens)'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............IDR forecasts (i.e. quantiles) based on qrf_train (n x n_q matrix)
  #......idr_ls..........Hyperparameters (list)
  #......runtime_est.....Estimation time (no sbg), incl. prediction (sbg) (numeric)
  #......runtime_pred....Prediction time (no sbg) (numeric)
  #......n_train.........Number of training samples (integer)
  #......n_test..........Number of test samples (integer)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank/pit.....Ranks of ensemble / PIT of IDR forecasts (n vector)
  #.........crps.........CRPS of ensemble/IDR forecasts (n vector)
  #.........logs.........Log-Score of ensemble/IDR forecasts (n vector)
  #.........lgt..........Ensemble range / Length of IDR prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(isodistrreg)
  library(scoringRules)
  
  # Define predictors if not given
  if(is.null(pred_vars)){ pred_vars <- paste0("ens_", 1:n_ens) }
  
  # Relevant variables for training and testing
  train_vars <- c("obs", pred_vars)
  test_vars <- pred_vars
  
  # Input check
  if(scores_ens){
    # Ensemble members
    ens_str <- paste0("ens_", 1:n_ens)
    
    # Check
    if(any(!is.element(ens_str, names(X)))){
      print(paste0("IDR: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("IDR: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("IDR: Test data does not include all of the relevant variables.") }
  
  # Cut training data to relevant variables
  train <- train[,train_vars]
  
  # Cut test data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(train))){
    print("IDR: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("IDR: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/(6*(n_ens + 1)),
                                         to = (6*(n_ens + 1) - 1)/(6*(n_ens + 1)),
                                         by = 1/(6*(n_ens + 1))) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(console = FALSE,
                  groups = rep(1, length(pred_vars)),
                  orders = c("sd" = 1),
                  n_sbg = 100,
                  n_sample = min(1000, floor(nrow(train)/2)),
                  max_iter = 1000L,
                  eps_abs = 1e-3,
                  eps_rel = 1e-3)
  
  # Update hyperparameters
  idr_ls <- update_hpar(hpar_ls = hpar_ls,
                        in_ls = idr_ls)
  
  # Set names for groups
  idr_ls$groups <- setNames(idr_ls$groups, pred_vars)

  # Define optimization settings
  pars_osqp <- list(verbose = FALSE,
                    eps_abs = idr_ls$eps_abs,
                    eps_rel = idr_ls$eps_rel,
                    max_iter = idr_ls$max_iter)
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(X)
  
  #### Estimation/Prediction ####
  # Differentiate subagging
  if(idr_ls$n_sbg == 0){ 
    # Take time
    start_tm <- Sys.time()
    
    # Estimate IDR fit
    est <- idr(y = train[["obs"]],
               X = train[,pred_vars],
               groups = idr_ls$groups,
               orders = idr_ls$orders,
               progress = idr_ls$console,
               pars = pars_osqp)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Take time
    start_tm <- Sys.time()
    
    # Calculate prediction object
    pred_idr <- predict(object = est,
                        data = X)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  }
  else{ 
    # Take time
    start_tm <- Sys.time()
    
    # Estimate and predict via subagging
    pred_idr <- idrbag(y = train[["obs"]],
                       X = train[,pred_vars],
                       groups = idr_ls$groups, 
                       orders = idr_ls$orders,
                       pars = pars_osqp, 
                       progress = idr_ls$console, 
                       newdata = X[,pred_vars], 
                       b = idr_ls$n_sbg, 
                       p = idr_ls$n_sample/nrow(train),
                       replace = (nrow(train) < idr_ls$n_sample*idr_ls$n_sbg)) 
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # No specific prediction runtime
    runtime_pred <- NA
  }
  
  #### Quantiles and Evaluation ####
  # Predict quantiles
  q <- qpred(predictions = pred_idr,
             quantiles = q_levels)
  
  # Initialize evaluation measure of IDR forecasts
  if(scores_pp){ 
    scores_pp <- data.frame(pit = numeric(length = n),
                            crps = numeric(length = n),
                            logs = numeric(length = n),
                            lgt = numeric(length = n),
                            e_me = numeric(length = n),
                            e_md = numeric(length = n))
    
    
    # Calculate PIT values
    scores_pp[["pit"]] <- pit(predictions = pred_idr,
                              y = X[["obs"]])
    
    
    # Calculate CRPS of IDR forecasts
    scores_pp[["crps"]] <- isodistrreg::crps(predictions = pred_idr,
                                             y = X[["obs"]])
    
    # Calculate bias of median forecast
    if((ncol(q) %% 2) == 1){ 
      scores_pp[["e_md"]] <- q[,(ncol(q) + 1)/2] - X[["obs"]] }
    else{ scores_pp[["e_md"]] <- qpred(predictions = pred_idr,
                                       quantiles = 0.5) - X[["obs"]] }
    
    # Calculate length of ~(n_ens-1)/(n_ens+1) % IDR prediction interval
    if(ncol(q) == 20){ scores_pp[["lgt"]] <- apply(t(apply(q, 1, range)), 1, diff) }
    # No COSMO ensemble size: Calculate corresponding range (via quantile fct.)
    else{ scores_pp[["lgt"]] <- apply(qpred(predictions = pred_idr,
                                            quantiles = c(1/(n_ens + 1), n_ens/(n_ens + 1))), 1, diff) }
    
    # Calculate Log-Score and bias of mean forecast as for ensemble
    temp_scores <- fn_scores_ens(ens = q,
                                 y = X[["obs"]],
                                 skip_evals = c("rank", "crps", "e_md", "lgt"),
                                 scores_ens = TRUE)
    
    # Read out Log-Score and bias of mean forecast
    scores_pp[["logs"]] <- temp_scores[["logs"]]
    scores_pp[["e_me"]] <- temp_scores[["e_me"]]
  }
  
  # Calculate evaluation measure of ensemble forecasts
  scores_ens <- fn_scores_ens(ens = as.matrix(X[,paste0("ens_", 1:n_ens)]),
                              y = X[["obs"]],
                              scores_ens = scores_ens)
  
  #### Output ####
  return(list(f = q, 
              pred_idr = pred_idr,
              idr_ls = idr_ls,
              scores_pp = scores_pp, 
              scores_ens = scores_ens,
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = nrow(train),
              n_test = nrow(X),
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}

