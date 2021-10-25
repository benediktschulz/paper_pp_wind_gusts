## File containing functions for postprocessing ensemble forecasts via MBM

# Comment: Incorporation of sub-ensemble structure
# Estimation: CRPS

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Prediction ####
# Function for prediction based on the MBM parameters #
mbm_pred <- function(X, par_mbm = NULL, 
                     n_sens = 4, n_ens = 20, 
                     scores_pp = TRUE, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #X...............Ensemble data for prediction including variables 'sens_mean', 'sens_spread' and ens (n x (min.) (2*n_sens + n_ens) df)
  #par_mbm.........MBM parameters a, b_i, c and d_i (b and d for each subsensemble) (2 + 2*n_sens vector)
  #................Default: a = d_i = 0, b_i = c = 1, i.e. no adjustments for all subensembles
  #n_sens..........Number of subensembles (integer)
  #................Default: 4 subensemble (COSMO)
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #scores_ens/pp...Should scores of ensemble and MBM forecasts, interval lengths, ranks be calculated? (logical)
  #................Data needs to include variables 'ens_1',...,'ens_(n_ens)', 'obs'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............MBM forecasts (i.e. ensemble) based on par_mbm (n x n_ens matrix)
  #......runtime.........Prediction time (numeric)
  #......n_test..........Number of test samples (integer)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank.........Ranks of ensemble/MBM forecasts (n vector)
  #.........crps.........CRPS of ensemble/MBM forecasts (n vector)
  #.........logs.........Log-Score of ensemble/MBM forecasts (n vector)
  #.........lgt..........Ensemble/MBM range (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Relevant variables for prediction
  test_vars <- c(paste0("ens_", 1:n_ens),
                 paste0("sens_mean_", 1:n_sens),
                 paste0("sens_spread_", 1:n_sens))
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(test_vars, names(X)))){ 
    print("MBM-pred: Data does not include all of the relevant variables.") }
  
  # Cut data to relevant variables
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(X))){
    print("MBM-pred: Data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(any(X[,paste0("sens_spread_", 1:n_sens)] <= 0)){ 
    print("MBM-pred: At least one subensemble spread is non-positive!") }
  
  # Default starting parameters
  if(is.null(par_mbm)){ par_mbm <- c(0, rep(1, n_sens), 1, rep(0, n_sens)) }
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(X)
  
  # Read out MBM parameters
  a <- par_mbm[1]
  b <- par_mbm[2 + (0:(n_sens-1))]
  c <- par_mbm[3 + (n_sens-1)]
  d <- par_mbm[4 + (n_sens-1) + (0:(n_sens-1))]
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Initiate matrix
  f <- matrix(nrow = n,
              ncol = n_ens)
  
  # Loop over sub-ensembles
  for(i_sens in 1:n_sens){
    # Indices of subensemble in ensemble
    i_ens <- 1:(n_ens/n_sens) + (i_sens - 1)*(n_ens/n_sens)
    
    # Calculate MBM forecasts (vectors operate rowwise on matrix!)
    f[,i_ens] <- as.matrix((a + b[i_sens]*X[[paste0("sens_mean_", i_sens)]]) + 
                             (c + d[i_sens]/X[[paste0("sens_spread_", i_sens)]])*
                             (X[,paste0("ens_", i_ens)] - X[[paste0("sens_mean_", i_sens)]]))
  }
  
  # Make sure calibrated ensemble is non-negative
  f[f < 0] <- 0
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Evaluation ####
  # Calculate evaluation measure of MBM forecasts
  scores_pp <- fn_scores_ens(ens = f,
                             y = X[["obs"]],
                             scores_ens = scores_pp)
  
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
# Function for estimating the MBM parameters #
mbm_est <- function(train, par_start = NULL,
                    n_sens = 4, n_ens = 20){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including variables 'obs', 'sens_mean', 'sens_spread' and ens (n x (min.) (1 + 2*n_sens + n_ens) df)
  #par_start.......Initial values of optimization (2 + 2*n_sens vector)
  #................Default: NULL -> a = d_i = 0, b_i = c = 1, i.e. no adjustments
  #n_sens..........Number of subensembles (integer)
  #................Default: 4 subensemble (COSMO)
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......par.............Estimated MBM subensemble parameters a, b_i, c and d_i (2 + 2*n_sens vector)
  #......n_train.........Number of training samples (integer)
  #......runtime.........Estimation time (numeric)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Relevant variables for training
  train_vars <- c("obs", paste0("ens_", 1:n_ens),
                  paste0("sens_mean_", 1:n_sens),
                  paste0("sens_spread_", 1:n_sens))
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("MBM-est: Training data does not include relevant variables.") }
  
  # Cut data to relevant variables
  train <- train[,train_vars]
  
  # Input check
  if(any(is.na(train))){
    print("MBM-est: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(train[,paste0("sens_spread_", 1:n_sens)] <= 0)){ 
    print("MBM-est: At least one subensemble spread is non-positive!") }
  
  # Threshold for Nelder-Mead improvement
  t_nelder <- 1e-3
  
  # Define loss function for optimization
  fn_sr <- function(y, ens, ens_spread){ mean( rowMeans( abs(ens - y) ) - ens_spread/2 ) }

  #### Data preparation ####
  # None
  
  #### Estimation ####
  # Set initial values (a = d_i = 0, b_i = c = 1, i.e. no adjustments)
  if(is.null(par_start)){ par_start <- c(0, rep(1, n_sens), 1, rep(0, n_sens)) }
  
  # Define wrapper function
  wrapper <- function(par_mbm){
    ###-----------------------------------------------------------------------------
    ###Input
    #par_mbm.....MBM parameters a, b_i, c and d_i (2 + 2*n_sens vector)
    ###-----------------------------------------------------------------------------
    ###Output
    #res...Mean score of MBM forecasts (of a, b, c, d) on training set (scalar)
    ###-----------------------------------------------------------------------------
    
    #### Initialization ####
    # Read out MBM parameters
    a <- par_mbm[1]
    b <- par_mbm[2 + (0:(n_sens-1))]
    c <- par_mbm[3 + (n_sens-1)]
    d <- par_mbm[4 + (n_sens-1) + (0:(n_sens-1))]
    
    #### Calculation ####
    # Initialize ensemble matrix
    f <- matrix(nrow = nrow(train),
                ncol = n_ens)
    
    # For-Loop over subensembles
    for(i_sens in 1:n_sens){
      # Indices of subensemble in ensemble
      i_ens <- 1:(n_ens/n_sens) + (i_sens - 1)*(n_ens/n_sens)
      
      # Calculate MBM forecasts based on parameters
      f[,i_ens] <- as.matrix((a + b[i_sens]*train[[paste0("sens_mean_", i_sens)]]) + 
                               (c + d[i_sens]/train[[paste0("sens_spread_", i_sens)]])*
                               (train[,paste0("ens_", i_ens)] - train[[paste0("sens_mean_", i_sens)]]))
    }
    
    # Make sure calibrated ensemble is non-negative
    f[f < 0] <- 0
    
    # Calculate ensemble loss (cut spread for MLE)
    res <- fn_sr(y = train[["obs"]],
                 ens = f,
                 ens_spread = apply(f, 1, fn_spread))
    
    # Output
    var_check(y = res, 
              check_type = "finite", 
              name_function = "wrapper in mbm_est")
    var_check(y = res, 
              name_function = "wrapper in mbm_est")
    return(res)
  }
  
  # Take time
  start_tm <- Sys.time()
  
  # Optimize ("L-BFGS-B" as for EMOS; here no indications of failures)
  est <- optim(par = par_start, 
               fn = wrapper, 
               method = "L-BFGS-B")
  
  # Check convergence
  if(est$convergence != 0){
    # Optimize
    temp <- optim(par = par_start,
                  fn = wrapper,
                  method = "Nelder-Mead")
    
    # Check if better fit
    if(temp$value < est$value - t_nelder){ est <- temp }
  }
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Output ####
  return(list(par = est$par,
              n_train = nrow(train),
              runtime = runtime))
}

