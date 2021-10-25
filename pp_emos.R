## File containing functions for postprocessing ensemble forecasts via EMOS

# Distribution: (zero-)truncated logistic distribution
# Estimation: CRPS

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Prediction ####
# Function for prediction based on the EMOS parameters #
emos_pred <- function(X, par_emos = NULL, t_sd = NULL,
                      n_ens = 20, scores_pp = TRUE, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #X...............Ensemble data for prediction including predictors (and obs.) (n x n_preds (+ 1) data.frame)
  #par_emos........EMOS parameters a, b, c and d (4 vector)
  #................Default: NULL -> a = c = 0, b = d = 1, i.e. ensemble mean and variance as location and scale
  #t_sd............Minimal threshold for ensemble standard deviations (positive scalar)
  #................Default: NULL -> No cut-off
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #scores_ens/pp...Should scores of ensemble and EMOS forecasts, interval lengths, PIT values be calculated? (logical)
  #................Data needs to include variables 'ens_1',...,'ens_(n_ens)', 'obs'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............EMOS forecasts (i.e. location, scale and shape) based on par_emos (n x 2/3 matrix)
  #......runtime.........Prediction time (numeric)
  #......n_test..........Number of test samples (integer)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank/pit.....Ranks of ensemble forecasts / PIT values of EMOS (n vector)
  #.........crps.........CRPS of ensemble/EMOS forecasts (n vector)
  #.........logs.........Log-Score of ensemble/EMOS forecasts (n vector)
  #.........lgt..........Ensemble range / Length of EMOS prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Relevant variables for prediction
  test_vars <- c("ens_mean", "ens_sd")
  
  # Input check
  if(scores_ens){
    # Ensemble members
    ens_str <- paste0("ens_", 1:n_ens)
    
    # Check
    if(any(!is.element(ens_str, names(X)))){
      print(paste0("EMOS-pred: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(test_vars, names(X)))){ 
    print("EMOS-pred: Data does not include all of the relevant variables.") }
  
  # Cut data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(X))){
    print("EMOS-pred: Data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(is.null(t_sd) & any(X[["ens_sd"]] <= 0)){ 
    print("EMOS-pred: At least one ensemble standard deviation is non-positive!") }
  
  # Initiate parameter vector
  if(is.null(par_emos)){ par_emos <- c(0, 0, 0, 1) }
  
  #### Data preparation ####
  # Cut ensemble standard deviations to t_sd
  if(!is.null(t_sd)){ X[["ens_sd"]] <- pmax(t_sd, X[["ens_sd"]]) }
  
  # Number of predictions
  n <- nrow(X)
  
  # Read out EMOS parameters
  a <- par_emos[1]
  b <- par_emos[2]
  c <- par_emos[3]
  d <- par_emos[4]
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Initiate forecast matrix
  f <- matrix(nrow = n, 
              ncol = 2)
  colnames(f) <- c("location", "scale")
  
  # Calculate location
  f[,1] <- a + exp(b)*X[["ens_mean"]]
  
  # Calculate scale
  f[,2] <- exp(c + d*log(X[["ens_sd"]]))
  
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
# Function for estimating the EMOS parameters #
emos_est <- function(train, n_ens = 20, par_start = NULL, t_sd = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #par_start.......Initial values of optimization (4 vector)
  #................Default: NULL -> a = c = 0, b = d = 1, i.e. ensemble mean and variance as location and scale
  #t_sd............Minimal threshold for ensemble standard deviations (positive scalar)
  #................Default: NULL -> No cut-off
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......par.............Estimated EMOS parameters a, b, c, d (4 vector)
  #......n_train.........Number of training samples (integer)
  #......runtime.........Estimation time (numeric)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Relevant variables for training
  train_vars <- c("obs", "ens_mean", "ens_sd")
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("EMOS-est: Training data does not include relevant variables.") }
  
  # Cut data to relevant variables
  train <- train[,train_vars]
  
  # Input check
  if(any(is.na(train))){
    print("EMOS-est: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(is.null(t_sd) & any(train[["ens_sd"]] <= 0)){ 
    print("EMOS-est: At least one ensemble standard deviation is non-positive!") }
  
  # Threshold for Nelder-Mead improvement
  t_nelder <- 1e-3
  
  # Minimum resp. maximum threshold for location and scale
  t_max <- 1e+3
  t_min <- 1e-3
  
  # Threshold for (almost) zero observations
  t_0 <- 1e-2
  
  # Location function
  fn_loc <- function(ens_mean, a, b){ pmax(-t_max, pmin(t_max, a + exp(b)*ens_mean)) }
  
  # Scale function
  fn_scale <- function(ens_sd, c, d){ pmax(t_min, pmin(t_max, exp(c + d*log(ens_sd)))) }
  
  # Define CRPS
  fn_sr <- function(y, location, scale){ crps_tlogis(y = y, 
                                                     location = location, 
                                                     scale = scale,
                                                     lower = 0) }
  
  # Define gradient for CRPS estimation
  fn_grad <- function(y, location, scale){ gradcrps_tlogis(y = y, 
                                                           location = location, 
                                                           scale = scale,
                                                           lower = 0) }
  
  #### Data preparation ####
  # Cut ensemble standard deviations to t_sd
  if(!is.null(t_sd)){ train[["ens_sd"]] <- pmax(t_sd, train[["ens_sd"]]) }
  
  # Set (almost) zero-observations to t_0
  obs_cut <- pmax(t_0, train[["obs"]])
  
  #### Estimation ####
  # Set initial values (a = c = 0, b/exp(b) = d = 1, i.e. ensemble mean and standard deviation as parameters)
  if(is.null(par_start)){ par_start <- c(0, 0, 0, 1) }
  
  # Define wrapper function
  wrapper <- function(par_emos){
    ###-----------------------------------------------------------------------------
    ###Input
    #par_emos.....EMOS parameters a, b, c, d (4 vector)
    ###-----------------------------------------------------------------------------
    ###Output
    #res...Mean score of EMOS forecasts (of a, b, c, d) on training set (scalar)
    ###-----------------------------------------------------------------------------
    
    #### Calculation ####
    # Calculate location and scale parameters
    loc_emos <- fn_loc(ens_mean = train[["ens_mean"]], 
                       a = par_emos[1], 
                       b = par_emos[2])
    scale_emos <- fn_scale(ens_sd = train[["ens_sd"]], 
                           c = par_emos[3], 
                           d = par_emos[4])
    
    # Calculate mean scores of training data
    res <- mean(fn_sr(y = obs_cut,
                      location = loc_emos,
                      scale = scale_emos))
    
    # Output
    var_check(y = res, 
              check_type = "finite", 
              name_function = "wrapper in emos_est")
    var_check(y = res, 
              name_function = "wrapper in emos_est")
    return(res)
  }
  
  # Define gradient (w.r.t. a, b, c, d)
  grad <- function(par_emos){
    ###-----------------------------------------------------------------------------
    ###Input
    #par_emos.....EMOS parameters a, b, c and d (4 vector)
    ###-----------------------------------------------------------------------------
    ###Output
    #res...Gradient of mean score w.r.t. EMOS parameters on training set (4 vector)
    ###-----------------------------------------------------------------------------
    
    #### Initialization ####
    # Initialize resulting gradient
    res <- vector(length = 4)
    
    #### Calculation ####
    # For each ensemble: Calculate location and scale parameters
    loc_emos <- fn_loc(ens_mean = train[["ens_mean"]], 
                       a = par_emos[1], 
                       b = par_emos[2])
    scale_emos <- fn_scale(ens_sd = train[["ens_sd"]], 
                           c = par_emos[3], 
                           d = par_emos[4])
    
    # Calculate gradient of crps
    s_grad <- fn_grad(y = obs_cut, 
                      location = loc_emos,
                      scale = scale_emos)
    
    # Derivatives w.r.t. a and b
    res[1] <- mean(s_grad[,"dloc"])
    res[2] <- mean(s_grad[,"dloc"]*exp(par_emos[2])*train[["ens_mean"]])
    
    # Derivatives w.r.t. c and d
    res[3] <- mean(s_grad[,"dscale"]*scale_emos)
    res[4] <- mean(s_grad[,"dscale"]*scale_emos*log(train[["ens_sd"]])) 
    
    # Output
    var_check(y = res, 
              name_function = "grad in emos_est")
    return(res)
  }
  
  # Take time
  start_tm <- Sys.time()
  
  # Try optimizing
  try_est <- try(expr = {
    # Optimize ("L-BFGS-B" outperforms "BFGS" (drastically) and "Nelder-Mead")
    est <- optim(par = par_start, 
                 fn = wrapper, 
                 gr = grad,
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
  },
  silent = TRUE)
  
  # If optimizing did not work, use Nelder-Mead
  if(class(try_est) == "try-error"){
    # Optimize
    est <- optim(par = par_start, 
                 fn = wrapper, 
                 gr = grad,
                 method = "Nelder-Mead")
  }
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  #### Output ####
  # Output
  return(list(par = est$par,
              n_train = nrow(train),
              runtime = runtime))
}
