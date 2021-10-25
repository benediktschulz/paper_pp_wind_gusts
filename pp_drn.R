## File containing functions for postprocessing ensemble forecasts via DRN

# Distribution: (zero-)truncated logistic distribution
# Estimation: CRPS
# Network combination: Parameter averaging
# Comment: Station embedding

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Estimation and Prediction ####
# Function for pp including estimation and prediction
drn_pp <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                   pred_vars = c("ens_mean", "ens_sd", "location"), 
                   nn_ls = list(), n_ens = 20, n_cores = NULL,
                   scores_ens = TRUE, scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (n_train x n_preds (+ 1) data.frame)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for NN (vector of strings)
  #................Default: c("ens_mean", "ens_sd", "location") -> Use only mean and variance (with embed.)
  #nn_ls...........List that may contain the following variables:
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...emb_dim......Dimension of station embedding (scalar)
  #................Default: 10
  #...lay1.........Number of nodes in first and second (use half) hidden layer (integer)
  #................Default: 64 -> 32 in second
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #scores_ens/pp...Should CRPS and Log-Score of ensemble and distributional forecasts be calculated (logical)
  #................Training data needs to include variables 'ens_1',...,'ens_(n_ens)'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............Distributional forecasts (i.e. parameters) (n x 2 matrix)
  #......nn_ls...........Hyperparameters (list)
  #......pred_vars.......Predictors (string vector)
  #......n_preds.........Number of predictors (integer)
  #......n_train.........Number of training samples (integer)
  #......n_valid.........Number of validation samples (integer)
  #......n_test..........Number of test samples (integer)
  #......runtime_est.....Estimation time (numeric)
  #......runtime_pred....Prediction time (numeric)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank/pit.....Ranks of ensemble / PIT of DRN forecasts (n vector)
  #.........crps.........CRPS of ensemble/DRN forecasts (n vector)
  #.........logs.........Log-Score of ensemble/DRN forecasts (n vector)
  #.........lgt..........Ensemble range / Length of DRN prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Relevant variables for training and testing
  train_vars <- unique(c("obs", "location", pred_vars))
  test_vars <- unique(c("location", pred_vars))
  
  # Input check
  if(scores_ens){
    # Ensemble members
    ens_str <- paste0("ens_", 1:n_ens)
    
    # Check
    if(any(!is.element(ens_str, names(X)))){
      print(paste0("DRN: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("DRN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("DRN: Test data does not include all of the relevant variables.") }
  
  # Cut training data to relevant variables
  train <- train[,train_vars]
  
  # Cut test data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(train))){
    print("DRN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("DRN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(is.element("ens_sd", pred_vars)){
    if(any(train[["ens_sd"]] < 0)){ print("DRN: At least one ensemble sd in training data is negative!") }
    if(any(X[["ens_sd"]] < 0)){ print("DRN: At least one ensemble sd in test data is negative!") }
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(n_sim = 10,
                  lr_adam = 5e-4, # previously 1e-3
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  emb_dim = 10,
                  lay1 = 64,
                  actv = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" }
  else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  # Get tensorflow probability
  tfp <- tf_probability()$distributions
  
  # Get standard logistic distribution
  tfd <- tfp$Logistic(loc = 0, scale = 1)
  
  # Custom loss function
  custom_loss <- function(y_true, y_pred){
    # Get location and scale
    mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
    sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
    
    # Truncated in zero
    z <- k_maximum(0, y_true)
    
    # Standardization
    z_0 <- -mu/sigma # q ~= -18.4 -> p = 1e-8
    z_y <- (z - mu)/sigma
    
    # Calculate CDFs
    p_0 <- tfd$cdf(z_0) # F(z_0)
    lp_0 <- tfd$log_cdf(z_0) # log( F(z_0) )
    p_m0 <- tfd$cdf(-z_0) # 1 - F(z_0) = F(-z_0)
    lp_m0 <- tfd$log_cdf(-z_0) # log( 1 - F(z_0)) = log( F(-z_0) )
    lp_my <- tfd$log_cdf(-z_y) # log( 1 - F(z_y)) = log( F(-z_y) )
    
    # Calculate sigma
    b <- lp_m0 - (1 + 2*lp_my)/p_m0 - k_square(p_0)*lp_0/k_square(p_m0)
    
    # Calculate CRPS
    res <- k_abs(z - y_true) - (z - mu)*(1 + p_0)/p_m0 + sigma*b
    
    # Calculate mean
    res <- k_mean(res)
    
    # Return mean CRPS
    return(res)
  }
  
  #### Data preparation ####
  # Divide data in training and validation set
  if(is.null(i_valid)){
    # Set ratio of validation data to 0.2
    r_valid <- 0.2
    
    # Take first n_train*r_valid samples for training
    i_train <- 1:floor(nrow(train)*(1 - r_valid))
    i_valid <- (max(i_train) + 1):nrow(train)
  }
  else{ i_train <- (1:nrow(train))[-i_valid] }
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(X)
  
  # Remove constant predictors
  pred_vars <- rm_const(data = train[i_train,],
                        cols = pred_vars,
                        t_c = t_c)
  
  # Predictors without location
  dir_preds <- pred_vars[pred_vars != "location"]
  
  # Get number of direct predictors
  n_dir_preds <- length(dir_preds)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train, dir_preds])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid, dir_preds],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(dir_input = as.matrix(X_train), 
                  id_input = sapply(train[i_train, "location"], function(x) which(x == loc_id_vec)))
  X_valid <- list(as.matrix(X_valid), 
                  sapply(train[i_valid, "location"], function(x) which(x == loc_id_vec)))
  
  # Get station IDs (for prediction)
  x_id <- sapply(X[,"location"], function(x) which(x == loc_id_vec))
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds], 
                                 center = tr_center, 
                                 scale = tr_scale)), x_id)
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate distribution parameter matrix
  f <- matrix(data = 0,
              nrow = nrow(X),
              ncol = 2)
  
  # For-Loop over ensemble size
  for(i_sim in 1:nn_ls$n_sim){
    #### Build network ####
    # Input
    id_input <- layer_input(shape = 1, name = "id_input")
    dir_input <- layer_input(shape = n_dir_preds, name = "dir_input")
    
    # Embedding part (see help for input_dim)
    station_embedding_part <- id_input %>%
      layer_embedding(input_dim = length(unique(train[i_train, "location"])) + 1, 
                      output_dim = nn_ls$emb_dim, input_length = 1) %>%
      layer_flatten()
    
    # Hidden layers
    hidden <- layer_concatenate(c(dir_input, station_embedding_part)) %>%
      layer_dense(units = nn_ls$lay1, activation = nn_ls$actv) %>%
      layer_dense(units = nn_ls$lay1/2, activation = nn_ls$actv)
    
    # Different activation functions for output (location positive for numerical reasons)
    loc_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
    scale_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
    
    # Concatenate output
    output <- layer_concatenate(c(loc_out, scale_out))
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = custom_loss
    )
    
    # Take time
    start_tm <- Sys.time()
    
    # Fit model
    history <- model %>% fit(
      x = X_train,
      y = train[i_train, "obs"],
      epochs = nn_ls$n_epochs,
      batch_size = nn_ls$n_batch,
      validation_data = list(X_valid, train[i_valid, "obs"]),
      verbose = nn_ls$nn_verbose,
      callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                          restore_best_weights = TRUE,
                                          monitor = "val_loss")
    )

    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- runtime_est + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Delete history
    rm(history)
    
    #### Prediction ####
    # Take time
    start_tm <- Sys.time()
    
    # Predict parameters of distributional forecasts (on scaled data)
    f <- f + predict(model, X_pred)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part, hidden, 
       loc_out, scale_out, output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average distribution parameters
  f <- f/nn_ls$n_sim
  
  #### Evaluation ####
  # Calculate evaluation measures of DRN forecasts
  if(scores_pp){ scores_pp <- fn_scores_distr(f = f,
                                              y = X[["obs"]]) }
  
  # Calculate evaluation measure of ensemble forecasts
  scores_ens <- fn_scores_ens(ens = as.matrix(X[,paste0("ens_", 1:n_ens)]),
                              y = X[["obs"]],
                              scores_ens = scores_ens)
  
  #### Output ####
  # Output
  return(list(f = f,
              nn_ls = nn_ls,
              scores_pp = scores_pp, 
              scores_ens = scores_ens,
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}
