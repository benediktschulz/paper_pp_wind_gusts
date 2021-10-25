## File containing functions for postprocessing ensemble forecasts via BQN

## Acknowledgments:
# We thank John B. Bremnes for supplying code of the original paper (https://doi.org/10.1175/MWR-D-19-0227.1),
# on which we built the following code.

# Network combination: Vincentization
# Comment: Station embedding

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_data.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Estimation and prediction ####
# Function for pp including estimation and prediction
bqn_pp <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                   pred_vars = c("ens_mean", "ens_sd", "location"), 
                   q_levels = NULL, nn_ls = list(), n_ens = 20, n_cores = NULL, 
                   scores_ens = TRUE, scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (n_train x n_preds (+ 1) data.frame)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for BQN (vector of strings)
  #................Default: c("ens_mean", "ens_sd", "location") -> Use only mean and variance (with embed.)
  #q_levels........Quantile levels used for output and evaluation (n_q probability vector)
  #................Default: NULL -> At least 100 member, incl. median and COSMO coverage
  #nn_ls...........List that may contain the following variables:
  #...p_degree.....Degree of Bernstein polynomials (integer)
  #................Default: 12
  #...n_q..........Number of equidistant quantile levels used in loss function (integer)
  #................Default: 99 (steps of 1%)
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience for early stopping (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...emb_dim......Dimension of station embedding (integer)
  #................Default: 10
  #...lay1.........Number of nodes in first hidden layer (integer)
  #................Default: 48
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...actv_out.....Activation function of output layer (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_ens...........Ensemble size (integer)
  #................Default: 20 member (COSMO)
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #scores_ens/pp...Should CRPS and Log-Score of ensemble and BQN forecasts be calculated (logical)
  #................Training data needs to include variables 'ens_1',...,'ens_(n_ens)'
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............BQN forecasts (i.e. quantiles) based on q_levels (n x n_q matrix)
  #......alpha...........BQN coefficients (n x p_degree matrix)
  #......nn_ls...........Hyperparameters (list)
  #......pred_vars.......Predictors (string vector)
  #......n_preds.........Number of predictors (integer)
  #......n_train.........Number of training samples (integer)
  #......n_valid.........Number of validation samples (integer)
  #......n_test..........Number of test samples (integer)
  #......runtime_est.....Estimation time (numeric)
  #......runtime_pred....Prediction time (numeric)
  #......scores_ens/pp...Data frames containing (n x 6 data frame):
  #.........rank.........Ranks of observations in ensemble/BQN forecasts (n vector)
  #.........crps.........CRPS of ensemble/BQN forecasts (n vector)
  #.........logs.........Log-Score of BQN/ensemble forecasts (n vector)
  #.........lgt..........Length of BQN prediction interval/Ensemble range (n vector)
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
      print(paste0("BQN: Data does not include all ensemble members!
            CRPS, Log-Score of raw ensemble and ranks can therefore not be calculated!"))
      scores_ens <- FALSE
    }
  }
  
  # Observations for scores
  if(scores_pp | scores_ens){ test_vars <- c(test_vars, "obs") }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("BQN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("BQN: Test data does not include all of the relevant variables.") }
  
  # Cut training data to relevant variables
  train <- train[,train_vars]
  
  # Cut test data to relevant variables
  if(scores_ens){ X <- X[,unique(c(test_vars, ens_str))] }
  else{ X <- X[,test_vars] }
  
  # Input check
  if(any(is.na(train))){
    print("BQN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("BQN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  if(is.element("ens_sd", pred_vars)){
    if(any(train[["ens_sd"]] < 0)){ print("BQN: At least one ensemble sd in training data is negative!") }
    if(any(X[["ens_sd"]] < 0)){ print("BQN: At least one ensemble sd in test data is negative!") }
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/(6*(n_ens + 1)),
                                         to = (6*(n_ens + 1) - 1)/(6*(n_ens + 1)),
                                         by = 1/(6*(n_ens + 1))) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(p_degree = 12,
                  n_q = 99,
                  n_sim = 10,
                  lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  emb_dim = 10,
                  lay1 = 48,
                  actv = "softplus",
                  actv_out = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Calculate equidistant quantile levels for loss function
  q_levels_loss <- seq(from = 1/(nn_ls$n_q + 1),
                       to = nn_ls$n_q/(nn_ls$n_q + 1),
                       by = 1/(nn_ls$n_q + 1))
  
  # Basis of Bernstein polynomials evaluated at quantile levels
  B <- sapply(0:nn_ls$p_degree, 
              dbinom, size = nn_ls$p_degree, prob = q_levels_loss)
  
  # Quantile loss functions (for neural network)
  qt_loss <- function(y_true, y_pred){
    # Quantiles calculated via basis and increments
    q  <- k_dot(k_cumsum(y_pred, axis = 0),
                k_constant(as.numeric(B), shape = c(nn_ls$p_degree + 1, nn_ls$n_q)))
    
    # Calculate individual quantile scores
    err  <- y_true - q
    e1   <- err * k_constant(q_levels_loss, shape = c(1, nn_ls$n_q))
    e2   <- err * k_constant(q_levels_loss - 1, shape = c(1, nn_ls$n_q))
    
    # Find correct values (max) and return mean
    return(k_mean( k_maximum(e1, e2), axis = 2 ))
  }
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" }
  else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
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
  
  # Get station IDs
  x_id <- sapply(X[,"location"], function(x) which(x == loc_id_vec))
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds],
                                 center = tr_center,
                                 scale = tr_scale)), x_id)
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate coefficient matrix
  coeff_bern <- matrix(data = 0,
                       nrow = nrow(X),
                       ncol = nn_ls$p_degree + 1)
  
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
    
    # Monotonicity in output (increments; strict via softplus, else relu)
    output <- hidden %>%
      layer_dense(units = nn_ls$p_degree + 1, activation = nn_ls$actv_out)
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = qt_loss
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
    
    # Predict coefficients of Bernstein polynomials
    coeff_bern <- coeff_bern + predict(model, X_pred)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part,
       hidden, output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average increments
  coeff_bern <- coeff_bern/nn_ls$n_sim
  
  # Accumulate increments
  coeff_bern <- t(apply(coeff_bern, 1, cumsum))
  
  #### Evaluation ####
  # Sum up calculated quantiles (Sum of basis at quantiles times coefficients)
  q <- bern_quants(alpha = coeff_bern,
                   q_levels = q_levels)
  
  # Calculate evaluation measure of BQN forecasts
  scores_pp <- fn_scores_ens(ens = q,
                             y = X[["obs"]],
                             skip_evals = c("e_me"),
                             scores_ens = scores_pp)
  
  # Transform ranks to n_(ens + 1) bins (for multiples of (n_ens + 1) exact)
  if(ncol(q) != n_ens){ scores_pp[["rank"]] <- ceiling(scores_pp[["rank"]]*(n_ens + 1)/(ncol(q) + 1)) }
  
  # Calculate bias of mean forecast (formula given)
  scores_pp[["e_me"]] <- rowMeans(coeff_bern) - X[["obs"]]
  
  # Calculate evaluation measure of ensemble forecasts
  scores_ens <- fn_scores_ens(ens = as.matrix(X[,paste0("ens_", 1:n_ens)]),
                              y = X[["obs"]],
                              scores_ens = scores_ens)
  
  #### Output ####
  return(list(f = q, 
              alpha = coeff_bern,
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

