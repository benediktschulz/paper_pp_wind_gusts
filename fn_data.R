## Functions for to work with data

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_basic.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### Remove constant columns ####
# Function that removes constant columns of data-frame
rm_const <- function(data, cols = NULL, t_c = 0){
  ###-----------------------------------------------------------------------------
  ###Input
  #data...Data to check (data frame)
  #cols...Columns to check (String or integer vector)
  #.......Default: NULL -> Check all
  #t_c....Threshold for (almost) constant column (non-negative scalar)
  #.......Default: 0 -> Constant
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Columns of data that are not constant (String or integer vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Set cols if not given
  if(is.null(cols)){ cols <- colnames(data) }
  
  # Use only data that is needed
  data <- data[,cols]
  
  #### Remove columns ####
  # Number of samples to check with
  n_check <- min(10, nrow(data))
  
  # Check on sample which rows are candidates (-> computational more feasible)
  bool_res <- (apply(data[sample(1:nrow(data), n_check),], 2, sd) <= t_c)
  
  # Check if any of candidates is constant (Special case: Apply only on matrices)
  if(sum(bool_res) == 1){ bool_res[bool_res] <- (sd(data[,bool_res]) <= t_c) }
  else if(any(bool_res)){ bool_res[bool_res] <- (apply(data[,bool_res], 2, sd) <= t_c) }
  
  #### Output ####
  # Return columns that are not (almost) constant
  return(cols[!bool_res])
}

#### Remove and replace NA ####
prep_na <- function(df, var_omit = NULL, var_replace = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #df............Data to mutate (data frame)
  #var_omit......Remove NA of these variables (string or integer vector)
  #..............Default: NULL -> None of the variables
  #var_replace...Replace NA of these variables with specific mean (string or integer vector)
  #..............Default: NULL -> None of the variables
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Mutated data frame (data frame)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(lubridate)
  
  #### Mutate data ####
  # Remove NA
  if(!is.null(var_omit)){ df <- df[!apply(df[,var_omit], 1, function(x) any(is.na(x))),] }
  
  # Replace NA
  if(!is.null(var_replace)){
    # Get variables including NA
    na_vars <- var_replace[apply(df[,var_replace], 2, function(x) any(is.na(x)))]
    
    # For-Loop over variables
    for(temp_var in na_vars){
      # Rows with NA
      na_rows <- which(is.na(df[[temp_var]]))
      
      # Get months including NA
      na_months <- unique(month(df[na_rows, "init_tm"]))
      
      # Get locations including NA
      na_locs <- unique(df[na_rows, "location"])
       
      # Calculate mean matrix for each month and location
      mean_mtx <- sapply(na_months, function(x_m){
        temp_sub <- subset(df, month(init_tm) == x_m)
        sapply(na_locs, function(x_l){
          res <- mean(temp_sub[temp_sub[["location"]] == x_l, 
                               temp_var], na.rm = TRUE)
        
          # Special case: Only NA's for month and location
          if(!is.na(res)){ return(res) }
          else{ return( mean(subset(df, location == x_l)[[temp_var]], 
                             na.rm = TRUE) )}
          })
        })
      
      # Replace NA by corresponding mean
      df[na_rows, temp_var] <- sapply(na_rows, function(i) mean_mtx[which(na_locs == df[i, "location"]),
                                                                    which(na_months == month(df[i, "init_tm"]))] )
    }
  }
  
  #### Output ####
  # Return mutated data frame
  return(df)
}

#### Prepare data ####
prep_data <- function(df, loc_path = NULL, pred_ens = c("ens_mean", "ens_sd"), 
                      return_ens = TRUE, sort_ens = FALSE, n_ens = 20, 
                      t_sd = NULL, t_spread = NULL,
                      pred_loc = NULL, pred_tm = NULL, pred_add = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #df............Data to prepare (data frame)
  #loc_path......Path to location data (string)
  #..............Default: NULL -> Only needed for spatial predictor
  #pred_ens......Predictors of wind gust ensemble (string vector)
  #..............Default: Ensemble mean and standard deviation
  #return_ens....Query if ensemble should be returned (logical)
  #..............Default: TRUE -> Return
  #..............Comment: If ensemble members are in pred_ens, they are returned anyways (!!)
  #ens_sort......Query, if ensemble should be sorted in case a ensemble member is predictor (logical)
  #..............Default: FALSE -> Do not order
  #n_ens.........Ensemble size (integer)
  #..............Default: 20 -> COSMO-DE-EPS
  #t_sd/spread...(Lower) threshold for ens_sd/ens_spread (non-negative scalar)
  #..............Default: NULL -> Ignore
  #pred_loc/tm...Spatial/Temporal predictors (string vector)
  #..............Default: NULL -> None
  #pred_add......Predictors of other meteorological variables (string vector)
  #..............Default: NULL -> None
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Mutated data frame (data frame)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(lubridate)
  
  # Ensemble members
  ens_str <- paste0("ens_", 1:n_ens)
  
  # Return the following columns of df
  return_cols <- c("init_tm", "fc_step", "location", "obs",
                   pred_ens, pred_loc, pred_tm, pred_add)
  
  # If ensemble members are not predictors but should be returned
  if(return_ens & any(!is.element(ens_str, return_cols))){
    return_cols <- c(return_cols, ens_str) }
  
  # Read out location data
  if(any(is.element(c("lat", "lon", "altitude", "orog_DE", "orog_D2"), pred_loc))){
    load(paste0(loc_path, "loc_data.RData")) }
  
  #### Data preparation ####
  # Get sub-ensemble means if needed
  if(any(grepl("sens_mean_", pred_ens, fixed = TRUE))){ for(i_sens in 1:4){ 
    if(is.element(paste0("sens_mean_", i_sens), pred_ens)){
      df[[paste0("sens_mean_", i_sens)]] <- 
        rowMeans(df[,paste0("ens_", (1:5) + 5*(i_sens - 1))], na.rm = TRUE)
    }
  }}
  
  # Get ensemble spread if needed
  if(is.element("ens_spread", pred_ens)){
    df[["ens_spread"]] <- apply(df[,paste0("ens_", 1:n_ens)], 1, fn_spread) }
  
  # Get sub-ensemble spreads if needed
  if(any(grepl("sens_spread_", pred_ens, fixed = TRUE))){ for(i_sens in 1:4){ 
    if(is.element(paste0("sens_spread_", i_sens), pred_ens)){
      df[[paste0("sens_spread_", i_sens)]] <- 
        apply(df[,paste0("ens_", (1:5) + 5*(i_sens - 1))], 1, fn_spread)
    }
  }}
  
  # Omit and replace NA
  df <- prep_na(df = df,
                var_omit = c(pred_ens, "obs"),
                var_replace = pred_add)
  
  # Cut ensemble standard deviation (if given as predictor)
  if(!is.null(t_sd) & is.element("ens_sd", pred_ens)){
    df[["ens_sd"]] <- pmax(df[["ens_sd"]], t_sd) }
  
  # Cut ensemble spread (if given as predictor)
  if(!is.null(t_spread) & is.element("ens_spread", pred_ens)){
    df[["ens_spread"]] <- pmax(df[["ens_spread"]], t_spread) }
  
  # Cut sub-ensemble spread (if given as predictor)
  if(!is.null(t_spread) & any(grepl("sens_spread_", pred_ens, fixed = TRUE))){ 
    for(i_sens in 1:4){ if(is.element(paste0("sens_spread_", i_sens), pred_ens)){
      df[[paste0("sens_spread_", i_sens)]] <- pmax(df[[paste0("sens_spread_", i_sens)]], t_spread) 
  }}}

  # Sort ensemble (if given as predictor)
  if(sort_ens & any(is.element(pred_ens, ens_str))){ df[,ens_str] <- t(apply(df[,ens_str], 1, sort)) }
  
  #### Temporal predictors ####  
  # Cosine transformed month
  if(is.element("month", pred_tm)){ 
    df[["month"]] <- cos(2*pi*(month(df[["init_tm"]])-1)/12) }
  
  # Cosine transformed day of the year
  if(is.element("yday", pred_tm)){
    df[["yday"]] <- cos(2*pi*(yday(df[["init_tm"]])-1)/
                          yday(str2time(paste0(year(df[["init_tm"]]), "123100")))) }
  
  #### Spatial predictors ####  
  # Latitude
  if(is.element("lat", pred_loc)){ df[["lat"]] <- sapply(df[["location"]], function(x_loc){ 
    loc_data[["latitude"]][which(loc_data[["station_id"]] == x_loc)] }) }
  
  # Longitude
  if(is.element("lon", pred_loc)){ df[["lon"]] <- sapply(df[["location"]], function(x_loc){ 
    loc_data[["longitude"]][which(loc_data[["station_id"]] == x_loc)] }) }
  
  # Altitude
  if(is.element("altitude", pred_loc)){ df[["altitude"]] <- sapply(df[["location"]], function(x_loc){ 
    loc_data[["height"]][which(loc_data[["station_id"]] == x_loc)] }) }
  
  # Model surface height (at nearest grid point) of COSMO-DE
  if(is.element("orog_DE", pred_loc)){ df[["orog_DE"]] <- sapply(df[["location"]], function(x_loc){ 
    loc_data[["orog_DE"]][which(loc_data[["station_id"]] == x_loc)] }) }
  
  # Difference of station altitude and model surface height (at nearest grid point) of COSMO-DE
  if(is.element("orog_diff", pred_loc)){ df[["orog_diff"]] <- sapply(df[["location"]], function(x_loc){ 
    loc_data[["orog_DE"]][which(loc_data[["station_id"]] == x_loc)] - 
      loc_data[["height"]][which(loc_data[["station_id"]] == x_loc)] }) }
  
  # Calculate bias of location
  if(is.element("loc_bias", pred_loc)){
    # Vector of station IDs
    id_vec <- unique(df[["location"]])
    
    # Calculate mean bias of ensemble median at station
    # Consider NA because ens_str might not be in var_omit
    loc_bias <- sapply(id_vec, function(x_loc){
      df_loc <- na.omit(subset(df, location == x_loc)[,c("obs", paste0("ens_", 1:n_ens))])
      return(mean(df_loc[["obs"]] - apply(df_loc[,paste0("ens_", 1:n_ens)], 1, median), 
                  na.rm = TRUE))
    })
    
    # Assign bias of location in training and test data
    df[["loc_bias"]] <- sapply(df[["location"]], function(x_loc) loc_bias[which(id_vec == x_loc)] )
  }
  
  # Calculate coverage of location
  if(is.element("loc_cover", pred_loc)){
    # Calculate coverage at station
    # Consider NA because ens_str might not be in var_omit
    loc_cover <- sapply(id_vec, function(x_loc){
      df_loc <- na.omit(subset(df, location == x_loc)[,c("obs", paste0("ens_", 1:n_ens))])
      ens_range <- t(apply(df_loc[,paste0("ens_", 1:n_ens)], 1, range, na.rm = TRUE))
      return( mean((ens_range[,1] <= df_loc[["obs"]]) & (df_loc[["obs"]] <= ens_range[,2]),
                   na.rm = TRUE) )
    })
    
    # Assign bias of location in training and test data
    df[["loc_cover"]] <- sapply(df[["location"]], function(x_loc) loc_cover[which(id_vec == x_loc)] )
  }
  
  #### Output ####
  # Return data
  return(df[,unique(return_cols)])
}

#### Assign station-specific predictors to data ####
assign_spatial_preds <- function(preds, df_train, df_test){
  ###-----------------------------------------------------------------------------
  ###Input
  #preds......Station-specific predictors to be assigned ("loc_bias", "loc_cover")
  #df_train...Training data including preds (data frame)
  #df_test....Test data to include preds (data frame)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Mutated data frame (data frame)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Only certain predictors
  preds <- preds[is.element(preds, c("loc_bias", "loc_cover"))]
  
  # Stop if no predictors given
  if(length(preds) == 0){ return(df_test) }
  
  # Read out locations of given data
  df_locs <- unique(df_test[["location"]])
  
  # Functions for finding first location row
  f_help <- function(f, b) function(a) f(a, b)
  g_help <- function(x) { Position(f_help(`==`, x), df_train[["location"]]) }
  
  #### Assign predictors ####
  # For-Loop over predictors
  for(temp_pred in preds){
    # Get value of each locations in test set (via training set)
    temp <- sapply(df_locs, function(x) df_train[[temp_pred]][g_help(x)] )
    
    # Assign values to test set via locations
    df_test[[temp_pred]] <- sapply(df_test[["location"]], function(x){
      return( temp[which(df_locs == x)]) })
  }
  
  #### Output ####
  # Return data
  return(df_test)
}

