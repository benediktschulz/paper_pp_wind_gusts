## Functions for data processing

#### Variable check ####
# Function to check variables
var_check <- function(y, check_type = "none", name_function = "var_check", name_y = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...............n variables (n vector)
  #check_type......Type of given variable ("probability", "binary", "none", "finite")
  #................Default: "none" -> check if NA or NaN
  #name_function...Name of function in which input is checked (string)
  #................Default: "var_check" -> function itself
  #name_y..........Name of variable (string)
  #................Default: NULL-> depending on check
  ###-----------------------------------------------------------------------------
  ###Output
  #...Warning if input accords to condition
  ###-----------------------------------------------------------------------------
  
  # Set name_y if not given
  if(is.null(name_y) & is.element(check_type, c("none", "finite"))){ name_y <- "variable" }
  else if(is.null(name_y) & (check_type == "probability")){ name_y <- "probability forecast" }
  else if(is.null(name_y) & (check_type == "binary")){ name_y <- "observation" }
  
  # Always check if NaN
  temp_bool <- is.nan(y)
  
  # Warning if nan is included
  if(any(temp_bool)){
    # Warning output
    warning(paste0("\n At least one ", name_y, " given to ", name_function, " is NaN!
                   \n Number and value: ", paste0(paste0(which(temp_bool), " (", y[temp_bool], ")"), collapse = ", ")))
  }
  
  # Always check if NA
  temp_bool <- is.na(y)
  
  # Warning if NA is included
  if(any(temp_bool)){
    # Warning output
    warning(paste0("\n At least one ", name_y, " given to ", name_function, " is NA!
                   \n Number and value: ", paste0(paste0(which(temp_bool), " (", y[temp_bool], ")"), collapse = ", ")))
  }
  
  # Check depending on condition
  if(check_type == "none"){ temp_bool <- FALSE }
  else if(check_type == "finite"){ temp_bool <- is.infinite(y) }
  else if(check_type == "probability"){ temp_bool <- ((y < 0) | (y > 1)) }
  else if(check_type == "binary"){ temp_bool <- !is.element(y, c(0, 1)) }
  
  # Warning output, if necessary
  if(any(temp_bool)){
    # Which condition is not satisfied?
    if(check_type == "probability"){ name_condition <- "in [0,1]" }
    else if(check_type == "binary"){ name_condition <- "binary" }
    else if(check_type == "finite"){ name_condition <- "finite" }
    
    # Warning output
    warning(paste0("\n At least one ", name_y, " given to ", name_function, " is not ", name_condition, "!
                   \n Number and value: ", paste0(paste0(which(temp_bool), " (", y[temp_bool], ")"), collapse = ", ")))
  }
}

#### Time and string conversions ####
# Function that converts string format of the data to time
str2time <- function(str_tm, tm_format = "%Y%m%d%H"){
  ###-----------------------------------------------------------------------------
  ###Input
  #str_tm......String representing the time (string)
  #tm_format...Date format of string (string)
  #............Default: "%Y%m%d%H" -> used in COSMO data set
  ###-----------------------------------------------------------------------------
  ###Output
  #tm...String transformed to time format/POSIXct (list; "time"/POSIXct type)
  ###-----------------------------------------------------------------------------
  
  # Convert time to string in given format
  return(strptime(str_tm, format = tm_format, tz = "UTC"))
}

# Function that converts time to the string format of the data
time2str <- function(tm, tm_format = "%Y%m%d%H"){
  ###-----------------------------------------------------------------------------
  ###Input
  #tm..........Time to be transformed (list; "time"/POSIXct type)
  #tm_format...Date format of string (string)
  #............Default: "%Y%m%d%H" -> used in COSMO data set
  ###-----------------------------------------------------------------------------
  ###Output
  #str...String representing the time to read out from data (string)
  ###-----------------------------------------------------------------------------
  
  # Convert time to string in given format
  return(format(tm, tm_format))
}

#### Update hyperparameters ####
update_hpar <- function(hpar_ls, in_ls){
  ###-----------------------------------------------------------------------------
  ###Input
  #hpar_ls...Default hyperparameter (list)
  #in_ls.....Selected hyperparameter given by user (list)
  ###-----------------------------------------------------------------------------
  ###Output
  #hpar_ls...All hyperparameters including users selection (list)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Names of hyperparameters
  hpar_names <- names(hpar_ls)
  
  # Names of hyperparameter to update
  in_names <- names(in_ls)
  
  #### Update ####
  # Loop over given names
  for(temp_hpar in in_names){
    # Update list if correct name is given
    if(is.element(temp_hpar, hpar_names)){ hpar_ls[[temp_hpar]] <- in_ls[[temp_hpar]] }
    else{ print(paste0("Wrong hyperparameter given: ", temp_hpar))}
  }
  
  #### Output ####
  # Return list
  return(hpar_ls)
}

#### Ensemble spread ####
# Calculate spread of ensemble (formula from Jordan et al. (2019))
fn_spread <- function(x){
  ###-----------------------------------------------------------------------------
  ###Input
  #x...Ensemble forecast (n_ens vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Ensemble spread
  ###-----------------------------------------------------------------------------
  
  #### Spread calculation ####
  # Get ensemble size
  n_ens <- length(x)
  
  # Matrix for double sum
  A <- matrix(data = x,
              ncol = n_ens,
              nrow = n_ens)
  
  # Calculate return spread
  return( sum(abs(A - t(A)), na.rm = TRUE) / (n_ens^2) )
}

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
