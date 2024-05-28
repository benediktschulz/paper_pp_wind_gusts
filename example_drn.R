## File with exemplary usage of DRN function

# Load package
library(lubridate)

# Path to station data
loc_path <- getwd()

#### Example with full data ####
### Initialization
# Load station data
load(file = paste0(loc_path, "/loc_data.RData"))

# Load training and test data
load(file = paste0(getwd(), "/ens_data_step18.RData"))

# Import postprocessing function
source(paste0(getwd(), "/pp_drn.R"))

# Wind gust ensemble predictors
pred_ens <- c("ens_mean", "ens_sd")

# Ensemble mean and standard deviation of additional variables as predictors
pred_add <- paste0(rep(x = add_vars, each = 2),
                   rep(x = c("_mean", "_sd"), times = length(add_vars)))

# Temporal predictors
pred_tm <- c("yday")

# Spatial predictors
pred_loc <- c("lat", "lon", "altitude", "orog_diff", "loc_cover", "loc_bias")

# Define predictor variables
pred_vars <- c(pred_ens, pred_add, pred_tm, pred_loc)

# Set validation set
i_valid <- which(year(df_train[["init_tm"]]) == 2015)


### Data preprocessing
# Prepare data (need to include ensemble mean for bias prediction)
df_train <- prep_data(df = df_train,
                      pred_ens = pred_ens,
                      pred_add = pred_add,
                      pred_tm = pred_tm,
                      pred_loc = pred_loc,
                      loc_path = loc_path)
df_test <- prep_data(df = df_test,
                     pred_ens = pred_ens,
                     pred_add = pred_add,
                     pred_tm = pred_tm,
                     pred_loc = pred_loc[!is.element(pred_loc, c("loc_bias", "loc_cover"))],
                     loc_path = loc_path)

# Assign spatial predictors to test set (bias and coverage are calculated on training set!)
df_test <- assign_spatial_preds(preds = pred_loc,
                                df_train = df_train,
                                df_test = df_test)

### Postprocessing application and evaluation
# Apply DRN function
pred <- drn_pp(train = df_train, 
               X = df_test, 
               i_valid = i_valid, 
               loc_id_vec = loc_data[["station_id"]],
               pred_vars = pred_vars,
               nn_ls = list(n_sim = 2,
                            nn_verbose = TRUE))

# Take a look at some forecasts
head(pred[["f"]])

# Take a look at evaluation measures of postprocessed forecasts
summary(pred[["scores_pp"]])

# PIT histogram
hist(x = pred[["scores_pp"]][["pit"]],
     freq = FALSE)
abline(h = 1,
       lty = 2)

# Calculate CRPSS w.r.t. ensemble forecasts (in %)
print(100*(1 - mean(pred[["scores_pp"]][["crps"]])/mean(pred[["scores_ens"]][["crps"]])))



#### Example with exemplary data from Github ####
# Load training and test data
load(file = paste0(getwd(), "/df_train.RData"))
load(file = paste0(getwd(), "/df_test.RData"))

# Set validation set
i_valid <- 800:984

# Import postprocessing function
source(paste0(getwd(), "/pp_drn.R"))

# Use variables besides observations and individual ensemble members as predictors
pred_vars <- names(df_train)[!is.element(names(df_train), c("obs", paste0("ens_", 1:20)))]

# Apply DRN function
pred <- drn_pp(train = df_train, 
               X = df_test, 
               i_valid = i_valid, 
               loc_id_vec = paste0(1:10),
               pred_vars = pred_vars,
               nn_ls = list(n_sim = 3,
                            nn_verbose = TRUE))

# Take a look at some forecasts
head(pred[["f"]])

# Take a look at evaluation measures of postprocessed forecasts
summary(pred[["scores_pp"]])

# PIT histogram
hist(x = pred[["scores_pp"]][["pit"]],
     freq = FALSE)
abline(h = 1,
       lty = 2)

# Calculate CRPSS w.r.t. ensemble forecasts (in %)
print(100*(1 - mean(pred[["scores_pp"]][["crps"]])/mean(pred[["scores_ens"]][["crps"]])))