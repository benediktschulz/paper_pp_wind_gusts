## File with exemplary usage of DRN function

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