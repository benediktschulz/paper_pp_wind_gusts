## File with exemplary usage of DRN function

# Load training and test data
load(file = paste0(getwd(), "/df_train.RData"))
load(file = paste0(getwd(), "/df_test.RData"))

# Consider only one location for local EMOS?
df_train <- subset(df_train, location == "2")
df_test <- subset(df_test, location == "2")

# Import postprocessing function
source(paste0(getwd(), "/pp_emos.R"))

# Estimate EMOS parameters
est <- emos_est(train = df_train)

# Take a look at estimated parameters
print(c("a" = est$par[1], 
        "b" = exp(est$par[2]),
        "c" = est$par[3],
        "d" = est$par[4]))

# Predict via EMOS
pred <- emos_pred(X = df_test, 
                  par_emos = est$par)

# Take a look at evaluation measures of postprocessed forecasts
summary(pred[["scores_pp"]])

# Calculate CRPSS
print(100*(1 - mean(pred[["scores_pp"]][["crps"]])/mean(pred[["scores_ens"]][["crps"]])))