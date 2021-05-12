## File containing functions for evaluation of ensemble and postprocessed forecasts

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_basic.R"))

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

#### Coverage ####
# Calculate coverage of a prediction interval of a probabilistic forecast
fn_cover <- function(x, alpha = NULL, n_ens = 20){
  ###-----------------------------------------------------------------------------
  ###Input
  #x.......PIT values / ranks (n vector)
  #alpha...Significance level (probability)
  #........Default: NULL -> Nominal n_ens-member coverage ((n_ens - 1)/(n_ens + 1))
  #n_ens...Size of ensemble (integer)
  #........Default: 20 -> COSMO-DE-EPS
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Coverage in percentage
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Nominal n_ens coverage
  if(is.null(alpha)){ alpha <- 2/(n_ens + 1) }
  
  #### Coverage calculation ####
  # PIT or rank?
  if(min(x) < 1){ res <- mean((alpha/2 <= x) & (x <= (1 - alpha/2))) }
  else{ res <- mean(is.element(x, 2:n_ens)) }
  
  # Output as percentage
  return(100*res)
}

#### Brier score ####
# Brier score for given distribution or ensemble
brier_score <- function(f, y, t = 0, distr = "ens", t_distr = 0){
  ###-----------------------------------------------------------------------------
  ###Input
  #f...........distr == "tlogis": Matrix with location and scale of forecast distribution (n x n_par matrix)
  #............distr == "ens": Ensemble forecasts (n x n_ens matrix)
  #............distr == "p" (or elsewise): Probability forecasts of exceeding t (n vector)
  #y...........Observations (n vector)
  #t...........Brier Score Threshold (non-negative scalar)
  #............Default: 0
  #distr.......Forecast type (specific distribution or ensemble) ("ens", "tlogis", "p")
  #............Default: "ens" -> Ensemble
  #t_distr.....Threshold for censored or truncated distribution (Scalar)
  #............Default: 0
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Brier Scores of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  ## Brier Score w.r.t. event of exceeding threshold t
  ## p_t = 1 - F(t). BS_t (F, y) = (p_t - 1(y > t))^2 = (F(t) - 1(y <= t))^2
  
  # Calculate F(t) depending on distribution
  if(distr == "ens"){ 
    if(is.vector(f)){ f <- mean(f <= t) }
    else{ f <- rowMeans(f <= t) }
  }
  # Truncated logistic
  else if(distr == "tlogis"){ f <- (t > t_distr)*crch::ptlogis(q = t, 
                                                               location = f[,1], 
                                                               scale = f[,2],
                                                               left = t_distr) }
  
  # Calculate Brier Score
  res <- (f - (y <= t))^2
  
  # Return score
  return(res)
}

#### BQN: Bernstein Quantile function ####
# Function that calculates quantiles for given coefficients
bern_quants <- function(alpha, q_levels){
  ###-----------------------------------------------------------------------------
  ###Input
  #alpha......Coefficients of Bernstein Basis (n x (p_degree + 1) matrix)
  #q_levels...Quantile levels (n_q vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Quantile forecasts for given coefficients (n x n_q matrix)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Get degree of polynomials from coefficients
  if(is.vector(alpha)){ p_degree <- length(alpha) - 1 }
  else{ p_degree <- ncol(alpha) - 1 }
  
  #### Calculation ####
  # Calculate quantiles (sum of coefficients times basis polynomials)
  res <- alpha %*% t(sapply(0:p_degree, dbinom, size = p_degree, prob = q_levels))
  
  # Return quantiles
  return(res)
}

#### HEN: Histogram distribution functions ####
# CDF of a histogram distribution
cdf_hd <- function(y, probs, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...........Observations (n vector)
  #probs.......Probabilities of corresponding bins (n_bins vector)
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...CDF evaluated at y (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Left/right edges of bins
  a <- bin_edges[-length(bin_edges)]
  b <- bin_edges[-1]
  
  #### Calculation ####
  # Calculate bin of observations
  k <- sapply(y, function(x) sum(x >= a) )
  
  # Calculate cumulative sums of bin probabilities
  # Entry k corresponds to sum of all k-1 bins below
  p_cum <- c(0, cumsum(probs)[-length(probs)])
  
  # Calculate CDF value via formula dependent on bin
  res <- p_cum[k] + probs[k]*(y - a[k])/(b[k] - a[k])
  
  # Output
  return(res)
}

# PDF of a histogram distribution
pdf_hd <- function(y, probs, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...........Observations (n vector)
  #probs.......Probabilities of corresponding bins (n_bins vector)
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...PDF evaluated at y (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Left/right edges of bins
  a <- bin_edges[-length(bin_edges)]
  b <- bin_edges[-1]
  
  #### Calculation ####
  # Transform probabilities
  p_tilde <- probs/(b - a)
  
  # Get transformed probability of corresponding bin
  res <- p_tilde[sapply(y, function(x) sum(x >= a) )]
  
  # Output
  return(res)
}

# Quantile function of a histogram distribution
quant_hd <- function(tau, probs, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #tau.........Quantile levels (n vector)
  #probs.......Probabilities of corresponding bins (n_bins vector)
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Quantiles at levels tau (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Left/right edges of bins
  a <- bin_edges[-length(bin_edges)]
  b <- bin_edges[-1]
  
  #### Calculation ####
  # Calculate cumulative sums of bin probabilities
  # Entry k corresponds to sum of all k-1 bins below
  p_cum <- c(0, cumsum(probs)[-length(probs)])
  
  # Calculate bin of quantiles (works for tau = 1)
  k <- sapply(tau, function(t) max(which(t >= p_cum)))
  
  # Calculate quantiles via formula dependent on bin
  res <- a[k] + (tau - p_cum[k])*(b[k] - a[k])/probs[k]
  
  # Output
  return(pmin(max(b), res))
}

#### HEN: Evaluation measures for a histogram distribution ####
# CRPS of a histogram distribution
crps_hd <- function(y, f, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...........Observations (n vector)
  #f...........Probabilities of corresponding bins (n x n_bins matrix)
  #............or: n list of (different) n_bins vectors
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  #............or: n x (n_bins + 1) matrix
  #............or: n list of (different) (n_bins + 1) vectors
  ###-----------------------------------------------------------------------------
  ###Output
  #res...CRPS of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # CRPS for a single forecast-observation pair
  fn_single <- function(obs, probs, bin_edges0){
    ## obs.........Wind gust observation (scalar)
    ## probs.......Bin probabilities (n_bins vector)
    ## bin_edges...Boundaries of bins ((n_bins + 1) vector)
    
    # Lower and upper edges of bins
    a <- bin_edges0[-length(bin_edges0)]
    b <- bin_edges0[-1]
    
    # Realizing values for individual uniform distributions
    z <- ((a <= obs) & (obs <= b))*obs + (obs < a)*a + (b < obs)*b
    if(obs < a[1]){ z[1] <- obs }
    if(obs > b[length(b)]){ z[length(b)] <- obs }
    
    # Lower and upper masses of individual uniform distributions
    L <- cumsum(c(0, probs[1:(length(probs)-1)]))
    U <- 1 - cumsum(probs)
    
    # Transform from bin to unit interval ([a_i, b_i] -> [0, 1])
    w <- (z - a)/(b - a)
    
    # Sum of standard uniform with lower mass L and upper mass U
    out <- sum((b - a)*( abs(w - punif(w)) 
                         + (1 - L - U)*punif(w)^2 
                         - punif(w)*(1 - 2*L) 
                         + ((1 - L - U)^2)/3 
                         + (1- L)*U ))
    
    # Output
    return(out)
  }
  
  # Function for apply (identical bin_edges?)
  if(is.list(bin_edges)){
    fn_apply <- function(i){ fn_single(obs = y[i],
                                       probs = f[[i]],
                                       bin_edges0 = bin_edges[[i]])} }
  else if(is.vector(bin_edges)){
    fn_apply <- function(i){ fn_single(obs = y[i],
                                       probs = f[i,],
                                       bin_edges0 = bin_edges)} }
  else if(is.matrix(bin_edges)){
    fn_apply <- function(i){ fn_single(obs = y[i],
                                       probs = f[i,],
                                       bin_edges0 = bin_edges[i,])} }
  
  #### Calculation ####
  # Apply function on all values
  res <- sapply(1:length(y), fn_apply)
  
  # Return
  return(res)
}

# Log-Score of a histogram distribution
logs_hd <- function(y, f, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...........Observations (n vector)
  #f...........Probabilities of corresponding bins (n x n_bins matrix)
  #............or: n list of (different) n_bins vectors
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  #............or: n x (n_bins + 1) matrix
  #............or: n list of (different) (n_bins + 1) vectors
  ###-----------------------------------------------------------------------------
  ###Output
  #res...PIT value of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Function for apply (identical bin_edges?)
  if(is.list(bin_edges)){
    fn_apply <- function(i){ -log(pdf_hd(y = y[i],
                                         probs = f[[i]],
                                         bin_edges = bin_edges[[i]]))} }
  else if(is.vector(bin_edges)){
    fn_apply <- function(i){ -log(pdf_hd(y = y[i],
                                         probs = f[i,],
                                         bin_edges = bin_edges))} }
  else if(is.matrix(bin_edges)){
    fn_apply <- function(i){ -log(pdf_hd(y = y[i],
                                         probs = f[i,],
                                         bin_edges = bin_edges[i,]))} }
  
  #### Calculation ####
  # Calculate logarithm of corresponding probability
  res <- sapply(1:length(y), fn_apply)
  
  # Output
  return(res)
}

# PIT of a histogram distribution
pit_hd <- function(y, f, bin_edges){
  ###-----------------------------------------------------------------------------
  ###Input
  #y...........Observations (n vector)
  #f...........Probabilities of corresponding bins (n x n_bins matrix)
  #............or: n list of (different) n_bins vectors
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  #............or: n x (n_bins + 1) matrix
  #............or: n list of (different) (n_bins + 1) vectors
  ###-----------------------------------------------------------------------------
  ###Output
  #res...PIT value of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Function for apply (identical bin_edges?)
  if(is.list(bin_edges)){
    fn_apply <- function(i){ cdf_hd(y = y[i],
                                    probs = f[[i]],
                                    bin_edges = bin_edges[[i]])} }
  else if(is.vector(bin_edges)){
    fn_apply <- function(i){ cdf_hd(y = y[i],
                                    probs = f[i,],
                                    bin_edges = bin_edges) } }
  else if(is.matrix(bin_edges)){
    fn_apply <- function(i){ cdf_hd(y = y[i],
                                    probs = f[i,],
                                    bin_edges = bin_edges[i,])} }
  
  #### Calculation ####
  # Plug in CDF
  res <- sapply(1:length(y), fn_apply)

  # Output (May be > 1 due to numerical reasons)
  return(pmin(1, res))
}

# Prediction interval length of a histogram distribution
lgt_hd <- function(f, bin_edges, n_ens = 20){
  ###-----------------------------------------------------------------------------
  ###Input
  #f...........Probabilities of corresponding bins (n x n_bins matrix)
  #............or: n list of (different) n_bins vectors
  #bin_edges...Boundaries of bins ((n_bins + 1) vector)
  #............or: n x (n_bins + 1) matrix
  #............or: n list of (different) (n_bins + 1) vectors
  #n_ens.......Ensemble size for significance level (integer)
  #............Default: 20 -> COSMO-DE-EPS
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Prediction interval length of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Half of significance level of prediction interval
  alpha <- 1/(n_ens + 1)
  
  #### Calculation ####
  # Calculate via quantile function(identical bin_edges?)
  if(is.list(bin_edges)){
    res <- sapply(1:length(f), function(i) diff(quant_hd(tau = c(alpha, 1 - alpha),
                                                         probs = f[[i]],
                                                         bin_edges = bin_edges[[i]]))) }
  else if(is.vector(bin_edges)){
    res <- apply(f, 1, function(probs) diff(quant_hd(tau = c(alpha, 1 - alpha),
                                                     probs = probs,
                                                     bin_edges = bin_edges))) }
  else if(is.matrix(bin_edges)){
    res <- sapply(1:nrow(f), function(i) diff(quant_hd(tau = c(alpha, 1 - alpha),
                                                       probs = f[i,],
                                                       bin_edges = bin_edges[i,]))) }
  # Output
  return(res)
}

# Function for prediction based on the bin probabilities #
fn_scores_hd <- function(f, y, bin_edges, n_ens = 20, skip_evals = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #f............Bin probabilities (n x n_bins matrix)
  #.............or: n list of (different) n_bins vectors
  #y............Observations (n vector)
  #bin_edges....Boundaries of bins ((n_bins + 1) vector)
  #.............or: n x (n_bins + 1) matrix
  #.............or: n list of (different) (n_bins + 1) vectors
  #n_ens........Ensemble size (integer)
  #.............Used for confidence level of prediction intervals
  #skip_evals...Skip the following evaluation measures (string vector)
  #.............Default: NULL -> Calculate all
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......scores_pp...Data frames containing (n x 4 data frame):
  #.........pit.........PIT values of forecasts (n vector)
  #.........crps........CRPS of forecasts (n vector)
  #.........logs........Log-Score of forecasts (n vector)
  #.........lgt.........Length of prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Number of predictions
  n <- length(y)
  
  #### Prediction and score calculation ####
  # Make data frame
  scores_pp <- data.frame(pit = numeric(length = n),
                          crps = numeric(length = n),
                          logs = numeric(length = n),
                          lgt = numeric(length = n),
                          e_me = numeric(length = n),
                          e_md = numeric(length = n))
  
  # Calculate PIT values
  if(is.element("pit", colnames(scores_pp))){
    scores_pp[["pit"]] <- pit_hd(y = y,
                               f = f,
                               bin_edges = bin_edges) }
  
  # Calculate CRPS
  if(is.element("crps", colnames(scores_pp))){
    scores_pp[["crps"]] <- crps_hd(y = y,
                                   f = f,
                                   bin_edges = bin_edges) }
  
  # Calculate Log-Score
  if(is.element("logs", colnames(scores_pp))){
    scores_pp[["logs"]] <- logs_hd(y = y,
                                   f = f,
                                   bin_edges = bin_edges) }
  
  # Calculate length of ~(n_ens-1)/(n_ens+1) % prediction interval
  if(is.element("lgt", colnames(scores_pp))){
    scores_pp[["lgt"]] <- lgt_hd(f = f,
                                 bin_edges = bin_edges,
                                 n_ens = n_ens) }
  
  # Calculate bias of median forecast (identical bins?)
  if(is.element("e_md", colnames(scores_pp))){
    if(is.list(bin_edges)){
      scores_pp[["e_md"]] <- sapply(1:length(y), function(i){
        quant_hd(tau = 0.5,
                 probs = f[[i]],
                 bin_edges = bin_edges[[i]]) }) - y }
    else if(is.vector(bin_edges)){
      scores_pp[["e_md"]] <- apply(f, 1, function(probs){
        quant_hd(tau = 0.5,
                 probs = probs,
                 bin_edges = bin_edges) }) - y }
    else if(is.matrix(bin_edges)){
      scores_pp[["e_md"]] <- sapply(1:nrow(f), function(i){
        quant_hd(tau = 0.5,
                 probs = f[i,],
                 bin_edges = bin_edges[i,]) }) - y }
  }
  
  # Calculate bias of mean forecast (t(t(..)..) is column-wise multiplication) (identical bins?)
  if(is.element("e_me", colnames(scores_pp))){
    if(is.list(bin_edges)){
      scores_pp[["e_me"]] <- ( sapply(1:length(y), function(i){
        sum( f[[i]]*((bin_edges[[i]])[-length(bin_edges[[i]])] + (bin_edges[[i]])[-1])/2 )
      }) - y ) }
    else if(is.vector(bin_edges)){
      scores_pp[["e_me"]] <- ( rowSums(t(t(f) * (bin_edges[-length(bin_edges)] + 
                                                   bin_edges[-1])/2)) - y
      ) }
    else if(is.matrix(bin_edges)){
      scores_pp[["e_me"]] <- ( sapply(1:nrow(f), function(i){
        sum( f[i,]*(bin_edges[i, -length(bin_edges)] + bin_edges[i, -1])/2 )
      }) - y ) }
  }
  
  #### Output ####
  # Skip evaluation measures
  scores_pp <- as.data.frame(scores_pp[,!is.element(colnames(scores_pp), skip_evals), drop = FALSE])
  
  # Return
  return(scores_pp)
}

#### ENS: Evaluation of ensemble ####
# Function to calculate evaluation measures of scores
fn_scores_ens <- function(ens, y, skip_evals = NULL, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #ens..........Ensemble data for prediction (n x n_ens matrix)
  #y............Observations for prediction (n vector)
  #skip_evals...Skip the following evaluation measures (string vector)
  #.............Default: NULL -> Calculate all
  #scores_ens...Should scores of ensemble forecasts, interval lengths, ranks be calculated? (logical)
  #.............Default: TRUE -> Calculate
  ###-----------------------------------------------------------------------------
  ###Output
  #...scores_ens...Data frames containing (n x 4 data frame):
  #......rank......Ranks of observations in ensemble forecasts (n vector)
  #......crps......CRPS of ensemble forecasts (n vector)
  #......logs......Log-Score of ensemble forecasts (n vector)
  #......lgt.......Ensemble range (n vector)
  #......e_md......Bias of median forecast (n vector)
  #......e_me......Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Calculate only if scores_ens is TRUE
  if(!scores_ens){ return(FALSE) }
  
  # Size of COSMO-DE-EPS
  n_cosmo <- 20
  
  # Check if vector is given
  if(is.vector(ens)){ ens <- matrix(data = ens,
                                    nrow = 1) }
  
  # Get number of ensembles
  n <- nrow(ens)
  
  # Get ensemble size
  n_ens <- ncol(ens)
  
  # Make data frame
  scores_ens <- data.frame(rank = numeric(length = n),
                           crps = numeric(length = n),
                           logs = numeric(length = n),
                           lgt = numeric(length = n),
                           e_me = numeric(length = n),
                           e_md = numeric(length = n))
  
  #### Calculation ####
  # Calculate observation ranks
  if(is.element("rank", colnames(scores_ens))){
    scores_ens[["rank"]] <- apply(cbind(y, ens), 1, function(x){ rank(x, ties = "random")[1] }) }
  
  # Calculate CRPS of raw ensemble
  if(is.element("crps", colnames(scores_ens))){
    scores_ens[["crps"]] <- crps_sample(y = y, 
                                        dat = ens) }
  
  # Calculate Log-Score of raw ensemble
  if(is.element("logs", colnames(scores_ens))){
    scores_ens[["logs"]] <- logs_sample(y = y, 
                                        dat = ens) }
  
  # Calculate ~(20 - 1)/(20 + 1)% prediction interval (corresponds to COSMO ensemble range)
  if(is.element("lgt", colnames(scores_ens))){
    # COSMO-ensemble size
    if(n_ens == n_cosmo){ scores_ens[["lgt"]] <- apply(t(apply(ens, 1, range)), 1, diff) }
    # Corresponding quantiles (1/21 and 20/21) are included
    else if(((n_ens + 1) %% (n_cosmo + 1)) == 0){ 
      # Indices of corresponding quantiles
      i_lgt <- (n_ens + 1)/(n_cosmo + 1)*c(1, n_cosmo)
      
      # Get quantiles
      q_lgt <- t(apply(ens, 1, sort))[,i_lgt]
      
      # Transform if vector
      if(n == 1){ q_lgt <- matrix(data =  q_lgt,
                                  nrow = 1) }
      
      # Calculate corresponding range
      scores_ens[["lgt"]] <- apply(t(apply(q_lgt, 1, range)), 1, diff) 
    }
    # Quantiles are not included: Calculate corresponding via quantile function
    else{ 
      scores_ens[["lgt"]] <- apply(ens, 1, function(x) 
        diff(quantile(x = x,
                      probs = c(1, 20)/(n_cosmo + 1))) ) 
    }
  }
  
  # Calculate bias of median forecast
  if(is.element("e_md", colnames(scores_ens))){
    scores_ens[["e_md"]] <- apply(ens, 1, median) - y }
  
  # Calculate bias of mean forecast
  if(is.element("e_me", colnames(scores_ens))){
    scores_ens[["e_me"]] <- rowMeans(ens) - y }
  
  #### Output ####
  # Skip evaluation measures
  scores_ens <- as.data.frame(scores_ens[,!is.element(colnames(scores_ens), skip_evals), drop = FALSE])
  
  # Return output
  return(scores_ens)
}

#### TLOGIS: Evaluation of distributional forecasts ####
# Function for prediction based on the distributional parameters #
fn_scores_distr <- function(f, y, n_ens = 20, skip_evals = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #f............Parameters of forecast distribution (n x n_par matrix)
  #y............Observations (n vector)
  #n_ens........Ensemble size (integer)
  #.............Used for confidence level of prediction intervals
  #.............Default: 20 -> COSMO-DE-EPS
  #skip_evals...Skip the following evaluation measures (string vector)
  #.............Default: NULL -> Calculate all
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......scores_pp...Data frames containing (n x 4 data frame):
  #.........pit.........PIT values of distributional forecasts (n vector)
  #.........crps........CRPS of forecasts (n vector)
  #.........logs........Log-Score of forecasts (n vector)
  #.........lgt.........Length of prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Input check
  if(any(f[,2] < 0)){ print("Non-positive scale forecast!") }
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(f)
  
  # Make data frame
  scores_pp <- data.frame(pit = numeric(length = n),
                          crps = numeric(length = n),
                          logs = numeric(length = n),
                          lgt = numeric(length = n),
                          e_me = numeric(length = n),
                          e_md = numeric(length = n))
  
  #### Prediction and score calculation ####
  # Calculate PIT values
  if(is.element("pit", colnames(scores_pp))){
    scores_pp[["pit"]] <- crch::ptlogis(q = y, 
                                        location = f[,1], 
                                        scale = f[,2],
                                        left = 0) }
  
  # Calculate CRPS of forecasts
  if(is.element("crps", colnames(scores_pp))){
    scores_pp[["crps"]] <- crps_tlogis(y = y, 
                                       location = f[,1], 
                                       scale = f[,2],
                                       lower = 0) }
  
  # Calculate Log-Score of forecasts
  if(is.element("logs", colnames(scores_pp))){
    scores_pp[["logs"]] <- logs_tlogis(y = y, 
                                       location = f[,1], 
                                       scale = f[,2],
                                       lower = 0) }
  
  # Calculate length of ~(n_ens-1)/(n_ens+1) % prediction interval
  if(is.element("lgt", colnames(scores_pp))){
    scores_pp[["lgt"]] <- crch::qtlogis(p = n_ens/(n_ens + 1), 
                                        location = f[,1], 
                                        scale = f[,2],
                                        left = 0) - crch::qtlogis(p = 1/(n_ens + 1), 
                                                                  location = f[,1], 
                                                                  scale = f[,2],
                                                                  left = 0) }
  
  # Calculate bias of median forecast
  if(is.element("e_md", colnames(scores_pp))){
    scores_pp[["e_md"]] <- crch::qtlogis(p = 0.5, 
                                         location = f[,1], 
                                         scale = f[,2],
                                         left = 0) - y }
  
  # Calculate bias of mean forecast
  if(is.element("e_me", colnames(scores_pp))){
    scores_pp[["e_me"]] <- (f[,1] - f[,2]*log(1 - plogis(- f[,1]/f[,2])))/(1 - plogis(- f[,1]/f[,2])) - y }

  #### Output ####
  # Skip evaluation measures
  scores_pp <- as.data.frame(scores_pp[,!is.element(colnames(scores_pp), skip_evals), drop = FALSE])
  
  # Return
  return(scores_pp)
}