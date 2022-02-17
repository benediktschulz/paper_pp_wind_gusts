# Machine learning methods for postprocessing ensemble forecasts of wind gusts: A systematic comparison

This repository provides R-code accompanying the paper

> Schulz, B. and Lerch, S. (2022). 
> Machine learning methods for postprocessing ensemble forecasts of wind gusts: A systematic comparison. 
> Monthly Weather Review, 150, 235-257, https://doi.org/10.1175/MWR-D-21-0150.1 (preprint version available at https://arxiv.org/abs/2106.09512).

In particular, code for the implementation and evaluation of the proposed postprocessing methods is available. In addition, two scripts exemplify the usage of the postprocessing functions.


## Data

The data was supplied by the German weather service (Deutscher Wetterdienst, DWD) and is not publicly available.

### Ensemble Forecasts

**COSMO-DE-EPS**

- Variables: Wind gust and various additional predictors
- Time period: 2010-12-08 - 2016-12-31
- Forecast initialization time: 00 UTC
- Members: 20
- Forecast lead times: 00-21 h
- Area: Germany
- Resolution: 2.8 km horizontal resolution

More information on the COSMO model can be found here: http://www.cosmo-model.org/.

### Observations

**DWD SYNOP stations**

- Forecasts: Taken from closest grid points
- Number of stations: 175
- Attributes: longitude, latitude, altitude, altitude of closest COSMO grid point


### Exemplary data set

We supply an additional training (`df_train.RData`) and test set (`df_test.RData`) that is derived from the data used in the paper together with the following comments:

- The forecasts in the training set are from the period of 2010-2015, those in the test set from 2016. The initialization times are not supplied.
- The forecasts have a common lead time, which is not supplied.
- Ten different stations are used and labeled by 1,...,10. This does not conform with the actual station IDs, which are not supplied.
- We are using only a small subset of the available predictor variables.
- We added random noise to all numeric variables besides the transformed day of the year.


## Post-processing

All models are estimated based on the period 2011-2015 and evaluated on 2016.

### Basic methods

The basic models make only use of the ensemble forecasts of wind gusts.

- Local Ensemble Model Output Statistics (EMOS) with CRPS estimation
- Local Member-by-Member (MBM) with CRPS estimation
- Local Isotonic Distributional Regression (IDR)

### Incorporating additional informtion

The machine learning-based methodsmake incorporation of additional features feasible.

- Local EMOS with gradient boosting (EMOS-GB) estimated via MLE
- Local Quantile Regression Forests (QRF)

### Locally adaptive networks

Based on neural networks and station embedding we built a locally adaptive joint model for all stations. The three variants are based on the same architecture and differ only in estimation and forecast type. All network models are estimated using stochastic gradient descent.

- Distributional Regression Network (DRN) estimated via CRPS and aggregated via parameter averaging
- Bernstein Quantile Network (BQN) estimated via Quantile Loss and aggregated via quantile averaging (Vincentization)
- Histogram Estimation Network (HEN) estimated via MLE and aggregated via quantile averaging (Vincentization)


## Code

Due to the fact that the data cannot be distributed publicly, we supply only the code used for postprocessing and evaluation. The scripts for data preprocessing and application of the provided functions are not supplied. However, two exemplary scripts for the usage of the postprocessing functions are given.
Each of the R-files includes functions that can be applied to data of the desired format.

| File | Description |
| ---- | ----------- | 
| `example_drn` | Example for the use of DRN. |
| `example_emos` | Example for the use of EMOS. |
| `fn_data` | Data processing. |
| `fn_eval` | Evaluation. |
| `pp_bqn` | Postprocessing via BQN. |
| `pp_drn` | Postprocessing via DRN. |
| `pp_emos` | Postprocessing via EMOS. |
| `pp_emos_bst` | Postprocessing via EMOS-GB. |
| `pp_hen` | Postprocessing via HEN. |
| `pp_idr` | Postprocessing via IDR. |
| `pp_mbm` | Postprocessing via MBM. |
| `pp_qrf` | Postprocessing via QRF. |


### Data structure used

The functions are based on the structure of the COSMO data, which is given by R-dataframes (generated from netCDF-files in a preprocessing step). We hope that the following comments will help the reader in understanding the code. 

The following table describes variable names that are referred to in the code (and the example data):

| Variable | Description |
| ---- | ----------- | 
| `obs` | Observed speed of a wind gust in m/s. |
| `location` | Station ID. |
| `ens_mean` | Mean of the wind gust ensemble in m/s. |
| `ens_sd` | Standard deviation of the wind gust ensemble in m/s. |
| `ens_i` | i-th Member of the wind gust ensemble in m/s, i = 1, ..., 20. |
| `sens_mean_i` | Mean of the i-th sub-ensemble of the wind gust ensemble in m/s, i = 1, ..., 4. The four sub-ensembles are given by the members 1-5 (i = 1), 6-10 (i = 2), 11-15 (i = 3) and 16-20 (i = 4). |
| `sens_spread_i` | Spread of the i-th sub-ensemble of the wind gust ensemble in m/s, i = 1, ..., 4. The four sub-ensembles are given by the members 1-5 (i = 1), 6-10 (i = 2), 11-15 (i = 3) and 16-20 (i = 4). |
