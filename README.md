# Statistical Postprocessing of ensemble forecasts of wind gusts

This repository provides R code accompanying the paper

> Schulz, B. and Lerch, S. (2021). 
> Name of the paper. 
> available at -

In particular, code for data processing as well as implementation and evaluation of the proposed postprocessing methods is available.

## Data

The data was supplied by the DWD and is not publicly available.

### Forecasts

**COSMO-DE-EPS**

Link to COSMO

- Variables: Wind gust and various additional predictors
- Time period: 2010-??-?? -- 2016-12-31
- Forecast initialization time: 00UTC
- Members: 20
- Forecast lead time: 00--21h (valid at 00--21UTC)
- Area: Germany
- Resolution: 2.8 kilometres

### Observations

**DWD SYNOP stations**

- Number of stations: 175
- Attributes: longitude, latitude, altitude, altitude of closest COSMO grid point


## Post-processing

All models are estimated based on the period 2011--2015 and evaluated on 2016.

### Basic methods

The basic models make only use of the ensemble forecasts of wind gusts

- Local Ensemble Model Output Statistics (EMOS) with CRPS estimation
- Local Member-by-Member (MBM) with CRPS estimation
- Local Isotonic Distributional Regression (IDR)

### Machine learning methods

The machine learning-based methods incorporate additional features

- Local EMOS with gradient boosting (EMOS-GB) estimated via MLE
- Local Quantile Regression Forests (QRF)
    
All network models are estimated by minimizing the CRPS using stochastic gradient descent.

### Locally adaptive networks

The locally adaptive networks are neural network-based postprocessing methods using station embedding to build one model for all stations. The three variants are based on the same architecture and differ only in estimation and forecast type.

- Distributional Regression Network (DRN) estimated via CRPS and aggregated via parameter averaging
- Bernstein Quantile Network (BQN) estimated via Quantile Loss and aggregated via quantile averaging (Vincentization)
- Histogram Estimation Network (HEN) estimated via MLE and aggregated via quantile averaging (Vincentization)
    
All network models are estimated using stochastic gradient descent.

## Code

Each R-file includes functions that can be applied to data of the desired format.

| File | Description |
| ---- | ----------- | 
| `fn_basic` | Basic function. |
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

## Comments on the data structure

The functions are based on the structure of the COSMO data, which was preprocessed to R-dataframes. We hope that the following comments will help the reader in understanding the code. 

The following table lists the variables names that we use in our data:

| Variable | Description |
| ---- | ----------- | 
| `obs` | Observed speed of a wind gust. |
| `location` | Station ID. |
| `ens_mean` | Mean of the wind gust ensemble. |
| `ens_sd` | Standard deviation of the wind gust ensemble. |
| `ens_i` | Member i of the wind gust ensemble, i = 1, ..., 20. |
| `sens_mean_i` | Mean of the i-th sub-ensemble (members 5*(i-1) + 1 to 5*i) of the wind gust ensemble, i = 1, ..., 4. |
| `sens_spread_i` | Spread of the i-th sub-ensemble (members 5*(i-1) + 1 to 5*i) of the wind gust ensemble, i = 1, ..., 4. |

The data processing further relies on the dataframe `loc_data` that includes the following variables:

| Variable | Description |
| ---- | ----------- | 
| `name` | Name of the station. |
| `station_id` | Station ID. |
| `latitude` | Latitude of the station in degree. |
| `longitude` | Longitude of the station in degree. |
| `height` | Altitude of the station in meters. |
| `orog_DE` | Model height of the closest grid point in meters. |
