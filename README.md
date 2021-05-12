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

The implementations are available in the directory: ?.

### Machine learning methods

The machine learning-based methods incorporate additional features

- Local EMOS with gradient boosting (EMOS-GB) estimated via MLE
- Local Quantile Regression Forests (QRF)
    
All network models are estimated by minimizing the CRPS using stochastic gradient descent. Implementations are available in the directory: ?.

### Locally adaptive networks

The locally adaptive networks are neural network-based postprocessing methods using station embedding to build one model for all stations. The three variants are based on the same architecture and differ only in estimation and forecast type.

- Distributional Regression Network (DRN) estimated via CRPS
- Bernstein Quantile Network (BQN) estimated via Quantile Loss
- Histogram Estimation Network (HEN) estimated via MLE
    
All network models are estimated using stochastic gradient descent. Implementations are available in the directory: ?.

## Overview of results (?)
