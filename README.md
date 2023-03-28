[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/LeonardoAlchieri/dcca/graphs/commit-activity)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://github.com/LeonardoAlchieri/dcca/blob/main/LICENSE)
[![PyPI pyversions](https://img.shields.io/badge/Python-3.11-informational)](https://github.com/LeonardoAlchieri/dcca)

# Time-Lagged Detrended Cross Correlation Analysis

This package allows to implement Detrended Cross Correlation Analysis [1] and its Time-Lagged version. Implementations are based on [1] and [2].

The code allows to easily calculate the **Detrended Cross-Correlation Coefficient**, which is a variation of the more famous Pearson's correlation coefficient, meant for cases in which the time series is non stationary. Indeed, in many domains, stationarity of a time series is not guaranteed and running traditional correlation coefficients can give skewed or untrue results. 

## Usage
The package allows to calculate the detrended version with both no time lag, or with a given time lag.
```python
from dcca import detrended_correlation
from numpy.random import rand

x: ndarray = rand(100)
y: ndarray = rand(100)

print(detrended_correlation(x=x, y=y, time_scale=3, time_lag=0))
```
```
-0.04848682827863634
```
The `time_lag` coefficient is meant to confront the time series when one is moved with respect to the other. At the moment, only the time series `y` can be moved using this. To run the full cross-correlation analysis, one may want to loop over different time lags:
```python
from dcca import detrended_correlation
from numpy.random import rand

x: ndarray = rand(100)
y: ndarray = rand(100)
time_lags: list[int] = range(0, 10)
dccas = [detrended_correlation(x=x, y=y, time_scale=3, time_lag=0) for time_lag in time_lags]
```
In the example provided, the coefficient is always going to be the same, since the arrays are randomly sampled.

The `time_scale` coefficient specifies the box size of the local detrending. In short, the Detrended Cross-Correlation Coefficient, in order to account for non-stationarity in the two series, calculates the correlation over smaller parts of the series (boxes), whose size is indeed given by the aforementioned parameter. For more detail, please look at the referenced material.

The package also allows to calculate the standard Pearson's correlation coefficient, even with different time lags. While the method can be easily implemented, I though it might be useful to provide it for confronttions.
```python
from dcca.cross_correlation import cross_correlation
from numpy.random import rand
x: ndarray = rand(100)
y: ndarray = rand(100)

print(cross_correlation(x=x, y=y, time_lag=0))
```

## Installation
At the moment, the packe is only installable via `pip`. For following releases, I will try to port it into `conda` as well.
```bash
pip install dcca
```

---
## References
[1] Podobnik, Boris, and H. Eugene Stanley. "Detrended cross-correlation analysis: a new method for analyzing two nonstationary time series." Physical review letters 100.8 (2008): 084102.

[2] Shen, Chenhua. "Analysis of detrended time-lagged cross-correlation between two nonstationary time series." Physics Letters A 379.7 (2015): 680-687.

## Contacts
For any information, contact me, Leonardo Alchieri, at leonardo@alchieri.eu. This package was developed as part of my PhD at USI (Università della Svizzera italiana), Switzerland.