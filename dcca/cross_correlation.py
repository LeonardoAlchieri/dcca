from numpy import mean, ndarray, sqrt

from dcca import integrate_series


def cross_covariance(x: ndarray, y: ndarray, time_lag: int) -> float:
    mean_x = mean(x)
    mean_y = mean(y)
    # NOTE: second series must be shifted by time_lag
    y = y[time_lag:]

    integrated_x = integrate_series(x, mean_x)
    # NOTE: the formula does not use the last N-τ elements of the first series
    integrated_x = integrated_x[: len(integrated_x) - time_lag]
    # FIXME: I believe mean should be calculated before removing the first
    # time_lag elements from the second series. However, I'm not 100% sure
    integrated_y = integrate_series(y, mean_y)

    return (integrated_x * integrated_y).sum()


def cross_correlation(x: ndarray, y: ndarray, time_lag: int) -> float:
    C_xy = cross_covariance(x, y, time_lag)
    C_xx = cross_covariance(x, x, 0)
    C_yy = cross_covariance(y, y, 0)

    return C_xy / (sqrt(C_xx) * sqrt(C_yy))
