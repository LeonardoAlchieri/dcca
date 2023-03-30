from numpy import apply_along_axis, arange, mean, ndarray, sqrt, stack, sum
from numpy.lib.stride_tricks import as_strided
from numpy.polynomial.polynomial import polyfit, polyval


def integrate_series(series: ndarray | list, mean_series: float) -> ndarray:
    """
    This method returns the integrated series. Given a series x, it will
    return the series x - mean(x).

    Parameters
    ----------
    series : ndarray | list
        timeseries over which to integrate
    mean_series: float
        mean of the series. Has to be given since it must be computed on
        the full series, and not the lagged ahead version

    Returns
    -------
    ndarray
        integrated time series, according to the formula reported above
    """
    return series - mean_series


def make_series_boxes(
    series: ndarray, time_scale: int, time_lag: int, N: int | None = None
) -> ndarray:
    """This method takes as an annary, the time scale (n parameter) and
    the time lag (τ) and gives an array with overlapping boxes. Each box

    Parameters
    ----------
    series : ndarray
        series to get overlapping boxes from
    time_scale : int
        time scale, n
    time_lag : int
        time lag, τ
    N : int
        number of points in the series. If None, it will be inferred from the `series` parameter

    Returns
    -------
    ndarray
        returns an array of size (N - time_scale - time_lag, time_scale + 1)
        containing the overlapping boxes
    """
    if N is None:
        N = len(series)
    return as_strided(
        series,
        shape=(N - time_scale - time_lag, time_scale + 1),
        strides=series.strides * 2,
    )


def get_local_detrend(integrated_data: ndarray, polynomial_degs: int = 1) -> ndarray:
    """This method takes a series box (subseries to be compared), computes
    a polynomial fit (by default linear regression) and then returns its
    values

    Parameters
    ----------
    series_box : ndarray
        series to be fittede
    polynomial_degs : int, optional
        degree of the polynomial fit, by default 1

    Returns
    -------
    ndarray
        polynomial values over the series range
    """

    x_in_box = arange(integrated_data.shape[1])
    coeffs = polyfit(x=x_in_box, y=integrated_data.T, deg=polynomial_degs)
    return polyval(x_in_box, coeffs)


def fddca_computation(
    boxes_residuals: ndarray,
    time_scale: int,
) -> float:
    """This method takes the residuals of the boxes (for both x and y) and computes the
    fdcca² value for that box. The input should be given with x and y residuals concatenated
    along the same dimension.

    This method assumes that the "time lag" has been already applied to the y residuals.

    Parameters
    ----------
    box residuals : ndarray
        residuals of the boxes, concatenated along the same dimension
    time_scale : int
        time scale, n

    Returns
    -------
    float
        the method returns the fdcca² value for the given box of x and y box
    """
    # NOTE: I am doing this reshaping because it is the only way to give
    # 2 dimensions to apply_along_axis

    # Reshape the boxes_residuals array to be 3D
    boxes_residuals = boxes_residuals.reshape(-1, 2, time_scale + 1)

    # Get the x_box_residuals and y_box_residuals arrays
    x_box_residuals = boxes_residuals[:, 0, :]
    y_box_residuals = boxes_residuals[:, 1, :]

    # Calculate the numerator of the FDDCA formula
    numerator = sum(x_box_residuals * y_box_residuals, axis=1)

    # Calculate the denominator of the FDDCA formula
    denominator = time_scale + 1

    # Calculate the FDCCA_squared value
    return sum(numerator) / denominator


def compute_FDCCA_squared(
    x: ndarray,
    y: ndarray,
    time_scale: int,
    time_lag: int,
    time_lag_y: int | None = None,
    N: int | None = None,
) -> float:
    """This method computes, for a given x and y time series, its FDCCA² value. Since
    the measure depends on a time scale and a time lag, these parameters are also required.

    As a suggestion, the time scale should be in the range of visible trends
    in the time series, since it is used to perform local detrends. The time
    lag should be investigated over different values, to gauge possible
    lagging effects between the series.

    Parameters
    ----------
    x : ndarray
        first array
    y : ndarray
        second array, to be lagged by time_lag
    time_scale : int
        time scale, i.e., the size of the boxes over which to perform the detrending
    time_lag : int
        time lag over which to shift the second array
    time_lag_y: int | None
        time lag over which to shift the second array. This should be used when x and y
        are the same signal and time_lag is given as 0. It allows to match the
        time lag of the convariance, and consider the same time series and
        performing analysis. Without its usage, results may be wrong.
    N : int, optional
        number of points in the series. If None, it will be inferred from the `series` parameter

    Returns
    -------
    float
        returns the FDCCA² value for the given time series and parameters
    """

    mean_x: float = mean(x)
    mean_y: float = mean(y)
    if N is None:
        N = len(x)
    # TODO: check that I am actually doing this correctly
    # NOTE: I rescale the second one since I have to move by tau
    # However, the rescalin will impact the calculation of the mean, as such
    # it is computed beforehand. This impact is larger the larger the
    # time lag is
    if time_lag_y is None:
        y = y[(time_lag):]
    elif time_lag_y > 0:
        y = y[(time_lag_y):]
        x = x[(time_lag_y):]
    # the size of the box is n+1 and the numbero of boxes is N - n - τ
    x_boxes = make_series_boxes(series=x, time_scale=time_scale, time_lag=time_lag, N=N)
    y_boxes = make_series_boxes(series=y, time_scale=time_scale, time_lag=time_lag, N=N)
    x_boxes_integrated: ndarray = x_boxes - mean_x.reshape(1, -1)
    y_boxes_integrated: ndarray = y_boxes - mean_y.reshape(1, -1)
    # x_boxes_integrated: ndarray = apply_along_axis(integrate_series, 1, x_boxes, mean_x)
    # y_boxes_integrated: ndarray = apply_along_axis(integrate_series, 1, y_boxes, mean_y)
    # x_local_detrend = apply_along_axis(get_local_detrend, 1, x_boxes_integrated)
    # y_local_detrend = apply_along_axis(get_local_detrend, 1, y_boxes_integrated)
    x_local_detrend = get_local_detrend(x_boxes_integrated)
    y_local_detrend = get_local_detrend(y_boxes_integrated)

    x_boxes_residuals = x_boxes_integrated - x_local_detrend
    y_boxes_residuals = y_boxes_integrated - y_local_detrend

    boxes_residuals = stack((x_boxes_residuals, y_boxes_residuals), axis=1)
    boxes_residuals = boxes_residuals.reshape(boxes_residuals.shape[0], -1)

    FDCCA_squared = fddca_computation(
        boxes_residuals=boxes_residuals, time_scale=time_scale
    ) / (N - time_scale - time_lag)
    # -0.0013385811760114692
    return FDCCA_squared


def detrended_correlation(
    x: ndarray, y: ndarray, time_scale: int, time_lag: int
) -> float:
    """This method computes the time lagged cross-correlation between two time series,
    x and y. The moethod is based on the paper "Analysis of detrended time-lagged
    cross-correlation between two nonstationary time series" by (Shen Chenhuaa, 2015).

    The output value, in range [-1,1], represents the detrended correlation, i.e.,
    a correlation without assuming any stationarity of the input series, between
    the two arrays. This value is dependent on the time scale chosen, which determines
    the size of the boxes over which to perform the detrending, and the time lag.
    In a Detrended Cross-Correlation Analysis, different time lags, and possibly
    different time scales, should be investigated.

    Parameters
    ----------
    x : ndarray
        first series
    y : ndarray
        second series
    time_scale : int
        time scale, over
    time_lag : int
        time lag, i.e., the shift to apply to the second series

    Returns
    -------
    float
        the detrended correlation between the two series
    """
    if len(x) != len(y):
        raise ValueError("The two series must have the same length")
    if not isinstance(x, ndarray) or not isinstance(y, ndarray):
        raise ValueError(
            "The two series must be numpy arrays. Got {} and {}".format(
                type(x), type(y)
            )
        )
    FDCCA_xy_squared = compute_FDCCA_squared(x, y, time_scale, time_lag)
    FDCCA_xx_squared = compute_FDCCA_squared(x, x, time_scale, 0, time_lag)
    FDCCA_yy_squared = compute_FDCCA_squared(y, y, time_scale, 0, time_lag)

    return FDCCA_xy_squared / (sqrt(FDCCA_xx_squared) * sqrt(FDCCA_yy_squared))
