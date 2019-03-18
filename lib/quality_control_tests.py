# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 18:54:16 2019

@author: r.realrangel
"""
import numpy as np
import xarray as xr


def normal(input_ts):
    return(xr.ufuncs.log(input_ts))


def standard(input_ts):
    return((input_ts - input_ts.mean()) / input_ts.std())


def zscore_check(input_ts, threshold, left_tail=True, right_tail=True):
    """ Test that data point exceeds min/max.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of the interest variable.
        threshold: float
            A threshold value to identifie outliers. It represents the
            deviation of the data point expressed as the number of
            times the standard deviation of the standardized time
            series.
            -----------------------------------------------------
            Freq (y)          Prob.     1T Z-score     2T Z-score
                                       (threshold)    (threshold)
            -----------------------------------------------------
                   2       0.500000        0.00000        0.67449
                   5       0.200000        0.84162        1.28155
                  10       0.100000        1.28155        1.64485
                  20       0.050000        1.64485        1.95996
                  50       0.020000        2.05375        2.32635
                 100       0.010000        2.32635        2.57583
                 500       0.002000        2.87816        3.09023
                1000       0.001000        3.09023        3.29053
               10000       0.000100        3.71902        3.89059
              100000       0.000010        4.26489        4.41717
             1000000       0.000001        4.75342        4.89164
            -----------------------------------------------------

    Reference
    ---------
        Gudmundsson, L., & Seneviratne, S. I. (2016). Observation-based
            gridded runoff estimates for Europe (E-RUN version 1.1).
            Earth System Science Data, 8(2), 279–295.
            https://doi.org/10.5194/essd-8-279-2016

        IOOS. (2018). Manual for Real-Time Quality Control of Stream
            Flow Observations: A Guide to Quality Control and Quality
            Assurance for Stream Flow Observations in Rivers and
            Streams.
    """
    lt = standard(normal(input_ts)) < -threshold
    rt = standard(normal(input_ts)) > threshold

    if left_tail and right_tail:
        return(lt | rt)

    elif left_tail:
        return(lt)

    elif right_tail:
        return(rt)


def range_test(input_ts, threshold=4.89164, climatology=True):
    """ Test that data point exceeds min/max.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
        threshold: float (default is 4.89164)
            A threshold value to identifie outliers. It represents the
            deviation of the data point expressed as the number of
            times the standard deviation of the standardized time
            series. The default value (4.89164) guarantees a
            probability of 0.999999 of "normal" values (and a probability
            of 0.0000001 of outliers).
        climatology: boolean
            Flag to specify wether to perform the test to the whole
            time series (gross test) or month by month (climatology test).

    Reference
    ---------
        Gudmundsson, L., & Seneviratne, S. I. (2016). Observation-based
            gridded runoff estimates for Europe (E-RUN version 1.1).
            Earth System Science Data, 8(2), 279–295.
            https://doi.org/10.5194/essd-8-279-2016

        IOOS. (2018). Manual for Real-Time Quality Control of Stream
            Flow Observations: A Guide to Quality Control and Quality
            Assurance for Stream Flow Observations in Rivers and
            Streams.
    """
    if climatology:
        return(input_ts.groupby('time.month').apply(func=zscore_check,
               threshold=threshold,
               left_tail=True,
               right_tail=True).drop('month'))

    else:
        return(zscore_check(input_ts=input_ts, threshold=threshold))


def spikes_data_test(input_ts, threshold=4.89164, climatology=True):
    """Test that data point n-1 exceeds a selected threshold relative
    to the adjacent data point.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
        threshold: float (default is 4.89164)
            A threshold value to identifie outliers. It represents the
            deviation of the data point expressed as the number of
            times the standard deviation of the standardized time
            series. The default value (4.89164) guarantees a
            probability of 0.999999 of "normal" values (and a probability
            of 0.0000001 of outliers).

    Reference
    ---------
        IOOS. (2018). Manual for Real-Time Quality Control of Stream
            Flow Observations: A Guide to Quality Control and Quality
            Assurance for Stream Flow Observations in Rivers and
            Streams.
    """
    spike_reference_weight = xr.DataArray([0.50, 0.0, 0.50], dims=['window'])
    spike_reference = input_ts.rolling(time=3, center=True).construct(
            'window').dot(spike_reference_weight)
    spikes = abs(input_ts - spike_reference)
    spikes[spikes == 0] = np.nan

    if climatology:
        return(spikes.groupby('time.month').apply(
                func=zscore_check,
                threshold=threshold,
                left_tail=False,
                right_tail=True).drop('month'))

    else:
        return(zscore_check(
                input_ts=spikes,
                threshold=threshold,
                left_tail=False,
                right_tail=True))


def change_rate_test(input_ts, threshold=4.89164, climatology=True):
    """Excessive rise/fall test.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
        threshold: float (default is 4.89164)
            A threshold value to identifie outliers. It represents the
            deviation of the data point expressed as the number of
            times the standard deviation of the standardized time
            series. The default value (4.89164) guarantees a
            probability of 0.999999 of "normal" values (and a probability
            of 0.0000001 of outliers).

    Reference
    ---------
        IOOS. (2018). Manual for Real-Time Quality Control of Stream
            Flow Observations: A Guide to Quality Control and Quality
            Assurance for Stream Flow Observations in Rivers and
            Streams.

        Krajewski, W. F., & Krajewski, K. L. (1989). Real-time quality
            control of streamflow data: A simulation study. Journal of
            the American Water Resources Association, 25(2), 391–399.
            https://doi.org/10.1111/j.1752-1688.1989.tb03076.x
    """
    change_rate = abs(input_ts.diff(dim='time'))
    change_rate = change_rate.reindex(time=sorted(input_ts.time.values))
    change_rate[change_rate == 0] = np.nan

    if climatology:
        return(change_rate.groupby('time.month').apply(
                func=zscore_check,
                threshold=threshold,
                left_tail=False,
                right_tail=True).drop('month'))

    else:
        return(zscore_check(
                input_ts=change_rate,
                threshold=threshold,
                left_tail=False,
                right_tail=True))


def flat_series_test(input_ts, value_tolerance=0.0, repetitions_tolerance=2,
                     skipzero=True):
    """ Invariant variable value.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
        tolerance: float (default is 0.0)
            Threshold of difference between two values to define if
            both are considered equivalent.
        repetitions_tolerance:
        skipzero:

    Reference
    ---------
        IOOS. (2018). Manual for Real-Time Quality Control of Stream
            Flow Observations: A Guide to Quality Control and Quality
            Assurance for Stream Flow Observations in Rivers and
            Streams.
    """
    invariant = (abs(input_ts.diff(dim='time')) <= value_tolerance)

    for recurse in range(1, repetitions_tolerance+1):
        invariant = invariant * (abs(
                input_ts.diff(
                        dim='time',
                        n=recurse)
                ) <= value_tolerance)

    invariant = invariant.reindex(
            time=sorted(input_ts.time.values)).astype('bool')

    for i in range(0, repetitions_tolerance + 1):
        input_ts.values[i] = False

    if skipzero:
        invariant[input_ts == 0] = False

    return(invariant)


def tmp_outlier_test(input_ts, c=7.5, threshold=4.89164):
    """ Applies the biweight mean and biweight standard deviation
    method to detect outliers where (i) data values are much larger (or
    smaller) than neighboring values but are not larger than the
    threshold to be detected by a consistency check (like the
    climatology test), known as 'spikes', and (2) data create a large
    step change from the previous daily value(s). Here, the method has
    been modified to compute mean and standard deviation of monthly
    groups, insteado of 3-daily groups as Lanzante (1996) originally
    proposed. This has been modified to consider the climatic
    variability.

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).

    Reference
    ---------
        Feng, S., Hu, Q., & Qian, W. (2004). Quality control of daily
            meteorological data in China, 1951–2000: a new dataset.
            International Journal of Climatology, 24(7), 853–870.
            https://doi.org/10.1002/joc.1047
        Gleason, B. E. (2002). Global Daily Climatology Network, V1.0
            (p. 26). Asheville, NC, USA: National Climatic Data Center.
        Lanzante, J. R. (1996). Resistant, robust and non-parametric
            techniques for the analysis of climate data: Theory and
            examples, including applications to historical radiosonde
            station data. International Journal of Climatology, 16(11),
            1197–1226.
            https://doi.org/10.1002/(SICI)1097-0088(199611)16:11<1197::
            AID-JOC89>3.0.CO;2-L


    """
    is_outlier = (input_ts.copy() * 0).astype('bool')
    dates = input_ts.time.values.astype('datetime64[D]')

    for date in dates:
        # Create a new time series (X_i) with the data from the day
        # before and after the day having the suspect value of the
        # current year, and the day before, after, and on the day from
        # all other years in the same station.
        X_i = xr.DataArray(data=[], coords={'time': []}, dims='time')

        for day in np.arange(date-15, date+16):
            is_the_month = input_ts['time.month'] == int(
                    str(day).split('-')[1])
            is_the_day = input_ts['time.day'] == int(str(day).split('-')[2])
            X_i = xr.concat(
                    [X_i, input_ts[(is_the_month) & (is_the_day)]], dim='time')

        X_i = X_i.reindex(time=sorted(X_i.time.values))
        X_i = X_i.drop(labels=date, dim='time')

        # After X_i series is obtained, the median (M) and absolute
        # deviatin from the median (MAD) are estimated. The MAD is the
        # median of the absolute deviations of the values from the
        # median.
        M = X_i.median()
        MAD = np.abs(X_i - M).median()

        # From the MAD, the weights u_i are calculated.
        u_i = (X_i - M) / (c * MAD)
        u_i[np.abs(u_i) > 1] = 1

        # With u_i, the biweight mean is estimated.
        upper = ((X_i - M) * ((1 - (u_i ** 2)) ** 2)).sum()
        lower = ((1 - (u_i ** 2)) ** 2).sum()
        Xmean_bi = M + (upper / lower)

        # And the biweight standard deviation.
        upper = np.sqrt(
                X_i.size * (((X_i - M) ** 2) * (1 - (u_i ** 2)) ** 4).sum())
        lower = np.abs(((1 - (u_i ** 2)) * (1 - (5 * (u_i ** 2)))).sum())
        s_bi = upper / lower

        # The Xmean_bi and s_bi are used to determine the Z-score of a
        # particular day's observation.
        Z = (input_ts.sel({'time': date}) - Xmean_bi) / s_bi
        is_outlier.loc[is_outlier.time == date] = abs(Z) >= threshold
        print("{} done".format(date))


def missd_ratio_test(input_ts, threshold=0.1):
    """ Test if the missing data in the time series is below (True) the
    defined threshold or if it is above (False).

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
        tolerance: float
            Threshold set up to missing data ratio (i. e.,
            missing_data / available_data). Values reported in the literature
            include 0.01, 0.05, 0.10, 0.15 and 0.20.

    Reference
    ---------
        Zhang, Y., & Post, D. (2018). How good are hydrological models
        for gap-filling streamflow data? Hydrology and Earth System
        Sciences, 22(8), 4593–4604.
        https://doi.org/10.5194/hess-22-4593-2018
    """
    gaps = float(input_ts.isnull().sum(dim='time'))
    return((gaps / input_ts.size) > threshold)


def minimlength_test(input_ts, threshold=10):
    """ Test if the length of the available records of the variable input_ts
    is above a defined threshold (in years).

    Parameters
    ----------
        input_ts: xarray.DataArray
            Time series of streamflow (or other variable).
    """
    pass
