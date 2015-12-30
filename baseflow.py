import numpy as np
import pandas as pd

def IHmethod(values, block_length=5, tp=0.9, interp_semilog=True):
    """Baseflow separation using the Institute of Hydrology method, as documented in
    Institute of Hydrology, 1980b, Low flow studies report no. 3--Research report: 
    Wallingford, Oxon, United Kingdom, Institute of Hydrology Report no. 3, p. 12-19,
    and
    Wahl, K.L and Wahl, T.L., 1988. Effects of regional ground-water level declines
    on streamflow in the Oklahoma Panhandle. In Proceedings of the Symposium on 
    Water-Use Data for Water Resources Management, American Water Resources Association. 
    
    Parameters
    ----------
    values : pandas Series
        Pandas time series (with datetime index) containing measured streamflow values.
    block_length : int
        N parameter in IH method. Streamflow is partitioned into N-day intervals;
        a minimum flow is recorded for each interval.
    tp : float
        f parameter in IH method. For each three N-day minima, if f * the central value
        is less than the adjacent two values, the central value is considered a 
        turning point. Baseflow is interpolated between the turning points.
    interp_semilog : boolean
        If False, linear interpolation is used to compute baseflow between turning points
        (as documented in the IH method). If True, the base-10 logs of the turning points
        are interpolated, and the interpolated values are transformed back to 
        linear space (producing a curved hydrograph). Semi-logarithmic interpolation
        as documented in Wahl and Wahl (1988), is used in the Base-Flow Index (BFI)
        fortran program.
    
    Returns
    -------
    Q : pandas DataFrame
        DataFrame containing the following columns:
        minima : N-day minima
        ordinate : selected turning points
        n : block number for each N-day minima
        QB : computed baseflow
        Q : discharge values
    
    Notes
    -----
    Whereas this program only selects turning points following the methodology above, 
    the BFI fortran program adds artificial turning points at the start and end of
    each calendar year. Therefore results for datasets consisting of multiple years
    will differ from those produced by the BFI program.
    
    """
    if values.dtype.name == 'object':
        values = values.convert_objects(convert_numeric=True)
    values = pd.DataFrame(values).resample('D')
    values.columns = ['discharge']

    # compute block numbers for grouping values on blocks
    #nblocks = int(len(values) / float(block_length) + 1)
    nblocks = int(np.floor(len(values) / float(block_length)))
    n = []
    for i in range(nblocks):
        n += [i+1] * block_length
    #values['n'] = n[:len(values)]
    n += [np.nan] * (len(values) - len(n))
    values['n'] = n
    
    # compute minima for block_length day blocks
    Q = [np.min(values.discharge.values[i-block_length:i]) 
         for i in np.arange(block_length, len(values))[::block_length]]

    #  compute the minimum for each block
    Q = values.groupby('n').min()
    Q['datetime'] = values[['discharge', 'n']].groupby('n').idxmin() # include dates of minimum values

    Q['ordinate'] = [np.nan] * len(Q)
    #for i in range(len(Q))[:-2]:
    #    end1, cv, end2 = Q.discharge[i:i+3]
    for i in np.arange(1, len(Q)-1):
        end1, cv, end2 = Q.discharge.values[i-1:i+2]
        if tp * cv < end1 and tp * cv < end2:
            Q.loc[Q.index[i], 'ordinate'] = cv
    Q['n'] = Q.index
    Q.index = Q.datetime

    Q = Q.dropna(subset=['datetime'], axis=0).resample('D')
    if interp_semilog:
        QB = 10**(np.log10(Q.ordinate).interpolate()).values
    else:
        QB = Q.ordinate.interpolate(limit=100).values
    Q['QB'] = QB
    Q['Q'] = values.discharge[Q.index]
    QBgreaterthanQ = Q.QB.values > Q.Q.values
    Q.loc[QBgreaterthanQ, 'QB'] = Q.ix[QBgreaterthanQ, 'Q']
    return Q

def WI_statewide_eqn(Qm, A, Qr, Q90):
    """Regression equation of Gebert and others (2007, 2011)
    for estimating average annual baseflow from a field measurement of streamflow
    during low-flow conditions.

    Parameters
    ----------
    Qm : float or 1-D array of floats
        Measured streamflow.
    A : float or 1-D array of floats
        Drainage area in watershed upstream of where Qm was taken.
    Qr : float or 1-D array of floats
        Recorded flow at index station when Qm was taken.
    Q90 : float or 1-D array of floats
        Q90 flow at index station.

    Returns
    -------
    Qb : float or 1-D array of floats, of length equal to input arrays
        Estimated average annual baseflow at point where Qm was taken.
    Bf : float or 1-D array of floats, of length equal to input arrays
        Baseflow factor. see Gebert and others (2007, 2011).

    Notes
    -----
    Gebert, W.A., Radloff, M.J., Considine, E.J., and Kennedy, J.L., 2007,
    Use of streamflow data to estimate base flow/ground-water recharge for Wisconsin:
    Journal of the American Water Resources Association,
    v. 43, no. 1, p. 220-236, http://dx.doi.org/10.1111/j.1752-1688.2007.00018.x

    Gebert, W.A., Walker, J.F., and Kennedy, J.L., 2011,
    Estimating 1970-99 average annual groundwater recharge in Wisconsin using streamflow data:
    U.S. Geological Survey Open-File Report 2009-1210, 14 p., plus appendixes,
    available at http://pubs.usgs.gov/ofr/2009/1210/.
    """
    Bf = (Qm / A) * (Q90 / Qr)
    Qb = 0.907 * A**1.02 * Bf**0.52
    return Qb.copy(), Bf.copy()