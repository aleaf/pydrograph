import numpy as np
import pandas as pd
import gisutils


def get_upstream_area(points, PlusFlow, NHDFlowlines, NHDCatchments, nearfield=None):
    """For each point in points, get upstream drainage area in km2, using
    NHDPlus PlusFlow routing table and NHDPlus Catchment areas. Upstream area
    within the containing catchment is estimated as a fraction of proportional
    to the distance of the measurment point along the NHDPlus Flowline associated with the catchment.

    Parameters
    ----------
    points : list of shapely Point objects
        Locations of streamflow measurements. Must be in same coordinate system as NHDCatchments
    PlusFlow : str or list of strings
        Path(s) to PlusFlow routing tables
    NHDFlowlines : str or list of strings
        Path(s) to Flowlines shapefiles
    NHDCatchments : str or list of strings
        Path(s) to Catchment shapefiles
    nearfield : shapefile or shapely Polygon
        Nearfield area of model. Used to filter NHDPlus flowlines and catchments to
        greatly speed reading them in and finding the COMIDs associated with points.
        Must be in same coordinate system as points and NHDPlus shapefiles.

    Returns
    -------
    upstream_area : list
        List of areas in km2, for each point in points.
    """
    try:
        import fiona
        from shapely.geometry import LineString, Polygon, shape
        from GISio import shp2df
    except ImportError:
        print('This method requires fiona, shapely and GIS_utils.')

    if isinstance(nearfield, Polygon):
        bbox = nearfield.bounds
    elif isinstance(nearfield, str):
        bbox = shape(fiona.open(nearfield).next()['geometry']).bounds()
    else:
        bbox = None

    # dialate the bounding box by half, so that features aren't missed.
    x = 0.5 * (bbox[2] - bbox[0])
    y = 0.5 * (bbox[3] - bbox[1])
    bbox = (bbox[0]-x, bbox[1]-y, bbox[2]+x, bbox[3]+y)

    pf = shp2df(PlusFlow)
    fl = shp2df(NHDFlowlines, index='COMID', filter=bbox)
    cmt = shp2df(NHDCatchments, index='FEATUREID', filter=bbox)

    # find the catchment containing each point in points
    comids = []
    for p in points:
        comids += cmt.FEATUREID[np.array([p.within(g) for g in cmt.geometry])].tolist()

    upstream_area = []
    for i, comid in enumerate(comids):
        comids = {comid}
        upstream = [comid]
        for j in range(1000):
            upstream = set(pf.ix[pf.TOCOMID.isin(upstream), 'FROMCOMID']).difference({0})
            if len(upstream) == 0:
                break
            comids.update(upstream)

        total_upstream_area = cmt.ix[comids, 'AreaSqKM'].sum()
        if comid == 11951607:
            j=2
        # estimate fraction of containing catchment that is upstream
        # by finding closest vertex on flowline,
        # and then dividing upstream length by downstream length
        #X = np.array(fl.ix[comid, 'geometry'].coords.xy[0])
        #Y = np.array(fl.ix[comid, 'geometry'].coords.xy[1])
        g = points[i] # misc measurement point
        #i = np.argmin(np.sqrt((X-g.x)**2 + (Y-g.y)**2)) # closest point on flowline

        # should be able to just project point onto flowline and divide by total length
        l = fl.ix[comid, 'geometry']
        frac = l.project(g)/l.length
        #frac = LineString(zip(X[:i+1], Y[:i+1])).length/LineString(zip(X[i:], Y[i:])).length
        upstream_in_catchment = cmt.ix[comid, 'AreaSqKM'] * frac
        total_upstream_area += upstream_in_catchment
        upstream_area.append(total_upstream_area)

    return upstream_area


def IHmethod(Qseries, block_length=5, tp=0.9, interp_semilog=True, freq='D', limit=100):
    """Baseflow separation using the Institute of Hydrology method, as documented in
    Institute of Hydrology, 1980b, Low flow studies report no. 3--Research report: 
    Wallingford, Oxon, United Kingdom, Institute of Hydrology Report no. 3, p. 12-19,
    and
    Wahl, K.L and Wahl, T.L., 1988. Effects of regional ground-water level declines
    on streamflow in the Oklahoma Panhandle. In Proceedings of the Symposium on 
    Water-Use Data for Water Resources Management, American Water Resources Association. 
    
    Parameters
    ----------
    Qseries : pandas Series
        Pandas time series (with datetime index) containing measured streamflow values.
    block_length : int
        N parameter in IH method. Streamflow is partitioned into N-day intervals;
        a minimum flow is recorded for each interval.
    tp : float
        f parameter in IH method. For each three N-day minima, if f * the central value
        is less than the adjacent two values, the central value is considered a 
        turning point. Baseflow is interpolated between the turning points.
    interp_semilog : boolean
        If False, linear interpolation is used to compute baseflow between  turning points
        (as documented in the IH method). If True, the base-10 logs of the turning points
        are interpolated, and the interpolated values are transformed back to 
        linear space (producing a curved hydrograph). Semi-logarithmic interpolation
        as documented in Wahl and Wahl (1988), is used in the Base-Flow Index (BFI)
        fortran program. This method reassigns zero values to -2 in log space (0.01)
        for the interpolation.
    freq : str or DateOffset, default ‘D’
        Any `pandas frequency alias <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases>`_
        Regular time interval that forms the basis for base-flow separation. Input data are
        resampled to this frequency, and block lengths represent the number of time increments
        of the frequency. By default, days ('D'), which is what all previous BFI methods
        are based on. Note that this is therefore an experimental option; it is up to the user t
        o verify any results produced by other frequencies.
    limit : int
        Maximum number of timesteps allowed during linear interploation between baseflow 
        ordinances. Must be greater than zero.

    
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
    if len(Qseries) < 2 * block_length:
        raise ValueError('Input Series must be at '
                         'least two block lengths\nblock_length: '
                         '{}\n{}'.format(block_length, Qseries))

    # convert flow values to numeric if they are objects
    # (pandas will cast column as objects if there are strings such as "ICE")
    # coerce any strings into np.nan values
    if Qseries.dtype.name == 'object':
        Qseries = pd.to_numeric(Qseries, errors='coerce')

    # convert the series to a dataframe; resample to daily values
    # missing days will be filled with nan values
    df = pd.DataFrame(Qseries).resample(freq).mean()
    df.columns = ['Q']

    # compute block numbers for grouping values on blocks
    nblocks = int(np.floor(len(df) / float(block_length)))

    # make list of ints, one per measurement, denoting the block
    # eg [1,1,1,1,1,2,2,2,2,2...] for block_length = 5
    n = []
    for i in range(nblocks):
        n += [i + 1] * block_length
    n += [np.nan] * (len(df) - len(n))  # pad any leftover values with nans
    df['n'] = n

    # compute the minimum for each block
    # create dataframe Q, which only has minimums for each block
    blocks = df[['Q', 'n']].reset_index(drop=True).dropna(axis=0).groupby('n')
    Q = blocks.min()
    Q = Q.rename(columns={'Q': 'block_Qmin'})
    Q['n'] = Q.index
    # get the index position of the minimum Q within each block
    idx_Qmins = blocks.idxmin()['Q'].values.astype(int)
    # get the date associated with each Q minimum
    Q['datetime'] = df.index[idx_Qmins]

    # compute baseflow ordinates
    Q['ordinate'] = [np.nan] * len(Q)
    Qlist = Q.block_Qmin.tolist()
    Q['Qi-1'] = [np.nan] + Qlist[:-2] + [np.nan]
    Q['Qi'] = [np.nan] + Qlist[1:-1] + [np.nan]
    Q['Qi+1'] = [np.nan] + Qlist[2:] + [np.nan]
    isordinate = tp * Q.Qi < Q[['Qi-1', 'Qi+1']].min(axis=1)
    Q.loc[isordinate, 'ordinate'] = Q.loc[isordinate, 'block_Qmin']

    # reset the index of Q to datetime
    Q.index = Q.datetime

    # expand Q dataframe back out to include row for each day
    Q = Q.dropna(subset=['datetime'], axis=0).resample(freq).mean()

    # interpolate between baseflow ordinates
    if interp_semilog:
        iszero = Q.ordinate.values == 0
	logQ = Q.ordinate.copy()
        logQ[iszero] = -2
        logQ = np.log10(Q.ordinate)
        QB = np.power(10.0, logQ.interpolate(limit=limit).values)
    else:
        QB = Q.ordinate.interpolate(limit=limit).values
    Q['QB'] = QB

    # reassign the original flow values back to Q
    Q['Q'] = df.Q.loc[Q.index]

    # ensure that no baseflow values are > Q measured
    QBgreaterthanQ = Q.QB.values > Q.Q.values
    Q.loc[QBgreaterthanQ, 'QB'] = Q.loc[QBgreaterthanQ, 'Q']
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


def baseflow_summary(self, field_sites, field_measurements, daily_values, q90_window=20, output_proj4=None):

    fm = field_measurements
    dvs = daily_values

    if fm['measurement_dt'].dtype != 'datetime64[ns]':
        fm['measurement_dt'] = pd.to_datetime(fm.measurement_dt)

    # reprojected the output X, Y coordinates
    print('reprojecting output from\n{}\nto\n{}...'.format(self.proj4, output_proj4))
    if output_proj4 is not None:
        field_sites['geometry'] = gisutils.project(field_sites, self.proj4, output_proj4)

    fm_site_no = []
    Qm = []
    measurement_dt = []
    measured_rating_diff = []
    drainage_area = []
    station_nm = []
    index_station = []
    indexQr = []
    indexQ90 = []
    X, Y = [], []
    for i in range(len(fm)):
        mdt = fm.measurement_dt.tolist()[i]
        Dt = dt.datetime(mdt.year, mdt.month, mdt.day)
        for site_no, data in list(dvs.items()):

            # check if index station covers measurement date
            try:
                dv = data.ix[Dt]
            except KeyError:
                continue
            dv = data.ix[Dt]
            site_no = dv.site_no
            DDcd = [k for k in list(data.keys()) if '00060' in k and not 'cd' in k][0]
            try:
                Qr = float(dv[DDcd]) # handle ice and other non numbers
            except:
                continue

            # get q90 values for window
            q90start = pd.Timestamp(Dt) - pd.Timedelta(0.5 * q90_window, unit='Y')
            q90end = pd.Timestamp(Dt) + pd.Timedelta(0.5 * q90_window, unit='Y')
            values = pd.to_numeric(data.ix[q90start:q90end, DDcd], errors='coerce')
            q90 = values.quantile(q=0.1)

            # append last to avoid mismatches in length
            site_info = field_sites.ix[fm.site_no.values[i]]
            fm_site_no.append(fm.site_no.values[i])
            station_nm.append(site_info['station_nm'])
            Qm.append(fm.discharge_va.values[i])
            measurement_dt.append(fm.measurement_dt.tolist()[i])
            measured_rating_diff.append(fm.measured_rating_diff.values[i])
            drainage_area.append(site_info['drain_area_va'])
            index_station.append(site_no)
            indexQr.append(Qr)
            indexQ90.append(q90)
            X.append(site_info['geometry'].xy[0][0])
            Y.append(site_info['geometry'].xy[1][0])

    df = pd.DataFrame({'site_no': fm_site_no,
                       'station_nm': station_nm,
                       'datetime': measurement_dt,
                       'Qm': Qm,
                       'quality': measured_rating_diff,
                       'drn_area': drainage_area,
                       'idx_station': index_station,
                       'indexQr': indexQr,
                       'indexQ90': indexQ90,
                       'X': X,
                       'Y': Y})
    df['est_error'] = [self.est_error.get(q.lower(), self.default_error) for q in df.quality]
    df = df[['site_no', 'datetime', 'Qm', 'quality', 'est_error',
             'idx_station', 'indexQr', 'indexQ90', 'drn_area', 'station_nm', 'X', 'Y']]
    return df