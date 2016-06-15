__author__ = 'aleaf'

import datetime as dt
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, shape
import pyproj
import GISio
from GISops import project, projectdf
from .attributes import streamflow_attributes, gw_attributes


coord_datums_epsg = {'NAD83': 4269,
                     'NAD27': 4267}

coord_datums_proj4 = {'NAD83': '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs',
                     'NAD27': '+proj=longlat +ellps=clrk66 +datum=NAD27 +no_defs'}

def WI_statewide_eqn(Qm, A, Qr, Q90):
    Bf = (Qm / A) * (Q90 / Qr)
    Qb = 0.907 * A**1.02 * Bf**0.52
    return Qb.copy(), Bf.copy()


class NWIS:
    """
    NWIS error codes:
    E     Excellent    The data is within 2% (percent) of the actual flow
    G     Good         The data is within 5% (percent) of the actual flow
    F     Fair         The data is within 8% (percent) of the actual flow
    P     Poor         The data are >8% (percent) of the actual flow
    """

    est_error = {'EXCELLENT': 0.02,
                 'GOOD': 0.05,
                 'FAIR': 0.08}
    default_error = 0.50

    urlbase = 'http://nwis.waterdata.usgs.gov/usa/nwis/'
    dtypes_dict = {'dv': 'dv?referred_module=sw&site_tp_cd=ST&',
                   'field_measurements': 'measurements?',
                   'gwlevels': 'gwlevels?',
                   'gwdv': 'dv?referred_module=gw&site_tp_cd=GW&',
                   'inventory': 'inventory?'}

    parameter_codes = {'discharge': '00060',
                       'gwlevels': '72019'}

    coordinate_format = 'decimal_degrees' #coordinate_format=decimal_degrees&
    group_key = 'NONE' #group_key=NONE&
    output_format = 'sitefile_output' #format=sitefile_output&
    sitefile_output_format = 'rdb' #sitefile_output_format=rdb&

    now = dt.datetime.now()
    range_selection = 'days' #range_selection=days
    period = 365 #period=365
    begin_date = '1880-01-01' #begin_date=2014-04-14
    end_date = '{:02d}-{:02d}-{:02d}'.format(now.year, now.month, now.day-1) #end_date=2015-04-13

    date_cols = ['measurement_dt', 'lev_dt']

    logscale = 1 #'set_logscale_y=1'
    channel_html_info = 0 #'channel_html_info=0'
    date_format = 'YYYY-MM-DD' #'date_format=YYYY-MM-DD'
    channel_rdb_info = 0 #'channel_rdb_info=0'
    rdb_compression = 'file' #'rdb_compression=file'
    list_of_search_criteria = 'lat_long_bounding_box' #'list_of_search_criteria=lat_long_bounding_box'

    def __init__(self, ll_bbox=None, extent=None, datum='NAD83',
                 get_sw_sites=True, get_gw_sites=True):
        """Class for retrieving data from NWIS.
        Currently only retrieves data within a specified lat/lon bounding box.

        Parameters
        ----------
        ll_bbox: list of floats
            List containing bounding latitudes and longitudes in decimal degrees, in the following order:
            [northwest longitude, northwest latitude, southeast longitude, southeast latitude]


        see the code for a list of default parameters
        """

        self.ll_bbox = ll_bbox
        self.datum = datum
        self.proj4 = coord_datums_proj4[self.datum]
        if extent is not None:
            if isinstance(extent, str):
                self.extent = self._read_extent_shapefile(extent)
            elif isinstance(extent, Polygon):
                self.extent = extent
            else:
                print('Warning: extent argument of unknown datatype!')
                self.extent = None
        else:
            self.extent = None
        if ll_bbox is None and self.extent is not None:
            self.ll_bbox = self.extent.bounds
        elif not ll_bbox and self.extent is None:
            print('Need bounding box or shapefile for NWIS queries.')
            return
        else:
            pass

        if get_gw_sites or get_sw_sites:
            print('Fetching site info...')
        if get_sw_sites:
            self.field_sites = self.get_siteinfo('field_measurements', streamflow_attributes)
            self.dv_sites = self.get_siteinfo('dv', streamflow_attributes)
        if get_gw_sites:
            self.gwfield_sites = self.get_siteinfo('gwlevels', gw_attributes)
            self.gwdv_sites = self.get_siteinfo('gwdv', gw_attributes)

        self.field_measurements = pd.DataFrame() # dataframe of all field measurements for area
        self.gwlevels = pd.DataFrame()
        self.dvs = {} # dictionary with dataframes of daily values for all dv sites, keyed by site no
        self.dv_q90 = {} # q90 flows for daily values stations, keyed by site no

    def _compute_geometries(self, df):

        geoms = []
        for i in range(len(df)):

            p = Point(df.dec_long_va[i], df.dec_lat_va[i])
            pr1 = "+init=EPSG:{}".format(coord_datums_epsg[df.dec_coord_datum_cd[i]])
            pr2 = self.proj4
            geom = project(p, pr1, pr2)
            geoms.append(geom)
        return geoms

    def _cull_to_extent(self, df):

        if not 'geometry' in df.columns:
            df['geometry'] = self._compute_geometries(df)

        within = np.array([g.within(self.extent) for g in df.geometry])
        return df[within].copy()

    def _read_extent_shapefile(self, shpfile, buffer=0):

        import fiona
        from fiona.crs import to_string, from_epsg

        print('reading extent from {}...'.format(shpfile))
        shp = fiona.open(shpfile)
        g = shape(shp.next()['geometry'])

        if to_string(from_epsg(coord_datums_epsg[self.datum])) != to_string(shp.crs):
            print('reprojecting extent from {} to {}'.format(to_string(shp.crs), self.proj4))
            return project(g, to_string(shp.crs), self.proj4)
        else:
            return g

    def make_site_url(self, data_type='inventory', attributes=None):
        """
        Parameters
        ----------
        data_type: str
            'dv' for Daily Values
            'field_measurements' for Field Measurements
            'inventory' for all measurements

        Returns
        -------
        url string
        """
        self.bbox_url = 'nw_longitude_va={:.3f}&'.format(self.ll_bbox[0]) +\
                        'nw_latitude_va={:.3f}&'.format(self.ll_bbox[3]) +\
                        'se_longitude_va={:.3f}&'.format(self.ll_bbox[2]) +\
                        'se_latitude_va={:.3f}&'.format(self.ll_bbox[1])

        self.stuff_at_beginning = 'coordinate_format={}&'.format(self.coordinate_format) +\
                                  'group_key={}&'.format(self.group_key) +\
                                  'format={}&'.format(self.output_format) +\
                                  'sitefile_output_format={}&'.format(self.sitefile_output_format)

        self.dv_info = 'range_selection={}&'.format(self.range_selection) +\
                        'period={}&'.format(self.period) +\
                        'begin_date={}&'.format(self.begin_date) +\
                        'end_date={}&'.format(self.end_date)

        self.stuff_at_end = 'date_format={}&'.format(self.date_format) +\
                            'rdb_compression={}&'.format(self.rdb_compression) +\
                            'list_of_search_criteria={}'.format(self.list_of_search_criteria)

        url = self.urlbase + self.dtypes_dict.get(data_type, data_type+'?')
        url += self.bbox_url
        url += self.stuff_at_beginning
        if attributes is not None:
            for a in attributes:
                url += 'column_name=' + a + '&'

        if data_type == 'dv':
            url += self.dv_info

        url += self.stuff_at_end
        #print '{}'.format(url)
        return url

    def make_dv_url(self, station_IDs, parameter_code='00060', start_date='1880-01-01', end_date=None):
        """Creates url to retrieve daily values for a site


        Parameters
        ----------
        stationIDs: int, str or list of ints or strings
            USGS station IDs

        parameter_code: (string)
            e.g. 00060 for discharge.
            See http://help.waterdata.usgs.gov/codes-and-parameters/parameters.

        start_date: (string) 'YYYY-DD-MM'
            To obtain the entire period-of-record use a start date of 1880-01-01 (default)...

        Notes
        -----
        A leading zero is added to the site number if the first digit is greater than 1
        (this can happend for basins 01 - 09 if the site number gets converted to an int).
        Note that this may cause site numbers for basin 01 (North Atlantic slope) to get confused with
        basins 10-16 (west coast and hawaii).
        See <http://help.waterdata.usgs.gov/faq/sites/do-station-numbers-have-any-particular-meaning>

        """

        if not isinstance(station_IDs, list):
            station_IDs = [str(station_IDs)]

        def add_leading_zero(station_ID):
            if 1 < int(str(station_ID)[0]) < 10:
                station_ID = '0{}'.format(station_IDs)
            return station_ID

        #station_IDs = ','.join(['0{}'.format(int(str(s))) for s in station_IDs])
        station_IDs = ','.join([self._correct_stationID(s) for s in station_IDs])

        url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb'

        url += '&sites={}'.format(station_IDs)
        url += '&startDT={}'.format(start_date)
        if end_date is not None:
            url += '&endDT={}'.format(end_date)
        url += '&parameterCd={}'.format(parameter_code)
        print('{}'.format(url))
        return url


    def make_measurements_url(self, station_ID, txt='measurements'):
        """Creates url to retrieve daily values for a site


        Parameters
        ----------
        stationID: (string)
            USGS station ID

        txt: (string)
            'measurements' for field measurements of streamflow
            'gwlevels' for field measurements of groundwater level

        """
        station_ID = self._correct_stationID(station_ID)

        url =  'http://nwis.waterdata.usgs.gov/nwis/{}?site_no={}&agency_cd=USGS&format=rdb'\
                .format(txt, station_ID)
        print('{}'.format(url))
        return url

    def get_header_length(self, sitefile_text, col0):
        knt = 0
        for line in sitefile_text:
            if '#' not in str(line) and col0 in str(line):
                knt += 2
                break
            else:
                knt += 1
        return knt

    def get_siteinfo(self, data_type, attributes):
        """Retrieves site information for the bounding box supplied to the NWIS class instance

        Parameters
        ----------
        data_type: str
            'dv' for Daily Values
            'field_measurements' for Field Measurements

        attributes: list of strings
            List of NWIS attributes to include (e.g. 'site_no', 'station_nm', etc.)
            Default sets of attributes for streamflow and groundwater levels data can
            be imported from the attributes.py file (work in progress)

        Returns
        -------
        the contents of an NWIS site information file in a dataframe format
        """
        url = self.make_site_url(data_type, attributes)
        sitefile_text = urlopen(url).readlines()
        skiprows = self.get_header_length(sitefile_text, attributes[0])
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=attributes)

        df['geometry'] = self._compute_geometries(df)
        df.index = df.site_no
        if self.extent is not None:
            within = np.array([g.within(self.extent) for g in df.geometry])
            return df[within].copy()
        return df

    def _correct_stationID(self, stationID):
        if 1 < int(str(stationID)[0]) < 10 and len(str(stationID)) < 15:
            return '0{}'.format(stationID)
        return str(stationID)

    @property
    def _get_dv_sites(self):
        print('Fetching info for sites with daily values...')
        self.dv_sites = self.get_siteinfo('dv', streamflow_attributes)

    def _get_date_col(self, df):
        return [d for d in self.date_cols if d in df.columns][0]

    def get_dvs(self, station_ID, parameter_code='00060', start_date='1880-01-01', end_date=None):
        """Retrieves daily values for a site.

        Parameters
        ----------
        stationID: (string)
            USGS station ID

        parameter_code: string, default is 00060 for discharge.
            See http://help.waterdata.usgs.gov/codes-and-parameters/parameters.

        start_date: (string) 'YYYY-DD-MM'
            To obtain the entire period-of-record use a start date of 1880-01-01 (default)...

        Returns
        -------
        dv: a datetime-index dataframe of daily discharge, with datagaps filled with NaNs
        """
        if parameter_code in list(self.parameter_codes.keys()):
            parameter_code = self.parameter_codes[parameter_code]

        url = self.make_dv_url(station_ID, parameter_code=parameter_code,
                               start_date=start_date, end_date=end_date)
        sitefile_text = urlopen(url).readlines()
        skiprows = self.get_header_length(sitefile_text, 'agency_cd')
        cols = sitefile_text[skiprows - 2].decode('utf-8').strip().split('\t')
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=cols)
        df.index = pd.to_datetime(df.datetime)
        return df

    def get_measurements(self, station_ID, txt='measurement'):
        """Retrieves field measurements for a site.

        Parameters
        ----------
        stationID: (string)
            USGS station ID

        Returns
        -------
        dv: a datetime-index dataframe of the measurements
        """

        url = self.make_measurements_url(station_ID, txt)
        sitefile_text = urlopen(url).readlines()
        skiprows = self.get_header_length(sitefile_text, 'agency_cd')
        cols = sitefile_text[skiprows - 2].decode('utf-8').strip().split('\t')
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=cols)
        if len(df) > 0:
            df.index = pd.to_datetime(df[self._get_date_col(df)])
        return df

    def get_all_measurements(self, site_numbers, txt='measurements'):
        """Get measurements for a list of site numbers.

        Parameters
        ----------
        site_numbers : list or 1D array
            USGS site numbers
        txt : str
            String in url specifying type of measurement
            measurements : field measurements
            dv : daily values
            gwlevels : groundwater levels
            qwdata : water quality data

        """
        all_measurements = pd.DataFrame()
        for s in site_numbers:
            print(s)
            df = self.get_measurements(s, txt=txt)
            df.index = pd.MultiIndex.from_product([[df.site_no.values[0]], df.index.values],
                                              names=['site_no', 'datetime'])
            df['measurement_dt'] = pd.to_datetime(df[self._get_date_col(df)])
            all_measurements = all_measurements.append(df)
        #self.field_measurements = all_measurements
        return all_measurements

    def get_all_dvs(self, stations, parameter_code='00060', start_date='1880-01-01', end_date=None):
        all_dvs = {}
        for station in stations:
            try:
                df = self.get_dvs(station, parameter_code=parameter_code, start_date=start_date, end_date=end_date)
            except Exception as e:
                print(e)
                continue
            all_dvs[station] = df
        self.dvs = all_dvs
        return all_dvs

    def q90(self, stations=None, start_date='1880-01-01', end_date=None):

        if len(self.dvs) == 0 and stations is not None:
            self.get_all_dvs(stations, start_date=start_date, end_date=end_date)
        Q90 = {}
        for site_no, dvs in list(self.dvs.items()):
            DDcd = [c for c in dvs.columns if '00060' in c and 'cd' not in c][0]
            DDvalues = dvs[DDcd].convert_objects(convert_numeric=True)
            Q90[site_no] = DDvalues.quantile(0.1)
        self.dv_q90 = Q90
        return Q90

    def number_of_sites_measured_by_year(self, df):
        """Computes the number of sites measured in each year. The dataframe is grouped by year,
        then by site, and the number of sites for each year is summed.

        Parameters
        ----------
        df:
            Dataframe of NWIS field measurement values indexed by datetime

        returns: nm
            Series of number of measurements, indexed by year
        """
        grouped = df.groupby(df.measurement_dt.dt.year)
        grouped = [(y, g.groupby('site_no').agg('mean')) for y, g in grouped]
        nmeasurements = [(y, len(g)) for y, g in grouped]
        nm = pd.DataFrame(nmeasurements, columns=['year', 'n'])
        nm.index = nm.year
        return nm['n']

    def baseflow_summary(self, q90_window=20, output_proj4=None):

        if self.field_measurements['measurement_dt'].dtype != 'datetime64[ns]':
            self.field_measurements['measurement_dt'] = \
                pd.to_datetime(self.field_measurements.measurement_dt)
        fm = self.field_measurements

        field_sites = self.field_sites.copy()

        # reprojected the output X, Y coordinates
        print('reprojecting output from\n{}\nto\n{}...'.format(self.proj4, output_proj4))
        if output_proj4 is not None:
            field_sites['geometry'] = projectdf(field_sites, self.proj4, output_proj4)

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
            for site_no, data in list(self.dvs.items()):

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
        df['est_error'] = [self.est_error.get(q, self.default_error) for q in df.quality]
        df = df[['site_no', 'datetime', 'Qm', 'quality', 'est_error',
                 'idx_station', 'indexQr', 'indexQ90', 'drn_area', 'station_nm', 'X', 'Y']]
        return df

    def write_shp(self, df, shpname='NWIS_export.shp', **kwargs):
        """Write a shapefile of points from NWIS site file

        Parameters
        ----------
        df: dataframe
            dataframe of site info, must have dec_long_va and dec_lat_va columns with lon/lat in DD

        shpname: string
            Name for output shapefile

        Notes
        -----
        NAD83 is assumed for dec_long_va and dec_lat_va.
        If some entries are in NAD27, a difference of ~5 to >15m will result for WI
        (see http://en.wikipedia.org/wiki/North_American_Datum#/media/File:Datum_Shift_Between_NAD27_and_NAD83.png)
        """
        shpdf = df.copy()
        shpdf['geometry'] = [Point(r.dec_long_va, r.dec_lat_va) for i, r in shpdf.iterrows()]
        GISio.df2shp(shpdf, shpname, epsg=4269)

'''
field measurements url:

url = "http://waterdata.usgs.gov/nwis/measurements? \
nw_longitude_va=-91.497& \
nw_latitude_va=46.748&
se_longitude_va=-90.228&
se_latitude_va=46.156&
coordinate_format=decimal_degrees&
group_key=NONE&
format=sitefile_output&
sitefile_output_format=rdb&
column_name=agency_cd&
column_name=site_no&
column_name=station_nm&
column_name=site_tp_cd&
column_name=lat_va&
column_name=long_va&
column_name=dec_lat_va&
column_name=dec_long_va&
column_name=coord_meth_cd&
column_name=coord_acy_cd&
column_name=coord_datum_cd&
column_name=dec_coord_datum_cd&
column_name=district_cd&
column_name=state_cd&
column_name=county_cd&
column_name=country_cd&
column_name=land_net_ds&
column_name=map_nm&
column_name=map_scale_fc&
column_name=alt_va&
column_name=alt_meth_cd&
column_name=alt_acy_va&
column_name=alt_datum_cd&
column_name=huc_cd&
column_name=basin_cd&
column_name=topo_cd&
column_name=data_types_cd&
column_name=instruments_cd&
column_name=construction_dt&
column_name=inventory_dt&
column_name=drain_area_va&
column_name=contrib_drain_area_va&
column_name=tz_cd&
column_name=local_time_fg&
column_name=reliability_cd&
column_name=gw_file_cd&
column_name=nat_aqfr_cd&
column_name=aqfr_cd&
column_name=aqfr_type_cd&
column_name=well_depth_va&
column_name=hole_depth_va&
column_name=depth_src_cd&
column_name=project_no&
column_name=rt_bol&
column_name=peak_begin_date&
column_name=peak_end_date&
column_name=peak_count_nu&
column_name=qw_begin_date&column_name=qw_end_date&column_name=qw_count_nu&column_name=gw_begin_date&column_name=gw_end_date&column_name=gw_count_nu&column_name=sv_begin_date&column_name=sv_end_date&column_name=sv_count_nu&set_logscale_y=1&channel_html_info=0&date_format=YYYY-MM-DD&channel_rdb_info=0&rdb_compression=file&list_of_search_criteria=lat_long_bounding_box"

Daily values url
http://waterdata.usgs.gov/nwis/dv?referred_module=sw&site_tp_cd=ST&nw_longitude_va=-91&nw_latitude_va=47&se_longitude_va=-90&se_latitude_va=46&coordinate_format=decimal_degrees&group_key=NONE&format=sitefile_output&sitefile_output_format=rdb&column_name=agency_cd&column_name=site_no&column_name=station_nm&range_selection=days&period=365&begin_date=2014-04-14&end_date=2015-04-13&date_format=YYYY-MM-DD&rdb_compression=file&list_of_search_criteria=lat_long_bounding_box
'''
