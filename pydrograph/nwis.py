import datetime as dt
from pathlib import Path
import time
from urllib.request import urlopen
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, shape, box
import gisutils
from .attributes import streamflow_attributes, gw_attributes, iv_attributes


coord_datums_epsg = {'NAD83': 4269,
                     'NAD27': 4267}


def WI_statewide_eqn(Qm, A, Qr, Q90):
    Bf = (Qm / A) * (Q90 / Qr)
    Qb = 0.907 * A**1.02 * Bf**0.52
    return Qb.copy(), Bf.copy()


class Nwis:
    """
    NWIS error codes:
    E     Excellent    The data is within 2% (percent) of the actual flow
    G     Good         The data is within 5% (percent) of the actual flow
    F     Fair         The data is within 8% (percent) of the actual flow
    P     Poor         The data are >8% (percent) of the actual flow
    """

    est_error = {'excellent': 0.02,
                 'good': 0.05,
                 'fair': 0.08
                 }
    default_error = 0.50

    urlbase = 'http://nwis.waterdata.usgs.gov/usa/nwis/'
    urlbase_iv = 'http://waterservices.usgs.gov/nwis/site/?format=rdb&bBox='
    dtypes_dict = {'dv': 'dv?referred_module=sw&site_tp_cd=ST&',
                   'daily_values': 'dv?referred_module=sw&site_tp_cd=ST&',
                   'field_measurements': 'measurements?',
                   'gwlevels': 'gwlevels?',
                   'gwdv': 'dv?referred_module=gw&site_tp_cd=GW&',
                   'gw_daily_values': 'dv?referred_module=gw&site_tp_cd=GW&',
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

    log_cols = ['site_no', 'url', 'retrieved', 'data_found']

    def __init__(self, bounds_latlon=None, extent=None, datum='NAD83',
                 log=False):
        """Class for retrieving data from NWIS.

        Parameters
        ----------
        bounds_latlon: sequence of floats
            Sequence containing bounding latitudes and longitudes of query area
            in decimal degrees, in the following order:
            [northwest longitude, northwest latitude, southeast longitude, southeast latitude]
        extent: filepath (str) or shapely polygon
            Polygon of area to query. A polygon from a shapefile
            will be automatically reprojected to lat/lon (epsg:4269);
            shapely polygons are assumed to be in geographic coordinates.

        see the code for a list of default parameters
        """

        self.bounds_latlon = bounds_latlon
        self.datum = datum
        self.crs = gisutils.get_authority_crs(coord_datums_epsg[self.datum])
        self.log = pd.DataFrame(columns=self.log_cols)
        self._extent = None
        self.extent = extent
        self._bounds_latlon = None
        self.bounds_latlon = bounds_latlon
        self.write_log_file = log

    @property
    def extent(self):
        """Polygon of query area in lat/lon (epsg:4269)"""
        return self._extent

    @extent.setter
    def extent(self, extent=None):
        if extent is not None:
            if isinstance(extent, str) or isinstance(extent, Path):
                # _read_extent_shapefile should automatically reproject to 4269
                self._extent = self._read_extent_shapefile(extent)
            elif isinstance(extent, Polygon):
                self._extent = extent
            else:
                try:
                    self._extent = box(*extent)
                except:
                    print('Warning: extent argument of unknown datatype!')
                    self._extent = None
        else:
            self._extent = None

    @property
    def bounds_latlon(self):
        """Bounding box of query area in lat/lon (epsg:4269)"""
        if self._bounds_latlon is None and self.extent is not None:
            self._bounds_latlon = self.extent.bounds
        return self._bounds_latlon

    @bounds_latlon.setter
    def bounds_latlon(self, bounds_latlon=None):
        self._bounds_latlon = bounds_latlon

    def _compute_geometries(self, df):

        datum = np.array([int(coord_datums_epsg[d]) for d in df.dec_coord_datum_cd])
        datums = np.unique(datum)
        x1, y1 = df.dec_long_va.values, df.dec_lat_va.values
        x2 = np.ones(len(df), dtype=float) * np.nan
        y2 = np.ones(len(df), dtype=float) * np.nan
        for dtm in datums:
            pr1 = gisutils.get_authority_crs(int(dtm))
            loc = datum == dtm
            x2[loc], y2[loc] = gisutils.project((x1[loc], y1[loc]), pr1, self.crs)
        geoms = [Point(x, y) for x, y in zip(x2, y2)]
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
        with fiona.open(shpfile) as src:
            g = shape(next(iter(src))['geometry'])
        shpfile_crs = gisutils.get_shapefile_crs(shpfile)
        if self.crs != shpfile_crs:
            print(f'reprojecting extent from\n{shpfile_crs}\nto\n{self.crs}')
            return gisutils.project(g, shpfile_crs, self.crs)
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
        self.bbox_url = 'nw_longitude_va={:.3f}&'.format(self.bounds_latlon[0]) +\
                        'nw_latitude_va={:.3f}&'.format(self.bounds_latlon[3]) +\
                        'se_longitude_va={:.3f}&'.format(self.bounds_latlon[2]) +\
                        'se_latitude_va={:.3f}&'.format(self.bounds_latlon[1])

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

        if data_type in {'dv', 'daily_values'}:
            url += self.dv_info

        url += self.stuff_at_end
        #print '{}'.format(url)
        return url
        
    def make_iv_site_url(self):
        """ This function makes a url to pull instantaneous values from the nwis webservice. lat/lon box must be defined earlier using a llbox initialization of the nwis class.
                
        Parameters
        ----------
        data_type: str

        Returns
        -------
        url string
        """
        
        self.bbox_url = '{:.3f},'.format(self.bounds_latlon[0]) +\
                '{:.3f},'.format(self.bounds_latlon[1]) +\
                '{:.3f},'.format(self.bounds_latlon[2]) +\
                '{:.3f}'.format(self.bounds_latlon[3])

        self.stuff_at_end = '&outputDataTypeCd=iv,id&siteStatus=all&hasDataTypeCd=iv'

        url = self.urlbase_iv
        url += self.bbox_url

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

        start_date: (string) 'YYYY-MM-DD'
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
        station_IDs = ','.join([Nwis.correct_stationID(s) for s in station_IDs])

        url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb'

        url += '&sites={}'.format(station_IDs)
        url += '&startDT={}'.format(start_date)
        if end_date is not None:
            url += '&endDT={}'.format(end_date)
        url += '&parameterCd={}'.format(parameter_code)
        print('{}'.format(url))
        return url
        
    def make_iv_url(self, station_IDs, parameter_code='00060', start_date='2000-01-01', end_date= '2000-12-31'):
        """Creates url to retrieve instantaneous values for a site given a list of station IDs.


        Parameters
        ----------
        stationIDs: int, str or list of ints or strings
            USGS station IDs

        parameter_code: (string)
            e.g. 00060 for discharge.
            See http://help.waterdata.usgs.gov/codes-and-parameters/parameters.

        start_date: (string) 'YYYY-MM-DD'
            default start and end dates are written to pull the year 2000, not the entire history of data,
            as these are instantaneous values, so there is a datapoint every ~15 minutes. pull smaller chunks of data
            to ensure that the url and code will run.

        Notes
        -----
        A leading zero is added to the site number if the first digit is greater than 1
        (this can happend for basins 01 - 09 if the site number gets converted to an int).
        Note that this may cause site numbers for basin 01 (North Atlantic slope) to get confused with
        basins 10-16 (west coast and hawaii).
        See <http://help.waterdata.usgs.gov/faq/sites/do-station-numbers-have-any-particular-meaning>

        """
        #need to code in a more systematic way of choosing dates, but for now, pulling a whole year is more than enough
        if not isinstance(station_IDs, list):
            station_IDs = [str(station_IDs)]

        def add_leading_zero(station_ID):
            if 1 < int(str(station_ID)[0]) < 10:
                station_ID = '0{}'.format(station_IDs)
            return station_ID

        #station_IDs = ','.join(['0{}'.format(int(str(s))) for s in station_IDs])
        station_IDs = ','.join([Nwis.correct_stationID(s) for s in station_IDs])

        url = 'http://waterservices.usgs.gov/nwis/iv/?format=rdb'

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
        station_ID = Nwis.correct_stationID(station_ID)

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

    def get_datetime_retrieved(self, sitefile_text):
        for line in sitefile_text:
            if 'retrieved' in str(line):
                return str(line).strip().split('retrieved:')[-1][:30].strip()
            elif '#' not in str(line):
                return None

    def get_siteinfo(self, data_type, attributes=None):
        """Retrieves site information for the bounding box supplied to the NWIS class instance

        Parameters
        ----------
        data_type: str
            'daily_values' for Daily Values
            'field_measurements' for Field Measurements
            'gwlevels' for groundwater field measurements
            'gwdv' for groundwater daily values

        attributes: list of strings
            List of NWIS attributes to include (e.g. 'site_no', 'station_nm', etc.)
            Default sets of attributes for streamflow and groundwater levels data can
            be imported from the attributes.py file (work in progress)

        Returns
        -------
        the contents of an NWIS site information file in a dataframe format
        """
        print('getting site inventory for {}...'.format(data_type))
        t0 = time.time()
        if attributes is None:
            if data_type in {'dv', 'daily_values', 'field_measurements'}:
                attributes = streamflow_attributes
            elif data_type in {'gwdv', 'gw_daily_values', 'gwlevels'}:
                attributes = gw_attributes
        url = self.make_site_url(data_type, attributes)
        print('url: {}'.format(url))
        sitefile_text = urlopen(url).readlines()
        if 'DOCTYPE html' in sitefile_text[0].decode():
            print(f"No sites found for {self.bounds_latlon}!")
            return None
            
        skiprows = self.get_header_length(sitefile_text, attributes[0])

        print('reading data with pandas...')
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=attributes,
                         dtype={'site_no': object})
        print("finished in {:.2f}s\n".format(time.time() - t0))
        df['geometry'] = self._compute_geometries(df)
        df.index = df.site_no
        n_sites = len(df)
        if self.extent is not None:
            print('culling {} sites to those within extent...'.format(n_sites))
            within = np.array([g.within(self.extent) for g in df.geometry])
            df = df[within].copy()
        print("finished inventory in {:.2f}s\n".format(time.time() - t0))
        return df

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

        start_date: (string) 'YYYY-MM-DD'
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
        loginfo = [str(station_ID), url, self.get_datetime_retrieved(sitefile_text)]
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=cols,
                         dtype={'site_no': object})
        if len(df) > 0:
            df.index = pd.to_datetime(df.datetime)
            loginfo.append(True)
        else:
            loginfo.append(False)
        self.log = pd.concat([self.log,
                              pd.DataFrame([loginfo], columns=self.log_cols)])
        return df
        
    def get_iv_siteinfo(self, attributes = 'iv_attributes'):
        """Retrieves site information for the bounding box supplied to the NWIS class instance.

        Parameters
        ----------

        attributes: preset, 'iv_attributes'

        Returns
        -------
        the contents of an NWIS site information file in a dataframe format
        """
        #print('getting site inventory for {}...'.format(data_type))
        t0 = time.time()
        attributes = iv_attributes
        url = self.make_iv_site_url()
        print('url: {}'.format(url))
        sitefile_text = urlopen(url).readlines()
        skiprows = self.get_header_length(sitefile_text, attributes[0])

        print('reading data with pandas...')
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=attributes, dtype={'site_no': object})
        print("finished in {:.2f}s\n".format(time.time() - t0))
        df['geometry'] = self._compute_geometries(df)
        df.index = df.site_no
        n_sites = len(df)
        if self.extent is not None:
            print('culling {} sites to those within extent...'.format(n_sites))
            within = np.array([g.within(self.extent) for g in df.geometry])
            df = df[within].copy()
        print("finished inventory in {:.2f}s\n".format(time.time() - t0))
        return df          

    
    def get_ivs(self, station_ID, parameter_code='00060', start_date='2000-01-01', end_date='2000-12-31',
        sample_period = 'D', agg_method = 'mean'):
        """Retrieves daily values for a site.

        Parameters
        ----------
        stationID: (string)
            USGS station ID

        parameter_code: string, default is 00060 for discharge.
            See http://help.waterdata.usgs.gov/codes-and-parameters/parameters.

        start_date: (string) 'YYYY-MM-DD'
            To obtain the entire period-of-record use a start date of 2000-01-01 (default)...
        end_date: (string) 'YYYY-MM-DD'
            preset to take the year 2000, don't set the range too long (longer than a year or so) or the code will be very slow and may not finish running
        sample_period: (string), default 'D' for daily
            change this string to how you would like the instantaneous values aggregated. this can be daily, weekly, monthly, etc.
            None will give just the raw data which is generally in 15 minute increments
        agg_method: (string), default 'mean'
            change this to chose how the aggregated instantaneous values are calculated. options include 'mean', 'median', 'max', etc.

        Returns
        -------
        df: a datetime-index dataframe of daily discharge, with datagaps filled with NaNs, aggregated daily or however specified with 'sample_period'
        """
        if parameter_code in list(self.parameter_codes.keys()):
            parameter_code = self.parameter_codes[parameter_code]

        url = self.make_iv_url(station_ID, parameter_code=parameter_code,
                               start_date=start_date, end_date=end_date)
        sitefile_text = urlopen(url).readlines()
        skiprows = self.get_header_length(sitefile_text, 'agency_cd')
        cols = sitefile_text[skiprows - 2].decode('utf-8').strip().split('\t')
        loginfo = [str(station_ID), url, self.get_datetime_retrieved(sitefile_text)]
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=cols, dtype={'site_no': object})

        if len(df) > 2:

            if len(df) > 0:
                df.index = pd.to_datetime(df.datetime)
                loginfo.append(True)
            else:
                loginfo.append(False)
            self.log = pd.concat([self.log,
                                  pd.DataFrame([loginfo], columns=self.log_cols)])
        
            if sample_period is not None:
                df = df.resample(sample_period).agg(agg_method)
                df = df.rename(columns = {df.columns[0]: 'discharge (cfs)'})
            else:
                df = df.rename(columns = {df.columns[4]: 'discharge (cfs)', df.columns[5]: 'code'})

        else:
            print('No data at this site during this timeframe.')

        if len(df) > 2:
            return df
        else:
            return None

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
        loginfo = [str(station_ID), url, self.get_datetime_retrieved(sitefile_text)]
        df = pd.read_csv(url, sep='\t', skiprows=skiprows, header=None, names=cols,
                         dtype={'site_no': object})
        if len(df) > 0:
            df.index = pd.to_datetime(df[self._get_date_col(df)])
            loginfo.append(True)
        else:
            loginfo.append(False)
        self.log = pd.concat([self.log,
                              pd.DataFrame([loginfo], columns=self.log_cols)])
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
            if len(df) == 0:
                print('no data returned.')
                continue
            df.index = pd.MultiIndex.from_product([[df.site_no.values[0]], df.index.values],
                                              names=['site_no', 'datetime'])
            df['measurement_dt'] = pd.to_datetime(df[self._get_date_col(df)])
            all_measurements = pd.concat([all_measurements, df])
        if self.write_log_file:
            out_logfile = 'retrieved_{}_log_{}.csv'.format(txt, time.strftime('%Y%m%d%H%M%S'))
            self.log.to_csv(out_logfile, index=False)
            print('Log of query saved to {}'.format(out_logfile))
            self.log = pd.DataFrame(columns=self.log_cols) # reset the log
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
        if self.write_log_file:
            out_logfile = 'retrieved_{}_dvs_log_{}.csv'.format(parameter_code,
                                                           time.strftime('%Y%m%d%H%M%S'))
            self.log.to_csv(out_logfile, index=False)
            print('Log of query saved to {}'.format(out_logfile))
        self.log = pd.DataFrame(columns=self.log_cols)  # reset the log
        return all_dvs

    def get_all_ivs(self, stations, parameter_code='00060', start_date='2000-01-01', end_date='2000-12-31'):
        ''' This function gets all instantaneous values for a list of station IDs, and places them in a dictionary of dataframes
        '''
        
        all_ivs = {}
        for station in stations:
            try:
                df = self.get_ivs(station, parameter_code=parameter_code, start_date=start_date, end_date=end_date)
            except Exception as e:
                print(e)
                continue
            all_ivs[station] = df
        if self.write_log_file:
            out_logfile = 'retrieved_{}_ivs_log_{}.csv'.format(parameter_code,
                                                           time.strftime('%Y%m%d%H%M%S'))
            self.log.to_csv(out_logfile, index=False)
            print('Log of query saved to {}'.format(out_logfile))
        self.log = pd.DataFrame(columns=self.log_cols)  # reset the log
        return all_ivs

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
        gisutils.df2shp(shpdf, shpname, epsg=4269)

    @staticmethod
    def correct_stationID(stationID):
        try:
            if 1 < int(str(stationID)[0]) < 10 and len(str(stationID)) < 15:
                return '0{}'.format(stationID)
        except:
            j=2
        return str(stationID)

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
