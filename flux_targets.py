__author__ = 'aleaf'

import sys
import datetime as dt
GISio_path = '../GIS_utils'
if GISio_path not in sys.path:
    sys.path.append(GISio_path)
import GISio
import urllib2
import numpy as np
import pandas as pd
import re

station_IDs_col = 'index_station'
start_dates_col = 'index_start'

measurements_file = '../../FishHatch/flux_target_info.csv'
CumulativeArea_file = None

# if True, get drainage areas for comids from NHDPlus, otherwise, supply in 'Area' column of measurements file
drainage_areas_from_NHD = False
area_units_mult = 0.386102 # multiplier to convert area units to mi2

output_file = 'adjusted_baseflows.csv'

# find start of data and column names in NWIS file
def NWIS_header(text):
    knt = 0
    for line in text:
        if line.strip().split('\t')[0] == 'agency_cd':
            columns = line.strip().split('\t')
            knt += 2
            break
        else:
            knt += 1
    return columns, knt

# get NWIS data for index sites
def get_nwis(station_ID, parameter_code, start_date='1880-01-01'):
    '''
    reads discharge values for an NWIS gage site into a pandas dataframe
    To obtain the entire period-of-record use a start date of 1880-01-01...

    stationID: (string) USGS station ID

    parameter_code: (string) e.g. 00060 for discharge

    start_date: (string) 'YYYY-DD-MM'
    '''

    print "getting data for site {0}...".format(station_ID)
    #url = 'http://nwis.waterdata.usgs.gov/wi/nwis/uv/?cb_00060=on&format=rdb&site_no={0}&begin_date={1}'\
    #    .format(station_ID, start_date)
    #url = 'http://waterservices.usgs.gov/nwis/dv?format=waterml,1.1&sites={0}&startDT={1}'\
    #    .format(station_ID, start_date)
    if len(str(station_ID)) == 7:
        url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&parameterCd={0}&sites=0{1}&startDT={2}'\
            .format(parameter_code, station_ID, start_date)
    else:
        url = 'http://waterservices.usgs.gov/nwis/dv/?format=rdb&parameterCd={0}&sites={1}&startDT={2}'\
            .format(parameter_code, station_ID, start_date)

    try:
        text = urllib2.urlopen(url).readlines()
    except:
        text = "{} ({})".format(sys.exc_info()[1], station_ID)
        print text

    return text


def NWIS2df(text, parameter='00060'):

    columns, knt = NWIS_header(text)

    try:
        parameter_cd = [c for c in columns if parameter in c and '_cd' not in c][0]
    except:
         print 'column for parameter {} not found!'.format(parameter)

    print text[knt]

    # convert tabular data to dataframe
    df = pd.DataFrame([s.strip().split('\t') for s in text[knt:]])
    df.columns = columns
    df.index = pd.to_datetime(df['datetime'])
    df = pd.DataFrame(df[parameter_cd])

    # remove 'ICE' or other text flags in discharge values
    non_decimal = re.compile(r'[^\d.]+') # expression for getting non-decimal characters (otherwise some measurements will be multiplied by 10!)

    df[parameter_cd] = df[parameter_cd].astype('str').apply(lambda x: non_decimal.sub("", x))

    # convert discharge values to floats
    df[parameter_cd] = df[parameter_cd].convert_objects(convert_numeric=True)

    return df


def statewide_eqn(Qm, A, Qr, Q90):
    Bf = (Qm / A) * (Q90 / Qr)
    Qb = 0.907 * A**1.02 * Bf**0.52
    return Qb, Bf


def get_drainage_areas(comids, CumulativeArea_file, units='mi2', areatype='TotDASqKM'):
    CA = GISio.shp2df(CumulativeArea_file)
    CA.index = CA.ComID
    areas_dict = {}
    for comid in comids:
        if units == 'mi2':
            areas_dict[comid] = CA.ix[comid, areatype] * 0.386102
        elif units == 'km2':
            areas_dict[comid] = CA.ix[comid, areatype]
    return areas_dict


def match_index_station(station, index_stations):
    for ID in index_stations:
        if str(station) in ID:
            station_ID = ID
            break
        else:
            station_ID = '9999'
    return station_ID


if __name__ == "__main__":

    # read in measurements
    M = pd.read_csv(measurements_file, parse_dates=True)



    if drainage_areas_from_NHD:
        # get drainage areas
        comids = list(np.unique(M['comid']))
        areas_dict = get_drainage_areas(comids, CumulativeArea_file)
        M['A'] = [areas_dict[c] for c in M.comid]
    else:
        area_col = [c for c in M.columns if c.lower() == 'area'][0]
        M['A'] = M[area_col] * area_units_mult

    # get data and calculate Q90 values for stations
    unique_IDs = np.unique(M[station_IDs_col])
    start_dates = [M.iloc[np.where(M[station_IDs_col] == uid)[0][0]][start_dates_col] for uid in unique_IDs]

    # format the dates correctly
    start_dates_f = []
    for s in start_dates:
        t = pd.to_datetime(s)
        start_dates_f.append('{:d}-{:02d}-{:02d}'.format(t.year, t.month, t.day))

    data = {}
    Q90 = {}
    # get Q90 values
    for i, u in enumerate(unique_IDs):
        rawtext = get_nwis(str(u), '00060', start_date=start_dates_f[i])
        data[u] = NWIS2df(rawtext)
        Q90[u] = data[u]['02_00060_00003'].quantile(q=0.1)

    # get recorded Q values at index stations, add index station info to measurements dataframe
    M['Qr'] = np.zeros(len(M))
    M['Q90'] = np.zeros(len(M))
    M['Bf'] = np.zeros(len(M))
    for i in M.index:
        datetime = pd.to_datetime(M.ix[i, 'datetime'])
        index_station = M.ix[i, station_IDs_col]
        dt_ind = data[index_station].index.searchsorted(datetime)
        M.ix[i, 'Qr'] = data[index_station].ix[data[index_station].index[dt_ind], '02_00060_00003']
        M.ix[i, 'Q90'] = Q90[index_station]

    # calculate baseflow, add to measurements dataframe
    M['Qb'], M['Bf'] = statewide_eqn(M['Qm'], M['A'], M['Qr'], M['Q90'])

    # include CFD
    M['Qb_cfd'] = M['Qb'] * 3600 * 24

    # write output
    M.to_csv(output_file, index=False)

