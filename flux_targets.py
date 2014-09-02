__author__ = 'aleaf'

import sys
GISio_path = 'D:\\ATLData\\Documents\\GitHub\\GIS_utils\\'
if GISio_path not in sys.path:
    sys.path.append(GISio_path)
import GISio
import urllib2
import numpy as np
import pandas as pd
import re

station_IDs = ['4027000', '4026561', '5362000', '5359500', '5332500']
start_dates = ['1970-01-01', '2011-05-20']


#start_date = '1970-01-01'
#station_ID = '04027000'
discharge_code = '02_00060_00003'

measurements_file = 'D:\\ATLData\\BadRiver\\Testpoints\\flux_targets\\20140924_seepage_run_04026561.csv'
CumulativeArea_file = 'D:\\ATLData\\BadRiver\\BCs\\NHDPlusGL\\NHDPlus04\\NHDPlusAttributes\\CumulativeArea.dbf'

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
    " To obtain the entire period-of-record use a start date of 1880-01-01..."
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

def NWIS2df(text):

    columns, knt = NWIS_header(text)

    if not start_date:
        start_date = header[knt].split()[2]

    print "downloading discharge data for {0} from {1} to present...".format(station_ID, start_date)

    try:
        parameter_cd = [c for c in columns if parameter in c and '_cd' not in c][0]
    except:
         print 'column for parameter {} not found!'.format(parameter)

    df = pd.read_csv(url, sep='\t', names=columns, skiprows=knt, index_col=2, parse_dates=True)
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

    # get drainage areas
    comids = list(np.unique(M['comid']))
    areas_dict = get_drainage_areas(comids, CumulativeArea_file)

    # get data and calculate Q90 values for stations
    data = {}
    Q90 = {}
    # get Q90 values
    for i in range(len(station_IDs)):
        data[station_IDs[i]] = nwis2df(station_IDs[i], start_dates[i])
        Q90[station_IDs[i]] = data[station_IDs[i]]['02_00060_00003'].quantile(q=0.1)

    # get recorded Q values at index stations, add index station info to measurements dataframe
    M['Qr'] = np.zeros(len(M))
    M['Q90'] = np.zeros(len(M))
    M['A'] = np.zeros(len(M))
    M['Bf'] = np.zeros(len(M))
    for i in M.index:
        datetime = pd.to_datetime(M.ix[i, 'datetime'])
        index_station = match_index_station(M.ix[i, 'index_station'], station_IDs)
        dt_ind = data[index_station].index.searchsorted(datetime)
        M.ix[i, 'Qr'] = data[index_station].ix[data[index_station].index[dt_ind], '02_00060_00003']
        M.ix[i, 'index_station'] = index_station
        M.ix[i, 'Q90'] = Q90[index_station]
        M.ix[i, 'A'] = areas_dict[M.ix[i, 'comid']]

    # calculate baseflow, add to measurements dataframe
    M['Qb'], M['Bf'] = statewide_eqn(M['Qm'], M['A'], M['Qr'], M['Q90'])

    # include CFD
    M['Qb_cfd'] = M['Qb'] * 3600 * 24

    # write output
    M.to_csv(output_file)

