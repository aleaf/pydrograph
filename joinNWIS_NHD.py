__author__ = 'aleaf'
'''
Match up NWIS stations with NHD flowlines (COMIDS)
'''

import sys
GISio_path = 'D:\\ATLData\\Documents\\GitHub\\GIS_utils\\'
if GISio_path not in sys.path:
    sys.path.append(GISio_path)
import GISio
import pandas as pd
from shapely.geometry import Point
import arcpy

# inputs for starting from raw NWIS files (requires arcpy)
model_domain_polygon = 'D:/ATLData/CNNF_Great_Divide/GIS/shps/GD_nearfield.shp' # must be in geographic coordinates!
NWIS_site_info_file = 'D:/ATLData/GFL files/Great_Divide/flux_targets/NWIS_sites.txt'
NWIS_measurements_file = 'D:/ATLData/GFL files/Great_Divide/flux_targets/NWIS_values.txt'
flowlines_clipped = 'D:/ATLData/GFL files/Great_Divide/flux_targets/NHDflowlines_clip.shp'

arcpy.env.overwriteOutput = True

# spatial join tolerance (in decimal degrees)
tol = .001 # 0.001 worked for the CNNF model (~ -46.5 latitude)

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


# read in NWIS site information and study area boundary
header_text = open(NWIS_site_info_file).readlines()
columns, header_rows = NWIS_header(header_text)

df = pd.read_csv(NWIS_site_info_file, sep='\t', names=columns, skiprows=header_rows)
bounds = GISio.shp2df(model_domain_polygon, geometry=True).geometry[0]

# make geomtries for each station, and drop stations not in the study area
df['geometry'] = df.apply(lambda x: Point(x['dec_long_va'], x['dec_lat_va']), axis=1)
GISio.df2shp(df, 'D:/ATLData/GFL files/Great_Divide/flux_targets/NWIS_sites_all.shp', prj='epsg:4269')

within = [p.within(bounds) for p in df.geometry]
df = df[within]
GISio.df2shp(df, NWIS_site_info_file[:-4]+'.shp', prj='epsg:4269')

# now do spatial join of NWIS locations to NHD comids
arcpy.SpatialJoin_analysis(NWIS_site_info_file[:-4]+'.shp', flowlines_clipped, NWIS_site_info_file[:-4]+'_joined.shp',
                           "JOIN_ONE_TO_ONE", "KEEP_ALL", '', "WITHIN_A_DISTANCE", .001)

# now read back in and make a csv file for input into flux_targets.py
df = GISio.shp2df(NWIS_site_info_file[:-4]+'_joined.shp')

site_info = df[['site_no', 'COMID']]

# read in NWIS measurements
header_text = open(NWIS_measurements_file).readlines()
columns, header_rows = NWIS_header(header_text)

df = pd.read_csv(NWIS_measurements_file, sep='\t', names=columns, skiprows=header_rows, index_col=3, parse_dates=True)
df['datetime'] = df.index

# join site info to measurements and save
# (inner only retains sites that are in both frames; sites outside of the study area dropped above will not be joined)
df = df.merge(site_info, how='inner', left_on='site_no', right_on='site_no', suffixes=['','1'])
df.index = df['datetime']
df.to_csv(NWIS_measurements_file[:-4]+'_NHD.csv', index_label='datetime')




