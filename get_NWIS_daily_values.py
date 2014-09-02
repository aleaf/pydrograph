__author__ = 'aleaf'
'''
retrieve daily values for a list of NWIS sites
'''
import os
import pandas as pd
from shapely.geometry import Point
import flux_targets as ft
import sys
sys.path.append('D:\\ATLData\\Documents\\GitHub\\GIS_utils\\')
import GISio

NWIS_daily_values_sites_file = 'D:/ATLData/GFL files/Great_Divide/flux_targets/NWIS_daily_values_sites.txt'
output_folder = 'D:/ATLData/GFL files/Great_Divide/flux_targets/daily'
model_domain_polygon = 'D:/ATLData/CNNF_Great_Divide/GIS/shps/GD_nearfield.shp' # must be in geographic coordinates!
bounds = GISio.shp2df(model_domain_polygon, geometry=True).geometry[0]

# read in info from daily values sites file
header_text = open(NWIS_daily_values_sites_file).readlines()
columns, header_rows = ft.NWIS_header(header_text)

df = pd.read_csv(NWIS_daily_values_sites_file, sep='\t', names=columns, skiprows=header_rows)

for n in df.site_no:

    # first test if the site is in the model domain
    site_location = Point(df[df.site_no == n][['dec_long_va', 'dec_lat_va']].get_values()[0])

    if site_location.within(bounds):
        text = ft.get_nwis(n, '00060')
        ofp = open(os.path.join(output_folder, '{}.txt'.format(n)), 'w')
        [ofp.write(line) for line in text]
        ofp.close()
