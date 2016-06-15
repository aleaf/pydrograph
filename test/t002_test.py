import sys
#sys.path.append('/Users/aleaf/Documents/GitHub/NWIS')
import os
import pandas as pd
import numpy as np
from nwis import NWIS
from GISio import shp2df


cpth = 'temp'

def test_write_shp():

    nwis = NWIS(extent='../examples/data/muk_nf.shp')

    nwis.write_shp(nwis.field_sites, os.path.join(cpth, 'sw_field_sites.shp'))
    nwis.write_shp(nwis.gwfield_sites, os.path.join(cpth, 'gw_field_sites.shp'))

if __name__ == '__main__':
    if not os.path.isdir(cpth):
        os.makedirs(cpth)
    test_write_shp()