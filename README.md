# NWIS
Classes for working with the USGS National Water Information System

###Example of retrieving site information for an area defined by a lat/lon bounding box:
```python
import nwis
from attributes import streamflow_attributes

ll_bbox = [-91.497, 46.748, -90.228, 46.156] # nw lon, nw lat, se lon, se lat

nwis = NWIS(ll_bbox)

# Generate a url to get field measurements for the bounding box
url = nwis.make_url('field_measurements', streamflow_attributes)

# Get a dataframe of site information for the bounding box (url is generated internally)
fm_siteinfo = nwis.get_siteinfo('field_measurements', streamflow_attributes)

# Write the site information out to a shapefile
nwis.write_shp(fm_siteinfo, 'shps/NWIS_field_measurements.shp')

# Get site information for daily values
dv_siteinfo = nwis.get_siteinfo('dv', streamflow_attributes)
```
