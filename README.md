# NWIS
Classes for working with the USGS National Water Information System. This repo is a work in progress. Currently the nwis.py module is limited to generating site information for a lat/lon bounding box (see below). Methods for pulling actual data from the results of the bounding box query will be added. 

The scripts folder contains some rough code examples invovling fetching and/or processing information from NWIS. The Notebooks folder has a Notebook from the USGS python class that illustrates fetching daily values using a url and then manipulating the data in pandas. You can view it here: <http://nbviewer.ipython.org/github/aleaf/NWIS/blob/master/Notebooks/NWIS_example.ipynb>



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
