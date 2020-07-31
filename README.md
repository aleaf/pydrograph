# pydrograph
Package for getting and processing stream flow and groundwater level measurements from the USGS National Water Information System (NWIS). Data can be queried by site number, lat/lon limits, or with a polygon shapefile. URLs are constructed for the rdb (tab-delimted) format. Data returned by the URLs are then parsed into pandas dataframes. Summary tables of available field measurements and daily values are generated automatically on instantiation of an `Nwis` object. Field measurements or daily values can then be fetched by site numbers. A baseflow module can perform hydrograph separation on daily streamflow values using the modified Base-Flow-Index (BFI) method (Wahl and Wahl, 1988; Institute of Hydrology, 1980). Annual base flows can also be estimated from miscellaneous measurements using the Wisconsin Statewide Equation method of Gebert and others (2007).

### Version 0

[![Build Status](https://travis-ci.com/aleaf/pydrograph.svg?branch=master)](https://travis-ci.com/aleaf/pydrograph.svg)
[![Build status](https://ci.appveyor.com/api/projects/status/u1mdn7enfel7u39n?svg=true)](https://ci.appveyor.com/project/aleaf/pydrograph)
[![codecov](https://codecov.io/gh/aleaf/pydrograph/branch/master/graph/badge.svg)](https://codecov.io/gh/aleaf/pydrograph)

Getting Started
-----------------------------------------------
### Example notebooks
<http://nbviewer.ipython.org/github/aleaf/pydrograph/blob/master/examples/Notebooks/NWIS_example.ipynb>
<https://github.com/aleaf/pydrograph/blob/master/examples/Notebooks/NWIS_gw_example.ipynb>  
[Demonstration/benchmark of the BFI method](https://github.com/aleaf/pydrograph/blob/master/examples/Notebooks/IHmethod_demo.ipynb)


### Get site information for an area defined by a lat/lon bounding box:
```python
import gisutils
import pydrograph
from pydrograph.attributes import streamflow_attributes

ll_bbox = [-91.497, 46.748, -90.228, 46.156] # nw lon, nw lat, se lon, se lat

nwis = pydrograph.Nwis(ll_bbox)

# Generate a url to get field measurements for the bounding box
url = nwis.make_site_url('field_measurements', streamflow_attributes)

# Get a dataframe of site information for the bounding box (url is generated internally)
field_sites = nwis.get_siteinfo('field_measurements')

# Write the site information out to a shapefile
gisutils.df2shp(field_sites, 'NWIS_field_measurements.shp')

# Get inventory of daily values sites
dv_sites = nwis.get_siteinfo('daily_values')

# Get daily values for a single site
df = nwis.get_dvs(4015475)
```
### Bugs

If you think you have discovered a bug in pydrograph in which you feel that the program does not work as intended, then we ask you to submit a [Github issue](https://github.com/aleaf/pydrograph/labels/bug).


Installation
-----------------------------------------------

**Python versions:**

pydrograph requires **Python** 3.7 (or higher)

**Dependencies:**  
numpy   
pandas  
fiona   
shapely  
pyproj  
gisutils  (available from pip or atleaf conda channel)

### Install python and dependency packages
Download and install the [Anaconda python distribution](https://www.anaconda.com/distribution/).
Open an Anaconda Command Prompt on Windows or a terminal window on OSX.
From the root folder for the package (that contains `requirements.yml`), install the above packages from `requirements.yml`:

```
conda env create -f requirements.yml
```
activate the environment:

```
conda activate pydrograph
```

### Install to site_packages folder
```
python setup.py install
```
### Install in current location (to current python path)
(i.e., for development)  

```  
pip install -e .
```



Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.


References
----------
Gebert, W.A., Radloff, M.J., Considine, E.J., and Kennedy, J.L., 2007,
Use of streamflow data to estimate base flow/ground-water recharge for Wisconsin:
Journal of the American Water Resources Association,
v. 43, no. 1, p. 220-236, http://dx.doi.org/10.1111/j.1752-1688.2007.00018.x

Gebert, W.A., Walker, J.F., and Kennedy, J.L., 2011,
Estimating 1970-99 average annual groundwater recharge in Wisconsin using streamflow data:
U.S. Geological Survey Open-File Report 2009-1210, 14 p., plus appendixes,
available at http://pubs.usgs.gov/ofr/2009/1210/.
    
Institute of Hydrology, 1980b, Low flow studies report no. 3--Research report: Wallingford, Oxon, United Kingdom, Institute of Hydrology Report no. 3, p. 12-19       
    
Wahl, K.L and Wahl, T.L., 1988. Effects of regional ground-water level declines
on streamflow in the Oklahoma Panhandle. In Proceedings of the Symposium on 
Water-Use Data for Water Resources Management, American Water Resources Association. 