import os
import pydrograph
from ..attributes import streamflow_attributes


def test_get_info_from_bounding_box(tmpdir):
    ll_bbox = [-91.497, 46.748, -90.228, 46.156] # nw lon, nw lat, se lon, se lat

    nwis = pydrograph.Nwis(ll_bbox)

    # Generate a url to get field measurements for the bounding box
    url = nwis.make_site_url('field_measurements', streamflow_attributes)

    # Get a dataframe of site information for the bounding box (url is generated internally)
    fm_siteinfo = nwis.get_siteinfo('field_measurements', streamflow_attributes)

    # Write the site information out to a shapefile
    nwis.write_shp(fm_siteinfo, '{}/NWIS_field_measurements.shp'.format(tmpdir))

    # Get site information for daily values
    dv_siteinfo = nwis.get_siteinfo('dv', streamflow_attributes)


def test_readme(tmpdir):

    # copy the python snippet in readme to a .py file in tests output folder
    dest_py_file = os.path.join(tmpdir, 'test_readme.py')
    with open('README.md') as src:
        with open(dest_py_file, 'w') as dest:
            for line in src:
                if "```python" in line:
                    for line in src:
                        if "```" in line:
                            break
                        dest.write(line)

    # execute the .py file
    wd = os.getcwd()
    os.chdir(tmpdir)
    ival = os.system('python test_readme.py')
    assert ival == 0, 'could not run python code in README.md'
    os.chdir(wd)