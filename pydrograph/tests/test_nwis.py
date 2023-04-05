from shapely.geometry import box, Polygon
import pytest
from gisutils import project
from pydrograph import Nwis
import numpy as np

@pytest.fixture(scope='session')
def extent_poly():
    extent_poly = box(390000, 1330000, 500000, 1455000)
    extent_poly_ll = project(extent_poly, "+init=epsg:{}".format(5070), "+init=epsg:4269")
    return extent_poly_ll

@pytest.fixture(scope='session')
def nwis_instance(extent_poly):
    nwis_instance = Nwis(extent=extent_poly)
    return nwis_instance

@pytest.fixture(scope='session')
def field_sites(nwis_instance):
    field_sites = nwis_instance.get_siteinfo('field_measurements')
    return field_sites

@pytest.fixture(scope='session')
def stations(nwis_instance):
    stations = nwis_instance.get_iv_siteinfo(attributes='iv_attributes')
    stations = stations.site_no.unique()
    #stations = stations.to_list()
    return stations

@pytest.mark.parametrize('shapefile', (
    # esri? shapefile that failed with old-style fiona CRS handling
    ('cacheeastcr_chd_perimeter_bufferinside50m.shp'),
    ))
def test_extent_shapefile(test_data_path, shapefile):
    shapefile = test_data_path / shapefile
    nwis_instance = Nwis(extent=shapefile)
    assert isinstance(nwis_instance.extent, Polygon)

def test_compute_geometries(extent_poly, nwis_instance, field_sites):
    geoms = nwis_instance._compute_geometries(field_sites)
    assert all([g.within(extent_poly) for g in geoms])
    
def test_instantaneous_value(nwis_instance):
    #make sure all inputs are available for this function
    #check that it creates a dataframe (use assert and exists)
    df = nwis_instance.get_iv_siteinfo(attributes = 'iv_attributes')
    assert len(df) > 0
    assert 'site_no' in df.columns
    assert df.site_no.dtype == object

def test_tuple_extent_no_data():

    bbox = (-91.45793026894977, 47.2, 
            -90.20509548401013, 47.3)
    nwis = Nwis(extent=bbox)
    assert np.allclose(bbox, nwis.extent.bounds)
    gwdv_sites = nwis.get_siteinfo('gwdv')
    assert gwdv_sites is None
    
def test_get_all_ivs(nwis_instance, stations):
    all_sites = nwis_instance.get_all_ivs(stations)
    site_one = list(all_sites.values())[0]
    assert len(all_sites) > 0
    assert len(site_one) > 2
    #assert all_sites.site_no.dtype == object