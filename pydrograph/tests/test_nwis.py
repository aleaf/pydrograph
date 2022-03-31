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
    assert df.site_no.dtype == np.object

def test_get_ivs_gw(nwis_instance):
    # check groundwater site (depth to water) 
    gw_site = '353831090502401'
    df = nwis_instance.get_ivs(gw_site, 
                               parameter_code='72019', 
                               start_date='2020-08-01',
                               end_date='2020-08-15')

    assert df.columns[0] == 'Depth to water level, feet below land surface' # check column header/units
    assert df.shape[0] == 15 # check daily resampling 

    # check without resampling
    df = nwis_instance.get_ivs(gw_site, 
                               parameter_code='72019', 
                               start_date='2020-08-01', 
                               end_date='2020-08-15', 
                               sample_period=None)
    assert df.shape[0] == 15 * 24 * 4 # 15 days at 15 min frequency                           
    
def test_get_ivs_sw(nwis_instance):
    # check surface water site (discharge)
    sw_site = '07047810'
    df = nwis_instance.get_ivs(sw_site, 
                               parameter_code='00060', 
                               start_date='2002-08-01',
                               end_date='2002-08-15')

    assert df.columns[0] == 'Discharge, cubic feet per second' # check column header/units
    assert df.shape[0] == 15 # check daily resampling 

    # check without resampling
    df = nwis_instance.get_ivs(sw_site, 
                               parameter_code='00060', 
                               start_date='2002-08-01', 
                               end_date='2002-08-15', 
                               sample_period=None)
    assert df.shape[0] == 15 * 24 # 15 days at 1 hr frequency 

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
    #assert all_sites.site_no.dtype == np.object

def test_get_all_ivs_sw(nwis_instance):
    sites = ['07040450', '07047800']

    # check set of surface water sites from nwis instnace
    df_dct = nwis_instance.get_all_ivs(stations=sites, 
                                       start_date='2000-12-01', 
                                       end_date='2000-12-31')

    for site, df in df_dct.items():
        assert site in sites
        assert df.shape[0] == 31 # aggregated daily
        assert df.columns == 'Discharge, cubic feet per second'
    
    # check without resampling
    df_dct = nwis_instance.get_all_ivs(stations=sites, 
                                       start_date='2000-12-01', 
                                       end_date='2000-12-31',
                                       sample_period=None)
    for site, df in df_dct.items():
        assert site in sites
        assert 'datetime' in df.columns and 'Discharge, cubic feet per second' in df.columns
        df.datetime[0] == '2000-12-01 00:00' # verify has hh:mm


def test_get_all_ivs_gw(nwis_instance):
    sites = ['353606090510701', '353831090502401']

    # check set of groundwater sites from nwis instnace
    df_dct = nwis_instance.get_all_ivs(stations=sites, 
                                       parameter_code='72019', 
                                       start_date='2020-12-01', 
                                       end_date='2020-12-31')
    for site, df in df_dct.items():
        assert site in sites
        assert df.shape[0] == 31 # aggregated daily
        assert df.columns == 'Depth to water level, feet below land surface'

    df_dct = nwis_instance.get_all_ivs(stations=sites, 
                                       parameter_code='72019', 
                                       start_date='2020-12-01', 
                                       end_date='2020-12-31',
                                       sample_period=None)
    for site, df in df_dct.items():
        assert site in sites
        assert 'datetime' in df.columns and 'Depth to water level, feet below land surface' in df.columns
        df.datetime[0] == '2020-12-01 00:00' # verify has hh:mm