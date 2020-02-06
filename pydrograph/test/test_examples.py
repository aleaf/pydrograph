import pandas as pd
import fiona
from shapely.geometry import shape, box
import pytest
from gisutils import project, df2shp
from pydrograph import Nwis


@pytest.fixture(scope='session')
def extent_poly():
    #extent_shp = 'examples/data/muk_nf.shp' # polygon of study area
    #epsg = 26915

    #extent_poly = shape(fiona.open(extent_shp).next()['geometry'])
    #extent_poly_ll = project(extent_poly, "+init=epsg:{}".format(epsg), "+init=epsg:4269")
    extent_poly_ll = box(-92.7, 46.7, -92.6, 46.8)

    extent_poly = project(extent_poly_ll, "+init=epsg:{}".format(4269), "+init=epsg:26915")
    df = pd.DataFrame({'geometry': [extent_poly],
                       'id': [0]})
    df2shp(df, 'examples/data/bbox.shp', epsg=26915)
    return extent_poly_ll


@pytest.fixture(scope='session')
def nwis_instance(extent_poly):
    nwis = Nwis(extent=extent_poly)
    return nwis


def test_instantiate_with_polygon_sw(nwis_instance):
    nwis = nwis_instance
    assert isinstance(nwis.field_sites, pd.DataFrame)
    assert len(nwis.field_sites) > 0
    assert isinstance(nwis.dv_sites, pd.DataFrame)
    assert len(nwis.dv_sites) > 0


def test_instantiate_with_polygon_gw(nwis_instance):
    nwis = nwis_instance
    assert isinstance(nwis.gwfield_sites, pd.DataFrame)
    assert len(nwis.gwfield_sites) > 0
    assert isinstance(nwis.gwdv_sites, pd.DataFrame)
    assert len(nwis.gwdv_sites) > 0


def test_get_daily_values_sw(nwis_instance):
    nwis = nwis_instance
    sites = nwis.dv_sites.site_no.tolist()[0:2]
    dvs = nwis.get_all_dvs(sites, start_date='1990-01-01')
    assert isinstance(dvs, dict)
    assert isinstance(dvs[4021520], pd.DataFrame)


def test_get_daily_values_gw(nwis_instance):
    nwis = nwis_instance
    sites = nwis.gwdv_sites.site_no.tolist()[0:2]
    dvs = nwis.get_all_dvs(sites, 'gwlevels', start_date='1990-01-01')
    assert isinstance(dvs, dict)
    assert isinstance(dvs[464222092403801], pd.DataFrame)


def test_get_single_site_sw(nwis_instance):
    nwis = nwis_instance
    df = nwis.get_dvs(4021520)
    assert isinstance(df, pd.DataFrame)


def test_get_single_site_gw(nwis_instance):
    nwis = nwis_instance
    df = nwis.get_dvs(464322092401401, 'gwlevels')
    assert isinstance(df, pd.DataFrame)


def test_make_url(nwis_instance):
    nwis = nwis_instance
    url = nwis.make_dv_url(4015475)
    assert isinstance(url, str)


def test_make_url_gw(nwis_instance):
    nwis = nwis_instance
    url = nwis.make_dv_url(464322092401401, parameter_code=72019)
    assert isinstance(url, str)


def test_get_field_measurements(nwis_instance):
    nwis = nwis_instance
    sites = nwis.field_sites.site_no.tolist()[:5]
    fm = nwis.get_all_measurements(sites)
    assert isinstance(fm, pd.DataFrame)


def test_get_gw_field_measurements(nwis_instance):
    nwis = nwis_instance
    sites = nwis.gwfield_sites.site_no.tolist()[:5]
    fm = nwis.get_all_measurements(sites, txt='gwlevels')
    assert isinstance(fm, pd.DataFrame)





