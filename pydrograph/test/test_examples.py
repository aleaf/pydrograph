import numpy as np
import pandas as pd
from shapely.geometry import box
import pytest
from gisutils import project, df2shp
from pydrograph import Nwis


@pytest.fixture(scope='session')
def extent_poly():
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


def test_nwis_from_shapefile(nwis_instance, extent_poly):
    """Check that supplying a shapefile directly in a different CRS
    results in the same extent as supplying a shapely polygon
    in the same CRS."""
    nwis = Nwis(extent='examples/data/bbox.shp')
    area_of_overlap = nwis.extent.intersection(nwis_instance.extent).area
    assert np.allclose(area_of_overlap, nwis_instance.extent.area)
    assert np.allclose(area_of_overlap, nwis.extent.area)
    assert np.allclose(area_of_overlap, extent_poly.area)
    return nwis


@pytest.fixture(scope='session')
def field_sites(nwis_instance):
    field_sites = nwis_instance.get_siteinfo('field_measurements')
    return field_sites


@pytest.fixture(scope='session')
def gw_field_sites(nwis_instance):
    field_sites = nwis_instance.get_siteinfo('gwlevels')
    return field_sites


@pytest.fixture(scope='session')
def dv_sites(nwis_instance):
    dv_sites = nwis_instance.get_siteinfo('daily_values')
    return dv_sites


@pytest.fixture(scope='session')
def gw_dv_sites(nwis_instance):
    dv_sites = nwis_instance.get_siteinfo('gw_daily_values')
    return dv_sites


def test_get_sw_sites(nwis_instance):
    nwis = nwis_instance
    field_sites = nwis.get_siteinfo('field_measurements')
    dv_sites = nwis.get_siteinfo('daily_values')
    assert isinstance(field_sites, pd.DataFrame)
    assert len(field_sites) > 0
    assert field_sites.site_no.dtype == np.object
    assert isinstance(dv_sites, pd.DataFrame)
    assert len(dv_sites) > 0
    assert dv_sites.site_no.dtype == np.object


def test_get_gw_sites(nwis_instance):
    nwis = nwis_instance
    gwfield_sites = nwis.get_siteinfo('gwlevels')
    gwdv_sites = nwis.get_siteinfo('gw_daily_values')
    assert isinstance(gwfield_sites, pd.DataFrame)
    assert len(gwfield_sites) > 0
    assert gwfield_sites.site_no.dtype == np.object
    assert isinstance(gwdv_sites, pd.DataFrame)
    assert len(gwdv_sites) > 0
    assert gwdv_sites.site_no.dtype == np.object


def test_get_daily_values_sw(nwis_instance, dv_sites):
    nwis = nwis_instance
    sites = dv_sites.site_no.tolist()[0:2]
    dvs = nwis.get_all_dvs(sites, start_date='1990-01-01')
    assert isinstance(dvs, dict)
    assert isinstance(dvs['04021520'], pd.DataFrame)


def test_get_daily_values_gw(nwis_instance, gw_dv_sites):
    nwis = nwis_instance
    sites = gw_dv_sites.site_no.tolist()[0:2]
    dvs = nwis.get_all_dvs(sites, 'gwlevels', start_date='1990-01-01')
    assert isinstance(dvs, dict)
    assert isinstance(dvs['464222092403801'], pd.DataFrame)


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


def test_get_field_measurements(nwis_instance, field_sites):
    nwis = nwis_instance
    sites = field_sites.site_no.tolist()[:5]
    fm = nwis.get_all_measurements(sites)
    assert isinstance(fm, pd.DataFrame)
    assert fm.site_no.dtype == np.object


def test_get_gw_field_measurements(nwis_instance, gw_field_sites):
    nwis = nwis_instance
    sites = gw_field_sites.site_no.tolist()[:5]
    fm = nwis.get_all_measurements(sites, txt='gwlevels')
    assert isinstance(fm, pd.DataFrame)
    assert fm.site_no.dtype == np.object





