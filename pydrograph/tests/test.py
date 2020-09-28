import fiona
from shapely.geometry import shape
import pytest
from gisutils import shp2df, project
from pydrograph.baseflow import get_upstream_area


@pytest.mark.skip(reason="need small test versions of input files")
def test_get_upstream_area():

    catchments = ['/Users/aleaf/Documents/NHDPlus/NHDPlusGL/NHDPlus04/NHDPlusCatchment/Catchment.shp',
                  '/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus07/NHDPlusCatchment/Catchment.shp']
    plusflow = ['/Users/aleaf/Documents/NHDPlus/NHDPlusGL/NHDPlus04/NHDPlusAttributes/PlusFlow.dbf',
                '/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus07/NHDPlusAttributes/PlusFlow.dbf']
    nodasites = '/Users/aleaf/Documents/USFS/Nicolet/targets/north/flux_field_no_da.shp'
    flowlines = ['/Users/aleaf/Documents/NHDPlus/NHDPlusGL/NHDPlus04/NHDSnapshot/Hydrography/NHDFlowline.shp',
                 '/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus07/NHDSnapshot/Hydrography/NHDFlowline.shp']
    nearfield = '/Users/aleaf/Documents/USFS/Nicolet/shps/Nicolet_north_NF.shp'

    nf = shape(fiona.open(nearfield).next()['geometry'])
    nf = project(nf, '+init=epsg:26716', '+init=epsg:4269')
    bbox = nf.bounds

    noda = shp2df(nodasites)

    get_upstream_area(noda.geometry.tolist(), plusflow, flowlines, catchments, nf)


if __name__ == '__main__':
    test_get_upstream_area()
