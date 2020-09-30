import numpy as np
import pandas as pd
import pytest
from pydrograph.baseflow import IHmethod


@pytest.fixture(scope='function')
def test_data(test_data_path):
    data = pd.read_csv(test_data_path / 'UV_04087088_Discharge_20071001_ab.csv')
    data.index = pd.to_datetime(data.date)
    return data


@pytest.mark.parametrize('freq', ('D', '6H'))
@pytest.mark.parametrize('interp_semilog', (False, True))
@pytest.mark.parametrize('block_length', [1, 2, 3])
def test_IHmethod(test_data, block_length, interp_semilog, freq):
    results = IHmethod(test_data.Q, block_length=block_length, tp=0.9,
                       freq=freq,
                       interp_semilog=interp_semilog)
    minimum_points = ~results.block_Qmin.isna()
    assert np.all(results.loc[minimum_points, 'block_Qmin'] <= results.loc[minimum_points, 'Q'])


@pytest.mark.parametrize('data', (pytest.param(pd.Series(),
                                               marks=pytest.mark.xfail(reason="index isn't datetime")),
                                  pytest.param(pd.Series(index=pd.to_datetime([])),
                                               marks=pytest.mark.xfail(reason="index isn't datetime")),
                                  pytest.param(pd.Series([15.2, 14.8, 14.5, 14.2],
                                                         index=pd.date_range('2020-09-30', '2020-10-03')),
                                               marks=pytest.mark.xfail(reason="index isn't datetime")),
                                  pd.Series([15.2, 14.8, 14.5, 14.2, 13.9, 13.8, 13.8, 13.8, 13.8, 13.8],
                                            index=pd.date_range('2020-01-01', '2020-01-10'))
                                  )
                         )
def test_IHmethod_with_not_enough_data(data):
    results = IHmethod(data)
    j=2




