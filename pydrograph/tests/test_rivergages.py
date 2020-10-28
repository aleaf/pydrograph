from pydrograph import rivergages

def test_pull_year_of_data(site_name = 'rcki2', begin_date = '2016-01-01'):
    data = rivergages.pull_year_of_data(site_name, begin_date)
    assert len(data) > 2
    assert 'stage' in data.columns
    assert 'datetime' in data.columns
