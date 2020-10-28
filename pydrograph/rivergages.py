__author__ = 'lschachter'

import datetime as dt
import time
from urllib.request import urlopen
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, shape
import pyproj
import gisutils
import urllib3
import xmltodict
import datetime

def make_url(site_code, method, begin_date = '2000-01-01', begin_time = '00:00'):
    ''' This function creates a url to interact with the Army Corps of Engineers webservices to access their rivergages data.
        Currently, the function only supports stage data, as that is what the code is written to pull; however, 
        it could be modified to pull in flow, precip, air/water temp, and other parameters as well that are collected by
        the Army Corps of Engineers.
        
        This function takes in a site code, which is usually alphanumeric such as "rcki2" and should be in a string format.
        The "method" that this function takes in should be one of the five methods that the webservice can interpret:
            getSites, getValues, getVariableInfo, getSiteInfo, and getFValues.
        
        The begin date is set to pull the January 2000, however these can be changed to start at any date and time. The function
        pulls one month of data after the start date and time, as the service seems to only return a month of data regardless of how
        far out the end date is set.
        
        Dates should be in 'YYYY-MM-DD' format, and time should be in 'HH:MM' 24 format. (data is hourly when available)
        
        Parameters
        ----------
        site_code = str
        
        method = str
          'getSites',
          'getValues',
          'getVariableInfo',
          'getSiteInfo',
          'getFValues'
          
        begin_date = str
            'YYYY-MM-DD'
        
        begin_time = str
            'HH:MM'
        
        Returns
        -------
        url string  
    '''
    begin = begin_date + 'T' + begin_time
    end = datetime.datetime.strptime(begin_date, '%Y-%m-%d')
    end = end + datetime.timedelta(days = 31)
    end = end.strftime('%Y-%m-%dT00:00')
    
    meth_list = ['getSites',
          'getValues',
          'getVariableInfo',
          'getSiteInfo',
          'getFValues']
    
    stuff_at_start = 'https://rivergages.mvr.usace.army.mil/watercontrol/webservices/rest/webserviceWaterML.cfc?method=RGWML&'
    method_url = 'meth={}'.format('getSiteInfo')
    stuff_in_middle = '&site={}&location={}&variable=HG&'.format(site_code, site_code)
    dates = 'beginDate={}&endDate={}'.format(begin, end)
    stuff_at_end = '&authtoken=RiverGages&authToken=RiverGages'
    url = stuff_at_start + method_url + stuff_in_middle + dates + stuff_at_end
    
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    data = xmltodict.parse(response.data)
    sitename = (data['sitesResponse']['site']['siteInfo']['siteName'])
    
    words = ['lock', 'Lock', 'dam', 'Dam']
    
    dam = 0
    for i in words:
        if i in sitename:
            dam = dam + 1
            #print('might be a Dam! [{}]'.format(i))
    
    method_url2 = 'meth={}'.format(method)
    
    if dam > 0:
        stuff_in_middle = '&site={}&location={}&variable=HT&'.format(site_code, site_code)
    
    url = stuff_at_start + method_url2 + stuff_in_middle + dates + stuff_at_end
    
    return url
    
    
def get_data(url):
    ''' This function takes a url made in the make_url function and returns the raw data, which is a ordered dictionary.
    
    Parameters
    ----------
    url = str
    
    Returns
    -------
    data = ordered dictionary
    
    '''

    url = url

    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    return data

def get_site_info(url, data_return = False):
    ''' This function prints site information about the site name, site code, available date range, geolocation, and coordinate reference system.
    This function will only run if the url was made using the method 'getSiteInfo'. With the correct url, this function pulls in variable information from the
    Army Corps of Engineers webservices. 
    
        Parameters
        ----------
        url = str
        data_return = optional, set to False, True/False
        
        Returns
        -------
        prints site information about the site name, site code, available date range, geolocation, and coordinate reference system. 
        If data_return is set to True, will return raw data (ordered dictionary).
    
    '''
    url = url

    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    sitename = (data['sitesResponse']['site']['siteInfo']['siteName'])
    sitecode = data['sitesResponse']['site']['siteInfo']['siteCode']['#text']
    geolocation = dict(data['sitesResponse']['site']['siteInfo']['geoLocation']['geogLocation'])
            
    print(sitename)
    print('---------------------------------------------------------------')
    print('Site Code: {}'.format(sitecode))
    print('Location: latitude {}, longitude {}'.format(geolocation['latitude'], geolocation['longitude']))
    print('Coordinate Reference System: {}'.format(geolocation['@srs']))
    
    
    if len(data['sitesResponse']['site']['seriesCatalog']) == 1:
            print('Variable: No variable info available at this site.')
            print('Available Date Range: No date info available at this site.')
            variable_info = None
    elif isinstance(data['sitesResponse']['site']['seriesCatalog']['series'], list):

        variable_info = dict(data['sitesResponse']['site']['seriesCatalog']['series'][1]['variable'])
        begin_date = data['sitesResponse']['site']['seriesCatalog']['series'][1]['variableTimeInterval']['beginDateTime']
        begin_date = begin_date.replace('T', ' ')
        end_date = data['sitesResponse']['site']['seriesCatalog']['series'][1]['variableTimeInterval']['endDateTime']
        end_date = end_date.replace('T', ' ')
        
    else:    
        
        variable_info = dict(data['sitesResponse']['site']['seriesCatalog']['series']['variable'])
        begin_date = data['sitesResponse']['site']['seriesCatalog']['series']['variableTimeInterval']['beginDateTime']
        begin_date = begin_date.replace('T', ' ')
        end_date = data['sitesResponse']['site']['seriesCatalog']['series']['variableTimeInterval']['endDateTime']
        end_date = end_date.replace('T', ' ')

    if variable_info is not None:
        print('Variable: {}'.format(variable_info['variableName']))
        print('Available Date Range: {} to {}'.format(begin_date, end_date))
    
    if data_return is True:
        return data
    
def get_variable_info(url, data_return = False):
    ''' This function prints information about the variable returned by 'getVariableInfo' from the Army Corps of Engineers webservice.
    This function will only run if the url was made using the method 'getVariableInfo'.
    
        Parameters
        ----------
        url = str
        data_return = optional, set to False, True/False
        
        Returns
        -------
        prints variable information about the variable, variable code, and no data value. 
        If data_return is set to True, will return raw data (ordered dictionary).
        
        variable_info - pandas dataframe of variable information
        variable_units - dictionary of variable units 
        variable_codes - dictionary of variable codes
        
    '''
    url = url

    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    variable_info = (data['variablesResponse']['variables'])
    variable_units = dict(data['variablesResponse']['variables']['variable']['units'])
    variable_codes = dict(data['variablesResponse']['variables']['variable']['variableCode'])
    
    print('Variable: {}'.format(variable_info['variable']['variableName']))
    print('Variable Code: {}'.format(variable_info['variable']['variableCode']['#text']))
    print('No Data Value: {}'.format(variable_info['variable']['NoDataValue']))
    
    if data_return is True:
        return data
    
def get_sites(url, data_return = False):
    ''' This function prints information about the site by 'getSites' from the Army Corps of Engineers webservice.
    This function will only run if the url was made using the method 'getSites'.
        
        Parameters
        ----------
        url = str
        data_return = optional, set to False, True/False
        
        Returns
        -------
        prints site information about site name, side code, location, and coordinate reference system.
        If data_return is set to True, will return raw data (ordered dictionary).
        
    '''
    url = url
    
    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    sitename = (data['sitesResponse']['site']['siteInfo']['siteName'])
    sitecode = dict(data['sitesResponse']['site']['siteInfo']['siteCode'])
    geolocation = dict(data['sitesResponse']['site']['siteInfo']['geoLocation']['geogLocation'])
    
    print(sitename)
    print('---------------------------------------------------------------')
    print('Site Code: {}, {}'.format(sitecode['#text'], sitecode['@siteID']))
    print('Location: latitude {}, longitude {}'.format(geolocation['latitude'], geolocation['longitude']))
    print('Coordinate Reference System: {}'.format(geolocation['@srs']))
    
    if data_return is True:
        return data
    
def get_values(url):
    ''' This function retrieves a month of stage data from the US Army Corps of Engineers webservices. Data is returned in a pandas
    dataframe with stage and datetime columns. This function will only run if the url was made using the method 'getValues'.
    
        Parameters
        ----------
        url = str
        
        Returns
        -------
        df - pandas dataframe of datetimes and stage data (as strings and as floats) in feet
        
    '''
    url = url

    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    if len(data['timeSeriesResponse']['timeSeries']['values']) == 2:
        df2 = pd.DataFrame(data['timeSeriesResponse']['timeSeries']['values']['value'])
        df2['datetime'] = pd.to_datetime(df2['@dateTime'], format= '%Y-%m-%dT%H:%M:%S')
        df2 = df2.drop(columns = {'@dateTime'})
        df2 = df2.rename(columns = {'#text' : 'stage (ft)'})
        df2 = df2[df2['stage (ft)'] >= '0']
        df2 = df2.assign(stage = df2['stage (ft)'].astype(float))
    
    else:
        print('No values available at this site in this date range')
        df2 = None
    
    return df2
    
def get_fvalues(url):
    ''' This function retrieves a month of stage data from the US Army Corps of Engineers webservices. Data is returned in a pandas
    dataframe with stage and datetime columns. This function will only run if the url was made using the method 'getValues'.
    
        Parameters
        ----------
        url = str
        
        Returns
        -------
        df - pandas dataframe of datetimes and stage data (as strings and as floats) in feet
       
    '''
    url = url

    http = urllib3.PoolManager()

    response = http.request('GET', url)

    data = xmltodict.parse(response.data)
    
    if len(data['timeSeriesResponse']['timeSeries']['values']) == 2:
        df2 = pd.DataFrame(data['timeSeriesResponse']['timeSeries']['values']['value'])
        df2['datetime'] = pd.to_datetime(df2['@dateTime'], format= '%Y-%m-%dT%H:%M:%S')
        df2 = df2.drop(columns = {'@dateTime'})
        df2 = df2.rename(columns = {'#text' : 'stage (ft)'})
        df2 = df2[df2['stage (ft)'] >= '0']
        df2 = df2.assign(stage = df2['stage (ft)'].astype(float))
    
    else:
        print('No values available at this site in this date range')
        df2 = None
    
    return df2
    
def pull_year_of_data(site_name, begin_date):
    '''This function uses the get_values function to pull one year of data (the Army site only gives a month at a time).
    
        Parameters
        ----------
        site_name = str
        begin_date = str 'YYYY-MM-DD'
        
        Returns
        -------
        data - pandas dataframe of datetimes and stage data (as strings and as floats) in feet
       
    '''
    begin_time = '00:00'
    begin = begin_date + 'T' + begin_time
    end = datetime.datetime.strptime(begin_date, '%Y-%m-%d')
    urls = [None]*12
    urls[0] = make_url(site_name, 'getValues', begin_date = begin_date, begin_time = '00:00')
    for i in range(0, 11):
        end = end + datetime.timedelta(days = 31)
        end_datetime = end.strftime('%Y-%m-%d')
        url = make_url(site_name, 'getValues', begin_date = end_datetime, begin_time = '00:00')
        urls[i+1] = url
    data = get_values(urls[0])
    for i in range(1, len(urls)):
        data2 = get_values(urls[i])
        data = pd.concat([data, data2])
    
    year = begin_date[0:4]
    next_year = str(int(year) + 1) + '-01-01'
    
    data = data[data['datetime'] < next_year]
    data = data.sort_values(by='datetime')
    
    return data