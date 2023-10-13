import sys
import pyproj
import numpy as np
import datetime as dt
from dateutil import tz

def utc_to_local(utc, zone='Europe/Paris'):
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz(zone)

    # Tell the datetime object that it is in UTC since
    # datetime object are 'naive' by default.
    utc = utc.replace(tzinfo=from_zone)

    # Convert time zone
    return utc.astimezone(to_zone)

# Requires holidays python library (https://pypi.org/project/holidays/)
def is_holiday(date, country_code):
    # Try import holidays
    try:
        import holidays
    except ImportError:
        print('WARNING: holidays library could not be imported.')
        return False
    # Return True if it is a holiday
    country_found = False
    for country in holidays.list_supported_countries():
        if country_code == country:
            country_found = True
            country_holidays = holidays.CountryHoliday(country_code)
            break
    if country_found == False:
        print('WARNING: country_code {} not supported.'.format(country_code))
        return False

    return date.date() in country_holidays


# From http://www.johndcook.com/blog/python_longitude_latitude/

import math

def distance_on_unit_sphere(lat1, long1, lat2, long2):
 
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.

    earth_radius = 6371000 # in m

    return arc * earth_radius
        
