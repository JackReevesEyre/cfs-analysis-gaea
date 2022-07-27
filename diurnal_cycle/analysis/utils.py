import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime

def rotation_xy(ds, u_name, v_name, angles, angle_convention, exp_dim=None):
    """Gives streamwise and cross-stream wind components in 'angles' coordinate.
    Inputs:
        ds [xarray dataset] -- dataset containing eastward and
                               northward wind components.
        u_name, v_name [string] -- the names of the u and v components
                                   in ds (e.g., 'UWND', 'VWND')
        angles [numpy array or xarray data array] -- the angles of the
                                                     new coordinate system
                                                     to rotate to. Must
                                                     broadcast with ds.
                                                     MUST BE IN RADIANS.
        angle_convention ['math' or 'met'] -- the convention used in the
                                              definition of angles.
            'math' -- towards, anticlockwise from E
            'met' -- from, clockwise from N
        exp_dim [int or None] -- the dimension to expand at to make angles
                                 broadcast-able with velocity components.
                                 If None, the method will try some options 
                                 to make it work automatically.
    Output:
        ds_out [xarray dataset] -- the same as ds but with two new varibles:
                                   'U_streamwise' and 'U_crossstream':
                                   the components in the new coordinates.
    """
    # Convert the orientation of the new coordinate system.
    if angle_convention == 'met':
        A = np.radians(270.0) - angles
    elif angle_convention == 'math':
        A = angles
    else:
        sys.exit('ACHTUNG: angle_convention must be \'math\' or \'met\' in rotation_xy')
    
    # Make sure that the angle and original dataset broadcast.
    if exp_dim is not None:
        cosA = np.expand_dims(np.cos(A), axis=exp_dim)
        sinA = np.expand_dims(np.sin(A), axis=exp_dim)
    else:
        if (angles.shape == ()):
            cosA = np.cos(A)
            sinA = np.sin(A)
        elif (ds[u_name].shape == angles.shape):
            cosA = np.cos(A)
            sinA = np.sin(A)
        elif ((ds[u_name].shape[0] == angles.shape[0]) &
              (len(angles.shape) == 1)):
            cosA = np.expand_dims(np.cos(A), axis=1)
            sinA = np.expand_dims(np.sin(A), axis=1)
        elif ((ds[u_name].shape[1] == angles.shape[0]) &
              (len(angles.shape) == 1)):
            cosA = np.expand_dims(np.cos(A), axis=0)
            sinA = np.expand_dims(np.sin(A), axis=0)
        else:
            sys.exit('ACHTUNG: ds and angles are not compatible shapes in rotation_xy')
    
    # Do the rotation.
    ds_out = ds.copy()
    ds_out['U_streamwise'] = ds[u_name]*cosA + ds[v_name]*sinA
    ds_out['U_crossstream'] = -1.0*ds[u_name]*sinA + ds[v_name]*cosA

    # Add metadata.
    ds_out['U_streamwise'].attrs = {
        'units':ds[u_name].attrs['units'],
        'long_name':'Streamwise wind speed: in direction of mean wind.',
        'standard_name':'streamwise_wind_speed'
    }
    ds_out['U_crossstream'].attrs = {
        'units':ds[u_name].attrs['units'],
        'long_name':'Cross-stream wind speed: normal to mean wind in right-handed coordinate system.',
        'standard_name':'cross_stream_wind_speed'
    }
    return ds_out


def wind_dir_twice(ds, u_name, v_name):
    """Calculates two wind directions: meteorological & mathematical.
    Meteorological means compass bearing from which the wind comes.
    Mathematical means angle anticlockwise from east of the direction
    TOWARDS which the wind blows.

    Inputs:
        ds [xarray dataset] -- dataset containing eastward and
                               northward wind components.
        u_name, v_name [string] -- the names of the u and v components
                                   in ds (e.g., 'UWND', 'VWND')
    Output:
        ds_out [xarray dataset] -- the same as ds but with two new varibles:
                                   'wind_dir' (meteorological)
                                   'arctan_VoverU' (mathematical)
    Both outputs are in DEGREES.
    """
    ds_out = ds.copy()
    ds_out['arctan_VoverU'] = \
        np.degrees(np.arctan2(ds[v_name],ds[u_name]))
    ds_out['wind_dir'] = np.mod(
        270.0 - np.degrees(np.arctan2(ds[v_name],ds[u_name])),
        360
    )
    # Add metadata.
    ds_out['arctan_VoverU'].attrs = {
        'units':'degrees',
        'long_name':'Mathematical convention wind direction: angle towards which the wind blows, anticlockwise from eastward.',
        'standard_name':'wind_towards_direction'
    }
    ds_out['wind_dir'].attrs = {
        'units':'degrees',
        'long_name':'Meteorological convention wind direction: compass bearing from which the wind comes.',
        'standard_name':'wind_from_direction'
    }
    #
    return ds_out


def switch_lon_lims(lon_list, min_lon=0.0):
    """Shifts longitudes to be in any chosen arbitrary range.

    Returns a list of longitudes in the interval [min_lon, min_lon+360).
    """
    result = (lon_list - min_lon) % 360.0 + min_lon
    return result
