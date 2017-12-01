""" This code corrects an observed KMOS exposure for barycentric velocity.

The core code was adapted from Jane Rigby's code for computing barycentric
corrections, which was in turn adapted from
https://gist.github.com/StuartLittlefair/5aaf476c5d7b52d20aa9544cfaa936a1
"""
from __future__ import print_function, division

import numpy as np
import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.coordinates import (SkyCoord, EarthLocation, solar_system, 
        UnitSphericalRepresentation, CartesianRepresentation)
from astropy.time import Time


#from ..utils.sinc_interpolation import sinc_interpolate

__all__ = ["rvcorr_cube", "rvcorr_exposure", "get_rvcorr", "calc_vbary"]


#def apply_barycentric_correction(df, barycor_vel,  colwav='wave', colwavnew='newwave'):
#    ''' Applies the barycentric velocity correction to the wavelength of a spectrum.  Assumes Pandas.
#    Inputs:
#    df:           the spectrum, as a pandas DataFrame
#    barycor_vel:  barycentric correction velocity, from compute_barycentric_correction()
#    colwav:       column name that holds the wavelength
#    colwavenew:   column name for new wavelength.  Have not tested what happens if you try to overwrite colwav'''
#    A_c = constants.c  # speed of light
#    df[colwavnew] = df[colwav] * (1.0 + (barycor_vel / A_c))
#    return(0)

#fhead = fits.getheader('/Users/jtmendel/astro/kmos_redux/new_dc/KMOS_SPEC_OBS304_0024.fits')
#print( get_velcorr(fhead))

def calc_vbary(time, skycoord, location):
    """Barycentric velocity correction.  was named 'velcorr'. Should return
    barycentric correction factor, in km/s, same sign convention as in jskycalc 
    
    Uses the ephemeris set with
    ``astropy.coordinates.solar_system_ephemeris.set`` for corrections.  For
    more information see `~astropy.coordinates.solar_system_ephemeris`.
    
    Parameters
    ----------
    mytime : `~astropy.time.Time`
      The time of observation.
    skycoord: `~astropy.coordinates.SkyCoord`
      The sky location to calculate the correction for.
    location: `~astropy.coordinates.EarthLocation`, optional
      The location of the observatory to calculate the correction for.
      If no location is given, the ``location`` attribute of the Time
      object is used
      
    Returns
    -------
    vel_corr : `~astropy.units.Quantity`
      The velocity correction to convert to Barycentric velocities. Should be
      added to the original velocity.
    """
    # ICRS position and velocity of Earth's geocenter     
    ep, ev = solar_system.get_body_barycentric_posvel('earth', time)
    
    # GCRS position and velocity of observatory
    op, ov = location.get_gcrs_posvel(time)
    
    # ICRS and GCRS are axes-aligned. Can add the velocities
    velocity = ev + ov # relies on PR5434 being merged
    
    # get unit ICRS vector in direction of SkyCoord
    sc_cartesian = skycoord.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
   
    return sc_cartesian.dot(velocity).to(u.km/u.s) # similarly requires PR5434


def get_rvcorr(header):
    """Computes barycentric velocity correction.

    Parameters
    ----------
    header : dict
        Fits header object (from, for example `~astropy.io.fits`) from
        KMOS observations

    Returns
    -------
    vbary: `~astropy.units.Quantity`
        Barycentric velocity correction for the input frame
    """
    date_obs = header['DATE-OBS']
    ra = header['RA']
    dec = header['DEC']
    
    lat = header['ESO TEL GEOLAT']
    lon = header['ESO TEL GEOLON']
    elev = header['ESO TEL GEOELEV']

    location = EarthLocation(lon=lon, lat=lat, height=elev)
    time = Time(date_obs, format='isot', scale='utc')
    skycoord = SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg), obstime=time)

    vbary = calc_vbary(time, skycoord, location)
    return vbary


def rvcorr_cube(wave_in, array_in, vbary):
    """ Apply barycentric velocity correction to input KMOS data cube

    Parameters
    ----------
    wave_in : ndarray
        Input wavelength array.
    array_in : ndarray
        Input datacube, oriented as (wavelength, x, y).
    vbary: `~astropy.units.Quantity`
        Barycentric velocity, computed from get_velcorr.  Should be
        added to the original velocity. 

    Returns
    -------
    array_out: ndarray
        Datacube corrected for barycentric velocity.
    """

    wave_rvcorr = wave_in * (1.0 + (vbary / c.c))

    #interpolate data from shifted wavelength grid back on to the input grid
    _, nx, ny = array_in.shape
    array_out = np.full_like(array_in, np.nan)
    goodpix_in = np.isfinite(array_in)
    for ix in range(nx):
        for iy in range(ny):
            old_spec = array_in[:, ix, iy]
            goodpix = goodpix_in[:, ix, iy]
            if np.count_nonzero(goodpix):
                new_spec = sinc_interpolate(wave_in, wave_rvcorr, old_spec)
                new_spec[~goodpix | ~np.isfinite(new_spec)] = np.nan
                array_out[:,ix,iy] = new_spec
    
    return array_out


def rvcorr_exposure(kmos_frame, suffix='heliocorr'):
    """ Apply barycentric velocity correction to every data cube in an
    a given KMOS exposure.  Writes output into the same path with the
    given suffix appended.

    Parameters
    ----------
    kmos_frame: string
        Full path to the KMOS exposure to be corrected

    """

    #set up path and file names
    file_path = '/'.join(kmos_frame.split('/')[:-1])
    file_name_in = kmos_frame.split('/')[-1]
    file_name_out = '{0}_{1}.fits'.format(file_name_in[:-5], suffix)

    #read in the primary fits header, add comments and initialize
    #output fits HDU



    return


