from __future__ import print_function
import os
import numpy as np
from astropy.io import fits

from ..utils.fileutils import (get_obs_exposures, duplicate)
from ..utils.obsutils import (get_cal_rotator, match_obs_rotator)

__all__ = ["illuminate"]


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

def illuminate(logger=None, subdir=None, obsband=None, **kwargs):

    #path for the output calibration directory
    caldir = os.path.join(kwargs['outputdir'], subdir, 'calib')

    #fetch rotator angles and filled extensions from the master flat file
    flat_file = fits.open(os.path.join(caldir, 'master_flat_{0}{0}{0}.fits'.format(obsband)))
    pro_naangles = get_cal_rotator(flat_file)
    nrot = len(pro_naangles)
    flat_file.close()

    # get list of frames to process for barycentric correction.  Assumes 
    # working with sky-subtracted frames in science directory
    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'science'), 
                                  obsband=obsband, suffix='s', **kwargs)

    
    # iterate over each exposure, applying closest angle illumination correction
    for exposure, _, _ in exposures:
        fpath = os.path.split(exposure)[0]

        #copy object to output
        exp_corr = duplicate(exposure, outputdir=fpath, suffix='i')
        
        _, match_rot, frot = match_obs_rotator(fits.getheader(exp_corr, ext=0), pro_naangles)
        illum_file = 'illum_corr_{0}{0}{0}_{1:03}.fits'.format(obsband.upper(), int(match_rot))

        with fits.open(exp_corr, mode='update') as ofile, \
                fits.open(os.path.join(caldir, illum_file)) as ifile:

            for ifu in range(1,25):
                cube = np.array(ofile[ifu].data)
                icorr = np.array(ifile[2*ifu-1].data)

                try:
                    oz, ox, oy = cube.shape
                except ValueError:
                    continue
                else:
                    ofile[ifu].data = cube/icorr

            ofile[0].header['HISTORY'] = "Applied illumination correction at {0:03} deg".format(int(match_rot))
            
            #for ifu in range(1,25):
            #    cube = np.array(ofile[ifu].data)
            #    try:
            #        oz, ox, oy = cube.shape
            #    except ValueError:
            #        continue
            #    else:
            #        wave = reconstruct_wave(ofile[ifu].header)
            #        corr_cube = rvcorr_cube(wave, cube, vbary)
            #        ofile[ifu].data = corr_cube



#def velocity_correct(logger=None, subdir=None, obsband=None, rvcorr=True,
#		     **kwargs):
#
#    # check that velocity correction is requested
#    if not rvcorr:
#        return
#
#    # get list of frames to process for barycentric correction.  Assumes 
#    # working with sky-subtracted frames in science directory
#    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'science'), 
#                                  obsband=obsband, suffix='s', **kwargs)
#
#    # iterate over each exposure, apply rv correction
#    for exposure, _, _ in exposures:
#        fpath = os.path.split(exposure)[0]
#
#        #copy object to output
#        exp_corr = duplicate(exposure, outputdir=fpath, suffix='rv')
#
#        with fits.open(exp_corr, mode='update') as ofile:
#            vbary = get_rvcorr(ofile[0].header) 
#            for ifu in range(1,25):
#                cube = np.array(ofile[ifu].data)
#                try:
#                    oz, ox, oy = cube.shape
#                except ValueError:
#                    continue
#                else:
#                    wave = reconstruct_wave(ofile[ifu].header)
#                    corr_cube = rvcorr_cube(wave, cube, vbary)
#                    ofile[ifu].data = corr_cube
#
#
#
#def rvcorr_exposure(kmos_frame, suffix='heliocorr'):
#    """ Apply barycentric velocity correction to every data cube in 
#    a given KMOS exposure.  Writes output into the same path with the
#    given suffix appended.
#
#    Parameters
#    ----------
#    kmos_frame: string
#        Full path to the KMOS exposure to be corrected
#
#    """
#    #set up path and file names
#    file_path = '/'.join(kmos_frame.split('/')[:-1])
#    file_name_in = kmos_frame.split('/')[-1]
#    file_name_out = '{0}_{1}.fits'.format(file_name_in[:-5], suffix)
#
#    #read in the primary fits header, add comments and initialize
#    #output fits HDU
#
#
#
#    return
#
#
