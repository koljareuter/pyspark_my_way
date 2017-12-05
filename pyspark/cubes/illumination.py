from __future__ import print_function
import os
import numpy as np
from astropy.io import fits

from ..utils.fileutils import (get_obs_exposures, duplicate)
from ..utils.obsutils import (get_cal_rotator, match_obs_rotator)

__all__ = ["illuminate"]


def illuminate(logger=None, subdir=None, obsband=None, 
               search_suffix=None, **kwargs):

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
                                  obsband=obsband, suffix=search_suffix, **kwargs)

    
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
            
