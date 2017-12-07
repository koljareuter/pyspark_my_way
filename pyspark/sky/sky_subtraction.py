from __future__ import print_function

import subprocess
import numpy as np
import os, sys
import time
import bottleneck as bn
import astropy.io.fits as fits

#from .zapv2 import (process, SVDoutput)
from .zapv2 import process as process_v2
from .zapv2 import SVDoutput as SVDoutput_v2
from .zapv1 import process as process_v1
from .zapv1 import SVDoutput as SVDoutput_v1

from ..utils.obsutils import find_lines
from ..utils.fileutils import (fetch_resource_path, duplicate, get_obs_exposures)

__all__ = ["subtract_sky", "match_obj_sky", "make_SVDs_v2"]

sky_methods = ['zapv1', 'zapv2']



def match_obj_sky(exposures):
    """ Given a list of exposures generate by get_obs_exposures, return
    a list of object,sky tuples, where the matching is done by closest time.
    """
    #pull out times
    time_dict = dict([(efile, etime) for efile, _, etime in exposures])

    #lookup table of times
    object_exposures = [efile for efile, etype, _ in exposures if etype == 'OBJ'] 
    sky_exposures = [efile for efile, etype, _ in exposures if etype == 'SKY'] 
    sky_times = [time_dict[efile] for efile in sky_exposures]

    obj_sky_pair = []
    for obj in object_exposures:
        sky_idx = np.argmin(np.fabs(np.array(sky_times) - time_dict[obj]))
        obj_sky_pair.append((obj, sky_exposures[sky_idx]))

    return obj_sky_pair


def make_SVDs_v2(spec_array, header, pcaoutfits=None, zlevel='none', cftype='med_smooth', 
              clean=False, cfwidth=101, boost=50, reject_factor=0.4, logger=None):
    """ Given an input (2D) array of spectra, generate PCA output following
    zap_v2.0 conventions.
    """
    #generate mask array for sky line rejection
    med_sky = bn.nanmedian(spec_array, axis=1)
    msky = find_lines(med_sky)

    #run zap instance
    zap_inst = SVDoutput_v2(None, cube='dummy', stack=spec_array, hdr=header, instrument='KMOS', 
                         zlevel=zlevel, cftype=cftype, clean=False, cfwidth=cfwidth, boost=boost, 
                         reject_factor=reject_factor, contmask=msky)

    if pcaoutfits is not None:
        zap_inst.writePCA(pcaoutputfits=pcaoutfits, overwrite=True)

   
def zap_sky_v2(matched_exposures, logger=None, obsband=None, subdir=None, 
            inputdir=None, suffix='s', **kwargs):

    #fetch the superdark reference frame for this detector
    svd_file = 'test_ZAPv2_SVD_{0}_n500004.fits'.format(obsband.upper())
    fpath = fetch_resource_path('data/{0}'.format(svd_file))
    if logger is not None:
        logger.debug("- Subtracting sky using {0}".format(svd_file))
    
    for obj_exp, sky_exp in matched_exposures:
        input_path = os.path.split(obj_exp)[0]
        output_path = os.path.join(inputdir, subdir, 'science')

        # copy object frame to output
        exp_sub = duplicate(obj_exp, outputdir=output_path, suffix=suffix) 

        with fits.open(exp_sub, mode='update') as ofile, fits.open(sky_exp) as sfile:

            for ifu in range(1,25):
                scube = np.array(sfile[ifu].data)
                ocube = np.array(ofile[ifu].data)
                header = ofile[ifu].header

                try:
                    oz, ox, oy = ocube.shape
                except ValueError:
                    continue
                else:
                    if logger is not None:
                        logger.info('-- Processing sky subtraction for IFU {0}'.format(ifu))
                    sub_cube = np.copy(ocube) - np.copy(scube)
                    nanmask = np.isfinite(sub_cube)
                    sub_cube[~nanmask] = 0.

                    #generate mask array for continuum
                    msky = find_lines(bn.nanmedian(scube[:,1:13,1:13].reshape(2048,-1), axis=1))

                    zap_obj = process_v2('dummy', cube=sub_cube, hdr=header, clean=False, 
                                        zlevel='none', cftype='med_smooth', cfwidthSP=101, 
                                        extSVD=fpath, interactive=True, overwrite=True,
                                        instrument='KMOS', boost=50, contmask=msky,
                                        reject_factor=0.4, logger=logger)
                
                    corr_cube = np.copy(zap_obj.cleancube)
                    corr_cube[~nanmask] = np.nan

                    ofile[ifu].data = corr_cube


def zap_sky_v1(matched_exposures, logger=None, obsband=None, subdir=None, 
            inputdir=None, suffix='s', **kwargs):

    #fetch the superdark reference frame for this detector
    svd_file = 'ZAPv1_SVD_{0}_n500004.fits'.format(obsband.upper())
    fpath = fetch_resource_path('data/{0}'.format(svd_file))
    if logger is not None:
        logger.debug("- Subtracting sky using {0}".format(svd_file))

    for obj_exp, sky_exp in matched_exposures:
        input_path = os.path.split(obj_exp)[0]
        output_path = os.path.join(inputdir, subdir, 'science')

        # copy object frame to output
        exp_sub = duplicate(obj_exp, outputdir=output_path, suffix=suffix) 

        with fits.open(exp_sub, mode='update') as ofile, fits.open(sky_exp) as sfile:

            for ifu in range(1,25):
                scube = np.array(sfile[ifu].data)
                ocube = np.array(ofile[ifu].data)
                header = ofile[ifu].header

                try:
                    oz, ox, oy = ocube.shape
                except ValueError:
                    continue
                else:
                    if logger is not None:
                        logger.info('-- Processing sky subtraction for IFU {0}'.format(ifu))
                    sub_cube = np.copy(ocube) - np.copy(scube)
                    nanmask = np.isfinite(sub_cube)
                    sub_cube[~nanmask] = 0.

                    #generate mask array for continuum
                    msky = find_lines(bn.nanmedian(scube[:,1:13,1:13].reshape(2048,-1), axis=1))

                    zap_obj = process_v1('dummy', cube=sub_cube, hdr=header, clean=False, 
                                        zlevel='none', cftype='med_smooth', cfwidthSP=51, 
                                        extSVD=fpath, interactive=True, overwrite=True,
                                        instrument='KMOS', boost=50, contmask=msky,
                                        reject_factor=0.4, logger=logger)
                
                    corr_cube = np.copy(zap_obj.cleancube)
                    corr_cube[~nanmask] = np.nan

                    ofile[ifu].data = corr_cube


def subtract_sky(logger=None, subdir=None, obsband=None, 
                 sky_method='zapv2', search_suffix=None, **kwargs):
    """ Driver function for sky subtraction.  Currently only supports
    subtraction using ZAP, but can be expanded to include other approaches.
    """
    #check that sky subtraction is handled by this routine
    if sky_method not in sky_methods:
        return 0

    if logger is not None:
        logger.debug('- Running sky subtraction')

    #get list of frames to process. Assumes location of un-sky-subtracted frames in "sky" directory
    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'sky'), obsband=obsband, 
                                  suffix=search_suffix, **kwargs) 

    if not len(exposures):
        if logger is not None:
            logger.warning('-- No reconstructed frames found, skipping sky subtraction')
        return 1

    #generate associated list of objects and skies
    obj_sky = match_obj_sky(exposures)
    
    if sky_method == 'zapv2':
        zap_sky_v2(obj_sky, logger=logger, obsband=obsband, subdir=subdir, **kwargs)    

    elif sky_method == 'zapv1':
        zap_sky_v1(obj_sky, logger=logger, obsband=obsband, subdir=subdir, **kwargs)    

    return 0
