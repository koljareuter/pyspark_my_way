from __future__ import print_function

import subprocess,copy
import numpy as np
import os, sys
import time
import bottleneck as bn
import astropy.io.fits as fits

import astropy.units as u
from scipy.optimize import minimize
#from astropy.modeling.blackbody import blackbody_lambda
from astropy.modeling import models

#from .zapv2 import (process, SVDoutput)
from .zapv2 import process as process_v2
from .zapv2 import Zap as ZapObj_v2
from .zapv2 import SVDoutput as SVDoutput_v2
from .zapv1 import process as process_v1
from .zapv1 import SVDoutput as SVDoutput_v1

from ..utils.obsutils import (find_lines, hsm, reconstruct_wave)
from ..utils.fileutils import (fetch_resource_path, duplicate, get_obs_exposures)


__all__ = ["subtract_sky", "match_obj_sky", "make_SVDs_v2"]

sky_methods = ['zapv1', 'zapv2']



def match_obj_sky(exposures):
    """ Given a list of exposures generate by get_obs_exposures, return
    a list of object,sky tuples, where the matching is done by closest time.
    """
    #pull out times
    time_dict = dict([(efile, etime) for efile, _, etime, _ in exposures])

    #lookup table of times
    object_exposures = [(efile, dit) for efile, etype, _, dit in exposures if etype == 'OBJ'] 
    sky_exposures = [(efile, dit) for efile, etype, _, dit in exposures if etype == 'SKY'] 
    #sky_times = [time_dict[efile] for efile, _ in sky_exposures]
    #sky_dits = np.array([dit for _, dit in sky_exposures])

    obj_sky_pair = []
    for obj, odit in object_exposures:
        #only consider skys matched by exposure time
        matched_skys = [efile for efile, dit in sky_exposures if dit == odit]
        matched_times = [time_dict[efile] for efile, dit in sky_exposures if dit == odit]

        sky_idx = np.argmin(np.fabs(np.array(matched_times) - time_dict[obj]))
        obj_sky_pair.append((obj, matched_skys[sky_idx], odit))

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


def sub_thermal(wave, cube, mask):
    #generate 1D spectrum from cube
    spec_1d = bn.nanmedian(cube[:,1:13,1:13].reshape(2048,-1), axis=1)

    #remove pixels with dodgy values in advance
    clean_mask = np.isfinite(spec_1d) & ~np.isnan(spec_1d) & (wave > 2.1) & (wave < 2.4)
    spec_1d = spec_1d[clean_mask]
    wave_init = wave[clean_mask]
    mask_init = ~mask[clean_mask]
   
    #iteratively fit for thermal component
    out_init = np.array([315., 3.30, 10.])
    for fiter in range(3):

        #initial fit using minizming the modal residual
        sol = minimize(res_bb, out_init, args=(wave_init, spec_1d, mask_init), method='L-BFGS-B', 
                            bounds=[(270,380), (0,10),(-200,300)])

        #store the output solution
        out_init = sol.x

        #generate best fit and modify the mask
        temp_bb = bb_model(out_init, wave_init)
        res = (spec_1d - temp_bb)[mask_init]
        med_offset = bn.nanmedian(res)
        var_offset = 1.4862 * bn.nanmedian(np.fabs(res-med_offset))
        
        #update mask
        mask_init = ((spec_1d - temp_bb) < 2.*var_offset + med_offset)

    ##remove subtract fitted thermal component and return
    fit_bb = bb_model([out_init[0], out_init[1], 0.], wave)
    corrected = np.copy(cube) - fit_bb[:,np.newaxis,np.newaxis]

    return corrected


def bb_model(params, wave):
    bb_spec = bb_lambda(wave, params[0]).value
    bb_spec /= bb_spec.max()
    return (10**params[1])*bb_spec + params[2]


def res_bb(p, x_in, y_in, mask_in):
    #generate blackbody
    model_spec = bb_model(p, x_in)

    #return rms residual
    residual = (y_in - model_spec)[mask_in]
    return np.sqrt(bn.nanmean(residual**2))


def bb_lambda(wave, temp):
    #generate blackbody spectrum given temperature.  In case we ever need it,
    #put it in physical units (similar to the calibrated spectrum)
    #temp_bb = (blackbody_lambda(wave*u.micron, temp) * u.sr).to(u.watt / u.meter**2 / u.micron) / 1e-17
    temp_bb = (modeling.BlackBody(wave*u.micron, temp*u.K) * u.sr).to(u.watt / u.meter**2 / u.micron) / 1e-17
    return temp_bb

def salva(cube, filename):
    hdu0=fits.PrimaryHDU(data=cube)
    hdul=fits.HDUList([hdu0])
    hdul.writeto(filename,overwrite=True)
    


def zap_sky_v2(matched_exposures, logger=None, obsband=None, subdir=None, 
            inputdir=None, suffix='s', tbsub=False, **kwargs):

    for obj_exp, sky_exp, dit in matched_exposures:
        input_path = os.path.split(obj_exp)[0]
        output_path = os.path.join(inputdir, subdir, 'science')

        #fetch the superdark reference frame for this detector
        if obsband.upper() == 'YJ':
            if int(dit) == 600:
                logger.warning("- Using 600s SVD file")
                svd_file = 'ZAPv2_SVD_{0}_n9504_e600.fits'.format(obsband.upper())
            else:
                svd_file = 'ZAPv2_SVD_{0}_n500004.fits'.format(obsband.upper())
            #if int(dit) == 300:
            #    svd_file = 'ZAPv2_SVD_{0}_n500004.fits'.format(obsband.upper())
            #elif int(dit) == 600:
            #    logger.warning("- Using 600s SVD file")
            #    svd_file = 'ZAPv2_SVD_{0}_n9504_e600.fits'.format(obsband.upper())
        else:
            svd_file = 'ZAPv2_SVD_{0}_n500004.fits'.format(obsband.upper())

        fpath = fetch_resource_path('data/{0}'.format(svd_file))
        # copy object frame to output
        exp_sub = duplicate(obj_exp, outputdir=output_path, suffix=suffix) 

        if logger is not None:
            logger.info('- Subtracting sky for {0}'.format(obj_exp.split('/')[-1])) 
            logger.debug('-- Using {0} for initial subtraction'.format(sky_exp.split('/')[-1]))
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

                    #generate mask array for continuum
                    msky = find_lines(bn.nanmedian(scube[:,1:13,1:13].reshape(2048,-1), axis=1))
                    
                    #apply thermal background correction if flagged
                    if tbsub:
                        if logger is not None:
                            logger.info('--- Removing thermal background from object and sky frames')
                        #reconstruct wavelength solution
                        #print('ricostruisci')
                        wave = reconstruct_wave(header)

                        #extract 1D spectra from the un-sky-subtracted cubes
                        scube_corr = sub_thermal(wave, scube, msky)
                        ocube_corr = sub_thermal(wave, ocube, msky)
                        sub_cube = np.copy(ocube_corr) - np.copy(scube_corr)
                    else:
                        sub_cube = np.copy(ocube) - np.copy(scube)
                    
                    nanmask = np.isfinite(sub_cube)
                    nanmask2=copy.deepcopy(nanmask)
                    sub_cube[~nanmask] = 0.
                    
                    zap_obj = process_v2('dummy', cube=sub_cube, hdr=header, clean=False, 
                                        zlevel='none', cftype='med_smooth', cfwidthSP=101, 
                                        extSVD=fpath, interactive=True, overwrite=True,
                                        instrument='KMOS', boost=50, contmask=msky,
                                        reject_factor=0.4, logger=logger)

                    #reprocess the reconstructed sky spectrum to correct for any residual continuum drift
                    if logger is not None:
                        logger.debug('---- running second-order continuum cleaning')

                    cleaner_obj = ZapObj_v2('dummy', pca_class=None, n_components=None, instrument='KMOS',
                                            cube=zap_obj.recon.reshape(sub_cube.shape), hdr=header, logger=None)
                    cleaner_obj._prepare(clean=False, zlevel='none', cftype='med_smooth',
                                         cfwidth=101, extzlevel=None, mask=None, boost=50,
                                         contmask=msky, reject_factor=0.4)

                    #reincorporate residual continuum array
                    corr_cube = np.copy(zap_obj.cleancube+cleaner_obj.contarray.reshape(sub_cube.shape))
                    corr_cube[~nanmask2] = np.nan
                    ofile[ifu].data = np.asarray(corr_cube, dtype=np.float32)


def zap_sky_v1(matched_exposures, logger=None, obsband=None, subdir=None, 
            inputdir=None, suffix='s', tbsub=None, **kwargs):

    #fetch the superdark reference frame for this detector
    svd_file = 'ZAPv1_SVD_{0}_n500004.fits'.format(obsband.upper())
    fpath = fetch_resource_path('data/{0}'.format(svd_file))
    
    #if logger is not None:
    #    logger.debug("- Subtracting sky using {0}".format(svd_file))

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

                    #generate mask array for continuum
                    msky = find_lines(bn.nanmedian(scube[:,1:13,1:13].reshape(2048,-1), axis=1))
                    
                    #apply thermal background correction if flagged
                    if tbsub:
                        if logger is not None:
                            logger.info('--- Removing thermal background from object and sky frames')
                        #reconstruct wavelength solution
                        wave = reconstruct_wave(header)

                        #extract 1D spectra from the un-sky-subtracted cubes
                        scube_corr = sub_thermal(wave, scube, msky)
                        ocube_corr = sub_thermal(wave, ocube, msky)
                        sub_cube = np.copy(ocube_corr) - np.copy(scube_corr)
                    else:
                        sub_cube = np.copy(ocube) - np.copy(scube)
                    
                    nanmask = np.isfinite(sub_cube)
                    sub_cube[~nanmask] = 0.

                    zap_obj = process_v1('dummy', cube=sub_cube*exp_scale, hdr=header, clean=False, 
                                        zlevel='none', cftype='med_smooth', cfwidthSP=51, 
                                        extSVD=fpath, interactive=True, overwrite=True,
                                        instrument='KMOS', boost=50, contmask=msky,
                                        reject_factor=0.4, logger=logger)
                
                    corr_cube = np.copy(zap_obj.cleancube)
                    corr_cube[~nanmask] = np.nan

                    ofile[ifu].data = np.asarray(corr_cube, dtype=np.float32)


def subtract_sky(logger=None, subdir=None, obsband=None, 
                 sky_method='zapv2', search_suffix=None, 
                 tbsub=True, segment_sky=False, **kwargs):
    """ Driver function for sky subtraction.  Currently only supports
    subtraction using ZAP, but can be expanded to include other approaches.
    """
    #check that sky subtraction is handled by this routine
    if sky_method not in sky_methods:
        return 0

    #get list of frames to process. Assumes location of un-sky-subtracted frames in "sky" directory
    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'sky'), obsband=obsband, 
                                  suffix=search_suffix, **kwargs) 
    if not len(exposures):
        if logger is not None:
            logger.warning('-- No reconstructed frames found, skipping sky subtraction')
        return 1

    #generate associated list of objects and skies
    obj_sky = match_obj_sky(exposures)
    #print(obj_sky)
    #only apply thermal background correction for k-band
    if tbsub and obsband == 'K':
        tbsub = True
    else:
        tbsub = False


    if sky_method == 'zapv2':
        zap_sky_v2(obj_sky, logger=logger, obsband=obsband, subdir=subdir, 
                   tbsub=tbsub, segment_sky=segment_sky, **kwargs)    

    elif sky_method == 'zapv1':
        zap_sky_v1(obj_sky, logger=logger, obsband=obsband, subdir=subdir, tbsub=tbsub, **kwargs)    

    return 0
