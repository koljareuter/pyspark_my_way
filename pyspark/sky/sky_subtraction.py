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


__all__ = ["subtract_sky", "match_obj_sky", "make_SVDs_v2", "build_custom_svd"]

sky_methods = ['zapv1', 'zapv2']



def match_obj_sky(exposures):
    """Given a list of exposures from get_obs_exposures, return
    (obj, sky_before, sky_after, weight_before, weight_after, dit) tuples.

    Each OBJ is bracketed by the closest DIT-matched SKY before and after it
    in time. Weights are linear in time so that
        w_before * sky_before + w_after * sky_after
    gives a time-interpolated sky frame. When only one bounding side exists,
    both slots point to the same file (w=1.0 on that side, 0.0 on the other).
    """
    time_dict       = {efile: etime for efile, _, etime, _ in exposures}
    object_exposures = [(efile, dit) for efile, etype, _, dit in exposures if etype == 'OBJ']
    sky_exposures    = [(efile, dit) for efile, etype, _, dit in exposures if etype == 'SKY']

    obj_sky_pair = []
    for obj, odit in object_exposures:
        t_obj   = time_dict[obj]
        matched = [(efile, time_dict[efile])
                   for efile, dit in sky_exposures if dit == odit]
        if not matched:
            continue

        before = [(f, t) for f, t in matched if t <= t_obj]
        after  = [(f, t) for f, t in matched if t >= t_obj]

        if before and after:
            sky_b, t_b = max(before, key=lambda x: x[1])  # latest before
            sky_a, t_a = min(after,  key=lambda x: x[1])  # earliest after
            if t_a == t_b:          # OBJ coincides with a SKY frame
                w_b, w_a = 1.0, 0.0
            else:
                dt  = t_a - t_b
                w_b = (t_a - t_obj) / dt
                w_a = 1.0 - w_b
        elif before:
            sky_b, _ = max(before, key=lambda x: x[1])
            sky_a, w_b, w_a = sky_b, 1.0, 0.0
        else:
            sky_a, _ = min(after, key=lambda x: x[1])
            sky_b, w_b, w_a = sky_a, 0.0, 1.0

        obj_sky_pair.append((obj, sky_b, sky_a, w_b, w_a, odit))
    return obj_sky_pair


def _oh_scale(science_cube, sky_cube, oh_mask):
    """Flux-scale sky_cube so its OH-line level matches science_cube.

    Corrects for OH variability faster than the bracketing cadence.
    Returns a scalar clamped to [0.5, 2.0].
    """
    if not oh_mask.any():
        return 1.0
    sci = bn.nanmedian(science_cube[oh_mask])
    sky = bn.nanmedian(sky_cube[oh_mask])
    if sky <= 0 or not np.isfinite(sky) or not np.isfinite(sci):
        return 1.0
    return float(np.clip(sci / sky, 0.5, 2.0))


def _make_object_mask(sub_cube, sigma_thresh=3.0):
    """Return a 2D bool mask marking galaxy spaxels in a sky-subtracted cube.

    Spaxels above background + sigma_thresh * noise are flagged so ZAP
    estimates its continuum from sky-only spaxels.
    """
    img  = bn.nanmedian(sub_cube, axis=0)
    flat = img[np.isfinite(img)].ravel()
    if flat.size < 9:
        return np.zeros(img.shape, dtype=bool)
    bkg   = np.nanpercentile(flat, 50)
    noise = 1.4826 * bn.nanmedian(np.abs(flat[flat <= bkg] - bkg))
    if noise <= 0 or not np.isfinite(noise):
        return np.zeros(img.shape, dtype=bool)
    return img > bkg + sigma_thresh * noise


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


def build_custom_svd(subdir=None, inputdir=None, obsband=None, pcaoutfits=None,
                     logger=None, **kwargs):
    """ Build a ZAPv2 eigenbasis (SVD) from this dataset's own SKY exposures,
    instead of using the generic bundled ZAPv2_SVD_*.fits file. Pools every
    spaxel's spectrum from every SKY-type exposure of `subdir`/`obsband`
    across all 24 IFUs, then hands that stack to make_SVDs_v2.

    Run this once per band before reduce_kmos.py, then point zap_sky_v2 at
    the result via svd_path=<pcaoutfits> (e.g. by setting svd_path in the
    [cubes] section of the .cfg -- subtract_sky forwards **kwargs straight
    through to zap_sky_v2).
    """
    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'sky'), obsband=obsband, inputdir=inputdir, **kwargs)
    sky_exposures = [efile for efile, etype, _, _ in exposures if etype == 'SKY']
    if not len(sky_exposures):
        if logger is not None:
            logger.warning('-- No SKY exposures found for {0}, cannot build custom SVD'.format(subdir))
        return None

    spectra, header = [], None
    for sky_exp in sky_exposures:
        with fits.open(sky_exp) as sfile:
            for ifu in range(1, 25):
                cube = np.array(sfile[ifu].data)
                if cube.ndim != 3:
                    continue
                if header is None:
                    header = sfile[ifu].header
                flat = cube.reshape(cube.shape[0], -1)
                # spaxels routinely have a few NaN samples near the edge of
                # their own valid-wavelength coverage (see CUBE_EDGE_MARGIN
                # in telluric_check.ipynb) -- requiring every pixel finite
                # throws away nearly everything. Keep spaxels that are
                # mostly finite and zero-fill the small remaining gaps.
                finite_frac = np.isfinite(flat).mean(axis=0)
                usable = flat[:, finite_frac >= 0.9]
                if usable.size:
                    spectra.append(np.where(np.isfinite(usable), usable, 0.0))

    if not spectra:
        if logger is not None:
            logger.warning('-- No finite sky spectra found for {0}, cannot build custom SVD'.format(subdir))
        return None

    spec_array = np.concatenate(spectra, axis=1)
    if logger is not None:
        logger.info('- Building custom SVD from {0} sky spaxel spectra across {1} exposure(s)'.format(
            spec_array.shape[1], len(sky_exposures)))

    make_SVDs_v2(spec_array, header, pcaoutfits=pcaoutfits, logger=logger)
    return pcaoutfits


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
    bb = models.BlackBody(temperature=temp*u.K)
    temp_bb = (bb(wave*u.micron) * u.sr).to(u.watt / u.meter**2 / u.micron, equivalencies=u.spectral_density(wave*u.micron)) / 1e-17
    return temp_bb

def salva(cube, filename):
    hdu0=fits.PrimaryHDU(data=cube)
    hdul=fits.HDUList([hdu0])
    hdul.writeto(filename,overwrite=True)
    


def zap_sky_v2(matched_exposures, logger=None, obsband=None, subdir=None,
            inputdir=None, suffix='s', tbsub=False, svd_path=None, **kwargs):

    for obj_exp, sky_b, sky_a, w_b, w_a, dit in matched_exposures:
        input_path = os.path.split(obj_exp)[0]
        output_path = os.path.join(inputdir, subdir, 'science')

        if svd_path is not None:
            # use a dataset-specific eigenbasis (see build_custom_svd below)
            # instead of the generic bundled one -- its OH-line residual
            # structure matches this run's own sky frames, not whatever
            # dataset ZAPv2_SVD_*_n500004.fits was originally built from.
            if logger is not None:
                logger.info('- Using custom SVD file {0}'.format(svd_path))
            fpath = svd_path
        else:
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

        bracketed = (sky_b != sky_a)
        if logger is not None:
            logger.info('- Subtracting sky for {0}'.format(obj_exp.split('/')[-1]))
            if bracketed:
                logger.debug('-- Bracketed sky: {0} (w={1:.2f}) + {2} (w={3:.2f})'.format(
                    sky_b.split('/')[-1], w_b, sky_a.split('/')[-1], w_a))
            else:
                logger.debug('-- Using {0} for initial subtraction'.format(sky_b.split('/')[-1]))
        with fits.open(exp_sub, mode='update') as ofile, \
             fits.open(sky_b) as sbfile, fits.open(sky_a) as safile:

            for ifu in range(1,25):
                scube_b = np.array(sbfile[ifu].data)
                scube_a = np.array(safile[ifu].data)
                ocube = np.array(ofile[ifu].data)
                header = ofile[ifu].header

                try:
                    oz, ox, oy = ocube.shape
                except ValueError:
                    continue
                else:
                    if logger is not None:
                        logger.info('-- Processing sky subtraction for IFU {0}'.format(ifu))

                    # Build interpolated sky; skip or fall back for dead sky IFUs
                    b_ok = np.ndim(scube_b) == 3
                    a_ok = np.ndim(scube_a) == 3
                    if not b_ok and not a_ok:
                        continue

                    # msky needed before interpolation so we can OH-scale each frame
                    _ref_for_msky = scube_b if b_ok else scube_a
                    msky = find_lines(bn.nanmedian(_ref_for_msky[:,1:13,1:13].reshape(2048,-1), axis=1))

                    # Scale each sky frame to match the science OH-line flux level before
                    # time-interpolation, correcting for variability faster than the cadence
                    if b_ok and a_ok:
                        scale_b = _oh_scale(ocube, scube_b, msky)
                        scale_a = _oh_scale(ocube, scube_a, msky)
                        if logger is not None:
                            logger.debug('--- OH scale: sky_b={0:.3f}  sky_a={1:.3f}'.format(scale_b, scale_a))
                        scube   = w_b * (scale_b * scube_b) + w_a * (scale_a * scube_a)
                        ref_sky = scube_b if w_b >= w_a else scube_a
                    elif b_ok:
                        scale_b = _oh_scale(ocube, scube_b, msky)
                        scube   = scale_b * scube_b
                        ref_sky = scube_b
                    else:
                        scale_a = _oh_scale(ocube, scube_a, msky)
                        scube   = scale_a * scube_a
                        ref_sky = scube_a

                    #apply thermal background correction if flagged
                    if tbsub:
                        if logger is not None:
                            logger.info('--- Removing thermal background from object and sky frames')
                        wave = reconstruct_wave(header)
                        scube_corr = sub_thermal(wave, scube, msky)
                        ocube_corr = sub_thermal(wave, ocube, msky)
                        sub_cube = np.copy(ocube_corr) - np.copy(scube_corr)
                    else:
                        sub_cube = np.copy(ocube) - np.copy(scube)

                    nanmask  = np.isfinite(sub_cube)
                    nanmask2 = copy.deepcopy(nanmask)
                    sub_cube[~nanmask] = 0.

                    # Mask bright galaxy spaxels so ZAP estimates sky from sky-only pixels
                    obj_mask_2d = _make_object_mask(sub_cube)
                    sub_cube_for_zap = np.copy(sub_cube)
                    if obj_mask_2d.any():
                        sky_only  = sub_cube[:, ~obj_mask_2d]
                        bkg_fill  = float(bn.nanmedian(sky_only)) if sky_only.size else 0.0
                        sub_cube_for_zap[:, obj_mask_2d] = bkg_fill if np.isfinite(bkg_fill) else 0.0
                        if logger is not None:
                            logger.debug('--- Masked {0} galaxy spaxels from ZAP'.format(obj_mask_2d.sum()))

                    zap_obj = process_v2('dummy', cube=sub_cube_for_zap, hdr=header, clean=True,
                                        zlevel='none', cftype='med_smooth', cfwidthSP=55,
                                        extSVD=fpath, interactive=True, overwrite=True,
                                        instrument='KMOS', boost=25, contmask=msky,
                                        reject_factor=0.4, logger=logger)

                    # Apply sky model (derived from masked cube) to the original unmasked cube
                    sky_model    = zap_obj.recon.reshape(sub_cube.shape)
                    sub_cube_corr = sub_cube - sky_model

                    #reprocess the reconstructed sky to correct for residual continuum drift
                    if logger is not None:
                        logger.debug('---- running second-order continuum cleaning')
                    cleaner_obj = ZapObj_v2('dummy', pca_class=None, n_components=None,
                                            instrument='KMOS', cube=sky_model,
                                            hdr=header, logger=None)
                    cleaner_obj._prepare(clean=True, zlevel='none', cftype='med_smooth',
                                         cfwidth=55, extzlevel=None, mask=None, boost=25,
                                         contmask=msky, reject_factor=0.4)

                    #reincorporate residual continuum array
                    corr_cube = np.copy(sub_cube_corr + cleaner_obj.contarray.reshape(sub_cube.shape))
                    corr_cube[~nanmask2] = np.nan
                    ofile[ifu].data = np.asarray(corr_cube, dtype=np.float32)


def zap_sky_v1(matched_exposures, logger=None, obsband=None, subdir=None, 
            inputdir=None, suffix='s', tbsub=None, **kwargs):

    #fetch the superdark reference frame for this detector
    svd_file = 'ZAPv1_SVD_{0}_n500004.fits'.format(obsband.upper())
    fpath = fetch_resource_path('data/{0}'.format(svd_file))
    
    #if logger is not None:
    #    logger.debug("- Subtracting sky using {0}".format(svd_file))

    for obj_exp, sky_b, sky_a, w_b, w_a, dit in matched_exposures:
        # v1 does not support time-interpolated bracketing; use the nearest sky
        sky_exp = sky_b if w_b >= w_a else sky_a
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
                        wave = reconstruct_wave(header)

                        scube_corr = sub_thermal(wave, scube, msky)
                        ocube_corr = sub_thermal(wave, ocube, msky)
                        sub_cube = np.copy(ocube_corr) - np.copy(scube_corr)
                    else:
                        sub_cube = np.copy(ocube) - np.copy(scube)

                    nanmask = np.isfinite(sub_cube)
                    sub_cube[~nanmask] = 0.

                    zap_obj = process_v1('dummy', cube=sub_cube, hdr=header, clean=False,
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
