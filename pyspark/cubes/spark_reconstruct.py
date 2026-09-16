from __future__ import print_function
import numpy as np
from astropy.io import fits
import os
import glob

from ..utils.fileutils import (get_sof_list, clean_move, patch_header)
from ..utils.obsutils import (get_cal_rotator, run_spark, reconstruct_wave)
from ..utils.setuputils import id_gen

__all__ = ["reconstruct"]


# Deep-absorption wavelength windows used to verify that the telluric
# correction was actually applied (kmos_sci_red NaN-masks T < 0.1 there).
_TELL_CHECK_BANDS = {
    'YJ': (1.34, 1.42),
    'H':  (1.79, 1.94),
    'K':  (1.95, 2.05),
}


def _verify_telluric(cube_path, obsband, logger):
    """Warn if the deep-absorption region contains no NaN channels.

    After telluric correction kmos_sci_red sets T<0.1 channels to NaN.
    If the median NaN fraction across IFUs is below 5% the correction
    likely was not applied, which warrants a loud warning.
    """
    lo, hi = _TELL_CHECK_BANDS.get(obsband.upper(), (None, None))
    if lo is None:
        return
    try:
        with fits.open(cube_path) as hdul:
            nan_fracs = []
            for ifu in range(1, 25):
                d = hdul[ifu].data
                if d is None or np.ndim(d) != 3:
                    continue
                wave = reconstruct_wave(hdul[ifu].header)
                in_band = (wave >= lo) & (wave <= hi)
                if not in_band.any():
                    continue
                nan_fracs.append(np.isnan(d[in_band]).mean())
            if not nan_fracs:
                return
            med = float(np.median(nan_fracs))
            if med < 0.05:
                if logger is not None:
                    logger.warning(
                        '!! Telluric check {0}: only {1:.1%} NaN in deep-band '
                        '{2:.2f}-{3:.2f} µm — telluric correction may not have '
                        'been applied'.format(
                            os.path.basename(cube_path), med, lo, hi))
            else:
                if logger is not None:
                    logger.info(
                        '-- Telluric check OK ({0}): {1:.1%} NaN in deep-band '
                        'across {2} IFUs'.format(
                            os.path.basename(cube_path), med, len(nan_fracs)))
    except Exception:
        pass

#static calibration directory based on SPARK_HOME environment variable
static_cal = os.environ['KMOS_CALIB']


def reconstruct(logger=None, subdir=None, obsband=None, sky_method='spark',
                combine=False, reconstruct_cubes=True, use_oh=True,
                telluric_in_sof=False, **kwargs):

    """ Reconstruct science cubes using spark
    """

    if not reconstruct_cubes:
        return 0

    #path for the calibration directory
    caldir = os.path.join(kwargs['inputdir'], subdir, 'calib')
    
    #get list of science frames for processing
    frames = get_sof_list(filetype='science', subdir=subdir, obsband=obsband, **kwargs)

    #write temporary .sof file to feed to
    outfile = '{0}.sof'.format(id_gen())
    
    if logger is not None:
        logger.debug('- Writing science frames to temporary file {0}'.format(outfile))

    with open(outfile, 'wt') as ofile:
        for frame, _ in frames:
            ofile.write('{0} SCIENCE\n'.format(frame))
            if logger is not None:
                logger.debug('-- Added {0} SCIENCE'.format(frame))

        #include additional files required for SPARK
        ofile.write('{0} XCAL\n'.format(os.path.join(caldir, 'XCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
        ofile.write('{0} YCAL\n'.format(os.path.join(caldir, 'YCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
        ofile.write('{0} LCAL\n'.format(os.path.join(caldir, 'LCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
        ofile.write('{0} MASTER_FLAT\n'.format(os.path.join(caldir, 'MASTER_FLAT_{0}{0}{0}.fits'.format(obsband.upper()))))
        ofile.write('{0} WAVE_BAND\n'.format(os.path.join(static_cal, 'kmos_wave_band.fits')))
        if use_oh:
            ofile.write('{0} OH_SPEC\n'.format(os.path.join(static_cal, 'kmos_oh_spec_{0}.fits'.format(obsband.lower()))))
        if telluric_in_sof:
            # TELLURIC files are raw data — look in rawdir, not the procdir
            # that inputdir points to for the cubes stage
            tel_dir = kwargs.get('rawdir', kwargs['inputdir'])
            tel_matches = glob.glob(os.path.join(tel_dir, subdir, '*TELLURIC*{0}*'.format(obsband.upper())))
            if tel_matches:
                ofile.write('{0} TELLURIC\n'.format(tel_matches[0]))
                if logger is not None:
                    logger.debug('-- Added {0} TELLURIC'.format(tel_matches[0]))
            elif logger is not None:
                logger.warning('-- telluric_in_sof=True but no TELLURIC file found in {0}'.format(
                    os.path.join(tel_dir, subdir)))

        if logger is not None:
            logger.debug('-- Added {0} XCAL'.format(os.path.join(caldir, 'XCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} YCAL'.format(os.path.join(caldir, 'YCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} LCAL'.format(os.path.join(caldir, 'LCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} MASTER_FLAT'.format(os.path.join(caldir, 'MASTER_FLAT_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} WAVE_BAND'.format(os.path.join(static_cal, 'kmos_wave_band.fits')))
            if use_oh:
                logger.debug('-- Added {0} OH_SPEC'.format(os.path.join(static_cal, 'kmos_oh_spec_{0}.fits'.format(obsband.lower()))))


    #control args
    routine_args = '-background=FALSE -tbsub=TRUE '
    
    #sky method
    if sky_method == 'spark':
        routine_args += '-no_subtract=FALSE -sky_tweak=TRUE '
    elif sky_method == 'edps':    
         routine_args += '-no_subtract=FALSE -sky_tweak=TRUE -flux=TRUE -lcmethod=NONE -method=header -collapse_reconstructed=TRUE '
    else:
        #routine_args += '-no_subtract=TRUE -sky_tweak=FALSE -lcmethod=SLICES_MEDIAN '
        routine_args += '-no_subtract=TRUE -sky_tweak=FALSE -lcmethod=SLICES_MEDIAN -skip_sky_oh_align=TRUE '

    #never combine at this stage
    routine_args += '-no_combine=TRUE '

    #execute spark script
    if logger is not None:
        logger.debug('- Executing reconstruction')
    ecode = run_spark('kmos_sci_red', args=routine_args, sof=outfile, logger=logger)
    if logger is not None:
        logger.debug('-- Complete')

    #check the exit status of the SPARK call
    if ecode != 0:
        return 1

    #move output to their final locations
    hdr_keys = ['DATE', 'DATE-OBS', 'OBJECT', 'RA', 'DEC', 'MJD-OBS', 'UTC', 'LST']
    for frame, _ in frames:
        frame_out = 'SCI_RECONSTRUCTED_{0}'.format(frame.split('/')[-1])
        frame_out2 = 'make_image_SCI_RECONSTRUCTED_{0}'.format(frame.split('/')[-1])
        if os.path.isfile(frame_out):
            patch_header(frame_out, hdr_keys, frame)
            if sky_method == 'spark':
                # kmos_sci_red already sky-subtracted; rename to *_s.fits and
                # send directly to science/ so illuminate() can find it
                frame_out_s = frame_out.replace('.fits', '_s.fits')
                os.rename(frame_out, frame_out_s)
                dest = os.path.join(kwargs['outputdir'], subdir, 'science')
                clean_move(frame_out_s, dest)
                if telluric_in_sof and tel_matches:
                    _verify_telluric(
                        os.path.join(dest, os.path.basename(frame_out_s)),
                        obsband, logger)
            else:
                # ZAP/zapv1 will read from sky/ and write *_s.fits to science/
                clean_move(frame_out, os.path.join(kwargs['outputdir'], subdir, 'sky'))
        if os.path.isfile(frame_out2):
            clean_move(frame_out2, os.path.join(kwargs['outputdir'], subdir, 'sky'))

    #move the .sof file for logging
    clean_move(outfile, os.path.join(kwargs['outputdir'], subdir, 'sof', '{0}_cubes_{1}.sof'\
            .format(subdir, obsband.lower())))
   
    #cleanup
    [os.remove(x) for x in glob.glob('obj_sky_table*.txt')]

    return 0


