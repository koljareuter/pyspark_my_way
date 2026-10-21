from __future__ import print_function
import numpy as np
from astropy.io import fits
import os
import glob

from ..utils.fileutils import (get_sof_list, clean_move, patch_header)
from ..utils.obsutils import (get_cal_rotator, run_spark)
from ..utils.setuputils import id_gen

__all__ = ["reconstruct"]

#static calibration directory based on SPARK_HOME environment variable
static_cal = os.environ['KMOS_CALIB']


def reconstruct(logger=None, subdir=None, obsband=None, sky_method='spark',
                combine=False, reconstruct_cubes=True, use_oh=True,
                telluric_in_sof=False, recipe_dir=None, lcorr_degree=None, **kwargs):

    """ Reconstruct science cubes using spark

    Optional [cubes] settings for the OH wavelength correction (docs/oh-lcorr-linear-rereduction.md):
      oh_spec_<band>  OH_SPEC file instead of $KMOS_CALIB/kmos_oh_spec_<band>.fits (e.g. oh_spec_iz)
      recipe_dir      esorex --recipe-dir, e.g. the KMOS DRP built with the KMOS_LCORR_DEGREE patch
      lcorr_degree    polynomial degree of the correction (DRP default 2); needs the patched recipe_dir
    """

    if not reconstruct_cubes:
        return 0

    #OH reference spectrum and polynomial degree of the OH wavelength correction
    oh_spec = kwargs.get('oh_spec_{0}'.format(obsband.lower()),
                         os.path.join(static_cal, 'kmos_oh_spec_{0}.fits'.format(obsband.lower())))
    if use_oh and not os.path.isfile(oh_spec):
        raise IOError('OH_SPEC file not found: {0}'.format(oh_spec))
    if recipe_dir and not os.path.isdir(recipe_dir):
        raise IOError('recipe_dir not found: {0}'.format(recipe_dir))
    if lcorr_degree is not None and not recipe_dir:
        raise ValueError('lcorr_degree needs recipe_dir (the patched KMOS DRP); the standard DRP ignores it')
    esorex_opts = '--recipe-dir={0}'.format(recipe_dir) if recipe_dir else ''
    lcorr_env = {'KMOS_LCORR_DEGREE': str(int(lcorr_degree))} if lcorr_degree is not None else None

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
            ofile.write('{0} OH_SPEC\n'.format(oh_spec))

        if logger is not None:
            logger.debug('-- Added {0} XCAL'.format(os.path.join(caldir, 'XCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} YCAL'.format(os.path.join(caldir, 'YCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} LCAL'.format(os.path.join(caldir, 'LCAL_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} MASTER_FLAT'.format(os.path.join(caldir, 'MASTER_FLAT_{0}{0}{0}.fits'.format(obsband.upper()))))
            logger.debug('-- Added {0} WAVE_BAND'.format(os.path.join(static_cal, 'kmos_wave_band.fits')))
            if use_oh:
                logger.debug('-- Added {0} OH_SPEC'.format(oh_spec))


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
    ecode = run_spark('kmos_sci_red', args=routine_args, sof=outfile, logger=logger,
                      esorex_opts=esorex_opts, env=lcorr_env)
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
                clean_move(frame_out_s, os.path.join(kwargs['outputdir'], subdir, 'science'))
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


