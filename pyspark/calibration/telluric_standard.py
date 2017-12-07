from __future__ import print_function
import numpy as np
from astropy.io import fits
from collections import Counter
import os

from ..utils.fileutils import (get_sof_list, clean_move)
from ..utils.obsutils import (get_cal_rotator, run_spark, match_obs_rotator)
from ..utils.setuputils import id_gen


__all__ = ["spark_standard"]

#static calibration directory based on SPARK_HOME environment variable
static_cal = os.environ['KMOS_CALIB']

# I've currently left this as a separate set of routines as its not yet
# clear if there's additional processing we'd like to do (e.g. update
# header information with a simbad query or extract our own spectra.


def spark_standard(logger=None, subdir=None, obsband=None, **kwargs):

    """ Reconstruct all standard stars for this band.
    """

    #path for the output calibration directory
    caldir = os.path.join(kwargs['outputdir'], subdir, 'calib')
    stddir = os.path.join(kwargs['outputdir'], subdir, 'std')

    #fetch rotator angles and filled extensions from the master flat file
    flat_file = fits.open(os.path.join(caldir, 'master_flat_{0}{0}{0}.fits'.format(obsband)))
    pro_naangles = get_cal_rotator(flat_file)
    nrot = len(pro_naangles)
    flat_file.close()


    #get list of flat frames for processing
    frames = get_sof_list(filetype='standard', subdir=subdir, obsband=obsband, **kwargs)

    #group frames by start time (these are unique OBs)
    group_dict = {}
    name_dict = {}
    dit_dict = {}
    rot_dict = {}
    for obs_file, obs_dit, obs_start, obs_name, obs_rot in frames:
        group_dict.setdefault(obs_start, []).append(obs_file)
        name_dict.setdefault(obs_name, []).append(obs_start)
        dit_dict[obs_start] = obs_dit #dits are always the same within the same OB
    
        if obs_rot < 0: obs_rot += 360
        rot_dict.setdefault(obs_start, []).append(obs_rot)
    
    name_dict = dict([(star, list(np.unique(start_list))) for star, start_list in name_dict.items()])
    rot_dict = dict([(ttime, np.mean(rot_list)) for ttime, rot_list in rot_dict.items()])

    # pass over each star, being careful to handle multiple observations
    for star, times in name_dict.items():
        dits = [dit_dict[ot] for ot in times]

        cnt_exp = 1
        for ttime in times:

            if len(list(np.unique(dits))) != len(times): #catch duplicate observations with same dit
                base_out = '{0}_{1:.1f}s_{2}_{3}'.format(star, float(dit_dict[times[0]]), 
                                                              obsband.lower(), cnt_exp) 
                cnt_exp += 1
            else:
                base_out = '{0}_{1:.1f}s_{2}'.format(star, float(dit_dict[times[0]]), obsband.lower())
      
            # get matched rotation for the illumination correctoin
            _, matched_rot, this_star_rot = match_obs_rotator('dummy', pro_naangles, 
                                                              input_rot=float(rot_dict[ttime])) 

            #write temporary .sof file to feed to kmos
            outfile = '{0}.sof'.format(id_gen())
            logger.debug('- Writing standard star frames to temporary file {0}'.format(outfile))

            with open(outfile, 'wt') as ofile:
                for frame in group_dict[ttime]:
                    ofile.write('{0} STD\n'.format(frame))
                    logger.debug('-- Added {0} STD'.format(frame))

                ofile.write('{0}/illum_corr_{1}{1}{1}_{2}.fits ILLUM_CORR\n'.format(caldir, obsband.upper(), matched_rot))
                logger.debug('-- Added {0}/illum_corr_{1}{1}{1}_{2}.fits ILLUM_CORR'.format(caldir, obsband.upper(), matched_rot))
                ofile.write('{0}/master_flat_{1}{1}{1}.fits MASTER_FLAT\n'.format(caldir, obsband.upper()))
                logger.debug('-- Added {0}/master_flat_{1}{1}{1}.fits MASTER_FLAT'.format(caldir, obsband.upper()))
                ofile.write('{0}/xcal_{1}{1}{1}.fits XCAL\n'.format(caldir, obsband.upper()))
                logger.debug('-- Added {0}/xcal_{1}{1}{1}.fits XCAL'.format(caldir, obsband.upper()))
                ofile.write('{0}/ycal_{1}{1}{1}.fits YCAL\n'.format(caldir, obsband.upper()))
                logger.debug('-- Added {0}/ycal_{1}{1}{1}.fits YCAL'.format(caldir, obsband.upper()))
                ofile.write('{0}/lcal_{1}{1}{1}.fits LCAL\n'.format(caldir, obsband.upper()))
                logger.debug('-- Added {0}/lcal_{1}{1}{1}.fits LCAL'.format(caldir, obsband.upper()))
                ofile.write('{0}/kmos_wave_band.fits WAVE_BAND\n'.format(static_cal))
                logger.debug('-- Added {0}/kmos_wave_band.fits WAVE_BAND'.format(static_cal))
                ofile.write('{0}/kmos_spec_type.fits SPEC_TYPE_LOOKUP\n'.format(static_cal))
                logger.debug('-- Added {0}/kmos_spec_type.fits SPEC_TYPE_LOOKUP'.format(static_cal))
                ofile.write('{0}/kmos_atmos_{1}.fits ATMOS_MODEL\n'.format(static_cal, obsband.lower()))
                logger.debug('-- Added {0}/kmos_atmos_{1}.fits ATMOS_MODEL'.format(static_cal, obsband.lower()))

       
            #control args
            routine_args = '--no_noise --save_cubes=TRUE'

            #execute spark script
            logger.debug('- Executing kmo_std_star')
            ecode = run_spark('kmo_std_star', args=routine_args, sof=outfile, logger=logger)
            logger.debug('-- Complete')
        
            #check the exit status of the SPARK call
            if ecode != 0:
                return 1

            #move files to their final location
            products = ['telluric_{0}{0}{0}{1}.fits', 'star_spec_{0}{0}{0}{1}.fits', 'std_mask_{0}{0}{0}{1}.fits',
                        'std_image_{0}{0}{0}{1}.fits', 'std_cube_{0}{0}{0}{1}.fits']
            for product in products:
                clean_move(product.format(obsband.upper(),''), 
                           os.path.join(stddir, product.format('', base_out)))
        
            #move the .sof file for logging
            clean_move(outfile, os.path.join(kwargs['outputdir'], subdir, 'sof', '{0}_stdstar_{1}.sof'\
                    .format(subdir, base_out)))
    return 0
