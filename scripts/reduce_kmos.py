#!/usr/local/bin/python
from __future__ import print_function
import numpy as np
import os, sys
import glob
import argparse

from pyspark.utils import (fileutils, setuputils, obsutils)
from pyspark.detector import fix_detector_levels
from pyspark.calibration import make_spark_calibs
from pyspark.cubes import (reconstruct, velocity_correct, illuminate)
from pyspark.sky import subtract_sky


def break_check(routine, logger):
    logger.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
    logger.error('!!!! An error occured during {0}, check the log'.format(routine))
    logger.error('!!!! Moving on to next night')
    logger.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')

    #clean up any dangling sof files
    [os.remove(x) for x in glob.glob('*.sof')]
    


#handle parsing of command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--parameter_file', type=str, help='paramter file controlling reduction')
parser.add_argument('--verbose', type=int, default=3, help='Logging level for stdout')
args = parser.parse_args()

param_file = args.parameter_file
verbosity = args.verbose


# --- initialize log for this reduction session ---
master_logger = setuputils.set_logging(logbase='master', log_level=verbosity)


# --- read input configuration file ------------------
run_params = setuputils.set_config(param_file, logger=master_logger)


#pull out information of interest for iteration
rawdir = run_params['rawdir']
procdir = run_params['procdir']
folders = run_params['subdirs']
bands = run_params['obsband']


# --- start iterating over bands ------------------
for band in bands:
    master_logger.info('######## Processing {0}-band observations ########'.format(band))

    #iterate over folders
    for folder in folders:

        # before anything, check if there are science files in this
        # folder observed with the given band
        iframes = fileutils.get_sof_list(filetype='science', subdir=folder,
                inputdir=rawdir, obsband=band)

        if not len(iframes):
            master_logger.info("No science exposures for this night+band, moving on")
            continue

        master_logger.info('Reducing {0} '.format(folder))

        #make directory structure in case it doesn't exist
        fileutils.verify_structure(outputdir=procdir, subdir=folder)

        #initialize night/band specific log
        dir_logger = setuputils.set_logging(logbase='{0}'.format(band), 
                outputdir=os.path.join(procdir, folder), log_level=verbosity, 
                master_log=master_logger)
       

        # -- detector corrections ----------------
        master_logger.info('- Running detector corrections')
        
        #get list of obs frames
        frames = fileutils.get_sof_list(filetype='science', subdir=folder, 
                inputdir=rawdir, obsband=band)

        #apply detector fix(es)
        if run_params['detector']['channel_correct'] or run_params['detector']['picture_correct']:
            try:
                fix_detector_levels([frame for frame, type in frames], subdir=folder, 
                        logger=dir_logger, **run_params['detector'])
            except:
                break_check('Det. correction', master_logger)
                continue
            else:
                master_logger.info('-- Finished detector corrections')
                master_logger.info(' ')
        else:
            master_logger.info('-- Finished detector corrections')
            master_logger.info(' ')


        # -- calibrations -----------------------
        master_logger.info('- Running calibrations')
        ecode, comp = make_spark_calibs(logger=dir_logger, subdir=folder, obsband=band, 
                          **run_params['calibration'])
        if ecode:
            break_check(comp, master_logger)
            continue
        master_logger.info('-- Finished calibrations')
        master_logger.info(' ')


        # -- reconstruct -----------------------
        master_logger.info('- Running reconstruction')
        ecode = reconstruct(logger=dir_logger, subdir=folder, obsband=band, 
                            **run_params['cubes'])
        if ecode:
            break_check(comp, master_logger)
            continue
        master_logger.info('-- Finished reconstruction')
        master_logger.info(' ')


        # -- sky subtraction -------------------
        master_logger.info('- Running sky subtraction')
        ecode = subtract_sky(logger=dir_logger, subdir=folder, obsband=band, 
                             **run_params['cubes'])
        if ecode:
            break_check('sky_subtraction', master_logger)
            continue
        master_logger.info('-- Finished sky subtraction')
        master_logger.info(' ')


        # -- illumination correction -----------
        master_logger.info('- Applying illumination correction')
        ecode = illuminate(logger=dir_logger, subdir=folder, obsband=band, 
                           search_suffix='s', **run_params['cubes'])
        if ecode:
            break_check('illumination correction', master_logger)
            continue
        master_logger.info('-- Finished illumination correction')
        master_logger.info(' ')


        # -- barycentric velocity correction -----
        master_logger.info('- Applying heliocentric correction')
        ecode = velocity_correct(logger=dir_logger, subdir=folder, obsband=band, 
                                 search_suffix='s_i', **run_params['cubes'])
        if ecode:
            break_check('velocity correction', master_logger)
            continue
        master_logger.info('-- Finished heliocentric correction')
        master_logger.info(' ')
 
        
        master_logger.info('- Finished reducing {0} '.format(folder))
