#!/home/lcoccato/pyspark/env/bin/python
from __future__ import print_function
import numpy as np
import os, sys
import glob
import argparse
import json
import tempfile
import shutil
import concurrent.futures

from pyspark.utils import (fileutils, setuputils, obsutils)
from pyspark.detector import fix_detector_levels
from pyspark.calibration import make_spark_calibs
from pyspark.cubes import (reconstruct, velocity_correct, illuminate, bkg_correct_2d)
from pyspark.sky import subtract_sky

from multiprocessing import freeze_support
freeze_support()

_CHECKPOINT_FILE = '.pyspark_checkpoint.json'


def _load_checkpoint(procdir):
    """Return the set of (band, folder) pairs already successfully reduced."""
    path = os.path.join(procdir, _CHECKPOINT_FILE)
    if os.path.isfile(path):
        try:
            with open(path) as f:
                return set(tuple(x) for x in json.load(f))
        except Exception:
            pass
    return set()


def _save_checkpoint(procdir, band, folder, done_set):
    """Record a successfully completed (band, folder) pair."""
    done_set.add((band, folder))
    path = os.path.join(procdir, _CHECKPOINT_FILE)
    with open(path, 'w') as f:
        json.dump(sorted(done_set), f, indent=2)


def reduce_one_ob(band, folder, run_params, verbosity):
    """Reduce a single (band, OB-folder) combination.

    Self-contained so it can be submitted to a ProcessPoolExecutor.
    EsoRex calls run inside a per-process temp directory, so parallel
    workers never clobber each other's output files.

    Returns (band, folder, status) where status is 'done', 'skipped',
    or 'failed:<step>'.
    """
    rawdir  = run_params['rawdir']
    procdir = run_params['procdir']

    iframes = fileutils.get_sof_list(filetype='science', subdir=folder,
                                     inputdir=rawdir, obsband=band)
    if not iframes:
        return band, folder, 'skipped'

    fileutils.verify_structure(outputdir=procdir, subdir=folder)
    dir_logger = setuputils.set_logging(
        logbase=band,
        outputdir=os.path.join(procdir, folder),
        log_level=verbosity)

    orig_cwd = os.getcwd()
    tmpdir   = tempfile.mkdtemp(prefix='pyspark_{0}_{1}_'.format(
        band, os.path.basename(folder)[:20]))
    try:
        # All EsoRex calls write output to CWD; isolate each worker in a
        # private temp directory so parallel runs never clobber each other.
        os.chdir(tmpdir)

        ecode, comp = make_spark_calibs(logger=dir_logger, subdir=folder,
                                        obsband=band, **run_params['calibration'])
        if ecode:
            return band, folder, 'failed:calibration'

        frames = fileutils.get_sof_list(filetype='science', subdir=folder,
                                        inputdir=rawdir, obsband=band)
        det = run_params['detector']
        if det['channel_correct'] or det['picture_correct'] or det['interchannel_correct']:
            try:
                fix_detector_levels([f for f, _ in frames], subdir=folder,
                                    logger=dir_logger, **det)
            except Exception:
                dir_logger.exception('Detector correction failed')
                return band, folder, 'failed:detector'

        ecode = reconstruct(logger=dir_logger, subdir=folder, obsband=band,
                            **run_params['cubes'])
        if ecode:
            return band, folder, 'failed:reconstruct'

        try:
            ecode = subtract_sky(logger=dir_logger, subdir=folder, obsband=band,
                                 **run_params['cubes'])
        except Exception:
            dir_logger.exception('Sky subtraction raised an exception')
            return band, folder, 'failed:sky'
        if ecode:
            return band, folder, 'failed:sky'

    finally:
        os.chdir(orig_cwd)
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Python-only steps: use absolute paths, no EsoRex, safe to run in any CWD
    ecode = illuminate(logger=dir_logger, subdir=folder, obsband=band,
                       search_suffix='s', **run_params['cubes'])
    if ecode:
        return band, folder, 'failed:illuminate'

    ecode = bkg_correct_2d(logger=dir_logger, subdir=folder, obsband=band,
                           search_suffix='s_i', **run_params['cubes'])
    if ecode:
        return band, folder, 'failed:bkgcorr'

    ecode = velocity_correct(logger=dir_logger, subdir=folder, obsband=band,
                             search_suffix='s_i_bkg2d', **run_params['cubes'])
    if ecode:
        return band, folder, 'failed:velocity'

    return band, folder, 'done'

def break_check(routine, logger):
    logger.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
    logger.error('!!!! An error occured during {0}, check the log'.format(routine))
    logger.error('!!!! Moving on to next night')
    logger.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')

    #clean up any dangling sof files
    [os.remove(x) for x in glob.glob('*.sof')]
    

if __name__ == "__main__":
 #handle parsing of command line arguments
 freeze_support()
 parser = argparse.ArgumentParser()
 parser.add_argument('-p', '--parameter_file', type=str, help='parameter file controlling reduction')
 parser.add_argument('--verbose', type=int, default=3, help='Logging level for stdout')
 parser.add_argument('--jobs', type=int, default=1,
                     help='Number of parallel OB workers (default 1 = serial). '
                          'Each worker runs EsoRex in its own temp directory.')
 args = parser.parse_args()

 param_file = args.parameter_file
 verbosity  = args.verbose
 n_jobs     = max(1, args.jobs)

 # --- initialize log for this reduction session ---
 #master_logger = setuputils.set_logging(logbase='master', log_level=verbosity)
 master_logger = setuputils.set_logging(logbase='master', log_level=3)


 # --- read input configuration file ------------------
 run_params = setuputils.set_config(param_file, logger=master_logger)


 #pull out information of interest for iteration
 rawdir = run_params['rawdir']
 procdir = run_params['procdir']
 folders = run_params['subdirs']
 bands = run_params['obsband']

 # --- checkpoint: skip already-completed (band, folder) pairs ---
 done_set = _load_checkpoint(procdir)
 if done_set:
     master_logger.info('Checkpoint: {0} (band, OB) pair(s) already done, will skip.'.format(
         len(done_set)))

 # --- parallel path ---
 if n_jobs > 1:
     tasks = [(band, folder) for band in bands for folder in folders
              if (band, folder) not in done_set]
     master_logger.info('Parallel mode: {0} task(s) across {1} worker(s)'.format(
         len(tasks), n_jobs))
     with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as pool:
         futs = {pool.submit(reduce_one_ob, b, f, run_params, verbosity): (b, f)
                 for b, f in tasks}
         for fut in concurrent.futures.as_completed(futs):
             try:
                 band_r, folder_r, status = fut.result()
             except Exception as exc:
                 band_r, folder_r = futs[fut]
                 master_logger.error('Worker raised exception for {0}/{1}: {2}'.format(
                     band_r, folder_r, exc))
                 continue
             if status == 'done':
                 master_logger.info('Finished {0}/{1}'.format(band_r, folder_r))
                 _save_checkpoint(procdir, band_r, folder_r, done_set)
             elif status == 'skipped':
                 master_logger.info('No {0}-band science for {1}'.format(band_r, folder_r))
             else:
                 master_logger.error('Failed {0}/{1}: {2}'.format(band_r, folder_r, status))
     sys.exit(0)

 # --- serial path (original) ---
 # --- start iterating over bands ------------------
 for band in bands:
    master_logger.info('######## Processing {0}-band observations ########'.format(band))

    #iterate over folders
    for folder in folders:

        if (band, folder) in done_set:
            master_logger.info('Skipping {0}/{1} (checkpoint)'.format(band, folder))
            continue

        # before anything, check if there are science files in this
        # folder observed with the given band
        iframes = fileutils.get_sof_list(filetype='science', subdir=folder,
                inputdir=rawdir, obsband=band)

        if not len(iframes):
            master_logger.info("No {0}-band science exposures for {1}, moving on".format(band, folder))
            continue

        master_logger.info('Reducing {0} '.format(folder))

        #make directory structure in case it doesn't exist
        fileutils.verify_structure(outputdir=procdir, subdir=folder)

        #initialize night/band specific log
        dir_logger = setuputils.set_logging(logbase='{0}'.format(band), 
                outputdir=os.path.join(procdir, folder), log_level=verbosity, 
                master_log=master_logger)
       

        # -- calibrations -----------------------
        master_logger.info('- Running calibrations')
        ecode, comp = make_spark_calibs(logger=dir_logger, subdir=folder, obsband=band, 
                          **run_params['calibration'])
        if ecode:
            break_check(comp, master_logger)
            continue
        master_logger.info('-- Finished calibrations')
        master_logger.info(' ')
   
   
        # -- detector corrections ----------------
        master_logger.info('- Running detector corrections')
        
        #get list of obs frames
        frames = fileutils.get_sof_list(filetype='science', subdir=folder, 
                inputdir=rawdir, obsband=band)
   
        #apply detector fix(es)
        if run_params['detector']['channel_correct'] or \
           run_params['detector']['picture_correct'] or \
           run_params['detector']['interchannel_correct']:
            try:
                fix_detector_levels([frame for frame, type in frames], subdir=folder, 
                        logger=dir_logger, **run_params['detector'])
            except:
                print('----detector correction failed')
                break_check('Det. correction', master_logger)
                continue
            else:
                master_logger.info('-- Finished detector corrections')
                master_logger.info(' ')
        else:
            master_logger.info('-- Detector corrections skipped')
            master_logger.info(' ')
   
   
   
        # -- reconstruct -----------------------
        master_logger.info('- Running reconstruction')
        ecode = reconstruct(logger=dir_logger, subdir=folder, obsband=band, 
                            **run_params['cubes'])
        if ecode:
            break_check('reconstruction', master_logger)
            continue
        master_logger.info('-- Finished reconstruction')
        master_logger.info(' ')
  
  
        # -- sky subtraction -------------------
        master_logger.info('- Running sky subtraction')
        try:
            ecode = subtract_sky(logger=dir_logger, subdir=folder, obsband=band,
                                 **run_params['cubes'])
        except Exception:
            master_logger.exception('sky subtraction raised an exception (likely no matching sky exposures for {0})'.format(folder))
            break_check('sky_subtraction', master_logger)
            continue
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
  
  
        # -- global background correction
        master_logger.info('- Applying global background correction')
        ecode = bkg_correct_2d(logger=dir_logger, subdir=folder, obsband=band,
                            search_suffix='s_i', **run_params['cubes'])
        if ecode:
            break_check('background correction', master_logger)
            continue
        master_logger.info('-- Finished background correction')
        master_logger.info(' ')
    
  
        # -- barycentric velocity correction -----
        master_logger.info('- Applying heliocentric correction')
        ecode = velocity_correct(logger=dir_logger, subdir=folder, obsband=band, 
                                 search_suffix='s_i_bkg2d', **run_params['cubes'])
        if ecode:
            break_check('velocity correction', master_logger)
            continue
        master_logger.info('-- Finished heliocentric correction')
        master_logger.info(' ')




        _save_checkpoint(procdir, band, folder, done_set)
        master_logger.info('- Finished reducing {0} '.format(folder))

