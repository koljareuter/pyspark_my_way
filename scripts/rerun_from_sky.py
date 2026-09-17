#!/usr/bin/env python
"""Re-run the Python-only pipeline steps from sky subtraction onwards.

Use this when EsoRex steps (calibrations, reconstruction) are already done
and you only want to test changes to sky subtraction, illumination, background
correction, or velocity correction.

The script removes existing science/ files for the specified OB and re-runs:
  subtract_sky → illuminate → bkg_correct_2d → velocity_correct

Run from the project root in WSL:
    PYTHONPATH=. python scripts/rerun_from_sky.py -p scripts/kmos_local.cfg
"""
from __future__ import print_function
import argparse, glob, os, sys

from pyspark.utils import fileutils, setuputils
from pyspark.cubes import illuminate, bkg_correct_2d, velocity_correct
from pyspark.sky import subtract_sky
from pyspark.sky.sky_subtraction import build_custom_svd


def _remove_science_files(procdir, subdir, band):
    """Delete existing science/ products so they are rebuilt cleanly."""
    sci_dir = os.path.join(procdir, subdir, 'science')
    pattern  = os.path.join(sci_dir, '*_s*.fits')
    removed  = 0
    for f in glob.glob(pattern):
        os.remove(f)
        removed += 1
    return removed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--parameter_file', required=True,
                        help='Config file (same one used for reduce_kmos.py)')
    parser.add_argument('--verbose', type=int, default=3)
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be removed/re-run without doing it')
    args = parser.parse_args()

    logger = setuputils.set_logging(logbase='rerun', log_level=args.verbose)
    params = setuputils.set_config(args.parameter_file, logger=logger)

    procdir = params['procdir']
    rawdir  = params['rawdir']
    folders = params['subdirs']
    bands   = params['obsband']

    for band in bands:
        for folder in folders:
            logger.info('=== {0} / {1} ==='.format(band, folder))

            # check pre-sky cubes exist
            sky_dir = os.path.join(procdir, folder, 'sky')
            pre_sky = glob.glob(os.path.join(sky_dir, 'SCI_RECONSTRUCTED_*.fits'))
            if not pre_sky:
                logger.warning('No reconstructed cubes in sky/ — run the full pipeline first')
                continue
            logger.info('Found {0} pre-sky cubes in sky/'.format(len(pre_sky)))

            # remove existing science/ outputs
            n = _remove_science_files(procdir, folder, band)
            if args.dry_run:
                logger.info('[dry-run] Would remove {0} science/ files'.format(n))
                continue
            logger.info('Removed {0} existing science/ files'.format(n))

            # Use the generic bundled SVD (purpose-built, wavelength-consistent).
            # Custom per-OB or multi-OB SVDs caused grating-drift artefacts across
            # nights and over-subtracted galaxy emission lines — abandoned.
            cubes_params = dict(params['cubes'])
            cubes_params.pop('svd_path', None)
            logger.info('Using generic ZAP SVD (+ bracketed sky + OH scaling)')

            # --- sky subtraction (with bracketing on pipeline-review) ---
            logger.info('Running sky subtraction ({0})...'.format(
                cubes_params.get('sky_method', 'zapv2')))
            try:
                ecode = subtract_sky(logger=logger, subdir=folder, obsband=band,
                                     **cubes_params)
            except Exception:
                logger.exception('Sky subtraction failed')
                continue
            if ecode:
                logger.error('Sky subtraction returned error code {0}'.format(ecode))
                continue
            logger.info('Sky subtraction done')

            # --- illumination correction ---
            logger.info('Running illumination correction...')
            ecode = illuminate(logger=logger, subdir=folder, obsband=band,
                               search_suffix='s', **params['cubes'])
            if ecode:
                logger.error('Illumination correction failed')
                continue
            logger.info('Illumination done')

            # --- 2D background correction ---
            logger.info('Running 2D background correction...')
            ecode = bkg_correct_2d(logger=logger, subdir=folder, obsband=band,
                                   search_suffix='s_i', **params['cubes'])
            if ecode:
                logger.error('Background correction failed')
                continue
            logger.info('Background correction done')

            # --- heliocentric velocity correction ---
            logger.info('Running velocity correction...')
            ecode = velocity_correct(logger=logger, subdir=folder, obsband=band,
                                     search_suffix='s_i_bkg2d', **params['cubes'])
            if ecode:
                logger.error('Velocity correction failed')
                continue
            logger.info('Velocity correction done')

            # final inventory
            final = glob.glob(os.path.join(procdir, folder, 'science', '*_rv.fits'))
            logger.info('Final products: {0} x _s_i_bkg2d_rv.fits'.format(len(final)))
            logger.info('=== Done ===')
