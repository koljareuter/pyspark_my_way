from __future__ import print_function

import subprocess
import numpy as np
import os, sys
import time
import bottleneck as bn
import astropy.io.fits as fits

from .zap import (process, SVDoutput)
from ..utils.obsutils import find_lines


__all__ = ["get_obs_exposures", "subtract_sky", "match_obj_sky", "make_SVDs"]

sky_methods = ['zap']

dateobsfmt = "%Y-%m-%dT%H:%M:%S.%f"



def get_obs_exposures(inputdir=None, subdir=None, obsband=None, prefix=None, 
		      suffix=None, **kwargs):
    """ Get list of reconstructed frames (assumes file naming convention from
    kmo_reconstruct, i.e. sci_reconstructed_***.fits.
    """
    if subdir is not None:
        fpath = os.path.join(inputdir, subdir)
    else:
        fpath = inputdir

    #keep track of the expected number of name elements to identify correct frames
    expected_split = 6
    if prefix is None:
        filepath = '{0}/sci_reconstructed_KMOS'.format(fpath)
    else:
        filepath = '{0}/sci_reconstructed{1}_KMOS'.format(fpath, '_{0}'.format(prefix))
        expected_split += len(prefix.split('_'))

    if suffix is None:
        dfits_string = 'dfits '+filepath+'*.fits | fitsort date-obs ins.grat1.id ocs.arm1.name ocs.arm2.name ocs.arm3.name | tail -n +2'
    else:
        dfits_string = 'dfits '+filepath+'*{0}.fits | fitsort date-obs ins.grat1.id ocs.arm1.name ocs.arm2.name ocs.arm3.name | tail -n +2'.format('_{0}'.format(suffix))
        expected_split += len(suffix.split('_'))
    
    #run the query, sort things into object and sky
    exposures = []
    p = subprocess.Popen(dfits_string, stdout=subprocess.PIPE, bufsize=1, shell=True)
    for line in iter(p.stdout.readline, b''):
        sline = line.decode('utf-8').strip().split(None)
        if len(os.path.split(sline[0])[1].split('_')) == expected_split:
            if all([True if ent not in ['ARM1_SKY','ARM2_SKY','ARM3_SKY'] else False
                    for ent in sline[-3:]]) and sline[2] == obsband:
                exposures.append((sline[0], 'OBJ', time.mktime(time.strptime(sline[1], dateobsfmt))))
            elif sline[2] == obsband:
                exposures.append((sline[0], 'SKY', time.mktime(time.strptime(sline[1], dateobsfmt))))
    p.wait()

    return exposures


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


def make_SVDs(spec_array, header, pcaoutfits=None, zlevel='none', cftype='med_smooth', 
              clean=False, cfwidth=101, boost=50, reject_factor=0.4, logger=None):
    """ Given an input (2D) array of spectra, generate PCA output following
    zap_v2.0 conventions.
    """

    #generate mask array for sky line rejection
    med_sky = bn.nanmedian(spec_array, axis=1)
    msky = find_lines(med_sky)

    #run zap instance
    zap_inst = SVDoutput(None, cube='dummy', stack=spec_array, hdr=header, instrument='KMOS', 
                         zlevel=zlevel, cftype=cftype, clean=False, cfwidth=cfwidth, boost=boost, 
                         reject_factor=reject_factor, contmask=msky)

    if pcaoutfits is not None:
        zap_inst.writePCA(pcaoutputfits=pcaoutfits, overwrite=True)

    

def subtract_sky(logger=None, subdir=None, obsband=None, 
                 sky_method='zap', **kwargs):

    #check that sky subtraction is handled by this routine
    if sky_method not in sky_methods:
        return

    #get list of frames to process. Assumes location of un-sky-subtracted frames in "sky" directory
    exposures = get_obs_exposures(subdir=os.path.join(subdir, 'sky'), obsband=obsband, **kwargs) 

    #generate associated list of objects and skies
    obj_sky = match_obj_sky(exposures)

    if sky_method == 'zap':
        pass

#def get_obs_exposures(inputdir=None, subdir=None, prefix='', suffix='', **kwargs):

    #            #get list of flat frames for processing
    #                frames = get_sof_list(filetype='science', subdir=subdir, obsband=obsband, **kwargs)

    pass

