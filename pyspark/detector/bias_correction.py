from __future__ import print_function
import numpy as np
import bottleneck as bn
import os
import astropy.io.fits as fits

from ..utils.obsutils import hl_diff
from ..utils.fileutils import (fetch_resource_path, clean_copy, duplicate)
from ..utils.setuputils import get_date

from .._version import __version__



__all__ = ["get_channel_correct", "get_picture_correct", "fix_detector_levels"]


def get_channel_correct(input_frame):
    """Given an input exposure, estimate and remove individual H2RG channel 
    drifts using reference pixels.  Applies both a row and column correction.
    
    Parameters
    ----------
    input_frame : ndarray
        Input (raw) detector image.

    Returns
    -------
    ref_correct : ndarray
        Detector image corrected for row and column shifts.
   """

    #first pixel is always a bit funny
    input_frame[2047,0] = np.nan

    #rip out the reference columns (skip the first and last)
    cols = np.zeros((2048,8), dtype=np.float32)
    cols[:,:4] = input_frame[:,:4]
    cols[:,4:] = input_frame[:,-4:]

    #compute the mean
    col = bn.nanmean(np.vstack((bn.nanmedian(cols[:,:4],axis=1),\
            bn.nanmedian(cols[:,4:],axis=1))), axis=0)

    #estimate mean of column pixels
    col -= bn.nanmean(col[4:-4])
    col[~np.isfinite(col)] = 0.

    #reconstruct the 2d column correction
    col[:4], col[-4:] = 0.0, 0.0
    ccorr = np.atleast_2d(col).T
    ccorr = np.repeat(ccorr, 2048, axis=1)
    ccorr[:,:4], ccorr[:,-4:] = 0., 0.
    ccorr[:4,:], ccorr[-4:,:] = 0., 0.

    ##produce the column-corrected cube
    rows = input_frame[-4:,:]

    oe_oned = np.zeros(2048, dtype=np.float32)
    for ii in range(1,32):
        this_channel = rows[:,ii*64:ii*64+64]
        prev_channel = rows[:,(ii-1)*64:(ii-1)*64+64]
        e_channel = this_channel[:,::2].flatten()
        o_channel = this_channel[:,1::2].flatten()

        oe_diff = hl_diff(e_channel, o_channel)
        oe_corr = np.r_[e_channel, o_channel+oe_diff]
        diff = hl_diff(prev_channel.flatten(), oe_corr)

        oe_oned[ii*64:ii*64+64:2] = -diff
        oe_oned[ii*64+1:ii*64+64:2] = -(diff+oe_diff)
        rows[:,ii*64:ii*64+64:2] += diff
        rows[:,ii*64+1:ii*64+64:2] += (diff+oe_diff)
        
    oe_oned -= bn.nanmean(oe_oned)
    oe_oned[:4], oe_oned[-4:] = 0.0, 0.0

    rcorr = np.atleast_2d(oe_oned)
    rcorr = np.repeat(rcorr, 2048, axis=0)

    rcorr[:,:4], rcorr[:,-4:] = 0., 0.
    rcorr[:4,:], rcorr[-4:,:] = 0., 0.

    #corrected = np.copy(input_frame) - rcorr + ccorr
    return rcorr-ccorr


def get_picture_correct(input_frame, detector=None):

    """ Apply correction for the H2RG "picture frame" effect based
    on the correlation between reference pixel values and residual
    bias signal.  Correction is based on a stack of a few thousand
    60s dark frames.
    """

    #fetch the superdark reference frame for this detector
    fpath = fetch_resource_path('data/superdark_coeffs_mean.fits.fz') 

    #pull the information for the picture frame correction
    sdfile = fits.open(fpath)
    picture_coeff = np.array(sdfile[detector].data, dtype=np.float32)
    sdfile.close() 

    #lookup values for reference classification
    idx, idy = np.meshgrid(range(2048), range(2048))
    ref_select = ((idx < 4) | (idx > 2043) | (idy > 2043)) & (idy > 3)
    ref_vals = np.array([-100, -4.5, -1.5, 1.5, 4.5, 100])

    #picture frame scaling based on median of input reference pixels
    reference_coeff = np.median(input_frame[ref_select & (input_frame != 0.) & np.isfinite(input_frame)])

    return picture_coeff*reference_coeff


def fix_detector_levels(frame_list, subdir=None, inputdir=None, outputdir=None, prefix=None, 
        channel_correct=True, picture_correct=True, logger=None, **kwargs):

    """ Given an input set of frames, apply detector-level 
    corrections and write out the new set of frames.
    """

    nframes = len(frame_list)

    #iterate over each frame and apply requested corrections
    for frame in frame_list[:3]:

        #make a copy of the input frame to outputdir/subdir
        if logger is not None:
            logger.info('-- Copying {0} to {1}'.format(frame.split('/')[-1], '_'.join([prefix, frame.split('/')[-1]])))
        outname = duplicate(frame, outputdir=os.path.join(outputdir, subdir), prefix=prefix)
    

        
        #perform any detector corrections (if requested)
        if channel_correct or picture_correct:
            fdata = fits.open(outname, mode='update')

            fdata[0].header['history'] = 'Processed by PySPARK {0} {1}'.format(__version__, get_date())

            #iterate over each detector (assumes KMOS)
            for det in range(1,4):
                channel_array = 0.
                picture_array = 0.

                orig_data = np.array(fdata[det].data)
                orig_header = fdata[det].header

                if channel_correct:
                    logger.debug('-- Applying channel correction to {0}'.format('_'.join([prefix, frame.split('/')[-1]])))
                    channel_array = get_channel_correct(orig_data)
                    orig_header['history'] = 'Applied channel correction {0}'.format(get_date()) 

                if picture_correct:
                    logger.debug('-- Applying picture frame correction to {0}'.format('_'.join([prefix, frame.split('/')[-1]])))
                    picture_array = get_picture_correct(orig_data, detector=det)
                    orig_header['history'] = 'Applied picture frame correction {0}'.format(get_date()) 

                #update the data in place
                fdata[det].data = orig_data - channel_array - picture_array
            
            fdata.flush()
    return
#
#                
#
#
##if __name__ == '__main__':
##    
##    #### OPTIONS ####
##    parser = argparse.ArgumentParser()
##    parser.add_argument('dirs', nargs='+', help='list of Object names to process')
##    parser.add_argument('--verbose', type=int, default=3, help='Set logging level for stdout')
##    parser.add_argument('--raw_dir', default=os.getcwd(), type=str, help='Absolute path to raw data')
##    parser.add_argument('--res_dir', default=os.getcwd(), type=str, help='Absolute path to results')
##    args = parser.parse_args()
##
##    #set up directories
##    dirs = args.dirs
##    raw_dir = args.raw_dir
##    res_dir = args.res_dir
##
##    #set up logging.
##    loglevels = {0:logging.CRITICAL, 1:logging.ERROR, 2:logging.WARNING, 3:logging.INFO, 4:logging.DEBUG}
##
##    #screen logging level is set by the verbose flag.  Default to file is always DEBUG
##    log_level = args.verbose
##    logger = logging.getLogger()
##    logger.setLevel(logging.DEBUG)
##    
##    #file output
##    if os.path.exists('make_lcorr.log'):
##        os.remove('make_lcorr.log')
##    fh = logging.FileHandler('make_lcorr.log')
##    fh.setLevel(logging.DEBUG)
##    
##    ch = logging.StreamHandler()
##    ch.setLevel(loglevels[log_level])
##    
##    formatter = logging.Formatter('%(asctime)s [%(levelname)-8s]: %(message)s', datefmt='%m/%d/%Y %I:%M:%S')
##    fh.setFormatter(formatter)
##    ch.setFormatter(formatter)
##    
##    logger.addHandler(ch)
##    logger.addHandler(fh)
##
##
##    #### END OPTIONS ####
##
##    ###########################################################
##    # FETCH FILES TO CORRECT
##    #
##    #fetch list of valid observation files and store a bit of information
##
##    obs_files, obs_night = [], []
##    for ddir in dirs:
##        cmd = 'dfits {0}/{1}/KMOS*OBS*.fits | fitsort ins.grat1.id dpr.type | tail -n +2'.format(raw_dir, ddir)
##        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, shell=True)
##        for line in iter(p.stdout.readline, b''):
##            temp = string.strip(line).split(None)
##            if temp[2] in ['OBJECT,SKY','DARK']:
##                obs_files.append(temp[0])
##                obs_night.append(ddir)
##                logger.debug('Found {0}/{1}'.format(ddir,temp[0].split('/')[-1])) 
##
##
##    #pixel lookup for superdark correction
##    idx, idy = np.meshgrid(range(2048), range(2048))
##    ref_select = ((idx < 4) | (idx > 2043) | (idy > 2043)) & (idy > 3)
##    
##    #load in the superdark file to apply correction 
##    sdfile = fits.open('{0}/data/superdark_coeffs_mean.fits.fz'.format(os.path.dirname(kutils.__file__)))
##
##    #go through the science frames one by one
##    for ii, sci_file in enumerate(obs_files[:1]):
##        this_file = sci_file.split('/')[-1]
##        logger.info('Estimating correction for {0}'.format(this_file))
##        out_file = 'ndc_{0}'.format(this_file)
##       
##        if not os.path.exists('{0}/{1}'.format(res_dir, obs_night[ii])):
##            os.makedirs('{0}/{1}'.format(res_dir, obs_night[ii])) 
##
##        file = fits.open(sci_file)
##        hdulist = fits.HDUList()
##        phdu = fits.PrimaryHDU(header=file[0].header)
##        hdulist.append(phdu)
##        for detector in range(1,4):
##            logger.debug('- processing detector {0}'.format(detector))
##            obj_data = np.asarray(file[detector].data, dtype=np.float32)
##
##            #reference pixel values for superdark correction
##            edge_vals = np.median(obj_data[ref_select & (obj_data != 0.) & np.isfinite(obj_data)])
##            sdark_corr = sdfile[detector].data * edge_vals
##
##            ref_out = _ref_corr(obj_data, det=detector) - sdark_corr
##
##            ##run second correction on image pixels to account for drift
##            #det_corr = _det_corr(ref_out, np.asarray(smfile[detector].data, dtype=bool)) 
##            
##            data_image = fits.ImageHDU(data=np.asarray(ref_out, dtype=np.float32), header=file[detector].header)
##            hdulist.append(data_image)
##       
##            #logger.debug('-- median reference value: {0}'.format(edge_vals)) 
##
##        hdulist.writeto('{0}/{1}/{2}'.format(res_dir, obs_night[ii], out_file), overwrite=True, output_verify='ignore')
##        file.close()
##    
##    sdfile.close()
##    smfile.close()
