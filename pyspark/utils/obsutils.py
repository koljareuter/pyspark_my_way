from __future__ import print_function
import bottleneck as bn
import numpy as np
import shutil
import subprocess

__all__ = ["nanMAD", "hl_diff", "hsm", "read_cube", "match_obs_rotator",
           "get_cal_rotator", "get_obs_extensions", "run_spark", "find_lines"]

def nanMAD(a, axis=None):
    """Compute the Median Absolute Deviation accounting
    for Inf and NaN values.
    """
    mu = bn.nanmedian(a, axis=axis)
    if axis is not None:
        mu = np.expand_dims(mu, axis=axis)    
    return bn.nanmedian(np.fabs(a-mu), axis=axis)

def hl_diff(a,b):
    """Hodges-Lehmann estimate of the difference between
    two input arrays.
    """
    pairs = a[:,None]-b[None,:]
    return bn.nanmedian(pairs.flatten())

def hsm(a, min_samp=5):
    """Estimate the mode of an input array using the
    the half-sample algorithm.
    """
    
    srt = np.argsort(a.ravel())
    cut_array = a.ravel()[srt]
    nsamp_full = cut_array.size

    power = 1
    while nsamp_full/float((2**(power+1))) > min_samp:
        interval = nsamp_full/(2**power)
        
        #find minimum interval
        cut_low, cut_high = cut_array[:interval], cut_array[interval:]
        #print cut_low, cut_high, nsamp_full
        if cut_low.size > cut_high.size:
            interval_size = cut_high.size
            cut_low = cut_low[:interval_size]
        elif cut_high.size > cut_low.size:
            interval_size = cut_low.size
            cut_high = cut_high[:interval_size]
        else:
            interval_size = interval

        min_interval = np.argmin(cut_high-cut_low)

        #cut down the sample
        cut_array = cut_array[min_interval:min_interval+interval_size]
        power += 1

    return np.mean(cut_array)


def read_cube(data_file, ext=None, fixzero=False):
    
    """ Read in a reconstructed data cube.  Return data, reconstructed wavelength
    array, size information, and header.
    """
    
    file = fits.open(data_file)
    data = file[ext].data
    data = np.asarray(data, dtype=np.float64)
    
    #pull header information for the wavelength array
    hdr = file[ext].header
    try:
        crval = hdr['CRVAL3']
        cdelt = hdr['CDELT3']
        crpix = hdr['CRPIX3']
        xsize = hdr['NAXIS1']
        ysize = hdr['NAXIS2']
        zsize = hdr['NAXIS3']

        file.close()

        #wavelength array in microns
        wl_array = (np.arange(zsize, dtype=np.float64)+1 - crpix)*cdelt + crval

        #set zero values to nan for proper handling
        findzero = np.zeros(data.shape, dtype=bool)
        if fixzero:
            findzero = (data == 0.0)
            if np.count_nonzero(findzero):
                data[findzero] = np.nan
     
        return wl_array, data, xsize, ysize, zsize, hdr, findzero
    except:
        return np.array([]), np.array([]), 0, 0, 0, hdr, None


def get_cal_rotator(cal_file):

    """Fetch the rotation angles from a KMOS calibration file
    returns an array with the observed angles.
    """
    
    pro_naangles = []
    for iter in range(1,100):
        try:
            head = cal_file[iter].header
            pro_naangles.append(head['HIERARCH ESO PRO ROT NAANGLE'])
        except:
            break
    pro_naangles = np.unique(np.array(pro_naangles))
    return pro_naangles


def match_obs_rotator(obs_file, cal_rot, get_second=False):
    
    """find closest match in a list of calibration rotation angles
    returns appropriate rotator index and closest match.
    """

    head = obs_file[0].header
    obs_rot = head['HIERARCH ESO OCS ROT NAANGLE']

    #find closest from calibration rotator angles
    offset = np.abs(((cal_rot-obs_rot + 540.)%360-180))
    srt = np.argsort(offset)
    rotator = srt[0]
    this_frame_match = int(cal_rot[rotator])

    if obs_rot < 0: obs_rot += 360.

    if get_second: #also get extension for the second closest frame
        return rotator, this_frame_match, (obs_rot), srt[1], int(cal_rot[srt[1]])
    else:
        return rotator, this_frame_match, (obs_rot)


def get_obs_extensions(obs_file):

    """find which extensions (arms for KMOS) contain objects
    returns list of filled IFUs
    """
    
    filled_ifus = []
    for ii in range(1,25):
        try:
            used = obs_file[0].header['HIERARCH ESO OCS ARM{0} NOTUSED'.format(ii)]
        except:
            filled_ifus.append(ii)

    return filled_ifus


def run_spark(routine, args='', sof=None, logger=None):
    
    """ Execute SPARK routine given by cmd.
    """

    cmd = 'esorex {0} {1} {2}'.format(routine, args, sof)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, shell=True)
    for line in iter(p.stdout.readline, b''):
        if logger is not None:
            logger.debug(line.decode('utf-8').strip())
        else:
            continue
    p.wait()
    
    return p.returncode


def find_lines(spec_in):
    """returns an boolean array which selects possible sky line pixels
    """

    #compute derivative spectrum
    derivatives = np.r_[np.array([0.]),np.diff(spec_in)]
    derivatives[~np.isfinite(spec_in) | ~np.isfinite(derivatives)] = 0.
    
    #array of line type flags
    dflags = np.zeros(len(derivatives), dtype=int)

    #identify pixels in the line wings by their large first derivatives
    dnon_zero = np.abs(derivatives[derivatives != 0])

    med = np.median(dnon_zero)
    sig = 1.4826*1.253*np.median(np.abs(dnon_zero-med))

    #work out relative amplitude to detect line given S/N and width
    if sig/med > 4.:
        der_lim = 3.*med
    else:
        der_lim = sig

    dlow, dhigh = -1.5*der_lim, 1.5*der_lim
    dflags[derivatives < dlow] = -1
    dflags[derivatives > dhigh] = 1

    #perform two searches of the flagged pixels, one from the left and one from the right.
    search_right = np.copy(dflags)
    search_left = np.copy(dflags)

    range_flag = 0
    for ii in range(0,len(derivatives)-1):
        if dflags[ii] == 1:
            range_flag = 1
        if (dflags[ii]==-1 and dflags[ii+1]==0) or (dflags[ii]==-1 and dflags[ii+1]==1): 
            range_flag = 0
        search_right[ii] = range_flag

    for ii in range(len(derivatives)-2,0,-1):
        if dflags[ii] == -1:
            range_flag = 1
        if (dflags[ii+1]==1 and dflags[ii]==0) or (dflags[ii+1]==1 and dflags[ii]==-1): 
            range_flag = 0
        search_left[ii] = range_flag
    
    lrange = (search_right == 1) & (search_left == 1)

    return lrange

