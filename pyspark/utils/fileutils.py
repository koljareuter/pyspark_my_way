from __future__ import print_function
import numpy as np
import shutil
import astropy.io.fits as fits
import os, sys
import subprocess


__all__ = ["clean_move", "read_cube", "get_cal_rotator", "match_obs_rotator", 
        "get_obs_extensions", "fetch_resource_path", "get_sof_list", "duplicate"]

def clean_move(src, dst):
    
    """Move a file, deleting the target if it already exists. If target is a directory, 
    preserve the original file name, otherwise move change name.

    Parameters
    ----------
    src: string
        Name of the source file.

    dst: string
        Name of the destination file or directory.

    """

    if os.path.isdir(dst): #assume copying same file to different directory
        if os.path.exists(os.path.join(dst,src)):
            os.remove(os.path.join(dst,src))
    else:
        if os.path.exists(dst):
            os.remove(dst)
    shutil.move(src, dst)
    return


def clean_copy(src, dst):
    
    """copy a file, deleting the target if it already exists. If target is a directory, 
    preserve the original file name, otherwise move change name.

    Parameters
    ----------
    src: string
        Name of the source file.

    dst: string
        Name of the destination file or directory.

    """

    if os.path.isdir(dst): #assume copying same file to different directory
        if os.path.exists(os.path.join(dst,src)):
            os.remove(os.path.join(dst,src))
    else:
        if os.path.exists(dst):
            os.remove(dst)
    shutil.copy(src, dst)
    return

def duplicate(full_path, outputdir=None, prefix=None):
    """ Copy a file from fullpath to outputdir, prepending 
    prefix if supplied.
    """

    splitfile = full_path.split('/')
    fpath, fname = '/'.join(splitfile[:-1]), splitfile[-1]

    if fpath == outputdir and prefix is None:
        raise ValueError("Change output directory or add prefix")
    
    if prefix is not None:
        outname = '{0}_{1}'.format(prefix, fname)
    else:
        outname = fname

    outpath = os.path.join(outputdir, outname)
    clean_copy(full_path, outpath)
    
    return outpath



def fetch_resource_path(resource):

    """ Fetch the path for a static resource included with the package.
    takes as input the relative path (usually data/XXX.fits) and returns 
    the full install path as a string.
    """
    
    d = os.path.dirname(sys.modules['pyspark'].__file__)
    #data = fits.open(os.path.join(d, resource))
    return os.path.join(d, resource)


def get_sof_list(filetype=None, subdir=None, inputdir=None, obsband=None, 
                 darkexp=None, prefix=None, rotangle=None, **kwargs):

    """ Return a list of frames matching the requested filetype and band
    in the given directory.

    Parameters
    ----------
    filetype: string
        File type to search for, should be one of arc, flat, dark, 
        standard, or science.

    input_dir: string
        Full path to the directory to search.

    obs_band: string
        Observing band to check for. Should be one of IZ, YJ, H, HK, K.

    exptime: int
        Exposure time constraint (only used for darks).  If no dark is 
        found with the given exposure time, this defaults to a 60s dark.

    """

    if subdir is not None:
        fpath = os.path.join(inputdir, subdir)
    else:
        fpath = inputdir

    if prefix is None:
        filepath = '{0}/KMOS'.format(fpath)
    else:
        filepath = '{0}/{1}_KMOS'.format(fpath, prefix)
    dfits_string = 'dfits '+filepath+'*{0}*.fits | fitsort tpl.id dpr.type ins.grat1.id det.seq1.dit ocs.rot.naangle date-obs | tail -n +2'
   

    sof_out = []
   
    # --- Flat field frames --- 
    if filetype.lower() == 'flat':
        p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if sline[1] == 'KMOS_spec_cal_calunitflat' and sline[3] == obsband.upper():
                if rotangle is None:
                    sof_out.append((sline[0], sline[2]))
                else:
                    if sline[2] == 'FLAT,OFF':
                        sof_out.append((sline[0], sline[2]))
                    else:
                        rangle = float(sline[5])
                        if rangle < 0: rangle += 360.
                        if float(rotangle) == rangle:
                            sof_out.append((sline[0], sline[2]))
        p.wait()

    # --- Arc line frames --- 
    elif filetype.lower() == 'arc':
        p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if sline[1] == 'KMOS_spec_cal_wave' and sline[3] == obsband.upper():
                #sof_out.append((sline[0], sline[2]))
                if rotangle is None:
                    sof_out.append((sline[0], sline[2]))
                else:
                    if sline[2] != 'WAVE,OFF':
                        rangle = float(sline[5])
                        if rangle < 0: rangle += 360.
                        if float(rotangle) == rangle:
                            sof_out.append((sline[0], sline[2]))

        p.wait()

    # --- Darks frames --- 
    elif filetype.lower() == 'dark':
        #track darks with two different exposure times
        exp_dark = []
        standard_dark = []

        p = subprocess.Popen(dfits_string.format('DARK'), stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if darkexp is not None and float(sline[4]) == darkexp:
                exp_dark.append((sline[0], None))
            elif float(sline[4]) == 60.:
                standard_dark.append((sline[0], None))
        p.wait()
        
        if len(exp_dark):
            sof_out = exp_dark
        else:
            sof_out = standard_dark

    # --- Standard star exposures --- 
    elif filetype.lower() == 'standard':
        p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if sline[1] in ['KMOS_spec_cal_stdstar','KMOS_spec_cal_stdstarscipatt'] and sline[3] == obs_band.upper():
                sof_out.append((sline[0], sline[2]))
        p.wait()

    # --- Science exposures --- 
    elif filetype.lower() == 'science':
        p = subprocess.Popen(dfits_string.format('OBS'), stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if sline[2] == 'OBJECT,SKY' and sline[3] == obsband.upper():
                sof_out.append((sline[0], sline[6]))
        p.wait()

    return sof_out


def verify_structure(outputdir=None, subdir=None):
    #checks if the necessary folder structure exists
    #if not, put it in place

    checkdirs = ['sky', 'science', 'calib', 'psf', 'sof', 'std']
    basedir = os.path.join(outputdir, subdir)

    if not os.path.exists(basedir):
        os.makedirs(basedir)
    for subdir in checkdirs:
        if not os.path.exists(os.path.join(basedir, subdir)):
            os.makedirs(os.path.join(basedir, subdir))

    return
