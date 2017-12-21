from __future__ import print_function
import numpy as np
import shutil
import astropy.io.fits as fits
import os, sys
import subprocess
import time

__all__ = ["clean_move", "get_obs_extensions", "fetch_resource_path", 
           "get_sof_list", "duplicate", "get_obs_exposures"]

dateobsfmt = "%Y-%m-%dT%H:%M:%S.%f"

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


def duplicate(full_path, outputdir=None, prefix=None, suffix=None):
    """ Copy a file from fullpath to outputdir, prepending/appending
    prefix or suffix if supplied.
    """
    splitfile = full_path.split('/')
    fpath, fname = '/'.join(splitfile[:-1]), splitfile[-1]
    ffile, fext = os.path.splitext(fname)

    if fpath == outputdir and prefix is None and suffix is None:
        raise ValueError("Change output directory or add prefix/suffix")
    
    if suffix is not None:
        soutname = '{0}_{1}{2}'.format(ffile, suffix, fext)
    else:
        soutname = fname

    if prefix is not None:
        outname = '{0}_{1}'.format(prefix, soutname)
    else:
        outname = soutname

    outpath = os.path.join(outputdir, outname)
    clean_copy(full_path, outpath)
    
    return outpath


def fetch_resource_path(resource):
    """ Fetch the path for a static resource included with the package.
    takes as input the relative path (usually data/XXX.fits) and returns 
    the full install path as a string.
    """
    d = os.path.dirname(sys.modules['pyspark'].__file__)
    return os.path.join(d, resource)


def get_sof_list(filetype=None, subdir=None, inputdir=None, obsband=None, 
                 darkexp=None, prefix=None, rotangle=None, logger=None, **kwargs):
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
        filepath = 'KMOS'
    else:
        filepath = '{0}_KMOS'.format(prefix)
    dfits_string = 'dfits '+filepath+'*{0}*.fits | fitsort tpl.id dpr.type ins.filt1.id det.seq1.dit ocs.rot.naangle date-obs obs.start obs.targ.name | tail -n +2'

    sof_out = []
   
    # --- Flat field frames --- 
    if filetype.lower() == 'flat':
        with cd(fpath):
            p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
            for line in iter(p.stdout.readline, b''):
                sline = line.decode('utf-8').strip().split(None)
                if sline[1] == 'KMOS_spec_cal_calunitflat' and sline[3].upper() in [obsband.upper(), 'BLOCK']:
                    if rotangle is None:
                        sof_out.append((os.path.join(fpath, sline[0]), sline[2]))
                    else:
                        if sline[2] == 'FLAT,OFF':
                            sof_out.append((os.path.join(fpath, sline[0]), sline[2]))
                        else:
                            rangle = float(sline[5])
                            if rangle < 0: rangle += 360.
                            if float(rotangle) == rangle:
                                sof_out.append((os.path.join(fpath, sline[0]), sline[2]))
            p.wait()

    # --- Arc line frames --- 
    elif filetype.lower() == 'arc':
        with cd(fpath):
            p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
            for line in iter(p.stdout.readline, b''):
                sline = line.decode('utf-8').strip().split(None)
                if sline[1] == 'KMOS_spec_cal_wave' and sline[3].upper() in [obsband.upper(), 'BLOCK']:
                    if rotangle is None:
                        sof_out.append((os.path.join(fpath, sline[0]), sline[2]))
                    else:
                        if sline[2] != 'WAVE,OFF':
                            rangle = float(sline[5])
                            if rangle < 0: rangle += 360.
                            if float(rotangle) == rangle:
                                sof_out.append((os.path.join(fpath, sline[0]), sline[2]))
    
            p.wait()

    # --- Darks frames --- 
    elif filetype.lower() == 'dark':
        #track darks with two different exposure times
        exp_dark = []
        standard_dark = []

        #first pass: identify darks and group by exposure time
        store_darks = {}
        with cd(fpath):
            p = subprocess.Popen(dfits_string.format('DARK'), stdout=subprocess.PIPE, bufsize=1, shell=True)
            for line in iter(p.stdout.readline, b''):
                sline = line.decode('utf-8').strip().split(None)
                store_darks.setdefault(int(float(sline[4])), []).append((os.path.join(fpath, sline[0]),None))
            p.wait()

        dark_times = list(store_darks.keys())
        if darkexp is not None and darkexp in dark_times:
            if logger is not None:
                logger.debug("- Found darks with {0}s".format(darkexp))
            sof_out = store_darks[darkexp]
        elif 60 in dark_times:
            if logger is not None:
                logger.debug("- {0}s darks not found, falling back to 60s".format(darkexp))
            sof_out = store_darks[60]
        else:
            #neither the provided time or 60s is in the dark folder
            if len(dark_times) == 1: #use the only time that's there
                if logger is not None:
                    logger.debug("- Only one set of darks found, using {0}s".format(dark_times[0]))
                sof_out = store_darks[dark_times[0]]
            else: #otherwise use the shortest time there
                if logger is not None:
                    logger.debug("- Multiple (non-specified) darks found, using {0}s".format(min(dark_times)))
                sof_out = store_darks[min(dark_times)]


    # --- Standard star exposures --- 
    elif filetype.lower() == 'standard':
        with cd(fpath):
            p = subprocess.Popen(dfits_string.format('CAL'), stdout=subprocess.PIPE, bufsize=1, shell=True)
            for line in iter(p.stdout.readline, b''):
                sline = line.decode('utf-8').strip().split(None)
                if obsband is not None:
                    if sline[1] in ['KMOS_spec_cal_stdstar','KMOS_spec_cal_stdstarscipatt'] and sline[3] == obsband.upper():
                        sof_out.append((os.path.join(fpath, sline[0]), sline[4], sline[7], sline[8], float(sline[5])))
                else:
                    if sline[1] in ['KMOS_spec_cal_stdstar','KMOS_spec_cal_stdstarscipatt']:
                        sof_out.append((os.path.join(fpath, sline[0]), sline[4], sline[7], sline[8], float(sline[5])))

            p.wait()

    # --- Science exposures --- 
    elif filetype.lower() == 'science':
        with cd(fpath):
            p = subprocess.Popen(dfits_string.format('OBS'), stdout=subprocess.PIPE, bufsize=1, shell=True)
            for line in iter(p.stdout.readline, b''):
                sline = line.decode('utf-8').strip().split(None)
                if sline[2] == 'OBJECT,SKY' and sline[3] == obsband.upper():
                    sof_out.append((os.path.join(fpath, sline[0]), sline[6]))
            p.wait()

    return sof_out


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
        filepath = 'sci_reconstructed_KMOS'
    else:
        filepath = 'sci_reconstructed{0}_KMOS'.format('_{0}'.format(prefix))
        expected_split += len(prefix.split('_'))

    if suffix is None:
        dfits_string = 'dfits '+filepath+'*.fits | fitsort date-obs ins.grat1.id ocs.arm1.name ocs.arm2.name ocs.arm3.name | tail -n +2'
    else:
        dfits_string = 'dfits '+filepath+'*{0}.fits | fitsort date-obs ins.grat1.id ocs.arm1.name ocs.arm2.name ocs.arm3.name | tail -n +2'.format('_{0}'.format(suffix))
        expected_split += len(suffix.split('_'))
    
    # run the query, sort things into object and sky
    # we actually temporarily change directories here
    # to deal with the length limit imposed by dfits
    exposures = []
    with cd(fpath):
        p = subprocess.Popen(dfits_string, stdout=subprocess.PIPE, bufsize=1, shell=True)
        for line in iter(p.stdout.readline, b''):
            sline = line.decode('utf-8').strip().split(None)
            if len(os.path.split(sline[0])[1].split('_')) == expected_split:
                if all([True if ent not in ['ARM1_SKY','ARM2_SKY','ARM3_SKY'] else False
                        for ent in sline[-3:]]) and sline[2] == obsband:
                    exposures.append((os.path.join(fpath, sline[0]), 'OBJ', time.mktime(time.strptime(sline[1], dateobsfmt))))
                elif sline[2] == obsband:
                    exposures.append((os.path.join(fpath, sline[0]), 'SKY', time.mktime(time.strptime(sline[1], dateobsfmt))))
        p.wait()

    return exposures


def verify_structure(outputdir=None, subdir=None):
    """ Routine which checks if the necessary folder structure exists
    if not, put it in place.
    """
    checkdirs = ['sky', 'science', 'calib', 'psf', 'sof', 'std']
    basedir = os.path.join(outputdir, subdir)

    if not os.path.exists(basedir):
        os.makedirs(basedir)
    for subdir in checkdirs:
        if not os.path.exists(os.path.join(basedir, subdir)):
            os.makedirs(os.path.join(basedir, subdir))

    return


class cd:
    """ A context manager class that can be used to temporarily change
    directories.
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


