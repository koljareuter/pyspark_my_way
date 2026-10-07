#!/usr/bin/env python
"""Preflight check before re-running sky subtraction with a corrected ZAP basis
(scripts/kmos_server_resky_<BAND>.cfg, see docs/zap-iz-svd-mean-bug.md).

For every OB folder in the cfg it checks that everything the re-run needs is in
PROCDIR, without touching any data:

  * the corrected eigenbasis (svd_path) exists and really has a zero mean spectrum
  * PROCDIR/<OB>/sky/ has un-sky-subtracted cubes of this band, incl. >= 1 sky frame
  * PROCDIR/<OB>/calib/ has MASTER_FLAT_<B><B><B>.fits and ILLUM_CORR_<B><B><B>_*.fits
  * raw science frames exist in RAWDIR/<OB>/ (reduce_kmos.py checks those first)
  * how many final *_s_i_bkg2d_rv.fits currently exist (will be overwritten)

  --backup   copy the current science/*_s_i_bkg2d_rv.fits of every OB to
             science_pre_zapfix/ (skipped where that copy already exists)

Usage (on hqempower1, from ~/pyspark_my_way):
    python scripts/check_resky_ready.py scripts/kmos_server_resky_IZ.cfg [--backup]
Exit code 0 if every OB is ready, 1 otherwise.
"""
import argparse
import configparser
import glob
import os
import re
import shutil
import sys

import numpy as np
from astropy.io import fits


def sky_frames(sky_dir, band):
    """(n_object, n_sky) un-sky-subtracted cubes of this band in sky/."""
    n_obj = n_sky = 0
    for f in glob.glob(os.path.join(sky_dir, "SCI_RECONSTRUCTED_*KMOS*.fits")):
        name = os.path.basename(f)
        if re.search(r"_s(_|\.fits)", name[len("SCI_RECONSTRUCTED_"):].split("KMOS", 1)[1]):
            continue                                   # not an un-subtracted cube
        try:
            h = fits.getheader(f, 0)
        except Exception:
            continue
        if (h.get("HIERARCH ESO INS GRAT1 ID") or "").strip() != band:
            continue
        arms = [str(h.get("HIERARCH ESO OCS ARM%d NAME" % i, "")) for i in (1, 2, 3)]
        if any(a.endswith("_SKY") for a in arms):
            n_sky += 1
        else:
            n_obj += 1
    return n_obj, n_sky


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cfg")
    ap.add_argument("--backup", action="store_true",
                    help="copy current science/*_s_i_bkg2d_rv.fits to science_pre_zapfix/")
    args = ap.parse_args()

    cp = configparser.ConfigParser()
    cp.read(args.cfg)
    bands = cp.get("instrument", "obsband").strip().split(",")
    if len(bands) != 1:
        sys.exit("ERROR: the re-sky cfg must contain exactly one band (svd_path applies to all bands)")
    band = bands[0]
    rawdir, procdir = cp.get("pathing", "rawdir"), cp.get("pathing", "procdir")
    subdirs = [s.strip() for s in re.split(",|\n", cp.get("pathing", "subdirs").strip()) if s.strip()]
    svd = cp.get("cubes", "svd_path")
    # rebuild mode: cubes are reconstructed from raw frames, so sky/ need not exist yet,
    # but the reconstruction calibrations must
    rebuild = cp.has_option("cubes", "reconstruct_cubes") and \
        cp.get("cubes", "reconstruct_cubes").strip() == "True"

    problems = []
    # ── corrected basis ──────────────────────────────────────────────────────
    if not os.path.isfile(svd):
        problems.append("svd_path not found: %s" % svd)
    else:
        with fits.open(svd) as h:
            if not np.all(np.asarray(h[1].data) == 0):
                problems.append("svd_path HDU1 (mean) is NOT zero: %s" % svd)
            if band not in os.path.basename(svd):
                problems.append("svd_path does not look like a %s basis: %s" % (band, svd))
    print("band %s   %d OB folders   basis %s%s" % (band, len(subdirs), svd,
          "   [rebuild: cubes reconstructed from raw, sky/ not required]" if rebuild else ""))
    print("%-30s %6s %6s %6s %6s %6s  %s" % ("OB folder", "raw", "obj", "sky", "illum", "rv", "status"))

    n_ok = 0
    for sub in subdirs:
        ob = os.path.join(procdir, sub)
        raw = len(glob.glob(os.path.join(rawdir, sub, "*.fits")))
        n_obj, n_sky = sky_frames(os.path.join(ob, "sky"), band)
        bbb = band * 3
        flat = os.path.isfile(os.path.join(ob, "calib", "MASTER_FLAT_%s.fits" % bbb))
        illum = len(glob.glob(os.path.join(ob, "calib", "ILLUM_CORR_%s_*.fits" % bbb)))
        rv = sorted(glob.glob(os.path.join(ob, "science", "*_s_i_bkg2d_rv.fits")))
        issues = []
        if raw == 0:
            issues.append("no raw frames")
        if rebuild:
            for cal in ("XCAL", "YCAL", "LCAL"):
                if not os.path.isfile(os.path.join(ob, "calib", "%s_%s.fits" % (cal, bbb))):
                    issues.append("%s missing" % cal)
        else:
            if n_obj == 0:
                issues.append("no object cubes in sky/")
            if n_sky == 0:
                issues.append("no sky cubes in sky/")
        if not flat:
            issues.append("MASTER_FLAT missing")
        if illum == 0:
            issues.append("ILLUM_CORR missing")
        status = "OK" if not issues else "; ".join(issues)
        print("%-30s %6d %6d %6d %6d %6d  %s" % (sub, raw, n_obj, n_sky, illum, len(rv), status))
        if issues:
            problems.append("%s: %s" % (sub, status))
        else:
            n_ok += 1
        if args.backup and rv:
            bdir = os.path.join(ob, "science_pre_zapfix")
            os.makedirs(bdir, exist_ok=True)
            for f in rv:
                dst = os.path.join(bdir, os.path.basename(f))
                if not os.path.exists(dst):
                    shutil.copy2(f, dst)

    print("\n%d/%d OB folders ready" % (n_ok, len(subdirs)))
    if args.backup:
        print("backup of *_s_i_bkg2d_rv.fits → <OB>/science_pre_zapfix/ done")
    if problems:
        print("\nPROBLEMS:")
        for p in problems:
            print("  - " + p)
        sys.exit(1)


if __name__ == "__main__":
    main()
