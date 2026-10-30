#!/usr/bin/env python
"""Set up a fresh PROCDIR for a re-reduction from the detector-corrected frames (links only) and report its progress.

For every OB folder of the given pyspark config(s) this creates <procdir of the cfg>/<OB>/ with
    calib/<file>         -> realpath(SRC/<OB>/calib/<file>)     every calibration file (XCAL, YCAL, LCAL, MASTER_FLAT, ILLUM_CORR ...)
    dc_KMOS.<...>.fits   -> realpath(SRC/<OB>/dc_KMOS.<...>.fits)
as absolute links straight to the real files (CPL follows only one link hop, docs/symlink-two-hop-kmos-flat-bugfix.md).
Nothing in SRC is created, changed or removed. sky/, science/, sof/ are created by reduce_kmos.py. Existing correct links
are kept; a path that exists but is not the expected link is reported and left alone.

--reuse DONE re-runs only the sky subtraction and the later steps on top of a finished run DONE (e.g. PROCDIR_HYB):
it also links DONE/<OB>/dcr_KMOS.*.fits, DONE/<OB>/sky/SCI_RECONSTRUCTED_*.fits and DONE/<OB>/skysub_det/SCI_RECONSTRUCTED_*.fits
(again straight to the real files). The cfg must then have reconstruct_cubes = False, no detector or row common-mode
step and no calibration step, so that nothing writes onto these links and the calibrations match the reused cubes; this
is checked. subtract_sky copies the linked sky/ cubes to the new science/ and reads the skysub_det/ cubes, so DONE is
only read (do not write-protect DONE: the copies inherit the file mode of the linked cubes).

Usage (on hqempower1, from ~/pyspark_my_way):
    python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg --dry-run    # check only
    python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg              # create the links
    python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg --status     # after/while reducing
    python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ.cfg --reuse /local_storage/KOLJA/PROCDIR_HYB
--status lists per OB the dc frames, sky/ cubes, final *_s_i_bkg2d_rv.fits cubes, the OH_SPEC in the saved sof and the
polynomial degree the patched DRP reported in the OB log (must equal lcorr_degree of the cfg); for the noise-fix
configs (kmos_server_hyb_*.cfg) also the row-corrected dcr_ frames and the skysub_det/ cubes of the hybrid sky
subtraction (docs/hybrid-sky.md; 0 for OBs whose SKY frames are at the OBJECT rotator angle), and in 'wing' the
final cubes made with the hybrid as <number>xW<hybrid_wing>, followed by /<oh_model>/<free_model> when those are
not the defaults (YJ). For a --reuse procdir, OH_SPEC and degree come from the reconstruction that made the linked
sky/ cubes (marked *).
Exit code 1 if an OB lacks dc frames or calibrations for its band, or a path conflict was found.
"""
import argparse
import configparser
import glob
import os
import re
import sys
from collections import Counter

from astropy.io import fits

REQUIRED = ("XCAL", "YCAL", "LCAL", "MASTER_FLAT")


def read_cfg(path, detector=False):
    c = configparser.ConfigParser()
    if not c.read(path):
        sys.exit("cannot read %s" % path)
    bands = [b.strip() for b in c.get("instrument", "obsband").split(",")]
    obs = [s.strip() for s in re.split(",|\n", c.get("pathing", "subdirs")) if s.strip()]
    cubes = dict(c.items("cubes")) if c.has_section("cubes") else {}
    if detector:
        return bands, c.get("pathing", "procdir").strip(), obs, cubes, \
            dict(c.items("detector")) if c.has_section("detector") else {}, \
            dict(c.items("calibration")) if c.has_section("calibration") else {}
    return bands, c.get("pathing", "procdir").strip(), obs, cubes


def link(src, dst, dry, problems):
    real = os.path.realpath(src)
    if os.path.islink(dst):
        if os.readlink(dst) == real:
            return "kept"
        problems.append("%s: links to %s, expected %s" % (dst, os.readlink(dst), real))
        return "conflict"
    if os.path.exists(dst):
        problems.append("%s: exists and is not a link" % dst)
        return "conflict"
    if not dry:
        os.symlink(real, dst)
    return "new"


def prepare(cfg, src_root, dry, problems, reuse_root=None):
    bands, dst_root, obs, cubes, det, calib = read_cfg(cfg, detector=True)
    if len(bands) != 1:
        sys.exit("%s: one band per config expected, got %s" % (cfg, bands))
    B3 = bands[0] * 3
    real_dst = os.path.realpath(dst_root)
    for root in [src_root] + ([reuse_root] if reuse_root else []):
        if real_dst == root or real_dst.startswith(root + os.sep):
            sys.exit("%s: procdir %s is (inside) the source PROCDIR %s - refused" % (cfg, dst_root, root))
    if reuse_root:
        # every step that writes dcr_ frames or sky/, skysub_det/ cubes must be off, or it would write onto the links
        # (calibrations default to on in make_spark_calibs: new ones would not match the reused sky/ cubes)
        if cubes.get("reconstruct_cubes", "True").strip() != "False" or \
                det.get("row_commonmode_correct", "False").strip() != "False" or \
                any(det.get(k, "False").strip() != "False" for k in ("channel_correct", "picture_correct", "interchannel_correct")) or \
                any(calib.get(k, "True").strip() != "False" for k in ("make_dark", "make_flat", "make_arc", "make_illum",
                                                                       "reconstruct_arc", "make_standard")):
            sys.exit("%s: --reuse needs reconstruct_cubes = False, all detector steps (row common mode, channel, "
                     "picture, interchannel) off and all [calibration] make_* / reconstruct_arc = False - refused" % cfg)
        prefix = det.get("prefix", "").strip()
    counts = Counter()
    for ob in obs:
        s, d = os.path.join(src_root, ob), os.path.join(dst_root, ob)
        dc = sorted(glob.glob(os.path.join(s, "dc_KMOS.*.fits")))
        cal = sorted(f for f in glob.glob(os.path.join(s, "calib", "*")) if os.path.isfile(f))
        missing = ["%s_%s.fits" % (c, B3) for c in REQUIRED if not os.path.isfile(os.path.join(s, "calib", "%s_%s.fits" % (c, B3)))]
        nill = len(glob.glob(os.path.join(s, "calib", "ILLUM_CORR_%s_*.fits" % B3)))
        if not dc or missing or nill == 0:
            problems.append("%s: dc frames %d, missing %s, ILLUM_CORR %d" % (ob, len(dc), missing, nill))
            continue
        reused = []
        if reuse_root:
            r = os.path.join(reuse_root, ob)
            pat = "SCI_RECONSTRUCTED_%sKMOS.*.fits" % (prefix + "_" if prefix else "")
            dcr = sorted(glob.glob(os.path.join(r, "%s_KMOS.*.fits" % prefix))) if prefix not in ("", "dc") else []
            sky = sorted(glob.glob(os.path.join(r, "sky", pat)))
            hyb = sorted(glob.glob(os.path.join(r, "skysub_det", pat)))
            if not sky:
                problems.append("%s: no sky/%s in %s (run not finished?)" % (ob, pat, r))
                continue
            if len(sky) != len(dc):
                print("NOTE %s: %d dc frames but %d sky/ cubes in %s; the new run uses the same %d"
                      % (ob, len(dc), len(sky), r, len(sky)))
            reused = [(f, d) for f in dcr] + [(f, os.path.join(d, "sky")) for f in sky] + \
                     [(f, os.path.join(d, "skysub_det")) for f in hyb]
        if not dry:
            for sub in ["calib"] + (["sky", "skysub_det"] if reuse_root else []):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
        for f in cal:
            counts[link(f, os.path.join(d, "calib", os.path.basename(f)), dry, problems)] += 1
        for f in dc:
            counts[link(f, os.path.join(d, os.path.basename(f)), dry, problems)] += 1
        for f, where in reused:
            counts[link(f, os.path.join(where, os.path.basename(f)), dry, problems)] += 1
        print("%s %s: %d dc frames, %d calibration files%s" % (
            ob, bands[0], len(dc), len(cal), "" if not reuse_root else
            "; reused %d dcr frames, %d sky/ and %d skysub_det/ cubes" % (
                sum(where == d for _, where in reused), sum(where.endswith("sky") for _, where in reused),
                sum(where.endswith("skysub_det") for _, where in reused))))
    print("%s: %d OBs -> %s; links new %d, kept %d, conflicts %d%s" % (
        cfg, len(obs), dst_root, counts["new"], counts["kept"], counts["conflict"], " (dry run, nothing created)" if dry else ""))


def status(cfg):
    bands, dst_root, obs, cubes = read_cfg(cfg)
    band = bands[0]
    want = cubes.get("lcorr_degree")
    print("%s (%s, lcorr_degree %s): %s" % (cfg, band, want, dst_root))
    print("%-30s %3s %4s %4s %4s %4s %-20s %-28s %s" % ("OB", "dc", "dcr", "sky", "det", "rv", "wing", "OH_SPEC",
                                                       "degree in log"))
    tot = Counter()
    for ob in obs:
        d = os.path.join(dst_root, ob)
        ndc = len(glob.glob(os.path.join(d, "dc_KMOS.*.fits")))
        ndcr = len(glob.glob(os.path.join(d, "dcr_KMOS.*.fits")))
        skyc = sorted(glob.glob(os.path.join(d, "sky", "SCI_RECONSTRUCTED_*.fits")))
        nsky = len(skyc)
        ndet = len(glob.glob(os.path.join(d, "skysub_det", "SCI_RECONSTRUCTED_*.fits")))
        rv = sorted(glob.glob(os.path.join(d, "science", "SCI_RECONSTRUCTED_*_s_i_bkg2d_rv.fits")))
        nrv = len(rv)
        wings = Counter()
        for f in rv:
            try:
                h = fits.getheader(f, 0)
            except Exception:
                continue
            if h.get("HIERARCH ESO PYSPARK SKY HYBRID"):
                wings[(h.get("HIERARCH ESO PYSPARK SKY HYBRID WING", 0),
                       h.get("HIERARCH ESO PYSPARK SKY HYBRID OHMODEL", "median"),
                       h.get("HIERARCH ESO PYSPARK SKY HYBRID FREEMODEL", "tmpl"))] += 1
        wing = ",".join("%dxW%d%s" % (n, w, "" if (om, fm) == ("median", "tmpl") else "/%s/%s" % (om, fm))
                        for (w, om, fm), n in sorted(wings.items())) or "-"
        # a --reuse procdir has no reconstruction of its own: OH_SPEC and degree of the run that made the sky/ cubes
        rec, mark = d, ""
        if skyc and os.path.islink(skyc[0]):
            rec, mark = os.path.dirname(os.path.dirname(os.path.realpath(skyc[0]))), "*"
        oh = "-"
        for sof in sorted(glob.glob(os.path.join(rec, "sof", "*_cubes_*.sof"))):
            with open(sof) as f:
                oh = ([os.path.basename(l.split()[0]) for l in f if l.strip().endswith("OH_SPEC")] or ["none"])[0]
        deg = set()
        for log in glob.glob(os.path.join(rec, "*_%s_*.log" % band)):     # reduce_kmos.py log <user>_<band>_<date>.log
            with open(log, errors="replace") as f:
                deg |= set(re.findall(r"OH lambda correction: polynomial degree (\d+)", f.read()))
        degs = ",".join(sorted(deg)) or "-"
        if mark:                                 # --reuse: started once the new run has written its OB log or cubes
            started = bool(glob.glob(os.path.join(d, "*_%s_*.log" % band))) or nrv > 0
        else:                                    # no degree line and no sky/ cubes: not reduced yet
            started = bool(deg) or nsky > 0
        flag = "" if (not started or want is None or degs == str(want)) else "   <-- expected %s" % want
        print("%-30s %3d %4d %4d %4d %4d %-20s %-28s %s%s" % (ob, ndc, ndcr, nsky, ndet, nrv, wing, oh + mark,
                                                             degs + mark if deg else degs, flag))
        tot["ob"] += 1; tot["started"] += started; tot["sky"] += nsky > 0; tot["rv"] += nrv > 0; tot["det"] += ndet
        tot["bad_degree"] += bool(flag)
    print("%d OBs: %d started, %d with sky/ cubes, %d with final cubes, %d started with a wrong or missing degree; "
          "%d skysub_det cubes" % (tot["ob"], tot["started"], tot["sky"], tot["rv"], tot["bad_degree"], tot["det"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cfg", nargs="+", help="pyspark config(s); their [pathing] procdir is the new PROCDIR")
    ap.add_argument("--src", default="/local_storage/KOLJA/PROCDIR", help="PROCDIR with dc_ frames and calib/ (read only)")
    ap.add_argument("--reuse", help="finished procdir whose dcr_ frames, sky/ and skysub_det/ cubes are linked as well "
                                    "(re-run from the sky subtraction on, read only)")
    ap.add_argument("--dry-run", action="store_true", help="check inputs, create nothing")
    ap.add_argument("--status", action="store_true", help="report the progress of the re-reduction instead")
    a = ap.parse_args()
    if a.status:
        for cfg in a.cfg:
            status(cfg)
        return 0
    problems = []
    for cfg in a.cfg:
        prepare(cfg, os.path.realpath(a.src), a.dry_run, problems,
                reuse_root=os.path.realpath(a.reuse) if a.reuse else None)
    for p in problems:
        print("PROBLEM:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
