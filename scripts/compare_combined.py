#!/usr/bin/env python
"""Old vs. new combined cubes, galaxy by galaxy (e.g. before/after the ZAP basis fix,
docs/zap-iz-svd-mean-bug.md).

For every (target list, target id, band) present in both directories:
  * line S/N (scripts/line_snr.py, Gaussian fit, correlated noise) of Hα, Hβ,
    [OIII]5007 where covered — same aperture in old and new (centred on the
    brightest continuum spaxel of the OLD cube, circular r = 2 px)
  * S/N per sqrt(on-source hour), since the new combination may contain fewer frames
  * coherent background residual: robust rms of the median background spectrum
    divided by the expected noise of that median (≈1 for pure noise)

Usage:
    python scripts/compare_combined.py OLD_DIR NEW_DIR [--band IZ] [--csv out.csv]
"""
import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmos_sample as ks                                   # noqa: E402
from line_snr import measure_line_snr, noise_from_cube    # noqa: E402

warnings.filterwarnings("ignore")
LINES = {"Ha": 0.65628, "Hb": 0.48627, "OIII": 0.50082}
MAD = 1.4826


def load(path):
    with fits.open(path) as h:
        ext = 1 if len(h) > 1 and h[1].data is not None else 0
        cube = np.asarray(h[ext].data, float)
        wave = ks.wavelength_axis(h[ext].header)
    return wave, cube


def centre(cube, margin=2):
    img = np.nanmedian(cube, axis=0)
    ny, nx = img.shape
    inner = np.full_like(img, np.nan)
    inner[margin:ny - margin, margin:nx - margin] = img[margin:ny - margin, margin:nx - margin]
    if not np.isfinite(inner).any():
        return ny // 2, nx // 2
    return np.unravel_index(np.nanargmax(inner), img.shape)


def coherent(cube, ap):
    bg = cube[:, ~ap & (np.isfinite(cube).mean(0) >= 0.9)]
    if bg.shape[1] < 8:
        return np.nan
    med = np.nanmedian(bg, 1)
    sd = MAD * np.nanmedian(np.abs(bg - med[:, None]), 1)
    core = slice(100, len(med) - 100)
    d = np.diff(med[core]) / np.sqrt(2)                       # remove continuum level
    expect = 1.2533 * np.nanmedian(sd[core]) / np.sqrt(bg.shape[1])
    return MAD * np.nanmedian(np.abs(d - np.nanmedian(d))) / expect


def measure(path, z, cy, cx, r=2):
    wave, cube = load(path)
    yy, xx = np.mgrid[:cube.shape[1], :cube.shape[2]]
    ap = np.hypot(yy - cy, xx - cx) <= r
    spec = np.nansum(cube[:, ap], axis=1)
    spec[~np.isfinite(cube[:, ap]).any(1)] = np.nan
    noise = noise_from_cube(cube, ap)
    out = {"coh": coherent(cube, ap)}
    for name, lam in LINES.items():
        res = measure_line_snr(wave, spec, z, lam, noise=noise)
        bad = {"no_coverage", "too_few_pixels", "fit_failed"} & set(res.flags)
        out[name] = np.nan if bad else float(res.snr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--band", default="IZ")
    ap.add_argument("--csv")
    a = ap.parse_args()

    for d in (a.old, a.new):
        if not os.path.isdir(d):
            sys.exit("ERROR: directory not found: %s" % os.path.abspath(d))
    zcat = ks.RedshiftCatalog(os.path.join(os.path.dirname(os.path.abspath(__file__)), "target_redshifts.csv"))
    def index(base):
        withz, _ = ks.scan_combined_cubes(base, zcat)
        d = {}
        for i in withz:
            if i["band"] == a.band and "_SKY" not in os.path.basename(i["path"]):
                d.setdefault((i["target_list"], i["tid"]), i)
        return d
    old, new = index(a.old), index(a.new)
    keys = sorted(set(old) & set(new))
    print("%s: %d galaxies old, %d new, %d in both" % (a.band, len(old), len(new), len(keys)))
    if not keys:
        sys.exit("nothing to compare: no %s DATACUBE_COMBINED_*.fits with a redshift found in both "
                 "directories" % a.band)

    rows = []
    for k in keys:
        o, n = old[k], new[k]
        cy, cx = centre(load(o["path"])[1])
        mo, mn = measure(o["path"], o["z"], cy, cx), measure(n["path"], n["z"], cy, cx)
        row = dict(target_list=k[0], tid=k[1], z=o["z"], t_old_h=o["exptime_s"] / 3600,
                   t_new_h=n["exptime_s"] / 3600, coh_old=mo["coh"], coh_new=mn["coh"])
        for L in LINES:
            row[L + "_old"], row[L + "_new"] = mo[L], mn[L]
        rows.append(row)

    hdr = "%-26s %8s %5s %9s %11s" % ("galaxy", "tid", "z", "t old/new", "coh old/new")
    hdr += "".join(" %15s" % ("S/N %s o->n" % L) for L in LINES)
    print(hdr)
    for r in rows:
        s = "%-26s %8s %5.3f %4.1f/%4.1f %5.1f/%5.1f" % (r["target_list"][:26], r["tid"], r["z"],
                                                          r["t_old_h"], r["t_new_h"], r["coh_old"], r["coh_new"])
        for L in LINES:
            s += " %7.1f->%6.1f" % (r[L + "_old"], r[L + "_new"]) if np.isfinite(r[L + "_old"]) else " %15s" % "-"
        print(s)

    print("\nsummary (median over galaxies; 'per sqrt(t)' corrects for different exposure)")
    co = np.array([[r["coh_old"], r["coh_new"]] for r in rows])
    print("  coherent background residual   old %.2f   new %.2f" % tuple(np.nanmedian(co, 0)))
    for L in LINES:
        p = np.array([[r[L + "_old"], r[L + "_new"], r["t_old_h"], r["t_new_h"]] for r in rows
                      if np.isfinite(r[L + "_old"]) and np.isfinite(r[L + "_new"])])
        if not len(p):
            continue
        ratio = p[:, 1] / p[:, 0]
        ratio_t = ratio * np.sqrt(p[:, 2] / p[:, 3])
        strong = p[:, 0] >= 5
        print("  %-5s N=%3d  S/N new/old %.2f  per sqrt(t) %.2f   (old S/N >= 5: N=%d, %.2f / %.2f)"
              % (L, len(p), np.median(ratio), np.median(ratio_t), strong.sum(),
                 np.median(ratio[strong]) if strong.any() else np.nan,
                 np.median(ratio_t[strong]) if strong.any() else np.nan))
    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print("written", a.csv)


if __name__ == "__main__":
    main()
