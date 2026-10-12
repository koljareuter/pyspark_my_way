#!/usr/bin/env python
"""Old vs. new combined spectra, side by side, one page per galaxy
(companion to scripts/compare_combined.py; same aperture: circular r = 2 px around
the brightest continuum spaxel of the OLD cube).

Each page: full spectrum (old above, new below, same flux scale) and zoom panels
around Hα, Hβ, [OIII]5007 where covered, with the line S/N of old and new.
Pages are sorted by the largest new line S/N.  The first N_SUMMARY galaxies are
also drawn into one overview PNG.

Usage:
    python scripts/plot_compare_spectra.py OLD_DIR NEW_DIR [--band IZ] [--out kmos_compare_IZ]
        → kmos_compare_IZ.pdf, kmos_compare_IZ_top.png
"""
import argparse
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmos_sample as ks                                   # noqa: E402
from line_snr import measure_line_snr, noise_from_cube    # noqa: E402
from compare_combined import load, centre                  # noqa: E402

warnings.filterwarnings("ignore")
LINES = {"Hα": 0.65628, "Hβ": 0.48627, "[OIII]5007": 0.50082}
EXTRA = {"[OIII]4959": 0.49589, "[NII]6583": 0.65835, "[NII]6548": 0.65480}
C_OLD, C_NEW = "#9a9a9a", "#1f5fa8"
ZOOM_KMS = 2500.0


def extract(path, cy, cx, r=2):
    wave, cube = load(path)
    yy, xx = np.mgrid[:cube.shape[1], :cube.shape[2]]
    ap = np.hypot(yy - cy, xx - cx) <= r
    spec = np.nansum(cube[:, ap], axis=1)
    spec[~np.isfinite(cube[:, ap]).any(1)] = np.nan
    return wave, spec, noise_from_cube(cube, ap)


def snrs(wave, spec, noise, z):
    out = {}
    for name, lam in LINES.items():
        r = measure_line_snr(wave, spec, z, lam, noise=noise)
        out[name] = np.nan if {"no_coverage", "too_few_pixels", "fit_failed"} & set(r.flags) else r.snr
    return out


def smooth(y, n=3):
    k = np.ones(n) / n
    m = np.isfinite(y)
    num = np.convolve(np.where(m, y, 0), k, "same")
    den = np.convolve(m.astype(float), k, "same")
    return np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)


def robust_ylim(*ys, lo=1, hi=99.5, pad=0.15):
    v = np.concatenate([y[np.isfinite(y)] for y in ys if np.isfinite(y).any()])
    if not v.size:
        return -1, 1
    a, b = np.percentile(v, [lo, hi])
    d = (b - a) or 1.0
    return a - pad * d, b + pad * d


def mark_lines(ax, z, wlo, whi):
    for name, lam in {**LINES, **EXTRA}.items():
        lo = lam * (1 + z)
        if wlo <= lo <= whi:
            ax.axvline(lo, color="#d62728" if name in LINES else "#e7a0a0", lw=0.8, ls="--", zorder=1)


def draw_galaxy(fig, gs_row, g, compact=False):
    """g: dict with name, z, t_old, t_new, old=(w,s,n), new=(w,s,n), snr_old, snr_new."""
    (wo, so, no), (wn, sn, nn) = g["old"], g["new"]
    covered = [L for L in LINES if np.isfinite(g["snr_old"][L]) or np.isfinite(g["snr_new"][L])]
    ncol = 1 + len(covered)
    sub = gs_row.subgridspec(2, ncol, width_ratios=[2.6] + [1] * len(covered), hspace=0.06, wspace=0.18)

    # full spectrum: same y scale for old and new
    so_s, sn_s = smooth(so), smooth(sn)
    yl = robust_ylim(so_s, sn_s)
    for row, (w, s, lab, col, t) in enumerate([(wo, so_s, "old", C_OLD, g["t_old"]),
                                               (wn, sn_s, "new", C_NEW, g["t_new"])]):
        ax = fig.add_subplot(sub[row, 0])
        ax.plot(w, s, color=col, lw=0.6)
        mark_lines(ax, g["z"], np.nanmin(w), np.nanmax(w))
        ax.set_ylim(*yl)
        ax.set_xlim(np.nanmin(wo), np.nanmax(wo))
        ax.text(0.01, 0.95, "%s  (%.1f h)" % (lab, t), transform=ax.transAxes, va="top", fontsize=8,
                color=col, fontweight="bold")
        if row == 0:
            ax.set_title("%s   z = %.3f" % (g["name"], g["z"]), fontsize=10, loc="left")
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("λ obs [µm]", fontsize=8)
        ax.tick_params(labelsize=7)

    # zoom panels
    for k, L in enumerate(covered):
        lam = LINES[L] * (1 + g["z"])
        half = lam * ZOOM_KMS / 2.998e5
        segs = []
        for (w, s, n) in (g["old"], g["new"]):
            m = (w > lam - half) & (w < lam + half)
            segs.append((w[m], s[m], n[m]))
        yl = robust_ylim(segs[0][1], segs[1][1], lo=0.5, hi=99.9)
        for row, ((w, s, n), col, key) in enumerate(zip(segs, (C_OLD, C_NEW), ("snr_old", "snr_new"))):
            ax = fig.add_subplot(sub[row, 1 + k])
            ax.fill_between(w, -n, n, color=col, alpha=0.18, lw=0, step="mid")
            ax.step(w, s, where="mid", color=col, lw=0.8)
            ax.axhline(0, color="k", lw=0.4)
            mark_lines(ax, g["z"], lam - half, lam + half)
            ax.set_xlim(lam - half, lam + half)
            ax.set_ylim(*yl)
            v = g[key][L]
            ax.text(0.03, 0.95, "S/N %s" % ("–" if not np.isfinite(v) else "%.1f" % v),
                    transform=ax.transAxes, va="top", fontsize=8, color=col, fontweight="bold")
            if row == 0:
                ax.set_title(L, fontsize=9)
                ax.tick_params(labelbottom=False)
            ax.tick_params(labelsize=6)
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--band", default="IZ")
    ap.add_argument("--out", default=None)
    ap.add_argument("--top", type=int, default=6, help="galaxies in the overview PNG")
    a = ap.parse_args()
    out = a.out or "kmos_compare_%s" % a.band

    zcat = ks.RedshiftCatalog(os.path.join(os.path.dirname(os.path.abspath(__file__)), "target_redshifts.csv"))
    def index(base):
        d = {}
        for i in ks.scan_combined_cubes(base, zcat)[0]:
            if i["band"] == a.band and "_SKY" not in os.path.basename(i["path"]):
                d.setdefault((i["target_list"], i["tid"]), i)
        return d
    old, new = index(a.old), index(a.new)
    keys = sorted(set(old) & set(new))
    print("%d galaxies in both" % len(keys))

    gals = []
    for k in keys:
        o, n = old[k], new[k]
        cy, cx = centre(load(o["path"])[1])
        go, gn = extract(o["path"], cy, cx), extract(n["path"], cy, cx)
        gals.append(dict(name="%s %s" % (ks.canonical_field(k[0]), k[1]), z=o["z"],
                         t_old=o["exptime_s"] / 3600, t_new=n["exptime_s"] / 3600,
                         old=go, new=gn, snr_old=snrs(*go, o["z"]), snr_new=snrs(*gn, n["z"])))
    best = lambda g: np.nanmax([*g["snr_new"].values(), *g["snr_old"].values(), -99])
    gals.sort(key=best, reverse=True)

    with PdfPages(out + ".pdf") as pdf:
        for g in gals:
            fig = plt.figure(figsize=(13, 4.2))
            draw_galaxy(fig, fig.add_gridspec(1, 1, left=0.05, right=0.99, top=0.9, bottom=0.13)[0], g)
            pdf.savefig(fig)
            plt.close(fig)
    print("written", out + ".pdf", "(%d pages, sorted by line S/N)" % len(gals))

    top = gals[:a.top]
    fig = plt.figure(figsize=(13, 3.6 * len(top)))
    gs = fig.add_gridspec(len(top), 1, left=0.05, right=0.99, top=0.97, bottom=0.04, hspace=0.42)
    for i, g in enumerate(top):
        draw_galaxy(fig, gs[i], g)
    fig.savefig(out + "_top.png", dpi=110)
    print("written", out + "_top.png")


if __name__ == "__main__":
    main()
