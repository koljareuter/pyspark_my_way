"""Hybrid sky subtraction for OBs whose SKY frames were taken at another rotator angle (docs/hybrid-sky.md).

pyspark subtracts the sky in cube space: O and S are reconstructed separately, each with its own calibrations
and its own per-frame OH wavelength correction. When the SKY frame's rotator angle differs from the OBJECT's
(37-90 deg for cosmos_0.6_0.8_h21.4/h21.7, cosmos_target_legac, goods_0.6_0.8_h21.6_lir), the same detector
pixels (warm/hot/RTS pixels) land in different voxels of O and S and do not cancel: 3-5 % outlier voxels and a
pattern shared by all O - S frames of a target, which does not average down in the combination. These lists
are 2-3x noisier per sqrt(t) than the matched ones.

The DRP's own detector-level subtraction (kmos_sci_red -no_subtract=FALSE) maps S with the OBJECT's
calibrations and lambda correction, so detector features cancel, but it leaves the OH-line misalignment between
the two frames (flexure) that the separate reconstruction corrects. The hybrid keeps both:
    R_ds = detector-level O - S (DRP)
    D    = (O - S) - R_ds = recon_O(S) - recon_S(S)       (linear reconstruction: S frame only, no target signal)
    R_hyb = R_ds + model(D),   model per slitlet: in OH-line planes the slitlet median of D, elsewhere a fit of
            the slitlet-median of D with (a0 + a1 x) g + (b0 + b1 x) S + (c0 + c1 x) g',  g = dS/dlambda,
            S = IFU median of the SKY cube, x = (lambda - 0.93)/0.14, all high-passed (101 planes)
Every term of the model comes from the SKY frame; R_hyb then goes through ZAP exactly like O - S.

Validated on IZ (2026-10-23): 6-frame Hbeta-window stacks noise x0.44-0.67 on two mismatched OBs, matched OB
unchanged (so use it only for mismatched OBs); injected detector lines keep their flux exactly when the
object's lambda correction is frozen (906 compact + 241 extended/edge sources, ratio 1.00000). Not for H (the
IZ model makes H 1.3-3.3x worse); YJ untested.
"""
from __future__ import print_function

import numpy as np
import bottleneck as bn

from ..utils.obsutils import find_lines

__all__ = ["hybrid_ifu", "highpass"]


def highpass(spec, width=101):
    """spec minus its centred running median (NaN-aware, at least width/4 valid values)."""
    s = np.asarray(spec, float)[:, None, None]
    nz, half = s.shape[0], width // 2
    pad = np.concatenate([np.repeat(s[:1], half, 0) * np.nan, s, np.repeat(s[-1:], half, 0) * np.nan], 0)
    med = bn.move_median(pad, width, axis=0, min_count=width // 4)[width - 1:width - 1 + nz]
    return (s - med)[:, 0, 0]


def hybrid_ifu(ocube, scube, rds, ifu, lam):
    """Hybrid sky-subtracted cube of one IFU.

    ocube, scube : un-subtracted OBJECT and nearest SKY cubes (pyspark sky/ reconstructions)
    rds          : the DRP detector-level sky-subtracted cube of the same OBJECT frame and IFU
    ifu          : IFU number (slitlets run along y for IFUs 1-16, along x for 17-24)
    lam          : wavelength of every plane (micron)
    Returns a float64 cube; NaN where rds is NaN.
    """
    oc = np.asarray(ocube, float)
    sc = np.asarray(scube, float)
    rds = np.asarray(rds, float)
    D = (oc - sc) - rds
    nz, ny, nx = D.shape
    x = (np.asarray(lam, float) - 0.93) / 0.14
    smed = np.nanmedian(sc.reshape(nz, -1), axis=1)
    smed = np.where(np.isfinite(smed), smed, 0.0)
    g = np.gradient(smed)
    g2 = np.gradient(g)
    T = np.column_stack([highpass(g), highpass(g) * x, highpass(smed), highpass(smed) * x,
                         highpass(g2), highpass(g2) * x])
    T = np.where(np.isfinite(T), T, 0.0)
    oh = find_lines(np.nanmedian(sc[:, 1:13, 1:13].reshape(nz, -1), axis=1)).astype(bool)
    oh = oh | np.roll(oh, 1) | np.roll(oh, -1)
    if ifu <= 16:
        groups = [(slice(None), slice(i, i + 1)) for i in range(nx)]
    else:
        groups = [(slice(i, i + 1), slice(None)) for i in range(ny)]
    model = np.zeros_like(D)
    for gy, gx in groups:
        sub = D[:, gy, gx].reshape(nz, -1)
        if np.isfinite(sub).sum() < nz:
            continue
        dmed = highpass(np.nanmedian(sub, axis=1))
        ok = np.isfinite(dmed)
        if ok.sum() < 200:
            continue
        coef = np.linalg.lstsq(T[ok], dmed[ok], rcond=None)[0]
        m = np.where(oh, np.nan_to_num(np.nanmedian(sub, axis=1)), T @ coef)
        model[:, gy, gx] = m[:, None, None]
    hyb = rds + model
    hyb[~np.isfinite(rds)] = np.nan
    return hyb
