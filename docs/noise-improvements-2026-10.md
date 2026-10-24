# Noise improvements, October 2026: row common mode, sky calibration angle, integer-shift combination

**Date:** 2026-10-23
**Code:** `pyspark/detector/row_commonmode.py`, `pyspark/cubes/sky_calangle.py`, `pyspark/cubes/wcssnap.py`,
`scripts/reduce_kmos.py` (three optional steps), `pyspark/sky/sky_combine.py` (multi-sky, see below)
**All options are off by default; with the defaults the output is bit-identical to before (verified).**
The largest gain, the hybrid sky subtraction for OBs with SKY frames at another rotator angle, has its own
document: docs/hybrid-sky.md; the server run with the adopted options: docs/procdir-hyb-server-run.md.

Context: after the OH wavelength fix (docs/oh-lcorr-linear-rereduction.md) the remaining blue-end "forest" is
noise. Single sky-subtracted frames have a 5x5-aperture noise 1.6x (YJ), 2.1x (IZ), 2.0x (H) above what
independent spaxels give. Four candidate causes were measured on local copies of one IZ, one YJ and one H OB
with the production-identical chain (patched DRP, KMOS_LCORR_DEGREE=1 and new OH_SPEC for IZ, then pyspark
O - S + ZAPv2); reproductions of the production products were bit-identical before any change.

## 1. Row common-mode correction  `[detector] row_commonmode_correct`

**Cause.** After the reference-column correction each frame still carries a per-row offset common to all
readout channels: 0.25-0.43 sigma_pix (IZ, YJ), 0.19-0.27 (H), white along the row index, independent between
frames and detectors. A detector row is one wavelength across the IFUs of that detector, so in the cubes it is
an offset coherent along the slitlets.

**Fix.** Per detector: subtract the OB-median frame; c(y) = median of the residual over the unilluminated
pixels of row y (NaN in XCAL at all six angles, dilated 1 px along x); subtract c(y) from the row. Corrected
copies are written with a new prefix; inputs are never modified.

**Measured** (worst IZ OB, 6-frame stacks of 5x5 apertures): noise x0.94-0.95 overall, x0.89-0.95 for IFUs on
detectors 2 and 3, no gain on detector 1 (only 40-75 unilluminated pixels per row) and none in H. About +5 %
S/N in IZ final spectra (+7-10 % for detector 2/3 IFUs); YJ expected similar (same amplitude), untested in cubes.
End to end (three IZ OBs, combined cubes, on top of the hybrid where it applies): 3x3-aperture noise
x0.92-0.95 for IFUs 9-24, IFUs 1-8 bit-identical.

```ini
[detector]
# frames the reconstruction reads
prefix = dcr
row_commonmode_correct = True
# frames the correction reads (existing dc_ frames)
row_commonmode_input = dc
row_commonmode_dets = 2,3
row_commonmode_grow = 1
row_commonmode_minpix = 150
```
With detector corrections also on, fix_detector_levels writes `row_commonmode_input` and the row step writes
`prefix`. Needs `calib/XCAL_<band>.fits`. Keep it off for H.

## 2. Sky frames at the object's calibration angle  `[cubes] sky_cal_angle = object`

**Cause.** kmos_sci_red selects MASTER_FLAT, LCAL and XCAL/YCAL by each frame's own rotator angle (6 angles,
nearest, no interpolation for flat/LCAL). The SKY frames are often 30-90 deg away from their OBJECT frames
(object time with the sky in another 60-deg calibration bin: IZ 41 %, H 81 %, YJ 84 %; since Feb 2026 almost
all), so O - S mixes two flats and two slitlet geometries. The DRP's own on-detector sky subtraction uses the
object's calibrations for both.

**Fix.** Each SKY frame is reconstructed from a copy whose `ESO OCS ROT NAANGLE` is the circular mean of the
angles of the OBJECT frames it is subtracted from (pyspark's nearest-sky rule); the true angle is written back
afterwards and the angle used is recorded (`ESO PYSPARK SKY CALANGLE`). Only the O - S + ZAP route.

**Measured.** IZ worst OB: 5x5 aperture noise x0.967-0.974 per frame (detector 2 x0.88-0.93), 6-frame stacks
x0.96-0.98. H OB: x0.915-0.962 per frame, x0.89-0.95 in stacks. Programme estimate (weighted by affected
object time): IZ +1.2 % S/N, H +5-6.5 %, YJ +2.5-6 % (unmeasured). Flat interpolation between angles gives
nothing; in IZ the slitlet geometry matters, in H the flat.

Verified: pyspark `reconstruct()` with the option reproduces the prototype bit for bit (object frames equal
production, sky frames equal the angle-forced reconstruction), headers restored, temporary copies removed.

## 3. Integer-shift combination  `[cubes] snap_wcs = True`

**Cause.** kmos_combine aligns frames by WCS and applies sub-pixel shifts with a hard-coded bicubic spline.
The O frames are dithered by half a spaxel, so 4 of 6 frames are interpolated, and the spline spreads
single-frame outliers into many voxels below the k-sigma threshold.

**Fix.** Header-only, after the velocity correction: CRPIX1/2 of every IFU move by at most 0.5 px so that the
target (ESO OCS ARMi ALPHA/DELTA) lies on the grid phase; kmos_combine then shifts by integers only. Original
CRPIX kept in `PYSPARK CRPIX1/2 ORIG`, frame flagged `ESO PYSPARK WCS SNAPPED`; repeated calls do nothing.

**Measured** with esorex kmos_combine on 6-frame OBs: matched-filter line S/N +5.0 % (IZ), +4.7 % (YJ),
+2.1 % (H); PSF FWHM +1.4-5.6 % (registration error up to 0.1"). Needs re-checking on full production sets
(35-128 frames). The spatial-noise calibration c(N) and PSF assumptions must be redone afterwards.

**Not used for now** (end-to-end test, docs/hybrid-sky.md): the gain came from the bicubic ringing of the many
single-frame outliers of the classic O - S. After the hybrid sky subtraction (and on the matched cdfs OB) the
empirical noise of 3x3 aperture sums was 5-14 % *higher* with snapped WCS in two of three OBs (neutral in the
third), because the spline also smooths the noise; with dithers of 1.5 spaxels (cdfs) four of six frames are
registered 0.5 px off in both axes. On the other hand, the spline turned one single-frame spike at an IFU edge
into a stripe that the line fit took for an [OIII] line at S/N 93 (cdfs 22737); the snapped cube keeps it in
one spaxel. Decide on full production sets with injected sources before switching it on.

## Not adopted / analysis side

* **Multi-sky subtraction** (`[cubes] sky_combine`, `sky_scale`, docs in `pyspark/sky/sky_combine.py`): mean of
  the OB's sky frames lowers single-frame noise by 13-19 % after ZAP, but the standard O S O O S O O S O
  sequence already uses every sky frame twice in the combined cube: gain 1-6 % for aperture spectra
  (8-15 % per spaxel). Per-band sky rescaling makes residuals worse after ZAP. Kept as an option, off.
* **WCS-only velocity correction**: +0.5-1.5 % S/N, -1-2 % line width; kmos_combine ignores CRVAL3, so it
  would need a post-combination step. Not implemented.
* **Analysis** (EMPOWER_Analysis_and_presentation/scripts): the Hbeta/[OIII] rest-frame stacks were already
  weighted per pixel by 1/noise^2 (x3.1 S/N over exposure-time weights). New `restframe_stack.py` adds
  `weighting="pergal_lines"` for line ratios (per-pixel weights bias stacked ratios by +0.2..+9 %, up to
  +25 %), and `kmos_sample.load_galaxy_band_stacks(normalise="after")` avoids the 3 % (up to 34 %) S/N loss
  of normalising cubes before stacking two-cube galaxies.

Investigation scripts, logs and products: WSL `~/gains/{cmnoise,calib-angle,resampling}`; reports in the
analysis session.
