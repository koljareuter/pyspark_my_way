# ZAP sky subtraction degrades IZ: bundled IZ eigenbasis has a wrong mean spectrum

**Date found:** 2026-10-07
**Affected stage:** sky subtraction, `pyspark/sky/sky_subtraction.py::zap_sky_v2` (`sky_method = zapv2`)
**Affected data:** every IZ OB (all used `pyspark/data/ZAPv2_SVD_IZ_n500004.fits`, confirmed in the reduction logs).
YJ and H are **not** affected by this mechanism.

---

## Symptom

[OIII]/Hβ in the IZ cubes are not detected even in rest-frame stacks of 24 galaxies;
IZ continuum S/N is ~3× lower than YJ for the same galaxies; IZ spectra show strong
structure common to all spaxels.

## Cause

`zap_sky_v2` projects the (object − sky) residual onto an *external* eigenbasis.
With `whiten=False` the cleaned spectrum is `(I − P)(X − mean_)`: the basis' stored
mean spectrum (`HDU1`, `mean_`) is subtracted from every spaxel.

| basis file | dtype | `mean_` range | `explained_variance_` max |
|---|---|---|---|
| ZAPv2_SVD_YJ_n500004 | float64 | −0.04 … +0.04 | 0.82 |
| ZAPv2_SVD_H_n500004  | float64 | −0.004 … +0.006 | 0.043 |
| ZAPv2_SVD_K_n500004  | float64 | −0.007 … +0.025 | 0.006 |
| **ZAPv2_SVD_IZ_n500004** | **float32** | **−104 … +31** | **1.7e5** |

The IZ file was produced differently (float32, `ZAPMASK=None`) and its mean spectrum is
~10³× larger than in the other bands, i.e. on a different normalisation than the data.
For YJ/H/K the `mean_` term is negligible; for IZ it injects a spurious coherent spectrum.

## Evidence

Re-running `zap_sky_v2` locally on `kmos_ob_proc_server_BACKUP_MEETINGS` (un-sky-subtracted
frames in `sky/`), identical code, reproduces the pipeline `_s.fits` to 1e-4 σ.

1. **Background statistics, one object/sky pair, all 24 IFUs** (coherent = rms of the median
   background spectrum):

   | | IZ coherent residual vs. plain (obj − sky) | YJ |
   |---|---|---|
   | pipeline (bundled basis) | **×1.5 worse** | ×0.85 (better) |
   | bundled basis with `mean_ = 0` | ×0.97 | identical to pipeline |

2. **Hα of CDFS 23420 (z = 0.605, Hα in IZ), 26 single 300 s frames, model-independent S/N:**
   median S/N pipeline **6.2**, plain obj − sky 5.8, ZAP with `mean_ = 0` **9.5**.
   (One frame, Gaussian fit: 9.5 vs 26.1; the line flux is preserved, the error triples.)

3. **All combined cubes** (`kmos_proc_and_comb_server`), coherent background residual /
   pure-noise expectation: **IZ 11.0**, YJ 3.8, H 4.3.

A per-OB basis built from raw sky frames (`build_custom_svd`) does **not** help: its
`mean_` is the mean *sky emission*, not a residual, so it has the same problem.

## Fix — prepared 2026-10-07, to be run on hqempower1

Prepared files (also applied to **H**, at the user's request — its mean is already tiny,
so only a small change is expected there; YJ is left unchanged):

| file | purpose |
|---|---|
| `pyspark/data/ZAPv2_SVD_IZ_n500004_nomean.fits` | IZ basis, HDU1 (mean) = 0, all else bit-identical |
| `pyspark/data/ZAPv2_SVD_H_n500004_nomean.fits`  | H basis, HDU1 (mean) = 0, all else bit-identical |
| `scripts/kmos_server_resky_IZ.cfg` | 76 IZ OB folders, sky subtraction onward, `svd_path` = IZ nomean basis |
| `scripts/kmos_server_resky_IZ_rebuild.cfg` | 7 IZ OBs (2026-01-13/14) whose un-sky-subtracted cubes are no longer in `sky/`: detector correction + cube reconstruction from raw, then ZAP with the nomean basis |
| `scripts/kmos_server_resky_H.cfg`  | 18 H OB folders, sky subtraction onward, `svd_path` = H nomean basis |
| `scripts/check_resky_ready.py` | preflight: basis, `sky/` cubes (or XCAL/YCAL/LCAL if `reconstruct_cubes = True`), calib files, raw frames; `--backup` |

Verified locally: the IZ nomean file reproduces the tested variant exactly (Hα S/N of
CDFS 23420 in the test frame 26.1 vs 9.5 with the bundled basis); `set_config` parses
both cfgs as intended (one band each, calibration/detector/reconstruction off, prefix `dc`).

One cfg per band on purpose: `svd_path` applies to every band of a run.

```bash
# on hqempower1, from ~/pyspark_my_way  (after the files are committed and pushed)
git pull
python scripts/check_resky_ready.py scripts/kmos_server_resky_IZ.cfg --backup
python scripts/check_resky_ready.py scripts/kmos_server_resky_IZ_rebuild.cfg --backup
python scripts/check_resky_ready.py scripts/kmos_server_resky_H.cfg  --backup
nohup env PYTHONPATH=. python scripts/reduce_kmos.py -p scripts/kmos_server_resky_IZ.cfg > resky_IZ.out 2>&1 &
# when that has finished:
nohup env PYTHONPATH=. python scripts/reduce_kmos.py -p scripts/kmos_server_resky_IZ_rebuild.cfg > resky_IZ_rebuild.out 2>&1 &
# when that has finished:
nohup env PYTHONPATH=. python scripts/reduce_kmos.py -p scripts/kmos_server_resky_H.cfg  > resky_H.out 2>&1 &
```

Then re-combine the IZ and H OBs exactly as before.  New IZ/H OBs (e.g. September 2026)
must also be reduced with `svd_path` pointing to the nomean basis of their band.

Original recipe for reference — zero the mean spectrum of the basis:

```python
from astropy.io import fits
import numpy as np
with fits.open('pyspark/data/ZAPv2_SVD_IZ_n500004.fits') as h:
    h[1].data = np.zeros_like(h[1].data)
    h.writeto('pyspark/data/ZAPv2_SVD_IZ_n500004_nomean.fits')
```

and use it for IZ (e.g. `svd_path` in the `[cubes]` section of the cfg, which
`zap_sky_v2` already supports).  Then re-run, for IZ OBs only, sky subtraction and all
later steps (`_s` → `_s_i` → `_s_i_bkg2d` → `_rv`) from the existing `sky/` frames, and
re-combine.  Calibrations and cube reconstruction do not need to be redone.

Expected gain: ~1.5× Hα S/N per frame in IZ (≈ 2.3× effective exposure time).

## Checked and found OK (no band-specific error)

* Calibration files per band (XCAL/YCAL/LCAL/MASTER_FLAT, OH_SPEC, WAVE_BAND) in every
  combined cube header match the band; all kmos_sci_red / kmos_combine parameters identical
  across bands.
* `make_image` YJ defaults (1.05–1.35 µm, mask 1.25–1.29 µm) — callers always pass the
  band limits; the mask is inert outside YJ.
* `bkg_correct_2d` (constant per IFU), `velocity_correct` (barycentric only) — band-independent.
* Number of ZAP eigenspectra chosen: YJ 53, IZ 24 — not pathological.
* Spatial (spaxel-to-spaxel) noise vs. STAT: IZ 0.91, YJ 1.01, H 1.02.
* H: no un-sky-subtracted frames available locally for the direct test; its basis `mean_`
  is smaller than YJ's, and its coherent residual (4.3) matches YJ (3.8).  The poorer look
  of H spectra is consistent with the much brighter/denser H-band OH sky, not a reduction error.
