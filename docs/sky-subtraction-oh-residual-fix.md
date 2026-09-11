## Sky-subtraction: OH residual spikes dominating the "before" spectra

**File:** `pyspark/sky/sky_subtraction.py`
**Function:** `zap_sky_v2`, plus new `build_custom_svd`
**Symptom:** in `telluric_check.ipynb`'s per-IFU flux plots, the tallest
feature in almost every panel didn't sit at the target's real (redshift
cross-matched) emission line — it was somewhere else entirely, and looked
different in every panel.

### Root cause

Cross-matching every KMOS target against public redshift catalogs
(`scripts/fetch_target_redshifts.py` /`target_redshifts.csv`) and checking
where each panel's dominant spike actually falls showed the same wavelengths
recurring across *unrelated* IFUs, OBs, and nights — e.g. 1.289 μm, 1.051 μm,
and 1.187 μm each showed up as "the loudest feature" in 3+ independent
(OB, IFU) combinations with different targets and different individual
redshifts. A real emission line tracks its own target's redshift and would
land at a different wavelength each time; a fixed-wavelength artifact
wouldn't move at all. These don't move — diagnostic of OH night-sky emission
residuals, not astrophysical signal.

Traced to how `zap_sky_v2` and `subtract_sky` do sky subtraction:

1. `match_obj_sky()` pairs each object exposure with only the single
   *closest-in-time* sky exposure (same DIT) — no interpolation between
   bracketing sky frames, no combining several. OH airglow strength drifts
   on a few-minutes timescale, so this leaves a systematic mismatch before
   ZAP even runs.
2. ZAP's residual-line cleaning depends entirely on a fixed, pre-baked PCA
   eigenbasis (`ZAPv2_SVD_YJ_n500004.fits`, bundled as a static package
   resource) — built from *some other* dataset, not this one. A
   `make_SVDs_v2()` function already existed in this file to build a custom
   eigenbasis from real data, but nothing in the pipeline ever called it.
3. `clean=False` was passed to both ZAP calls inside `zap_sky_v2`, disabling
   ZAP's own pixel-level cleaning pass.

### Fix

- Added `build_custom_svd()`: pools every finite spaxel spectrum from this
  dataset's own SKY-type exposures (found via the existing
  `get_obs_exposures`) across all 24 IFUs, and hands the stack to the
  already-existing `make_SVDs_v2()` to write a dataset-specific eigenbasis
  FITS file.
- `zap_sky_v2` now takes an optional `svd_path` kwarg; if set, it's used
  instead of the bundled generic SVD file. `subtract_sky()` already forwards
  `**kwargs` straight through, so this only needs a `svd_path = ...` line in
  a config's `[cubes]` section — no signature changes needed upstream.
  Unset, behavior is unchanged (falls back to the bundled file).
- Flipped `clean=False` → `clean=True` in both ZAP calls in `zap_sky_v2`.

### How to use

```python
from pyspark.sky.sky_subtraction import build_custom_svd
from pyspark.utils import setuputils

logger = setuputils.set_logging(logbase='svd')
build_custom_svd(
    subdir="KMOS.2026-01-13T03:42:28.829",   # any OB folder for the target band
    inputdir="/home/kreuter/kmos_ob_sorted",
    obsband="YJ",
    pcaoutfits="/home/kreuter/kmos_ob_proc/ZAPv2_SVD_YJ_custom.fits",
    logger=logger,
)
```

then in the reduction config's `[cubes]` section:

```
svd_path = /home/kreuter/kmos_ob_proc/ZAPv2_SVD_YJ_custom.fits
```

### Update: tested on the server — not the actual cause

Built a custom YJ eigenbasis from `KMOS.2026-01-12T05:36:55.243`'s own sky
exposures (`build_custom_svd`, see `scripts/test_custom_sky_svd.py`) and
re-ran sky subtraction with it. Result: essentially no change. Peak |flux|
per IFU before/after was within ~1-6% everywhere
(`scripts/compare_sky_svd.py`), and a targeted check at the three specific
"shared" wavelengths (1.289/1.051/1.187 μm) that motivated this whole
investigation showed tiny, noise-level values (0.01-0.03) at a representative
non-source spaxel, both before and after (`scripts/compare_sky_svd_at_lines.py`)
— nowhere near the 10-124 flux spikes actually dominating the flux plots. A
genuine sky-emission residual should be elevated across most of the field at
that wavelength (sky background is roughly uniform spatially); this isn't.

Traced the biggest offender directly (IFU 7, spaxel (10,5), exposure
`06:08:54.225`, aperture-summed peak 124.03): the raw single-spaxel spectrum
neighborhood is `[..., -6.58, 1.67, -0.17, 0.02, 56.36, 78.37, -13.22, 3.49,
-0.56, ...]` — a sharp rise/fall over ~4 pixels, not the single isolated
sample a hot pixel would leave, and not the field-wide bump a sky line would
leave. Checked whether it survives sky subtraction at all: pulled the same
spaxel/wavelength from the **pre**-sky-subtraction reconstructed cube
(`.../sky/SCI_RECONSTRUCTED_dc_..._06:08:54.225.fits`, no `_s` suffix) —
values there are 56.50/78.47/-13.22, essentially identical. **The spike is
already present before any sky subtraction happens.** It cannot be a sky
artifact.

Checked the same spaxel across every other exposure of this OB: 2.2-5.8
(normal noise) in 8 of 9 exposures, 78.5 in exactly one. Not a chronically
bad detector pixel (which would repeat every exposure) — a one-off transient
hit, i.e. a cosmic ray, in that single exposure.

**Conclusion:** the custom-SVD / `clean=True` sky-subtraction changes above
are a real improvement to keep (marginal gains elsewhere, no downside seen),
but they are not the fix for the dominant spikes seen in
`telluric_check.ipynb`'s flux plots. Those are cosmic rays, and there is
currently no cosmic-ray or bad-pixel rejection anywhere in this pipeline
(`pyspark/detector/` only does bias/reference-channel/interchannel
corrections). That's the actual next thing to build.

### Not yet done

- `match_obj_sky()`'s single-nearest-sky-frame matching (point 1 above)
  hasn't been touched. Lower priority now that sky subtraction isn't the
  dominant issue, but still worth doing at some point.
- **Cosmic-ray / per-exposure outlier rejection** — the real fix for the
  spikes in the flux plots. Nothing in this codebase does this yet.
