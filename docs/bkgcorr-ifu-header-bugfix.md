## Bug: `bkg_correct_2d` reads the wrong FITS extension's header

**File:** `pyspark/cubes/bkgcorr.py`
**Function:** `bkg_correct_2d`
**Symptom:** `KeyError: "Keyword 'NAXIS3' not found."` during the global background
correction step of `reduce_kmos.py`, raised from `reconstruct_wave()` in
`pyspark/utils/obsutils.py`.

### Root cause

The function loops over all 24 KMOS IFU extensions correctly:

```python
for ifu in range(1, 25):
    cube = np.array(ofile[ifu].data)
    ...
```

but then reconstructs the wavelength array from a **hardcoded extension index**
instead of the current loop variable:

```python
wave_array = reconstruct_wave(ofile[1].header)   # always extension 1
```

This works only as long as IFU 1 happens to hold a valid 3D reconstructed cube
for that exposure. Not every IFU is allocated to a target on every OB — when
IFU 1 itself is inactive/empty, its header has no `NAXIS3` (it isn't a cube),
and `reconstruct_wave()` raises `KeyError`. This aborts the *entire* exposure,
even for otherwise-valid IFUs later in the loop, since the traceback happens
before the current IFU's correction is applied.

### Fix

```diff
-                    wave_array = reconstruct_wave(ofile[1].header)
+                    wave_array = reconstruct_wave(ofile[ifu].header)
```

Read the wavelength solution from the header of the IFU actually being
processed, not a fixed extension.

### Where deployed

- Source: `pyspark/cubes/bkgcorr.py` (this repo)
- Runtime copy on `empower.hq.eso.org`: `~/miniforge3/envs/pyspark/lib/python3.9/site-packages/pyspark/cubes/bkgcorr.py`
  (installed as a standalone copy, not editable-linked back to `~/pyspark_my_way` —
  patched there directly via `sed` rather than through a package reinstall)

### Observed trigger

`SCI_RECONSTRUCTED_dc_KMOS.2026-01-14T08:44:59.864_s.fits`, IZ band,
night `2026-01-14`, during a 3-night (`2026-01-12`–`2026-01-14`) test run of
`reduce_kmos.py` against a night-sorted symlink tree built from
`kmos_nights.json` via `lodo_sorting_algrithm.py`.
