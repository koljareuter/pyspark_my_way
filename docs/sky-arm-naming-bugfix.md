# Sky-subtraction: ARM4_SKY / ARM5_SKY not recognized as sky frames

**Date:** 2026-09-30  
**Affected file:** `pyspark/utils/fileutils.py`, function `get_obs_exposures`  
**Affected stage:** sky subtraction (`subtract_sky` → `match_obj_sky`)  
**Symptom:** 18 OBs fail at sky subtraction with:

```
sky subtraction raised an exception (likely no matching sky exposures for <OB>)
Traceback: ... match_obj_sky ... np.argmin(np.fabs(np.array(matched_times) - ...))
ValueError: zero-size array to reduction operation minimum which has no identity
```

All 18 have complete calibrations (5/5 key files) and fully reconstructed
cubes in `sky/` — only sky subtraction fails.

---

## Background

The KMOS observing strategy for these OBs uses a dedicated-sky pattern:
alternating object (`O`) and sky (`S`) exposures. kmos_sci_red reconstructs
one cube per exposure and `spark_reconstruct.py` moves all of them to the OB's
`sky/` subdirectory. `subtract_sky` then calls `get_obs_exposures` to sort
those cubes into object and sky lists before running ZAP.

The OBJ/SKY classification in `get_obs_exposures` reads the FITS header
keywords `OCS.ARM1.NAME`, `OCS.ARM2.NAME`, `OCS.ARM3.NAME` (via
`dfits | fitsort`) and checks whether any of them equals one of a hard-coded
list:

```python
if all([True if ent not in ['ARM1_SKY','ARM2_SKY','ARM3_SKY'] else False
        for ent in sline[-3:]]) and sline[2] == obsband:
    exposures.append((..., 'OBJ', ...))
elif sline[2] == obsband:
    exposures.append((..., 'SKY', ...))
```

The failing OBs were observed with sky IFU arms named **`ARM4_SKY`** and
**`ARM5_SKY`** (confirmed via `dfits | fitsort ocs.arm1.name ocs.arm2.name
ocs.arm3.name` on the sky/ cubes). These names are not in the hard-coded list,
so every exposure — including the true sky frames — is classified as `OBJ`.
`sky_exposures` is empty, `match_obj_sky` tries to call `np.argmin` on an
empty array, and raises `ValueError`.

Example fitsort output for `KMOS.2026-01-14T05:23:13.608`:

```
FILE                                                    OCS.ARM1.NAME   OCS.ARM2.NAME   OCS.ARM3.NAME
SCI_RECONSTRUCTED_dc_KMOS.2026-01-14T05:23:13.608.fits                  3671            3662   ← OBJ
SCI_RECONSTRUCTED_dc_KMOS.2026-01-14T05:28:50.550.fits                  ARM4_SKY        ARM5_SKY  ← SKY (missed)
SCI_RECONSTRUCTED_dc_KMOS.2026-01-14T05:34:43.375.fits                  3671            3662   ← OBJ
...
```

---

## Root cause

`get_obs_exposures` hard-codes exactly three arm-name strings
(`ARM1_SKY`, `ARM2_SKY`, `ARM3_SKY`). KMOS OBs can allocate sky arms at any
IFU number; any OB using arms 4, 5, 6 … as sky positions is silently broken.

---

## Fix

**`pyspark/utils/fileutils.py`** — replace the fixed-list membership test with
an `endswith('_SKY')` check that matches any arm name in the conventional
`ARM<N>_SKY` pattern:

```python
# before
if all([True if ent not in ['ARM1_SKY','ARM2_SKY','ARM3_SKY'] else False
        for ent in sline[-3:]]) and sline[2] == obsband:

# after
if all([True if not ent.endswith('_SKY') else False
        for ent in sline[-3:]]) and sline[2] == obsband:
```

The `elif` branch that classifies the frame as SKY is unchanged — it triggers
whenever the `if` condition is false and the band matches.

---

## Deployment

The fix must be applied to the **server's running pyspark installation**
(not just the local venv), since the pipeline runs under
`/home/kreuter/miniforge3/envs/pyspark/`:

The fix lives in the project's own `pyspark/` package, which Python loads in preference
to the conda-env installation.  Commit the change and pull on the server:

```bash
# on hqempower1, from ~/pyspark_my_way
git pull
python3 -c "import pyspark.utils.fileutils as f; import inspect; src=inspect.getsource(f.get_obs_exposures); print('PATCHED' if 'endswith' in src else 'NOT PATCHED')"
```

Or patch the file directly without a git pull:

```bash
# on hqempower1
python3 -c "
f='/home/kreuter/pyspark_my_way/pyspark/utils/fileutils.py'
with open(f) as g: txt=g.read()
old=\"ent not in ['ARM1_SKY','ARM2_SKY','ARM3_SKY']\"
new=\"not ent.endswith('_SKY')\"
assert old in txt, 'Pattern not found -- check the file'
txt=txt.replace(old, new)
with open(f,'w') as g: g.write(txt)
print('Patched OK')
"
```

After patching, re-run the pipeline for the 18 previously failing OBs.

---

## Note: why only 18 OBs, not all 97?

The 79 OBs that already have `_rv.fits` in `science/` were reduced in a
**previous run** (before September 2026). `check_progress.py` sees those
pre-existing files and counts them as done. In the September run, sky
subtraction failed for all OBs — but since the 79 already had output, only
the 18 new OBs (never successfully reduced before) were left without results.
