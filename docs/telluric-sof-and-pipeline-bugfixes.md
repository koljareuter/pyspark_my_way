## Telluric correction via `kmos_sci_red` SOF + pipeline bug fixes

### 1. Telluric correction: SOF-based approach

**Files:** `pyspark/cubes/spark_reconstruct.py`, `pyspark/cubes/__init__.py`,
`scripts/reduce_kmos.py`, `scripts/kmos_test.cfg`

The KMOS DRS recipe `kmos_sci_red` natively supports a `TELLURIC` entry in its
SOF file. When present, it reconstructs the standard star, derives the telluric
absorption and flux response, and applies both corrections to the science cubes
during reconstruction — the intended, pipeline-internal approach.

A previous attempt had implemented telluric correction as a post-processing
Python step (`pyspark/cubes/telluric.py`, dividing reconstructed cubes by a
per-IFU transmission spectrum). That file has been removed; it was not part of
the original pyspark codebase and would have double-corrected if left alongside
the SOF approach.

#### Change in `spark_reconstruct.py`

A new parameter `telluric_in_sof=False` is added to `reconstruct()`. When
`True`, the function searches for a band-matched TELLURIC file in the OB's raw
directory and appends it to the SOF:

```python
if telluric_in_sof:
    tel_matches = glob.glob(os.path.join(kwargs['inputdir'], subdir,
                            '*TELLURIC*{0}*'.format(obsband.upper())))
    if tel_matches:
        ofile.write('{0} TELLURIC\n'.format(tel_matches[0]))
```

The glob filters by band (`obsband.upper()`) so that multi-band runs (e.g.
`IZ,YJ,H`) do not accidentally pass the wrong band's standard star.

#### Config

Enable with `telluric_in_sof = True` in the `[cubes]` section of the config.
Do **not** combine with a separate post-processing telluric step.

---

### 2. Bug: wrong label passed to `break_check` after reconstruction failure

**File:** `scripts/reduce_kmos.py`  
**Symptom:** When `reconstruct()` returns a non-zero exit code, the error log
reads as if a *calibration* step failed, because the variable `comp` (the
completion label from `make_spark_calibs`) was reused:

```python
# before
if ecode:
    break_check(comp, master_logger)   # comp = calibration label
```

**Fix:**

```diff
-    break_check(comp, master_logger)
+    break_check('reconstruction', master_logger)
```

---

### 3. Bug: unconditional `logger.debug` call in `spark_reconstruct.py`

**File:** `pyspark/cubes/spark_reconstruct.py`  
**Function:** `reconstruct`  
**Symptom:** `AttributeError: 'NoneType' object has no attribute 'debug'` when
`logger=None` (the default). Every other logger call in the function is guarded
by `if logger is not None`, but the call inside the science-frame loop was not:

```python
# before
for frame, _ in frames:
    ofile.write('{0} SCIENCE\n'.format(frame))
    logger.debug(...)   # crashes if logger is None
```

**Fix:**

```diff
+    if logger is not None:
         logger.debug('-- Added {0} SCIENCE'.format(frame))
```

---

### 4. Config: missing `interchannel_correct` key

**File:** `scripts/kmos_test.cfg`  
**Symptom:** `KeyError: 'interchannel_correct'` at `reduce_kmos.py:104`, since
the key was absent from the `[detector]` section while the code reads it
unconditionally.

`fix_detector_levels()` defaults `interchannel_correct=False` (it is an
optional interslice bias correction, off by default). The config now makes that
explicit:

```ini
[detector]
channel_correct = True
picture_correct = True
interchannel_correct = False
```
