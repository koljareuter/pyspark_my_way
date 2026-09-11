## Implementing telluric correction in `reduce_kmos.py`

### Starting question

Whether telluric correction just needed to be *activated* (an existing but
disabled code path) rather than newly written.

**Answer: it needed to be written.** Confirmed by grepping the whole repo for
`get_instrument_response` and `extract_spectrum` (the two functions in
`pyspark/calibration/telluric_standard.py` that would derive/apply a telluric
correction) — both only appear in their own definitions, never called
anywhere else. `reduce_kmos.py` doesn't import either of them, and neither
function is actually finished (`extract_spectrum` is a no-op;
`get_instrument_response` calls `trace()`, which itself just builds up
wavelength chunks and returns without doing anything with them). Nothing in
the existing pipeline consumes the `TELLURIC_KMOS.<timestamp>.fits_<band>.fits`
files either.

### What made this tractable rather than a from-scratch derivation

Inspected the actual structure of one telluric file, from a folder in the
new OB-sorted tree (`kmos_ob_sorted/KMOS.2026-01-15T01:00:33.152/TELLURIC_KMOS.2026-01-15T01:00:33.152.fits_IZ.fits`):

- Per-IFU extensions named `IFU.<n>.DATA` and `IFU.<n>.NOISE` (n = 1-24),
  each a 1D array of 2048 points.
- Wavelength solution: `CRVAL1 = 0.779999971389771`, `CRPIX1 = 1.0`.

Compared against a real final reduced science product
(`kmos_test_proc/2026-01-13/science/SCI_RECONSTRUCTED_dc_KMOS.2026-01-13T02:54:58.606_s_i_bkg2d_rv.fits`):

- Per-IFU extension `IFU.<n>.DATA`, shape `(2048, 14, 14)` — wavelength ×
  14×14 spaxels.
- Wavelength solution: `CRVAL3 = 0.779999971389771`, `CDELT3 = 0.000151367217767984`,
  `CRPIX3 = 1.0`.

Same `CRVAL`, same axis length (2048) on the same band/night. This is an
already-computed, per-IFU transmission spectrum that lines up directly with
the science cube's own wavelength axis — so the actual work is "read it and
divide," not "derive a transmission curve from a standard star spectrum,"
which is what the unfinished `get_instrument_response`/`extract_spectrum`
stubs would have required.

Also relevant: `kmos_ob_structure.json` (as opposed to the original
`kmos_nights.json`) puts exactly one OB's worth of files, and exactly one
matching telluric file, in each sorted folder — this removes the ambiguity
that existed when a whole night (with several OBs, several telluric files)
was sorted into one folder.

### Files created / changed

- **New:** [pyspark/cubes/telluric.py](../pyspark/cubes/telluric.py) —
  `telluric_correct()`, following the same structure as the existing
  `bkg_correct_2d`/`velocity_correct` steps (locate exposures via
  `get_obs_exposures`, `duplicate()` each into a new suffixed file, open with
  `mode='update'`, loop over the 24 IFU extensions).
- **Changed:** [pyspark/cubes/__init__.py](../pyspark/cubes/__init__.py) —
  export `telluric_correct`.
- **Changed:** [scripts/reduce_kmos.py](../scripts/reduce_kmos.py) — import
  `telluric_correct` and call it after `velocity_correct`, as the last step
  before "Finished reducing {night}".
- **Changed:** [scripts/kmos_test.cfg](../scripts/kmos_test.cfg) — `rawdir`/
  `procdir` repointed at the new OB-sorted tree (`kmos_ob_sorted` /
  `kmos_ob_proc`), `subdirs` rewritten as the 20 individual OB folder names
  for the same three test nights used before, and `telluric_corr = True`
  added under `[cubes]`.

### Design decisions made (and why)

- **Where to find the telluric file.** `run_params['cubes']` (what
  `reduce_kmos.py` passes into every cube-stage function) sets `inputdir` to
  `procdir`, not `rawdir` — the cube steps operate on already-reconstructed
  products in `procdir`. But the telluric file lives in `rawdir/<OB>/`
  alongside the raw data. So `telluric_correct()` takes an explicit `rawdir`
  argument, passed separately at the call site in `reduce_kmos.py`
  (`rawdir=rawdir`, alongside `**run_params['cubes']`), rather than relying
  on `kwargs['inputdir']`.
- **Applied last in the chain**, after heliocentric velocity correction
  (`search_suffix='s_i_bkg2d_rv'`), producing a final `..._s_i_bkg2d_rv_tell.fits`.
  Nothing in the existing code or config said where in the order telluric
  correction "should" go, so it was added as an additional final step rather
  than reordering the existing chain — the safer choice given the existing
  order wasn't something this task was asked to change.
- **Guard against dividing by near-zero transmission.** Deep telluric
  absorption cores can approach zero transmission; dividing by that would
  blow up the noise arbitrarily. Values ≤ 0.1 are set to NaN before
  dividing, rather than corrected. The 0.1 floor is a judgement call, not
  something derived from the data — worth revisiting once real corrected
  cubes can be inspected.
- **Wavelength-grid consistency check.** Before applying, each IFU's
  telluric header (`NAXIS1`, `CRVAL1`) is compared against the science
  cube's (`NAXIS3`, `CRVAL3`); a mismatch is logged and that IFU is skipped
  rather than silently dividing misaligned arrays.
- **Per-IFU spectrum logged.** The function logs
  `min`/`max`/`median` of each IFU's telluric spectrum before applying it —
  see "What's not yet verified" below for why.

### Verified against a real run

The first real run (3 test nights, 20 OBs, OB-sorted tree, `telluric_corr =
True`) completed the telluric step cleanly for multiple OBs — `- Applying
telluric correction` / `-- Finished telluric correction` with no traceback.

1. **Whether the `DATA` extension is actually a 0-1 transmission fraction —
   confirmed.** The per-IFU `min`/`max`/`median` debug lines (found in each
   OB's own `kreuter_<band>_<timestamp>.log` under `kmos_ob_proc/<OB>/` —
   not the master log, since the file handler logs at DEBUG regardless of
   `--verbose`, while `--verbose 3` only raises the console/stdout stream to
   INFO) came back in the range expected for a transmission spectrum. The
   load-bearing assumption behind the whole implementation holds.

### Still not independently verified

1. **The `CDELT1` fallback.** The diagnostic dump didn't show a `CDELT1`
   key in the telluric header (only `CRVAL1`/`CRPIX1` were printed), which
   suggests it may be stored as `CD1_1` instead — matching the
   `CDELT3`/`CD3_3` fallback pattern already used elsewhere in this codebase
   (`reconstruct_wave()` in `pyspark/utils/obsutils.py`). The consistency
   check implemented here only compares `NAXIS1`/`CRVAL1`, not the actual
   wavelength step size, so this wasn't blocking — but it means the step
   size itself is still an unverified assumption.
2. **Whether the corrected output is scientifically sound**, beyond "ran
   without crashing and stayed in a plausible value range" — e.g. whether
   known telluric absorption features actually disappear from the corrected
   spectra when inspected directly, and whether the 0.1 transmission floor
   (see "Design decisions" above) is the right threshold rather than just a
   reasonable-sounding guess.

### How to run it

```bash
python3 lodo_sorting_algrithm.py --json kmos_ob_structure.json --limit 0 \
    --dest /home/kreuter/kmos_ob_sorted
cd ~/pyspark_my_way/scripts
nohup python3 reduce_kmos.py -p kmos_test.cfg --verbose 3 > reduce_run_telluric.log 2>&1 &
disown
```

Watch for the new `- Applying telluric correction` / `-- Finished telluric
correction` log lines, and check the per-IFU spectrum range logged at debug
level to confirm assumption (1) above before trusting the output.
