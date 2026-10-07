# Two-hop Symlink Fix: kmos_flat "File not found"

**Date:** 2026-09-28  
**Affected stage:** calibration (`kmos_flat` DRS recipe, called from `make_spark_calibs`)  
**Symptom:** Every OB failed at calibration with `Error: File not found` in the
band-specific log, despite the flat FITS files visibly existing in
`kmos_ob_sorted/<OB>/`. `astropy.io.fits.open()` could open the same files
fine — only the ESO CPL layer inside `kmos_flat` could not.

---

## Background

`kmos_ob_sorted` is built by `scripts/lodo_sorting_algrithm.py`, which reads
a JSON dataset description and symlinks each raw frame into a per-OB directory.
The source paths in the JSON point into `.../EMPOWER_raw_data/input_all/`, but
the actual files live under `.../EMPOWER_raw_data/cosmos/<night>/`.
`resolve_source_path()` in the sorting script falls back to the cosmos path,
so each link is created as:

```
kmos_ob_sorted/<OB>/KMOS.*.fits  →  ../input_all/KMOS.*.fits
```

where `input_all/KMOS.*.fits` is itself a symlink into `input3/`:

```
input_all/KMOS.*.fits  →  ../input3/KMOS.*.fits
```

This produces a **two-hop chain**: `kmos_ob_sorted/<OB>/X.fits` →
`input_all/X.fits` → `input3/X.fits`.

`astropy` (and CFITSIO generally) follows symlink chains fine.
The CPL layer used by the ESO DRS (`kmos_flat`, `kmos_dark`, etc.) resolves
only the **first** symlink hop, lands at `input_all/X.fits`, and then attempts
to open that path literally — which is itself a symlink file, not a regular
file. CPL treats it as non-existent and raises "File not found".

---

## Fix

### Immediate: `scripts/relink_direct.py` (new file)

A one-shot script that walks every symlink in `kmos_ob_sorted` and replaces
two-hop chains with a direct single-hop link using `os.path.realpath()`:

```python
import os

OB_DIR = '/home/kreuter/kmos_ob_sorted'
updated = 0; already_direct = 0; broken = 0

for ob in sorted(os.listdir(OB_DIR)):
    ob_path = os.path.join(OB_DIR, ob)
    if not os.path.isdir(ob_path):
        continue
    for fname in sorted(os.listdir(ob_path)):
        fpath = os.path.join(ob_path, fname)
        if not os.path.islink(fpath):
            continue
        current_target = os.readlink(fpath)
        real_target = os.path.realpath(fpath)
        if not os.path.exists(real_target):
            broken += 1
            continue
        if current_target == real_target:
            already_direct += 1
            continue
        os.unlink(fpath)
        os.symlink(real_target, fpath)
        updated += 1

print(f'Updated to direct: {updated}')
print(f'Already direct:    {already_direct}')
print(f'Broken (skipped):  {broken}')
```

Run on hqempower1 (result: `Updated: 6069, Already direct: 275, Broken: 0`).

### Permanent: `scripts/lodo_sorting_algrithm.py` line 131

Changed `os.path.abspath` → `os.path.realpath` so future sorts always create
direct symlinks to the canonical target:

```python
# before
os.symlink(os.path.abspath(real_path), dest_path)

# after
os.symlink(os.path.realpath(real_path), dest_path)
```

`abspath` normalizes `../` components but does **not** follow symlinks.
`realpath` follows all symlinks and returns the canonical on-disk path,
guaranteeing a single-hop link regardless of how many intermediate symlinks
exist in the raw-data tree.

---

## Deployment

```bash
# on hqempower1, from ~/pyspark_my_way
git pull
python scripts/relink_direct.py      # one-time repair of existing links
# restart reduction as normal
```

No PROCDIR cleanup needed — calibrations had not completed for any affected OB.
