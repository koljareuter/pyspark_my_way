#!/usr/bin/env python3
"""Add reduced DATACUBE targets missing from target_redshifts.csv.

Scans kmos_proc_and_comb_server/ for all DATACUBE_COMBINED_*.fits,
finds stems with no matching entry in target_redshifts.csv, reads
RA/Dec from FITS headers, and appends skeleton rows.

Run update_redshifts.py afterwards to fill redshifts from catalogs.
"""
import csv
import glob
import math
import os
import re

from astropy.io import fits

COMB_DIR   = "kmos_proc_and_comb_server"
TARGET_CSV = os.path.join("scripts", "target_redshifts.csv")

SKIP_RE      = re.compile(r"recenter[0-9]+", re.I)
_BAND_SUFFIX = re.compile(r"_(IZ|YJ|HK|H|K)$")

CSV_FIELDNAMES = [
    "field", "arm_target_id", "ra_deg", "dec_deg",
    "redshift", "redshift_type", "source_catalog",
    "source_id", "quality", "separation_arcsec",
]


def _parse_stem(stem):
    """Return (catalog, source_id) from a DATACUBE stem.

    Splits on the LAST underscore, so for 'XMM2215_EMPOWER.cat_R_28'
    the split gives catalog='XMM2215_EMPOWER.cat_R' and source_id='28'.
    The caller must handle compound source_ids separately.
    """
    parts = stem.rsplit("_", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (stem, stem)


def _get_radec(fp):
    try:
        with fits.open(fp, memmap=False) as hdul:
            h = hdul[0].header
            ra  = float(h.get("RA",  "nan"))
            dec = float(h.get("DEC", "nan"))
        return ra, dec
    except Exception:
        return math.nan, math.nan


# ── Load existing CSV ─────────────────────────────────────────────────────────
with open(TARGET_CSV, newline="", encoding="utf-8") as f:
    existing_rows = list(csv.DictReader(f))

# Build lookup: arm_target_id → all known variants (for flexible matching)
# Compound IDs like "R_28" are indexed under "R_28" AND "28".
existing_sids = set()
for row in existing_rows:
    tid = row["arm_target_id"].strip()
    existing_sids.add(tid)
    last = tid.rsplit("_", 1)[-1]
    if last != tid:
        existing_sids.add(last)

print(f"Existing CSV rows : {len(existing_rows)}")
print(f"Known arm_target_ids (+ last-parts): {len(existing_sids)}")


# ── Scan FITS files ───────────────────────────────────────────────────────────
new_rows  = []
stem_seen = set()

for fp in sorted(
    glob.glob(os.path.join(COMB_DIR, "**", "DATACUBE_COMBINED_*.fits"), recursive=True)
):
    fn   = os.path.basename(fp)
    stem = fn[len("DATACUBE_COMBINED_"):].replace(".fits", "")

    m = _BAND_SUFFIX.search(stem)
    if m:
        stem = stem[:m.start()]

    if SKIP_RE.search(stem) or stem in stem_seen:
        continue
    stem_seen.add(stem)

    catalog, source_id = _parse_stem(stem)

    # Already in CSV (by source_id or compound last-part)?
    if source_id in existing_sids:
        continue
    # Also try compound: e.g. "R_" + source_id → handles case where CSV has
    # "R_28" and source_id="28" (we DON'T want to add a duplicate "28" row).
    if ("R_" + source_id) in existing_sids:
        continue

    ra, dec = _get_radec(fp)
    new_rows.append({
        "field":             catalog,   # exact catalog name → update_redshifts.py needs this
        "arm_target_id":     source_id,
        "ra_deg":            f"{ra:.6f}"  if math.isfinite(ra)  else "",
        "dec_deg":           f"{dec:.6f}" if math.isfinite(dec) else "",
        "redshift":          "",
        "redshift_type":     "",
        "source_catalog":    "",
        "source_id":         "",
        "quality":           "",
        "separation_arcsec": "",
    })

print(f"\nNew rows to add  : {len(new_rows)}")
for r in new_rows:
    print(f"  {r['field']:<38} {r['arm_target_id']:<12}  "
          f"RA={r['ra_deg']}  Dec={r['dec_deg']}")

if not new_rows:
    print("Nothing to add.")
else:
    all_rows = existing_rows + new_rows
    with open(TARGET_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWritten {len(all_rows)} total rows → {TARGET_CSV}")
    print("Next: run   python scripts/update_redshifts.py")
