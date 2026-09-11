#!/usr/bin/env python3
"""Cross-match this dataset's KMOS arm targets (COSMOS + CDFS/GOODS-S fields)
against public redshift catalogs hosted on IRSA, by sky position.

The KMOS OB headers carry a per-arm target ID (HIERARCH ESO OCS ARM{n} NAME)
and RA/Dec, but no redshift -- the parent target catalog that assigned those
IDs isn't part of these FITS files (see the note in telluric_check.ipynb's
first cell). Matching by ID against a guessed parent catalog is unreliable
(checked: COSMOS2020's own `id` column collides with our arm IDs by
coincidence, at a position half a degree away) -- position is the only
robust cross-match key here.

For each field, this pulls every catalog row inside a bounding box around
that field's targets (one query per catalog, not per target) and picks the
nearest match within MATCH_RADIUS_ARCSEC for each arm target, preferring
spectroscopic redshifts over photometric ones. Writes
scripts/target_redshifts.csv.

Run standalone: python scripts/fetch_target_redshifts.py
"""
import csv
import glob
import os
import ssl
import urllib.parse
import urllib.request

from astropy.io import fits
import numpy as np

try:
    import certifi
    SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CONTEXT = ssl.create_default_context()

BASE = "kmos_ob_proc_done"
IRSA_TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
MATCH_RADIUS_ARCSEC = 1.5
BOX_MARGIN_DEG = 0.05  # ~3 arcmin pad around each field's target spread
OUT_PATH = os.path.join("scripts", "target_redshifts.csv")

# One or more catalogs to try per field, in priority order (spec-z sources
# first, photo-z as a fallback), as
# (table, id_col, ra_col, dec_col, z_col, quality_col, "spec"|"photo").
# All are on IRSA and each is pulled with a single bounding-box query per
# field rather than one query per target.
FIELD_CATALOGS = {
    "cosmos": [
        ("cosmos.zcosmos_bright_v3_5_spec", "name", "ra", "dec", "z", "cc", "spec"),
        ("cosmos.cosmos_deimos_redshift", "id", "ra", "dec", "zspec", "qf", "spec"),
        ("cosmos2020c", "id", "ra", "dec", "lp_zbest", "lp_type", "photo"),
    ],
    "cdfs_targets_2.cat": [
        # MUSYC (ECDFS) carries both a compiled best spec-z and an EAzY
        # photo-z in the same table -- checked twice, spec-z first.
        ("musyc_photz", "num", "ra", "dec", "z_spec", "z_qual", "spec"),
        ("musyc_photz", "num", "ra", "dec", "z_peak", "q_z", "photo"),
    ],
}


def sexagesimal_to_deg(value, is_ra):
    """ESO-style packed sexagesimal header value (+/-[D]DMMSS.sss, or
    HHMMSS.sss for RA) to decimal degrees."""
    sign = -1 if value < 0 else 1
    value = abs(value)
    whole_deg_or_hr = int(value // 10000)
    minutes = int((value % 10000) // 100)
    seconds = value % 100
    deg_or_hr = whole_deg_or_hr + minutes / 60 + seconds / 3600
    return sign * (deg_or_hr * 15 if is_ra else deg_or_hr)


def collect_targets(base):
    """field -> {arm_target_id: (ra_deg, dec_deg)} for every active-IFU
    target found across every OB under base."""
    targets = {}
    for after in sorted(glob.glob(os.path.join(base, "*", "*_s_i_bkg2d_rv_tell.fits"))):
        with fits.open(after) as hdul:
            hdr = hdul[0].header
            field = hdr.get("ESO OBS TARG NAME")
            if field is None:
                continue
            for ifu in range(1, 25):
                name_key = f"ESO OCS ARM{ifu} NAME"
                alpha_key = f"ESO OCS ARM{ifu} ALPHA"
                delta_key = f"ESO OCS ARM{ifu} DELTA"
                if name_key not in hdr or alpha_key not in hdr or delta_key not in hdr:
                    continue
                alpha, delta = hdr[alpha_key], hdr[delta_key]
                if alpha is None or delta is None:
                    continue
                ra = sexagesimal_to_deg(alpha, is_ra=True)
                dec = sexagesimal_to_deg(delta, is_ra=False)
                targets.setdefault(field, {})[str(hdr[name_key])] = (ra, dec)
    return targets


def tap_query_csv(adql):
    url = IRSA_TAP + "?" + urllib.parse.urlencode({"QUERY": adql, "LANG": "ADQL", "FORMAT": "csv"})
    req = urllib.request.Request(url, headers={"User-Agent": "pyspark_my_way redshift lookup"})
    with urllib.request.urlopen(req, timeout=60, context=SSL_CONTEXT) as r:
        return r.read().decode("utf-8", errors="replace")


def parse_csv(text):
    lines = text.strip().splitlines()
    return list(csv.DictReader(lines)) if lines else []


def angular_sep_arcsec(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    cos_sep = (np.sin(dec1) * np.sin(dec2)
               + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
    return np.degrees(np.arccos(np.clip(cos_sep, -1, 1))) * 3600


def nearest_row(ra0, dec0, rows, ra_col, dec_col, radius_arcsec):
    best = None
    for row in rows:
        try:
            ra, dec = float(row[ra_col]), float(row[dec_col])
        except (TypeError, ValueError, KeyError):
            continue
        sep = angular_sep_arcsec(ra0, dec0, ra, dec)
        if sep <= radius_arcsec and (best is None or sep < best[0]):
            best = (sep, row)
    return best


def usable_redshift(z_raw, kind):
    if z_raw in ("", None):
        return None
    try:
        z = float(z_raw)
    except ValueError:
        return None
    if kind == "spec" and (z <= 0 or z >= 9):  # e.g. zCOSMOS's 9.9999 = unknown
        return None
    return z


def main():
    targets = collect_targets(BASE)
    results = []

    for field, arm_targets in targets.items():
        catalogs = FIELD_CATALOGS.get(field)
        if not catalogs:
            print(f"no catalog configured for field {field!r}, skipping {len(arm_targets)} targets")
            continue

        ras = [ra for ra, dec in arm_targets.values()]
        decs = [dec for ra, dec in arm_targets.values()]
        ra_lo, ra_hi = min(ras) - BOX_MARGIN_DEG, max(ras) + BOX_MARGIN_DEG
        dec_lo, dec_hi = min(decs) - BOX_MARGIN_DEG, max(decs) + BOX_MARGIN_DEG

        table_rows = {}
        for table, id_col, ra_col, dec_col, z_col, qual_col, kind in catalogs:
            if table in table_rows:
                continue
            adql = (
                f"SELECT * FROM {table} WHERE "
                f"{ra_col} BETWEEN {ra_lo} AND {ra_hi} AND "
                f"{dec_col} BETWEEN {dec_lo} AND {dec_hi}"
            )
            print(f"querying {table} for field {field!r} ({len(arm_targets)} targets)...")
            table_rows[table] = parse_csv(tap_query_csv(adql))
            print(f"  -> {len(table_rows[table])} candidate rows in box")

        for target_id, (ra0, dec0) in sorted(arm_targets.items()):
            match_result = None
            for table, id_col, ra_col, dec_col, z_col, qual_col, kind in catalogs:
                nearest = nearest_row(ra0, dec0, table_rows[table], ra_col, dec_col, MATCH_RADIUS_ARCSEC)
                if nearest is None:
                    continue
                sep, row = nearest
                z = usable_redshift(row.get(z_col), kind)
                if z is None:
                    continue
                match_result = {
                    "field": field, "arm_target_id": target_id,
                    "ra_deg": round(ra0, 6), "dec_deg": round(dec0, 6),
                    "redshift": z, "redshift_type": kind,
                    "source_catalog": table, "source_id": row.get(id_col, ""),
                    "quality": row.get(qual_col, ""),
                    "separation_arcsec": round(sep, 3),
                }
                break  # first catalog in priority order with a usable z wins
            results.append(match_result or {
                "field": field, "arm_target_id": target_id,
                "ra_deg": round(ra0, 6), "dec_deg": round(dec0, 6),
                "redshift": "", "redshift_type": "", "source_catalog": "",
                "source_id": "", "quality": "", "separation_arcsec": "",
            })

    fieldnames = ["field", "arm_target_id", "ra_deg", "dec_deg", "redshift",
                  "redshift_type", "source_catalog", "source_id", "quality",
                  "separation_arcsec"]
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    n_matched = sum(1 for r in results if r["redshift"] != "")
    print(f"\n{n_matched}/{len(results)} targets matched -> wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
