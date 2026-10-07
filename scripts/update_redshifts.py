#!/usr/bin/env python3
"""Update target_redshifts.csv with better redshift sources.

COSMOS fields  : cross-match against cosmos_specz_good.csv (local file,
                 spec-z with quality flag).  Only entries with flag >= MIN_FLAG
                 and within MATCH_RADIUS_ARCSEC are used.

CDFS / GOODS-S : query ASTRODEEP GS43 (Merlin+2021, J/A+A/649/A22) via
                 VizieR TAP — zBest is spec-z when available, else photo-z.

Cluster / other fields are left untouched.

Usage:
    python scripts/update_redshifts.py
    python scripts/update_redshifts.py --min-flag 3   # only secure COSMOS spec-z
    python scripts/update_redshifts.py --dry-run       # print stats, no write
    python scripts/update_redshifts.py --no-astrodeep  # skip VizieR download
"""
import argparse
import csv
import math
import os
import ssl
import urllib.parse
import urllib.request

import numpy as np

try:
    import certifi
    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    _SSL_CTX = ssl.create_default_context()

TARGET_CSV        = os.path.join('scripts', 'target_redshifts.csv')
COSMOS_SPECZ_CSV  = os.path.join('scripts', 'cosmos_specz_good.csv')

MATCH_RADIUS_ARCSEC = 1.5
BOX_MARGIN_DEG      = 0.05

# All field names (case-insensitive substring match) treated as COSMOS
COSMOS_FIELD_TAG  = 'cosmos'

# Substrings in field name (lower-case) that trigger ASTRODEEP GS43 query
ASTRODEEP_TAGS = ('cdfs', 'goods', 'ecdfs')

VIZIER_TAP = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ang_sep_arcsec(ra0: float, dec0: float,
                    ras: np.ndarray, decs: np.ndarray) -> np.ndarray:
    ra0r  = math.radians(ra0)
    dec0r = math.radians(dec0)
    rasr  = np.radians(ras)
    decsr = np.radians(decs)
    cos_sep = (np.sin(dec0r) * np.sin(decsr)
               + np.cos(dec0r) * np.cos(decsr) * np.cos(ra0r - rasr))
    return np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0))) * 3600.0


def _nearest(ra0, dec0, ras, decs, radius):
    """Return index and separation of nearest entry within radius, or None."""
    seps = _ang_sep_arcsec(ra0, dec0, ras, decs)
    within = np.where(seps <= radius)[0]
    if len(within) == 0:
        return None, None
    best = within[np.argmin(seps[within])]
    return int(best), float(seps[best])


# ---------------------------------------------------------------------------
# COSMOS loader
# ---------------------------------------------------------------------------

def load_cosmos_specz(csv_path: str, min_flag: int = 2):
    """Return parallel arrays (ra, dec, specz, flag) for flag >= min_flag."""
    ras, decs, speczs, flags = [], [], [], []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                flag  = int(row['flag'])
                specz = float(row['specz'])
            except (ValueError, KeyError):
                continue
            if flag < min_flag or specz <= 0 or specz >= 9:
                continue
            ra = dec = None
            for prefix in ('COS20_Classic', 'COS20_Farmer', 'COSMOS25'):
                try:
                    r = float(row[f'ra_{prefix}'])
                    d = float(row[f'dec_{prefix}'])
                    if r != -999.0:
                        ra, dec = r, d
                        break
                except (ValueError, KeyError):
                    continue
            if ra is not None:
                ras.append(ra); decs.append(dec)
                speczs.append(specz); flags.append(flag)
    print(f'  COSMOS specz: loaded {len(ras)} entries '
          f'(flag >= {min_flag}) from {os.path.basename(csv_path)}')
    return (np.array(ras), np.array(decs),
            np.array(speczs), np.array(flags, dtype=int))


# ---------------------------------------------------------------------------
# ASTRODEEP GS43 (Merlin+2021) via VizieR TAP
# ---------------------------------------------------------------------------

def _tap_csv(adql: str) -> list:
    url = (VIZIER_TAP + '?'
           + urllib.parse.urlencode({'QUERY': adql, 'LANG': 'ADQL',
                                     'FORMAT': 'csv', 'REQUEST': 'doQuery'}))
    req = urllib.request.Request(url, headers={'User-Agent': 'pyspark_my_way'})
    with urllib.request.urlopen(req, timeout=60, context=_SSL_CTX) as r:
        text = r.read().decode('utf-8', errors='replace')
    lines = [l for l in text.strip().splitlines() if not l.startswith('#')]
    return list(csv.DictReader(lines)) if lines else []


def fetch_astrodeep_gs43(ra_lo: float, ra_hi: float,
                         dec_lo: float, dec_hi: float):
    """Download ASTRODEEP GS43 positions + best redshifts via VizieR TAP.

    zBest is the spectroscopic redshift when the source has one
    (zspec-survey non-empty), otherwise the best photometric estimate.
    Returns parallel arrays (ra, dec, redshift, redshift_type).
    """
    # Positions for sources in the bounding box
    q_pos = (
        f'SELECT ID, RAJ2000, DEJ2000 '
        f'FROM "J/A+A/649/A22/gs43phot" '
        f'WHERE RAJ2000 BETWEEN {ra_lo} AND {ra_hi} '
        f'AND DEJ2000 BETWEEN {dec_lo} AND {dec_hi}'
    )
    print(f'  ASTRODEEP: querying gs43phot in box '
          f'({ra_lo:.3f}-{ra_hi:.3f}, {dec_lo:.3f}-{dec_hi:.3f})...')
    pos_rows = _tap_csv(q_pos)
    if not pos_rows:
        print('  ASTRODEEP: no sources found in box')
        return np.array([]), np.array([]), np.array([]), np.array([])
    print(f'  ASTRODEEP: {len(pos_rows)} sources in box')

    id_set = {r['ID'] for r in pos_rows}
    pos_map = {r['ID']: (float(r['RAJ2000']), float(r['DEJ2000']))
               for r in pos_rows}

    # Redshifts for those IDs
    q_z = (
        'SELECT ID, zBest, "zspec-survey" '
        'FROM "J/A+A/649/A22/gs43phys" '
        'WHERE zBest IS NOT NULL'
    )
    print('  ASTRODEEP: querying gs43phys for redshifts...')
    z_rows = _tap_csv(q_z)
    print(f'  ASTRODEEP: {len(z_rows)} redshift entries downloaded')

    ras, decs, zs, ztypes = [], [], [], []
    for row in z_rows:
        if row.get('ID') not in id_set:
            continue
        try:
            z = float(row['zBest'])
        except (ValueError, TypeError):
            continue
        if not (0 < z < 15):
            continue
        ra, dec = pos_map[row['ID']]
        zspec = row.get('zspec-survey', '').strip()
        ztype = 'spec' if zspec and zspec != '-' else 'photo'
        ras.append(ra); decs.append(dec); zs.append(z); ztypes.append(ztype)

    n_spec = sum(1 for t in ztypes if t == 'spec')
    print(f'  ASTRODEEP: {len(ras)} matched (spec: {n_spec}, '
          f'photo: {len(ras)-n_spec})')
    return (np.array(ras), np.array(decs),
            np.array(zs), np.array(ztypes))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-flag', type=int, default=2,
                        help='Minimum quality flag for cosmos_specz_good (default 2)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print statistics without writing the CSV')
    parser.add_argument('--no-astrodeep', action='store_true',
                        help='Skip ASTRODEEP VizieR query (offline mode)')
    args = parser.parse_args()

    # ── Load CSV ──────────────────────────────────────────────────────────
    with open(TARGET_CSV, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys())

    # ── COSMOS catalog ────────────────────────────────────────────────────
    cos_ras, cos_decs, cos_z, cos_flag = load_cosmos_specz(
        COSMOS_SPECZ_CSV, min_flag=args.min_flag)

    # ── ASTRODEEP GS43 for CDFS  ──────────────────────────────────────────
    ad_ras = ad_decs = ad_z = ad_zt = None
    if not args.no_astrodeep:
        cdfs_rows = [r for r in rows
                     if any(tag in r.get('field', '').lower() for tag in ASTRODEEP_TAGS)]
        if cdfs_rows:
            ras_c  = [float(r['ra_deg'])  for r in cdfs_rows]
            decs_c = [float(r['dec_deg']) for r in cdfs_rows]
            ra_lo  = min(ras_c)  - BOX_MARGIN_DEG
            ra_hi  = max(ras_c)  + BOX_MARGIN_DEG
            dec_lo = min(decs_c) - BOX_MARGIN_DEG
            dec_hi = max(decs_c) + BOX_MARGIN_DEG
            try:
                ad_ras, ad_decs, ad_z, ad_zt = fetch_astrodeep_gs43(
                    ra_lo, ra_hi, dec_lo, dec_hi)
            except Exception as exc:
                print(f'  ASTRODEEP query failed: {exc}')

    # ── Update rows ───────────────────────────────────────────────────────
    n_cos_up = n_cos_kept = n_ad_up = n_ad_kept = 0

    for row in rows:
        field = row.get('field', '')
        try:
            ra0  = float(row['ra_deg'])
            dec0 = float(row['dec_deg'])
        except (ValueError, KeyError):
            continue

        if COSMOS_FIELD_TAG in field.lower():
            idx, sep = _nearest(ra0, dec0, cos_ras, cos_decs, MATCH_RADIUS_ARCSEC)
            if idx is not None:
                row['redshift']          = round(float(cos_z[idx]), 6)
                row['redshift_type']     = 'spec'
                row['source_catalog']    = 'cosmos_specz_good'
                row['source_id']         = ''
                row['quality']           = str(cos_flag[idx])
                row['separation_arcsec'] = round(sep, 3)
                n_cos_up += 1
            else:
                n_cos_kept += 1

        elif any(tag in field.lower() for tag in ASTRODEEP_TAGS) and ad_ras is not None and len(ad_ras):
            idx, sep = _nearest(ra0, dec0, ad_ras, ad_decs, MATCH_RADIUS_ARCSEC)
            if idx is not None:
                row['redshift']          = round(float(ad_z[idx]), 6)
                row['redshift_type']     = str(ad_zt[idx])
                row['source_catalog']    = 'astrodeep_gs43'
                row['source_id']         = ''
                row['quality']           = ''
                row['separation_arcsec'] = round(sep, 3)
                n_ad_up += 1
            else:
                n_ad_kept += 1

    print(f'\nCOSMOS    : {n_cos_up} updated, {n_cos_kept} no match (kept existing)')
    print(f'ASTRODEEP : {n_ad_up} updated, {n_ad_kept} no match (kept existing)')

    if args.dry_run:
        print('Dry run — nothing written.')
        return

    with open(TARGET_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f'Written -> {TARGET_CSV}')


if __name__ == '__main__':
    main()
