#!/usr/bin/env python3
"""Download a night of KMOS science data plus its correctly-associated raw
calibrations from the ESO archive, into a single flat folder as expected by
pyspark's reduce_kmos.py (rawdir/<night>/*.fits, no science/flat/dark
subfolders).

Usage
-----
    python3 fetch_kmos_night.py --night 2013-06-30 --band K \\
        --rawdir /mnt/c/Users/kolja/kmos_data

    # restrict science frames to a specific target name (substring match)
    python3 fetch_kmos_night.py --night 2013-06-30 --band K \\
        --rawdir /mnt/c/Users/kolja/kmos_data --target NGC1234

    # proprietary data: log in first so retrieve_data can fetch it
    python3 fetch_kmos_night.py --night 2013-06-30 --band K \\
        --rawdir /mnt/c/Users/kolja/kmos_data --eso-user myUsername

Why calSelector, not a "same night" CALIB query
-------------------------------------------------
An earlier version of this script queried dp_cat=CALIB for the same
calendar night and downloaded everything found. That is wrong on two
counts, confirmed directly against the archive for 2013-06-30:

1. It pulls in calibrations for every band observed that night, not just
   the one you want -- for 2013-06-30 that's IZ/YJ/H/K/HK mixed together
   (364 frames), when the K-band reduction only actually needs ~90.
2. Worse, it can simply miss what you need: 2013-06-30 has zero K-band (and
   zero IZ-band) FLAT,LAMP frames at all -- only HK/H/YJ were taken that
   night. Calibration validity windows span multiple days by design (KMOS
   flats/darks are often taken the following midday), so "same calendar
   night" is the wrong search entirely. That's why kmos_flat failed with
   "Null input data" / "Input frameset is not consistent" -- the SOF file
   pyspark built had no flat frames to give it, because none existed in
   the folder.

ESO's archive has a dedicated service for this: calSelector. Given a
science dataset ID, it returns exactly the raw calibration files that
recipe validity rules associate with it (wherever in time they actually
are). This script uses that via astroquery.eso.Eso.retrieve_data(...,
with_calib='raw'), and get_associated_files(...) for the --dry-run count.
For the 30 K-band science frames on 2013-06-30, calSelector finds 63
associated calibration frames (some from 2013-07-01 midday) -- a fifth of
what the old per-night approach would have downloaded, and it actually
includes the K-band flats.

Notes
-----
* Uses the ESO archive's 'night' filter (evening-to-morning) to select
  science frames -- simpler and more reliable than building a datetime
  stime/etime range by hand (which silently returned zero rows when
  tested against astroquery.eso 0.4.11).
* astroquery.eso.Eso.ROW_LIMIT defaults to 50 and silently truncates
  results above that; this script sets ROW_LIMIT = -1 (unlimited).
* Safe to re-run: retrieve_data() skips files already present in the
  destination directory, and any half-finished ".part" downloads from an
  interrupted run are removed first so they get retried instead of left
  as gaps.
* Only run this in WSL/Linux -- the reduction itself needs dfits/fitsort/
  esorex on PATH, even though this download step doesn't.
"""
from __future__ import print_function
import argparse
import os

from astroquery.eso import Eso


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--night', required=True,
                    help='Night start date, YYYY-MM-DD (local evening the run started)')
    p.add_argument('--band', required=True,
                    help='Observing band to restrict science frames to, e.g. K, H, YJ, IZ, HK')
    p.add_argument('--rawdir', required=True,
                    help='Parent raw data directory (matches "rawdir" in reduce.cfg); '
                         'a subfolder named <night> is created/used inside it')
    p.add_argument('--target', default=None,
                    help='Optional substring to filter science target names on')
    p.add_argument('--eso-user', default=None,
                    help='ESO archive username, only needed for proprietary (non-public) data. '
                         'You will be prompted for your password interactively.')
    p.add_argument('--sci-limit', type=int, default=None,
                    help='Maximum number of science frames to download (and to fetch '
                         'calSelector-associated calibrations for). Omit for no limit. '
                         'Note calibrations are NOT separately limited -- calSelector '
                         'returns whatever a valid reduction of the requested science '
                         'frames actually needs, and trimming that would just break '
                         'the reduction the way the old --limit did.')
    p.add_argument('--dry-run', action='store_true',
                    help='Only query and print how many science/calibration frames would be '
                         'downloaded; download nothing.')
    p.add_argument('--clean', action='store_true',
                    help='Delete all existing .fits/.part files in the night folder before '
                         'downloading, for a fresh start. Asks for confirmation unless --yes '
                         'is also given. Note: this is normally unnecessary -- already-downloaded '
                         'files are skipped automatically, so re-running is safe without it.')
    p.add_argument('--yes', action='store_true',
                    help='Skip the confirmation prompt for --clean.')
    return p.parse_args()


def dp_id_column(table):
    return 'DP.ID' if 'DP.ID' in table.colnames else table.colnames[0]


def target_name_column(table):
    for name in table.colnames:
        if 'TARG' in name.upper() and 'NAME' in name.upper():
            return name
    return None


def clean_night_folder(target_dir, assume_yes=False):
    """Delete all .fits/.part files already in target_dir, for a fresh
    start. Asks for confirmation first unless assume_yes is set."""
    victims = [f for f in os.listdir(target_dir) if f.endswith('.fits') or f.endswith('.part')]
    if not victims:
        print(f"--clean: nothing to delete in {target_dir}")
        return
    print(f"--clean: about to permanently delete {len(victims)} file(s) from {target_dir}")
    if not assume_yes:
        reply = input(f"Type 'yes' to confirm deleting these {len(victims)} file(s): ")
        if reply.strip().lower() != 'yes':
            print("Aborted -- nothing deleted.")
            raise SystemExit(1)
    for fname in victims:
        os.remove(os.path.join(target_dir, fname))
    print(f"Deleted {len(victims)} file(s).")


def clean_stale_partial_downloads(target_dir):
    """Remove any *.part left behind by a previous interrupted download
    so retrieve_data will fetch them again instead of leaving a gap."""
    for fname in os.listdir(target_dir):
        if fname.endswith('.part'):
            print(f" -> Removing stale partial download: {fname}")
            os.remove(os.path.join(target_dir, fname))


def main():
    args = parse_args()
    target_dir = os.path.join(args.rawdir, args.night)
    os.makedirs(target_dir, exist_ok=True)

    if args.clean:
        clean_night_folder(target_dir, assume_yes=args.yes)

    eso = Eso()
    eso.ROW_LIMIT = -1  # unlimited -- default of 50 silently truncates results
    if args.eso_user:
        eso.login(args.eso_user)

    clean_stale_partial_downloads(target_dir)

    # --- science frames -------------------------------------------------
    print(f"Querying science frames for {args.night} ({args.band}-band)...")
    sci_table = eso.query_instrument('kmos', column_filters={
        'night': args.night, 'dp_cat': 'SCIENCE', 'ins_filt1_name': args.band,
    })
    sci_ids = []
    if sci_table is not None:
        idcol = dp_id_column(sci_table)
        namecol = target_name_column(sci_table)
        for row in sci_table:
            if args.target and namecol and args.target not in str(row[namecol]):
                continue
            sci_ids.append(row[idcol])
    print(f"Found {len(sci_ids)} science frame(s)")

    if args.sci_limit is not None and len(sci_ids) > args.sci_limit:
        print(f"Limiting to {args.sci_limit} of {len(sci_ids)} science frame(s)")
        sci_ids = sci_ids[:args.sci_limit]

    if not sci_ids:
        print("No science frames found -- nothing to do.")
        return

    # --- calSelector: find the raw calibrations actually associated with
    # these science frames (may span multiple days) -------------------
    print("Querying calSelector for associated raw calibrations...")
    calib_ids = eso.get_associated_files(sci_ids)
    print(f"calSelector found {len(calib_ids)} associated calibration frame(s)")

    if args.dry_run:
        print(f"Dry run: would download {len(sci_ids)} science + {len(calib_ids)} "
              f"calibration frame(s) into {target_dir}; downloading nothing.")
        return

    all_ids = sci_ids + calib_ids
    print(f"Downloading {len(all_ids)} total frame(s) into {target_dir}")
    # retrieve_data skips files already present in destination and
    # gunzips compressed downloads automatically (unzip=True default).
    eso.retrieve_data(all_ids, destination=target_dir)

    print(f"Done. {target_dir} now has "
          f"{len([f for f in os.listdir(target_dir) if f.endswith('.fits')])} FITS files.")


if __name__ == '__main__':
    main()
