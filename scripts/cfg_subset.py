#!/usr/bin/env python
"""Write a copy of a pyspark config that reduces only some of its OBs.

The OBs are those of the config's [pathing] subdirs whose folder in its procdir has a regular file
(not a link) modified since --modified-since, or the ones given with --obs. Everything else in the
config is copied unchanged.  Used to redo the OBs of an interrupted or doubled run.

Usage (on hqempower1, from ~/pyspark_my_way):
    python scripts/cfg_subset.py scripts/kmos_server_hyb_IZ.cfg /local_storage/KOLJA/repair_IZ.cfg \
           --modified-since "2026-11-04 00:00"
    python scripts/cfg_subset.py scripts/kmos_server_hyb_YJ.cfg /local_storage/KOLJA/repair_YJ.cfg \
           --obs KMOS.2026-02-11T06:31:26.207,KMOS.2026-06-08T05:24:25.056
"""
import argparse
import configparser
import os
import re
import sys
from datetime import datetime


def modified_since(folder, t0):
    for root, _, files in os.walk(folder):
        for f in files:
            p = os.path.join(root, f)
            if not os.path.islink(p) and os.path.getmtime(p) >= t0:
                return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cfg")
    ap.add_argument("out")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--modified-since", help="'YYYY-MM-DD HH:MM' (local time)")
    g.add_argument("--obs", help="comma-separated OB folder names")
    a = ap.parse_args()

    text = open(a.cfg).read()
    c = configparser.ConfigParser()
    c.read_string(text)
    procdir = c.get("pathing", "procdir").strip()
    obs = [s.strip() for s in re.split(",|\n", c.get("pathing", "subdirs")) if s.strip()]
    if a.obs:
        want = [s.strip() for s in a.obs.split(",") if s.strip()]
        unknown = [o for o in want if o not in obs]
        if unknown:
            sys.exit("not in %s: %s" % (a.cfg, ", ".join(unknown)))
        keep = [o for o in obs if o in want]
    else:
        t0 = datetime.strptime(a.modified_since, "%Y-%m-%d %H:%M").timestamp()
        keep = [o for o in obs if os.path.isdir(os.path.join(procdir, o)) and modified_since(os.path.join(procdir, o), t0)]
    if not keep:
        sys.exit("no OB selected; nothing written")

    # replace the subdirs block (first line 'subdirs = ...' plus its indented continuation lines)
    new = re.sub(r"(?m)^subdirs\s*=.*(?:\n[ \t]+\S.*)*", "subdirs = " + "\n          ".join(keep), text, count=1)
    with open(a.out, "w") as f:
        f.write(new)
    print("%s: %d of %d OBs -> %s (procdir %s)" % (a.cfg, len(keep), len(obs), a.out, procdir))
    for o in keep:
        print("  " + o)


if __name__ == "__main__":
    main()
