"""Per-galaxy archival KMOS time for the SPT2040 and XMM2215 cluster fields.

The local IFU tables only cover COSMOS/ECDFS, so here we query the ESO archive
(TAP dbo.raw) for KMOS science frames on the cluster fields, fetch one raw
header per OB to read the 24 arm positions + band, and match arms to our
targets. Output: kmos_council_cluster_archive.csv
"""
import csv, glob, html, os, re, json, time, urllib.parse, urllib.request
import numpy as np
from collections import defaultdict

OUR_PROG = "116.29EG"
CACHE = os.path.join(os.environ.get("TEMP", "."), "claude", "hdr")
os.makedirs(CACHE, exist_ok=True)
# field -> (ra, dec, only_ours). COSMOS/ECDFS archival time comes from the local
# IFU tables (council_numbers.py); here we only add our own not-yet-stacked OBs.
FIELDS = {"SPT2040_Hbeta_OIII": (310.24, -44.86, False), "XMM2215_EMPOWER.cat": (333.995, -17.635, False),
          "cosmos_0.6_0.8_h21.4.cat": (149.925, 2.53, True), "cdfs_targets_2.cat": (53.095, -27.735, True)}


def tap(q):
    u = "https://archive.eso.org/tap_obs/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
    return list(csv.DictReader(urllib.request.urlopen(u, timeout=300).read().decode().splitlines()))


def header(dpid):
    fn = os.path.join(CACHE, dpid.replace(":", "_") + ".txt")
    if not os.path.exists(fn):
        u = "https://archive.eso.org/hdr?DpId=" + dpid
        txt = urllib.request.urlopen(u, timeout=120).read().decode("utf-8", "replace")
        open(fn, "w", encoding="utf-8").write(txt)
        time.sleep(0.3)
    return html.unescape(open(fn, encoding="utf-8").read())


def sexa(v, ra):
    # ESO OCS ARMi ALPHA = HHMMSS.sss ; DELTA = +/-DDMMSS.sss
    v = float(v); s = -1 if v < 0 else 1; v = abs(v)
    d, m = int(v // 10000), int((v % 10000) // 100); sec = v % 100
    x = d + m / 60 + sec / 3600
    return x * 15 if ra else s * x


def arms(txt):
    kv = dict(re.findall(r"HIERARCH ESO (OCS ARM\d+ (?:ALPHA|DELTA|NAME|TYPE))\s*=\s*'?([^'/]+?)'?\s*/", txt))
    band = re.search(r"HIERARCH ESO INS FILT1 ID\s*=\s*'([^']+)'", txt)
    out = []
    for i in range(1, 25):
        a, d = kv.get(f"OCS ARM{i} ALPHA"), kv.get(f"OCS ARM{i} DELTA")
        if a and d:
            out.append((i, sexa(a, True), sexa(d, False), kv.get(f"OCS ARM{i} NAME", "").strip(),
                        kv.get(f"OCS ARM{i} TYPE", "").strip()))
    mode = re.search(r"HIERARCH ESO TPL MODE OBS\s*=\s*'([^']+)'", txt)
    return out, band.group(1).strip() if band else "?", mode.group(1).strip() if mode else "?"


# our targets in these fields
targets = [r for r in csv.DictReader(open("kmos_council_numbers.csv", encoding="utf-8"))
           if r["field"] in FIELDS]

frames = []
for f, (ra, dec, only_ours) in FIELDS.items():
    rows = tap(f"""SELECT dp_id, prog_id, ob_id, exposure, dp_type, tpl_start, filter_path
FROM dbo.raw WHERE instrument='KMOS' AND dp_cat='SCIENCE'
AND ra BETWEEN {ra-0.12} AND {ra+0.12} AND dec BETWEEN {dec-0.09} AND {dec+0.09}""")
    if only_ours:
        rows = [r for r in rows if r["prog_id"].startswith(OUR_PROG)]
    for r in rows:
        r["field"] = f
    frames += rows
print(f"{len(frames)} science frames; dp_types: {sorted(set(r['dp_type'] for r in frames))}")

by_ob = defaultdict(list)
for r in frames:
    by_ob[(r["field"], r["prog_id"], r["ob_id"], r["tpl_start"])].append(r)

res = defaultdict(lambda: defaultdict(float))  # (field,tid) -> {(prog,band): t_on}
for (f, prog, ob, tpl), fr in sorted(by_ob.items()):
    fr.sort(key=lambda r: r["dp_id"])
    try:
        al, band, mode = arms(header(fr[0]["dp_id"]))
    except Exception as e:
        print("  header fail", fr[0]["dp_id"], e); continue
    # NOD_TO_SKY runs O-S-O sequences: 2/3 of the frames are on source
    # (checked against our own OBs: 12 of 18 frames combined)
    t_on = sum(float(r["exposure"]) for r in fr) * (2 / 3 if mode == "NOD_TO_SKY" else 1.0)
    for t in targets:
        if t["field"] != f:
            continue
        tra, tde = float(t["ra"]), float(t["dec"])
        for i, ara, ade, nm, ty in al:
            if ty != "O":
                continue
            sep = np.hypot((ara - tra) * np.cos(np.radians(tde)), ade - tde) * 3600
            if sep < 1.0:
                res[(f, t["tid"])][(prog, band)] += t_on
print(f"{len(by_ob)} OB/template blocks scanned")

with open("kmos_council_cluster_archive.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["field", "tid", "prog_id", "band", "t_on_s", "is_ours"])
    for (f, tid), d in sorted(res.items()):
        for (prog, band), t in sorted(d.items()):
            w.writerow([f, tid, prog, band, t, prog.startswith(OUR_PROG)])
for (f, tid), d in sorted(res.items()):
    print(f[:8], tid, {f"{p}/{b}": round(t / 3600, 2) for (p, b), t in d.items()})
