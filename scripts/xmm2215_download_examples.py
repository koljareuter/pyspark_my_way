"""Download archival (non-116.29EG) KMOS Phase-3 cubes for example galaxies.
Usage: python scripts/xmm2215_download_examples.py 32400 29284"""
import csv, os, sys, urllib.request

tids = sys.argv[1:] or ["32400", "29284"]
rows = [p for p in csv.DictReader(open("xmm2215_archive_products.csv", encoding="utf-8"))
        if p["obs_collection"] == "KMOS" and p["tid"] in tids and not p["proposal_id"].startswith("116.29EG")]
with open(os.path.join("archive_xmm2215", "index.csv"), "w", newline="", encoding="utf-8") if os.makedirs("archive_xmm2215", exist_ok=True) is None else None as fh:
    w = csv.writer(fh)
    w.writerow(["tid", "dp_id", "proposal_id", "em_min", "em_max", "t_exptime", "file"])
    for i, p in enumerate(rows, 1):
        d = os.path.join("archive_xmm2215", p["tid"]); os.makedirs(d, exist_ok=True)
        fn = os.path.join(d, p["dp_id"].replace(":", "_") + ".fits")
        if not (os.path.exists(fn) and os.path.getsize(fn) > 0):
            urllib.request.urlretrieve("https://dataportal.eso.org/dataPortal/file/" + p["dp_id"], fn + ".part")
            os.replace(fn + ".part", fn)
        w.writerow([p["tid"], p["dp_id"], p["proposal_id"], p["em_min"], p["em_max"], p["t_exptime"], fn])
        print(f"[{i}/{len(rows)}] {fn} {os.path.getsize(fn)/1e6:.1f} MB", flush=True)
