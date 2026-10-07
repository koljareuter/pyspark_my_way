"""Programme-wide time accounting for 116.29EG from the ESO archive (dbo.raw).
Writes kmos_council_programme.csv (one row per OB execution)."""
import csv, urllib.parse, urllib.request
from collections import defaultdict
from datetime import datetime

rows = []
for sub in range(1, 10):  # 116.29EG.001 ... .009 (equality queries hit the index)
    q = f"""SELECT dp_id, prog_id, ob_id, ob_name, object, dp_cat, exposure, exp_start, filter_path, tpl_start
FROM dbo.raw WHERE prog_id='116.29EG.{sub:03d}'"""
    u = "https://archive.eso.org/tap_obs/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q, "MAXREC": 100000})
    got = list(csv.DictReader(urllib.request.urlopen(u, timeout=300).read().decode().splitlines()))
    print(f"116.29EG.{sub:03d}: {len(got)} frames")
    rows += got
print(len(rows), "frames for 116.29EG")

t = lambda s: datetime.fromisoformat(s.replace("Z", "")[:23])
# one OB execution = consecutive frames of the same ob_id with gaps < 20 min
ex = defaultdict(list)
byob = defaultdict(list)
for r in rows:
    byob[(r["prog_id"], r["ob_id"])].append(r)
for (prog, ob), fr in byob.items():
    fr.sort(key=lambda r: r["exp_start"])
    n, last = 0, None
    for r in fr:
        ts = t(r["exp_start"]).timestamp()
        if last is not None and ts - last > 20 * 60:
            n += 1
        last = ts + float(r["exposure"])
        ex[(prog, ob, n)].append(r)

out = []
for (prog, ob, night), fr in ex.items():
    sci = [r for r in fr if r["dp_cat"] == "SCIENCE"]
    if not sci:
        continue
    t0 = min(t(r["exp_start"]) for r in fr)
    t1 = max(t(r["exp_start"]).timestamp() + float(r["exposure"]) for r in fr)
    band = sci[0]["filter_path"].split(",")[0]
    out.append({"prog_id": prog, "ob_id": ob, "ob_name": sci[0]["ob_name"], "target": sci[0]["object"],
                "band": band, "date": t0.date().isoformat(), "n_sci": len(sci),
                "shutter_s": sum(float(r["exposure"]) for r in sci),
                "elapsed_s": round(t1 - t0.timestamp())})
out.sort(key=lambda r: r["date"])
with open("kmos_council_programme.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)

tot = defaultdict(lambda: [0, 0, 0.0, 0.0, set()])
for r in out:
    k = (r["prog_id"], r["band"])
    tot[k][0] += 1; tot[k][1] += r["n_sci"]; tot[k][2] += r["shutter_s"]; tot[k][3] += r["elapsed_s"]
    tot[k][4].add(r["date"])
for k, v in sorted(tot.items()):
    print(f"{k[0]:14s} {k[1]:3s} OBexec={v[0]:3d} sci={v[1]:4d} shutter={v[2]/3600:6.1f}h elapsed={v[3]/3600:6.1f}h nights={len(v[4])}")
print("TOTAL OBexec", len(out), "shutter %.1fh" % (sum(r['shutter_s'] for r in out)/3600),
      "elapsed %.1fh" % (sum(r['elapsed_s'] for r in out)/3600), "nights", len(set(r['date'] for r in out)))
