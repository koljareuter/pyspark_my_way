"""Per-galaxy ESO archive search for the 20 XMM2215 targets.

Phase-3 products (ivoa.ObsCore): footprint contains the target (MUSE, HAWK-I,
ALMA) or product centre within 1" (KMOS per-arm IDP cubes).
Raw frames (dbo.raw): SINFONI / X-shooter / MUSE pointing within a 30" box;
FORS field-level only (MXU slit positions are not in dbo.raw).
Output: xmm2215_archive_products.csv, xmm2215_archive_raw.csv
"""
import csv, io, urllib.error, urllib.parse, urllib.request

def tap(q, ep="tap_obs"):
    u = f"https://archive.eso.org/{ep}/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q, "MAXREC": 200000})
    try:
        return list(csv.DictReader(io.StringIO(urllib.request.urlopen(u, timeout=300).read().decode())))
    except urllib.error.HTTPError as e:
        raise RuntimeError(e.read().decode()[:600])

targets = [r for r in csv.DictReader(open("scripts/target_redshifts.csv", encoding="utf-8"))
           if r["field"] == "XMM2215_EMPOWER.cat"]
COLS = ("obs_collection, instrument_name, dp_id, obs_id, proposal_id, target_name, s_ra, s_dec, "
        "em_min, em_max, t_exptime, abmaglim, s_resolution, obs_release_date, access_estsize, "
        "access_url, dataproduct_type")

prod, raw = [], []
for t in targets:
    ra, dec = float(t["ra_deg"]), float(t["dec_deg"])
    pt = f"POINT('ICRS',{ra},{dec})"
    rows = tap(f"""SELECT {COLS} FROM ivoa.ObsCore WHERE
        (obs_collection='KMOS' AND CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',{ra},{dec},{1/3600}))=1)
     OR (obs_collection<>'KMOS' AND INTERSECTS({pt}, s_region)=1)""")
    for r in rows:
        r["tid"] = t["arm_target_id"]; r["z"] = t["redshift"]
    prod += rows
    rows = tap(f"""SELECT instrument, dp_id, prog_id, ob_name, object, ra, dec, exposure, filter_path,
        grat_path, ins_mode, dp_tech, date_obs, release_date FROM dbo.raw
        WHERE dp_cat='SCIENCE' AND instrument IN ('SINFONI','XSHOOTER','MUSE')
        AND ra BETWEEN {ra - 30/3600/0.953} AND {ra + 30/3600/0.953}
        AND dec BETWEEN {dec - 30/3600} AND {dec + 30/3600}""")  # box: dbo.raw has invalid dec rows
    for r in rows:
        r["tid"] = t["arm_target_id"]
    raw += rows
    print(f"{t['arm_target_id']:>6s}: {len([p for p in prod if p['tid']==t['arm_target_id']])} products, "
          f"{len(rows)} raw SINFONI/XSHOOTER/MUSE frames within 30\"")

for fn, data in (("xmm2215_archive_products.csv", prod), ("xmm2215_archive_raw.csv", raw)):
    if data:
        with open(fn, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
