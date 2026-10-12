"""Per-galaxy summary of xmm2215_archive_products.csv -> xmm2215_archive_summary.md
KMOS hours = sum of t_exptime of the per-target Phase-3 cubes (on-source)."""
import csv
from collections import defaultdict

P = list(csv.DictReader(open("xmm2215_archive_products.csv", encoding="utf-8")))
Z = {r["arm_target_id"]: r for r in csv.DictReader(open("scripts/target_redshifts.csv", encoding="utf-8"))
     if r["field"] == "XMM2215_EMPOWER.cat"}
S = defaultdict(lambda: defaultdict(float))
progs = defaultdict(lambda: defaultdict(set))
for p in P:
    tid, col = p["tid"], p["obs_collection"]
    if col == "KMOS":
        band = "H" if float(p["em_min"]) > 1.4e-6 else "YJ"
        who = "own" if p["proposal_id"].startswith("116.29EG") else "arch"
        S[tid][f"{who}_{band}"] += float(p["t_exptime"]) / 3600
        if who == "arch":
            progs[tid][band].add(p["proposal_id"])
    elif col == "MUSE-DEEP":
        S[tid]["muse"] = max(S[tid]["muse"], float(p["t_exptime"]) / 3600)
    elif col == "ALMA":
        S[tid]["alma"] += 1
    elif col == "HAWKI" and p["dataproduct_type"] == "image":
        S[tid]["hawki"] += 1

L = ["# XMM2215 — ESO archive per galaxy\n",
     "KMOS: Hα (1.61 µm) in H, [OIII]/Hβ in YJ. MUSE: [OII]3727 at ≈0.92 µm (red end). Hours on source.\n",
     "| Target | z | z type | Own YJ | Own H | Archival KMOS YJ ([OIII]) | Archival KMOS H (Hα) | MUSE-DEEP [OII] | ALMA products | HAWK-I images |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for tid, z in Z.items():
    s = S[tid]
    f = lambda k: f"{s[k]:.1f} h" if s[k] else "—"
    pr = lambda b: f" ({', '.join(sorted(progs[tid][b]))})" if progs[tid][b] else ""
    L.append(f"| {tid} | {z['redshift']} | {z['redshift_type']} | {f('own_YJ')} | {f('own_H')} | "
             f"{f('arch_YJ')}{pr('YJ')} | {f('arch_H')}{pr('H')} | {f('muse')} | {int(s['alma']) or '—'} | {int(s['hawki']) or '—'} |")
n = lambda k, x=0: sum(S[t][k] > x for t in Z)
L.append(f"\n- Archival KMOS H band (Hα): **{n('arch_H')}/20** targets; YJ ([OIII]): **{n('arch_YJ')}/20**; "
         f"MUSE: **{n('muse')}/20**; ALMA: {n('alma')}/20.")
L.append("- SINFONI: 091.A-0364(A), 3 BCG OBs (XCS2215BCG1-3), H band, 0.25″ scale, 2.7 h, Sept 2013. "
         "The field centre is 9″ from 99991, so none of the 20 targets is in the 8″ field of view.")
L.append("- X-shooter: no data. FORS1/FORS2 (2003–08): field-level only; MXU slit positions would need the headers.")
open("xmm2215_archive_summary.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
