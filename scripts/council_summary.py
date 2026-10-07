"""Merge per-galaxy S/N + own/archival time into the council summary.

Inputs (run these first):
  python scripts/council_numbers.py           -> kmos_council_numbers.csv
  python scripts/council_cluster_archive.py   -> kmos_council_cluster_archive.csv
  python scripts/council_programme_time.py    -> kmos_council_programme.csv
Outputs: kmos_council_galaxies.csv, kmos_council_summary.md, kmos_council_snr.png
"""
import csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINES = {"Ha": 0.65628, "OIII": 0.50069, "Hb": 0.48613}
LABEL = {"Ha": "Hα", "OIII": "[OIII]5007", "Hb": "Hβ"}
BAND_RANGES = {"IZ": (0.800, 1.080), "YJ": (1.000, 1.359), "H": (1.456, 1.850),
               "K": (1.934, 2.450), "HK": (1.456, 2.450)}
SHORT = {"cdfs_targets_2.cat": "CDFS", "cosmos_0.6_0.8_h21.4.cat": "COSMOS",
         "SPT2040_Hbeta_OIII": "SPT2040", "XMM2215_EMPOWER.cat": "XMM2215"}


def covers(band, line, z, pad=0.01):
    lo, hi = BAND_RANGES.get(band, (0, 0))
    return lo + pad < LINES[line] * (1 + z) < hi - pad


f = lambda v: float(v) if v not in ("", "nan", None) else np.nan
gal = list(csv.DictReader(open("kmos_council_numbers.csv", encoding="utf-8")))

# time per (field, tid) -> list of (source, band, seconds); source in {own, arch}
times = defaultdict(list)
for g in gal:  # archival (other programmes) from the local COSMOS/ECDFS IFU tables
    for item in filter(None, g["arch_bands"].split(";")):
        b, h = item.split(":")
        times[(g["field"], g["tid"])].append(("arch", b, float(h[:-1]) * 3600))
for r in csv.DictReader(open("kmos_council_cluster_archive.csv", encoding="utf-8")):
    src = "own" if r["is_ours"] == "True" else "arch"
    if src == "arch" and r["field"] not in ("SPT2040_Hbeta_OIII", "XMM2215_EMPOWER.cat"):
        continue  # COSMOS/CDFS archival already taken from the local tables
    times[(r["field"], r["tid"])].append((src, r["band"], float(r["t_on_s"])))

out = []
for g in gal:
    z = f(g["z"]); key = (g["field"], g["tid"])
    t_now = f(g["t_ours_s"])
    snr = {ln: f(g[f"snr_{ln}"]) for ln in LINES}
    # primary line: Ha if in the stacked cube, else the stronger of [OIII]/Hb
    line = next((ln for ln in LINES if np.isfinite(snr[ln])), None)
    if line is None:  # no line measurable (all at noisy band edges)
        print(f"  no measurable line: {g['field']} {g['tid']} z={z}")
        continue
    rec = {"field": SHORT.get(g["field"], g["field"]), "tid": g["tid"], "z": z,
           "z_type": g["z_type"], "band_stacked": g["band"], "line": line,
           "t_stacked_h": t_now / 3600, "snr_now": snr[line]}
    for ln in LINES:
        own = sum(t for s, b, t in times[key] if s == "own" and covers(b, ln, z))
        arch = sum(t for s, b, t in times[key] if s == "arch" and covers(b, ln, z))
        rec[f"t_own_{ln}_h"] = own / 3600
        rec[f"t_arch_{ln}_h"] = arch / 3600
    # own data covering the line that is not in the stack yet (never below 0)
    own_line = max(rec[f"t_own_{line}_h"], rec["t_stacked_h"])
    t_all = own_line + rec[f"t_arch_{line}_h"]
    rec["gain_own"] = np.sqrt(own_line / rec["t_stacked_h"])
    rec["gain_all"] = np.sqrt(t_all / rec["t_stacked_h"])
    rec["snr_own_all"] = rec["snr_now"] * rec["gain_own"]
    rec["snr_projected"] = rec["snr_now"] * rec["gain_all"]
    rec["t_Ha_total_h"] = max(rec["t_own_Ha_h"], rec["t_stacked_h"] if line == "Ha" else 0) + rec["t_arch_Ha_h"]
    out.append(rec)

with open("kmos_council_galaxies.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    for r in out:
        w.writerow({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})

# ── programme totals ─────────────────────────────────────────────────────────
prog = list(csv.DictReader(open("kmos_council_programme.csv", encoding="utf-8")))
P = {"n_exec": len(prog), "nights": len({p["date"] for p in prog}),
     "shutter_h": sum(float(p["shutter_s"]) for p in prog) / 3600,
     "elapsed_h": sum(float(p["elapsed_s"]) for p in prog) / 3600,
     "last": max(p["date"] for p in prog)}
by_band = defaultdict(float)
for p in prog:
    by_band[p["band"]] += float(p["elapsed_s"]) / 3600

# ── markdown summary ─────────────────────────────────────────────────────────
L = []
L.append("# KMOS programme 116.29EG — numbers for the council\n")
L.append(f"## 1. Telescope time used (ESO archive, as of {P['last']})\n")
L.append(f"- **{P['elapsed_h']:.0f} h** of telescope time (first acquisition to end of last exposure) "
         f"in **{P['n_exec']} OB executions** on **{P['nights']} nights**")
L.append(f"- {P['shutter_h']:.0f} h open shutter, about {P['shutter_h']*2/3:.0f} h of it on source "
         f"(O-S-O nod-to-sky: 2/3 of frames are on target)")
L.append("- By band: " + ", ".join(f"{b} {h:.0f} h" for b, h in sorted(by_band.items())) + "\n")

L.append("## 2. Stacked data set (kmos_proc_server_comb)\n")
L.append("| Field | z | Band | Galaxies | On-source / galaxy | Line | S/N ≥ 3 | S/N ≥ 5 | Best S/N |")
L.append("|---|---|---|---|---|---|---|---|---|")
for fld in ["COSMOS", "CDFS", "SPT2040", "XMM2215"]:
    rs = [r for r in out if r["field"] == fld]
    if not rs:
        continue
    for ln in sorted({r["line"] for r in rs}, key=list(LINES).index):
        rr = [r for r in rs if r["line"] == ln]
        s = np.array([r["snr_now"] for r in rr])
        L.append(f"| {fld} | {np.nanmedian([r['z'] for r in rr]):.2f} | {rr[0]['band_stacked']} | {len(rr)} | "
                 f"{rr[0]['t_stacked_h']:.1f} h | {LABEL[ln]} | {(s >= 3).sum()} | {(s >= 5).sum()} | {np.nanmax(s):.1f} |")
n = len(gal)
ha = [r for r in out if r["line"] == "Ha"]
L.append(f"\n- **{n} galaxies** in the stacked cubes; Hα falls in the stacked band for **{len(ha)}**, "
         f"of which **{sum(r['snr_now'] >= 5 for r in ha)} have Hα S/N ≥ 5** "
         f"({sum(r['snr_now'] >= 3 for r in ha)} with S/N ≥ 3).")
det_o = [r for r in out if r["line"] != "Ha" and r["snr_now"] >= 5]
L.append(f"- Another **{len(det_o)}** galaxies have [OIII]/Hβ at S/N ≥ 5 (z≈1.47 clusters and CDFS).")
top = sorted(ha, key=lambda r: -r["snr_now"])[:6]
L.append("- Brightest Hα: " + ", ".join(f"{r['field']} {r['tid']} (S/N {r['snr_now']:.0f})" for r in top) + "\n")

L.append("## 3. Extra data on the same galaxies\n")
L.append("| Field | Galaxies | Line | Stacked now | Own, not yet stacked | Archival (other programmes) | Depth gain √(t_all/t_now) | Expected S/N ≥ 5 (now → all) |")
L.append("|---|---|---|---|---|---|---|---|")
for fld in ["COSMOS", "CDFS", "SPT2040", "XMM2215"]:
    for ln in LINES:
        rr = [r for r in out if r["field"] == fld and r["line"] == ln]
        if not rr:
            continue
        own_extra = np.median([max(r[f"t_own_{ln}_h"] - r["t_stacked_h"], 0) for r in rr])
        arch = [r[f"t_arch_{ln}_h"] for r in rr]
        gain = np.median([r["gain_all"] for r in rr])
        L.append(f"| {fld} | {len(rr)} | {LABEL[ln]} | {rr[0]['t_stacked_h']:.1f} h | {own_extra:.1f} h (median) | "
                 f"median {np.median(arch):.1f} h, max {np.max(arch):.1f} h ({sum(a > 0 for a in arch)}/{len(rr)} galaxies) | "
                 f"×{gain:.1f} | {sum(r['snr_now'] >= 5 for r in rr)} → {sum(r['snr_projected'] >= 5 for r in rr)} |")
L.append("")
# Ha in other bands (CDFS IZ misses Ha at z>0.63, clusters need H band)
noha = [r for r in out if r["line"] != "Ha"]
L.append("**Hα outside the stacked band.** For the "
         f"{len(noha)} galaxies where Hα is not in the stacked cube (CDFS at z>0.63 in IZ; clusters at z≈1.47, where Hα is in H):")
for fld in ["CDFS", "SPT2040", "XMM2215"]:
    rr = [r for r in noha if r["field"] == fld]
    if not rr:
        continue
    own = [r["t_own_Ha_h"] for r in rr]; arch = [r["t_arch_Ha_h"] for r in rr]
    L.append(f"- {fld}: {len(rr)} galaxies; own Hα-band data median {np.median(own):.1f} h; archival Hα-band data "
             f"median {np.median(arch):.1f} h (max {np.max(arch):.1f} h, {sum(a > 0 for a in arch)} galaxies with some)")
L.append("")
L.append("## Method / caveats\n")
L.append("- S/N = integrated line flux in a 1.1″-diameter aperture centred on the line peak (±350 km/s, re-centred within "
         "±400 km/s of the catalogue z), divided by the scatter of same-width windows across the line-free spectrum. "
         "Sky-line residuals are therefore included in the noise.")
L.append("- Projection: S/N ∝ √t (background-limited, similar conditions). This assumes the archival data are of similar "
         "quality and combine optimally, so treat it as an upper estimate. For non-detections only the depth gain is meaningful.")
L.append("- Archival time = on-source exposure from other programmes with an IFU within 1″ of the target "
         "(COSMOS/ECDFS: local ESO IFU tables; SPT2040/XMM2215: ESO TAP + raw-header arm positions, nod-to-sky ×2/3).")
open("kmos_council_summary.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))

# ── figure: S/N now vs. with all available data ──────────────────────────────
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
cols = {"COSMOS": "#1f77b4", "CDFS": "#ff7f0e", "SPT2040": "#2ca02c", "XMM2215": "#d62728"}
for r in out:
    m = "o" if r["line"] == "Ha" else "s"
    ax[0].scatter(max(r["snr_now"], 0.1), max(r["snr_projected"], 0.1), c=cols[r["field"]], marker=m,
                  s=40, alpha=.8, edgecolor="k", lw=.4)
for x in (3, 5):
    ax[0].axhline(x, color="grey", ls=":", lw=.8); ax[0].axvline(x, color="grey", ls=":", lw=.8)
ax[0].plot([0.1, 60], [0.1, 60], "k-", lw=.6)
ax[0].set(xscale="log", yscale="log", xlim=(0.3, 40), ylim=(0.3, 80),
          xlabel="line S/N, stacked cube now", ylabel="expected S/N, adding own + archival data")
for k, c in cols.items():
    ax[0].scatter([], [], c=c, label=k)
ax[0].scatter([], [], c="w", edgecolor="k", marker="o", label="Hα"); ax[0].scatter([], [], c="w", edgecolor="k", marker="s", label="[OIII]/Hβ")
ax[0].legend(fontsize=8, loc="upper left")
flds = ["COSMOS", "CDFS", "SPT2040", "XMM2215"]
now = [np.median([r["t_stacked_h"] for r in out if r["field"] == f_]) for f_ in flds]
own = [np.median([max(r[f"t_own_{r['line']}_h"] - r["t_stacked_h"], 0) for r in out if r["field"] == f_]) for f_ in flds]
arc = [np.median([r[f"t_arch_{r['line']}_h"] for r in out if r["field"] == f_]) for f_ in flds]
ax[1].bar(flds, now, label="stacked now", color="#4c72b0")
ax[1].bar(flds, own, bottom=now, label="own, not yet stacked", color="#8fb0de")
ax[1].bar(flds, arc, bottom=np.add(now, own), label="archival (other programmes)", color="#c4c4c4")
ax[1].set_ylabel("median on-source hours per galaxy (same line)"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("kmos_council_snr.png", dpi=150)
print("\nwrote kmos_council_galaxies.csv, kmos_council_summary.md, kmos_council_snr.png")
