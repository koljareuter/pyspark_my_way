"""Headline numbers for the council presentation.

Uses the combined (stacked) cubes in kmos_proc_server_comb/ and the ESO-archive
IFU tables (COSMOS, ECDFS) to compute, per galaxy:
  * on-source integration time in our programme
  * integrated line S/N (Ha where in band, else [OIII]5007 / Hb)
  * archival KMOS time on the same position (other programmes), split by
    whether the archival band covers Ha at the galaxy redshift
  * projected Ha S/N after adding archival data (S/N scales as sqrt(t))

Run from the repo root:  python scripts/council_numbers.py
"""
import csv, glob, os, re, warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

COMB_BASE    = "kmos_proc_server_comb"
OB_BASE      = "kmos_ob_proc_server"
REDSHIFT_CSV = os.path.join("scripts", "target_redshifts.csv")
ARCH_TABLES  = [os.path.join("scripts", "table_ifus_COSMOS_Field_exptime.fits"),
                os.path.join("scripts", "table_ifus_ECDFS 2.fits")]
OUR_PROG     = "116.29EG"
OUT_CSV      = "kmos_council_numbers.csv"
OUT_MD       = "kmos_council_summary.md"
OUT_PNG      = "kmos_council_snr.png"

LINES = {"Ha": 0.65628, "OIII": 0.50069, "Hb": 0.48613}
# rest-frame continuum windows (um), chosen to avoid [NII], [SII], [OIII]4959, Hb
CONT_REST = {"Ha":   [(0.6440, 0.6530), (0.6610, 0.6700)],
             "OIII": [(0.5030, 0.5150)],
             "Hb":   [(0.4780, 0.4840), (0.4890, 0.4940)]}
# rest wavelengths masked when estimating the empirical noise
MASK_REST = [0.48613, 0.49590, 0.50069, 0.65480, 0.65628, 0.65835, 0.67164, 0.67308]
BAND_RANGES = {"IZ": (0.780, 1.090), "YJ": (1.000, 1.359), "H": (1.456, 1.850),
               "K": (1.934, 2.450), "HK": (1.456, 2.450)}
VWIN_KMS   = 350.0   # +/- line integration window
ZSEARCH_KMS = 400.0  # allowed offset of the line peak from the catalogue z
AP_RADIUS  = 3.5     # px (0.159"/px -> ~1.1" diameter aperture)
MAX_OFFSET = 4.0     # px; line peak farther than this from centre -> use centre
MATCH_ARCSEC = 1.0
C_KMS = 299792.458


def wave_axis(h):
    n = h["NAXIS3"]
    d = h.get("CDELT3") or h.get("CD3_3")
    return h["CRVAL3"] + (np.arange(n) - (h.get("CRPIX3", 1.0) - 1)) * d


def in_band(band, obs_wl, pad=0.003):
    lo, hi = BAND_RANGES.get(band, (0, 0))
    return lo + pad < obs_wl < hi - pad


def line_snr(cube, wave, z, line):
    """Integrated line S/N in a circular aperture, empirical noise."""
    lam = LINES[line] * (1 + z)
    dl = lam * VWIN_KMS / C_KMS
    sig = np.abs(wave - lam) <= dl
    if sig.sum() < 3:
        return np.nan, np.nan, None

    # aperture centre: peak of the line image (smoothed), else IFU centre
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        img = np.nansum(cube[sig], axis=0) - sig.sum() * np.nanmedian(cube, axis=0)
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    sm = np.nan_to_num(img)
    sm = (sm + np.roll(sm, 1, 0) + np.roll(sm, -1, 0) + np.roll(sm, 1, 1) + np.roll(sm, -1, 1)) / 5
    py, px = np.unravel_index(np.argmax(sm), sm.shape)
    if np.hypot(py - cy, px - cx) > MAX_OFFSET:
        py, px = cy, cx
    ap = np.hypot(yy - py, xx - px) <= AP_RADIUS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = np.nansum(cube[:, ap], axis=1)
        nval = np.isfinite(cube[:, ap]).sum(axis=1)
    spec[nval < 0.5 * ap.sum()] = np.nan

    # local continuum (linear fit in rest-frame windows)
    cmask = np.zeros_like(wave, bool)
    for a, b in CONT_REST[line]:
        cmask |= (wave >= a * (1 + z)) & (wave <= b * (1 + z))
    cmask &= np.isfinite(spec)
    if cmask.sum() < 10:
        return np.nan, np.nan, (py, px)
    p = np.polyfit(wave[cmask], spec[cmask], 1)
    resid = spec - np.polyval(p, wave)
    # re-centre on the line peak (spec-z can be off by a few 100 km/s)
    srch = np.abs(wave - lam) <= lam * ZSEARCH_KMS / C_KMS
    smr = np.convolve(np.nan_to_num(resid), np.ones(3) / 3, mode="same")
    lam_pk = wave[srch][np.argmax(smr[srch])]
    sig = np.abs(wave - lam_pk) <= dl
    flux = np.nansum(resid[sig])

    # empirical noise: scatter of same-width window sums over line-free spectrum
    # (includes sky-line residuals and spatial/spectral correlated noise)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        k = 151
        pad = np.pad(spec, k // 2, mode="edge")
        run_med = np.array([np.nanmedian(pad[i:i + k]) for i in range(len(spec))])
    r2 = spec - run_med
    bad = ~np.isfinite(r2)
    for lr in MASK_REST:
        bad |= np.abs(wave - lr * (1 + z)) <= 2 * dl
    w = int(sig.sum())
    sums = [np.sum(r2[i:i + w]) for i in range(0, len(r2) - w, max(w // 2, 1))
            if not bad[i:i + w].any()]
    if len(sums) < 20:
        return np.nan, np.nan, (py, px)
    sums = np.asarray(sums)
    noise = 1.4826 * np.median(np.abs(sums - np.median(sums)))
    return flux / noise, flux, (py, px)


# ── redshifts ────────────────────────────────────────────────────────────────
ZCAT = {}
with open(REDSHIFT_CSV, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        try:
            ZCAT[(r["field"].strip(), r["arm_target_id"].strip())] = (
                float(r["redshift"]), r["redshift_type"].strip())
        except ValueError:
            pass

# ── archive (other programmes only, science-length frames) ───────────────────
arch = vstack([Table.read(t) for t in ARCH_TABLES], metadata_conflicts="silent")
arch["s_run_id"] = [str(s).strip() for s in arch["s_run_id"]]
arch["s_filt_name"] = [str(s).strip() for s in arch["s_filt_name"]]
arch = arch[[not s.startswith(OUR_PROG) for s in arch["s_run_id"]]]
arch = arch[np.asarray(arch["s_exptime"], float) >= 60]
a_ra = np.radians(np.asarray(arch["s_ra"], float))
a_de = np.radians(np.asarray(arch["s_dec"], float))


def archive_for(ra, dec, z):
    r, d = np.radians(ra), np.radians(dec)
    cosd = np.sin(d) * np.sin(a_de) + np.cos(d) * np.cos(a_de) * np.cos(a_ra - r)
    sep = np.degrees(np.arccos(np.clip(cosd, -1, 1))) * 3600
    m = arch[sep < MATCH_ARCSEC]
    out = {"t_any": float(np.sum(m["s_exptime"])), "t_ha": 0.0,
           "progs_any": sorted(set(m["s_run_id"])), "progs_ha": set(), "bands": {}}
    for row in m:
        b = row["s_filt_name"]
        out["bands"][b] = out["bands"].get(b, 0.0) + float(row["s_exptime"])
        if in_band(b, LINES["Ha"] * (1 + z)):
            out["t_ha"] += float(row["s_exptime"])
            out["progs_ha"].add(row["s_run_id"])
    out["progs_ha"] = sorted(out["progs_ha"])
    return out


# ── per-galaxy measurements from the combined cubes ──────────────────────────
rows = []
for path in sorted(glob.glob(os.path.join(COMB_BASE, "*", "*", "DATACUBE_COMBINED_*.fits"))):
    m = re.match(r"DATACUBE_COMBINED_(.+)_([^_]+)_(IZ|YJ|HK|H|K)\.fits$", os.path.basename(path))
    if not m:
        continue
    field, tid, band = m.groups()
    if field.endswith("_R") or tid.startswith("recenter"):
        continue  # reference / acquisition stars
    with fits.open(path) as hd:
        h0 = hd[0].header
        cube = np.asarray(hd["DATA"].data, float)
        wave = wave_axis(hd["DATA"].header)
    z, ztype = ZCAT.get((field, tid), (np.nan, ""))
    rec = {"field": field, "tid": tid, "band": band, "z": z, "z_type": ztype,
           "ra": h0["RA"], "dec": h0["DEC"], "ncombine": h0["NCOMBINE"],
           "t_ours_s": float(h0["TEXPTIME"]), "seeing": h0.get("ESO TEL IA FWHM"),
           "obs_name": h0.get("ESO OBS NAME")}
    # usable range of this cube (band edges are noisy -> generous pad)
    wlo, whi = h0["WAVELMIN"] / 1e3 + 0.030, h0["WAVELMAX"] / 1e3 - 0.012
    for ln in LINES:
        snr = flux = np.nan
        if np.isfinite(z) and wlo < LINES[ln] * (1 + z) < whi:
            snr, flux, _ = line_snr(cube, wave, z, ln)
        rec[f"snr_{ln}"] = snr
    rec["ha_in_band"] = bool(np.isfinite(z) and wlo < LINES["Ha"] * (1 + z) < whi)

    a = archive_for(rec["ra"], rec["dec"], z if np.isfinite(z) else 0)
    rec["t_arch_any_s"] = a["t_any"]
    rec["t_arch_ha_s"] = a["t_ha"]
    rec["arch_progs_ha"] = ";".join(a["progs_ha"])
    rec["arch_bands"] = ";".join(f"{b}:{t/3600:.1f}h" for b, t in sorted(a["bands"].items()))
    if rec["ha_in_band"] and np.isfinite(rec["snr_Ha"]):
        rec["snr_Ha_proj"] = rec["snr_Ha"] * np.sqrt((rec["t_ours_s"] + a["t_ha"]) / rec["t_ours_s"])
    else:
        rec["snr_Ha_proj"] = np.nan
    rows.append(rec)
    print(f"{field[:20]:20s} {tid:>8s} {band} z={z:.3f}  t={rec['t_ours_s']/3600:5.2f}h "
          f"SNR Ha={rec['snr_Ha']:6.1f} OIII={rec['snr_OIII']:6.1f} Hb={rec['snr_Hb']:6.1f} "
          f"| arch Ha-band {a['t_ha']/3600:5.1f}h (any {a['t_any']/3600:5.1f}h) -> {rec['snr_Ha_proj']:6.1f}")

keys = list(rows[0].keys())
with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6f}" if k in ("ra", "dec") else f"{v:.3f}" if isinstance(v, float) else v)
                    for k, v in r.items()})

# ── programme-wide telescope time from the per-OB folders ────────────────────
ob_rows = []
for ob in sorted(glob.glob(os.path.join(OB_BASE, "KMOS.*"))):
    sci = sorted(glob.glob(os.path.join(ob, "science", "*_rv_tell.fits"))) or \
          sorted(glob.glob(os.path.join(ob, "science", "*.fits")))
    if not sci:
        continue
    try:
        h = fits.getheader(sci[0])
    except Exception:
        continue
    ob_rows.append((os.path.basename(ob), h.get("ESO OBS NAME"), h.get("ESO OBS TARG NAME"),
                    float(h.get("ESO OBS EXECTIME", np.nan)), len(sci)))
with open("kmos_council_obs.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["ob_folder", "obs_name", "targ_name", "exectime_s", "n_sci_frames"])
    w.writerows(ob_rows)
print(f"\n{len(rows)} galaxies written to {OUT_CSV}; {len(ob_rows)} OB folders scanned")
