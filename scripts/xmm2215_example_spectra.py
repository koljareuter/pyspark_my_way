"""Stack archival KMOS Phase-3 cubes of the XMM2215 example galaxies and compare
with our own stacked cube.

Per cube: 1.1"-diameter aperture spectrum + error from the STAT extension
(sigma per spaxel, rescaled so the line-free residuals have unit scatter).
Aperture centre: continuum peak of the cube (3-px border eroded, smoothed) if it
is significant (S/N >= 5) and within 1" of the WCS target position, else WCS.
Stack: inverse-variance weighted per pixel on a common grid (sky lines are
down-weighted automatically). Lines: weighted Gaussian fit with a linear
continuum -- Ha + [NII]6548,6583 (1:3) or [OIII]5007 + 4959 (1/3) or Hb;
common z (+/-400 km/s) and width (>= 50 km/s); strong sky-line pixels excluded.
Archival cubes are flux calibrated; ours is ADU/s, so plots are in S/N units.
Inputs: archive_xmm2215/index.csv (xmm2215_download_examples.py)
Outputs: xmm2215_example_spectra.pdf/.png, xmm2215_example_snr.csv
"""
import csv, glob, os, warnings
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_erosion, uniform_filter
from scipy.optimize import curve_fit

warnings.simplefilter("ignore")
C_KMS = 299792.458
AP_DIAM = 1.1        # arcsec
CENTRE_MAX = 1.0     # arcsec: max shift from WCS to continuum peak
EDGE_ERODE = 3       # px: IFU border masked when searching for the continuum peak (edge artefacts)
DV_MAX = 400.0       # km/s allowed offset from catalogue z
FIT_HALF = 0.020     # um half-width of the fit window
SKY_CLIP = 2.0       # pixels with error > SKY_CLIP x median are excluded from fits (and plots)
SIG_MIN = 50.0       # km/s: minimum intrinsic line width (stops fits locking onto single-pixel spikes)
REST = {"Ha": 0.65628, "NII83": 0.65835, "NII48": 0.65480, "OIII": 0.50069, "OIII59": 0.49590, "Hb": 0.48613}
GRID = {"YJ": np.arange(1.00, 1.36, 0.000175), "H": np.arange(1.45, 1.86, 0.0002)}
COMB_BASE = "kmos_proc_server_comb"
Z = {r["arm_target_id"]: r for r in csv.DictReader(open(os.path.join("scripts", "target_redshifts.csv"), encoding="utf-8"))
     if r["field"] == "XMM2215_EMPOWER.cat"}


def wave_axis(h):
    return h["CRVAL3"] + (np.arange(h["NAXIS3"]) - (h.get("CRPIX3", 1.0) - 1)) * (h.get("CDELT3") or h["CD3_3"])


def running_resid(f, k=101):
    pad = np.pad(f, k // 2, mode="edge")
    return f - np.array([np.nanmedian(pad[i:i + k]) for i in range(len(f))])


def extract(path, ra, dec):
    """Aperture spectrum, rescaled 1-sigma error, and how the centre was chosen."""
    with fits.open(path) as hd:
        cube = np.asarray(hd["DATA"].data, float)
        stat = np.asarray(hd["STAT"].data, float)
        h = hd["DATA"].header
    pix = abs(h.get("CD1_1") or h.get("CDELT1")) * 3600
    ny, nx = cube.shape[1:]
    x, y = WCS(h).celestial.world_to_pixel_values(ra, dec)
    mode = "wcs"
    img = np.nanmedian(cube[int(0.15 * len(cube)):int(0.85 * len(cube))], axis=0)
    inner = binary_erosion(np.isfinite(img), iterations=EDGE_ERODE)
    if inner.sum() > 10:
        sm = uniform_filter(np.nan_to_num(img), 3)
        v = img[inner]; noise = 1.4826 * np.median(np.abs(v - np.median(v))) / 3
        sm[~inner] = -np.inf
        py, px = np.unravel_index(np.argmax(sm), sm.shape)
        if (sm[py, px] - np.median(v)) / noise >= 5 and np.hypot(px - x, py - y) * pix <= CENTRE_MAX:
            x, y, mode = px, py, "continuum"
    if not (0 <= x < nx and 0 <= y < ny):
        x, y, mode = (nx - 1) / 2, (ny - 1) / 2, "centre"
    yy, xx = np.mgrid[:ny, :nx]
    ap = np.hypot(yy - y, xx - x) <= AP_DIAM / 2 / pix
    f = np.nansum(cube[:, ap], axis=1)
    e = np.sqrt(np.nansum(stat[:, ap] ** 2, axis=1))
    bad = np.isfinite(cube[:, ap]).sum(axis=1) < 0.5 * ap.sum()
    f[bad] = np.nan; e[bad] = np.nan
    e[e <= 0] = np.nan
    r = running_resid(f) / e
    k = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
    return wave_axis(h), f, e * (k if np.isfinite(k) and k > 0 else 1.0), mode


def stack(entries, band):
    g = GRID[band]
    num = np.zeros_like(g); den = np.zeros_like(g)
    for w, f, e, *_ in entries:
        fi = np.interp(g, w, f, left=np.nan, right=np.nan)
        ei = np.interp(g, w, e, left=np.nan, right=np.nan)
        ok = np.isfinite(fi) & np.isfinite(ei) & (ei > 0)
        num[ok] += fi[ok] / ei[ok] ** 2; den[ok] += 1 / ei[ok] ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return g, np.where(den > 0, num / den, np.nan), np.where(den > 0, 1 / np.sqrt(den), np.nan)


def fit_line(w, f, e, z, line):
    """Weighted Gaussian fit; returns (S/N of the main line flux, model, z_fit, line flux, sigma_kms)."""
    comps = {"Ha": [("Ha", 1.0), ("NII83", None), ("NII48", None)],
             "OIII": [("OIII", 1.0), ("OIII59", 1 / 3)], "Hb": [("Hb", 1.0)]}[line]
    lam0 = REST[line] * (1 + z)
    m = (np.abs(w - lam0) < FIT_HALF) & np.isfinite(f) & np.isfinite(e)
    m &= e <= SKY_CLIP * np.nanmedian(e[m])  # drop strong sky-line pixels: their residuals are systematic
    if m.sum() < 20:
        return np.nan, None, np.nan, np.nan, np.nan
    ww, ff, ee = w[m], f[m], e[m]
    inst = lam0 / 3500 / 2.355  # KMOS instrumental sigma (R~3500), um

    def model(x, c0, c1, zz, sig_kms, a, a_nii):
        s = np.hypot(REST[line] * (1 + zz) * sig_kms / C_KMS, inst)
        out = c0 + c1 * (x - lam0)
        for name, ratio in comps:
            amp = a * ratio if ratio is not None else (a_nii if name == "NII83" else a_nii / 3)
            out = out + amp * np.exp(-0.5 * ((x - REST[name] * (1 + zz)) / s) ** 2)
        return out

    dz = DV_MAX / C_KMS * (1 + z)
    lo = [-np.inf, -np.inf, z - dz, SIG_MIN, 0, 0]  # emission only: no negative (absorption-like) fits to sky residuals
    hi = [np.inf, np.inf, z + dz, 300, np.inf, np.inf]
    # several starting redshifts across the allowed range; keep the best chi^2
    # (a single start can lock onto a sky residual or a local minimum)
    best = None
    for z0 in np.linspace(z - 0.9 * dz, z + 0.9 * dz, 7):
        p0 = [np.nanmedian(ff), 0, z0, 100, np.nanmedian(ee), 0.1 * np.nanmedian(ee)]
        try:
            p_, cov_ = curve_fit(model, ww, ff, p0=p0, sigma=ee, absolute_sigma=True, bounds=(lo, hi), maxfev=20000)
        except Exception:
            continue
        c2 = np.sum(((ff - model(ww, *p_)) / ee) ** 2)
        if best is None or c2 < best[0]:
            best = (c2, p_, cov_)
    if best is None:
        return np.nan, None, np.nan, np.nan, np.nan
    _, p, cov = best
    chi2r = np.sum(((ff - model(ww, *p)) / ee) ** 2) / max(m.sum() - 6, 1)
    s_um = np.hypot(REST[line] * (1 + p[2]) * p[3] / C_KMS, inst)
    flux = p[4] * s_um * np.sqrt(2 * np.pi)
    # flux error from amplitude and width covariance
    J = np.array([0, 0, 0, p[4] * np.sqrt(2 * np.pi) * (REST[line] * (1 + p[2]) / C_KMS) ** 2 * p[3] / s_um,
                  s_um * np.sqrt(2 * np.pi), 0])
    err = np.sqrt(J @ cov @ J) * np.sqrt(max(chi2r, 1))
    return flux / err, (ww, model(ww, *p)), p[2], flux, p[3]


idx = list(csv.DictReader(open(os.path.join("archive_xmm2215", "index.csv"), encoding="utf-8")))
groups = defaultdict(list)
for r in idx:
    band = "H" if float(r["em_min"]) > 1.4e-6 else "YJ"
    t = Z[r["tid"]]
    groups[(r["tid"], band)].append((*extract(r["file"], float(t["ra_deg"]), float(t["dec_deg"])),
                                     r["proposal_id"], float(r["t_exptime"])))

rows, figs = [], []
for tid in sorted({k[0] for k in groups}):
    z = float(Z[tid]["redshift"]); t = Z[tid]
    own = glob.glob(os.path.join(COMB_BASE, "*", "*", f"DATACUBE_COMBINED_XMM2215_EMPOWER.cat_{tid}_YJ.fits"))
    own_spec = extract(own[0], float(t["ra_deg"]), float(t["dec_deg"])) if own else None
    t_own = fits.getheader(own[0])["TEXPTIME"] / 3600 if own else np.nan
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), gridspec_kw={"width_ratios": [2, 1]})
    fig.suptitle(f"XMM2215 {tid}   z={z}", fontweight="bold")
    for row, (band, lines) in enumerate((("YJ", ["OIII", "Hb"]), ("H", ["Ha"]))):
        ent = groups.get((tid, band), [])
        if not ent:
            axes[row, 0].set_axis_off(); axes[row, 1].set_axis_off(); continue
        g, s, se = stack(ent, band)
        t_arch = sum(x[5] for x in ent) / 3600
        progs = sorted({x[4] for x in ent})
        n_cont = sum(x[3] == "continuum" for x in ent)
        res = {ln: fit_line(g, s, se, z, ln) for ln in lines}
        res_own = {ln: fit_line(*own_spec[:3], z, ln) for ln in lines} if band == "YJ" and own_spec else {}
        for ln in lines:
            rows.append({"tid": tid, "band": band, "line": ln, "t_arch_h": round(t_arch, 2), "n_cubes": len(ent),
                         "n_continuum_centred": n_cont, "programmes": ";".join(progs),
                         "snr_arch": round(float(res[ln][0]), 1),
                         "snr_own": round(float(res_own[ln][0]), 1) if ln in res_own else "",
                         "t_own_h": round(t_own, 2) if ln in res_own else ""})
        ax = axes[row, 0]
        ax.plot(g, s / se, lw=.5, color="k", label=f"archival {t_arch:.1f} h, {len(ent)} cubes")
        for ln in lines:
            ax.axvline(REST[ln] * (1 + z), color="r", lw=.7, ls="--")
        ax.set(xlim=(g[0], g[-1]), ylim=(-5, 10), ylabel="flux / σ (per pixel)", xlabel="λ (µm)",
               title=f"{band} archival stack: " + ", ".join(progs))
        ax.legend(fontsize=8, loc="upper right")
        main = lines[0]; lam = REST[main] * (1 + z)
        ax = axes[row, 1]
        mz = np.abs(g - lam) < FIT_HALF
        sky = se > SKY_CLIP * np.nanmedian(se[mz])  # strong sky lines: masked in the plot
        ax.step(g[mz], np.where(sky, np.nan, s)[mz] * 1e19, where="mid", color="k", lw=.9,
                label=f"archival  S/N={res[main][0]:.1f}")
        if res[main][1] is not None:
            ax.plot(*res[main][1][0:1], res[main][1][1] * 1e19, color="r", lw=1.2)
        if main in res_own:
            ow, of, oe, _ = own_spec; mo = np.abs(ow - lam) < FIT_HALF
            # ours is ADU/s: scale to the archival noise level, for display only
            k = np.nanmedian(se[mz]) / np.nanmedian(oe[mo])
            osky = oe > SKY_CLIP * np.nanmedian(oe[mo])
            off = -6 * np.nanmedian(se[mz]) * 1e19
            ax.step(ow[mo], np.where(osky, np.nan, of * k)[mo] * 1e19 + off, where="mid", color="C1", lw=.9,
                    label=f"ours {t_own:.1f} h  S/N={res_own[main][0]:.1f} (scaled, offset)")
            if res_own[main][1] is not None:
                ax.plot(res_own[main][1][0], res_own[main][1][1] * k * 1e19 + off, color="C3", lw=1.2)
        for ln in (["Ha", "NII83", "NII48"] if main == "Ha" else ["OIII", "OIII59"]):
            ax.axvline(REST[ln] * (1 + z), color="grey", lw=.6, ls=":")
        ax.set(xlabel="λ (µm)", ylabel="flux (1e-19 erg/s/cm²/Å)",
               title=("Hα + [NII]" if main == "Ha" else "[OIII]4959, 5007") + " — sky lines masked")
        ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    figs.append(fig)

with PdfPages("xmm2215_example_spectra.pdf") as pdf:
    for f in figs:
        pdf.savefig(f)
for i, f in enumerate(figs):
    f.savefig(f"xmm2215_example_spectra_{i}.png", dpi=110)
with open("xmm2215_example_snr.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
