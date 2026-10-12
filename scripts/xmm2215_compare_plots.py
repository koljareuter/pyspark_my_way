"""'Ours' vs 'ours + archival' detail pages for XMM2215 example galaxies
(layout after the kmos_stacked_detail pages).

Ours: local combined cubes in kmos_proc_and_comb_server (YJ and H, ADU/s).
Archival: ESO Phase-3 KMOS cubes in archive_xmm2215/ (erg/s/cm2/A).
Our cubes are not flux calibrated and their continuum is largely removed by
the 2-D background step, so they are put on the archival flux scale by
assuming equal sensitivity per hour (noise in flux units ~ 1/sqrt(t)); the
Ha line-flux ratio is printed as a cross-check. Centring: the two pipelines' astrometry differs by ~0.5", so each data set is
centred on its own peak: our H cube on its Ha peak (our YJ cube uses the same
offset from the catalogue position); archival cubes per programme on the
peak of their stacked Ha (H) or continuum (YJ) map. Everything is then
stacked with inverse-variance weights. Line S/N: the profile (z within +/-400 km/s of the catalogue z, width) is found
once by a matched-filter scan of the deepest spectrum (ours + archival); [OIII]
uses the Ha redshift and only its width is scanned. Both
'ours' and 'ours + archival' are then measured with that fixed profile.
Maps: per-spaxel S/N of line / continuum images, every cube resampled onto a
common 0.2" north-up grid centred on its own aperture centre.

Usage: python scripts/xmm2215_compare_plots.py [tid ...]   (default 32400 29284)
Outputs: xmm2215_compare_<tid>.pdf/.png, xmm2215_compare_snr.csv
"""
import csv, glob, os, sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from scipy.ndimage import map_coordinates, median_filter

# shared helpers: extract, stack, fit_line, wave_axis, REST, GRID, Z, constants
_src = open(os.path.join("scripts", "xmm2215_example_spectra.py"), encoding="utf-8").read()
exec(_src.split("idx = list(")[0])

OWN_BASE = "kmos_proc_and_comb_server"
MAP_PIX, MAP_HALF = 0.2, 1.6          # arcsec: common map grid
CONT_GAP, CONT_WIDTH = 0.004, 0.010   # um: continuum windows either side of the line
LINE_HALF_KMS = 350.0                 # half-width of the line image window
COL_OURS, COL_ALL = "#1f77b4", "#d62728"
LINE_COL = {"Ha": "#d62728", "OIII": "#17becf", "Hb": "#1f77b4"}
LINE_LAB = {"Ha": r"H$\alpha$", "OIII": r"[OIII]$\lambda$5007", "Hb": r"H$\beta$"}
COMPANIONS = {"Ha": [("NII48", "[NII]6548"), ("NII83", "[NII]6583")], "OIII": [("OIII59", "[OIII]4959")]}


def load(path, ra, dec, line_lam=None, fixed=False):
    """Aperture spectrum (+ rescaled error), cube, stat, aperture centre, CD matrix.

    Centre: archival cubes as in extract() (continuum peak or WCS). For our own
    cubes pass line_lam: the 2-D background step removes most of the continuum,
    so we centre on the peak of the line S/N map within CENTRE_MAX of the WCS.
    """
    w, f, e, mode = extract(path, ra, dec)
    with fits.open(path) as hd:
        cube = np.asarray(hd["DATA"].data, float)
        stat = np.asarray(hd["STAT"].data, float)
        h = hd["DATA"].header
    from astropy.wcs import WCS
    from scipy.ndimage import binary_erosion, uniform_filter
    cd = np.array([[h.get("CD1_1", 0.0), h.get("CD1_2", 0.0)], [h.get("CD2_1", 0.0), h.get("CD2_2", 0.0)]])
    pix = abs(h["CD1_1"]) * 3600 if h.get("CD1_1") else abs(h["CDELT1"]) * 3600
    x, y = WCS(h).celestial.world_to_pixel_values(ra, dec)
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    inner = binary_erosion(np.isfinite(np.nanmedian(cube, axis=0)), iterations=EDGE_ERODE)
    near = inner & (np.hypot(xx - x, yy - y) * pix <= CENTRE_MAX)
    c = dict(w=wave_axis(h), cube=cube, stat=stat, cd=cd)
    if line_lam is not None:
        img, err = line_image(c, line_lam, line_lam * LINE_HALF_KMS / C_KMS)
        if img is not None:
            snr = uniform_filter(np.nan_to_num(img), 3) / np.sqrt(uniform_filter(np.nan_to_num(err ** 2), 3) / 9)
            snr[~near] = -np.inf
            py, px = np.unravel_index(np.argmax(snr), snr.shape)
            if snr[py, px] >= 3:
                x, y, mode = float(px), float(py), "line"
            else:
                mode = "wcs"
    elif fixed:
        mode = "wcs"  # (ra, dec) already is the Ha centre measured in our H cube
    elif mode == "continuum":
        img = np.nanmedian(cube[int(0.15 * len(cube)):int(0.85 * len(cube))], axis=0)
        sm = uniform_filter(np.nan_to_num(img), 3); sm[~inner] = -np.inf
        y, x = np.unravel_index(np.argmax(sm), sm.shape)
    elif mode == "centre":
        y, x = (ny - 1) / 2, (nx - 1) / 2
    if mode in ("line", "wcs"):  # re-extract at the chosen centre (extract() only knows WCS/continuum)
        ap = np.hypot(yy - y, xx - x) <= AP_DIAM / 2 / pix
        f = np.nansum(cube[:, ap], axis=1)
        e = np.sqrt(np.nansum(stat[:, ap] ** 2, axis=1))
        bad = np.isfinite(cube[:, ap]).sum(axis=1) < 0.5 * ap.sum()
        f[bad] = np.nan; e[bad] = np.nan; e[e <= 0] = np.nan
        r = running_resid(f) / e
        k = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        e = e * (k if np.isfinite(k) and k > 0 else 1.0)
    c.update(w=w, f=f, e=e, x=float(x), y=float(y), mode=mode,
             radec=tuple(float(v) for v in WCS(h).celestial.pixel_to_world_values(x, y)))
    return c


def line_image(c, lam, dl):
    """Continuum-subtracted line image and its 1-sigma map (both in cube units)."""
    w, cube, stat = c["w"], c["cube"], c["stat"]
    sig = np.abs(w - lam) <= dl
    cont = (np.abs(w - lam) >= dl + CONT_GAP) & (np.abs(w - lam) <= dl + CONT_GAP + CONT_WIDTH)
    # drop strong sky-line channels (high STAT) from both windows
    s_ch = np.nanmedian(stat, axis=(1, 2))
    ok = s_ch <= SKY_CLIP * np.nanmedian(s_ch[sig | cont])
    sig &= ok; cont &= ok
    if sig.sum() < 2 or cont.sum() < 4:
        return None, None
    cmed = np.nanmedian(cube[cont], axis=0)
    img = np.nansum(cube[sig] - cmed, axis=0)
    err = np.sqrt(np.nansum(stat[sig] ** 2, axis=0))
    return img, err


def cont_image(c):
    w, cube, stat = c["w"], c["cube"], c["stat"]
    s_ch = np.nanmedian(stat, axis=(1, 2))
    ok = (s_ch <= SKY_CLIP * np.nanmedian(s_ch)) & (np.arange(len(w)) > 0.1 * len(w)) & (np.arange(len(w)) < 0.9 * len(w))
    img = np.nanmean(cube[ok], axis=0)
    err = np.sqrt(np.nansum(stat[ok] ** 2, axis=0)) / ok.sum()
    return img, err


def resample(c, img, err):
    """Map an image + error onto the common grid (offsets from the aperture centre).
    The IFU border (EDGE_ERODE px) is masked first: it carries edge artefacts."""
    from scipy.ndimage import binary_erosion
    inner = binary_erosion(np.isfinite(np.nanmedian(c["cube"], axis=0)), iterations=EDGE_ERODE)
    img = np.where(inner, img, np.nan); err = np.where(inner, err, np.nan)
    n = int(round(2 * MAP_HALF / MAP_PIX)) + 1
    off = (np.arange(n) - (n - 1) / 2) * MAP_PIX / 3600  # deg
    east, north = np.meshgrid(off, off)
    # CD maps pixel offsets to (east, north) intermediate coordinates in degrees
    inv = np.linalg.inv(c["cd"])
    dx = inv[0, 0] * east + inv[0, 1] * north
    dy = inv[1, 0] * east + inv[1, 1] * north
    xs, ys = c["x"] + dx, c["y"] + dy
    out = map_coordinates(np.nan_to_num(img, nan=np.nan), [ys, xs], order=1, cval=np.nan)
    oute = map_coordinates(np.nan_to_num(err, nan=np.nan), [ys, xs], order=1, cval=np.nan)
    return out, oute


def combine_maps(items):
    num = den = None
    for im, er in items:
        ok = np.isfinite(im) & np.isfinite(er) & (er > 0)
        if num is None:
            num = np.zeros_like(im); den = np.zeros_like(im)
        num[ok] += im[ok] / er[ok] ** 2; den[ok] += 1 / er[ok] ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(den > 0, num / den, np.nan); s = np.where(den > 0, 1 / np.sqrt(den), np.nan)
    return m / s  # S/N map


def flux_scale(ours, arch_stack, t_ours, t_arch, lam):
    """Factor putting our ADU/s onto the archival flux scale, assuming the same
    instrument delivers the same noise per hour: sigma_flux ~ 1/sqrt(t).
    Uses the median per-pixel error of clean (non-sky) pixels within 0.05 um of the line."""
    g, s, se = arch_stack
    oe = np.interp(g, ours["w"], ours["e"], left=np.nan, right=np.nan)
    m = (np.abs(g - lam) < 0.05) & np.isfinite(se) & np.isfinite(oe)
    m &= (se <= SKY_CLIP * np.nanmedian(se[m])) & (oe <= SKY_CLIP * np.nanmedian(oe[m]))
    return np.nanmedian(se[m]) * np.sqrt(t_arch / t_ours) / np.nanmedian(oe[m])


def matched_snr(g, s, se, z_l, sig_kms, line):
    """Line flux and S/N for a FIXED profile (z, width): optimal (matched-filter)
    estimate with per-pixel errors, after a linear continuum fit to line-free pixels.
    Errors are inflated by sqrt(chi2_red) of the line-free residuals if > 1."""
    lam0 = REST[line] * (1 + z_l)
    m = (np.abs(g - lam0) < FIT_HALF) & np.isfinite(s) & np.isfinite(se)
    m &= se <= SKY_CLIP * np.nanmedian(se[m])
    inst = lam0 / 3500 / 2.355
    names = {"Ha": ["Ha", "NII83", "NII48"], "OIII": ["OIII", "OIII59"], "Hb": ["Hb"]}[line]
    width = np.hypot(lam0 * sig_kms / C_KMS, inst)
    free = m.copy()
    for n_ in names:
        free &= np.abs(g - REST[n_] * (1 + z_l)) > 3 * width
    if free.sum() < 10:
        return np.nan, np.nan, None
    cf = np.polyfit(g[free] - lam0, s[free], 1, w=1 / se[free])
    r = s - np.polyval(cf, g - lam0)
    chi2r = np.sum((r[free] / se[free]) ** 2) / max(free.sum() - 2, 1)
    P = np.exp(-0.5 * ((g - lam0) / width) ** 2)
    if line == "OIII":
        P = P + np.exp(-0.5 * ((g - REST["OIII59"] * (1 + z_l)) / width) ** 2) / 3
    P = P / (width * np.sqrt(2 * np.pi))  # unit integrated flux of the main line
    mm = m & (P > 1e-3 * P.max())
    F = np.sum(r[mm] * P[mm] / se[mm] ** 2) / np.sum(P[mm] ** 2 / se[mm] ** 2)
    eF = 1 / np.sqrt(np.sum(P[mm] ** 2 / se[mm] ** 2)) * np.sqrt(max(chi2r, 1))
    wm = np.abs(g - lam0) < FIT_HALF
    return F / eF, F, (g[wm], np.polyval(cf, g[wm] - lam0) + F * P[wm])


def profile_scan(spec, z, line, dv=20.0, widths=(50, 75, 100, 150, 200, 250)):
    """Matched-filter search for the line profile: z on a 20 km/s grid within
    +/-ZSEARCH of the catalogue z (sigma = 100 km/s), then the width at that z."""
    zs = z + np.arange(-DV_MAX, DV_MAX + dv, dv) / C_KMS * (1 + z) if dv < DV_MAX else np.array([z])
    sn = [matched_snr(*spec, zz, 100.0, line)[0] for zz in zs]
    z_best = zs[int(np.nanargmax(sn))] if np.isfinite(sn).any() else z
    sw = [matched_snr(*spec, z_best, w, line)[0] for w in widths]
    return z_best, widths[int(np.nanargmax(sw))] if np.isfinite(sw).any() else 100.0


def sky_offset(ra0, dec0, ra1, dec1):
    """(east, north) arcsec of (ra1, dec1) relative to (ra0, dec0)."""
    return ((ra1 - ra0) * np.cos(np.radians(dec0)) * 3600, (dec1 - dec0) * 3600)


def shift(ra, dec, east, north):
    return ra + east / 3600 / np.cos(np.radians(dec)), dec + north / 3600


def stacked_peak(cubes, kind, lam, min_snr=4.0):
    """Offset (east, north; arcsec) of the peak of the stacked line/continuum S/N map,
    all cubes aligned on their WCS. Zero if nothing significant within CENTRE_MAX."""
    items = []
    for c in cubes:
        im, er = cont_image(c) if kind == "cont" else line_image(c, lam, lam * LINE_HALF_KMS / C_KMS)
        if im is not None:
            items.append(resample(c, im, er))
    if not items:
        return (0.0, 0.0)
    m = combine_maps(items)
    from scipy.ndimage import uniform_filter
    sm = uniform_filter(np.nan_to_num(m), 3) * 3  # ~S/N of a 3x3 box
    n = m.shape[0]; off = (np.arange(n) - (n - 1) / 2) * MAP_PIX
    E, N = np.meshgrid(off, off)
    sm[np.hypot(E, N) > CENTRE_MAX] = -np.inf
    r, c_ = np.unravel_index(np.argmax(sm), sm.shape)
    return (float(E[r, c_]), float(N[r, c_])) if sm[r, c_] >= min_snr else (0.0, 0.0)


if os.environ.get("AP_DIAM"):
    AP_DIAM = float(os.environ["AP_DIAM"])
tids = sys.argv[1:] or ["32400", "29284"]
idx = [r for r in csv.DictReader(open(os.path.join("archive_xmm2215", "index.csv"), encoding="utf-8"))
       if r["tid"] in tids]
rows = []
for tid in tids:
    t = Z[tid]; ra, dec, z = float(t["ra_deg"]), float(t["dec_deg"]), float(t["redshift"])
    data = {}
    for band in ("H", "YJ"):
        own_f = [p for p in glob.glob(os.path.join(OWN_BASE, "*", "*", f"DATACUBE_COMBINED_XMM2215_EMPOWER.cat_{tid}.fits"))
                 if (fits.getheader(p).get("ESO INS GRAT1 ID") or "").strip() == band]
        arch = [r for r in idx if r["tid"] == tid and (float(r["em_min"]) > 1.4e-6) == (band == "H")]
        if not own_f or not arch:
            continue
        main = "Ha" if band == "H" else "OIII"
        lam = REST[main] * (1 + z)
        # ours: Ha peak (H); in YJ the same offset from the catalogue position as found in our H cube
        if band == "H":
            ours = load(own_f[0], ra, dec, line_lam=lam)
            own_off = sky_offset(ra, dec, *ours["radec"]) if ours["mode"] == "line" else (0.0, 0.0)
        else:
            ours = load(own_f[0], *shift(ra, dec, *own_off), fixed=True)
        t_ours = fits.getheader(own_f[0])["TEXPTIME"] / 3600
        # archival: per programme, peak of the stacked Ha (H) or continuum (YJ) map near the catalogue position
        arch_c, offsets = [], {}
        for prog in sorted({r["proposal_id"] for r in arch}):
            rs = [r for r in arch if r["proposal_id"] == prog]
            cs = [load(r["file"], ra, dec, fixed=True) for r in rs]
            off = stacked_peak(cs, "line" if band == "H" else "cont", lam)
            offsets[prog] = off
            cs = [load(r["file"], *shift(ra, dec, *off), fixed=True) for r in rs]
            for c, r in zip(cs, rs):
                c["t"] = float(r["t_exptime"]); c["prog"] = prog
            arch_c += cs
        t_arch = sum(c["t"] for c in arch_c) / 3600
        a_stack = stack([(c["w"], c["f"], c["e"]) for c in arch_c], band)
        k = flux_scale(ours, a_stack, t_ours, t_arch, lam)
        ours_s = (ours["w"], ours["f"] * k, ours["e"] * k)
        data[band] = dict(
            ours=stack([ours_s], band), all=stack([ours_s] + [(c["w"], c["f"], c["e"]) for c in arch_c], band),
            arch=a_stack, ours_c=ours, arch_c=arch_c, k=k, t_ours=t_ours, t_arch=t_arch,
            progs=sorted(offsets))
        print(f"{tid} {band}: ours centred on {ours['mode']} (offset E,N = {own_off[0]:+.2f},{own_off[1]:+.2f}\"); "
              "archival offsets " + ", ".join(f"{p_} {o[0]:+.2f},{o[1]:+.2f}" for p_, o in offsets.items()))
        print(f"    flux scale {k:.3e} (noise per hour); ours {t_ours:.1f} h + archival {t_arch:.1f} h ({len(arch_c)} cubes)")

    # line measurements: profile (z, width) from the deepest spectrum, then the
    # same fixed profile measured in both 'ours' and 'all'
    res = defaultdict(dict)
    z_ha = None
    for band, ln in (("H", "Ha"), ("YJ", "OIII")):
        if band not in data:
            continue
        d = data[band]
        if ln == "Ha" or z_ha is None:
            z_l, sig = profile_scan(d["all"], z, ln)
            if ln == "Ha":
                z_ha = z_l
        else:  # [OIII]: same galaxy -> Ha redshift; only the width is scanned
            z_l, sig = profile_scan(d["all"], z_ha, ln, dv=1e9)
        for which in ("ours", "all"):
            snr, flux, mod = matched_snr(*d[which], z_l, sig, ln)
            res[which][ln] = (snr, mod, z_l, flux, sig)
        # cross-check of the flux scale: line-flux ratio archival / ours (same profile; only if both detect it)
        sa, fa, _ = matched_snr(*d["arch"], z_l, sig, ln)
        so, fo = res["ours"][ln][0], res["ours"][ln][3]
        d["k_line"] = d["k"] * fa / fo if sa >= 5 and so >= 5 else np.nan
        d["snr_arch"] = sa
        print(f"{tid} {ln}: S/N ours {so:.1f}, archival only {sa:.1f}, combined {res['all'][ln][0]:.1f}; "
              f"flux scale {d['k']:.3e} (noise per hour) vs {d['k_line']:.3e} (line-flux ratio)")
        rows.append({"tid": tid, "band": band, "line": ln, "t_ours_h": round(d["t_ours"], 2),
                     "t_arch_h": round(d["t_arch"], 2), "snr_ours": round(float(res["ours"][ln][0]), 1),
                     "snr_all": round(float(res["all"][ln][0]), 1), "snr_arch_only": round(float(d["snr_arch"]), 1),
                     "gain_measured": round(float(res["all"][ln][0] / res["ours"][ln][0]), 2),
                     "gain_sqrt_t": round(float(np.sqrt((d["t_ours"] + d["t_arch"]) / d["t_ours"])), 2),
                     "z_line": round(float(z_l), 5), "sigma_kms": round(float(sig), 0),
                     "flux_scale": f"{d['k']:.3e}", "flux_scale_line_check": f"{d['k_line']:.3e}",
                     "archival_programmes": ";".join(d["progs"])})

    # ── figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(17, 11.5), constrained_layout=True)
    tt = {b: (data[b]["t_ours"], data[b]["t_arch"]) for b in data}
    fig.suptitle(
        f"XMM2215 {tid}  ·  z = {z:.4f}  ·  left: our data  ·  right: ours + archival  ·  "
        + "  ·  ".join(f"{b}: {o:.1f} h + {a:.1f} h" for b, (o, a) in tt.items()),
        fontsize=11, fontweight="bold")
    outer = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.6, 1.5, 1.4])
    cols = (("ours", "Our data", COL_OURS), ("all", "Ours + archival", COL_ALL))

    # row 0: H-band overview (falls back to YJ)
    ob = "H" if "H" in data else "YJ"
    ylims = None
    for ci, (which, lab, colr) in enumerate(cols):
        ax = fig.add_subplot(outer[0, ci])
        g, s, se = data[ob][which]
        sky = se > SKY_CLIP * np.nanmedian(se)
        y = median_filter(np.where(sky, np.nan, s), 3, mode="nearest") * 1e19
        ax.plot(g * 1e4, y, color=colr, lw=.6)
        if ylims is None:
            v = y[np.isfinite(y)]; ylims = (np.nanpercentile(v, 1) - 5, max(np.nanpercentile(v, 99.7), 5) * 1.3)
        ax.set_ylim(ylims)
        for ln, (key, nm) in {"Ha": ("Ha", "Hα"), "NII83": ("NII83", "[NII]"), "OIII": ("OIII", "[OIII]"), "Hb": ("Hb", "Hβ")}.items():
            lam = REST[key] * (1 + z)
            if g[0] < lam < g[-1]:
                ax.axvline(lam * 1e4, color=LINE_COL.get(ln, "#ff7f0e"), lw=.9, ls="--", alpha=.7)
                ax.text(lam * 1e4, .97, nm, transform=ax.get_xaxis_transform(), rotation=90, ha="right", va="top",
                        fontsize=7, color=LINE_COL.get(ln, "#ff7f0e"))
        tot = data[ob]["t_ours"] + (data[ob]["t_arch"] if which == "all" else 0)
        ax.set_title(f"{lab} — {ob} band, {tot:.1f} h on source (sky lines masked, 3-px median)",
                     fontsize=9, color=colr, fontweight="bold")
        ax.set_xlabel(r"$\lambda_{\rm obs}$ (Å)", fontsize=8); ax.set_ylabel(r"flux ($10^{-19}$ erg/s/cm$^2$/Å)", fontsize=8)
        ax.tick_params(labelsize=7); ax.set_xlim(g[0] * 1e4, g[-1] * 1e4)

    # row 1: line zooms (Ha, [OIII]) with fits
    zoom_y = {}
    for ci, (which, lab, colr) in enumerate(cols):
        inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1, ci])
        for j, (band, ln) in enumerate((("H", "Ha"), ("YJ", "OIII"))):
            ax = fig.add_subplot(inner[j])
            if band not in data:
                ax.set_axis_off(); continue
            g, s, se = data[band][which]
            lam = REST[ln] * (1 + res[which][ln][2])
            m = np.abs(g - lam) < FIT_HALF
            sky = se > SKY_CLIP * np.nanmedian(se[m])
            ax.step(g[m] * 1e4, np.where(sky, np.nan, s)[m] * 1e19, where="mid", color=colr, lw=1.1)
            ax.fill_between(g[m] * 1e4, -se[m] * 1e19, se[m] * 1e19, step="mid", color="grey", alpha=.15, lw=0)
            snr, mod, zf = res[which][ln][:3]
            if mod is not None:
                ax.plot(mod[0] * 1e4, mod[1] * 1e19, color="k", lw=1.2)
            ax.axvline(lam * 1e4, color=LINE_COL[ln], lw=1.3, ls="--", alpha=.8)
            for key, nm in COMPANIONS.get(ln, []):
                ax.axvline(REST[key] * (1 + res[which][ln][2]) * 1e4, color="#ff7f0e", lw=.8, ls=":")
            if (band, ln) not in zoom_y:
                v = (np.where(sky, np.nan, s)[m] * 1e19); v = v[np.isfinite(v)]
                lo, hi = np.nanpercentile(v, 1), np.nanmax(v)
                zoom_y[(band, ln)] = (lo - 0.25 * (hi - lo), hi + 0.25 * (hi - lo))
            ax.set_ylim(zoom_y[(band, ln)])
            ax.set_xlim((lam - FIT_HALF) * 1e4, (lam + FIT_HALF) * 1e4)
            tot = data[band]["t_ours"] + (data[band]["t_arch"] if which == "all" else 0)
            ax.set_title(f"{LINE_LAB[ln]}   S/N = {snr:.1f}   ({tot:.1f} h)", fontsize=9, color=LINE_COL[ln])
            ax.set_xlabel(r"$\lambda_{\rm obs}$ (Å)", fontsize=7); ax.tick_params(labelsize=7)

    # row 2: S/N maps (white light, Ha, [OIII]) on a common grid
    vmax = {}
    for ci, (which, lab, colr) in enumerate(cols):
        inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2, ci])
        for j, (kind, band, ln, cmap) in enumerate((("cont", "H" if "H" in data else "YJ", None, "viridis"),
                                                    ("line", "H", "Ha", "inferno"), ("line", "YJ", "OIII", "cividis"))):
            ax = fig.add_subplot(inner[j])
            if band not in data:
                ax.set_axis_off(); continue
            d = data[band]
            cubes = [(d["ours_c"], d["k"])] + ([(c, 1.0) for c in d["arch_c"]] if which == "all" else [])
            items = []
            for c, k in cubes:
                if kind == "cont":
                    im, er = cont_image(c)
                else:
                    zf = res["all"][ln][2] if np.isfinite(res["all"][ln][2]) else z
                    lam = REST[ln] * (1 + zf)
                    im, er = line_image(c, lam, lam * LINE_HALF_KMS / C_KMS)
                if im is None:
                    continue
                items.append(resample(c, im * k, er * k))
            if not items:
                ax.set_axis_off(); continue
            snr_map = combine_maps(items)
            key = (kind, band)
            if key not in vmax:
                vmax[key] = max(np.nanpercentile(snr_map, 99), 3)
            ext = [MAP_HALF + MAP_PIX / 2, -MAP_HALF - MAP_PIX / 2, -MAP_HALF - MAP_PIX / 2, MAP_HALF + MAP_PIX / 2]
            # grid columns run west->east; flip so east is left (sky convention)
            imh = ax.imshow(snr_map[:, ::-1], origin="lower", cmap=cmap, vmin=-1, vmax=vmax[key], extent=ext)
            ax.add_patch(plt.Circle((0, 0), AP_DIAM / 2, fill=False, color="w", lw=1.2, ls="--"))
            ax.set_title(("continuum" if kind == "cont" else LINE_LAB[ln]) + " — S/N per spaxel", fontsize=8)
            ax.set_xlabel("ΔRA (″)", fontsize=7); ax.set_ylabel("ΔDec (″)", fontsize=7); ax.tick_params(labelsize=6)
            plt.colorbar(imh, ax=ax, fraction=0.046, pad=0.02).ax.tick_params(labelsize=6)

    # improvement box
    txt = []
    for r in [r for r in rows if r["tid"] == tid]:
        txt.append(f"{LINE_LAB[r['line']]}: S/N {r['snr_ours']:.1f} → {r['snr_all']:.1f}  "
                   f"(×{r['gain_measured']:.1f} measured, ×{r['gain_sqrt_t']:.1f} expected from √t)   "
                   f"+{r['t_arch_h']:.1f} h archival from {r['archival_programmes'].replace(';', ', ')}")
    fig.text(0.5, -0.01, "\n".join(txt), ha="center", va="top", fontsize=10,
             bbox=dict(fc="#f5f5f5", ec="#999", boxstyle="round,pad=0.5"))
    fig.savefig(f"xmm2215_compare_{tid}.pdf", bbox_inches="tight")
    fig.savefig(f"xmm2215_compare_{tid}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)

with open("xmm2215_compare_snr.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
