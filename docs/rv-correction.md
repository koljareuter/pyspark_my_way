# Barycentric (rv) correction and other existing steps: measured (2026-11-03)

Local test data only: 7 OBs with 6 OBJECT frames each (IZ 2026-01-13T02:54, 01-14T05:23, 02-12T06:02; YJ 01-13T04:40,
02-11T06:31, 06-08T05:24; H 06-07T07:24), local esorex `kmos_combine` with the production settings, line fluxes from
forced fits at fixed kinematics. Every result was checked by an independent second analysis.

## rv correction

`velocity_correct` resamples every spaxel of every frame onto the barycentric grid with linear interpolation
(`rvcorr_cube`, np.interp; before 2026-09-24 a damped sinc). Because |v_bary| is below one pixel, the shift fraction is
f = |v|/pixel and the interpolation smooths by sqrt(f(1-f)) pixels: 19-26 km/s for the OBs with |v| = 18-20 km/s
(larger at the blue end), about 14 km/s for |v| = 4.5 km/s (pixel 45-49 km/s); in the server combinations
(kmos_comb_hybrid, 275 cubes) a median of 23 km/s (IZ) and 19 km/s (YJ).

What it costs, and what not:
- flux is conserved, and for broad lines (sigma >= 60 km/s) the true flux S/N too;
- neighbouring planes become correlated: lag-1 correlation about 0.5 instead of 0.15 and per-pixel sigma about 1.28x
  lower for |v| = 19 km/s (plane-to-plane difference noise 1.65-1.75x lower; only 1.13-1.17x for |v| = 4.5 km/s). An S/N
  that treats planes as independent, or uses a per-pixel STAT, is too high by about 1.3x; `line_snr.measure_lines`
  (noise_from_cube, side-band rescaling, autocorrelation factor) is calibrated in both cases (reported/true error
  1.05-1.13 in injections);
- lines are broadened by 19-26 km/s in quadrature: the peak drops by 7-11 % (sigma 50 km/s) or 10-16 % (40 km/s), and
  forced fits with a width taken from elsewhere (another band, LSF only) lose up to 6-8 % of the flux;
- for lines at the instrumental width (sigma_obs 35-50 km/s) the production fit loses 2-5 % true S/N (injections in the
  |v| = 19 km/s OBs: raw 4-8 % minus a 0-4 % method floor); an optimal, correlation-aware estimator would lose nothing,
  since the interpolation kernel is invertible for f != 0.5;
- 0.05 % of the voxels turn NaN (np.interp spreads NaNs by one plane).
The frame-to-frame wavelength alignment itself is fine: O frames agree to 0.1-3 km/s after the OH correction (worst
IFU 9 km/s; the old +-60-90 km/s per-IFU offsets are gone), v_bary of production agrees with astropy to 0.005 km/s.

Alternative, `[cubes] rv_mode = header` (new option; default `interp` unchanged and bit-identical, other values refused):
the data stay in the telescope frame and only VHELIO (and RV MODE = header) are stored. The frames are combined in the
telescope frame (kmos_combine ignores the spectral WCS and combines plane by plane) and the combined cube is relabelled
lambda_bary = lambda (1 + <v>/c) afterwards. The frames of one combination differ in v_bary by only 0.2 km/s (YJ) /
0.4 km/s (IZ) rms in the median, at most 4-12 km/s, so this is sharper than the interpolation in all 275 server
combinations (for multi-night combinations with std(v) > 2 km/s the relabelling is exact only on average).
Tested on the 7 OBs (same frames, same combination; T = header mode + relabelling, I = interpolation): line flux T/I
0.995 (IZ) and 0.999 (YJ); true S/N (null test, 21-42 line-free offsets) T/I 0.99 (95 % 0.94-1.11) over the 12 lines of
the |v| = 19 km/s OBs and 0.96 (0.92-1.07) over all 18 lines (the two |v| = 4.5 km/s OBs are not informative), i.e.
unchanged to about +-5-10 %, which cannot resolve the 2-5 % narrow-line effect; widths of the anchored lines 1.3-1.4
km/s narrower (quadrature expectation 1.5-3 km/s for 60-130 km/s lines). The two H lines (S/N 4-5, width pinned in T)
are not used. So the gain is honest per-pixel noise, sharper narrow lines and peaks, a few per cent S/N on lines at the
instrumental width, and no interpolation of either data or telluric transmission; broad integrated lines are not
detected more often.

Adoption (2026-11-03): `rv_mode = header` is set in the PROCDIR_HYBE cfgs of all three bands
(docs/procdir-hyb-server-run.md 'Slitlet-end + header-mode re-run'). Done:
- analysis project `scripts/telluric_combined.py`: the index stores RV MODE per frame (indexes must be rebuilt from the
  new frames); v = 0 for header-mode frames in all shifted arrays of method best and in frame_transmission (standard
  alignment d kept); cubes whose frames mix modes, header cubes with an index without RV MODE, interp cubes with a
  header index and already relabelled cubes are refused; after dividing out T_eff, CRVAL3 and CDELT3/CD3_3 of DATA and
  STAT (and the primary WAVELMIN/WAVELMAX) are scaled by (1 + <v>/c), <v> = plain mean of the index VHELIO over the
  PROV frames of that IFU; new keys RV MODE = relabelled, VBARY MEAN/STD/NOFRAME, SPECSYS = BARYCENT (extensions and
  primary); a warning when the VHELIO spread exceeds 2 km/s. Interp cubes: output bit-identical to the old script
  (84/84 files).
- end-to-end test (H OB 2026-06-07, YJ OB 2026-06-08; velocity_correct I/T, esorex kmos_combine, telluric apply):
  VBARY MEAN equals the frame mean, relabelled grid = interp grid x (1 + <v>/c), star spectra T vs I shifted by
  0-1.3 km/s, star continua T/I 0.997-1.000, galaxy forced fluxes T/I 0.96 (H) / 1.04 (YJ) for the 3 lines with
  S/N >= 4 (median (F_T - F_I)/error over all 76 forced lines 0.01), fit errors 1.4-1.5 % larger in T (honest per-pixel
  noise). Telluric residuals equal at matched resolution; on T's own unsmoothed grid the residual rms in absorption
  bands is 1.15-1.3x larger relative to the noise (the interpolation no longer smooths them).
- analysis code: kmos_sample.load_galaxy_band_stacks / load_stack_with_stat resample members whose grids differ by
  0.01-70 km/s onto the first member's grid (with their TELLEFF arrays; equal grids untouched, existing stacks
  bit-identical on 459 stacks); new helpers is_relabelled, telescope_wavelength_axis, align_to_grid. Adapted:
  telluric_eval.py, kmos_hbeta_detection, kmos_iz_blue_end_checks, kmos_best3_single_gauss, kmos_xmm2215_31179,
  pipeline explainer steps 8 and 10 (not re-run).
Still open:
- telluric constants: in-field k in YJ 1.116 (T) vs 1.044 (I) on the same stars, so the band default 1.075 (and IZ
  1.13) under-corrects header-mode cubes by about 5 % in line depth. Re-derive K_HF_DEFAULT from the in-field run
  values of the first PROCDIR_HYBE apply run, re-check SCI_REF (16-star validation) and SYS_A/SUBRES_C, before
  science use;
- EDPS must keep 'HIERARCH ESO PYSPARK RV MODE' in the combined products (it keeps PYSPARK VHELIO; check on the first
  server product, the telluric script refuses header cubes without it);
- notebooks and scripts that read single _rv frames (kmos_statistics, kmos_snr_annuli*, the frame part of
  kmos_showcase, pyspark_my_way scripts/telluric_compare.py) take the frame axis as barycentric; for header-mode
  frames it is wave (1 + VHELIO/c);
- noise calibrations (null tests) on the new cubes; multi-night relabelling (VBARY STD > 2 km/s) untested end to end.
Scripts: WSL `~/impl_test/rv` (variants, combinations) and `~/impl_test/rv/scripts` (rv_variants.py, rv_measure.py,
multi_night_v.py, rv_mode_test.py); end-to-end test WSL `~/impl_test/rvh/e2e` with scripts in the session scratchpad
rvh/e2e/; review scripts in the session scratchpad review_rv/.

## Other existing steps

- Spatial alignment (`kmos_combine --method=header`): the header offsets leave 0.02-0.05" rms misalignment within an
  OB (a drift common to all IFUs, partly atmospheric dispersion). It costs a point source about 1.7 % of its flux at
  r = 0.34" and compact galaxies < 1 %; re-registering on stars recovers about 1 % with no S/N gain. Keep `header`.
- Rejection in `kmos_combine` (ksigma 3/3/3): no measurable line or continuum flux is clipped (milder or no rejection
  changes forced fluxes by +0.3 % median, star continua by 0.000 %); median and min_max bias fluxes low by 1-5 %. Keep
  it. With about 6 frames it removes only a small part of the spikes (one clipping rule at n = 6).
- Row common mode on detector 1 (production: dets 2,3, minpix 150, which leaves det 1 practically uncorrected):
  dets 1,2,3 with minpix 30 lowers the combined-cube aperture noise of IFUs 1-8 to x0.92 (3x3) / x0.90 (5x5) in the IZ
  OB 2026-01-13T02:54 (about +9 % line S/N), but only to x0.99 in the YJ OB 2026-02-11T06:31 (det 1 has about 52 usable
  pixels per row there, too few). IFUs 9-24 bit-identical, no line flux lost.
  More OBs (2026-11-03; IZ 2026-02-12 and 2026-01-14 with hybrid sky, YJ 2026-06-08, H 2026-06-07 with the bright star
  H129R on det 1; WSL ~/impl_test/rvh/rowcm_det1): det-1 aperture noise IZ x0.87-0.96 (5-13 % depending on OB and
  estimator; x0.85 at 0.81 um, about 1.0 at 1.08 um, OH planes unchanged), YJ x0.97-0.99, H x1.00. IFUs 9-24
  bit-identical in every OB. At grow 1 det 1 has at least 49-54 (IZ), 38 (YJ), 30 (H) usable pixels per row in the 7
  local XCALs, so minpix 30 corrects every row (minpix 150 only 4-10 rows per XCAL); grow 0 leaks the bright star into
  the row estimate (keep grow 1). At grow 1 no significant leakage (star flux x0.9994-1.0005; other det-1 targets'
  continuum after local background within 1.1 sigma; a few 4-sigma shape projections occur for det 2/3 under the
  production settings too). The old -5 sigma forced-fit drop of 1170582 Halpha (YJ 2026-02-11) is a fit artefact (direct
  window flux +0.1 sigma). Adopt for IZ (cfg kmos_server_hybr_IZ.cfg, full chain), optional for YJ, not for H.
  Possible further gain: with an odd number of frames the OB-median residuals contain exact zeros, so 5-45 % of rows
  get c = 0 (also on det 2/3 in production); a tie-robust row estimator would correct more.
