# Hybrid sky subtraction for OBs with SKY frames at another rotator angle (IZ)

**Date:** 2026-10-23
**Code:** `pyspark/sky/hybrid.py` (model), `pyspark/cubes/skysub_det.py` (DRP detector-level run),
`pyspark/cubes/spark_reconstruct.py` and `pyspark/sky/sky_subtraction.py` (option `[cubes] hybrid_sky`)
**Off by default; with `hybrid_sky = False` nothing changes (no extra DRP run, same code path).**

## Problem

Four IZ target lists are 2-3x noisier per sqrt(t) than cdfs_targets_2 and cosmos: cosmos_0.6_0.8_h21.4,
goods_0.6_0.8_h21.6_lir, cosmos_target_legac and cosmos_0.6_0.8_h21.7. They are exactly the lists whose SKY
frames were taken 37-90 deg away in rotator angle from the OBJECT frames (cdfs/cosmos: 0.5-3 deg). Conditions,
sky brightness, exposure bookkeeping and the det3 flat bug were ruled out.

pyspark reconstructs O and S separately (each with its own calibrations and its own OH wavelength correction)
and subtracts in cube space. With a rotator mismatch the same detector pixel lands in different voxels of O
and S, so warm, hot and noisy pixels no longer cancel: 3-5 % outlier voxels instead of 1 %, plus a pattern
shared by all O - S frames of a target (voxel correlation +0.54 between frames with different SKY frames, 0.00
in matched OBs). The shared pattern does not average down in the combination. In the same OBs, O - O
differences (same angle) are as clean as matched O - S.

## Method

The DRP's own detector-level subtraction (kmos_sci_red -no_subtract=FALSE -sky_tweak=FALSE) maps S with the
OBJECT's calibrations and OH wavelength correction, so detector features cancel. It leaves the OH-line
misalignment between the two exposures (flexure), which the separate reconstruction corrects. The hybrid keeps
both, per OBJECT frame and IFU:

    R_ds  = DRP detector-level O - S                     (skysub_det/, object's calibrations + lcorr)
    D     = (O - S) - R_ds = recon_O(S) - recon_S(S)     (linear reconstruction: S frame only, no target)
    R_hyb = R_ds + model(D)

`model(D)` per slitlet (along y for IFUs 1-16, along x for 17-24): the slitlet median of D in the OH-line
planes (pyspark's `find_lines` mask of the SKY cube, grown by one plane), elsewhere a least-squares fit of the
slitlet median of D with (a0 + a1 x) g + (b0 + b1 x) S + (c0 + c1 x) g', where S is the IFU median of the SKY
cube, g = dS/dlambda, x = (lambda - 0.93)/0.14, all high-passed with a 101-plane running median. Every term
comes from the SKY frame. R_hyb then goes through ZAPv2 exactly like the classic O - S (same SVD, same
continuum mask from the SKY cube).

## Pipeline

1. `reconstruct()` (unchanged main run), then, with `hybrid_sky = True`: pair every OBJECT frame with its
   nearest SKY frame (same DIT, closest DATE-OBS, as `match_obj_sky`) and run kmos_sci_red with
   `-no_subtract=FALSE -sky_tweak=FALSE -skip_sky_oh_align=TRUE -lcmethod=SLICES_MEDIAN` on the OBJECT frames
   whose SKY frame is more than `hybrid_min_dangle` (default 20 deg) away, plus all SKY frames, with the same
   calibrations, OH_SPEC, recipe_dir and lcorr_degree as the main run, in its own working directory. The DRP's
   object/sky table is checked against the pyspark pairing; each accepted cube goes to `<OB>/skysub_det/` with
   `ESO PYSPARK SKYSUB DET SKY` (the SKY frame) and `DANGLE`. Original frames are used even with
   `sky_cal_angle = object`.
2. `subtract_sky()`: for every OBJECT frame with a `skysub_det/` cube made with the same SKY frame, the hybrid
   residual replaces O - S before ZAP; the output is flagged `ESO PYSPARK SKY HYBRID = T`. All other frames
   (matched OBs, IFUs without a detector-level cube) get the classic O - S.

Only with the classic nearest-sky zapv2 path (no multi-sky, no K tbsub); `subtract_sky` raises otherwise.

```ini
[cubes]
hybrid_sky = True
hybrid_min_dangle = 20
# bands the hybrid is used for (default IZ)
hybrid_bands = IZ
# OH-line wings: IFU-wide median of D in the planes 2..3 from a line core (default 0 = as first validated)
hybrid_wing = 3
```
Server run: docs/procdir-hyb-server-run.md (scripts/kmos_server_hyb_IZ.cfg).

Cost: one extra kmos_sci_red run per mismatched OB (about the time of the main reconstruction, 4-5 min for a
9-frame IZ OB locally) and one skysub_det cube (34 MB) per mismatched OBJECT frame.

## Validation (local, IZ)

* **Implementation = prototype, bit for bit**: `subtract_sky(hybrid_sky=True)` on the worst IZ OB
  (KMOS.2026-02-12T06-02-50.540) gives the same `_s` cubes as the prototype chain (6/6 frames, 0 differing
  voxels, same NaN pattern); the full chain `reconstruct(hybrid_sky=True)` -> `subtract_sky` was checked
  stage by stage (sky/ = production reconstruction, skysub_det/ = prototype detector-level cubes).
* **Noise**, two mismatched OBs, 6-frame stacks after ZAP, hybrid/classic in the OH-free blue window (Hbeta):
  target k-sigma stack 0.60 / 0.53, signal-free jackknife 0.53 / 0.67, background apertures 0.44 / 0.53;
  red and OH windows 0.70-0.83. Matched control OB unchanged (0.96-1.01), so the option is applied only to
  mismatched frames. IFU 16 at 0.80-0.83 um: from 3.4-4.2x the other IFUs to 1.2-1.6x.
* **Flux**: the detector-level subtraction and the model are linear in the target signal. With the object's
  OH wavelength correction frozen (it is computed from the object frame alone and is the same in both
  products in production), 906 injected compact bumps (379 on edge spaxels, 122 below 0.833 um) and 241
  extended / point / edge / corner sources keep exactly the classic flux (ratio 1.00000, 0 voxels lost); the
  model changes by at most 6e-7 per voxel when signal is added. Cube-level lines survive ZAP identically
  (0.995). An earlier apparent loss for 12 bumps was a test artefact: the injected bumps changed the object
  frame's own wavelength correction.

* **End to end** (three local IZ OBs, pyspark chain -> esorex kmos_combine with the production settings ->
  the measurements of the OH-linear comparison notebook; WSL `~/impl_test/e2e`). Empirical noise of 3x3
  aperture sums on full-depth background spaxels, ratio to the classic chain:

  | OB (list, rotator offset) | hybrid, blue OH-free | blue OH | red OH-free | red OH |
  |---|---|---|---|---|
  | KMOS.2026-02-12T06:02:50.540 (cosmos_0.6_0.8_h21.4, 87-93 deg) | 0.39 | 0.55 | 0.69 | 0.78 |
  | KMOS.2026-01-14T05:23:13.608 (cosmos_target_legac, 37 deg) | 0.48 | 0.61 | 0.61 | 0.68 |
  | KMOS.2026-01-13T02:54:58.606 (cdfs_targets_2, matched; hybrid forced on) | 0.98 | 0.95 | 0.99 | 1.00 |

  With `hybrid_min_dangle = 20` the matched OB gets no detector-level cube at all (classic path). Two
  classic "[OIII] detections" of cosmos_0.6_0.8_h21.4 (S/N 5.7 and 10.7, 1224797 and 1170582) vanish with the
  hybrid: their line maps in the classic cube are positive/negative dipoles and streaks across the IFU, not a
  galaxy; in the hybrid cube the maps are flat. Real lines are kept: cdfs 23420 Halpha flux ratio 1.00
  (S/N 40 -> 46), [OIII] 1.01; 23840 Halpha 0.92 (within 1 sigma); the IFU-16 target 1168518 [OIII] S/N 4.1 -> 6.9.
  Row common mode on top: noise x0.92-0.95 for IFUs 9-24 in all three OBs.
* **Guards** (verified): a skysub_det cube is used only if it was made with the same SKY frame, is newer than
  the sky/ OBJECT cube and has a rotator offset above `hybrid_min_dangle` (with 95 deg all frames of the test
  OB fall back to classic, bit-identical to the classic chain); relative procdirs work.

Expected gain: per-galaxy Hbeta S/N about x2 (x1.6-2.6) for the 79 of 123 IZ galaxies in mismatched lists
(6-frame combinations; full combined cubes may gain more because the shared pattern grows with the number of
frames). The inverse-variance Hbeta stack of all galaxies gains x1.08-1.24 (already dominated by cdfs/cosmos).

## OH wings (`hybrid_wing`, 2026-10-27)

The first server run (PROCDIR_HYB, IZ) showed the remaining weakness: in the wings of OH lines and on faint OH
lines the hybrid left a residual that is the same in all slitlets and shared by the IFUs of a detector
(correlation of the per-plane background between target cubes 0.36-0.69 in OH-wing planes, ~0 classic), so red
target-aperture noise in those planes got 16-20 % worse in cosmos_0.6_0.8_h21.4 / h21.7 while everything else
improved.

**Cause.** model(D) uses the exact slitlet median of D only in the find_lines core grown by one plane; further out
it uses the smooth template fit (g, S, g' times a linear function of lambda, one set per slitlet). The fit is
dominated by the strongest lines and leaves most of the flank structure: from 2 planes from a core outward the
residual is 85-90 % of the slitlet median of D. On a single frame the flank residual is within the noise of a
slitlet median (a wider median mask or piecewise-linear coefficients only add noise or overfit), but it is common to
all slitlets and frames, so it survives in the combined cube.

**Fix.** In the planes 2..W from a core, use the high-passed IFU-wide median of D (~200 spaxels instead of the 14
of a slitlet: noise about /3.7, sqrt(14)) instead of the template fit; core planes and all other planes as before.
Every term still comes from D (SKY frame only), so the model stays independent of the target signal. The wing mask
is grown only from finite core planes and without wrap-around at the band ends.

**Choice of W = 3.** W = 5 also replaces planes 4-5 from a core; there the IFU median gains only above ~0.95 um
(red window: noise W5/W3 x0.995) and costs +0.6-0.7 % in the Hbeta/[OIII] windows of the test galaxies (pooled 5x5
apertures W5/W3 = 1.007 [0.983, 1.022] h21.4, 1.006 [1.002, 1.009] legac). W was chosen on the same two OBs it was
validated on (in-sample); the first server run should be checked with the same correlation test.

**Validation** (local, the two mismatched test OBs, row CM + hybrid, 6-frame esorex combinations; WSL
`~/impl_test/wing`, `~/impl_test/e2e/comb/<OB>/{C,W0,W3,W5}`; numbers current -> W3 unless noted):
- `hybrid_wing = 0` reproduces the previous output bit for bit (6/6 frames per OB). The hardened wing mask (finite
  cores, no wrap-around) gives the same W = 3 cubes bit for bit on both OBs (no such band-end case there).
- Cross-validated on D (fit on half of each slitlet's spaxels, scored on the other half, residual averaged over the
  slitlets of an IFU): between-IFU correlation of the wing residual 0.45 -> 0.01 (h21.4) and 0.32 -> 0.02 (legac),
  wing residual -17 % / -30 %.
- Combined cubes, between-IFU correlation of the per-plane background (same detector): h21.4 wing mid
  0.43 -> 0.03, wing red 0.36 -> 0.18 (W5 0.06), OH red 0.44 -> 0.15; legac wing mid/red 0.14/0.15 -> 0.00/0.05,
  OH red 0.34 -> 0.03 (classic ~0 everywhere).
- Target 5x5 aperture noise in wing planes: h21.4 x0.85 (mid), x0.80-0.88 (red); legac x0.98-1.00; OH-free
  planes and edge spaxels unchanged (x0.99-1.03 / x1.000), 5-sigma tails x0.92-0.99 in the red OH and wing planes.
  Hbeta/[OIII] windows of the test galaxies: x0.96 (h21.4), x0.99 (legac).
- Per detector (planes 2-3 from a core only): detectors 1 and 3 carry the shared residual and gain strongly
  (h21.4 x0.28-0.58 mid/red, x0.82-0.92 blue; legac detector 1 x0.47-0.81 mid/red); detector 2 has no shared
  residual and gets nothing (x0.98-1.02), in the blue wing planes of h21.4 slightly worse (x1.04, 5 cubes). Kept
  for all detectors: the loss is small against the gain and a per-detector switch would need a server-side test.
- Flux: in-field star continua x0.9986-0.9999; 1168518 [OIII] x1.04 at S/N 7.7.

## Not for H

H has the same O - S excess (O - S / O - O = 1.70 blue), but only about a quarter of its shared pattern comes
from SKY-frame detector features; the IZ model makes H 1.3-3.3x worse and every full-hybrid variant increases
the H common mode. Keep `hybrid_sky = False` for H (`sky_cal_angle = object` gives +5-11 % there).

Separate open issue found on the way: production ZAP removes part of the flux of H-band emission lines that
fall on strong OH lines (classic and hybrid alike; 16 % of such lines lose 16-24 % or more). Needs its own
investigation.

## YJ (2026-10-27): own model options

Test data: two mismatched YJ OBs (KMOS.2026-02-11T06:31:26.207 cosmos h21.4, KMOS.2026-06-08T05:24:25.056
SPT2040; scripts/fetch_yj_test.sh), local chain classic / row CM / row CM + hybrid (W = 3), 6-frame esorex
combinations (WSL `~/impl_test/e2e/comb/<OB>/YJ*`, variants of the sky subtraction in `~/impl_test/yjoh`).

**With the IZ model the hybrid helps off the OH lines but fails on them.** Interior noise in OH-free planes x0.23-0.62
of classic, but bright OH planes up to x1.7 in the interior and up to x3.1 / x3.6 on the outermost / second spaxel
from the cube edge, and a per-plane offset shared by the
IFUs of a detector (between-IFU covariance x2.6-130 of classic) that costs S/N in apertures (fits forced at the
catalogue z: 1147707 Halpha 10.9 -> 6.2, SPT2040 H45 [OIII] 6.7 -> 4.3).

**Cause.** In YJ the detector-level cube R_ds leaves OH residuals about 10x larger than in IZ (flexure at the
~90 deg rotator offset): in the line cores as large as the sky itself, sign-changing from plane to plane, shared by
all IFUs of a detector (correlation 0.9-1.0), and present at a lower level in the "free" planes. D = (O - S) - R_ds
carries the same with the opposite sign, and the model has to put it back. Two parts of the IZ model are not good
enough for that:
- outside the OH cores the model is high-passed (101-plane running median removed); with dense, large OH residuals
  in D that running median is not smooth, and its structure stays in the result as the shared per-plane offset;
- in the OH cores the slitlet median misses that D changes towards the slitlet ends (first/last 1-3 spaxels along
  the slitlet, the same for all slitlets of an IFU: a reconstruction edge effect), and the same pattern comes back
  in every frame (frame-to-frame correlation of the residual 0.73-0.86 in bright OH planes against 0.30-0.40
  classic).
The target cancels in D in YJ as in IZ (in-field stars: continuum in D +-0.02 %, wavelength offset R_ds vs O - S
0.002 planes), so every model of D stays target-independent.

**Options** (`[cubes]`, defaults = the IZ model, bit-identical to before):
- `hybrid_free_model = ifumed`: outside the OH cores the exact per-plane IFU median of D (no high-pass, ~200
  spaxels), plus the per-slitlet template fit of the slitlet deviations from it; in the wing planes the exact IFU
  median.
- `hybrid_oh_model = pos3`: in the OH-core planes the slitlet median plus a position profile along the slitlet at
  the 3 end positions: per plane, the median over the IFU's ~14 slitlets of (D - slitlet median) at each position.
  The median over slitlets keeps the common end pattern and rejects single-spaxel features of D (SKY-frame detector
  artefacts, i.e. the classic mismatch problem).
  Also tested and rejected: per-spaxel models along the slitlet (`run3`, `run5`, `lin`, `full` = classic in the cores,
  `run3b<T>` = run3 on bright lines only). They put D's per-spaxel structure back: `run3b30` made the central 5x5
  apertures x1.05-1.28 noisier in OH planes than the slitlet median (86-90 % of targets in the blue) and re-added
  1755 / 5523 voxels above 10 sigma (129 / 314 in the central 5x5); `full` 12052 / 26296. `pos3` re-adds 60 / 1991
  (none in the central 5x5) and removes 4636 / 41521 relative to the slitlet median. `pos2` behaves like `pos3`
  within 1-2 % except for a somewhat larger shared OH-plane covariance; `pos` (profile at all positions) is closest
  to classic in that covariance but adds noise at the target (central-aperture OH noise x0.96-1.05 of `pos3`, higher
  in 5 of 6 OB/window cases). Against the plain slitlet median (`median` + `ifumed`) `pos3` has the same central
  aperture noise but about 5 % lower forced line S/N at the target (11 of 13 lines with S/N >= 3; H45 [OIII] 8.7 ->
  7.5); its gain is at the cube edges, in the interior, in the tails, the artefact count and the shared covariance.

**Result** (median ratios to classic over the target cubes of both OBs; `ifumed` + `pos3`, W = 3; WSL
`~/impl_test/e2e/comb/<OB>/YJhP3I`):
- interior spaxels (>= 3 from the cube edge): free planes x0.22-0.54, wing x0.34-0.65, OH x0.54-0.87 (IZ model: OH up
  to x1.70);
- outermost spaxels (distance 0): free x0.52-0.83, OH x0.51-1.01 (IZ model: OH up to x3.1); second spaxel
  (distance 1): free x0.47-0.69, OH x0.59-1.08 (IZ model: OH up to x3.6);
- central 5x5 target apertures (robust): free x0.47-0.71, wing x0.65-0.84, OH x0.74-1.00, the same as with the
  slitlet median;
- 5-sigma tails, normalised with classic's per-spaxel sigma for all variants: OH planes x0.18-0.94, wing x0.05-0.64;
- between-IFU shared part (mean product of the per-plane background medians of IFU pairs on one detector): x0.6-2.2
  of classic where the classic value is clearly positive (OB 1 blue OH planes x4.5 of a small classic value);
  IZ model x2.6-130;
- lines, all fitted forced at the catalogue z with the default width (like for like): 1147707 Halpha S/N 10.9 ->
  11.2 (flux x1.06), H45 [OIII] 6.7 -> 7.5 (x0.90), 1123917 Halpha 5.2 -> 10.8 (x1.11), 1168518 Halpha 27 -> 44
  (x1.01), 1216183 Halpha 2.1 -> 4.8; lower: 1195224 Halpha 3.2 -> 2.8 and 1220318 3.7 -> 3.0 (flux x0.65 / x0.79,
  within the errors at S/N ~3, the same in every hybrid variant), and 1180013 Halpha 2.2 -> 0.7 (flux x0.18; in every
  hybrid variant 0.6-1.6, unchanged with row CM alone, so either a classic residual or a faint line the hybrid
  removes; undecided at 2 sigma); in-field star continua x0.996-1.003 (blue -0.35 % in all ifumed variants). (Free
  'anchored' fits are not comparable between variants: classic H45 was a broad, poor fit at sigma 345 km/s, and
  1120854 crossed the anchor threshold.)
- Row CM alone: cube noise x0.98-1.00, but forced line errors x0.75-0.89 on the h21.4 OB, so part of the line gains
  above comes from it (1168518: classic 27, row CM 34, row CM + hybrid 44).
- Chosen on these two OBs (in-sample) at 80-93 deg rotator offset. Not tested at 36-54 deg (XMM2215, 9 YJ OBs;
  cosmos_target_legac, cosmos_target_1), where R_ds should be closer to IZ.

**IZ keeps the IZ model.** On the two IZ test OBs `ifumed` changes little (interior and edges x1.00, remaining
red-wing correlation 0.18 -> 0.01) but raises the central-aperture noise in free planes by 3-8 %; `run3b30` makes
the IZ OH planes worse (x1.17-1.33); `pos3` is neutral in IZ (all metrics x0.97-1.03, OH-plane between-IFU
correlation of h21.4 0.13-0.15 -> 0.08-0.11). So the IZ configs keep the defaults (tmpl, median); `pos3` would be
harmless there but is not worth a re-run.

Server run: scripts/kmos_server_hyb_YJ.cfg (and _test), docs/procdir-hyb-server-run.md 'YJ'.

Investigation products: WSL `~/gains/{depth,hybrid_fix,hybrid_h}`, `~/impl_test/{hyb,hyb2,yjoh}`.
