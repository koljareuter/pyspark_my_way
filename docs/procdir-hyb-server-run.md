# Server re-reduction with the October 2026 noise fixes (PROCDIR_HYB)

**What changes against PROCDIR_OHLIN** (same OBs, same patched DRP, same OH_SPEC and lcorr degrees):

| band | new | expected gain | doc |
|---|---|---|---|
| IZ | hybrid sky subtraction for OBJECT frames whose SKY frame is > 20 deg away in rotator angle | 3x3-aperture noise of the combined cubes x0.39-0.48 (OH-free Hbeta window), x0.61-0.69 (red) for the mismatched lists; matched OBs untouched | hybrid-sky.md |
| IZ | row common-mode correction, detectors 2 and 3 (dcr_ frames) | noise x0.92-0.95 for IFUs 9-24 | noise-improvements-2026-10.md |
| H | SKY frames with the OBJECT frames' calibration angle (`sky_cal_angle = object`) | +5-11 % per frame | noise-improvements-2026-10.md |
| YJ | hybrid with the YJ options `ifumed` / `pos3` + row common mode (added 2026-10-27, section 'YJ' at the end) | interior noise x0.22-0.65 off the OH lines and x0.54-0.87 on them, central apertures x0.47-0.84 / x0.74-1.00, mismatched lists | hybrid-sky.md, 'YJ' |

`snap_wcs` (integer-shift combination) is not used: after the hybrid it raised the aperture noise by 5-14 %
in the end-to-end test. The H results of the first run were not adopted (edge strips); the IZ OH-wing fix is the
re-run PROCDIR_HYBW (section 'OH-wing re-run').
PROCDIR, PROCDIR_OHLIN and kmos_ob_sorted are only read; everything new goes to PROCDIR_HYB(TEST) and PROCDIR_HYBW.

Order (2026-10-28): the first run (IZ + H, next section) is done. Pending are the OH-wing re-run and YJ; start with
'Pending runs: commit and pull first'.

## First run IZ + H on hqempower1 (done 2026-10-25/26, do not re-run)

Results checked on 2026-10-26; downloaded to kmos_comb_hyb, telluric-corrected into kmos_comb_hyb_tell (that folder
also holds byte copies of the old OHLIN YJ cubes). Running these steps again would rewrite the dcr_ frames and the sky/
and skysub_det/ cubes of PROCDIR_HYB, which PROCDIR_HYBW links to.

```bash
cd ~/pyspark_my_way && git pull            # after the changes are committed and pushed

# 0. smoke test on one OB (about 20-30 min); compare with the expectations in the cfg header
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_IZ_test.cfg
mkdir -p /local_storage/KOLJA/hyb_work/IZtest && cd /local_storage/KOLJA/hyb_work/IZtest
env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
    -p $HOME/pyspark_my_way/scripts/kmos_server_hyb_IZ_test.cfg > hyb_IZtest.out 2>&1
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_IZ_test.cfg --status
grep -h "Hybrid sky\|classic O - S" /local_storage/KOLJA/PROCDIR_HYBTEST/KMOS.2026-02-12T06:02:50.540/*_IZ_*.log

# 1. links for the full run (dry run first)
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_IZ.cfg scripts/kmos_server_hyb_H.cfg --dry-run
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_IZ.cfg scripts/kmos_server_hyb_H.cfg

# 2. one working directory per band (kmos_sci_red writes into the cwd)
mkdir -p /local_storage/KOLJA/hyb_work/IZ && cd /local_storage/KOLJA/hyb_work/IZ
nohup env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
      -p $HOME/pyspark_my_way/scripts/kmos_server_hyb_IZ.cfg > hyb_IZ.out 2>&1 &
mkdir -p /local_storage/KOLJA/hyb_work/H && cd /local_storage/KOLJA/hyb_work/H
nohup env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
      -p $HOME/pyspark_my_way/scripts/kmos_server_hyb_H.cfg > hyb_H.out 2>&1 &

# 3. progress (dcr = row-corrected frames, det = skysub_det cubes of the hybrid; 0 for matched OBs)
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_IZ.cfg scripts/kmos_server_hyb_H.cfg --status

# 4. combination exactly as for PROCDIR_OHLIN (docs/oh-lcorr-linear-rereduction.md, step 6), new folders:
for B in IZ H; do
  python scripts/make_edps_input.py $B --procdir /local_storage/KOLJA/PROCDIR_HYB --out /local_storage/KOLJA/edps_in_hyb_$B
done
mkdir -p /home/kreuter/EDPS_data_hyb
sed -i.bak 's#^base_dir=.*#base_dir=/home/kreuter/EDPS_data_hyb#' ~/.edps/application.properties
#    restart edps-gui --port 5007, run cubes_combination on edps_in_hyb_<BAND>, export to /local_storage/KOLJA/res_hyb_<BAND>
python scripts/edps_pick_jobs.py --base /home/kreuter/EDPS_data_hyb/KMOS/cubes_combination --band IZ    # H alike
```

Time: the hybrid adds a second kmos_sci_red run per mismatched IZ OB (about the length of the reconstruction), so
IZ takes roughly 1.5x the PROCDIR_OHLIN run (≈ 20 h), H as before (≈ 3 h). Disk: as PROCDIR_OHLIN plus
dcr_ frames (≈ 50 MB per IZ frame, ≈ 40 GB) and skysub_det cubes (34 MB per mismatched IZ OBJECT frame).

## First run: locally afterwards (done)

As written for the first run; it was done into kmos_comb_hyb and kmos_comb_hyb_tell. For the pending runs: the download
in WSL, the telluric step with the Windows Python (the WSL python has no astropy).

```bash
REMOTE_BASE=/local_storage/KOLJA/res_hyb_IZ LOCAL_BASE=kmos_comb_hyb_IZ bash scripts/download_combined.sh    # H alike
python scripts/telluric_combined.py apply --comb kmos_comb_hyb_IZ --index telluric_frame_index_ohlin.json \
       --tell telluric_all --out kmos_comb_hyb_IZ_tell
```

The telluric frame index is keyed by frame timestamps, so the PROCDIR_OHLIN index also covers the dcr_ frames.
Since 2026-10-25, `telluric_combined.py` defaults to method `best`: per-night wavelength alignment of the standards
(IZ, YJ), an in-field line-depth factor, a 0.02 floor and flux-scale keywords. See the analysis project's
docs/telluric-correction-v2.md. Its constants were measured on PROCDIR_OHLIN cubes and stay valid here, because the
wavelength calibration is the same.

## First run: what to check in the results (checked 2026-10-26)

* Mismatched IZ lists (cosmos_0.6_0.8_h21.4, goods_0.6_0.8_h21.6_lir, cosmos_target_legac, cosmos_0.6_0.8_h21.7):
  noise in the Hbeta window about x0.4-0.6 of PROCDIR_OHLIN at equal frames; spurious "lines" made of
  positive/negative stripes in the classic cubes disappear (locally two [OIII] "detections" at S/N 5.7 and
  10.7 of cosmos_0.6_0.8_h21.4 were such artefacts; the line maps show dipoles and streaks, no galaxy).
* cdfs_targets_2 / cosmos IZ: only the row common mode, a few per cent for IFUs 9-24.
* Real lines keep their flux (locally cdfs 23420: Halpha flux ratio 1.00, [OIII] within 1 sigma).
* The single frames of the smoke-test OB should equal the local ones (WSL ~/impl_test/e2e/procdir_C) up to
  the sky subtraction if the dc frames and calibrations are the same.

## Pending runs: commit and pull first (2026-10-28)

The OH-wing re-run and YJ need the same uncommitted changes. In the local checkout, commit and push together:
pyspark/sky/hybrid.py, pyspark/sky/sky_subtraction.py, scripts/prepare_procdir_ohlin.py, scripts/make_edps_input.py
(never committed so far), scripts/fetch_yj_test.sh, scripts/kmos_server_hyb_YJ.cfg, scripts/kmos_server_hyb_YJ_test.cfg,
scripts/kmos_server_hybw_IZ.cfg, scripts/kmos_server_hybw_IZ_test.cfg, docs/hybrid-sky.md, docs/procdir-hyb-server-run.md.

The server has an untracked hand copy of scripts/make_edps_input.py (step 4 of the first run used it). `git pull` stops
on it once the file is in git, so move it aside before this one pull (and diff it against the pulled version if you
edited it there):
```bash
cd ~/pyspark_my_way
mv scripts/make_edps_input.py ~/make_edps_input.py.handcopy     # only before this first pull
git pull && git log -1 --oneline                                 # must show the new commit
```

Then:
- The two reductions are independent and may run at the same time or one after the other, in either order: the OH-wing
  re-run only reads the IZ OB folders of PROCDIR_HYB (links) and writes PROCDIR_HYBW; YJ adds new YJ OB folders to
  PROCDIR_HYB; each has its own working directory.
- The two EDPS combinations share edps-gui (port 5007) and ~/.edps/application.properties: run them one after the
  other. Finish one cubes_combination, its export and edps_pick_jobs before switching base_dir for the other.
- Combine only once a full run has finished. The --status summary counts an OB 'with final cubes' as soon as it has
  one, so wait until the nohup job has ended (pgrep prints nothing) and check its .out file, as written in each section.
- Locally, download in WSL; run telluric_combined.py with the Windows Python in the same folder (the WSL python has no
  astropy).

## OH-wing re-run (PROCDIR_HYBW, `hybrid_wing = 3`, 2026-10-27)

The first PROCDIR_HYB run left a residual in the OH-line wings that is shared by the IFUs of a detector
(hybrid-sky.md, 'OH wings'). The fix changes only the sky subtraction, so PROCDIR_HYBW reuses the dcr_ frames and
the sky/ and skysub_det/ cubes of PROCDIR_HYB (links straight to the real files, only read) and runs
subtract_sky -> illumination -> 2D background -> velocity correction again. scripts/kmos_server_hybw_IZ.cfg is
kmos_server_hyb_IZ.cfg with procdir PROCDIR_HYBW, `reconstruct_cubes = False` and `row_commonmode_correct = False`;
`prepare_procdir_ohlin.py --reuse` refuses a cfg that would run the reconstruction, a detector step or a calibration
step, because those would write onto the links or mix calibrations. The IZ/H OB folders of PROCDIR_HYB stay as they
are, as the reference (the YJ run only adds new YJ OB folders there). Do
not write-protect it: the science/ copies inherit the file mode of the linked sky/ cubes and are then updated in place.

```bash
cd ~/pyspark_my_way && git log -1 --oneline     # must show the commit from 'Pending runs' (else commit, push, pull)

# 0. smoke test on the worst mismatched OB (cosmos_0.6_0.8_h21.4), same procdir as the full run, a few minutes
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ_test.cfg --reuse /local_storage/KOLJA/PROCDIR_HYB
mkdir -p /local_storage/KOLJA/hyb_work/IZwtest && cd /local_storage/KOLJA/hyb_work/IZwtest
env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
    -p $HOME/pyspark_my_way/scripts/kmos_server_hybw_IZ_test.cfg > hybw_IZtest.out 2>&1
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ_test.cfg --status
#    expected: dc 9, dcr 9, sky 9, det 6, rv 6, wing 6xW3, OH_SPEC and degree 1 marked * (from PROCDIR_HYB)
grep -h "Hybrid sky subtraction\|classic O - S" $(ls -t /local_storage/KOLJA/PROCDIR_HYBW/KMOS.2026-02-12T06:02:50.540/*_IZ_*.log | head -1)
#    expected: six '-- Hybrid sky subtraction with skysub_det/...' lines, no 'classic O - S' (newest OB log only: every
#    run adds a log, and the full run redoes this OB with the same files)
grep -h "reconstruction" /local_storage/KOLJA/hyb_work/IZwtest/hybw_IZtest.out
#    expected: '- Running reconstruction' and '-- Finished reconstruction' in the same second (no kmos_sci_red)

# 1. all 83 IZ OBs (dry run first; a NOTE line means PROCDIR_HYB has a different number of sky/ cubes than dc frames
#    there, normally fewer; the new run then uses the frames that have sky/ cubes)
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ.cfg --reuse /local_storage/KOLJA/PROCDIR_HYB --dry-run
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ.cfg --reuse /local_storage/KOLJA/PROCDIR_HYB
mkdir -p /local_storage/KOLJA/hyb_work/IZw && cd /local_storage/KOLJA/hyb_work/IZw
nohup env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
      -p $HOME/pyspark_my_way/scripts/kmos_server_hybw_IZ.cfg > hybw_IZ.out 2>&1 &

# 2. progress: rv = final cubes, wing = how many of them were made with the hybrid and which hybrid_wing
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ.cfg --status

# 3. finished? only once the job has ended (pgrep -f kmos_server_hybw_IZ.cfg prints nothing)
grep -c "Finished reducing" /local_storage/KOLJA/hyb_work/IZw/hybw_IZ.out      # expected 83
grep -c "An error occ" /local_storage/KOLJA/hyb_work/IZw/hybw_IZ.out           # expected 0
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybw_IZ.cfg --status
#    expected summary '83 OBs: 83 started, 83 with sky/ cubes, 83 with final cubes, 0 started with a wrong or missing
#    degree'; per OB rv = number of OBJECT frames, wing NxW3 for the mismatched lists

# 4. combination as for PROCDIR_HYB, new folders (EDPS: not while the YJ combination is running)
python scripts/make_edps_input.py IZ --procdir /local_storage/KOLJA/PROCDIR_HYBW --out /local_storage/KOLJA/edps_in_hybw_IZ
#    expected 'OB folders: 83'; fewer means an OB failed: do not combine
mkdir -p /home/kreuter/EDPS_data_hybw
sed -i.bak 's#^base_dir=.*#base_dir=/home/kreuter/EDPS_data_hybw#' ~/.edps/application.properties
#    restart edps-gui --port 5007, run cubes_combination on edps_in_hybw_IZ, export to /local_storage/KOLJA/res_hybw_IZ
python scripts/edps_pick_jobs.py --base /home/kreuter/EDPS_data_hybw/KMOS/cubes_combination --band IZ
#    one job per IZ target list, each '<- clean'; no jobs: edps-gui was not restarted after the base_dir change
```

Locally afterwards, in EMPOWER_Analysis_and_presentation, new folders: the download in WSL, the telluric step with the
Windows Python (the WSL python has no astropy). `--method best` is the default of the analysis copy of
telluric_combined.py and makes an old copy elsewhere fail instead of applying v1:
```bash
REMOTE_BASE=/local_storage/KOLJA/res_hybw_IZ LOCAL_BASE=kmos_comb_hybw_IZ bash scripts/download_combined.sh
python scripts/telluric_combined.py apply --comb kmos_comb_hybw_IZ --index telluric_frame_index_ohlin.json \
       --tell telluric_all --out kmos_comb_hybw_IZ_tell --method best
```

Expected: matched lists (no skysub_det cubes) take the unchanged classic path and should equal PROCDIR_HYB (local:
kmos_comb_hyb, kmos_comb_hyb_tell). In the
mismatched lists the change is largest in the OH-wing planes; ZAP (SVD projection, 101-plane continuum) moves the
other planes at the per-cent level (x0.96-1.03 locally). Locally the between-IFU correlation of the background in
the wing planes went 0.43 -> 0.03 (mid) and 0.36 -> 0.18 (red) for h21.4 and to <= 0.05 for legac, and the red
wing-plane noise of h21.4 to x0.8-0.9. The first check is that correlation test on the combined cubes (W = 3 was
chosen on the two local test OBs). Disk: the new science/ folders only (four stages per OBJECT frame, as in
PROCDIR_HYB). Time: about 4 min per OB on the server (the sky-subtraction-to-velocity part of the PROCDIR_OHLIN run,
docs/oh-lcorr-linear-rereduction.md), so 5-6 h for the 83 OBs; 2-2.5 min per OB locally.

Tested locally on a mock of this layout (WSL ~/impl_test/hybw: raw, PROCDIR and PROCDIR_HYB as links into the local
hybrid run of three IZ OBs, colon folder names): `--reuse` refuses a cfg with the reconstruction on; final _s cubes
bit-identical to the validated W = 3 frames (matched OB: to the classic frames); all four stages written; WING = 3
on the 12 hybrid frames; every reused entry still a link and every reused file unchanged (md5, mtime).

## YJ (PROCDIR_HYB, 2026-10-27)

About 80 % of the YJ object time is in OBs with mismatched SKY frames (the census's 84 % counts SKY frames in another
calibration-angle bin instead; 33 of the 45 YJ OBs: SPT2040 15 OBs, mostly 80-83 deg,
cosmos_0.6_0.8_h21.4 7 at 83-90 deg, XMM2215 9, mostly 36-39 deg with pairs up to 54 deg, cosmos_target_legac and
cosmos_target_1 1 each at 38 deg; the 11 cosmos OBs and the cdfs_targets_2 OB are matched; the hybrid decides per
frame from the rotator angles anyway, det 0 in --status = classic). The hybrid works for YJ with two YJ-specific model options
(`hybrid_free_model = ifumed`, `hybrid_oh_model = pos3`; hybrid-sky.md, 'YJ'), tested locally on one h21.4 and one
SPT2040 OB fetched with scripts/fetch_yj_test.sh (ssh + tar + md5; the server has no rsync). Not tested at 36-54 deg:
XMM2215 (the YJ BPT target) is the list to look at first. scripts/kmos_server_hyb_YJ.cfg is kmos_server_ohlin_YJ.cfg
(same 45 OBs, patched DRP, lcorr degree 1, standard YJ OH_SPEC, default ZAP basis) with procdir PROCDIR_HYB, the row
common mode of detectors 2 and 3 (dcr_ frames) and the hybrid. This is a full chain (reconstruction + detector-level
cubes), so it takes about as long per OB as the first IZ hybrid run; the YJ OB folders are new in PROCDIR_HYB (no
overlap with the IZ/H folders there).

Commit, move-aside and pull: section 'Pending runs' (once for both pending runs). With older server code the YJ options
are silently ignored (the IZ model is applied, which is harmful in YJ). The smoke test shows it: with the new code the
log lines end '(wing 3, oh pos3, free ifumed)' and --status has a 'wing' column reading '6xW3/pos3/ifumed'; if the
suffix is missing, or --status has no 'wing' column or shows only '6xW3', the server code is older: git pull, check
`git log -1` and run the smoke test again.

```bash
cd ~/pyspark_my_way && git log -1 --oneline     # must show the commit from 'Pending runs'

# 0. smoke test on the h21.4 YJ OB (same procdir as the full run), about 30-40 min; expectations in the cfg header
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ_test.cfg
mkdir -p /local_storage/KOLJA/hyb_work/YJtest && cd /local_storage/KOLJA/hyb_work/YJtest
env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
    -p $HOME/pyspark_my_way/scripts/kmos_server_hyb_YJ_test.cfg > hyb_YJtest.out 2>&1
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ_test.cfg --status
#    expected: dc 9, dcr 9, sky 9, det 6, rv 6, wing 6xW3/pos3/ifumed, degree 1
grep -h "Hybrid sky\|classic O - S" $(ls -t /local_storage/KOLJA/PROCDIR_HYB/KMOS.2026-02-11T06:31:26.207/*_YJ_*.log | head -1)
#    expected: 'Hybrid sky: 6 of 6 ...' and six lines ending '(wing 3, oh pos3, free ifumed)', no 'classic O - S'
#    (newest OB log only: a re-run of the smoke test adds a log and keeps the old one)

# 1. all 45 YJ OBs
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ.cfg --dry-run
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ.cfg
mkdir -p /local_storage/KOLJA/hyb_work/YJ && cd /local_storage/KOLJA/hyb_work/YJ
nohup env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
      -p $HOME/pyspark_my_way/scripts/kmos_server_hyb_YJ.cfg > hyb_YJ.out 2>&1 &
# progress: rv = final cubes, det = skysub_det cubes (0 for matched OBs), wing = hybrid frames and options
cd ~/pyspark_my_way && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ.cfg --status

# 2. finished? only once the job has ended (pgrep -f kmos_server_hyb_YJ.cfg prints nothing)
grep -c "Finished reducing" /local_storage/KOLJA/hyb_work/YJ/hyb_YJ.out      # expected 45
grep -c "An error occ" /local_storage/KOLJA/hyb_work/YJ/hyb_YJ.out           # expected 0
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hyb_YJ.cfg --status
#    expected summary '45 OBs: 45 started, 45 with sky/ cubes, 45 with final cubes, 0 started with a wrong or missing
#    degree'; per OB rv = number of OBJECT frames, wing NxW3/pos3/ifumed for the mismatched OBs, '-' for matched ones

# 3. combination as for IZ, own EDPS store (make_edps_input.py takes only the YJ final cubes of PROCDIR_HYB);
#    EDPS: not while the OH-wing combination is running (finish one combination, its export and edps_pick_jobs first)
python scripts/make_edps_input.py YJ --procdir /local_storage/KOLJA/PROCDIR_HYB --out /local_storage/KOLJA/edps_in_hyb_YJ
#    expected 'OB folders: 45'; fewer means an OB failed: do not combine
mkdir -p /home/kreuter/EDPS_data_hyb_YJ
sed -i.bak 's#^base_dir=.*#base_dir=/home/kreuter/EDPS_data_hyb_YJ#' ~/.edps/application.properties
#    restart edps-gui --port 5007, run cubes_combination on edps_in_hyb_YJ, export to /local_storage/KOLJA/res_hyb_YJ
python scripts/edps_pick_jobs.py --base /home/kreuter/EDPS_data_hyb_YJ/KMOS/cubes_combination --band YJ
#    one job per YJ target list, each '<- clean'; no jobs: edps-gui was not restarted after the base_dir change
```

Locally afterwards, in EMPOWER_Analysis_and_presentation: the download in WSL, the telluric step with the Windows
Python (the WSL python has no astropy):
```bash
REMOTE_BASE=/local_storage/KOLJA/res_hyb_YJ LOCAL_BASE=kmos_comb_hyb_YJ bash scripts/download_combined.sh
python scripts/telluric_combined.py apply --comb kmos_comb_hyb_YJ --index telluric_frame_index_ohlin.json \
       --tell telluric_all --out kmos_comb_hyb_YJ_tell --method best
```

What to check, against the old OHLIN YJ cubes (local metrics: scratchpad wing/yj_metrics.py style, edge-inclusive):
- cosmos_0.6_0.8_h21.4 and SPT2040 (as tested): interior noise (>= 3 spaxels from the cube edge) x0.2-0.65 in the
  free/wing planes and x0.5-0.9 in the OH planes, OH planes on the two outermost spaxels x0.5-1.1 (test: up to x1.08
  on the second spaxel, red, SPT2040), between-IFU covariance of the per-plane background within ~x2 of classic,
  classic streak/dipole residuals gone, line fluxes at fixed kinematics unchanged within the errors (part of the line
  S/N gain comes from the row common mode alone: forced line errors x0.75-0.89 on the h21.4 test OB);
- XMM2215, legac, cosmos_target_1 (36-54 deg, untested): the same checks first; if the OH planes or central
  apertures come out worse than classic there, re-run only the sky subtraction for those OBs with other options
  (PROCDIR_HYBW-style `--reuse`, a few minutes per OB);
- matched YJ OBs (cosmos, cdfs_targets_2): row common mode only, x0.99-1.00.
The options were chosen on two OBs.
