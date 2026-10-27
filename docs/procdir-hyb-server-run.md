# Server re-reduction with the October 2026 noise fixes (PROCDIR_HYB)

**What changes against PROCDIR_OHLIN** (same OBs, same patched DRP, same OH_SPEC and lcorr degrees):

| band | new | expected gain | doc |
|---|---|---|---|
| IZ | hybrid sky subtraction for OBJECT frames whose SKY frame is > 20 deg away in rotator angle | 3x3-aperture noise of the combined cubes x0.39-0.48 (OH-free Hbeta window), x0.61-0.69 (red) for the mismatched lists; matched OBs untouched | hybrid-sky.md |
| IZ | row common-mode correction, detectors 2 and 3 (dcr_ frames) | noise x0.92-0.95 for IFUs 9-24 | noise-improvements-2026-10.md |
| H | SKY frames with the OBJECT frames' calibration angle (`sky_cal_angle = object`) | +5-11 % per frame | noise-improvements-2026-10.md |

`snap_wcs` (integer-shift combination) is not used: after the hybrid it raised the aperture noise by 5-14 %
in the end-to-end test. YJ stays as in PROCDIR_OHLIN until the hybrid has been tested on a mismatched YJ OB
(see the end).
PROCDIR, PROCDIR_OHLIN and kmos_ob_sorted are only read; everything new goes to PROCDIR_HYB(TEST).

## Steps on hqempower1

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

## Locally afterwards (WSL, in EMPOWER_Analysis_and_presentation, new folders only)

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

## What to check in the results

* Mismatched IZ lists (cosmos_0.6_0.8_h21.4, goods_0.6_0.8_h21.6_lir, cosmos_target_legac, cosmos_0.6_0.8_h21.7):
  noise in the Hbeta window about x0.4-0.6 of PROCDIR_OHLIN at equal frames; spurious "lines" made of
  positive/negative stripes in the classic cubes disappear (locally two [OIII] "detections" at S/N 5.7 and
  10.7 of cosmos_0.6_0.8_h21.4 were such artefacts; the line maps show dipoles and streaks, no galaxy).
* cdfs_targets_2 / cosmos IZ: only the row common mode, a few per cent for IFUs 9-24.
* Real lines keep their flux (locally cdfs 23420: Halpha flux ratio 1.00, [OIII] within 1 sigma).
* The single frames of the smoke-test OB should equal the local ones (WSL ~/impl_test/e2e/procdir_C) up to
  the sky subtraction if the dc frames and calibrations are the same.

## YJ (pending)

84 % of the YJ object time has mismatched SKY frames, but there is no mismatched YJ OB locally. To test the
hybrid on YJ, copy two mismatched YJ OBs (detector-corrected frames + calibrations, about 1 GB) into WSL:

```bash
mkdir -p ~/yj_test && cd ~/yj_test
for OB in KMOS.2026-02-11T06:31:26.207 KMOS.2026-06-08T05:24:25.056; do
  D=$(echo $OB | tr : -); mkdir -p $D/in $D/calib
  rsync -avL "kreuter@hqempower1:/local_storage/KOLJA/PROCDIR/$OB/dc_KMOS.*.fits" $D/in/
  rsync -avL "kreuter@hqempower1:/local_storage/KOLJA/PROCDIR/$OB/calib/{XCAL,YCAL,LCAL,MASTER_FLAT}_YJYJYJ.fits" $D/calib/
done
```
