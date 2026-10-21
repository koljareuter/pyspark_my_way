# OH-line wavelength correction: linear fit and new IZ OH spectrum (re-reduction of all OBs, H stays quadratic)

**Date:** 2026-10-21
**Affected stage:** cube reconstruction (`kmos_sci_red` with `OH_SPEC`, called from `pyspark/cubes/spark_reconstruct.py`)
**Affected data:** every OB, worst in IZ below ~0.9 µm (the Hβ range of the z ≈ 0.6–0.8 targets)
**Analysis:** EMPOWER_Analysis_and_presentation, notebooks `kmos_sky_residuals`, `kmos_iz_wavelength_fix`,
`kmos_wavelength_census`, `kmos_iz_new_ohspec_test`

---

## Symptom

Towards the blue end of IZ, every spectrum shows a "forest" of OH residuals. The census of all 149 OBs found
the object and sky frames misaligned by 7–29 km/s (median 13 km/s) at 0.79–0.87 µm in all 82 IZ OBs, and an
absolute offset scatter between IFUs of ~60 km/s. YJ and H are fine (≤ 1.3 and 0.5 km/s).

## Cause

`kmos_sci_red` corrects the arc wavelength solution per exposure and IFU with the OH lines of the
reference spectrum `kmos_oh_spec_<band>.fits` (`kmos/kmo_priv_lcorr.c`, `kmo_lcorr_crosscorrelate_spectra`):
lines brighter than 20 % of the brightest one, Gaussian centroids, outlier clipping, then an **unweighted
2nd-order polynomial** Δλ(λ) (`cpl_size degree[1] = {2};`, hard-coded).

The standard IZ reference spectrum gives only 10 usable lines, all at 0.979–1.090 µm. The quadratic is
extrapolated down to 0.78 µm: a centroid error is amplified ×12 at 0.83 µm and ×16 at 0.80 µm. Object and
sky frames get independent, wrong curvatures, so they no longer line up in the blue.

## Fix

1. **New IZ OH spectrum** `/local_storage/EMPOWER_raw_data/input_all/static/kmos_oh_spec_iz_NEW.fits`:
   27 usable lines, the bluest at 0.835 µm (amplification ×0.8 instead of ×12 at 0.83 µm).
2. **Linear instead of quadratic fit**: a patched DRP reads the degree from the environment variable
   `KMOS_LCORR_DEGREE` (unset = 2 = original behaviour).

### DRP patch

`scripts/drp_patch/kmos_lcorr_degree_patch.py` edits `kmos/kmo_priv_lcorr.c` of KMOS 4.5.11 (adds
`#include <stdlib.h>` and reads `KMOS_LCORR_DEGREE`, allowed 0–4, after `cpl_size degree[1] = {2};`). The
recipe log reports the degree once per run:

    [ INFO  ] kmos_sci_red: OH lambda correction: polynomial degree 1 (KMOS_LCORR_DEGREE 1)

`scripts/drp_patch/build_kmos_lcorr_drp.sh build` unpacks the same source tarball as the existing pipeline,
applies the patch and installs into `/home/kreuter/KMOS_PIPELINE_4.5.11_lcorr/install`, linked against the
CPL and telluriccorr of the existing install. The existing pipeline is not touched.

### pyspark

New optional `[cubes]` settings; without them nothing changes:

| setting | effect |
|---|---|
| `oh_spec_<band>` (e.g. `oh_spec_iz`) | OH_SPEC file instead of `$KMOS_CALIB/kmos_oh_spec_<band>.fits` |
| `recipe_dir` | passed to esorex as `--recipe-dir` (the patched DRP) |
| `lcorr_degree` | sets `KMOS_LCORR_DEGREE` for the `kmos_sci_red` call; refused without `recipe_dir` |

Code: `pyspark/cubes/spark_reconstruct.py::reconstruct`, `pyspark/utils/obsutils.py::run_spark` (new
arguments `esorex_opts`, `env`).

## Verification (local, WSL, KMOS 4.5.11)

All versions reconstructed from the same dc_ frames and calibrations with the same `kmos_sci_red` options as
pyspark; measured with the analysis tools (OH centroids vs the UVES/Oliva line lists, O − S of the pairs pyspark
subtracts, frame jitter). Tables: `EMPOWER_Analysis_and_presentation/kmos_drp_linear_verify/`, script
`scripts/kmos_drp_linear_verify.py` there.

**Patched DRP without the variable = standard DRP**, bit by bit (data of all extensions): 9 of 9 frames of the
worst IZ OB (new OH_SPEC), 9 of 9 of the H OB, and 9 of 9 against the production sky/ cubes of the worst IZ OB
(`build_kmos_lcorr_drp.sh smoke`, run locally with the server copies). pyspark `reconstruct()` with `oh_spec_iz`, `recipe_dir`,
`lcorr_degree = 1` gives sky/ cubes bit-identical to the direct DRP run (9 of 9), and the OB log carries the
degree line.

**IZ, worst OB `KMOS.2026-02-12T06:02:50.540`** (Δrot 93° between object and sky frames), median in km/s:

| version | O − S, Hβ range 0.79–0.87 | O − S, red 0.95–1.07 | IFU scatter of the absolute offset, Hβ / red | absolute offset, Hβ | jitter, Hβ |
|---|---|---|---|---|---|
| uncorrected (no OH_SPEC) | 17.0 | 17.6 | 8.3 / 1.5 | +3.8 | 10.0 |
| quadratic, old OH_SPEC (production) | 33.6 | 1.4 | 54.0 / 1.0 | +9.0 | 23.5 |
| quadratic, new OH_SPEC | 7.1 | 1.0 | 4.3 / 0.5 | −8.3 | 5.5 |
| **linear, new OH_SPEC (patched DRP)** | **5.0** | **0.9** | **1.8 / 0.3** | **−1.5** | **3.8** |
| replica prediction, linear + new | 4.7 | 0.8 | 2.4 / 0.5 | −1.2 | 3.5 |

The linear fit also removes the common −8 km/s bend that the quadratic leaves at Hβ, and the red end does not
get worse. The DRP agrees with the Python replica's prediction. The 90th percentile of O − S at Hβ stays at
≈ 30 km/s in every version, including uncorrected; those are 0.01 µm chunks with too few or too faint OH lines to
measure, and the residual they leave is small (90th percentile 0.34 of the box noise, production 2.1).

**YJ, `KMOS.2026-01-13T04:40:23.303`** (standard OH_SPEC), median in km/s:

| version | O − S, blue 1.00–1.08 / red 1.10–1.33 | IFU scatter, blue / red | jitter, blue / red |
|---|---|---|---|
| uncorrected | 0.30 / 0.20 | 8.1 / 6.6 | 0.39 / 0.29 |
| quadratic (production) | 0.42 / 0.20 | 1.0 / 0.3 | 0.50 / 0.24 |
| **linear (patched DRP)** | **0.38 / 0.19** | **0.9 / 0.2** | **0.44 / 0.21** |

YJ was already fine; linear is as good as quadratic (all differences < 0.2 km/s). The OH correction is still
needed for the per-IFU arc zero points (8 km/s scatter without it).

**H, `KMOS.2026-06-07T07:24:38.627`** (standard OH_SPEC, Δrot 89°), median in km/s; residual = O − S residual
in units of the box noise:

| version | O − S, blue 1.45–1.52 / red 1.55–1.78 | IFU scatter, blue / red | residual, blue / red |
|---|---|---|---|
| uncorrected | 9.7 / 8.9 | 2.0 / 1.4 | 9.1 / 15.2 |
| **quadratic (production)** | **0.58 / 0.21** | **0.41 / 0.17** | **0.53 / 0.38** |
| linear (patched DRP) | 0.66 / 0.29 | 0.84 / 0.16 | 0.59 / 0.54 |

In H the linear fit is slightly **worse** than the quadratic: with 38 lines spread over the band the quadratic
term is constrained and apparently real (blue IFU scatter 0.4 → 0.8 km/s, red residual 0.38 → 0.54 of the box
noise, which matters with H's bright lines). The replica predicts the same (0.70 in red).

**Conclusion:** IZ needs the new OH_SPEC and gains clearly from the linear fit; YJ is indifferent; H is better
off with the DRP's quadratic. A degree per band is one line in each config (`lcorr_degree`); without it a band
runs exactly as the standard DRP.

**Choice for the re-reduction (2026-10-21):** IZ linear + new OH_SPEC, YJ linear, H quadratic (`lcorr_degree = 2`,
so the H sky/ cubes equal the current ones; H is re-run only to have all bands in one folder tree).

## Re-reduction on hqempower1

Everything goes to new folders; `PROCDIR`, `kmos_ob_sorted`, `edps_in_<BAND>` and the existing results are only read.

| what | where |
|---|---|
| patched DRP | `/home/kreuter/KMOS_PIPELINE_4.5.11_lcorr/` |
| new PROCDIR (links to dc_ frames + calib/, new sky/ science/) | `/local_storage/KOLJA/PROCDIR_OHLIN/<OB>/` |
| pyspark working directories (one per band) | `/local_storage/KOLJA/ohlin_work/<BAND>/` |
| EDPS input links | `/local_storage/KOLJA/edps_in_ohlin_<BAND>/` |
| combined cubes | `/local_storage/KOLJA/res_ohlin_<BAND>/` |
| census of the new sky/ cubes | `/local_storage/KOLJA/ohspec_census/ohlin/` |

Configs: `scripts/kmos_server_ohlin_{IZ,YJ,H}.cfg` (83 + 45 + 18 OBs; the three single-exposure OBs without a
sky frame are left out, as they have no science cubes in the current reduction either). Calibrations and
detector corrections are not redone: the configs start at the reconstruction, from the dc_ frames and calib/
of PROCDIR (current state, including the det3 IZ flat fix of 2026-10-15).

### Steps (on hqempower1, from `~/pyspark_my_way`, after `git pull`)

```bash
# 1. patched DRP (writes only /home/kreuter/KMOS_PIPELINE_4.5.11_lcorr)
bash scripts/drp_patch/build_kmos_lcorr_drp.sh build
#    expect: "patch present in .../libkmos.so.0", ldd shows libkmos from KMOS_PIPELINE_4.5.11_lcorr
bash scripts/drp_patch/build_kmos_lcorr_drp.sh smoke
#    expect: default 9 cubes "degree 2 (KMOS_LCORR_DEGREE unset)", lin 9 cubes "degree 1",
#            "9 of 9 frames bit-identical" (patched DRP without the variable = production sky/ cubes)

# 2. new PROCDIR with links (writes only PROCDIR_OHLIN)
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg --dry-run    # 146 OBs, no PROBLEM lines
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg

# 3. reduction, one working directory per band (kmos_sci_red writes into the cwd)
for B in IZ YJ H; do
  mkdir -p /local_storage/KOLJA/ohlin_work/$B && cd /local_storage/KOLJA/ohlin_work/$B &&
  nohup env PYTHONPATH=$HOME/pyspark_my_way python $HOME/pyspark_my_way/scripts/reduce_kmos.py \
        -p $HOME/pyspark_my_way/scripts/kmos_server_ohlin_$B.cfg > ohlin_$B.out 2>&1 &
done
cd ~/pyspark_my_way

# 4. progress and check that every OB ran with degree 1 and the right OH_SPEC
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_ohlin_*.cfg --status
#    IZ: OH_SPEC kmos_oh_spec_iz_NEW.fits, YJ/H: kmos_oh_spec_yj/h.fits; "degree in log" = 1 for IZ and YJ, 2 for H

# 5. wavelength census of the new sky/ cubes (EMPOWER_Analysis_and_presentation tools, new output folder)
OMP_NUM_THREADS=1 nice -n 10 /local_storage/KOLJA/ohspec_census/venv/bin/python \
    /local_storage/KOLJA/ohspec_census/tools/kmos_census_measure.py \
    --sky-root /local_storage/KOLJA/PROCDIR_OHLIN --out /local_storage/KOLJA/ohspec_census/ohlin --workers 8

# 6. combination: NEW EDPS input folders (--out is mandatory, the default edps_in_<BAND> holds the current links)
for B in IZ YJ H; do
  python scripts/make_edps_input.py $B --procdir /local_storage/KOLJA/PROCDIR_OHLIN --out /local_storage/KOLJA/edps_in_ohlin_$B
done
#    run EDPS cubes_combination on edps_in_ohlin_<BAND> as before, then collect the clean jobs of this run:
python scripts/edps_pick_jobs.py --band IZ --since <date of step 6> --copy /local_storage/KOLJA/res_ohlin_IZ   # same for YJ, H
#    (do not run other combinations meanwhile, or check in the report that the inputs are PROCDIR_OHLIN cubes)
```

Locally (Windows/WSL, `pyspark_my_way`): download into new folders and correct the tellurics into a new mirror
folder; the current combined cubes (COMB_BASE) stay as they are.

```bash
REMOTE_BASE=/local_storage/KOLJA/res_ohlin_IZ LOCAL_BASE=kmos_comb_ohlin_IZ bash scripts/download_combined.sh   # YJ, H alike
python scripts/telluric_combined.py apply --comb kmos_comb_ohlin_IZ --index telluric_frame_index.json \
       --tell telluric_all --out kmos_comb_ohlin_IZ_tell
```

The frame index can be reused: the single frames keep their names and barycentric velocities.

### Disk and time

Links only for the inputs; new per 9-frame OB as the current folders: sky/ 0.29 GB + science/ 0.78 GB, so
≈ 160 GB for 146 OBs (/local_storage had 5.4 TB free on 2026-10-21). From the det3-fix log of the worst IZ OB:
≈ 6 min reconstruction + ≈ 4 min sky subtraction to velocity correction per OB, run one OB after the other per band,
so ≈ 14 h for IZ (83 OBs), ≈ 8 h YJ, ≈ 3 h H, with the three bands in parallel.
