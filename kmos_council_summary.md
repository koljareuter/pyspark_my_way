# KMOS programme 116.29EG — numbers for the council

## 1. Telescope time used (ESO archive, as of 2026-10-01)

- **149 h** of telescope time (first acquisition to end of last exposure) in **80 OB executions** on **36 nights**
- 127 h open shutter, about 85 h of it on source (O-S-O nod-to-sky: 2/3 of frames are on target)
- By band: H 26 h, IZ 73 h, YJ 51 h

## 2. Stacked data set (kmos_proc_server_comb)

| Field | z | Band | Galaxies | On-source / galaxy | Line | S/N ≥ 3 | S/N ≥ 5 | Best S/N |
|---|---|---|---|---|---|---|---|---|
| COSMOS | 0.70 | YJ | 19 | 2.9 h | Hα | 8 | 6 | 19.4 |
| CDFS | 0.61 | IZ | 4 | 1.0 h | Hα | 1 | 1 | 17.7 |
| CDFS | 0.68 | IZ | 16 | 1.0 h | [OIII]5007 | 1 | 0 | 5.0 |
| SPT2040 | 1.48 | YJ | 21 | 11.2 h | [OIII]5007 | 3 | 1 | 6.1 |
| XMM2215 | 1.46 | YJ | 18 | 6.7 h | [OIII]5007 | 6 | 2 | 6.8 |

- **80 galaxies** in the stacked cubes; Hα falls in the stacked band for **23**, of which **7 have Hα S/N ≥ 5** (9 with S/N ≥ 3).
- Another **3** galaxies have [OIII]/Hβ at S/N ≥ 5 (z≈1.47 clusters and CDFS).
- Brightest Hα: COSMOS 1147707 (S/N 19), CDFS 23420 (S/N 18), COSMOS 1168518 (S/N 14), COSMOS 1123917 (S/N 8), COSMOS 1137103 (S/N 6), COSMOS 1195224 (S/N 5)

## 3. Extra data on the same galaxies

| Field | Galaxies | Line | Stacked now | Own, not yet stacked | Archival (other programmes) | Depth gain √(t_all/t_now) | Expected S/N ≥ 5 (now → all) |
|---|---|---|---|---|---|---|---|
| COSMOS | 19 | Hα | 2.9 h | 0.0 h (median) | median 0.0 h, max 0.0 h (0/19 galaxies) | ×1.0 | 6 → 6 |
| CDFS | 4 | Hα | 1.0 h | 10.6 h (median) | median 3.2 h, max 6.0 h (3/4 galaxies) | ×3.8 | 1 → 1 |
| CDFS | 16 | [OIII]5007 | 1.0 h | 9.7 h (median) | median 0.0 h, max 0.0 h (0/16 galaxies) | ×3.3 | 0 → 4 |
| SPT2040 | 21 | [OIII]5007 | 11.2 h | 0.9 h (median) | median 0.0 h, max 3.4 h (2/21 galaxies) | ×1.0 | 1 → 1 |
| XMM2215 | 18 | [OIII]5007 | 6.7 h | 0.8 h (median) | median 6.5 h, max 27.9 h (13/18 galaxies) | ×1.5 | 2 → 7 |

**Hα outside the stacked band.** For the 55 galaxies where Hα is not in the stacked cube (CDFS at z>0.63 in IZ; clusters at z≈1.47, where Hα is in H):
- CDFS: 16 galaxies; own Hα-band data median 0.8 h; archival Hα-band data median 2.9 h (max 6.0 h, 15 galaxies with some)
- SPT2040: 21 galaxies; own Hα-band data median 6.7 h; archival Hα-band data median 0.0 h (max 9.1 h, 8 galaxies with some)
- XMM2215: 18 galaxies; own Hα-band data median 8.5 h; archival Hα-band data median 8.9 h (max 14.9 h, 15 galaxies with some)

## Method / caveats

- S/N = integrated line flux in a 1.1″-diameter aperture centred on the line peak (±350 km/s, re-centred within ±400 km/s of the catalogue z), divided by the scatter of same-width windows across the line-free spectrum. Sky-line residuals are therefore included in the noise.
- Projection: S/N ∝ √t (background-limited, similar conditions). This assumes the archival data are of similar quality and combine optimally, so treat it as an upper estimate. For non-detections only the depth gain is meaningful.
- Archival time = on-source exposure from other programmes with an IFU within 1″ of the target (COSMOS/ECDFS: local ESO IFU tables; SPT2040/XMM2215: ESO TAP + raw-header arm positions, nod-to-sky ×2/3).
