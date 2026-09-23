"""Compare _rv.fits white-light maps: new test run vs old reference.

Usage:
    python scripts/compare_rv.py

Saves compare_rv.png in the current directory.
"""
import glob, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

OLD_DIR = "/home/kreuter/kmos_ob_proc/KMOS.2026-01-13T02:54:58.606/science"
NEW_DIR = "/local_storage/KOLJA/PROCDIR_TEST/KMOS.2026-01-13T02:54:58.606/science"
OB_TAG  = "KMOS.2026-01-13T02:54:58.606"


def find_rv(directory, suffix="_rv.fits"):
    matches = sorted(glob.glob(f"{directory}/*{suffix}"))
    return matches


def whitelight(cube):
    """Collapse a (nwave, nx, ny) cube to a 2D image, ignoring NaNs."""
    return np.nanmedian(cube, axis=0)


def plot_ifu(ax, img, title, vmin=None, vmax=None):
    if vmin is None:
        vmin = np.nanpercentile(img, 5)
    if vmax is None:
        vmax = np.nanpercentile(img, 95)
    ax.imshow(img, origin="lower", vmin=vmin, vmax=vmax, cmap="viridis")
    ax.set_title(title, fontsize=7)
    ax.axis("off")


old_files = find_rv(OLD_DIR, "_rv.fits")
new_files = find_rv(NEW_DIR, "_rv.fits")

if not old_files:
    sys.exit(f"No _rv.fits found in OLD dir: {OLD_DIR}")
if not new_files:
    sys.exit(f"No _rv.fits found in NEW dir: {NEW_DIR}")

print(f"Old: {len(old_files)} files   New: {len(new_files)} files")

# Match by filename stem (dc_KMOS.<timestamp>.fits part)
old_map = {f.split("/")[-1].replace("SCI_RECONSTRUCTED_", "").replace("_s_i_bkg2d_rv.fits", ""): f
           for f in old_files}
new_map = {f.split("/")[-1].replace("SCI_RECONSTRUCTED_", "").replace("_s_i_bkg2d_rv.fits", ""): f
           for f in new_files}

common = sorted(set(old_map) & set(new_map))
print(f"Matched frames: {len(common)}")
if not common:
    print("Old stems:", list(old_map)[:3])
    print("New stems:", list(new_map)[:3])
    sys.exit("No matching frames — check filename format")

# Pick up to 3 frames and up to 6 IFUs each
n_frames = min(3, len(common))
n_ifus   = 6

fig, axes = plt.subplots(n_frames * 3, n_ifus,
                          figsize=(n_ifus * 1.8, n_frames * 3 * 1.8),
                          constrained_layout=True)
if n_frames == 1:
    axes = axes[np.newaxis]

fig.suptitle(f"_rv.fits comparison: OLD (kmos_ob_proc) vs NEW (PROCDIR_TEST)\n{OB_TAG}",
             fontsize=9, fontweight="bold")

for fi, stem in enumerate(common[:n_frames]):
    old_f, new_f = old_map[stem], new_map[stem]
    with fits.open(old_f) as ho, fits.open(new_f) as hn:
        ifu_count = 0
        col = 0
        for ifu in range(1, 25):
            if col >= n_ifus:
                break
            try:
                od = np.array(ho[ifu].data, dtype=float)
                nd = np.array(hn[ifu].data, dtype=float)
            except (IndexError, TypeError):
                continue
            if od.ndim != 3 or nd.ndim != 3:
                continue
            if np.all(np.isnan(od)) or np.all(np.isnan(nd)):
                continue

            ow = whitelight(od)
            nw = whitelight(nd)
            diff = nw - ow

            vmin = np.nanpercentile(np.concatenate([ow.ravel(), nw.ravel()]), 5)
            vmax = np.nanpercentile(np.concatenate([ow.ravel(), nw.ravel()]), 95)
            dmax = np.nanpercentile(np.abs(diff.ravel()), 99)

            plot_ifu(axes[fi * 3 + 0, col], ow, f"OLD IFU{ifu}\n{stem[-12:]}", vmin, vmax)
            plot_ifu(axes[fi * 3 + 1, col], nw, f"NEW IFU{ifu}", vmin, vmax)
            plot_ifu(axes[fi * 3 + 2, col], diff, f"NEW-OLD", -dmax, dmax)
            axes[fi * 3 + 2, col].images[0].set_cmap("RdBu_r")

            col += 1
            ifu_count += 1

    row_labels = [f"Frame {fi+1}  OLD", f"Frame {fi+1}  NEW", f"Frame {fi+1}  DIFF (new-old)"]
    for ri, lbl in enumerate(row_labels):
        axes[fi * 3 + ri, 0].set_ylabel(lbl, fontsize=7, rotation=90, labelpad=2)

outfile = "compare_rv.png"
plt.savefig(outfile, dpi=120, bbox_inches="tight")
print(f"Saved → {outfile}")
