#!/usr/bin/env python3
"""
Run the full zapv2+telluric pipeline on every OB in kmos_ob_sorted.

Usage (on hqempower1, inside screen/tmux):
    cd ~/pyspark_my_way
    git pull
    PYTHONPATH=. python scripts/run_server_full.py [--dry-run] [--proc-dir NAME]

The script auto-discovers all OB subdirs in RAW_DIR, writes a temporary
config, and calls reduce_kmos.py.  Results land in PROC_DIR (default:
~/kmos_ob_proc_telluric) — a separate directory so existing kmos_ob_proc
is untouched.

Options:
  --dry-run       Print config and OB list, do not run the pipeline.
  --proc-dir DIR  Override the output directory name (relative to HOME).
"""

import argparse
import os
import subprocess
import sys
import tempfile
from datetime import datetime

# ── paths ─────────────────────────────────────────────────────────────────────
HOME     = os.path.expanduser("~")
RAW_DIR  = os.path.join(HOME, "kmos_ob_sorted")
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CONFIG_TEMPLATE = """\
[instrument]
obsband = IZ,YJ,H

[pathing]
rawdir = {rawdir}
procdir = {procdir}
{subdirs_block}

[detector]
channel_correct = True
picture_correct = True
interchannel_correct = False
prefix = dc

[calibration]
make_dark = True
make_flat = True
make_arc = True
make_illum = True
reconstruct_arc = True
make_standard = True
darkexp = 300

[cubes]
reconstruct_cubes = True
sky_method = zapv2
rvcorr = True
use_oh = True
illum_corr = True
bkg_corr = True
telluric_in_sof = True
get_star_params = False
star_prefix = None
"""


def make_subdirs_block(subdirs):
    lines = [f"subdirs = {subdirs[0]}"]
    for s in subdirs[1:]:
        lines.append(f"          {s}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print config and exit without running")
    parser.add_argument("--proc-dir", default="kmos_ob_proc_telluric",
                        help="Output directory name under HOME (default: kmos_ob_proc_telluric)")
    args = parser.parse_args()

    proc_dir = os.path.join(HOME, args.proc_dir)

    # ── discover OBs ──────────────────────────────────────────────────────────
    if not os.path.isdir(RAW_DIR):
        sys.exit(f"ERROR: rawdir not found: {RAW_DIR}")

    subdirs = sorted(
        d for d in os.listdir(RAW_DIR)
        if os.path.isdir(os.path.join(RAW_DIR, d)) and d.startswith("KMOS.")
    )
    if not subdirs:
        sys.exit(f"ERROR: no KMOS.* subdirs found in {RAW_DIR}")

    print(f"Found {len(subdirs)} OB subdirs in {RAW_DIR}")
    print(f"Output → {proc_dir}")

    # ── write config ──────────────────────────────────────────────────────────
    cfg_text = CONFIG_TEMPLATE.format(
        rawdir=RAW_DIR,
        procdir=proc_dir,
        subdirs_block=make_subdirs_block(subdirs),
    )

    if args.dry_run:
        print("\n── generated config ──────────────────────────────────────────\n")
        print(cfg_text)
        print("── dry-run: pipeline NOT started ────────────────────────────")
        return

    cfg_file = tempfile.NamedTemporaryFile(
        mode="w", suffix=".cfg", prefix="kmos_server_",
        delete=False, dir=HOME
    )
    cfg_file.write(cfg_text)
    cfg_file.flush()
    cfg_name = cfg_file.name
    cfg_file.close()
    print(f"Config written to {cfg_name}")

    # ── run ───────────────────────────────────────────────────────────────────
    env = dict(os.environ, PYTHONPATH=REPO_DIR)
    reduce_script = os.path.join(REPO_DIR, "scripts", "reduce_kmos.py")
    cmd = [sys.executable, reduce_script, "-p", cfg_name]

    print(f"Starting: {' '.join(cmd)}")
    print(f"Started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    try:
        subprocess.run(cmd, env=env, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Pipeline exited with code {e.returncode}")
    finally:
        print(f"\nFinished at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
