#!/bin/bash
# KMOS DRP 4.5.11 with the KMOS_LCORR_DEGREE patch (degree of the OH-line wavelength correction in kmos_sci_red, see
# docs/oh-lcorr-linear-rereduction.md), built next to the existing pipeline. The existing pipeline, PROCDIR and the raw
# data are only read; everything is written under $DEST.
#
# On hqempower1, from ~/pyspark_my_way:
#   bash scripts/drp_patch/build_kmos_lcorr_drp.sh build          # unpack, patch, configure, make, install (~10 min)
#   bash scripts/drp_patch/build_kmos_lcorr_drp.sh smoke [OB]     # patched DRP on one OB (default the worst IZ OB):
#        default = degree unset, standard OH_SPEC  → must equal PROCDIR/<OB>/sky/ bit by bit
#        lin     = KMOS_LCORR_DEGREE=1, new IZ OH_SPEC (standard file for YJ/H) → log must report degree 1
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
PIPE=${PIPE:-/home/kreuter/KMOS_PIPELINE_4.5.11/install}      # existing pipeline: CPL, telluriccorr, esorex, static calibrations
TARBALL=${TARBALL:-/home/kreuter/KMOS_PIPELINE_4.5.11/pipeline-kit-2186009.tmp/tmp/kmos-4.5.11.tar.gz}   # md5 f8d1e09b7516834bb7edbbc69e3ec9e7
DEST=${DEST:-/home/kreuter/KMOS_PIPELINE_4.5.11_lcorr}
PROCDIR=${PROCDIR:-/local_storage/KOLJA/PROCDIR}
OHNEW=${OHNEW:-/local_storage/EMPOWER_raw_data/input_all/static/kmos_oh_spec_iz_NEW.fits}
PY=${PY:-python3}
RECIPES=$DEST/install/lib/esopipes-plugins
STATIC=${STATIC:-$PIPE/calib/kmos-4.5.11}

case ${1:-} in
  build)
    [ -e "$DEST/install" ] && { echo "$DEST/install exists already (remove it or set DEST)"; exit 1; }
    md5sum "$TARBALL"
    mkdir -p "$DEST/src" && tar xzf "$TARBALL" -C "$DEST/src" || exit 1
    $PY "$HERE/kmos_lcorr_degree_patch.py" "$DEST/src/kmos-4.5.11" || exit 1
    cd "$DEST/src/kmos-4.5.11" || exit 1
    { ./configure --prefix="$DEST/install" --with-cpl="$PIPE" --with-telluriccorr="$PIPE" &&
      nice -n 10 make -j8 && make install; } > "$DEST/build.log" 2>&1
    echo "build exit $? (log $DEST/build.log)"
    ls "$RECIPES/kmos-4.5.11/kmos_sci_red.so" || exit 1
    lib=$(find "$DEST/install" -name libkmos.so.0 | head -1)
    ldd "$RECIPES/kmos-4.5.11/kmos_sci_red.so" | grep -E "libkmos|libcplcore"   # libkmos must come from $DEST
    strings "$lib" | grep -q KMOS_LCORR_DEGREE && echo "patch present in $lib" || { echo "PATCH MISSING in $lib"; exit 1; }
    echo "recipe_dir for the pyspark configs: $RECIPES"
    ;;
  smoke)
    OB=${2:-KMOS.2026-02-12T06:02:50.540}; src=$PROCDIR/$OB; out=$DEST/smoke/$(echo "$OB" | tr ':' '-')
    band=$($PY -c "import sys; from astropy.io import fits; print(fits.getheader(sys.argv[1])['HIERARCH ESO INS GRAT1 ID'].strip())" \
           "$(ls "$src"/dc_KMOS.*.fits | head -1)")
    B3=$band$band$band; b=$(echo "$band" | tr 'A-Z' 'a-z')
    oh_lin=$STATIC/kmos_oh_spec_$b.fits; [ "$band" = IZ ] && oh_lin=$OHNEW
    echo "$OB  band $band  → $out"
    for v in default lin; do
      d=$out/$v; mkdir -p "$d"
      { for f in "$src"/dc_KMOS.*.fits; do echo "$f SCIENCE"; done
        for c in XCAL YCAL LCAL MASTER_FLAT; do echo "$src/calib/${c}_$B3.fits $c"; done
        echo "$STATIC/kmos_wave_band.fits WAVE_BAND"
        if [ $v = default ]; then echo "$STATIC/kmos_oh_spec_$b.fits OH_SPEC"; else echo "$oh_lin OH_SPEC"; fi
      } > "$d/in.sof"
      (cd "$d" || exit 1
       if [ $v = lin ]; then export KMOS_LCORR_DEGREE=1; else unset KMOS_LCORR_DEGREE; fi
       nice -n 10 "$PIPE/bin/esorex" --recipe-dir="$RECIPES" --output-dir="$d" --log-dir="$d" kmos_sci_red \
         -background=FALSE -tbsub=TRUE -no_subtract=TRUE -sky_tweak=FALSE -lcmethod=SLICES_MEDIAN -skip_sky_oh_align=TRUE \
         -no_combine=TRUE "$d/in.sof" > "$d/run.log" 2>&1)
      echo "$v: exit $?, $(ls "$d"/SCI_RECONSTRUCTED_*.fits 2>/dev/null | wc -l) cubes; $(grep -h 'polynomial degree' "$d/run.log" | head -1)"
    done
    $PY - "$out/default" "$src/sky" <<'EOF'
import glob, os, sys
import numpy as np
from astropy.io import fits
a, b = sys.argv[1:3]
n = same = 0
for p in sorted(glob.glob(os.path.join(a, "SCI_RECONSTRUCTED_*.fits"))):
    q = os.path.join(b, os.path.basename(p))
    if not os.path.exists(q):
        print("not in %s: %s" % (b, os.path.basename(p)))
        continue
    n += 1
    with fits.open(p) as ha, fits.open(q) as hb:
        same += len(ha) == len(hb) and all(
            (ea.data is None and eb.data is None) or
            (ea.data is not None and eb.data is not None and ea.data.shape == eb.data.shape
             and np.array_equal(ea.data, eb.data, equal_nan=True)) for ea, eb in zip(ha, hb))
print("patched DRP, degree unset vs PROCDIR sky/: %d of %d frames bit-identical" % (same, n))
EOF
    ;;
  *) echo "usage: $0 build | smoke [OB]"; exit 1 ;;
esac
