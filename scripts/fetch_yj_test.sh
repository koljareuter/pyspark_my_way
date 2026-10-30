#!/usr/bin/env bash
# Copy the two mismatched YJ test OBs (dc_ frames + YJ calibrations) from the server PROCDIR into
#   ~/yj_test/<OB, ':' -> '-'>/{in,calib}
# for the local hybrid sky test (chain_band.py YJ classic|row|hyb). Needs only ssh and tar on the server (no
# rsync): one tar stream per OB, then an md5 comparison of every file. About 2 GB per OB.
#
# Run from WSL:
#   bash "/mnt/c/Users/kreuter/OneDrive - ESO/Dokumente/ESO_rep/pyspark_my_way/scripts/fetch_yj_test.sh"
# Other OBs / target:  OBS="KMOS.2026-...:..." DST=$HOME/other bash scripts/fetch_yj_test.sh
set -euo pipefail

SRV=${SRV:-kreuter@hqempower1}
SRC=${SRC:-/local_storage/KOLJA/PROCDIR}
DST=${DST:-$HOME/yj_test}
OBS=${OBS:-"KMOS.2026-02-11T06:31:26.207 KMOS.2026-06-08T05:24:25.056"}
if [ -z "${SSH:-}" ]; then
    if command -v ssh.exe &>/dev/null; then SSH=ssh.exe; else SSH=ssh; fi
fi
CAL="XCAL_YJYJYJ.fits YCAL_YJYJYJ.fits LCAL_YJYJYJ.fits MASTER_FLAT_YJYJYJ.fits"
CALP=$(for c in $CAL; do printf 'calib/%s ' "$c"; done)

for OB in $OBS; do
    D="$DST/${OB//:/-}"
    mkdir -p "$D/in" "$D/calib"
    echo "== $OB -> $D"
    # -h: store the files behind links; dc_ frames go to in/, calibrations stay in calib/
    $SSH "$SRV" "cd '$SRC/$OB' && tar -chf - --transform='s,^dc_,in/dc_,' dc_KMOS.*.fits $CALP" | tar -xf - -C "$D"
    $SSH "$SRV" "cd '$SRC/$OB' && md5sum dc_KMOS.*.fits | sed 's,  ,  in/,' && md5sum $CALP" > "$D/md5.remote"
    N=$(grep -c . "$D/md5.remote")
    if (cd "$D" && md5sum -c --quiet md5.remote); then
        echo "   $N files copied, md5 identical ($(du -sh "$D" | cut -f1))"
    else
        echo "   md5 MISMATCH in $D - run the script again" >&2
        exit 1
    fi
done
echo "Done. Next: python chain_band.py <OB> YJ classic|row|hyb  (OB with '-' instead of ':')"
