#!/usr/bin/env bash
# Download DATACUBE_COMBINED_*.fits files produced by the stacking pipeline
# from /local_storage/tmp_results/ on hqempower1.
#
# Run from WSL in the pyspark_my_way directory:
#   bash scripts/download_combined.sh
#
# Files land in:
#   kmos_proc_server_comb/<OB_name>/<timestamp>/DATACUBE_COMBINED_*.fits
# (colons replaced by dashes for Windows/NTFS compatibility)
#
# Other EDPS export folder / local target (keep re-reduced results separate,
# otherwise old and new combinations of the same OB end up side by side):
#   REMOTE_BASE=/local_storage/KOLJA/res_new_IZ_full_2 LOCAL_BASE=kmos_comb_new_IZ_full \
#       bash scripts/download_combined.sh

SERVER=kreuter@hqempower1
REMOTE_BASE=${REMOTE_BASE:-/local_storage/KOLJA/tmp_results}
LOCAL_BASE=${LOCAL_BASE:-kmos_proc_and_comb_server}


if command -v ssh.exe &>/dev/null; then SSH=ssh.exe; else SSH=ssh; fi
if command -v scp.exe &>/dev/null; then SCP=scp.exe; else SCP=scp; fi

mkdir -p "$LOCAL_BASE"

echo "=== Querying server for stacked OBs ==="
OB_DIRS=$($SSH "$SERVER" \
    "ls -1d $REMOTE_BASE/SCI_RECONSTRUCTED_dc*_s_i_bkg2d_rv 2>/dev/null | xargs -I{} basename {}")   # dc_ and dcr_ (row-corrected IZ)

if [ -z "$OB_DIRS" ]; then
    echo "No stacked OB directories found under $REMOTE_BASE."
    exit 0
fi

echo "Found: $(echo "$OB_DIRS" | wc -l | tr -d ' ') OB directories"
echo ""

NEW=0
SKIP=0

# Load into array to avoid ssh consuming the while-read stdin
mapfile -t OB_ARRAY <<< "$OB_DIRS"

for OB_DIR in "${OB_ARRAY[@]}"; do
    [ -z "$OB_DIR" ] && continue

    # Pick the most recent timestamp subdirectory that holds DATACUBE_COMBINED files:
    # every EDPS export also writes a twin folder with EDPS-internal product names,
    # and taking the newest folder blindly skipped whole datasets.
    LATEST_TS=$($SSH "$SERVER" \
        "cd '$REMOTE_BASE/$OB_DIR/' 2>/dev/null && for ts in \$(ls -1 | sort -r); do ls \"\$ts\"/DATACUBE_COMBINED_*.fits >/dev/null 2>&1 && { echo \"\$ts\"; break; }; done")
    if [ -z "$LATEST_TS" ]; then
        echo "  Skip  $OB_DIR (no timestamp subdir with DATACUBE_COMBINED files)"
        continue
    fi

    REMOTE_DIR="$REMOTE_BASE/$OB_DIR/$LATEST_TS"

    N_REMOTE=$($SSH "$SERVER" \
        "ls '$REMOTE_DIR'/DATACUBE_COMBINED_*.fits 2>/dev/null | wc -l | tr -d ' '")
    if [ "${N_REMOTE:-0}" -eq 0 ]; then
        echo "  Skip  $OB_DIR/$LATEST_TS (no DATACUBE_COMBINED files)"
        continue
    fi

    # Convert colons → dashes for local NTFS paths
    LOCAL_OB="${OB_DIR//:/-}"
    LOCAL_TS="${LATEST_TS//:/-}"
    LOCAL_DIR="$LOCAL_BASE/$LOCAL_OB/$LOCAL_TS"

    N_LOCAL=$(ls "$LOCAL_DIR"/DATACUBE_COMBINED_*.fits 2>/dev/null | wc -l | tr -d ' ')
    if [ "${N_LOCAL:-0}" -ge "$N_REMOTE" ]; then
        echo "  Skip  $LOCAL_OB  ($N_LOCAL/$N_REMOTE files already present)"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "  Fetch $OB_DIR/$LATEST_TS  ($N_REMOTE files) ..."
    mkdir -p "$LOCAL_DIR"

    mapfile -t FILE_ARRAY <<< "$($SSH "$SERVER" \
        "ls '$REMOTE_DIR'/DATACUBE_COMBINED_*.fits 2>/dev/null")"

    for REMOTE_FILE in "${FILE_ARRAY[@]}"; do
        [ -z "$REMOTE_FILE" ] && continue
        FNAME=$(basename "$REMOTE_FILE" | tr ':' '-')
        LOCAL_FILE="$LOCAL_DIR/$FNAME"
        [ -f "$LOCAL_FILE" ] && continue
        $SCP "$SERVER:$REMOTE_FILE" "$LOCAL_FILE"
    done

    echo "  Done  → $LOCAL_OB/$LOCAL_TS/"
    NEW=$((NEW + 1))

done

echo ""
echo "=== Downloaded $NEW new OB(s), skipped $SKIP ==="
