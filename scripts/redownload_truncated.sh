#!/usr/bin/env bash
# Re-download the 3 truncated _rv.fits files detected on 2026-09-24.
# Run from WSL in the pyspark_my_way directory:
#   bash scripts/redownload_truncated.sh

SERVER=kreuter@hqempower1
REMOTE_BASE=/local_storage/KOLJA/PROCDIR
LOCAL_BASE=kmos_ob_proc_server

if command -v scp.exe &>/dev/null; then SCP=scp.exe; else SCP=scp; fi

# Helper: convert local dash-timestamp name back to server colon-timestamp name
# Only the time part (THH-MM-SS) gets colons; date dashes stay.
to_server_name() {
    echo "$1" | sed 's/T\([0-9][0-9]\)-\([0-9][0-9]\)-\([0-9][0-9]\)/T\1:\2:\3/'
}

fetch() {
    local ob_local="$1"   # e.g. KMOS.2026-07-13T04-28-03.961
    local file_local="$2" # e.g. SCI_RECONSTRUCTED_dc_..._rv.fits
    local ob_remote
    ob_remote=$(to_server_name "$ob_local")
    local file_remote
    file_remote=$(to_server_name "$file_local")
    local dest="$LOCAL_BASE/$ob_local/science/$file_local"
    echo "Fetching $file_local ..."
    $SCP "$SERVER:$REMOTE_BASE/$ob_remote/science/$file_remote" "$dest"
}

fetch "KMOS.2026-07-13T04-28-03.961" \
      "SCI_RECONSTRUCTED_dc_KMOS.2026-07-13T04-28-03.961_s_i_bkg2d_rv.fits"

fetch "KMOS.2026-07-13T04-28-03.961" \
      "SCI_RECONSTRUCTED_dc_KMOS.2026-07-13T04-46-31.826_s_i_bkg2d_rv.fits"

fetch "KMOS.2026-07-13T07-18-37.747" \
      "SCI_RECONSTRUCTED_dc_KMOS.2026-07-13T07-36-36.574_s_i_bkg2d_rv.fits"

echo "Done."
