#!/usr/bin/env bash
# Download _rv_tell.fits files for finished OBs from hqempower1 into
# kmos_ob_proc_server/<OB>/science/, on top of existing _rv.fits files.
#
# Run from WSL in the pyspark_my_way directory:
#   bash scripts/download_tellur.sh
#
# To migrate February OBs that were downloaded to the OB root before this fix:
#   for d in kmos_ob_proc_server/KMOS.*; do
#       if ls "$d"/SCI_RECONSTRUCTED_dc_*_s_i_bkg2d_rv.fits &>/dev/null; then
#           mkdir -p "$d/science"
#           mv "$d"/SCI_RECONSTRUCTED_dc_*_s_i_bkg2d_rv.fits "$d/science/"
#       fi
#   done

SERVER=kreuter@hqempower1
REMOTE_PROC=/local_storage/KOLJA/PROCDIR
LOCAL_DIR="kmos_ob_proc_server"

# On WSL2 the VPN only works through the Windows network stack,
# so prefer Windows ssh.exe/scp.exe when available.
if command -v ssh.exe &>/dev/null; then
    SSH=ssh.exe
    SCP=scp.exe
else
    SSH=ssh
    SCP=scp
fi

mkdir -p "$LOCAL_DIR"

echo "=== Querying server for finished OBs ==="
FINISHED=$($SSH "$SERVER" \
    'for d in /local_storage/KOLJA/PROCDIR/KMOS.*; do
         n=$(ls "$d/science/"*_s_i_bkg2d_rv_tell.fits 2>/dev/null | wc -l)
         [ "$n" -gt 0 ] && basename "$d"
     done')

if [ -z "$FINISHED" ]; then
    echo "No finished OBs found on server."
    exit 0
fi

echo "Finished OBs on server: $(echo "$FINISHED" | wc -l)"
echo ""

NEW=0
for OB in $FINISHED; do
    # Colons → dashes for NTFS/Windows compatibility
    LOCAL_NAME="${OB//:/-}"
    LOCAL_OB="$LOCAL_DIR/$LOCAL_NAME"
    LOCAL_SCI="$LOCAL_OB/science"

    if [ -d "$LOCAL_SCI" ] && [ "$(ls "$LOCAL_SCI"/*_rv_tell.fits 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "  Skip  $LOCAL_NAME (tell files already downloaded)"
        continue
    fi

    echo "  Fetch $OB ..."
    mkdir -p "$LOCAL_SCI"
    # Get the list of files from the server, then download each one individually
    # renaming colons to dashes (Windows/NTFS does not allow colons in filenames).
    REMOTE_FILES=$($SSH "$SERVER" "ls $REMOTE_PROC/$OB/science/*_s_i_bkg2d_rv_tell.fits 2>/dev/null")
    for REMOTE_FILE in $REMOTE_FILES; do
        LOCAL_FILE="${LOCAL_SCI}/$(basename "$REMOTE_FILE" | tr ':' '-')"
        $SCP "$SERVER:$REMOTE_FILE" "$LOCAL_FILE"
    done
    echo "  Done  → $LOCAL_NAME/science/"
    NEW=$((NEW + 1))
done

echo ""
echo "=== Downloaded $NEW new OB(s) to $LOCAL_DIR ==="
