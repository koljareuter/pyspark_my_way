#!/bin/bash
# PROCDIR_HYBE run of one band (IZ or YJ) after the doubled first-run jobs of 2026-11-04 (docs/procdir-hyb-server-run.md,
# 'Run of 2026-11-04'):
#   1. repair: redo, with the original first-run cfg, the PROCDIR_HYB OB folders that got files modified since SINCE
#      (cfg_subset.py), in a working directory of its own; stops if that run reports an error
#   2. links for PROCDIR_HYBE from the (repaired) PROCDIR_HYB (prepare_procdir_ohlin.py --reuse)
#   3. the PROCDIR_HYBE reduction (hybrid_edge = 3, sky subtraction and later steps only)
# Usage on hqempower1 (one call per band; both bands may run at the same time):
#   cd ~/pyspark_my_way && nohup bash scripts/run_hybe_chain.sh IZ > /local_storage/KOLJA/hyb_work/chain_IZ.log 2>&1 &
set -u
B=$1
SINCE=${2:-"2026-11-04 00:00"}
K=/local_storage/KOLJA
P=$HOME/pyspark_my_way
export PYTHONPATH=$P
case $B in IZ|YJ) ;; *) echo "band must be IZ or YJ"; exit 1 ;; esac
if pgrep -f "kmos_server_hyb_${B}.cfg|kmos_server_hybe_${B}.cfg|repair_${B}.cfg" > /dev/null; then
  echo "a $B reduction is still running:"; pgrep -af "kmos_server_hyb_${B}.cfg|kmos_server_hybe_${B}.cfg|repair_${B}.cfg"; exit 1
fi

echo "== $(date '+%F %T') 1. repair of the $B OB folders modified since $SINCE"
cd $P
if python scripts/cfg_subset.py scripts/kmos_server_hyb_${B}.cfg $K/repair_${B}.cfg --modified-since "$SINCE"; then
  mkdir -p $K/hyb_work/${B}rep && cd $K/hyb_work/${B}rep
  python $P/scripts/reduce_kmos.py -p $K/repair_${B}.cfg > repair_${B}.out 2>&1
  n=$(grep -c "An error occ" repair_${B}.out)
  echo "   repair finished: $(grep -c 'Reducing ' repair_${B}.out) OBs, $n errors"
  if [ "$n" -gt 0 ]; then echo "   STOP: errors in $K/hyb_work/${B}rep/repair_${B}.out"; exit 1; fi
else
  echo "   nothing to repair"
fi

echo "== $(date '+%F %T') 2. links for PROCDIR_HYBE"
cd $P
python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybe_${B}.cfg --reuse $K/PROCDIR_HYB || { echo "   STOP: prepare failed"; exit 1; }

echo "== $(date '+%F %T') 3. PROCDIR_HYBE reduction"
mkdir -p $K/hyb_work/${B}e && cd $K/hyb_work/${B}e
python $P/scripts/reduce_kmos.py -p $P/scripts/kmos_server_hybe_${B}.cfg > hybe_${B}.out 2>&1
echo "   finished: $(grep -c 'Finished reducing' hybe_${B}.out) OBs, $(grep -c 'An error occ' hybe_${B}.out) errors"
cd $P && python scripts/prepare_procdir_ohlin.py scripts/kmos_server_hybe_${B}.cfg --status | tail -1
echo "== $(date '+%F %T') done"
