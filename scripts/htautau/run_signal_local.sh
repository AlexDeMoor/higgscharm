#!/bin/bash
# Local (lxplus) production of the H->tautau workflows over the private signal samples.
#
# Columnar reads straight from EOS are latency bound (~10 MB/s per node for the whole
# process pool) and the local /tmp disk of the lxplus VMs sustains only ~100 random
# reads/s, so each sample is processed in parts of NPART files copied to /dev/shm (RAM):
# copy a part (~130 MB/s), run all workflows on it in parallel (CPU bound, ~13k events/s
# per workflow with 5 workers), remove it, next part. One part = ~17 GB of tmpfs; keep
# NPART x 2 GB well below the free memory (check `df -h /dev/shm`, `free -g`).
#
#   bash scripts/htautau/run_signal_local.sh [workflows] [datasets]
#   bash scripts/htautau/run_signal_local.sh "htautau_mutau htautau_etau htautau_tautau" "HplusCharm_HtoTauTau HplusBottom_HtoTauTau"
#
# Outputs: outputs/<workflow>/2024/<dataset>/<dataset>_<part>.coffea (condor layout) and
# outputs/<workflow>/2024/<dataset>/local_parts.txt (expected part tags, checked by
# run_postprocess.py), then
#   python3 run_postprocess.py -w <workflow> -y 2024 --postprocess --output_format coffea
#
# Any staging failure or failed part is fatal for that dataset (exit code 1 at the end);
# stale <dataset>_<part>.coffea files of a previous run are removed before processing.
set -u -o pipefail
WORKFLOWS=${1:-"htautau_mutau htautau_etau htautau_tautau"}
DATASETS=${2:-"HplusCharm_HtoTauTau HplusBottom_HtoTauTau"}
YEAR=${YEAR:-2024}
NPART=${NPART:-9}
WORKERS=${WORKERS:-5}
# python with coffea 0.7.30: $PY if set (see README.md), else the condor image through apptainer
CONDA_PY=/afs/cern.ch/work/a/ademoor/conda/conda_envs/btv_coffea/bin/python3
IMG=${IMG:-/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base-almalinux9:0.7.30-py3.10}
if [ -z "${PY:-}" ]; then
  if [ -x $CONDA_PY ]; then PY=$CONDA_PY
  else PY="apptainer exec -B /afs -B /eos -B /cvmfs -B /run/user/$(id -u) $IMG python3"; fi
fi
SHM=${SHM:-/dev/shm/$USER/stage}
LOGDIR=${LOGDIR:-logs/htautau_local}
OUTROOT=${OUTROOT:-/eos/user/${USER:0:1}/$USER/higgscharm/outputs}
FILESET=analysis/filesets/fileset_${YEAR}_nanov15_lxplus.json
mkdir -p $LOGDIR
LOG=$LOGDIR/run_signal_local.log
export XRD_NETWORKSTACK=IPv4
export XRD_RUNFORKHANDLER=1
status=0

stage_file() {  # stage_file <url> <dst>: 3 attempts, never leaves a truncated file behind
  local url=$1 dst=$2 attempt
  [ -f "$dst" ] && return 0
  for attempt in 1 2 3; do
    xrdcp -f -s "$url" "$dst.part" && mv "$dst.part" "$dst" && return 0
    rm -f "$dst.part"
    sleep 5
  done
  return 1
}

for d in $DATASETS; do
  mapfile -t FILES < <($PY -c "import json; print('\n'.join(json.load(open('$FILESET'))['$d']))")
  n=${#FILES[@]}
  if [ $n -eq 0 ]; then echo "[$(date)] no files for $d in $FILESET" | tee -a $LOG; status=1; continue; fi
  nparts=$(( (n + NPART - 1) / NPART ))
  dataset_ok=1
  for wf in $WORKFLOWS; do
    outdir=$OUTROOT/$wf/$YEAR/$d
    mkdir -p $outdir
    # stale parts of a previous run would be re-counted by the postprocess
    rm -f $outdir/${d}_*.coffea
    seq 1 $nparts | sed "s/^/${d}_/" > $outdir/local_parts.txt
  done
  echo "[$(date)] $d: $n files -> $nparts parts of $NPART" | tee -a $LOG
  part=0
  for ((i=0; i<n; i+=NPART)); do
    part=$((part+1)); mkdir -p $SHM/$d; local_files=""; staged_ok=1
    echo "[$(date)] $d part $part/$nparts: staging files $i-$((i+NPART-1)) of $n to $SHM" | tee -a $LOG
    for f in "${FILES[@]:i:NPART}"; do
      dst=$SHM/$d/$(basename $f)
      if stage_file "$f" "$dst"; then
        local_files="$local_files $dst"
      else
        echo "[$(date)] STAGING FAILED $f" | tee -a $LOG; staged_ok=0
      fi
    done
    if [ $staged_ok -eq 0 ]; then
      echo "[$(date)] $d part $part skipped (staging failure)" | tee -a $LOG
      dataset_ok=0; rm -rf $SHM/$d; continue
    fi
    pids=""; wfs=""
    for wf in $WORKFLOWS; do
      $PY run_local.py -w $wf -y $YEAR -d $d --files $local_files --eos --workers $WORKERS --chunksize 50000 --tag _$part > $LOGDIR/${wf}_${d}_$part.log 2>&1 &
      pids="$pids $!"; wfs="$wfs $wf"
    done
    set -- $wfs
    for pid in $pids; do
      wf=$1; shift
      if wait $pid; then
        echo "[$(date)] OK $wf $d part $part" | tee -a $LOG
      else
        echo "[$(date)] FAILED $wf $d part $part (see $LOGDIR/${wf}_${d}_$part.log)" | tee -a $LOG
        dataset_ok=0
      fi
    done
    rm -rf $SHM/$d
  done
  if [ $dataset_ok -eq 1 ]; then
    echo "[$(date)] finished $d ($nparts parts)" | tee -a $LOG
  else
    echo "[$(date)] finished $d WITH FAILURES" | tee -a $LOG; status=1
  fi
done
if [ $status -eq 0 ]; then echo "[$(date)] ALL DONE" | tee -a $LOG; else echo "[$(date)] DONE WITH FAILURES" | tee -a $LOG; fi
exit $status
