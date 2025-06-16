#!/bin/bash
echo "[INFO] Starting job on $(hostname) at $(date)"
echo "[INFO] Arguments: $1 $2 $3"

cd /afs/cern.ch/user/m/mmanoni/FiducialXS/CMSSW_14_1_0_pre4/src/
cmsenv
source setup.sh
cd FiducialXSFWK/coefficients/

echo "[INFO] Environment set up."

python3 -u run_pipeline_22_23.py --obsName "$1" --obsBins "$2" --year "$3"


