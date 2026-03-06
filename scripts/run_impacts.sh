#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/fit

obsName="${1//_/' vs '}"
year="$2"

python3 impacts.py --obsName "$obsName" --year "$year" --interpolation
python3 impacts.py --obsName "$obsName" --year "$year" --interpolation --unblind
