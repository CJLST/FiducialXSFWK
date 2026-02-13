#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/coefficients/JES

obsName="${1//_/' vs '}"
year="$2"

python3 RunJES.py --obsName "$obsName" --year "$year"
