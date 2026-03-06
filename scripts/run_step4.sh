#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/fit

#obsName="${1//_/ }"
obsName="${1//_/' vs '}"
#obsName="$1"
year="$2"

python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --split #--interpolation
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --nnlops --split #--interpolation

cd ../coefficients

python3 pdfUncertainties.py --obsName "$obsName" --year "$year" --split
