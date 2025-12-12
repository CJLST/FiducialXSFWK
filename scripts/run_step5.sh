#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/fit

obsName="${1//_/' vs '}"
year="$2"

python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year"
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --nnlops

cd ../coefficients

python3 pdfUncertainties.py --obsName "$obsName" --year "$year" --split --merge
