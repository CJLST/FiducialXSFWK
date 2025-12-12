#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/coefficients

obsName="${1//_/' vs '}"
year="$2"

python3 pdfUncertainties.py --obsName "$obsName" --year "$year" --merge

cd ..

obsName="$1"

python3 py_to_json.py $obsName $year

python3 differential_plotting.py --variables "$obsName" --config-json jsons/${obsName}_results_${year}.json --year "$year"
