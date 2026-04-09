#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv

obsName="${1//_/' vs '}"
#obsName="$1"
year="$2"

cd FiducialXSFWK/fit

#python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --interpolation --ZZfloating
#python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --nnlops --interpolation --ZZfloating

cd ../coefficients

#python3 pdfUncertainties.py --obsName "$obsName" --year "$year" --merge --ZZfloating

cd ..

obsName="$1"

#python3 py_to_json.py $obsName $year 

python3 differential_plotting.py --variables "$obsName" --config-json jsons/${obsName}_results_${year}.json --year "$year" #--no-preliminary
