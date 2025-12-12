#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/coefficients

obsName="${1//_/' vs '}"
year="$2"

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --nnlops

cd ../templates

python3 RunTemplates.py --obsName "$obsName" --year "$year"
python3 plot_templates.py --obsName "$obsName" --year "$year"

cd ../fit

python3 RunFiducialXS.py --obsName "$obsName" --year "$year"
#python3 RunFiducialXS.py --obsName "$obsName" --year "$year" --eff_unc
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year"
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --nnlops

cd ../LHScans

python3 plot_LLScan.py --obsName "$obsName" --year "$year"
