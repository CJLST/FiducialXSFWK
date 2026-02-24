#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/coefficients

obsName="${1//_/' vs '}"
year="$2"

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --nnlops

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --interpolation --hypothesis "24"
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --nnlops --interpolation --hypothesis "24"

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --interpolation --hypothesis "26"
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --merge --nnlops --interpolation --hypothesis "26"

python3 RunInterpolation.py --obsName "$obsName" --year "$year"
python3 RunInterpolation.py --obsName "$obsName" --year "$year" --nnlops

#python3 RunPlotCoefficients.py --obsName "$obsName" --year "$year" --interpolation
#python3 RunPlotCoefficients.py --obsName "$obsName" --year "$year" --nnlops --interpolation

cd ../templates

python3 RunTemplates.py --obsName "$obsName" --year "$year"
#python3 plot_templates.py --obsName "$obsName" --year "$year"

cd ../fit

python3 RunFiducialXS.py --obsName "$obsName" --year "$year" --eff_unc --interpolation --NOK1K2
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --interpolation
python3 expected_xsec_allPmodes.py --obsName "$obsName" --year "$year" --nnlops --interpolation
python3 impacts.py --obsName "$obsName" --year "$year" --interpolation --NOK1K2

python3 RunFiducialXS.py --obsName "$obsName" --year "$year" --eff_unc --interpolation --unblind --NOK1K2
python3 impacts.py --obsName "$obsName" --year "$year" --interpolation --unblind --NOK1K2

cd ../LHScans

python3 plot_LLScan.py --obsName "$obsName" --year "$year" --interpolation --unblind
