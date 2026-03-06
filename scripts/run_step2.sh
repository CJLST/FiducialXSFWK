#!/bin/bash

cd /afs/cern.ch/user/s/sellissp/public/HZZ/CMSSW_14_1_0_pre4/src
source setup.sh
cmsenv
cd FiducialXSFWK/coefficients

obsName="${1//_/' vs '}"
#obsName="$1"
year="$2"

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge --nnlops

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge --interpolation --hypothesis "24"
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge --nnlops --interpolation --hypothesis "24"

python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge --interpolation --hypothesis "26"
python3 RunCoefficients.py --obsName "$obsName" --year "$year" --split --merge --nnlops --interpolation --hypothesis "26"

python3 RunInterpolation.py --obsName "$obsName" --year "$year"
python3 RunInterpolation.py --obsName "$obsName" --year "$year" --nnlops
