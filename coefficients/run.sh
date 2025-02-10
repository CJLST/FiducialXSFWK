#!/bin/bash
cd /afs/cern.ch/work/l/lurda/HZZ/CMSSW_14_1_0_pre4/src/FiducialXSFWK/coefficients
cmsenv

python3 RunCoefficients.py --obsName "$1" --year "$2"
