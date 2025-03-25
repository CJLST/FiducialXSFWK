#!/bin/bash
cd /afs/cern.ch/work/l/lurda/HZZ/StatisticalAnalysis/CMSSW_14_1_0_pre4/src/FiducialXSFWK/coefficients
cmsenv

#python3 RunCoefficients.py --obsName '$1' --year '$2' --interpolation --hypothesis '$3' ZZFLOATING
python3 RunCoefficients_.py --obsName 'pT4l' --year 'Run3'
