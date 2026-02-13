import ROOT
#import sys, os, pwd, commands
import sys, os, pwd # spencer
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json
import importlib.util

sys.path.append('../helperstuff/')
from higgs_xsbr_13TeV import *
from binning import binning
from paths import path

sys.path.append('../inputs/')

def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',      dest='SOURCEDIR',type='string',default='./', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--doHIG', action='store_true', dest='DOHIG', default=False, help='use HIG 19 001 acceptances')
    parser.add_option('',   '--split', action='store_true', dest='SPLIT', default=False, help='')
    parser.add_option('',   '--interpolation', action='store_true', dest='INTER', default=False, help='Calculate acceptances at 124 and 126 GeV')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()



# parse the arguments and options
global opt, args
parseOptions()

def exp_xsec():
    # prepare the set of bin boundaries to run over, only 1 bin in case of the inclusive measurement
    observableBins, doubleDiff = binning(opt.OBSNAME)
    if doubleDiff:
        obs_name = opt.OBSNAME.split(' vs ')[0]
        obs_name_2nd = opt.OBSNAME.split(' vs ')[1]
        obs_name_2d = opt.OBSNAME
        doubleDiff = True
    else:
        obs_name = opt.OBSNAME
        doubleDiff = False
    ## Run for the given observable
    # obsName = opt.OBSNAME
    if doubleDiff: obsName = obs_name+'_'+obs_name_2nd
    else: obsName = obs_name
    _th_MH = opt.THEORYMASS
    print('Running Fiducial XS computation - '+obsName+' - bin boundaries: ', observableBins, '\n')
    print('Theory xsec and BR at MH = '+_th_MH)

    #_temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs_xs_136TeV','higgs4l_br'], -1)
    _temp = __import__('higgs_xsbr_13TeV', globals(), locals(), ['higgs_xs','higgs_xs_136TeV','higgs4l_br'], 0) # spencer
    if(opt.YEAR=='Run3'):
        higgs_xs = _temp.higgs_xs_136TeV
    else:
        higgs_xs = _temp.higgs_xs
    higgs4l_br = _temp.higgs4l_br

    #if opt.SPLIT:
    #    fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_'+opt.YEAR+".py"
    #else:
    if opt.INTER:
        fname = path['eos_path']+'inputs/inputs_sig_extrap_'+obsName+'_'+opt.YEAR+".py"
    else:
        fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_'+opt.YEAR+".py"

    if opt.DOHIG: fname = fname + '_HIG19001'
    #_temp = __import__(fname, globals(), locals(), ['acc'], -1)
    #_temp = __import__(fname, globals(), locals(), ['acc'], 0) # spencer
    
    module_name = os.path.splitext(os.path.basename(fname))[0]
    spec = importlib.util.spec_from_file_location(module_name, fname)
    _temp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_temp)
    
    acc = _temp.acc
    XH = []
    nBins = len(observableBins)
    if not doubleDiff: nBins = nBins - 1
    xs = {}
    for obsBin in range(nBins):
        XH.append(0.0)
        # if('mass4l' not in obsName):
        for channel in ['4e','4mu','2e2mu']:
            XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            print('Bin ', obsBin, '\t SigmaBin', obsBin, channel, ' = ', XH_fs)
            print('XH',  XH_fs - higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ggH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])
            XH[obsBin]+=XH_fs

        _obsxsec = XH[obsBin]

        print('\n')
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' = ', _obsxsec)
        xs['SigmaBin'+str(obsBin)] = _obsxsec

    with open('../inputs/xsec_'+obsName+'_'+opt.YEAR+'.py', 'w') as f:
        f.write('xsec = '+str(xs)+' \n')

exp_xsec()
