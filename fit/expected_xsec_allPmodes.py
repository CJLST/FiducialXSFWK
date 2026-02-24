import ROOT
#import sys, os, pwd, commands
import sys, os, pwd, subprocess
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
    parser.add_option('',   '--nnlops', action='store_true', dest='NNLOPS', default=False, help='nnlops prediction')
    parser.add_option('',   '--split', action='store_true', dest='SPLIT', default=False, help='use HIG 19 001 acceptances')
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
    #    fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_'+opt.YEAR+'.py'
    #else:
    if opt.INTER:
        fname = path['eos_path']+'inputs/inputs_sig_extrap_'+obsName+'_'+opt.YEAR+".py"
    else:
        fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_'+opt.YEAR+".py"

    if opt.DOHIG: fname = fname + '_HIG19001'

    module_name = os.path.splitext(os.path.basename(fname))[0]
    spec = importlib.util.spec_from_file_location(module_name, fname)
    _temp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_temp)
    
    acc = _temp.acc
    if opt.NNLOPS:

        #if opt.SPLIT:
        #    fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_NNLOPS_'+opt.YEAR+'.py'
        #else:
        if opt.INTER:
            fname = path['eos_path']+'inputs/inputs_sig_extrap_'+obsName+'_NNLOPS_'+opt.YEAR+".py"
        else:
            fname = path['eos_path']+'inputs/inputs_sig_'+obsName+'_NNLOPS_'+opt.YEAR+".py"

        module_name = os.path.splitext(os.path.basename(fname))[0]
        spec = importlib.util.spec_from_file_location(module_name, fname)
        _temp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_temp)
        
        acc_ggh = _temp.acc
    else:
        acc_ggh = acc

    if opt.NNLOPS:
        suffix = 'NNLOPS_'
    else:
        suffix = ''

    XH = []
    XH_ggh = []
    XH_vbf = []
    XH_zh = []
    XH_wh = []
    XH_ttH = []
    nBins = len(observableBins)
    if not doubleDiff: nBins = nBins - 1
    xs = {}
    xs_ggh = {}
    xs_vbf = {}
    xs_zh = {}
    xs_wh = {}
    xs_vh = {}
    xs_tth = {}
    xs_xh = {}

    # --- NEW: per-final-state storage (only used/written for mass4l) ---
    channels = ['4e','4mu','2e2mu']
    xs_ggh_fs = {ch:{} for ch in channels}
    xs_vbf_fs = {ch:{} for ch in channels}
    xs_vh_fs  = {ch:{} for ch in channels}
    xs_tth_fs = {ch:{} for ch in channels}
    xs_xh_fs  = {ch:{} for ch in channels}

    combo = '4e4mu'
    xs_ggh_fs[combo] = {}
    xs_vbf_fs[combo] = {}
    xs_vh_fs[combo]  = {}
    xs_tth_fs[combo] = {}
    xs_xh_fs[combo]  = {}

    # ---------------------------------------------------------------

    for obsBin in range(nBins):
        XH.append(0.0)
        XH_ggh.append(0.0)
        XH_vbf.append(0.0)
        XH_zh.append(0.0)
        XH_wh.append(0.0)
        XH_ttH.append(0.0)

        # --- NEW: per-channel accumulators for this bin ---
        bin_ggh_fs = {ch:0.0 for ch in channels}
        bin_vbf_fs = {ch:0.0 for ch in channels}
        bin_vh_fs  = {ch:0.0 for ch in channels}
        bin_tth_fs = {ch:0.0 for ch in channels}
        bin_xh_fs  = {ch:0.0 for ch in channels}
        # --------------------------------------------------

        for channel in ['4e','4mu','2e2mu']:
            print(channel)
            print(acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])
            xxs_ggh = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_vbf = higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_zh = higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_tth = higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            xxs_wh = higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]

            XH_fs = higgs_xs['ggH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc_ggh['ggH125_'+suffix+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['VBF_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['VBFH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['WH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['WH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ZH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ZH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH_fs += higgs_xs['ttH_'+opt.THEORYMASS]*higgs4l_br[opt.THEORYMASS+'_'+channel]*acc['ttH125_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)]
            XH[obsBin]+=XH_fs

            XH_ggh[obsBin]+=xxs_ggh
            XH_vbf[obsBin]+=xxs_vbf
            XH_zh[obsBin]+=xxs_zh
            XH_wh[obsBin]+=xxs_wh
            XH_ttH[obsBin]+=xxs_tth

            # --- NEW: fill per-channel bin sums (used only for mass4l outputs) ---
            bin_ggh_fs[channel] += xxs_ggh
            bin_vbf_fs[channel] += xxs_vbf
            bin_vh_fs[channel]  += (xxs_wh + xxs_zh)
            bin_tth_fs[channel] += xxs_tth
            bin_xh_fs[channel]  += (xxs_vbf + xxs_wh + xxs_zh + xxs_tth)
            # --------------------------------------------------------------------

        _obsxsec = XH[obsBin]

        print('\n')
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' = ', _obsxsec)
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' (ggh) = ', XH_ggh[obsBin])
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' (vbf) = ', XH_vbf[obsBin])
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' (zh) = ', XH_zh[obsBin])
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' (wh) = ', XH_wh[obsBin])
        print('Bin ', obsBin, '\t SigmaBin', obsBin, ' (tth) = ', XH_ttH[obsBin])
        print('(xcheck : ', XH_ggh[obsBin]+XH_vbf[obsBin]+XH_zh[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin],')')
        print('(XH : ', XH_vbf[obsBin]+XH_zh[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin],')')
        print('\n\n')
        xs['SigmaBin'+str(obsBin)] = _obsxsec
        xs_ggh['SigmaBin'+str(obsBin)] = XH_ggh[obsBin]
        xs_vbf['SigmaBin'+str(obsBin)] = XH_vbf[obsBin]
        xs_vh['SigmaBin'+str(obsBin)] = XH_wh[obsBin]+XH_zh[obsBin]
        xs_wh['SigmaBin'+str(obsBin)] = XH_wh[obsBin]
        xs_tth['SigmaBin'+str(obsBin)] = XH_ttH[obsBin]
        xs_xh['SigmaBin'+str(obsBin)] = XH_vbf[obsBin]+XH_wh[obsBin]+XH_ttH[obsBin]+XH_zh[obsBin]

        # --- NEW: store per-channel dict entries for this bin ---
        for ch in channels:
            xs_ggh_fs[ch]['SigmaBin'+str(obsBin)] = bin_ggh_fs[ch]
            xs_vbf_fs[ch]['SigmaBin'+str(obsBin)] = bin_vbf_fs[ch]
            xs_vh_fs[ch]['SigmaBin'+str(obsBin)]  = bin_vh_fs[ch]
            xs_tth_fs[ch]['SigmaBin'+str(obsBin)] = bin_tth_fs[ch]
            xs_xh_fs[ch]['SigmaBin'+str(obsBin)]  = bin_xh_fs[ch]

        k = 'SigmaBin'+str(obsBin)
        xs_ggh_fs['4e4mu'][k] = bin_ggh_fs['4e'] + bin_ggh_fs['4mu']
        xs_vbf_fs['4e4mu'][k] = bin_vbf_fs['4e'] + bin_vbf_fs['4mu']
        xs_vh_fs ['4e4mu'][k] = bin_vh_fs ['4e'] + bin_vh_fs ['4mu']
        xs_tth_fs['4e4mu'][k] = bin_tth_fs['4e'] + bin_tth_fs['4mu']
        xs_xh_fs ['4e4mu'][k] = bin_xh_fs ['4e'] + bin_xh_fs ['4mu']
        # --------------------------------------------------------

    with open('../inputs/xsec_'+obsName+'_'+opt.YEAR+'.py', 'w') as f:
        f.write('xsec = '+str(xs)+' \n')

    obsFull = obsName

    # helper to write fidXS files (keeps your exact content structure)
    import math

    def _is_bad(x):
        try:
            x = float(x)
        except Exception:
            return True
        return (not math.isfinite(x))  # catches nan/inf

    def write_fid_file(outname, boundaries, vals):
        # sanitize: if nominal is 0 or nan/inf -> set nominal and all variations to 0
        clean = []
        for v in vals:
            if _is_bad(v) or float(v) == 0.0:
                clean.append(0.0)
            else:
                clean.append(float(v))

        with open(outname, 'w') as f:
            f.write('Boundaries = '+str(boundaries)+'\n')
            f.write('fidXS = '+str(clean)+'\n')
            f.write('fidXS_scale_up = '+str(clean)+'\n')
            f.write('fidXS_scale_dn = '+str(clean)+'\n')
            f.write('fidXS_pdf_up = '+str(clean)+'\n')
            f.write('fidXS_pdf_dn = '+str(clean)+'\n')
            f.write('fidXS_alpha_up = '+str(clean)+'\n')
            f.write('fidXS_alpha_dn = '+str(clean)+'\n')
    
    if obsName in special_obs:

        # --- existing 4l-summed outputs (UNCHANGED behavior) ---
        write_fid_file('../inputs/fidXS_'+suffix+obsFull+'_ggH_'+opt.YEAR+'.py', observableBins, list(xs_ggh.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_VBFH_'+opt.YEAR+'.py',          observableBins, list(xs_vbf.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_VH_'+opt.YEAR+'.py',            observableBins, list(xs_vh.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_ttH_'+opt.YEAR+'.py',           observableBins, list(xs_tth.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_xH_'+opt.YEAR+'.py',            observableBins, list(xs_xh.values()))

        # --- NEW: additional per-final-state outputs (4e/4mu/2e2mu) ---
        for ch in ['4e','4mu','2e2mu','4e4mu']:
            # ggH keeps the same suffix convention as your 4l file
            write_fid_file('../inputs/fidXS_'+suffix+obsFull+'_ggH_'+ch+'_'+opt.YEAR+'.py',
                           observableBins, list(xs_ggh_fs[ch].values()))
            write_fid_file('../inputs/fidXS_'+obsFull+'_VBFH_'+ch+'_'+opt.YEAR+'.py',
                           observableBins, list(xs_vbf_fs[ch].values()))
            write_fid_file('../inputs/fidXS_'+obsFull+'_VH_'+ch+'_'+opt.YEAR+'.py',
                           observableBins, list(xs_vh_fs[ch].values()))
            write_fid_file('../inputs/fidXS_'+obsFull+'_ttH_'+ch+'_'+opt.YEAR+'.py',
                           observableBins, list(xs_tth_fs[ch].values()))
            write_fid_file('../inputs/fidXS_'+obsFull+'_xH_'+ch+'_'+opt.YEAR+'.py',
                           observableBins, list(xs_xh_fs[ch].values()))
        # --------------------------------------------------------------

    else:
    
        # --- existing non-mass4l behavior (UNCHANGED) ---
        write_fid_file('../inputs/fidXS_'+suffix+obsFull+'_ggH_'+opt.YEAR+'.py', observableBins, list(xs_ggh.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_VBFH_'+opt.YEAR+'.py',       observableBins, list(xs_vbf.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_VH_'+opt.YEAR+'.py',         observableBins, list(xs_vh.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_ttH_'+opt.YEAR+'.py',        observableBins, list(xs_tth.values()))
        write_fid_file('../inputs/fidXS_'+obsFull+'_xH_'+opt.YEAR+'.py',         observableBins, list(xs_xh.values()))

    
special_obs = [] # ['mass4l', 'massZ1', 'massZ2', 'costhetaZ1', 'costhetaZ2', 'costhetastar', 'phi', 'phi1']
exp_xsec()