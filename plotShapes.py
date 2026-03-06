import sys, os, string, re, pwd,  ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
from paths import path
import ROOT
ROOT.gROOT.SetBatch(True)

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

### Define function for parsing options
def parseOptions():

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-d', '--dir',    dest='SOURCEDIR',  type='string',default='./combine_files/', help='run from the SOURCEDIR as working area, skip if SOURCEDIR is an empty string')
    parser.add_option('',   '--asimovModel',dest='ASIMOV',type='string',default='SM_125', help='Name of the asimov data mode')
    parser.add_option('',   '--asimovMass',dest='ASIMOVMASS',type='string',default='125.38', help='Asimov Mass')
    parser.add_option('',   '--unfoldModel',dest='UNFOLD',type='string',default='SM_125', help='Name of the unfolding model')
    parser.add_option('',   '--obsName',dest='OBSNAME',    type='string',default='mass4l',   help='Name of the observalbe, supported: "inclusive", "pT", "eta", "Njets"')
    parser.add_option('',   '--obsBins',dest='OBSBINS',    type='string',default='|105|160|',   help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
    parser.add_option('',   '--fixFrac', action='store_true', dest='FIXFRAC', default=False, help='Use results from fixed fraction fit, default is False')
    parser.add_option('',   '--unblind', action='store_true', dest='UNBLIND', default=False, help='Use real data')
    parser.add_option('',   '--prefit', action='store_true', dest='PREFIT', default=False, help='Prefit plots')
    parser.add_option('',   '--bonly', action='store_true', dest='BONLY', default=False, help='Bonly postfit plots')
    parser.add_option('',   '--theoryMass',dest='THEORYMASS',    type='string',default='125.38',   help='Mass value for theory prediction')
    parser.add_option('',   '--year',  dest='YEAR',  type='string',default='Run3',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--m4lLower',  dest='LOWER_BOUND',  type='int',default=105.0,   help='Lower bound for m4l')
    parser.add_option('',   '--m4lUpper',  dest='UPPER_BOUND',  type='int',default=160.0,   help='Upper bound for m4l')
    parser.add_option('',   '--interpolation', action='store_true', dest='INTER', default=False, help='Calculate acceptances at 124 and 126 GeV')

    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# parse the arguments and options
global opt, args, runAllSteps
parseOptions()
sys.argv = grootargs

from ROOT import *
from tdrStyle import *
setTDRStyle()

def checkDir(folder_path):
    os.makedirs(folder_path, exist_ok=True)

def generateName(_year, _fStateNumber, _recobin, _fState, _bin, _physicalModel, _observableBins, _obsName):
    #_years = {"1":"2022", "2":"2022EE", "3":"2018", "2018":"2018", "2022":"2022", "2022EE":"2022EE"}
    print("_physicalModel", _physicalModel)
    print("_year", _year)
    
    if _physicalModel == 'v4':
        print("generateName v4 model")
        binName = "ch"+_fStateNumber+"_ch"+str(_recobin+1)
        print("binName", binName)
        procName = _fState+"Bin"+str(_bin)
        
        return binName, procName

    elif _physicalModel != 'v3':
        # In case of mass4l we have only one bin, so the third 'chan' part is not included in the name of the function
        if (_obsName!='mass4l' and _obsName!='mass4l_zzfloating'):
            binName = "ch"+_year+"_ch"+_fStateNumber+"_ch"+str(_recobin+1)
            procName = _fState+"Bin"+str(_bin)
        else:
            binName = "ch"+_year+"_ch"+_fStateNumber
            procName = _fState+"Bin"+str(_bin)
        return binName, procName

    else:
        _obsName_v3 = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'PTJET', 'njets_pt30_eta4p7': 'NJ'}
        if _obsName not in _obsName_v3:
            _obsName_v3[_obsName] = _obsName

        if '_' in obsName and not 'floating' in obsName and not 'kL' in obsName and obsName != 'njets_pt30_eta4p7':
            _recobin_final = str(_observableBins[_recobin][0]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][1]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][2]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin][3]).replace('.', 'p').replace('-','m')
            _genbin_final = str(_observableBins[_bin][0]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][1]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][2]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin][3]).replace('.', 'p').replace('-','m')
        else:
            _recobin_final = str(_observableBins[_recobin]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_recobin+1]).replace('.', 'p').replace('-','m')
            if int(_observableBins[_recobin+1]) > 1000:
                _recobin_final = 'GT'+str(int(_observableBins[_recobin]))

            _genbin_final = str(_observableBins[_bin]).replace('.', 'p').replace('-','m')+'_'+str(_observableBins[_bin+1]).replace('.', 'p').replace('-','m')
            if int(_observableBins[_bin+1]) > 1000:
                _genbin_final = 'GT'+str(int(_observableBins[_bin]))

        # print(_obsName, _year, _years[_year])
        binName = "hzz_" + _obsName_v3[_obsName] + "_" + _recobin_final + "_cat" + _fState + "_" + _year
        procName = _obsName_v3[_obsName] + "_" + _genbin_final
        return binName, procName



def plotAsimov_sim(modelName, physicalModel, obsName, fstate, observableBins, recobin):
    print("Entering in plotAsimov_sim......")
    eos_inputs_path = path["eos_path"]
    sourcedir = opt.SOURCEDIR
    theorymass = opt.THEORYMASS
    year = opt.YEAR
    if year == '2016':
        lumi = '36.3'
        years = [""]
    elif year == '2017':
        lumi = '41.5'
        years = [""]
    elif year == '2018':
        lumi = '59.7'
        years = ["2018"]
    elif year == 'Full':
        lumi = '138'
        years = ["1", "2", "3"]
    elif year == '2022':
        lumi = '7.98'
        years = ["2022"]
    elif year == '2022EE':
        lumi = '26.67'
        years = ["2022EE"]
    elif year == '2022full':
        lumi = '34.65'
        years = ["2022", "2022EE"]
    elif year == '2023preBPix':
        lumi = '18.06'
        years = ["2023preBPix"]
    elif year == '2023postBPix':
        lumi = '9.93'
        years = ["2023postBPix"]
    elif year == '2024':
        lumi = '108.822'
        years = ["2024"]
    elif year == 'Run3':
        lumi = '171.22'
        years = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]

    print("year", year, "lumi", lumi)
#    nBins = len(observableBins)
    # if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries

    channel = {"4mu":"1", "4e":"2", "2e2mu":"3", "4l":"2"} # 4l is dummy, won't be used
    run = {"2016":"1", "2017":"2", "2018":"3", "Full": "2", "Run3": "2", "2022": "2"}
    SignalNames = {"v3":"smH_", "v2":"trueH", "v4":"trueH", "kLambda": "trueH"}
    CombNames = {"v3":"nonResH", "v2":"fakeH", "v4":"fakeH", "kLambda":"fakeH"}
    OutNames = {"v3":"OutsideAcceptance", "v2":"out_trueH", "v4":"out_trueH", "kLambda":"out_trueH"}

    # Load some libraries
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

    if(not opt.UNBLIND):
        theorymass = theorymass + '.123456'
    if physicalModel == 'v3':
        fname = 'higgsCombine_'+obsName+'_r_smH_0.MultiDimFit.mH'+theorymass+'.root'
    elif physicalModel == 'kLambda':
        fname = 'higgsCombine_'+obsName+'.MultiDimFit.mH'+theorymass+'.root'
    else:
        fname = 'higgsCombine_'+obsName+'_r2e2muBin0.MultiDimFit.mH'+theorymass+'.root'
    # print('combine file: ', fname)

    print("fname", fname)

    f_asimov = TFile(sourcedir + fname, "READ")


    #f_asimov = TFile(os.path.join(eos_inputs_path, fname), "READ")

    print("f_asimov", f_asimov)

    if (not opt.UNBLIND):
        data = f_asimov.Get("toys/toy_asimov")
    w_asimov = f_asimov.Get("w")
    w_asimov.loadSnapshot("clean")
    if (opt.UNBLIND):
        data = w_asimov.data("data_obs")
        w_asimov.loadSnapshot("MultiDimFit")
    w_asimov.loadSnapshot("clean")

    print("-------Workspace----------")
    #w_asimov.Print("v")

    ##### ------------------------ Normalisation values for Asimov dataset ------------------------ #####
    print("-------Normalisation values for Asimov dataset----------")
    trueH_asimov = {}
    zjets_asimov = {}
    ggzz_asimov = {}
    fakeH_asimov = {}
    out_trueH_asimov = {}
    qqzz_asimov = {}
    n_trueH_asimov = {}
    for year in years:
        print("year", year)
        n_trueH_asimov["4l_"+year] = 0.0
        print("n_trueH_asimov", n_trueH_asimov)
    n_trueH_otherfid_asimov = {}
    for year in years:
        n_trueH_otherfid_asimov["4l_"+year] = 0.0
    n_zjets_asimov = {}
    for year in years:
        n_zjets_asimov["4l_"+year] = 0.0
    n_ggzz_asimov = {}
    for year in years:
        n_ggzz_asimov["4l_"+year] = 0.0
    n_fakeH_asimov = {}
    for year in years:
        n_fakeH_asimov["4l_"+year] = 0.0
    n_out_trueH_asimov = {}
    for year in years:
        n_out_trueH_asimov["4l_"+year] = 0.0
    n_qqzz_asimov = {}
    for year in years:
        n_qqzz_asimov["4l_"+year] = 0.0
    n_zz_asimov = {}
    for year in years:
        n_zz_asimov["4l_"+year] = 0.0

    fStates = ['4mu','4e','2e2mu']
    for year in years:
        for fState in fStates:
            print("fState", fState)
            for bin in range(nBins):
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                print("bin_name", bin_name)
                print("process_name", process_name)
                print("--->", "n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                trueH_asimov[fState+"_"+year+"Bin"+str(bin)] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)
                print("trueH_asimov", trueH_asimov)

            zjets_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_zjets")
            print("zjets_asimov[fState+year]", zjets_asimov[fState+"_"+year])
            ggzz_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_ggzz")
            fakeH_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+CombNames[physicalModel])
            out_trueH_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_"+OutNames[physicalModel])
            qqzz_asimov[fState+"_"+year] = w_asimov.function("n_exp_final_bin"+bin_name+"_proc_bkg_qqzz")

            n_trueH_otherfid_asimov[fState+"_"+year] = 0.0
            for bin in range(nBins):
                if (bin==recobin): 
                    print("recobin:",recobin, "  fState:",fState,"   year", year)
                    n_trueH_asimov[fState+"_"+year] = trueH_asimov[fState+"_"+year+"Bin"+str(bin)].getVal()

                else: n_trueH_otherfid_asimov[fState+"_"+year] += trueH_asimov[fState+"_"+year+"Bin"+str(bin)].getVal()
            n_zjets_asimov[fState+"_"+year] = zjets_asimov[fState+"_"+year].getVal()
            n_ggzz_asimov[fState+"_"+year] = ggzz_asimov[fState+"_"+year].getVal()
            n_fakeH_asimov[fState+"_"+year] = fakeH_asimov[fState+"_"+year].getVal()
            n_out_trueH_asimov[fState+"_"+year] = out_trueH_asimov[fState+"_"+year].getVal()
            n_qqzz_asimov[fState+"_"+year] = qqzz_asimov[fState+"_"+year].getVal()
            n_zz_asimov[fState+"_"+year] = n_ggzz_asimov[fState+"_"+year]+n_qqzz_asimov[fState+"_"+year]
            n_trueH_asimov["4l_"+year] += n_trueH_asimov[fState+"_"+year]
            n_trueH_otherfid_asimov["4l_"+year] += n_trueH_otherfid_asimov[fState+"_"+year]
            n_zjets_asimov["4l_"+year] += zjets_asimov[fState+"_"+year].getVal()
            n_ggzz_asimov["4l_"+year] += ggzz_asimov[fState+"_"+year].getVal()
            n_fakeH_asimov["4l_"+year] += fakeH_asimov[fState+"_"+year].getVal()
            n_out_trueH_asimov["4l_"+year] += out_trueH_asimov[fState+"_"+year].getVal()
            n_qqzz_asimov["4l_"+year] += qqzz_asimov[fState+"_"+year].getVal()
            n_zz_asimov["4l_"+year] += n_ggzz_asimov[fState+"_"+year]+n_qqzz_asimov[fState+"_"+year]

    # Sum over the years
    for fState in fStates:
        n_trueH_asimov[fState] = 0.0
        n_trueH_otherfid_asimov[fState] = 0.0
        n_fakeH_asimov[fState] = 0.0
        n_out_trueH_asimov[fState] = 0.0
        n_zz_asimov[fState] = 0.0
        n_qqzz_asimov[fState] = 0.0
        n_zjets_asimov[fState] = 0.0
        for year in years:
            n_trueH_asimov[fState] += n_trueH_asimov[fState+"_"+year]
            n_trueH_otherfid_asimov[fState] += n_trueH_otherfid_asimov[fState+"_"+year]
            n_fakeH_asimov[fState] += n_fakeH_asimov[fState+"_"+year]
            n_out_trueH_asimov[fState] += n_out_trueH_asimov[fState+"_"+year]
            n_zz_asimov[fState] += n_zz_asimov[fState+"_"+year]
            n_qqzz_asimov[fState] += n_qqzz_asimov[fState+"_"+year]
            n_zjets_asimov[fState] += n_zjets_asimov[fState+"_"+year]

    # Sum over the final states
    n_trueH_asimov["4l"] = 0.0
    n_trueH_otherfid_asimov["4l"] = 0.0
    n_fakeH_asimov["4l"] = 0.0
    n_out_trueH_asimov["4l"] = 0.0
    n_zjets_asimov["4l"] = 0.0
    n_zz_asimov["4l"] = 0.0
    n_qqzz_asimov["4l"] = 0.0
    for year in years:
        n_trueH_asimov["4l"] += n_trueH_asimov["4l_"+year]
        n_trueH_otherfid_asimov["4l"] += n_trueH_otherfid_asimov["4l_"+year]
        n_fakeH_asimov["4l"] += n_fakeH_asimov["4l_"+year]
        n_out_trueH_asimov["4l"] += n_out_trueH_asimov["4l_"+year]
        n_zjets_asimov["4l"] += n_zjets_asimov["4l_"+year]
        n_qqzz_asimov["4l"] += n_zz_asimov["4l_"+year]
        n_zz_asimov["4l"] += n_zz_asimov["4l_"+year]

    ##### ------------------------ Normalisation values for modelfit ------------------------ #####
    print("--------Normalisation values for modelfit------------")
    f_modelfit = TFile(sourcedir + fname, "READ")
    w_modelfit = f_modelfit.Get("w")
    w_prefit = f_modelfit.Get("w")
    sim_prefit = w_modelfit.pdf("model_s")
    w_prefit.loadSnapshot("clean")
    if opt.BONLY:
        sim = w_modelfit.pdf("model_b")
    else:
        sim = w_modelfit.pdf("model_s")
    
    print("sim", sim)
    CMS_zz4l_mass = w_modelfit.var("CMS_zz4l_mass")
    if not opt.PREFIT:
        w_modelfit.loadSnapshot("MultiDimFit")
    # w_modelfit.loadSnapshot("clean")

    trueH_modelfit = {}
    zjets_modelfit = {}
    ggzz_modelfit = {}
    fakeH_modelfit = {}
    out_trueH_modelfit = {}
    qqzz_modelfit = {}
    n_trueH_modelfit = {}
    for year in years:
        n_trueH_modelfit["4l_"+year] = 0.0
    n_trueH_otherfid_modelfit = {}
    for year in years:
        n_trueH_otherfid_modelfit["4l_"+year] = 0.0
    n_zjets_modelfit = {}
    for year in years:
        n_zjets_modelfit["4l_"+year] = 0.0
    n_ggzz_modelfit = {}
    for year in years:
        n_ggzz_modelfit["4l_"+year] = 0.0
    n_fakeH_modelfit = {}
    for year in years:
        n_fakeH_modelfit["4l_"+year] = 0.0
    n_out_trueH_modelfit = {}
    for year in years:
        n_out_trueH_modelfit["4l_"+year] = 0.0
    n_qqzz_modelfit = {}
    for year in years:
        n_qqzz_modelfit["4l_"+year] = 0.0
    n_zz_modelfit = {}
    for year in years:
        n_zz_modelfit["4l_"+year] = 0.0

    fStates = ['4mu','4e','2e2mu']
    for year in years:
        for fState in fStates:
            print("fState", fState)
            for bin in range(nBins):
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                trueH_modelfit[fState+"_"+year+"Bin"+str(bin)] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+SignalNames[physicalModel]+process_name)

            zjets_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_zjets")
            print("zjets_modelfit[fState+year]", zjets_modelfit[fState+"_"+year])
            ggzz_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_ggzz")
            fakeH_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+CombNames[physicalModel])
            out_trueH_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_"+OutNames[physicalModel])
            qqzz_modelfit[fState+"_"+year] = w_modelfit.function("n_exp_final_bin"+bin_name+"_proc_bkg_qqzz")

            n_trueH_otherfid_modelfit[fState+"_"+year] = 0.0
            for bin in range(nBins):
                if (bin==recobin): n_trueH_modelfit[fState+"_"+year] = trueH_modelfit[fState+"_"+year+"Bin"+str(bin)].getVal()
                else: n_trueH_otherfid_modelfit[fState+"_"+year] += trueH_modelfit[fState+"_"+year+"Bin"+str(bin)].getVal()
            n_zjets_modelfit[fState+"_"+year] = zjets_modelfit[fState+"_"+year].getVal()
            n_ggzz_modelfit[fState+"_"+year] = ggzz_modelfit[fState+"_"+year].getVal()
            n_fakeH_modelfit[fState+"_"+year] = fakeH_modelfit[fState+"_"+year].getVal()
            n_out_trueH_modelfit[fState+"_"+year] = out_trueH_modelfit[fState+"_"+year].getVal()
            n_qqzz_modelfit[fState+"_"+year] = qqzz_modelfit[fState+"_"+year].getVal()
            n_zz_modelfit[fState+"_"+year] = n_ggzz_modelfit[fState+"_"+year]+n_qqzz_modelfit[fState+"_"+year]
            n_trueH_modelfit["4l_"+year] += n_trueH_modelfit[fState+"_"+year]
            n_trueH_otherfid_modelfit["4l_"+year] += n_trueH_otherfid_modelfit[fState+"_"+year]
            n_zjets_modelfit["4l_"+year] += zjets_modelfit[fState+"_"+year].getVal()
            n_ggzz_modelfit["4l_"+year] += ggzz_modelfit[fState+"_"+year].getVal()
            n_fakeH_modelfit["4l_"+year] += fakeH_modelfit[fState+"_"+year].getVal()
            n_out_trueH_modelfit["4l_"+year] += out_trueH_modelfit[fState+"_"+year].getVal()
            n_qqzz_modelfit["4l_"+year] += qqzz_modelfit[fState+"_"+year].getVal()
            n_zz_modelfit["4l_"+year] += n_ggzz_modelfit[fState+"_"+year]+n_qqzz_modelfit[fState+"_"+year]

    # Sum over the years
    for fState in fStates:
        n_trueH_modelfit[fState] = 0.0
        n_trueH_otherfid_modelfit[fState] = 0.0
        n_fakeH_modelfit[fState] = 0.0
        n_out_trueH_modelfit[fState] = 0.0
        n_zz_modelfit[fState] = 0.0
        n_qqzz_modelfit[fState] = 0.0
        n_zjets_modelfit[fState] = 0.0
        for year in years:
            n_trueH_modelfit[fState] += n_trueH_modelfit[fState+"_"+year]
            n_trueH_otherfid_modelfit[fState] += n_trueH_otherfid_modelfit[fState+"_"+year]
            n_fakeH_modelfit[fState] += n_fakeH_modelfit[fState+"_"+year]
            n_out_trueH_modelfit[fState] += n_out_trueH_modelfit[fState+"_"+year]
            n_zz_modelfit[fState] += n_zz_modelfit[fState+"_"+year]
            n_qqzz_modelfit[fState] += n_qqzz_modelfit[fState+"_"+year]
            n_zjets_modelfit[fState] += n_zjets_modelfit[fState+"_"+year]

    # Sum over the final states
    n_trueH_modelfit["4l"] = 0.0
    n_trueH_otherfid_modelfit["4l"] = 0.0
    n_fakeH_modelfit["4l"] = 0.0
    n_out_trueH_modelfit["4l"] = 0.0
    n_zjets_modelfit["4l"] = 0.0
    n_zz_modelfit["4l"] = 0.0
    n_qqzz_modelfit["4l"] = 0.0
    for year in years:
        n_trueH_modelfit["4l"] += n_trueH_modelfit["4l_"+year]
        n_trueH_otherfid_modelfit["4l"] += n_trueH_otherfid_modelfit["4l_"+year]
        n_fakeH_modelfit["4l"] += n_fakeH_modelfit["4l_"+year]
        n_out_trueH_modelfit["4l"] += n_out_trueH_modelfit["4l_"+year]
        n_zjets_modelfit["4l"] += n_zjets_modelfit["4l_"+year]
        n_qqzz_modelfit["4l"] += n_zz_modelfit["4l_"+year]
        n_zz_modelfit["4l"] += n_zz_modelfit["4l_"+year]


    ##### ------------------------ Data ------------------------ #####
    print("--------Data---------")
    
    CMS_channel = w_asimov.cat("CMS_channel")
    print("1---CMS_channel", CMS_channel)
    mass = w_asimov.var("ZZMass")
    mass= mass.frame(RooFit.Bins(30))

    massVar = w_asimov.var("ZZMass")
    print("massVar", massVar)
    if massVar is None:
        raise RuntimeError("ZZMass not found in workspace")

    print("fstate", fstate)
    if (fstate=="4l"):
        datacut = ''
        comp_otherfid = ''
        for bin in range(nBins):
            if bin!=recobin: continue
            for year in years:
                for fState in fStates:
                    bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                    comp_otherfid += "shapeSig_"+SignalNames[physicalModel]+process_name+"_"+bin_name+","

        for year in years:
            for fState in fStates:
                bin_name, process_name = generateName(year, channel[fState], recobin, fState, bin, physicalModel, observableBins, obsName)
                ch_idx = str(CMS_channel.lookupIndex(bin_name))
                #print("1---ch_idx", ch_idx)
                datacut += "CMS_channel=="+ch_idx+" || "
                #print("datacut", datacut)
        datacut = datacut.rstrip(" || ")
        comp_otherfid = comp_otherfid.rstrip(',')
        data = data.reduce(RooFit.Cut(datacut))
        data.plotOn(mass)
        # Plot the combined 4l categories simultaneously using datacut, NO Slice
        sim.plotOn(
        mass,                     # Select all 4l categories combined
        RooFit.ProjWData(CMS_channel, data, True),
        RooFit.LineColor(kOrange-3)
        )
        
        if opt.BONLY:
            sim_prefit.plotOn(
            mass,
            RooFit.ProjWData(CMS_channel, data, True),
            RooFit.Components(comp_otherfid),
            RooFit.LineColor(kRed)
            )

    else:
        datacut = ''
        comp_otherfid = ''
        for bin in range(nBins):
            if bin!=recobin: continue
            for year in years:
                print("year", year)
                print("quiiii")
                bin_name, process_name = generateName(year, channel[fstate], recobin, fstate, bin, physicalModel, observableBins, obsName)
                comp_otherfid += "shapeSig_"+SignalNames[physicalModel]+process_name+"_"+bin_name+","
        for year in years:
           bin_name, process_name = generateName(year, channel[fstate], recobin, fstate, bin, physicalModel, observableBins, obsName)
           ch_idx = str(CMS_channel.lookupIndex(bin_name))
           print("1---ch_idx", ch_idx)
           datacut += "CMS_channel=="+ch_idx+" || "
           print("datacut", datacut)
           cat_val = CMS_channel.getLabel()

           idx_to_label = {}
           for i in range(CMS_channel.numBins()):
                CMS_channel.setIndex(i)
                idx_to_label[i] = CMS_channel.getLabel()
           cat_val = idx_to_label[int(ch_idx)]
           print("cat_val", cat_val)

        comp_otherfid = comp_otherfid.rstrip(',')
        print("comp_otherfid", comp_otherfid)
        datacut = datacut.rstrip(" || ")
        print("datacut", datacut)
        data = data.reduce(RooFit.Cut(datacut))
        print("Variables in 'data':")
        data.get().Print("v")
        print("data", data)
        print("mass", mass)
        data.plotOn(mass)
        print("data enetries", data.numEntries())
        print()

        print("mass", mass)
        print("RooFit.ProjWData(data,True)", RooFit.ProjWData(data,True))
        sim.plotOn(
        mass,
        #RooFit.Slice(CMS_channel, cat_val),                    # tell the sim‐PDF which channel to draw
        RooFit.ProjWData(CMS_channel, data, True),             # <— category first, dataset second
        RooFit.LineColor(kOrange-3)
        )

        if opt.BONLY:
            sim_prefit.plotOn(
            mass,
            RooFit.ProjWData(CMS_channel, data, True),
            RooFit.Components(comp_otherfid),
            RooFit.LineColor(kRed)
            )

    ##### ------------------------ Shapes ------------------------ #####
    print("-------Shapes----------") 

    # --- helper to build comma‑separated component lists --------------------------
    def make_comp_list(exprs):
        """Take a list of strings, drop empties, return a comma‑joined string."""
        return ",".join([e for e in exprs if e])

    # --- signal in other fiducial bins --------------------------------------------
    comp_otherfid_parts = []
    for bin in range(nBins):
        if bin == recobin:
            continue
        for year in years:
            if fstate != "4l":                                         # single final state
                bin_name, proc = generateName(year, channel[fstate], recobin,
                                            fstate, bin, physicalModel,
                                            observableBins, obsName)
                comp_otherfid_parts.append(f"shapeSig_{SignalNames[physicalModel]}{proc}_{bin_name}")
            else:                                                      # loop all final states
                for fState in fStates:
                    bin_name, proc = generateName(year, channel[fState], recobin,
                                                fState, bin, physicalModel,
                                                observableBins, obsName)
                    print("    bin_name:", bin_name,"    proc:",proc)
                    comp_otherfid_parts.append(f"shapeSig_{SignalNames[physicalModel]}{proc}_{bin_name}")

    comp_otherfid = make_comp_list(comp_otherfid_parts)

    # --- background components ----------------------------------------------------
    def collect_bkg(expr_tpl_list):
        parts = []
        for year in years:
            #print("year", year)
            #print("year_cat", year_cat)
            #if year != year_cat:
                #continue
            if fstate != "4l":
                print("---------------------------------------")
                print("fstate", fstate)
                #print("fstate_cat", fstate_cat)
                print("---------------------------------------")
                #if fstate != fstate_cat:
                    #continue
                bin_name, _ = generateName(year, channel[fstate], recobin, fstate, 0, physicalModel, observableBins, obsName)
                for tpl in expr_tpl_list:
                    parts.append(tpl.format(bin_name=bin_name))
            else:
                for fState in fStates:
                    #if fState != fstate_cat:
                        #continue
                    bin_name, _ = generateName(year, channel[fState], recobin, fState, 0, physicalModel, observableBins, obsName)
                    print("----bin_name", bin_name)
                    for tpl in expr_tpl_list:
                        parts.append(tpl.format(bin_name=bin_name))
        print("parts", parts)
        return make_comp_list(parts)

    comp_out  = collect_bkg([f"shapeBkg_{OutNames[physicalModel]}_{{bin_name}}"])
    comp_fake = collect_bkg([f"shapeBkg_{CombNames[physicalModel]}_{{bin_name}}"])
    comp_zz   = collect_bkg(["shapeBkg_bkg_ggzz_{bin_name}","shapeBkg_bkg_qqzz_{bin_name}"])
    comp_zx   = collect_bkg(["shapeBkg_bkg_zjets_{bin_name}"])

    print("comp_out", comp_out)
    print("comp_fake", comp_fake)
    print("comp_zz", comp_zz)
    print("comp_zx", comp_zx)

    # ------------------------ plotting --------------------------------------------
    def overlay(line_color, line_style, *components, slice_label=None):
        """
        Plot `components` with optional category slice.

        If `slice_label` is None → draw the PDF projected over *all*
        categories present in `data` (works for the 4l combined case).

        If `slice_label` is a string → draw only that single category.
        """
        # common options
        opts = [
            RooFit.Components(make_comp_list(components)),
            RooFit.ProjWData(CMS_channel, data, True),
            RooFit.LineColor(line_color),
            RooFit.LineStyle(line_style)
        ]

        # add the slice only when requested
        if slice_label is not None:
            opts.insert(0, RooFit.Slice(CMS_channel, slice_label))

        sim.plotOn(mass, *opts)

    if fstate == "4l":          # ← combined 4l: draw everything together
        overlay(kGreen+2, 1, comp_zx, comp_zz, comp_fake, comp_otherfid, comp_out)
        overlay(kOrange-3, 2, comp_zx, comp_zz, comp_fake, comp_otherfid)
        overlay(kAzure-3, 1, comp_zx, comp_zz, comp_fake)
        overlay(kViolet,   1, comp_zx, comp_zz)
        overlay(kViolet+2, 1, comp_zx)
    else:                       # ← single final state: keep the slice
        # build cat_val once, e.g. "hzz_mass4l_105p0_160p0_cat4e_2022"
        cat_val = CMS_channel.getLabel()
        print("CAT LABEL", cat_val)
        overlay(kGreen+2, 1, comp_zx, comp_zz, comp_fake, comp_otherfid, comp_out)
                #slice_label=cat_val)
        overlay(kOrange-3, 2, comp_zx, comp_zz, comp_fake, comp_otherfid)
                #slice_label=cat_val)
        overlay(kAzure-3, 1, comp_zx, comp_zz, comp_fake)
                #slice_label=cat_val)
        overlay(kViolet,   1, comp_zx, comp_zz)
                #slice_label=cat_val)
        overlay(kViolet+2, 1, comp_zx)
                #slice_label=cat_val)

    # Finally the data points
    print("mass", mass)
    data.plotOn(mass)

    ##### ------------------------ Plot ------------------------ #####
    print("-------Plot----------")
    gStyle.SetOptStat(0)

    c = TCanvas("c","c",1000,800)
    c.cd()

    #dummy = TH1D("","",1,105.6,160.6)
    dummy = TH1D("","",1,opt.LOWER_BOUND,opt.UPPER_BOUND)
    dummy.SetBinContent(1,2)
    dummy.SetFillColor(0)
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.SetMarkerColor(0)
    dummy.GetYaxis().SetTitle("Events / (1.83 GeV)")
    dummy.GetXaxis().SetTitle("m_{"+fstate.replace("mu","#mu")+"} [GeV]")
    if (opt.UNBLIND):
        dummy.SetMaximum(max(1*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        # dummy.SetMaximum(max(0.8*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        
        '''if fstate == "4l":
            dummy.SetMaximum(42)
        elif fstate == "2e2mu":
            dummy.SetMaximum(25)
        elif fstate == "4mu":
            dummy.SetMaximum(18)
        elif fstate == "4e":
            dummy.SetMaximum(8)'''
        
    else:
        dummy.SetMaximum(max(1*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),1.0))
        if (obsName=="massZ2" and recobin==0): dummy.SetMaximum(max(3.0*max(n_trueH_asimov[fstate],n_trueH_modelfit[fstate]),3.5))

    dummy.SetMinimum(0.0)
    dummy.Draw()

    dummy_data = TH1D()
    dummy_data.SetMarkerColor(kBlack)
    dummy_data.SetMarkerStyle(20)
    dummy_fid = TH1D()
    if opt.BONLY:
        dummy_fid.SetLineColor(kRed)
    else:
        dummy_fid.SetLineColor(kOrange-3)
    dummy_fid.SetLineWidth(2)
    dummy_other = TH1D()
    dummy_other.SetLineColor(kOrange-3)
    dummy_other.SetLineWidth(2)
    dummy_other.SetLineStyle(2)
    dummy_out = TH1D()
    dummy_out.SetLineColor(kGreen+2)
    dummy_out.SetLineWidth(2)
    dummy_fake = TH1D()
    dummy_fake.SetLineColor(kAzure-3)
    dummy_fake.SetLineWidth(2)
    dummy_zz = TH1D()
    dummy_zz.SetLineColor(kViolet)
    dummy_zz.SetLineWidth(2)
    dummy_zx = TH1D()
    dummy_zx.SetLineColor(kViolet+2)
    dummy_zx.SetLineWidth(2)

    legend = TLegend(.20,.6,.93,.89)
    if (not opt.UNBLIND):
       legend.AddEntry(dummy_data,"Asimov Data","ep")
    else:
       legend.AddEntry(dummy_data,"Data","ep")

    if opt.BONLY:
        legend.AddEntry(dummy_fid,"N_{fid} (prefit)", "l")
    else:
        legend.AddEntry(dummy_fid,"N_{fid}", "l")
    legend.AddEntry(dummy_out, "N_{nonfid}", "l")
    legend.AddEntry(dummy_fake, "N_{nonres}", "l")
    legend.AddEntry(dummy_zz, "N_{ZZ}", "l")
    legend.AddEntry(dummy_zx, "N_{ZX}", "l")

    legend.SetShadowColor(0)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetNColumns(2)
    legend.SetTextSize(0.045)
    legend.Draw()

    mass.Draw("same")

    if (obsName=="pT4l"):
        label="p_{T}^{H}"
        unit="GeV"
    elif(obsName=="pT4lj"):
        label="p_{T}^{Hj}"
        unit="GeV"
    elif (obsName=="massZ2"):
        label = "m(Z_{2})"
        unit = "GeV"
    elif (obsName=="massZ1"):
        label = "m(Z_{1})"
        unit = "GeV"
    elif (obsName=="nJets" or obsName=="njets_reco_pt30_eta4p7"):
        label = "N(jets) |#eta|<4.7"
        unit = ""
    elif (obsName=="njets_pt30_eta4p7"):
        label = "N(jets) |#eta|<4.7"
        unit = ""
    elif (obsName=='pt_leadingjet_pt30_eta4p7'):
        label = "p_{T}(jet)"
        unit = "GeV"
    elif (obsName=='pt_leadingjet_pt30_eta2p5'):
        label = "p_{T}(jet) |#eta|<2.5"
        unit = "GeV"
    elif (obsName=='absdeltarapidity_hleadingjet_pt30_eta4p7'):
        label = "|y(H)-y(jet)|"
        unit = ""
    elif (obsName=='absdeltarapidity_hleadingjet_pt30_eta2p5'):
        label = "|y(H)-y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=='absrapidity_leadingjet_pt30_eta4p7'):
        label = "|y(jet)|"
        unit = ""
    elif (obsName=='absrapidity_leadingjet_pt30_eta2p5'):
        label = "|y(jet)| |#eta|<2.5"
        unit = ""
    elif (obsName=="rapidity4l"):
        label = "|y^{H}|"
        unit = ""
    elif (obsName=="costhetastar"):
        label = "|cos#theta*|"
        unit = ""
    elif (obsName=="costhetaZ1"):
        label = "|cos#theta_{1}|"
        unit = ""
    elif (obsName=="cosTheta2"):
        label = "|cos#theta_{2}|"
        unit = ""
    elif (obsName=="phi"):
        label = "|#Phi|"
        unit = ""
    elif (obsName=="Phi"):
        label = "|#Phi^{#star}|"
        unit = ""
    elif (obsName=="mass4l") or (obsName=='mass4l_zzfloating'):
        label = "inclusive"
        unit = ""
    else:
        label = obsName
        unit = ""

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
#    latex2.DrawLatex(0.87, 0.95,"35.9 fb^{-1} at #sqrt{s} = 13 TeV")
#    latex2.DrawLatex(0.87, 0.95,"41.4 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
#    latex2.DrawLatex(0.87, 0.95,"59.7 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
#    latex2.DrawLatex(0.87, 0.95,"136.0 fb^{-1} at #sqrt{s} = 13 TeV") # 2017
    latex2.DrawLatex(0.95, 0.95,lumi+" fb^{-1} (#sqrt{s} = 13.6 TeV)") # 2017
    latex2.SetTextSize(0.8*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.16, 0.94, "CMS")
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    # latex2.DrawLatex(0.30, 0.95, "Preliminary")
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.45*c.GetTopMargin())
    # if (obsName!='mass4l' and obsName!='mass4l_zzfloating'): latex2.DrawLatex(0.65,0.85, str(observableBins[recobin])+" "+unit+" < "+label+" < "+str(observableBins[recobin+1])+" "+unit)

    # Base plots directory from paths.py
    plots_base = path["plots_path"]

    # Create base directory if needed
    checkDir(plots_base)
    plots_base = os.path.join(plots_base, "SHAPES")

    checkDir(os.path.join(plots_base, obsName))

    suffix = ""

    if opt.PREFIT:
        suffix = "_prefit"
    if opt.BONLY:
        suffix = "_bonlyprefit"
    if opt.UNBLIND:
        suffix = "_sbpostfit"
    if opt.BONLY and opt.UNBLIND:
        suffix = "_bonlypostfit"


    if (not opt.UNBLIND):
        asimov_dir = os.path.join(plots_base, obsName, "asimov", "model")
        os.makedirs(asimov_dir, exist_ok=True)

        outname = f"asimovdata_{physicalModel}_{opt.YEAR}_{obsName}_{fstate}_recobin{recobin}{suffix}"

        c.SaveAs(os.path.join(asimov_dir, outname + ".pdf"))
        c.SaveAs(os.path.join(asimov_dir, outname + ".png"))
        c.SaveAs(os.path.join(asimov_dir, outname + ".root"))
    else:
        data_dir = os.path.join(plots_base, obsName, "data", "model")
        os.makedirs(data_dir, exist_ok=True)

        outname = f"data_unfoldwith_{modelName}_{physicalModel}_{opt.YEAR}_{obsName}_{fstate}_recobin{recobin}{suffix}"

        c.SaveAs(os.path.join(data_dir, outname + ".pdf"))
        c.SaveAs(os.path.join(data_dir, outname + ".png"))
        c.SaveAs(os.path.join(data_dir, outname + ".root"))




######### ---------------------------------------- #########
######### ----------------- Main ----------------- #########
######### ---------------------------------------- #########
print("-------Main----------")
modelName = opt.UNFOLD

if 'vs' in opt.OBSNAME:
    obsName_tmp = opt.OBSNAME.split(' vs ')
    obsName = obsName_tmp[0]+'_'+obsName_tmp[1]
    doubleDiff = True
else:
    obsName = opt.OBSNAME
    doubleDiff = False

# Load observableBins from eos_path
print("LOADING EOS PATH")
eos_inputs_path = path["eos_path"]

sys.path.append(eos_inputs_path)

if opt.INTER:
    module_name = 'inputs.inputs_sig_extrap_' + obsName + '_' + opt.YEAR
else:
    module_name = 'inputs.inputs_sig_' + obsName + '_' + opt.YEAR
print("MODULE NAME", module_name)
_temp = __import__(module_name, globals(), locals(), ['observableBins'])
observableBins = _temp.observableBins

sys.path.remove(eos_inputs_path)

if obsName.startswith("mass4l"):
    PhysicalModels = ['v3']
elif obsName == 'D0m' or obsName == 'Dcp' or obsName == 'D0hp' or obsName == 'Dint' or obsName == 'DL1' or obsName == 'DL1Zg' or obsName == 'costhetaZ1' or obsName == 'costhetaZ2'or obsName == 'costhetastar' or obsName == 'phi' or obsName == 'phistar' or obsName == 'massZ1' or obsName == 'massZ2':
    PhysicalModels = ['v3']
elif 'kL' in obsName:
    PhysicalModels = ['kLambda']
else:
    PhysicalModels = ['v3']


nBins = len(observableBins)
if not doubleDiff: nBins = nBins-1 #in case of 1D measurement the number of bins is -1 the length of the list of bin boundaries
# print(nBins)
fStates = ["4e","4mu","2e2mu","4l"]
for fState in fStates:
    print("fState", fState)
    for recobin in range(nBins):
        print("recobin", recobin)
        for physicalModel in PhysicalModels:
            plotAsimov_sim(opt.UNFOLD, physicalModel, obsName, fState, observableBins, recobin)