import os, sys
import numpy as np
import pandas as pd
import uproot as uproot
import ROOT
import argparse
from array import array
import random
import math
from tqdm import tqdm

from observables import observables
from binning import binning
from paths import path

np.seterr(invalid='ignore', divide='ignore')

np.set_printoptions(linewidth=np.inf)

print()
print('Welcome to RunBinning!')
print()

def parseOptions():
    global opt, args, runAllSteps

    parser = argparse.ArgumentParser(description="Run binning for fiducial cross-section measurement")

    parser.add_argument('--obsName', required=True, type=str, help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_argument('--obsBins', required=False, type=str, help='Bin boundaries for the differential measurement separated by "|", e.g. "|0|50|100|". Use empty string for default.')
    parser.add_argument('--year', required=True, type=str, help='Year (e.g., "2022", "2023preBPix", etc.)')
    parser.add_argument('--m4lLower', dest='m4lLower', type=int, default=110, help='Lower bound for m4l (default: 105)')
    parser.add_argument('--m4lUpper', dest='m4lUpper', type=int, default=130, help='Upper bound for m4l (default: 160)')
    parser.add_argument('--finstate', dest='finstate', type=str, default="4l", help='4mu, 4e, 2e2mu, 2mu2e, 4l')
    parser.add_argument('--do2024', dest='do2024', action='store_true', help='add 2024 lumi (109.08 invfb) to 2023postBPix')
    parser.add_argument('--PRINT', dest='PRINT', action='store_true', help='print output of each binning')
    parser.add_argument('--AUTO', dest='AUTO', action='store_true', help='do binning')
    parser.add_argument('--NBINS', required=False, type=int, help='number of bins to try')

    opt = parser.parse_args()
    return opt

opt = parseOptions()

# ---------------------------------------------------- FUNCTIONS TO GENERATE DATAFRAMES ----------------------------------------------------

# Weights for histogram
def weight(df, xsec, gen, lumi, type, additional = None):

    # overallEventWeight[0] = PUWeight[0]*tree.genHEPMCweight*tree.dataMCWeight

    if type == "sig":

        coeff = (lumi * 1000 * xsec) / gen

        #weight_reco = np.sign(df.genHEPMCweight) * df.PUWeight * df.dataMCWeight
        weight_reco = df.genHEPMCweight * df.PUWeight * df.dataMCWeight
        #weight_gen = np.sign(df.genHEPMCweight)
        #eight_gen = df.genHEPMCweight
        
        if additional == "ggH":
           weight_reco *= df.ggH_NNLOPS_weight

        weight_event = weight_reco * coeff
        #weight_total = weight_reco.sum()

        weight = weight_event # / weight_total

        df['weight'] = weight

    if type =="bkg":

        weight = (lumi * 1000 * df.overallEventWeight * df.dataMCWeight)/gen
        df['weight'] = weight

    return df

# Uproot to generate pandas
def prepareTrees(samples, year, treename):
    df = {}
    for sample in samples:

        if year == "2024" and sample == "ZZTo4l":
            year = "2023postBPix" # use 2023postBPix for 2024 ZZTo4l sample
        elif year == "2024" and sample.startswith("ggTo"):
            year = "2023postBPix" # use 2023postBPix for 2024 ggZZ samples
        elif year == "2024" and sample == "ZH125":
            year = "2023postBPix" # use 2023postBPix for 2024 ZH samples

        fname = path['eos_path_sig']+year+"_MC/"+sample+"/ZZ4lAnalysis_SKIMMED.root"
        df[sample] = uproot.open(fname)[treename]
    return df

# Calculate cross sections
def xsecs(samples, year, treename):

    xsec_ = {}
    if treename == treepass:
        df = prepareTrees(samples, year, treename)
        for sample in samples:
            total_weight = df[sample].arrays("overallEventWeight",library="np")["overallEventWeight"]
            puweight = df[sample].arrays("PUWeight", library="np")["PUWeight"]
            genweight = df[sample].arrays("genHEPMCweight", library="np")["genHEPMCweight"]
            if 'ZZTo' in sample:
                # TODO: Add EW KFactor once in the samples
                KFactor_QCD_qqZZ_M_Weight = df[sample].arrays("KFactor_QCD_qqZZ_M_Weight", library="np")["KFactor_QCD_qqZZ_M_Weight"]
                xsec = total_weight/(puweight*genweight*KFactor_QCD_qqZZ_M_Weight)
            elif 'ggTo' in sample:
                KFactor_QCD_ggZZ_Nominal_Weight = df[sample].arrays("KFactor_QCD_ggZZ_Nominal_Weight", library="np")["KFactor_QCD_ggZZ_Nominal_Weight"]
                xsec = total_weight/(puweight*genweight*KFactor_QCD_ggZZ_Nominal_Weight)
            else:
                xsec = total_weight/(puweight*genweight)

            xsec_[sample] = xsec[0]
    else:
        for sample in samples:
            xsec_[sample] = 1

    return xsec_

# Get the "number" of MC events to divide the weights
def generators(samples, year, df):
    gen = {}
    
    for sample in samples:

        if year == "2024" and sample == "ZZTo4l":
            year = "2023postBPix" # use 2023postBPix for 2024 ZZTo4l sample
        elif year == "2024" and sample.startswith("ggTo"):
            year = "2023postBPix" # use 2023postBPix for 2024 ggZZ samples
        elif year == "2024" and sample == "ZH125":
            year = "2023postBPix" # use 2023postBPix for 2024 ZH samples


        fname = path['eos_path_sig']+year+"_MC/"+sample+"/ZZ4lAnalysis_SKIMMED.root"
        gen[sample] = uproot.open(fname)["Counters"].values()[39] 
    return gen

# Define the final state
def add_fin_state(i, j):
    if abs(i) == 121 and abs(j) == 121:
        fin = '4e'
    elif abs(i) == 169 and abs(j) == 169:
        fin = '4mu'
    elif abs(i) == 121 and abs(j) == 169:
        fin = '2e2mu'
    elif abs(i) == 169 and abs(j) == 121:
        fin = '2mu2e'
    else: print('Problem with add_fin_state')
    return fin

def dataframes(samples, year, year_mc, treename, type):
    if year_mc == '2022':
        lumi = 7.9804
    elif year_mc == '2022EE':
        lumi = 26.6728
    elif year_mc == '2023preBPix':
        lumi = 17.794
    elif year_mc == '2023postBPix':
        lumi = 9.451
        if opt.do2024:
            lumi += 109.08
    elif year_mc == '2024':
        lumi = 109.08   

    d = {}
    df_all = prepareTrees(samples, year_mc, treename)
    gen_df = generators(samples, year_mc, df_all)
    xsec_df = xsecs(samples, year_mc, treename)

    for sample in samples:
        if sample == "ggH125":
            vars.append("ggH_NNLOPS_weight")
        gen = gen_df[sample]
        xsec = xsec_df[sample]
        df_np = df_all[sample].arrays(vars, library="np")
        df = pd.DataFrame({var: df_np[var] for var in vars})
        df['FinState'] = [add_fin_state(i, j) for i, j in zip(df.Z1Flav, df.Z2Flav)]
        if type == "sig":
            df['cuth4l_gen'] = [add_cuth4l_gen(i,j) for i,j in zip(df.GENlep_MomMomId,df.GENlep_Hindex)]
            df['cuth4l_reco'] = [add_cuth4l_reco(list(row[0]), list(row[1]), list(row[2]), list(row[3])) for row in zip(df['lep_Hindex'], df['lep_genindex'], df['GENlep_MomMomId'], df['GENlep_MomId'])]
        if sample == "ggH125":
            d[sample] = weight(df, xsec, gen, lumi, type, "ggH")
            vars.remove("ggH_NNLOPS_weight")
        else:
            d[sample] = weight(df, xsec, gen, lumi, type)

    return d

def skim_df(samples, year, year_mc, treename, type):

    df = dataframes(samples, year, year_mc, treename, type)
    df_skim = {}
    frames = []
    
    ggZZ_frames = []

    for sample in samples:
        if sample.startswith('ggTo'):
            ggZZ_frames.append(df[sample])
        elif sample == 'ZZTo4l':
            df_skim['qqZZ'] = df[sample]
            frames.append(df[sample])
        else:
            df_skim[sample] = df[sample]
            frames.append(df[sample])

    if ggZZ_frames:
        df_skim['ggZZ'] = pd.concat(ggZZ_frames)
        frames.append(df_skim['ggZZ'])

    if type == 'bkg':
        df_skim['ZX'] = do_ZX(df_skim, year, year_mc)

    print(f"{year} {type} skimmed df created")
    return df_skim

def add_cuth4l_gen(momMomId,Hindex):
    if (int(Hindex[0])==99) | (int(Hindex[1])==99) | (int(Hindex[2])==99) | (int(Hindex[3])==99):
        return False
    if int(momMomId[int(Hindex[0])])==25 and int(momMomId[int(Hindex[1])])==25 and int(momMomId[int(Hindex[2])])==25 and int(momMomId[int(Hindex[3])])==25:
        return True
    else:
        return False
#Checks if all four generator-level leptons are descendants of a Higgs boson (i.e., if it’s a true H→ZZ→4l decay).

def add_cuth4l_reco(Hindex,genIndex,momMomId,momId): #(Hindex, momMomId,momId):
    if (Hindex[0]==99) | (Hindex[1]==99) | (Hindex[2]==99) | (Hindex[3]==99) | (int(Hindex[0])==-1) | (int(Hindex[1])==-1) | (int(Hindex[2])==-1) | (int(Hindex[3])==-1):
        return False
    if ((genIndex[Hindex[0]]>-0.5)*momMomId[max(0,genIndex[Hindex[0]])]==25) and ((genIndex[Hindex[0]]>-0.5)*momId[max(0,genIndex[Hindex[0]])]==23) and ((genIndex[Hindex[1]]>-0.5)*momMomId[max(0,genIndex[Hindex[1]])]==25) and ((genIndex[Hindex[1]]>-0.5)*momId[max(0,genIndex[Hindex[1]])]==23) and ((genIndex[Hindex[2]]>-0.5)*momMomId[max(0,genIndex[Hindex[2]])]==25) and ((genIndex[Hindex[2]]>-0.5)*momId[max(0,genIndex[Hindex[2]])]==23) and ((genIndex[Hindex[3]]>-0.5)*momMomId[max(0,genIndex[Hindex[3]])]==25) and ((genIndex[Hindex[3]]>-0.5)*momId[max(0,genIndex[Hindex[3]])]==23):
        return True
    else:
        return False
# ------------------------------------------------- FUNCTIONS TO GENERATE DATAFRAME FOR ZX ---------------------------------------------------

def FindFinalState(z1_flav, z2_flav):
    if(z1_flav == -121):
        if(z2_flav == +121): return 0 # 4e
        if(z2_flav == +169): return 2 # 2e2mu
    if(z1_flav == -169):
        if(z2_flav == +121): return 3 # 2mu2e
        if(z2_flav == +169): return 1 # 4mu

def findFSZX(df):
    df['FinState'] = [FindFinalState(x,y) for x,y in zip(df['Z1Flav'], df['Z2Flav'])]
    return df

def openFR(year):
    
    if (year == "2022" or year == "2022EE" or year == "2023preBPix" or year == "2023postBPix"):
        fnameFR = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250303/FAKERATES/%s/FakeRates_SS_%s.root" % (year, year)
    elif (year == "2024"):
        fnameFR = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250303/FAKERATES/2023postBPix/FakeRates_SS_2023postBPix.root" # using 2023postBPix for 2024
    else:
        raise ValueError(f"ERROR: Unsupported year")

    if not os.path.exists(fnameFR):
        raise FileNotFoundError(f"Fake rate file not found: {fnameFR}")

    file = uproot.open(fnameFR)

    input_file_FR = ROOT.TFile(fnameFR)
    FR_mu_EB = input_file_FR.Get("FR_SS_muon_EB")
    FR_mu_EE = input_file_FR.Get("FR_SS_muon_EE")
    FR_e_EB  = input_file_FR.Get("FR_SS_electron_EB")
    FR_e_EE  = input_file_FR.Get("FR_SS_electron_EE")
    return FR_mu_EB, FR_mu_EE, FR_e_EB, FR_e_EE

def GetFakeRate(lep_Pt, lep_eta, lep_ID):
    if(lep_Pt >= 80):
        my_lep_Pt = 79
    else:
        my_lep_Pt = lep_Pt
    my_lep_ID = abs(lep_ID)
    if((my_lep_Pt > 5) & (my_lep_Pt <= 7)): bin = 0
    if((my_lep_Pt > 7) & (my_lep_Pt <= 10)): bin = 1
    if((my_lep_Pt > 10) & (my_lep_Pt <= 20)): bin = 2
    if((my_lep_Pt > 20) & (my_lep_Pt <= 30)): bin = 3
    if((my_lep_Pt > 30) & (my_lep_Pt <= 40)): bin = 4
    if((my_lep_Pt > 40) & (my_lep_Pt <= 50)): bin = 5
    if((my_lep_Pt > 50) & (my_lep_Pt <= 80)): bin = 6
    if(abs(my_lep_ID) == 11): bin = bin-1 # There is no [5, 7] bin in the electron fake rate # CHECK THIS - SPENCER
    if(my_lep_ID == 11):
        if(abs(lep_eta) < 1.479): return FR_e_EB.GetY()[bin]
        else: return FR_e_EE.GetY()[bin]
    if(my_lep_ID == 13):
        if(abs(lep_eta) < 1.2): return FR_mu_EB.GetY()[bin]
        else: return FR_mu_EE.GetY()[bin]

def comb(year):
    if year == "2022": # 2022 from HIG 24 13, 2023 from SPENCER
        cb_SS = np.array([
            1.239, # 4e
            1.093, # 4mu
            1.057, # 2e2mu
            1.254, # 2mu2e
        ])
    elif year == "2022EE":
        cb_SS = np.array([
            1.067, # 4e
            1.015, # 4mu
            1.049, # 2e2mu
            0.905, # 2mu2e
        ])
    elif year == "2023preBPix":
        cb_SS = np.array([
            1.116, # 4e
            1.036, # 4mu
            0.989, # 2e2mu
            1.141, # 2mu2e
        ])
    elif year == "2023postBPix":
        cb_SS = np.array([
            0.795, # 4e
            1.025, # 4mu
            1.074, # 2e2mu
            1.078, # 2mu2e
        ])
    elif year == "2024": # using 2023postBPix
        cb_SS = np.array([
            0.795, # 4e
            1.025, # 4mu
            1.074, # 2e2mu
            1.078, # 2mu2e
        ])
    return cb_SS

def ratio(year): # 2022 from HIG 24 13, 2023 from SPENCER
    if year == "2022":
        OS_SS = np.array([
            1.030,   # 4e
            1.165,  # 4mu
            0.966,   # 2e2mu
            1.041,  # 2mu2e
            ])
    elif year == "2022EE":
        OS_SS = np.array([
            0.990,   # 4e
            0.997,  # 4mu
            1.039,   # 2e2mu
            1.016,  # 2mu2e
            ])
    elif year == "2023preBPix":
        OS_SS = np.array([
            0.992,   # 4e
            1.024,  # 4mu
            1.102,   # 2e2mu
            1.024,  # 2mu2e
            ])
    elif year == "2023postBPix":
        OS_SS = np.array([
            1.006,   # 4e
            1.040,  # 4mu
            1.078,   # 2e2mu
            1.025,  # 2mu2e
            ])
    elif year == "2024": # using 2023postBPix
        OS_SS = np.array([
            1.006,   # 4e
            1.040,  # 4mu
            1.078,   # 2e2mu
            1.025,  # 2mu2e
            ])
    return OS_SS

def ZXYield(df, year, year_mc):
    cb_SS = comb(year_mc)
    OS_SS = ratio(year_mc)
    vec = df.to_numpy()
    Yield = np.zeros(len(vec), float)
    for i in range(len(vec)):
        finSt  = vec[i][len(vars)] 
        lepPt  = vec[i][5]
        lepEta = vec[i][4]
        lepID  = vec[i][3]
        Yield[i] = cb_SS[finSt] * OS_SS[finSt] * GetFakeRate(lepPt[2], lepEta[2], lepID[2]) * GetFakeRate(lepPt[3], lepEta[3], lepID[3])
    return Yield

def do_ZX(df_skim, year, year_mc):

    keyZX = 'CRZLLTree/candTree'

    #if (year=="2022"): data = '/eos/home-s/sellissp/HZZ/SAMPLES/062025/2022_Data/Data_eraCD_preEE_SKIMMED.root'
    #if (year=="2022EE"): data = '/eos/home-s/sellissp/HZZ/SAMPLES/062025/2022_Data/Data_eraEFG_postEE_SKIMMED.root'
    #if (year=="2023preBPix"): data = '/eos/home-s/sellissp/HZZ/SAMPLES/062025/2023_Data/Data_eraC_preBPix_SKIMMED.root'
    #if (year=="2023postBPix"): data = '/eos/home-s/sellissp/HZZ/SAMPLES/062025/2023_Data/Data_eraD_postBPix_SKIMMED.root'
    if (year=="2022"): data = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2022_Data/Data_eraCD_preEE_SKIMMED.root'
    if (year=="2022EE"): data = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2022_Data/Data_eraEFG_postEE_SKIMMED.root'
    if (year=="2023preBPix"): data = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2023_Data/Data_eraC_preBPix_SKIMMED.root'
    if (year=="2023postBPix"): data = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2023_Data/Data_eraD_postBPix_SKIMMED.root'
    if (year=="2024"): data = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/062025/2024_Data/2024_Data.root'

    ttreeZX = uproot.open(data)[keyZX]
    ttreeZX = ttreeZX.arrays(zx_vars, library="np")
    dfZX = pd.DataFrame(columns=zx_vars)
    for var in zx_vars:
        dfZX[var] = ttreeZX[var]
    dfZX['overallEventWeight'] = 1
    dfZX['dataMCWeight'] = 1
    dfZX = dfZX[dfZX.Z2Flav > 0] # keep just same-sign events
    dfZX = findFSZX(dfZX)
    dfZX['ZX_yield'] = ZXYield(dfZX, year, year_mc)

    return dfZX


# ---------------------------------------------------- FUNCTIONS FOR COMPUTING BINNING RESULTS ----------------------------------------------------

def do_AMS(sig_counts, bkg_counts, bkg_errors):

    qqZZ, ggZZ, ZX = bkg_counts['qqZZ'], bkg_counts['ggZZ'], bkg_counts['ZX']
    qqZZ_err, ggZZ_err, ZX_errors = bkg_errors['qqZZ'], bkg_errors['ggZZ'], bkg_errors['ZX']

    if len(sig_counts) != len(qqZZ) or len(sig_counts) != len(ggZZ) or len(sig_counts) != len(ZX):
        raise ValueError("sig and bkg counts must have the same length")

    AMS = []
    for i in range(len(sig_counts)):
        s = sig_counts[i]
        b = qqZZ[i] + ggZZ[i] + ZX[i]

        sigma_b = np.sqrt(qqZZ_err[i]**2 + ggZZ_err[i]**2 + ZX_errors[i]**2)

        with np.errstate(divide='ignore', invalid='ignore'):  # suppress div by zero warnings
            term1 = (s+b) * np.log( ((s+b)*(b+sigma_b**2))/(b**2+(s+b)*sigma_b**2) )
            term2 = (b**2/sigma_b**2) * np.log(1 + (s*sigma_b**2)/(b*(b+sigma_b**2)))

            AMS_val = np.sqrt(2*term1 - 2*term2)
            AMS.append(AMS_val)

    return AMS

def do_AMS_simple(sig_counts, bkg_counts, bkg_errors):

    s = sig_counts
    b = bkg_counts
    sigma_b = bkg_errors

    with np.errstate(divide='ignore', invalid='ignore'):  # suppress div by zero warnings
        term1 = (s+b) * np.log( ((s+b)*(b+sigma_b**2))/(b**2+(s+b)*sigma_b**2) )
        term2 = (b**2/sigma_b**2) * np.log(1 + (s*sigma_b**2)/(b*(b+sigma_b**2)))

        AMS_val = np.sqrt(2*term1 - 2*term2)

    return AMS_val

def get_bin_index(x, bins):
    for i in range(len(bins) - 1):
        if bins[i] <= x < bins[i + 1]:
            return i
    if x == bins[-1]:
        return len(bins) - 2
    return -1  # out of rangef

def do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff, doprint):

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins, dtype=float)
        hist = ROOT.TH2F("tmp_hist", "tmp_hist", nbins, bin_array, nbins, bin_array)
        hist_fidonly = ROOT.TH1F("tmp1_hist", "tmp1_hist", nbins, bin_array)
        hist_recoonly = ROOT.TH1F("tmp2_hist", "tmp2_hist", nbins, bin_array)

    else:
        region_strs = [r.strip() for r in obsBins.split('/')]

        all_xbins = []
        all_ybins = []
        nbins_x_list = []
        nbins_y_list = []

        for region in region_strs:
            x_str, y_str = [s.strip() for s in region.split('vs')]
            xbins = [float(x) for x in x_str.strip('|').split('|') if x]
            ybins = [float(y) for y in y_str.strip('|').split('|') if y]

            all_xbins.append(xbins)
            all_ybins.append(ybins)
            nbins_x_list.append(len(xbins) - 1)
            nbins_y_list.append(len(ybins) - 1)

        nbins = sum(nx * ny for nx, ny in zip(nbins_x_list, nbins_y_list))

        hist = ROOT.TH2F("tmp_hist", "tmp_hist", nbins, 0, nbins, nbins, 0, nbins)

    for year in years:
        for sig in dfs_sig[year].keys():

            df_sel = dfs_sig[year][sig][(dfs_sig[year][sig]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['ZZMass'] <= opt.m4lUpper) &
                                        (dfs_sig[year][sig]['GENmass4l'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['GENmass4l'] <= opt.m4lUpper)]

            df_fidonly = df_sel[(df_sel['passedFullSelection'] == 0) &
                                (df_sel['passedFiducial'] == 1)]

            df_recoonly = df_sel[(df_sel['passedFullSelection'] == 1) &
                                 (df_sel['passedFiducial'] == 0)]

            df_sel = df_sel[(df_sel['passedFullSelection'] == 1) &
                            (df_sel['passedFiducial'] == 1)]

            if not doubleDiff:
                for x, y, w in zip(df_sel[observables[obs_name]['obs_gen']], 
                                   df_sel[observables[obs_name]['obs_reco']], 
                                   df_sel['weight']):
                    hist.Fill(x, y, w)
                for x, w in zip(df_fidonly[observables[obs_name]['obs_gen']], 
                                df_fidonly['weight']):
                    hist_fidonly.Fill(x, w)
                for x, w in zip(df_recoonly[observables[obs_name]['obs_reco']], 
                                df_recoonly['weight']):
                    hist_recoonly.Fill(x, w)
            else:
                obs_name_1st = obs_name.split('vs')[0].strip()
                obs_name_2nd = obs_name.split('vs')[1].strip()

                for x1, x2, y1, y2, w in zip(
                    df_sel[observables[obs_name_1st]['obs_gen']],
                    df_sel[observables[obs_name_2nd]['obs_gen']],
                    df_sel[observables[obs_name_1st]['obs_reco']],
                    df_sel[observables[obs_name_2nd]['obs_reco']],
                    df_sel['weight']):

                    # Find gen bin and gen offset region
                    gen_flat_bin = -1
                    gen_bin_offset = 0
                    gen_found = False
                    for xbins, ybins, nbins_x, nbins_y in zip(all_xbins, all_ybins, nbins_x_list, nbins_y_list):
                        gen_bin1 = get_bin_index(x1, xbins)
                        gen_bin2 = get_bin_index(x2, ybins)
                        if -1 not in (gen_bin1, gen_bin2):
                            gen_flat_bin = gen_bin_offset + gen_bin1 * nbins_y + gen_bin2
                            gen_found = True
                            break
                        gen_bin_offset += nbins_x * nbins_y

                    reco_flat_bin = -1
                    reco_bin_offset = 0
                    reco_found = False
                    for xbins, ybins, nbins_x, nbins_y in zip(all_xbins, all_ybins, nbins_x_list, nbins_y_list):
                        reco_bin1 = get_bin_index(y1, xbins)
                        reco_bin2 = get_bin_index(y2, ybins)
                        if -1 not in (reco_bin1, reco_bin2):
                            reco_flat_bin = reco_bin_offset + reco_bin1 * nbins_y + reco_bin2
                            reco_found = True
                            break
                        reco_bin_offset += nbins_x * nbins_y

                    if not (gen_found and reco_found):
                        continue

                    hist.Fill(gen_flat_bin + 0.5, reco_flat_bin + 0.5, w)

    matrix = np.zeros((nbins, nbins), dtype=float)
    for i in range(0, nbins + 1):
        for j in range(0, nbins + 1):
            matrix[i - 1, j - 1] = max(0.0, hist.GetBinContent(i, j))

    #recoonly = []
    #for i in range(1, nbins + 1):
    #    recoonly.append(max(0.0, hist_recoonly.GetBinContent(i, j)))
    #recoonly_ = recoonly / np.sum(recoonly)

    matrix = np.flipud(matrix)

    # normalize the matrix by gen count
    column_sums = matrix.sum(axis=0)
    matrix = np.divide( matrix, column_sums )

    np.set_printoptions(formatter={'float_kind': lambda x: f"{x:.2f}"})

    mig_score = do_migration_score(matrix)
    try:
        con_number = np.linalg.cond(matrix)
    except np.linalg.LinAlgError:
        con_number = 99

    if doprint:
        print(f"migration matrix:")
        print(matrix)
        #print(f"{np.round(recoonly_, 2)} - nonfid = {np.sum(recoonly):.2f}")
        print()
        print(f"migration score: {mig_score}")
        print(f"condition number: {con_number}")

    return mig_score, con_number


def do_migration_score(matrix):

    matrix = np.flipud(matrix)

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    total = np.sum(np.abs(matrix))
    diag_sum = np.sum(np.abs(np.diag(matrix)))

    return diag_sum / total

def do_binning(dfs_sig, dfs_bkg, years, obs_name, obsBins, doubleDiff, doprint):
    
    if doprint: print()

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins)

        if doprint:
            print("bins: ", bin_array)

    if doubleDiff:
        obs_name_1st = obs_name.split('vs')[0].strip()
        obs_name_2nd = obs_name.split('vs')[1].strip()

        regions = [region.strip() for region in obsBins.split('/') if region.strip()]
        
        all_xbins = []
        all_ybins = []
        
        if doprint: print("bins:")

        n = 0
        for region in regions:
            n += 1
            if 'vs' not in region:
                raise ValueError(f"Missing 'vs' in region: {region}")
            
            x_str, y_str = [s.strip() for s in region.split('vs')]
            xbins = [float(x) for x in x_str.strip('|').split('|') if x]
            ybins = [float(y) for y in y_str.strip('|').split('|') if y]

            all_xbins.append(np.array(xbins, dtype=float))
            all_ybins.append(np.array(ybins, dtype=float))

            if doprint: print(f"bin {n}: ", end="")
            for i in range(len(xbins)-1):
                for j in range(len(ybins)-1):
                    if doprint: print(f"[{xbins[i]}, {xbins[i+1]}]x[{ybins[j]}, {ybins[j+1]}]", end="  ")
            if doprint: print()

        nbins = sum((len(xb)-1)*(len(yb)-1) for xb, yb in zip(all_xbins, all_ybins))
        
        if doprint:
            print(f"total bins: {nbins}")
            print()

    sum_sig_weighted_counts = {sig: np.zeros(nbins) for sig in dfs_sig[years[0]].keys()}
    total_sig_weighted_counts, total_sig_bin_error2 = np.zeros(nbins), np.zeros(nbins)
    sum_bkg_weighted_counts = {bkg: np.zeros(nbins) for bkg in dfs_bkg[years[0]].keys()}
    total_bkg_weighted_counts, total_bkg_bin_error2 = np.zeros(nbins), np.zeros(nbins)

    for year in years:
        for sig in dfs_sig[year].keys():

            df_sel = dfs_sig[year][sig][(dfs_sig[year][sig]['passedFullSelection'] == 1) &
                                        (dfs_sig[year][sig]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['ZZMass'] <= opt.m4lUpper)
                                        & (dfs_sig[year][sig]['cuth4l_reco'] == True)]

            if opt.finstate != "4l":
                df_sel = df_sel[(df_sel['FinState'] == opt.finstate)]

            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_sig_{sig}_{year}", "", nbins, bin_array) 
                if obs_name == 'rapidity4l': 
                    col = observables[obs_name]['obs_reco']
                    df_sel.loc[:, col] = df_sel[col].abs()
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_sel['weight']):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]
                if obs_name_1st == 'rapidity4l': df_obs_x = np.abs(df_obs_x)
                if obs_name_2nd == 'rapidity4l': df_obs_y = np.abs(df_obs_y)

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_sel['weight'].values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_sig_{sig}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)

                    for x, y, w in zip(obs_x_vals, obs_y_vals, weights):
                        if x >= xbins[0] and x < xbins[-1] and y >= ybins[0] and y < ybins[-1]:
                            hist.Fill(x, y, w)

                    region_counts = np.array([
                        hist.GetBinContent(i + 1, j + 1)
                        for i in range(nbins_x)
                        for j in range(nbins_y)
                    ])
                    weighted_counts_list.append(region_counts)
                    bin_errors_list.append(np.sqrt(region_counts))

                weighted_counts = np.concatenate(weighted_counts_list)
                bin_errors = np.concatenate(bin_errors_list)


            sum_sig_weighted_counts[sig] += weighted_counts
            total_sig_weighted_counts += weighted_counts
            total_sig_bin_error2 += bin_errors ** 2

    if doprint:
        for sample, counts in sum_sig_weighted_counts.items():
            print(f"{sample} weighted counts per bin: {counts}")
            #print(f"    total: {np.sum(counts)}")
        print(f"sig total: {np.sum(total_sig_weighted_counts)}")
        print(f"total sig weighted counts per bin: {total_sig_weighted_counts}")
        print(f"total sig errors per bin: {np.sqrt(total_sig_bin_error2)}")
        print()

    bkg_weighted_counts = {}
    bkg_binerrors = {}

    for year in years:
        for bkg in dfs_bkg[year].keys():

            df_sel = dfs_bkg[year][bkg][(dfs_bkg[year][bkg]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_bkg[year][bkg]['ZZMass'] <= opt.m4lUpper)]

            if opt.finstate != "4l":
                if bkg == 'ZX': 
                    if opt.finstate == "4e": df_sel = df_sel[(df_sel['FinState'] == 0)]
                    if opt.finstate == "4mu": df_sel = df_sel[(df_sel['FinState'] == 1)]
                    if opt.finstate == "2e2mu": df_sel = df_sel[(df_sel['FinState'] == 2)]
                    if opt.finstate == "2mu2e": df_sel = df_sel[(df_sel['FinState'] == 3)]
                else:
                    df_sel = df_sel[(df_sel['FinState'] == opt.finstate)]

            if bkg == 'ZX': 
                df_weights = df_sel['ZX_yield']
            else: 
                df_weights = df_sel['weight']
                
            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_bkg_{bkg}_{year}", "", nbins, bin_array)
                if obs_name == 'rapidity4l': 
                    col = observables[obs_name]['obs_reco']
                    df_sel.loc[:, col] = df_sel[col].abs()
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_weights):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]
                if obs_name_1st == 'rapidity4l': df_obs_x = np.abs(df_obs_x)
                if obs_name_2nd == 'rapidity4l': df_obs_y = np.abs(df_obs_y)

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_weights.values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_bkg_{bkg}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)

                    # Loop over events and fill if in range of current region
                    for x, y, w in zip(obs_x_vals, obs_y_vals, weights):
                        if x >= xbins[0] and x < xbins[-1] and y >= ybins[0] and y < ybins[-1]:
                            hist.Fill(x, y, w)

                    # Collect bin contents and errors from this region
                    region_counts = np.array([
                        hist.GetBinContent(i + 1, j + 1)
                        for i in range(nbins_x)
                        for j in range(nbins_y)
                    ])
                    weighted_counts_list.append(region_counts)
                    bin_errors_list.append(np.sqrt(region_counts))

                weighted_counts = np.concatenate(weighted_counts_list)
                bin_errors = np.concatenate(bin_errors_list)

            if bkg == 'ZX' and opt.year == 'Run3' and year == "2023postBPix" and opt.do2024 == True:
                weighted_counts *= (109.08+9.451)/9.451
                bin_errors = np.sqrt(weighted_counts)
                
            sum_bkg_weighted_counts[bkg] += weighted_counts
            total_bkg_weighted_counts += weighted_counts
            total_bkg_bin_error2 += bin_errors**2
            bkg_weighted_counts[bkg] = weighted_counts
            bkg_binerrors[bkg] = bin_errors
            

    if doprint:
        for bkg, counts in sum_bkg_weighted_counts.items():
            print(f"{bkg} weighted counts per bin: {counts}")
        for bkg, counts in sum_bkg_weighted_counts.items():
            print(f"{bkg} total: {np.sum(counts)}")
        print(f"total bkg weighted counts per bin: {total_bkg_weighted_counts}")
        print(f"total bkg errors per bin: {np.sqrt(total_bkg_bin_error2)}")
        print()

    with np.errstate(divide='ignore', invalid='ignore'):
        ssqrtb = np.divide(total_sig_weighted_counts, np.sqrt(total_bkg_weighted_counts), out=np.zeros_like(total_sig_weighted_counts), where=total_bkg_weighted_counts!=0)
    
    AMS = do_AMS(total_sig_weighted_counts, bkg_weighted_counts, bkg_binerrors)

    if doprint:
        print(f"AMS per bin: {AMS}")
        print()

    migration_score, condition_number = do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff, doprint)

    return AMS, migration_score, condition_number

# ------------------------------------------------- FUNCTIONS TO AUTOMATICALLY OPTIMIZE BINNING ---------------------------------------------------



def parse_bins(obs_bins):
    return [float(x) for x in obs_bins.strip('|').split('|') if x]

def bins_to_string(bins):
    return "|" + "|".join(str(x) for x in bins) + "|"

def flatten_dfs(dfs_sig, obs_name):
    vals = []
    wts  = []
    for year_dict in dfs_sig.values():
        for name, df in year_dict.items():
            df = df[(df['ZZMass'] >= opt.m4lLower) & (df['ZZMass'] <= opt.m4lUpper)]
            obs_val = df[observables[obs_name]['obs_reco']]
            if name == 'ZX':
                weight = df['ZX_yield']
            else:
                weight = df['weight']
            vals.append(obs_val)
            wts.append(weight)
    return pd.concat(vals, ignore_index=True), pd.concat(wts, ignore_index=True)

def get_init_binning_soversqrtb(obs_name, dfs_sig, dfs_bkg, step_size, targetams, targetsig, eps=1e-3):

    # flatten signal and background
    sig_values, sig_weights = flatten_dfs(dfs_sig, obs_name)
    bkg_values, bkg_weights = flatten_dfs(dfs_bkg, obs_name)

    sig_values = sig_values.to_numpy()
    sig_weights = sig_weights.to_numpy()
    bkg_values = bkg_values.to_numpy()
    bkg_weights = bkg_weights.to_numpy()

    # restrict to range [xmin, xmax]
    if obs_name == 'pTj1' or obs_name == 'pTj2':
        xmin = 30.0  # jets must have pT > 30 GeV
    elif obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        xmin = 0.0   # mjj, detajj, mHj, pTHj, pTHjj must be positive
    else:
        xmin = min(sig_values.min(), bkg_values.min())
    xmax = max(sig_values.max(), bkg_values.max())

    # scan from xmin to xmax
    scan_points = np.arange(xmin, xmax + step_size, step_size)

    edges = [xmin]   # always start at xmin

    S, B = 0.0, 0.0  # running sums

    for i in range(len(scan_points) - 1):
        lo, hi = scan_points[i], scan_points[i + 1]

        mask_sig = (sig_values >= lo) & (sig_values < hi)
        mask_bkg = (bkg_values >= lo) & (bkg_values < hi)

        S += sig_weights[mask_sig].sum()
        B += bkg_weights[mask_bkg].sum()

        current_ams = do_AMS_simple(S, B, np.sqrt(B)) if B > 0 else 0.0
        #print(f"scanning [{lo:.2f}, {hi:.2f}]: S = {S:.2f}, B = {B:.2f}, AMS = {current_ams:.2f}")

        if B > 0 and current_ams >= targetams and S>targetsig:
            edges.append(hi)   # cut bin here
            #print(f"--> Bin closed at {hi:.2f} (AMS reached {current_ams:.2f})")
            S, B = 0.0, 0.0    # reset accumulators

    # make sure last bin closes at xmax
    if abs(edges[-1] - xmax) > eps:
        edges.append(xmax)

    return edges


def perturb_bins(bins, n_bins, var_min, var_max, max_shift):
    """
    Perturb a random number of bin edges while keeping order and number of bins fixed.
    
    bins      : current list of bin edges
    var_min   : minimum allowed value
    var_max   : maximum allowed value
    n_bins    : desired number of bins
    max_shift : maximum fraction of bin width to shift an edge
    """
    new_bins = bins[:]
    n_edges = len(new_bins) - 2  # exclude first and last edges

    # Randomly decide how many edges to shift (at least 1)
    num_edges_to_shift = random.randint(1, n_edges)

    # Randomly pick edges to shift
    edges_to_shift = random.sample(range(1, len(new_bins) - 1), num_edges_to_shift)

    for idx in edges_to_shift:
        left = new_bins[idx - 1]
        right = new_bins[idx + 1]
        
        # fractional shift relative to local bin width
        shift = random.uniform(-max_shift, max_shift) * (right - new_bins[idx])
        new_edge = new_bins[idx] + shift

        # Clamp to maintain ascending order
        new_edge = max(left + 1e-6, min(right - 1e-6, new_edge))
        new_bins[idx] = new_edge

    #new_bins = [var_min] + list(new_bins) + [var_max]

    return new_bins


def get_init_binning(obs_name, dfs_sig, doubleDiff, nbins, xmin=None, xmax=None, eps=1e-3):
    values, weights = flatten_dfs(dfs_sig, obs_name)
    values = values.to_numpy()
    weights = weights.to_numpy()

    # restrict to range [xmin, xmax]
    if xmin is not None:
        mask_min = values >= xmin
    else:
        mask_min = np.ones_like(values, dtype=bool)

    if xmax is not None:
        mask_max = values <= xmax
    else:
        mask_max = np.ones_like(values, dtype=bool)

    mask = mask_min & mask_max
    values = values[mask]
    weights = weights[mask]

    # sort and compute weighted CDF
    sorter = np.argsort(values)
    values_sorted = values[sorter]
    weights_sorted = weights[sorter]

    cdf = np.cumsum(weights_sorted)
    cdf /= cdf[-1]  # normalize to 1

    # compute percentiles between xmin and xmax
    percentiles = np.linspace(0, 1, nbins)
    array = np.interp(percentiles, cdf, values_sorted)[1:-1]  # exclude first & last edges


    return array

def random_bins(obs_name, n_bins_min, n_bins_max, step_size, targetams, targetsig):

    sig_values, sig_weights = flatten_dfs(dfs_sig, obs_name)
    bkg_values, bkg_weights = flatten_dfs(dfs_bkg, obs_name)

    sig_values = sig_values.to_numpy()
    sig_weights = sig_weights.to_numpy()
    bkg_values = bkg_values.to_numpy()
    bkg_weights = bkg_weights.to_numpy()

    # restrict to range [xmin, xmax]
    xmin = min(sig_values.min(), bkg_values.min())
    xmax = max(sig_values.max(), bkg_values.max())

    # Randomly pick number of bins
    n_bins = random.randint(n_bins_min, n_bins_max)

    #array = get_init_binning(obs_name, dfs_sig, doubleDiff, n_bins, xmin, xmax)
    array = get_init_binning_soversqrtb(obs_name, dfs_sig, dfs_bkg, step_size, targetams, targetsig, eps=1e-3)
    
    array = perturb_bins(array, n_bins, xmin, xmax, max_shift=2)

    if obs_name == 'pTj1' or obs_name == 'pTj2' or obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        array = [-100] + array

    return array


def meets_significance(obs_name, sig, low=3, high=4.5):
    """
    Check significance conditions:
    - All bins except the last must be >= low
    - All bins must be <= high
    """
    cond1 = all(x >= low for x in sig[:-1])
    if obs_name == 'pTj1' or obs_name == 'pTj2' or obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        cond2 = all(x <= high for x in sig[1:])  # all but first bin
    else:
        cond2 = all(x <= high for x in sig)       # all bins
    return cond1 and cond2

'''
def optimize_binning(initial_bins_str, dfs_sig, dfs_bkg, years, obs_name, doubleDiff, steps=1000):

    solutions = []

    initial_bins = parse_bins(initial_bins_str)
    n_bins = len(initial_bins) - 1

    print("INITIAL BINNING")
    sig, mig, con = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(initial_bins), doubleDiff, True)
    ams_score = np.sqrt(np.sum(np.array(sig)**2))
    solutions.append(((ams_score, con, mig), bins_to_string(initial_bins), sig))

    if obs_name == 'pTj1' or obs_name == 'pTj2':
        var_min, var_max = 30, initial_bins[-1]
        bins = initial_bins[1:]
    elif obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        var_min, var_max = 0, initial_bins[-1]
        bins = initial_bins[1:]
    else:
        var_min, var_max = initial_bins[0], initial_bins[-1]
        bins = initial_bins

    for step in range(steps):

        # propose new bins
        if obs_name == 'pTj1' or obs_name == 'pTj2' or obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
            new_bins = perturb_bins(bins, var_min, var_max, n_bins-1)
            new_bins = [-100] + new_bins  # keep first bin fixed
        else:
            new_bins = perturb_bins(bins, var_min, var_max, n_bins)

        new_bins_str = bins_to_string(new_bins)
        sig, mig, con = do_binning(dfs_sig, dfs_bkg, years, obs_name, new_bins_str, doubleDiff, opt.PRINT)

        # reject invalid proposals
        if not is_valid_sig(sig) or not meets_significance(sig):
            continue
        else:
            ams_score = np.sqrt(np.sum(np.array(sig)**2))
            # accept if not dominated
            if not any(pareto_dominates(other[0], (ams_score, con, mig)) for other in solutions):
                bins = new_bins
                solutions.append(((ams_score, con, mig), new_bins_str, sig))

    # extract Pareto front
    front = pareto_front(solutions)

    # pick best by tie-breaker
    best_metrics, best_bins, best_sig = max(
        front, key=lambda x: (x[0][2]) / (x[0][1] + 1e-6)
    )

    best_ams, best_cond, best_mig = best_metrics
    return best_bins, best_sig, best_cond, best_mig, front
'''
# -----------------------------------------------------------------------------------------                                                                                                                                                                                                                                                 
# ------------------------------- MAIN ----------------------------------------------------                                                                                                                                                                                                                                                 
# -----------------------------------------------------------------------------------------

sigs = ['ggH125', 'VBFH125', 'WminusH125', 'WplusH125', 'ZH125', 'ttH125']
bkgs = ['ZZTo4l', 'ggTo2e2mu_Contin_MCFM701','ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']

#eos_path_FR = path['eos_path_FR']
eos_path = path['eos_path']

treepass = 'ZZTree/candTree'
treefail = 'ZZTree/candTree_failed'

obs_name = opt.obsName

if opt.obsBins is not None:
    obs_bins = opt.obsBins
else:
    opt.AUTO = True


doubleDiff = 'vs' in obs_name

vars = ['overallEventWeight', 'dataMCWeight', 'ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt']

# TO MAKE CSV #
#obsvars = ['pT4l', 'rapidity4l', 'massZ1', 'massZ2', 'pTj1', 'pTj2', 'Nj', 'mjj', 'absdetajj', 'dphijj', 'pTHj', 'pTHjj', 'mHj', 'TBjmax', 'TCjmax']
#for obs in obsvars:
#    vars.append(observables[obs]['obs_gen'])
# TO MAKE CSV #

if not doubleDiff: 
    vars.append(observables[obs_name]['obs_reco'])
if doubleDiff:
    obs_name_1st = obs_name.split('vs')[0].strip()
    obs_name_2nd = obs_name.split('vs')[1].strip()
    vars.append(observables[obs_name_1st]['obs_reco'])
    vars.append(observables[obs_name_2nd]['obs_reco'])

zx_vars = vars.copy()
for var in ['overallEventWeight', 'dataMCWeight']:
    zx_vars.remove(var)

if (opt.year == '2022'):
    years_MC = ['2022']
    years = ["2022"]
if (opt.year == '2022EE'):
    years_MC = ['2022EE']
    years = ["2022EE"]
if (opt.year == '2023preBPix'):
    years_MC = ['2023preBPix']
    years = ["2023preBPix"]
if (opt.year == '2023postBPix'):
    years_MC = ['2023postBPix']
    years = ["2023postBPix"]
if (opt.year == '2024'):
    years_MC = ['2024']
    years = ["2024"]

if (opt.year == '2022full'):
    years_MC = ['2022', '2022EE']
    years = ["2022", "2022EE"]
if (opt.year == '2023full'):
    years_MC = ['2023preBPix', '2023postBPix']
    years = ["2023preBPix", "2023postBPix"]

if (opt.year == 'Run3'):
    years_MC = ['2022', '2022EE', '2023preBPix', '2023postBPix', '2024']
    years = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]


# --------------------------- RUN --------------------------------------------


dfs_sig = {}
dfs_bkg = {}

#save_dir = "/eos/home-s/sellissp/HZZ/SAMPLES/CSV/"
#os.makedirs(save_dir, exist_ok=True)
                                                       
for year, year_mc in zip(years, years_MC):
    FR_mu_EB, FR_mu_EE, FR_e_EB, FR_e_EE = openFR(year_mc)
    df_bkg = skim_df(bkgs, year, year_mc, treepass, 'bkg')
    dfs_bkg[year] = df_bkg
    #for key, df in dfs_bkg[year].items(): # TO MAKE CSV #
    #    df.to_csv(os.path.join(save_dir, f"bkg_{year}_{key}.csv"), index=False) # TO MAKE CSV #

for var in ['GENmass4l', 'genHEPMCweight', 'PUWeight', 'passedFiducial', 'passedFullSelection', 'GENlep_Hindex', 'lep_Hindex','lep_genindex','GENlep_MomMomId','GENlep_MomId', 'GENlep_id', 'EventNumber', 'GENZ_DaughtersId', 'GENZ_MomId']:
    vars.append(var)

if not doubleDiff: 
    vars.append(observables[obs_name]['obs_gen'])
if doubleDiff:
    vars.append(observables[obs_name_1st]['obs_gen'])
    vars.append(observables[obs_name_2nd]['obs_gen'])

for year, year_mc in zip(years, years_MC):
    df_sig_pass = skim_df(sigs, year, year_mc, treepass, 'sig')
    #df_sig_fail = skim_df(sigs, year, year_mc, treefail, 'sig')
    #df_sig = pd.concat([df_sig_pass, df_sig_fail], ignore_index=True)
    dfs_sig[year] = df_sig_pass
    #for key, df in dfs_sig[year].items(): # TO MAKE CSV #
    #    df.to_csv(os.path.join(save_dir, f"sig_{year}_{key}.csv"), index=False) # TO MAKE CSV #

'''
# TO USE CSV #
for year in years:
    dfs_bkg[year] = {}
    for filename in os.listdir(save_dir):
        if filename.startswith(f"bkg_{year}_") and filename.endswith(".csv"):
            key = filename[len(f"bkg_{year}_"):-4]
            filepath = os.path.join(save_dir, filename)
            dfs_bkg[year][key] = pd.read_csv(filepath)

for year in years:
    dfs_sig[year] = {}
    for filename in os.listdir(save_dir):
        if filename.startswith(f"sig_{year}_") and filename.endswith(".csv"):
            key = filename[len(f"sig_{year}_"):-4]
            filepath = os.path.join(save_dir, filename)
            dfs_sig[year][key] = pd.read_csv(filepath)
# TO USE CSV #
'''


if opt.AUTO:

    print(opt.obsName, " - AUTO BINNING")

    best_bins = None

    '''
    if obs_name in ("pTj1", "pTj2"):
        xmin, xmax = 30, 2000
        step_size = 1
    elif obs_name in ("mjj", "absdetajj", "pTHj", "pTHjj", "mHj"):
        xmin, xmax = 0, 10000
        step_size = 1
    elif obs_name in ("rapidity4l"):
        xmin, xmax = 0, 2.5
        step_size = 0.01
    elif obs_name in ("massZ1"):
        xmin, xmax = 40, 120
        step_size = 1
    elif obs_name in ("massZ2"):
        xmin, xmax = 12, 65
        step_size = 1
    else:
        step_size = 0.1


    # PERCENTILES
    idx = 0
    binnings = []

    for n_bins in range(3, 20+1):

        init_binning = get_init_binning(obs_name, dfs_sig, doubleDiff, n_bins+1, xmin, xmax, eps=1e-3)

        if obs_name in ("pTj1", "pTj2", "mjj", "absdetajj", "pTHj", "pTHjj", "mHj"):
            print("NBINS = ", n_bins+1)
            init_binning = [-100] + [xmin] + init_binning.tolist() + [xmax]
        else:
            print("NBINS = ", n_bins)
            init_binning = [xmin] + init_binning.tolist() + [xmax]

        binnings.append(init_binning)
        sig, mig, con = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(init_binning), doubleDiff, opt.PRINT)

        if obs_name in ("pTj1", "pTj2"):
            sig = np.array(sig[1:], dtype=float)  # ignore first bin
        else:
            sig = np.array(sig, dtype=float)

        if np.sum(sig) > 3.5*n_bins:
            idx += 1
            continue

        else:

            if obs_name in ("pTj1", "pTj2", "mjj", "absdetajj", "pTHj", "pTHjj", "mHj"):
                print("IDEAL NUMBER OF BINS = ", n_bins)
                best_bins, sig, con, mig, front = optimize_binning(bins_to_string(binnings[idx-1]), dfs_sig, dfs_bkg, years, obs_name, doubleDiff, steps=100)
            else:
                print("IDEAL NUMBER OF BINS = ", n_bins-1)
                best_bins, sig, con, mig, front = optimize_binning(bins_to_string(binnings[idx-1]), dfs_sig, dfs_bkg, years, obs_name, doubleDiff, steps=100)

            #try:
            #best_bins, score = optimize_binning(dfs_sig, dfs_bkg, years, obs_name, doubleDiff, xmin, xmax, bins_to_string(init_binning))

            #except Exception as e:
            #    print(f"n_bins={n_bins-1}: optimization exception: {e}")
            

            break
        
    
    # AMS
    init_binning = get_init_binning_soversqrtb(obs_name, dfs_sig, dfs_bkg, step_size, eps=0)

    if obs_name in ("pTj1", "pTj2", "mjj", "absdetajj", "pTHj", "pTHjj", "mHj"):
            print("IDEAL NUMBER OF BINS = ", len(parse_bins(bins_to_string(init_binning))))
            init_binning = [-100] + init_binning.tolist()
    else:
        print("IDEAL NUMBER OF BINS = ", len(parse_bins(bins_to_string(init_binning)))-1)
        init_binning = init_binning.tolist()


    best_bins, sig, con, mig, front = optimize_binning(bins_to_string(init_binning), dfs_sig, dfs_bkg, years, obs_name, doubleDiff, steps=5)
    '''

    n_bins_min = 99
    n_bins_max = 99

    step_size = 0.1

    n_trials = 1000

    if obs_name in ("pTj1", "pTj2", "mjj", "absdetajj", "pTHj", "pTHjj", "mHj"):
        target_ams = 3.5
        target_sig = 25
    else:
        target_ams = 3.5
        target_sig = 20

    all_binnings = []
    with tqdm(total=n_trials, desc="testing binnings (0 valid)") as pbar:
        for i in range(n_trials):
            array = random_bins(obs_name, n_bins_min, n_bins_max, step_size, target_ams, target_sig)
            sig, mig, con = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(array), doubleDiff, opt.PRINT)
            if meets_significance(obs_name, sig, low=3, high=4.5):
                all_binnings.append((array, mig/con))
                pbar.set_description(f"testing binnings ({len(all_binnings)} valid)")

            pbar.update(1)

    if len(all_binnings) == 0:
        print("No valid binning found, try more iterations")
    else:
        max_len = max(len(a[0]) for a in all_binnings) 
        max_nbins_binnings = [a for a in all_binnings if len(a[0]) == max_len]
        best_bins = max(max_nbins_binnings, key=lambda x: x[1])[0]

        print("\n" + "-"*100)
        print("OPTIMAL BINNING:", best_bins)
        do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(best_bins), doubleDiff, True)
        print("-"*100)

else:
    do_binning(dfs_sig, dfs_bkg, years, obs_name, obs_bins, doubleDiff, True)