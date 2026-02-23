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
from math import erf, sqrt


from observables import observables
from binning import binning
from paths import path

ROOT.TH1.AddDirectory(False) 

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
    parser.add_argument('--PRINT', dest='PRINT', action='store_true', help='print output of each binning')
    parser.add_argument('--AUTO', dest='AUTO', action='store_true', help='do binning')
    parser.add_argument('--doRun3', dest='doRun3', action='store_true', help='add 22/23 lumi (62 invfb) to 2024')

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
                KFactor_QCD_qqZZ_M_Weight = df[sample].arrays("KFactor_QCD_qqZZ_M_weight", library="np")["KFactor_QCD_qqZZ_M_weight"]
                xsec = total_weight/(puweight*genweight*KFactor_QCD_qqZZ_M_Weight)
            elif 'ggTo' in sample:
                KFactor_QCD_ggZZ_Nominal_Weight = df[sample].arrays("KFactor_QCD_ggZZ_Nominal_weight", library="np")["KFactor_QCD_ggZZ_Nominal_weight"]
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
        lumi = 18.06265
    elif year_mc == '2023postBPix':
        lumi = 9.693
    elif year_mc == '2024':
        lumi = 108.822
        if opt.doRun3:
            lumi += 7.9804+26.6728+18.06265+9.693

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
    
    if (year == "2022" or year == "2022EE" or year == "2023preBPix" or year == "2023postBPix" or year == "2024"):
        fnameFR = path['eos_path_FR']+"FAKERATES/%s/FakeRates_SS_%s.root" % (year, year)
    else:
        raise ValueError(f"ERROR: Unsupported year")

    if not os.path.exists(fnameFR):
        raise FileNotFoundError(f"Fake rate file not found: {fnameFR}")

    f = ROOT.TFile.Open(fnameFR)
    FR_mu_EB = f.Get("FR_SS_muon_EB").Clone()
    FR_mu_EE = f.Get("FR_SS_muon_EE").Clone()
    FR_e_EB  = f.Get("FR_SS_electron_EB").Clone()
    FR_e_EE  = f.Get("FR_SS_electron_EE").Clone()

    f.Close()
    del f
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
    elif year == "2024": 
        cb_SS = np.array([
            0.787, # 4e
            0.960, # 4mu
            0.958, # 2e2mu
            0.749, # 2mu2e
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
    elif year == "2024": 
        OS_SS = np.array([
            0.997,   # 4e
            1.051,  # 4mu
            1.051,   # 2e2mu
            1.024,  # 2mu2e
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

    if opt.doRun3:
        Yield *= 171/109
    return Yield

def do_ZX(df_skim, year, year_mc):

    keyZX = 'CRZLLTree/candTree'

    PATH = fname = path['eos_path_sig']

    if (year=="2022"): data = PATH + '2022_Data/Data_eraCD_preEE_SKIMMED.root'
    if (year=="2022EE"): data = PATH + '2022_Data/Data_eraEFG_postEE_SKIMMED.root'
    if (year=="2023preBPix"): data = PATH + '2023_Data/Data_eraC_preBPix_SKIMMED.root'
    if (year=="2023postBPix"): data = PATH + '2023_Data/Data_eraD_postBPix_SKIMMED.root'
    if (year=="2024"): data = PATH + '2024_Data/ZZ4lAnalysis_SKIMMED.root'

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


# ----------------------------------------------------------- HELPER FUNCTIONS  -------------------------------------------------------------------

def get_bin_index(x, bins):
    for i in range(len(bins) - 1):
        if bins[i] <= x < bins[i + 1]:
            return i
    if x == bins[-1]:
        return len(bins) - 2
    return -1  # out of range

def bins_to_array(bins):
    if isinstance(bins, str):
        return [float(x) for x in bins.split('|') if x]
    else:
        return bins

def bins_to_string(bins):
    if isinstance(bins, list):
        return "|" + "|".join(str(x) for x in bins) + "|"
    else:
        return bins

def bins_to_array_2d(bins_2d):

    if isinstance(bins_2d, str):

        rects = []
        rect_strs = bins_2d.split(" / ")

        for r in rect_strs:
            x_str, y_str = r.split(" vs ")

            # Extract values between pipes, ignoring empty tokens
            x_vals = [float(v) for v in x_str.split("|") if v]
            y_vals = [float(v) for v in y_str.split("|") if v]

            if len(x_vals) != 2 or len(y_vals) != 2:
                raise ValueError(f"Bad 2D bin format in rectangle: {r}")

            rects.append((x_vals[0], x_vals[1], y_vals[0], y_vals[1]))

        return rects

    else:
        return bins_2d


def bins_to_string_2d(rects):

    if isinstance(rects, (list, tuple)):

        rect_strs = []
        for (x1, x2, y1, y2) in rects:
            x_str = f"|{x1}|{x2}|"
            y_str = f"|{y1}|{y2}|"
            rect_strs.append(f"{x_str} vs {y_str}")

        return " / ".join(rect_strs)

    else:
        return rects

def flatten_dfs(dfs_sig, obs_name, doubleDiff_var=None, doubleDiff_constraint=None):

    vals = []
    wts  = []

    for year_dict in dfs_sig.values():
        for name, df in year_dict.items():

            df = df[(df['ZZMass'] >= opt.m4lLower) & (df['ZZMass'] <= opt.m4lUpper)]
                # & (df['passedFullSelection'] == 1) & (df['cuth4l_reco'] == True)]

            if doubleDiff_var is not None:
                df = df[(df[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                        (df[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])
                ]

            obs_val = df[observables[obs_name]['obs_reco']]

            if obs_name == 'rapidity4l':
                obs_val = abs(obs_val)

            if name == 'ZX':
                weight = df['ZX_yield']
            else:
                weight = df['weight']

            vals.append(obs_val)
            wts.append(weight)

    return pd.concat(vals, ignore_index=True), pd.concat(wts, ignore_index=True)

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

def do_migration_score(matrix):

    matrix = np.flipud(matrix)

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    total = np.sum(np.abs(matrix))
    diag_sum = np.sum(np.abs(np.diag(matrix)))

    return diag_sum / total

def meets_sig_and_counts(obs_name, sig, events, lowsig, lowevents):

    cond1 = all(x >= 0.9*lowsig for x in sig[:-1])
    cond2 = all(x >= lowevents for x in events[:-1])

    return cond1 and cond2

# ---------------------------------------------------- FUNCTIONS FOR COMPUTING BINNING RESULTS ----------------------------------------------------

def do_binning(dfs_sig, dfs_bkg, years, obs_name, obsBins, doubleDiff, doprint, doubleDiff_var=None, doubleDiff_constraint=None):
    
    if doprint: print()

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins)

        if doprint:
            print("bins:")
            print(bin_array)
            print(f"total bins: {nbins}")

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

            if doubleDiff_constraint is not None: 
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])]

            if opt.finstate != "4l":
                df_sel = df_sel[(df_sel['FinState'] == opt.finstate)]

            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_sig_{sig}_{year}", "", nbins, bin_array) 
                hist.SetDirectory(0)
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_sel['weight']):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)
                del hist

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_sel['weight'].values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_sig_{sig}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)
                    hist.SetDirectory(0)

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

                    del hist

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

            if doubleDiff_constraint is not None: 
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])]

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
                hist.SetDirectory(0)
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_weights):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

                del hist

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_weights.values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_bkg_{bkg}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)
                    hist.SetDirectory(0)

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

                    del hist

                weighted_counts = np.concatenate(weighted_counts_list)
                bin_errors = np.concatenate(bin_errors_list)

                
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

    migration_score, condition_number, matrix = do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff, doprint, doubleDiff_var=None, doubleDiff_constraint=None)

    return total_sig_weighted_counts, AMS, migration_score, condition_number, matrix

def do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff, doprint, doubleDiff_var=None, doubleDiff_constraint=None):

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins, dtype=float)
        hist = ROOT.TH2F("tmp_hist", "tmp_hist", nbins, bin_array, nbins, bin_array)
        hist.SetDirectory(0)

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
        hist.SetDirectory(0)

    for year in years:
        for sig in dfs_sig[year].keys():

            df_sel = dfs_sig[year][sig][(dfs_sig[year][sig]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['ZZMass'] <= opt.m4lUpper) &
                                        (dfs_sig[year][sig]['GENmass4l'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['GENmass4l'] <= opt.m4lUpper) &
                                        (dfs_sig[year][sig]['passedFullSelection'] == 1) &
                                        (dfs_sig[year][sig]['passedFiducial'] == 1)
                                        ]
            
            if doubleDiff_constraint is not None: 
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1]) &
                                (df_sel[observables[doubleDiff_var]['obs_gen']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_gen']] < doubleDiff_constraint[1])]

            if not doubleDiff:
                for x, y, w in zip(df_sel[observables[obs_name]['obs_gen']], 
                                   df_sel[observables[obs_name]['obs_reco']], 
                                   df_sel['weight']):
                    hist.Fill(x, y, w)

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

    del hist

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
        np.set_printoptions(
        threshold=np.inf,
        suppress=True,
        linewidth=200,
        formatter={'float_kind': '{:.2f}'.format}
        )
        print(f"migration matrix:")
        print(matrix)
        #print(f"{np.round(recoonly_, 2)} - nonfid = {np.sum(recoonly):.2f}")
        print()
        print(f"migration score: {mig_score}")
        print(f"condition number: {con_number}")

    return mig_score, con_number, matrix

# ------------------------------------------------- FUNCTIONS TO AUTOMATICALLY OPTIMIZE 1D BINNING ---------------------------------------------------

def get_init_binning_AMS(obs_name, dfs_sig, dfs_bkg, step_size, targetams, targetsig, doubleDiff, doubleDiff_var=None, doubleDiff_constraint=None, eps=1e-3):

    sig_values, sig_weights = flatten_dfs(dfs_sig, obs_name, doubleDiff_var, doubleDiff_constraint)
    bkg_values, bkg_weights = flatten_dfs(dfs_bkg, obs_name, doubleDiff_var, doubleDiff_constraint)

    sig_values = sig_values.to_numpy()
    sig_weights = sig_weights.to_numpy()
    bkg_values = bkg_values.to_numpy()
    bkg_weights = bkg_weights.to_numpy()

    # restrict range [xmin, xmax]
    if obs_name == 'pTj1' or obs_name == 'pTj2':
        xmin = 30.0  # jets must have pT > 30 GeV
        xmax = max(sig_values.max(), bkg_values.max())
    elif obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        xmin = 0.0 # mjj, detajj, mHj, pTHj, pTHjj must be positive
        xmax = max(sig_values.max(), bkg_values.max())
    elif obs_name == 'TBjmax' or obs_name == 'TCjmax':
        xmin = 0.0 
        xmax = max(sig_values.max(), bkg_values.max())
    elif obs_name == 'massZ1':
        xmin = 40
        xmax = 120
    elif obs_name == 'massZ2':
        xmin = 12
        xmax = 65
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
            S, B = 0.0, 0.0    # reset accumulators

    # make sure last bin closes at xmax
    if abs(edges[-1] - xmax) > eps:
        edges.append(xmax)

    if obs_name == 'pTj1' or obs_name == 'pTj2' or obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj' or obs_name == 'TBjmax' or obs_name == 'TCjmax':
        if edges[0] != -100:
            edges = [-100] + edges

    edges = bins_to_string(edges)

    return edges

def bifurcate_gaussian(mu, sigma1, sigma2):

    if (opt.PRINT): print("MU: ", mu, "S1: ", sigma1, "S2: ", sigma2)

    if np.random.rand() < sigma1 / (sigma1 + sigma2):
        return mu - abs(np.random.randn()) * sigma1
    else:
        return mu + abs(np.random.randn()) * sigma2

def perturb_bins(bins, var_min, var_max):

    new_bins = bins[:]
    n_edges = len(new_bins) - 2  # exclude first and last edges (min/max)

    num_edges_to_shift = random.randint(1, n_edges)
    edges_to_shift = random.sample(range(1, len(new_bins) - 1), num_edges_to_shift)

    for idx in edges_to_shift:
        left = new_bins[idx - 1] + 1e-2
        right = new_bins[idx + 1] - 1e-2

        sigmaL = (new_bins[idx] - left) / 2
        sigmaR = (right - new_bins[idx]) / 2
        mu = new_bins[idx]

        new_edge = bifurcate_gaussian(mu, sigmaL, sigmaR)

        new_edge = max(left + 1e-6, min(right - 1e-6, new_edge))
        new_edge = max(var_min, min(var_max, new_edge))

        new_bins[idx] = new_edge

    return new_bins

def randomize_bins(obs_name, array_init):

    bins = bins_to_array(array_init)

    sig_values, sig_weights = flatten_dfs(dfs_sig, obs_name)
    bkg_values, bkg_weights = flatten_dfs(dfs_bkg, obs_name)

    sig_values = sig_values.to_numpy()
    sig_weights = sig_weights.to_numpy()
    bkg_values = bkg_values.to_numpy()
    bkg_weights = bkg_weights.to_numpy()

    # restrict to range [xmin, xmax]
    xmin = min(sig_values.min(), bkg_values.min())
    xmax = max(sig_values.max(), bkg_values.max())

    if obs_name == 'pTj1' or obs_name == 'pTj2':
        xmin = 30.0  # jets must have pT > 30 GeV
    elif obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj':
        xmin = 0.0 # mjj, detajj, mHj, pTHj, pTHjj must be positive
    elif obs_name == 'TBjmax' or obs_name == 'TCjmax':
        xmin = 0.0
    elif obs_name == 'massZ1':
        xmin = 40
        xmax = 120
    elif obs_name == 'massZ2':
        xmin = 12
        xmax = 65

    array = perturb_bins(bins, xmin, xmax)

    if obs_name == 'pTj1' or obs_name == 'pTj2' or obs_name == 'mjj' or obs_name == 'absdetajj' or obs_name == 'mHj' or obs_name == 'pTHj' or obs_name == 'pTHjj' or obs_name == 'TBjmax' or obs_name == 'TCjmax':
        if array[0] != -100:
            array = [-100] + array

    return array

def check_matrix(obs_name, array, matrix, eff, res, doubleDiff, doubleDiff_var=None, doubleDiff_constraint=None, threshold_base=0.0):

    array = bins_to_array(array)

    matrix = np.flipud(matrix)  # keep your convention

    new_array = [array[0]]
    i = 0
    merged = False

    while i < len(matrix):

        delta = array[i+1] - array[i]

        if res[i] == 1e-8:
            P_expected = 0
        else:
            P_expected = eff * erf(delta / (2 * sqrt(2) * res[i]))

        if (opt.PRINT): print(f"Bin [{array[i]:.2f}, {array[i+1]:.2f}]: eff = {eff:.2f} res = {res[i]:.2f} delta = {delta:.2f} => P_expected = {P_expected:.2f}")

        # impose a floor
        P_target = max(P_expected, threshold_base)

        if matrix[i, i] > P_target:
            new_array.append(array[i+1])
            i += 1
        else:
            # merge this bin with the next by skipping appending upper edge
            if i + 1 < len(array) - 1:
                array = array[:i+1] + array[i+2:]  # merge with next bin
                merged = True
                break  # stop loop and recompute matrix
            else:
                # last bin, just keep it
                new_array.append(array[i+1])
                i += 1

    # if any merge happened, recompute matrix and re-run recursively
    if merged:
        events, sig, mig, con, new_matrix = do_binning(
            dfs_sig, dfs_bkg, years, obs_name, bins_to_string(array), doubleDiff, opt.PRINT, doubleDiff_var, doubleDiff_constraint
        )
        res = get_res(obs_name, array, dfs_sig, doubleDiff_var, doubleDiff_constraint)
        return check_matrix(obs_name, array, new_matrix, eff, res, doubleDiff, doubleDiff_var, doubleDiff_constraint, threshold_base)

    return bins_to_string(new_array)

def get_res(obs_name, obs_bins, dfs_sig, doubleDiff_var = None, doubleDiff_constraint = None):

    obs_bins = bins_to_array(obs_bins)  # make sure it's float array
    if (opt.PRINT): ("Bins:", obs_bins)

    obs_gen, obs_reco = [], []

    # flatten over years
    for year_dict in dfs_sig.values():
        for name, df in year_dict.items():
            df_sel = df[
                (df['ZZMass'] >= opt.m4lLower) & (df['ZZMass'] <= opt.m4lUpper) &
                (df['GENmass4l'] >= opt.m4lLower) & (df['GENmass4l'] <= opt.m4lUpper) &
                (df['passedFullSelection'] == 1) & (df['cuth4l_reco'] == True)
            ]

            if doubleDiff_constraint is not None:
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])]

            obs_gen_val = df_sel[observables[obs_name]['obs_gen']]
            obs_reco_val = df_sel[observables[obs_name]['obs_reco']]

            obs_gen.append(obs_gen_val)
            obs_reco.append(obs_reco_val)

    obs_gen = pd.concat(obs_gen, ignore_index=True).astype(float)
    obs_reco = pd.concat(obs_reco, ignore_index=True).astype(float)

    res = []

    # jet-related variables
    jet_vars = {"pTj1", "pTj2", "mjj", "absdetajj", "mHj", "pTHj", "pTHjj", "TBjmax", "TCjmax"}

    for i in range(len(obs_bins)-1):
        # if first bin and jet variable, assign tiny value
        if i == 0 and obs_name in jet_vars:
            res.append(1e-8)
            continue

        mask = (obs_reco >= obs_bins[i]) & (obs_reco < obs_bins[i+1])

        if obs_name in jet_vars:
            mask &= (obs_reco != -99) & (obs_gen != -99)

        x = obs_gen[mask] - obs_reco[mask]

        #---------------------------------------
        # REMOVE OUTLIERS
        #---------------------------------------
        # drop NaN/inf first
        x = x[np.isfinite(x)]

        # ---- outlier mask (sigma clip) ----
        if x.size >= 5:  # avoid doing this on tiny samples
            mu = np.mean(x)
            sig = np.std(x, ddof=1)
            if sig > 0:
                keep = np.abs(x - mu) < 1.0 * sig 
                x = x[keep]
        #---------------------------------------
        #---------------------------------------

        #rms = np.sqrt(np.nanmean(x**2)) if len(x) > 0 else np.nan
        rms = np.sqrt(np.nanmean((x-np.mean(x))**2)) if len(x) > 0 else np.nan


        res.append(rms)

    return res


def doublediff_quantiles(obs_name, values, N):

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if obs_name in {"pTj1","pTj2"}:#,"mjj","absdetajj","mHj","pTHj","pTHjj","TBjmax","TCjmax"}:
        vq = values[values >= 30]
        if vq.size == 0:
            return []
        qs = np.linspace(0, 1, N + 1)
        edges = np.quantile(vq, qs, method="lower")
        edges[0] = 30
        edges[-1] = vq.max()
    elif obs_name in {"mjj","absdetajj","mHj","pTHj","pTHjj","TBjmax","TCjmax"}:
        vq = values[values >= 0]
        if vq.size == 0:
            return []
        qs = np.linspace(0, 1, N + 1)
        edges = np.quantile(vq, qs, method="lower")
        edges[0] = 0
        edges[-1] = vq.max()
    else:
        qs = np.linspace(0, 1, N + 1)
        edges = np.quantile(values, qs, method="lower")
        edges[0] = values.min()
        edges[-1] = values.max()

    edges = np.maximum.accumulate(edges)

    return list(zip(edges[:-1], edges[1:]))

def make2dbinning(obs_name, array, bins, swap_xy):

    rects = []

    if obs_name == "pTj1 vs pTj2":
        rects.append((-100.0, 30.0, -100.0, 30.0))
        rects.append((30.0, 10000.0, -100.0, 30.0)) # 1j
    elif obs_name == "pT4l vs pTHj":
        rects.append((0.0, 10000.0, -100.0, 0.0)) # flip
    elif obs_name == "TCjmax vs pT4l":
        rects.append((-100.0, 0.0, 0.0, 10000.0))

    for xdef, (y1, y2) in zip(array, bins):

        if len(xdef) > 0 and not isinstance(xdef[0], (list, tuple)):
            x_pairs = zip(xdef[:-1], xdef[1:])
        else:
            x_pairs = xdef

        for x1, x2 in x_pairs:

            if obs_name == "pTj1 vs pTj2":
                # force everything else to start at 30
                x1c = max(x1, 30.0)
                y1c = max(y1, 30.0)

                # skip regions fully below 30
                if x2 <= 30.0 or y2 <= 30.0:
                    continue
                if x2 <= x1c or y2 <= y1c:
                    continue

                rect = (y1c, y2, x1c, x2) if swap_xy else (x1c, x2, y1c, y2)

            if obs_name == "TCjmax vs pT4l" or obs_name == "pT4l vs pTHj":
                # force everything else to start at 30
                x1c = max(x1, 0.0)
                y1c = max(y1, 0.0)

                # skip regions fully below 30
                if x2 <= 0.0 or y2 <= 0.0:
                    continue
                if x2 <= x1c or y2 <= y1c:
                    continue

                rect = (y1c, y2, x1c, x2) if swap_xy else (x1c, x2, y1c, y2)
            else:
                rect = (y1, y2, x1, x2) if swap_xy else (x1, x2, y1, y2)

            rects.append(rect)

    return rects





def do_binning(dfs_sig, dfs_bkg, years, obs_name, obsBins, doubleDiff, doprint, doubleDiff_var=None, doubleDiff_constraint=None):
    
    if doprint: print()

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins)

        if doprint:
            print("bins:")
            print(bin_array)
            print(f"total bins: {nbins}")

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

            if doubleDiff_constraint is not None: 
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])]

            if opt.finstate != "4l":
                df_sel = df_sel[(df_sel['FinState'] == opt.finstate)]

            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_sig_{sig}_{year}", "", nbins, bin_array) 
                hist.SetDirectory(0)
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_sel['weight']):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)
                del hist

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_sel['weight'].values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_sig_{sig}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)
                    hist.SetDirectory(0)

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

                    del hist

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

            if doubleDiff_constraint is not None: 
                df_sel = df_sel[(df_sel[observables[doubleDiff_var]['obs_reco']] >= doubleDiff_constraint[0]) &
                                (df_sel[observables[doubleDiff_var]['obs_reco']] < doubleDiff_constraint[1])]

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
                hist.SetDirectory(0)
                for x, w in zip(df_sel[observables[obs_name]['obs_reco']], df_weights):
                    hist.Fill(x, w)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

                del hist

            if doubleDiff:
                weighted_counts_list = []
                bin_errors_list = []

                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]

                obs_x_vals = df_obs_x.values
                obs_y_vals = df_obs_y.values
                weights = df_weights.values

                for idx, (xbins, ybins) in enumerate(zip(all_xbins, all_ybins)):
                    nbins_x = len(xbins) - 1
                    nbins_y = len(ybins) - 1

                    hist = ROOT.TH2F(f"hist_bkg_{bkg}_{year}_region{idx}", "", nbins_x, xbins, nbins_y, ybins)
                    hist.SetDirectory(0)

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

                    del hist

                weighted_counts = np.concatenate(weighted_counts_list)
                bin_errors = np.concatenate(bin_errors_list)

                
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

    migration_score, condition_number, matrix = do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff, doprint, doubleDiff_var=None, doubleDiff_constraint=None)

    return total_sig_weighted_counts, AMS, migration_score, condition_number, matrix

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
    df_bkg = skim_df(bkgs, year, year_mc, treepass, 'bkg').copy()

    if obs_name == 'rapidity4l':
        for key in df_bkg.keys():
            df_bkg[key]['ZZy'] = df_bkg[key]['ZZy'].abs()

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
    df_sig_pass = skim_df(sigs, year, year_mc, treepass, 'sig').copy()

    if obs_name == 'rapidity4l':
        for key in df_sig_pass.keys():
            df_sig_pass[key]['ZZy'] = df_sig_pass[key]['ZZy'].abs()
            df_sig_pass[key]['GENrapidity4l'] = df_sig_pass[key]['GENrapidity4l'].abs()

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

    n_trials = 1

    eff_4lep = 0.95**4/4 + 0.88**4/4 + 0.95**2*0.88**2/2
    eff_jet = 0.8

    step_size = 0.01

    target_sig_count = 0.0

    if obs_name in ("pTj1", "pTHj", "mHj","TBjmax","TCjmax"):  
        eff = eff_4lep * eff_jet
        target_ams = 3
    elif obs_name in ("pTj2", "mjj", "absdetajj", "pTHjj"): 
        eff = eff_4lep * eff_jet**2
        target_ams = 2
    elif obs_name in ("massZ1", "massZ2"):  
        eff = eff_4lep
        target_ams = 3.75
    else: 
        eff = eff_4lep
        target_ams = 3

    if doubleDiff:
        if obs_name_1st in ("pTj1", "pTHj", "mHj","TBjmax","TCjmax"):  
            eff1 = eff_4lep * eff_jet
            target_ams1 = 3
        elif obs_name_1st in ("pTj2", "mjj", "absdetajj", "pTHjj"): 
            eff1 = eff_4lep * eff_jet**2
            target_ams1 = 2
        elif obs_name_1st in ("massZ1", "massZ2"):  
            eff1 = eff_4lep
            target_ams1 = 3.75
        else: 
            eff1 = eff_4lep
            target_ams1 = 3

        if obs_name_2nd in ("pTj1", "pTHj", "mHj","TBjmax","TCjmax"):  
            eff2 = eff_4lep * eff_jet
            target_ams2 = 3
        elif obs_name_2nd in ("pTj2", "mjj", "absdetajj", "pTHjj"): 
            eff2 = eff_4lep * eff_jet**2
            target_ams2 = 2
        elif obs_name_2nd in ("massZ1", "massZ2"):
            eff2 = eff_4lep
            target_ams2 = 3.75
        else: 
            eff2 = eff_4lep
            target_ams2 = 3

        if obs_name == "rapidity4l vs pT4l":
            eff = eff_4lep
            target_ams = 3
        elif obs_name == "massZ1 vs massZ2":
            eff = eff_4lep
            target_ams = 3.75
        elif obs_name == "pTj1 vs pTj2":
            eff = eff_4lep * eff_jet**2
            target_ams = 2.5
        elif obs_name == "pT4l vs pTHj":
            eff = eff_4lep * eff_jet
            target_ams = 3
        elif obs_name == "Nj vs pT4l":
            eff = eff_4lep * eff_jet
            target_ams = 3
        elif obs_name == "TCjmax vs pT4l":
            eff = eff_4lep * eff_jet
            target_ams = 2.75
        


    all_binnings = []

    if doubleDiff:

        if obs_name == "Nj vs pT4l":

            # bin var 2 fix var 1
            bins_var2 = []
            var1_edges = [(0,0.99),(1,1.99), (2,100)]
            for constraint in var1_edges:
                print(f"---------------------------------------BINNING ({obs_name_2nd}) W/ {obs_name_1st} FIXED {constraint} ---------------------------------------")
                array_init = get_init_binning_AMS(obs_name_2nd, dfs_sig, dfs_bkg, step_size, target_ams, target_sig_count, False, obs_name_1st, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array_init, False, opt.PRINT, obs_name_1st, constraint)

                print(f"----------------------------------------AFTER MERGING ({obs_name_2nd}) W/ {obs_name_1st} FIXED {constraint} ----------------------------------------")
                res = get_res(obs_name_2nd, array_init, dfs_sig, obs_name_1st, constraint)
                array2 = check_matrix(obs_name_2nd, array_init, matrix, eff, res, False, obs_name_1st, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array2, False, opt.PRINT, obs_name_1st, constraint)

                bins_var2.append(bins_to_array(array2))

            binning2 = make2dbinning(obs_name, bins_var2, var1_edges, True)

            print()
            print(bins_to_string_2d(binning2))
            print()

            events2, sig2, mig2, con2, matrix2 = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning2), True, True)

            all_binnings.append((bins_to_string_2d(binning2), mig2/con2))
            print(f"BINNING {obs_name_2nd} W/ {obs_name_1st} FIXED")

            new_array2 = bins_to_array(bins_var2).copy()

            with tqdm(total=n_trials, desc=f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )----------------------------------") as pbar:
                
                for i in range(n_trials):
                    
                    n = math.isqrt(len(events))
                    numregions_to_shift = random.randint(1, n)

                    for i in range(0, numregions_to_shift):

                        j = random.randint(0, n-1)

                        array = bins_var2[j]
                        constraint = var1_edges[j]

                        array = randomize_bins(obs_name_2nd, array)
                        events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, bins_to_string(array), False, opt.PRINT, obs_name_1st, constraint)
                        new_array2[j] = array
                        
                    binning2 = make2dbinning(obs_name,new_array2, var1_edges, False)
                    events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning2), doubleDiff, opt.PRINT)

                    if meets_sig_and_counts(obs_name, sig, events, target_ams, target_sig_count):

                        all_binnings.append((binning2, mig/con))

                    pbar.set_description(f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )---------------------------------")
                    pbar.update(1)

        else:

            print(f"---------------------------------------INITIAL BINNING ({obs_name_1st})---------------------------------------")
            array_init_1st = get_init_binning_AMS(obs_name_1st, dfs_sig, dfs_bkg, step_size, target_ams1, target_sig_count, False)
            events_1st, sig_1st, mig_1st, con_1st, matrix_1st = do_binning(dfs_sig, dfs_bkg, years, obs_name_1st, array_init_1st, False, opt.PRINT)

            print(f"----------------------------------------AFTER MERGING ({obs_name_1st})----------------------------------------")
            res_1st = get_res(obs_name_1st, array_init_1st, dfs_sig)

            array_1st = check_matrix(obs_name_1st, array_init_1st, matrix_1st, eff1, res_1st, False)
            events_1st, sig_1st, mig_1st, con_1st, matrix_1st = do_binning(dfs_sig, dfs_bkg, years, obs_name_1st, array_1st, False, opt.PRINT)

            print(f"---------------------------------------INITIAL BINNING ({obs_name_2nd})---------------------------------------")
            array_init_2nd = get_init_binning_AMS(obs_name_2nd, dfs_sig, dfs_bkg, step_size, target_ams2, target_sig_count, False)
            events_2nd, sig_2nd, mig_2nd, con_2nd, matrix_2nd = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array_init_2nd, False, opt.PRINT)

            print(f"----------------------------------------AFTER MERGING ({obs_name_2nd})----------------------------------------")
            res_2nd = get_res(obs_name_2nd, array_init_2nd, dfs_sig)
            array_2nd = check_matrix(obs_name_2nd, array_init_2nd, matrix_2nd, eff2, res_2nd, False)
            events_2nd, sig_2nd, mig_2nd, con_2nd, matrix_2nd = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array_2nd, False, opt.PRINT)

            obs1_vals, obs1_weights = flatten_dfs(dfs_sig, obs_name_1st)
            var1_edges = doublediff_quantiles(obs_name_1st, obs1_vals, math.isqrt(len(events_1st)))

            obs2_vals, obs2_weights = flatten_dfs(dfs_sig, obs_name_2nd)
            var2_edges = doublediff_quantiles(obs_name_2nd, obs2_vals, math.isqrt(len(events_2nd)))

            # bin var 1 fix var 2
            bins_var1 = []
            for constraint in var2_edges:
                print(f"---------------------------------------BINNING ({obs_name_1st}) W/ {obs_name_2nd} FIXED {constraint} ---------------------------------------")
                array_init = get_init_binning_AMS(obs_name_1st, dfs_sig, dfs_bkg, step_size, target_ams, target_sig_count, False, obs_name_2nd, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_1st, array_init, False, opt.PRINT, obs_name_2nd, constraint)

                print(f"----------------------------------------AFTER MERGING ({obs_name_1st}) W/ {obs_name_2nd} FIXED {constraint} ----------------------------------------")
                res = get_res(obs_name_1st, array_init, dfs_sig, obs_name_2nd, constraint)
                array1 = check_matrix(obs_name_1st, array_init, matrix, eff, res, False, obs_name_2nd, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_1st, array1, False, opt.PRINT, obs_name_2nd, constraint)

                bins_var1.append(bins_to_array(array1))

            print(f"{obs_name_1st}: BINS {bins_var1}")
            print(f"{obs_name_2nd}: BINS {var2_edges}")

            binning1 = make2dbinning(obs_name, bins_var1, var2_edges, False)

            print()
            print(bins_to_string_2d(binning1))
            print()

            events1, sig1, mig1, con1, matrix1 = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning1), True, True)


            # bin var 2 fix var 1
            bins_var2 = []
            for constraint in var1_edges:
                print(f"---------------------------------------BINNING ({obs_name_2nd}) W/ {obs_name_1st} FIXED {constraint} ---------------------------------------")
                array_init = get_init_binning_AMS(obs_name_2nd, dfs_sig, dfs_bkg, step_size, target_ams, target_sig_count, False, obs_name_1st, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array_init, False, opt.PRINT, obs_name_1st, constraint)

                print(f"----------------------------------------AFTER MERGING ({obs_name_2nd}) W/ {obs_name_1st} FIXED {constraint} ----------------------------------------")
                res = get_res(obs_name_2nd, array_init, dfs_sig, obs_name_1st, constraint)
                array2 = check_matrix(obs_name_2nd, array_init, matrix, eff, res, False, obs_name_1st, constraint)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, array2, False, opt.PRINT, obs_name_1st, constraint)

                bins_var2.append(bins_to_array(array2))

            print(f"{obs_name_2nd}: BINS {bins_var2}")
            print(f"{obs_name_1st}: BINS {var1_edges}")

            binning2 = make2dbinning(obs_name, bins_var2, var1_edges, True)

            print()
            print(bins_to_string_2d(binning2))
            print()

            events2, sig2, mig2, con2, matrix2 = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning2), True, True)


            # toys


            if ( mig1/con1 ) > ( mig2/con2 ):

                all_binnings.append((bins_to_string_2d(binning1), mig1/con1))
                print(f"BINNING {obs_name_1st} W/ {obs_name_2nd} FIXED")

                new_array1 = bins_to_array(bins_var1).copy()

                with tqdm(total=n_trials, desc=f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )----------------------------------") as pbar:
                    
                    for i in range(n_trials):

                        n = math.isqrt(len(events_2nd))
                        numregions_to_shift = random.randint(1, n)

                        for i in range(0, numregions_to_shift):

                            j = random.randint(0, n-1)

                            array = bins_var1[j]
                            constraint = var2_edges[j]

                            array = randomize_bins(obs_name_1st, array)
                            events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_1st, bins_to_string(array), False, opt.PRINT, obs_name_2nd, constraint)

                            new_array1[j] = array

                        binning1 = make2dbinning(obs_name, new_array1, var2_edges, False)
                        events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning1), doubleDiff, opt.PRINT)

                        if meets_sig_and_counts(obs_name, sig, events, target_ams, target_sig_count):

                            all_binnings.append((binning1, mig/con))

                        pbar.set_description(f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )---------------------------------")
                        pbar.update(1)

            else:
                all_binnings.append((bins_to_string_2d(binning2), mig2/con2))
                print(f"BINNING {obs_name_2nd} W/ {obs_name_1st} FIXED")

                new_array2 = bins_to_array(bins_var2).copy()

                with tqdm(total=n_trials, desc=f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )----------------------------------") as pbar:
                    
                    for i in range(n_trials):
                        
                        n = math.isqrt(len(events_1st))
                        numregions_to_shift = random.randint(1, n)

                        for i in range(0, numregions_to_shift):

                            j = random.randint(0, n-1)

                            array = bins_var2[j]
                            constraint = var1_edges[j]

                            array = randomize_bins(obs_name_2nd, array)
                            events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name_2nd, bins_to_string(array), False, opt.PRINT, obs_name_1st, constraint)
                            new_array2[j] = array
                            
                        binning2 = make2dbinning(obs_name,new_array2, var1_edges, False)
                        events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(binning2), doubleDiff, opt.PRINT)

                        if meets_sig_and_counts(obs_name, sig, events, target_ams, target_sig_count):

                            all_binnings.append((binning2, mig/con))

                        pbar.set_description(f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )---------------------------------")
                        pbar.update(1)

                
    else:
        
        print("---------------------------------------INITIAL BINNING---------------------------------------")
        array_init = get_init_binning_AMS(obs_name, dfs_sig, dfs_bkg, step_size, target_ams, target_sig_count, doubleDiff)
        events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, array_init, doubleDiff, True)

        print("----------------------------------------AFTER MERGING----------------------------------------")
        res = get_res(obs_name, array_init, dfs_sig)
        array = check_matrix(obs_name, array_init, matrix, eff, res, doubleDiff)
        events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, array, doubleDiff, True)

        if meets_sig_and_counts(obs_name, sig, events, target_ams, target_sig_count) and con < 99 and not math.isnan(mig):
            all_binnings.append((array, mig/con))

        with tqdm(total=n_trials, desc=f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )----------------------------------") as pbar:
            for i in range(n_trials):

                array = randomize_bins(obs_name, array)
                events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(array), doubleDiff, opt.PRINT)
                if meets_sig_and_counts(obs_name, sig, events, target_ams, target_sig_count):
                    all_binnings.append((array, mig/con))
                    pbar.set_description(f"--------------------------------TESTING BINNINGS ( {len(all_binnings)} VALID )---------------------------------")
                pbar.update(1)

    if len(all_binnings) == 0:
        print("NO VALID BINNING FOUND")
    else:
        max_len = max(len(a[0]) for a in all_binnings) 
        max_nbins_binnings = [a for a in all_binnings if len(a[0]) == max_len]
        best_bins = max(max_nbins_binnings, key=lambda x: x[1])[0]

    print("\n" + "-"*100)
    if best_bins is not None:
        if doubleDiff:
            print("OPTIMAL BINNING:", bins_to_string_2d(best_bins))
            do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string_2d(best_bins), doubleDiff, True)
        else:
            print("OPTIMAL BINNING:", bins_to_string(best_bins))
            do_binning(dfs_sig, dfs_bkg, years, obs_name, bins_to_string(best_bins), doubleDiff, True)
    print("-"*100)

else:
    do_binning(dfs_sig, dfs_bkg, years, obs_name, obs_bins, doubleDiff, True)



'''
# ------------------------------------------------- FUNCTIONS TO AUTOMATICALLY OPTIMIZE 2D BINNING ---------------------------------------------------


def get_init_binning_2D_AMS(xmin, xmax, ymin, ymax,
                     obs1_values, obs2_values, obs1_weights,
                     bkg1_values, bkg2_values, bkg1_weights,
                     targetams, targetsig, step_size_x, step_size_y,
                     rng_seed=None):

    if rng_seed is not None:
        random.seed(rng_seed)
        np.random.seed(rng_seed)

    x_vals = np.arange(xmin, xmax + step_size_x, step_size_x)
    y_vals = np.arange(ymin, ymax + step_size_y, step_size_y)

    Nx = len(x_vals) - 1
    Ny = len(y_vals) - 1

    if Nx <= 0 or Ny <= 0:
        return []

    occupied = np.zeros((Nx, Ny), dtype=bool)
    rectangles = []

    def cell_to_coords(x0_idx, x1_idx, y0_idx, y1_idx):
        lo_x = x_vals[x0_idx]
        hi_x = x_vals[x1_idx + 1] 
        lo_y = y_vals[y0_idx]
        hi_y = y_vals[y1_idx + 1]
        return lo_x, hi_x, lo_y, hi_y

    def build_integral(x, y, w):
        
        ix = np.floor((x - xmin) / step_size_x).astype(np.int32)
        iy = np.floor((y - ymin) / step_size_y).astype(np.int32)

        # keep only in-range
        m = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Ny) & np.isfinite(w)
        ix = ix[m]; iy = iy[m]; w = w[m].astype(np.float64)

        # weighted grid
        grid = np.zeros((Nx, Ny), dtype=np.float64)
        np.add.at(grid, (ix, iy), w)

        # integral image with padding row/col at 0
        integ = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
        integ[1:, 1:] = grid.cumsum(axis=0).cumsum(axis=1)
        return integ

    def rect_sum(integ, x0i, x1i, y0i, y1i):
        # inclusive cell indices -> half-open corners in integral coords
        x0 = x0i
        x1 = x1i + 1
        y0 = y0i
        y1 = y1i + 1
        return float(integ[x1, y1] - integ[x0, y1] - integ[x1, y0] + integ[x0, y0])

    def mark_occupied(x0i, x1i, y0i, y1i):
        occupied[x0i:x1i+1, y0i:y1i+1] = True

    def any_free():
        return not occupied.all()

    def find_free_adjacent():
        # search for any free cell that has a neighbor occupied (4-neighborhood)
        for xi in range(Nx):
            for yi in range(Ny):
                if occupied[xi, yi]:
                    continue
                neighbors = []
                if xi > 0:
                    neighbors.append(occupied[xi-1, yi])
                if xi < Nx-1:
                    neighbors.append(occupied[xi+1, yi])
                if yi > 0:
                    neighbors.append(occupied[xi, yi-1])
                if yi < Ny-1:
                    neighbors.append(occupied[xi, yi+1])
                if any(neighbors):
                    return xi, yi
        # fallback: first free cell
        for xi in range(Nx):
            for yi in range(Ny):
                if not occupied[xi, yi]:
                    return xi, yi
        return None

    seed = (0, 0)
    first = True
    last_rect_cell = None 

    S_integ = build_integral(obs1_values, obs2_values, obs1_weights)
    B_integ = build_integral(bkg1_values, bkg2_values, bkg1_weights)

    while any_free():

        if first:
            x0i, y0i = seed
            first = False
        else:
            # try to follow requested "start at X+1,0 or 0,Y+1"
            # prefer valid and free; if chosen invalid, fallback to other; if both invalid, find adjacent free cell
            if last_rect_cell is None:
                pos = find_free_adjacent()
                if pos is None:
                    break
                x0i, y0i = pos
            else:
                last_x1i, last_y1i = last_rect_cell
                opt1 = (last_x1i + 1, 0)   # right-of-last at y=0
                opt2 = (0, last_y1i + 1)   # above-last at x=0
                choices = [opt1, opt2]
                random.shuffle(choices)
                chosen = None
                for cx, cy in choices:
                    if 0 <= cx < Nx and 0 <= cy < Ny and not occupied[cx, cy]:
                        chosen = (cx, cy)
                        break
                if chosen is None:
                    # find any free cell adjacent to existing rectangles
                    pos = find_free_adjacent()
                    if pos is None:
                        break
                    x0i, y0i = pos
                else:
                    x0i, y0i = chosen

        if occupied[x0i, y0i]:
            pos = find_free_adjacent()
            if pos is None:
                break
            x0i, y0i = pos

        x1i, y1i = x0i, y0i

        while True:
            S_total = rect_sum(S_integ, x0i, x1i, y0i, y1i)
            B_total = rect_sum(B_integ, x0i, x1i, y0i, y1i)
            current_ams = do_AMS_simple(S_total, B_total, np.sqrt(B_total) if B_total > 0 else 0.0) if B_total > 0 else 0.0

            # check stopping condition:
            # - if rectangle meets both AMS and signal targets -> stop
            # - otherwise, if no growth direction is possible -> accept (last bin allowed)
            if (B_total > 0 and current_ams >= targetams and S_total >= targetsig):
                accepted = True
            else:
                accepted = False

            allowed_dirs = []

            if x1i + 1 < Nx:
                if not occupied[x1i + 1, y0i:y1i+1].any():
                    allowed_dirs.append("+x")
            if x0i - 1 >= 0:
                if not occupied[x0i - 1, y0i:y1i+1].any():
                    allowed_dirs.append("-x")
            if y1i + 1 < Ny:
                if not occupied[x0i:x1i+1, y1i + 1].any():
                    allowed_dirs.append("+y")
            if y0i - 1 >= 0:
                if not occupied[x0i:x1i+1, y0i - 1].any():
                    allowed_dirs.append("-y")

            # enforce origin constraint: if this rectangle started at (0,0) we should not expand negative
            if x0i == 0 and y0i == 0:
                # remove negative options if present (they shouldn't be)
                allowed_dirs = [d for d in allowed_dirs if d in ("+x", "+y")]

            # If we already accepted because AMS+S satisfied, stop without further growth.
            if accepted:
                break

            # If no allowed directions -> accept rectangle (last-bin rule), even if AMS not met
            if not allowed_dirs:
                #print("No allowed expansions left for this rectangle -> accepting (may be last bin).")
                break

            # choose a direction by flipping a coin among allowed directions
            dir_choice = random.choice(allowed_dirs)

            # apply the chosen expansion (grow rectangle by 1 cell)
            if dir_choice == "+x":
                x1i += 1
            elif dir_choice == "-x":
                x0i -= 1
            elif dir_choice == "+y":
                y1i += 1
            elif dir_choice == "-y":
                y0i -= 1
            else:
                # shouldn't happen
                pass

        mark_occupied(x0i, x1i, y0i, y1i)
        lo_x, hi_x, lo_y, hi_y = cell_to_coords(x0i, x1i, y0i, y1i)
        rectangles.append((lo_x, hi_x, lo_y, hi_y))
        
        if (opt.PRINT): print("ACCEPTED rectangle:", (lo_x, hi_x, lo_y, hi_y))

        last_rect_cell = (x1i, y1i)

    if (opt.PRINT):
        if not occupied.all():
            print("Warning: tiling did not fully cover the grid. Remaining free cells:",
                np.argwhere(~occupied).tolist())
        else:
            print("Tiling complete: fully covered.")


    return rectangles   

def resolve_overlaps(rects):

    changed = True
    rects = rects[:]

    while changed:
        changed = False
        new_rects = []
        skip = set()

        for i in range(len(rects)):
            if i in skip:
                continue
            x0a, x1a, y0a, y1a = rects[i]
            merged_rect = (x0a, x1a, y0a, y1a)

            for j in range(i+1, len(rects)):
                if j in skip:
                    continue
                x0b, x1b, y0b, y1b = rects[j]

                # overlap check
                if not (x1a <= x0b or x1b <= x0a or y1a <= y0b or y1b <= y0a):
                    # merge them
                    merged_rect = (
                        min(merged_rect[0], x0b),
                        max(merged_rect[1], x1b),
                        min(merged_rect[2], y0b),
                        max(merged_rect[3], y1b)
                    )
                    skip.add(j)
                    changed = True

            new_rects.append(merged_rect)

        rects = new_rects

    return rects

def check_matrix_2D(obs_name, binning_string, matrix, eff, res, dfs_sig, dfs_bkg, doubleDiff, threshold_base=0.0):

    rects = bins_to_array_2d(binning_string)

    while True:

        matrix = np.flipud(matrix)
        merged_any = False
        rects = bins_to_array_2d(binning_string)

        for i_diag in range(len(matrix)):
            
            # Diagonal rectangle
            rect = rects[i_diag]
            x0, x1, y0, y1 = rect

            if res[i_diag][0] == 1e-8 and res[i_diag][1] == 1e-8:
                delta_x = x1 - x0
                delta_y = y1 - y0
                arg_x = delta_x / (2 * sqrt(2) * res[i_diag][0])
                arg_y = delta_y / (2 * sqrt(2) * res[i_diag][1])
                P_expected = 0
            else:
                delta_x = x1 - x0
                delta_y = y1 - y0
                arg_x = delta_x / (2 * sqrt(2) * res[i_diag][0])
                arg_y = delta_y / (2 * sqrt(2) * res[i_diag][1])
                P_expected = eff * erf(arg_x) * erf(arg_y)

            P_target = max(P_expected, threshold_base)
            diag = matrix[i_diag, i_diag]
            
            if (opt.PRINT): print(f"Bin {rects[i_diag]}: eff = {eff:.2f} res_x={res[i_diag][0]:.2f} res_y={res[i_diag][1]:.2f} delta_x={delta_x:.2f} delta_y={delta_y:.2f} => P_expected={P_expected:.3f} diag={diag:.3f}")

            if diag >= P_target:
                continue  # bin is OK

            # merge with next bin
            if i_diag + 1 >= len(matrix):
                continue

            # merged rectangle
            x0_new = x0
            x1_new = rects[i_diag + 1][1]
            y0_new = y0
            y1_new = rects[i_diag + 1][3]
            rects[i_diag] = (x0_new, x1_new, y0_new, y1_new)
            rects.pop(i_diag + 1)

            rects = resolve_overlaps(rects)

            binning_string = bins_to_string_2d(rects)

            events, sig, mig, con, matrix = do_binning(dfs_sig, dfs_bkg, years, obs_name, binning_string, doubleDiff, False)
            
            res = get_res_2D(obs_name, binning_string, dfs_sig)

            merged_any = True
            break  # restart from first bin

        if not merged_any:
            break

    rects = resolve_overlaps(rects)
    binning_string = bins_to_string_2d(rects)
    return binning_string

def get_res_2D(obs_name, binning_string, dfs_sig):
    rects = bins_to_array_2d(binning_string)
    obs_name_1st, obs_name_2nd = [o.strip() for o in obs_name.split("vs")]

    reco_x_list, reco_y_list, gen_x_list, gen_y_list = [], [], [], []

    for year_dict in dfs_sig.values():
        for df in year_dict.values():
            df_sel = df[
                (df['ZZMass'] >= opt.m4lLower) & (df['ZZMass'] <= opt.m4lUpper) &
                (df['GENmass4l'] >= opt.m4lLower) & (df['GENmass4l'] <= opt.m4lUpper) &
                (df['passedFullSelection'] == 1) & (df['cuth4l_reco'] == True)
            ]

            rx = df_sel[observables[obs_name_1st]['obs_reco']].astype(float)
            ry = df_sel[observables[obs_name_2nd]['obs_reco']].astype(float)
            gx = df_sel[observables[obs_name_1st]['obs_gen']].astype(float)
            gy = df_sel[observables[obs_name_2nd]['obs_gen']].astype(float)

            if obs_name_1st == "rapidity4l":
                rx = rx.abs(); gx = gx.abs()
            if obs_name_2nd == "rapidity4l":
                ry = ry.abs(); gy = gy.abs()

            reco_x_list.append(rx); reco_y_list.append(ry)
            gen_x_list.append(gx); gen_y_list.append(gy)

    reco_x = pd.concat(reco_x_list, ignore_index=True)
    reco_y = pd.concat(reco_y_list, ignore_index=True)
    gen_x  = pd.concat(gen_x_list,  ignore_index=True)
    gen_y  = pd.concat(gen_y_list,  ignore_index=True)

    jet_vars = {"pTj1","pTj2","mjj","absdetajj","mHj","pTHj","pTHjj"}

    res = []
    for x1, x2, y1, y2 in rects:
        mask = (reco_x >= x1) & (reco_x < x2) & (reco_y >= y1) & (reco_y < y2)

        if obs_name_1st in jet_vars:
            mask &= (reco_x != -99) & (gen_x != -99)
        if obs_name_2nd in jet_vars:
            mask &= (reco_y != -99) & (gen_y != -99)

        dx = gen_x[mask] - reco_x[mask]
        dy = gen_y[mask] - reco_y[mask]

        #---------------------------------------
        # REMOVE OUTLIERS
        #---------------------------------------
        # drop NaN/inf first
        dx = dx[np.isfinite(dx)]
        dy = dy[np.isfinite(dy)]

        # ---- outlier mask (sigma clip) ----
        if dx.size >= 5:  # avoid doing this on tiny samples
            mu_x = np.mean(dx)
            sig_x = np.std(dx, ddof=1)
            if sig_x > 0:
                keep_x = np.abs(dx - mu_x) < 1.0 * sig_x
                dx = dx[keep_x]
        if dy.size >= 5:  # avoid doing this on tiny samples
            mu_y = np.mean(dy)
            sig_y = np.std(dy, ddof=1)
            if sig_y > 0:
                keep_y = np.abs(dy - mu_y) < 1.0 * sig_y
                dy = dy[keep_y]
        #---------------------------------------
        #---------------------------------------
    
        with np.errstate(all='ignore'):
            res_x = np.sqrt(np.nanmean(dx**2)) if len(dx) > 0 else 1e-8
            res_y = np.sqrt(np.nanmean(dy**2)) if len(dy) > 0 else 1e-8
        res.append((res_x, res_y))


    return res
'''