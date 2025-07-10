import os, sys
import numpy as np
import pandas as pd
import uproot as uproot
import ROOT
import argparse
from array import array

from observables import observables
from binning import binning
from paths import path

print()
print('Welcome in RunBinning!')
print()

def parseOptions():
    global opt, args, runAllSteps

    parser = argparse.ArgumentParser(description="Run binning for fiducial cross-section measurement")

    parser.add_argument('--obsName', required=True, type=str, help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_argument('--obsBins', required=True, type=str, help='Bin boundaries for the differential measurement separated by "|", e.g. "|0|50|100|". Use empty string for default.')
    parser.add_argument('--year', required=True, type=str, help='Year (e.g., "2022", "2023preBPix", etc.)')
    parser.add_argument('--m4lLower', dest='m4lLower', type=int, default=105, help='Lower bound for m4l (default: 105)')
    parser.add_argument('--m4lUpper', dest='m4lUpper', type=int, default=160, help='Upper bound for m4l (default: 160)')

    opt = parser.parse_args()
    return opt

opt = parseOptions()

# ---------------------------------------------------- FUNCTIONS TO GENERATE DATAFRAMES ----------------------------------------------------

# Weights for histogram
def weight(df, xsec, gen, lumi):
    # xsec is in overallEventWeight now
    weight = (lumi * 1000 * df.overallEventWeight * df.dataMCWeight)/gen
    df['weight'] = weight
    return df

# Uproot to generate pandas
def prepareTrees(samples, year):
    df = {}
    for sample in samples:
        fname = path['eos_path_sig']+year+"_MC/"+sample+"/ZZ4lAnalysis_SKIMMED.root"
        df[sample] = uproot.open(fname)[key]
    return df

# Calculate cross sections
def xsecs(samples, year):
    xsec_ = {}
    df = prepareTrees(samples, year)
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

    return xsec_

# Get the "number" of MC events to divide the weights
def generators(samples, year):
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
    elif (abs(i) == 121 and abs(j) == 169):
        fin = '2e2mu'
    elif (abs(i) == 169 and abs(j) == 121):
        fin = '2mu2e'
    else: print('Problem with add_fin_state')
    return fin

def dataframes(samples, year, year_mc):
    if year_mc == '2022':
        lumi = 7.9804
    elif year_mc == '2022EE':
        lumi = 26.6728
    elif year_mc == '2023preBPix':
        lumi = 17.794
    elif year_mc == '2023postBPix':
        lumi = 9.451

    d = {}
    df_all = prepareTrees(samples, year_mc)
    gen_df = generators(samples, year_mc)
    xsec_df = xsecs(samples, year_mc)

    for sample in samples:
        gen = gen_df[sample]
        xsec = xsec_df[sample]
        df_np = df_all[sample].arrays(vars, library="np")
        df = pd.DataFrame({var: df_np[var] for var in vars})
        df['FinState'] = [add_fin_state(i, j) for i, j in zip(df.Z1Flav, df.Z2Flav)]
        d[sample] = weight(df, xsec, gen, lumi)

    return d

def skim_df(samples, year, year_mc, type):
    df = dataframes(samples, year, year_mc)
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

    if (year=="2022"): data = '/eos/user/s/sellissp/HZZ/SAMPLES/062025/2022_Data/Data_eraCD_preEE_SKIMMED.root'
    if (year=="2022EE"): data = '/eos/user/s/sellissp/HZZ/SAMPLES/062025/2022_Data/Data_eraEFG_postEE_SKIMMED.root'
    if (year=="2023preBPix"): data = '/eos/user/s/sellissp/HZZ/SAMPLES/062025/2023_Data/Data_eraC_preBPix_SKIMMED.root'
    if (year=="2023postBPix"): data = '/eos/user/s/sellissp/HZZ/SAMPLES/062025/2023_Data/Data_eraD_postBPix_SKIMMED.root'

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


# ---------------------------------------------------- FUNCTIONS FOR BINNING OPTIMIZATION ----------------------------------------------------

def do_AMS(sig_counts, bkg_counts, bkg_errors):

    qqZZ, ggZZ, ZX = bkg_counts['qqZZ'], bkg_counts['ggZZ'], bkg_counts['ZX']
    qqZZ_err, ggZZ_err, ZX_errors = bkg_errors['qqZZ'], bkg_errors['ggZZ'], bkg_errors['ZX']

    if len(sig_counts) != len(qqZZ) or len(sig_counts) != len(ggZZ) or len(sig_counts) != len(ZX):
        raise ValueError("sig and bkg counts must have the same length")

    AMS = []
    for i in range(len(sig_counts)):

        s = sig_counts[i]
        b = qqZZ[i] + ggZZ[i] + ZX[i]

        #sigma_b = np.sqrt( (qqZZ_err[i]**2/qqZZ[i]**2 + ggZZ_err[i]**2/ggZZ[i]**2 + ZX_err[i]**2/ZX[i]**2 ) / (qqZZ[i]**2+ggZZ[i]**2+ZX[i]**2) ) # FROM HIG 21 09 AN
        sigma_b = np.sqrt(qqZZ_err[i]**2 + ggZZ_err[i]**2 + ZX_errors[i]**2)

        term1 = (s+b) * np.log( ((s+b)*(b+sigma_b**2))/(b**2+(s+b)*sigma_b**2) )
        term2 = (b**2/sigma_b**2) * np.log(1 + (s*sigma_b**2)/(b*(b+sigma_b**2)))

        AMS.append( np.sqrt(2*term1 - 2*term2) )
    
    return AMS

def get_bin_index(x, bins):
    for i in range(len(bins) - 1):
        if bins[i] <= x < bins[i + 1]:
            return i
    if x == bins[-1]:
        return len(bins) - 2
    return -1  # out of range

def do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff):

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins, dtype=float)
        hist = ROOT.TH2F("tmp_hist", "tmp_hist", nbins, bin_array, nbins, bin_array)

    else:
        obs_name_1st = obs_name.split('vs')[0].strip()
        obs_name_2nd = obs_name.split('vs')[1].strip()
        
        x_str, y_str = [s.strip() for s in obsBins.split('vs')]
        xbins = [float(x) for x in x_str.strip('|').split('|') if x]
        ybins = [float(y) for y in y_str.strip('|').split('|') if y]
        
        nbins_x = len(xbins) - 1
        nbins_y = len(ybins) - 1
        xbin_array = np.array(xbins, dtype=float)
        ybin_array = np.array(ybins, dtype=float)

        nbins = nbins_x * nbins_y
        hist = ROOT.TH2F("tmp_hist", "tmp_hist", nbins, 0, nbins, nbins, 0, nbins)

    for year in years:
        for sig in dfs_sig[year].keys():
            df_sel = dfs_sig[year][sig][(dfs_sig[year][sig]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['ZZMass'] <= opt.m4lUpper) &
                                        (dfs_sig[year][sig]['GENmass4l'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['GENmass4l'] <= opt.m4lUpper)]

            if not doubleDiff:
                for x, y, w in zip(df_sel[observables[obs_name]['obs_gen']], 
                                   df_sel[observables[obs_name]['obs_reco']], 
                                   df_sel['weight']):
                    hist.Fill(x, y, w)
            if doubleDiff:
                for x1, x2, y1, y2, w in zip(df_sel[observables[obs_name_1st]['obs_gen']],
                                             df_sel[observables[obs_name_2nd]['obs_gen']],
                                             df_sel[observables[obs_name_1st]['obs_reco']],
                                             df_sel[observables[obs_name_2nd]['obs_reco']],
                                             df_sel['weight']):

                    gen_bin1 = get_bin_index(x1, xbins)
                    gen_bin2 = get_bin_index(x2, ybins)
                    reco_bin1 = get_bin_index(y1, xbins)
                    reco_bin2 = get_bin_index(y2, ybins)

                    if -1 in (gen_bin1, gen_bin2, reco_bin1, reco_bin2):
                        continue  # skip out-of-range

                    gen_flat_bin = gen_bin1 * nbins_y + gen_bin2
                    reco_flat_bin = reco_bin1 * nbins_y + reco_bin2

                    hist.Fill(gen_flat_bin + 0.5, reco_flat_bin + 0.5, w)  # center of bin

    matrix = np.zeros((nbins, nbins), dtype=float)
    for i in range(1, nbins + 1):
        for j in range(1, nbins + 1):
            matrix[i - 1, j - 1] = max(0.0, hist.GetBinContent(i, j))

    matrix = np.flipud(matrix)

    # normalize the matrix by gen count
    column_sums = matrix.sum(axis=0)
    matrix = np.divide( matrix, column_sums )

    np.set_printoptions(formatter={'float_kind': lambda x: f"{x:.2f}"})

    return matrix


def do_migration_score(matrix):

    matrix = np.flipud(matrix)

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    total = np.sum(np.abs(matrix))
    diag_sum = np.sum(np.abs(np.diag(matrix)))

    return diag_sum / total

def do_binning(dfs_sig, dfs_bkg, years, obs_name, obsBins, doubleDiff):
    
    print()

    if not doubleDiff:
        bins = [float(x) for x in obsBins.strip('|').split('|') if x]
        nbins = len(bins) - 1
        bin_array = np.array(bins)
    if doubleDiff:
        obs_name_1st = obs_name.split('vs')[0].strip()
        obs_name_2nd = obs_name.split('vs')[1].strip()
        x_str, y_str = [s.strip() for s in obsBins.split('vs')]
        xbins = [float(x) for x in x_str.strip('|').split('|') if x]
        ybins = [float(y) for y in y_str.strip('|').split('|') if y]
        nbins_x, nbins_y = len(xbins) - 1, len(ybins) - 1
        xbin_array, ybin_array = np.array(xbins, dtype=float), np.array(ybins, dtype=float)
        nbins = nbins_x * nbins_y
        print("binning:", "(" + ", ".join(f"[{xbins[i]}, {xbins[i+1]}] [{ybins[j]}, {ybins[j+1]}]" for i in range(len(xbins)-1) for j in range(len(ybins)-1)) + ")")

    sum_sig_weighted_counts = {sig: np.zeros(nbins) for sig in dfs_sig[years[0]].keys()}
    total_sig_weighted_counts, total_sig_bin_error2 = np.zeros(nbins), np.zeros(nbins)
    sum_bkg_weighted_counts = {bkg: np.zeros(nbins) for bkg in dfs_bkg[years[0]].keys()}
    total_bkg_weighted_counts, total_bkg_bin_error2 = np.zeros(nbins), np.zeros(nbins)

    for year in years:
        for sig in dfs_sig[year].keys():

            df_sel = dfs_sig[year][sig][(dfs_sig[year][sig]['ZZMass'] >= opt.m4lLower) & 
                                        (dfs_sig[year][sig]['ZZMass'] <= opt.m4lUpper)]

            #df_sel = df_sel[df_sel['FinState'] == '4e' ] # keep only 4e events
            #df_sel = df_sel[df_sel['FinState'] == '4mu' ] # keep only 4mu events
            #df_sel = df_sel[df_sel['FinState'] == '2e2mu'] # keep only 2e2mu events
            #df_sel = df_sel[df_sel['FinState'] == '2mu2e'] # keep only 2mu2e events

            df_weights = df_sel['weight']
        
            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_sig_{sig}_{year}", "", nbins, bin_array) 
                df_obs = df_sel[observables[obs_name]['obs_reco']]
                if obs_name == 'rapidity4l': df_obs = np.abs(df_obs)
                for val, wgt in zip(df_obs, df_weights):
                    hist.Fill(val, wgt)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])

                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

            if doubleDiff:
                hist = ROOT.TH2F(f"hist_sig_{sig}_{year}", "", nbins_x, xbin_array, nbins_y, ybin_array)
                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]
                if obs_name_1st == 'rapidity4l': df_obs_x = np.abs(df_obs_x)
                if obs_name_2nd == 'rapidity4l': df_obs_y = np.abs(df_obs_y)
                

                for val_x, val_y, wgt in zip(df_obs_x, df_obs_y, df_weights):
                    hist.Fill(val_x, val_y, wgt)

                weighted_counts = np.array([
                    hist.GetBinContent(i + 1, j + 1)
                    for i in range(nbins_x)
                    for j in range(nbins_y)
                ])

                #bin_errors = np.array([
                #    hist.GetBinError(i + 1, j + 1)
                #    for i in range(nbins_x)
                #    for j in range(nbins_y)
                #])
                bin_errors = np.sqrt(weighted_counts)


            sum_sig_weighted_counts[sig] += weighted_counts
            total_sig_weighted_counts += weighted_counts
            total_sig_bin_error2 += bin_errors ** 2

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

            df_sel = dfs_bkg[year][bkg][(dfs_bkg[year][bkg]['ZZMass'] >= opt.m4lLower) & (dfs_bkg[year][bkg]['ZZMass'] <= opt.m4lUpper)]

            #if bkg == 'ZX':
            #    df_sel = df_sel[df_sel['FinState'] == 3]
            #else:
            #    df_sel = df_sel[df_sel['FinState'] == '2mu2e']

            if bkg == 'ZX': 
                df_weights = df_sel['ZX_yield']
            else: 
                df_weights = df_sel['weight']
                
            if not doubleDiff:
                hist = ROOT.TH1F(f"hist_bkg_{bkg}_{year}", "", nbins, bin_array)
                df_obs = df_sel[observables[obs_name]['obs_reco']]
                if obs_name == 'rapidity4l': df_obs = np.abs(df_obs)
                for val, wgt in zip(df_obs, df_weights):
                    hist.Fill(val, wgt)

                weighted_counts = np.array([hist.GetBinContent(i + 1) for i in range(nbins)])
                
                #bin_errors = np.array([hist.GetBinError(i + 1) for i in range(nbins)])
                bin_errors = np.sqrt(weighted_counts)

            if doubleDiff:
                hist = ROOT.TH2F(f"hist_bkg_{bkg}_{year}", "", nbins_x, xbin_array, nbins_y, ybin_array)
                df_obs_x = df_sel[observables[obs_name_1st]['obs_reco']]
                df_obs_y = df_sel[observables[obs_name_2nd]['obs_reco']]
                if obs_name_1st == 'rapidity4l': df_obs_x = np.abs(df_obs_x)
                if obs_name_2nd == 'rapidity4l': df_obs_y = np.abs(df_obs_y)
                for val_x, val_y, wgt in zip(df_obs_x, df_obs_y, df_weights):
                    hist.Fill(val_x, val_y, wgt)

                weighted_counts = np.array([
                    hist.GetBinContent(i + 1, j + 1)
                    for i in range(nbins_x)
                    for j in range(nbins_y)
                ])
                
                #bin_errors = np.array([
                #    hist.GetBinError(i + 1, j + 1)
                #    for i in range(nbins_x)
                #    for j in range(nbins_y)
                #])
                bin_errors = np.sqrt(weighted_counts)

            sum_bkg_weighted_counts[bkg] += weighted_counts
            total_bkg_weighted_counts += weighted_counts
            total_bkg_bin_error2 += bin_errors**2
            bkg_weighted_counts[bkg] = weighted_counts
            bkg_binerrors[bkg] = bin_errors
            

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
    print(f"AMS per bin: {AMS}")
    print()

    gen_reco = do_migration(dfs_sig, years, obs_name, obsBins, doubleDiff)
    print(f"migration matrix:")
    print(gen_reco)
    print()
    print(f"migration score: {do_migration_score(gen_reco)}")
    print(f"condition number: {np.linalg.cond(gen_reco)}")

# -----------------------------------------------------------------------------------------                                                                                                                                                                                                                                                 
# ------------------------------- MAIN ----------------------------------------------------                                                                                                                                                                                                                                                 
# -----------------------------------------------------------------------------------------

sigs = ['ggH125', 'VBFH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
bkgs = ['ZZTo4l', 'ggTo2e2mu_Contin_MCFM701','ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']

#sigs = ['VBFH125']
#bkgs = ['ggH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125', 'ZZTo4l', 'ggTo2e2mu_Contin_MCFM701','ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701']


#eos_path_FR = path['eos_path_FR']
eos_path = path['eos_path']
key = 'ZZTree/candTree'

obs_name = opt.obsName
obs_bins = opt.obsBins

doubleDiff = 'vs' in obs_name

vars = ['overallEventWeight', 'dataMCWeight', 'ZZMass', 'Z1Flav', 'Z2Flav', 'LepLepId', 'LepEta', 'LepPt']

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

if (opt.year == '2022full'):
    years_MC = ['2022', '2022EE']
    years = ["2022", "2022EE"]
if (opt.year == '2023full'):
    years_MC = ['2023preBPix', '2023postBPix']
    years = ["2023preBPix", "2023postBPix"]

if (opt.year == 'Run3'):
    years_MC = ['2022', '2022EE', '2023preBPix', '2023postBPix']
    years = ["2022", "2022EE", "2023preBPix", "2023postBPix"]


# --------------------------- RUN --------------------------------------------


dfs_sig = {}
dfs_bkg = {}
                                                                                                                                                                                                                                                                                                   
for year, year_mc in zip(years, years_MC):
    FR_mu_EB, FR_mu_EE, FR_e_EB, FR_e_EE = openFR(year_mc)
    df_bkg = skim_df(bkgs, year, year_mc, 'bkg')
    dfs_bkg[year] = df_bkg

vars.append('GENmass4l')
if not doubleDiff: 
    vars.append(observables[obs_name]['obs_gen'])
if doubleDiff:
    vars.append(observables[obs_name_1st]['obs_gen'])
    vars.append(observables[obs_name_2nd]['obs_gen'])

for year, year_mc in zip(years, years_MC):
    df_sig = skim_df(sigs, year, year_mc, 'sig')
    dfs_sig[year] = df_sig

do_binning(dfs_sig, dfs_bkg, years, obs_name, obs_bins, doubleDiff)