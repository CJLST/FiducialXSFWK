import numpy as np
import pandas as pd
import uproot #3 as uproot
from math import sqrt, log
import sys,os
import optparse
import itertools
import math
import json
import ROOT
import awkward as ak

from paths import path


signals_original = ['ggH125', 'VBFH125', 'ttH125', 'WminusH125', 'WplusH125', 'ZH125']
bkgs = ['ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701', 'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701', 'ggTo4tau_Contin_MCFM701', 'ZZTo4l'] 

key = 'ZZTree/candTree'
#key = 'Events'


# ------------------------------- FUNCTIONS TO GENERATE DATAFRAMES ----------------------------------------------------
# Weights for histogram
def weight(df, xsec, gen, lumi, additional = None):
    #Coefficient to calculate weights for histograms
    coeff = (lumi * 1000 * xsec) / gen
    #Reco
    weight_reco = (df.overallEventWeight) # * df.L1prefiringWeight * df.SFcorr)
    if additional == 'ggH':
        weight_reco *= df.ggH_NNLOPS_weight
    #elif additional == 'qqzz':
    #    weight_reco *= df.KFactor_EW_qqZZ*df.KFactor_QCD_qqZZ_M
    #elif additional == 'ggzz':
    #    weight_reco *= df.KFactor_QCD_ggZZ_Nominal
    weight_histo_reco = weight_reco * coeff
    #Columns in pandas
    df['weight_reco'] = weight_reco #Powheg
    df['weight_histo_reco'] = weight_histo_reco #Powheg
    # if additional == 'ggH':
    #     weight_reco_NNLOPS = weight_reco * df.ggH_NNLOPS_weight
    #     weight_histo_reco_NNLOPS = weight_histo_reco * df.ggH_NNLOPS_weight
    #     df['weight_reco_NNLOPS'] = weight_reco_NNLOPS #NNLOPS (only ggH)
    #     df['weight_histo_reco_NNLOPS'] = weight_histo_reco_NNLOPS #NNLOPS (only ggH)
    # else:
    #     df['weight_reco_NNLOPS'] = -1
    #     df['weight_histo_reco_NNLOPS'] = -1
    return df

# Uproot to generate pandas
def prepareTrees(year):
    d_sig = {}
    d_bkg = {}

    for bkg in bkgs:
        #fname = eos_path_sig + 'MC_samples_UL/%s_MELA' %year
        #fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'

        fname = path['eos_path_sig']+year+"_MC/"+bkg+"/ZZ4lAnalysis_SKIMMED.root"
        #fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/{year}_MC/{bkg}/ZZ4lAnalysis_SKIMMED_addJES.root"

        print(fname)
#         fname = eos_path_sig + '%i_MELA' %year
# #         if year == 2016:
# #             fname += '_CorrectBTag'
#         if (year == 2018) & (bkg == 'ZZTo4lext'):
#             fname += '/'+bkg+'1/'+bkg+'1_reducedTree_MC_'+str(year)+'.root'
#         else:
#             fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'
        d_bkg[bkg] = uproot.open(fname)[key]

    if year!='2016pre':
        for signal in signals_original:
            # if(opt.AC==False):
            #fname = eos_path_sig + 'MC_samples_UL/%s_MELA' %year
            # else:
    	        # fname = eos_path_sig + 'AC%i' %year
            # if (year == 2017) & (signal == 'VBFH125'):
            #     fname += '/'+signal+'ext/'+signal+'ext_reducedTree_MC_'+str(year)+'.root'
            # else:
            #     fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
            #fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'

            fname = path['eos_path_sig']+year+"_MC/"+signal+"/ZZ4lAnalysis_SKIMMED.root"
            #fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/{year}_MC/{signal}/ZZ4lAnalysis_SKIMMED_addJES.root"

            print(fname)
            d_sig[signal] = uproot.open(fname)[key]


    return d_sig, d_bkg


# Calculate cross sections
def xsecs(year):
    xsec_sig = {}
    xsec_bkg = {}
    #Hard-coded values
    xsec_sig['ggH125'] = 0.0133352
    xsec_sig['VBFH125'] = 0.0010381
    xsec_sig['ttH125'] = 0.0003639
    xsec_sig['WminusH125'] = 0.0001462
    xsec_sig['WplusH125'] = 0.0002305
    xsec_sig['ZH125'] = 0.0005321

    xsec_bkg['ZZTo4l'] = 1.2560000
    xsec_bkg['ggTo2e2mu_Contin_MCFM701'] = 0.0031914
    xsec_bkg['ggTo2e2tau_Contin_MCFM701'] = 0.0031914
    xsec_bkg['ggTo2mu2tau_Contin_MCFM701'] = 0.0031914
    xsec_bkg['ggTo4e_Contin_MCFM701'] = 0.0015854
    xsec_bkg['ggTo4mu_Contin_MCFM701'] = 0.0015854
    xsec_bkg['ggTo4tau_Contin_MCFM701'] = 0.0015854

    # d_sig, d_bkg = prepareTrees(year)
    # for signal in signals_original:
    #     xsec_sig[signal] = d_sig[signal].pandas.df('xsec').xsec[0]
    # for bkg in bkgs:
    #     xsec_bkg[bkg] = d_bkg[bkg].pandas.df('xsec').xsec[0]
    return xsec_sig, xsec_bkg


def add_fin_state_reco(i, j):
    if abs(i) == 121 and abs(j) == 121:
        fin = '4e'
    elif abs(i) == 169 and abs(j) == 169:
        fin = '4mu'
    elif (abs(i) == 121 and abs(j) == 169) or (abs(i) == 169 and abs(j) == 121):
        fin = '2e2mu'
    elif (abs(i) == 225 and abs(j) == 169) or (abs(i) == 169 and abs(j) == 225):
        fin = '2tau2mu'
    elif (abs(i) == 225 and abs(j) == 121) or (abs(i) == 121 and abs(j) == 225):
        fin = '2tau2e'
    elif abs(i) == 225 and abs(j) == 225:
        fin = '4tau'
    #elif abs(i) == 0 and abs(j) == 0:
    else:
        fin = 'other'
    return fin


# Get the "number" of MC events to divide the weights
def generators(year):
    gen_sig = {}
    gen_bkg = {}

    for bkg in bkgs:
        #fname = eos_path_sig + 'MC_samples_UL/%s_MELA' %year
        #fname += '/'+bkg+'/'+bkg+'_reducedTree_MC_'+str(year)+'.root'

        fname = path['eos_path_sig']+year+"_MC/"+bkg+"/ZZ4lAnalysis_SKIMMED.root"
        #fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/{year}_MC/{bkg}/ZZ4lAnalysis_SKIMMED_addJES.root"

        gen_bkg[bkg] = uproot.open(fname)["Counters"].values()[39] # spencer
        
    if year!='2016pre':
        for signal in signals_original:
            # if(opt.AC==False):
            #fname = eos_path_sig + 'MC_samples_UL/%s_MELA' %year
            # else:
    	        # fname = eos_path_sig + 'AC%i' %year
            # if (year == 2017) & (signal == 'VBFH125'):
            #     fname += '/'+signal+'ext/'+signal+'ext_reducedTree_MC_'+str(year)+'.root'
            # else:
            #     fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'
            #fname += '/'+signal+'/'+signal+'_reducedTree_MC_'+str(year)+'.root'

            fname = path['eos_path_sig']+year+"_MC/"+signal+"/ZZ4lAnalysis_SKIMMED.root"
            #fname = f"/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/{year}_MC/{signal}/ZZ4lAnalysis_SKIMMED_addJES.root"

            gen_sig[signal] = uproot.open(fname)["Counters"].values()[39] # spencer 

    return gen_sig, gen_bkg

def createDataframe(jesNames, year, dataFrame,isBkg,gen,xsec,signal,lumi,obs_reco,obs_reco_2nd='None'):

    b_sig = ['overallEventWeight', 'Z1Flav', 'Z2Flav', 'ZZMass' ]


    if signal == 'ggH125':
        b_sig.append('ggH_NNLOPS_weight') #Additional entry for the weight in case of ggH

    #elif signal == 'ZZTo4l':
        #b_sig.append('KFactor_EW_qqZZ')
        #b_sig.append('KFactor_QCD_qqZZ_M')

    #elif 'gg' in signal:
    #    b_sig.append('KFactor_QCD_ggZZ_Nominal')

    if obs_reco == 'ZZMass':
        b_sig.pop('ZZMass')
    if obs_reco_2nd == 'ZZMass':
        b_sig.pop('ZZMass')

    b_sig.append(obs_reco)
    if obs_reco_2nd != 'None': b_sig.append(obs_reco_2nd)

    for name in jesNames:

        if "year" in name:
            name1st = name.split("_")[0]
            name = name1st + "_" + year

        b_sig.append(obs_reco+"_"+name+"_ScaleUp")
        b_sig.append(obs_reco+"_"+name+"_ScaleDn")

        if obs_reco_2nd != 'None':
            b_sig.append(obs_reco_2nd+"_"+name+"_ScaleUp")
            b_sig.append(obs_reco_2nd+"_"+name+"_ScaleDn")

    df_np = dataFrame.arrays(library="np")
    df = pd.DataFrame({var: df_np[var] for var in b_sig})

    df['gen'] = gen
    df['xsec'] = xsec
    
    df["FinState_reco"] = [add_fin_state_reco(i, j) for i, j in zip(df.Z1Flav, df.Z2Flav)]

    if signal != 'ggH125':
        df = weight(df, xsec, gen, lumi)

    elif 'ZZTo' in signal and isBkg:
        df = weight(df, xsec, gen, lumi, 'qqzz')
        #df = df.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M'])

    elif 'ggTo' in signal and isBkg:
        df = weight(df, xsec, gen, lumi, 'ggzz')
        #df = df.drop(columns=['KFactor_EW_qqZZ','KFactor_QCD_qqZZ_M'])

    else:
        df = weight(df, xsec, gen, lumi, 'ggH')
        df = df.drop(columns=['ggH_NNLOPS_weight'])

    return df


# Set up data frames
def dataframes(jesNames, year, doubleDiff, obs_reco, obs_reco_2nd):
    if year == '2016pre':
        lumi_bkg = 19.52
    elif year == '2016post':
        lumi_bkg = 16.81
        lumi_sig = 35.9
    elif year == '2017':
        lumi_sig = 41.5
        lumi_bkg = 41.5
    elif year == '2018':
        lumi_sig = 59.7
        lumi_bkg = 59.7
    elif year == '2022':
        lumi_sig = 7.9804
        lumi_bkg = 7.9804
    elif year == '2022EE':
        lumi_sig = 26.6728
        lumi_bkg = 26.6728
    elif year == '2023preBPix':
        lumi_sig = 17.794
        lumi_bkg = 17.794
    elif year == '2023postBPix':
        lumi_sig = 9.451
        lumi_bkg = 9.451
    elif year == '2024':
        lumi_sig = 109.08
        lumi_bkg = 109.08
        
    d_df_sig = {}
    d_df_bkg = {}
    d_sig, d_bkg = prepareTrees(year)
    gen_sig, gen_bkg = generators(year)
    xsec_sig, xsec_bkg = xsecs(year)

    for bkg in bkgs:
        print ('Processing', bkg, year)
        if doubleDiff:
            d_df_bkg[bkg] = createDataframe(jesNames, year, d_bkg[bkg],True,gen_bkg[bkg],xsec_bkg[bkg],bkg,lumi_bkg,obs_reco,obs_reco_2nd)
        else:
            d_df_bkg[bkg] = createDataframe(jesNames, year, d_bkg[bkg],True,gen_bkg[bkg],xsec_bkg[bkg],bkg,lumi_bkg,obs_reco)
        print ('Background created')

    if year!='2016pre':
        for signal in signals_original:
            print ('Processing', signal, year)
            if doubleDiff:
                d_df_sig[signal] = createDataframe(jesNames, year, d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi_sig,obs_reco,obs_reco_2nd)
            else:
                d_df_sig[signal] = createDataframe(jesNames, year, d_sig[signal],False,gen_sig[signal],xsec_sig[signal],signal,lumi_sig,obs_reco)
            print ('Signal created')


    return d_df_sig, d_df_bkg


# Merge WplusH125 and WminusH125
def skim_df(jesNames, year, doubleDiff, obs_reco, obs_reco_2nd = ''):
    d_df_sig, d_df_bkg = dataframes(jesNames, year, doubleDiff, obs_reco, obs_reco_2nd)
    d_skim_sig = {}
    d_skim_bkg = {}

    frames = []
    for bkg in bkgs:
        if (bkg == 'ZZTo4l'):
            d_skim_bkg['qqzz'] = d_df_bkg[bkg]
        else:
            frames.append(d_df_bkg[bkg])
    d_skim_bkg['ggzz'] = pd.concat(frames)

    if year!='2016pre':
        frames = []
        for signal in signals_original:
            if (signal == 'WplusH125') or (signal == 'WminusH125'):
                frames.append(d_df_sig[signal])
            else:
                d_skim_sig[signal] = d_df_sig[signal]
        d_skim_sig['WH125'] = pd.concat(frames)


    print ('%s SKIMMED df CREATED' %year)
    return d_skim_sig, d_skim_bkg
