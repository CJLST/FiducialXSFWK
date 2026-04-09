#!/usr/bin/env python
# coding: utf-8

import ROOT

import sys
import os
import importlib.util
import uproot 
import awkward as ak
from collections import defaultdict

import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

import optparse # spencer

sys.path.append('../helperstuff/') # spencer
from paths import path
from split import split
from observables import observables # spencer
from binning import binning # spencer

#sys.path.append('/eos/user/m/mbonanom/run3_trees/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/FiducialXSFWK/inputs')
sys.path.append('../inputs/')

# spencer
def parseOptions(): 

    global opt, args, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('',   '--obsName',  dest='OBSNAME',  type='string',default='pT4l',   help='Name of the observable, supported: "inclusive", "pT4l", "eta4l", "massZ2", "nJets"')
    parser.add_option('',   '--year',  dest='YEAR',  type='string', default='2022',   help='Year -> 2016 or 2017 or 2018 or Full')
    parser.add_option('',   '--split', action='store_true', dest='SPLIT', default=False,   help='split eras')
    parser.add_option('',   '--merge', action='store_true', dest='MERGE', default=False,   help='2022full - 2023full - Run3')
    parser.add_option('',   '--ZZfloating', action='store_true', dest='ZZ', default=False,   help='ZZ floating')

    global opt, args
    (opt, args) = parser.parse_args()

global opt, args, runAllSteps
parseOptions()
# spencer

def get_zh_chan_cuts(channel, tree):

    genz_idx0 = abs(tree['GENZ_DaughtersId'][:,0])
    genz_idx1 = abs(tree['GENZ_DaughtersId'][:,1])
    genz_idx2 = abs(tree['GENZ_DaughtersId'][:,2])
    
    mom_id0 = tree['GENZ_MomId'][:,0]
    mom_id1 = tree['GENZ_MomId'][:,1]
    mom_id2 = tree['GENZ_MomId'][:,2]

    if channel == "4l":
        sel_l0l1 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (genz_idx0 != genz_idx1))
        sel_l0l2 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx0 != genz_idx2))
        sel_l1l2 = ((mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx1 != genz_idx2))
        cutchan_gen_out_4e = ((genz_idx0 == 11) & (genz_idx1 == 11) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 11) & (genz_idx2 == 11) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 11) & (genz_idx2 == 11) & (mom_id1 == 25) & (mom_id2 == 25))
        cutchan_gen_out_4mu = ((genz_idx0 == 13) & (genz_idx1 == 13) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 13) & (genz_idx2 == 13) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 13) & (genz_idx2 == 13) & (mom_id1 == 25) & (mom_id2 == 25))
        cutchan_gen_out_2e2mu = sel_l0l1 | sel_l0l2 | sel_l1l2
        cutchan_gen_out = cutchan_gen_out_4e | cutchan_gen_out_4mu | cutchan_gen_out_2e2mu
    elif channel == "4e":
        cutchan_gen_out = ((genz_idx0 == 11) & (genz_idx1 == 11) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 11) & (genz_idx2 == 11) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 11) & (genz_idx2 == 11) & (mom_id1 == 25) & (mom_id2 == 25))
    elif channel == "4mu":
        cutchan_gen_out = ((genz_idx0 == 13) & (genz_idx1 == 13) & (mom_id0 == 25) & (mom_id1 == 25)) | ((genz_idx0 == 13) & (genz_idx2 == 13) & (mom_id0 == 25) & (mom_id2 == 25)) | ((genz_idx1 == 13) & (genz_idx2 == 13) & (mom_id1 == 25) & (mom_id2 == 25))
    elif channel == "2e2mu":
        sel_l0l1 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (genz_idx0 != genz_idx1))
        sel_l0l2 = ((mom_id0 == 25) & ((genz_idx0 == 11) | (genz_idx0 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx0 != genz_idx2))
        sel_l1l2 = ((mom_id1 == 25) & ((genz_idx1 == 11) | (genz_idx1 == 13)) & (mom_id2 == 25) & ((genz_idx2 == 11) | (genz_idx2 == 13)) & (genz_idx1 != genz_idx2))
        cutchan_gen_out = sel_l0l1 | sel_l0l2 | sel_l1l2
        
    return cutchan_gen_out

def get_chan_cuts(channel, tree):

    genz_idx0 = abs(tree['GENZ_DaughtersId'][:,0])
    genz_idx1 = abs(tree['GENZ_DaughtersId'][:,1])

    if channel == "4l":
        cutchan_gen_out = (((genz_idx0 == 11) | (genz_idx0 == 13)) & ((genz_idx1 == 11) | (genz_idx1 == 13)))
    elif channel == "4e":
        cutchan_gen_out = (genz_idx0 == 11) & (genz_idx1 == 11)
    elif channel == "4mu":
        cutchan_gen_out = (genz_idx0 == 13) & (genz_idx1 == 13)
    elif channel == "2e2mu":
        cutchan_gen_out = (((genz_idx0 == 11) & (genz_idx1 == 13)) | ((genz_idx0 == 13) & (genz_idx1 == 11)))

    return cutchan_gen_out

def m4l_cuts(tree, obs_gen, obs_gen_low, obs_gen_high, doubleDiff, obs_gen_2nd=None, obs_gen_2nd_low=None, obs_gen_2nd_high=None):

    cutm4l_gen = (tree['GENmass4l'] > m4l_low) & (tree['GENmass4l'] < m4l_high)
    
    if doubleDiff:

        obs_val = abs(tree[obs_gen]) if 'rapidity4l' in obs_gen else tree[obs_gen]
        if obs_gen_2nd is not None:
            obs_val_2nd = abs(tree[obs_gen_2nd]) if 'rapidity4l' in obs_gen_2nd else tree[obs_gen_2nd]
            cutobs_gen = ((obs_val >= obs_gen_low) & (obs_val < obs_gen_high) &
                          (obs_val_2nd >= obs_gen_2nd_low) & (obs_val_2nd < obs_gen_2nd_high))
        else:
            raise ValueError("doubleDiff=True requires obs_gen_2nd and its bounds")
    else:
        obs_val = abs(tree[obs_gen]) if 'rapidity4l' in obs_gen else tree[obs_gen]
        cutobs_gen = (obs_val >= obs_gen_low) & (obs_val < obs_gen_high)

    return cutm4l_gen, cutobs_gen

def build_histos(tree, sel, w_gen, w_nom, w_scale):

    h_nom = ak.sum(w_gen[sel] * w_nom[sel])

    h_scale = {}

    for i in w_scale:
        h_scale[i] = ak.sum(w_gen[sel] * w_scale[i][sel])

    return h_nom, h_scale


def get_qcd_weights(qcd_weights):
    w_nom = qcd_weights[:, 4]
    w_scale = {}
    for i in range(9):
        if ((i == 2) | (i == 4) | (i == 6)): continue
        w_scale[i] = qcd_weights[:, i]

    return w_nom, w_scale

def get_pdf_weights(pdf_weights):
    w_nom = pdf_weights[:, 0]
    w_scale = {}
    for i in range(1,101):
        w_scale[i] = pdf_weights[:, i]
        
    w_as = {
            0: pdf_weights[:, -1],
            1: pdf_weights[:, -2]
           }

    return w_nom, w_scale, w_as

def get_scale_unc(process, channel, tree, scale, observable, nnlops = False):

    if doubleDiff:
        obs_gen, obs_gen_low, obs_gen_high, obs_gen_2nd, obs_gen_2nd_low, obs_gen_2nd_high = observable
        cutm4l_gen, cutobs_gen = m4l_cuts(tree, obs_gen, obs_gen_low, obs_gen_high, True, obs_gen_2nd, obs_gen_2nd_low, obs_gen_2nd_high)
    else:
        obs_gen, obs_gen_low, obs_gen_high = observable
        cutm4l_gen, cutobs_gen = m4l_cuts(tree, obs_gen, obs_gen_low, obs_gen_high, False)

    if process == "ZH125":
        cutchan_gen = get_zh_chan_cuts(channel, tree)
    else:
        cutchan_gen = get_chan_cuts(channel, tree)

    full_sel = cutm4l_gen & cutobs_gen & cutchan_gen & tree['passedFiducial'] == 1
    
    fname = path['eos_path_sig']+year+"_MC/"+process+"/ZZ4lAnalysis.root"

    with uproot.open(f'{fname}') as f: # spencer
        tree = f['AllEvents'].arrays() # spencer

    #w_gen = tree['genHEPMCweight']
    w_gen = tree['Generator_weight'] # spencer
    
    if nnlops:
        #w_nnlops = tree['ggH_NNLOPS_weight']
        w_nnlops = tree['ggH_NNLOPS_Weight'] # spencer
        w_gen = w_gen * w_nnlops
        
    if scale == "qcd":
        w_name = "LHEScaleWeight"
        scale_weights = tree[w_name]
        w_nom, w_scale = get_qcd_weights(scale_weights)
        w_var = w_scale
    elif ((scale == "pdf") or (scale == "as")):
        w_name = "LHEPdfWeight"
        scale_weights = tree[w_name]
        w_nom, w_scale, w_as = get_pdf_weights(scale_weights)
        if scale == "pdf":
            w_var = w_scale
        else:
            w_var = w_as
    else:
        assert '... Unsupported type of uncertainty! Supported types are "qcd", "pdf", or "as" ...'

    h_nom, h_scale = build_histos(tree, full_sel, w_gen, w_nom, w_var)
    h_tot, h_tot_scale = build_histos(tree, cutchan_gen, w_gen, w_nom, w_var)
        
    return h_nom, h_scale, h_tot, h_tot_scale

def load_tree(fname):
    with uproot.open(f'{fname}') as f:
        #tree = f['ZZTree/candTree'].arrays()
        tree = f['ZZTree/candTree'].arrays() # spencer
        #tree_failed = f['ZZTree/candTree_failed'].arrays()
        tree_failed = f['ZZTree/candTree_failed'].arrays() # spencer
        
    tree_tot = ak.concatenate([tree, tree_failed])
    return tree_tot
'''
def get_th_xsec(process, obs_gen, suffix, year):

    # Angular variables added by Martina
    obs_name_dict = {"GENmass4l": "mass4l", "GENpT4l": "pT4l", "GENrapidity4l": "rapidity4l", "GENcostheta1": "costhetaZ1", "GENcostheta2": "costhetaZ2", "GENPhi": "phi", "GENPhi1": "phi1", "GENcosthetastar": "costhetastarZZ", "GENmassZ1": "massZ1", "GENmassZ2": "massZ2", "GENpTj1": "pTj1", "GENpTj2": "pTj2", "GENmjj": "mjj", "GENabsdetajj": "absdetajj", "GENdphijj": "dphijj", "GENmHj":  "mHj", "GENpTHj": "pTHj", "GENpTHjj": "pTHjj", "GENNj": "Nj", "GENTBMax": "TBMax", "GENTCMax": "TCMax"}
    obs_name = obs_name_dict[obs_gen]
    
    # Determine module name
    if ('ZH' not in process) and ('W' not in process):
        module_name = f'fidXS_{suffix}{obs_name}_{process.split("125")[0]}_{year}'
    else:
        module_name = f'fidXS_{suffix}{obs_name}_VH_{year}'
    
    # Import module dynamically
    th_xs = __import__(module_name, fromlist=[''])
    return th_xs
'''
def get_th_xsec(process, obs_gen, suffix, year, doubleDiff, channel, obs_gen_2nd=None):
    import sys, os

    # Add the inputs folder to sys.path so Python can find the modules
    inputs_dir = os.path.abspath("../inputs")
    if inputs_dir not in sys.path:
        sys.path.append(inputs_dir)

    # Map GEN variable names
    obs_name_dict = {
        "GENmass4l": "mass4l", "GENpT4l": "pT4l", "GENrapidity4l": "rapidity4l",
        "GENcostheta1": "costhetaZ1", "GENcostheta2": "costhetaZ2", "GENPhi": "phi",
        "GENPhi1": "phi1", "GENcosthetastar": "costhetastar",
        "GENmassZ1": "massZ1", "GENmassZ2": "massZ2", "GENpTj1": "pTj1",
        "GENpTj2": "pTj2", "GENmjj": "mjj", "GENabsdetajj": "absdetajj",
        "GENdphijj": "dphijj", "GENmHj":  "mHj", "GENpTHj": "pTHj",
        "GENpTHjj": "pTHjj", "GENNj": "Nj", "GENTBjMax": "TBjmax", "GENTCjMax": "TCjmax",
    }

    obs_name = obs_name_dict[obs_gen]

    chan_tag = f'_{channel}' if channel != '4l' else ''

    # Determine module name
    if doubleDiff:

        obs_name_2nd = obs_name_dict[obs_gen_2nd]

        if ('ZH' not in process) and ('W' not in process):
            pname = process.replace("NNLOPS_", "").split("125")[0]
            module_name = f'fidXS_{suffix}{obs_name}_{obs_name_2nd}_{pname}{chan_tag}_{year}'
        else:
            module_name = f'fidXS_{suffix}{obs_name}_{obs_name_2nd}_VH{chan_tag}_{year}'
    
    else:

        if ('ZH' not in process) and ('W' not in process):
            pname = process.replace("NNLOPS_", "").split("125")[0]
            if opt.ZZ: module_name = f'fidXS_{suffix}{obs_name}_zzfloating_{pname}{chan_tag}_{year}'
            else: module_name = f'fidXS_{suffix}{obs_name}_{pname}{chan_tag}_{year}'
        else:
            if opt.ZZ: module_name = f'fidXS_{suffix}{obs_name}_zzfloating_VH{chan_tag}_{year}'
            else: module_name = f'fidXS_{suffix}{obs_name}_VH{chan_tag}_{year}'
 
    print(module_name)

    # Import module dynamically
    th_xs = __import__(module_name, fromlist=[''])

    return th_xs

def get_unc_dict(obs_gen, process, year, tree, th_xs, channel, bins, nnlops, doubleDiff, obs_gen_2nd = None):

    unc_var_dn = {}
    unc_var_up = {}
    h_nom_var = {}
    h_scale_var = {}
    h_tot_var = {}
    h_tot_scale_var = {}

    for var in ["qcd", "pdf", "as"]: # SPENCER

        unc_bin_dn = {}
        unc_bin_up = {}
        h_nom_var[var] = {}
        h_scale_var[var] = {}
        h_tot_var[var] = {}
        h_tot_scale_var[var] = {}

        if doubleDiff:
            length = len(bins)
        else:
            length = len(bins)-1

        for idx in range(length):

            if doubleDiff:
                obs_gen_low = bins[idx][0]
                obs_gen_high = bins[idx][1]
                obs_gen_2nd_low = bins[idx][2]
                obs_gen_2nd_high = bins[idx][3]
                observable = [obs_gen, obs_gen_low, obs_gen_high, obs_gen_2nd, obs_gen_2nd_low, obs_gen_2nd_high]
            else:
                obs_gen_low = bins[idx]
                obs_gen_high = bins[idx+1]
                #if "mass4l" in obs_gen:
                #    obs_gen_low = 105
                #    obs_gen_high = 160
                observable = [obs_gen, obs_gen_low, obs_gen_high]

            h_nom, h_scale, h_tot, h_tot_scale = get_scale_unc(process, channel, tree, var, observable, nnlops)

            acc_nom = h_nom/h_tot

            scale_vars = {}
            for i in h_scale:
                acc_scale = h_scale[i]/h_tot_scale[i]
                scale_vars[i] = acc_scale/acc_nom
                
            if var == "pdf":
                scale_vars = np.array(list(scale_vars.values()))*acc_nom
                scale_unc = np.sqrt(np.sum(np.square(scale_vars - acc_nom)))
                unc_up = (acc_nom + scale_unc)/acc_nom
                unc_dn = (acc_nom - scale_unc)/acc_nom
            else:
                unc_up = max(scale_vars.values())
                unc_dn = min(scale_vars.values())
            
            up = 1 + (unc_up-1) + YR4_UNC[process][var]['up']
            # Correlate
            # np.sqrt((var_up-1)**2 + YR4_UNC[process][var]['up']**2)
            dn = 1 - ((1-unc_dn) + YR4_UNC[process][var]['dn'])
            # Correlate
            # np.sqrt((1-var_dn)**2 + YR4_UNC[process][var]['dn']**2)

            print(f'( {year} - {var}) {obs_gen}, bin_{idx} ({channel}): {up*th_xs.fidXS[idx]}/{dn*th_xs.fidXS[idx]}')

            unc_bin_up[idx] = up*th_xs.fidXS[idx]
            unc_bin_dn[idx] = dn*th_xs.fidXS[idx]

            h_nom_var[var][idx] = h_nom
            h_scale_var[var][idx] = h_scale
            h_tot_var[var][idx] = h_tot
            h_tot_scale_var[var][idx] = h_tot_scale
        
        unc_var_up[var] = unc_bin_up
        unc_var_dn[var] = unc_bin_dn


    return unc_var_up, unc_var_dn, h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var

def get_uncerainties(obs_gen, year, process, channel, bins, nnlops, doubleDiff, obs_gen_2nd = None):

    unc = {}
    h = {}

    fname = path['eos_path_sig']+year+"_MC/"+process+"/ZZ4lAnalysis_SKIMMED.root"

    if (process == "ggH125") and nnlops:
        suffix = "NNLOPS_"
    else:
        suffix = ""

    tree = load_tree(fname)
    
    th_xs = get_th_xsec(process, obs_gen, suffix, year, doubleDiff, channel, obs_gen_2nd)
    
    if doubleDiff:
        unc_var_up, unc_var_dn, h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var = get_unc_dict(obs_gen, process, year, tree, th_xs, channel, bins, nnlops, doubleDiff, obs_gen_2nd)
    else:
        unc_var_up, unc_var_dn, h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var = get_unc_dict(obs_gen, process, year, tree, th_xs, channel, bins, nnlops, doubleDiff)

    unc[process] = (unc_var_dn, unc_var_up)
    h[process] = (h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var)
    
    return unc, h

def save_uncertainties(process, obs_gen, nnlops, year, years, unc, h, channel, doubleDiff, obs_gen_2nd=None):

    pname = process.split("125")[0]

    if (pname == "ZH"): pname = "VH"

    if doubleDiff:
        obs_1 = opt.OBSNAME.split(' vs ')[0]
        obs_2 = opt.OBSNAME.split(' vs ')[1]
        obs_name = obs_1 + '_' + obs_2
    else:
        obs_name = opt.OBSNAME

    if (process == "ggH125") and nnlops:
        suffix = "NNLOPS_"
    else:
        suffix = ""


    th_xs = get_th_xsec(process, obs_gen, suffix, year, doubleDiff, channel, obs_gen_2nd)

    channel_tag = "" if channel == "4l" else f"_{channel}"

    fname = f'../inputs/fidXS_{suffix}{obs_name}_{pname}{channel_tag}_{opt.YEAR}.py'
    if opt.ZZ: fname = f'../inputs/fidXS_{suffix}{obs_name}_zzfloating_{pname}{channel_tag}_{opt.YEAR}.py'

    print(f'WRITING - {fname}')
         
    with open(fname, mode = 'w') as f:
        f.write(f'Boundaries = {bins}\n')
        f.write(f'fidXS = {th_xs.fidXS}\n')
        f.write(f'fidXS_scale_up = {list(unc[process][1]["qcd"].values())}\n')
        f.write(f'fidXS_scale_dn = {list(unc[process][0]["qcd"].values())}\n')
        f.write(f'fidXS_pdf_up = {list(unc[process][1]["pdf"].values())}\n')
        f.write(f'fidXS_pdf_dn = {list(unc[process][0]["pdf"].values())}\n')
        f.write(f'fidXS_alpha_up = {list(unc[process][1]["as"].values())}\n')
        f.write(f'fidXS_alpha_dn = {list(unc[process][0]["as"].values())}\n')

        # save these to compute ratios for merging eras
        f.write(f'h_nom_qcd = {h[process][0]["qcd"]}\n')
        f.write(f'h_scale_qcd = {list(h[process][1]["qcd"].values())}\n')
        f.write(f'h_tot_qcd = {h[process][2]["qcd"]}\n')
        f.write(f'h_tot_scale_qcd = {list(h[process][3]["qcd"].values())}\n')
        f.write(f'h_nom_pdf = {h[process][0]["pdf"]}\n')
        f.write(f'h_scale_pdf = {list(h[process][1]["pdf"].values())}\n')
        f.write(f'h_tot_pdf = {h[process][2]["pdf"]}\n')
        f.write(f'h_tot_scale_pdf = {list(h[process][3]["pdf"].values())}\n')
        f.write(f'h_nom_as = {h[process][0]["as"]}\n')
        f.write(f'h_scale_as = {list(h[process][1]["as"].values())}\n')
        f.write(f'h_tot_as = {h[process][2]["as"]}\n')
        f.write(f'h_tot_scale_as = {list(h[process][3]["as"].values())}\n')

YR4_UNC = {'ggH125': {'qcd': {'up': 0.077, 'dn': 0.088},
                     'pdf': {'up': 0.018, 'dn': 0.018},
                     'as':{'up': 0.025, 'dn': 0.025}},
           
           'VBFH125': {'qcd': {'up': 0.004, 'dn': 0.003},
                     'pdf': {'up': 0.021, 'dn': 0.021},
                     'as':{'up': 0.005, 'dn': 0.005}},
           
           'ZH125': {'qcd': {'up': 0.038, 'dn': 0.031},
                     'pdf': {'up': 0.013, 'dn': 0.013},
                     'as':{'up': 0.009, 'dn': 0.009}},
           
           'ttH125': {'qcd': {'up': 0.058, 'dn': 0.092},
                     'pdf': {'up': 0.03, 'dn': 0.03},
                     'as':{'up': 0.02, 'dn': 0.02}},
          }

def _write_combined_4e4mu_module(process, obs_name, suffix, out_year, era_to_read, bins, doubleDiff):
    """
    Build fidXS_*_{process}_4e4mu_<out_year>.py by combining the already-written
    4e and 4mu modules for the same era_to_read (e.g. 2022, 2022EE, ... or merged era).

    It recomputes qcd/pdf/as up/down the same way as your MERGE block:
      - QCD/as: envelope of acceptance variations
      - PDF: RMS of acceptance variations
      - then adds correlated YR4_UNC and multiplies by fidXS
    """

    # ---- figure out pname used in filenames ----
    # save_uncertainties uses: pname = process.split("125")[0], except ZH->VH
    pname = process.split("125")[0]
    if pname == "ZH":
        pname = "VH"

    # ---- module names to load (4e and 4mu) ----
    # matches your save_uncertainties naming:
    # ../inputs/fidXS_{suffix}{obs_name}_{pname}_{channel}_{era}.py
    mod4e_fname  = f'../inputs/fidXS_{suffix}{obs_name}_{pname}_4e_{era_to_read}.py'
    mod4mu_fname = f'../inputs/fidXS_{suffix}{obs_name}_{pname}_4mu_{era_to_read}.py'

    print(mod4e_fname)
    print(mod4mu_fname)

    if os.path.exists(mod4e_fname):

        module_name_4e = f'fidXS_{suffix}{obs_name}_{pname}_4e_{era_to_read}'
        spec_4e = importlib.util.spec_from_file_location(module_name_4e, mod4e_fname)
        m4e = importlib.util.module_from_spec(spec_4e)
        spec_4e.loader.exec_module(m4e)

    if os.path.exists(mod4mu_fname):

        module_name_4mu = f'fidXS_{suffix}{obs_name}_{pname}_4mu_{era_to_read}'
        spec_4mu = importlib.util.spec_from_file_location(module_name_4mu, mod4mu_fname)
        m4mu = importlib.util.module_from_spec(spec_4mu)
        spec_4mu.loader.exec_module(m4mu)

    # ---- helpers ----
    def sum_dict(a, b):
        return {k: a[k] + b[k] for k in a}

    def sum_scale_by_bin(scale4e, scale4mu):
        # scale4e, scale4mu: list[dict] with same length (one dict per bin)
        if len(scale4e) != len(scale4mu):
            raise ValueError(f"Different number of bins: 4e={len(scale4e)} 4mu={len(scale4mu)}")

        out = []
        for d4e, d4mu in zip(scale4e, scale4mu):
            # make a new dict with summed values per key
            dsum = dict(d4e)  # copy
            for k, v in d4mu.items():
                dsum[k] = dsum.get(k, 0.0) + v
            out.append(dsum)

        return out

    # determine number of bins
    if doubleDiff:
        length = len(bins)
    else:
        length = len(bins) - 1

    # ---- combine fidXS (the nominal theory cross section) ----
    fidXS = (np.array(m4e.fidXS) + np.array(m4mu.fidXS)).tolist()

    # ---- combine the stored hist-integral pieces ----
    h_nom_qcd = sum_dict(m4e.h_nom_qcd, m4mu.h_nom_qcd)
    h_tot_qcd = sum_dict(m4e.h_tot_qcd, m4mu.h_tot_qcd)
    h_scale_qcd = sum_scale_by_bin(m4e.h_scale_qcd, m4mu.h_scale_qcd)
    h_tot_scale_qcd = sum_scale_by_bin(m4e.h_tot_scale_qcd, m4mu.h_tot_scale_qcd)

    h_nom_as = sum_dict(m4e.h_nom_as, m4mu.h_nom_as)
    h_tot_as = sum_dict(m4e.h_tot_as, m4mu.h_tot_as)
    h_scale_as = sum_scale_by_bin(m4e.h_scale_as, m4mu.h_scale_as)
    h_tot_scale_as = sum_scale_by_bin(m4e.h_tot_scale_as, m4mu.h_tot_scale_as)

    h_nom_pdf = sum_dict(m4e.h_nom_pdf, m4mu.h_nom_pdf)
    h_tot_pdf = sum_dict(m4e.h_tot_pdf, m4mu.h_tot_pdf)
    h_scale_pdf = sum_scale_by_bin(m4e.h_scale_pdf, m4mu.h_scale_pdf)
    h_tot_scale_pdf = sum_scale_by_bin(m4e.h_tot_scale_pdf, m4mu.h_tot_scale_pdf)

    # ---- recompute acceptance variation factors (like your MERGE block) ----
    up_qcd_TOT, dn_qcd_TOT = {}, {}
    up_as_TOT,  dn_as_TOT  = {}, {}
    up_pdf_TOT, dn_pdf_TOT = {}, {}

    for idx in range(length):
        # QCD
        acc_nom_qcd = h_nom_qcd[idx] / h_tot_qcd[idx]
        acc_scale_qcd = np.array(list(h_scale_qcd[idx].values())) / np.array(list(h_tot_scale_qcd[idx].values()))
        scale_vars_qcd = acc_scale_qcd / acc_nom_qcd
        up_qcd_TOT[idx] = float(np.max(scale_vars_qcd))
        dn_qcd_TOT[idx] = float(np.min(scale_vars_qcd))

        # alpha_s
        acc_nom_as = h_nom_as[idx] / h_tot_as[idx]
        acc_scale_as = np.array(list(h_scale_as[idx].values())) / np.array(list(h_tot_scale_as[idx].values()))
        scale_vars_as = acc_scale_as / acc_nom_as
        up_as_TOT[idx] = float(np.max(scale_vars_as))
        dn_as_TOT[idx] = float(np.min(scale_vars_as))

    for idx in range(length):
        # PDF (RMS)
        acc_nom_pdf = h_nom_pdf[idx] / h_tot_pdf[idx]
        acc_scale_pdf = np.array(list(h_scale_pdf[idx].values())) / np.array(list(h_tot_scale_pdf[idx].values()))
        scale_unc_pdf = np.sqrt(np.sum(np.square(acc_scale_pdf - acc_nom_pdf)))
        up_pdf_TOT[idx] = float((acc_nom_pdf + scale_unc_pdf) / acc_nom_pdf)
        dn_pdf_TOT[idx] = float((acc_nom_pdf - scale_unc_pdf) / acc_nom_pdf)

    # process key for YR4_UNC dict
    # your YR4_UNC keys are ggH125, VBFH125, ttH125, ZH125
    p_for_unc = process
    if p_for_unc.startswith("NNLOPS_"):
        p_for_unc = p_for_unc.replace("NNLOPS_", "")

    # ---- turn factors into absolute fidXS_up/dn arrays ----
    fidXS_scale_up, fidXS_scale_dn = [], []
    fidXS_pdf_up,   fidXS_pdf_dn   = [], []
    fidXS_alpha_up, fidXS_alpha_dn = [], []

    for idx in range(length):
        # PDF
        up_pdf = 1 + (up_pdf_TOT[idx] - 1) + YR4_UNC[p_for_unc]['pdf']['up']
        dn_pdf = 1 - ((1 - dn_pdf_TOT[idx]) + YR4_UNC[p_for_unc]['pdf']['dn'])

        # QCD
        up_qcd = 1 + (up_qcd_TOT[idx] - 1) + YR4_UNC[p_for_unc]['qcd']['up']
        dn_qcd = 1 - ((1 - dn_qcd_TOT[idx]) + YR4_UNC[p_for_unc]['qcd']['dn'])

        # alpha_s
        up_as = 1 + (up_as_TOT[idx] - 1) + YR4_UNC[p_for_unc]['as']['up']
        dn_as = 1 - ((1 - dn_as_TOT[idx]) + YR4_UNC[p_for_unc]['as']['dn'])

        fidXS_scale_up.append(up_qcd * fidXS[idx])
        fidXS_scale_dn.append(dn_qcd * fidXS[idx])
        fidXS_pdf_up.append(up_pdf * fidXS[idx])
        fidXS_pdf_dn.append(dn_pdf * fidXS[idx])
        fidXS_alpha_up.append(up_as * fidXS[idx])
        fidXS_alpha_dn.append(dn_as * fidXS[idx])

    # ---- write output module ----
    outname = f'../inputs/fidXS_{suffix}{obs_name}_{pname}_4e4mu_{out_year}.py'
    print(f'WRITING COMBINED 4e4mu - {outname}')

    with open(outname, 'w') as f:
        f.write(f'Boundaries = {bins}\n')
        f.write(f'fidXS = {fidXS}\n')
        f.write(f'fidXS_scale_up = {fidXS_scale_up}\n')
        f.write(f'fidXS_scale_dn = {fidXS_scale_dn}\n')
        f.write(f'fidXS_pdf_up = {fidXS_pdf_up}\n')
        f.write(f'fidXS_pdf_dn = {fidXS_pdf_dn}\n')
        f.write(f'fidXS_alpha_up = {fidXS_alpha_up}\n')
        f.write(f'fidXS_alpha_dn = {fidXS_alpha_dn}\n')

        # keep the merging payload too (same as other channels)
        f.write(f'h_nom_qcd = {h_nom_qcd}\n')
        f.write(f'h_scale_qcd = {h_scale_qcd}\n')
        f.write(f'h_tot_qcd = {h_tot_qcd}\n')
        f.write(f'h_tot_scale_qcd = {h_tot_scale_qcd}\n')

        f.write(f'h_nom_pdf = {h_nom_pdf}\n')
        f.write(f'h_scale_pdf = {h_scale_pdf}\n')
        f.write(f'h_tot_pdf = {h_tot_pdf}\n')
        f.write(f'h_tot_scale_pdf = {h_tot_scale_pdf}\n')

        f.write(f'h_nom_as = {h_nom_as}\n')
        f.write(f'h_scale_as = {h_scale_as}\n')
        f.write(f'h_tot_as = {h_tot_as}\n')
        f.write(f'h_tot_scale_as = {h_tot_scale_as}\n')

if __name__ == '__main__':

    SPECIAL_OBS = { } #{'massZ1','massZ2','costhetaZ1','costhetaZ2','costhetastar','phi','phi1'}

    # spencer

    m4l_low = 105
    m4l_high = 160
    #year = opt.YEAR

    if "zzfloating" in opt.OBSNAME:
        OBSNAME = opt.OBSNAME.rsplit('_', 1)[0]
    else:
        OBSNAME = opt.OBSNAME

    bins, doubleDiff = binning(OBSNAME)
    if doubleDiff:
        obs_name = OBSNAME.split(' vs ')[0]
        obs_name_2nd = OBSNAME.split(' vs ')[1]
        obs_name_2d = OBSNAME
    else:
        obs_name = OBSNAME

    
    _temp = __import__('observables', globals(), locals(), ['observables'], 0) # spencer
    observables = _temp.observables

    if doubleDiff:
        obs_gen = observables[obs_name_2d]['obs_gen']
        obs_gen_2nd = observables[obs_name_2d]['obs_gen_2nd']
        length = len(bins)
    else:
        obs_gen = observables[obs_name]['obs_gen']
        obs_gen_2nd = None
        length = len(bins)-1

    if obs_name == "TCjMax": 
        obs_name = "TCjmax"
        obs_gen = "GENTCjmax"
    if obs_name == "TBjMax": 
        obs_name = "TBjmax"
        obs_gen = "GENTBjmax"
    if doubleDiff:
        if obs_name_2nd == "TCjMax": 
            obs_name_2nd = "TCjmax"
            obs_gen_2nd = "GENTCjmax"
        if obs_name_2nd == "TBjMax": 
            obs_name_2nd = "TBjmax"
            obs_gen_2nd = "GENTBjmax"

    channels = ['4l']

    if opt.OBSNAME == 'mass4l' or opt.OBSNAME == 'mass4l_zzfloating':
        channels = ['2e2mu', '4e', '4mu', '4l']
    elif opt.OBSNAME in SPECIAL_OBS:
        # we must compute 4e and 4mu so we can build 4e4mu later
        channels = ['2e2mu', '4e', '4mu', '4l']

    if (opt.YEAR == '2022'): years = ['2022']
    if (opt.YEAR == '2022EE'): years = ['2022EE']
    if (opt.YEAR == '2023preBPix'): years = ['2023preBPix']
    if (opt.YEAR == '2023postBPix'): years = ['2023postBPix']
    if (opt.YEAR == '2024'): years = ['2024']

    if (opt.YEAR == '2022full'): years = ['2022', '2022EE']
    if (opt.YEAR == '2023full'): years = ['2023preBPix', '2023postBPix']
    if (opt.YEAR == '2022_2023'): years = ['2022', '2022EE', '2023preBPix', '2023postBPix']
    if (opt.YEAR == 'Run3'): years = ['2022', '2022EE', '2023preBPix', '2023postBPix', '2024']


    if opt.SPLIT:
        YEAR_SPLIT = opt.YEAR.split("_")[0]
        nSplit = split[YEAR_SPLIT]
        years = [opt.YEAR]
    else:
        nSplit = 0

    if opt.MERGE:
        if nSplit != 0:
            years = []
            for i in range(0,nSplit):
                years.append(opt.YEAR+"_"+str(i+1))
        else:
            years = years


    print(years)


    if opt.MERGE:

        for channel in channels:

            channel_tag = "" if channel == "4l" else f"_{channel}"

            unc = {}
            h = {}

            processes = ["NNLOPS_ggH125", "ggH125", "VBFH125", "ttH125", "ZH125"]

            for process in processes:

                # --- initialize per-bin totals ---
                h_nom_qcd_TOT = {idx: 0.0 for idx in range(length)}
                h_tot_qcd_TOT = {idx: 0.0 for idx in range(length)}
                h_scale_qcd_TOT = {idx: [] for idx in range(length)}
                h_tot_scale_qcd_TOT = {idx: [] for idx in range(length)}

                h_nom_as_TOT = {idx: 0.0 for idx in range(length)}
                h_tot_as_TOT = {idx: 0.0 for idx in range(length)}
                h_scale_as_TOT = {idx: [] for idx in range(length)}
                h_tot_scale_as_TOT = {idx: [] for idx in range(length)}

                h_nom_pdf_TOT = {idx: 0.0 for idx in range(length)}
                h_tot_pdf_TOT = {idx: 0.0 for idx in range(length)}
                h_scale_pdf_TOT = {idx: [] for idx in range(length)}
                h_tot_scale_pdf_TOT = {idx: [] for idx in range(length)}

                suffix = "NNLOPS_" if process.startswith("NNLOPS_") else ""
                if process.startswith("NNLOPS_"):
                    pname = process.replace("NNLOPS_", "").split("125")[0]
                else:
                    pname = process.split("125")[0]

                print(f'Processing {process}...')

                if doubleDiff:
                    if opt.ZZ: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{obs_name_2nd}_zzfloating_{pname}{channel_tag}_{opt.YEAR}.py'
                    else: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{obs_name_2nd}_{pname}{channel_tag}_{opt.YEAR}.py'
                else:
                    if opt.ZZ: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_zzfloating_{pname}{channel_tag}_{opt.YEAR}.py'
                    else: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{pname}{channel_tag}_{opt.YEAR}.py'


                th_xs = get_th_xsec(pname, obs_gen, suffix, opt.YEAR, doubleDiff, channel, obs_gen_2nd)

                #if os.path.exists(fname_ALL):
                #    th_xs = get_th_xsec(pname, obs_gen, suffix, opt.YEAR, doubleDiff, obs_gen_2nd)
                #else:
                #    print(f'{fname_ALL} does not exist!')
                #    continue

                for year in years:

                    if (pname == "ZH"): pname = "VH"

                    if doubleDiff:
                        if opt.ZZ: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{obs_name_2nd}_zzfloating_{pname}{channel_tag}_{year}.py'
                        else: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{obs_name_2nd}_{pname}{channel_tag}_{year}.py'
                    else:
                        if opt.ZZ: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_zzfloating_{pname}{channel_tag}_{year}.py'
                        else: fname_ALL = f'../inputs/fidXS_{suffix}{obs_name}_{pname}{channel_tag}_{year}.py'

                    if os.path.exists(fname_ALL):

                        module_name = f'fidXS_{suffix}{obs_name}_{pname}{channel_tag}_{year}'
                        if opt.ZZ: module_name = f'fidXS_{suffix}{obs_name}_zzfloating_{pname}{channel_tag}_{year}'

                        print(module_name)
                        spec = importlib.util.spec_from_file_location(module_name, fname_ALL)
                        module = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(module)

                        for idx in range(length):

                            # QCD
                            print(h_nom_qcd_TOT)
                            print(module.h_nom_qcd)
                            h_nom_qcd_TOT[idx] += module.h_nom_qcd[idx]
                            h_tot_qcd_TOT[idx] += module.h_tot_qcd[idx]

                            if not h_scale_qcd_TOT[idx]:
                                h_scale_qcd_TOT[idx] = module.h_scale_qcd[idx].copy()
                                h_tot_scale_qcd_TOT[idx] = module.h_tot_scale_qcd[idx].copy()
                            else:
                                # sum per variation (dict value by key)
                                for key, val in module.h_scale_qcd[idx].items():
                                    h_scale_qcd_TOT[idx][key] += val
                                for key, val in module.h_tot_scale_qcd[idx].items():
                                    h_tot_scale_qcd_TOT[idx][key] += val

                            # alpha_s
                            h_nom_as_TOT[idx] += module.h_nom_as[idx]
                            h_tot_as_TOT[idx] += module.h_tot_as[idx]

                            if not h_scale_as_TOT[idx]:
                                h_scale_as_TOT[idx] = module.h_scale_as[idx].copy()
                                h_tot_scale_as_TOT[idx] = module.h_tot_scale_as[idx].copy()
                            else:
                                # sum per variation (dict value by key)
                                for key, val in module.h_scale_as[idx].items():
                                    h_scale_as_TOT[idx][key] += val
                                for key, val in module.h_tot_scale_as[idx].items():
                                    h_tot_scale_as_TOT[idx][key] += val

                            # PDF
                            h_nom_pdf_TOT[idx] += module.h_nom_pdf[idx]
                            h_tot_pdf_TOT[idx] += module.h_tot_pdf[idx]
                            
                            if not h_scale_pdf_TOT[idx]:
                                h_scale_pdf_TOT[idx] = module.h_scale_pdf[idx].copy()
                                h_tot_scale_pdf_TOT[idx] = module.h_tot_scale_pdf[idx].copy()
                            else:
                                # sum per variation (dict value by key)
                                for key, val in module.h_scale_pdf[idx].items():
                                    h_scale_pdf_TOT[idx][key] += val
                                for key, val in module.h_tot_scale_pdf[idx].items():
                                    h_tot_scale_pdf_TOT[idx][key] += val

                    else:
                        print(f'{fname_ALL}' + ' does not exist!')
                        continue
                

                up_qcd_TOT, dn_qcd_TOT, dn_pdf_TOT = {}, {}, {}
                up_as_TOT, dn_as_TOT, up_pdf_TOT = {}, {}, {}

                for idx in range(length):

                    acc_nom_qcd_TOT = h_nom_qcd_TOT[idx] / h_tot_qcd_TOT[idx]
                    acc_nom_as_TOT  = h_nom_as_TOT[idx] / h_tot_as_TOT[idx]

                    acc_scale_qcd_TOT = np.array(list(h_scale_qcd_TOT[idx].values()))/np.array(list(h_tot_scale_qcd_TOT[idx].values()))
                    acc_scale_as_TOT =np.array(list(h_scale_as_TOT[idx].values()))/np.array(list(h_tot_scale_as_TOT[idx].values()))

                    scale_vars_qcd_TOT = acc_scale_qcd_TOT/acc_nom_qcd_TOT
                    scale_vars_as_TOT = acc_scale_as_TOT/acc_nom_as_TOT

                    up_qcd_TOT[idx] = max(scale_vars_qcd_TOT)
                    dn_qcd_TOT[idx] = min(scale_vars_qcd_TOT)

                    up_as_TOT[idx] = max(scale_vars_as_TOT)
                    dn_as_TOT[idx] = min(scale_vars_as_TOT)
            
                
        
                for idx in range(length):

                    acc_nom_pdf_TOT = h_nom_pdf_TOT[idx] / h_tot_pdf_TOT[idx]
                    acc_scale_pdf_TOT = np.array(list(h_scale_pdf_TOT[idx].values()))/np.array(list(h_tot_scale_pdf_TOT[idx].values()))
                    scale_unc_pdf_TOT = np.sqrt(np.sum(np.square(acc_scale_pdf_TOT -  acc_nom_pdf_TOT)))
                    up_pdf_TOT[idx] = (acc_nom_pdf_TOT + scale_unc_pdf_TOT)/acc_nom_pdf_TOT
                    dn_pdf_TOT[idx] = (acc_nom_pdf_TOT - scale_unc_pdf_TOT)/acc_nom_pdf_TOT

                '''
                print(f'--- Final results for {process} (last bin) ---')
                print('qcd all', scale_vars_qcd_TOT)
                print('qcd max', max(scale_vars_qcd_TOT))
                print('qcd min', min(scale_vars_qcd_TOT))
                print('as all', scale_vars_as_TOT)
                print('as max', max(scale_vars_as_TOT))
                print('as min', min(scale_vars_as_TOT))
                print('pdf all', acc_scale_pdf_TOT)
                print('pdf up', up_pdf_TOT)
                print('pdf dn', dn_pdf_TOT)
                print('-------------------------------')
                '''

                print('qcd up:', up_qcd_TOT)
                print('qcd dn:', dn_qcd_TOT)
                print('as up:', up_as_TOT)
                print('as dn:', dn_as_TOT) 
                print('pdf up:', up_pdf_TOT)
                print('pdf dn:', dn_pdf_TOT)

                unc_bin_up_pdf, unc_bin_dn_pdf = {}, {}
                unc_bin_up_qcd, unc_bin_dn_qcd = {}, {}
                unc_bin_up_as, unc_bin_dn_as = {}, {}

                for idx in range(length):

                    if process.startswith("NNLOPS_"):
                        p = process.replace("NNLOPS_", "")
                    else:
                        p = process

                    # PDF: symmetric RMS + correlated uncertainty
                    up_pdf = 1 + (up_pdf_TOT[idx]-1) + YR4_UNC[p]['pdf']['up']
                    dn_pdf = 1 - ((1 - dn_pdf_TOT[idx]) + YR4_UNC[p]['pdf']['dn'])
                    unc_bin_up_pdf[idx] = up_pdf * th_xs.fidXS[idx]
                    unc_bin_dn_pdf[idx] = dn_pdf * th_xs.fidXS[idx]

                    # QCD: envelope + correlated uncertainty
                    up_qcd = 1 + (up_qcd_TOT[idx]-1) + YR4_UNC[p]['qcd']['up']
                    dn_qcd = 1 - ((1 - dn_qcd_TOT[idx]) + YR4_UNC[p]['qcd']['dn'])
                    unc_bin_up_qcd[idx] = up_qcd * th_xs.fidXS[idx]
                    unc_bin_dn_qcd[idx] = dn_qcd * th_xs.fidXS[idx]

                    # alpha_s: envelope + correlated uncertainty
                    up_as = 1 + (up_as_TOT[idx]-1) + YR4_UNC[p]['as']['up']
                    dn_as = 1 - ((1 - dn_as_TOT[idx]) + YR4_UNC[p]['as']['dn'])
                    unc_bin_up_as[idx] = up_as * th_xs.fidXS[idx]
                    unc_bin_dn_as[idx] = dn_as * th_xs.fidXS[idx]

                print(up_qcd)
                print(dn_qcd)

                unc_var_up, unc_var_dn = {}, {}
                h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var = {}, {}, {}, {}

                unc_var_up['pdf'] = unc_bin_up_pdf
                unc_var_dn['pdf'] = unc_bin_dn_pdf
                unc_var_up['qcd'] = unc_bin_up_qcd
                unc_var_dn['qcd'] = unc_bin_dn_qcd
                unc_var_up['as'] = unc_bin_up_as
                unc_var_dn['as'] = unc_bin_dn_as

                h_nom_var['pdf'] = h_nom_pdf_TOT
                h_scale_var['pdf'] = h_scale_pdf_TOT
                h_tot_var['pdf'] = h_tot_pdf_TOT
                h_tot_scale_var['pdf'] = h_tot_scale_pdf_TOT
                h_nom_var['qcd'] = h_nom_qcd_TOT
                h_scale_var['qcd'] = h_scale_qcd_TOT
                h_tot_var['qcd'] = h_tot_qcd_TOT
                h_tot_scale_var['qcd'] = h_tot_scale_qcd_TOT
                h_nom_var['as'] = h_nom_as_TOT
                h_scale_var['as'] = h_scale_as_TOT
                h_tot_var['as'] = h_tot_as_TOT
                h_tot_scale_var['as'] = h_tot_scale_as_TOT

                unc[p] = (unc_var_dn, unc_var_up)
                h[p] = (h_nom_var, h_scale_var, h_tot_var, h_tot_scale_var)

                if doubleDiff:
                    if process.startswith("NNLOPS_"):
                        save_uncertainties("ggH125", obs_gen, True,  opt.YEAR, years, unc, h, channel, True,  obs_gen_2nd)
                    else:
                        save_uncertainties(process,  obs_gen, False, opt.YEAR, years, unc, h, channel, True,  obs_gen_2nd)
                else:
                    if process.startswith("NNLOPS_"):
                        save_uncertainties("ggH125", obs_gen, True,  opt.YEAR, years, unc, h, channel, False)
                    else:
                        save_uncertainties(process,  obs_gen, False, opt.YEAR, years, unc, h, channel, False)

                # NEW: after producing merged 4e and 4mu modules, make 4e4mu
                if (opt.OBSNAME in SPECIAL_OBS) and (opt.OBSNAME != 'mass4l' and opt.OBSNAME != 'mass4l_zzfloating') and (channel == '4mu'):

                    # process name + suffix must match what save_uncertainties() wrote
                    if process.startswith("NNLOPS_"):
                        proc_nom = process.replace("NNLOPS_", "")   # "ggH125"
                        suffix = "NNLOPS_"
                    else:
                        proc_nom = process                           # e.g. "VBFH125"
                        suffix = ""

                    obs_tag = opt.OBSNAME if not doubleDiff else f"{obs_name}_{obs_name_2nd}"

                    _write_combined_4e4mu_module(
                        process=proc_nom,
                        obs_name=obs_tag,
                        suffix=suffix,
                        out_year=opt.YEAR,     # merged tag like 2022full/Run3
                        era_to_read=opt.YEAR,  # read the merged 4e + 4mu modules you just wrote
                        bins=bins,
                        doubleDiff=doubleDiff
                    )

    else:

        processes = ["NNLOPS_ggH125", "ggH125", "VBFH125", "ttH125", "ZH125"]

        for process in processes:

            use_NNLOPS_flag = process.startswith("NNLOPS_")
            proc_nom = process.replace("NNLOPS_", "") if use_NNLOPS_flag else process

            for year in years:

                # ---- run the needed channels ----
                if opt.OBSNAME == 'mass4l' or opt.OBSNAME == 'mass4l_zzfloating':
                    run_channels = ['2e2mu','4e','4mu','4l']  # unchanged
                elif opt.OBSNAME in SPECIAL_OBS:
                    run_channels = ['2e2mu','4e','4mu','4l']  # needed to build 4e4mu
                else:
                    run_channels = ['4l']  # unchanged

                # compute+write each channel as before
                for channel in run_channels:
                    if doubleDiff:
                        unc, h = get_uncerainties(obs_gen, year, proc_nom, channel, bins, use_NNLOPS_flag, True, obs_gen_2nd)
                        save_uncertainties(proc_nom, obs_gen, use_NNLOPS_flag, opt.YEAR, year, unc, h, channel, True, obs_gen_2nd)
                    else:
                        unc, h = get_uncerainties(obs_gen, year, proc_nom, channel, bins, use_NNLOPS_flag, False)
                        save_uncertainties(proc_nom, obs_gen, use_NNLOPS_flag, opt.YEAR, year, unc, h, channel, False)

                # ---- NEW: build combined 4e4mu only for special obs (NOT mass4l) ----
                if (opt.OBSNAME in SPECIAL_OBS) and (opt.OBSNAME != 'mass4l' and opt.OBSNAME != 'mass4l_zzfloating'):
                    suffix = "NNLOPS_" if (proc_nom == "ggH125" and use_NNLOPS_flag) else ""
                    _write_combined_4e4mu_module(
                        process=proc_nom,
                        obs_name=obs_name if not doubleDiff else f"{obs_name}_{obs_name_2nd}",
                        suffix=suffix,
                        out_year=year,          # this is the era module name you just wrote (e.g. 2022)
                        era_to_read=year,       # read 4e and 4mu for the same era
                        bins=bins,
                        doubleDiff=doubleDiff
                    )

            
    # spencer
            
    '''
    m4l_low = 105
    m4l_high = 160

    obs_gen = 'GENmass4l'
    bins = [0, 2, 4, 6, 8]

    channel = '4l'
    process = 'ggH125'

    for process in ["ggH125"]:#, "VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, True)
        save_uncertainties(process, obs_gen, True)

    obs_gen = 'GENrapidity4l'
    bins = [0.0, 0.15, 0.30, 0.6, 0.9, 2.5]

    for process in ["ggH125", "VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, False)
        save_uncertainties(process, obs_gen, False)
        
    for process in ["ggH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, True)
        save_uncertainties(process, obs_gen, True)

    obs_gen = 'GENpT4l'
    bins = [0, 30, 80, 10000]
    for process in ["VBFH125", "ttH125", "ZH125"]:
        unc = get_uncerainties(obs_gen, process, channel, bins, False)
        save_uncertainties(process, obs_gen, False)
    '''
