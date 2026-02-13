import os, sys
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import json
from argparse import ArgumentParser

sys.path.append('helperstuff/')
from paths import path

############################## mass4l helper functions ##############################

def _mass4l_json_paths(config_json_path: str, era: str):
    """
    Given the inclusive mass4l json path (the one with NO final state in the name),
    build the 4 paths in the order: 4l, 2e2mu, 4mu, 4e.

    Example:
      mass4l_results_2022.json
        -> mass4l_results_2e2mu_2022.json
        -> mass4l_results_4mu_2022.json
        -> mass4l_results_4e_2022.json
    """
    d = os.path.dirname(config_json_path)
    base = os.path.basename(config_json_path)

    suffix = f"_{era}.json"
    if not base.endswith(suffix):
        # fallback: just append _<fs> before .json
        stem, ext = os.path.splitext(base)
        return [
            os.path.join(d, base),                    # 4l inclusive
            os.path.join(d, f"{stem}_2e2mu{ext}"),
            os.path.join(d, f"{stem}_4mu{ext}"),
            os.path.join(d, f"{stem}_4e{ext}"),
        ]

    stem = base[:-len(suffix)]  # e.g. "mass4l_results"
    return [
        os.path.join(d, base),                       # 4l inclusive
        os.path.join(d, f"{stem}_2e2mu{suffix}"),
        os.path.join(d, f"{stem}_4mu{suffix}"),
        os.path.join(d, f"{stem}_4e{suffix}"),
    ]


def load_mass4l_final_states(inclusive_mass4l_json_path: str, era: str):
    """
    Read 4 JSON files (4l, 2e2mu, 4mu, 4e) and return:
      - a combined config dict where exp_xs/uncs are arrays length 4
      - a list of per-final-state configs (each one is the mass4l block)
      - categorical binning appropriate for final states (0..4 with centers 0.5..3.5)
    """
    paths = _mass4l_json_paths(inclusive_mass4l_json_path, era)
    final_state_labels = ["4l", "2e2mu", "4mu", "4e"]

    per_fs = []
    for p in paths:
        with open(p, "r") as f:
            js = json.load(f)
        if "mass4l" not in js:
            raise KeyError(f"'mass4l' key not found in {p}")
        per_fs.append(js["mass4l"])

    # Build arrays of length 4 from the single-entry lists in each JSON
    def one(x):  # each x is like [3.02]
        return float(x[0])

    combined = dict(per_fs[0])  # start from inclusive as the template

    combined["exp_xs"]    = [one(cfg["exp_xs"])    for cfg in per_fs]
    combined["err_up"]    = [one(cfg["err_up"])    for cfg in per_fs]
    combined["err_down"]  = [one(cfg["err_down"])  for cfg in per_fs]
    combined["stat_up"]   = [one(cfg["stat_up"])   for cfg in per_fs]
    combined["stat_down"] = [one(cfg["stat_down"]) for cfg in per_fs]

    # Keep module names per final state so we can import theory bin-by-bin
    combined["_per_fs_configs"] = per_fs
    combined["_final_state_labels"] = final_state_labels

    # Categorical bins: 4 categories -> edges [0,1,2,3,4], centers [0.5,1.5,2.5,3.5]
    bins_plot = np.arange(0, 5, 1, dtype=float)
    bins_c = (bins_plot[1:] + bins_plot[:-1]) * 0.5
    bin_w = np.ones(4, dtype=float)  # IMPORTANT: do NOT divide by 25 for categories

    return combined, bins_plot, bins_c, bin_w


#####################################################################################

#SPECIAL_OBS = {'massZ1','massZ2','costhetaZ1','costhetaZ2','costhetastar','phi','phi1'}
SPECIAL_OBS = {}

def load_var_config(json_path: str, variable: str):
    with open(json_path, "r") as f:
        js = json.load(f)
    if variable not in js:
        raise KeyError(f"'{variable}' key not found in {json_path}")
    return js[variable]

def _specialobs_json_paths(variable: str, era: str):
    # assumes your per-variable jsons live in jsons/
    return [
        f"jsons/{variable}_results_{era}.json",
        f"jsons/{variable}_results_2e2mu_{era}.json",
        f"jsons/{variable}_results_4e4mu_{era}.json",
    ]

#####################################################################################

# Use CMS style from mplhep for plotting
plt.style.use(hep.style.CMS)

# Append custom directory 'fidXS' to the system path for importing modules
sys.path.append(os.path.join(os.path.dirname(__file__), 'fidXS'))

# Function to parse command-line arguments
def get_options():
    parser = ArgumentParser()
    parser.add_argument("--variables", dest="variables", required=True, default='', help="Variable(s) to be plotted. If multiple variables are asked, separate them with commas (e.g. PTH,rapidity,Njets2p5). One can also plot all available variables.")
    parser.add_argument("--config-json", dest="config_json", required=False, default='./config_unblinded.json', help="Path to the config JSON file which contains the information for the plot (LL values and uncertainties, plotting informations, etc.)")
    parser.add_argument("--no-preliminary", dest="no_preliminary", required=False, default=False, action='store_true', help="Flag that determines, if Preliminary label should be printed.")
    parser.add_argument("--year", dest="YEAR", required=True, default='', help="Flag that determines, if Preliminary label should be printed.")

    return parser.parse_args()
args = get_options()

# Load configuration from JSON file
with open(args.config_json, 'r') as f:
    config = json.load(f)

# Get list of variables to be plotted
variable_list = args.variables.split(",")

# If 'all' is specified, plot all available variables from the config
if 'all' in variable_list:
    variable_list = list(config.keys())

# Loop through each variable in the variable list
for variable in variable_list:

    # build a list of (label, current_config, out_suffix) to plot
    plot_jobs = []

    if variable == "mass4l":
        current_config, bins_plot, bins_c, bin_w = load_mass4l_final_states(args.config_json, args.YEAR)
        plot_jobs = [("mass4l", current_config, "")]  # keep your mass4l special treatment

    elif variable in SPECIAL_OBS:
        # derive the 3 json paths from the provided 4l json
        paths = _specialobs_json_paths(variable, args.YEAR)
        labels = ["4l", "2e2mu", "4e4mu"]
        for p, lab in zip(paths, labels):
            cfg = load_var_config(p, variable)
            plot_jobs.append((lab, cfg, f"_{lab}"))

    else:
        current_config = config[variable]
        plot_jobs = [("4l", current_config, "")]

    # Now loop over plot_jobs and run your existing plotting code
    for ch_label, current_config, out_suffix in plot_jobs:


        if current_config["first_bin_center"] == -99 and current_config["last_bin_center"] == -99:
            doubleDiff = True
        else:
            doubleDiff = False

        vars = ['fidXS', 'fidXS_scale_up', 'fidXS_scale_dn', 'fidXS_pdf_up', 'fidXS_pdf_dn', 'fidXS_alpha_up', 'fidXS_alpha_dn', 'Boundaries']

        if variable != "mass4l":
            ggh_xs = __import__(current_config['ggh_xs'], globals(), locals(), vars)
            if doubleDiff:
                bin_number = len(ggh_xs.Boundaries)
                bins = [25*i for i in range(0, bin_number+1)]
                bins_plot = np.array(bins)
                bins_c = (bins_plot[1:]+bins_plot[:-1])*0.5
            else:
                bins = ggh_xs.Boundaries
                bins_plot = np.array(bins)
                bins_plot[0]  = current_config['x_lim'][0]
                bins_plot[-1] = current_config['x_lim'][-1]
                bins_c = (bins_plot[1:]+bins_plot[:-1])*0.5
                bins_c[0]  = current_config['first_bin_center']
                bins_c[-1] = current_config['last_bin_center']

            bin_w = np.array([bins_plot[k+1]-bins_plot[k] for k in range(len(bins_plot)-1)])
        # else: mass4l already set bins_plot/bins_c/bin_w (categorical)
            
        if variable == "mass4l":
            per_fs = current_config["_per_fs_configs"]

            def import_mode_array(mode_key):
                out = []
                for cfg in per_fs:
                    m = __import__(cfg[mode_key], globals(), locals(), vars)
                    out.append(float(m.fidXS[0]))
                return np.array(out, dtype=float)

            def import_unc_array(mode_key, source_key):
                out = []
                for cfg in per_fs:
                    m = __import__(cfg[mode_key], globals(), locals(), vars)
                    out.append(abs(float(getattr(m, source_key)[0]) - float(m.fidXS[0])))
                return np.array(out, dtype=float)

            xs = {}
            xs['ggh'] = import_mode_array('ggh_xs')
            xs['vbf'] = import_mode_array('vbf_xs')
            xs['vh']  = import_mode_array('vh_xs')
            xs['tth'] = import_mode_array('tth_xs')

            ggh_powheg_vals = import_mode_array('ggh_powheg_xs')

            # categorical: bin_w = 1, so "_norm" is identical
            ggh_xs_norm = xs['ggh'] / bin_w
            vbf_xs_norm = xs['vbf'] / bin_w
            vh_xs_norm  = xs['vh']  / bin_w
            tth_xs_norm = xs['tth'] / bin_w
            ggh_powheg_xs_norm = ggh_powheg_vals / bin_w

            xh_xs_norm = vbf_xs_norm + vh_xs_norm + tth_xs_norm

            # build per-source up/down arrays (shape: (3, 4)) matching your existing logic
            sources_up = ["fidXS_scale_up", "fidXS_pdf_up", "fidXS_alpha_up"]
            sources_dn = ["fidXS_scale_dn", "fidXS_pdf_dn", "fidXS_alpha_dn"]

            ggh_up = np.stack([import_unc_array('ggh_xs', s) for s in sources_up], axis=0)
            vbf_up = np.stack([import_unc_array('vbf_xs', s) for s in sources_up], axis=0)
            vh_up  = np.stack([import_unc_array('vh_xs',  s) for s in sources_up], axis=0)
            tth_up = np.stack([import_unc_array('tth_xs', s) for s in sources_up], axis=0)

            ggh_dn = np.stack([import_unc_array('ggh_xs', s) for s in sources_dn], axis=0)
            vbf_dn = np.stack([import_unc_array('vbf_xs', s) for s in sources_dn], axis=0)
            vh_dn  = np.stack([import_unc_array('vh_xs',  s) for s in sources_dn], axis=0)
            tth_dn = np.stack([import_unc_array('tth_xs', s) for s in sources_dn], axis=0)

            ggh_powheg_up = np.stack([import_unc_array('ggh_powheg_xs', s) for s in sources_up], axis=0)
            ggh_powheg_dn = np.stack([import_unc_array('ggh_powheg_xs', s) for s in sources_dn], axis=0)

            # Then continue with your existing "madgraph_up = sum([...])" style,
            # but note these are now arrays, so sum along production modes:
            madgraph_up = ggh_up + vbf_up + vh_up + tth_up
            madgraph_dn = ggh_dn + vbf_dn + vh_dn + tth_dn

            powheg_up = ggh_powheg_up + vbf_up + vh_up + tth_up
            powheg_dn = ggh_powheg_dn + vbf_dn + vh_dn + tth_dn

            total_xs = (xs['ggh'] + xs['vbf'] + xs['vh'] + xs['tth'])

            # Combine sources in quadrature + 2% BR, then divide by bin_w (==1)
            unc_th_up = (np.sqrt((np.sqrt(np.sum(madgraph_up**2, axis=0)) / total_xs)**2 + 0.02**2) * total_xs) / bin_w
            unc_th_dn = (np.sqrt((np.sqrt(np.sum(madgraph_dn**2, axis=0)) / total_xs)**2 + 0.02**2) * total_xs) / bin_w

            unc_th_powheg_up = (np.sqrt((np.sqrt(np.sum(powheg_up**2, axis=0)) / total_xs)**2 + 0.02**2) * total_xs) / bin_w
            unc_th_powheg_dn = (np.sqrt((np.sqrt(np.sum(powheg_dn**2, axis=0)) / total_xs)**2 + 0.02**2) * total_xs) / bin_w

        else:

            xs = {}

            # Dynamically import the required module for ggH cross-section
            ggh_xs = __import__(current_config['ggh_xs'], globals(), locals(), vars)
            xs['ggh'] = np.array(ggh_xs.fidXS)
            ggh_xs_norm = xs['ggh'] / bin_w

            # Dynamically import the required module for xH cross-section
            # xh_xs = __import__(current_config['xh_xs'], globals(), locals(), vars)
            # xh_xs_norm = np.array(xh_xs.fidXS) / bin_w

            # Dynamically import the required module for VBF, VH, and ttH cross-section
            vbf_xs = __import__(current_config['vbf_xs'], globals(), locals(), vars)
            xs['vbf'] = np.array(vbf_xs.fidXS)
            vbf_xs_norm = xs['vbf'] / bin_w

            vh_xs = __import__(current_config['vh_xs'], globals(), locals(), vars)
            xs['vh'] = np.array(vh_xs.fidXS)
            vh_xs_norm = xs['vh'] / bin_w
            
            tth_xs = __import__(current_config['tth_xs'], globals(), locals(), vars)
            xs['tth'] = np.array(tth_xs.fidXS)
            tth_xs_norm = xs['tth'] / bin_w

            # Dynamically import the required module for ggH POWHEG cross-section
            ggh_powheg_xs = __import__(current_config['ggh_powheg_xs'], globals(), locals(), vars)
            ggh_powheg_xs_norm = np.array(ggh_powheg_xs.fidXS) / bin_w

            # Dynamically import the required module for ggH MadGraph (w/o NNLOPS reweighting) cross-section
            #ggh_no_nnlops_xs = __import__(current_config['ggh_no_nnlops_xs'], globals(), locals(), vars)
            #ggh_no_nnlops_xs_norm = np.array(ggh_no_nnlops_xs.fidXS) / bin_w
            
            if variable == "pTj1" or variable == "pTj2" or variable == "mjj" or variable == "absdetajj" or variable == "dphijj" or variable == "pTHj" or variable == "pTHjj" or variable == "mHj":
                ggh_xs_norm[0] = ggh_xs_norm[0] * bin_w[0]
                vbf_xs_norm[0] = vbf_xs_norm[0] * bin_w[0]
                vh_xs_norm[0] = vh_xs_norm[0] * bin_w[0]
                tth_xs_norm[0] = tth_xs_norm[0] * bin_w[0]
                ggh_powheg_xs_norm[0] = ggh_powheg_xs_norm[0] * bin_w[0]
                #ggh_no_nnlops_xs_norm[0] = ggh_no_nnlops_xs_norm[0] * bin_w[0]
                
            xh_xs_norm = vbf_xs_norm + vh_xs_norm + tth_xs_norm

            # Theoretical uncertainty
            sources_up = ["fidXS_scale_up", "fidXS_pdf_up", "fidXS_alpha_up"]

            ## [[fidXS_scale_up unc per bin], [fidXS_pdf_up unc per bin], [fidXS_alpha_up unc per bin]]
            ggh_up = [abs(np.array(getattr(ggh_xs, source)) - np.array(ggh_xs.fidXS)) for source in sources_up]
            ggh_powheg_up = [abs(np.array(getattr(ggh_powheg_xs, source)) - np.array(ggh_powheg_xs.fidXS)) for source in sources_up]
            #ggh_no_nnlops_up = [abs(np.array(getattr(ggh_no_nnlops_xs, source)) - np.array(ggh_no_nnlops_xs.fidXS)) for source in sources_up]
            vbf_up = [abs(np.array(getattr(vbf_xs, source)) - np.array(vbf_xs.fidXS)) for source in sources_up]
            vh_up = [abs(np.array(getattr(vh_xs, source)) - np.array(vh_xs.fidXS)) for source in sources_up]
            tth_up = [abs(np.array(getattr(tth_xs, source)) - np.array(tth_xs.fidXS)) for source in sources_up]

            sources_dn = ["fidXS_scale_dn", "fidXS_pdf_dn", "fidXS_alpha_dn"]

            ## [[fidXS_scale_dn unc per bin], [fidXS_pdf_dn unc per bin], [fidXS_alpha_dn unc per bin]]
            ggh_dn = [abs(np.array(getattr(ggh_xs, source)) - np.array(ggh_xs.fidXS)) for source in sources_dn]
            ggh_powheg_dn = [abs(np.array(getattr(ggh_powheg_xs, source)) - np.array(ggh_powheg_xs.fidXS)) for source in sources_dn]
            #ggh_no_nnlops_dn = [abs(np.array(getattr(ggh_no_nnlops_xs, source)) - np.array(ggh_no_nnlops_xs.fidXS)) for source in sources_dn]
            vbf_dn = [abs(np.array(getattr(vbf_xs, source)) - np.array(vbf_xs.fidXS)) for source in sources_dn]
            vh_dn = [abs(np.array(getattr(vh_xs, source)) - np.array(vh_xs.fidXS)) for source in sources_dn]
            tth_dn = [abs(np.array(getattr(tth_xs, source)) - np.array(tth_xs.fidXS)) for source in sources_dn]

            ## sum of the contributions for each production mode
            ## Linear sum as each source of uncertainty is fully correlated across the production modes
            madgraph_up = np.sum([ggh_up,vbf_up,vh_up,tth_up], axis=0)
            powheg_up = np.sum([ggh_powheg_up,vbf_up,vh_up,tth_up], axis=0)
            #no_nnlops_up = np.sum([ggh_no_nnlops_up,vbf_up,vh_up,tth_up], axis=0)

            madgraph_dn = np.sum([ggh_dn,vbf_dn,vh_dn,tth_dn], axis=0)
            powheg_dn = np.sum([ggh_powheg_dn,vbf_dn,vh_dn,tth_dn], axis=0)
            #no_nnlops_dn = np.sum([ggh_no_nnlops_dn,vbf_dn,vh_dn,tth_dn], axis=0)

            ## Sum in quadrature of the different sources of uncertainties + uncertainty on the BR
            unc_th_up = (np.sqrt((np.sqrt(np.sum(np.square(madgraph_up), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth'])) / bin_w
            unc_th_dn = (np.sqrt((np.sqrt(np.sum(np.square(madgraph_dn), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth'])) / bin_w

            unc_th_powheg_up = np.sqrt((np.sqrt(np.sum(np.square(powheg_up), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']) / bin_w
            unc_th_powheg_dn = np.sqrt((np.sqrt(np.sum(np.square(powheg_dn), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']) / bin_w

            #unc_th_no_nnlops_up = np.sqrt((np.sqrt(np.sum(np.square(no_nnlops_up), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']) / bin_w
            #unc_th_no_nnlops_dn = np.sqrt((np.sqrt(np.sum(np.square(no_nnlops_dn), axis=0)) / (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']))**2 + 0.02**2 ) * (xs['ggh']+xs['vbf']+xs['vh']+xs['tth']) / bin_w
        
        if variable == "pTj1" or variable == "pTj2" or variable == "mjj" or variable == "absdetajj" or variable == "dphijj" or variable == "pTHj" or variable == "pTHjj" or variable == "mHj":
            unc_th_up[0] = unc_th_up[0] * bin_w[0]
            unc_th_dn[0] = unc_th_dn[0] * bin_w[0]
            unc_th_powheg_up[0] = unc_th_powheg_up[0] * bin_w[0]
            unc_th_powheg_dn[0] = unc_th_powheg_dn[0] * bin_w[0]
            #unc_th_no_nnlops_up[0] = unc_th_no_nnlops_up[0] * bin_w[0]
            #unc_th_no_nnlops_dn[0] = unc_th_no_nnlops_dn[0] * bin_w[0]
            #bin_w[0] = 30
            #bin_w[-1] = 40



        # Compute expected cross-section and uncertainties
        #exp_xs = np.array(current_config['exp_xs']) * (ggh_xs_norm + xh_xs_norm)
        #err_up = np.array(current_config['err_up']) * (ggh_xs_norm + xh_xs_norm)
        #err_down = np.array(current_config['err_down']) * (ggh_xs_norm + xh_xs_norm)
        #stat_up = np.array(current_config['stat_up']) * (ggh_xs_norm + xh_xs_norm)
        #stat_down = np.array(current_config['stat_down']) * (ggh_xs_norm + xh_xs_norm)

        # SPENCER
        exp_xs = np.array(current_config['exp_xs']) / bin_w
        err_up = np.array(current_config['err_up']) / bin_w
        err_down = np.array(current_config['err_down']) / bin_w
        stat_up = np.array(current_config['stat_up']) / bin_w
        stat_down = np.array(current_config['stat_down']) / bin_w
        
        sys_up = np.sqrt(np.array(err_up)**2 - np.array(stat_up)**2)
        sys_down = np.sqrt(np.array(err_down)**2 - np.array(stat_down)**2)


        if variable == "pTj1" or variable == "pTj2" or variable == "mjj" or variable == "absdetajj" or variable == "dphijj" or variable == "pTHj" or variable == "pTHjj" or variable == "mHj": # or variable == "TBjmax" or variable == "TCjmax":
            exp_xs[0] = current_config['exp_xs'][0]
            
        
        #print("Sys_up:", sys_up)
        #print("Sys_down:", sys_down)

        # #######################################
        # ##### S T A R T   P L O T T I N G #####
        # #######################################

        fig = plt.figure(figsize=(10,8), dpi=120)
        frame1 = fig.add_axes((.1, .35, .8, .6)) # frame1 = fig.add_axes((.1, .35, .8, .8))
        
        if args.no_preliminary:
            cms_label = ""
        else:
            cms_label = "Preliminary"

        era = args.YEAR

        if era == "2022":
            lumi = 7.9804
        elif era == "2022EE":
            lumi = 26.6728
        elif era == '2023preBPix':
            lumi = 17.794
        elif era == '2023postBPix':
            lumi = 9.451
        elif era == '2022full':
            lumi = 34.6
        elif era == '2023full':
            lumi = 27.2
        elif era == '2022_2023':
            lumi = 62
        elif era == 'Run3':
            lumi = 171

        lumi = round(lumi, 2)
        hep.cms.label(cms_label, data=current_config['is_data'], lumi=lumi, fontsize=20, com=13.6)

        # Plot theoretical predictions and experimental data
        plt.stairs((ggh_xs_norm+xh_xs_norm), bins_plot, linewidth=2, label='ggH (NNLOPS + JHUGen + Pythia) + xH', color='brown')
        plt.stairs((ggh_powheg_xs_norm+xh_xs_norm), bins_plot, linewidth=2, label='ggH (POWHEG + JHUGen + Pythia) + xH', color='tab:blue') 
        #plt.stairs((ggh_no_nnlops_xs_norm+xh_xs_norm), bins_plot, linewidth=2, label='ggH (MadGraph5_aMC@NLO + Pythia) + xH', color='tab:purple')
        plt.stairs(xh_xs_norm, bins_plot, linewidth=2, color='green')
        plt.stairs(xh_xs_norm, bins_plot, linewidth=2, label='xH = ttH + VH + VBF (POWHEG + JHUGen + Pythia)', alpha=0.2, color='green', fill=True)
        
        plt.rcParams['hatch.linewidth'] = 2
        # NNLOPS
        for center, value, err_low, err_high, width in zip(bins_c, ggh_xs_norm+xh_xs_norm, unc_th_dn, unc_th_up, bin_w):
            plt.gca().add_patch(plt.Rectangle((center - width/4, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='brown', hatch='///'))
        # POWHEG
        for center, value, err_low, err_high, width in zip(bins_c, ggh_powheg_xs_norm+xh_xs_norm, unc_th_powheg_dn, unc_th_powheg_up, bin_w):
            plt.gca().add_patch(plt.Rectangle((center + width/10, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='tab:blue', hatch='////'))
        # Madgraph w/o NNLOPS
        #for center, value, err_low, err_high, width in zip(bins_c, ggh_no_nnlops_xs_norm+xh_xs_norm, unc_th_no_nnlops_dn, unc_th_no_nnlops_up, bin_w):
        #    plt.gca().add_patch(plt.Rectangle((center + width/3.6, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='tab:purple', hatch='////'))

        # Plot data
        #plt.errorbar(bins_c, exp_xs , yerr=[err_down,err_up], marker = 'o', linestyle = 'None', color = 'k', linewidth = 2, ms=5, capsize=4, label=r'Data (stat $\oplus$ sys unc.)')
        plt.errorbar(bins_c, exp_xs , yerr=[err_down,err_up], marker = 'o', linestyle = 'None', color = 'k', linewidth = 2, ms=5, capsize=4, label=r' Toy Data (stat $\oplus$ sys unc.)') 
        plt.errorbar(bins_c, exp_xs , yerr=[sys_down,sys_up], marker = 'None', linestyle = 'None', color = 'red', linewidth = 6, ms=5, capsize=4, label='Systematic Uncertainty')
        
        # custom settings

        var_map = {
            "mass4l": r"m_{4\ell}",
            "pT4l": r"p^T_{H}",
            "rapidity4l": r"|y_{H}|",
            "massZ1": r"m_{Z1}",
            "massZ2": r"m_{Z2}",
            "Nj": r"N_{jets}",
            "pTj1": r"p^T_{j1}",
            "pTj2": r"p^T_{j2}",
            "mjj": r"m_{jj}",
            "absdetajj": r"|\Delta\eta_{jj}|",
            "dphijj": r"\Delta\Phi_{jj}",
            "pTHj": r"p^T_{Hj}",
            "pTHjj": r"p^T_{Hjj}",
            "mHj": r"m_{Hj}",
            "TCjmax": r"\mathcal{T}_C^{max}",
            "TBjmax": r"\mathcal{T}_B^{max}",
            "rapidity4l_pT4l": r"|y_{H}|p^T_H",
            "massZ1_massZ2":  r"m_{Z1}m_{Z2}",
            "Nj_pT4l": r"N_{jets}p^T_H",
            "pTj1_pTj2": r"p^T_{j1}p^T_{j2}",
            "pTHj_pT4l": r"p^T_{Hj}p^T_{H}",
            "TCjmax_pT4l": r"\mathcal{T}_C^{max}p^T_{H}",
            "pT4l_pTHj": r"p^T_{H} p^T_{Hj}",
            "absdetajj_mjj": r"|\Delta\eta_{jj}|m_{jj}",
            "costhetaZ1": r"cos \theta_{Z_1}",
            "costhetaZ2": r"cos \theta_{Z_2}",
            "costhetastar": r"cos \theta^{*}_{ZZ}",
            "phi": r"\Phi",
            "phi1": r"\Phi_1",
        }

        var = var_map.get(current_config["variable"], current_config["variable"])  # fallback to variable name if not found
        
        if variable == "mass4l":
            plt.ylabel(r"$\sigma_{\text{fid}}$ (fb)", fontsize=20)
        elif doubleDiff:
            plt.ylabel(r"$d^2\sigma_{\text{fid}} / d" + var + r"\;" + current_config["y_unit"] + "$", fontsize=20)
        else:
            plt.ylabel(r"$d\sigma_{\text{fid}} / d" + var + r"\;" + current_config["y_unit"] + "$", fontsize=20)

        #if variable == "pT4l":
            #plt.figtext(.875, .50, r"$\frac{\sigma_{fid}(\it{p}_{\rm T}^{\rm H} > 350.0\ \rm{GeV})}{150.0\ \rm{GeV}}$", horizontalalignment='right', rotation=90, fontsize=16)
        
        if current_config['plot_log']:
            plt.yscale('log')
            
        if "y_lim_top" in current_config.keys(): plt.ylim(top=current_config["y_lim_top"])
        if "y_lim_bottom" in current_config.keys(): plt.ylim(bottom=current_config["y_lim_bottom"])
        plt.xlim(bins_plot[0], bins_plot[-1])
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        legend_location = current_config.get('legend_location', 'upper right')
        fontsize=14
        title_fontsize=14
        #legend = plt.legend(fontsize=fontsize, title='p-value (MadGraph NNLOPS) = '+ str(current_config["pvalue"]), loc= legend_location, title_fontsize=title_fontsize)
        legend = plt.legend(fontsize=fontsize, loc= legend_location, title_fontsize=title_fontsize)
        legend._legend_box.align = "left"
        
        frame1.set_xticklabels([])

        # Plot ratio (Data/Prediction) in a separate frame
        
        frame2 = fig.add_axes((.1,.05,.8,.25))

        ratio_xs = exp_xs / (ggh_xs_norm+xh_xs_norm)
        ratio_powheg = (ggh_powheg_xs_norm+xh_xs_norm) / (ggh_xs_norm+xh_xs_norm)
        #ratio_no_nnlops = (ggh_no_nnlops_xs_norm+xh_xs_norm) / (ggh_xs_norm+xh_xs_norm)
        ratio_madgraph = (ggh_xs_norm+xh_xs_norm) / (ggh_xs_norm+xh_xs_norm) # Dummy, it is always 1
        ratio_err_up = err_up / (ggh_xs_norm+xh_xs_norm)
        ratio_err_down = err_down / (ggh_xs_norm+xh_xs_norm)
        ratio_sys_up = sys_up / (ggh_xs_norm+xh_xs_norm)
        ratio_sys_down = sys_down / (ggh_xs_norm+xh_xs_norm)
        ratio_unc_up = unc_th_up / (ggh_xs_norm+xh_xs_norm)
        ratio_unc_dn = unc_th_dn / (ggh_xs_norm+xh_xs_norm)
        ratio_unc_powheg_up = unc_th_powheg_up / (ggh_xs_norm+xh_xs_norm)
        ratio_unc_powheg_dn = unc_th_powheg_dn / (ggh_xs_norm+xh_xs_norm)
        #ratio_unc_no_nnlops_up = unc_th_no_nnlops_up / (ggh_xs_norm+xh_xs_norm)
        #ratio_unc_no_nnlops_dn = unc_th_no_nnlops_dn / (ggh_xs_norm+xh_xs_norm)

        plt.errorbar(bins_c, ratio_xs , yerr=[ratio_err_down,ratio_err_up], marker = 'o', linestyle = 'None', color = 'k', linewidth = 2, ms=5, capsize=4, label='Toy Data (Stat + Syst)')
        plt.errorbar(bins_c, ratio_xs , yerr=[ratio_sys_down,ratio_sys_up], marker = 'None', linestyle = 'None', color = 'red', linewidth = 6, ms=5, capsize=4, label='Systematic error')
        
        plt.hlines(1, -500,500, color='tab:brown')

        # NNLOPS
        for center, value, err_low, err_high, width in zip(bins_c, ratio_madgraph, ratio_unc_dn, ratio_unc_up, bin_w):
            plt.gca().add_patch(plt.Rectangle((center - width/4, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='brown', hatch='/////')) #/2
        # POWHEG
        for center, value, err_low, err_high, width in zip(bins_c, ratio_powheg, ratio_unc_powheg_dn, ratio_unc_powheg_up, bin_w):
            plt.gca().add_patch(plt.Rectangle((center + width/10, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='tab:blue', hatch='/////'))
        # Madgraph w/o NNLOPS
        #for center, value, err_low, err_high, width in zip(bins_c, ratio_no_nnlops, ratio_unc_no_nnlops_dn, ratio_unc_no_nnlops_up, bin_w):
        #    plt.gca().add_patch(plt.Rectangle((center + width/3.6, value - err_low), width/8, err_low + err_high, fill=False, lw=0, color='tab:blue', hatch='////'))

        plt.stairs(ratio_powheg, bins_plot, linewidth=2, color='tab:blue')
        #plt.stairs(ratio_no_nnlops, bins_plot, linewidth=2, color='tab:blue')

        '''
        if variable == "Nj":
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            if current_config['is_data']:
                custom_ytick_labels = ['0.0', '1.0', '2.0', '3.0']
                custom_yticks = [0, 1, 2, 3]
            else:
                custom_ytick_labels = ['0.0', '0.5', '1.0', '1.5', '2.0']
                custom_yticks = [0, 0.5, 1, 1.5, 2]
            frame2.set_yticklabels(custom_ytick_labels)
            frame2.set_yticks(custom_yticks)
        
        if variable == "pTj1":
            frame2.set_xticklabels(custom_xtick_labels)
            if current_config['is_data']:
                additional_labels = [plt.Text(-25, -2.8, r'$\sigma(\it{N}_{\rm{Jets}}=0)$')] # Unblinded
            else:
                additional_labels = [plt.Text(-25, -1, r'$\sigma(\it{N}_{\rm{Jets}}=0)$')] # Blinded
            for label in additional_labels:
                frame2.text(label.get_position()[0], label.get_position()[1], label.get_text(), rotation=45, fontsize=18)
            # frame2.set_xticklabels(labels)
            frame2.set_xticks(custom_xticks)
            plt.axvline(x=30, color='black', ls='dashed', lw=1, alpha=1)
        '''

        if variable == "mass4l":
            custom_xtick_labels = ['$4\ell$', '$2e2\mu$', '$4\mu$', '$4e$']
            custom_xticks = [0.5, 1.5, 2.5, 3.5]
            plt.ylabel(r'$\sigma_{{fid}}$ (fb)', fontsize=20)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())

        if variable == "Nj":
            custom_xtick_labels = ['0', '1', '2', '3', '$\geq4$']
            custom_xticks = [0.5, 1.5, 2.5, 3.5, 4.5]
            plt.ylabel(r'$\sigma_{{fid}}$ (fb)', fontsize=20)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())

        if variable == "pTj1":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}=0)$', '30', '60', '90', '120', '150', '180', '210', '240']
            custom_xticks = [10, 30, 60, 90, 120, 150, 180, 210, 240]
            plt.axvline(x=30, color='black', ls='dashed', lw=1, alpha=1)
            #plt.figtext(.875, .52, r"$\frac{\sigma_{fid}(\it{p}_{\rm T}^{\rm j_{1}} > 200.0\ \rm{GeV})}{200.0\ \rm{GeV}}$", horizontalalignment='right', rotation=90, fontsize=16)  
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "pTj2":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}} \leq 1)$', '40', '60', '80', '100', '120', '140']
            custom_xticks = [15, 40, 60, 80, 100, 120, 140]
            plt.axvline(x=30, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "mjj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}} \leq 1)$', '0', '75', '150', '225', '300', '375', '450', '525']
            custom_xticks = [-50, 0, 75, 150, 225, 300, 375, 450, 525]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "absdetajj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}} \leq 1)$', '0', '2', '4', '6', '8', '10']
            custom_xticks = [-1, 0,2,4,6,8,10]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "dphijj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}} \leq 1)$', '-3', '-2', '-1', '0', '1', '2', '3']
            custom_xticks = [-3.57, -3,-2,-1,0,1,2,3]
            plt.axvline(x=-3.14, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "mHj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}=0)$', '0', '110', '220', '330', '440', '550', '660', '770', '880']
            custom_xticks = [-50, 0, 110, 220, 330,440,550,660,770,880]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "pTHj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}=0)$', '0', '25', '50', '75', '100', '125', '150', '175', '200']
            custom_xticks = [-12.5,0,25,50,75,100,125,150,175,200]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "pTHjj":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}\leq 1)$', '0', '25', '50', '75', '100']
            custom_xticks = [-12.5,0,25,50,75,100]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')
            
        if variable == "TCjmax":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}=0)$', '0', '20', '40', '60', '80']
            custom_xticks = [-10,0,20,40,60,80]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if variable == "TBjmax":
            custom_xtick_labels = [r'$\sigma(\it{N}_{\rm{Jets}}=0)$', '0', '20', '40', '60', '80']
            custom_xticks = [-10,0,20,40,60,80]
            plt.axvline(x=0, color='black', ls='dashed', lw=1, alpha=1)
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
            frame2.get_xticklabels()[0].set_rotation(30)
            frame2.get_xticklabels()[0].set_ha('right')

        if doubleDiff:
            custom_xtick_labels = []
            for i in range(1, bin_number+1):
                custom_xtick_labels.append(f'Bin {i}')
            custom_xticks = bins_c
            frame2.set_xticklabels(custom_xtick_labels)
            frame2.set_xticks(custom_xticks)
            frame2.xaxis.set_minor_locator(plt.NullLocator())
        
        for b in bins_plot:
            plt.axvline(x=b, color='gray', ls='dashed', lw=1, alpha=0.5)

        plt.ylabel(r'Ratio to MG5+NNLOPS', fontsize=15) #Not MC, since NNLOPS is based on a calculation

        fs_label = "" #if ch_label == "4l" else f" ({ch_label})"
        plt.xlabel(r"$" + var + r"$" + current_config["x_unit"] + fs_label, fontsize=20)

        fontsizex=20
        if variable == "Nj_pT4l":
            fontsizex=16
            
        plt.xticks(fontsize=fontsizex)
        plt.yticks(fontsize=20)

        plt.xlim(bins_plot[0], bins_plot[-1])
        plt.ylim(current_config['y_lim'])

        plt.ylim(0, 2)
        
        # Save the plot to a file
        #plt.savefig(f"./Plots/{current_config['output_name']}.pdf", bbox_inches='tight', dpi=120)
        plt.savefig(f"{path['plots_path']}PLOTS/{current_config['output_name']}{out_suffix}.pdf", bbox_inches='tight', dpi=120)
            
    
