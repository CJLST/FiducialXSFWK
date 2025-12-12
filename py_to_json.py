import importlib.util
import json
import sys
import os
import re
import numpy as np

from observables import observables
from binning import binning
from paths import path

def load_results(file_path):
    # Dynamically import the Python file to get the 'resultsXS' dictionary
    spec = importlib.util.spec_from_file_location("results_module", file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.resultsXS

def parse_results(resultsXS, variable_name, year):
    num_bins = len([k for k in resultsXS if variable_name in k and "_statOnly" not in k])

    # Load acceptance uncertainties
    acc_file = f"LHScans/resultsXS_LHScan_expected_{variable_name}_v3.py"
    if not os.path.exists(acc_file):
        print(f"Warning: Acceptance uncertainty file {acc_file} not found.")
        acc_unc = {}
    else:
        acc_unc = load_results(acc_file)

    if "_" in variable_name:
        bins, doubleDiff = binning(variable_name.split("_")[0] + " vs " + variable_name.split("_")[1])
    else:
        bins, doubleDiff = binning(variable_name)
    
    exp_xs, err_up, err_down = [], [], []
    stat_up, stat_down = [], []

    for i in range(num_bins):
        key = f"SM_125_{variable_name}_genbin{i}"
        stat_key = f"{key}_statOnly"

        if key in resultsXS and stat_key in resultsXS:
            central = resultsXS[key]['central']
            main_up = abs(resultsXS[key]['uncerUp'])
            main_dn = abs(resultsXS[key]['uncerDn'])
            statUp  = abs(resultsXS[stat_key]['uncerUp'])
            statDn  = abs(resultsXS[stat_key]['uncerDn'])
            
            # Append for every valid bin
            exp_xs.append(round(central, 3))
            err_up.append(round(main_up, 3))
            err_down.append(round(abs(main_dn), 3))
            stat_up.append(round(statUp, 3))
            stat_down.append(round(abs(statDn), 3))

        else:
            print(f"Warning: Missing keys for bin {i}")


    do_log = 1
    x_unit = " (GeV)"
    y_unit = " (fb/GeV)" 
    y_lim_bottom = 10e-5
    y_lim_top = 10
        
    if doubleDiff:
        last_bin_center = -99
        first_bin_center = -99
        x_lim = [-99, -99]
        y_lim_bottom = 10e-6
        y_lim_top = 10e4
    else: 
        if variable_name == "pT4l":
            x_lim = [0, 280]
            y_lim_bottom = 1e-4
            y_lim_top = 5
        elif variable_name == "rapidity4l":
            do_log = 0
            x_unit = ""
            y_unit = " (fb)"
            x_lim = [0, 2.5]
            y_lim_bottom = 0
            y_lim_top = 5
        elif variable_name == "massZ1":
            x_lim = [40, 120]
            y_lim_bottom = 10e-4
            y_lim_top = 100
        elif variable_name == "massZ2":
            x_lim = [0, 65]
            y_lim_bottom = 10e-4
            y_lim_top = 100
        elif variable_name == "Nj":
            x_unit = ""
            y_unit = " (fb)"
            x_lim = [0, 5]
            y_lim_bottom = 10e-4
            y_lim_top = 100
        elif variable_name == "pTj1":
            x_lim = [-30, 240]
        elif variable_name == "pTj2":
            x_lim = [-10, 140]
        elif variable_name == "mjj":
            x_lim = [-100, 525]
        elif variable_name == "absdetajj":
            x_unit = ""
            y_unit = " (fb)"
            x_lim = [-2, 10]
            y_lim_bottom = 10e-4
        elif variable_name == "dphijj":
            x_unit = ""
            y_unit = " (fb)"
            x_lim = [-5, 3.14159]
            y_lim_bottom = 10e-3
        elif variable_name == "mHj":
            x_lim = [-110, 880]
        elif variable_name == "pTHj":
            x_lim = [-25, 200]
        elif variable_name == "pTHjj":
            x_lim = [-25, 100]
        elif variable_name == "costhetastar" or variable_name == "costhetaZ1" or variable_name == "costhetaZ2":
            x_lim = [-1, 1]
            y_lim_bottom = 10e-3
            y_lim_top = 1000
        elif variable_name == "phi" or variable_name == "phi1":
            x_lim = [-3.14159, 3.14159]
            y_lim_bottom = 10e-3
            y_lim_top = 1000
        else:
            x_lim = [-1000, 1000]

        first_bin_center =  x_lim[0] + (bins[1] - x_lim[0])/2
        last_bin_center = bins[-2] + (x_lim[1] - bins[-2])/2

    return {
        variable_name: {
            "ggh_xs": f"fidXS_NNLOPS_{variable_name}_ggH_{year}",
            "vbf_xs": f"fidXS_{variable_name}_VBFH_{year}",
            "vh_xs": f"fidXS_{variable_name}_VH_{year}",
            "tth_xs": f"fidXS_{variable_name}_ttH_{year}",
            "xh_xs": f"fidXS_{variable_name}_xH_{year}",
            "ggh_powheg_xs": f"fidXS_{variable_name}_ggH_{year}",
            "exp_xs": exp_xs,
            "err_up": err_up,
            "err_down": err_down,
            "stat_up": stat_up,
            "stat_down": stat_down,
            "output_name": variable_name,
            "is_data": 0,
            "last_bin_center": last_bin_center,
            "first_bin_center": first_bin_center,
            "plot_log": do_log,
            "variable": variable_name,
            "x_unit": x_unit,
            "y_unit": y_unit,
            "x_lim": x_lim,
            "y_lim_bottom": y_lim_bottom,
            "y_lim_top": y_lim_top,
            "y_lim": [-1,3],
            "pvalue": -99
        }
    }

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_results_to_json.py <variable_name> <year>")
        sys.exit(1)

    input_file = 'LHScans/resultsXS_LHScan_expected_'+sys.argv[1]+'_v3.py'
    variable_name = sys.argv[1]
    year = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)

    resultsXS = load_results(input_file)
    result_dict = parse_results(resultsXS, variable_name, year)

    json_file = f"jsons/{variable_name}_results_{year}.json"
    with open(json_file, "w") as f:
        json.dump(result_dict, f, indent=4)

    print(f"JSON saved to {json_file}")

if __name__ == "__main__":
    main()
