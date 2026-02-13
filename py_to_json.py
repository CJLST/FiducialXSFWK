import importlib.util
import json
import sys
import os
import re
import numpy as np

from observables import observables
from binning import binning
from paths import path

#SPECIAL_OBS = {"massZ1","massZ2","costhetaZ1","costhetaZ2","costhetastar","phi","phi1"}
SPECIAL_OBS = {}

def load_results(file_path):
    # Dynamically import the Python file to get the 'resultsXS' dictionary
    spec = importlib.util.spec_from_file_location("results_module", file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.resultsXS

def _get_key(var, ch, i):
    # v4 has SM_125_<var>_<ch>_genbin<i> for ch in {4l,2e2mu}
    # and SM_125_<var>_genbin<i> is NOT what you showed for v4 special obs
    return f"SM_125_{var}_{ch}_genbin{i}"

def _get_stat_key(key):
    return f"{key}_statOnly"

def _quad(a, b):
    return float(np.sqrt(float(a)**2 + float(b)**2))

def _get_key_v4(var, ch, i):
    return f"SM_125_{var}_{ch}_genbin{i}"

def _get_key_v3_inclusive(var, i):
    # v3 inclusive is usually old style (no channel)
    return f"SM_125_{var}_genbin{i}"

def _theory_suffix(variable_name, channel):
    if variable_name in SPECIAL_OBS:
        if channel == "4l":
            return ""
        if channel == "2e2mu":
            return "_2e2mu"
        if channel == "4e4mu":
            return "_4e4mu"
    return "" if channel == "4l" else f"_{channel}"

def parse_results(channel, variable_name, year):

    # --- Choose scan output model by observable/channel ---
    if variable_name in SPECIAL_OBS:
        # per your rule:
        #   inclusive (4e+4mu+2e2mu) -> v3 4l
        #   2e2mu -> v4 2e2mu
        #   4e+4mu -> v4 4l
        if channel == "4l":
            input_file = f'LHScans/resultsXS_LHScan_expected_{variable_name}_v3.py'
        elif channel in ("2e2mu", "4e4mu"):
            input_file = f'LHScans/resultsXS_LHScan_expected_{variable_name}_v4.py'
        else:
            raise ValueError(f"Unsupported channel for SPECIAL_OBS: {channel}")

    elif variable_name == "mass4l":
        input_file = f'LHScans/resultsXS_LHScan_expected_{variable_name}_v3.py' if channel == "4l" else f'LHScans/resultsXS_LHScan_expected_{variable_name}_v2.py'
    else:
        input_file = f'LHScans/resultsXS_LHScan_expected_{variable_name}_v3.py' if channel == "4l" else f'LHScans/resultsXS_LHScan_expected_{variable_name}_v2.py'

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)

    resultsXS = load_results(input_file)

    # --- Count bins correctly for each format ---
    if variable_name in SPECIAL_OBS:
        if channel == "4l":
            # v3 inclusive: old-style keys
            num_bins = len([k for k in resultsXS
                            if k.startswith(f"SM_125_{variable_name}_genbin")
                            and "_statOnly" not in k])
        else:
            # v4: channelized keys; count using the actual v4 channel you're reading
            v4ch = "2e2mu" if channel == "2e2mu" else "4l"  # 4e4mu uses v4 "4l"
            num_bins = len([k for k in resultsXS
                            if k.startswith(f"SM_125_{variable_name}_{v4ch}_genbin")
                            and "_statOnly" not in k])

    elif variable_name == "mass4l" and channel != "4l":
        num_bins = 1
    else:
        num_bins = len([k for k in resultsXS if k.startswith(f"SM_125_{variable_name}_genbin") and "_statOnly" not in k])

    if "_" in variable_name:
        bins, doubleDiff = binning(variable_name.split("_")[0] + " vs " + variable_name.split("_")[1])
    else:
        bins, doubleDiff = binning(variable_name)

    exp_xs, err_up, err_down = [], [], []
    stat_up, stat_down = [], []

    if variable_name == "mass4l" and channel != "4l":
        # unchanged mass4l special handling
        key = f"SM_125_{variable_name}_{channel}_genbin0"
        stat_key = f"{key}_statOnly"
        if key in resultsXS and stat_key in resultsXS:
            central = resultsXS[key]['central']
            main_up = abs(resultsXS[key]['uncerUp'])
            main_dn = abs(resultsXS[key]['uncerDn'])
            statUp  = abs(resultsXS[stat_key]['uncerUp'])
            statDn  = abs(resultsXS[stat_key]['uncerDn'])

            exp_xs.append(round(central, 3))
            err_up.append(round(main_up, 3))
            err_down.append(round(abs(main_dn), 3))
            stat_up.append(round(statUp, 3))
            stat_down.append(round(abs(statDn), 3))
        else:
            print(f"Warning: Missing keys for channel {channel} in mass4l")

    else:
        for i in range(num_bins):

            if variable_name in SPECIAL_OBS:
                if channel == "4l":
                    # v3 inclusive
                    key = _get_key_v3_inclusive(variable_name, i)
                elif channel == "2e2mu":
                    # v4 2e2mu
                    key = _get_key_v4(variable_name, "2e2mu", i)
                elif channel == "4e4mu":
                    # v4 4l == 4e+4mu
                    key = _get_key_v4(variable_name, "4l", i)
                else:
                    raise ValueError(f"Unsupported channel for SPECIAL_OBS: {channel}")

                stat_key = f"{key}_statOnly"

                if key in resultsXS and stat_key in resultsXS:
                    central = resultsXS[key]['central']
                    main_up = abs(resultsXS[key]['uncerUp'])
                    main_dn = abs(resultsXS[key]['uncerDn'])
                    statUp  = abs(resultsXS[stat_key]['uncerUp'])
                    statDn  = abs(resultsXS[stat_key]['uncerDn'])

                    exp_xs.append(round(central, 3))
                    err_up.append(round(main_up, 3))
                    err_down.append(round(main_dn, 3))
                    stat_up.append(round(statUp, 3))
                    stat_down.append(round(statDn, 3))
                else:
                    print(f"Warning: Missing keys for {channel}, bin {i}: {key} / {stat_key}")

            else:
                # your existing non-SPECIAL_OBS logic (unchanged)
                if channel in ["4l", "2e2mu"]:
                    key = f"SM_125_{variable_name}_genbin{i}"
                    stat_key = f"{key}_statOnly"
                    if key in resultsXS and stat_key in resultsXS:
                        central = resultsXS[key]['central']
                        main_up = abs(resultsXS[key]['uncerUp'])
                        main_dn = abs(resultsXS[key]['uncerDn'])
                        statUp  = abs(resultsXS[stat_key]['uncerUp'])
                        statDn  = abs(resultsXS[stat_key]['uncerDn'])
                        exp_xs.append(round(central, 3))
                        err_up.append(round(main_up, 3))
                        err_down.append(round(main_dn, 3))
                        stat_up.append(round(statUp, 3))
                        stat_down.append(round(statDn, 3))
                    else:
                        print(f"Warning: Missing keys for {channel}, bin {i}: {key} / {stat_key}")
                else:
                    raise ValueError(f"Unsupported channel: {channel}")

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
        if variable_name == "mass4l":
            do_log = 0
            x_lim = [0, 4]
            y_lim_bottom = 0
            y_lim_top = 6
            x_unit = ""
            y_unit = " (fb)"
        elif variable_name == "pT4l":
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
            x_unit = ""
            y_unit = " (fb)"
        elif variable_name == "phi" or variable_name == "phi1":
            x_lim = [-3.14159, 3.14159]
            y_lim_bottom = 10e-3
            y_lim_top = 1000
            x_unit = ""
            y_unit = " (fb)"
        elif variable_name == "TCjmax":
            x_lim = [-20, 80]
        elif variable_name == "TBjmax":
            x_lim = [-20, 80]
        else:
            x_lim = [-1000, 1000]

        first_bin_center =  x_lim[0] + (bins[1] - x_lim[0])/2
        last_bin_center = bins[-2] + (x_lim[1] - bins[-2])/2

    channel_tag = _theory_suffix(variable_name, channel)

    return {
        variable_name: {
            "ggh_xs": f"fidXS_NNLOPS_{variable_name}_ggH{channel_tag}_{year}",
            "vbf_xs": f"fidXS_{variable_name}_VBFH{channel_tag}_{year}",
            "vh_xs": f"fidXS_{variable_name}_VH{channel_tag}_{year}",
            "tth_xs": f"fidXS_{variable_name}_ttH{channel_tag}_{year}",
            "xh_xs": f"fidXS_{variable_name}_xH{channel_tag}_{year}",
            "ggh_powheg_xs": f"fidXS_{variable_name}_ggH{channel_tag}_{year}",
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
        print("Usage: python py_to_json.py <variable_name> <year>")
        sys.exit(1)

    variable_name = sys.argv[1]
    year = sys.argv[2]

    if variable_name == "mass4l":
        channels = ["4l", "4e", "4mu", "2e2mu"]
    elif variable_name in SPECIAL_OBS:
        channels = ["4l", "2e2mu", "4e4mu"]
    else:
        channels = ["4l"]

    for channel in channels:

        channel_tag = "" if channel == "4l" else f"_{channel}"

        result_dict = parse_results(channel, variable_name, year)
        json_file = f"jsons/{variable_name}_results{channel_tag}_{year}.json"

        with open(json_file, "w") as f:
            json.dump(result_dict, f, indent=4)

        print(f"JSON saved to {json_file}")

if __name__ == "__main__":
    main()
