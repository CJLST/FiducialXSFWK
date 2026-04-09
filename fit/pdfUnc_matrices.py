import numpy as np
import importlib.util
import sys
import os
import matplotlib.pyplot as plt

sys.path.append('../helperstuff/')
from paths import path

def load_module_from_file(filepath):
    """Load a Python module from a given .py file path."""
    module_name = os.path.splitext(os.path.basename(filepath))[0]
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def build_all_matrices(eff_num_var, eff_den_var, NNLOPS):
    import numpy as np

    variations = ['pdf', 'qcd', 'as']

    if NNLOPS:
        prod_modes = ['ggH125_NNLOPS']
    else:
        prod_modes = sorted({k.split('_')[0] for k in eff_num_var.keys()})

    keys = list(eff_num_var.keys())

    genbins = sorted(set(
        int(part.replace('genbin',''))
        for k in keys
        for part in k.split('_')
        if 'genbin' in part
    ))

    recobins = sorted(set(
        int(part.replace('recobin',''))
        for k in keys
        for part in k.split('_')
        if 'recobin' in part
    ))

    n_gen = len(genbins)
    n_reco = len(recobins)

    matrices_all = {}  # final structure with all prod modes

    # Loop over production modes
    for prod in prod_modes:
        matrices_all[prod] = {}
        # Filter keys for this production mode
        prod_keys = [k for k in keys if k.startswith(prod + "_")]

        # Extract final states from keys
        final_states = sorted({k.split('_')[1] for k in prod_keys})

        for fs in final_states:
            matrices = {var: {} for var in variations}
            variation_indices = {var: set() for var in variations}

            # Keys for this prod + final state
            fs_keys = [k for k in prod_keys if f"_{fs}_" in k]

            # Identify variation indices
            for key in fs_keys:
                for var in variations:
                    if var in eff_num_var[key]:
                        variation_indices[var].update(eff_num_var[key][var].keys())

            # Initialize matrices
            for var in variations:
                for idx in variation_indices[var]:
                    matrices[var][idx] = np.full((n_gen, n_reco), np.nan)

            # Fill matrices
            for key in fs_keys:
                parts = key.split('_')
                genbin_part = next((p for p in parts if 'genbin' in p), None)
                recobin_part = next((p for p in parts if 'recobin' in p), None)

                genbin_idx = int(genbin_part.replace('genbin',''))
                recobin_idx = int(recobin_part.replace('recobin',''))

                genbin_pos = genbins.index(genbin_idx)
                recobin_pos = recobins.index(recobin_idx)

                for var in variations:
                    if var not in eff_num_var[key] or var not in eff_den_var[key]:
                        continue
                    for var_idx in eff_num_var[key][var]:
                        num = eff_num_var[key][var][var_idx]
                        den = eff_den_var[key][var].get(var_idx, None)
                        if den is None or den == 0:
                            matrices[var][var_idx][genbin_pos, recobin_pos] = np.nan
                        else:
                            matrices[var][var_idx][genbin_pos, recobin_pos] = num / den

            matrices_all[prod][fs] = matrices

    # Build "allH125" as sum over all prod modes
    matrices_all["allH125"] = {}
    # Use final states from first production mode as reference
    ref_prod = prod_modes[0]
    for fs in matrices_all[ref_prod]:
        matrices = {var: {} for var in variations}
        for var in variations:
            # Get all variation indices
            var_indices = set()
            for prod in prod_modes:
                var_indices.update(matrices_all[prod][fs][var].keys())
            for idx in var_indices:
                sum_num = np.zeros((n_gen, n_reco))
                sum_den = np.zeros((n_gen, n_reco))
                for prod in prod_modes:
                    mat = matrices_all[prod][fs][var].get(idx, np.zeros((n_gen, n_reco)))
                    # Retrieve numerator and denominator to reconstruct sum? We'll approximate by multiplying mat * den?
                    # But we don't have original numerator/denominator now; we need to recompute sum
                    # So better: sum directly from eff_num_var and eff_den_var
                    for key in keys:
                        if key.startswith(prod + f"_{fs}_") and var in eff_num_var[key] and idx in eff_num_var[key][var]:
                            parts = key.split('_')
                            genbin_part = next((p for p in parts if 'genbin' in p), None)
                            recobin_part = next((p for p in parts if 'recobin' in p), None)
                            genbin_idx = int(genbin_part.replace('genbin',''))
                            recobin_idx = int(recobin_part.replace('recobin',''))
                            genbin_pos = genbins.index(genbin_idx)
                            recobin_pos = recobins.index(recobin_idx)
                            num = eff_num_var[key][var][idx]
                            den = eff_den_var[key][var].get(idx, 0)
                            sum_num[genbin_pos, recobin_pos] += num
                            sum_den[genbin_pos, recobin_pos] += den
                # Compute ratio
                with np.errstate(invalid='ignore', divide='ignore'):
                    mat_all = sum_num / sum_den
                    mat_all[sum_den == 0] = np.nan
                matrices[var][idx] = mat_all
        matrices_all["allH125"][fs] = matrices

    return matrices_all, genbins, recobins, keys


import numpy as np

def build_matrices(matrices_all_input, genbins, recobins):
    """
    Compute RMS (for pdf) and envelope (max/min for qcd, as) 
    using the matrices_all dict returned by build_all_matrices.
    
    Parameters:
        matrices_all_input : dict
            Output of build_all_matrices
        genbins : list of int
        recobins : list of int
    
    Returns:
        matrices_variation : dict
            Nested dictionary with same structure as matrices_all, 
            but containing RMS / max-min values
    """
    import numpy as np

    variations = ['pdf', 'qcd', 'as']

    matrices_variation = {}

    # Loop over all prod modes in the input (including allH125)
    for prod, fs_dict in matrices_all_input.items():
        matrices_variation[prod] = {}

        for fs, var_dict in fs_dict.items():
            matrices = {
                'pdf': {'rms': np.full((len(genbins), len(recobins)), np.nan)},
                'qcd': {'max': np.full((len(genbins), len(recobins)), np.nan),
                        'min': np.full((len(genbins), len(recobins)), np.nan)},
                'as':  {'max': np.full((len(genbins), len(recobins)), np.nan),
                        'min': np.full((len(genbins), len(recobins)), np.nan)}
            }

            # Loop over variations
            for var in variations:
                # Get all available variation indices
                var_indices = list(var_dict[var].keys())

                if not var_indices:
                    continue

                # Create a 3D array: (n_indices, n_gen, n_reco)
                mats = np.array([var_dict[var][idx] for idx in var_indices])

                if var == "pdf":
                    # Compute RMS along the variation axis
                    rms = np.sqrt(np.nanmean((mats - np.nanmean(mats, axis=0))**2, axis=0))
                    matrices['pdf']['rms'] = rms
                else:
                    # Compute max / min along the variation axis
                    matrices[var]['max'] = np.nanmax(mats, axis=0)
                    matrices[var]['min'] = np.nanmin(mats, axis=0)

            matrices_variation[prod][fs] = matrices

    return matrices_variation


def compute_percent_variations(matrices, matrices_all):
    """
    Compute percent variation wrt nominal for each production mode and final state.

    - For PDF: (RMS / nominal) * 100
    - For QCD and AS: ((max - nominal) / nominal) * 100 and ((min - nominal) / nominal) * 100

    Args:
        matrices : dict
            Output from build_matrices_from_eff (contains 'pdf'->'rms', 'qcd'->'max/min', 'as'->'max/min')
        matrices_all : dict
            Output from build_all_matrices (contains all variation indices, including nominal)

    Returns:
        dict with the same structure as matrices but containing percent differences.
    """
    import numpy as np

    percent_diffs_all = {}

    for prod in matrices.keys():
        percent_diffs_all[prod] = {}

        for fs in matrices[prod].keys():
            percent_diffs_all[prod][fs] = {}

            # --- Retrieve nominal matrix ---
            # always take nominal as "pdf" index 1 from matrices_all
            try:
                nominal_matrix = matrices_all[prod][fs]['pdf'][1]
            except KeyError:
                raise KeyError(f"Nominal ('pdf' index 1) missing for {prod}, {fs}")

            # --- Loop over the variations available in 'matrices' ---
            for var in matrices[prod][fs].keys():
                percent_diffs_all[prod][fs][var] = {}

                for key, mat in matrices[prod][fs][var].items():
                    if mat is None or np.all(np.isnan(mat)):
                        percent_diffs_all[prod][fs][var][key] = np.full_like(nominal_matrix, np.nan)
                        continue

                    with np.errstate(divide='ignore', invalid='ignore'):
                        if var == "pdf":
                            diff = (mat / nominal_matrix) * 100.0
                        else:
                            diff = ((mat - nominal_matrix) / nominal_matrix) * 100.0

                        diff[np.isnan(nominal_matrix)] = np.nan

                    percent_diffs_all[prod][fs][var][key] = np.abs(diff)

    return percent_diffs_all

def plot_and_save_matrices(matrices, obsName, year, genbins, recobins, input_file, dir):
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    os.makedirs("pdfUnc_matrices", exist_ok=True)

    colors = {
        "pdf": "Reds",
        "qcd": "Blues",
        "as": "Greens"
    }

    v_min, v_max = (0, 100) if dir == "percent_diffs" else (0, 1)

    # Extract base name for labeling
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    name = obsName

    # Loop over production modes
    for prod in sorted(matrices.keys()):
        for fs in sorted(matrices[prod].keys()):  # loop over final states
            for var in matrices[prod][fs]:
                for idx, mat in matrices[prod][fs][var].items():
                    plt.figure(figsize=(8, 6))

                    # Choose colormap
                    if var == "pdf" and idx == 1:
                        cmap = plt.get_cmap("Greys")  # Nominal PDF RMS
                    else:
                        cmap = plt.get_cmap(colors[var])

                    im = plt.imshow(mat, interpolation='nearest', origin='lower', aspect='auto',
                                    cmap=cmap, vmin=v_min, vmax=v_max)

                    # Print values on cells
                    for i in range(mat.shape[0]):
                        for j in range(mat.shape[1]):
                            val = mat[i, j]
                            if not np.isnan(val):
                                text_fmt = f"{val:.5f}%" if dir == "percent_diffs" else f"{val:.5f}"
                                plt.text(j, i, text_fmt, ha="center", va="center", color="black", fontsize=8)

                    plt.colorbar(im, label=f"{var} variation index {idx}")
                    plt.title(f"{base_name} — {prod} — {fs} — {var} idx {idx}")
                    plt.xlabel("Reco bin index")
                    plt.ylabel("Gen bin index")
                    plt.xticks(ticks=np.arange(len(recobins)), labels=recobins)
                    plt.yticks(ticks=np.arange(len(genbins)), labels=genbins)
                    plt.tight_layout()

                    filename = f"{path['plots_path']}UNC_MATRICES/{name}/{year}/{dir}/{base_name}_{prod}_{fs}_{var}_varIdx{idx}.png"
                    os.makedirs(os.path.dirname(filename), exist_ok=True)
                    plt.savefig(filename)
                    plt.close()

    #print(f"Saved matrices to directory: pdfUnc_matrices/{name}/{dir}")

import numpy as np

def transpose(mat):
    """
    Return a copy of mat where the top-left square block of size n = min(rows, cols)
    is replaced by its transpose. Values outside that square block are preserved.
    """
    mat = np.array(mat, copy=True)
    rows, cols = mat.shape
    n = min(rows, cols)
    # Extract top-left square block and transpose it
    block = mat[:n, :n].T
    out = mat.copy()
    out[:n, :n] = block
    return out

def transpose_all(obj):
    """Recursively transpose the top-left square block of all numpy arrays in nested dict/list."""
    if isinstance(obj, dict):
        return {k: transpose_all(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [transpose_all(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return transpose(obj)
    return obj

import importlib.util
import numpy as np
from pathlib import Path

def safe_get(x, j):
    """Return x[j] unless missing/NaN/None, else return 0."""
    try:
        val = x[j]
    except (IndexError, TypeError):
        return 0  # j out of range or x not indexable
    
    # If it's a list with one element, unwrap it
    if isinstance(val, (list, tuple)) and len(val) == 1:
        val = val[0]
    
    # Convert to float and handle NaN/None
    try:
        val = float(val)
    except (ValueError, TypeError):
        return 0
    
    if np.isnan(val):
        return 0
    
    return val

def append_uncertainties(input_file, obsName, year, physicalModel, percent_diffs_all):
 
    variable = obsName

    if variable == "pT4l":
        variable = "PTH"
    if variable == "pT4l_zzfloating":
        variable = "PTH_zzfloating"
    if variable == "rapidity4l":
        variable = "YH"
    if variable == "rapidity4l_zzfloating":
        variable = "YH_zzfloating"

    out_dir = Path(f"../datacard/datacard_{year}")

    for fs, values in percent_diffs_all.items():

        if fs == "4l":
            continue

        recobins = len(values["pdf"]["rms"])
        genbins = len(values["pdf"]["rms"][0])

        for i in range(recobins):

            # Prepare target file path
            datacard_path = out_dir / f"hzz4l_{fs}S_13TeV_xs_{variable}_bin{i}_{physicalModel}.txt"
            print(f"Processing datacard: {datacard_path}")

            # Build the lines
            pdf_rms = values["pdf"]["rms"][i]
            qcd_max, qcd_min = values["qcd"]["max"][i], values["qcd"]["min"][i]
            as_max, as_min = values["as"]["max"][i], values["as"]["min"][i]

            pdf_, qcd_, as_ = [], [], []

            for j in range(genbins):

                pdf_rms_scalar = safe_get(pdf_rms, j)
                qcd_max_scalar = safe_get(qcd_max, j)
                qcd_min_scalar = safe_get(qcd_min, j)
                as_max_scalar  = safe_get(as_max, j)
                as_min_scalar  = safe_get(as_min, j)

                if pdf_rms_scalar == 100: 
                    pdf_rms_scalar = 0
                if qcd_max_scalar == 100: 
                    qcd_max_scalar = 0
                if qcd_min_scalar == 100:
                    qcd_min_scalar = 0
                if as_max_scalar == 100: 
                    as_max_scalar = 0
                if as_min_scalar == 100: 
                    as_min_scalar = 0
                    
                pdf_.append(f"{1+(pdf_rms_scalar/100):.7f}/{1-(pdf_rms_scalar/100):.7f}")
                qcd_.append(f"{1+(qcd_max_scalar/100):.7f}/{1-(qcd_min_scalar/100):.7f}")
                as_.append(f"{1+(as_max_scalar/100):.7f}/{1-(as_min_scalar/100):.7f}")

                pdf_.append(" ")
                qcd_.append(" ")
                as_.append(" ")

            # Join into strings
            pdf_str = "CMS_HIG25015_pdf_effMatrix lnN " + "".join(pdf_) + "- - - - -\n"
            qcd_str = "CMS_HIG25015_QCDscale_effMatrix lnN " + "".join(qcd_) + "- - - - -\n"
            as_str  = "CMS_HIG25015_alphaS_effMatrix lnN " + "".join(as_) + "- - - - -\n"

            new_lines = [pdf_str, qcd_str, as_str]

            #print(new_lines)

            # Read current file (if exists)
            if datacard_path.exists():
                with open(datacard_path, "r") as f:
                    lines = f.readlines()
            else:
                print(f"file {datacard_path} does not exist")
                
            filtered = [line for line in lines if not line.strip().startswith(("CMS_HIG25015_pdf_effMatrix", "CMS_HIG25015_QCDscale_effMatrix", "CMS_HIG25015_alphaS_effMatrix"))]
            filtered.extend(new_lines)

            # Append to datacard
            with open(datacard_path, "w") as f:
                f.writelines(filtered)

            print(f"Appended uncertainties to {datacard_path}")


def run_pdf_unc_matrices(input_file, obsName, year, physicalModel):
    if "NNLOPS" in input_file:
        NNLOPS = True
    else:
        NNLOPS = False

    if 'zzfloating' in input_file:
        input_file = input_file.replace('_zzfloating', '')

    module = load_module_from_file(input_file)
    eff_num_var = module.eff_num_var
    eff_den_var = module.eff_den_var

    matrices_all, genbins, recobins, keys = build_all_matrices(eff_num_var, eff_den_var, NNLOPS)
    matrices = build_matrices(matrices_all, genbins, recobins)
    percent_diffs = compute_percent_variations(matrices, matrices_all)
    append_uncertainties(input_file, obsName, year, physicalModel, transpose_all(percent_diffs)['allH125'])

    #plot_and_save_matrices(matrices_all, obsName, year, genbins, recobins, input_file, "all")
    #plot_and_save_matrices(matrices, obsName, year, genbins, recobins, input_file, "variations")
    #plot_and_save_matrices(percent_diffs, obsName, year, genbins, recobins, input_file, "percent_diffs")


    print("Processing complete.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python pdfUnc_matrices.py <input_file.py>")
        sys.exit(1)
    run_pdf_unc_matrices(sys.argv[1])