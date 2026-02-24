#!/bin/bash

# Usage check
# Usage: 
# ./Run3Skimmer.sh MC 2022EE --retry
# ./Run3Skimmer.sh MC 2023 [--retry]
# ./Run3Skimmer.sh Data 2022

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <MC|Data> <Subdirectory> [--retry]"
    exit 1
fi

data_type="$1"
subdir="$2"
retry_mode=false

if [ "$3" == "--retry" ]; then
    retry_mode=true
fi

# Validate input
if [[ "$data_type" != "MC" && "$data_type" != "Data" ]]; then
    echo "First argument must be 'MC' or 'Data'"
    exit 2
fi

# Define paths
base_path="/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/"
#base_path="/eos/home-s/sellissp/HZZ/SAMPLES/"
full_path="$base_path/${subdir}_${data_type}"

# Check directory
if [ ! -d "$full_path" ]; then
    echo "Directory not found: $full_path"
    exit 3
fi

# ----- CASE 1: MC -----
if [ "$data_type" == "MC" ]; then
    sample_list="$full_path/sampleList.txt"

    echo "Generating sample list for MC (excluding *Chunk*)..."
    #find "$full_path" -mindepth 1 -maxdepth 1 -type d ! -name '*Chunk*' -exec basename {} \; | sort > "$sample_list"
    find "$full_path" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort > "$sample_list"
    echo "Sample list has $(wc -l < "$sample_list") entries."
    echo "Starting Run3Skimmer processing for MC..."

    while IFS= read -r sample; do
        sample=$(echo "$sample" | xargs)
        if [ -n "$sample" ]; then
            input_file="$full_path/$sample/ZZ4lAnalysis.root"
            output_file="$full_path/$sample/ZZ4lAnalysis_SKIMMED.root"
	    #output_file="$full_path/$sample/ZZ4lAnalysis_FR.root" 
	    
            if [ -f "$input_file" ]; then
                if [ "$retry_mode" = true ] && [ -f "$output_file" ]; then
                    echo "-- Skipped (already skimmed): $sample"
                    continue
                fi
                echo "-> Skimming $sample"
                python3 Run3Skimmer.py --input "$input_file" --output "$output_file" --mc
		#python3 NanoConverter.py --input "$input_file" --output "$output_file" --mc --skipZL
            else
                echo "!! Missing input file for sample: $sample"
            fi
        fi
    done < "$sample_list"

# ----- CASE 2: Data -----
else
    echo "Processing Data files directly in: $full_path"

    year=""
    if [[ "$subdir" == *"2022"* ]]; then
        year="2022"
        data_files=("Data_eraCD_preEE.root" "Data_eraEFG_postEE.root")
    elif [[ "$subdir" == *"2023"* ]]; then
        year="2023"
        data_files=("Data_eraC_preBPix.root" "Data_eraD_postBPix.root")
    elif [[ "$subdir" == *"2024"* ]]; then
	sample_list="$full_path/sampleList1.txt"
	#find "$full_path" -mindepth 1 -maxdepth 1 -type d ! -name '*Chunk*' -exec basename {} \; | sort > "$sample_list"
	echo "Sample list has $(wc -l < "$sample_list") entries."
	echo "Starting Run3Skimmer processing for 2024 Data..."
	while IFS= read -r sample; do
	    input_file="$full_path/$sample/ZZ4lAnalysis.root"
	    output_file="$full_path/$sample/ZZ4lAnalysis_SKIMMED.root"
	    if [ -f "$input_file" ]; then
		if [ "$retry_mode" = true ] && [ -f "$output_file" ]; then
		    echo "-- Skipped (already skimmed): $sample"
		    continue
		fi
		echo "-> Skimming $sample"
		python3 Run3Skimmer.py --input "$input_file" --output "$output_file"
	    else
                echo "!! Missing input file for sample: $sample"
            fi
	done < "$sample_list"
	
    else
        echo "Could not determine year (2022/2023) from directory name: $subdir"
        exit 4
    fi

    for file in "${data_files[@]}"; do
        input_file="$full_path/$file"
        output_file="${input_file/.root/_SKIMMED.root}"

        if [ -f "$input_file" ]; then
            echo "-> Skimming $file"
            python3 Run3Skimmer.py --input "$input_file" --output "$output_file"
        else
            echo "!! Missing file: $input_file"
        fi
    done
fi

echo "All done."

