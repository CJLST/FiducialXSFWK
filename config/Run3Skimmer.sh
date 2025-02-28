#!/bin/bash

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Please provide a directory as an argument."
    exit 1
fi

# Define the main path
path="/eos/user/l/lurda/CMS/HZZ/XS_analysis/250226"

# Get the directory from the command line argument
directory="$1"
echo "Processing directory: $directory"

# Loop over the sample list (sampleList.txt)
while IFS= read -r sample; do
    # Ensure the line is not empty and trim leading/trailing spaces
    sample=$(echo "$sample" | xargs)

    if [ -n "$sample" ]; then
        echo "Processing sample: $sample"
        
        # Execute the python command with the corresponding path and directory
        python3 Run3Skimmer.py --input "$path/$directory/$sample/ZZ4lAnalysis.root" --output "$path/$directory/$sample/ZZ4lAnalysis_SKIMMED.root" --mc
    fi
done < "$path/$directory/sampleList.txt"
