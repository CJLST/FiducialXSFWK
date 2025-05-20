#!/bin/bash
# How to Run: ./copy_to_www.sh costhetaZ1 2022EE
# Works on lxplus 

# Check if observable and year are provided
if [[ $# -lt 2 ]] ; then
    echo "Usage: $0 <observable> <year>"
    echo "Example: $0 costhetaZ1 2022"
    echo "Optional: Provide 'overwrite' as the third argument to overwrite only the year-specific directory."
    exit 1
fi

observable=$1
year=$2
overwrite_year=$3 

path="/eos/user/m/mmanoni/www/HZZ_Fiducial/" #Run3

# Check if observable directory exists
if [ -d "$path/${observable}" ]; then
  echo "Observable directory exists: $observable"
else
  # Observable directory does not exist, create it
  mkdir $path/${observable}
  cp $path/index.php $path/${observable}
fi

mkdir $path/${observable}/datacard_${year}
cp $path/index.php $path/${observable}
cp $path/index.php $path/${observable}/datacard_${year}

mkdir $path/${observable}/templatesBkgs_${year}
cp $path/index.php $path/${observable}/templatesBkgs_${year}/.
# Move tempaltes plots
cp templates/plots/${year}/${observable}/* $path/${observable}/templatesBkgs_${year}/.

mkdir $path/${observable}/eff_${year}
mkdir $path/${observable}/nonfid_${year}

cp $path/index.php $path/${observable}/eff_${year}
cp $path/index.php $path/${observable}/nonfid_${year}
# Move coefficients plots
cp coefficients/matrix_eff/${year}/eff_${year}_${observable}_* $path/${observable}/eff_${year}/.
cp coefficients/matrix_nonfid/${year}/nonFid_${year}_${observable}_* $path/${observable}/nonfid_${year}/.

# Move datacards
if [ ${observable} = "njets_pt30_eta4p7" ]; then
  obs_datacards="NJ"
elif [ ${observable} = "pT4l" ]; then
  obs_datacards="PTH"
elif [ ${observable} = "rapidity4l" ]; then
  obs_datacards="YH"
elif [ ${observable} = "pTj1" ]; then
  obs_datacards="PTJET"
else
  obs_datacards=${observable}
fi

cp datacard/datacard_${year}/hzz4l_*_13TeV_xs_${obs_datacards}_bin*_* $path/${obs_datacards}/datacard_${year}/.
if [ ${1} = "pT4l_kL" ]; then
  cp datacard/datacard_${year}/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_kLambda* $path/${observable}/datacard_${year}/.
  cp datacard/hzz4l_all_13TeV_xs_${observable}_bin_kLambda* $path/${observable}/datacard_${year}/.
else
  cp datacard/datacard_${year}/hzz4l_*_13TeV_xs_SM_125_${obs_datacards}_v* $path/${observable}/datacard_${year}/.
  cp datacard/hzz4l_all_13TeV_xs_${observable}_bin_v* $path/${observable}/datacard_${year}/.
fi

mkdir $path/${observable}/LHScans_${year}
cp $path/index.php $path/${observable}/LHScans_${year}
if [ ${1} = "pT4l_kL" ]; then
  cp LHScans/plots/${year}_lhscan_compare_${observable}_kappa* $path/${observable}/LHScans_${year}/.
else
  cp LHScans/plots/${year}_lhscan_compare_${observable}_r* $path/${observable}/LHScans_${year}/.
fi

mkdir $path/${observable}/impacts_${year}
cp $path/index.php $path/${observable}/impacts_${year}
cp impacts/impacts_${year}_*_${observable}*  $path/${observable}/impacts_${year}/.

mkdir $path/${observable}/combine_files_${year}
cp $path/index.php $path/${observable}/combine_files_${year}
cp combine_files/higgsCombine_${observable}* $path/${observable}/combine_files_${year}/.
cp combine_files/SM_125_all_13TeV_xs_${observable}_*_${year}*  $path/${observable}/combine_files_${year}/.

#cp plots/${1}/asimov/${1}_unfoldwith* $path/${1}/.
#cp plots/${1}/asimov/corr_${1}_*.png $path/${1}/.
#cp plots/${1}/asimov/model/* $path/${1}/.
#if [ ${1} = "pT4l_kL" ]; then
# cp impacts/${1}/impacts_v3_kappa_lambda* $path/${1}/.
# cp LHScans/plots/lhscan_compare_${1}_kappa* $path/${1}/.
#else
# cp impacts/${1}/impacts_*_${1}_*_asimov* $path/${1}/.
# cp LHScans/plots/lhscan_compare_${1}_r* $path/${1}/.
#fi

#cp fit/commands_${1}.py $path/${1}/.
#cp impacts/${1}/commands_impacts_${1}_v* $path/${1}/.
