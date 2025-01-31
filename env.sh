export base_path=$PWD
echo $base_path;
cd LHScans; mkdir plots;
cd $base_path; 
mkdir combine_files;
mkdir datacard; cd datacard; mkdir datacard_2022; mkdir datacard_2022EE;
mkdir datacard_2023postBPix; mkdir datacard_2023preBPix;
cd $base_path;
mkdir plots;
cd templates; mkdir 2022; mkdir 2023;
cd $base_path;
mkdir scripts;
cd $base_path;
cp ./model/*.py $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/python
