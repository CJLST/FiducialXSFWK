import sys
sys.path.append('../helperstuff/')
from split import split

DO_SPLIT = False

#obsNames=["mass4l", "pT4l", "rapidity4l", "massZ1", "massZ2", "pTj1", "pTj2", "Nj", "mjj", "absdetajj", "dphijj", "pTHj", "pTHjj", "mHj", "massZ1_massZ2", "rapidity4l_pT4l", "pTj1_pTj2", "Nj_pT4l", "pT4l_pTHj", "absdetajj_mjj", "phi", "phi1", "costhetaZ1", "costhetaZ2","costhetastar"]   
#obsNames=["pT4l", "rapidity4l", "massZ1", "massZ2", "pTj1", "pTj2", "Nj", "mjj", "absdetajj", "dphijj", "pTHj", "pTHjj", "mHj", "phi", "phi1", "costhetaZ1", "costhetaZ2", "costhetastar"]

#obsNames=[ "massZ1", "massZ2", "phi", "phi1", "costhetaZ1", "costhetaZ2", "costhetastar"]

#obsNames=["massZ1_massZ2", "rapidity4l_pT4l", "pTj1_pTj2", "Nj_pT4l", "pT4l_pTHj", "absdetajj_mjj"]
#obsNames=["mass4l"]

obsNames=["mass4l", "pT4l", "rapidity4l", "pTj1"]

YEARS = ["2022", "2022EE"]#, "2023preBPix", "2023postBPix", "2024"]
#YEARS = ["Run3"]

#YEARS = ["2022", "2022EE", "2023preBPix", "2023postBPix"]
#YEARS = ["2022_2023"] 

#YEARS = ["2022"]
#YEARS = ["2022EE"]
#YEARS = ["2023preBPix"]
#YEARS = ["2023postBPix"]

years = []


if DO_SPLIT:
    for YEAR in YEARS:
        if split[YEAR]:
            nSplit = split[YEAR]
            for i in range(0,nSplit):
                years.append(YEAR+"_"+str(i+1))
else:
    years = YEARS
    
with open("input_args.txt", "w") as f:
    for obs in obsNames:
        for year in years:
            f.write(f'{obs} {year}\n')
