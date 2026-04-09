
import sys
sys.path.append('../helperstuff/')
from split import split

DO_SPLIT = False

#obsNames=["mass4l", "pT4l", "rapidity4l", "massZ1", "massZ2", "pTj1", "pTj2", "Nj", "mjj", "absdetajj", "dphijj", "pTHj", "pTHjj", "mHj", "TCjmax", "TBjmax", "massZ1_massZ2", "rapidity4l_pT4l", "pTj1_pTj2", "Nj_pT4l", "pT4l_pTHj", "absdetajj_mjj", "TCjmax_pT4l", "phi", "phi1", "costhetaZ1", "costhetaZ2", "costhetastar"]   

#obsNames=["mass4l"]
obsNames=["mass4l_zzfloating"]
#obsNames=["pT4l", "rapidity4l", "massZ1", "massZ2"]
#obsNames=["phi", "phi1", "costhetaZ1", "costhetaZ2", "costhetastar"]
#obsNames=["pTj1", "pTj2", "mjj", "absdetajj", "dphijj", "pTHj", "pTHjj", "mHj", "TCjmax", "TBjmax", "Nj"]
#obsNames=["massZ1_massZ2", "rapidity4l_pT4l", "pTj1_pTj2", "pT4l_pTHj", "absdetajj_mjj", "TCjmax_pT4l"] #, "Nj_pT4l"


#obsNames=["pTj2", "mjj", "absdetajj", "dphijj", "pTHjj", "pTj1_pTj2", "absdetajj_mjj", "Nj", "Nj_pT4l"]
#obsNames=["pTj2", "mjj", "absdetajj", "dphijj", "pTHjj", "pTj1_pTj2", "absdetajj_mjj", "Nj"]
#obsNames=["Nj_pT4l"]

#YEARS = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]
#YEARS = ["2024"]
YEARS = ["Run3"]

#YEARS = ["2022", "2022EE"]
#YEARS = ["2022full"] 

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
