## python2 Run3Skimmer.py --input ZZ4lAnalysis.root --output ZZ4lAnalysis_SKIMMED.root --mc
## Drop `--mc` option when running on Data
## Latest files: /eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/240820/
##  /eos/user/m/mmanoni/HZZ_prod_300425_angles/MC/PROD_samplesNano_2022EE_MC_7cebd27b/ggH125/ZZ4lAnalysis.root
#ciao
import ROOT
import os,sys
import optparse

#from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw      

## spencer                                                                                                                                                                                           
def get_genEventSumw(input_file, maxEntriesPerSample=None):

    f = input_file

    runs  = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun) :
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun +=1
    print ("gen=", genEventCount, "sumw=", genEventSumw)

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries>maxEntriesPerSample :
            genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
            nEntries=maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw
## spencer        

usage = ('usage: %prog [options]\n' + '%prog -h for help')
parser = optparse.OptionParser(usage)
parser.add_option('',   '--mc', action='store_true', dest='MC', default=False, help='MC samples')
parser.add_option('',   '--input',dest='INPUT',type='string',default='', help='Name and path to the input file')
parser.add_option('',   '--output',dest='OUTPUT',type='string',default='', help='Name and path to the output file')
(opt, args) = parser.parse_args()


######################################################################## MAKE DATAFRAMES ########################################################################


def makeCR(_df, _flag, year):

    # CRZLLos_3P1F 23 --> 8388608
    # CRZLLos_2P2F 22 --> 4194304
    # CRZLLss 21 --> 2097152

    if _flag == '3P1F':
        bit = '8388608'
        filter = 'ZLLbest3P1FIdx>-1'
        key = "ZLLbest"+_flag+"Idx"
        cand ="ZLLCand"
    elif _flag == '2P2F':
        bit = '4194304'
        filter = 'ZLLbest2P2FIdx>-1'
        key = "ZLLbest"+_flag+"Idx"
        cand ="ZLLCand"
    elif _flag == 'SS':
        bit = '2097152'
        filter = 'ZLLbestSSIdx>-1'
        key = "ZLLbest"+_flag+"Idx"
        cand ="ZLLCand"
    elif _flag == 'SIPCR':
        bit = '21'
        filter = 'ZLLbestSIPCRIdx>-1'
        key = "ZLLbest"+_flag+"Idx"
        cand ="ZLLCand"
    elif _flag == 'SR':
        bit = '0' ## dummy
        filter = 'bestCandIdx>=0'
        key = "bestCandIdx"
        cand ="ZZCand"
    else:
        raise Exception("The CR "+_flag+" is not known")

    df_out = ( _df.Filter(filter)
        .Define('RunNumber', "run")
        .Define('EventNumber', "event")
        .Define('LumiNumber', "luminosityBlock")
        .Define("CRflag", bit)
        .Define("ZZMass", f"{cand}_mass[{key}]")
        .Define("ZZPt",   f"{cand}_pt[{key}]")
        .Define("ZZy",    f"TMath::Abs({cand}_rapidity[{key}])")
        .Define("ZZEta",  f"{cand}_eta[{key}]")
        .Define("ZZPhi",  f"{cand}_phi[{key}]")
        .Define("Z1Flav", f"{cand}_Z1flav[{key}]")
        .Define("Z2Flav", f"{cand}_Z2flav[{key}]")
        .Define("Z1Mass", f"{cand}_Z1mass[{key}]")
        .Define("Z2Mass", f"{cand}_Z2mass[{key}]")
        .Define('Leptons_pt',  "concatenate(Electron_pt,Muon_pt)")
        .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
        .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
        .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
        .Define('Leptons_dz',  "concatenate(Electron_dz,Muon_dz)")
        .Define('Leptons_id',  "concatenate(Electron_pdgId,Muon_pdgId)")
        .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
        .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
        .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_ZZFullId)")
        .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
        .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
        .Define('LepPt',  f"std::vector<float> LepPt{{Leptons_pt[{cand}_Z1l1Idx[{key}]], Leptons_pt[{cand}_Z1l2Idx[{key}]], Leptons_pt[{cand}_Z2l1Idx[{key}]], Leptons_pt[{cand}_Z2l2Idx[{key}]]}}; return LepPt;")
        .Define('LepEta', f"std::vector<float> LepEta{{Leptons_eta[{cand}_Z1l1Idx[{key}]], Leptons_eta[{cand}_Z1l2Idx[{key}]], Leptons_eta[{cand}_Z2l1Idx[{key}]], Leptons_eta[{cand}_Z2l2Idx[{key}]]}}; return LepEta;")
        .Define('LepPhi', f"std::vector<float> LepPhi{{Leptons_phi[{cand}_Z1l1Idx[{key}]], Leptons_phi[{cand}_Z1l2Idx[{key}]], Leptons_phi[{cand}_Z2l1Idx[{key}]], Leptons_phi[{cand}_Z2l2Idx[{key}]]}}; return LepPhi;")
        .Define('Lepdxy', f"std::vector<float> Lepdxy{{Leptons_dxy[{cand}_Z1l1Idx[{key}]], Leptons_dxy[{cand}_Z1l2Idx[{key}]], Leptons_dxy[{cand}_Z2l1Idx[{key}]], Leptons_dxy[{cand}_Z2l2Idx[{key}]]}}; return Lepdxy;")
        .Define('Lepdz',  f"std::vector<float> Lepdz{{Leptons_dz[{cand}_Z1l1Idx[{key}]],  Leptons_dz[{cand}_Z1l2Idx[{key}]],  Leptons_dz[{cand}_Z2l1Idx[{key}]],  Leptons_dz[{cand}_Z2l2Idx[{key}]]}}; return Lepdz;")
        .Define('LepLepId', f"std::vector<short> LepLepId{{Leptons_id[{cand}_Z1l1Idx[{key}]], Leptons_id[{cand}_Z1l2Idx[{key}]], Leptons_id[{cand}_Z2l1Idx[{key}]], Leptons_id[{cand}_Z2l2Idx[{key}]]}}; return LepLepId;")
        .Define('LepisID', f"std::vector<bool> LepisID{{Leptons_isid[{cand}_Z1l1Idx[{key}]], Leptons_isid[{cand}_Z1l2Idx[{key}]], Leptons_isid[{cand}_Z2l1Idx[{key}]], Leptons_isid[{cand}_Z2l2Idx[{key}]]}}; return LepisID;")
        .Define('LepSIP', f"std::vector<float> LepSIP{{Leptons_sip[{cand}_Z1l1Idx[{key}]], Leptons_sip[{cand}_Z1l2Idx[{key}]], Leptons_sip[{cand}_Z2l1Idx[{key}]], Leptons_sip[{cand}_Z2l2Idx[{key}]]}}; return LepSIP;")
        .Define('LepCombRelIsoPF', f"std::vector<float> LepCombRelIsoPF{{Leptons_iso[{cand}_Z1l1Idx[{key}]], Leptons_iso[{cand}_Z1l2Idx[{key}]], Leptons_iso[{cand}_Z2l1Idx[{key}]], Leptons_iso[{cand}_Z2l2Idx[{key}]]}}; return LepCombRelIsoPF;")
        .Define('LepMissingHit', f"std::vector<unsigned char> LepMissingHit{{Leptons_missinghit[{cand}_Z1l1Idx[{key}]], Leptons_missinghit[{cand}_Z1l2Idx[{key}]], Leptons_missinghit[{cand}_Z2l1Idx[{key}]], Leptons_missinghit[{cand}_Z2l2Idx[{key}]]}}; return LepMissingHit;")
        .Define('lep_Hindex', "getHindex(LepPt)")
        .Define('passedFullSelection', "1")
        #.Define('xsec', '1') ## Dummy
        .Define('Leptons_ZZFullSel', "concatenate(Electron_ZZFullSel,Muon_ZZFullSel)")
        .Define('jet_idx', f'get_jet_idx("{year}", Jet_jetId, Jet_ZZMask, Jet_pt, Jet_eta, Jet_phi, Leptons_ZZFullSel, Leptons_eta, Leptons_phi)')
        .Define('pTj1', "jet_idx.size()>0 ? Jet_pt[jet_idx[0]] : -99")
        .Define('Mj1',  "jet_idx.size()>0 ? Jet_mass[jet_idx[0]] : -99")
        .Define('Etaj1',"jet_idx.size()>0 ? Jet_eta[jet_idx[0]] : -99")
        .Define('Phij1',"jet_idx.size()>0 ? Jet_phi[jet_idx[0]] : -99")
        .Define('pTj2', "jet_idx.size()>1 ? Jet_pt[jet_idx[1]] : -99")
        .Define('Mj2',  "jet_idx.size()>1 ? Jet_mass[jet_idx[1]] : -99")
        .Define('Etaj2',"jet_idx.size()>1 ? Jet_eta[jet_idx[1]] : -99")
        .Define('Phij2',"jet_idx.size()>1 ? Jet_phi[jet_idx[1]] : -99")
        .Define('Nj',   "jet_idx.size()")
        .Define('mjj',      "get_mjj(jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('absdetajj',"get_absdetajj(jet_idx, Jet_eta)")
        .Define('dphijj',   "get_dphijj(jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('pTHj',  f"get_pTHj({cand}_pt[{key}], {cand}_eta[{key}], {cand}_phi[{key}], {cand}_mass[{key}], jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('pTHjj', f"get_pTHjj({cand}_pt[{key}], {cand}_eta[{key}], {cand}_phi[{key}], {cand}_mass[{key}], jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('mHj',   f"get_mHj({cand}_pt[{key}], {cand}_eta[{key}], {cand}_phi[{key}], {cand}_mass[{key}], jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('TBjMax', f"get_TBjmax({cand}_rapidity[{key}], jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('TCjMax', f"get_TCjmax({cand}_rapidity[{key}], jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define('costheta1',    f'{cand}_costheta1[{key}]')
        .Define('costheta2',    f'{cand}_costheta2[{key}]')
        .Define('Phi',          f'{cand}_Phi[{key}]')
        .Define('costhetastar', f'{cand}_costhetastar[{key}]')
        .Define('Phi1',         f'{cand}_Phi1[{key}]')
    )

    if  year == "2024":
        df_out = ( df_out.Define('MET', "PFMET_pt") )
    else:
        df_out = ( df_out.Define('MET', "MET_pt") )

    return df_out

def makeGenCR(df, isfail):

    if isfail:
        df = df.Filter("isNotInEvents(event)")

    df_out = (df.Define('GENlep_pt', "getGENlep_vector(FidDressedLeps_pt)")
                    .Define('GENlep_eta', "getGENlep_vector(FidDressedLeps_eta)")
                    .Define('GENlep_phi', "getGENlep_vector(FidDressedLeps_phi)")
                    .Define('GENlep_mass', "getGENlep_vector(FidDressedLeps_mass)")
                    .Define('GENlep_id', "getGENlep_vector(FidDressedLeps_id)")
                    .Define('GENlep_MomId', "getGENlep_vector(FidDressedLeps_momid)")
                    .Define('GENlep_MomMomId', "getGENlep_vector(FidDressedLeps_mommomid)")
                    .Define('GENlep_RelIso', "getGENlep_vector(FidDressedLeps_RelIso)")
                    .Define('GENmass4l', "FidZZ_mass")
                    .Define('GENpT4l', "FidZZ_pt")
                    .Define('GENeta4l', "FidZZ_eta")
                    .Define('GENphi4l', "FidZZ_phi")
                    .Define('GENrapidity4l', "FidZZ_rapidity")
                    .Define('GENmassZ1', "FidZ1_mass")
                    .Define('GENmassZ2', "FidZ2_mass")
                    .Define('GENZ_DaughtersId', "getGENlep_int(FidZ_DauPdgId)")
                    .Define('GENZ_MomId', "getGENlep_int(FidZ_MomPdgId)")
                    .Define('GENlep_Hindex', "getGENHindex(FidZZ_Z1l1Idx, FidZZ_Z1l2Idx, FidZZ_Z2l1Idx, FidZZ_Z2l2Idx)")
                    .Define('genjet_idx', f'get_genjet_idx("{year}", GenJet_pt, GenJet_eta, GenJet_phi, GENlep_eta, GENlep_phi)')
                    .Define('GENpTj1', "genjet_idx.size() > 0 ? GenJet_pt[genjet_idx[0]] : -99")
                    .Define('GENMj1', "genjet_idx.size() > 0 ? GenJet_mass[genjet_idx[0]] : -99")
                    .Define('GENEtaj1', "genjet_idx.size() > 0 ? GenJet_eta[genjet_idx[0]] : -99")
                    .Define('GENPhij1', "genjet_idx.size() > 0 ? GenJet_phi[genjet_idx[0]] : -99")
                    .Define('GENpTj2', "genjet_idx.size() > 1 ? GenJet_pt[genjet_idx[1]] : -99")
                    .Define('GENMj2', "genjet_idx.size() > 1 ? GenJet_mass[genjet_idx[1]] : -99")
                    .Define('GENEtaj2', "genjet_idx.size() > 1 ? GenJet_eta[genjet_idx[1]] : -99")
                    .Define('GENPhij2', "genjet_idx.size() > 1 ? GenJet_phi[genjet_idx[1]] : -99")
                    .Define('GENNj', "genjet_idx.size()")
                    .Define('GENmjj', "get_mjj(genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENabsdetajj', "get_absdetajj(genjet_idx, GenJet_eta)")
                    .Define('GENdphijj', "get_dphijj(genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENpTHj', "get_pTHj(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass, genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENpTHjj', "get_pTHjj(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass, genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENmHj', "get_mHj(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass, genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENTBjMax', "get_TBjmax(FidZZ_rapidity, genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENTCjMax', "get_TCjmax(FidZZ_rapidity, genjet_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")
                    .Define('GENcostheta1', 'FidZZ_costheta1')
                    .Define('GENcostheta2', 'FidZZ_costheta2')
                    .Define('GENPhi', 'FidZZ_Phi')
                    .Define('GENcosthetastar', 'FidZZ_costhetastar')
                    .Define('GENPhi1', 'FidZZ_Phi1')
        )
    
    return df_out

def makeZLCR(df, year):
    
    df_ZL = ( df.Filter('ZLCand_lepIdx>-1').Define("ZZMass", "0") ## Dummy
                                            .Define("CRflag", "0") ## Dummy
                                            .Define("Z1Flav", "ZCand_flav[bestZIdx]")
                                            .Define("Z2Flav", "0") ## Dummy
                                            .Define("Z1Mass", "ZCand_mass[bestZIdx]")
                                            .Define("Z2Mass", "0") ## Dummy
                                            .Define("RunNumber", "run")
                                            .Define("EventNumber", "event")
                                            .Define("LumiNumber", "luminosityBlock")
                                            .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                            .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                            .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                            .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                            .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                            .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                            .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                            .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                            .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_ZZFullId)")
                                            ## Need to add the LepMissingHit branch for SS FR method
                                            ## First create a dummy branch for muons filled with zeroes
                                            .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                            .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                            .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZCand_l1Idx[bestZIdx]], Leptons_pt[ZCand_l2Idx[bestZIdx]], Leptons_pt[ZLCand_lepIdx]}; return LepPt")
                                            .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZCand_l1Idx[bestZIdx]], Leptons_eta[ZCand_l2Idx[bestZIdx]], Leptons_eta[ZLCand_lepIdx]}; return LepEta")
                                            .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZCand_l1Idx[bestZIdx]], Leptons_phi[ZCand_l2Idx[bestZIdx]], Leptons_phi[ZLCand_lepIdx]}; return LepPhi")
                                            .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZCand_l1Idx[bestZIdx]], Leptons_dxy[ZCand_l2Idx[bestZIdx]], Leptons_dxy[ZLCand_lepIdx]}; return Lepdxy")
                                            .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZCand_l1Idx[bestZIdx]], Leptons_dz[ZCand_l2Idx[bestZIdx]], Leptons_dz[ZLCand_lepIdx]}; return Lepdz")
                                            .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZCand_l1Idx[bestZIdx]], Leptons_id[ZCand_l2Idx[bestZIdx]], Leptons_id[ZLCand_lepIdx]}; return LepLepId")
                                            .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZCand_l1Idx[bestZIdx]], Leptons_sip[ZCand_l2Idx[bestZIdx]], Leptons_sip[ZLCand_lepIdx]}; return LepSIP")
                                            .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZCand_l1Idx[bestZIdx]], Leptons_iso[ZCand_l2Idx[bestZIdx]], Leptons_iso[ZLCand_lepIdx]}; return LepCombRelIsoPF")
                                            .Define('LepisID', "std::vector<bool> LepisID{Leptons_isid[ZCand_l1Idx[bestZIdx]], Leptons_isid[ZCand_l2Idx[bestZIdx]], Leptons_isid[ZLCand_lepIdx]}; return LepisID")
                                            .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZCand_l1Idx[bestZIdx]], Leptons_missinghit[ZCand_l2Idx[bestZIdx]], Leptons_missinghit[ZLCand_lepIdx]}; return LepMissingHit")
                                            .Define('L1prefiringWeight', "1") ## Dummy
                                            .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                            .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                            .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                            .Define('xsec', '1') ## Dummy
                                            # not sure if these are necessary
                                            .Define('Leptons_ZZFullSel', "concatenate(Electron_ZZFullSel,Muon_ZZFullSel)")
                                            .Define('jet_idx', f'get_jet_idx("{year}", Jet_jetId, Jet_ZZMask, Jet_pt, Jet_eta, Jet_phi, Leptons_ZZFullSel, Leptons_eta, Leptons_phi)')
                                            .Define('pTj1', "jet_idx.size()>0 ? Jet_pt[jet_idx[0]] : -99")
                                            .Define('Mj1', "jet_idx.size()>0 ? Jet_mass[jet_idx[0]] : -99")
                                            .Define('Etaj1', "jet_idx.size()>0 ? Jet_eta[jet_idx[0]] : -99")
                                            .Define('Phij1', "jet_idx.size()>0 ? Jet_phi[jet_idx[0]] : -99")
                                            .Define('pTj2',  "jet_idx.size()>1 ? Jet_pt[jet_idx[1]] : -99")
                                            .Define('Mj2', "jet_idx.size()>1 ? Jet_mass[jet_idx[1]] : -99")
                                            .Define('Etaj2', "jet_idx.size()>1 ? Jet_eta[jet_idx[1]] : -99")
                                            .Define('Phij2', "jet_idx.size()>1 ? Jet_phi[jet_idx[1]] : -99")
                                            .Define('Nj', "jet_idx.size()")
                                            .Define('mjj', "get_mjj(jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('absdetajj', "get_absdetajj(jet_idx, Jet_eta)")
                                            .Define('dphijj', "get_dphijj(jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('pTHj', "get_pTHj(0, 0, 0, 0, jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('pTHjj', "get_pTHjj(0, 0, 0, 0, jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('mHj', "get_mHj(0, 0, 0, 0, jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('TBjMax', "get_TBjmax(0, jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('TCjMax', "get_TCjmax(0, jet_idx, Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
                                            .Define('costheta1', '0')
                                            .Define('costheta2', '0')
                                            .Define('Phi', '0')
                                            .Define('costhetastar', '0')
                                            .Define('Phi1', '0')
            )

    if  year == "2024":
        df_ZL = ( df_ZL.Define('MET', "PFMET_pt") )
    else:
        df_ZL = ( df_ZL.Define('MET', "MET_pt") )

    return df_ZL

def add_jes(df, jes, year):

    up_idx = f"jet_{jes}_ScaleUp_idx"
    dn_idx = f"jet_{jes}_ScaleDn_idx"

    df = (df
        .Define(up_idx, f'get_jet_idx("{year}", Jet_jetId, Jet_ZZMask, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Leptons_ZZFullSel, Leptons_eta, Leptons_phi)')
        .Define(dn_idx, f'get_jet_idx("{year}", Jet_jetId, Jet_ZZMask, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Leptons_ZZFullSel, Leptons_eta, Leptons_phi)')
        .Define(f"pTj1_{jes}_ScaleUp", f"{up_idx}.size()>0 ? Jet_{jes}_ScaleUp_pt[{up_idx}[0]] : -99")
        .Define(f"pTj1_{jes}_ScaleDn", f"{dn_idx}.size()>0 ? Jet_{jes}_ScaleDn_pt[{dn_idx}[0]] : -99")
        .Define(f"Mj1_{jes}_ScaleUp",  f"{up_idx}.size()>0 ? Jet_{jes}_ScaleUp_mass[{up_idx}[0]] : -99")
        .Define(f"Mj1_{jes}_ScaleDn",  f"{dn_idx}.size()>0 ? Jet_{jes}_ScaleDn_mass[{dn_idx}[0]] : -99")
        .Define(f"pTj2_{jes}_ScaleUp", f"{up_idx}.size()>1 ? Jet_{jes}_ScaleUp_pt[{up_idx}[1]] : -99")
        .Define(f"pTj2_{jes}_ScaleDn", f"{dn_idx}.size()>1 ? Jet_{jes}_ScaleDn_pt[{dn_idx}[1]] : -99")
        .Define(f"Mj2_{jes}_ScaleUp",  f"{up_idx}.size()>1 ? Jet_{jes}_ScaleUp_mass[{up_idx}[1]] : -99")
        .Define(f"Mj2_{jes}_ScaleDn",  f"{dn_idx}.size()>1 ? Jet_{jes}_ScaleDn_mass[{dn_idx}[1]] : -99")
        .Define(f"Nj_{jes}_ScaleUp",   f"{up_idx}.size()")
        .Define(f"Nj_{jes}_ScaleDn",   f"{dn_idx}.size()")
        .Define(f"mjj_{jes}_ScaleUp",  f"get_mjj({up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"mjj_{jes}_ScaleDn",  f"get_mjj({dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"absdetajj_{jes}_ScaleUp", f"get_absdetajj({up_idx}, Jet_eta)")
        .Define(f"absdetajj_{jes}_ScaleDn", f"get_absdetajj({dn_idx}, Jet_eta)")
        .Define(f"dphijj_{jes}_ScaleUp", f"get_dphijj({up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"dphijj_{jes}_ScaleDn", f"get_dphijj({dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"pTHj_{jes}_ScaleUp", f"get_pTHj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"pTHj_{jes}_ScaleDn", f"get_pTHj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"pTHjj_{jes}_ScaleUp", f"get_pTHjj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"pTHjj_{jes}_ScaleDn", f"get_pTHjj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"mHj_{jes}_ScaleUp", f"get_mHj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"mHj_{jes}_ScaleDn", f"get_mHj(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx], {dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"TBjMax_{jes}_ScaleUp", f"get_TBjmax(ZZCand_rapidity[bestCandIdx], {up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"TBjMax_{jes}_ScaleDn", f"get_TBjmax(ZZCand_rapidity[bestCandIdx], {dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
        .Define(f"TCjMax_{jes}_ScaleUp", f"get_TCjmax(ZZCand_rapidity[bestCandIdx], {up_idx}, Jet_{jes}_ScaleUp_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleUp_mass)")
        .Define(f"TCjMax_{jes}_ScaleDn", f"get_TCjmax(ZZCand_rapidity[bestCandIdx], {dn_idx}, Jet_{jes}_ScaleDn_pt, Jet_eta, Jet_phi, Jet_{jes}_ScaleDn_mass)")
    )

    return df


####################################################################### HELPER FUNCTIONS ########################################################################


ROOT.gInterpreter.Declare("""
ROOT::RVec<float> concatenate(ROOT::RVec<float> &A, ROOT::RVec<float> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<short> concatenate(ROOT::RVec<int> &A, ROOT::RVec<int> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<bool> concatenate(ROOT::RVec<bool> &A, ROOT::RVec<bool> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> concatenate(ROOT::RVec<unsigned char> &A, ROOT::RVec<unsigned char> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> addDummyBranch(ROOT::RVec<float> &A){
    int sizeA = A.size();
    ROOT::RVec<unsigned char> B(sizeA);

    for (int i = 0; i < sizeA; ++i) {
        B[i] = 0;
    }

    return B;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getHindex(std::vector<float> LepPt) {
      std::vector<int> _lep_Hindex;
      for(unsigned int i=0;i<LepPt.size();i++){
        _lep_Hindex.push_back(-1);
      }
      float lead_Z1 = max(LepPt.at(0),LepPt.at(1));
      if(lead_Z1 == LepPt.at(0)){
        _lep_Hindex[0] = 0;
        _lep_Hindex[1] = 1;
      }
      else if(lead_Z1 == LepPt.at(1)){
        _lep_Hindex[0] = 1;
        _lep_Hindex[1] = 0;
      }
      float lead_Z2 = max(LepPt.at(2),LepPt.at(3));
      if(lead_Z2 == LepPt.at(2)){
        _lep_Hindex[2] = 2;
        _lep_Hindex[3] = 3;
      }
      else if(lead_Z2 == LepPt.at(3)){
        _lep_Hindex[2] = 3;
        _lep_Hindex[3] = 2;
      }

      return _lep_Hindex;
}

""")

ROOT.gInterpreter.Declare("""
std::vector<int> getGENHindex(int FidZZ_Z1l1Idx, int FidZZ_Z1l2Idx, int FidZZ_Z2l1Idx, int FidZZ_Z2l2Idx) {
    std::vector<int> GENlep_Hindex;
    GENlep_Hindex.push_back(FidZZ_Z1l1Idx);
    GENlep_Hindex.push_back(FidZZ_Z1l2Idx);
    GENlep_Hindex.push_back(FidZZ_Z2l1Idx);
    GENlep_Hindex.push_back(FidZZ_Z2l2Idx);
    return GENlep_Hindex;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getLepGENindex(std::vector<float> LepPt, std::vector<float> LepEta, std::vector<float> LepPhi, ROOT::VecOps::RVec<float> GENlep_id, std::vector<short> LepLepId, ROOT::VecOps::RVec<float> GENlep_pt, ROOT::VecOps::RVec<float> GENlep_eta, ROOT::VecOps::RVec<float> GENlep_phi){
    std::vector<int> lep_genindex;
    for(unsigned int i = 0; i < LepPt.size(); i++){
        lep_genindex.push_back(-1);
    }
    for(unsigned int i = 0; i < LepPt.size(); i++){
        double minDr = 9999.0;
        TLorentzVector reco, gen;
        reco.SetPtEtaPhiM(LepPt.at(i), LepEta.at(i), LepPhi.at(i), 0.0001); 
        for(unsigned int j = 0; j < GENlep_id.size(); j++){
            if (GENlep_id.at(j) != LepLepId.at(i)) continue;
            gen.SetPtEtaPhiM(GENlep_pt.at(j), GENlep_eta.at(j), GENlep_phi.at(j), 0.0001);
            double deta = reco.Eta() - gen.Eta();
            double dphi = abs(reco.Phi() - gen.Phi());
            if (dphi > 3.1415){
                dphi = dphi - 2*3.1415;
            }
            double thisDr = sqrt(deta*deta + dphi*dphi);
            if (thisDr<minDr && thisDr<0.5) {
                lep_genindex[i] = j;
                minDr = thisDr;
            }
        }
    }
    return lep_genindex;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<float> getGENlep_vector(ROOT::VecOps::RVec<float> GENlep_input){
    std::vector<float> GENlep_vector;
    for(unsigned int i = 0; i < GENlep_input.size(); i++){
        GENlep_vector.push_back(GENlep_input.at(i));
    }
    return GENlep_vector;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<int> getGENlep_int(ROOT::VecOps::RVec<int> GENlep_input){
    std::vector<int> GENlep_vector;
    for(unsigned int i = 0; i < GENlep_input.size(); i++){
        GENlep_vector.push_back(GENlep_input.at(i));
    }
    return GENlep_vector;
}
""")

ROOT.gInterpreter.Declare("""
#include <cmath>

bool pass_dR(float eta_1, float phi_1,
             float eta_2, float phi_2,
             float min_dR = 0.4f)
{
    float deta = eta_1 - eta_2;

    float dphi = std::fabs(phi_1 - phi_2);
    if (dphi > M_PI) dphi = 2.0f*M_PI - dphi;

    float dr = std::sqrt(deta*deta + dphi*dphi);
    return dr > min_dR;
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <cmath>
#include <cstddef>

ROOT::VecOps::RVec<int> get_jet_idx(const std::string& year,
                                        const ROOT::VecOps::RVec<UChar_t>& jetId,
                                        const ROOT::VecOps::RVec<Bool_t>& zzMask,
                                        const ROOT::VecOps::RVec<Float_t>& jet_pt,
                                        const ROOT::VecOps::RVec<Float_t>& jet_eta,
                                        const ROOT::VecOps::RVec<Float_t>& jet_phi,
                                        const ROOT::VecOps::RVec<Float_t>& lep_zzfullsel,
                                        const ROOT::VecOps::RVec<Float_t>& lep_eta,
                                        const ROOT::VecOps::RVec<Float_t>& lep_phi,
                                        int requiredId = 6,
                                        float minPt = 30.0,
                                        float maxAbsEta = 5.0)
{
    
    ROOT::VecOps::RVec<int> jet_indices_final;

    const std::size_t njet = jet_eta.size();
    const std::size_t nlep = lep_eta.size();

    for (std::size_t i = 0; i < njet; ++i) {

        if ((int)jetId[i] != requiredId) continue;
        if ((bool)zzMask[i]) continue;
        if (std::abs(jet_eta[i]) >= maxAbsEta) continue;

        // mask checks jep lep energy fraction so dr cut not necessary
        //bool passdR = true;
        //for (std::size_t j = 0; j < nlep; ++j) {
        //    if (!lep_zzfullsel[j]) continue;
        //    if (!pass_dR(jet_eta[i], jet_phi[i], lep_eta[j], lep_phi[j])) {
        //        passdR = false;
        //        break;
        //    }
        //}
        //if (!passdR) continue;

        float pt_cut = minPt;

        // Example: raise pT threshold in endcap region; and in forward region for Run-3 years
        if ((std::abs(jet_eta[i]) >= 2.5f && std::abs(jet_eta[i]) < 3.0f) ||
            (std::abs(jet_eta[i]) >= 3.0f && (year == "2022" || year == "2022EE" || year == "2023" || year == "2023BPix" || year == "2024"))) // TESTING 2024
        {
            pt_cut = 50.0f;
        }

        if (jet_pt[i] > pt_cut) {
            jet_indices_final.push_back(i);
        }
    }

    return jet_indices_final;
}
""")


ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <cmath>
#include <cstddef>

ROOT::VecOps::RVec<int> get_genjet_idx(const std::string& year,
                                        const ROOT::VecOps::RVec<Float_t>& jet_pt,
                                        const ROOT::VecOps::RVec<Float_t>& jet_eta,
                                        const ROOT::VecOps::RVec<Float_t>& jet_phi,
                                        std::vector<float>& lep_eta,
                                        std::vector<float>& lep_phi,
                                        float minPt = 30.0,
                                        float maxAbsEta = 5.0)
{
    
    ROOT::VecOps::RVec<int> genjet_indices_final;

    // Be defensive about sizes
    const std::size_t njet = jet_eta.size();
    const std::size_t nlep = lep_eta.size();

    for (std::size_t i = 0; i < njet; ++i) {

        if (std::abs(jet_eta[i]) >= maxAbsEta) continue;

        bool passdR = true;
        for (std::size_t j = 0; j < nlep; ++j) {
            if (!pass_dR(jet_eta[i], jet_phi[i], lep_eta[j], lep_phi[j])) {
                passdR = false;
                break;
            }
        }
        if (!passdR) continue;

        float pt_cut = minPt;

        // Example: raise pT threshold in endcap region; and in forward region for Run-3 years
        if ((std::abs(jet_eta[i]) >= 2.5f && std::abs(jet_eta[i]) < 3.0f) ||
            (std::abs(jet_eta[i]) >= 3.0f && (year == "2022" || year == "2022EE" || year == "2023" || year == "2023BPix" || year == "2024"))) // TESTING 2024
        {
            pt_cut = 50.0f;
        }

        if (jet_pt[i] > pt_cut) {
            genjet_indices_final.push_back(i);
        }
    }

    return genjet_indices_final;
}
""")


ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <cmath>

float get_mjj(const ROOT::VecOps::RVec<int>& jet_indices,
                    const ROOT::VecOps::RVec<float>& jet_pt,
                    const ROOT::VecOps::RVec<float>& jet_eta,
                    const ROOT::VecOps::RVec<float>& jet_phi,
                    const ROOT::VecOps::RVec<float>& jet_mass,
                    float default_value = -99.0f)
{
    if (jet_indices.size() < 2) {
        return default_value;
    }
    
    int idx1 = jet_indices[0];
    int idx2 = jet_indices[1];
    
    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jet_pt[idx1], jet_eta[idx1], jet_phi[idx1], jet_mass[idx1]);
    jet2.SetPtEtaPhiM(jet_pt[idx2], jet_eta[idx2], jet_phi[idx2], jet_mass[idx2]);
    
    return (jet1 + jet2).M();
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <cmath>

float get_absdetajj(const ROOT::VecOps::RVec<int>& jet_indices,
                         const ROOT::VecOps::RVec<float>& jet_eta,
                         float default_value = -99.0f)
{
    if (jet_indices.size() < 2) {
        return default_value;
    }
    
    int idx1 = jet_indices[0];
    int idx2 = jet_indices[1];
    
    return std::abs(jet_eta[idx1] - jet_eta[idx2]);
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TVector3.h>
#include <cmath>

float get_dphijj(const ROOT::VecOps::RVec<int>& jet_indices,
               const ROOT::VecOps::RVec<float>& jet_pt,
               const ROOT::VecOps::RVec<float>& jet_eta,
               const ROOT::VecOps::RVec<float>& jet_phi,
               const ROOT::VecOps::RVec<float>& jet_mass,
               float default_value = -99.0f)
{
    if (jet_indices.size() < 2) {
        return default_value;
    }

    int idx1 = jet_indices[0];
    int idx2 = jet_indices[1];

    // Build 4-vectors (mass doesn't matter for directions)
    ROOT::Math::PtEtaPhiMVector j1(jet_pt[idx1], jet_eta[idx1], jet_phi[idx1], jet_mass[idx1]);
    ROOT::Math::PtEtaPhiMVector j2(jet_pt[idx2], jet_eta[idx2], jet_phi[idx2], jet_mass[idx2]);

    // Direction of the two jets
    TVector3 j1dir(j1.Px(), j1.Py(), j1.Pz());
    TVector3 j2dir(j2.Px(), j2.Py(), j2.Pz());

    // Transverse components
    TVector3 jt1(j1.Px(), j1.Py(), 0);
    TVector3 jt2(j2.Px(), j2.Py(), 0);

    // Unit vectors
    TVector3 jt1_norm = jt1.Unit();
    TVector3 jt2_norm = jt2.Unit();

    TVector3 z(0,0,1);

    // Cross product sign
    double cross = jt1_norm.Cross(jt2_norm) * z;
    double cross_norm = cross / std::abs(cross);

    // Dot product
    double dot = jt1_norm * jt2_norm;

    // Direction sign
    double diff = (j1dir - j2dir) * z;
    double diff_norm = diff / std::abs(diff);

    return std::acos(dot) * diff_norm * cross_norm;
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TVector3.h>
#include <cmath>

float get_pTHj(float higgs_pt,
                float higgs_eta,
                float higgs_phi,
                float higgs_m,
                const ROOT::VecOps::RVec<int>& jet_indices,
                const ROOT::VecOps::RVec<float>& jet_pt,
                const ROOT::VecOps::RVec<float>& jet_eta,
                const ROOT::VecOps::RVec<float>& jet_phi,
                const ROOT::VecOps::RVec<float>& jet_mass,
                float default_value = -99.0f)
{
    if (jet_indices.size() < 1) {
        return default_value;
    }

    int idx1 = jet_indices[0];

    TLorentzVector H, j1;
    H.SetPtEtaPhiM(higgs_pt, higgs_eta, higgs_phi, higgs_m);
    j1.SetPtEtaPhiM(jet_pt[idx1], jet_eta[idx1], jet_phi[idx1], jet_mass[idx1]);


    return (H+j1).Pt();
}
""")


ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TVector3.h>
#include <cmath>

float get_pTHjj(float higgs_pt,
                float higgs_eta,
                float higgs_phi,
                float higgs_m,
                const ROOT::VecOps::RVec<int>& jet_indices,
                const ROOT::VecOps::RVec<float>& jet_pt,
                const ROOT::VecOps::RVec<float>& jet_eta,
                const ROOT::VecOps::RVec<float>& jet_phi,
                const ROOT::VecOps::RVec<float>& jet_mass,
                float default_value = -99.0f)
{
    if (jet_indices.size() < 2) {
        return default_value;
    }

    int idx1 = jet_indices[0];
    int idx2 = jet_indices[1];

    TLorentzVector H, j1, j2;
    H.SetPtEtaPhiM(higgs_pt, higgs_eta, higgs_phi, higgs_m);
    j1.SetPtEtaPhiM(jet_pt[idx1], jet_eta[idx1], jet_phi[idx1], jet_mass[idx1]);
    j2.SetPtEtaPhiM(jet_pt[idx2], jet_eta[idx2], jet_phi[idx2], jet_mass[idx2]);

    return (H+j1+j2).Pt();
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TVector3.h>
#include <cmath>

float get_mHj(float higgs_pt,
                float higgs_eta,
                float higgs_phi,
                float higgs_m,
                const ROOT::VecOps::RVec<int>& jet_indices,
                const ROOT::VecOps::RVec<float>& jet_pt,
                const ROOT::VecOps::RVec<float>& jet_eta,
                const ROOT::VecOps::RVec<float>& jet_phi,
                const ROOT::VecOps::RVec<float>& jet_mass,
                float default_value = -99.0f)
{
    if (jet_indices.size() < 1) {
        return default_value;
    }

    int idx1 = jet_indices[0];

    TLorentzVector H, j1;
    H.SetPtEtaPhiM(higgs_pt, higgs_eta, higgs_phi, higgs_m);
    j1.SetPtEtaPhiM(jet_pt[idx1], jet_eta[idx1], jet_phi[idx1], jet_mass[idx1]);


    return (H+j1).M();
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <cmath>

float get_TBjmax(float higgs_y,
         const ROOT::VecOps::RVec<int>& jet_indices,
         const ROOT::VecOps::RVec<Float_t>& jet_pt,
         const ROOT::VecOps::RVec<Float_t>& jet_eta,
         const ROOT::VecOps::RVec<Float_t>& jet_phi,
         const ROOT::VecOps::RVec<Float_t>& jet_mass,
         float default_value = -99.0f)
{
    if (jet_indices.size() < 1) {
        return default_value;
    }

    float TBjmax = default_value;

    for (int idx : jet_indices) {

        TLorentzVector j;
        j.SetPtEtaPhiM(jet_pt[idx], jet_eta[idx], jet_phi[idx], jet_mass[idx]);

        float TBj = j.Mt()*std::exp(-std::abs(j.Rapidity() - higgs_y));

        if (TBj > TBjmax) TBjmax = TBj;
    }

    return TBjmax;
}
""")

ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <cmath>

float get_TCjmax(float higgs_y,
         const ROOT::VecOps::RVec<int>& jet_indices,
         const ROOT::VecOps::RVec<Float_t>& jet_pt,
         const ROOT::VecOps::RVec<Float_t>& jet_eta,
         const ROOT::VecOps::RVec<Float_t>& jet_phi,
         const ROOT::VecOps::RVec<Float_t>& jet_mass,
         float default_value = -99.0f)
{
    if (jet_indices.size() < 1) {
        return default_value;
    }

    float TCjmax = default_value;

    for (int idx : jet_indices) {

        TLorentzVector j;
        j.SetPtEtaPhiM(jet_pt[idx], jet_eta[idx], jet_phi[idx], jet_mass[idx]);

        float TCj = std::sqrt(j.E()*j.E() - j.Pz()*j.Pz()) / (2*std::cosh(j.Rapidity() - higgs_y));

        if (TCj > TCjmax) TCjmax = TCj;
    }

    return TCjmax;
}
""")


############################################################################# MAIN ##############################################################################


MC = opt.MC
if MC: 
    print("MC TRUE")

inFileName = opt.INPUT
outFileName = opt.OUTPUT

skip_ZLL = True

if "ZZTo4l" in inFileName: skip_ZL = True
else: skip_ZL = False

do_CR = False
if not opt.MC: do_CR = True
if "DYJets" in inFileName: do_CR = True
if "WZto3LNu" in inFileName: do_CR = True
if "TTto2L2Nu" in inFileName: do_CR = True
if "ZZTo4l" in inFileName: do_CR = True

year = None
if "2022EE" in inFileName: year = "2022EE"
elif "2022" in inFileName: year = "2022"
elif "2023postBPix" in inFileName: year = "2023BPix"
elif "2023" in inFileName: year = "2023"
elif "2024" in inFileName: year = "2024"

jesNames = ['Absolute', f'Absolute_{year}','BBEC1', f'BBEC1_{year}','EC2', f'EC2_{year}','FlavorQCD','HF', f'HF_{year}','RelativeBal',f'RelativeSample_{year}']

df = ROOT.RDataFrame('Events', inFileName)
# df = df.Range(0, 100)

print(f"df nEntries: {df.Count().GetValue()}")

opts = ROOT.RDF.RSnapshotOptions()

df_SR = makeCR(df, "SR", year)

vars_keep = []

if MC:

    vars_keep = ["overallEventWeight", "LHEPdfWeight", "LHEScaleWeight"]
    
    df_SR = ( df_SR.Define("dataMCWeight", "ZZCand_dataMCWeight[bestCandIdx]")
                   .Define('genHEPMCweight', "Generator_weight")
                   .Define('PUWeight', "puWeight"))
    
    if 'ggH' in inFileName: df_SR = df_SR.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")
    if 'ZZTo' in inFileName: df_SR = df_SR.Define('KFactor_QCD_qqZZ_M_weight', "KFactor_QCD_qqZZ_M_Weight")
    if 'ggTo' in inFileName: df_SR = df_SR.Define('KFactor_QCD_ggZZ_Nominal_weight', "KFactor_QCD_ggZZ_Nominal_Weight")

    #for jes in jesNames:
    #    df_SR = add_jes(df_SR, jes, year)

    if "H12" in inFileName:

        vars_keep = vars_keep + ["passedFiducial"]

        events = df.Take[df.GetColumnType('event')]('event')
        event_values = df.AsNumpy(["event"])["event"]
        event_list = list(set(event_values))

        
        ROOT.gInterpreter.Declare(f"""
            #include <vector>
            #include <set>

            std::vector<unsigned long long> getEventVector() {{
                return std::vector<unsigned long long>({{ {', '.join(map(str, event_list))} }});
            }}

            auto eventVector = getEventVector();
            auto eventSet = std::set<unsigned long long>(eventVector.begin(), eventVector.end());
            auto isNotInEvents = [](unsigned long long event) {{
                return eventSet.find(event) == eventSet.end();
            }};
        """)

        df_SR = makeGenCR(df_SR, False)
        df_SR = ( df_SR.Define('lep_genindex', "getLepGENindex(LepPt, LepEta, LepPhi, GENlep_id, LepLepId, GENlep_pt, GENlep_eta, GENlep_phi)") )

        df_all = ROOT.RDataFrame('AllEvents', inFileName)
        print(f"df_all nEntries: {df_all.Count().GetValue()}")
        df_fail = makeGenCR(df_all, True)

        df_fail = ( df_fail.Define('genHEPMCweight', "Generator_weight")
                           .Define('PUWeight', "puWeight")
                           .Define('passedFullSelection', "0")
                           .Define('RunNumber', "run")
                           .Define('EventNumber', "event")
                           .Define('LumiNumber', "luminosityBlock") )

        if 'ggH' in inFileName: df_fail = df_fail.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")

        df_fail.Snapshot('ZZTree/candTree_failed', outFileName, list(df_fail.GetDefinedColumnNames())+vars_keep, opts)

if do_CR:

    if MC:
        vars_keep = ["overallEventWeight", "LHEPdfWeight", "LHEScaleWeight"]

    path = outFileName.rsplit('/', 1)[0] + '/'

    df_3P1F = makeCR(df, "3P1F", year)
    df_2P2F = makeCR(df, "2P2F", year)
    df_2P2Lss = makeCR(df, "SS", year)
    df_SIP = makeCR(df, "SIPCR", year)

    
    opts.fMode = 'RECREATE'
    df_3P1F.Snapshot('CRZLLTree/candTree', path+"test_3P1F.root", list(df_3P1F.GetDefinedColumnNames())+vars_keep, opts)
    df_2P2F.Snapshot('CRZLLTree/candTree', path+"test_2P2F.root", list(df_2P2F.GetDefinedColumnNames())+vars_keep, opts)
    df_2P2Lss.Snapshot('CRZLLTree/candTree', path+"test_2P2Lss.root", list(df_2P2Lss.GetDefinedColumnNames())+vars_keep, opts)
    df_SIP.Snapshot('CRZLLTree/candTree', path+"test_SIP.root", list(df_SIP.GetDefinedColumnNames())+vars_keep, opts)
    

    ## RooDataFrames cannot be concatenated.
    ## Solution: save each CR in a different tree and then merge them through another RooDataFrame
    ## Not fancy, but it works
    df_bis = ROOT.RDataFrame('CRZLLTree/candTree', f"{path}test_*.root")
    print("df_bis.Count().GetValue(): ", df_bis.Count().GetValue())
    if df_bis.Count().GetValue() != 0:
        ## If the number of entries in df_bis is zero (when debugging on few events)
        ## you get an error when trying to take the snapshot
        opts.fMode = 'UPDATE'
        df_bis.Snapshot('CRZLLTree/candTree', outFileName, list(df_3P1F.GetDefinedColumnNames())+vars_keep, opts)
        skip_ZLL = False
    ## Remove the intermediate files, we don't need them anymore
    os.system(f'rm {path}test_*.root')

    if not skip_ZL:
        df_ZL = makeZLCR(df, year)
        opts.fMode = 'UPDATE'
        df_ZL.Snapshot('CRZLTree/candTree', outFileName, list(df_ZL.GetDefinedColumnNames())+vars_keep, opts)

opts.fMode = 'UPDATE'
df_SR.Snapshot('ZZTree/candTree', outFileName, list(df_SR.GetDefinedColumnNames())+vars_keep, opts)

## Add counter only with the 40th entry
counters = ROOT.TH1F("Counters", "Counters", 50, 0, 100)
if opt.MC:
    root = ROOT.TFile.Open(inFileName)
    genEventSumw = get_genEventSumw(root)
    counters.SetBinContent(40, genEventSumw)
    root.Close()

root_file = ROOT.TFile(outFileName, "UPDATE")

if not skip_ZLL:
    sub_dir = root_file.Get("CRZLLTree")
    if sub_dir:
        sub_dir.cd()
        counters.Clone("Counters").Write("", ROOT.TObject.kOverwrite)

if not skip_ZL:
    sub_dir = root_file.Get("CRZLTree")
    if sub_dir:
        sub_dir.cd()
        counters.Clone("Counters").Write("", ROOT.TObject.kOverwrite)

root_file.cd()
counters.Clone("Counters").Write("", ROOT.TObject.kOverwrite)
root_file.Close()