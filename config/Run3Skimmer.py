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

usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
parser = optparse.OptionParser(usage)
parser.add_option('',   '--mc', action='store_true', dest='MC', default=False, help='MC samples')
parser.add_option('',   '--input',dest='INPUT',type='string',default='', help='Name and path to the input file')
parser.add_option('',   '--output',dest='OUTPUT',type='string',default='', help='Name and path to the output file')
parser.add_option('',   '--skipZL', action='store_true', dest='SKIPZL', default=False, help='Skip the ZL tree (e.g., necessary for ZZto4l samples)')

(opt, args) = parser.parse_args()

def makeCR(_df, _flag):

    # CRZLLss 21 --> 2097152
    # CRZLLos_2P2F 22 --> 4194304
    # CRZLLos_3P1F 23 --> 8388608
    if _flag == '3P1F':
        bit = '8388608'
    elif _flag == '2P2F':
        bit = '4194304'
    elif _flag == 'SS':
        bit = '2097152'
    elif _flag == 'SIPCR':
        bit = '21'
    else:
        raise Exception("The CR "+_flag+" is not known")

    df_out = ( _df.Filter('ZLLbest'+_flag+'Idx>-1').Define("ZZMass", "ZLLCand_mass[ZLLbest"+_flag+"Idx]")
                                                    .Define("RunNumber", "run")
                                                    .Define("EventNumber", "event")
                                                    .Define("LumiNumber", "luminosityBlock")
                                                    .Define("CRflag", bit)
                                                    .Define("Z1Mass", "ZLLCand_Z1mass[ZLLbest"+_flag+"Idx]")
                                                    .Define("Z2Mass", "ZLLCand_Z2mass[ZLLbest"+_flag+"Idx]")
                                                    .Define("Z1Flav", "ZLLCand_Z1flav[ZLLbest"+_flag+"Idx]")
                                                    .Define("Z2Flav", "ZLLCand_Z2flav[ZLLbest"+_flag+"Idx]")
                                                    .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                                    .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                                    .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                                    .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                                    .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                                    .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                                    .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                                    .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                                    .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_passID)")
                                                    ## Need to add the LepMissingHit branch for SS FR method
                                                    ## First create a dummy branch for muons filled with zeroes
                                                    .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                                    .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                                    ## Variable miniAOD-style
                                                    .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepPt")
                                                    .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepEta")
                                                    .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepPhi")
                                                    .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Lepdxy")
                                                    .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Lepdz")
                                                    .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepLepId")
                                                    .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepSIP")
                                                    .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepCombRelIsoPF")
                                                    .Define('LepisID', "std::vector<bool> LepisID{Leptons_isid[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Leptons_isid")
                                                    .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepMissingHit")
                                                    #.Define('PFMET', "PFMET_pt")
                                                    ## overallEventWeight contains everything in NanoAODs
                                                    .Define('L1prefiringWeight', "1") ## Dummy
                                                    .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                                    .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                                    .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                                    .Define('xsec', '1') ## Dummy
                                                    # spencer
                                                    .Define('ZZPt', "ZLLCand_pt[ZLLbest"+_flag+"Idx]")
                                                    .Define('ZZy', "ZLLCand_rapidity[ZLLbest"+_flag+"Idx]")
                                                    .Define('pTj1', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? Jet_pt[JetLeadingIdx]: -99")
                                                    .Define('pTj2',  "JetSubleadingIdx >= 0 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetSubleadingIdx]==6 ? Jet_pt[JetSubleadingIdx]: -99")
                                                    .Define('Nj', "nCleanedJetsPt30")
                                                    .Define('mjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).M() : -99")
                                                    .Define('absdetajj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? TMath::Abs(Jet_eta[JetLeadingIdx] - Jet_eta[JetSubleadingIdx]) : -99")
                                                    .Define('dphijj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? deltaphi(ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]), ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])) : -99")
                                                    .Define('pTHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[ZLLbest"+_flag+"Idx], ZZCand_eta[ZLLbest"+_flag+"Idx], ZZCand_phi[ZLLbest"+_flag+"Idx], ZZCand_mass[ZLLbest"+_flag+"Idx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).Pt() : -99")
                                                    .Define('pTHjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[ZLLbest"+_flag+"Idx], ZZCand_eta[ZLLbest"+_flag+"Idx], ZZCand_phi[ZLLbest"+_flag+"Idx], ZZCand_mass[ZLLbest"+_flag+"Idx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).Pt() : -99")
                                                    .Define('mHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[ZLLbest"+_flag+"Idx], ZZCand_eta[ZLLbest"+_flag+"Idx], ZZCand_phi[ZLLbest"+_flag+"Idx], ZZCand_mass[ZLLbest"+_flag+"Idx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).M() : -99")
                                                    .Define('Nj_JESUP', 'nCleanedJetsPt30_jesUp')
                                                    .Define('Nj_JESDOWN', 'nCleanedJetsPt30_jesDn')
                                                    #.Define('TBMax', "Jet_pt.size() > 0 ? TB(ROOT::Math::PtEtaPhiMVector(ZZCand_pt[ZLLbest"+_flag+"Idx], ZZCand_eta[ZLLbest"+_flag+"Idx], ZZCand_phi[ZLLbest"+_flag+"Idx], ZZCand_mass[ZLLbest"+_flag+"Idx]) , Jet_pt, Jet_eta, Jet_phi, Jet_mass) : -99")
                                                    #.Define('TCMax', "Jet_pt.size() > 0 ? TC(ROOT::Math::PtEtaPhiMVector(ZZCand_pt[ZLLbest"+_flag+"Idx], ZZCand_eta[ZLLbest"+_flag+"Idx], ZZCand_phi[ZLLbest"+_flag+"Idx], ZZCand_mass[ZLLbest"+_flag+"Idx]) , Jet_pt, Jet_eta, Jet_phi, Jet_mass) : -99")
                                                    #.Define('Jet_pt_scaleUp', 'Jet_scaleUp_pt')
                                                    #.Define('Jet_pt_scaleDn', 'Jet_scaleDn_pt')
                                                    #.Define('Jet_mass_scaleUp', 'Jet_scaleUp_mass')
                                                    #.Define('Jet_mass_scaleDn', 'Jet_scaleDn_mass')
                                                    #.Define('Jet_pt_smearUp', 'Jet_smearUp_pt')
                                                    #.Define('Jet_pt_smearDn', 'Jet_smearDn_pt')
                                                    #.Define('Jet_mass_smearUp', 'Jet_smearUp_mass')
                                                    #.Define('Jet_mass_smearDn', 'Jet_smearDn_mass')
                                                    # Addition of angluar variables by Martina
                                                    .Define('costheta1', "ZZCand_costheta1[ZLLbest"+_flag+"Idx]")
                                                    .Define('costheta2', "ZZCand_costheta2[ZLLbest"+_flag+"Idx]")
                                                    .Define('Phi', "ZZCand_Phi[ZLLbest"+_flag+"Idx]")
                                                    .Define('costhetastar', "ZZCand_costhetastar[ZLLbest"+_flag+"Idx]")
                                                    .Define('Phi1', "ZZCand_Phi1[ZLLbest"+_flag+"Idx]")
                                                    # tagging
                                                    #.Define('Jet_btagPNetB_', 'Jet_btagPNetB')
                                                )
    return df_out

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

# spencer
ROOT.gInterpreter.Declare("""
float deltaphi (ROOT::Math::PtEtaPhiMVector j1, ROOT::Math::PtEtaPhiMVector j2){ 

  //Direction of the two jets - vectors in the lab frame
  TVector3 j1dir(j1.Px(), j1.Py(), j1.Pz());
  TVector3 j2dir(j2.Px(), j2.Py(), j2.Pz());

  //Transverse component in the xy plane
  TVector3 jt1(j1.Px(), j1.Py(), 0);
  TVector3 jt2(j2.Px(), j2.Py(), 0);

  //Unit vectors of the transverse components
  TVector3 jt1_norm = jt1.Unit();
  TVector3 jt2_norm = jt2.Unit();

  //Unit vector of the z axis
  TVector3 z(0,0,1);

  //Cross product between transverse components
  Double_t cross      = jt1_norm.Cross(jt2_norm) * z;
  Double_t cross_norm = cross / std::abs(cross);

  //Dot product between transverse components
  Double_t dot = jt1_norm * jt2_norm;

  //Difference between the direction of the two jets
  Double_t diff      = (j1dir - j2dir) * z;
  Double_t diff_norm = diff / std::abs(diff);

  return std::acos(dot) * diff_norm * cross_norm;
}
""")

ROOT.gInterpreter.Declare("""
float TB(ROOT::Math::PtEtaPhiMVector H,
         ROOT::VecOps::RVec<Float_t> JetPt,
         ROOT::VecOps::RVec<Float_t> JetEta,
         ROOT::VecOps::RVec<Float_t> JetPhi,
         ROOT::VecOps::RVec<Float_t> JetMass)
{
    float _TBj = -1;
    float _TBjmax = -1;

    for (unsigned int i = 0; i < JetPt.size(); ++i) {
        if (JetPt[i] > 30 && abs(JetEta[i]) < 4.7) {
            TLorentzVector theJet;
            theJet.SetPtEtaPhiM(JetPt[i], JetEta[i], JetPhi[i], JetMass[i]);
            _TBj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2)) * exp(-1 * abs(theJet.Rapidity() - H.Rapidity()));
            if (_TBj > _TBjmax) _TBjmax = _TBj;
        }
    }
    return _TBjmax;
}
""")

ROOT.gInterpreter.Declare("""
float TC(ROOT::Math::PtEtaPhiMVector H,
         ROOT::VecOps::RVec<Float_t> JetPt,
         ROOT::VecOps::RVec<Float_t> JetEta,
         ROOT::VecOps::RVec<Float_t> JetPhi,
         ROOT::VecOps::RVec<Float_t> JetMass)
{
    float _TCj = -1;
    float _TCjmax = -1;

    for (unsigned int i = 0; i < JetPt.size(); ++i) {
        if (JetPt[i] > 30 && abs(JetEta[i]) < 4.7) {
            TLorentzVector theJet;
            theJet.SetPtEtaPhiM(JetPt[i], JetEta[i], JetPhi[i], JetMass[i]);
            _TCj = sqrt(pow(theJet.Pt(), 2) + pow(theJet.M(), 2)) / (2 * cosh(theJet.Rapidity() - H.Rapidity()));
            if (_TCj > _TCjmax) _TCjmax = _TCj;
        }
    }
    return _TCjmax;
}
""")

ROOT.gInterpreter.Declare("""
TLorentzVector makeP4(float pt, float eta, float phi, float mass) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, mass);
    return v;
}
""")

# spencer


##################################### MAIN #####################################

_SKIPZL = opt.SKIPZL
MC = opt.MC
if MC: 
    print("MC TRUE")
    _SKIPZL = True

inFileName = opt.INPUT
outFileName = opt.OUTPUT

df = ROOT.RDataFrame('Events', inFileName)
# df = df.Range(0, 100)

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

if not MC:
    df_3P1F = makeCR(df, "3P1F")
    df_2P2F = makeCR(df, "2P2F")
    df_2P2Lss = makeCR(df, "SS")
    df_SIP = makeCR(df, "SIPCR")
    
#print(f"df nEntries: {df.Count().GetValue()}")

if "H12" in inFileName:
    df_all = ROOT.RDataFrame('AllEvents', inFileName)
    print(f"df_all nEntries: {df_all.Count().GetValue()}")

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = 'RECREATE'

## Variables to store in the output root file
vars = {'RunNumber',
        'EventNumber',
        'LumiNumber',
        'ZZMass',
        'ZZPt',
        'ZZy',
        # 'CRflag',
        'Z1Flav',
        'Z2Flav',
        'Z1Mass',
        'Z2Mass',
        #'dataMCWeight',
        #'PFMET',
        'LepPt',
        'LepEta',
        'LepPhi',
        'Lepdxy',
        'Lepdz',
        'LepLepId',
        'LepSIP',
        'LepCombRelIsoPF',
        # 'LepMissingHit',
        # spencer
        'pTj1', 
        'pTj2',
        'Nj',
        'mjj',
        'absdetajj',
        'dphijj',  
        'pTHj',
        'pTHjj',
        'mHj', 
        'Nj_JESUP',
        'Nj_JESDOWN',
        #'TBMax',
        #'TCMax',
        #'Jet_pt_scaleUp',
        #'Jet_pt_scaleDn',
        #'Jet_mass_scaleUp',
        #'Jet_mass_scaleDn',
        #'Jet_pt_smearUp',
        #'Jet_pt_smearDn',
        #'Jet_mass_smearUp',
        #'Jet_mass_smearDn',
        # Addition of angluar variables by Martina
        'costheta1',
        'costheta2',
        'Phi',
        'costhetastar',
        'Phi1',
        #'Jet_btagPNetB_',
        }
if MC:
    vars.add('overallEventWeight')
    vars.add('LHEPdfWeight_') # spencer
    vars.add('LHEScaleWeight_') # spencer
    vars.add('PUWeight')
    # vars.add('LHEWeight_originalXWGTUP')
    vars.add('lep_Hindex') # SHOULD THIS BE GENLEP_HINDEX ? - SPENCER
    vars.add('dataMCWeight') # spencer
    if "H12" in inFileName:
        vars.add('GENlep_pt')
        vars.add('GENlep_eta')
        vars.add('GENlep_phi')
        vars.add('GENlep_mass')
        vars.add('GENlep_id')
        vars.add('GENlep_MomId')
        vars.add('GENlep_MomMomId')
        vars.add('GENlep_RelIso')
        vars.add('GENmass4l')
        vars.add('GENpT4l')
        vars.add('GENeta4l')
        vars.add('GENphi4l')
        vars.add('GENrapidity4l')
        vars.add('GENZ_DaughtersId')
        vars.add('GENZ_MomId')
        vars.add('GENlep_Hindex')
        # vars.add('GENlep_index')
        vars.add('lep_genindex')
        vars.add('passedFiducial')
        # spencer
        vars.add('GENpTj1')
        vars.add('GENpTj2')
        vars.add('GENNj')
        vars.add('GENmjj')
        vars.add('GENabsdetajj')
        vars.add('GENdphijj')
        vars.add('GENpTHj')
        vars.add('GENpTHjj')
        vars.add('GENmHj')
        vars.add('GENmassZ1')
        vars.add('GENmassZ2')
        #ars.add('GENTBMax')
        #vars.add('GENTCMax')
        # Addition of angluar varaibles by Martina
        vars.add('GENcostheta1')
        vars.add('GENcostheta2')
        vars.add('GENPhi')
        vars.add('GENcosthetastar')
        vars.add('GENPhi1')
        vars.add('Jet_hadronFlavour_')

    vars.add('passedFullSelection')
    vars.add('genHEPMCweight')
    
    if 'ggH' in inFileName:
        vars.add('ggH_NNLOPS_weight')
    if 'ZZTo' in inFileName:
        # vars.add('KFactor_EW_qqZZ_Weight')
        vars.add('KFactor_QCD_qqZZ_M_Weight')
    if 'ggTo' in inFileName:
        vars.add('KFactor_QCD_ggZZ_Nominal_Weight')

skip_ZLL = True
if not MC:
    df_3P1F.Snapshot('CRZLLTree/candTree', "test_3P1F.root", vars, opts)
    df_2P2F.Snapshot('CRZLLTree/candTree', "test_2P2F.root", vars, opts)
    df_2P2Lss.Snapshot('CRZLLTree/candTree', "test_2P2Lss.root", vars, opts)
    df_SIP.Snapshot('CRZLLTree/candTree', "test_SIP.root", vars, opts)

    ## RooDataFrames cannot be concatenated.
    ## Solution: save each CR in a different tree and then merge them through another RooDataFrame
    ## Not fancy, but it works
    df_bis = ROOT.RDataFrame('CRZLLTree/candTree', "test_*.root")
    skip_ZLL = True
    if df_bis.Count().GetValue() != 0:
        ## If the number of entries in df_bis is zero (when debugging on few events)
        ## you get an error when trying to take the snapshot
        df_bis.Snapshot('CRZLLTree/candTree', outFileName, vars, opts)
        skip_ZLL = False
    ## Remove the intermediate files, we don't need them anymore
    os.system('rm test_*.root')
        
## SR
df_SR = ( df.Filter('bestCandIdx>=0').Define("ZZMass", "ZZCand_mass[bestCandIdx]") ## Dummy
                                    .Define("ZZPt", "ZZCand_pt[bestCandIdx]")
                                    .Define("ZZy", "abs(ZZCand_rapidity[bestCandIdx])")
                                    #.Define("CRflag", "0") ## Dummy
                                    .Define("Z1Flav", "ZZCand_Z1flav[bestCandIdx]")
                                    .Define("Z2Flav", "ZZCand_Z2flav[bestCandIdx]") ## Dummy
                                    .Define("Z1Mass", "ZZCand_Z1mass[bestCandIdx]")
                                    .Define("Z2Mass", "ZZCand_Z2mass[bestCandIdx]") ## Dummy
                                    #.Define("dataMCWeight", "ZZCand_dataMCWeight[bestCandIdx]")
                                    .Define('RunNumber', "run")
                                    .Define('EventNumber', "event")
                                    .Define('LumiNumber', "luminosityBlock")
                                    .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                    .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                    .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                    .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                    .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                    .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                    .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                    .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                    ## Need to add the LepMissingHit branch for SS FR method
                                    ## First create a dummy branch for muons filled with zeroes
                                    # .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                    # .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                    ## Variable miniAOD-style
                                    .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPt")
                                    .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepEta")
                                    .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPhi")
                                    .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdxy")
                                    .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdz")
                                    .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepLepId")
                                    .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepSIP")
                                    .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepCombRelIsoPF")
                                    # .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepMissingHit")
                                    #.Define('PFMET', "PFMET_pt")
                                    .Define('lep_Hindex', "getHindex(LepPt)")
                                    .Define('passedFullSelection', "1")
                                    # spencer
                                    .Define('pTj1', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? Jet_pt[JetLeadingIdx]: -99")
                                    .Define('pTj2',  "JetSubleadingIdx >= 0 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetSubleadingIdx]==6 ? Jet_pt[JetSubleadingIdx]: -99")
                                    .Define('Nj', "nCleanedJetsPt30")
                                    .Define('mjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).M() : -99")
                                    .Define('absdetajj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? TMath::Abs(Jet_eta[JetLeadingIdx] - Jet_eta[JetSubleadingIdx]) : -99")
                                    .Define('dphijj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? deltaphi(ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]), ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])) : -99")
                                    .Define('pTHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).Pt() : -99")
                                    .Define('pTHjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).Pt() : -99")
                                    .Define('mHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).M() : -99")
                                    .Define('Nj_JESUP', 'nCleanedJetsPt30_jesUp')
                                    .Define('Nj_JESDOWN', 'nCleanedJetsPt30_jesDn')
                                    #.Define('TBMax', "Jet_pt.size() > 0 ? TB(ROOT::Math::PtEtaPhiMVector(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx]) , Jet_pt, Jet_eta, Jet_phi, Jet_mass) : -99")
                                    #.Define('TCMax', "Jet_pt.size() > 0 ? TC(ROOT::Math::PtEtaPhiMVector(ZZCand_pt[bestCandIdx], ZZCand_eta[bestCandIdx], ZZCand_phi[bestCandIdx], ZZCand_mass[bestCandIdx]) , Jet_pt, Jet_eta, Jet_phi, Jet_mass) : -99")
                                    #.Define('Jet_pt_scaleUp', 'Jet_scaleUp_pt')
                                    #.Define('Jet_pt_scaleDn', 'Jet_scaleDn_pt')
                                    #.Define('Jet_mass_scaleUp', 'Jet_scaleUp_mass')
                                    #.Define('Jet_mass_scaleDn', 'Jet_scaleDn_mass')
                                    #.Define('Jet_pt_smearUp', 'Jet_smearUp_pt')
                                    #.Define('Jet_pt_smearDn', 'Jet_smearDn_pt')
                                    #.Define('Jet_mass_smearUp', 'Jet_smearUp_mass')
                                    #.Define('Jet_mass_smearDn', 'Jet_smearDn_mass')
                                    # Addition of angluar variables by Martina
                                    .Define('costheta1', 'ZZCand_costheta1[bestCandIdx]')
                                    .Define('costheta2', 'ZZCand_costheta2[bestCandIdx]')
                                    .Define('Phi', 'ZZCand_Phi[bestCandIdx]')
                                    .Define('costhetastar', 'ZZCand_costhetastar[bestCandIdx]')
                                    .Define('Phi1', 'ZZCand_Phi1[bestCandIdx]')
                                    # tagging
                                    #.Define('Jet_btagPNetB_', 'Jet_btagPNetB')
    )

if MC:
    df_SR = (df_SR.Define('LHEPdfWeight_', "LHEPdfWeight") # spencer
         .Define('LHEScaleWeight_', "LHEScaleWeight") # spencer
         .Define("dataMCWeight", "ZZCand_dataMCWeight[bestCandIdx]") # spencer
         .Define('genHEPMCweight', "Generator_weight") # spencer
         .Define('PUWeight', "puWeight")) # spencer

    
if 'ggH' in inFileName: df_SR = df_SR.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")
if 'ZZTo' in inFileName: df_SR = df_SR.Define('KFactor_QCD_qqZZ_M_weight', "KFactor_QCD_qqZZ_M_Weight") # spencer
if 'ggTo' in inFileName: df_SR = df_SR.Define('KFactor_QCD_ggZZ_Nominal_weight', "KFactor_QCD_ggZZ_Nominal_Weight") # spencer

if "H12" in inFileName:
    df_SR = (df_SR.Define('GENlep_pt', "getGENlep_vector(FidDressedLeps_pt)")
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
                    .Define('GENZ_DaughtersId', "getGENlep_int(FidZ_DauPdgId)")
                    .Define('GENZ_MomId', "getGENlep_int(FidZ_MomPdgId)")
                    .Define('GENlep_Hindex', "getGENHindex(FidZZ_Z1l1Idx, FidZZ_Z1l2Idx, FidZZ_Z2l1Idx, FidZZ_Z2l2Idx)")
                    .Define('lep_genindex', "getLepGENindex(LepPt, LepEta, LepPhi, GENlep_id, LepLepId, GENlep_pt, GENlep_eta, GENlep_phi)")
                    # spencer
                    .Define('GENpTj1', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? GenJet_pt[FidZZ_GenJetLeadingIdx] : -99")
                    .Define('GENpTj2', "FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? GenJet_pt[FidZZ_GenJetSubleadingIdx] : -99")
                    .Define('GENNj', "FidZZ_nCleanedGenJetsPt30")
                    .Define('GENmjj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])).M() : -99")
                    .Define('GENabsdetajj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0  && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx] - GenJet_eta[FidZZ_GenJetSubleadingIdx]) : -99")
                    .Define('GENdphijj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? deltaphi(ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]), ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])) : -99")
                    .Define('GENpTHj', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx])).Pt() : -99")
                    .Define('GENpTHjj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])).Pt() : -99")
                    .Define('GENmHj', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx])).M() : -99")
                    .Define('GENmassZ1', "FidZ1_mass")
                    .Define('GENmassZ2', "FidZ2_mass")
                    #.Define('GENTBMax', "GenJet_pt.size() > 0 ? TB(ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) , GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass) : -99")
                    #Define('GENTCMax', "GenJet_pt.size() > 0 ? TC(ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) , GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass) : -99")
                    # Addition of angluar variables by Martina
                    .Define('GENcostheta1', 'FidZZ_costheta1')
                    .Define('GENcostheta2', 'FidZZ_costheta2')
                    .Define('GENPhi', 'FidZZ_Phi')
                    .Define('GENcosthetastar', 'FidZZ_costhetastar')
                    .Define('GENPhi1', 'FidZZ_Phi1')
                    .Define('Jet_hadronFlavour_', 'Jet_hadronFlavour')
            )
    
    df_filter = df_all.Filter("isNotInEvents(event)")
    print(f"df_filter nEntries: {df_filter.Count().GetValue()}")
    vars_fail = {'GENlep_pt',
    'GENlep_eta',
    'GENlep_phi',
    'GENlep_mass',
    'GENlep_id',
    'GENlep_MomId',
    'GENlep_MomMomId',
    'GENlep_RelIso',
    'GENmass4l',
    'GENpT4l',
    'GENeta4l',
    'GENphi4l',
    'GENrapidity4l',
    'GENZ_DaughtersId',
    'GENZ_MomId',
    'GENlep_Hindex',
    'RunNumber',
    'EventNumber',
    'LumiNumber',
    'passedFiducial',
    'passedFullSelection',
    'LHEPdfWeight_', # spencer
    'LHEScaleWeight_', # spencer
    'genHEPMCweight',
    'PUWeight',
    'GENpTj1', # spencer
    'GENpTj2', 
    'GENNj',
    'GENmjj',
    'GENabsdetajj',
    'GENdphijj',
    'GENpTHj',
    'GENpTHjj', 
    'GENmHj', 
    'GENmassZ1',
    'GENmassZ2',
    #'GENTBMax',
    #'GENTCMax',
    'GENcostheta1', # Martina
    'GENcostheta2',
    'GENPhi',
    'GENcosthetastar',
    'GENPhi1',}

    df_fail = (df_filter.Define('GENlep_pt', "getGENlep_vector(FidDressedLeps_pt)")
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
                        .Define('GENZ_DaughtersId', "getGENlep_int(FidZ_DauPdgId)")
                        .Define('GENZ_MomId', "getGENlep_int(FidZ_MomPdgId)")
                        .Define('GENlep_Hindex', "getGENHindex(FidZZ_Z1l1Idx, FidZZ_Z1l2Idx, FidZZ_Z2l1Idx, FidZZ_Z2l2Idx)")
                        .Define('RunNumber', "run")
                        .Define('EventNumber', "event")
                        .Define('LumiNumber', "luminosityBlock")
                        .Define('passedFullSelection', "0")
                        .Define('genHEPMCweight', "Generator_weight")  #LHEWeight_originalXWGTUP")
                        .Define('PUWeight', "puWeight")
                        .Define('LHEPdfWeight_', "LHEPdfWeight") # spencer
                        .Define('LHEScaleWeight_', "LHEScaleWeight") # spencer    
                        # spencer
                        .Define('GENpTj1', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? GenJet_pt[FidZZ_GenJetLeadingIdx] : -99")
                        .Define('GENpTj2', "FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? GenJet_pt[FidZZ_GenJetSubleadingIdx] : -99")
                        .Define('GENNj', "FidZZ_nCleanedGenJetsPt30")
                        .Define('GENmjj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])).M() : -99")
                        .Define('GENabsdetajj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0  && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx] - GenJet_eta[FidZZ_GenJetSubleadingIdx]) : -99")
                        .Define('GENdphijj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? deltaphi(ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]), ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])) : -99")
                        .Define('GENpTHj', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx])).Pt() : -99")
                        .Define('GENpTHjj', "FidZZ_GenJetLeadingIdx >= 0 && FidZZ_GenJetSubleadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 && GenJet_pt[FidZZ_GenJetSubleadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetSubleadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetSubleadingIdx], GenJet_eta[FidZZ_GenJetSubleadingIdx], GenJet_phi[FidZZ_GenJetSubleadingIdx], GenJet_mass[FidZZ_GenJetSubleadingIdx])).Pt() : -99")
                        .Define('GENmHj', "FidZZ_GenJetLeadingIdx >= 0 && GenJet_pt[FidZZ_GenJetLeadingIdx]>30 && TMath::Abs(GenJet_eta[FidZZ_GenJetLeadingIdx])<5 ? (ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) + ROOT::Math::PtEtaPhiMVector(GenJet_pt[FidZZ_GenJetLeadingIdx], GenJet_eta[FidZZ_GenJetLeadingIdx], GenJet_phi[FidZZ_GenJetLeadingIdx], GenJet_mass[FidZZ_GenJetLeadingIdx])).M() : -99")
                        #.Define('GENTBMax', "GenJet_pt.size() > 0 ? TB(ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) , GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass) : -99")
                        #.Define('GENTCMax', "GenJet_pt.size() > 0 ? TC(ROOT::Math::PtEtaPhiMVector(FidZZ_pt, FidZZ_eta, FidZZ_phi, FidZZ_mass) , GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass) : -99")
                        .Define('GENmassZ1', "FidZ1_mass")
                        .Define('GENmassZ2', "FidZ2_mass")
                        # Addition of angluar variables by Martina
                        .Define('GENcostheta1', 'FidZZ_costheta1')
                        .Define('GENcostheta2', 'FidZZ_costheta2')
                        .Define('GENPhi', 'FidZZ_Phi')
                        .Define('GENcosthetastar', 'FidZZ_costhetastar')
                        .Define('GENPhi1', 'FidZZ_Phi1')
    )

    if 'ggH' in inFileName:
        vars_fail.add('ggH_NNLOPS_weight')
        df_fail = df_fail.Define('ggH_NNLOPS_weight', "ggH_NNLOPS_Weight")

opts.fMode = 'UPDATE'
df_SR.Snapshot('ZZTree/candTree', outFileName, vars, opts)


## ZL CR for the computation of fake rates
if not _SKIPZL:
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
                                           .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_passID)")
                                           ## Need to add the LepMissingHit branch for SS FR method
                                           ## First create a dummy branch for muons filled with zeroes
                                           .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                           .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                           ## Variable miniAOD-style
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
                                           #.Define('PFMET', "PFMET_pt")
                                           ## overallEventWeight contains everything in NanoAODs
                                           .Define('L1prefiringWeight', "1") ## Dummy
                                           .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                           .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                           .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                           .Define('xsec', '1') ## Dummy
                                            # spencer
                                            .Define('ZZPt', "0")
                                            .Define('ZZy', "0")
                                            .Define('pTj1', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? Jet_pt[JetLeadingIdx]: -99")
                                            .Define('pTj2',  "JetSubleadingIdx >= 0 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetSubleadingIdx]==6 ? Jet_pt[JetSubleadingIdx]: -99")
                                            .Define('Nj', "nCleanedJetsPt30")
                                            .Define('mjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).M() : -99")
                                            .Define('absdetajj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? TMath::Abs(Jet_eta[JetLeadingIdx] - Jet_eta[JetSubleadingIdx]) : -99")
                                            .Define('dphijj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? deltaphi(ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]), ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])) : -99")
                                            .Define('pTHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).Pt() : -99")
                                            .Define('pTHjj', "JetLeadingIdx >= 0 && JetSubleadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_pt[JetSubleadingIdx]>30 && TMath::Abs(Jet_eta[JetSubleadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 && Jet_jetId[JetSubleadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx]) + ROOT::Math::PtEtaPhiMVector(Jet_pt[JetSubleadingIdx], Jet_eta[JetSubleadingIdx], Jet_phi[JetSubleadingIdx], Jet_mass[JetSubleadingIdx])).Pt() : -99")
                                            .Define('mHj', "JetLeadingIdx >= 0 && Jet_pt[JetLeadingIdx]>30 && TMath::Abs(Jet_eta[JetLeadingIdx])<5 && Jet_jetId[JetLeadingIdx]==6 ? (ROOT::Math::PtEtaPhiMVector(Jet_pt[JetLeadingIdx], Jet_eta[JetLeadingIdx], Jet_phi[JetLeadingIdx], Jet_mass[JetLeadingIdx])).M() : -99")
                                            .Define('Nj_JESUP', 'nCleanedJetsPt30_jesUp')
                                            .Define('Nj_JESDOWN', 'nCleanedJetsPt30_jesDn')
                                            #.Define('TBMax', "Jet_pt.size() > 0 ? TB(0 , Jet_pt, Jet_eta, Jet_phi, Jet_mass): -99")
                                            #.Define('TCMax', "Jet_pt.size() > 0 ? TC(0 , Jet_pt, Jet_eta, Jet_phi, Jet_mass): -99")
                                            #.Define('Jet_pt_scaleUp', 'Jet_scaleUp_pt')
                                            #.Define('Jet_pt_scaleDn', 'Jet_scaleDn_pt')
                                            #.Define('Jet_mass_scaleUp', 'Jet_scaleUp_mass')
                                            #.Define('Jet_mass_scaleDn', 'Jet_scaleDn_mass')
                                            #.Define('Jet_pt_smearUp', 'Jet_smearUp_pt')
                                            #.Define('Jet_pt_smearDn', 'Jet_smearDn_pt')
                                            #.Define('Jet_mass_smearUp', 'Jet_smearUp_mass')
                                            #.Define('Jet_mass_smearDn', 'Jet_smearDn_mass')
                                            # Addition of angluar variables by Martina
                                            .Define('costheta1', 'ZZCand_costheta1[bestCandIdx]')
                                            .Define('costheta2', 'ZZCand_costheta2[bestCandIdx]')
                                            .Define('Phi', 'ZZCand_Phi[bestCandIdx]')
                                            .Define('costhetastar', 'ZZCand_costhetastar[bestCandIdx]')
                                            .Define('Phi1', 'ZZCand_Phi1[bestCandIdx]')
                                            # tagging
                                            #.Define('Jet_btagPNetB_', 'Jet_btagPNetB')
                                            )
    opts.fMode = 'UPDATE'
    df_ZL.Snapshot('CRZLTree/candTree', outFileName, vars, opts)


opts.fMode = 'UPDATE'
print(f"df nEntries: {df_SR.Count().GetValue()}")

df_SR.Snapshot('ZZTree/candTree', outFileName, vars, opts)
if "H12" in inFileName:
    opts.fMode = 'UPDATE'
    df_fail.Snapshot('ZZTree/candTree_failed', outFileName, vars_fail, opts)
#'''
## Add counter only with the 40th entry
counters = ROOT.TH1F("Counters", "Counters", 50, 0, 100)
if opt.MC:
    root = ROOT.TFile.Open(inFileName)
    genEventSumw = get_genEventSumw(root)
    counters.SetBinContent(40, genEventSumw)
    counters.Write()
    print("Counter belongs to:", counters.GetDirectory().GetName())
    root.Close()
#'''

if opt.MC:
    root = ROOT.TFile.Open(inFileName)
    genEventSumw = get_genEventSumw(root)
    root.Close()

    # Open output file in UPDATE mode
    root_file = ROOT.TFile(outFileName, "UPDATE")
    root_file.cd()  # Ensure you're in the file's directory

    # Create and write histogram
    counters = ROOT.TH1F("Counters", "Counters", 50, 0, 100)
    counters.SetBinContent(40, genEventSumw)
    counters.SetDirectory(root_file)  # Attach it to the file
    counters.Write("", ROOT.TObject.kOverwrite)  # Force write

    print("Counter belongs to:", counters.GetDirectory().GetName())
    root_file.Close()

root_file = ROOT.TFile(outFileName, "UPDATE")
#'''
if not skip_ZLL:
    sub_dir = root_file.Get("CRZLLTree")
    sub_dir.cd()
    counters.Write()

if not _SKIPZL:
    sub_dir = root_file.Get("CRZLTree")
    sub_dir.cd()
    counters.Write()
#'''
root_file.Close()
