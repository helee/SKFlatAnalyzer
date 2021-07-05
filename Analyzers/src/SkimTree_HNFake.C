#include "SkimTree_HNFake.h"

void SkimTree_HNFake::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNFake::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  //triggers.clear();
  triggers_mu.clear();
  triggers_el.clear();
  if(DataYear==2016){
    triggers_mu = {
      "HLT_Mu3_PFJet40_v",                             // DoubleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
    };
    triggers_el = {
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // DoubleEG
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
      "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // DoubleEG
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // DoubleEG
    };
  }else if(DataYear==2017){
    triggers_mu = {
      "HLT_Mu3_PFJet40_v",                             // SingleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
    };
    triggers_el = {
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // SingleElectron
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // SingleElectron
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // SingleElectron
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // SingleElectron
    };
  }else if(DataYear==2018){
    triggers_mu = {
      "HLT_Mu3_PFJet40_v",                             // SingleMuon
      "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
      "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
    };
    triggers_el = {
      "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // EGamma
      "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // EGamma
      "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",          // EGamma
      "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // EGamma
    };
  }else{
    cout<<"[SkimTree_HNFake::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNFake::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HNFake::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

}

void SkimTree_HNFake::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if(!(ev.PassTrigger(triggers_mu) || ev.PassTrigger(triggers_el))) return;

  //std::vector<Muon> muons_veto = GetMuons("ISRVeto", 3., 2.4);
  //std::vector<Electron> electrons_veto = GetElectrons("ISRVeto", 8., 2.5);
  //std::vector<Electron> electrons_veto_IsoUp = GetElectrons("ISRVetoIsoUp", 8., 2.5);
  std::vector<Muon> muons_loose = GetMuons("ISRLooseIsoUp", 3., 2.4);
  std::vector<Electron> electrons_loose = GetElectrons("ISRLooseIsoUp", 8., 2.5);
  std::vector<Jet> jets = GetJets("HNTight", 15., 2.7); // Do not require lepveto to consider b tagging for the away jet

  bool muonSel = false, electronSel = false;
  //if(muons_veto.size()==1 && electrons_veto.size()==0) muonSel = true;
  //if(muons_veto.size()==0 && electrons_veto_IsoUp.size()==1) electronSel = true;
  if(muons_loose.size()==1 && electrons_loose.size()==0) muonSel = true;
  if(muons_loose.size()==0 && electrons_loose.size()==1) electronSel = true;

  if(!(muonSel || electronSel)) return;

  //std::vector<Jet> jets;
  //jets.clear();
  //jets = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);

  std::vector<Jet> jets_awayFromMuon;
  std::vector<Jet> jets_awayFromElectron;
  jets_awayFromMuon.clear();
  jets_awayFromElectron.clear();

  bool passMuon = false, passElectron = false;

  Particle METv = ev.GetMETVector();
  //double MET = METv.Pt();
  double Mt = 0.;

  if(muonSel){
    if(ev.PassTrigger(triggers_mu)){
      Mt = MT(muons_loose.at(0), METv);
      jets_awayFromMuon = JetsAwayFromLepton(jets, muons_loose.at(0), 1.4);
      if(Mt<40. && jets_awayFromMuon.size()>0) passMuon = true;
    }
  }

  if(electronSel){
    if(ev.PassTrigger(triggers_el)){
      Mt = MT(electrons_loose.at(0), METv);
      jets_awayFromElectron = JetsAwayFromLepton(jets, electrons_loose.at(0), 1.4);
      if(Mt<40. && jets_awayFromElectron.size()>0) passElectron = true;
    }
  }

  if(passMuon || passElectron){
    newtree->Fill();
  }

}

void SkimTree_HNFake::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNFake::SkimTree_HNFake(){
  newtree=NULL;
}

SkimTree_HNFake::~SkimTree_HNFake(){

}

void SkimTree_HNFake::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
