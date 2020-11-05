#ifndef HNtypeI_Cutflow_h
#define HNtypeI_Cutflow_h

#include "HNAnalyzerCore.h"

class HNtypeI_Cutflow : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;
  bool RunMuon, RunTightIP, RunElectron, RunEMu, RunTightEMu, RunFake, RunCF;

  // Trigger
  vector<TString> MuonTriggers;
  vector<TString> MuonTriggersH;
  vector<TString> ElectronTriggers;
  vector<TString> EMuTriggers;
  vector<TString> EMuTriggersH;
  vector<TString> Mu8Ele23Triggers;
  vector<TString> Mu23Ele12Triggers;

  // Lepton ID
  vector<TString> MuonVetoIDs;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronVetoIDs;
  vector<TString> ElectronLooseIDs;
  vector<TString> ElectronTightIDs;

  // Fake rate file
  vector<TString> MuonFRNames;
  vector<TString> ElectronFRNames;

  // Lepton pT cut
  double MuonPtCut1;
  double MuonPtCut2;
  double ElectronPtCut1;
  double ElectronPtCut2;
  double EMuPtCut1;
  double EMuPtCut2;

  //vector<TString> EleIDs, EleIDSFKeys, MuonIDs, MuonIDSFKeys;
  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<FatJet> AllFatJets;

  //double weight_Prefire;

  HNtypeI_Cutflow();
  ~HNtypeI_Cutflow();

};

#endif
