#ifndef HNtypeI_DY_CR_h
#define HNtypeI_DY_CR_h

#include "HNAnalyzerCore.h"

class HNtypeI_DY_CR : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  //bool RunNewPDF;
  //bool RunXSecSyst;
  bool RunFake;

  //==== Trigger
  vector<TString> MuonTriggers;
  vector<TString> MuonTriggersTight;
  vector<TString> ElectronTriggers;
  vector<TString> EMuTriggers;
  vector<TString> EMuTriggersTight;
  vector<TString> EMuTriggersMu8;
  vector<TString> EMuTriggersMu23;

  //==== Lepton ID
  vector<TString> MuonVetoIDs;
  vector<TString> MuonLooseIDs;
  vector<TString> MuonTightIDs;
  vector<TString> ElectronVetoIDs;
  vector<TString> ElectronLooseIDs;
  vector<TString> ElectronTightIDs;

  //==== Fake rate
  vector<TString> MuonFRNames;
  vector<TString> ElectronFRNames;

  //==== Lepton pT cut
  double MuonPtCut1;
  double MuonPtCut2;
  double ElectronPtCut1;
  double ElectronPtCut2;
  double EMuPtCut1;
  double EMuPtCut2;

  //vector<TString> MuonIDs, MuonIDSFKeys;
  //double weight_Prefire;

  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<FatJet> AllFatJets;

  HNtypeI_DY_CR();
  ~HNtypeI_DY_CR();

};



#endif

