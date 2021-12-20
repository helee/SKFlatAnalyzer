#ifndef HNtypeI_VV_CR_h
#define HNtypeI_VV_CR_h

#include "HNAnalyzerCore.h"

class HNtypeI_VV_CR : public HNAnalyzerCore {

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
  vector<TString> MuonTriggersNoDZ;
  vector<TString> ElectronTriggers;
  vector<TString> EMuTriggers;
  vector<TString> EMuTriggersNoDZ;

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

  //TString IsoMuTriggerName;
  //double TriggerSafePtCut;
  //vector<TString> MuonIDs, MuonIDSFKeys;

  vector<Electron> AllElectrons;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<FatJet> AllFatJets;

  //double weight_Prefire;

  HNtypeI_VV_CR();
  ~HNtypeI_VV_CR();

};



#endif

