#ifndef HNtypeI_SR_h
#define HNtypeI_SR_h

#include "HNAnalyzerCore.h"

class HNtypeI_SR : public HNAnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  double weight_Prefire;

  HNtypeI_SR();
  ~HNtypeI_SR();

};



#endif

