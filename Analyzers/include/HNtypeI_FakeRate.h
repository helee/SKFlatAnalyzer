#ifndef HNtypeI_FakeRate_h
#define HNtypeI_FakeRate_h

#include "HNAnalyzerCore.h"

class HNtypeI_FakeRate : public HNAnalyzerCore {

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

  HNtypeI_FakeRate();
  ~HNtypeI_FakeRate();

};



#endif

