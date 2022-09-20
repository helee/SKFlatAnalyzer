#ifndef SkimTree_HNFake_h
#define SkimTree_HNFake_h

#include "AnalyzerCore.h"

class SkimTree_HNFake : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNFake();
  ~SkimTree_HNFake();

  TTree *newtree;

  double TriggerSafePt_Electron;
  double TriggerSafePt_Muon;

  vector<TString> triggers;
  vector<TString> validation_muon_triggers;
  vector<TString> validation_electron_triggers;
  void WriteHist();

};



#endif
