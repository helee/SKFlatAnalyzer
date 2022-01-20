#ifndef SkimTree_HNFake_UL_h
#define SkimTree_HNFake_UL_h

#include "HNAnalyzerCore.h"

class SkimTree_HNFake_UL : public HNAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNFake_UL();
  ~SkimTree_HNFake_UL();

  TTree *newtree;

  double TriggerSafePt_Electron;
  double TriggerSafePt_Muon;

  vector<TString> triggers;
  vector<TString> validation_muon_triggers;
  vector<TString> validation_electron_triggers;
  void WriteHist();

};



#endif
