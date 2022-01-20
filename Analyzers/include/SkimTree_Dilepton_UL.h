#ifndef SkimTree_Dilepton_UL_h
#define SkimTree_Dilepton_UL_h

#include "AnalyzerCore.h"

class SkimTree_Dilepton_UL : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_Dilepton_UL();
  ~SkimTree_Dilepton_UL();

  TTree *newtree;

  vector<TString> double_triggers;
  vector<TString> single_muon_triggers;
  vector<TString> single_electron_triggers;
  void WriteHist();

};



#endif
