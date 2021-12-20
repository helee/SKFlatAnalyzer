#ifndef SkimTree_Dilepton_v1_h
#define SkimTree_Dilepton_v1_h

#include "AnalyzerCore.h"

class SkimTree_Dilepton_v1 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_Dilepton_v1();
  ~SkimTree_Dilepton_v1();

  TTree *newtree;

  vector<TString> double_triggers;
  vector<TString> single_muon_triggers;
  vector<TString> single_electron_triggers;
  void WriteHist();

};



#endif
