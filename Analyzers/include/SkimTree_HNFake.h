#ifndef SkimTree_HNFake_h
#define SkimTree_HNFake_h

#include "HNAnalyzerCore.h"

class SkimTree_HNFake : public HNAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNFake();
  ~SkimTree_HNFake();

  TTree *newtree;

  vector<TString> triggers, triggers_mu, triggers_el;
  void WriteHist();

};



#endif
