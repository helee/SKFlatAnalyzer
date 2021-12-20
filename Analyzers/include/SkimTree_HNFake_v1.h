#ifndef SkimTree_HNFake_v1_h
#define SkimTree_HNFake_v1_h

#include "AnalyzerCore.h"

class SkimTree_HNFake_v1 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNFake_v1();
  ~SkimTree_HNFake_v1();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif
