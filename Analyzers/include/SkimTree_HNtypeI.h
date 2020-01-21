#ifndef SkimTree_HNtypeI_h
#define SkimTree_HNtypeI_h

#include "AnalyzerCore.h"

class SkimTree_HNtypeI : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNtypeI();
  ~SkimTree_HNtypeI();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif

