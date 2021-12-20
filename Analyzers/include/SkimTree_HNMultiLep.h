#ifndef SkimTree_HNMultiLep_h
#define SkimTree_HNMultiLep_h

#include "HNAnalyzerCore.h"

class SkimTree_HNMultiLep : public HNAnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HNMultiLep();
  ~SkimTree_HNMultiLep();

  TTree *newtree;

  bool isSingleMu;

  vector<TString> triggers;
  vector<TString> triggers_dimu;
  void WriteHist();

  double LeptonPtCut, AK4JetPtCut, AK8JetPtCut;

};



#endif

