#ifndef HNL_ControlRegionPlotter_h
#define HNL_ControlRegionPlotter_h

#include "HNL_RegionDefinitions.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class HNL_ControlRegionPlotter : public HNL_RegionDefinitions {

 public:


  void initializeAnalyzer();
  void executeEvent();

  HNL_ControlRegionPlotter();
  ~HNL_ControlRegionPlotter();

  void RunControlRegions(AnalyzerParameter param);


};



#endif
