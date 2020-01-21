#include "SkimTree_HNtypeI.h"

void SkimTree_HNtypeI::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNtypeI::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     // B-G
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    // B-G
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  // H
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"  // H
    };
  }else if(DataYear==2017){
    triggers = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }else if(DataYear==2018){
    triggers = {
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
  }else{
    cout<<"[SkimTree_HNtypeI::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNtypeI::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HNtypeI::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

}

void SkimTree_HNtypeI::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(triggers) ){
    newtree->Fill();
  }

}

void SkimTree_HNtypeI::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNtypeI::SkimTree_HNtypeI(){
  newtree=NULL;
}

SkimTree_HNtypeI::~SkimTree_HNtypeI(){

}

void SkimTree_HNtypeI::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


