#include "Tau.h"

ClassImp(Tau)

Tau::Tau(){

  j_IDBit = 0;
  j_decaymode=-1;
  j_idDecayModeNewDMs=false;

  this->SetLeptonFlavour(TAU);
}

Tau::~Tau(){

}

void Tau::SetIDBit(unsigned int idbit){
  j_IDBit = idbit;
}


void Tau::SetDecayMode(int decaymode){

  j_decaymode= decaymode;
}


void Tau::SetDecayModeNewDM(bool DecayModeNewDMs){

  j_idDecayModeNewDMs= DecayModeNewDMs;
}

   
bool Tau::PassID(TString ID) const{

  //==== list of IDs for analyis
  if(ID=="NoCut") return true;

  //==== HNVeto IDs
  int id = 1;
  vector<bool> IDvJet = {passVVLIDvJet(), passVLIDvJet(), passLIDvJet(), passMIDvJet(), passTIDvJet(), passVTIDvJet(), passVVTIDvJet()};
  vector<bool> IDvEl  = {passVVLIDvEl(), passVLIDvEl(), passLIDvEl(), passMIDvEl(), passTIDvEl(), passVTIDvEl(), passVVTIDvEl()};
  vector<bool> IDvMu  = {passVLIDvMu(), passLIDvMu(), passMIDvMu(), passTIDvMu()};

  for(unsigned int i=0; i<IDvJet.size(); i++){

    for(unsigned int j=0; j<IDvEl.size(); j++){

      for(unsigned int k=0; k<IDvMu.size(); k++){

        if(ID=="HNVetoV"+TString::Itoa(id, 10)){

          if(j_decaymode==0 || j_decaymode==1 || j_decaymode==10 || j_decaymode==11){
      
            if(DecayModeNewDM() && IDvJet.at(i) && IDvEl.at(j) && IDvMu.at(k)) return true;

          }

        }

        id++;

      }

    }

  }

  return false;

}
