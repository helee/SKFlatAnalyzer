#include "SkimTree_HNFake_UL.h"

void SkimTree_HNFake_UL::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNFake_UL::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);


  triggers.clear();
  validation_muon_triggers.clear();
  validation_electron_triggers.clear();

  if(DataYear==2016)  cout << "2016  " << endl;
  if(DataYear==2017)  cout << "2017  " << endl;
  if(DataYear==2018)  cout << "2018  " << endl;

  //bool IsDATA= (this->DataStream != "");
  if(IsDATA) cout << this->DataStream << endl;
  else    cout <<this->MCSample << endl;

  cout <<  " IsDATA = " << IsDATA << endl;

  if(DataYear==2016){

    if(IsDATA){
      
      if (this->DataStream == "SingleMuon"){
	
	triggers = { "HLT_Mu50_v"};
	validation_muon_triggers = {
          "HLT_IsoMu24_v",
          "HLT_IsoTkMu24_v",
        };
	
      } // 2016 SMu DATA
      if (this->DataStream == "SingleElectron"){    
	validation_electron_triggers = {
	  "HLT_Ele27_WPTight_Gsf_v",
	};
      } // 2016 SE DATA    

      if (this->DataStream == "DoubleMuon"){    
	
	triggers = {
	  "HLT_Mu3_PFJet40_v",                             // DoubleMuon
	  "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
	  "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
	};
      } // // 2016 DMu DATA    

      if (this->DataStream == "DoubleEG"){
        triggers = {
	  "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // DoubleEG
	  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
	  "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
          //"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",        // DoubleEG
	  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // DoubleEG
	};
      } // 2016EGMu DATA    
    } // 2016 DATA

    else{


      validation_muon_triggers = {
	"HLT_IsoMu24_v",
	"HLT_IsoTkMu24_v",
      };

      validation_electron_triggers = {
	"HLT_Ele27_WPTight_Gsf_v",
      };

      triggers = {
	"HLT_Mu50_v",
        "HLT_Mu3_PFJet40_v",                             // DoubleMuon
        "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
        "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // DoubleEG
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // DoubleEG
        //"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",        // DoubleEG
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // DoubleEG
      };

    } // 2016 MC
    
    TriggerSafePt_Electron= 30.;
    TriggerSafePt_Muon = 26.;
  } // 2016
  else if(DataYear==2017){
    

    if(IsDATA){
      
      if (this->DataStream == "SingleMuon"){
	triggers = { "HLT_Mu3_PFJet40_v","HLT_Mu50_v"};
	validation_muon_triggers = {
	  "HLT_IsoMu27_v",
	};
	
      } // 2017 DATA SMu

      if (this->DataStream == "SingleElectron"){
	
	triggers = {
	  "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // SingleElectron
	  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // SingleElectron
	  //"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",        // SingleElectron
	  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // SingleElectron
	};
	validation_electron_triggers = {
	  "HLT_Ele35_WPTight_Gsf_v",
	};
      }// 2017 DATA SEl       

      if (this->DataStream == "DoubleMuon"){

	triggers = {
          "HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon                                                                
          "HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon                                                                
        };

      } // 2017 DATA DMu       
    } // 2017
    
    else {
      
      validation_muon_triggers = {
	"HLT_IsoMu27_v",
      };
      validation_electron_triggers = {
	"HLT_Ele35_WPTight_Gsf_v",
      };
      
      triggers = {
	"HLT_Mu50_v",
	"HLT_Mu3_PFJet40_v",                             // SingleMuon
	"HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
	"HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
	"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // SingleElectron
	"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // SingleElectron
	//"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",        // SingleElectron
	"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // SingleElectron
      };
      cout << "Filling 2017 trigger " << endl;

    } // 2017 MC
    TriggerSafePt_Electron= 38.;
    TriggerSafePt_Muon = 29.;
    
  }
  else if(DataYear==2018){
    if(IsDATA){
   
      if (this->DataStream == "SingleMuon"){
        validation_muon_triggers = {
	  "HLT_IsoMu24_v",
	};
      
	triggers = {"HLT_Mu3_PFJet40_v","HLT_Mu50_v"};
	
      }
      if (this->DataStream == "EGamma"){

        validation_electron_triggers = {
	  "HLT_Ele32_WPTight_Gsf_v",

	};
        triggers = {
          "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v", 
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
          //"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",      
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v" 
        };
      }
      if (this->DataStream == "DoubleMuon"){
	triggers = {
	  "HLT_Mu8_TrkIsoVVL_v", 
	  "HLT_Mu17_TrkIsoVVL_v",
	};
	
      }
    } // 2018 DATA
    else {
      
      validation_muon_triggers = {
	"HLT_IsoMu24_v",
      };
      validation_electron_triggers = {
	"HLT_Ele32_WPTight_Gsf_v",

      };

      triggers = {
        "HLT_Mu50_v",
	"HLT_Mu3_PFJet40_v",                             // SingleMuon
	"HLT_Mu8_TrkIsoVVL_v",                           // DoubleMuon
	"HLT_Mu17_TrkIsoVVL_v",                          // DoubleMuon
	"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",     // EGamma
	"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",    // EGamma
	//"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",        // EGamma
	"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"     // EGamma
      };
      TriggerSafePt_Electron= 35.;
      TriggerSafePt_Muon = 26.;

    }

  }
  else{
    cout<<"[SkimTree_HNFake_UL::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNFake_UL::initializeAnalyzer] triggers to skim = " << endl;
  for(auto i :triggers)  cout << "[SkimTree_HNFake_UL::initializeAnalyzer]   " <<  i << endl;
  for(auto i :validation_electron_triggers)  cout << "[SkimTree_HNFake_UL::initializeAnalyzer]   " <<  i << endl;
  for(auto i :validation_muon_triggers)  cout << "[SkimTree_HNFake_UL::initializeAnalyzer]   " <<  i << endl;

}

void SkimTree_HNFake_UL::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);


  if(ev.PassTrigger(validation_electron_triggers)){
    
    vector<Electron> allel = GetElectrons("HNLoosest", 8., 2.4);
    std::sort(allel.begin(),allel.end(),PtComparing);

    if(allel.size() > 0) {
      if (allel[0].Pt() > TriggerSafePt_Electron) {
	newtree->Fill();
	return;
      }
    }
  }

  if(ev.PassTrigger(validation_muon_triggers)){

    vector<Muon> allmuons = GetMuons("HNLoosest", 5., 2.4);
    std::sort(allmuons.begin(),allmuons.end(),PtComparing);

    if(allmuons.size() > 0) {
      if (allmuons[0].Pt() > TriggerSafePt_Muon) {
        newtree->Fill();
        return;
      }
    }

  }

  //==== Skim 1 ) trigger
  if(! (ev.PassTrigger(triggers)) ) return;


  //==== Skim 2) at least one loose leptons (e or mu) 

  vector<Muon> allmuons = GetMuons("HNLoosest", 5., 2.4);
  vector<Electron> allel = GetElectrons("HNLoosest", 8., 2.4);

  int NLep = allmuons.size() + allel.size();
  
  if( NLep == 0 ) return;

  vector<Jet> alljet = GetJets("HNTight", 30., 2.7);
  
  bool dphi_lj(false);
  for(unsigned int imu=0; imu < allmuons.size(); imu++){
    for(unsigned int ij=0; ij <alljet.size(); ij++){
      float dphi =fabs(TVector2::Phi_mpi_pi(allmuons[imu].Phi()- alljet.at(ij).Phi()));
      if(dphi > 2.5) dphi_lj=true;
    }
  }
  
  for(unsigned int iel=0; iel <allel.size(); iel++){
    for(unsigned int ij=0; ij <alljet.size(); ij++){
      float dphi =fabs(TVector2::Phi_mpi_pi(allel[iel].Phi()- alljet.at(ij).Phi()));
      if(dphi >2.5) dphi_lj=true;
    }
  }

  if(IsDATA){
    if( !dphi_lj ) return;
  }
  else if(!MCSample.Contains("QCD")){
    if( !dphi_lj ) return;
  }
  
  newtree->Fill();
  return;

}

void SkimTree_HNFake_UL::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNFake_UL::SkimTree_HNFake_UL(){
  newtree=NULL;
}

SkimTree_HNFake_UL::~SkimTree_HNFake_UL(){

}

void SkimTree_HNFake_UL::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
