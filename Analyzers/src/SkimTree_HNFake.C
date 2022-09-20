#include "SkimTree_HNFake.h"

void SkimTree_HNFake::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNFake::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
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
	//triggers = { "HLT_Mu3_PFJet40_v","HLT_Mu50_v"};
        triggers = {"HLT_Mu3_PFJet40_v"};
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
      
	//triggers = {"HLT_Mu3_PFJet40_v","HLT_Mu50_v"};
	triggers = {"HLT_Mu3_PFJet40_v"};

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
    cout<<"[SkimTree_HNFake::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNFake::initializeAnalyzer] triggers to skim = " << endl;
  for(auto i :triggers)  cout << "[SkimTree_HNFake::initializeAnalyzer]   " <<  i << endl;
  for(auto i :validation_electron_triggers)  cout << "[SkimTree_HNFake::initializeAnalyzer]   " <<  i << endl;
  for(auto i :validation_muon_triggers)  cout << "[SkimTree_HNFake::initializeAnalyzer]   " <<  i << endl;

}

void SkimTree_HNFake::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  /*if(ev.PassTrigger(validation_electron_triggers)){
    
    vector<Electron> allel = GetElectrons("HNLoosest", 8., 2.5);
    std::sort(allel.begin(),allel.end(),PtComparing);

    if(allel.size() > 0) {
      if(allel[0].Pt() > TriggerSafePt_Electron){
	newtree->Fill();
	return;
      }
    }

  }

  if(ev.PassTrigger(validation_muon_triggers)){

    vector<Muon> allmuons = GetMuons("HNLoosest", 4., 2.4);
    std::sort(allmuons.begin(),allmuons.end(),PtComparing);

    if(allmuons.size() > 0) {
      if(allmuons[0].Pt() > TriggerSafePt_Muon){
        newtree->Fill();
        return;
      }
    }

  }*/

  //==== Skim 1 ) trigger
  if(!(ev.PassTrigger(triggers))) return;

  if(IsDATA){
    newtree->Fill();
    return;
  }

  //==== Skim 2) only one loose leptons (e or mu)  //==== TODO : To be updated (only for MC) 

  //vector<Muon> allmuons = GetMuons("HNLoosest", 4., 2.4);
  //vector<Electron> allel = GetElectrons("HNLoosest", 8., 2.5);

  vector<Muon> allmuons0p4  = GetMuons("HNLooseIso0p4", 4., 2.4);
  vector<Muon> allmuons0p5  = GetMuons("HNLooseIso0p5", 4., 2.4);
  vector<Electron> allel0p6 = GetElectrons("HNLooseIso0p6", 8., 2.5);
  vector<Electron> allel0p7 = GetElectrons("HNLooseIso0p7", 8., 2.5);

  std::sort(allmuons0p4.begin(),allmuons0p4.end(),PtComparing);
  std::sort(allel0p6.begin(),allel0p6.end(),PtComparing);

  int NLep = allmuons0p4.size() + allel0p6.size();
  int NLepSyst = allmuons0p5.size() + allel0p7.size();
 
  if(!(NLep==1 || NLepSyst==1 || NLep==2)) return;

  if(NLep==1 || NLepSyst==1){

    vector<Jet> alljet = GetJets("tight", 20., 2.7);
  
    bool dphi_lj(false);
    double dphi_cut = 1.5; // 3.0, 2.5 (nominal), 2.0, 1.5

    if(NLep == 1){

      for(unsigned int imu=0; imu<allmuons0p4.size(); imu++){
        for(unsigned int ij=0; ij<alljet.size(); ij++){
          float dphi = fabs(TVector2::Phi_mpi_pi(allmuons0p4[imu].Phi() - alljet.at(ij).Phi()));
          if(dphi > dphi_cut) dphi_lj = true;
        }
      }
  
      for(unsigned int iel=0; iel<allel0p6.size(); iel++){
        for(unsigned int ij=0; ij<alljet.size(); ij++){
          float dphi = fabs(TVector2::Phi_mpi_pi(allel0p6[iel].Phi() - alljet.at(ij).Phi()));
          if(dphi > dphi_cut) dphi_lj = true;
        }
      }

    }
    else{

      for(unsigned int imu=0; imu<allmuons0p5.size(); imu++){
        for(unsigned int ij=0; ij<alljet.size(); ij++){
          float dphi = fabs(TVector2::Phi_mpi_pi(allmuons0p5[imu].Phi() - alljet.at(ij).Phi()));
          if(dphi > dphi_cut) dphi_lj = true;
        }
      }

      for(unsigned int iel=0; iel<allel0p7.size(); iel++){
        for(unsigned int ij=0; ij<alljet.size(); ij++){
          float dphi = fabs(TVector2::Phi_mpi_pi(allel0p7[iel].Phi() - alljet.at(ij).Phi()));
          if(dphi > dphi_cut) dphi_lj = true;
        }
      }

    }

    if(IsDATA){
      if(!dphi_lj) return;
    }
    else if(!MCSample.Contains("QCD")){
      if(!dphi_lj) return;
    }

  }

  if(NLep == 2){

    if(allmuons0p4.size() == 2){
      if(!(allmuons0p4.at(0).Pt()>17. && allmuons0p4.at(1).Pt()>8.)) return;
    }
    else if(allel0p6.size() == 2){
      if(!(allel0p6.at(0).Pt()>23. && allel0p6.at(1).Pt()>12.)) return;
    }
    else return;

  }

  newtree->Fill();
  return;

}

void SkimTree_HNFake::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNFake::SkimTree_HNFake(){
  newtree=NULL;
}

SkimTree_HNFake::~SkimTree_HNFake(){

}

void SkimTree_HNFake::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
