#include "HNL_JetIDOpt.h"

void HNL_JetIDOpt::initializeAnalyzer(){

  // All default settings like trigger/ PD/ BJet are decalred in HNL_LeptonCore::initializeAnalyzer to make them consistent for all HNL codes

  HNL_LeptonCore::initializeAnalyzer();
  cout << "SetupMVAReader " << endl;
  SetupMVAReader();
}


void HNL_JetIDOpt::executeEvent(){
  
  if((_jentry==0)){ 
    // Print out trigger info in HNL_LeptonCore::initializeAnalyzer
    TriggerPrintOut(GetEvent());
  }
  
  if(!IsData)  gens = GetGens();

  //AnalyzerParameter param_signal = HNL_LeptonCore::InitialiseHNLParameter("MVAUL","_UL");
  AnalyzerParameter param_signal = HNL_LeptonCore::InitialiseHNLParameter("HNL","_UL");
  RunULAnalysis(param_signal);

  if(!IsData) RunSyst=true;
  if(RunSyst){
    TString param_signal_name = param_signal.Name;
    vector<AnalyzerParameter::Syst> SystList;// = GetSystList("Initial");

    for(auto isyst : SystList){
      param_signal.syst_ = AnalyzerParameter::Syst(isyst);
      
      param_signal.Name = "Syst_"+param_signal.GetSystType()+param_signal_name;
      param_signal.DefName = "Syst_"+param_signal.GetSystType()+param_signal_name;
      RunULAnalysis(param_signal);
    }
  }    


  return ;
}

void HNL_JetIDOpt::RunULAnalysis(AnalyzerParameter param){

  if(run_Debug) cout << "HNL_JetIDOpt::executeEvent " << endl;

  Event ev = GetEvent();
  double weight =SetupWeight(ev,param);


  ///// MERGE WJet samples for more stats                                                                                                                                      
  //if(MCSample.Contains("WJet")){
  //  vector<TString> vec = {"WJet"};
  //  double merge_weight = MergeMultiMC( vec, "" );
  //  weight*= merge_weight;
  // }

  /// Merge DY samples for more stats                                                                                                                                        
  //if(MCSample.Contains("DYJets_MG")){
  //  vector<TString> vec = {"DYMG"};
  //  double merge_weight = MergeMultiMC( vec, "" );
  //  weight*= merge_weight;
  // }


  
  // HL ID
  std::vector<Electron>   ElectronCollV = GetElectrons(param.Electron_Veto_ID, 10., 2.5); 
  std::vector<Muon>       MuonCollV     = GetMuons    (param.Muon_Veto_ID, 5., 2.4);

  TString el_ID = (RunFake) ?  param.Electron_FR_ID : param.Electron_Tight_ID ;
  TString mu_ID = (RunFake) ?  param.Muon_FR_ID :  param.Muon_Tight_ID ;

  double Min_Muon_Pt     = (RunFake) ? 3. : 5.;
  double Min_Electron_Pt = (RunFake) ? 7. : 10.;

  std::vector<Muon>       MuonCollTInit = GetMuons    ( param,mu_ID, Min_Muon_Pt, 2.4, false);
  std::vector<Electron>   ElectronCollTInit = GetElectrons( param,el_ID, Min_Electron_Pt, 2.5, false)  ;

  std::vector<Muon>       MuonCollT     = GetLepCollByRunType    ( MuonCollTInit,gens,param);
  std::vector<Electron>   ElectronCollT  =  GetLepCollByRunType   ( ElectronCollTInit,gens,param);

  std::vector<Lepton *> leps_veto  = MakeLeptonPointerVector(MuonCollV,ElectronCollV);

  std::vector<Tau>        TauColl        = GetTaus     (leps_veto,param.Tau_Veto_ID,20., 2.3);

  //==== Creat Lepton vector to have lepton blind codes 

  //==== MET
  Particle METv = GetvMET("PuppiT1xyCorr",param); // returns MET with systematic correction

  //==== Jet map
  map<TString, std::vector<Jet> > jet_map;
  map<TString, std::vector<Jet> > vbfjet_map;
  map<TString, std::vector<FatJet> > fatjet_map;
  map<TString, std::vector<Jet> > bjet_map;
  
  std::vector<FatJet> fatjets_tmp    = GetFatJets(param, param.FatJet_ID, 200., 5.);
  std::vector<Jet> jets_tmp          = GetJets(param, param.Jet_ID, 15., 5.);
  std::vector<Jet> All_JetColl       = GetJets("NoID", 10., 3.0);
  
  //std::vector<FatJet> AK8_JetColl                 = SelectAK8Jets(fatjets_tmp, 200., 5., true,  1., false, -999, false, 0., 20000., ElectronCollV, MuonCollV);
  //std::vector<FatJet> AK8_JetColl                 = SelectAK8Jets(fatjets_tmp, 200., 2.7, true,  1., false, -999, false, 40., 130., ElectronCollV, MuonCollV);

  if(HasFlag("JetOptAK8")){

    vector<double> pT_AK8 = {200., 225., 250.};
    vector<double> eta_AK8 = {2.4, 2.7};

    vector<TString> s_pT_AK8 = {"200", "225", "250"};
    vector<TString> s_eta_AK8 = {"2p4", "2p7"};

    //TString idtag = "";

    for(unsigned int ipt=0; ipt<pT_AK8.size(); ipt++){
      for(unsigned int ieta=0; ieta<eta_AK8.size(); ieta++){

        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_NoMassNoTau21"] = SelectAK8Jets(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., false, -999., -999., ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_NoMassPNLoose"] = SelectAK8Jetsv2(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., false, -999., -999., "Loose", ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_NoTau21"]       = SelectAK8Jets(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., true,  40., 130., ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_Tau21HP"]       = SelectAK8Jets(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., true, -1., true,  40., 130., ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_Tau21LP"]       = SelectAK8Jets(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., true,  1., true,  40., 130., ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_PNLoose"]       = SelectAK8Jetsv2(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., true, 40., 130., "Loose", ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_PNMedium"]      = SelectAK8Jetsv2(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., true, 40., 130., "Medium", ElectronCollV, MuonCollV);
        fatjet_map["AK8_Pt"+s_pT_AK8.at(ipt)+"_Eta"+s_eta_AK8.at(ieta)+"_PNTight"]       = SelectAK8Jetsv2(fatjets_tmp, pT_AK8.at(ipt), eta_AK8.at(ieta), true,  1., false, 0., true, 40., 130., "Tight", ElectronCollV, MuonCollV);

      }
    }

    for(std::map<TString, std::vector<FatJet> >::iterator ak8_it = fatjet_map.begin(); ak8_it != fatjet_map.end(); ak8_it++){

      double weight_optjet = weight;

      param.Name = param.DefName + "_" + ak8_it->first;

      //==== AK4 Jet selection
 
      TString PUIDWP = "";
 
      std::vector<Jet> JetCollLoose  = SelectAK4Jets(jets_tmp, 15., 4.7, true, 0.4, 0.8, "",  ElectronCollV, MuonCollV, ak8_it->second);

      std::vector<Jet> JetColl       = SelectAK4Jets(jets_tmp, 20., 2.7, true, 0.4, 0.8, PUIDWP, ElectronCollV, MuonCollV, ak8_it->second);
      std::vector<Jet> VBF_JetColl   = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, PUIDWP, ElectronCollV, MuonCollV, ak8_it->second);    // High Eta jets                 
      std::vector<Jet> BJetColltmp   = SelectAK4Jets(jets_tmp, 20., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, ak8_it->second);
  
      //double PJet_PUID_weight = GetJetPileupIDSF(JetColl, PUIDWP, param);
      //    weight*= PJet_PUID_weight;
      //FillWeightHist("PJet_PUID_weight_" ,PJet_PUID_weight);
  
  
      // select B jets
      JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  
      // Get BJets  and EV weight to corr BTag Eff
      std::vector<Jet> BJetColl    = SelectBJets(param, BJetColltmp, param_jets);
      double sf_btag               = GetBJetSF(param, BJetColltmp, param_jets);

      //if(!IsData )weight*= sf_btag;

      //JetTagging::Parameters param_jetsT = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::mujets);
      //std::vector<Jet> BJetCollSR1    = SelectBJets(param,  BJetColltmp, param_jets);
      //double sf_btagSR1               = GetBJetSF(param, BJetColltmp, param_jets);

      //if(!IsData && AK8_JetColl.size()==0)weight = weight*sf_btag;
      //if(!IsData && AK8_JetColl.size()>0)weight = weight*sf_btagSR1;

      //weight = weight*sf_btag;

      weight_optjet *= sf_btag;

      if(ak8_it->first.Contains("Tau21HP"))  weight_optjet *= GetEventFatJetSF(ak8_it->second, "HP", 0);
      if(ak8_it->first.Contains("Tau21LP"))  weight_optjet *= GetEventFatJetSF(ak8_it->second, "LP", 0);
      if(ak8_it->first.Contains("PNLoose"))  weight_optjet *= GetEventFatJetSFPN(ak8_it->second, "Loose", 0);
      if(ak8_it->first.Contains("PNMedium")) weight_optjet *= GetEventFatJetSFPN(ak8_it->second, "Medium", 0);
      if(ak8_it->first.Contains("PNTight"))  weight_optjet *= GetEventFatJetSFPN(ak8_it->second, "Tight", 0);

      RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose, All_JetColl, JetColl, VBF_JetColl, ak8_it->second, BJetColl, BJetColl, ev, METv, param, weight_optjet);

      param.Name = param.DefName;

    }

  }

  if(HasFlag("JetOptAK4")){

    vector<double> pT_AK4 = {20., 25., 30.};
    //vector<double> eta_AK4 = {2.4, 2.7};
    vector<double> eta_AK4 = {2.7};
    vector<double> pT_VBF = {15., 20., 25., 30.};

    vector<TString> s_pT_AK4 = {"20", "25", "30"};
    //vector<TString> s_eta_AK4 = {"2p4", "2p7"};
    vector<TString> s_eta_AK4 = {"2p7"};
    vector<TString> s_pT_VBF = {"15", "20", "25", "30"};
    //vector<TString> s_PUWP = {"NoCut", "Loose", "Medium", "Tight"};

    std::vector<FatJet> FatJetColl_1 = SelectAK8Jetsv2(fatjets_tmp, 200., 2.7, true,  1., false, 0., true, 40., 130., "Loose", ElectronCollV, MuonCollV);
    std::vector<Jet> BJetColltmp_1   = SelectAK4Jets(jets_tmp, 20., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_1);

    JetTagging::Parameters param_jets_1 = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
    double sf_fatjet_1           = GetEventFatJetSFPN(FatJetColl_1, "Loose", 0);
    std::vector<Jet> BJetColl_1  = SelectBJets(param, BJetColltmp_1, param_jets_1);
    double sf_btag_1             = GetBJetSF(param, BJetColltmp_1, param_jets_1);

    for(unsigned int ipt=0; ipt<pT_AK4.size(); ipt++){
      for(unsigned int ieta=0; ieta<eta_AK4.size(); ieta++){

        jet_map["AK4_Pt"+s_pT_AK4.at(ipt)+"_Eta"+s_eta_AK4.at(ieta)+"_PUNoCut"] = SelectAK4Jets(jets_tmp, pT_AK4.at(ipt), eta_AK4.at(ieta), true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_1);
        jet_map["AK4_Pt"+s_pT_AK4.at(ipt)+"_Eta"+s_eta_AK4.at(ieta)+"_PULoose"] = SelectAK4Jets(jets_tmp, pT_AK4.at(ipt), eta_AK4.at(ieta), true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_1);
        jet_map["AK4_Pt"+s_pT_AK4.at(ipt)+"_Eta"+s_eta_AK4.at(ieta)+"_PUMedium"] = SelectAK4Jets(jets_tmp, pT_AK4.at(ipt), eta_AK4.at(ieta), true, 0.4, 0.8, "Medium", ElectronCollV, MuonCollV, FatJetColl_1);
        jet_map["AK4_Pt"+s_pT_AK4.at(ipt)+"_Eta"+s_eta_AK4.at(ieta)+"_PUTight"] = SelectAK4Jets(jets_tmp, pT_AK4.at(ipt), eta_AK4.at(ieta), true, 0.4, 0.8, "Tight", ElectronCollV, MuonCollV, FatJetColl_1);

      }
    }

    for(unsigned int ivbf=0; ivbf<pT_VBF.size(); ivbf++){

      vbfjet_map["HTJetPt"+s_pT_VBF.at(ivbf)] = SelectAK4Jets(jets_tmp, pT_VBF.at(ivbf), 4.7, true, 0.4, 0.8, "",  ElectronCollV, MuonCollV, FatJetColl_1);

    }

    for(std::map<TString, std::vector<Jet> >::iterator ak4_it = jet_map.begin(); ak4_it != jet_map.end(); ak4_it++){

      double weight_optjet = weight;
      double weight_SR2 = 1., weight_SR3 = 1.;

      std::vector<Jet> VBF_JetColl_1, Forward_JetColl;
      VBF_JetColl_1.clear();
      Forward_JetColl.clear();
      if(ak4_it->first.Contains("PUNoCut")){
        VBF_JetColl_1 = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p4")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.4, 4.7, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p7")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.7, 4.7, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_1);
      }
      else if(ak4_it->first.Contains("PULoose")){
        VBF_JetColl_1 = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p4")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.4, 4.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p7")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.7, 4.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_1);
      }
      else if(ak4_it->first.Contains("PUMedium")){
        VBF_JetColl_1 = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, "Medium", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p4")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.4, 4.7, true, 0.4, 0.8, "Medium", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p7")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.7, 4.7, true, 0.4, 0.8, "Medium", ElectronCollV, MuonCollV, FatJetColl_1);
      }
      else if(ak4_it->first.Contains("PUTight")){
        VBF_JetColl_1 = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, "Tight", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p4")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.4, 4.7, true, 0.4, 0.8, "Tight", ElectronCollV, MuonCollV, FatJetColl_1);
        if(ak4_it->first.Contains("2p7")) Forward_JetColl = SelectAK4Jets(jets_tmp, 30., 2.7, 4.7, true, 0.4, 0.8, "Tight", ElectronCollV, MuonCollV, FatJetColl_1);
      }

      //std::vector<Jet> BJetColl_1  = SelectBJets(param, BJetColltmp_1, param_jets_1);
      //double sf_btag_1             = GetBJetSF(param, BJetColltmp_1, param_jets_1);
      weight_optjet *= sf_fatjet_1*sf_btag_1;

      if(ak4_it->first.Contains("PULoose")){
        weight_SR3 = GetJetPileupIDSF(ak4_it->second, "Loose", param);
        weight_SR2 = weight_SR3*GetJetPileupIDSF(Forward_JetColl, "Loose", param);
      }
      if(ak4_it->first.Contains("PUMedium")){
        weight_SR3 = GetJetPileupIDSF(ak4_it->second, "Medium", param);
        weight_SR2 = weight_SR3*GetJetPileupIDSF(Forward_JetColl, "Medium", param);
      }
      if(ak4_it->first.Contains("PUTight")){
        weight_SR3 = GetJetPileupIDSF(ak4_it->second, "Tight", param);
        weight_SR2 = weight_SR3*GetJetPileupIDSF(Forward_JetColl, "Tight", param);
      }

      for(std::map<TString, std::vector<Jet> >::iterator vbf_it = vbfjet_map.begin(); vbf_it != vbfjet_map.end(); vbf_it++){

        param.Name = param.DefName + "_" + ak4_it->first + "_" + vbf_it->first;

        RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, vbf_it->second, All_JetColl, ak4_it->second, VBF_JetColl_1, FatJetColl_1, BJetColl_1, BJetColl_1, ev, METv, param, weight_optjet, weight_SR2, weight_SR3, true);

        param.Name = param.DefName;

      }
    }

  }

  if(HasFlag("JetOptBTag")){

    std::vector<FatJet> FatJetColl_2   = SelectAK8Jetsv2(fatjets_tmp, 200., 2.7, true,  1., false, 0., true, 40., 130., "Loose", ElectronCollV, MuonCollV);;
    std::vector<Jet> BJetColltmp_2     = SelectAK4Jets(jets_tmp, 15., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_2);
    std::vector<Jet> BJetColltmp_3     = SelectAK4Jets(jets_tmp, 20., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_2);
    std::vector<Jet> BJetColltmp_4     = SelectAK4Jets(jets_tmp, 25., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_2);
    std::vector<Jet> BJetColltmp_5     = SelectAK4Jets(jets_tmp, 30., 2.4, true, 0.4, 0.8, "", ElectronCollV, MuonCollV, FatJetColl_2);

    std::vector<Jet> JetCollLoose_2    = SelectAK4Jets(jets_tmp, 15., 4.7, true, 0.4, 0.8, "",  ElectronCollV, MuonCollV, FatJetColl_2);

    std::vector<Jet> JetColl_2         = SelectAK4Jets(jets_tmp, 20., 2.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_2);
    std::vector<Jet> VBF_JetColl_2     = SelectAK4Jets(jets_tmp, 30., 4.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_2);
    std::vector<Jet> Forward_JetColl_2 = SelectAK4Jets(jets_tmp, 30., 2.7, 4.7, true, 0.4, 0.8, "Loose", ElectronCollV, MuonCollV, FatJetColl_2);

    double sf_fatjet_2 = GetEventFatJetSFPN(FatJetColl_2, "Loose", 0);
    double weightSR2 = 1., weightSR3 = 1.;
    weightSR3 = GetJetPileupIDSF(JetColl_2, "Loose", param);
    weightSR2 = weightSR3*GetJetPileupIDSF(Forward_JetColl_2, "Loose", param);

    //==== B tagging WP

    JetTagging::Parameters param_jets_2 = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::mujets);
    JetTagging::Parameters param_jets_3 = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
    JetTagging::Parameters param_jets_4 = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Tight, JetTagging::incl, JetTagging::mujets);
    
    std::vector<Jet> BJetColl_22   = SelectBJets(param, BJetColltmp_2, param_jets_2);
    double sf_btag_22              = GetBJetSF(param, BJetColltmp_2, param_jets_2);

    std::vector<Jet> BJetColl_23   = SelectBJets(param, BJetColltmp_2, param_jets_3);
    double sf_btag_23              = GetBJetSF(param, BJetColltmp_2, param_jets_3);

    std::vector<Jet> BJetColl_24   = SelectBJets(param, BJetColltmp_2, param_jets_4);
    double sf_btag_24              = GetBJetSF(param, BJetColltmp_2, param_jets_4);

    std::vector<Jet> BJetColl_32   = SelectBJets(param, BJetColltmp_3, param_jets_2);
    double sf_btag_32              = GetBJetSF(param, BJetColltmp_3, param_jets_2);

    std::vector<Jet> BJetColl_33   = SelectBJets(param, BJetColltmp_3, param_jets_3);
    double sf_btag_33              = GetBJetSF(param, BJetColltmp_3, param_jets_3);

    std::vector<Jet> BJetColl_34   = SelectBJets(param, BJetColltmp_3, param_jets_4);
    double sf_btag_34              = GetBJetSF(param, BJetColltmp_3, param_jets_4);

    std::vector<Jet> BJetColl_42   = SelectBJets(param, BJetColltmp_4, param_jets_2);
    double sf_btag_42              = GetBJetSF(param, BJetColltmp_4, param_jets_2);

    std::vector<Jet> BJetColl_43   = SelectBJets(param, BJetColltmp_4, param_jets_3);
    double sf_btag_43              = GetBJetSF(param, BJetColltmp_4, param_jets_3);

    std::vector<Jet> BJetColl_44   = SelectBJets(param, BJetColltmp_4, param_jets_4);
    double sf_btag_44              = GetBJetSF(param, BJetColltmp_4, param_jets_4);

    std::vector<Jet> BJetColl_52   = SelectBJets(param, BJetColltmp_5, param_jets_2);
    double sf_btag_52              = GetBJetSF(param, BJetColltmp_5, param_jets_2);

    std::vector<Jet> BJetColl_53   = SelectBJets(param, BJetColltmp_5, param_jets_3);
    double sf_btag_53              = GetBJetSF(param, BJetColltmp_5, param_jets_3);

    std::vector<Jet> BJetColl_54   = SelectBJets(param, BJetColltmp_5, param_jets_4);
    double sf_btag_54              = GetBJetSF(param, BJetColltmp_5, param_jets_4);

    //==== RunAllSR

    param.Name = param.DefName + "_BTagPt15Loose";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_22, BJetColl_22, ev, METv, param, weight*sf_fatjet_2*sf_btag_22, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt15Medium";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_23, BJetColl_23, ev, METv, param, weight*sf_fatjet_2*sf_btag_23, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt15Tight";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_24, BJetColl_24, ev, METv, param, weight*sf_fatjet_2*sf_btag_24, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt20Loose";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_32, BJetColl_32, ev, METv, param, weight*sf_fatjet_2*sf_btag_32, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt20Medium";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_33, BJetColl_33, ev, METv, param, weight*sf_fatjet_2*sf_btag_33, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt20Tight";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_34, BJetColl_34, ev, METv, param, weight*sf_fatjet_2*sf_btag_34, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt25Loose";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_42, BJetColl_42, ev, METv, param, weight*sf_fatjet_2*sf_btag_42, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt25Medium";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_43, BJetColl_43, ev, METv, param, weight*sf_fatjet_2*sf_btag_43, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt25Tight";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_44, BJetColl_44, ev, METv, param, weight*sf_fatjet_2*sf_btag_44, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt30Loose";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_52, BJetColl_52, ev, METv, param, weight*sf_fatjet_2*sf_btag_52, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt30Medium";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_53, BJetColl_53, ev, METv, param, weight*sf_fatjet_2*sf_btag_53, weightSR2, weightSR3, true);

    param.Name = param.DefName + "_BTagPt30Tight";
    RunAllSignalRegions(Inclusive, ElectronCollT, ElectronCollV, MuonCollT, MuonCollV, TauColl, JetCollLoose_2, All_JetColl, JetColl_2, VBF_JetColl_2, FatJetColl_2, BJetColl_54, BJetColl_54, ev, METv, param, weight*sf_fatjet_2*sf_btag_54, weightSR2, weightSR3, true);

    param.Name = param.DefName;

  }

}



 


HNL_JetIDOpt::HNL_JetIDOpt(){

  cout << "HNL_JetIDOpt::HNL_JetIDOpt  TMVA::Tools::Instance() " << endl;
  TMVA::Tools::Instance();
  cout << "Create Reader class " << endl;
  //MVAReader = new TMVA::Reader();
  //MVAReaderMM = new TMVA::Reader();
  MVAReaderEE = new TMVA::Reader();
  MVAReaderEM = new TMVA::Reader();
  
}
 
HNL_JetIDOpt::~HNL_JetIDOpt(){

  //delete MVAReader;
  //delete MVAReaderMM;
  delete MVAReaderEE;
  delete MVAReaderEM;

}




