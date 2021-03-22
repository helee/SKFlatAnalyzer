#include "HNtypeI_QCDFakeRateISR.h"

HNtypeI_QCDFakeRateISR::HNtypeI_QCDFakeRateISR(){

}

void HNtypeI_QCDFakeRateISR::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");

  cout << "[HNtypeI_QCDFakeRateISR::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  MuonTightIDs     = {"ISRTight"};
  MuonLooseIDs     = {"ISRLoose"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"ISRTight"};
  ElectronLooseIDs = {"ISRLoose"};
  ElectronVetoIDs  = {"ISRVeto"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_QCDFakeRateISR.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
  MuonTriggers.clear();
  ElectronTriggers.clear();

  MuonTrig1 = "HLT_Mu3_PFJet40_v";       // DoubleMuon(2016), SingleMuon(2017,2018)
  MuonTrig2 = "HLT_Mu8_TrkIsoVVL_v";     // DoubleMuon
  MuonTrig3 = "HLT_Mu17_TrkIsoVVL_v";    // DoubleMuon

  MuonTriggers.push_back(MuonTrig1);      
  MuonTriggers.push_back(MuonTrig2);    
  MuonTriggers.push_back(MuonTrig3);   
  MuonPtCut1 = 5., MuonPtCut2 = 10., MuonPtCut3 = 20.;
  MuonPtconeCut1 = 5., MuonPtconeCut2 = 15., MuonPtconeCut3 = 30.;

  // DoubleEG (2016), SingleElectron (2017), EGamma (2018)
  if(DataYear==2016){
    ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig3 = "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

    ElectronTriggers.push_back(ElectronTrig1);
    ElectronTriggers.push_back(ElectronTrig2);
    ElectronTriggers.push_back(ElectronTrig3);
    ElectronTriggers.push_back(ElectronTrig4);
    ElectronPtCut1 = 10., ElectronPtCut2 = 15., ElectronPtCut3 = 20., ElectronPtCut4 = 25.;
    ElectronPtconeCut1 = 15., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 45.;
  }
  else if(DataYear==2017){
    ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig3 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

    ElectronTriggers.push_back(ElectronTrig1);
    ElectronTriggers.push_back(ElectronTrig2);
    ElectronTriggers.push_back(ElectronTrig3);
    ElectronTriggers.push_back(ElectronTrig4);
    ElectronPtCut1 = 10., ElectronPtCut2 = 15., ElectronPtCut3 = 15., ElectronPtCut4 = 25.;
    ElectronPtconeCut1 = 15., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 45.;
  }
  else if(DataYear==2018){
    ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig3 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

    ElectronTriggers.push_back(ElectronTrig1);
    ElectronTriggers.push_back(ElectronTrig2);
    ElectronTriggers.push_back(ElectronTrig3);
    ElectronTriggers.push_back(ElectronTrig4);
    ElectronPtCut1 = 10., ElectronPtCut2 = 15., ElectronPtCut3 = 15., ElectronPtCut4 = 25.;
    ElectronPtconeCut1 = 15., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 45.;
  }

  //cout << "[HNtypeI_QCDFakeRateISR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_QCDFakeRateISR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== b tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_QCDFakeRateISR::~HNtypeI_QCDFakeRateISR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_QCDFakeRateISR::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_QCDFakeRateISR.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_QCDFakeRateISR.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){

    //TString MuonID = "HNTight2016";
    //TString MuonIDSFKey = "NUM_TightID_DEN_genTracks";
    TString MuonTightID = MuonTightIDs.at(it_id);
    TString MuonLooseID = MuonLooseIDs.at(it_id);
    TString MuonVetoID  = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);

    param.Clear();

    param.fakesyst_ = AnalyzerParameter::FakeCentral;

    //param.Name = MuonID+"_"+"Central";
    param.Name = "FakeCentral";

    // Muon ID
    param.Muon_Tight_ID = MuonTightID;
    param.Muon_Loose_ID = MuonLooseID;
    param.Muon_Veto_ID  = MuonVetoID;
    param.Muon_ID_SF_Key = "";
    param.Muon_ISO_SF_Key = "";

    // Electron ID
    param.Electron_Tight_ID = ElectronTightID;
    param.Electron_Loose_ID = ElectronLooseID;
    param.Electron_Veto_ID  = ElectronVetoID;
    param.Electron_ID_SF_Key = "";

    // Jet ID
    param.Jet_ID = "HNTight";

    executeEventFromParameter(param);

    /*if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NFakeSyst; it_syst++){
        param.fakesyst_ = AnalyzerParameter::FakeSyst(it_syst);
        //param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        param.Name  = "FakeSyst_"+param.GetFakeSystType();
        executeEventFromParameter(param);
      }
    }*/

  }

}

void HNtypeI_QCDFakeRateISR::executeEventFromParameter(AnalyzerParameter param){

  TString MuonIDname = "MuonISRRun2";

  TString ElectronIDname = "ElectronISRRun2";

  vector<TString> regions = {"FR", "DY", "Wjet"};
  TString btagDirName = "";

  // Boolean : QCD MC
  bool isMuon = false, isElectron = false;
  if(MCSample.Contains("MuEnriched")) isMuon = true;
  if(MCSample.Contains("EMEnriched") || MCSample.Contains("bcToE")) isElectron = true;

  TString systName = param.Name;

  //========================================================
  //==== Luminosity of prescaled triggers
  //========================================================

  // Luminosity
  if(DataYear==2016){
    MuonLumi1 = 7.408, MuonLumi2 = 7.801, MuonLumi3 = 216.748;
    ElectronLumi1 = 6.988, ElectronLumi2 = 14.851 , ElectronLumi3 = 58.639, ElectronLumi4 = 62.808;
  }
  if(DataYear==2017){
    MuonLumi1 = 4.612, MuonLumi2 = 2.903, MuonLumi3 = 65.943;
    ElectronLumi1 = 3.973, ElectronLumi2 = 27.698, ElectronLumi3 = 35.594, ElectronLumi4 = 43.468;
  }
  if(DataYear==2018){
    MuonLumi1 = 2.696, MuonLumi2 = 8.561, MuonLumi3 = 45.781;
    ElectronLumi1 = 6.412, ElectronLumi2 = 38.849, ElectronLumi3 = 38.861, ElectronLumi4 = 38.906;
  }

  // Muon : ISR
  if(param.Muon_Tight_ID.Contains("ISRTight")){
    if(DataYear==2016){
      SFMuonLumi1 = 0.727537, SFMuonLumi2 = 1.30377, SFMuonLumi3 = 0.96938;
    }
    if(DataYear==2017){
      SFMuonLumi1 = 1.13727, SFMuonLumi2 = 1.34793, SFMuonLumi3 = 0.999073;
    }
    if(DataYear==2018){
      SFMuonLumi1 = 1.98061, SFMuonLumi2 = 1.10568, SFMuonLumi3 = 0.965261;
    }
  }

  // Electron : ISR
  if(param.Electron_Tight_ID.Contains("ISRTight")){
    if(DataYear==2016){
      SFElectronLumi1 = 1.17675, SFElectronLumi2 = 1.09899, SFElectronLumi3 = 1.04746, SFElectronLumi4 = 1.03955;
    }
    if(DataYear==2017){
      SFElectronLumi1 = 1.18292, SFElectronLumi2 = 1.04972, SFElectronLumi3 = 1.04972, SFElectronLumi4 = 0.954682;
    }
    if(DataYear==2018){
      SFElectronLumi1 = 1.07043, SFElectronLumi2 = 1.07253, SFElectronLumi3 = 1.07253, SFElectronLumi4 = 0.915935;
    }
  }

  Event ev = GetEvent();

  //========================================================
  //==== No Cut
  //========================================================

  //JSFillHist(param.Name, "NoCut_"+param.Name, 0., 1., 1, 0., 1.);

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //========================================================
  //==== Trigger
  //========================================================

  //if(! (ev.PassTrigger(MuonTriggers) )) return;

  //========================================================
  //==== Copy AllObjects
  //========================================================

  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;
  vector<Gen> gens = GetGens();

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  double dphiCut = 2.5;
  //double jetPtCut_syst = 40.;
  double PtRatioCut = 1.;

  /*if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnUp){
    this_AllMuons = ScaleMuons( this_AllMuons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnDown){
    this_AllMuons = ScaleMuons( this_AllMuons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResUp){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResDown){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnUp){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnDown){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
  }
  else{
    cout << "[HNtypeI_QCDFakeRateISR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  if(param.fakesyst_ == AnalyzerParameter::FakeCentral){

  }

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  //==== Leptons
  vector<Muon> muons_tight, muons_loose, muons_veto, muons_POGTight, muons_POGLoose;
  muons_tight.clear();
  muons_loose.clear();
  muons_veto.clear();
  muons_POGTight.clear();
  muons_POGLoose.clear();

  vector<Electron> electrons_tight, electrons_loose, electrons_veto;
  electrons_tight.clear();
  electrons_loose.clear();
  electrons_veto.clear();

  muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, MuonPtCut1, 2.4);
  muons_veto  = SelectMuons(this_AllMuons, param.Muon_Veto_ID, MuonPtCut1, 2.4);
  muons_POGTight = SelectMuons(this_AllMuons, "POGTight", 10., 2.4);
  muons_POGLoose = SelectMuons(this_AllMuons, "POGLoose", 5., 2.4);

  electrons_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  electrons_loose = SelectElectrons(this_AllElectrons, param.Electron_Loose_ID, ElectronPtCut1, 2.5);
  electrons_veto  = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, ElectronPtCut1, 2.5);

  //==== Truth matching
  vector<Muon> muons_fake;
  vector<Electron> electrons_fake;
  muons_fake.clear();
  electrons_fake.clear();

  //==== Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);
  vector<Jet> jets;
  jets.clear();

  if(muons_veto.size()+electrons_veto.size() == 1) jets = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto); // We use jets only for FR measurement
  //vector<Jet> jets_PUveto = JetsPassPileupMVA(jets_lepveto);
  
  vector<Jet> jets_awayFromMuon;
  vector<Jet> jets_awayFromElectron;
  jets_awayFromMuon.clear();
  jets_awayFromElectron.clear();

  //========================================================
  //==== Sort in pt-order
  //========================================================

  std::sort(muons_tight.begin(), muons_tight.end(), PtComparing);
  std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(muons_POGTight.begin(), muons_POGTight.end(), PtComparing);
  std::sort(muons_POGLoose.begin(), muons_POGLoose.end(), PtComparing);
  std::sort(electrons_tight.begin(), electrons_tight.end(), PtComparing);
  std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  if(muons_veto.size()+electrons_veto.size() == 1) std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);

  //========================================================
  //==== b tagging
  //========================================================

  int Nbjet_medium = 0, Nbjet_lepveto_medium = 0;
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
  }

  if(muons_veto.size()+electrons_veto.size() == 1){
    for(unsigned int ij=0; ij<jets.size(); ij++){
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets.at(ij))) Nbjet_lepveto_medium++;
    }
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv_central = ev.GetMETVector();

  /*if(muons_loose.size()==1 || electrons_loose.size()==1){
    METv = UpdateMETMuon(METv, muons_loose);
    METv = UpdateMETElectron(METv, electrons_loose);
  }*/

  double MET = METv_central.Pt();
  double METPhi = METv_central.Phi();

  //========================================================    
  //==== Define particles, variables
  //========================================================

  double muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1.;
  double mu_tight_iso = 0.05;

  double el_tight_iso = 0.;   

  // POG cut-based Medium  
  // barrel : 0.0478+0.506/pT, endcap : 0.0658+0.963/pT
  // POG cut-based Tight
  // barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT

  //double pi = 3.14159265358979323846;
  //double MZ = 91.1876;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;
  double jet_emfraction = 0.;

  double trigLumi = 1.; 
  double jetPtCut = 40.;

  double ptcone_mu = 0.;
  double ptcone_el = 0.;
  //double trkiso_Pt = 0.;
  double relTrkIso_MiniAODPt = 0.;
  //double ptcone_mu1 = 0.;
  TString PtConeRange = "";
  Particle ZCand, METv;

  /*Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2);*/

  //========================================================
  //==== Muon Fake Rate Measurement
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muonIDSF = 1., muonIsoSF = 1.;
    if(!(muons_loose.size()>0)) break;
    if(!ev.PassTrigger(MuonTriggers)) break;
    if(!isMuon) break;

    // Fake rate measurement region
    if(it_rg == 0){

      if(!(muons_loose.size()==1 && electrons_loose.size()==0)) continue;
      if(!(muons_veto.size()==1 && electrons_veto.size()==0)) continue;
      if(!(jets.size() >= 1)) continue;

      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium", Nbjet_medium, weight, 10, 0., 10.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_LepVeto", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      if(Nbjet_medium > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_gt0", Nbjet_medium, weight, 10, 0., 10.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_LepVeto_gt0", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      }

      // MET
      METv = UpdateMETMuon(METv_central, muons_loose);
      MET = METv.Pt();
      METPhi = METv.Phi();

      // Set up pTcone
      ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
      relTrkIso_MiniAODPt = muons_loose.at(0).TrkIso()/muons_loose.at(0).MiniAODPt();
      //ptcone_mu1 = muons_loose.at(0).Pt()*(1.+std::max(0., muons_loose.at(0).RelIso()-mu_tight_iso));
      //FillHist("PtCone_ratio", ptcone_mu1/ptcone_mu, weight, 20, 0., 2.);

      // Truth matching
      muons_fake.clear();
      muons_fake = MuonFakeOnly(muons_loose, gens);
      if(!(muons_fake.size() == 1)) continue;

      // Event weights except trigger luminosity
      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        if(param.Muon_Tight_ID.Contains("ISRTight")){
          muonIDSF = mcCorr->MuonID_SF("NUM_TightID_DEN_genTracks", muons_loose.at(0).Eta(), muons_loose.at(0).MiniAODPt(), 0);
          if(muons_tight.size() > 0) muonIsoSF = mcCorr->MuonISO_SF("NUM_TightRelIso_DEN_TightIDandIPCut", muons_tight.at(0).Eta(), muons_tight.at(0).MiniAODPt(), 0);
          else muonIsoSF = 1.;
        }
        else{
          muonIDSF = 1.;
          muonIsoSF = 1.;
        }

        weight *= muonIDSF*muonIsoSF;

      }

      /*if(IsDATA){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIso", muons_loose.at(0).RelIso(), weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIsoTrkIso2D", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
        if(muons_tight.size() > 0) FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PFIsoTrkIso2D", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
      }*/

      trigLumi = 1.;
      // One prescaled trigger for each pTcone range, setup lumi
      if(!(ptcone_mu >= MuonPtconeCut1)) continue;
      if(ptcone_mu >= MuonPtconeCut1 && ptcone_mu < MuonPtconeCut2){
        if(!(muons_loose.at(0).Pt() > MuonPtCut1)) continue;
        if(!ev.PassTrigger(MuonTrig1)) continue;
        if(!IsDATA) trigLumi = MuonLumi1*SFMuonLumi1;
        jetPtCut = 50.;
        PtConeRange = "Range0";
      }
      if(ptcone_mu >= MuonPtconeCut2 && ptcone_mu < MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut2)) continue;
        if(!ev.PassTrigger(MuonTrig2)) continue;
        if(!IsDATA) trigLumi = MuonLumi2*SFMuonLumi2;
        PtConeRange = "Range1";
      }
      if(ptcone_mu >= MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut3)) continue;
        if(!ev.PassTrigger(MuonTrig3)) continue;
        if(!IsDATA) trigLumi = MuonLumi3*SFMuonLumi3;
        PtConeRange = "Range2";
      }

      weight *= trigLumi;

      // For checking TrkIso and RelIso
      if(ev.PassTrigger(MuonTrig1)){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_PassMu3", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIso_PassMu3", muons_loose.at(0).RelIso(), weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIsoTrkIso2D_PassMu3", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
        if(muons_tight.size() > 0) FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PFIsoTrkIso2D_PassMu3", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
      }
      if(ev.PassTrigger(MuonTrig2)){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_PassMu8", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIso_PassMu8", muons_loose.at(0).RelIso(), weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIsoTrkIso2D_PassMu8", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
        if(muons_tight.size() > 0) FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PFIsoTrkIso2D_PassMu8", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
      }
      if(ev.PassTrigger(MuonTrig3)){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_PassMu17", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIso_PassMu17", muons_loose.at(0).RelIso(), weight, 100, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PFIsoTrkIso2D_PassMu17", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
        if(muons_tight.size() > 0) FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PFIsoTrkIso2D_PassMu17", muons_loose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 100, 0., 1., 100, 0., 1.);
      }


      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoDijet", muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoDijet_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      // Away jet selection
      jets_awayFromMuon.clear();
      jets_awayFromMuon = JetsAwayFromLepton(jets, muons_loose.at(0), dphiCut);
      std::sort(jets_awayFromMuon.begin(), jets_awayFromMuon.end(), PtComparing);

      //for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        // define dphi between a jet and the loose lepton
      //  dphi = fabs(jets.at(ijet).Phi() - muons_loose.at(0).Phi());
      //  if(dphi > pi) dphi = 2.*pi-dphi;
        //dphi = fabs(muons_loose.at(0).DeltaPhi(jets.at(ijet)));
      //  FillHist(regions.at(it_rg)+"/dphi_"+PtConeRange, dphi, weight, 32, 0., 3.2);

      //  if(dphi > 2.5) awayjet++; 
      //  if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      //}

      if(!(jets_awayFromMuon.size() > 0)) continue;
      if(!(jets_awayFromMuon.at(0).Pt() > jetPtCut)) continue;
 
      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets_awayFromMuon.at(0).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/MET_NoCut", MET, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/METPhi_NoCut", METPhi, weight, 64, -3.2, 3.2);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/METPhi_NoCut_"+PtConeRange, METPhi, weight, 64, -3.2, 3.2);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Mt_NoCut", Mt, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Ptratio_NoCut", Pt_ratio, weight, 50, 0., 5.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_NoCut", ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Pt_NoCut", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoCut", muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_NoCut", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_NoCut_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_NoCut_"+PtConeRange, 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_NoCut", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Pt_NoCut", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Eta_NoCut", muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Eta_NoCut_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt_NoCut", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt_NoCut_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_NoCut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;

      // Histograms after applying cuts
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone", ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta", muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_TrkIso_MiniAODPt_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_"+PtConeRange, 0.5, weight, 2, 0., 2.);

      if(muons_tight.size() > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Eta", muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Eta_"+PtConeRange, muons_tight.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_TrkIso_MiniAODPt_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        // Loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_IB", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_IB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Pt_IB", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_IB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_IB_NoBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_IB_NoBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_IB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_IB_WithBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_IB_WithBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_IB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_IB", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_IB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Pt_IB", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_IB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_IB_NoBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_IB_NoBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_IB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_IB_WithBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_IB_WithBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_IB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }

          }

        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){
        // Loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_OB", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_OB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Pt_OB", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_OB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_OB_NoBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_OB_NoBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_OB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_OB_WithBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_OB_WithBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_OB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_OB", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_OB_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Pt_OB", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_OB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_OB_NoBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_OB_NoBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_OB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_OB_WithBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_OB_WithBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_OB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }

          }

        }

      }

      // Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.5){
        // Loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_EC", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_PtCone_EC_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Pt_EC", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_EC_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_EC_NoBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_EC_NoBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_EC_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_EC_WithBJet", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_LooseID_Pt_EC_WithBJet", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_EC_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(muons_tight.size() > 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_EC", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_PtCone_EC_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_TightID_Pt_EC", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_Events_EC_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_EC_NoBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_EC_NoBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_EC_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_EC_WithBJet", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Muon_TightID_Pt_EC_WithBJet", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/"+btagDirName+"/Number_Events_EC_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
            }

          }

        }

      }

    }

    if(it_rg > 0) continue;

  }



  //========================================================
  //==== Electron Fake Rate Measurement
  //========================================================
 
  for(unsigned int it_rg2=0; it_rg2<regions.size(); it_rg2++){
    weight = 1.;
    if(!(electrons_loose.size()>0)) break;
    if(!ev.PassTrigger(ElectronTriggers)) break;
    if(!isElectron) break;

    // Fake rate measurement region 
    if(it_rg2 == 0){

      if(!(muons_loose.size()==0 && electrons_loose.size()==1)) continue;
      if(!(muons_veto.size()==0 && electrons_veto.size()==1)) continue;
      if(!(jets.size() >= 1)) continue;

      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium", Nbjet_medium, weight, 10, 0., 10.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_LepVeto", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      if(Nbjet_medium > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_gt0", Nbjet_medium, weight, 10, 0., 10.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_LepVeto_gt0", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      }

      // MET
      METv = UpdateMETElectron(METv_central, electrons_loose);
      MET = METv.Pt();
      METPhi = METv.Phi();

      // Set up pTcone
      el_tight_iso = 0.08; // 2016 or PGO MVA

      if(param.Electron_Tight_ID.Contains("ISR")){ // POG CB Medium
        el_tight_iso = 0.0478+0.506/electrons_loose.at(0).UncorrPt();
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0658+0.963/electrons_loose.at(0).UncorrPt();
      }

      if(param.Electron_Tight_ID.Contains("HNTight")){ // POG CB Tight
        el_tight_iso = 0.0287+0.506/electrons_loose.at(0).UncorrPt();
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons_loose.at(0).UncorrPt();
      }

      // For the same bin in 2016 analysis
      //if(param.Electron_Tight_ID.Contains("2016")) ElectronPtconeCut2 = 23.;
      //else ElectronPtconeCut2 = 25.;

      ptcone_el = electrons_loose.at(0).CalcPtCone(electrons_loose.at(0).RelIso(), el_tight_iso);
      
      // Truth matching
      electrons_fake.clear();
      electrons_fake = ElectronFakeOnly(electrons_loose, gens);
      if(!(electrons_fake.size() == 1)) continue;

      // Event weights except trigger luminosity
      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        electronRecoSF = mcCorr->ElectronReco_SF(electrons_loose.at(0).scEta(), electrons_loose.at(0).UncorrPt(), 0);

        if(param.Electron_Tight_ID.Contains("ISRTight")){
          if(electrons_tight.size() > 0) electronIDSF = mcCorr->ElectronID_SF("passMediumID", electrons_tight.at(0).scEta(), electrons_tight.at(0).UncorrPt(), 0);
          else electronIDSF = 1.;
        }
        else{
          electronIDSF = 1.;
        }

        weight *= electronRecoSF*electronIDSF;

      }

      trigLumi = 1.;
      // One prescaled trigger for each PtCone range, setup lumi
      if(!(ptcone_el >= ElectronPtconeCut1)) continue;
      if(ptcone_el >= ElectronPtconeCut1 && ptcone_el < ElectronPtconeCut2){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut1)) continue;
        if(!ev.PassTrigger(ElectronTrig1)) continue;
        if(!IsDATA) trigLumi = ElectronLumi1*SFElectronLumi1;
        PtConeRange = "Range0";
      }
      if(ptcone_el >= ElectronPtconeCut2 && ptcone_el < ElectronPtconeCut3){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut2)) continue;      
        if(!ev.PassTrigger(ElectronTrig2)) continue;
        if(!IsDATA) trigLumi = ElectronLumi2*SFElectronLumi2;
        PtConeRange = "Range1";
      }
      if(ptcone_el >= ElectronPtconeCut3 && ptcone_el < ElectronPtconeCut4){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut3)) continue;
        if(!ev.PassTrigger(ElectronTrig3)) continue;
        if(!IsDATA) trigLumi = ElectronLumi3*SFElectronLumi3;
        PtConeRange = "Range2";
      }
      if(ptcone_el >= ElectronPtconeCut4){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut4)) continue;
        if(!ev.PassTrigger(ElectronTrig4)) continue;
        if(!IsDATA) trigLumi = ElectronLumi4*SFElectronLumi4;
        PtConeRange = "Range3";
      }

      weight *= trigLumi;

      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoDijet", electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoDijet_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);        
      
      // Away jet selection
      jets_awayFromElectron.clear();
      jets_awayFromElectron = JetsAwayFromLepton(jets, electrons_loose.at(0), dphiCut);
      std::sort(jets_awayFromElectron.begin(), jets_awayFromElectron.end(), PtComparing);

      //for(unsigned int ijet=0; ijet<jets.size(); ijet++){
        //dphi = fabs(jets.at(ijet).Phi() - electrons_loose.at(0).Phi());
        //if(dphi > pi) dphi = 2.*pi-dphi;
      //  dphi = fabs(electrons_loose.at(0).DeltaPhi(jets.at(ijet)));
      //  FillHist("dphi_"+PtConeRange+"_"+regions.at(it_rg2), dphi, weight, 32, 0., 3.2);

      //  if(dphi > 2.5) awayjet++;
      //  if(dphi > 2.5 && awayjet == 1) leadingjet = ijet;
      //}

      if(!(jets_awayFromElectron.size() > 0)) continue;
      if(!(jets_awayFromElectron.at(0).Pt() > jetPtCut)) continue;

      Mt = MT(electrons_loose.at(0), METv);
      Pt_ratio = jets_awayFromElectron.at(0).Pt()/electrons_loose.at(0).Pt();
      jet_emfraction = jets_awayFromElectron.at(0).ChargedEmEnergyFraction();

      // Histograms before applying cuts
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/MET_NoCut", MET, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/METPhi_NoCut", METPhi, weight, 64, -3.2, 3.2);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/METPhi_NoCut_"+PtConeRange, METPhi, weight, 64, -3.2, 3.2);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Mt_NoCut", Mt, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Ptratio_NoCut", Pt_ratio, weight, 50, 0., 5.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Jet_ChargedEmEnergyFraction", jet_emfraction, weight, 100, 0., 1.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Jet_ChargedEmEnergyFraction_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_NoCut", ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Pt_NoCut", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoCut", electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoCut_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_NoCut_"+PtConeRange, 0.5, weight, 2, 0., 2.); 

      if(electrons_tight.size() > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_NoCut", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Pt_NoCut", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Eta_NoCut", electrons_tight.at(0).scEta(), weight, 50, -2.5, 2.5);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Eta_NoCut_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 50, -2.5, 2.5);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_NoCut_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;

      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Jet_ChargedEmEnergyFraction_WithCuts", jet_emfraction, weight, 100, 0., 1.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Jet_ChargedEmEnergyFraction_WithCuts_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);

      if(!(jet_emfraction < 0.65)) continue;

      // Histograms after applying cuts
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone", ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta", electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_"+PtConeRange, 0.5, weight, 2, 0., 2.);
      if(electrons_tight.size() > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Eta", electrons_tight.at(0).scEta(), weight, 50, -2.5, 2.5);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Eta_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 50, -2.5, 2.5);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_"+PtConeRange, 1.5, weight, 2, 0., 2.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(electrons_loose.at(0).scEta()) < 0.8){
        // Loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_IB", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_IB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Pt_IB", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_IB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_IB_NoBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_IB_NoBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_IB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_IB_WithBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_IB_WithBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_IB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_IB", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_IB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Pt_IB", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_IB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_IB_NoBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_IB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_IB_NoBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_IB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_IB_WithBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_IB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_IB_WithBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_IB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }

          }

        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(electrons_loose.at(0).scEta()) >= 0.8 && fabs(electrons_loose.at(0).scEta()) < 1.479){
        // Loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_OB", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_OB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Pt_OB", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_OB_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_OB_NoBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_OB_NoBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_OB_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_OB_WithBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_OB_WithBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_OB_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_OB", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_OB_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Pt_OB", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_OB_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_OB_NoBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_OB_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_OB_NoBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_OB_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_OB_WithBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_OB_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_OB_WithBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_OB_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }

          }

        }

      }

      // Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(electrons_loose.at(0).scEta()) >= 1.479 && fabs(electrons_loose.at(0).scEta()) < 2.5){
        // Loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_EC", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_PtCone_EC_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Pt_EC", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_EC_"+PtConeRange, 0.5, weight, 2, 0., 2.);

        if(systName == "FakeCentral"){

          if(Nbjet_medium == 0){
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_EC_NoBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_EC_NoBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_EC_NoBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }
          else{
            btagDirName = "btagEvent";
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_EC_WithBJet", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_LooseID_Pt_EC_WithBJet", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_EC_WithBJet_"+PtConeRange, 0.5, weight, 2, 0., 2.);
          }

        }

        // Tight ID
        if(electrons_tight.size() > 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_EC", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_PtCone_EC_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_TightID_Pt_EC", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_Events_EC_"+PtConeRange, 1.5, weight, 2, 0., 2.);

          if(systName == "FakeCentral"){

            if(Nbjet_medium == 0){
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_EC_NoBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_EC_NoBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_EC_NoBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_EC_NoBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }
            else{
              btagDirName = "btagEvent";
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_EC_WithBJet", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_PtCone_EC_WithBJet_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Electron_TightID_Pt_EC_WithBJet", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/"+btagDirName+"/Number_Events_EC_WithBJet_"+PtConeRange, 1.5, weight, 2, 0., 2.);
            }

          }

        }

      }

    }

    if(it_rg2 > 0) continue;

  }

}

