#include "HNtypeI_FakeRate.h"

HNtypeI_FakeRate::HNtypeI_FakeRate(){

}

void HNtypeI_FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  RunIso  = HasFlag("RunIso");
  RunNorm = HasFlag("RunNorm");
  RunSF   = HasFlag("RunSF");
  RunPt   = HasFlag("RunPt");

  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunIso = " << RunIso << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunNorm = " << RunNorm << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSF = " << RunSF << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunPt = " << RunPt << endl;

  MuonTightIDs     = {"HNTightV2"};
  MuonLooseIDs     = {"HNLooseV2"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_FakeRate.h 
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
  if(RunPt) MuonPtconeCut1 = MuonPtCut1, MuonPtconeCut2 = MuonPtCut2, MuonPtconeCut3 = MuonPtCut3; 
  else MuonPtconeCut1 = 5., MuonPtconeCut2 = 15., MuonPtconeCut3 = 30.;

  // DoubleEG (2016), SingleElectron (2017), EGamma (2018)
  ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; 
  ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  ElectronTrig3 = "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v";
  ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  if(DataYear==2016) ElectronTrig17L = "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  else ElectronTrig17L = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

  ElectronTriggers.push_back(ElectronTrig1);
  ElectronTriggers.push_back(ElectronTrig2);
  ElectronTriggers.push_back(ElectronTrig3);
  ElectronTriggers.push_back(ElectronTrig4);
  ElectronTriggers.push_back(ElectronTrig17L);
  ElectronPtCut1 = 10., ElectronPtCut2 = 15., ElectronPtCut3 = 20., ElectronPtCut4 = 25.;
  if(RunPt) ElectronPtconeCut1 = ElectronPtCut1, ElectronPtconeCut2 = ElectronPtCut2, ElectronPtconeCut3 = ElectronPtCut3, ElectronPtconeCut4 = ElectronPtCut4;
  else ElectronPtconeCut1 = 15., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 45.;

  //cout << "[HNtypeI_FakeRate::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_FakeRate::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== b tagging
  //==== Qdd taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_FakeRate::~HNtypeI_FakeRate(){

  //==== Destructor of this Analyzer

}

void HNtypeI_FakeRate::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_FakeRate.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllElectrons = GetAllElectrons();
  AllMuons = GetAllMuons();
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_FakeRate.h
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

    if(RunSyst){
      for(int it_syst=3; it_syst<AnalyzerParameter::NFakeSyst; it_syst++){
        param.fakesyst_ = AnalyzerParameter::FakeSyst(it_syst);
        //param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        param.Name  = "FakeSyst_"+param.GetFakeSystType();
        executeEventFromParameter(param);
      }
    }

  }

}

void HNtypeI_FakeRate::executeEventFromParameter(AnalyzerParameter param){

  TString MuonIDname = "MuonHNRun2";

  TString ElectronIDname = "ElectronHNRun2";

  vector<TString> regions = {"FR", "DY", "Wjet"};
  TString btagDirName = "";

  // Boolean : primary datasets
  bool isMuon = false, isElectron = false;
  if(IsDATA){
    if(DataStream.Contains("SingleMuon") || DataStream.Contains("DoubleMuon")) isMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")) isElectron = true;
  }

  TString systName = param.Name;

  //========================================================
  //==== Luminosity of prescaled triggers
  //========================================================

  // Luminosity
  if(DataYear==2016){
    MuonLumi1 = 7.408, MuonLumi2 = 7.801, MuonLumi3 = 216.748;
    ElectronLumi1 = 6.988, ElectronLumi2 = 14.851 , ElectronLumi3 = 62.761, ElectronLumi4 = 62.808, ElectronLumi17L = 58.639;  // Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v : 58.639
  }
  if(DataYear==2017){
    MuonLumi1 = 4.612, MuonLumi2 = 2.903, MuonLumi3 = 65.943;
    ElectronLumi1 = 3.973, ElectronLumi2 = 27.698, ElectronLumi3 = 35.594, ElectronLumi4 = 43.468, ElectronLumi17L = ElectronLumi2;
  }
  if(DataYear==2018){
    MuonLumi1 = 2.696, MuonLumi2 = 8.561, MuonLumi3 = 45.781;
    ElectronLumi1 = 6.412, ElectronLumi2 = 38.849, ElectronLumi3 = 38.861, ElectronLumi4 = 38.906, ElectronLumi17L = ElectronLumi2;
  }

  // Muon, Electron : 2016 HN
  /*if(param.Muon_Tight_ID.Contains("2016")){
    if(DataYear==2016){
      MuonLumi1 = 7.408*0.679317, MuonLumi2 = 7.801*1.26668, MuonLumi3 = 216.748*0.940637;
      ElectronLumi1 = 6.988*1.04596, ElectronLumi2 = 14.851*0.968804, ElectronLumi3 = 62.761*0.933706, ElectronLumi4 = 62.808*0.933263;
    }
    if(DataYear==2017){
      MuonLumi1 = 4.612*1.14633, MuonLumi2 = 2.903*1.32612, MuonLumi3 = 65.943*0.981435;
      ElectronLumi1 = 3.973*1.05833, ElectronLumi2 = 27.698*0.913181, ElectronLumi3 = 35.594*0.854491, ElectronLumi4 = 43.468*0.821393;
    }
    if(DataYear==2018){
      MuonLumi1 = 2.696*1.89483, MuonLumi2 = 8.561*1.0494, MuonLumi3 = 45.781*0.915734;
      ElectronLumi1 = 6.412*0.932546, ElectronLumi2 = 38.849*0.928263, ElectronLumi3 = 38.861*0.797685, ElectronLumi4 = 38.906*0.793358;
    }
  }*/

  // Muon : HNTight
  if(param.Muon_Tight_ID.Contains("HNTightV2")){
    if(DataYear==2016){
      SFMuonLumi1 = 0.713927, SFMuonLumi2 = 1.31533, SFMuonLumi3 = 0.977825;
    }
    if(DataYear==2017){
      SFMuonLumi1 = 1.17044, SFMuonLumi2 = 1.35542, SFMuonLumi3 = 1.00362;
    }
    if(DataYear==2018){
      SFMuonLumi1 = 2.0007, SFMuonLumi2 = 1.09757, SFMuonLumi3 = 0.954471;
    }
  }

  // Electron : HNTight
  if(param.Electron_Tight_ID.Contains("HNTightV2")){
    if(DataYear==2016){
      SFElectronLumi1 = 1.20126, SFElectronLumi2 = 1.10488, SFElectronLumi3 = 1.06424, SFElectronLumi4 = 1.06333, SFElectronLumi17L = 1.06763;
    }
    if(DataYear==2017){
      SFElectronLumi1 = 1.23103, SFElectronLumi2 = 1.06039, SFElectronLumi3 = 0.990935, SFElectronLumi4 = 0.958795, SFElectronLumi17L = SFElectronLumi2;
    }
    if(DataYear==2018){
      SFElectronLumi1 = 1.05831, SFElectronLumi2 = 1.0892, SFElectronLumi3 = 0.939894, SFElectronLumi4 = 0.933807, SFElectronLumi17L = SFElectronLumi2;
    }
  }

  // Without normalization scale factors
  if(RunNorm && !RunSF){
    if(DataYear==2016){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1., SFElectronLumi17L = 1.;
    }
    if(DataYear==2017){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1., SFElectronLumi17L = SFElectronLumi2;
    }
    if(DataYear==2018){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1., SFElectronLumi17L = SFElectronLumi2;
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

  double jetPtCut_syst = 40.;
  double dPhiCut = 2.5;
  double PtRatioCut = 1.;

  if(param.fakesyst_ == AnalyzerParameter::FakeCentral){

  }
  else if(param.fakesyst_ == AnalyzerParameter::AwayJetPt20){
    jetPtCut_syst = 20.;
  }
  else if(param.fakesyst_ == AnalyzerParameter::AwayJetPt30){
    jetPtCut_syst = 30.;
  }
  else if(param.fakesyst_ == AnalyzerParameter::AwayJetPt60){
    jetPtCut_syst = 60.;
  }
  else if(param.fakesyst_ == AnalyzerParameter::AwayJetPt100){
    jetPtCut_syst = 100.;
  }
  else if(param.fakesyst_ == AnalyzerParameter::AwayJetPt200){
    jetPtCut_syst = 200.;
  }
  else if(param.fakesyst_ == AnalyzerParameter::dPhi1){
    dPhiCut = 1.5;
  }
  else if(param.fakesyst_ == AnalyzerParameter::dPhi1){
    dPhiCut = 2.0;
  }
  else if(param.fakesyst_ == AnalyzerParameter::dPhi1){
    dPhiCut = 3.0;
  }
  else if(param.fakesyst_ == AnalyzerParameter::PtRatioUp){
    PtRatioCut = 1.2;
  }
  else if(param.fakesyst_ == AnalyzerParameter::PtRatioDown){
    PtRatioCut = 0.8;
  }
  /*else if(param.syst_ == AnalyzerParameter::JetResUp){
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
  }*/
  else{
    cout << "[HNtypeI_FakeRate::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
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
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

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

  double MET = METv_central.Pt();
  double METPhi = METv_central.Phi();

  //========================================================    
  //==== Define particles, variables
  //========================================================

  double muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1.;
  double mu_tight_iso = 0.07;
  double el_tight_iso = 0.;   

  // POG cut-based Medium  
  // barrel : 0.0478+0.506/pT, endcap : 0.0658+0.963/pT
  // POG cut-based Tight
  // barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT

  //double pi = 3.14159265358979323846;
  double MZ = 91.1876;
  double weight = 1.;
  double Mt = 0.;
  double Pt_ratio = 0.;
  double jet_emfraction = 0.;

  double trigLumi = 1.;
  double jetPtCut = jetPtCut_syst;

  bool IsAwayJetBtag = false, IsNearbyJetBtag = false;

  double ptcone_mu = 0., ptcone_el = 0.;
  //double relIsoLoose_mu = 0., relIsoLoose_el = 0.;
  //double trkiso_Pt = 0.;
  double relTrkIso_MiniAODPt = 0.;
  //double ptcone_mu1 = 0.;
  TString PtConeRange = "";
  Particle ZCand, METv;

  // Track Iso vs PF Iso (Use DoubleMuon w/o skim)
  if(RunIso){

    if(muons_POGLoose.size()==1 && muons_POGTight.size()==1){

      FillHist("Muon_MiniAODPt", muons_POGLoose.at(0).MiniAODPt(), weight, 500, 0., 500.);

      relTrkIso_MiniAODPt = muons_POGLoose.at(0).TrkIso()/muons_POGLoose.at(0).MiniAODPt();
      FillHist("Muon_TrkIso", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
      FillHist("Muon_PFIso", muons_POGLoose.at(0).RelIso(), weight, 100, 0., 1.);
      FillHist("Muon_PFIsoTrkIso2D", muons_POGLoose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 1000, 0., 1., 1000, 0., 1.);

      if(ev.PassTrigger(MuonTrig1)){
        if(muons_POGLoose.at(0).Pt() > 5.){
          FillHist("Muon_TrkIso_Mu3", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
          FillHist("Muon_PFIso_Mu3", muons_POGLoose.at(0).RelIso(), weight, 100, 0., 1.);
          FillHist("Muon_PFIsoTrkIso2D_Mu3", muons_POGLoose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 1000, 0., 1., 1000, 0., 1.);
        }
      }
 
      if(ev.PassTrigger(MuonTrig2)){
        if(muons_POGLoose.at(0).Pt() > 10.){
          FillHist("Muon_TrkIso_Mu8", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
          FillHist("Muon_PFIso_Mu8", muons_POGLoose.at(0).RelIso(), weight, 100, 0., 1.);
          FillHist("Muon_PFIsoTrkIso2D_Mu8", muons_POGLoose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 1000, 0., 1., 1000, 0., 1.);
        }
      }

      if(ev.PassTrigger(MuonTrig3)){
        if(muons_POGLoose.at(0).Pt() > 20.){
          FillHist("Muon_TrkIso_Mu17", relTrkIso_MiniAODPt, weight, 100, 0., 1.);
          FillHist("Muon_PFIso_Mu17", muons_POGLoose.at(0).RelIso(), weight, 100, 0., 1.);
          FillHist("Muon_PFIsoTrkIso2D_Mu17", muons_POGLoose.at(0).RelIso(), relTrkIso_MiniAODPt, weight, 1000, 0., 1., 1000, 0., 1.);
        }
      }

    }

  }

  /*Gen gen_test;
  FillHist("gen_mother", gen_test.MotherIndex(), weight, 4, -2, 2);
  FillHist("gen_pid", gen_test.PID(), weight, 4, -2, 2);
  FillHist("gen_status", gen_test.Status(), weight, 4, -2, 2);*/

  FillHist("MET_NoCut", MET, weight, 500, 0., 500.);
  FillHist("METPhi_NoCut", METPhi, weight, 64, -3.2, 3.2);

  //========================================================
  //==== Muon Fake Rate Measurement
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muonIDSF = 1., muonIsoSF = 1.;
    if(!(muons_loose.size()>0)) break;
    if(!ev.PassTrigger(MuonTriggers)) break;
    if(IsDATA){ if(!isMuon) break; }

    // Fake rate measurement region
    if(it_rg == 0){

      if(RunNorm) continue;

      if(!(muons_loose.size()==1 && electrons_loose.size()==0)) continue;
      if(!(muons_veto.size()==1 && electrons_veto.size()==0)) continue;
      if(!(jets.size() >= 1)) continue;

      /*FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium", Nbjet_medium, weight, 10, 0., 10.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_LepVeto", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      if(Nbjet_medium > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_gt0", Nbjet_medium, weight, 10, 0., 10.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Number_BJets_Medium_LepVeto_gt0", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      }*/

      // MET
      METv = UpdateMETMuon(METv_central, muons_loose);
      MET = METv.Pt();
      METPhi = METv.Phi();

      // Set up pTcone
      if(RunPt) ptcone_mu = muons_loose.at(0).Pt();
      else ptcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);

      //relIsoLoose_mu = std::max(0., muons_loose.at(0).RelIso() - mu_tight_iso);
      relTrkIso_MiniAODPt = muons_loose.at(0).TrkIso()/muons_loose.at(0).MiniAODPt();
      //ptcone_mu1 = muons_loose.at(0).Pt()*(1.+std::max(0., muons_loose.at(0).RelIso()-mu_tight_iso));
      //FillHist("PtCone_ratio", ptcone_mu1/ptcone_mu, weight, 20, 0., 2.);

      // Truth matching
      muons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons_loose, gens);
      if(!(muons_prompt.size() == 1)) continue;

      /*if(muons_tight.size() == 0){
        if(ev.PassTrigger(MuonTrig1)){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt5to10_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 5, 5., 10.);
        }
        if(ev.PassTrigger(MuonTrig2)){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt10to15_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 5, 10., 15.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt15to20_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 5, 15., 20.);
        }
        if(ev.PassTrigger(MuonTrig3)){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt20to30_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 20., 30.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt30to40_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 30., 40.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt40to50_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 40., 50.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt50to60_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 50., 60.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt60to70_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 60., 70.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_RelIsoLoose_Pt70to80_2D", relIsoLoose_mu, muons_loose.at(0).Pt(), weight, 45, 0., 0.45, 10, 70., 80.);
        }
      }*/

      // Event weights except trigger luminosity
      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        if(param.Muon_Tight_ID.Contains("HNTight")){
          if(muons_tight.size() > 0) muonIDSF = mcCorr->MuonID_SF_HNtypeI(param.Muon_Tight_ID, muons_tight.at(0).Eta(), muons_tight.at(0).MiniAODPt(), 0);
          else muonIDSF = 1.;
          muonIsoSF = 1.;
        }
        else{
          muonIDSF  = 1.;
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
        if(jetPtCut_syst == 40.) jetPtCut = 50.;
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
      /*if(ev.PassTrigger(MuonTrig1)){
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
      }*/

      //FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoDijet", muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);
      //FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"/Muon_LooseID_Eta_NoDijet_"+PtConeRange, muons_loose.at(0).Eta(), weight, 50, -2.5, 2.5);

      // Away jet selection
      jets_awayFromMuon.clear();
      jets_awayFromMuon = JetsAwayFromLepton(jets, muons_loose.at(0), dPhiCut);
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

      // B tagging for the away/nearby jet 
      IsAwayJetBtag = false, IsNearbyJetBtag = false;

      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_awayFromMuon.at(0))) IsAwayJetBtag = true;

      for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
        if(jets_nolepveto.at(ij).DeltaR(muons_loose.at(0)) < 0.4){
          if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) IsNearbyJetBtag = true;
        }
        if(IsNearbyJetBtag) break;
      }

      Mt = MT(muons_loose.at(0), METv);
      Pt_ratio = jets_awayFromMuon.at(0).Pt()/muons_loose.at(0).Pt();

      // Histograms before applying cuts
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_MET_NoCut", MET, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_METPhi_NoCut", METPhi, weight, 64, -3.2, 3.2);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_METPhi_NoCut_"+PtConeRange, METPhi, weight, 64, -3.2, 3.2);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mt_NoCut", Mt, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Ptratio_NoCut", Pt_ratio, weight, 50, 0., 5.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);

      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_PtCone_NoCut", ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Pt_NoCut", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Eta_NoCut", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_TrkIso_MiniAODPt_NoCut", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_TrkIso_MiniAODPt_NoCut_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 0.5, weight, 3, 0., 3.);

      /*if(muons_tight.size() == 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_PtCone_NoCut", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_Pt_NoCut", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_Eta_NoCut", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_TrkIso_MiniAODPt_NoCut", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_NoTight_TrkIso_MiniAODPt_NoCut_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 1.5, weight, 3, 0., 3.);
      }*/

      if(muons_tight.size() > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_PtCone_NoCut", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_PtCone_NoCut_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Pt_NoCut", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Eta_NoCut", muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Eta_NoCut_"+PtConeRange, muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_TrkIso_MiniAODPt_NoCut", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_TrkIso_MiniAODPt_NoCut_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;

      // Histograms after applying cuts
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Eta", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_TrkIso_MiniAODPt", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Loose_TrkIso_MiniAODPt_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
      FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

      if(muons_tight.size() > 0){
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Eta", muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_Eta_"+PtConeRange, muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_TrkIso_MiniAODPt", relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Muon_Tight_TrkIso_MiniAODPt_"+PtConeRange, relTrkIso_MiniAODPt, weight, 20, 0., 1.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){

        // Passing loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_IB_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){

        // Passing loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_OB_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      // Endcap ( 1.479 < |eta| < 2.4 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.4){

        // Passing loose ID
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Loose_PtCone", ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Loose_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_NoTight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_NoTight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Tight_PtCone", ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Tight_PtCone_"+PtConeRange, ptcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_EC_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

    }

    // DY control region
    if(it_rg == 1){

      if(RunSyst) continue;

      if(!(muons_tight.size()==2 && electrons_tight.size()==0)) continue;
      //if(!(muons_veto.size()==2 && electrons_veto.size()==0)) continue;

      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<muons_tight.size(); i++){

          if(param.Muon_Tight_ID.Contains("HNTight")){
            muonIDSF  = mcCorr->MuonID_SF_HNtypeI(param.Muon_Tight_ID, muons_tight.at(i).Eta(), muons_tight.at(i).MiniAODPt(), 0);
            muonIsoSF = 1.;
          }
          else{
            muonIDSF  = 1.;
            muonIsoSF = 1.;
          }

          weight *= muonIDSF*muonIsoSF;

        }

      }

      ZCand = muons_tight.at(0) + muons_tight.at(1);

      if(ev.PassTrigger(MuonTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }
      if(ev.PassTrigger(MuonTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3*SFMuonLumi3;
 
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      // Event selection
      if(!(muons_tight.at(0).Pt()>20. && muons_tight.at(1).Pt()>10.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }
        
        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }
        
        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

    }

    // W+jets control region
    if(it_rg == 2){

      if(RunSyst) continue;

      if(!(muons_tight.size()==1 && electrons_tight.size()==0)) continue;
      //if(!(muons_veto.size()==1 && electrons_veto.size()==0)) continue;

      // MET
      METv = UpdateMETMuon(METv_central, muons_tight);
      MET = METv.Pt();
      METPhi = METv.Phi();

      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        if(param.Muon_Tight_ID.Contains("HNTight")){
          muonIDSF  = mcCorr->MuonID_SF_HNtypeI(param.Muon_Tight_ID, muons_tight.at(0).Eta(), muons_tight.at(0).MiniAODPt(), 0);
          muonIsoSF = 1.;
        }
        else{
          muonIDSF  = 1.;
          muonIsoSF = 1.;
        }

        weight *= muonIDSF*muonIsoSF;

      }

      Mt = MT(muons_tight.at(0), METv);

      if(ev.PassTrigger(MuonTrig1)){
        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      // Event selection
      if(!(muons_tight.at(0).Pt() > 20.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;
      
      // Histograms for each trigger
      if(ev.PassTrigger(MuonTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(MuonIDname+"/"+systName+"/"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

    }

  }



  //========================================================
  //==== Electron Fake Rate Measurement
  //========================================================
 
  for(unsigned int it_rg2=0; it_rg2<regions.size(); it_rg2++){
    weight = 1.;
    if(!(electrons_loose.size()>0)) break;
    if(!ev.PassTrigger(ElectronTriggers)) break;
    if(IsDATA){ if(!isElectron) break; }

    // Fake rate measurement region 
    if(it_rg2 == 0){

      if(RunNorm) continue;

      if(!(muons_loose.size()==0 && electrons_loose.size()==1)) continue;
      if(!(muons_veto.size()==0 && electrons_veto.size()==1)) continue;
      if(!(jets.size() >= 1)) continue;

      /*FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium", Nbjet_medium, weight, 10, 0., 10.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_LepVeto", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      if(Nbjet_medium > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_gt0", Nbjet_medium, weight, 10, 0., 10.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Number_BJets_Medium_LepVeto_gt0", Nbjet_lepveto_medium, weight, 10, 0., 10.);
      }*/

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

      if(RunPt) ptcone_el = electrons_loose.at(0).Pt();
      else ptcone_el = electrons_loose.at(0).CalcPtCone(electrons_loose.at(0).RelIso(), el_tight_iso);
      //relIsoLoose_el = std::max(0., electrons_loose.at(0).RelIso() - el_tight_iso);

      // Truth matching
      electrons_prompt.clear();
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons_loose, gens);
      if(!(electrons_prompt.size() == 1)) continue;

      /*if(electrons_tight.size() == 0){
        if(ev.PassTrigger(ElectronTrig1)){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt9to15_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 6, 9., 15.);
        }
        if(ev.PassTrigger(ElectronTrig2)){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt15to20_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 5, 15., 20.);
        }
        if(ev.PassTrigger(ElectronTrig3)){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt20to25_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 5, 20., 25.);
        }
        if(ev.PassTrigger(ElectronTrig4)){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt25to30_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 5, 25., 30.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt30to40_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 30., 40.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt40to50_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 40., 50.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt50to60_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 50., 60.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt60to70_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 60., 70.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt70to80_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 70., 80.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_RelIsoLoose_Pt80to90_2D", relIsoLoose_el, electrons_loose.at(0).Pt(), weight, 65, 0., 0.65, 10, 80., 90.);
        }
      }*/

      // Event weights except trigger luminosity
      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        electronRecoSF = mcCorr->ElectronReco_SF(electrons_loose.at(0).scEta(), electrons_loose.at(0).UncorrPt(), 0);

        if(param.Electron_Tight_ID.Contains("HNTight")){
          if(electrons_tight.size() > 0) electronIDSF = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons_tight.at(0).scEta(), electrons_tight.at(0).UncorrPt(), 0);
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

      //FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoDijet", electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);
      //FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"/Electron_LooseID_Eta_NoDijet_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 50, -2.5, 2.5);        
      
      // Away jet selection
      jets_awayFromElectron.clear();
      jets_awayFromElectron = JetsAwayFromLepton(jets, electrons_loose.at(0), dPhiCut);
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

      // B tagging for the away/nearby jet
      IsAwayJetBtag = false, IsNearbyJetBtag = false;

      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_awayFromElectron.at(0))) IsAwayJetBtag = true;

      for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
        if(jets_nolepveto.at(ij).DeltaR(electrons_loose.at(0)) < 0.4){
          if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) IsNearbyJetBtag = true;
        }
        if(IsNearbyJetBtag) break;
      }

      Mt = MT(electrons_loose.at(0), METv);
      Pt_ratio = jets_awayFromElectron.at(0).Pt()/electrons_loose.at(0).Pt();
      jet_emfraction = jets_awayFromElectron.at(0).ChargedEmEnergyFraction();

      // Histograms before applying cuts
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_MET_NoCut", MET, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_METPhi_NoCut", METPhi, weight, 64, -3.2, 3.2);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_METPhi_NoCut_"+PtConeRange, METPhi, weight, 64, -3.2, 3.2);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Mt_NoCut", Mt, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ptratio_NoCut", Pt_ratio, weight, 50, 0., 5.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ptratio_NoCut_"+PtConeRange, Pt_ratio, weight, 50, 0., 5.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction", jet_emfraction, weight, 100, 0., 1.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);

      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_PtCone_NoCut", ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Pt_NoCut", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Eta_NoCut", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Eta_NoCut_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Number_Events_NoCut_"+PtConeRange, 0.5, weight, 3, 0., 3.); 

      if(electrons_tight.size() > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_PtCone_NoCut", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_PtCone_NoCut_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Pt_NoCut", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Eta_NoCut", electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Eta_NoCut_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Number_Events_NoCut_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      // Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(Pt_ratio > PtRatioCut)) continue;

      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_WithCuts", jet_emfraction, weight, 100, 0., 1.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_WithCuts_"+PtConeRange, jet_emfraction, weight, 100, 0., 1.);

      if(!(jet_emfraction < 0.65)) continue;

      // Histograms after applying cuts
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Eta", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Loose_Eta_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

      if(electrons_tight.size() > 0){
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Eta", electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Electron_Tight_Eta_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      // Inner barrel ( |eta| < 0.8 )
      if(fabs(electrons_loose.at(0).scEta()) < 0.8){

        // Passing loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(electrons_tight.size() > 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_IB_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      // Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(electrons_loose.at(0).scEta()) >= 0.8 && fabs(electrons_loose.at(0).scEta()) < 1.479){

        // Passing loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(electrons_tight.size() > 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_OB_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      // Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(electrons_loose.at(0).scEta()) >= 1.479 && fabs(electrons_loose.at(0).scEta()) < 2.5){

        // Passing loose ID
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBtag && IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

          if(IsAwayJetBtag || IsNearbyJetBtag){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Loose_PtCone", ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Loose_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        // Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_NoTight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_NoTight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        // Passing tight ID
        if(electrons_tight.size() > 0){

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBtag && IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_2Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

            if(IsAwayJetBtag || IsNearbyJetBtag){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_1Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Tight_PtCone", ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Tight_PtCone_"+PtConeRange, ptcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_EC_0Btag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

    }

    // DY control region
    if(it_rg2 == 1){

      if(RunSyst) continue;

      if(!(muons_tight.size()==0 && electrons_tight.size()==2)) continue;
      //if(!(muons_veto.size()==0 && electrons_veto.size()==2)) continue;

      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        for(unsigned int i=0; i<electrons_tight.size(); i++){

          electronRecoSF = mcCorr->ElectronReco_SF(electrons_tight.at(i).scEta(), electrons_tight.at(i).UncorrPt(), 0);

          if(param.Electron_Tight_ID.Contains("HNTight")){
            electronIDSF = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons_tight.at(i).scEta(), electrons_tight.at(i).UncorrPt(), 0);
          }
          else{
            electronIDSF = 1.;
          }

          weight *= electronRecoSF*electronIDSF;

        }

      }

      ZCand = electrons_tight.at(0) + electrons_tight.at(1);

      if(ev.PassTrigger(ElectronTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
        }

      }

      if(DataYear==2016){

        if(ev.PassTrigger(ElectronTrig17L)){

          trigLumi = 1.;
          if(!IsDATA) trigLumi = ElectronLumi17L*SFElectronLumi17L;

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass_NoCut", ZCand.M(), weight*trigLumi, 80, 50., 130.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigLumi, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigLumi, 80, 50., 130.);
          }

        }

      }

      // Event selection
      if(!(electrons_tight.at(0).Pt()>25. && electrons_tight.at(1).Pt()>15.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }        

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
          }
        }

      }

      if(DataYear==2016){

        if(ev.PassTrigger(ElectronTrig17L)){

          trigLumi = 1.;
          if(!IsDATA) trigLumi = ElectronLumi17L*SFElectronLumi17L;

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass_Inclusive", ZCand.M(), weight*trigLumi, 50, 70., 120.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
          }

          if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass", ZCand.M(), weight*trigLumi, 50, 70., 120.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events", 1.5, weight*trigLumi, 2, 0., 2.);

            if(muons_veto.size()==0 && electrons_veto.size()==2){
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigLumi, 50, 70., 120.);
              FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events_OnlyTight", 1.5, weight*trigLumi, 2, 0., 2.);
            }
          }

        }

      }

    }


    // W+jets control region
    if(it_rg2 == 2){

      if(RunSyst) continue;

      if(!(muons_tight.size()==0 && electrons_tight.size()==1)) continue;
      //if(!(muons_veto.size()==0 && electrons_veto.size()==1)) continue;

      // MET
      METv = UpdateMETElectron(METv_central, electrons_tight);
      MET = METv.Pt();
      METPhi = METv.Phi();

      if(!IsDATA){

        weight *= weight_norm_1invpb;
        weight *= ev.MCweight();
        weight *= GetPrefireWeight(0);
        weight *= GetPileUpWeight(nPileUp,0);

        electronRecoSF = mcCorr->ElectronReco_SF(electrons_tight.at(0).scEta(), electrons_tight.at(0).UncorrPt(), 0);

        if(param.Electron_Tight_ID.Contains("HNTight")){
          electronIDSF = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons_tight.at(0).scEta(), electrons_tight.at(0).UncorrPt(), 0);
        }
        else{
          electronIDSF = 1.;
        }

        weight *= electronRecoSF*electronIDSF;

      }

      Mt = MT(electrons_tight.at(0), METv);

      if(ev.PassTrigger(ElectronTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
        }

      }

      if(DataYear==2016){

        if(ev.PassTrigger(ElectronTrig17L)){

          trigLumi = 1.;
          if(!IsDATA) trigLumi = ElectronLumi17L*SFElectronLumi17L;

          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_MET_NoCut", MET, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_METPhi_NoCut", METPhi, weight*trigLumi, 64, -3.2, 3.2);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Mt_NoCut", Mt, weight*trigLumi, 500, 0., 500.);

          if(muons_veto.size()==0 && electrons_veto.size()==1){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigLumi, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_MET_NoCut_OnlyTight", MET, weight*trigLumi, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_METPhi_NoCut_OnlyTight", METPhi, weight*trigLumi, 64, -3.2, 3.2);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Mt_NoCut_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          }

        }

      }

      // Event selection
      if(!(electrons_tight.at(0).Pt() > 25.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;

      // Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17M_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigLumi = 1.;
        if(!IsDATA) trigLumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Mt", Mt, weight*trigLumi, 500, 0., 500.);
        FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
        }

      }

      if(DataYear==2016){

        if(ev.PassTrigger(ElectronTrig17L)){

          trigLumi = 1.;
          if(!IsDATA) trigLumi = ElectronLumi17L*SFElectronLumi17L;
  
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Mt", Mt, weight*trigLumi, 500, 0., 500.);
          FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events", 0.5, weight*trigLumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==1){
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Mt_OnlyTight", Mt, weight*trigLumi, 500, 0., 500.);
            FillHist(ElectronIDname+"/"+systName+"/"+regions.at(it_rg2)+"_Ele17L_Number_Events_OnlyTight", 0.5, weight*trigLumi, 2, 0., 2.);
          }

        }

      }

    }

  }

}

