#include "HNtypeI_DY_CR.h"

HNtypeI_DY_CR::HNtypeI_DY_CR(){

}

void HNtypeI_DY_CR::initializeAnalyzer(){

  //==== If you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  RunFake = HasFlag("RunFake");

  cout << "[HNtypeI_DY_CR::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_DY_CR::initializeAnalyzer] RunFake = " << RunFake << endl;

  MuonTightIDs     = {"HNTightV2"};
  MuonLooseIDs     = {"HNLooseV2"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};
  MuonFRNames      = {"HNTightV2"};
  ElectronFRNames  = {"HNTightV2"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_DY_CR.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
  //==== Run number : ~280385 (2016G), 281613~ (2016H)

  MuonTriggers.clear();
  MuonTriggersNoDZ.clear();
  ElectronTriggers.clear();
  EMuTriggers.clear();
  EMuTriggersNoDZ.clear();
  Mu8Ele23Trigger = "";
  Mu23Ele12Trigger = "";
  Mu8Ele23TriggerNoDZ = "";
  Mu23Ele12TriggerNoDZ = "";

  if(DataEra == "2016preVFP"){                                                        // Lumi values of triggers in 2016 data (/pb)

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                     // 19517.523849710 
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                   // 19517.523849710
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                  // 19517.523849710
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                // 19517.523849710
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");        // 19517.523849710
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");        // 19517.523849710
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");       // 19517.523849710

    //Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }
  else if(DataEra == "2016postVFP"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                  // 16812.151722311
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                // 16812.151722311
    MuonTriggersNoDZ.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                 // 8072.032418212  (F+G)
    MuonTriggersNoDZ.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");               // 8072.032418212  (F+G)
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");        // 16812.151722311
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 16812.151722311
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");    // 16812.151722311
    EMuTriggersNoDZ.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");    // 8072.032418212  (F+G)
    EMuTriggersNoDZ.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");   // 8072.032418212  (F+G)

  }
  else if(DataEra == "2017"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }
  else if(DataEra == "2018"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }

  //cout << "[HNtypeI_DY_CR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_DY_CR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B Tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_DY_CR::~HNtypeI_DY_CR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_DY_CR::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_DY_CR.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets());

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_DY_CR.h
  //weight_Prefire = GetPrefireWeight(0);

  //==== Declare AnalyzerParameter

  AnalyzerParameter param;

  //==== Loop over muon IDs

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){

    TString MuonTightID     = MuonTightIDs.at(it_id);
    TString MuonLooseID     = MuonLooseIDs.at(it_id);
    TString MuonVetoID      = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);
    TString MuonFRName      = MuonFRNames.at(it_id);
    TString ElectronFRName  = ElectronFRNames.at(it_id);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;
    param.Name = "Central";

    //==== Muon ID
    param.Muon_Tight_ID       = MuonTightID;
    param.Muon_Loose_ID       = MuonLooseID;
    param.Muon_Veto_ID        = MuonVetoID;
    param.Muon_FR_ID          = MuonFRName;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key         = "AwayJetPt40";  // histname
    param.Muon_ID_SF_Key      = "";
    param.Muon_ISO_SF_Key     = "";
    param.Muon_Trigger_SF_Key = "";
    param.Muon_UsePtCone      = true;

    //==== Electron ID
    param.Electron_Tight_ID       = ElectronTightID;
    param.Electron_Loose_ID       = ElectronLooseID;
    param.Electron_Veto_ID        = ElectronVetoID;
    param.Electron_FR_ID          = ElectronFRName; // ID name in histmap_Electron.txt
    param.Electron_FR_Key         = "AwayJetPt40";  // histname
    param.Electron_ID_SF_Key      = "";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone      = true;

    //==== Jet ID
    param.Jet_ID = "HNTight";
    if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45";

    executeEventFromParameter(param);

    /*if(RunSyst){

      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){

        //==== Everything else remains same, but only change syst_ and parameter name

        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }

    }*/

  }

}

void HNtypeI_DY_CR::executeEventFromParameter(AnalyzerParameter param){

  TString IDName = "HNTightV2";

  vector<TString> regions = {"WZ", "ZG", "Fake", "WG", "ZZ", "WW"};
  vector<TString> channels3L = {"mmm", "mme", "mee", "eee"};
  vector<TString> channels4L = {"mmmm", "mmee", "eeee"};
  TString tight_leptons = "";

  TString systName = param.Name;

  double cutflow_max = 12.;
  int cutflow_bin = 12;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trigger_weight = 0., diel_trigger_weight = 0.;
  int tight_muons = 0, tight_electrons = 0;

  Event ev = GetEvent();

  bool isDoubleMuon = false, isDoubleEG = false;

  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) isDoubleMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) isDoubleEG = true;
  }

  bool passMuMu = false, passEE = false;

  if(DataEra == "2016postVFP"){
    if(IsDATA){
      if(run < 281000) passMuMu = ev.PassTrigger(MuonTriggers) ||  ev.PassTrigger(MuonTriggersNoDZ);
      else passMuMu = ev.PassTrigger(MuonTriggers);
    }
    else passMuMu = ev.PassTrigger(MuonTriggers) ||  ev.PassTrigger(MuonTriggersNoDZ);
  }
  else passMuMu = ev.PassTrigger(MuonTriggers);

  passEE = ev.PassTrigger(ElectronTriggers);

/*
  //=============
  //==== No Cut
  //=============

  FillHist(param.Name+"/NoCut_"+param.Name, 0., 1., 1, 0., 1.);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  //==============
  //==== Trigger
  //==============
  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;



  //======================
  //==== Copy AllObjects
  //======================

  vector<Muon> this_AllMuons = AllMuons;
  vector<Jet> this_AllJets = AllJets;

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  if(param.syst_ == AnalyzerParameter::Central){

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
    cout << "[HNtypeI_DY_CR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 20., 2.4);
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);

  //=======================
  //==== Sort in pt-order
  //=======================

  //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
  std::sort(muons.begin(), muons.end(), PtComparing);
  //==== 2) jets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
  std::sort(jets.begin(), jets.end(), PtComparing);

  int NBJets_NoSF(0), NBJets_WithSF_2a(0);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV,
                                                                     JetTagging::Medium,
                                                                     JetTagging::incl, JetTagging::comb);

  //==== b tagging

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij = 0 ; ij < jets.size(); ij++){

    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::DeepCSV);
    //==== No SF
    if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
    //==== 2a
    if( mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets.at(ij)) ) NBJets_WithSF_2a++;

  }

  
  //=========================
  //==== Event selections..
  //=========================

  //==== dimuon
  if(muons.size() != 2) return;

  //==== leading muon has trigger-safe pt
  if( muons.at(0).Pt() <= TriggerSafePtCut ) return;

  //==== On-Z
  Particle ZCand = muons.at(0) + muons.at(1);
  if(!IsOnZ(ZCand.M(), 15.)) return;

  //===================
  //==== Event weight
  //===================

  double weight = 1.;
  //==== If MC
  if(!IsDATA){

    //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
    //==== So, you have to multiply trigger luminosity
    //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"

    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");

    //==== MCweight is +1 or -1. Should be multiplied if you are using e.g., aMC@NLO NLO samples
    weight *= ev.MCweight();

    //==== L1Prefire reweight
    weight *= weight_Prefire;

    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){

      double this_idsf = 1.;
      //double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      weight *= this_idsf*this_isosf;

    }

  }

  //==========================
  //==== Now fill histograms
  //==========================

  FillHist(param.Name+"/ZCand_Mass_"+param.Name, ZCand.M(), weight, 40, 70., 110.);
*/

}



