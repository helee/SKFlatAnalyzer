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
  MuonVetoIDs      = {"HNVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"HNVeto"};
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
  MuonTriggersTight.clear();
  ElectronTriggers.clear();
  EMuTriggers.clear();
  EMuTriggersTight.clear();

  EMuTriggersMu8.clear();
  EMuTriggersMu23.clear();

  if(DataEra == "2016preVFP"){                                                             // Lumi values of triggers in 2016 data (/pb)

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                          // 19517.523849710 
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                        // 19517.523849710
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                       // 19517.523849710
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                     // 19517.523849710
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");             // 19517.523849710
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");             // 19517.523849710
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");            // 19517.523849710

    EMuTriggersMu8.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");

    //Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");

    MuonPtCut1 = 20., MuonPtCut2 = 15.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }
  else if(DataEra == "2016postVFP"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                       // 8072.032418212  (FG)
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                     // 8072.032418212  (FG)
    MuonTriggersTight.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                  // 16812.151722311
    MuonTriggersTight.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                // 16812.151722311
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");             // 16812.151722311
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");             // 8072.032418212  (FG)
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");            // 8072.032418212  (FG)
    EMuTriggersTight.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 16812.151722311
    EMuTriggersTight.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");    // 16812.151722311

    EMuTriggersMu8.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggersMu8.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  }
  else if(DataEra == "2017"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                       // 4803.366325775  (B)
    MuonTriggersTight.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");          // 36674.511073518 (CDEF)
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");                // 41477.877399293
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");            // 36674.511073518 (CDEF)
    EMuTriggersTight.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 41477.877399293
    EMuTriggersTight.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");    // 41477.877399293

    EMuTriggersMu8.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 15.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }
  else if(DataEra == "2018"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");               // 59827.879505586
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");                // 59827.879505586
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");            // 59827.879505586
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 59827.879505586
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");         // 59827.879505586

    EMuTriggersMu8.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggersMu23.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 15.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }

  //cout << "[HNtypeI_DY_CR::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_DY_CR::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B Tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_DY_CR::~HNtypeI_DY_CR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_DY_CR::executeEvent(){

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
    param.Muon_Tight_ID           = MuonTightID;
    param.Muon_Loose_ID           = MuonLooseID;
    param.Muon_Veto_ID            = MuonVetoID;
    param.Muon_FR_ID              = MuonFRName;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key             = "AwayJetPt40";  // histname
    param.Muon_ID_SF_Key          = "";
    param.Muon_ISO_SF_Key         = "";
    param.Muon_Trigger_SF_Key     = "";
    param.Muon_UsePtCone          = true;

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

    //==== Systematics (JES, JER, L1Prefire, PU, Lepton ID/trigger SF, etc.)
    if(RunSyst){
      for(int it_syst=1; it_syst<23; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name  = "Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }

  }

}

void HNtypeI_DY_CR::executeEventFromParameter(AnalyzerParameter param){

  TString IDName = "HNTightV2";

  TString channel = "";
  vector<TString> regions = {"DY", "TT"};
  TString tight_leptons = "";  

  TString systName = param.Name;

  double cutflow_max = 12.;
  int cutflow_bin = 12;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trigger_weight = 0., diel_trigger_weight = 0., emu_trigger_weight = 0.;
  int tight_muons = 0, tight_electrons = 0;

  Event ev = GetEvent();

  //==== Boolean : primary datasets
  bool isDoubleMuon = false, isDoubleEG = false, isMuonEG = false;

  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) isDoubleMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) isDoubleEG = true;
    if(DataStream.Contains("MuonEG")) isMuonEG = true;
  }

  //==== Boolean : passing dilepton triggers
  //==== Run Numbers (DoubleMuon)
  //==== 2016postVBF : G (278820-280385), H (281613-284044)
  //==== 2017 : B (297047-299329), C (299368-302029)
  bool passMuMu = false, passEE = false, passEMu = false, passEMu8 = false, passEMu23 = false;

  if(DataEra == "2016postVFP"){

    if(IsDATA){

      if(run < 281000){
        passMuMu  = ev.PassTrigger(MuonTriggers) || ev.PassTrigger(MuonTriggersTight);
        passEMu   = ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersTight);
        passEMu8  = ev.PassTrigger(EMuTriggersMu8);
        passEMu23 = ev.PassTrigger(EMuTriggersMu23);
      }
      else{
        passMuMu  = ev.PassTrigger(MuonTriggers);
        passEMu   = ev.PassTrigger(EMuTriggers);
        passEMu8  = ev.PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        passEMu23 = ev.PassTrigger("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      }

    }
    else{
      passMuMu  = ev.PassTrigger(MuonTriggers) || ev.PassTrigger(MuonTriggersTight);
      passEMu   = ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersTight);
      passEMu8  = ev.PassTrigger(EMuTriggersMu8);
      passEMu23 = ev.PassTrigger(EMuTriggersMu23);
    }

  }
  else if(DataEra == "2017"){

    if(IsDATA){

      if(run < 299350){
        passMuMu  = ev.PassTrigger(MuonTriggers);
        passEMu   = ev.PassTrigger(EMuTriggersTight);
        passEMu8  = ev.PassTrigger(EMuTriggersMu8);
        passEMu23 = ev.PassTrigger("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      }
      else{
        passMuMu  = ev.PassTrigger(MuonTriggersTight);
        passEMu   = ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersTight);
        passEMu8  = ev.PassTrigger(EMuTriggersMu8);
        passEMu23 = ev.PassTrigger(EMuTriggersMu23);
      }

    }
    else{
      passMuMu  = ev.PassTrigger(MuonTriggers) || ev.PassTrigger(MuonTriggersTight);
      passEMu   = ev.PassTrigger(EMuTriggers) || ev.PassTrigger(EMuTriggersTight);
      passEMu8  = ev.PassTrigger(EMuTriggersMu8);
      passEMu23 = ev.PassTrigger(EMuTriggersMu23);
    }

  }
  else{
    passMuMu  = ev.PassTrigger(MuonTriggers);
    passEMu   = ev.PassTrigger(EMuTriggers);
    passEMu8  = ev.PassTrigger(EMuTriggersMu8);
    passEMu23 = ev.PassTrigger(EMuTriggersMu23);
  }

  passEE = ev.PassTrigger(ElectronTriggers);

  //==== Period-dependent trigger weights

  if(!IsDATA){

    if(DataEra == "2016postVFP"){
      if(ev.PassTrigger(MuonTriggers)) dimu_trigger_weight = 8072.032418212;
      if(ev.PassTrigger(MuonTriggersTight)) dimu_trigger_weight = ev.GetTriggerLumi("Full");
      if(ev.PassTrigger(EMuTriggers)) emu_trigger_weight = 8072.032418212;
      if(ev.PassTrigger(EMuTriggersTight)) emu_trigger_weight = ev.GetTriggerLumi("Full");
    }
    else if(DataEra == "2017"){
      if(ev.PassTrigger(MuonTriggers)) dimu_trigger_weight = 4803.366325775;
      if(ev.PassTrigger(MuonTriggersTight)) dimu_trigger_weight = ev.GetTriggerLumi("Full");
      if(ev.PassTrigger(EMuTriggers)) emu_trigger_weight = 36674.511073518;
      if(ev.PassTrigger(EMuTriggersTight)) emu_trigger_weight = ev.GetTriggerLumi("Full");
    }
    else{
      dimu_trigger_weight = ev.GetTriggerLumi("Full");
      emu_trigger_weight  = ev.GetTriggerLumi("Full");
    }

    diel_trigger_weight = ev.GetTriggerLumi("Full");

  }


  //========================================================
  //==== No Cut
  //========================================================

  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp, 0);
  }

  int Nvtx = nPV;
  FillHist(systName+"_Number_Vertices_NoCut", Nvtx, weight, 100, 0., 100.);

  //==== Cutflow 1
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //==== Cutflow 2
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 1.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_"+IDName, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 1.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== Trigger
  //========================================================

  if(!(passMuMu || passEE || passEMu)) return;

  //========================================================
  //==== Copy AllObjects
  //========================================================

  vector<Muon> this_AllMuons;
  if(param.Muon_Tight_ID.Contains("HighPt")) this_AllMuons = UseTunePMuon(AllMuons);
  else this_AllMuons = AllMuons;
  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Jet> this_AllJets = AllJets;
  vector<FatJet> this_AllFatJets = AllFatJets;
  vector<Gen> gens = GetGens();

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  string systBtag = "central";
  int systL1 = 0, systPU = 0, systMuonID = 0, systElectronReco = 0, systElectronID = 0, systMuonTrigger = 0, systElectronTrigger = 0, systTau21 = 0;

  if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(0), pfMET_Type1_phi_shifts->at(0));
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(1), pfMET_Type1_phi_shifts->at(1));
  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(2), pfMET_Type1_phi_shifts->at(2));
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(3), pfMET_Type1_phi_shifts->at(3));
  }
  else if(param.syst_ == AnalyzerParameter::UnclusteredEnUp){
    ev.SetMET(pfMET_Type1_pt_shifts->at(10), pfMET_Type1_phi_shifts->at(10));
  }
  else if(param.syst_ == AnalyzerParameter::UnclusteredEnDown){
    ev.SetMET(pfMET_Type1_pt_shifts->at(11), pfMET_Type1_phi_shifts->at(11));
  }
  /*else if(param.syst_ == AnalyzerParameter::BtagSFUp){
    systBtag = "up";
  }
  else if(param.syst_ == AnalyzerParameter::BtagSFDown){
    systBtag = "down";
  }*/
  else if(param.syst_ == AnalyzerParameter::L1PrefireUp){
    systL1 = 1;
  }
  else if(param.syst_ == AnalyzerParameter::L1PrefireDown){
    systL1 = -1;
  }
  else if(param.syst_ == AnalyzerParameter::PileupUp){
    systPU = 1;
  }
  else if(param.syst_ == AnalyzerParameter::PileupDown){
    systPU = -1;
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnUp){
    this_AllMuons = ScaleMuons( this_AllMuons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnDown){
    this_AllMuons = ScaleMuons( this_AllMuons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResUp){
    this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResDown){
    this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnUp){
    this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnDown){
    this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonIDSFUp){
    systMuonID = 1;
  }
  else if(param.syst_ == AnalyzerParameter::MuonIDSFDown){
    systMuonID = -1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronRecoSFUp){
    systElectronReco = 1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronRecoSFDown){
    systElectronReco = -1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronIDSFUp){
    systElectronID = 1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronIDSFDown){
    systElectronID = -1;
  }
  /*else if(param.syst_ == AnalyzerParameter::MuonTriggerSFUp){
    systMuonTrigger = 1;
  }
  else if(param.syst_ == AnalyzerParameter::MuonTriggerSFDown){
    systMuonTrigger = -1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronTriggerSFUp){
    systElectronTrigger = 1;
  }
  else if(param.syst_ == AnalyzerParameter::ElectronTriggerSFDown){
    systElectronTrigger = -1;
  }*/
  else{
    //cout << "[HNtypeI_DY_CR::executeEventFromParameter] Wrong syst" << endl;
    cerr << "[HNtypeI_DY_CR::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  //==== Leptons
  TString MuonID = param.Muon_Tight_ID;
  TString ElectronID = param.Electron_Tight_ID;
  if(RunFake){
    MuonID = param.Muon_Loose_ID;
    ElectronID = param.Electron_Loose_ID;
  }

  vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 5., 2.4);
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);
  vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 10., 2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);

  //==== Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  //==== Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  //==== Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  //==== Fatjet selection in CATanalyzer (see the links)
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  //vector<Jet> jets = JetsPassPileupMVA(jets_lepveto);
  vector<Jet> jets = JetsAwayFromFatJet(jets_lepveto, fatjets);  // AK4jets used in SR, CR

  vector<Jet> jets_Pt30;
  jets_Pt30.clear();

  for(unsigned int i=0; i<jets.size(); i++){
    if(jets.at(i).Pt() > 30.) jets_Pt30.push_back(jets.at(i));
  }

  std::vector<Lepton*> leptons, leptons_minus, leptons_plus, leptons_veto;

  //==================================================
  //==== Sort in pT-order
  //==================================================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);
  std::sort(jets_Pt30.begin(), jets_Pt30.end(), PtComparing);

  //========================================================
  //==== B tagging
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0, Nbjet_Pt30_loose = 0., Nbjet_Pt30_medium = 0, Nbjet_medium_lepveto = 0;
  JetTagging::Parameters jtp_DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepJet_Medium = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){

    if(jets_nolepveto.at(ij).Pt() > 20.){
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Loose, jets_nolepveto.at(ij), systBtag)) Nbjet_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_nolepveto.at(ij), systBtag)) Nbjet_medium++;
    }

    if(jets_nolepveto.at(ij).Pt() > 30.){
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Loose, jets_nolepveto.at(ij), systBtag)) Nbjet_Pt30_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_nolepveto.at(ij), systBtag)) Nbjet_Pt30_medium++;
    }

  }

  for(unsigned int ij=0; ij<jets.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets.at(ij), systBtag)) Nbjet_medium_lepveto++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  if((muons.size()+electrons.size() > 2) && (muons.size()+electrons.size() < 5)){
    METv = UpdateMETMuon(METv, muons);
    METv = UpdateMETElectron(METv, electrons);
  }

  double MET = METv.Pt();
  double METPhi = METv.Phi();

  //========================================================
  //==== Define particles, variables
  //========================================================

  double ST = 0., MET2ST = 0.;
  double Mt = 0., Mt3l = 0.;
  double dRll = 0., dPhill = 0., PtDiff = 0.;
  double MZ = 91.1876;
  double mllCut = 55.;  // 10 GeV cut in EXO-17-028
  double muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1., fatjetTau21SF = 1.;
  int lepton_veto_size = 0;
  double lepton1_eta = 0., lepton2_eta = 0.;

  bool passPtCut = false;
  Particle ZCand;

  //==== Set up pTcone if RunFake=true
  double tightIsoCut_muon = 0.07, tightIsoCut_electron = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(RunFake){

    if((muons.size()+electrons.size() > 2) && (muons.size()+electrons.size() < 5)){

      for(unsigned int i=0; i<muons.size(); i++){
        this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), tightIsoCut_muon);
        muons.at(i).SetPtCone(this_ptcone_muon);
      }

      for(unsigned int i=0; i<electrons.size(); i++){

        if(param.Electron_Tight_ID.Contains("HNTight")){ // POG cut-based tight WP
          tightIsoCut_electron = 0.0287+0.506/electrons.at(i).UncorrPt();
          if(fabs(electrons.at(i).scEta()) > 1.479) tightIsoCut_electron = 0.0445+0.963/electrons.at(i).UncorrPt();
        }

        this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), tightIsoCut_electron);
        electrons.at(i).SetPtCone(this_ptcone_electron);

      }

      //==== Correct MET if RunFake=true, because pT was replaced by pTcone
      METv = UpdateMETFake(METv, electrons, muons);

      muons = MuonUsePtCone(muons);
      electrons = ElectronUsePtCone(electrons);
      std::sort(muons.begin(), muons.end(), PtComparing);
      std::sort(electrons.begin(), electrons.end(), PtComparing);

    }

  }

  //==== Define leptons (pT order)
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));
  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  //==== Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  //==== Leptons (minus, plus charge)
  for(unsigned int i=0; i<muons.size(); i++){
    if(muons.at(i).Charge() < 0) leptons_minus.push_back(&muons.at(i));
    if(muons.at(i).Charge() > 0) leptons_plus.push_back(&muons.at(i));
  }

  for(unsigned int i=0; i<electrons.size(); i++){
    if(electrons.at(i).Charge() < 0) leptons_minus.push_back(&electrons.at(i));
    if(electrons.at(i).Charge() > 0) leptons_plus.push_back(&electrons.at(i));
  }

  lepton_veto_size = leptons_veto.size() - leptons.size();

  //==== Define ST, MET^2/ST
  MET = METv.Pt();
  METPhi = METv.Phi();

  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;

  //==== Number of tight leptons
  if(RunFake){

    if(leptons.size() == 2){

      tight_muons = 0, tight_electrons = 0;

      for(unsigned int i=0; i<muons.size(); i++){
        if(muons.at(i).PassID(param.Muon_Tight_ID)) tight_muons++;
      }
      for(unsigned int i=0; i<electrons.size(); i++){
        if(electrons.at(i).PassID(param.Electron_Tight_ID)) tight_electrons++;
      }

      if(tight_muons + tight_electrons == 0) tight_leptons = "Tight0";
      if(tight_muons + tight_electrons == 1) tight_leptons = "Tight1";
      if(tight_muons + tight_electrons == 2) tight_leptons = "Tight2";

    }

  }

  //========================================================
  //==== Event selection
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

    weight = 1., muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1.;

    if(!IsDATA){
      weight *= weight_norm_1invpb;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(systL1);
      weight *= GetPileUpWeight(nPileUp, systPU);
    }

    //==== Cutflow 3
    //==== Passing dilepton triggers
    if(passMuMu){
      if(!IsDATA || isDoubleMuon){
        if(!IsDATA) trigger_lumi = dimu_trigger_weight;
        FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }

    if(passEE){
      if(!IsDATA || isDoubleEG){
        if(!IsDATA) trigger_lumi = diel_trigger_weight;
        FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_"+IDName, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }

    if(passEMu){
       if(!IsDATA || isMuonEG){
        if(!IsDATA) trigger_lumi = emu_trigger_weight;
        FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_"+IDName, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }

  }

  //========================================================
  //==== DY, ttbar control region
  //========================================================

  //==== Two leptons

  if(leptons.size() == 2){

    //==== pT > trigger thresholds
    passPtCut = false;

    if(muons.size()==2 && electrons.size()==0){
      if(!passMuMu) return;
      if(!IsDATA) trigger_lumi = dimu_trigger_weight;
      if(IsDATA){ if(!isDoubleMuon) return; }
      if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      channel = "dimu";
    }
    
    if(muons.size()==0 && electrons.size()==2){
      if(!passEE) return;
      if(!IsDATA) trigger_lumi = diel_trigger_weight;
      if(IsDATA){ if(!isDoubleEG) return; }
      if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      channel = "diel";
    }

    if(muons.size()==1 && electrons.size()==1){
      if(!passEMu) return;
      if(!IsDATA) trigger_lumi = emu_trigger_weight;
      if(IsDATA){ if(!isMuonEG) return; }
      if(passEMu8){
        if(electrons.at(0).Pt()>EMuPtCut1 && muons.at(0).Pt()>EMuPtCut2) passPtCut = true;
      }
      if(passEMu23){
        if(muons.at(0).Pt()>EMuPtCut1 && electrons.at(0).Pt()>EMuPtCut2) passPtCut = true;;
      }
      channel = "emu";
    }

     if(!passPtCut) return;

    //==== Truth matching
    muons_prompt.clear();
    electrons_prompt.clear();
    muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
    electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

    if(channel == "dimu"){
      if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) return;
    }
    if(channel == "diel"){
      if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) return;
    }
    if(channel == "emu"){
      if(!(muons_prompt.size()==1 && electrons_prompt.size()==1)) return;
    }

    //==== Event weights for MC
    weight = 1.;

    if(!IsDATA){

      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(systL1);
      weight *= GetPileUpWeight(nPileUp, systPU);

      //==== Muons
      for(unsigned int i=0; i<muons.size(); i++){

        if(param.Muon_Tight_ID.Contains("HNTight")){

          muonRecoSF = 1.;
          muonIDSF   = mcCorr->MuonID_SF_HNtypeI(param.Muon_Tight_ID, muons.at(i).Eta(), muons.at(i).MiniAODPt(), systMuonID);
          muonIsoSF  = 1.;  // HNTight ID contains both ID and Iso. For POG ID muons, ID/Iso SFs are measured separately.

          if(RunFake){  // When subtracting prompt contribution from fake contribution, we apply ID SF only for muons passing the tight ID
            if(!muons.at(i).PassID(param.Muon_Tight_ID)) muonIDSF = 1.;
          }

        }
        else{
          muonRecoSF = 1.;
          muonIDSF   = 1.;
          muonIsoSF  = 1.;
        }

        weight *= muonRecoSF*muonIDSF*muonIsoSF;

      }

      //==== Electrons
      for(unsigned int i=0; i<electrons.size(); i++){

        electronRecoSF = mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), systElectronReco);

        if(param.Electron_Tight_ID.Contains("HNTight")){

          electronIDSF = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons.at(i).scEta(), electrons.at(i).UncorrPt(), systElectronID);

          if(RunFake){  // When subtracting prompt contribution from fake contribution, we apply ID SF only for electrons passing the tight ID
            if(!electrons.at(i).PassID(param.Electron_Tight_ID)) electronIDSF = 1.;
          }

        }
        else electronIDSF = 1.;

        weight *= electronRecoSF*electronIDSF;

      }

      //==== Triggers

      //==== AK8 jets
      /*for(unsigned int i=0; i<fatjets.size(); i++){

        fatjetTau21SF = mcCorr->FatJetWTagSF(param.FatJet_ID, systTau21);

        weight *= fatjetTau21SF;

      }*/

    }

    if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

    ZCand  = *leptons.at(0) + *leptons.at(1);
    dRll   = leptons.at(0)->DeltaR(*leptons.at(1));
    dPhill = fabs(leptons.at(0)->DeltaPhi(*leptons.at(1)));
    PtDiff = fabs(leptons.at(0)->Pt() - leptons.at(1)->Pt())/(leptons.at(0)->Pt() + leptons.at(1)->Pt());

    if(channel=="dimu"){
      lepton1_eta = muons.at(0).Eta();
      lepton2_eta = muons.at(1).Eta();
    }
    if(channel=="diel"){
      lepton1_eta = electrons.at(0).scEta();
      lepton2_eta = electrons.at(1).scEta();
    }
    if(channel=="emu"){
      lepton1_eta = muons.at(0).Eta();
      lepton2_eta = electrons.at(0).scEta();
    }

    //==== Cutflow 4
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 3.5, 1., cutflow_bin, 0., cutflow_max);
    }

    //==== OS events
    if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) return;

    //==== Cutflow 5
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 4.5, 1., cutflow_bin, 0., cutflow_max);
    }

    //==== No more leptons
    if(!(lepton_veto_size == 0)) return;

    //==== Cutflow 6
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 5.5, 1., cutflow_bin, 0., cutflow_max);
    }

    //==== m(ll) > 55 GeV (Not using DY10to50)
    if(!(ZCand.M() > mllCut)) return;

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

      //==== Cutflow 7
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 6.5, 1., cutflow_bin, 0., cutflow_max);

      if(it_rg == 0){ // DY CR

        //==== No b-tagged jets
        if(Nbjet_medium == 0){

          //==== Cutflow 8
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 7.5, 1., cutflow_bin, 0., cutflow_max);

          //==== Histograms
          if(systName == "Central"){
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_NoBJet_"+IDName, 0.5, weight, 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_NoBJet_"+IDName, 0.5, 1., 2, 0., 2.);
          }
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoBJet_"+IDName, Nvtx, weight, 200, 0., 200.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoBJet_"+IDName, jets.size(), weight, 10, 0., 10.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoBJet_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoBJet_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoBJet_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoBJet_"+IDName, dRll, weight, 60, 0., 6.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoBJet_"+IDName, dPhill, weight, 32, 0., 3.2);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoBJet_"+IDName, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoBJet_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoBJet_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoBJet_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoBJet_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_NoBJet_"+IDName, MET, weight, 2000, 0., 2000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoBJet_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          //==== Z peak region
          if(fabs(ZCand.M() - MZ) < 10.){

            //==== Cutflow 9
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);

            //==== Histograms
            if(systName == "Central"){
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
            }
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDName, Nvtx, weight, 200, 0., 200.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_"+IDName, MET, weight, 2000, 0., 2000.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          }

        }

      }
      else{ // TT CR

        //==== N(j) >= 2 && N(b) >= 1 
        if(jets.size()>=2 && Nbjet_medium_lepveto>=1){

          //==== Cutflow 8
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 7.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 7.5, 1., cutflow_bin, 0., cutflow_max);

          //==== Histograms
          if(systName == "Central"){
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_BJet_"+IDName, 0.5, weight, 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_BJet_"+IDName, 0.5, 1., 2, 0., 2.);
          }
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_BJet_"+IDName, Nvtx, weight, 200, 0., 200.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_BJet_"+IDName, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_BJet_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_BJet_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_BJet_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_BJet_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_BJet_"+IDName, dRll, weight, 60, 0., 6.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_BJet_"+IDName, dPhill, weight, 32, 0., 3.2);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_BJet_"+IDName, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_BJet_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_BJet_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_BJet_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_BJet_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_BJet_"+IDName, MET, weight, 2000, 0., 2000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_BJet_"+IDName, MET2ST, weight, 2000, 0., 2000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_BJet_"+IDName, fatjets.at(0).Pt(), weight, 2000, 0., 2000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_BJet_"+IDName, fatjets.at(0).SDMass(), weight, 2000, 0., 2000.);

          //==== MET cut
          if(MET > 40.){

            //==== Cutflow 9
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);

            //==== Histograms
            if(systName == "Central"){
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
            }
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDName, Nvtx, weight, 200, 0., 200.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_"+IDName, MET, weight, 2000, 0., 2000.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_"+IDName, fatjets.at(0).Pt(), weight, 2000, 0., 2000.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_"+IDName, fatjets.at(0).SDMass(), weight, 2000, 0., 2000.);

          }

        }

      }

    }

  }

}



