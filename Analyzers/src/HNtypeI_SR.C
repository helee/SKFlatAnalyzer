#include "HNtypeI_SR.h"

HNtypeI_SR::HNtypeI_SR(){

}

void HNtypeI_SR::initializeAnalyzer(){

  //==== If you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  RunFake = HasFlag("RunFake");
  RunCF   = HasFlag("RunCF");
  RunOS   = HasFlag("RunOS");
  RunAK8  = HasFlag("RunAK8");

  cout << "[HNtypeI_SR::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_SR::initializeAnalyzer] RunFake = " << RunFake << endl;
  cout << "[HNtypeI_SR::initializeAnalyzer] RunCF = " << RunCF << endl;
  cout << "[HNtypeI_SR::initializeAnalyzer] RunOS = " << RunOS << endl;
  cout << "[HNtypeI_SR::initializeAnalyzer] RunAK8 = " << RunAK8 << endl;

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

    MuonPtCut1 = 20., MuonPtCut2 = 15.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;

  }
  else if(DataEra == "2016postVFP"){

    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                          // 8072.032418212  (FG)
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                        // 8072.032418212  (FG)
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

  //==== B tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_SR::~HNtypeI_SR(){

  //==== Destructor of this Analyzer

}

void HNtypeI_SR::executeEvent(){

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_SR.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SR.h
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
    if(DataEra.Contains("2016")) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45";

    executeEventFromParameter(param);

    //==== Systematics (JES, JER, L1Prefire, PU, Lepton ID/trigger SF, etc.)
    if(RunSyst){

      /*for(int it_syst=1; it_syst<23; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name  = "Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }*/

    }   

  }

}

void HNtypeI_SR::executeEventFromParameter(AnalyzerParameter param){

  TString IDName = "HNTightV2";

  TString channel = "";
  vector<TString> regions = {"SR1", "CR1", "SR2", "CR2", "CR2b", "SR3", "CR3"};
  vector<double> Lep1PtCutSR1, Lep1PtCutSR3, Lep2PtCutSR1, Lep2PtCutSR3, mlljjCut, mljjCut1, mljjCut2, mlJCut1, mlJCut2, MET2STCut;
  vector<int> mass = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1500, 2000, 2500};
  TString tight_leptons = "";

  TString systName = param.Name;

  double cutflow_max = 17.;
  int cutflow_bin = 17;
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
    //cout << "[HNtypeI_SR::executeEventFromParameter] Wrong syst" << endl;
    cerr << "[HNtypeI_SR::executeEventFromParameter] Wrong syst" << endl;
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

  //==== For charge flip
  vector<Electron> electrons_beforeShift;
  vector<Electron> electrons_afterShift;
  electrons_beforeShift.clear();
  electrons_afterShift.clear();

  //==== Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 4.7);
  vector<Jet> jets_bcand = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 4.7);

  //==== Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  //==== Fatjet selection in CATanalyzer (see the links)
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  vector<Jet> jets_pv = JetsPassPileupMVA(jets_lepveto, "loose");
  vector<Jet> jets = JetsAwayFromFatJet(jets_pv, fatjets);  // AK4jets used in SR, CR

  vector<Jet> jets_eta2p7;
  jets_eta2p7.clear();

  for(unsigned int i=0; i<jets.size(); i++){
    if(fabs(jets.at(i).Eta()) < 2.7) jets_eta2p7.push_back(jets.at(i));
  }

  std::vector<Lepton*> leptons, leptons_minus, leptons_plus, leptons_veto;

  vector<Jet> jets_WCand;
  FatJet fatjets_WCand;
  jets_WCand.clear();

  //==================================================
  //==== Sort in pT-order
  //==================================================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_bcand.begin(), jets_bcand.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);
  std::sort(jets_eta2p7.begin(), jets_eta2p7.end(), PtComparing);

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
  for(unsigned int ij=0; ij<jets_bcand.size(); ij++){

    if(jets_bcand.at(ij).Pt() > 20.){
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Loose, jets_bcand.at(ij), systBtag)) Nbjet_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_bcand.at(ij), systBtag)) Nbjet_medium++;
    }

    /*if(jets_bcand.at(ij).Pt() > 30.){
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Loose, jets_bcand.at(ij), systBtag)) Nbjet_Pt30_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_bcand.at(ij), systBtag)) Nbjet_Pt30_medium++;
    }*/

  }

  for(unsigned int ij=0; ij<jets.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets.at(ij), systBtag)) Nbjet_medium_lepveto++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  if(muons.size()+electrons.size() == 2){
    METv = UpdateMETMuon(METv, muons);
    METv = UpdateMETElectron(METv, electrons);
    METv = UpdateMETSmearedJet(METv, jets);
  }

  double MET = METv.Pt();
  double METPhi = METv.Phi();

  

  //========================================================
  //==== Define particles, variables
  //========================================================
 
  double ST = 0., HT = 0., MET2ST = 0., HTPt1 = 0.;
  //double Mt = 0., Mt3l = 0.;
  double dRll = 0., dPhill = 0., PtDiff = 0.;
  double dRl1jj = 0., dRl2jj = 0., dRlCjj = 0., dRlAjj = 0., dRjj = 0., dRl1J = 0., dRl2J = 0., dRlCJ = 0., dRlAJ = 0.;
  double avgEta = 0., dEtajj = 0., zep = 0.;
  double MZ = 91.1876, MW = 80.379;
  double mllCut1 = 10.;   // 10 GeV cut in EXO-17-028
  double mllCut2 = 55.;   // For OS events with AK8 jets
  double muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1., fatjetTau21SF = 1.;
  int lepton_veto_size = 0;
  double lepton1_eta = 0., lepton2_eta = 0.;

  bool passPtCut = false;
  Particle ZCand, WCand, jjVBF;
  Particle lCloseSR1, lAwaySR1, lCJ, lAJ, l1J, l2J, llJ;
  Particle lCloseSR3, lAwaySR3, lCjj, lAjj, l1jj, l2jj, lljj;

  //==== Set up pTcone when RunFake is true
  double tightIsoCut_muon = 0.07, tightIsoCut_electron = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(RunFake){

    if(muons.size()+electrons.size() == 2){

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

      //==== Correct MET when RunFake is true, because pT was replaced by pTcone
      METv = UpdateMETFake(METv, electrons, muons);

      muons = MuonUsePtCone(muons);
      electrons = ElectronUsePtCone(electrons);
      std::sort(muons.begin(), muons.end(), PtComparing);
      std::sort(electrons.begin(), electrons.end(), PtComparing);

    }

  }

  //==== Shift electron energy and MET wuen RunCF is true
  if(RunCF){

    if(muons.size()==0 && electrons.size()==2){

      electrons_beforeShift.push_back(electrons.at(0));
      electrons_beforeShift.push_back(electrons.at(1));
      electrons = ShiftElectronEnergy(param.Electron_Tight_ID, electrons, true);
      electrons_afterShift.push_back(electrons.at(0));
      electrons_afterShift.push_back(electrons.at(1));
      METv = UpdateMETElectronCF(METv, electrons_beforeShift, electrons_afterShift);

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

  for(unsigned int i=0; i<jets.size(); i++){
    ST += jets.at(i).Pt();
    HT += jets.at(i).Pt();
  }
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
      for(unsigned int i=0; i<fatjets.size(); i++){

        fatjetTau21SF = mcCorr->FatJetWTagSF(param.FatJet_ID, systTau21);

        weight *= fatjetTau21SF;

      }

    }

    if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

    if(RunCF) weight *= GetCFWeight(param.Electron_Tight_ID, leptons, true, 0);

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
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 3.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 3.5, 1., cutflow_bin, 0., cutflow_max);

    //==== Same-sign events (opposite-sign when RunCF is true)
    if(RunCF || RunOS || RunAK8){ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) return; }
    else{ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) return; }

    //==== Cutflow 5
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 4.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 4.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 4.5, 1., cutflow_bin, 0., cutflow_max);

    //==== No more leptons
    if(!(lepton_veto_size == 0)) return;

    //==== Cutflow 6
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 5.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 5.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 5.5, 1., cutflow_bin, 0., cutflow_max);

    //==== m(ll) cut
    if(!(ZCand.M() > mllCut1)) return;
    if(RunOS || RunAK8){
      if(!(ZCand.M() > mllCut2)) return;
    }
    else{
      if(channel=="diel" && IsOnZ(ZCand.M(), 10.)) return;
    }

    //==== Cutflow 7
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 6.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 6.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 6.5, 1., cutflow_bin, 0., cutflow_max);

    //==== Non-prompt CRs in EXO-17-028

    //==== CR1
    if(Nbjet_medium > 0){

      FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 7.5, 1., cutflow_bin, 0., cutflow_max);

      FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
      FillHist(systName+"_"+channel+"_FakeCR1_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
      FillHist(systName+"_"+channel+"_FakeCR1_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
      FillHist(systName+"_"+channel+"_FakeCR1_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
      FillHist(systName+"_"+channel+"_FakeCR1_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
      //FillHist(systName+"_"+channel+"_FakeCR1_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
      FillHist(systName+"_"+channel+"_FakeCR1_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
      FillHist(systName+"_"+channel+"_FakeCR1_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
      FillHist(systName+"_"+channel+"_FakeCR1_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
      FillHist(systName+"_"+channel+"_FakeCR1_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
      //FillHist(systName+"_"+channel+"_FakeCR1_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
      FillHist(systName+"_"+channel+"_FakeCR1_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"_"+channel+"_FakeCR1_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"_"+channel+"_FakeCR1_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
      FillHist(systName+"_"+channel+"_FakeCR1_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
      FillHist(systName+"_"+channel+"_FakeCR1_MET_"+IDName, MET, weight, 2000, 0., 2000.);
      FillHist(systName+"_"+channel+"_FakeCR1_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

    }

    //==== CR2
    if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){

      FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 7.5, 1., cutflow_bin, 0., cutflow_max);

      if(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5){

        FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);

        FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
        FillHist(systName+"_"+channel+"_FakeCR2_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
        FillHist(systName+"_"+channel+"_FakeCR2_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
        //FillHist(systName+"_"+channel+"_FakeCR2_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
        //FillHist(systName+"_"+channel+"_FakeCR2_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
        //FillHist(systName+"_"+channel+"_FakeCR2_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"_"+channel+"_FakeCR2_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_FakeCR2_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_FakeCR2_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
        FillHist(systName+"_"+channel+"_FakeCR2_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
        //FillHist(systName+"_"+channel+"_FakeCR2_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"_"+channel+"_FakeCR2_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_FakeCR2_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_FakeCR2_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"_"+channel+"_FakeCR2_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"_"+channel+"_FakeCR2_MET_"+IDName, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_FakeCR2_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

      }

    }

    //==== Jet requirement
    if(RunAK8){
      if(!(fatjets.size() >= 1)) return;
    }
    else{
      if(!(jets.size()>=2 || fatjets.size()>=1)) return;
    }

    //==== Preselection
    FillHist(systName+"_"+channel+"_Pre_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
    FillHist(systName+"_"+channel+"_Pre_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
    FillHist(systName+"_"+channel+"_Pre_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
    FillHist(systName+"_"+channel+"_Pre_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
    FillHist(systName+"_"+channel+"_Pre_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
    FillHist(systName+"_"+channel+"_Pre_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
    FillHist(systName+"_"+channel+"_Pre_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
    FillHist(systName+"_"+channel+"_Pre_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
    FillHist(systName+"_"+channel+"_Pre_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
    FillHist(systName+"_"+channel+"_Pre_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
    FillHist(systName+"_"+channel+"_Pre_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
    FillHist(systName+"_"+channel+"_Pre_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
    FillHist(systName+"_"+channel+"_Pre_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
    FillHist(systName+"_"+channel+"_Pre_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
    FillHist(systName+"_"+channel+"_Pre_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
    FillHist(systName+"_"+channel+"_Pre_MET_"+IDName, MET, weight, 2000, 0., 2000.);
    FillHist(systName+"_"+channel+"_Pre_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);


    //==== Event selection for SRs and CRs
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

      jets_WCand.clear();

      //==== Cutflow 8
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 7.5, 1., cutflow_bin, 0., cutflow_max);

      //==== SR1 (with AK8 jets)
      if(it_rg < 2){

        if(!(fatjets.size() >= 1)) continue;

        fatjets_WCand = FatJetWCand(fatjets, MW);

        l1J   = *leptons.at(0) + fatjets_WCand;
        l2J   = *leptons.at(1) + fatjets_WCand;
        llJ   = *leptons.at(0) + *leptons.at(1) + fatjets_WCand;
        dRl1J = leptons.at(0)->DeltaR(fatjets_WCand);
        dRl2J = leptons.at(1)->DeltaR(fatjets_WCand);

        if(dRl1J < dRl2J){
          lCloseSR1 = *leptons.at(0);
          lAwaySR1  = *leptons.at(1);
          dRlCJ = dRl1J;
          dRlAJ = dRl2J;
        }
        else{
          lCloseSR1 = *leptons.at(1);
          lAwaySR1  = *leptons.at(0);
          dRlCJ = dRl2J;
          dRlAJ = dRl1J;
        }

        lCJ = lCloseSR1 + fatjets_WCand;
        lAJ = lAwaySR1 + fatjets_WCand;

        //==== Cutflow 9
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);

        //==== Histograms

        if(it_rg==0 && systName=="Central"){

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nocut_"+IDName, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nocut_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nocut_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nocut_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nocut_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nocut_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nocut_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_nocut_"+IDName, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_nocut_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_nocut_"+IDName, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);

          if(Nbjet_medium == 0){

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nobjet_"+IDName, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nobjet_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nobjet_"+IDName, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nobjet_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nobjet_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nobjet_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nobjet_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_nobjet_"+IDName, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_nobjet_"+IDName, MET2ST, weight, 2000, 0., 2000.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_nobjet_"+IDName, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_nobjet_"+IDName, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_nobjet_"+IDName, llJ.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_nobjet_"+IDName, l1J.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_nobjet_"+IDName, l2J.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_nobjet_"+IDName, dRl1J, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_nobjet_"+IDName, dRl2J, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_Mass_nobjet_"+IDName, lCJ.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_Mass_nobjet_"+IDName, lAJ.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_DeltaR_nobjet_"+IDName, dRlCJ, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_DeltaR_nobjet_"+IDName, dRlAJ, weight, 60, 0., 6.);

            if(RunAK8 && MET2ST<15.){

              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_MET2ST15_"+IDName, 0.5, weight, 1, 0., 1.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_MET2ST15_"+IDName, 0.5, 1., 1, 0., 1.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_MET2ST15_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);

            }

          }

        }

        //==== OS SR
        if(RunAK8){
          if(!(ZCand.M() > 110.)) continue;
        }

        //==== SR1
        if(it_rg == 0){
          if(!(Nbjet_medium == 0)) continue;
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(MET2ST < 15.)) continue;
        }

        //==== CR1
        if(it_rg == 1){
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
        }

        //==== Cutflow 10
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 9.5, 1., cutflow_bin, 0., cutflow_max);

        //==== Histograms
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
        //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDName, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDName, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_"+IDName, MET, weight, 2000, 0., 2000.);
        //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDName, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_"+IDName, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_"+IDName, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_"+IDName, llJ.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_"+IDName, l1J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_"+IDName, l2J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_"+IDName, dRl1J, weight, 60, 0., 6.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_"+IDName, dRl2J, weight, 60, 0., 6.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_Mass_"+IDName, lCJ.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_Mass_"+IDName, lAJ.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_DeltaR_"+IDName, dRlCJ, weight, 60, 0., 6.);
        FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_DeltaR_"+IDName, dRlAJ, weight, 60, 0., 6.);

        if(RunOS){

          if(IsOnZ(ZCand.M(), 10.)) continue;

          //==== Cutflow 11
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 10.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 10.5, 1., cutflow_bin, 0., cutflow_max);

          //==== Histograms
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_NoMZ_"+IDName, 0.5, weight, 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_NoMZ_"+IDName, 0.5, 1., 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoMZ_"+IDName, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoMZ_"+IDName, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_NoMZ_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoMZ_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoMZ_"+IDName, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoMZ_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoMZ_"+IDName, dRll, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoMZ_"+IDName, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoMZ_"+IDName, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoMZ_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoMZ_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoMZ_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoMZ_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_NoMZ_"+IDName, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoMZ_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_NoMZ_"+IDName, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_NoMZ_"+IDName, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_NoMZ_"+IDName, llJ.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_NoMZ_"+IDName, l1J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_NoMZ_"+IDName, l2J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_NoMZ_"+IDName, dRl1J, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_NoMZ_"+IDName, dRl2J, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_Mass_NoMZ_"+IDName, lCJ.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_Mass_NoMZ_"+IDName, lAJ.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCJ_DeltaR_NoMZ_"+IDName, dRlCJ, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAJ_DeltaR_NoMZ_"+IDName, dRlAJ, weight, 60, 0., 6.);

        }

        //==== Optimized cuts in EXO-17-028
        if(!(it_rg == 0)) continue;

        Lep1PtCutSR1.clear();
        Lep2PtCutSR1.clear();
        mlJCut1.clear();
        mlJCut2.clear();

        if(!RunAK8){

          if(channel=="dimu"){  // For these cuts, see page 22-24 of arXiv:1806.10905
            Lep1PtCutSR1 = {25.,  100., 140., 140., 140., 140., 140., 140., 140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.};
            Lep2PtCutSR1 = {15.,  20.,  40.,  65.,  65.,  10.,  10.,  10.,  10.,   10.,   10.,   10.,   10.,   10.,   10.,   10.};
            mlJCut1      = {98.,  175., 280., 340., 445., 560., 635., 755., 840.,  900.,  990.,  1035., 1100., 1330., 1700., 2200.};
            mlJCut2      = {145., 235., 340., 445., 560., 685., 825., 960., 1055., 1205., 1250., 1430., 1595., 1800., 2300., 2800.};
          }
          if(channel=="diel"){
            Lep1PtCutSR1 = {25.,  100., 100., 100., 120., 120., 140., 140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.};
            Lep2PtCutSR1 = {15.,  20.,  30.,  35.,  35.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
            mlJCut1      = {100., 173., 270., 330., 440., 565., 635., 740.,  865.,  890.,  1035., 1085., 1140., 1300., 1700., 2200.};
            mlJCut2      = {220., 220., 330., 440., 565., 675., 775., 1005., 1030., 1185., 1395., 1460., 1590., 1800., 2300., 2800.};
          }
          if(channel=="emu"){
            Lep1PtCutSR1 = {30.,  70.,  95.,  125., 145., 160., 170., 170., 180.,  180.,  180.,  180.,  180.,  180.,  180.,  180.};
            Lep2PtCutSR1 = {15.,  30.,  55.,  55.,  60.,  15.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
            mlJCut1      = {100., 180., 280., 340., 460., 555., 610., 730., 845.,  930.,  1020., 1080., 1155., 1345., 1800., 2300.};
            mlJCut2      = {335., 225., 340., 475., 555., 645., 780., 895., 1015., 1075., 1340., 1340., 1595., 1615., 2200., 2700.};
          }

          for(unsigned int it_m=0; it_m<mass.size(); it_m++){

            if(!(leptons.at(0)->Pt() > Lep1PtCutSR1.at(it_m))) continue;
            if(!(leptons.at(1)->Pt() > Lep2PtCutSR1.at(it_m))) continue;

            if(it_m < 2){
              if(!(l2J.M()>mlJCut1.at(it_m) && l2J.M()<mlJCut2.at(it_m))) continue;
            }
            else{
              if(!(l1J.M()>mlJCut1.at(it_m) && l1J.M()<mlJCut2.at(it_m)) && !(l2J.M()>mlJCut1.at(it_m) && l2J.M()<mlJCut2.at(it_m))) continue;
            }

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);

            /*if(channel=="diel"){
              if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
                FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_"+IDName, 0.5, weight, 1, 0., 1.);
                FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_unweighted_"+IDName, 0.5, 1., 1, 0., 1.);
              }
            }*/

          }

        }

      }

      //==== SR2, SR3 (without AK8 jets)
      if(it_rg >= 2){

        if(RunAK8) continue;
        if(!(jets.size()>=2 && fatjets.size()==0)) continue;

        jjVBF  = jets.at(0) + jets.at(1);
        avgEta = 0.5*(jets.at(0).Eta() + jets.at(1).Eta());
        dEtajj = fabs(jets.at(0).Eta() - jets.at(1).Eta());
        zep    = std::max(fabs(lepton1_eta - avgEta), fabs(lepton2_eta - avgEta))/dEtajj;
        HTPt1  = HT/leptons.at(0)->Pt();

        //==== SR2, CR2
        //==== See https://github.com/Hooooon12/SKFlatAnalyzer/blob/UL_HNtype1/Analyzers/src/SSWW.C
        if(it_rg < 5){

          //==== Cutflow 9
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(leptons.at(0)->Pt()>30. && leptons.at(1)->Pt()>30.)) continue;

          //==== Cutflow 10
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 9.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 9.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(ZCand.M() > 20.)) continue;

          //==== Cutflow 11
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 10.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 10.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(jets.at(1).Pt() > 30.)) continue;

          //==== Cutflow 12
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 11.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 11.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(jjVBF.M() > 750.)) continue;

          //==== Cutflow 13
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 12.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 12.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(dEtajj > 2.5)) continue;

          //==== Cutflow 14
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 13.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 13.5, 1., cutflow_bin, 0., cutflow_max);

          if(!(zep < 0.75)) continue;

          //==== Cutflow 15
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 14.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 14.5, 1., cutflow_bin, 0., cutflow_max);

          //==== SR2
          if(it_rg == 2){
            if(!(Nbjet_medium == 0)) continue;
            if(!(dPhill > 2.)) continue;
          }

          //==== CR2
          if(it_rg == 3){
            if(!(Nbjet_medium == 0)) continue;
            if(!(dPhill <= 2.)) continue;
          }

          //==== CR2b
          if(it_rg == 4){
            if(!(Nbjet_medium > 0)) continue;
          }

          //==== Cutflow 16
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 15.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 15.5, 1., cutflow_bin, 0., cutflow_max);

          //==== Histograms
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_"+IDName, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet1_Pt_"+IDName, jets.at(0).Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet2_Pt_"+IDName, jets.at(1).Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet1_Eta_"+IDName, jets.at(0).Eta(), weight, 100, -5., 5.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet2_Eta_"+IDName, jets.at(1).Eta(), weight, 100, -5., 5.);
          
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Zeppenfeld_"+IDName, zep, weight, 150, 0., 1.5);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_jjVBF_Mass_"+IDName, jjVBF.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_jjVBF_DeltaEta_"+IDName, dEtajj, weight, 100, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_HTPt1_"+IDName, HTPt1, weight, 10, 0., 10.);

          if(RunOS){

            if(IsOnZ(ZCand.M(), 10.)) continue;

            //==== Cutflow 17
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 16.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 16.5, 1., cutflow_bin, 0., cutflow_max);

            //==== Histograms
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_NoMZ_"+IDName, 0.5, weight, 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_NoMZ_"+IDName, 0.5, 1., 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoMZ_"+IDName, Nvtx, weight, 100, 0., 100.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoMZ_"+IDName, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_NoMZ_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoMZ_"+IDName, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoMZ_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoMZ_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoMZ_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoMZ_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_NoMZ_"+IDName, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoMZ_"+IDName, MET2ST, weight, 2000, 0., 2000.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoMZ_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoMZ_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoMZ_"+IDName, dRll, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoMZ_"+IDName, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoMZ_"+IDName, PtDiff, weight, 100, 0., 1.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet1_Pt_NoMZ_"+IDName, jets.at(0).Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet2_Pt_NoMZ_"+IDName, jets.at(1).Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet1_Eta_NoMZ_"+IDName, jets.at(0).Eta(), weight, 100, -5., 5.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Jet2_Eta_NoMZ_"+IDName, jets.at(1).Eta(), weight, 100, -5., 5.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Zeppenfeld_NoMZ_"+IDName, zep, weight, 150, 0., 1.5);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_jjVBF_Mass_NoMZ_"+IDName, jjVBF.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_jjVBF_DeltaEta_NoMZ_"+IDName, dEtajj, weight, 100, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_HTPt1_NoMZ_"+IDName, HTPt1, weight, 10, 0., 10.);

          }

        }

        //==== SR3, CR3
        if(it_rg >= 5){

          //==== SSWW veto
          if(leptons.at(1)->Pt()>30. && ZCand.M()>20. && jets.at(1).Pt()>30. && jjVBF.M()>750. && dEtajj>2.5 && zep<0.75) continue;

          jets_WCand = JetsWCandHighMass(jets, MW);
          WCand      = jets_WCand.at(0) + jets_WCand.at(1);
          lljj       = *leptons.at(0) + *leptons.at(1) + jets_WCand.at(0) + jets_WCand.at(1);
          l1jj       = *leptons.at(0) + jets_WCand.at(0) + jets_WCand.at(1);
          l2jj       = *leptons.at(1) + jets_WCand.at(0) + jets_WCand.at(1);
          dRl1jj     = leptons.at(0)->DeltaR(WCand);
          dRl2jj     = leptons.at(1)->DeltaR(WCand);
          dRjj       = jets_WCand.at(0).DeltaR(jets_WCand.at(1));

          if(dRl1jj < dRl2jj){
            lCloseSR3 = *leptons.at(0);
            lAwaySR3  = *leptons.at(1);
            dRlCjj = dRl1jj;
            dRlAjj = dRl2jj;
          }
          else{
            lCloseSR3 = *leptons.at(1);
            lAwaySR3  = *leptons.at(0);
            dRlCjj = dRl2jj;
            dRlAjj = dRl1jj;
          }

          lCjj = lCloseSR3 + jets_WCand.at(0) + jets_WCand.at(1);
          lAjj = lAwaySR3 + jets_WCand.at(0) + jets_WCand.at(1);

          //==== Cutflow 9
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 8.5, 1., cutflow_bin, 0., cutflow_max);
 
          //==== Histograms

          if(it_rg==5 && systName=="Central"){

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nocut_"+IDName, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nocut_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nocut_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nocut_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nocut_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nocut_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_nocut_"+IDName, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_nocut_"+IDName, MET2ST, weight, 2000, 0., 2000.);

            if(Nbjet_medium == 0){

              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_nobjet_"+IDName, 0.5, weight, 2, 0., 2.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_nobjet_"+IDName, 0.5, 1., 2, 0., 2.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nobjet_"+IDName, jets.size(), weight, 10, 0., 10.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nobjet_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nobjet_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nobjet_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nobjet_"+IDName, lepton1_eta, weight, 60, -3., 3.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nobjet_"+IDName, lepton2_eta, weight, 60, -3., 3.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_nobjet_"+IDName, MET, weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_nobjet_"+IDName, MET2ST, weight, 2000, 0., 2000.);

              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nobjet_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nobjet_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nobjet_"+IDName, dRll, weight, 60, 0., 6.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nobjet_"+IDName, dPhill, weight, 32, 0., 3.2);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nobjet_"+IDName, PtDiff, weight, 100, 0., 1.);

              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_nobjet_"+IDName, WCand.M(), weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_nobjet_"+IDName, WCand.Pt(), weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_nobjet_"+IDName, dRjj, weight, 60, 0., 6.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_nobjet_"+IDName, jets_WCand.at(0).Pt(), weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_nobjet_"+IDName, jets_WCand.at(1).Pt(), weight, 2000, 0., 2000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lljj_Mass_nobjet_"+IDName, lljj.M(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_Mass_nobjet_"+IDName, l1jj.M(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_Mass_nobjet_"+IDName, l2jj.M(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_DeltaR_nobjet_"+IDName, dRl1jj, weight, 60, 0., 6.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_DeltaR_nobjet_"+IDName, dRl2jj, weight, 60, 0., 6.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_Mass_nobjet_"+IDName, lCjj.M(), weight, 4000, 0., 4000.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_Mass_nobjet_"+IDName, lAjj.M(), weight, 4000, 0., 4000.);      
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_DeltaR_nobjet_"+IDName, dRlCjj, weight, 60, 0., 6.);
              FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_DeltaR_nobjet_"+IDName, dRlAjj, weight, 60, 0., 6.);

            }

          }

          //==== SR3
          if(it_rg == 5){
            if(!(Nbjet_medium == 0)) continue;
            if(!(WCand.M()>30. && WCand.M()<150.)) continue;
            if(!(MET2ST < 15.)) continue;
          }

          //==== CR3
          if(it_rg == 6){
            if(!(WCand.M()>30. && WCand.M()<150.)) continue;
            if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
          }

          //==== Cutflow 10
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 9.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 9.5, 1., cutflow_bin, 0., cutflow_max);

          //==== Histograms
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDName, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDName, jets.size(), weight, 10, 0., 10.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDName, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDName, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDName, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDName, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_"+IDName, MET, weight, 2000, 0., 2000.);
          //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDName, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDName, MET2ST, weight, 2000, 0., 2000.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDName, dRll, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDName, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDName, PtDiff, weight, 100, 0., 1.);

          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_"+IDName, WCand.M(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_"+IDName, WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_"+IDName, dRjj, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_"+IDName, jets_WCand.at(0).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_"+IDName, jets_WCand.at(1).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lljj_Mass_"+IDName, lljj.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_Mass_"+IDName, l1jj.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_Mass_"+IDName, l2jj.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_DeltaR_"+IDName, dRl1jj, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_DeltaR_"+IDName, dRl2jj, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_Mass_"+IDName, lCjj.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_Mass_"+IDName, lAjj.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_DeltaR_"+IDName, dRlCjj, weight, 60, 0., 6.);
          FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_DeltaR_"+IDName, dRlAjj, weight, 60, 0., 6.);

          if(RunOS){

            if(IsOnZ(ZCand.M(), 10.)) continue;

            //==== Cutflow 11
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDName, 10.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDName, 10.5, 1., cutflow_bin, 0., cutflow_max);

            //==== Histograms
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_NoMZ_"+IDName, 0.5, weight, 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_NoMZ_"+IDName, 0.5, 1., 2, 0., 2.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoMZ_"+IDName, Nvtx, weight, 100, 0., 100.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoMZ_"+IDName, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_NoMZ_"+IDName, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoMZ_"+IDName, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoMZ_"+IDName, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoMZ_"+IDName, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoMZ_"+IDName, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoMZ_"+IDName, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET_NoMZ_"+IDName, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoMZ_"+IDName, MET2ST, weight, 2000, 0., 2000.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoMZ_"+IDName, ZCand.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoMZ_"+IDName, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoMZ_"+IDName, dRll, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoMZ_"+IDName, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoMZ_"+IDName, PtDiff, weight, 100, 0., 1.);

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_NoMZ_"+IDName, WCand.M(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_NoMZ_"+IDName, WCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_NoMZ_"+IDName, dRjj, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_NoMZ_"+IDName, jets_WCand.at(0).Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_NoMZ_"+IDName, jets_WCand.at(1).Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lljj_Mass_NoMZ_"+IDName, lljj.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_Mass_NoMZ_"+IDName, l1jj.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_Mass_NoMZ_"+IDName, l2jj.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l1jj_DeltaR_NoMZ_"+IDName, dRl1jj, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_l2jj_DeltaR_NoMZ_"+IDName, dRl2jj, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_Mass_NoMZ_"+IDName, lCjj.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_Mass_NoMZ_"+IDName, lAjj.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lCjj_DeltaR_NoMZ_"+IDName, dRlCjj, weight, 60, 0., 6.);
            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_lAjj_DeltaR_NoMZ_"+IDName, dRlAjj, weight, 60, 0., 6.);

          }
 
          //==== Optimized cuts in EXO-17-028
          if(!(it_rg == 5)) continue;

          Lep1PtCutSR3.clear();
          Lep2PtCutSR3.clear();
          mlljjCut.clear();
          mljjCut1.clear();
          mljjCut2.clear();
          MET2STCut.clear();

          if(channel=="dimu"){  // For these cuts, see page 22-24 of arXiv:1806.10905
            Lep1PtCutSR3 = {25.,  50.,  100., 110., 110., 110., 110., 110., 110.,  110.,  110.,  110.,  110.,  110.,  110.,  110.};
            Lep2PtCutSR3 = {15.,  40.,  50.,  60.,  60.,  10.,  10.,  10.,  10.,   10.,   10.,   10.,   10.,   10.,   10.,   10.};
            mlljjCut     = {110., 250., 370., 490., 610., 680., 800., 800., 800.,  800.,  800.,  800.,  800.,  800.,  800.,  800.};
            mljjCut1     = {55.,  160., 225., 295., 370., 370., 370., 370., 370.,  370.,  370.,  370.,  370.,  370.,  370.,  370.};
            mljjCut2     = {115., 215., 340., 490., 550., 630., 885., 890., 1225., 1230., 1245., 1690., 1890., 2220., 2700., 3200.};
            MET2STCut    = {9.,   7.,   7.,   7.,   7.,   7.,   7.,   7.,   7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.};
          }
          if(channel=="diel"){
            Lep1PtCutSR3 = {25.,  55.,  80.,  100., 125., 125., 125., 125.,  125.,  125.,  125.,  125.,  125.,  125.,  125.,  125.};
            Lep2PtCutSR3 = {15.,  40.,  60.,  65.,  65.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
            mlljjCut     = {120., 220., 370., 450., 560., 760., 760., 760.,  760.,  760.,  760.,  760.,  760.,  760.,  760.,  760.};
            mljjCut1     = {50.,  160., 235., 335., 400., 400., 400., 400.,  400.,  400.,  400.,  400.,  400.,  400.,  400.,  400.};
            mljjCut2     = {110., 225., 335., 450., 555., 690., 955., 1130., 1300., 1490., 1490., 1600., 1930., 1930., 2400., 2900.};
            MET2STCut    = {6.,   6.,   6.,   6.,   6.,   6.,   6.,   6.,    6.,    6.,    6.,    6.,    6.,    6.,    6.,    6.};
          }
          if(channel=="emu"){
            Lep1PtCutSR3 = {25.,  65.,  95.,  120., 150., 175., 180.,  180.,  185.,  185.,  185.,  185.,  185.,  185.,  185.,  185.};
            Lep2PtCutSR3 = {20.,  35.,  60.,  60.,  60.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
            mlljjCut     = {110., 270., 340., 530., 580., 670., 720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.};
            mljjCut1     = {60.,  170., 255., 325., 315., 315., 350.,  400.,  450.,  500.,  550.,  600.,  650.,  650.,  650.,  650.};
            mljjCut2     = {115., 230., 325., 450., 530., 740., 1030., 1030., 1040., 1415., 1640., 1780., 1880., 1885., 2400., 2900.};
            MET2STCut    = {7.,   7.,   7.,   7.,   7.,   7.,   7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.};
          }

          for(unsigned int it_m=0; it_m<mass.size(); it_m++){

            if(!(jets_eta2p7.size() < 4)) continue;
            if(!(jets_WCand.at(0).Pt() > 25.)) continue;
            if(!(leptons.at(0)->Pt() > Lep1PtCutSR3.at(it_m))) continue;
            if(!(leptons.at(1)->Pt() > Lep2PtCutSR3.at(it_m))) continue;
            if(!(WCand.M()>50. && WCand.M()<120.)) continue;
            if(!(lljj.M() > mlljjCut.at(it_m))) continue;

            if(it_m < 2){
              if(!(dRl2jj < 3.1)) continue;
              if(!(l2jj.M()>mljjCut1.at(it_m) && l2jj.M()<mljjCut2.at(it_m))) continue;
            }
            else{
              if(!(l1jj.M()>mljjCut1.at(it_m) && l1jj.M()<mljjCut2.at(it_m)) && !(l2jj.M()>mljjCut1.at(it_m) && l2jj.M()<mljjCut2.at(it_m))) continue;
            }
            if(!(MET2ST < MET2STCut.at(it_m))) continue;

            FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_"+IDName, 0.5, weight, 2, 0., 2.);
            //FillHist(systName+"_"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_unweighted_"+IDName, 0.5, 1., 2, 0., 2.);

            /*if(channel=="diel"){
              if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
                FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_"+IDName, 0.5, weight, 1, 0., 1.);
                FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_unweighted_"+IDName, 0.5, 1., 1, 0., 1.);
              }
            }*/

          }

        }

      }

    }

    //jj = jets.at(0) + jets.at(1);
    //double avgEta = 0.5*(jets.at(0).Eta() + jets.at(1).Eta());
    //double dEta = fabs(jets.at(0).Eta() - jets.at(1).Eta());
    //double max_zep = std::max(fabs(muons.at(0).Eta()-avgEta), fabs(muons.at(1).Eta()-avgEta))/dEta;

    //==== SSWW selection
    //==== pT(l) > 30/30 GeV, 3rd lepton veto, m(ll) > 20 GeV, pT(j) > 30/30 GeV, max_zep < 0.75, dEta > 2.5, m(jj) > 750 GeV, N(b) = 0

  }

}
