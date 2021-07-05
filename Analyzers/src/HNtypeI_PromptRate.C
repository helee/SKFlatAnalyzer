#include "HNtypeI_PromptRate.h"

HNtypeI_PromptRate::HNtypeI_PromptRate(){

}

void HNtypeI_PromptRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst     = HasFlag("RunSyst");
  //RunMET      = HasFlag("RunMET");
  RunFake     = HasFlag("RunFake");

  cout << "[HNtypeI_PromptRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  //cout << "[HNtypeI_PromptRate::initializeAnalyzer] RunMET = " << RunMET << endl;
  cout << "[HNtypeI_PromptRate::initializeAnalyzer] RunFake = " << RunFake << endl;

  MuonTightIDs     = {"HNTightV2"};
  MuonLooseIDs     = {"HNLooseV2"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};
  MuonFRNames      = {"HNRun2"};
  ElectronFRNames  = {"HNRun2"};
  //MuonFRNames      = {"HNRun2METPhi"};
  //ElectronFRNames  = {"HNRun2METPhi"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_PromptRate.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  MuonTriggers.clear();
  MuonTriggersH.clear();
  ElectronTriggers.clear();
  EMuTriggers.clear();
  EMuTriggersH.clear();
  Mu8Ele23Triggers.clear();
  Mu23Ele12Triggers.clear();
  Mu8Ele23TriggersH.clear();
  Mu23Ele12TriggersH.clear();

  if(DataYear==2016){                                                                   // Lumi values for trigger weight (/pb)
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");                       // 27267.591112919
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");                     // 27267.591112919
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                    // 35918.219492947
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                  // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");                   // 35918.219492947
    MuonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");                 // 35918.219492947
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 27267.591112919
    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
    EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028

    //==== These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    Mu8Ele23TriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12TriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 17., MuonPtCut2 = 8.;
    ElectronPtCut1 = 23., ElectronPtCut2 = 12.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    //==== These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 17., MuonPtCut2 = 8.;
    ElectronPtCut1 = 23., ElectronPtCut2 = 12.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }
  else if(DataYear==2018){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    //==== These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 17., MuonPtCut2 = 8.;
    ElectronPtCut1 = 23., ElectronPtCut2 = 12.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

  //cout << "[HNtypeI_PromptRate::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_PromptRate::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== b tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_PromptRate::~HNtypeI_PromptRate(){

  //==== Destructor of this Analyzer

}

void HNtypeI_PromptRate::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_PromptRate.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_PromptRate.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

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

    //param.Name = MuonID+"_"+"Central";
    param.Name = "Central";

    //==== Muon ID
    param.Muon_Tight_ID       = MuonTightID;
    param.Muon_Loose_ID       = MuonLooseID;
    param.Muon_Veto_ID        = MuonVetoID;
    param.Muon_FR_ID          = MuonFRName;     // ID name in histmap_Muon.txt
    param.Muon_FR_Key         = "AwayJetPt40";  // histname
    //param.Muon_ID_SF_Key      = "NUM_TightID_DEN_genTracks";
    //param.Muon_ISO_SF_Key     = "NUM_TightRelIso_DEN_TightIDandIPCut";
    param.Muon_ID_SF_Key      = "";
    param.Muon_ISO_SF_Key     = "";
    param.Muon_Trigger_SF_Key = "";
    param.Muon_UsePtCone      = true;

    //==== Electron ID
    param.Electron_Tight_ID       = ElectronTightID;
    //if(DataYear==2018 && it_id==1) param.Electron_Tight_ID = "HEEP2018_dZ";
    param.Electron_Loose_ID       = ElectronLooseID;
    param.Electron_Veto_ID        = ElectronVetoID;
    param.Electron_FR_ID          = ElectronFRName; // ID name in histmap_Electron.txt
    param.Electron_FR_Key         = "AwayJetPt40";  // histname
    //param.Electron_ID_SF_Key      = "passTightID";
    param.Electron_ID_SF_Key      = "";
    param.Electron_Trigger_SF_Key = "";
    param.Electron_UsePtCone      = true;

    //==== Jet ID
    //param.Jet_ID = "tightLepVeto";
    param.Jet_ID = "HNTight";
    if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45";

    executeEventFromParameter(param);

    if(RunSyst){
      //for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
      for(int it_syst=1; it_syst<31; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name  = "Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }

  }

}

void HNtypeI_PromptRate::executeEventFromParameter(AnalyzerParameter param){

  TString IDsuffix = "HNRun2";

  TString channel = "";
  //vector<TString> regions = {"DY", "TT"};

  TString systName = param.Name;

  TString NjetBin = "", VtxBin = "";
  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., diel_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;
 
  Event ev = GetEvent();

  //==== Boolean : primary datasets
  bool isDoubleMuon = false, isDoubleEG = false, isMuonEG = false;
  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) isDoubleMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) isDoubleEG = true;
    if(DataStream.Contains("MuonEG")) isMuonEG = true;
  }

  //==== Boolean : passTrigger
  bool passMuMu  = ev.PassTrigger(MuonTriggers);
  bool passEE    = ev.PassTrigger(ElectronTriggers);
  bool passEMu   = ev.PassTrigger(EMuTriggers);
  bool passE23Mu = ev.PassTrigger(Mu8Ele23Triggers);
  bool passEMu23 = ev.PassTrigger(Mu23Ele12Triggers);

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
  //if(!IsDATA) Nvtx = nPileUp+1;
  FillHist(systName+"/Number_Vertices_NoCut", Nvtx, weight, 100, 0., 100.);

  //FillHist(systName+"/"+"nPileUp", nPileUp, weight, 200., 0., 200.);
  //FillHist(systName+"/"+"Nvtx", Nvtx, weight, 200., 0., 200.);

  //==== Cutflow : No Cuts
  /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    if(!IsDATA || isDoubleMuon){
      FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isDoubleEG){
      FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isMuonEG){
      FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }*/

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //==== Cutflow : MET filter
  /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    if(!IsDATA || isDoubleMuon){
      FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isDoubleEG){
      FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isMuonEG){
      FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }*/

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
  int systL1 = 0, systPU = 0, systMuonID = 0, systElectronReco = 0, systElectronID = 0, systMuonTrigger = 0, systElectronTrigger = 0, systEMuTrigger = 0;

  if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(2), pfMET_Type1_phi_shifts->at(2));
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(3), pfMET_Type1_phi_shifts->at(3));
  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(0), pfMET_Type1_phi_shifts->at(0));
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(1), pfMET_Type1_phi_shifts->at(1));
  }
  else if(param.syst_ == AnalyzerParameter::UnclusteredEnUp){
    ev.SetMET(pfMET_Type1_pt_shifts->at(10), pfMET_Type1_phi_shifts->at(10));
  }
  else if(param.syst_ == AnalyzerParameter::UnclusteredEnDown){
    ev.SetMET(pfMET_Type1_pt_shifts->at(11), pfMET_Type1_phi_shifts->at(11));
  }
  else if(param.syst_ == AnalyzerParameter::BtagSFUp){
    systBtag = "up";
  }
  else if(param.syst_ == AnalyzerParameter::BtagSFDown){
    systBtag = "down";
  }
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
  else if(param.syst_ == AnalyzerParameter::MuonTriggerSFUp){
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
  }
  else if(param.syst_ == AnalyzerParameter::EMuTriggerSFUp){
    systEMuTrigger = 1;
  }
  else if(param.syst_ == AnalyzerParameter::EMuTriggerSFDown){
    systEMuTrigger = -1;
  }
  else{
    cout << "[HNtypeI_PromptRate::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  //==== Leptons
  vector<Muon> muons_tight, muons_loose, muons_veto, muons_POGTight, muons_POGLoose;
  muons_tight.clear();
  muons_loose.clear();

  vector<Electron> electrons_tight, electrons_loose, electrons_veto;
  electrons_tight.clear();
  electrons_loose.clear();
  electrons_veto.clear();

  muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, 5., 2.4);
  muons_veto  = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 5., 2.4);

  electrons_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  electrons_loose = SelectElectrons(this_AllElectrons, param.Electron_Loose_ID, 10., 2.5);
  electrons_veto  = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 10., 2.5);

  //==== Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  //==== Charge flip
  /*vector<Electron> electrons_beforeShift;
  vector<Electron> electrons_afterShift;
  electrons_beforeShift.clear();
  electrons_afterShift.clear();*/

  //==== Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 10., 2.7);  // AK4jets used for b tag
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

  //vector<Jet> jets_Pt10to18;
  vector<Jet> jets_Pt20;
  vector<Jet> jets_Pt30;
  //jets_Pt10to18.clear();
  jets_Pt20.clear();
  jets_Pt30.clear();

  /*for(unsigned int i=0; i<jets_lepveto.size(); i++){
    if(jets_lepveto.at(i).Pt()>10. && jets_lepveto.at(i).Pt()<18.) jets_Pt10to18.push_back(jets_lepveto.at(i));
  }*/

  for(unsigned int i=0; i<jets.size(); i++){
    if(jets.at(i).Pt() > 20.) jets_Pt20.push_back(jets.at(i));
    if(jets.at(i).Pt() > 30.) jets_Pt30.push_back(jets.at(i));
  }
  //vector<Jet> jets = JetsAwayFromFatJet(jets_PUveto, fatjets);
  
  std::vector<Lepton*> leptons, leptons_veto;

  //========================================================
  //==== Sort in pT-order
  //========================================================

  std::sort(muons_tight.begin(), muons_tight.end(), PtComparing);
  std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons_tight.begin(), electrons_tight.end(), PtComparing);
  std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);
  //std::sort(jets_Pt10to18.begin(), jets_Pt10to18.end(), PtComparing);
  std::sort(jets_Pt20.begin(), jets_Pt20.end(), PtComparing);
  std::sort(jets_Pt30.begin(), jets_Pt30.end(), PtComparing);

  //========================================================
  //==== b tagging
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0, Nbjet_Pt30_loose = 0., Nbjet_Pt30_medium = 0, Nbjet_medium_lepveto = 0;
  JetTagging::Parameters jtp_DeepCSV_Loose = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_nolepveto.size(); ij++){
    //double this_discr = jets_nolepveto.at(ij).GetTaggerResult(JetTagging::DeepCSV);
    //==== No SF
    //if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;
    if(jets_nolepveto.at(ij).Pt() > 20.){
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij), systBtag)) Nbjet_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij), systBtag)) Nbjet_medium++;
    }
    if(jets_nolepveto.at(ij).Pt() > 30.){
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij), systBtag)) Nbjet_Pt30_loose++;
      if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij), systBtag)) Nbjet_Pt30_medium++;
    }
  }

  for(unsigned int ij=0; ij<jets_Pt20.size(); ij++){
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_Pt20.at(ij), systBtag)) Nbjet_medium_lepveto++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  //ev.SetMET(pfMET_Type1_PhiCor_pt, pfMET_Type1_PhiCor_phi);
  //Particle METv_dxy = ev.GetMETVector();

  if(muons_loose.size()+electrons_loose.size() == 2){
    METv = UpdateMETMuon(METv, muons_loose);
    METv = UpdateMETElectron(METv, electrons_loose);
    //METv_dxy = UpdateMETMuon(METv_dxy, muons);
    //METv_dxy = UpdateMETElectron(METv_dxy, electrons);
  }

  double MET = METv.Pt();
  double METPhi = METv.Phi();
  //double MET_dxy = METv_dxy.Pt();

  //========================================================
  //==== Define particles, variables
  //========================================================

  //double ST = 0., MET2ST = 0.;
  //double dRll = 0., dPhill = 0., PtDiff = 0.;
  double MZ = 91.1876;
  //double MW = 80.379;
  //double mllCut = 55.;
  double muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1.;
  int lepton_veto_size = 0;
  //double lepton1_eta = 0., lepton2_eta = 0.;

  bool passPtCut = false;
  Particle ZCand;
  //Particle llj, l1j, l2j,  lljj, l1jj, l2jj, l1J, l2J;
  
  //==== Set up pTcone
  double mu_tight_iso = 0.07, el_tight_iso = 0.;
  double ptcone_mu1 = 0., ptcone_mu2 = 0., ptcone_el1 = 0., ptcone_el2 = 0.;
  //double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(muons_loose.size() == 2){
    ptcone_mu1 = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);
    ptcone_mu2 = muons_loose.at(1).CalcPtCone(muons_loose.at(1).RelIso(), mu_tight_iso);
  }

  if(electrons_loose.size() == 2){

    if(param.Electron_Tight_ID.Contains("HNTight")){ // POG CB Tight
      el_tight_iso = 0.0287+0.506/electrons_loose.at(0).UncorrPt();
      if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons_loose.at(0).UncorrPt();
    }

    ptcone_el1 = electrons_loose.at(0).CalcPtCone(electrons_loose.at(0).RelIso(), el_tight_iso);
    ptcone_el2 = electrons_loose.at(1).CalcPtCone(electrons_loose.at(1).RelIso(), el_tight_iso);

  }

  //==== Define leptons (pT order)
  for(unsigned int i=0; i<muons_loose.size(); i++) leptons.push_back(&muons_loose.at(i));
  for(unsigned int i=0; i<electrons_loose.size(); i++) leptons.push_back(&electrons_loose.at(i));
  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  //==== Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  //==== Leptons (minus, plus charge)
  /*for(unsigned int i=0; i<muons.size(); i++){
    if(muons.at(i).Charge() < 0) leptons_minus.push_back(&muons.at(i));
    if(muons.at(i).Charge() > 0) leptons_plus.push_back(&muons.at(i));
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    if(electrons.at(i).Charge() < 0) leptons_minus.push_back(&electrons.at(i));
    if(electrons.at(i).Charge() > 0) leptons_plus.push_back(&electrons.at(i));
  }*/

  lepton_veto_size = leptons_veto.size() - leptons.size();

  //==== Define ST, MET^2/ST
  MET = METv.Pt();
  METPhi = METv.Phi();
  //MET_dxy = METv_dxy.Pt();

  /*for(unsigned int i=0; i<jets_Pt20.size(); i++) ST += jets_Pt20.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;*/

  //=====================================================================================
  //==== SM background CR (DY, TT)
  //=====================================================================================

  //==== Period-dependent trigger weight (only for 2016 MC)
  trigger_lumi = 1., dimu_trig_weight = 0., diel_trig_weight = 0., emu_trig_weight = 0.;
  if(!IsDATA){
    if(DataYear==2016){
      if(ev.PassTrigger(MuonTriggers))  dimu_trig_weight += 27267.591;
      if(ev.PassTrigger(MuonTriggersH)) dimu_trig_weight += 8650.628;
      diel_trig_weight = ev.GetTriggerLumi("Full");
      if(ev.PassTrigger(EMuTriggers))  emu_trig_weight += 27267.591;
      if(ev.PassTrigger(EMuTriggersH)) emu_trig_weight += 8650.628;
    }
    else{
      dimu_trig_weight = ev.GetTriggerLumi("Full");
      diel_trig_weight = ev.GetTriggerLumi("Full");
      emu_trig_weight  = ev.GetTriggerLumi("Full");
    }
  }

  //========================================================
  //==== Event selection
  //========================================================

  /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    weight = 1., muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1.;

    if(!IsDATA){
      weight *= weight_norm_1invpb;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(systL1);
      weight *= GetPileUpWeight(nPileUp, systPU);
    }
    //==== Cutflow : passing dilepton triggers (dimu || diel || emu)
    if(passMuMu){
      if(!IsDATA || isDoubleMuon){
        if(!IsDATA) trigger_lumi = dimu_trig_weight;
        FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/dimu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
    if(passEE){
      if(!IsDATA || isDoubleEG){
        if(!IsDATA) trigger_lumi = diel_trig_weight;
        FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/diel_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
    if(passEMu){
      if(!IsDATA || isMuonEG){
        if(!IsDATA) trigger_lumi = emu_trig_weight;
        FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/emu_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
  }*/

  //========================================================
  //==== Prompt rate measurement region
  //========================================================

  //==== Cutflow : 2 tight leptons (gen-matched, pT > trigger thresholds)

  if(leptons.size() == 2){

    //==== pTcone,pT > trigger thresholds
    passPtCut = false;

    if(muons_loose.size()==2 && electrons_loose.size()==0){
      if(!passMuMu) return;
      if(!IsDATA) trigger_lumi = dimu_trig_weight;
      if(IsDATA){ if(!isDoubleMuon) return; }
      if(ptcone_mu1>20. && ptcone_mu2>10. && muons_loose.at(0).Pt()>MuonPtCut1 && muons_loose.at(1).Pt()>MuonPtCut2) passPtCut = true;
      channel = "dimu";
    }
    if(muons_loose.size()==0 && electrons_loose.size()==2){
      if(!passEE) return;
      if(!IsDATA) trigger_lumi = diel_trig_weight;
      if(IsDATA){ if(!isDoubleEG) return; }
      if(ptcone_el1>25. && ptcone_el2>15. && electrons_loose.at(0).Pt()>ElectronPtCut1 && electrons_loose.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      channel = "diel";
    }
    if(muons_loose.size()==1 && electrons_loose.size()==1){
      return;
      /*if(!passEMu) return;
      if(!IsDATA) trigger_lumi = emu_trig_weight;
      if(IsDATA){ if(!isMuonEG) return; }
      if(passE23Mu){
        if(electrons.at(0).Pt()>EMuPtCut1 && muons.at(0).Pt()>EMuPtCut2) passPtCut = true;
      }
      if(passEMu23){
        if(muons.at(0).Pt()>EMuPtCut1 && electrons.at(0).Pt()>EMuPtCut2) passPtCut = true;;
      }
      channel = "emu";*/
    }

    if(!passPtCut) return;

    //==== Truth matching
    muons_prompt.clear();
    electrons_prompt.clear();
    muons_prompt = MuonPromptOnlyHNtypeI(muons_loose, gens);
    electrons_prompt = ElectronPromptOnlyHNtypeI(electrons_loose, gens);

    if(channel=="dimu"){
      if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) return;
    }
    if(channel=="diel"){
      if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) return;
    }
    /*if(channel=="emu"){
      if(!(muons_prompt.size()==1 && electrons_prompt.size()==1)) return;
    }*/

    //==== Event weights for MC
    weight = 1.;

    if(!IsDATA){

      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(systL1);
      weight *= GetPileUpWeight(nPileUp, systPU);

      //PUweight_up = GetPileUpWeight(nPileUp, 1);
      //PUweight_down = GetPileUpWeight(nPileUp, -1);

      /*for(unsigned int i=0; i<muons.size(); i++){

        if(param.Muon_Tight_ID.Contains("HighPt")){
          muon_miniaodP = sqrt( muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz() );
          muonRecoSF    = mcCorr->MuonReco_SF("HighPtMuonRecoSF", muons.at(i).Eta(), muon_miniaodP, 0);
          muonIDSF      = mcCorr->MuonID_SF("NUM_HighPtID_DEN_genTracks",  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          muonIsoSF     = mcCorr->MuonISO_SF("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut", muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        else if(param.Muon_Tight_ID.Contains("HNTight")){
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

      for(unsigned int i=0; i<electrons.size(); i++){

        electronRecoSF = mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), systElectronReco);

        if(param.Electron_Tight_ID.Contains("HEEP")){
          electronIDSF = mcCorr->ElectronID_SF("HEEP", electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
        }
        else if(param.Electron_Tight_ID.Contains("HNTight")){
          electronIDSF = mcCorr->ElectronID_SF(param.Electron_Tight_ID, electrons.at(i).scEta(), electrons.at(i).UncorrPt(), systElectronID);
          if(RunFake){  // When subtracting prompt contribution from fake contribution, we apply ID SF only for electrons passing the tight ID
            if(!electrons.at(i).PassID(param.Electron_Tight_ID)) electronIDSF = 1.;
          }
        }
        else electronIDSF = 1.;

        weight *= electronRecoSF*electronIDSF;

      }

      if(channel=="dimu") triggerSF = mcCorr->MuonTrigger_SF_HNtypeI(param.Muon_Tight_ID, muons, RunFake, systMuonTrigger);
      if(channel=="diel") triggerSF = mcCorr->ElectronTrigger_SF_HNtypeI(param.Electron_Tight_ID, electrons, RunFake, systElectronTrigger);
      //if(channel=="emu")  triggerSF = mcCorr->EMuTrigger_SF_HNtypeI(param.Muon_Tight_ID, param.Electron_Tight_ID, muons, electrons, RunFake, systEMuTrigger);

      weight *= triggerSF;*/

    }
    if(RunFake) weight *= fakeEst->GetWeight(leptons, param);

    ZCand = *leptons.at(0) + *leptons.at(1);
    //dRll   = leptons.at(0)->DeltaR(*leptons.at(1));
    //dPhill = fabs(leptons.at(0)->DeltaPhi(*leptons.at(1)));
    //PtDiff = fabs(leptons.at(0)->Pt() - leptons.at(1)->Pt())/(leptons.at(0)->Pt() + leptons.at(1)->Pt());

    /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
    }*/

    //==== Cutflow : OS event
    if(leptons.at(0)->Charge()*leptons.at(1)->Charge() > 0) return;

    /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
    }*/

    //==== Cutflow : veto 3rd leptons using veto ID
    if(lepton_veto_size > 0) return;

    /*for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
    }*/

    //==== Cutflow : |m(ll) - m(Z)| < 10 GeV
    if(!IsOnZ(ZCand.M(), 10.)) return;

    //==== Cutflow : veto b jets
    if(!(Nbjet_medium == 0)) return;

    FillHist(systName+"/"+channel+"_2L_ZCand_Mass", ZCand.M(), weight, 40, 70., 110.);
    if(muons_tight.size()>0 || electrons_tight.size()>0){
      FillHist(systName+"/"+channel+"_1T_ZCand_Mass", ZCand.M(), weight, 40, 70., 110.);
    }

    if(channel == "dimu"){

      //==== Two loose muons
      if(fabs(muons_loose.at(0).Eta()) < 0.8){
        FillHist(systName+"/"+channel+"_2L_IB_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
      }
      if(fabs(muons_loose.at(0).Eta())>=0.8 && fabs(muons_loose.at(0).Eta())<1.479){
        FillHist(systName+"/"+channel+"_2L_OB_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
      }
      if(fabs(muons_loose.at(0).Eta()) > 1.479){
        FillHist(systName+"/"+channel+"_2L_EC_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
        }
      }

      if(fabs(muons_loose.at(1).Eta()) < 0.8){
        FillHist(systName+"/"+channel+"_2L_IB_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
      }
      if(fabs(muons_loose.at(1).Eta())>=0.8 && fabs(muons_loose.at(1).Eta())<1.479){
        FillHist(systName+"/"+channel+"_2L_OB_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
      }
      if(fabs(muons_loose.at(1).Eta()) > 1.479){
        FillHist(systName+"/"+channel+"_2L_EC_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
        }
      }

      //==== At least one tight muon (like TnP)
      if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){

        if(fabs(muons_loose.at(1).Eta()) < 0.8){
          FillHist(systName+"/"+channel+"_1T_IB_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
        }
        if(fabs(muons_loose.at(1).Eta())>=0.8 && fabs(muons_loose.at(1).Eta())<1.479){
          FillHist(systName+"/"+channel+"_1T_OB_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
        }
        if(fabs(muons_loose.at(1).Eta()) > 1.479){
          FillHist(systName+"/"+channel+"_1T_EC_Lepton_Loose_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_Tight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_NoTight_PtCone", ptcone_mu2, weight, 500, 0., 500.);
          }
        }

      }

      if(muons_loose.at(1).PassID(param.Muon_Tight_ID)){

        if(fabs(muons_loose.at(0).Eta()) < 0.8){
          FillHist(systName+"/"+channel+"_1T_IB_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
        }
        if(fabs(muons_loose.at(0).Eta())>=0.8 && fabs(muons_loose.at(0).Eta())<1.479){
          FillHist(systName+"/"+channel+"_1T_OB_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
        }
        if(fabs(muons_loose.at(0).Eta()) > 1.479){
          FillHist(systName+"/"+channel+"_1T_EC_Lepton_Loose_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          if(muons_loose.at(0).PassID(param.Muon_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_Tight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_NoTight_PtCone", ptcone_mu1, weight, 500, 0., 500.);
          }
        }

      }

    }

    if(channel == "diel"){

      //==== Two loose electrons
      if(fabs(electrons_loose.at(0).scEta()) < 0.8){
        FillHist(systName+"/"+channel+"_2L_IB_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
        if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        } 
      } 
      if(fabs(electrons_loose.at(0).scEta())>=0.8 && fabs(electrons_loose.at(0).scEta())<1.479){
        FillHist(systName+"/"+channel+"_2L_OB_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
        if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        }
      } 
      if(fabs(electrons_loose.at(0).scEta()) > 1.479){
        FillHist(systName+"/"+channel+"_2L_EC_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
        if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
        } 
      }

      if(fabs(electrons_loose.at(1).scEta()) < 0.8){
        FillHist(systName+"/"+channel+"_2L_IB_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
        if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_IB_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
      }
      if(fabs(electrons_loose.at(1).scEta())>=0.8 && fabs(electrons_loose.at(1).scEta())<1.479){
        FillHist(systName+"/"+channel+"_2L_OB_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
        if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_OB_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
      }
      if(fabs(electrons_loose.at(1).scEta()) > 1.479){
        FillHist(systName+"/"+channel+"_2L_EC_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
        if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
        else{
          FillHist(systName+"/"+channel+"_2L_EC_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
        }
      }

      //==== At least one tight electron (like TnP)
      if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){

        if(fabs(electrons_loose.at(1).scEta()) < 0.8){
          FillHist(systName+"/"+channel+"_1T_IB_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
          if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
        }
        if(fabs(electrons_loose.at(1).scEta())>=0.8 && fabs(electrons_loose.at(1).scEta())<1.479){
          FillHist(systName+"/"+channel+"_1T_OB_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
          if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
        }
        if(fabs(electrons_loose.at(1).scEta()) > 1.479){
          FillHist(systName+"/"+channel+"_1T_EC_Lepton_Loose_PtCone", ptcone_el2, weight, 500, 0., 500.);
          if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_Tight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_NoTight_PtCone", ptcone_el2, weight, 500, 0., 500.);
          }
        }

      }

      if(electrons_loose.at(1).PassID(param.Electron_Tight_ID)){

        if(fabs(electrons_loose.at(0).scEta()) < 0.8){
          FillHist(systName+"/"+channel+"_1T_IB_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
          if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_IB_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
        }
        if(fabs(electrons_loose.at(0).scEta())>=0.8 && fabs(electrons_loose.at(0).scEta())<1.479){
          FillHist(systName+"/"+channel+"_1T_OB_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
          if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_OB_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
        }
        if(fabs(electrons_loose.at(0).scEta()) > 1.479){
          FillHist(systName+"/"+channel+"_1T_EC_Lepton_Loose_PtCone", ptcone_el1, weight, 500, 0., 500.);
          if(electrons_loose.at(0).PassID(param.Electron_Tight_ID)){
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_Tight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
          else{
            FillHist(systName+"/"+channel+"_1T_EC_Lepton_NoTight_PtCone", ptcone_el1, weight, 500, 0., 500.);
          }
        }

      }

    }

  }

}



