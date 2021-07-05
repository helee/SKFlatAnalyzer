#include "HNtypeI_SR_2016H.h"

HNtypeI_SR_2016H::HNtypeI_SR_2016H(){

}

void HNtypeI_SR_2016H::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  RunFake = HasFlag("RunFake");
  RunCF   = HasFlag("RunCF");
  RunOS   = HasFlag("RunOS");

  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunFake = " << RunFake << endl;
  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunCF = " << RunCF << endl;
  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunOS = " << RunOS << endl;

  MuonTightIDs     = {"HNTightV2"};
  MuonLooseIDs     = {"HNLooseV2"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};
  MuonFRNames      = {"HNTightV2"};
  ElectronFRNames  = {"HNTightV2"};
  //MuonFRNames      = {"HNRun2METPhi"};
  //ElectronFRNames  = {"HNRun2METPhi"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_SR_2016H.h 
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
    //MuonTriggers.push_back("HLT_Mu17_Mu8_SameSign_DZ_v");
    //MuonTriggersH.push_back("HLT_Mu17_Mu8_SameSign_DZ_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");          // 35918.219492947
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");         // 27267.591112919
    //EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");          // 27267.591112919
    EMuTriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028
    EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 8650.628380028
    //EMuTriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");      // 8650.628380028

    //==== These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    Mu8Ele23TriggersH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12TriggersH.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    //Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
    //EMuPtCut1 = 25., EMuPtCut2 = 10.;
  }
  else if(DataYear==2017){
    MuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    ElectronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    EMuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    EMuTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    //==== These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
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

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

  //cout << "[HNtypeI_SR_2016H::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_SR_2016H::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== b tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== Set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_SR_2016H::~HNtypeI_SR_2016H(){

  //==== Destructor of this Analyzer

}

void HNtypeI_SR_2016H::executeEvent(){

  //================================================================
  //====  Example 1
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/HNtypeI_SR_2016H.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();
  //AllFatJets = GetAllFatJets();
  AllFatJets = puppiCorr->Correct(GetAllFatJets());

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/HNtypeI_SR_2016H.h
  //weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<MuonTightIDs.size(); it_id++){
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
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name  = "Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }

  }

}

void HNtypeI_SR_2016H::executeEventFromParameter(AnalyzerParameter param){

  TString IDsuffix = "HNTightV2";

  TString channel = "";
  vector<TString> regions = {"LowSR1", "LowCR1", "HighSR1", "HighCR1", "LowSR2", "LowCR2", "HighSR2", "HighCR2"};
  vector<double> Lep1PtCutSR1, Lep1PtCutSR2, Lep2PtCutSR1, Lep2PtCutSR2, mlljjCut, mljjCut1, mljjCut2, mlJCut1, mlJCut2, MET2STCut;
  vector<int> mass = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1500, 2000, 2500};
  TString tight_leptons = "";

  TString systName = param.Name;

  double cutflow_max = 12.;
  int cutflow_bin = 12;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., diel_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;
  int tight_muons = 0, tight_electrons = 0;

  Event ev = GetEvent();
 
  //==== Boolean : primary datasets
  bool isDoubleMuon = false, isDoubleEG = false, isMuonEG = false;
  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) isDoubleMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) isDoubleEG = true;
    if(DataStream.Contains("MuonEG")) isMuonEG = true;
  }

  //==== Boolean : passTrigger
  bool passMuMu  = ev.PassTrigger(MuonTriggersH);      // NOTE : Change for 2016H
  bool passEE    = ev.PassTrigger(ElectronTriggers);
  bool passEMu   = ev.PassTrigger(EMuTriggersH);       // NOTE : Change for 2016H
  bool passEMu8  = ev.PassTrigger(Mu8Ele23TriggersH);  // NOTE : Change for 2016H
  bool passEMu23 = ev.PassTrigger(Mu23Ele12TriggersH); // NOTE : Change for 2016H

  //========================================================
  //==== No Cut
  //========================================================

  weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp, 0);
  }

  int Nvtx = nPV;
  FillHist(systName+"/Number_Vertices_NoCut", Nvtx, weight, 100, 0., 100.);

  //==== Cutflow : No Cuts
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
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
  }
  if(!IsDATA || isDoubleMuon){
    FillHist(systName+"/dimu_FakeCR1_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR2_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isDoubleEG){
    FillHist(systName+"/diel_FakeCR1_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR1_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR2_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR2_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isMuonEG){
    FillHist(systName+"/emu_FakeCR1_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR2_Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //==== Cutflow : MET filter
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
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
  }
  if(!IsDATA || isDoubleMuon){
    FillHist(systName+"/dimu_FakeCR1_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR2_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isDoubleEG){
    FillHist(systName+"/diel_FakeCR1_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR1_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR2_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel_FakeCR2_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isMuonEG){
    FillHist(systName+"/emu_FakeCR1_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR2_Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
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
  int systL1 = 0, systPU = 0, systMuonID = 0, systElectronReco = 0, systElectronID = 0, systMuonTrigger = 0, systElectronTrigger = 0, systEMuTrigger = 0, systTau21 = 0;

  if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    this_AllFatJets = ScaleFatJets( this_AllFatJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(2), pfMET_Type1_phi_shifts->at(2));
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    this_AllFatJets = ScaleFatJets( this_AllFatJets, -1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(3), pfMET_Type1_phi_shifts->at(3));
  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    this_AllFatJets = SmearFatJets( this_AllFatJets, +1 );
    ev.SetMET(pfMET_Type1_pt_shifts->at(0), pfMET_Type1_phi_shifts->at(0));
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    this_AllFatJets = SmearFatJets( this_AllFatJets, -1 );
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
  else if(param.syst_ == AnalyzerParameter::SDMassScaleUp){
    this_AllFatJets = ScaleSDMassFatJets( this_AllFatJets, +1);
  }
  else if(param.syst_ == AnalyzerParameter::SDMassScaleDown){
    this_AllFatJets = ScaleSDMassFatJets( this_AllFatJets, -1);
  }
  else if(param.syst_ == AnalyzerParameter::SDMassResUp){
    this_AllFatJets = SmearSDMassFatJets( this_AllFatJets, +1);
  }
  else if(param.syst_ == AnalyzerParameter::SDMassResDown){
    this_AllFatJets = SmearSDMassFatJets( this_AllFatJets, -1);
  }
  else if(param.syst_ == AnalyzerParameter::Tau21SFUp){
    systTau21 = 1;
  }
  else if(param.syst_ == AnalyzerParameter::Tau21SFDown){
    systTau21 = -1;
  }
  else{
    cerr << "[HNtypeI_SR_2016H::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

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
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  //==== Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  //==== Fatjet selection in CATanalyzer (see the links)
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  //==== https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  //vector<Jet> jets_PUveto = JetsPassPileupMVA(jets_lepveto);
  vector<Jet> jets = JetsAwayFromFatJet(jets_lepveto, fatjets);  // AK4jets used in SR, CR

  //vector<Jet> jets_WCandLowMass;
  vector<Jet> jets_WCandHighMass;
  FatJet fatjets_WCand;
  //jets_WCandLowMass.clear();
  jets_WCandHighMass.clear();

  std::vector<Lepton*> leptons, leptons_veto;

  //========================================================
  //==== Sort in pT-order
  //========================================================

  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(electrons.begin(), electrons.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  std::sort(jets.begin(), jets.end(), PtComparing);
  std::sort(jets_nolepveto.begin(), jets_nolepveto.end(), PtComparing);
  std::sort(fatjets.begin(), fatjets.end(), PtComparing);

  //========================================================
  //==== b tagging 
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0;
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
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij), systBtag)) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij), systBtag)) Nbjet_medium++;
  }

  //========================================================
  //==== Set up MET
  //========================================================

  Particle METv = ev.GetMETVector();

  if(muons.size()+electrons.size() == 2){
    METv = UpdateMETMuon(METv, muons);
    METv = UpdateMETElectron(METv, electrons);
  }

  double MET = METv.Pt();
  double METPhi = METv.Phi();

  //==== Correct MET with jets inside a fatjet
  if(systName == "Syst_JetResUp") METv = UpdateMETJet(METv, jets_insideFatjets, 1);
  if(systName == "Syst_JetResDown") METv = UpdateMETJet(METv, jets_insideFatjets, -1);

  //========================================================
  //==== Define particles, variables
  //========================================================

  double ST = 0., MET2ST = 0.;
  double dRll = 0., dPhill = 0., dRl1jj = 0., dRl2jj = 0., dRl1J = 0., dRl2J = 0., PtDiff = 0., dRjj = 0.;
  //double MZ = 91.1876;
  double MW = 80.379;
  double mllCut1 = 10., mllCut2 = 55.;  // mllCut1 : 10 GeV cut in EXO-17-028, mllCut2 : No DY10to50 MC samples
  double muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., triggerSF = 1., fatjetTau21SF = 1.;
  int lepton_veto_size = 0;
  double lepton1_eta = 0., lepton2_eta = 0.;

  bool passPtCut = false;
  Particle ZCand, Ztemp;
  Particle WCand, lljjHigh, l1jjHigh, l2jjHigh;  // High Mass SR1
  //Particle lljjLow, l1jjLow, l2jjLow;          // Low Mass SR1
  Particle llj, l1j, l2j;                        // Low Mass SR2
  Particle l1J, l2J, llJ;                        // High Mass SR2
  
  /*if(muons.size()==2 && electrons.size()==0){
    FillHist("Pt_muon1", muons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_muon2", muons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon1", muons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_muon2", muons.at(1).PtCone(), weight, 1000, 0., 1000.);
  }
  if(muons.size()==0 && electrons.size()==2){
    FillHist("Pt_electron1", electrons.at(0).Pt(), weight, 1000, 0., 1000.);
    FillHist("Pt_electron2", electrons.at(1).Pt(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron1", electrons.at(0).PtCone(), weight, 1000, 0., 1000.);
    FillHist("PtCone_electron2", electrons.at(1).PtCone(), weight, 1000, 0., 1000.);
  }

  if(electrons.size() > 0) cout << electrons.at(0).PtCone() << endl;*/

  //==== Set up pTcone if RunFake=true
  double mu_tight_iso = 0.07, el_tight_iso = 0.;
  double this_ptcone_muon = 0., this_ptcone_electron = 0.;

  if(RunFake){
    if(muons.size()+electrons.size() == 2){

      for(unsigned int i=0; i<muons.size(); i++){
        this_ptcone_muon = muons.at(i).CalcPtCone(muons.at(i).RelIso(), mu_tight_iso);
        muons.at(i).SetPtCone(this_ptcone_muon);
      }

      for(unsigned int i=0; i<electrons.size(); i++){

        if(param.Electron_Tight_ID.Contains("HNTight")){ // POG cut-based tight WP
          el_tight_iso = 0.0287+0.506/electrons.at(i).UncorrPt();
          if(fabs(electrons.at(i).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons.at(i).UncorrPt();
        }

        //if(param.Electron_Tight_ID.Contains("HNMVA")) el_tight_iso = 0.08;

        this_ptcone_electron = electrons.at(i).CalcPtCone(electrons.at(i).RelIso(), el_tight_iso);
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

  //==== Shift electron energy and MET if RunCF=true
  if(RunCF){
    if(muons.size()==0 && electrons.size()==2){

      electrons_beforeShift.push_back(electrons.at(0));
      electrons_beforeShift.push_back(electrons.at(1));
      electrons = ShiftElectronEnergy(param.Muon_Tight_ID, electrons, true);
      electrons_afterShift.push_back(electrons.at(0));
      electrons_afterShift.push_back(electrons.at(1));
      METv = UpdateMETElectronCF(METv, electrons_beforeShift, electrons_afterShift);

    }
  }

  //==== Define leptons
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));

  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  //==== Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

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

  trigger_lumi = 1.;

  //==== Cutflow : passing dilepton triggers
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

  weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb;
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(systL1);
    weight *= GetPileUpWeight(nPileUp, systPU);
  }

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
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
  }
  if(passMuMu){
    if(!IsDATA || isDoubleMuon){
      if(!IsDATA) trigger_lumi = dimu_trig_weight;
      FillHist(systName+"/dimu_FakeCR1_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu_FakeCR2_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(passEE){
    if(!IsDATA || isDoubleEG){
      if(!IsDATA) trigger_lumi = diel_trig_weight;
      FillHist(systName+"/diel_FakeCR1_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel_FakeCR1_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel_FakeCR2_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel_FakeCR2_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(passEMu){
    if(!IsDATA || isMuonEG){
      if(!IsDATA) trigger_lumi = emu_trig_weight;
      FillHist(systName+"/emu_FakeCR1_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu_FakeCR1_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu_FakeCR2_Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu_FakeCR2_Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }

  //========================================================
  //==== Preselection
  //========================================================

  if(leptons.size() == 2){ 

    //==== Cutflow : 2 tight leptons (truth-matched, pT > trigger thresholds)

    passPtCut = false;

    //==== Pass triggers, pT cuts
    if(muons.size()==2 && electrons.size()==0){
      if(!passMuMu) return;
      if(!IsDATA) trigger_lumi = dimu_trig_weight;
      if(IsDATA){ if(!isDoubleMuon) return; }
      if(muons.at(0).Pt()>MuonPtCut1 && muons.at(1).Pt()>MuonPtCut2) passPtCut = true;
      channel = "dimu";
    }

    if(muons.size()==0 && electrons.size()==2){
      if(!passEE) return;
      if(!IsDATA) trigger_lumi = diel_trig_weight;
      if(IsDATA){ if(!isDoubleEG) return; }
      if(electrons.at(0).Pt()>ElectronPtCut1 && electrons.at(1).Pt()>ElectronPtCut2) passPtCut = true;
      channel = "diel";
    }

    if(muons.size()==1 && electrons.size()==1){
      if(!passEMu) return;
      if(!IsDATA) trigger_lumi = emu_trig_weight;
      if(IsDATA){ if(!isMuonEG) return; }
      if(passEMu8){
        if(electrons.at(0).Pt()>EMuPtCut1 && muons.at(0).Pt()>EMuPtCut2) passPtCut = true;
      }
      if(passEMu23){
        if(muons.at(0).Pt()>EMuPtCut1 && electrons.at(0).Pt()>EMuPtCut2) passPtCut = true;
      }
      channel = "emu";
    }
    
    if(!passPtCut) return;

    //==== Truth matching
    muons_prompt.clear();
    electrons_prompt.clear();
    muons_prompt = MuonPromptOnlyHNtypeI(muons, gens);
    electrons_prompt = ElectronPromptOnlyHNtypeI(electrons, gens);

    if(channel=="dimu"){
      if(!(muons_prompt.size()==2 && electrons_prompt.size()==0)) return;
    }
    if(channel=="diel"){
      if(!(muons_prompt.size()==0 && electrons_prompt.size()==2)) return;
    }
    if(channel=="emu"){
      if(!(muons_prompt.size()==1 && electrons_prompt.size()==1)) return;
    }

    //==== Event weights for MC
    weight = 1., muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1., fatjetTau21SF = 1.;

    if(!IsDATA){

      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(systL1);
      weight *= GetPileUpWeight(nPileUp, systPU);

      for(unsigned int i=0; i<muons.size(); i++){

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
      if(channel=="emu")  triggerSF = mcCorr->EMuTrigger_SF_HNtypeI(param.Muon_Tight_ID, param.Electron_Tight_ID, muons, electrons, RunFake, systEMuTrigger);

      weight *= triggerSF;

      for(unsigned int i=0; i<fatjets.size(); i++){

        fatjetTau21SF = mcCorr->FatJetWtagSF(param.FatJet_ID, systTau21);

        weight *= fatjetTau21SF;

      }

    }

    //==== Event weights for fake, charge flip
    if(RunFake) weight *= fakeEst->GetWeight(leptons, param);
    if(RunCF) weight *= GetCFWeight(param.Muon_Tight_ID, leptons, true, 0);
    //if(RunCF) weight *= GetCFWeight2D("HNTightV1", leptons);

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

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

    //==== Cutflow : same-sign (oppsite-sign when RunCF=true) 

    if(RunCF || RunOS){ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) return; }
    else{ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) return; }

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

    //==== Cutflow : veto 3rd leptons using veto ID

    if(lepton_veto_size > 0) return;

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

    //==== Cutflow : m(ll) > 10 GeV (|m(ll)-m(Z)| > 10 GeV in OS/SS(only ee) channel)

    if(!(ZCand.M() > mllCut1)) return;
    if(RunOS){
      if(IsOnZ(ZCand.M(), mllCut2)) return;
    }
    else{
      if(!RunOS && channel=="diel" && IsOnZ(ZCand.M(), 10.)) return;
    }

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);     

    /*FillHist(systName+"/"+channel+"/PreNoJetCut/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/PreNoJetCut/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/PreNoJetCut/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/PreNoJetCut/METPhi_"+IDsuffix, METPhi, weight, 64, -3.2, 3.2);
    FillHist(systName+"/"+channel+"/PreNoJetCut/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);*/

    //==== Non-prompt CR1 : same-sign 2 leptons + b jets
    if(Nbjet_medium > 0){

      FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_FakeCR1_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_FakeCR1_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_FakeCR1_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
      FillHist(systName+"/"+channel+"_FakeCR1_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
      FillHist(systName+"/"+channel+"_FakeCR1_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
      FillHist(systName+"/"+channel+"_FakeCR1_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
      FillHist(systName+"/"+channel+"_FakeCR1_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_FakeCR1_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_FakeCR1_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
      FillHist(systName+"/"+channel+"_FakeCR1_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
      FillHist(systName+"/"+channel+"_FakeCR1_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
      FillHist(systName+"/"+channel+"_FakeCR1_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
      FillHist(systName+"/"+channel+"_FakeCR1_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

      if(RunFake){
        FillHist(systName+"/"+channel+"_FakeCR1_"+tight_leptons+"_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_FakeCR1_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
      }

      if(leptons.at(1)->Pt() > 15.){

        FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_Vertices_Pt15_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_Jets_Pt15_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_BJets_Loose_Pt15_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_BJets_Medium_Pt15_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR1_Number_FatJets_Pt15_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR1_ZCand_Mass_Pt15_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR1_ZCand_Pt_Pt15_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_FakeCR1_ZCand_DeltaR_Pt15_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_FakeCR1_ZCand_DeltaPhi_Pt15_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_FakeCR1_ZCand_PtDiff_Pt15_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_FakeCR1_Lep1_Pt_Pt15_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR1_Lep2_Pt_Pt15_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR1_Lep1_Eta_Pt15_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_FakeCR1_Lep2_Eta_Pt15_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_FakeCR1_MET_Pt15_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_FakeCR1_METPhi_Pt15_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_FakeCR1_MET2ST_Pt15_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

        if(RunFake){
          FillHist(systName+"/"+channel+"_FakeCR1_"+tight_leptons+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_FakeCR1_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        }

      }

    }

    //==== Non-prompt CR2 : no jets && same-sign back-to-back 2 leptons
    if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){
     
      // Cutflow : jet requirement for non-prompt CR2 
      FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

      if(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5){

        FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR2_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_FakeCR2_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR2_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_FakeCR2_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_FakeCR2_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_FakeCR2_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_FakeCR2_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR2_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_FakeCR2_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_FakeCR2_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_FakeCR2_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_FakeCR2_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_FakeCR2_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

        if(RunFake){
          FillHist(systName+"/"+channel+"_FakeCR2_"+tight_leptons+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_FakeCR2_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        }

        if(leptons.at(1)->Pt() > 15.){
        
          FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_Vertices_Pt15_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_Jets_Pt15_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_BJets_Loose_Pt15_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_BJets_Medium_Pt15_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_FakeCR2_Number_FatJets_Pt15_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_FakeCR2_ZCand_Mass_Pt15_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_FakeCR2_ZCand_Pt_Pt15_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_FakeCR2_ZCand_DeltaR_Pt15_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_FakeCR2_ZCand_DeltaPhi_Pt15_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_FakeCR2_ZCand_PtDiff_Pt15_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"/"+channel+"_FakeCR2_Lep1_Pt_Pt15_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_FakeCR2_Lep2_Pt_Pt15_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_FakeCR2_Lep1_Eta_Pt15_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_FakeCR2_Lep2_Eta_Pt15_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_FakeCR2_MET_Pt15_"+IDsuffix, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_FakeCR2_METPhi_Pt15_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_FakeCR2_MET2ST_Pt15_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

          if(RunFake){
            FillHist(systName+"/"+channel+"_FakeCR2_"+tight_leptons+"_Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_FakeCR2_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
          }

        }

      }

    }

    //==== Cutflow : jet requirement
    if(!(jets.size()>1 || fatjets.size()>0)) return; 
    //if(!(fatjets.size()>0) && !(jets.size()>1 && fatjets.size()==0) && !(jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.)) continue;  // Preselection jet cut in EXO-17-028

    FillHist(systName+"/"+channel+"_Pre_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_Pre_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"_Pre_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
    FillHist(systName+"/"+channel+"_Pre_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"_Pre_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"_Pre_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"_Pre_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"_Pre_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
    FillHist(systName+"/"+channel+"_Pre_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
    FillHist(systName+"/"+channel+"_Pre_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
    FillHist(systName+"/"+channel+"_Pre_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
    FillHist(systName+"/"+channel+"_Pre_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
    FillHist(systName+"/"+channel+"_Pre_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
    FillHist(systName+"/"+channel+"_Pre_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
    FillHist(systName+"/"+channel+"_Pre_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
    FillHist(systName+"/"+channel+"_Pre_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
    FillHist(systName+"/"+channel+"_Pre_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
    FillHist(systName+"/"+channel+"_Pre_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
    FillHist(systName+"/"+channel+"_Pre_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

    if(RunFake){
      FillHist(systName+"/"+channel+"_Pre_"+tight_leptons+"_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_Pre_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);
    }

    if(leptons.at(1)->Pt() > 15.){

      FillHist(systName+"/"+channel+"_Pre_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_Pre_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_Pre_Number_Vertices_Pt15_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
      FillHist(systName+"/"+channel+"_Pre_Number_Jets_Pt15_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_Pre_Number_BJets_Loose_Pt15_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_Pre_Number_BJets_Medium_Pt15_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_Pre_Number_FatJets_Pt15_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
      FillHist(systName+"/"+channel+"_Pre_ZCand_Mass_Pt15_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_Pre_ZCand_Pt_Pt15_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
      FillHist(systName+"/"+channel+"_Pre_ZCand_DeltaR_Pt15_"+IDsuffix, dRll, weight, 60, 0., 6.);
      FillHist(systName+"/"+channel+"_Pre_ZCand_DeltaPhi_Pt15_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
      FillHist(systName+"/"+channel+"_Pre_ZCand_PtDiff_Pt15_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
      FillHist(systName+"/"+channel+"_Pre_Lep1_Pt_Pt15_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_Pre_Lep2_Pt_Pt15_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
      FillHist(systName+"/"+channel+"_Pre_Lep1_Eta_Pt15_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
      FillHist(systName+"/"+channel+"_Pre_Lep2_Eta_Pt15_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
      FillHist(systName+"/"+channel+"_Pre_MET_Pt15_"+IDsuffix, MET, weight, 2000, 0., 2000.);
      FillHist(systName+"/"+channel+"_Pre_METPhi_Pt15_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
      FillHist(systName+"/"+channel+"_Pre_MET2ST_Pt15_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);

      if(RunFake){
        FillHist(systName+"/"+channel+"_Pre_"+tight_leptons+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_Pre_"+tight_leptons+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
      }

    }

    //==== Event selections for each CR
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      //jets_WCandLowMass.clear();
      jets_WCandHighMass.clear();

      //==== This is the number or events at preselection
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

      //==== Non-prompt CR0 : Preselection + b jets
      /*if(it_rg == 0){

        if(!(Nbjet_medium > 0)) continue;

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDsuffix, METPhi, weight, 64, -3.2, 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      }*/
    
      //==== Low mass SR1, CR1 & High mass SR1, CR1
      if(it_rg < 4){

        if(it_rg < 2) continue; // No Low mass for Run2        

        if(!(jets.size()>=2 && fatjets.size()==0)) continue;

        //jets_WCandLowMass  = JetsWCandLowMass(*leptons.at(0), *leptons.at(1), jets, MW);
        jets_WCandHighMass = JetsWCandHighMass(jets, MW);

        WCand    = jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        //lljjLow  = *leptons.at(0) + *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        //l1jjLow  = *leptons.at(0) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        //l2jjLow  = *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        lljjHigh = *leptons.at(0) + *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        l1jjHigh = *leptons.at(0) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        l2jjHigh = *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        dRl1jj   = leptons.at(0)->DeltaR(WCand);
        dRl2jj   = leptons.at(1)->DeltaR(WCand);
        dRjj     = jets_WCandHighMass.at(0).DeltaR(jets_WCandHighMass.at(1));

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_nocut_"+IDsuffix, Nvtx, weight, 100, 0., 100.); 
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_nocut_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_nocut_"+IDsuffix, WCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_nocut_"+IDsuffix, dRjj, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nocut_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nocut_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nocut_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nocut_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_nocut_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_nocut_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_nocut_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjLow_Mass_nocut_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjLow_Mass_nocut_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjLow_Mass_nocut_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjHigh_Mass_nocut_"+IDsuffix, lljjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_Mass_nocut_"+IDsuffix, l1jjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_Mass_nocut_"+IDsuffix, l2jjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_DeltaR_nocut_"+IDsuffix, dRl1jj, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_DeltaR_nocut_"+IDsuffix, dRl2jj, weight, 60, 0., 6.);

        if(systName=="Central" && Nbjet_medium==0){
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_nobjet_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_nobjet_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_nobjet_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nobjet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_nobjet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nobjet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nobjet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_nobjet_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_nobjet_"+IDsuffix, WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_nobjet_"+IDsuffix, dRjj, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_nobjet_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_nobjet_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nobjet_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nobjet_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nobjet_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nobjet_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nobjet_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nobjet_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nobjet_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nobjet_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nobjet_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_nobjet_"+IDsuffix, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_nobjet_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_nobjet_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjHigh_Mass_nobjet_"+IDsuffix, lljjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_Mass_nobjet_"+IDsuffix, l1jjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_Mass_nobjet_"+IDsuffix, l2jjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_DeltaR_nobjet_"+IDsuffix, dRl1jj, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_DeltaR_nobjet_"+IDsuffix, dRl2jj, weight, 60, 0., 6.);
        }

        //==== Low mass SR1 
        /*if(it_rg == 0){
          if(!(Nbjet_medium == 0)) continue;
          if(!(lljjLow.M() < 300.)) continue;
          if(!(MET < 80.)) continue;
        }*/

        //==== Low mass CR1
        /*if(it_rg == 1){
          if(!(lljjLow.M() < 300.)) continue;
          if(!(Nbjet_medium>0 || MET>100.)) continue;
        }*/

        //==== High mass SR1
        if(it_rg == 2){
          if(!(Nbjet_medium == 0)) continue;
          if(!(WCand.M()>30. && WCand.M()<150.)) continue;
          if(!(MET2ST < 15.)) continue;
        }

        //==== High mass CR1
        if(it_rg == 3){
          if(!(WCand.M()>30. && WCand.M()<150.)) continue;
          if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
        }

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_"+IDsuffix, WCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_"+IDsuffix, dRjj, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjLow_Mass_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjLow_Mass_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
        //FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjLow_Mass_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjHigh_Mass_"+IDsuffix, lljjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_Mass_"+IDsuffix, l1jjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_Mass_"+IDsuffix, l2jjHigh.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_DeltaR_"+IDsuffix, dRl1jj, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_DeltaR_"+IDsuffix, dRl2jj, weight, 60, 0., 6.);

        if(leptons.at(1)->Pt() > 15.){
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_Pt15_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_Pt15_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_Pt15_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_Pt15_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_Pt15_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_Pt15_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_Pt15_"+IDsuffix, WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_Pt15_"+IDsuffix, dRjj, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_Pt15_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_Pt15_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_Pt15_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_Pt15_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_Pt15_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_Pt15_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_Pt15_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_Pt15_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_Pt15_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_Pt15_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_Pt15_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_Pt15_"+IDsuffix, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_Pt15_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_Pt15_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjHigh_Mass_Pt15_"+IDsuffix, lljjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_Mass_Pt15_"+IDsuffix, l1jjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_Mass_Pt15_"+IDsuffix, l2jjHigh.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_DeltaR_Pt15_"+IDsuffix, dRl1jj, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_DeltaR_Pt15_"+IDsuffix, dRl2jj, weight, 60, 0., 6.);
        }

        if(channel=="diel"){
          if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoEC_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoEC_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_NoEC_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_NoEC_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoEC_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Mass_NoEC_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_Pt_NoEC_"+IDsuffix, WCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCand_DeltaR_NoEC_"+IDsuffix, dRjj, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet1_Pt_NoEC_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_WCandJet2_Pt_NoEC_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoEC_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoEC_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoEC_"+IDsuffix, dRll, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoEC_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoEC_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoEC_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoEC_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoEC_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoEC_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_NoEC_"+IDsuffix, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_NoEC_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoEC_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_lljjHigh_Mass_NoEC_"+IDsuffix, lljjHigh.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_Mass_NoEC_"+IDsuffix, l1jjHigh.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_Mass_NoEC_"+IDsuffix, l2jjHigh.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1jjHigh_DeltaR_NoEC_"+IDsuffix, dRl1jj, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2jjHigh_DeltaR_NoEC_"+IDsuffix, dRl2jj, weight, 60, 0., 6.);
          }
        }

        //==== Apply optimized cuts for each mass point and count the number of events
        if(!(it_rg == 2)) continue;

        Lep1PtCutSR1.clear();
        Lep2PtCutSR1.clear();
        mlljjCut.clear();
        mljjCut1.clear();
        mljjCut2.clear();
        MET2STCut.clear();

        if(channel=="dimu"){  // For these cuts, see page 22-24 of arXiv:1806.10905
          Lep1PtCutSR1 = {25.,  50.,  100., 110., 110., 110., 110., 110., 110.,  110.,  110.,  110.,  110.,  110.,  110.,  110.};
          Lep2PtCutSR1 = {15.,  40.,  50.,  60.,  60.,  10.,  10.,  10.,  10.,   10.,   10.,   10.,   10.,   10.,   10.,   10.};
          mlljjCut     = {110., 250., 370., 490., 610., 680., 800., 800., 800.,  800.,  800.,  800.,  800.,  800.,  800.,  800.};
          mljjCut1     = {55.,  160., 225., 295., 370., 370., 370., 370., 370.,  370.,  370.,  370.,  370.,  370.,  370.,  370.};
          mljjCut2     = {115., 215., 340., 490., 550., 630., 885., 890., 1225., 1230., 1245., 1690., 1890., 2220., 2700., 3200.};
          MET2STCut    = {9.,   7.,   7.,   7.,   7.,   7.,   7.,   7.,   7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.};
        }
        if(channel=="diel"){
          Lep1PtCutSR1 = {25.,  55.,  80.,  100., 125., 125., 125., 125.,  125.,  125.,  125.,  125.,  125.,  125.,  125.,  125.};
          Lep2PtCutSR1 = {15.,  40.,  60.,  65.,  65.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
          mlljjCut     = {120., 220., 370., 450., 560., 760., 760., 760.,  760.,  760.,  760.,  760.,  760.,  760.,  760.,  760.};
          mljjCut1     = {50.,  160., 235., 335., 400., 400., 400., 400.,  400.,  400.,  400.,  400.,  400.,  400.,  400.,  400.};
          mljjCut2     = {110., 225., 335., 450., 555., 690., 955., 1130., 1300., 1490., 1490., 1600., 1930., 1930., 2400., 2900.};
          MET2STCut    = {6.,   6.,   6.,   6.,   6.,   6.,   6.,   6.,    6.,    6.,    6.,    6.,    6.,    6.,    6.,    6.};
        }
        if(channel=="emu"){
          Lep1PtCutSR1 = {25.,  65.,  95.,  120., 150., 175., 180.,  180.,  185.,  185.,  185.,  185.,  185.,  185.,  185.,  185.};
          Lep2PtCutSR1 = {20.,  35.,  60.,  60.,  60.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
          mlljjCut     = {110., 270., 340., 530., 580., 670., 720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.,  720.};
          mljjCut1     = {60.,  170., 255., 325., 315., 315., 350.,  400.,  450.,  500.,  550.,  600.,  650.,  650.,  650.,  650.};
          mljjCut2     = {115., 230., 325., 450., 530., 740., 1030., 1030., 1040., 1415., 1640., 1780., 1880., 1885., 2400., 2900.};
          MET2STCut    = {7.,   7.,   7.,   7.,   7.,   7.,   7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.,    7.};
        }

        for(unsigned int it_m=0; it_m<mass.size(); it_m++){
          if(!(jets.size() < 4)) continue;
          if(!(jets_WCandHighMass.at(0).Pt() > 25.)) continue;
          if(!(leptons.at(0)->Pt() > Lep1PtCutSR1.at(it_m))) continue;
          if(!(leptons.at(1)->Pt() > Lep2PtCutSR1.at(it_m))) continue;
          if(!(WCand.M()>50. && WCand.M()<120.)) continue;
          if(!(lljjHigh.M() > mlljjCut.at(it_m))) continue;
          if(it_m < 2){
            if(!(dRl2jj < 3.1)) continue;
            if(!(l2jjHigh.M()>mljjCut1.at(it_m) && l2jjHigh.M()<mljjCut2.at(it_m))) continue;
          }
          else{
            if(!(l1jjHigh.M()>mljjCut1.at(it_m) && l1jjHigh.M()<mljjCut2.at(it_m)) && !(l2jjHigh.M()>mljjCut1.at(it_m) && l2jjHigh.M()<mljjCut2.at(it_m))) continue;
          }
          if(!(MET2ST < MET2STCut.at(it_m))) continue;

          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
          if(leptons.at(1)->Pt() > 15.){
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Pt15_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Pt15_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
          }
          if(channel=="diel"){
            if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
              FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
              FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
            }
          }
        } 

        /*if(it_rg == 3){
          if(!(jets.size() < 4)) continue;
          if(!(jets.at(0).Pt() > 25.)) continue;
          if(!(leptons.at(0)->Pt() > 110.)) continue;
          if(!(WCand.M()>50. && WCand.M()<120.)) continue;
          if(!(lljjHigh.M() > 800.)) continue;
          if(!(l1jjHigh.M() > 370.)) continue;
          if(!(MET2ST < 7.)) continue;

          if(MCSample.Contains("M700")){
            if(!(l1jjHigh.M() < 885.)) continue;
          }
          if(MCSample.Contains("M1000")){
            if(!(l1jjHigh.M() < 1230.)) continue;
          }
          if(MCSample.Contains("M1500")){
            if(!(l1jjHigh.M() < 2220.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        }*/

        /*if(it_rg==1) cout << "In LowSR1, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;
        if(it_rg==3) cout << "In HighSR1, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;
        if(it_ch==0 && it_rg==3){
          if(run==276283 && event==1252562683){
            cout << "Muon1 : (" << muons.at(0).Pt() << ", " << muons.at(0).Eta() << ", " << muons.at(0).Phi() << ", " << muons.at(0).E() << ")" << endl;
            cout << "Muon2 : (" << muons.at(1).Pt() << ", " << muons.at(1).Eta() << ", " << muons.at(1).Phi() << ", " << muons.at(1).E() << ")" << endl;
            for(unsigned int i=0; i<jets.size(); i++){
              cout << TString::Itoa(i, 10)+"th Jets : (" << jets.at(i).Pt() << ", " << jets.at(i).Eta() << ", " << jets.at(i).Phi() << ", " << jets.at(i).E() << ")" << endl;
            }
          }
        }*/

      }

      //==== Low mass SR2, CR2
      if(it_rg>=4 && it_rg<6){

        continue; // No Low mass SR for Run2

        if(!(jets.size()==1 && fatjets.size()==0)) continue;

        llj = *leptons.at(0) + *leptons.at(1) + jets.at(0);
        l1j = *leptons.at(0) + jets.at(0);
        l2j = *leptons.at(1) + jets.at(0);

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_nocut_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 21000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nocut_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nocut_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nocut_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nocut_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_nocut_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_nocut_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_nocut_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llj_Mass_nocut_"+IDsuffix, llj.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1j_Mass_nocut_"+IDsuffix, l1j.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2j_Mass_nocut_"+IDsuffix, l2j.M(), weight, 4000, 0., 4000.);

        //==== Low mass SR2
        if(it_rg == 4){
          if(!(Nbjet_medium == 0)) continue;
          if(!(llj.M() < 300.)) continue;
          if(!(MET < 80.)) continue;
        }

        //==== Low mass CR2
        if(it_rg == 5){
          if(!(llj.M() < 300.)) continue;
          if(!(Nbjet_medium>0 || MET>100.)) continue;
        }

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llj_Mass_"+IDsuffix, llj.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1j_Mass_"+IDsuffix, l1j.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2j_Mass_"+IDsuffix, l2j.M(), weight, 4000, 0., 4000.);

        //if(it_rg==5) cout << "In LowSR2, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;

      }

      //==== High mass SR2, CR2
      if(it_rg >= 6){

        if(!(fatjets.size() > 0)) continue;

        fatjets_WCand = FatJetWCand(fatjets, MW);

        l1J   = *leptons.at(0) + fatjets_WCand;
        l2J   = *leptons.at(1) + fatjets_WCand;
        llJ   = *leptons.at(0) + *leptons.at(1) + fatjets_WCand;
        dRl1J = leptons.at(0)->DeltaR(fatjets_WCand);
        dRl2J = leptons.at(1)->DeltaR(fatjets_WCand);

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_nocut_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nocut_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nocut_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nocut_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nocut_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_nocut_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_nocut_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_nocut_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_nocut_"+IDsuffix, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_nocut_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_nocut_"+IDsuffix, l1J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_nocut_"+IDsuffix, l2J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_nocut_"+IDsuffix, dRl1J, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_nocut_"+IDsuffix, dRl2J, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_nocut_"+IDsuffix, llJ.M(), weight, 4000, 0., 4000.);

        if(systName=="Central" && Nbjet_medium == 0){

          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_nobjet_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_nobjet_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_nobjet_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_nobjet_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_nobjet_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_nobjet_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_nobjet_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nobjet_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_nobjet_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_nobjet_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_nobjet_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_nobjet_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_nobjet_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_nobjet_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_nobjet_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_nobjet_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_nobjet_"+IDsuffix, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_nobjet_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_nobjet_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_nobjet_"+IDsuffix, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_nobjet_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_nobjet_"+IDsuffix, l1J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_nobjet_"+IDsuffix, l2J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_nobjet_"+IDsuffix, dRl1J, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_nobjet_"+IDsuffix, dRl2J, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_nobjet_"+IDsuffix, llJ.M(), weight, 4000, 0., 4000.);

          if(RunOS && MET2ST<15.){
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_nobjet_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_nobjet_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_nobjet_MET2ST15_"+IDsuffix, ZCand.M(), weight, 4000, 0., 4000.);
          }

        }

        if(RunOS){
          if(!(ZCand.M() > 120.)) continue;
        }

        //==== High mass SR2
        if(it_rg == 6){
          if(!(Nbjet_medium == 0)) continue;
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(MET2ST < 15.)) continue;
        }

        //==== High mass CR2
        if(it_rg == 7){
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
        }

        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_"+IDsuffix, MET, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_"+IDsuffix, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_"+IDsuffix, l1J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_"+IDsuffix, l2J.M(), weight, 4000, 0., 4000.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_"+IDsuffix, dRl1J, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_"+IDsuffix, dRl2J, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_"+IDsuffix, llJ.M(), weight, 4000, 0., 4000.);

        if(leptons.at(1)->Pt() > 15.){
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 10.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 10.5, 1., cutflow_bin, 0., cutflow_max);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_Pt15_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_Pt15_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_Pt15_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_Pt15_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_Pt15_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_Pt15_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_Pt15_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_Pt15_"+IDsuffix, dRll, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_Pt15_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_Pt15_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_Pt15_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_Pt15_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_Pt15_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_Pt15_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_Pt15_"+IDsuffix, MET, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_Pt15_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_Pt15_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_Pt15_"+IDsuffix, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_Pt15_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_Pt15_"+IDsuffix, l1J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_Pt15_"+IDsuffix, l2J.M(), weight, 4000, 0., 4000.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_Pt15_"+IDsuffix, dRl1J, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_Pt15_"+IDsuffix, dRl2J, weight, 60, 0., 6.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_Pt15_"+IDsuffix, llJ.M(), weight, 4000, 0., 4000.);
        }

        if(channel=="diel"){
          if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_"+IDsuffix, 11.5, weight, cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Events_unweighted_"+IDsuffix, 11.5, 1., cutflow_bin, 0., cutflow_max);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Vertices_NoEC_"+IDsuffix, Nvtx, weight, 100, 0., 100.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_Jets_NoEC_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Loose_NoEC_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_BJets_Medium_NoEC_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Number_FatJets_NoEC_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Mass_NoEC_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_Pt_NoEC_"+IDsuffix, ZCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaR_NoEC_"+IDsuffix, dRll, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_DeltaPhi_NoEC_"+IDsuffix, dPhill, weight, 32, 0., 3.2);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_ZCand_PtDiff_NoEC_"+IDsuffix, PtDiff, weight, 100, 0., 1.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Pt_NoEC_"+IDsuffix, leptons.at(0)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Pt_NoEC_"+IDsuffix, leptons.at(1)->Pt(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep1_Eta_NoEC_"+IDsuffix, lepton1_eta, weight, 60, -3., 3.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Lep2_Eta_NoEC_"+IDsuffix, lepton2_eta, weight, 60, -3., 3.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET_NoEC_"+IDsuffix, MET, weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_METPhi_NoEC_"+IDsuffix, METPhi, weight, 32, 0., 3.2);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_MET2ST_NoEC_"+IDsuffix, MET2ST, weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Pt_NoEC_"+IDsuffix, fatjets_WCand.Pt(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_Fatjet_Mass_NoEC_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_Mass_NoEC_"+IDsuffix, l1J.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_Mass_NoEC_"+IDsuffix, l2J.M(), weight, 4000, 0., 4000.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l1J_DeltaR_NoEC_"+IDsuffix, dRl1J, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_l2J_DeltaR_NoEC_"+IDsuffix, dRl2J, weight, 60, 0., 6.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_llJ_Mass_NoEC_"+IDsuffix, llJ.M(), weight, 4000, 0., 4000.);
          }
        }

        //==== Apply optimized cuts for each mass point and count the number of events
        if(!(it_rg == 6)) continue;

        Lep1PtCutSR2.clear();
        Lep2PtCutSR2.clear();
        mlJCut1.clear();
        mlJCut2.clear();

        if(channel=="dimu"){  // For these cuts, see page 22-24 of arXiv:1806.10905
          Lep1PtCutSR2 = {25.,  100., 140., 140., 140., 140., 140., 140., 140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.};
          Lep2PtCutSR2 = {15.,  20.,  40.,  65.,  65.,  10.,  10.,  10.,  10.,   10.,   10.,   10.,   10.,   10.,   10.,   10.};
          mlJCut1      = {98.,  175., 280., 340., 445., 560., 635., 755., 840.,  900.,  990.,  1035., 1100., 1330., 1700., 2200.};
          mlJCut2      = {145., 235., 340., 445., 560., 685., 825., 960., 1055., 1205., 1250., 1430., 1595., 1800., 2300., 2800.};
        }
        if(channel=="diel"){
          Lep1PtCutSR2 = {25.,  100., 100., 100., 120., 120., 140., 140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.};
          Lep2PtCutSR2 = {15.,  20.,  30.,  35.,  35.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
          mlJCut1      = {100., 173., 270., 330., 440., 565., 635., 740.,  865.,  890.,  1035., 1085., 1140., 1300., 1700., 2200.};
          mlJCut2      = {220., 220., 330., 440., 565., 675., 775., 1005., 1030., 1185., 1395., 1460., 1590., 1800., 2300., 2800.};
        }
        if(channel=="emu"){
          Lep1PtCutSR2 = {30.,  70.,  95.,  125., 145., 160., 170., 170., 180.,  180.,  180.,  180.,  180.,  180.,  180.,  180.};
          Lep2PtCutSR2 = {15.,  30.,  55.,  55.,  60.,  15.,  15.,  15.,  15.,   15.,   15.,   15.,   15.,   15.,   15.,   15.};
          mlJCut1      = {100., 180., 280., 340., 460., 555., 610., 730., 845.,  930.,  1020., 1080., 1155., 1345., 1800., 2300.};
          mlJCut2      = {335., 225., 340., 475., 555., 645., 780., 895., 1015., 1075., 1340., 1340., 1595., 1615., 2200., 2700.};
        }

        for(unsigned int it_m=0; it_m<mass.size(); it_m++){
          if(!(leptons.at(0)->Pt() > Lep1PtCutSR2.at(it_m))) continue;
          if(!(leptons.at(1)->Pt() > Lep2PtCutSR2.at(it_m))) continue;
          if(it_m < 2){ 
            if(!(l2J.M()>mlJCut1.at(it_m) && l2J.M()<mlJCut2.at(it_m))) continue;
          }
          else{
            if(!(l1J.M()>mlJCut1.at(it_m) && l1J.M()<mlJCut2.at(it_m)) && !(l2J.M()>mlJCut1.at(it_m) && l2J.M()<mlJCut2.at(it_m))) continue;
          }

          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
          FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
          if(leptons.at(1)->Pt() > 15.){
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Pt15_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
            FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_Pt15_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
          }
          if(channel=="diel"){
            if(fabs(electrons.at(0).scEta())<1.479 && fabs(electrons.at(1).scEta())<1.479){
              FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_"+IDsuffix, 0.5, weight, 1, 0., 1.);
              FillHist(systName+"/"+channel+"_"+regions.at(it_rg)+"_M"+TString::Itoa(mass.at(it_m), 10)+"_NoEC_Number_Events_unweighted_"+IDsuffix, 0.5, 1., 1, 0., 1.);
            }
          }
        }

        /*if(it_rg == 7){
          if(!(leptons.at(0)->Pt() > 140.)) continue;
          if(!(fatjets_WCand.SDMass()>40. && fatjets_WCand.SDMass()<130.)) continue;
          if(!(MET2ST < 15.)) continue;

          if(MCSample.Contains("M700")){
            if(!(l1J.M()>635. && l1J.M()<825.)) continue;
          }
          if(MCSample.Contains("M1000")){
            if(!(l1J.M()>900. && l1J.M()<1205.)) continue;
          }
          if(MCSample.Contains("M1500")){
            if(!(l1J.M()>1330. && l1J.M()<1800.)) continue;
          }

          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
          FillHist(channels.at(it_ch)+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        }*/
        //if(it_rg==7) cout << "In HighSR2, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;

      }

    } // Region Loop 

  } // Dilepton selection

}

