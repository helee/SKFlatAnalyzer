#include "HNtypeI_SR_2016H.h"

HNtypeI_SR_2016H::HNtypeI_SR_2016H(){

}

void HNtypeI_SR_2016H::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst     = HasFlag("RunSyst");
  RunFake     = HasFlag("RunFake");
  RunCF       = HasFlag("RunCF");

  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunFake = " << RunFake << endl;
  cout << "[HNtypeI_SR_2016H::initializeAnalyzer] RunCF = " << RunCF << endl;

  /*if(RunMuon){
    MuonTightIDs     = {"ISRTightV1", "ISRTightV2", "HNTightV1", "HNTightV2"};
    MuonLooseIDs     = {"ISRLoose", "ISRLoose", "HNLooseV1", "HNLooseV1"};
    MuonVetoIDs      = {"ISRVeto", "ISRVeto", "ISRVeto", "ISRVeto"};
    ElectronTightIDs = {"HNTightV1", "HNTightV1", "HNTightV1", "HNTightV1"};  // Not used
    ElectronLooseIDs = {"HNLooseV1", "HNLooseV1", "HNLooseV1", "HNLooseV1"};  // Not used
    ElectronVetoIDs  = {"ISRVeto", "ISRVeto", "ISRVeto", "ISRVeto"};
    MuonFRNames      = {"POGCBV1", "POGCBV2", "POGCBV3", "POGCBV4"};
    ElectronFRNames  = {"POGCB", "POGCB", "POGCB", "POGCB"};
  }

  if(RunTightIP){
    MuonTightIDs     = {"HNTightV3", "HNTightV4", "HNTightV5", "HNTightV6"};
    MuonLooseIDs     = {"HNLooseV2", "HNLooseV2", "HNLooseV2", "HNLooseV2"};
    MuonVetoIDs      = {"ISRVeto", "ISRVeto", "ISRVeto", "ISRVeto"};
    ElectronTightIDs = {"HNTightV1", "HNTightV1", "HNTightV1", "HNTightV1"};  // Not used
    ElectronLooseIDs = {"HNLooseV1", "HNLooseV1", "HNLooseV1", "HNLooseV1"};  // Not used
    ElectronVetoIDs  = {"ISRVeto", "ISRVeto", "ISRVeto", "ISRVeto"};
    MuonFRNames      = {"POGCBV5", "POGCBV6", "POGCBV7", "POGCBV8"};
    ElectronFRNames  = {"POGCB", "POGCB", "POGCB", "POGCB"};
  }

  if(RunElectron){
    MuonTightIDs     = {"ISRTightV1", "ISRTightV1"};  // Not used
    MuonLooseIDs     = {"ISRLoose", "ISRLoose"};      // Not used
    MuonVetoIDs      = {"ISRVeto", "ISRVeto"};
    ElectronTightIDs = {"HNTightV1", "HNMVATight"};
    ElectronLooseIDs = {"HNLooseV1", "HNMVALoose"};
    ElectronVetoIDs  = {"ISRVeto", "HNMVAVeto"};
    MuonFRNames      = {"POGCBV1", "POGCBV1"};
    ElectronFRNames  = {"POGCB", "POGMVA"};
  }

  if(RunTightEMu){
    MuonTightIDs     = {"HNTightV5", "HNTightV5", "HNTightV6", "HNTightV6"};
    MuonLooseIDs     = {"HNLooseV2", "HNLooseV2", "HNLooseV2", "HNLooseV2"};
    MuonVetoIDs      = {"ISRVeto", "ISRVeto", "ISRVeto", "ISRVeto"};
    ElectronTightIDs = {"HNTightV1", "HNMVATight", "HNTightV1", "HNMVATight"};
    ElectronLooseIDs = {"HNLooseV1", "HNMVALoose", "HNLooseV1", "HNMVALoose"};
    ElectronVetoIDs  = {"ISRVeto", "HNMVAVeto", "ISRVeto", "HNMVAVeto"};
    MuonFRNames      = {"POGCBV7", "POGCBV7", "POGCBV8", "POGCBV8"};
    ElectronFRNames  = {"POGCB", "POGMVA", "POGCB", "POGMVA"};
  }

  if(RunMuon || RunEMu){
    MuonTightIDs     = {"HNTightV1", "HNTightV1", "HNTightV1"};
    MuonLooseIDs     = {"HNLooseV1", "HNLooseV2", "HNLooseV3"};
    MuonVetoIDs      = {"ISRVeto", "ISRVeto", "ISRVeto"};
    ElectronTightIDs = {"HNTightV1", "HNTightV1", "HNTightV1"};
    ElectronLooseIDs = {"HNLooseV1", "HNLooseV1", "HNLooseV1"};
    ElectronVetoIDs  = {"ISRVeto", "ISRVeto", "ISRVeto"};
    MuonFRNames      = {"HNV11", "HNV12", "HNV13"};
    ElectronFRNames  = {"HNV11", "HNV11", "HNV11"};
  }*/

  MuonTightIDs     = {"HNTightV1"};
  MuonLooseIDs     = {"HNLooseV3"};
  MuonVetoIDs      = {"ISRVeto"};
  ElectronTightIDs = {"HNTightV1"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"ISRVeto"};
  MuonFRNames      = {"HNV13"};
  ElectronFRNames  = {"HNV11"};

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

    // These are needed for applying lepton pT cuts 
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

    // These are needed for applying lepton pT cuts 
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

    // These are needed for applying lepton pT cuts 
    Mu8Ele23Triggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    Mu23Ele12Triggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    MuonPtCut1 = 20., MuonPtCut2 = 10.;
    ElectronPtCut1 = 25., ElectronPtCut2 = 15.;
    EMuPtCut1 = 25., EMuPtCut2 = 15.;
  }

  //cout << "[HNtypeI_SR_2016H::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  //cout << "[HNtypeI_SR_2016H::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
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
    param.Jet_ID    = "HNTight";
    if(DataYear==2016) param.FatJet_ID = "HNTight0p55";
    else param.FatJet_ID = "HNTight0p45";

    executeEventFromParameter(param);

    /*if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }*/
  }
}

void HNtypeI_SR_2016H::executeEventFromParameter(AnalyzerParameter param){

  TString IDsuffix = "HNRun2";
  /*if(RunMuon){
    if(param.Muon_Tight_ID.Contains("ISRTightV1")) IDsuffix = "POGCBV1";
    if(param.Muon_Tight_ID.Contains("ISRTightV2")) IDsuffix = "POGCBV2";
    if(param.Muon_Tight_ID.Contains("HNTightV1"))  IDsuffix = "POGCBV3";
    if(param.Muon_Tight_ID.Contains("HNTightV2"))  IDsuffix = "POGCBV4";
  }
  if(RunTightIP){
    if(param.Muon_Tight_ID.Contains("HNTightV3"))  IDsuffix = "POGCBV5";
    if(param.Muon_Tight_ID.Contains("HNTightV4"))  IDsuffix = "POGCBV6";
    if(param.Muon_Tight_ID.Contains("HNTightV5"))  IDsuffix = "POGCBV7";
    if(param.Muon_Tight_ID.Contains("HNTightV6"))  IDsuffix = "POGCBV8";
  }
  if(RunElectron || RunEMu){
    if(param.Electron_Tight_ID.Contains("HNTight")) IDsuffix = "POGCB";
    if(param.Electron_Tight_ID.Contains("HNMVA"))  IDsuffix = "POGMVA";
  }
  if(RunTightEMu){
    if(param.Muon_Tight_ID.Contains("HNTightV5") && param.Electron_Tight_ID.Contains("HNTight")) IDsuffix = "POGCBV1";
    if(param.Muon_Tight_ID.Contains("HNTightV5") && param.Electron_Tight_ID.Contains("MVA")) IDsuffix = "POGMVAV1";
    if(param.Muon_Tight_ID.Contains("HNTightV6") && param.Electron_Tight_ID.Contains("HNTight")) IDsuffix = "POGCBV2";
    if(param.Muon_Tight_ID.Contains("HNTightV6") && param.Electron_Tight_ID.Contains("MVA")) IDsuffix = "POGMVAV2";
  }

  if(RunMuon || RunEMu){
    if(param.Muon_Tight_ID.Contains("V1") && param.Muon_Loose_ID.Contains("V1")) IDsuffix = "HNV11";
    if(param.Muon_Tight_ID.Contains("V1") && param.Muon_Loose_ID.Contains("V2")) IDsuffix = "HNV12";
    if(param.Muon_Tight_ID.Contains("V1") && param.Muon_Loose_ID.Contains("V3")) IDsuffix = "HNV13";
    if(param.Muon_Tight_ID.Contains("V2") && param.Muon_Loose_ID.Contains("V2")) IDsuffix = "HNV22";
  }

  if(RunElectron){
    if(param.Electron_Tight_ID.Contains("V1") && param.Electron_Loose_ID.Contains("V1")) IDsuffix = "HNV11";
  }*/

  //vector<TString> channels = {"dimu", "diel", "emu"};
  TString channel = "";
  vector<TString> regions = {"fakeCR1", "lowSR1", "lowCR1", "highSR1", "highCR1", "lowSR2", "lowCR2", "highSR2", "highCR2"};

  TString systName = param.Name;

  double cutflow_max = 10.;
  int cutflow_bin = 10;
  double weight = 1.;
  double trigger_lumi = 1., dimu_trig_weight = 0., diel_trig_weight = 0., emu_trig_weight = 0.;
  double muon_miniaodP = 0.;

  Event ev = GetEvent();
 
  // Boolean : primary datasets
  bool isDoubleMuon = false, isDoubleEG = false, isMuonEG = false;
  if(IsDATA){
    if(DataStream.Contains("DoubleMuon")) isDoubleMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("EGamma")) isDoubleEG = true;
    if(DataStream.Contains("MuonEG")) isMuonEG = true;
  }

  // Boolean : passTrigger
  bool passMuMu  = ev.PassTrigger(MuonTriggersH);      // NOTE : Change for 2016H 
  bool passEE    = ev.PassTrigger(ElectronTriggers);
  bool passEMu   = ev.PassTrigger(EMuTriggersH);       // NOTE : Change for 2016H 
  bool passE23Mu = ev.PassTrigger(Mu8Ele23TriggersH);  // NOTE : Change for 2016H 
  bool passEMu23 = ev.PassTrigger(Mu23Ele12TriggersH); // NOTE : Change for 2016H 

  //========================================================
  //==== No Cut
  //========================================================

  weight = 1.;
  if(!IsDATA){
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
    weight *= ev.MCweight();
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);
  }

  // Cutflow : No Cuts
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    if(!IsDATA || isDoubleMuon){
      FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isDoubleEG){
      FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isMuonEG){
      FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(!IsDATA || isDoubleMuon){
    FillHist(systName+"/dimu/fakeCR2/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isDoubleEG){
    FillHist(systName+"/diel/fakeCR2/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel/fakeCR2/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isMuonEG){
    FillHist(systName+"/emu/fakeCR2/Number_Events_"+IDsuffix, 0.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 0.5, 1., cutflow_bin, 0., cutflow_max);
  }

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  // Cutflow : MET filter
  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    if(!IsDATA || isDoubleMuon){
      FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isDoubleEG){
      FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
    if(!IsDATA || isMuonEG){
      FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(!IsDATA || isDoubleMuon){
    FillHist(systName+"/dimu/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/dimu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isDoubleEG){
    FillHist(systName+"/diel/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/diel/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
  }
  if(!IsDATA || isMuonEG){
    FillHist(systName+"/emu/fakeCR2/Number_Events_"+IDsuffix, 1.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/emu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 1.5, 1., cutflow_bin, 0., cutflow_max);
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
    cout << "[HNtypeI_SR_2016H::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }*/

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  // Leptons
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

  // Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  // For charge flip
  vector<Electron> electrons_beforeShift;
  vector<Electron> electrons_afterShift;
  electrons_beforeShift.clear();
  electrons_afterShift.clear();

  // Jets
  vector<Jet> jets_nolepveto = SelectJets(this_AllJets, param.Jet_ID, 20., 2.7);  // AK4jets used for b tag
  vector<FatJet> fatjets_nolepveto = SelectFatJets(this_AllFatJets, param.FatJet_ID, 200., 2.7);

  // Jet, FatJet selection to avoid double counting due to jets matched geometrically with a lepton
  // Fatjet selection in CATanalyzer (see the links)
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/CATConfig/SelectionConfig/user_fatjets.sel
  // https://github.com/jedori0228/LQanalyzer/blob/CatAnalyzer_13TeV_v8-0-7.36_HNAnalyzer/LQCore/Selection/src/FatJetSelection.cc#L113-L124

  vector<FatJet> fatjets = FatJetsVetoLeptonInside(fatjets_nolepveto, electrons_veto, muons_veto);  // AK8jets used in SR, CR
  vector<Jet> jets_lepveto = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto);
  vector<Jet> jets_insideFatjets = JetsInsideFatJet(jets_lepveto, fatjets);  // For jets inside a fatjet, remove their smearing from MET. Because FatJet smearing is already propagted to MET.
  vector<Jet> jets_PUveto = JetsPassPileupMVA(jets_lepveto);
  vector<Jet> jets = JetsAwayFromFatJet(jets_PUveto, fatjets);  // AK4jets used in SR, CR

  vector<Jet> jets_WCandLowMass;
  vector<Jet> jets_WCandHighMass;
  FatJet fatjets_WCand;
  jets_WCandLowMass.clear();
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
  //==== B-Tagging 
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
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Loose, jets_nolepveto.at(ij))) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets_nolepveto.at(ij))) Nbjet_medium++;
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

  //========================================================
  //==== Define particles, variables
  //========================================================

  double ST = 0., MET2ST = 0.;
  double dRll = 0., dRl2WCand = 0.;
  //double MZ = 91.1876;
  double MW = 80.379;
  double muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1.;
  int lepton_veto_size = 0;

  bool passPtCut = false;
  Particle ZCand, Ztemp;
  Particle WCand, lljjLow, l1jjLow, l2jjLow, lljjHigh, l1jjHigh, l2jjHigh;  // Low Mass SR1, High Mass SR1
  Particle llj, l1j, l2j;                                                   // Low Mass SR2
  Particle l1J, l2J;                                                        // High Mass SR2
  
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

  // Set up pTcone if RunFake=true
  double mu_tight_iso = 0.05, el_tight_iso = 0.;
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

      // Correct MET if RunFake=true, because pT was replaced by pTcone
      METv = UpdateMETFake(METv, electrons, muons);

      muons = MuonUsePtCone(muons);
      electrons = ElectronUsePtCone(electrons);
      std::sort(muons.begin(), muons.end(), PtComparing);
      std::sort(electrons.begin(), electrons.end(), PtComparing);

    }
  }

  // Shift electron energy and MET if RunCF=true
  if(RunCF){
    if(muons.size()==0 && electrons.size()==2){

      electrons_beforeShift.push_back(electrons.at(0));
      electrons_beforeShift.push_back(electrons.at(1));
      electrons = ShiftElectronEnergy(electrons, param, true);
      electrons_afterShift.push_back(electrons.at(0));
      electrons_afterShift.push_back(electrons.at(1));
      METv = UpdateMETElectronCF(METv, electrons_beforeShift, electrons_afterShift);

    }
  }

  // Define leptons
  for(unsigned int i=0; i<muons.size(); i++) leptons.push_back(&muons.at(i));
  for(unsigned int i=0; i<electrons.size(); i++) leptons.push_back(&electrons.at(i));

  std::sort(leptons.begin(), leptons.end(), PtComparingPtr);

  // Define leptons passing veto IDs
  for(unsigned int i=0; i<muons_veto.size(); i++) leptons_veto.push_back(&muons_veto.at(i));
  for(unsigned int i=0; i<electrons_veto.size(); i++) leptons_veto.push_back(&electrons_veto.at(i));

  lepton_veto_size = leptons_veto.size() - leptons.size();

  // Define ST, MET^2/ST
  MET = METv.Pt();

  for(unsigned int i=0; i<jets.size(); i++) ST += jets.at(i).Pt();
  for(unsigned int i=0; i<fatjets.size(); i++) ST += fatjets.at(i).Pt();
  for(unsigned int i=0; i<leptons.size(); i++) ST += leptons.at(i)->Pt();

  ST += MET;
  MET2ST = MET*MET/ST;

  //========================================================
  //==== Event selection
  //========================================================

  trigger_lumi = 1.;

  // Cutflow : passing dilepton triggers
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
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp,0);
  }

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
    if(passMuMu){
      if(!IsDATA || isDoubleMuon){
        if(!IsDATA) trigger_lumi = dimu_trig_weight;
        FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/dimu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
    if(passEE){
      if(!IsDATA || isDoubleEG){
        if(!IsDATA) trigger_lumi = diel_trig_weight;
        FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/diel/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
    if(passEMu){
      if(!IsDATA || isMuonEG){
        if(!IsDATA) trigger_lumi = emu_trig_weight;
        FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/emu/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
      }
    }
  }
  if(passMuMu){
    if(!IsDATA || isDoubleMuon){
      if(!IsDATA) trigger_lumi = dimu_trig_weight;
      FillHist(systName+"/dimu/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/dimu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(passEE){
    if(!IsDATA || isDoubleEG){
      if(!IsDATA) trigger_lumi = diel_trig_weight;
      FillHist(systName+"/diel/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/diel/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }
  if(passEMu){
    if(!IsDATA || isMuonEG){
      if(!IsDATA) trigger_lumi = emu_trig_weight;
      FillHist(systName+"/emu/fakeCR2/Number_Events_"+IDsuffix, 2.5, weight*trigger_lumi, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/emu/fakeCR2/Number_Events_unweighted_"+IDsuffix, 2.5, 1., cutflow_bin, 0., cutflow_max);
    }
  }

  //========================================================
  //==== Preselection
  //========================================================

  if(leptons.size() == 2){ 

    // Cutflow : 2 tight leptons (truth-matched, pT > trigger thresholds)

    passPtCut = false;

    // Pass triggers, pT cuts
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
      if(passE23Mu){
        if(electrons.at(0).Pt()>EMuPtCut1 && muons.at(0).Pt()>EMuPtCut2) passPtCut = true;
      }
      if(passEMu23){
        if(muons.at(0).Pt()>EMuPtCut1 && electrons.at(0).Pt()>EMuPtCut2) passPtCut = true;;
      }
      channel = "emu";
    }
    
    if(!passPtCut) return;

    // Truth matching
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

    // Event weights for MC
    weight = 1., muonRecoSF = 1., muonIDSF = 1., muonIsoSF = 1., electronRecoSF = 1., electronIDSF = 1.;

    if(!IsDATA){

      weight *= weight_norm_1invpb*trigger_lumi;
      weight *= ev.MCweight();
      weight *= GetPrefireWeight(0);
      weight *= GetPileUpWeight(nPileUp,0);

      for(unsigned int i=0; i<muons.size(); i++){
        if(param.Muon_Tight_ID.Contains("HighPt")){
          muon_miniaodP = sqrt( muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz() );
          muonRecoSF   *= mcCorr->MuonReco_SF("HighPtMuonRecoSF", muons.at(i).Eta(), muon_miniaodP, 0);
          muonIDSF     *= mcCorr->MuonID_SF("NUM_HighPtID_DEN_genTracks",  muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
          muonIsoSF    *= mcCorr->MuonISO_SF("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut", muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0);
        }
        else{
          muonRecoSF *= 1.;
          muonIDSF   *= 1.;
          muonIsoSF  *= 1.;
        }
        weight *= muonRecoSF*muonIDSF*muonIsoSF;
      }
 
      for(unsigned int i=0; i<electrons.size(); i++){
        electronRecoSF *= mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
        if(param.Electron_Tight_ID.Contains("HEEP")){
          electronIDSF *= mcCorr->ElectronID_SF("HEEP", electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0);
        }
        else electronIDSF *= 1.;
        weight *= electronRecoSF*electronIDSF;
      }

    }

    // Event weights for fake, charge flip
    if(RunFake) weight *= fakeEst->GetWeight(leptons, param);
    //if(RunCF) weight *= GetCFweight(leptons, param, true, 0);
    if(RunCF) weight *= GetCFWeightElectron(leptons, param.Electron_Tight_ID);

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 3.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 3.5, 1., cutflow_bin, 0., cutflow_max);

    ZCand = *leptons.at(0) + *leptons.at(1);
    dRll  = leptons.at(0)->DeltaR(*leptons.at(1));

    // Cutflow : same-sign (oppsite-sign when RunCF=true) 

    if(!RunCF){ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()<0) return; }
    if(RunCF){ if(leptons.at(0)->Charge()*leptons.at(1)->Charge()>0) return; }
    /*if(MCSample.Contains("M1500")){
      if(muons.at(0).Charge() < 0) continue;
    }*/

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 4.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 4.5, 1., cutflow_bin, 0., cutflow_max);

    // Cutflow : veto 3rd leptons using veto ID

    if(lepton_veto_size > 0) return;

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 5.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 5.5, 1., cutflow_bin, 0., cutflow_max);

    // Cutflow : m(ll) > 10 GeV (|m(ll)-m(Z)| > 10 GeV in ee channel)

    if(!(ZCand.M() > 10.)) return;
    if(channel=="diel" && IsOnZ(ZCand.M(), 10.)) return;

    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);
    }
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 6.5, weight, cutflow_bin, 0., cutflow_max);
    FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 6.5, 1., cutflow_bin, 0., cutflow_max);     

    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/"+"PreNoJetCut/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

    // non-prompt CR2 : no jets && same-sign back-to-back 2 leptons
    if(jets.size()+fatjets.size()==0 && Nbjet_medium==0){
     
      // Cutflow : jet requirement for non-prompt CR2 
      FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

      if(leptons.at(0)->DeltaR(*leptons.at(1)) > 2.5){
        FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/fakeCR2/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/fakeCR2/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/fakeCR2/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/fakeCR2/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/fakeCR2/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/fakeCR2/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/fakeCR2/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/fakeCR2/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/fakeCR2/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/fakeCR2/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/fakeCR2/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/fakeCR2/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      }

    }

    // Cutflow : jet requirement
    if(!(jets.size()>0 || fatjets.size()>0)) return; 
    //if(!(fatjets.size()>0) && !(jets.size()>1 && fatjets.size()==0) && !(jets.size()==1 && fatjets.size()==0 && ZCand.M()<80.)) continue;
    
    FillHist(systName+"/"+channel+"/Pre/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/Pre/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/Pre/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/Pre/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
    FillHist(systName+"/"+channel+"/Pre/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
    FillHist(systName+"/"+channel+"/Pre/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/Pre/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/Pre/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/Pre/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/Pre/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
    FillHist(systName+"/"+channel+"/Pre/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
    FillHist(systName+"/"+channel+"/Pre/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

    /*if(muons.at(0).Charge() > 0){
      FillHist(channels.at(it_ch)+"/Pre/MET2ST_MupMup_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Pt_MupMup_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
      FillHist(channels.at(it_ch)+"/Pre/Lep2_Eta_MupMup_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
    }*/

    // Event selections for each CR
    for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){
      jets_WCandLowMass.clear();
      jets_WCandHighMass.clear();

      // This is the number or events at preselection
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 7.5, weight, cutflow_bin, 0., cutflow_max);
      FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 7.5, 1., cutflow_bin, 0., cutflow_max);

      // non-prompt CR1 : SS 2 leptons with b-tagged jets
      if(it_rg == 0){

        if(!(Nbjet_medium > 0)) continue;

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);

      }
    
      // Low mass SR1, CR1 & High mass SR1, CR1
      if(it_rg>=1 && it_rg<5){

        if(!(jets.size()>=2 && fatjets.size()==0)) continue;

        jets_WCandLowMass  = JetsWCandLowMass(*leptons.at(0), *leptons.at(1), jets, MW);
        jets_WCandHighMass = JetsWCandHighMass(jets, MW);

        WCand     = jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        lljjLow   = *leptons.at(0) + *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        l1jjLow   = *leptons.at(0) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        l2jjLow   = *leptons.at(1) + jets_WCandLowMass.at(0) + jets_WCandLowMass.at(1);
        lljjHigh  = *leptons.at(0) + *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        l1jjHigh  = *leptons.at(0) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        l2jjHigh  = *leptons.at(1) + jets_WCandHighMass.at(0) + jets_WCandHighMass.at(1);
        dRl2WCand = leptons.at(0)->DeltaR(WCand);

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max); 
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCand_Mass_nocut_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCandJet1_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCandJet2_Pt_nocut_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/lljjLow_Mass_nocut_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1jjLow_Mass_nocut_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2jjLow_Mass_nocut_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/lljjHigh_Mass_nocut_"+IDsuffix, lljjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1jjHigh_Mass_nocut_"+IDsuffix, l1jjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2jjHigh_Mass_nocut_"+IDsuffix, l2jjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2WCand_DeltaR_nocut_"+IDsuffix, dRl2WCand, weight, 60, 0., 6.);

        // Low mass SR1 
        if(it_rg == 1){
          if(!(Nbjet_medium == 0)) continue;
          if(!(lljjLow.M() < 300.)) continue;
          if(!(MET < 80.)) continue;
        }

        // Low mass CR1
        if(it_rg == 2){
          if(!(lljjLow.M() < 300.)) continue;
          if(!(Nbjet_medium>0 || MET>100.)) continue;
        }

        // High mass SR1
        if(it_rg == 3){
          if(!(Nbjet_medium == 0)) continue;
          if(!(WCand.M() < 150.)) continue;
          if(!(MET2ST < 15.)) continue;
        }

        // High mass CR1
        if(it_rg == 4){
          if(!(WCand.M() < 150.)) continue;
          if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
        }

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCand_Mass_"+IDsuffix, WCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCandJet1_Pt_"+IDsuffix, jets_WCandHighMass.at(0).Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/WCandJet2_Pt_"+IDsuffix, jets_WCandHighMass.at(1).Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/lljjLow_Mass_"+IDsuffix, lljjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1jjLow_Mass_"+IDsuffix, l1jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2jjLow_Mass_"+IDsuffix, l2jjLow.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/lljjHigh_Mass_"+IDsuffix, lljjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1jjHigh_Mass_"+IDsuffix, l1jjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2jjHigh_Mass_"+IDsuffix, l2jjHigh.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2WCand_DeltaR_"+IDsuffix, dRl2WCand, weight, 60, 0., 6.);

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

      // Low mass SR2, CR2
      if(it_rg>=5 && it_rg<7){

        if(!(jets.size()==1 && fatjets.size()==0)) continue;

        llj = *leptons.at(0) + *leptons.at(1) + jets.at(0);
        l1j = *leptons.at(0) + jets.at(0);
        l2j = *leptons.at(1) + jets.at(0);

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_DeltaR_nocut_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/llj_Mass_nocut_"+IDsuffix, llj.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1j_Mass_nocut_"+IDsuffix, l1j.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2j_Mass_nocut_"+IDsuffix, l2j.M(), weight, 2000, 0., 2000.);

        // Low mass SR2
        if(it_rg == 5){
          if(!(Nbjet_medium == 0)) continue;
          if(!(llj.M() < 300.)) continue;
          if(!(MET < 80.)) continue;
        }

        // Low mass CR2
        if(it_rg == 6){
          if(!(llj.M() < 300.)) continue;
          if(!(Nbjet_medium>0 || MET>100.)) continue;
        }

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_DeltaR_"+IDsuffix, dRll, weight, 60, 0., 6.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/llj_Mass_"+IDsuffix, llj.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1j_Mass_"+IDsuffix, l1j.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2j_Mass_"+IDsuffix, l2j.M(), weight, 2000, 0., 2000.);

        //if(it_rg==5) cout << "In LowSR2, RunNumber:Lumi:EventNumber = " << run << ":" << lumi << ":" << event << endl;

      }

      // High mass SR2, CR2
      if(it_rg >= 7){

        if(!(fatjets.size() > 0)) continue;

        fatjets_WCand = FatJetWCand(fatjets, MW);

        l1J = *leptons.at(0) + fatjets_WCand;
        l2J = *leptons.at(1) + fatjets_WCand;

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_nocut_"+IDsuffix, 8.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_nocut_"+IDsuffix, 8.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_nocut_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_nocut_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_nocut_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_nocut_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_nocut_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_nocut_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_nocut_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_nocut_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_nocut_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_nocut_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_nocut_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_nocut_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Fatjet_Pt_nocut_"+IDsuffix, fatjets_WCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Fatjet_Mass_nocut_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1J_Mass_nocut_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2J_Mass_nocut_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

        // High mass SR2
        if(it_rg == 7){
          if(!(Nbjet_medium == 0)) continue;
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(MET2ST < 15.)) continue;
        }

        // High mass CR2
        if(it_rg == 8){
          if(!(fatjets_WCand.SDMass() < 150.)) continue;
          if(!(Nbjet_medium>0 || MET2ST>20.)) continue;
        }

        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_"+IDsuffix, 9.5, weight, cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Events_unweighted_"+IDsuffix, 9.5, 1., cutflow_bin, 0., cutflow_max);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_Jets_"+IDsuffix, jets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Loose_"+IDsuffix, Nbjet_loose, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_BJets_Medium_"+IDsuffix, Nbjet_medium, weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Number_FatJets_"+IDsuffix, fatjets.size(), weight, 10, 0., 10.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Mass_"+IDsuffix, ZCand.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/ZCand_Pt_"+IDsuffix, ZCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Pt_"+IDsuffix, leptons.at(0)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Pt_"+IDsuffix, leptons.at(1)->Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep1_Eta_"+IDsuffix, leptons.at(0)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Lep2_Eta_"+IDsuffix, leptons.at(1)->Eta(), weight, 50, -2.5, 2.5);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET_"+IDsuffix, MET, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/MET2ST_"+IDsuffix, MET2ST, weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Fatjet_Pt_"+IDsuffix, fatjets_WCand.Pt(), weight, 1000, 0., 1000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/Fatjet_Mass_"+IDsuffix, fatjets_WCand.SDMass(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l1J_Mass_"+IDsuffix, l1J.M(), weight, 2000, 0., 2000.);
        FillHist(systName+"/"+channel+"/"+regions.at(it_rg)+"/l2J_Mass_"+IDsuffix, l2J.M(), weight, 2000, 0., 2000.);

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

