#include "HNtypeI_FakeRate.h"

HNtypeI_FakeRate::HNtypeI_FakeRate(){

}

void HNtypeI_FakeRate::initializeAnalyzer(){

  //==== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");
  RunNorm = HasFlag("RunNorm");
  RunSF   = HasFlag("RunSF");

  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSyst = " << RunSyst << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunNorm = " << RunNorm << endl;
  cout << "[HNtypeI_FakeRate::initializeAnalyzer] RunSF = " << RunSF << endl;

  MuonTightIDs     = {"HNTightV2"};
  MuonLooseIDs     = {"HNLooseV2"};
  MuonVetoIDs      = {"HNVeto"};
  ElectronTightIDs = {"HNTightV2"};
  ElectronLooseIDs = {"HNLooseV1"};
  ElectronVetoIDs  = {"HNVeto"};

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/HNtypeI_FakeRate.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  //==== Muon triggers
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

  //==== Electron triggers
  //==== DoubleEG (2016), SingleElectron (2017), EGamma (2018)
  ElectronTrig1 = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  ElectronTrig2 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  if(DataEra.Contains("2016")) ElectronTrig3 = "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  else ElectronTrig3 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  ElectronTrig4 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

  ElectronTriggers.push_back(ElectronTrig1);
  ElectronTriggers.push_back(ElectronTrig2);
  ElectronTriggers.push_back(ElectronTrig3);
  ElectronTriggers.push_back(ElectronTrig4);
  ElectronPtCut1 = 10., ElectronPtCut2 = 15., ElectronPtCut3 = 20., ElectronPtCut4 = 25.;
  ElectronPtconeCut1 = 15., ElectronPtconeCut2 = 25., ElectronPtconeCut3 = 35., ElectronPtconeCut4 = 45.;

  //==== B tagging
  //==== Add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb) );
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

}

HNtypeI_FakeRate::~HNtypeI_FakeRate(){

  //==== Destructor of this Analyzer

}

void HNtypeI_FakeRate::executeEvent(){

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

  //==== Declare AnalyzerParameter

  AnalyzerParameter param;

  for(unsigned int it_id=0; it_id<ElectronTightIDs.size(); it_id++){

    TString MuonTightID = MuonTightIDs.at(it_id);
    TString MuonLooseID = MuonLooseIDs.at(it_id);
    TString MuonVetoID  = MuonVetoIDs.at(it_id);
    TString ElectronTightID = ElectronTightIDs.at(it_id);
    TString ElectronLooseID = ElectronLooseIDs.at(it_id);
    TString ElectronVetoID  = ElectronVetoIDs.at(it_id);

    param.Clear();

    param.fakesyst_ = AnalyzerParameter::FakeCentral;

    param.Name = "FakeCentral";

    //==== Muon ID
    param.Muon_Tight_ID = MuonTightID;
    param.Muon_Loose_ID = MuonLooseID;
    param.Muon_Veto_ID  = MuonVetoID;
    param.Muon_ID_SF_Key = "";
    param.Muon_ISO_SF_Key = "";

    //==== Electron ID
    param.Electron_Tight_ID = ElectronTightID;
    param.Electron_Loose_ID = ElectronLooseID;
    param.Electron_Veto_ID  = ElectronVetoID;
    param.Electron_ID_SF_Key = "";

    //==== Jet ID
    param.Jet_ID = "HNTight";

    executeEventFromParameter(param);

    if(RunSyst){

      /*for(int it_syst=1; it_syst<AnalyzerParameter::NFakeSyst; it_syst++){
        param.fakesyst_ = AnalyzerParameter::FakeSyst(it_syst);
        param.Name  = "FakeSyst_"+param.GetFakeSystType();
        executeEventFromParameter(param);
      }*/

    }

  }

}

void HNtypeI_FakeRate::executeEventFromParameter(AnalyzerParameter param){

  TString MuonIDName = "HNTightV2";

  TString ElectronIDName = "HNTightV2";

  vector<TString> regions = {"FR", "DY", "WJ"};

  //==== Boolean : primary datasets
  bool isMuon = false, isElectron = false;
  if(IsDATA){
    if(DataStream.Contains("SingleMuon") || DataStream.Contains("DoubleMuon")) isMuon = true;
    if(DataStream.Contains("DoubleEG") || DataStream.Contains("SingleElectron") || DataStream.Contains("EGamma")) isElectron = true;
  }

  TString systName = param.Name;

  //========================================================
  //==== Luminosity of prescaled triggers
  //========================================================

  //==== Lumi values from brilcalc
  if(DataEra == "2016preVFP"){
    MuonLumi1 = 3.911703648, MuonLumi2 = 6.585395243, MuonLumi3 = 192.124599592;
    ElectronLumi1 = 4.127188014, ElectronLumi2 = 11.033890219, ElectronLumi3 = 52.443756312, ElectronLumi4 = 52.790026618;
  }

  if(DataEra == "2016postVFP"){
    MuonLumi1 = 3.580526616, MuonLumi2 = 1.296837118, MuonLumi3 = 26.872793405;
    ElectronLumi1 = 2.938006053, ElectronLumi2 = 3.980060594, ElectronLumi3 = 6.834870773, ElectronLumi4 = 10.699483420;
  }

  if(DataEra == "2017"){
    MuonLumi1 = 4.607782551, MuonLumi2 = 2.899706509, MuonLumi3 = 65.898920913;
    ElectronLumi1 = 3.970041451, ElectronLumi2 = 27.683553584, ElectronLumi3 = 27.683553584, ElectronLumi4 = 43.453223348;
  }

  if(DataEra == "2018"){
    MuonLumi1 = 2.704239929, MuonLumi2 = 8.581579807, MuonLumi3 = 45.852032058;
    ElectronLumi1 = 6.424543223, ElectronLumi2 = 38.917235487, ElectronLumi3 = 38.917235487, ElectronLumi4 = 38.973910597;
  }

  //==== Lumi SF (Muon)
  if(param.Muon_Tight_ID.Contains("HNTightV2")){

    if(DataEra == "2016preVFP"){
      SFMuonLumi1 = 1.02112, SFMuonLumi2 = 1.52457, SFMuonLumi3 = 1.11109;
    }

    if(DataEra == "2016postVFP"){
      SFMuonLumi1 = 0.635706, SFMuonLumi2 = 1.26902, SFMuonLumi3 = 1.07973;
    }

    if(DataEra == "2017"){
      SFMuonLumi1 = 1.27997, SFMuonLumi2 = 1.50086, SFMuonLumi3 = 1.12999;
    }

    if(DataEra == "2018"){
      SFMuonLumi1 = 2.10869, SFMuonLumi2 = 1.15663, SFMuonLumi3 = 1.00586;
    }

  }

  //==== Lumi SF (Electron)
  if(param.Electron_Tight_ID.Contains("HNTightV2")){

    if(DataEra == "2016preVFP"){
      SFElectronLumi1 = 1.15521, SFElectronLumi2 = 1.06986, SFElectronLumi3 = 1.06932, SFElectronLumi4 = 1.06766;
    }

    if(DataEra == "2016postVFP"){
      SFElectronLumi1 = 1.17585, SFElectronLumi2 = 1.15739, SFElectronLumi3 = 0.957963, SFElectronLumi4 = 0.968462;
    }

    if(DataEra == "2017"){
      SFElectronLumi1 = 1.28124, SFElectronLumi2 = 1.12829, SFElectronLumi3 = 1.12829, SFElectronLumi4 = 1.03026;
    }

    if(DataEra == "2018"){
      SFElectronLumi1 = 1.11554, SFElectronLumi2 = 1.13261, SFElectronLumi3 = 1.13261, SFElectronLumi4 = 0.971968;
    }

  }

  //==== No lumi SF
  if(RunNorm && !RunSF){

    if(DataEra == "2016preVFP"){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1.;
    }

    if(DataEra == "2016postVFP"){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1.;
    }

    if(DataEra == "2017"){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1.;
    }

    if(DataEra == "2018"){
      SFMuonLumi1 = 1., SFMuonLumi2 = 1., SFMuonLumi3 = 1.;
      SFElectronLumi1 = 1., SFElectronLumi2 = 1., SFElectronLumi3 = 1., SFElectronLumi4 = 1.;
    }

  }

  Event ev = GetEvent();

  //========================================================
  //==== No Cut
  //========================================================

  //========================================================
  //==== MET Filter
  //========================================================

  if(!PassMETFilter()) return;

  //========================================================
  //==== Trigger
  //========================================================

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
  double jetEtaCut = 4.7;

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
  else if(param.fakesyst_ == AnalyzerParameter::dPhi1){
    dPhiCut = 1.7;
  }
  else if(param.fakesyst_ == AnalyzerParameter::dPhi2){
    dPhiCut = 2.1;
  }
  else if(param.fakesyst_ == AnalyzerParameter::dPhi3){
    dPhiCut = 2.9;
  }
  else if(param.fakesyst_ == AnalyzerParameter::PtRatioUp){
    PtRatioCut = 1.2;
  }
  else if(param.fakesyst_ == AnalyzerParameter::PtRatioDown){
    PtRatioCut = 0.8;
  }
  else{
    cerr << "[HNtypeI_FakeRate::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //========================================================
  //==== Then, apply ID selections using this_AllXXX
  //========================================================

  //==== Leptons
  vector<Muon> muons_tight, muons_loose, muons_veto;
  muons_tight.clear();
  muons_loose.clear();
  muons_veto.clear();

  vector<Electron> electrons_tight, electrons_loose, electrons_veto;
  electrons_tight.clear();
  electrons_loose.clear();
  electrons_veto.clear();

  muons_tight = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 10., 2.4);
  muons_loose = SelectMuons(this_AllMuons, param.Muon_Loose_ID, MuonPtCut1, 2.4);
  muons_veto  = SelectMuons(this_AllMuons, param.Muon_Veto_ID, MuonPtCut1, 2.4);

  electrons_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 10., 2.5);
  electrons_loose = SelectElectrons(this_AllElectrons, param.Electron_Loose_ID, ElectronPtCut1, 2.5);
  electrons_veto  = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, ElectronPtCut1, 2.5);

  //==== Truth matching
  vector<Muon> muons_prompt;
  vector<Electron> electrons_prompt;
  muons_prompt.clear();
  electrons_prompt.clear();

  //==== Jets
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 20., jetEtaCut);
  vector<Jet> jets_bcand = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);

  /*if(muons_veto.size()+electrons_veto.size() == 1){
    jets = JetsVetoLeptonInside(jets_nolepveto, electrons_veto, muons_veto); // We use jets only for FR measurement
  }*/

  vector<Jet> jets_awayFromMuon;
  vector<Jet> jets_awayFromElectron;
  jets_awayFromMuon.clear();
  jets_awayFromElectron.clear();

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
  std::sort(jets_bcand.begin(), jets_bcand.end(), PtComparing);

  //========================================================
  //==== B tagging
  //========================================================

  int Nbjet_loose = 0, Nbjet_medium = 0;
  JetTagging::Parameters jtp_DeepJet_Loose = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepJet_Medium = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

  //==== method 1a)
  //==== multiply "btagWeight" to the event weight
  //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

  //==== method 2a)
  for(unsigned int ij=0; ij<jets_bcand.size(); ij++){

    if(mcCorr->IsBTagged_2a(jtp_DeepJet_Loose, jets_bcand.at(ij), "central")) Nbjet_loose++;
    if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_bcand.at(ij), "central")) Nbjet_medium++;

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

  //==== POG cut-based Medium  
  //==== barrel : 0.0478+0.506/pT, endcap : 0.0658+0.963/pT
  //==== POG cut-based Tight
  //==== barrel : 0.0287+0.506/pT, endcap : 0.0445+0.963/pT

  double MZ = 91.1876;
  double weight = 1.;
  double Mt = 0.;
  double pT_ratio = 0.;
  double jet_EmFraction = 0.;

  double trigger_lumi = 1.;
  double jetPtCut = jetPtCut_syst;

  bool IsAwayJetBTag = false, IsCloseJetBTag = false;

  double pTcone_mu = 0., pTcone_el = 0.;
  TString PtConeRange = "";

  Particle ZCand, METv;

  if(systName == "FakeCentral"){
    FillHist("MET_NoCut", MET, weight, 500, 0., 500.);
    FillHist("METPhi_NoCut", METPhi, weight, 32, 0., 3.2);
  }

  //========================================================
  //==== Muon
  //========================================================

  for(unsigned int it_rg=0; it_rg<regions.size(); it_rg++){

    weight = 1., muonIDSF = 1., muonIsoSF = 1.;

    if(!(muons_loose.size() > 0)) break;
    if(!ev.PassTrigger(MuonTriggers)) break;
    if(IsDATA){ if(!isMuon) break; }

    //==== Fake rate measurement region
    if(it_rg == 0){

      if(RunNorm) continue;

      if(!(muons_loose.size()==1 && electrons_loose.size()==0)) continue;
      if(!(muons_veto.size()==1 && electrons_veto.size()==0)) continue;
      if(!(jets.size() >= 1)) continue;

      //==== MET
      METv = UpdateMETMuon(METv_central, muons_loose);
      METv = UpdateMETSmearedJet(METv, jets);
      MET = METv.Pt();
      METPhi = METv.Phi();

      //==== pTcone
      //if(RunPt) pTcone_mu = muons_loose.at(0).Pt();
      pTcone_mu = muons_loose.at(0).CalcPtCone(muons_loose.at(0).RelIso(), mu_tight_iso);

      //==== Truth matching
      muons_prompt.clear();
      muons_prompt = MuonPromptOnlyHNtypeI(muons_loose, gens);
      if(!(muons_prompt.size() == 1)) continue;

      //==== Event weights except trigger luminosity
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

      //==== Passing one prescaled trigger for each pTcone range, applying trigger lumi

      trigger_lumi = 1.;

      if(!(pTcone_mu >= MuonPtconeCut1)) continue;

      if(pTcone_mu >= MuonPtconeCut1 && pTcone_mu < MuonPtconeCut2){
        if(!(muons_loose.at(0).Pt() > MuonPtCut1)) continue;
        if(!ev.PassTrigger(MuonTrig1)) continue;
        if(!IsDATA) trigger_lumi = MuonLumi1*SFMuonLumi1;
        if(jetPtCut_syst == 40.) jetPtCut = 50.;
        PtConeRange = "Range0";
      }

      if(pTcone_mu >= MuonPtconeCut2 && pTcone_mu < MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut2)) continue;
        if(!ev.PassTrigger(MuonTrig2)) continue;
        if(!IsDATA) trigger_lumi = MuonLumi2*SFMuonLumi2;
        PtConeRange = "Range1";
      }

      if(pTcone_mu >= MuonPtconeCut3){
        if(!(muons_loose.at(0).Pt() > MuonPtCut3)) continue;
        if(!ev.PassTrigger(MuonTrig3)) continue;
        if(!IsDATA) trigger_lumi = MuonLumi3*SFMuonLumi3;
        PtConeRange = "Range2";
      }

      weight *= trigger_lumi;

      //==== Away jet selection
      jets_awayFromMuon.clear();
      jets_awayFromMuon = JetsAwayFromLepton(jets, muons_loose.at(0), dPhiCut);
      std::sort(jets_awayFromMuon.begin(), jets_awayFromMuon.end(), PtComparing);

      if(!(jets_awayFromMuon.size() > 0)) continue;
      if(!(jets_awayFromMuon.at(0).Pt() > jetPtCut)) continue;

      //==== B tagging for the away/close jet 
      IsAwayJetBTag = false, IsCloseJetBTag = false;

      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_awayFromMuon.at(0))) IsAwayJetBTag = true;

      for(unsigned int ij=0; ij<jets_bcand.size(); ij++){
        if(jets_bcand.at(ij).DeltaR(muons_loose.at(0)) < 0.4){
          if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_bcand.at(ij))) IsCloseJetBTag = true;
        }
        if(IsCloseJetBTag) break;
      }

      Mt = MT(muons_loose.at(0), METv);
      pT_ratio = jets_awayFromMuon.at(0).Pt()/muons_loose.at(0).Pt();

      //==== Histograms before applying cuts
      if(systName=="FakeCentral"){

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_MET_NoCut", MET, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_METPhi_NoCut", METPhi, weight, 32, 0., 3.2);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_METPhi_NoCut_"+PtConeRange, METPhi, weight, 32, 0., 3.2);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mt_NoCut", Mt, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Ptratio_NoCut", pT_ratio, weight, 50, 0., 5.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Ptratio_NoCut_"+PtConeRange, pT_ratio, weight, 50, 0., 5.);

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_PtCone_NoCut", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_PtCone_NoCut_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Pt_NoCut", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Eta_NoCut", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(muons_tight.size() == 0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_PtCone_NoCut", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_PtCone_NoCut_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Pt_NoCut", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Eta_NoCut", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Eta_NoCut_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 1.5, weight, 3, 0., 3.);
        }

        if(muons_tight.size() > 0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_PtCone_NoCut", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_PtCone_NoCut_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Pt_NoCut", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Eta_NoCut", muons_tight.at(0).Eta(), weight, 60, -3., 3.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Eta_NoCut_"+PtConeRange, muons_tight.at(0).Eta(), weight, 60, -3., 3.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_NoCut_"+PtConeRange, 2.5, weight, 3, 0., 3.);
        }

      }

      //==== Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(pT_ratio > PtRatioCut)) continue;

      //==== Histograms after applying cuts
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Eta", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Loose_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
      FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

      if(muons_tight.size() == 0){
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Eta", muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_NoTight_Eta_"+PtConeRange, muons_loose.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
      }

      if(muons_tight.size() > 0){
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Eta", muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Muon_Tight_Eta_"+PtConeRange, muons_tight.at(0).Eta(), weight, 60, -3., 3.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      //==== Inner barrel ( |eta| < 0.8 )
      if(fabs(muons_loose.at(0).Eta()) < 0.8){

        //==== Passing loose ID
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        //==== Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){

          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        //==== Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_IB_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      //==== Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(muons_loose.at(0).Eta()) >= 0.8 && fabs(muons_loose.at(0).Eta()) < 1.479){

        //==== Passing loose ID
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
              
        }   

        //==== Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){
        
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
          
          if(systName == "FakeCentral"){
          
            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            } 
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            } 
              
          }   
              
        }   

        //==== Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_OB_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      //==== Endcap ( 1.479 < |eta| < 2.4 )
      if(fabs(muons_loose.at(0).Eta()) >= 1.479 && fabs(muons_loose.at(0).Eta()) < 2.4){

        //==== Passing loose ID
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Loose_PtCone", pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Loose_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Loose_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        //==== Passing loose ID but not passing tight ID
        if(muons_tight.size() == 0){

          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_NoTight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_NoTight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_NoTight_Pt", muons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        //==== Passing tight ID
        if(muons_tight.size() > 0){

          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Tight_PtCone", pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Tight_PtCone_"+PtConeRange, pTcone_mu, weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Muon_Tight_Pt", muons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_EC_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

    }

    //==== DY CR

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

      //==== Histograms
      if(ev.PassTrigger(MuonTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep1_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep2_Pt_NoCut", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep1_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep2_Pt_NoCut_OnlyTight", muons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      //==== Event selection criteria
      if(!(muons_tight.at(0).Pt()>20. && muons_tight.at(1).Pt()>15.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      //==== Histograms
      if(ev.PassTrigger(MuonTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==2 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==2 && electrons_veto.size()==0){
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

    }

    //==== W+jets CR
    if(it_rg == 2){

      if(RunSyst) continue;

      if(!(muons_tight.size()==1 && electrons_tight.size()==0)) continue;
      //if(!(muons_veto.size()==1 && electrons_veto.size()==0)) continue;

      //==== MET
      METv = UpdateMETMuon(METv_central, muons_tight);
      METv = UpdateMETSmearedJet(METv, jets);
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

      //==== Histograms
      if(ev.PassTrigger(MuonTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep_Pt_NoCut", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Lep_Pt_NoCut_OnlyTight", muons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      //==== Event selection
      if(!(muons_tight.at(0).Pt() > 20.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;

      //==== Histograms
      if(ev.PassTrigger(MuonTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi1*SFMuonLumi1;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu3_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(MuonTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi2*SFMuonLumi2;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu8_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(MuonTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = MuonLumi3*SFMuonLumi3;

        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==1 && electrons_veto.size()==0){
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(MuonIDName+"_"+systName+"_"+regions.at(it_rg)+"_Mu17_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

    }

  }

  //========================================================
  //==== Electron
  //========================================================

  for(unsigned int it_rg2=0; it_rg2<regions.size(); it_rg2++){

    weight = 1.;

    if(!(electrons_loose.size() > 0)) break;
    if(!ev.PassTrigger(ElectronTriggers)) break;
    if(IsDATA){ if(!isElectron) break; }

    //==== Fake rate measurement region
    if(it_rg2 == 0){

      if(RunNorm) continue;

      if(!(muons_loose.size()==0 && electrons_loose.size()==1)) continue;
      if(!(muons_veto.size()==0 && electrons_veto.size()==1)) continue;
      if(!(jets.size() >= 1)) continue;

      //==== MET
      METv = UpdateMETElectron(METv_central, electrons_loose);
      METv = UpdateMETSmearedJet(METv, jets);
      MET = METv.Pt();
      METPhi = METv.Phi();

      //==== pTcone
      el_tight_iso = 0.08; // 2016 or POG MVA

      if(param.Electron_Tight_ID.Contains("ISR")){ // POG cut-based Medium
        el_tight_iso = 0.0478+0.506/electrons_loose.at(0).UncorrPt();
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0658+0.963/electrons_loose.at(0).UncorrPt();
      }

      if(param.Electron_Tight_ID.Contains("HNTight")){ // POG cut-based Tight
        el_tight_iso = 0.0287+0.506/electrons_loose.at(0).UncorrPt();
        if(fabs(electrons_loose.at(0).scEta()) > 1.479) el_tight_iso = 0.0445+0.963/electrons_loose.at(0).UncorrPt();
      }

      //if(RunPt) pTcone_el = electrons_loose.at(0).Pt();
      pTcone_el = electrons_loose.at(0).CalcPtCone(electrons_loose.at(0).RelIso(), el_tight_iso);

      //==== Truth matching
      electrons_prompt.clear();
      electrons_prompt = ElectronPromptOnlyHNtypeI(electrons_loose, gens);
      if(!(electrons_prompt.size() == 1)) continue;

      //==== Event weights except trigger luminosity
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

      //==== Passing one prescaled trigger for each PtCone range, applying trigger lumi

      trigger_lumi = 1.;

      if(!(pTcone_el >= ElectronPtconeCut1)) continue;

      if(pTcone_el >= ElectronPtconeCut1 && pTcone_el < ElectronPtconeCut2){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut1)) continue;
        if(!ev.PassTrigger(ElectronTrig1)) continue;
        if(!IsDATA) trigger_lumi = ElectronLumi1*SFElectronLumi1;
        PtConeRange = "Range0";
      }

      if(pTcone_el >= ElectronPtconeCut2 && pTcone_el < ElectronPtconeCut3){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut2)) continue;
        if(!ev.PassTrigger(ElectronTrig2)) continue;
        if(!IsDATA) trigger_lumi = ElectronLumi2*SFElectronLumi2;
        PtConeRange = "Range1";
      }

      if(pTcone_el >= ElectronPtconeCut3 && pTcone_el < ElectronPtconeCut4){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut3)) continue;
        if(!ev.PassTrigger(ElectronTrig3)) continue;
        if(!IsDATA) trigger_lumi = ElectronLumi3*SFElectronLumi3;
        PtConeRange = "Range2";
      }

      if(pTcone_el >= ElectronPtconeCut4){
        if(!(electrons_loose.at(0).Pt() > ElectronPtCut4)) continue;
        if(!ev.PassTrigger(ElectronTrig4)) continue;
        if(!IsDATA) trigger_lumi = ElectronLumi4*SFElectronLumi4;
        PtConeRange = "Range3";
      }

      weight *= trigger_lumi;

      //==== Away jet selection
      jets_awayFromElectron.clear();
      jets_awayFromElectron = JetsAwayFromLepton(jets, electrons_loose.at(0), dPhiCut);
      std::sort(jets_awayFromElectron.begin(), jets_awayFromElectron.end(), PtComparing);

      if(!(jets_awayFromElectron.size() > 0)) continue;
      if(!(jets_awayFromElectron.at(0).Pt() > jetPtCut)) continue;

      //==== B tagging for the away/close jet
      IsAwayJetBTag = false, IsCloseJetBTag = false;

      if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_awayFromElectron.at(0))) IsAwayJetBTag = true;

      for(unsigned int ij=0; ij<jets_bcand.size(); ij++){
        if(jets_bcand.at(ij).DeltaR(electrons_loose.at(0)) < 0.4){
          if(mcCorr->IsBTagged_2a(jtp_DeepJet_Medium, jets_bcand.at(ij))) IsCloseJetBTag = true;
        }
        if(IsCloseJetBTag) break;
      }

      Mt = MT(electrons_loose.at(0), METv);
      pT_ratio = jets_awayFromElectron.at(0).Pt()/electrons_loose.at(0).Pt();
      jet_EmFraction = jets_awayFromElectron.at(0).ChargedEmEnergyFraction();

      //==== Histograms before applying cuts
      if(systName=="FakeCentral"){

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_MET_NoCut", MET, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_METPhi_NoCut", METPhi, weight, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_MET_NoCut_"+PtConeRange, MET, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_METPhi_NoCut_"+PtConeRange, METPhi, weight, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Mt_NoCut", Mt, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Mt_NoCut_"+PtConeRange, Mt, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ptratio_NoCut", pT_ratio, weight, 50, 0., 5.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ptratio_NoCut_"+PtConeRange, pT_ratio, weight, 50, 0., 5.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction", jet_EmFraction, weight, 100, 0., 1.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_"+PtConeRange, jet_EmFraction, weight, 100, 0., 1.);

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_PtCone_NoCut", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_PtCone_NoCut_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Pt_NoCut", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Eta_NoCut", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Eta_NoCut_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_NoCut_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(electrons_tight.size() == 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_PtCone_NoCut", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_PtCone_NoCut_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Pt_NoCut", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Eta_NoCut", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Eta_NoCut_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_NoCut_"+PtConeRange, 1.5, weight, 3, 0., 3.);
        }

        if(electrons_tight.size() > 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_PtCone_NoCut", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_PtCone_NoCut_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Pt_NoCut", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Eta_NoCut", electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Eta_NoCut_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_NoCut_"+PtConeRange, 2.5, weight, 3, 0., 3.);
        }

      }

      //==== Additional cuts to reduce prompt contribution
      if(!(MET < 80.)) continue;
      if(!(Mt < 25.)) continue;
      if(!(pT_ratio > PtRatioCut)) continue;

      if(systName=="FakeCentral"){
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_WithCuts", jet_EmFraction, weight, 100, 0., 1.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Jet_ChargedEmEnergyFraction_WithCuts_"+PtConeRange, jet_EmFraction, weight, 100, 0., 1.);
      }

      if(!(jet_EmFraction < 0.65)) continue;

      //==== Histograms after applying cuts
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Eta", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Loose_Eta_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
      FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

      if(electrons_tight.size() == 0){
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Eta", electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_NoTight_Eta_"+PtConeRange, electrons_loose.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
      }

      if(electrons_tight.size() > 0){
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Eta", electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Electron_Tight_Eta_"+PtConeRange, electrons_tight.at(0).scEta(), weight, 60, -3., 3.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
      }

      //==== Inner barrel ( |eta| < 0.8 )
      if(fabs(electrons_loose.at(0).scEta()) < 0.8){

        //==== Passing loose ID
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        //==== Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }

        //==== Passing tight ID
        if(electrons_tight.size() > 0){

          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_IB_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      //==== Outer barrel ( 0.8 < |eta| < 1.479 )
      if(fabs(electrons_loose.at(0).scEta()) >= 0.8 && fabs(electrons_loose.at(0).scEta()) < 1.479){

        //==== Passing loose ID
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        //==== Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }
            
        //==== Passing tight ID
        if(electrons_tight.size() > 0){

          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.); 
              
          if(systName == "FakeCentral"){
              
            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_OB_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

      //==== Endcap ( 1.479 < |eta| < 2.5 )
      if(fabs(electrons_loose.at(0).scEta()) >= 1.479 && fabs(electrons_loose.at(0).scEta()) < 2.5){

        //==== Passing loose ID
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);

        if(systName == "FakeCentral"){

          if(IsAwayJetBTag){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }
          else{
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Loose_PtCone", pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Loose_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Loose_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Number_Events_"+PtConeRange, 0.5, weight, 3, 0., 3.);
          }

        }

        //==== Passing loose ID but not passing tight ID
        if(electrons_tight.size() == 0){

          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);

          if(systName == "FakeCentral"){

            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_NoTight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_NoTight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_NoTight_Pt", electrons_loose.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Number_Events_"+PtConeRange, 1.5, weight, 3, 0., 3.);
            }

          }

        }
            
        //==== Passing tight ID
        if(electrons_tight.size() > 0){
              
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.); 
              
          if(systName == "FakeCentral"){
              
            if(IsAwayJetBTag){
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_BTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }
            else{
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Tight_PtCone", pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Tight_PtCone_"+PtConeRange, pTcone_el, weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Electron_Tight_Pt", electrons_tight.at(0).Pt(), weight, 500, 0., 500.);
              FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_EC_NoBTag_Number_Events_"+PtConeRange, 2.5, weight, 3, 0., 3.);
            }

          }

        }

      }

    }

    //==== DY CR
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

      //==== Histograms
      if(ev.PassTrigger(ElectronTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep1_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep2_Pt_NoCut", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_NoCut", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep1_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep2_Pt_NoCut_OnlyTight", electrons_tight.at(1).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_NoCut_OnlyTight", ZCand.M(), weight*trigger_lumi, 80, 50., 130.);
        }

      }

      //==== Event selection criteria
      if(!(electrons_tight.at(0).Pt()>25. && electrons_tight.at(1).Pt()>15.)) continue;
      if(!(fabs(ZCand.M() - MZ) < 10.)) continue;

      //==== Histograms
      if(ev.PassTrigger(ElectronTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_Inclusive", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==2){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_Inclusive_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

        if(electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events", 1.5, weight*trigger_lumi, 2, 0., 2.);

          if(muons_veto.size()==0 && electrons_veto.size()==2){
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_ZCand_Mass_OnlyTight", ZCand.M(), weight*trigger_lumi, 50, 70., 120.);
            FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 1.5, weight*trigger_lumi, 2, 0., 2.);
          }
        }

      }

    }


    //==== W+jets CR
    if(it_rg2 == 2){

      if(RunSyst) continue;

      if(!(muons_tight.size()==0 && electrons_tight.size()==1)) continue;
      //if(!(muons_veto.size()==0 && electrons_veto.size()==1)) continue;

      //==== MET
      METv = UpdateMETElectron(METv_central, electrons_tight);
      METv = UpdateMETSmearedJet(METv, jets);
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

      //==== Histograms
      if(ev.PassTrigger(ElectronTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep_Pt_NoCut", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_MET_NoCut", MET, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_METPhi_NoCut", METPhi, weight*trigger_lumi, 32, 0., 3.2);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Mt_NoCut", Mt, weight*trigger_lumi, 500, 0., 500.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Lep_Pt_NoCut_OnlyTight", electrons_tight.at(0).Pt(), weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_MET_NoCut_OnlyTight", MET, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_METPhi_NoCut_OnlyTight", METPhi, weight*trigger_lumi, 32, 0., 3.2);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Mt_NoCut_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
        }

      }

      //==== Event selection criteria
      if(!(electrons_tight.at(0).Pt() > 25.)) continue;
      if(!(MET > 40.)) continue;
      if(!(Mt>60. && Mt<100.)) continue;

      //==== Histograms for each trigger
      if(ev.PassTrigger(ElectronTrig1)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi1*SFElectronLumi1;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele8_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig2)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi2*SFElectronLumi2;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele12_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig3)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi3*SFElectronLumi3;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele17_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

      if(ev.PassTrigger(ElectronTrig4)){

        trigger_lumi = 1.;
        if(!IsDATA) trigger_lumi = ElectronLumi4*SFElectronLumi4;

        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Mt", Mt, weight*trigger_lumi, 500, 0., 500.);
        FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events", 0.5, weight*trigger_lumi, 2, 0., 2.);

        if(muons_veto.size()==0 && electrons_veto.size()==1){
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Mt_OnlyTight", Mt, weight*trigger_lumi, 500, 0., 500.);
          FillHist(ElectronIDName+"_"+systName+"_"+regions.at(it_rg2)+"_Ele23_Number_Events_OnlyTight", 0.5, weight*trigger_lumi, 2, 0., 2.);
        }

      }

    }

  }

}
