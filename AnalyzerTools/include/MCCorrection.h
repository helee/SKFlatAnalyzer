#ifndef MCCorrection_h
#define MCCorrection_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TString.h"
#include "TRegexp.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "Lepton.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "Jet.h"
#include "Gen.h"

#include "TDirectoryHelper.h"
#include "GenFinderForDY.h"
#include "BTagCalibrationStandalone.h"
#include "TRandom3.h"

class MCCorrection{

public:

  MCCorrection();
  ~MCCorrection();

  TDirectory *histDir;
  void ReadHistograms();

  TString MCSample;
  void SetMCSample(TString s);

  int DataYear;
  void SetDataYear(int i);

  bool IsDATA;
  void SetIsDATA(bool b);

  int run,lumi,event;
  void SetEventInfo(int r, int l, int e);

  bool IsFastSim;
  void SetIsFastSim(bool b);

  bool IgnoreNoHist;

  //=========================================================
  //==== Muon
  //=========================================================

  double MuonReco_SF(TString key, double eta, double p, int sys=0);
  double MuonID_SF(TString ID, double eta, double pt, int sys=0);
  double MuonISO_SF(TString ID, double eta, double pt, int sys=0);
  double MuonTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon>& muons, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon *>& muons, int sys=0);

  //==== For HNtypeI
  double MuonID_Eff_Period(TString ID, TString period, TString DataMC, double eta, double pt, int sys=0);
  double MuonID_SF_HNtypeI(TString ID, double eta, double pt, int sys=0);
  double MuonTrigger_Eff_Period(TString leg, TString ID, TString period, TString DataMC, double eta, double pt, double ptcone, bool RunFake, int sys=0);
  double MuonTrigger_Eff_HNtypeI(TString leg, TString ID, TString DataMC, double eta, double pt, double ptcone, bool RunFake, int sys=0);
  double MuonTrigger_SF_HNtypeI(TString ID, const std::vector<Muon>& muons, bool RunFake, int sys=0);

  std::map< TString, TH2F* > map_hist_Muon;

  //=========================================================
  //==== Electron
  //=========================================================

  double ElectronReco_SF(double sceta, double pt, int sys=0);
  double ElectronID_SF(TString ID, double sceta, double pt, int sys=0);
  double ElectronTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron>& electrons, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron *>& electrons, int sys=0);

  //==== For HNtypeI
  double ElectronTrigger_Eff_HNtypeI(TString leg, TString ID, TString DataMC, double sceta, double pt, double ptcone, bool RunFake, int sys=0);
  double ElectronTrigger_SF_HNtypeI(TString ID, const std::vector<Electron>& electrons, bool RunFake, int sys=0);
  double EMuTrigger_SF_HNtypeI(TString MuonID, TString ElectronID, const std::vector<Muon>& muons, const std::vector<Electron>& electrons, bool RunFake, int sys=0);

  std::map< TString, TH2F* > map_hist_Electron;
  std::map< TString, TGraphAsymmErrors* > map_graph_Electron;

  //=========================================================
  //==== FatJet
  //=========================================================

  double FatJetWtagSF(TString ID, int sys=0);  // For HNtypeI

  //=========================================================
  //==== Event reweight
  //=========================================================

  std::map< TString, TH2F* > map_hist_prefire;
  double GetPrefireWeight(const std::vector<Photon>& photons, const std::vector<Jet>& jets, int sys);

  std::map< TString, TH1D* > map_hist_pileup;
  double GetPileUpWeightBySampleName(int N_pileup, int syst);
  double GetPileUpWeight(int N_pileup, int syst);
  double GetPileUpWeight2017(int N_pileup, int syst); // For HNtypeI

  std::map< TString, TH1D* > map_hist_vertex;  // For HNtypeI
  double GetVertexWeight(int Nvtx, TString channel);

  std::map< TString, TH1D* > map_hist_rho;  // For HNtypeI
  double GetRhoWeight(double rho, TString channel);

  double GetTopPtReweight(const std::vector<Gen>& gens);

  TH2D *hist_DYPtReweight_2D;
  GenFinderForDY *genFinderDY;
  double GetOfficialDYReweight(const vector<Gen>& gens, int sys=0);

  //=========================================================
  //==== b tagging
  //=========================================================

  double GetJetTaggingCutValue(JetTagging::Tagger tagger, JetTagging::WP wp);

  std::vector<JetTagging::Parameters> jetTaggingPars;
  void SetJetTaggingParameters(std::vector<JetTagging::Parameters> v);

  std::map< std::string, BTagCalibrationReader* > map_BTagCalibrationReader;
  void SetupJetTagging();

  TH2D *hist_JetTagEff_B;
  TH2D *hist_JetTagEff_C;
  TH2D *hist_JetTagEff_Light;

  double GetMCJetTagEff(JetTagging::Tagger tagger, JetTagging::WP wp, int JetFlavor, double JetPt, double JetEta);
  double GetJetTaggingSF(JetTagging::Parameters jtp, int JetFlavor, double JetPt, double JetEta, double Jetdiscr, string Syst="central");

  //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

  //==== 1a) Event reweighting using scale factors and MC b-tagging efficiencies
  double GetBTaggingReweight_1a(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst="central");
  //==== 1d) Event reweighting using discriminator-dependent scale factors
  double GetBTaggingReweight_1d(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst="central");

  //==== 2a) Jet-by-jet updating of the b-tagging status
  bool IsBTagged_2a(JetTagging::Parameters jtp, const Jet& jet, string Syst="central");

};

#endif
