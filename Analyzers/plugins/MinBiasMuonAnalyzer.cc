#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"


// _____________________________________________________________________________
class MinBiasMuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MinBiasMuonAnalyzer(const edm::ParameterSet&);
  ~MinBiasMuonAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // Main functions
  void process();

  // Aux functions
  void getHandles(const edm::Event& iEvent);

  void bookHistograms();
  void writeHistograms();

  // Configurables
  const edm::InputTag   simTrackTag_;
  const edm::InputTag   trkPartTag_;
  const std::string     outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<edm::SimTrackContainer>      simTrackToken_;
  edm::EDGetTokenT<TrackingParticleCollection>  trkPartToken_;

  edm::SimTrackContainer      simTracks_;
  TrackingParticleCollection  trkParts_;

  std::map<TString, TH1F*> histograms_;
  std::map<TString, TH2F*> histogram2Ds_;
};

// _____________________________________________________________________________
//static std::random_device rd;
//static std::mt19937 genrd(rd());
static std::mt19937 genrd0(20230);
static std::mt19937 genrd1(20231);


MinBiasMuonAnalyzer::MinBiasMuonAnalyzer(const edm::ParameterSet& iConfig) :
    simTrackTag_  (iConfig.getParameter<edm::InputTag>("simTrackTag")),
    trkPartTag_   (iConfig.getParameter<edm::InputTag>("trkPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");

  simTrackToken_ = consumes<edm::SimTrackContainer>     (simTrackTag_);
  trkPartToken_  = consumes<TrackingParticleCollection> (trkPartTag_);
}

MinBiasMuonAnalyzer::~MinBiasMuonAnalyzer() {}

void MinBiasMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  getHandles(iEvent);
  process();
}

// _____________________________________________________________________________
void MinBiasMuonAnalyzer::getHandles(const edm::Event& iEvent) {

  // Sim Tracks
  edm::Handle<decltype(simTracks_)> simTracks_handle;

  if (!iEvent.isRealData()) {
    if (!simTrackToken_.isUninitialized()) {
      iEvent.getByToken(simTrackToken_, simTracks_handle);
    }
    if (!simTracks_handle.isValid()) {
      edm::LogError("MinBiasMuonAnalyzer") << "Cannot get the product: " << simTrackTag_;
    }
  }

  // Tracking particles
  edm::Handle<decltype(trkParts_)> trkParts_handle;

  if (!iEvent.isRealData()) {
    if (!trkPartToken_.isUninitialized()) {
      iEvent.getByToken(trkPartToken_, trkParts_handle);
    }
    if (!trkParts_handle.isValid()) {
      edm::LogError("MinBiasMuonAnalyzer") << "Cannot get the product: " << trkPartTag_;
    }
  }

  // ___________________________________________________________________________
  // Object filters

  //simTracks_.clear();
  //for (const auto& trk : (*simTracks_handle)) {
  //  simTracks_.push_back(trk);
  //}

  trkParts_.clear();
  for (const auto& part : (*trkParts_handle)) {
    if (!(part.pt() >= 2.))     continue;  // only pT > 2
    //if (!(1.2 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap
    if (!(std::abs(part.pdgId()) == 13))  continue;  // only muons

    // Plot dxy, dz before selection (pT > 2 GeV, z0 < 50 cm)
    if (part.pt() >= 2. && std::abs(part.vz()) < 50.0) {
      double dxy = std::sqrt(part.vx() * part.vx() + part.vy() * part.vy());
      double dz = std::abs(part.vz());
      histograms_.at("muon_ptmin2_absEtaMin0_absEtaMax2.5_dxy")->Fill(dxy);
      histograms_.at("muon_ptmin2_absEtaMin0_absEtaMax2.5_dz")->Fill(dz);
      histogram2Ds_.at("muon_pt_vs_dxy")->Fill(dxy, part.pt());
      histogram2Ds_.at("muon_pt_vs_dz")->Fill(dz, part.pt());
    }

    {
      // Selection

      // Signal event
      //bool signal = (part.eventId().event() == 0);

      // In time bunch-crossing
      bool intime = (part.eventId().bunchCrossing() == 0);

      // Primary+charged: pT > 0.2 GeV, |eta| < 2.5, |rho0| < 0.5 cm, |z0| < 30 cm
      //bool primary = (part.charge() != 0 && part.pt() > 0.2 && std::abs(part.eta()) < 2.5 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 0.5 && std::abs(part.vz()) < 30.0);

      // Primary+secondary
      bool secondary = (part.charge() != 0 && part.pt() > 0.2 && std::abs(part.eta()) < 2.5 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 15.0 && std::abs(part.vz()) < 50.0);

      //if (!signal)  continue;
      if (!intime)  continue;
      //if (!primary) continue;
      if (!secondary) continue;
    }

    trkParts_.push_back(part);
  }

}

// _____________________________________________________________________________
void MinBiasMuonAnalyzer::process() {
  TString hname;
  TH1F* h;

  // number of events
  {
    hname = "nevents";
    histograms_.at(hname)->Fill(1.0);
  }

  //// Loop over sim tracks
  //for (const auto& trk : simTracks_) {
  //
  //}

  // Loop over tracking particles
  for (const auto& part : trkParts_) {
    double pt = part.pt();
    double absEta = std::abs(part.eta());

    double invPt = 1.0 / pt;
    double invPt2 = 1.0 / pt / pt;
    double invPt3 = 1.0 / pt / pt / pt;
    double invPt4 = 1.0 / pt / pt / pt / pt;
    double invPt5 = 1.0 / pt / pt / pt / pt / pt;
    double logPt = std::log2(pt);

    if (1.24 < absEta && absEta <= 2.5) {
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_pt")->Fill(pt);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_invPt")->Fill(invPt);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_invPt2")->Fill(invPt2);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_invPt3")->Fill(invPt3);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_invPt4")->Fill(invPt4);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_invPt5")->Fill(invPt5);
      histograms_.at("muon_absEtaMin1.24_absEtaMax2.5_logPt")->Fill(logPt);
    }

    if (1.24 < absEta && absEta <= 1.64) {
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_pt")->Fill(pt);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_invPt")->Fill(invPt);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_invPt2")->Fill(invPt2);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_invPt3")->Fill(invPt3);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_invPt4")->Fill(invPt4);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_invPt5")->Fill(invPt5);
      histograms_.at("muon_absEtaMin1.24_absEtaMax1.64_logPt")->Fill(logPt);
    } else if (1.64 < absEta && absEta <= 2.14) {
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_pt")->Fill(pt);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_invPt")->Fill(invPt);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_invPt2")->Fill(invPt2);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_invPt3")->Fill(invPt3);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_invPt4")->Fill(invPt4);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_invPt5")->Fill(invPt5);
      histograms_.at("muon_absEtaMin1.64_absEtaMax2.14_logPt")->Fill(logPt);
    } else if (2.14 < absEta && absEta <= 2.5) {
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_pt")->Fill(pt);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_invPt")->Fill(invPt);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_invPt2")->Fill(invPt2);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_invPt3")->Fill(invPt3);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_invPt4")->Fill(invPt4);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_invPt5")->Fill(invPt5);
      histograms_.at("muon_absEtaMin2.14_absEtaMax2.5_logPt")->Fill(logPt);
    }

    histogram2Ds_.at("muon_pt_vs_eta")->Fill(absEta, pt);
    histogram2Ds_.at("muon_invPt_vs_eta")->Fill(absEta, invPt);
    histogram2Ds_.at("muon_invPt2_vs_eta")->Fill(absEta, invPt2);
    histogram2Ds_.at("muon_invPt3_vs_eta")->Fill(absEta, invPt3);
    histogram2Ds_.at("muon_invPt4_vs_eta")->Fill(absEta, invPt4);
    histogram2Ds_.at("muon_invPt5_vs_eta")->Fill(absEta, invPt5);
    histogram2Ds_.at("muon_logPt_vs_eta")->Fill(absEta, logPt);
  }

  // Rate studies
  {
    decltype(trkParts_) myparts;
    //auto get_pt = [](const auto& x) { return x->pt(); };
    auto get_pt = [](const auto& x) { return std::min(100.-1e-3, (double) x->pt()); };
    auto cmp_pt = [](const auto& lhs, const auto& rhs) { return (lhs.pt() < rhs.pt()); };

    auto fill_pt = [&](auto&& myparts, auto&& h, auto get_pt, auto cmp_pt, auto selection) {
      double highest_pt = 0.;
      myparts.clear();
      std::copy_if(trkParts_.begin(), trkParts_.end(), std::back_inserter(myparts), selection);
      if (!myparts.empty())  highest_pt = get_pt(std::max_element(myparts.begin(), myparts.end(), cmp_pt) );
      if (!myparts.empty())  h->Fill((highest_pt));
    };

    hname = "highest_muon_absEtaMin0_absEtaMax2.5_mc_pt";
    fill_pt(myparts, histograms_.at(hname), get_pt, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5); });

    hname = "highest_muon_absEtaMin0_absEtaMax2.1_mc_pt";
    fill_pt(myparts, histograms_.at(hname), get_pt, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.1); });

    hname = "highest_muon_absEtaMin0_absEtaMax0.83_mc_pt";
    fill_pt(myparts, histograms_.at(hname), get_pt, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 0.83); });

    hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_mc_pt";
    fill_pt(myparts, histograms_.at(hname), get_pt, cmp_pt, [](const auto& trk) { return (0.83 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 1.24); });

    hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_mc_pt";
    fill_pt(myparts, histograms_.at(hname), get_pt, cmp_pt, [](const auto& trk) { return (1.24 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5); });


    std::normal_distribution<double> rng(0.,0.2);  // 20% pT resolution

    auto get_pt_fake = [](const auto& x, const auto& res) { return std::min(100.-1e-3, (double) (x->pt() * (1.0+res))); };

    auto fill_pt_fake = [&](auto&& myparts, auto&& h, auto get_pt_fake, auto cmp_pt, auto selection) {
      double highest_pt = 0.;
      double res = std::max(0., rng(genrd0));
      myparts.clear();
      std::copy_if(trkParts_.begin(), trkParts_.end(), std::back_inserter(myparts), selection);
      if (!myparts.empty())  highest_pt = get_pt_fake(std::max_element(myparts.begin(), myparts.end(), cmp_pt), res);  // CUIDADO: only smear with pT resolution after sorting
      if (!myparts.empty())  h->Fill((highest_pt));
    };

    hname = "highest_muon_absEtaMin0_absEtaMax2.5_fake_pt";
    fill_pt_fake(myparts, histograms_.at(hname), get_pt_fake, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5); });

    hname = "highest_muon_absEtaMin0_absEtaMax2.1_fake_pt";
    fill_pt_fake(myparts, histograms_.at(hname), get_pt_fake, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.1); });

    hname = "highest_muon_absEtaMin0_absEtaMax0.83_fake_pt";
    fill_pt_fake(myparts, histograms_.at(hname), get_pt_fake, cmp_pt, [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 0.83); });

    hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_fake_pt";
    fill_pt_fake(myparts, histograms_.at(hname), get_pt_fake, cmp_pt, [](const auto& trk) { return (0.83 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 1.24); });

    hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_fake_pt";
    fill_pt_fake(myparts, histograms_.at(hname), get_pt_fake, cmp_pt, [](const auto& trk) { return (1.24 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5); });

  }

  if (h) {}

}

// _____________________________________________________________________________
void MinBiasMuonAnalyzer::beginJob() {
  bookHistograms();
}

void MinBiasMuonAnalyzer::endJob() {
  writeHistograms();
}

// _____________________________________________________________________________
void MinBiasMuonAnalyzer::bookHistograms() {
  TString hname;

  // TFileService
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2();

  // TH1F
  hname = "muon_ptmin2_absEtaMin0_absEtaMax2.5_dxy";
  histograms_[hname] = fs->make<TH1F>(hname, "; |d_{xy}| [cm]", 200, 0, 100);
  hname = "muon_ptmin2_absEtaMin0_absEtaMax2.5_dz";
  histograms_[hname] = fs->make<TH1F>(hname, "; |d_{z}| [cm]", 200, 0, 100);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]", 200, 0, 200);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]", 200, 0, 200);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]", 200, 0, 200);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]", 200, 0, 200);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_invPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T} [1/GeV]", 200, 0, 0.5);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_invPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T} [1/GeV]", 200, 0, 0.5);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T} [1/GeV]", 200, 0, 0.5);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_invPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T} [1/GeV]", 200, 0, 0.5);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_invPt2";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{2} [1/GeV^{2}]", 200, 0, 0.25);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_invPt2";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{2} [1/GeV^{2}]", 200, 0, 0.25);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt2";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{2} [1/GeV^{2}]", 200, 0, 0.25);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_invPt2";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{2} [1/GeV^{2}]", 200, 0, 0.25);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_invPt3";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{3} [1/GeV^{3}]", 200, 0, 0.125);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_invPt3";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{3} [1/GeV^{3}]", 200, 0, 0.125);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt3";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{3} [1/GeV^{3}]", 200, 0, 0.125);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_invPt3";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{3} [1/GeV^{3}]", 200, 0, 0.125);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_invPt4";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{4} [1/GeV^{4}]", 200, 0, 0.0625);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_invPt4";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{4} [1/GeV^{4}]", 200, 0, 0.0625);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt4";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{4} [1/GeV^{4}]", 200, 0, 0.0625);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_invPt4";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{4} [1/GeV^{4}]", 200, 0, 0.0625);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_invPt5";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{5} [1/GeV^{5}]", 200, 0, 0.03125);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_invPt5";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{5} [1/GeV^{5}]", 200, 0, 0.03125);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt5";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{5} [1/GeV^{5}]", 200, 0, 0.03125);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_invPt5";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T}^{5} [1/GeV^{5}]", 200, 0, 0.03125);

  hname = "muon_absEtaMin1.24_absEtaMax2.5_logPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; log_{2}(p_{T}) [log_{2} GeV]", 200, 1, 8);
  hname = "muon_absEtaMin1.24_absEtaMax1.64_logPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; log_{2}(p_{T}) [log_{2} GeV]", 200, 1, 8);
  hname = "muon_absEtaMin1.64_absEtaMax2.14_logPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; log_{2}(p_{T}) [log_{2} GeV]", 200, 1, 8);
  hname = "muon_absEtaMin2.14_absEtaMax2.5_logPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; log_{2}(p_{T}) [log_{2} GeV]", 200, 1, 8);

  hname = "highest_muon_absEtaMin0_absEtaMax2.5_mc_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax2.1_mc_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax0.83_mc_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_mc_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_mc_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);

  hname = "highest_muon_absEtaMin0_absEtaMax2.5_fake_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax2.1_fake_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax0.83_fake_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_fake_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_fake_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);


  // TH2F
  hname = "muon_pt_vs_dxy";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |d_{xy}| [cm]; p_{T} [GeV]", 40, 0., 100, 200, 0, 200);
  hname = "muon_pt_vs_dz";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |d_{z}| [cm]; p_{T} [GeV]", 40, 0., 100, 200, 0, 200);

  hname = "muon_pt_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; p_{T} [GeV]", 30, 0., 3.0, 200, 0, 200);

  hname = "muon_invPt_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; 1/p_{T} [1/GeV]", 30, 0., 3.0, 200, 0, 0.5);

  hname = "muon_invPt2_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; 1/p_{T}^{2} [1/GeV^{2}]", 30, 0., 3.0, 200, 0, 0.25);

  hname = "muon_invPt3_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; 1/p_{T}^{3} [1/GeV^{3}]", 30, 0., 3.0, 200, 0, 0.125);

  hname = "muon_invPt4_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; 1/p_{T}^{4} [1/GeV^{4}]", 30, 0., 3.0, 200, 0, 0.0625);

  hname = "muon_invPt5_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; 1/p_{T}^{5} [1/GeV^{5}]", 30, 0., 3.0, 200, 0, 0.03125);

  hname = "muon_logPt_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |#eta|; log_{2}(p_{T}) [log_{2} GeV]", 30, 0., 3.0, 200, 1, 8);


  // bookkeeping
  hname = "nevents";
  histograms_[hname] = fs->make<TH1F>(hname, "; count", 5, 0, 5);
}

void MinBiasMuonAnalyzer::writeHistograms() {
  // do nothing
}

// _____________________________________________________________________________
void MinBiasMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(MinBiasMuonAnalyzer);
