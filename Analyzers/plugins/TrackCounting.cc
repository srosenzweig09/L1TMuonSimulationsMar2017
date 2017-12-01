#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <cstdint>

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "L1Trigger/L1TMuonEndCap/interface/TrackTools.h"


// From L1Trigger/L1TMuonEndCap/interface/MuonTriggerPrimitive.h
class TriggerPrimitive {
public:
  enum subsystem_type{kDT,kCSC,kRPC,kGEM,kNSubsystems};
};


// _____________________________________________________________________________
class TrackCounting : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackCounting(const edm::ParameterSet& iConfig);
  ~TrackCounting();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  virtual void endJob() override;

  // Main functions
  void process();

  // Aux functions
  void getHandles(const edm::Event& iEvent);

  void bookHistograms();
  void writeHistograms();

  // Configurables
  const edm::InputTag emuHitTag_;
  const edm::InputTag emuTrackTag_;
  const edm::InputTag gmtMuonTag_;
  const edm::InputTag genPartTag_;
  const std::string   outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<l1t::EMTFHitCollection>          emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection>        emuTrackToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection>           gmtMuonToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>     genPartToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;

  l1t::EMTFHitCollection      emuHits_;
  l1t::EMTFTrackCollection    emuTracks_;
  //l1t::MuonBxCollection       gmtMuons_;
  std::vector<l1t::Muon>      gmtMuons_;
  reco::GenParticleCollection genParts_;
  int nPV1_, nPV2_;

  std::map<TString, TH1F*> histograms_;
  std::map<TString, TH2F*> histogram2Ds_;
};

// _____________________________________________________________________________
TrackCounting::TrackCounting(const edm::ParameterSet& iConfig) :
    emuHitTag_    (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_  (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    gmtMuonTag_   (iConfig.getParameter<edm::InputTag>("gmtMuonTag")),
    genPartTag_   (iConfig.getParameter<edm::InputTag>("genPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");

  emuHitToken_   = consumes<l1t::EMTFHitCollection>         (emuHitTag_);
  emuTrackToken_ = consumes<l1t::EMTFTrackCollection>       (emuTrackTag_);
  gmtMuonToken_  = consumes<l1t::MuonBxCollection>          (gmtMuonTag_);
  genPartToken_  = consumes<reco::GenParticleCollection>    (genPartTag_);
  puInfoToken_   = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));
}

TrackCounting::~TrackCounting() {}

void TrackCounting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent);
  process();
}

// _____________________________________________________________________________
void TrackCounting::getHandles(const edm::Event& iEvent) {

  // EMTF hits and tracks
  edm::Handle<decltype(emuHits_)>   emuHits_handle;
  edm::Handle<decltype(emuTracks_)> emuTracks_handle;

  if (!emuHitToken_.isUninitialized()) {
    iEvent.getByToken(emuHitToken_, emuHits_handle);
  }
  if (!emuHits_handle.isValid()) {
    edm::LogError("TrackCounting") << "Cannot get the product: " << emuHitTag_;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("TrackCounting") << "Cannot get the product: " << emuTrackTag_;
  }

  // GMT muons
  //edm::Handle<decltype(gmtMuons_)> gmtMuons_handle;
  edm::Handle<l1t::MuonBxCollection> gmtMuons_handle;

  if (!gmtMuonToken_.isUninitialized()) {
    iEvent.getByToken(gmtMuonToken_, gmtMuons_handle);
  }
  if (!gmtMuons_handle.isValid()) {
    edm::LogError("TrackCounting") << "Cannot get the product: " << gmtMuonTag_;
  }

  // Pileup summary info
  edm::Handle<std::vector<PileupSummaryInfo> > puInfos_handle;

  if (!iEvent.isRealData()) {  // MC
    nPV1_ = -1;
    nPV2_ = -1;

    if (!puInfoToken_.isUninitialized()) {
      iEvent.getByToken(puInfoToken_, puInfos_handle);
    }
    if (!puInfos_handle.isValid()) {
      edm::LogError("TrackCounting") << "Cannot get the product: " << edm::InputTag("addPileupInfo");
    }

    for (std::vector<PileupSummaryInfo>::const_iterator it = puInfos_handle->begin(); it != puInfos_handle->end(); ++it) {
      int bx = it->getBunchCrossing();
      if (bx == 0) {
        nPV1_ = it->getTrueNumInteractions();
        nPV2_ = it->getPU_NumInteractions();
      }
    }
    assert(nPV1_ != -1);
    assert(nPV2_ != -1);

  } else {  // Data
    nPV1_ = 0;  //FIXME
    nPV2_ = 0;  //FIXME
  }


  // ___________________________________________________________________________
  // Object filters

  //emuHits_.clear();
  //for (const auto& hit : (*emuHits_handle)) {
  //  if (!(-1 <= hit.BX() && hit.BX() <= 1))  continue;  // only BX=[-1,+1]
  //  //if (hit.Endcap() != 1)  continue;  // only positive endcap
  //  emuHits_.push_back(hit);
  //}

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    //if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }

  gmtMuons_.clear();
  for (int ibx = gmtMuons_handle->getFirstBX(); ibx <= gmtMuons_handle->getLastBX(); ++ibx) {
    if (ibx != 0)  continue;  // only BX=0
    for (l1t::MuonBxCollection::const_iterator it = gmtMuons_handle->begin(ibx); it != gmtMuons_handle->end(ibx); ++it) {
      if (!(it->pt() > 0.))  continue;  // only valid muons
      gmtMuons_.push_back(*it);
    }
  }
}

// _____________________________________________________________________________
// Functions

auto modify_ptcut = [](const auto& h) {
  int nbins = h->GetNbinsX();
  int entries = 0;
  for (int i = nbins; i > 0; i--) {
    entries += h->GetBinContent(i);
    h->SetBinContent(i, entries);
  }
  h->SetEntries(h->GetEntries() - nbins);
};


// _____________________________________________________________________________
void TrackCounting::process() {
  TString hname;
  TH1F* h;

  // number of events
  {
    hname = "nevents";
    histograms_.at(hname)->Fill(1.0);
  }

  // number of vertices
  {
    hname = "nvertices";
    histograms_.at(hname)->Fill(nPV1_);
    hname = "nvertices_a";
    histograms_.at(hname)->Fill(nPV2_);
  }

  // uGMT
  {
    decltype(gmtMuons_) mytracks;
    std::vector<bool> eta_bins(10+2, false);
    double highest_pt = 0.;
    //auto get_pt = [](const auto& x) { return x->pt(); };
    auto get_pt = [](const auto& x) { return std::min(100.-1e-3, (double) x->pt()); };
    auto cmp_pt = [](const auto& lhs, const auto& rhs) { return (lhs.pt() < rhs.pt()); };

    hname = "highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5) && (trk.hwQual() >= 12); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.1) && (trk.hwQual() >= 12); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.eta()) && std::abs(trk.eta()) <= 0.83) && (trk.hwQual() >= 12); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0.83 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 1.24) && (trk.hwQual() >= 12); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (1.24 < std::abs(trk.eta()) && std::abs(trk.eta()) <= 2.5) && (trk.hwQual() >= 12); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "muon_ptmin10_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.pt() > 10.) && (trk.hwQual() >= 12); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "muon_ptmin14_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.pt() > 14.) && (trk.hwQual() >= 12); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "muon_ptmin20_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.pt() > 20.) && (trk.hwQual() >= 12); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "muon_ptmin22_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(gmtMuons_.begin(), gmtMuons_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.pt() > 22.) && (trk.hwQual() >= 12); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }
  }

  // EMTF
  {
    decltype(emuTracks_) mytracks;
    std::vector<bool> eta_bins(10+2, false);
    double highest_pt = 0.;
    auto get_pt = [](const auto& x) { return std::min(100.-1e-3, (double) x->Pt()); };
    auto cmp_pt = [](const auto& lhs, const auto& rhs) { return (lhs.Pt() < rhs.Pt()); };

    hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 2.5) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_emtf_absEtaMin1.24_absEtaMax1.64_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (1.24 < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 1.64) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_emtf_absEtaMin1.64_absEtaMax2.14_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (1.64 < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 2.14) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_emtf_absEtaMin2.14_absEtaMax2.5_qmin12_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (2.14 < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 2.5) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin8_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 2.5) && (trk.Mode() == 7 || trk.Mode() == 10 || trk.Mode() == 12 || trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin4_pt";
    h = histograms_.at(hname);
    mytracks.clear();
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (0. < std::abs(trk.Eta()) && std::abs(trk.Eta()) <= 2.5) && (trk.Mode() == 3 || trk.Mode() == 5 || trk.Mode() == 6 || trk.Mode() == 9 || trk.Mode() == 7 || trk.Mode() == 10 || trk.Mode() == 12 || trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    if (!mytracks.empty())  highest_pt = get_pt(std::max_element(mytracks.begin(), mytracks.end(), cmp_pt) );
    if (!mytracks.empty())  h->Fill((highest_pt));

    hname = "emtf_ptmin10_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.Pt() > 10.) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.Eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "emtf_ptmin14_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.Pt() > 14.) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.Eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "emtf_ptmin20_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.Pt() > 20.) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.Eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }

    hname = "emtf_ptmin22_qmin12_eta";
    h = histograms_.at(hname);
    mytracks.clear();
    std::fill(eta_bins.begin(), eta_bins.end(), false);
    std::copy_if(emuTracks_.begin(), emuTracks_.end(), std::back_inserter(mytracks), [](const auto& trk) { return (trk.Pt() > 22.) && (trk.Mode() == 11 || trk.Mode() == 13 || trk.Mode() == 14 || trk.Mode() == 15); });
    for (const auto& trk : mytracks) { int bin = h->FindFixBin(std::abs(trk.Eta())); eta_bins.at(bin) = true; }
    for (unsigned bin=0; bin < eta_bins.size(); ++bin) { if (eta_bins.at(bin))  h->Fill(h->GetBinCenter(bin)); }
  }

}

// _____________________________________________________________________________
void TrackCounting::beginJob() {
  bookHistograms();
}

void TrackCounting::endJob() {
  writeHistograms();
}

// _____________________________________________________________________________
void TrackCounting::bookHistograms() {
  TString hname;

  TH1::SetDefaultSumw2();

  // number of events
  hname = "nevents";
  histograms_[hname] = new TH1F(hname, "; count", 5, 0, 5);

  // number of vertices
  hname = "nvertices";
  histograms_[hname] = new TH1F(hname, "; # of vertices", 300, 0, 300);
  hname = "nvertices_a";
  histograms_[hname] = new TH1F(hname, "; # of vertices", 300, 0, 300);

  // uGMT
  hname = "highest_muon_absEtaMin0_absEtaMax2.5_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax2.1_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0_absEtaMax0.83_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin0.83_absEtaMax1.24_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_muon_absEtaMin1.24_absEtaMax2.5_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "muon_ptmin10_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "muon_ptmin14_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "muon_ptmin20_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "muon_ptmin22_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);

  // emtf
  hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_emtf_absEtaMin1.24_absEtaMax1.64_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_emtf_absEtaMin1.64_absEtaMax2.14_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_emtf_absEtaMin2.14_absEtaMax2.5_qmin12_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin8_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "highest_emtf_absEtaMin0_absEtaMax2.5_qmin4_pt";
  histograms_[hname] = new TH1F(hname, "; p_{T} [GeV]; entries", 100, 0., 100.);
  hname = "emtf_ptmin10_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "emtf_ptmin14_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "emtf_ptmin20_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
  hname = "emtf_ptmin22_qmin12_eta";
  histograms_[hname] = new TH1F(hname, "; |#eta|; entries", 10, 1.55, 2.55);
}

void TrackCounting::writeHistograms() {
  // TFileService
  edm::Service<TFileService> fs;

  TH1F* h;
  TH2F* h2;
  TString ts_hname;

  for (const auto& kv : histograms_) {
    h = fs->make<TH1F>(*kv.second);
    if (h) {}

    // Modify cutoff pT histograms
    ts_hname = kv.first;
    if (ts_hname.EndsWith("_pt")) {
      ts_hname = ts_hname(0, ts_hname.Sizeof()-4) + "_ptcut";
      h = kv.second;
      h = (TH1F*) h->Clone(ts_hname.Data());
      h = fs->make<TH1F>(*h);
      modify_ptcut(h);
    }
  }

  for (const auto& kv : histogram2Ds_) {
    h2 = fs->make<TH2F>(*kv.second);
    if (h2) {}
  }
}

// _____________________________________________________________________________
void TrackCounting::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(TrackCounting);
