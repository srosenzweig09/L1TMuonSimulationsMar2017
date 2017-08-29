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
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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


  //// Loop over sim tracks
  //for (const auto& part : (*simTracks_handle)) {
  //
  //}

  // Loop over tracking particles
  for (const auto& part : (*trkParts_handle)) {
    if (!(part.pt() >= 2.))     continue;  // only pT > 2
    //if (!(1.2 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap

    // Signal event
    //bool signal = (part.eventId().event() == 0);
    // In time bunch-crossing
    bool intime = (part.eventId().bunchCrossing() == 0);
    // Primary+charged: pT > 0.2 GeV, |eta| < 2.5, |rho0| < 0.5 cm, |z0| < 30 cm
    bool primary = (part.charge() != 0 && part.pt() > 0.2 && std::abs(part.eta()) < 2.5 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 0.5 && std::abs(part.vz()) < 30.0);
    bool is_muon = (std::abs(part.pdgId()) == 13);
    //if (!signal)  continue;
    if (!intime)  continue;
    if (!primary) continue;
    if (!is_muon) continue;

    double pt = part.pt();
    double invPt = 1.0 / pt;
    double absEta = std::abs(part.eta());

    if (1.24 <= absEta && absEta <= 1.64) {
      histograms_["muon_absEtaMin1.64_absEtaMax2.14_pt"]->Fill(pt);
      histograms_["muon_absEtaMin1.64_absEtaMax2.14_invPt"]->Fill(invPt);
    }

    histogram2Ds_["muon_pt_vs_eta"]->Fill(absEta, pt);
    histogram2Ds_["muon_invPt_vs_eta"]->Fill(absEta, invPt);
  }

}

void MinBiasMuonAnalyzer::beginJob() {
  bookHistograms();
}

void MinBiasMuonAnalyzer::endJob() {
  writeHistograms();
}

void MinBiasMuonAnalyzer::bookHistograms() {
  TString hname;

  // TFileService
  edm::Service<TFileService> fs;

  // TH1F
  hname = "muon_absEtaMin1.64_absEtaMax2.14_pt";
  histograms_[hname] = fs->make<TH1F>(hname, "; p_{T} [GeV]", 200, 0, 200);

  hname = "muon_absEtaMin1.64_absEtaMax2.14_invPt";
  histograms_[hname] = fs->make<TH1F>(hname, "; 1/p_{T} [1/GeV]", 200, 0, 0.5);

  // TH2F
  hname = "muon_pt_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |eta|; p_{T} [GeV]", 30, 0., 3.0, 200, 0, 200);

  hname = "muon_invPt_vs_eta";
  histogram2Ds_[hname] = fs->make<TH2F>(hname, "; |eta|; 1/p_{T} [1/GeV]", 30, 0., 3.0, 200, 0, 0.5);
}

void MinBiasMuonAnalyzer::writeHistograms() {
  // do nothing
}

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
