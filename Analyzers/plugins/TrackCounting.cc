#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <cstdint>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"

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
  const edm::InputTag genPartTag_;
  const std::string outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<l1t::EMTFHitCollection>   emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> emuTrackToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartToken_;

  l1t::EMTFHitCollection    emuHits_;
  l1t::EMTFTrackCollection  emuTracks_;
  reco::GenParticleCollection genParts_;

  std::map<TString, TH1F*> histograms_;
  std::map<TString, TH2F*> histogram2Ds_;
};

// _____________________________________________________________________________
TrackCounting::TrackCounting(const edm::ParameterSet& iConfig) :
    emuHitTag_    (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_  (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    genPartTag_   (iConfig.getParameter<edm::InputTag>("genPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");

  emuHitToken_   = consumes<l1t::EMTFHitCollection>     (emuHitTag_);
  emuTrackToken_ = consumes<l1t::EMTFTrackCollection>   (emuTrackTag_);
  genPartToken_  = consumes<reco::GenParticleCollection>(genPartTag_);
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
    return;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("TrackCounting") << "Cannot get the product: " << emuTrackTag_;
    return;
  }

  // Object filters
  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    if (!(-1 <= hit.BX() && hit.BX() <= 1))  continue;  // only BX=[-1,+1]
    //if (hit.Endcap() != 1)  continue;  // only positive endcap
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    //if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }
}

// _____________________________________________________________________________
// Functions

// from RecoMuon/DetLayers/src/MuonCSCDetLayerGeometryBuilder.cc
//      RecoMuon/DetLayers/src/MuonRPCDetLayerGeometryBuilder.cc
//      RecoMuon/DetLayers/src/MuonGEMDetLayerGeometryBuilder.cc
auto isFront_detail = [](int subsystem, int station, int ring, int chamber, int subsector) {
  bool result = false;

  if (subsystem == TriggerPrimitive::kCSC) {
    bool isOverlapping = !(station == 1 && ring == 3);
    // not overlapping means back
    if(isOverlapping)
    {
      bool isEven = (chamber % 2 == 0);
      // odd chambers are bolted to the iron, which faces
      // forward in 1&2, backward in 3&4, so...
      result = (station < 3) ? isEven : !isEven;
    }
  } else if (subsystem == TriggerPrimitive::kRPC) {
    // 10 degree rings have even subsectors in front
    // 20 degree rings have odd subsectors in front
    bool is_10degree = !((station == 3 || station == 4) && (ring == 1));
    bool isEven = (subsector % 2 == 0);
    result = (is_10degree) ? isEven : !isEven;
  } else if (subsystem == TriggerPrimitive::kGEM) {
    //
    result = (chamber % 2 == 0);
  }
  return result;
};

auto isFront = [](const auto& hit) {
  return isFront_detail(hit.Subsystem(), hit.Station(), hit.Ring(), hit.Chamber(), (hit.Subsystem() == TriggerPrimitive::kRPC ? hit.Subsector_RPC() : hit.Subsector()));
};

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

  std::random_device rd;
  std::mt19937 genrd(rd());

  auto get_mode_bin = [](const auto& trk) {
    int mode      = trk.Mode();
    assert(0 < mode && mode <= 15);
    if (mode == 15)                              return 3;
    if (mode == 11 || mode == 13 || mode == 14)  return 2;
    if (mode ==  7 || mode == 10 || mode == 12)  return 1;
    if (mode ==  3 || mode ==  5 || mode ==  6 || mode == 9)  return 0;
    return -1;
  };

  auto is_in_eta_range = [](const auto& trk) {
    double absEta = std::abs(trk.Eta());
    return (1.64 <= absEta && absEta <= 2.14);
  };

  auto is_in_bx_range =  [](const auto& trk) {
    return (trk.BX() == 0);
  };

  // ___________________________________________________________________________
  // Hits
  //for (const auto& hit : emuHits_) {
  //
  //}

  // Tracks
  bool found_leading_muon = false;

  for (const auto& trk : emuTracks_) {
    bool in_eta_range = is_in_eta_range(trk);
    bool in_bx_range  = is_in_bx_range(trk);
    int mode_bin      = get_mode_bin(trk);

    // Select SingleMu quality BX=0 in 1.64 < eta < 2.14 (only the leading muon)
    if (in_eta_range && in_bx_range && mode_bin >= 2 && !found_leading_muon) {
      found_leading_muon = true;

      hname = "trk_pt_csc";
      h = histograms_.at(hname);
      h->Fill(trk.Pt());

      hname = "trk_ptcut_csc";
      h = histograms_.at(hname);
      h->Fill(trk.Pt());

      bool pass = true;
      {
        // Apply gem-csc bending angle cut
        std::vector<l1t::EMTFHit> myhits0;
        const int sector = trk.Sector();
        std::copy_if(emuHits_.begin(), emuHits_.end(), std::back_inserter(myhits0), [sector](const auto& hit) { return (hit.PC_sector() == sector); });

        std::vector<l1t::EMTFHit> myhits1;  // GEM
        std::vector<l1t::EMTFHit> myhits2;  // CSC

        std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits1), [](const auto& hit) { return (hit.Station() == 1 && (hit.Ring() == 1 || hit.Ring() == 4) && hit.Subsystem() == TriggerPrimitive::kGEM); });
        std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits2), [](const auto& hit) { return (hit.Station() == 1 && (hit.Ring() == 1 || hit.Ring() == 4) && hit.Subsystem() == TriggerPrimitive::kCSC); });

        int idx = 0;
        int min_abs_dphi = 9999;
        int min_abs_dphi_idx = -1;
        for (const auto& hit : myhits2) {  // CSC
          int abs_dphi = std::abs(hit.Phi_fp() - trk.Phi_fp());
          if (min_abs_dphi > abs_dphi) {
            min_abs_dphi = abs_dphi;
            min_abs_dphi_idx = idx;
          }
          idx++;
        }

        if (min_abs_dphi_idx != -1) {
          const l1t::EMTFHit& myhit2 = myhits2.at(min_abs_dphi_idx);

          idx = 0;
          min_abs_dphi = 9999;
          min_abs_dphi_idx = -1;
          for (const auto& hit : myhits1) {  // GEM
            int abs_dphi = std::abs(hit.Phi_fp() - myhit2.Phi_fp());
            if (min_abs_dphi > abs_dphi) {
              min_abs_dphi = abs_dphi;
              min_abs_dphi_idx = idx;
            }
            idx++;
          }

          if (min_abs_dphi_idx != -1) {
            const l1t::EMTFHit& myhit1 = myhits1.at(min_abs_dphi_idx);
            float gem_csc_abs_bend = std::abs(emtf::range_phi_deg(myhit1.Phi_sim() - myhit2.Phi_sim()));
            bool is_front = isFront(myhit2);  // CSC

            // The cut values from Sven (even is front)
            // https://raw.githubusercontent.com/dildick/MuJetAnalysis/for-DisplacedMuonL1-PtAssignment-CMSSW-91X/DisplacedL1MuFilter/test/GEMCSCdPhiDict_wholeChamber.py
            //   'Pt5 ' : { 'odd' :  0.02112152, 'even' : 0.00948039 },
            //   'Pt7 ' : { 'odd' :  0.01460424, 'even' : 0.00664357 },
            //   'Pt10' : { 'odd' :  0.01001365, 'even' : 0.00463343 },
            //   'Pt15' : { 'odd' :  0.00666509, 'even' : 0.00317360 },
            //   'Pt20' : { 'odd' :  0.00506147, 'even' : 0.00251524 },
            //   'Pt30' : { 'odd' :  0.00352464, 'even' : 0.00193779 },
            //   'Pt40' : { 'odd' :  0.00281599, 'even' : 0.00169830 }
            double the_cut = !is_front ? 0.03 : 0.015;
            if (trk.Pt() > 5. )  the_cut = !is_front ? 0.02112152 : 0.00948039;
            if (trk.Pt() > 7. )  the_cut = !is_front ? 0.01460424 : 0.00664357;
            if (trk.Pt() > 10.)  the_cut = !is_front ? 0.01001365 : 0.00463343;
            if (trk.Pt() > 15.)  the_cut = !is_front ? 0.00666509 : 0.00317360;
            if (trk.Pt() > 20.)  the_cut = !is_front ? 0.00506147 : 0.00251524;
            if (trk.Pt() > 30.)  the_cut = !is_front ? 0.00352464 : 0.00193779;
            if (trk.Pt() > 40.)  the_cut = !is_front ? 0.00281599 : 0.00169830;
            the_cut = emtf::rad_to_deg(the_cut);  // convert to degrees

            pass = (gem_csc_abs_bend <= the_cut);
          } else {
            pass = false;
          }

        } else {
          pass = false;
        }
      }  // end gem-csc bending angle cut

      if (pass) {
        hname = "trk_pt_gem";
        h = histograms_.at(hname);
        h->Fill(trk.Pt());

        hname = "trk_ptcut_gem";
        h = histograms_.at(hname);
        h->Fill(trk.Pt());
      }
    }  // end track selection
  }  // end loop over tracks

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
  TH1F* h;

  hname = "trk_pt_csc";
  h = new TH1F(hname, "; EMTFv5 p_{T} [GeV]; entries", 400, 5., 205.);
  histograms_[hname] = h;

  hname = "trk_ptcut_csc";
  h = new TH1F(hname, "; EMTFv5 cutoff p_{T} [GeV]; entries", 400, 5., 205.);
  histograms_[hname] = h;

  hname = "trk_pt_gem";
  h = new TH1F(hname, "; EMTFv5 p_{T} [GeV]; entries", 400, 5., 205.);
  histograms_[hname] = h;

  hname = "trk_ptcut_gem";
  h = new TH1F(hname, "; EMTFv5 cutoff p_{T} [GeV]; entries", 400, 5., 205.);
  histograms_[hname] = h;
}

void TrackCounting::writeHistograms() {
  // TFileService
  edm::Service<TFileService> fs;

  TH1F* h;
  TH2F* h2;

  for (const auto& kv : histograms_) {
    // Modify cutoff pT histograms
    if (kv.first == "trk_ptcut_csc" || kv.first == "trk_ptcut_gem") {
      modify_ptcut(kv.second);
    }

    h = fs->make<TH1F>(*kv.second);
    if (h) {}
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
