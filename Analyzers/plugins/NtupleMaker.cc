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
//#include "TH1F.h"
//#include "TH2F.h"
//#include "TEfficiency.h"
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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "L1Trigger/L1TMuonEndCap/interface/TrackTools.h"

#include "L1TMuonSimulations/Analyzers/interface/EMTFMCTruth.h"


// From L1Trigger/L1TMuonEndCap/interface/MuonTriggerPrimitive.h
class TriggerPrimitive {
public:
  enum subsystem_type{kDT,kCSC,kRPC,kGEM,kNSubsystems};
};

// From L1Trigger/L1TMuonEndCap/interface/TTMuonTriggerPrimitive.h
class TTTriggerPrimitive {
public:
  enum subsystem_type{kTT = 20, kNSubsystems};
};


// _____________________________________________________________________________
class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit NtupleMaker(const edm::ParameterSet& iConfig);
  ~NtupleMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  virtual void endJob() override;

  // Main functions
  void process();

  // Aux functions
  void getHandles(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  void makeTree();
  void writeTree();

  // Configurables
  EMTFMCTruth  truth_;  // MC truth for EMTFHit

  const edm::InputTag   emuHitTag_;
  const edm::InputTag   emuTrackTag_;
  const edm::InputTag   genPartTag_;
  const edm::InputTag   trkPartTag_;
  const std::string     outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<l1t::EMTFHitCollection>      emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection>    emuTrackToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartToken_;
  edm::EDGetTokenT<TrackingParticleCollection>  trkPartToken_;

  l1t::EMTFHitCollection      emuHits_;
  l1t::EMTFTrackCollection    emuTracks_;
  reco::GenParticleCollection genParts_;
  TrackingParticleCollection  trkParts_;

  // TTree
  TTree* tree;

  // Hits
  std::unique_ptr<std::vector<int16_t> >  vh_endcap;
  std::unique_ptr<std::vector<int16_t> >  vh_station;
  std::unique_ptr<std::vector<int16_t> >  vh_ring;
  std::unique_ptr<std::vector<int16_t> >  vh_sector;
  std::unique_ptr<std::vector<int16_t> >  vh_subsector;
  std::unique_ptr<std::vector<int16_t> >  vh_chamber;
  std::unique_ptr<std::vector<int16_t> >  vh_cscid;
  std::unique_ptr<std::vector<int16_t> >  vh_bx;
  std::unique_ptr<std::vector<int16_t> >  vh_type;  // subsystem: DT=0,CSC=1,RPC=2,GEM=3
  std::unique_ptr<std::vector<int16_t> >  vh_neighbor;
  //
  std::unique_ptr<std::vector<int16_t> >  vh_strip;
  std::unique_ptr<std::vector<int16_t> >  vh_wire;
  std::unique_ptr<std::vector<int16_t> >  vh_roll;
  std::unique_ptr<std::vector<int16_t> >  vh_pattern;
  std::unique_ptr<std::vector<int16_t> >  vh_quality;
  std::unique_ptr<std::vector<int16_t> >  vh_bend;
  std::unique_ptr<std::vector<int16_t> >  vh_time;
  std::unique_ptr<std::vector<int16_t> >  vh_fr;
  //
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_phi;
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_theta;
  //
  std::unique_ptr<std::vector<float  > >  vh_sim_phi;
  std::unique_ptr<std::vector<float  > >  vh_sim_theta;
  std::unique_ptr<std::vector<float  > >  vh_sim_eta;
  std::unique_ptr<std::vector<float  > >  vh_sim_r;
  std::unique_ptr<std::vector<float  > >  vh_sim_z;
  std::unique_ptr<std::vector<int32_t> >  vh_sim_tp1;
  std::unique_ptr<std::vector<int32_t> >  vh_sim_tp2;
  //
  std::unique_ptr<int32_t              >  vh_size;

  // Tracks
  std::unique_ptr<std::vector<float  > >  vt_pt;
  std::unique_ptr<std::vector<float  > >  vt_xml_pt;
  std::unique_ptr<std::vector<float  > >  vt_phi;
  std::unique_ptr<std::vector<float  > >  vt_eta;
  std::unique_ptr<std::vector<float  > >  vt_theta;
  std::unique_ptr<std::vector<int16_t> >  vt_q;  // charge
  //
  std::unique_ptr<std::vector<int16_t> >  vt_mode;
  std::unique_ptr<std::vector<int16_t> >  vt_endcap;
  std::unique_ptr<std::vector<int16_t> >  vt_sector;
  std::unique_ptr<std::vector<int16_t> >  vt_bx;
  //
  std::unique_ptr<int32_t              >  vt_size;

  // Gen particles
  std::unique_ptr<std::vector<float  > >  vp_pt;
  std::unique_ptr<std::vector<float  > >  vp_phi;
  std::unique_ptr<std::vector<float  > >  vp_eta;
  std::unique_ptr<std::vector<float  > >  vp_theta;
  std::unique_ptr<std::vector<float  > >  vp_vx;
  std::unique_ptr<std::vector<float  > >  vp_vy;
  std::unique_ptr<std::vector<float  > >  vp_vz;
  std::unique_ptr<std::vector<int16_t> >  vp_q;  // charge
  std::unique_ptr<std::vector<int32_t> >  vp_pdgid;
  std::unique_ptr<int32_t              >  vp_size;
};


// _____________________________________________________________________________
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) :
    truth_        (iConfig, consumesCollector()),
    emuHitTag_    (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_  (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    genPartTag_   (iConfig.getParameter<edm::InputTag>("genPartTag")),
    trkPartTag_   (iConfig.getParameter<edm::InputTag>("trkPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");

  emuHitToken_   = consumes<l1t::EMTFHitCollection>     (emuHitTag_);
  emuTrackToken_ = consumes<l1t::EMTFTrackCollection>   (emuTrackTag_);
  genPartToken_  = consumes<reco::GenParticleCollection>(genPartTag_);
  trkPartToken_  = consumes<TrackingParticleCollection> (trkPartTag_);
}

NtupleMaker::~NtupleMaker() {}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent, iSetup);
  process();
}

// _____________________________________________________________________________
void NtupleMaker::getHandles(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // MC truth for EMTFHit
  truth_.initEvent(iEvent, iSetup);

  // EMTF hits and tracks
  edm::Handle<decltype(emuHits_)>   emuHits_handle;
  edm::Handle<decltype(emuTracks_)> emuTracks_handle;

  if (!emuHitToken_.isUninitialized()) {
    iEvent.getByToken(emuHitToken_, emuHits_handle);
  }
  if (!emuHits_handle.isValid()) {
    edm::LogError("NtupleMaker") << "Cannot get the product: " << emuHitTag_;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("NtupleMaker") << "Cannot get the product: " << emuTrackTag_;
  }

  // Gen particles
  edm::Handle<decltype(genParts_)> genParts_handle;

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      edm::LogError("NtupleMaker") << "Cannot get the product: " << genPartTag_;
    }
  }

  // Tracking particles
  edm::Handle<decltype(trkParts_)> trkParts_handle;

  if (!iEvent.isRealData()) {
    if (!trkPartToken_.isUninitialized()) {
      iEvent.getByToken(trkPartToken_, trkParts_handle);
    }
    if (!trkParts_handle.isValid()) {
      edm::LogError("NtupleMaker") << "Cannot get the product: " << trkPartTag_;
    }
  }


  // ___________________________________________________________________________
  // Object filters
  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    if (!(-1 <= hit.BX() && hit.BX() <= 1))  continue;  // only BX=[-1,+1]
    //if (hit.Endcap() != 1)  continue;  // only positive endcap
    if (hit.Subsystem() == TTTriggerPrimitive::kTT)  continue;  // ignore TTStubs
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    //if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }

  genParts_.clear();
  for (const auto& part : (*genParts_handle)) {
    if (!(part.pt() >= 2.))     continue;  // only pT > 2
    //if (!(1.2 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap
    genParts_.push_back(part);
  }

  trkParts_.clear();
  for (const auto& part : (*trkParts_handle)) {
    //if (!(part.pt() >= 2.))     continue;  // only pT > 2
    //if (!(1.2 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap
    if (!(std::abs(part.pdgId()) == 13))  continue;  // only muons

    // Tracking particle selection
    {
      // Signal event
      //bool signal = (part.eventId().event() == 0);

      // In time bunch-crossing
      //bool intime = (part.eventId().bunchCrossing() == 0);

      // In time + out of time bunch-crossing (-2 <= BX <= +2)
      bool outoftime = (-2 <= part.eventId().bunchCrossing() && part.eventId().bunchCrossing() <= +2);

      // Primary+charged: pT > 0.2 GeV, |eta| < 2.5, |rho0| < 0.5 cm, |z0| < 30 cm
      //bool primary = (part.charge() != 0 && part.pt() > 0.2 && std::abs(part.eta()) < 2.5 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 0.5 && std::abs(part.vz()) < 30.0);

      // Primary+secondary pT > 0.5 GeV, |eta| < 2.5, |rho0| < 120 cm, |z0| < 300 cm (tracker volume)
      bool secondary = (part.charge() != 0 && part.pt() > 0.5 && std::abs(part.eta()) < 2.5 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 120.0 && std::abs(part.vz()) < 300.0);

      //if (!signal)  continue;
      //if (!intime)  continue;
      if (!outoftime) continue;
      //if (!primary) continue;
      if (!secondary) continue;
    }

    trkParts_.push_back(part);
  }
}


// _____________________________________________________________________________
void NtupleMaker::process() {

  // In-place functions

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

  auto isFront = [&](const auto& hit) {
    return isFront_detail(hit.Subsystem(), hit.Station(), hit.Ring(), hit.Chamber(), (hit.Subsystem() == TriggerPrimitive::kRPC ? hit.Subsector_RPC() : hit.Subsector()));
  };

  auto prepare_sim_tp = [&]() {
    truth_.makeTrackingParticleLinks(trkParts_);
    return;
  };

  auto get_sim_tp1 = [&](const auto& hit) {
    if (hit.Subsystem() == TriggerPrimitive::kCSC) {
      return truth_.findCSCStripSimLink(hit, trkParts_);
    } else if (hit.Subsystem() == TriggerPrimitive::kRPC) {
      return truth_.findRPCDigiSimLink(hit, trkParts_);
    } else if (hit.Subsystem() == TriggerPrimitive::kGEM) {
      return truth_.findGEMDigiSimLink(hit, trkParts_);
    }
    return -1;
  };

  auto get_sim_tp2 = [&](const auto& hit) {
    if (hit.Subsystem() == TriggerPrimitive::kCSC) {
      return truth_.findCSCWireSimLink(hit, trkParts_);
    } else if (hit.Subsystem() == TriggerPrimitive::kRPC) {
      return truth_.findRPCDigiSimLink(hit, trkParts_);
    } else if (hit.Subsystem() == TriggerPrimitive::kGEM) {
      return truth_.findGEMDigiSimLink(hit, trkParts_);
    }
    return -1;
  };

  auto get_time = [](double time) {
    return static_cast<int>(std::round(time * 100));  // to integer unit of 0.01 ns (arbitrary)
  };


  // ___________________________________________________________________________
  bool please_use_trkParts = true;

  if (verbose_ > 0) {
    std::cout << "[DEBUG] # hits: " << emuHits_.size() << " #  tracks: " << emuTracks_.size() << " # gen parts: " << genParts_.size() << " # trk parts: " << trkParts_.size() << std::endl;
  }

  // Hits
  prepare_sim_tp();  // must be called before calling get_sim_tp1() and get_sim_tp2()

  for (const auto& hit : emuHits_) {
    vh_endcap     ->push_back(hit.Endcap());
    vh_station    ->push_back(hit.Station());
    vh_ring       ->push_back(hit.Ring());
    vh_sector     ->push_back(hit.PC_sector());
    vh_subsector  ->push_back(hit.Subsector());
    vh_chamber    ->push_back(hit.Chamber());
    vh_cscid      ->push_back(hit.CSC_ID());
    vh_bx         ->push_back(hit.BX());
    vh_type       ->push_back(hit.Subsystem());
    vh_neighbor   ->push_back(hit.Neighbor());
    //
    vh_strip      ->push_back(hit.Strip());
    vh_wire       ->push_back(hit.Wire());
    vh_roll       ->push_back(hit.Roll());
    vh_pattern    ->push_back(hit.Pattern());
    vh_quality    ->push_back(hit.Quality());
    vh_bend       ->push_back(hit.Bend());
    vh_time       ->push_back(get_time(hit.Time()));
    vh_fr         ->push_back(isFront(hit));
    //
    vh_emtf_phi   ->push_back(hit.Phi_fp());
    vh_emtf_theta ->push_back(hit.Theta_fp());
    //
    vh_sim_phi    ->push_back(hit.Phi_sim());
    vh_sim_theta  ->push_back(hit.Theta_sim());
    vh_sim_eta    ->push_back(hit.Eta_sim());
    vh_sim_r      ->push_back(hit.Rho_sim());
    vh_sim_z      ->push_back(hit.Z_sim());
    vh_sim_tp1    ->push_back(get_sim_tp1(hit));
    vh_sim_tp2    ->push_back(get_sim_tp2(hit));
  }
  (*vh_size) = emuHits_.size();

  // Tracks
  for (const auto& trk : emuTracks_) {
    vt_pt         ->push_back(trk.Pt());
    vt_xml_pt     ->push_back(trk.Pt_XML());
    vt_phi        ->push_back(trk.Phi_glob());
    vt_eta        ->push_back(trk.Eta());
    vt_theta      ->push_back(trk.Theta());
    vt_q          ->push_back(trk.Charge());
    //
    vt_mode       ->push_back(trk.Mode());
    vt_endcap     ->push_back(trk.Endcap());
    vt_sector     ->push_back(trk.Sector());
    vt_bx         ->push_back(trk.BX());
  }
  (*vt_size) = emuTracks_.size();

  // Gen particles
  if (!please_use_trkParts) {
    for (const auto& part : genParts_) {
      vp_pt         ->push_back(part.pt());
      vp_phi        ->push_back(part.phi());
      vp_eta        ->push_back(part.eta());
      vp_theta      ->push_back(part.theta());
      vp_vx         ->push_back(part.vx());
      vp_vy         ->push_back(part.vy());
      vp_vz         ->push_back(part.vz());
      vp_q          ->push_back(part.charge());
      vp_pdgid      ->push_back(part.pdgId());
    }
    (*vp_size) = genParts_.size();
    assert(static_cast<size_t>(*vp_size) == vp_pt->size());
    assert((*vp_size) <= 1);  // expect 0 or 1 gen particle
  }

  // Tracking particles
  if (please_use_trkParts) {
    for (const auto& part : trkParts_) {
      vp_pt         ->push_back(part.pt());
      vp_phi        ->push_back(part.phi());
      vp_eta        ->push_back(part.eta());
      vp_theta      ->push_back(part.theta());
      vp_vx         ->push_back(part.vx());
      vp_vy         ->push_back(part.vy());
      vp_vz         ->push_back(part.vz());
      vp_q          ->push_back(part.charge());
      vp_pdgid      ->push_back(part.pdgId());
    }
    (*vp_size) = trkParts_.size();
    assert(static_cast<size_t>(*vp_size) == vp_pt->size());
    //assert((*vp_size) <= 1);  // expect 0 or 1 gen particle
  }

  // Fill
  tree->Fill();

  // Clear
  vh_endcap     ->clear();
  vh_station    ->clear();
  vh_ring       ->clear();
  vh_sector     ->clear();
  vh_subsector  ->clear();
  vh_chamber    ->clear();
  vh_cscid      ->clear();
  vh_bx         ->clear();
  vh_type       ->clear();
  vh_neighbor   ->clear();
  //
  vh_strip      ->clear();
  vh_wire       ->clear();
  vh_roll       ->clear();
  vh_pattern    ->clear();
  vh_quality    ->clear();
  vh_bend       ->clear();
  vh_time       ->clear();
  vh_fr         ->clear();
  //
  vh_emtf_phi   ->clear();
  vh_emtf_theta ->clear();
  //
  vh_sim_phi    ->clear();
  vh_sim_theta  ->clear();
  vh_sim_eta    ->clear();
  vh_sim_r      ->clear();
  vh_sim_z      ->clear();
  vh_sim_tp1    ->clear();
  vh_sim_tp2    ->clear();
  //
  (*vh_size)    = 0;

  // Tracks
  vt_pt         ->clear();
  vt_xml_pt     ->clear();
  vt_phi        ->clear();
  vt_eta        ->clear();
  vt_theta      ->clear();
  vt_q          ->clear();
  //
  vt_mode       ->clear();
  vt_endcap     ->clear();
  vt_sector     ->clear();
  vt_bx         ->clear();
  //
  (*vt_size)    = 0;

  // Gen particles
  vp_pt         ->clear();
  vp_phi        ->clear();
  vp_eta        ->clear();
  vp_theta      ->clear();
  vp_vx         ->clear();
  vp_vy         ->clear();
  vp_vz         ->clear();
  vp_q          ->clear();
  vp_pdgid      ->clear();
  (*vp_size)    = 0;
}

// _____________________________________________________________________________
void NtupleMaker::beginJob() {
  makeTree();
}

void NtupleMaker::endJob() {
  writeTree();
}

// _____________________________________________________________________________
void NtupleMaker::makeTree() {

  // TFileService
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  // Hits
  vh_endcap     .reset(new std::vector<int16_t>());
  vh_station    .reset(new std::vector<int16_t>());
  vh_ring       .reset(new std::vector<int16_t>());
  vh_sector     .reset(new std::vector<int16_t>());
  vh_subsector  .reset(new std::vector<int16_t>());
  vh_chamber    .reset(new std::vector<int16_t>());
  vh_cscid      .reset(new std::vector<int16_t>());
  vh_bx         .reset(new std::vector<int16_t>());
  vh_type       .reset(new std::vector<int16_t>());
  vh_neighbor   .reset(new std::vector<int16_t>());
  //
  vh_strip      .reset(new std::vector<int16_t>());
  vh_wire       .reset(new std::vector<int16_t>());
  vh_roll       .reset(new std::vector<int16_t>());
  vh_pattern    .reset(new std::vector<int16_t>());
  vh_quality    .reset(new std::vector<int16_t>());
  vh_bend       .reset(new std::vector<int16_t>());
  vh_time       .reset(new std::vector<int16_t>());
  vh_fr         .reset(new std::vector<int16_t>());
  //
  vh_emtf_phi   .reset(new std::vector<int32_t>());
  vh_emtf_theta .reset(new std::vector<int32_t>());
  //
  vh_sim_phi    .reset(new std::vector<float  >());
  vh_sim_theta  .reset(new std::vector<float  >());
  vh_sim_eta    .reset(new std::vector<float  >());
  vh_sim_r      .reset(new std::vector<float  >());
  vh_sim_z      .reset(new std::vector<float  >());
  vh_sim_tp1    .reset(new std::vector<int32_t>());
  vh_sim_tp2    .reset(new std::vector<int32_t>());
  //
  vh_size       .reset(new int32_t(0)            );

  // Tracks
  vt_pt         .reset(new std::vector<float  >());
  vt_xml_pt     .reset(new std::vector<float  >());
  vt_phi        .reset(new std::vector<float  >());
  vt_eta        .reset(new std::vector<float  >());
  vt_theta      .reset(new std::vector<float  >());
  vt_q          .reset(new std::vector<int16_t>());
  //
  vt_mode       .reset(new std::vector<int16_t>());
  vt_endcap     .reset(new std::vector<int16_t>());
  vt_sector     .reset(new std::vector<int16_t>());
  vt_bx         .reset(new std::vector<int16_t>());
  //
  vt_size       .reset(new int32_t(0)            );

  // Gen particles
  vp_pt         .reset(new std::vector<float  >());
  vp_phi        .reset(new std::vector<float  >());
  vp_eta        .reset(new std::vector<float  >());
  vp_theta      .reset(new std::vector<float  >());
  vp_vx         .reset(new std::vector<float  >());
  vp_vy         .reset(new std::vector<float  >());
  vp_vz         .reset(new std::vector<float  >());
  vp_q          .reset(new std::vector<int16_t>());
  vp_pdgid      .reset(new std::vector<int32_t>());
  vp_size       .reset(new int32_t(0)            );

  // Set branches
  // Hits
  tree->Branch("vh_endcap"    , &(*vh_endcap    ));
  tree->Branch("vh_station"   , &(*vh_station   ));
  tree->Branch("vh_ring"      , &(*vh_ring      ));
  tree->Branch("vh_sector"    , &(*vh_sector    ));
  tree->Branch("vh_subsector" , &(*vh_subsector ));
  tree->Branch("vh_chamber"   , &(*vh_chamber   ));
  tree->Branch("vh_cscid"     , &(*vh_cscid     ));
  tree->Branch("vh_bx"        , &(*vh_bx        ));
  tree->Branch("vh_type"      , &(*vh_type      ));
  tree->Branch("vh_neighbor"  , &(*vh_neighbor  ));
  //
  tree->Branch("vh_strip"     , &(*vh_strip     ));
  tree->Branch("vh_wire"      , &(*vh_wire      ));
  tree->Branch("vh_roll"      , &(*vh_roll      ));
  tree->Branch("vh_pattern"   , &(*vh_pattern   ));
  tree->Branch("vh_quality"   , &(*vh_quality   ));
  tree->Branch("vh_bend"      , &(*vh_bend      ));
  tree->Branch("vh_time"      , &(*vh_time      ));
  tree->Branch("vh_fr"        , &(*vh_fr        ));
  //
  tree->Branch("vh_emtf_phi"  , &(*vh_emtf_phi  ));
  tree->Branch("vh_emtf_theta", &(*vh_emtf_theta));
  //
  tree->Branch("vh_sim_phi"   , &(*vh_sim_phi   ));
  tree->Branch("vh_sim_theta" , &(*vh_sim_theta ));
  tree->Branch("vh_sim_eta"   , &(*vh_sim_eta   ));
  tree->Branch("vh_sim_r"     , &(*vh_sim_r     ));
  tree->Branch("vh_sim_z"     , &(*vh_sim_z     ));
  tree->Branch("vh_sim_tp1"   , &(*vh_sim_tp1   ));
  tree->Branch("vh_sim_tp2"   , &(*vh_sim_tp2   ));
  //
  tree->Branch("vh_size"      , &(*vh_size      ));

  // Tracks
  tree->Branch("vt_pt"        , &(*vt_pt        ));
  tree->Branch("vt_xml_pt"    , &(*vt_xml_pt    ));
  tree->Branch("vt_phi"       , &(*vt_phi       ));
  tree->Branch("vt_eta"       , &(*vt_eta       ));
  tree->Branch("vt_theta"     , &(*vt_theta     ));
  tree->Branch("vt_q"         , &(*vt_q         ));
  //
  tree->Branch("vt_mode"      , &(*vt_mode      ));
  tree->Branch("vt_endcap"    , &(*vt_endcap    ));
  tree->Branch("vt_sector"    , &(*vt_sector    ));
  tree->Branch("vt_bx"        , &(*vt_bx        ));
  //
  tree->Branch("vt_size"      , &(*vt_size      ));

  // Gen particles
  tree->Branch("vp_pt"        , &(*vp_pt        ));
  tree->Branch("vp_phi"       , &(*vp_phi       ));
  tree->Branch("vp_eta"       , &(*vp_eta       ));
  tree->Branch("vp_theta"     , &(*vp_theta     ));
  tree->Branch("vp_vx"        , &(*vp_vx        ));
  tree->Branch("vp_vy"        , &(*vp_vy        ));
  tree->Branch("vp_vz"        , &(*vp_vz        ));
  tree->Branch("vp_q"         , &(*vp_q         ));
  tree->Branch("vp_pdgid"     , &(*vp_pdgid     ));
  tree->Branch("vp_size"      , &(*vp_size      ));
}

void NtupleMaker::writeTree() {
  // Handled by TFileService
}

// _____________________________________________________________________________
void NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
