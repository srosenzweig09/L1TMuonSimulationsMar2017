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
#include "FWCore/Framework/interface/EDAnalyzer.h"
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
class NtupleMaker : public edm::EDAnalyzer {
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
  void getHandles(const edm::Event& iEvent);

  void makeTree();
  void writeTree();

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
  std::unique_ptr<std::vector<int16_t> >  vh_type;  // subsystem: DT=1,CSC=2,RPC=3,GEM=4
  //
  std::unique_ptr<std::vector<int16_t> >  vh_strip;
  std::unique_ptr<std::vector<int16_t> >  vh_wire;
  std::unique_ptr<std::vector<int16_t> >  vh_roll;
  std::unique_ptr<std::vector<int16_t> >  vh_pattern;
  std::unique_ptr<std::vector<int16_t> >  vh_quality;
  std::unique_ptr<std::vector<int16_t> >  vh_bend;
  //
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_phi;
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_theta;
  //
  std::unique_ptr<std::vector<float  > >  vh_sim_phi;
  std::unique_ptr<std::vector<float  > >  vh_sim_theta;
  std::unique_ptr<std::vector<float  > >  vh_sim_eta;

  // Tracks
  std::unique_ptr<std::vector<float  > >  vp_pt;
  std::unique_ptr<std::vector<float  > >  vp_eta;
  std::unique_ptr<std::vector<float  > >  vp_phi;
  std::unique_ptr<std::vector<float  > >  vp_p;  // |p|
  std::unique_ptr<std::vector<float  > >  vp_q;  // charge
};


// _____________________________________________________________________________
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) :
    emuHitTag_    (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_  (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    genPartTag_   (iConfig.getParameter<edm::InputTag>("genPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  emuHitToken_   = consumes<l1t::EMTFHitCollection>     (emuHitTag_);
  emuTrackToken_ = consumes<l1t::EMTFTrackCollection>   (emuTrackTag_);
  genPartToken_  = consumes<reco::GenParticleCollection>(genPartTag_);
}

NtupleMaker::~NtupleMaker() {}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent);
  process();
}

// _____________________________________________________________________________
void NtupleMaker::getHandles(const edm::Event& iEvent) {

  // EMTF hits and tracks
  edm::Handle<decltype(emuHits_)>   emuHits_handle;
  edm::Handle<decltype(emuTracks_)> emuTracks_handle;

  if (!emuHitToken_.isUninitialized()) {
    iEvent.getByToken(emuHitToken_, emuHits_handle);
  }
  if (!emuHits_handle.isValid()) {
    edm::LogError("RPCIntegration") << "Cannot get the product: " << emuHitTag_;
    return;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("RPCIntegration") << "Cannot get the product: " << emuTrackTag_;
    return;
  }

  // Gen particles
  edm::Handle<decltype(genParts_)> genParts_handle;

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      edm::LogError("RPCIntegration") << "Cannot get the product: " << genPartTag_;
      return;
    }
  }

  // Object filters
  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    if (!(-1 <= hit.BX() && hit.BX() <= 1))  continue;  // only BX=[-1,+1]
    if (hit.Endcap() != 1)  continue;  // only positive endcap
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }

  genParts_.clear();
  for (const auto& part : (*genParts_handle)) {
    //if (!(part.pt() >= 2.))     continue;  // only pT > 2
    if (!(1.24 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap
    genParts_.push_back(part);
  }
}

// _____________________________________________________________________________
void NtupleMaker::process() {

  // Hits
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
    //
    vh_strip      ->push_back(hit.Strip());
    vh_wire       ->push_back(hit.Wire());
    vh_roll       ->push_back(hit.Roll());
    vh_pattern    ->push_back(hit.Pattern());
    vh_quality    ->push_back(hit.Quality());
    vh_bend       ->push_back(hit.Bend());
    //
    vh_emtf_phi   ->push_back(hit.Phi_fp());
    vh_emtf_theta ->push_back(hit.Theta_fp());
    //
    vh_sim_phi    ->push_back(hit.Phi_sim());
    vh_sim_theta  ->push_back(hit.Theta_sim());
    vh_sim_eta    ->push_back(hit.Eta_sim());
  }

  // Tracks
  for (const auto& part : genParts_) {
    vp_pt         ->push_back(part.pt());
    vp_eta        ->push_back(part.eta());
    vp_phi        ->push_back(part.phi());
    vp_p          ->push_back(part.p());
    vp_q          ->push_back(part.charge());
  }

  assert(vp_pt->size() <= 1);  // expect 0 or 1 gen particle


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
  //
  vh_strip      ->clear();
  vh_wire       ->clear();
  vh_roll       ->clear();
  vh_pattern    ->clear();
  vh_quality    ->clear();
  vh_bend       ->clear();
  //
  vh_emtf_phi   ->clear();
  vh_emtf_theta ->clear();
  //
  vh_sim_phi    ->clear();
  vh_sim_theta  ->clear();
  vh_sim_eta    ->clear();

  // Tracks
  vp_pt         ->clear();
  vp_eta        ->clear();
  vp_phi        ->clear();
  vp_p          ->clear();
  vp_q          ->clear();
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

  // Create TTree
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
  //
  vh_strip      .reset(new std::vector<int16_t>());
  vh_wire       .reset(new std::vector<int16_t>());
  vh_roll       .reset(new std::vector<int16_t>());
  vh_pattern    .reset(new std::vector<int16_t>());
  vh_quality    .reset(new std::vector<int16_t>());
  vh_bend       .reset(new std::vector<int16_t>());
  //
  vh_emtf_phi   .reset(new std::vector<int32_t>());
  vh_emtf_theta .reset(new std::vector<int32_t>());
  //
  vh_sim_phi    .reset(new std::vector<float  >());
  vh_sim_theta  .reset(new std::vector<float  >());
  vh_sim_eta    .reset(new std::vector<float  >());

  // Tracks
  vp_pt         .reset(new std::vector<float  >());
  vp_eta        .reset(new std::vector<float  >());
  vp_phi        .reset(new std::vector<float  >());
  vp_p          .reset(new std::vector<float  >());
  vp_q          .reset(new std::vector<float  >());

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
  //
  tree->Branch("vh_strip"     , &(*vh_strip     ));
  tree->Branch("vh_wire"      , &(*vh_wire      ));
  tree->Branch("vh_roll"      , &(*vh_roll      ));
  tree->Branch("vh_pattern"   , &(*vh_pattern   ));
  tree->Branch("vh_quality"   , &(*vh_quality   ));
  tree->Branch("vh_bend"      , &(*vh_bend      ));
  //
  tree->Branch("vh_emtf_phi"  , &(*vh_emtf_phi  ));
  tree->Branch("vh_emtf_theta", &(*vh_emtf_theta));
  //
  tree->Branch("vh_sim_phi"   , &(*vh_sim_phi   ));
  tree->Branch("vh_sim_theta" , &(*vh_sim_theta ));
  tree->Branch("vh_sim_eta"   , &(*vh_sim_eta   ));

  // Tracks
  tree->Branch("vp_pt"        , &(*vp_pt        ));
  tree->Branch("vp_eta"       , &(*vp_eta       ));
  tree->Branch("vp_phi"       , &(*vp_phi       ));
  tree->Branch("vp_p"         , &(*vp_p         ));
  tree->Branch("vp_q"         , &(*vp_q         ));

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
