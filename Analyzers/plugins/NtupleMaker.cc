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

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"


// From L1Trigger/L1TMuonEndCap/interface/MuonTriggerPrimitive.h
class TriggerPrimitive {
public:
  enum subsystem_type{kDT,kCSC,kRPC,kGEM,kME0,kNSubsystems};
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

  // Typedefs
  typedef TTTrack<Ref_Phase2TrackerDigi_>               L1TrackTriggerTrack;
  typedef std::vector<L1TrackTriggerTrack>              L1TrackTriggerTrackCollection;
  typedef TTTrackAssociationMap<Ref_Phase2TrackerDigi_> L1TrackTriggerTrackAssociator;

  template<typename T>
  edm::Handle<T> make_handle(T t)
  {
    return edm::Handle<T>();
  }

  const edm::InputTag   emuHitTag_;
  const edm::InputTag   emuTrackTag_;
  const edm::InputTag   tkTrackTag_;
  const edm::InputTag   tkTrackAssocTag_;
  const edm::InputTag   genPartTag_;
  const edm::InputTag   trkPartTag_;
  const std::string     outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<l1t::EMTFHitCollection>        emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection>      emuTrackToken_;
  edm::EDGetTokenT<L1TrackTriggerTrackCollection> tkTrackToken_;
  edm::EDGetTokenT<L1TrackTriggerTrackAssociator> tkTrackAssocToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>   genPartToken_;
  edm::EDGetTokenT<TrackingParticleCollection>    trkPartToken_;

  l1t::EMTFHitCollection        emuHits_;
  l1t::EMTFTrackCollection      emuTracks_;
  L1TrackTriggerTrackCollection tkTracks_;
  L1TrackTriggerTrackAssociator tkTrackAssoc_;
  reco::GenParticleCollection   genParts_;
  TrackingParticleCollection    trkParts_;

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
  std::unique_ptr<std::vector<int16_t> >  vh_type;  // subsystem: DT=0,CSC=1,RPC=2,GEM=3,ME0=4
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
  std::unique_ptr<std::vector<uint64_t> > vt_address;
  std::unique_ptr<std::vector<int16_t> >  vt_mode;
  std::unique_ptr<std::vector<int16_t> >  vt_endcap;
  std::unique_ptr<std::vector<int16_t> >  vt_sector;
  std::unique_ptr<std::vector<int16_t> >  vt_bx;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref1;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref2;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref3;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref4;
  //
  std::unique_ptr<int32_t              >  vt_size;

  // Track trigger tracks
  std::unique_ptr<std::vector<float  > >  vu_pt;
  std::unique_ptr<std::vector<float  > >  vu_phi;
  std::unique_ptr<std::vector<float  > >  vu_eta;
  std::unique_ptr<std::vector<float  > >  vu_theta;
  std::unique_ptr<std::vector<float  > >  vu_vx;
  std::unique_ptr<std::vector<float  > >  vu_vy;
  std::unique_ptr<std::vector<float  > >  vu_vz;
  std::unique_ptr<std::vector<int16_t> >  vu_q;  // charge
  std::unique_ptr<std::vector<float  > >  vu_rinv;
  std::unique_ptr<std::vector<float  > >  vu_chi2;
  std::unique_ptr<std::vector<int16_t> >  vu_ndof;
  std::unique_ptr<std::vector<int16_t> >  vu_sector;
  //
  std::unique_ptr<std::vector<float  > >  vu_sim_pt;
  std::unique_ptr<std::vector<float  > >  vu_sim_phi;
  std::unique_ptr<std::vector<float  > >  vu_sim_eta;
  std::unique_ptr<std::vector<int    > >  vu_sim_tp;
  std::unique_ptr<std::vector<int    > >  vu_sim_pdgid;
  std::unique_ptr<std::vector<int16_t> >  vu_sim_assoc;  // isGenuine, isLooselyGenuine, isCombinatoric, isUnknown
  //
  std::unique_ptr<int32_t>                vu_size;

  // Gen particles
  std::unique_ptr<std::vector<float  > >  vp_pt;
  std::unique_ptr<std::vector<float  > >  vp_phi;
  std::unique_ptr<std::vector<float  > >  vp_eta;
  std::unique_ptr<std::vector<float  > >  vp_theta;
  std::unique_ptr<std::vector<float  > >  vp_vx;
  std::unique_ptr<std::vector<float  > >  vp_vy;
  std::unique_ptr<std::vector<float  > >  vp_vz;
  std::unique_ptr<std::vector<int16_t> >  vp_q;  // charge
  std::unique_ptr<std::vector<int16_t> >  vp_bx;
  std::unique_ptr<std::vector<int16_t> >  vp_event;
  std::unique_ptr<std::vector<int32_t> >  vp_pdgid;
  std::unique_ptr<int32_t              >  vp_size;
};


// _____________________________________________________________________________
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) :
    truth_          (iConfig, consumesCollector()),
    emuHitTag_      (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_    (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    tkTrackTag_     (iConfig.getParameter<edm::InputTag>("tkTrackTag")),
    tkTrackAssocTag_(iConfig.getParameter<edm::InputTag>("tkTrackAssocTag")),
    genPartTag_     (iConfig.getParameter<edm::InputTag>("genPartTag")),
    trkPartTag_     (iConfig.getParameter<edm::InputTag>("trkPartTag")),
    outFileName_    (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_        (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");  // shared resources

  emuHitToken_       = consumes<l1t::EMTFHitCollection>       (emuHitTag_);
  emuTrackToken_     = consumes<l1t::EMTFTrackCollection>     (emuTrackTag_);
  tkTrackToken_      = consumes<L1TrackTriggerTrackCollection>(tkTrackTag_);
  tkTrackAssocToken_ = consumes<L1TrackTriggerTrackAssociator>(tkTrackAssocTag_);
  genPartToken_      = consumes<reco::GenParticleCollection>  (genPartTag_);
  trkPartToken_      = consumes<TrackingParticleCollection>   (trkPartTag_);
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
  auto emuHits_handle = make_handle(emuHits_);
  auto emuTracks_handle = make_handle(emuTracks_);

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

  // Track trigger tracks
  auto tkTracks_handle = make_handle(tkTracks_);
  auto tkTrackAssoc_handle = make_handle(tkTrackAssoc_);

  if (!tkTrackToken_.isUninitialized()) {
    iEvent.getByToken(tkTrackToken_, tkTracks_handle);
  }
  if (!tkTracks_handle.isValid()) {
    edm::LogError("NtupleMaker") << "Cannot get the product: " << tkTrackTag_;
  }

  if (!tkTrackAssocToken_.isUninitialized()) {
    iEvent.getByToken(tkTrackAssocToken_, tkTrackAssoc_handle);
  }
  if (!tkTrackAssoc_handle.isValid()) {
    edm::LogError("NtupleMaker") << "Cannot get the product: " << tkTrackAssocTag_;
  }

  // Gen particles
  auto genParts_handle = make_handle(genParts_);

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      edm::LogError("NtupleMaker") << "Cannot get the product: " << genPartTag_;
    }
  }

  // Tracking particles
  auto trkParts_handle = make_handle(trkParts_);

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
    if (!(-2 <= hit.BX() && hit.BX() <= 2))  continue;  // only BX=[-2,+2]
    //if (hit.Endcap() != 1)  continue;  // only positive endcap
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    //if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }

  if (tkTracks_handle.isValid()) {
    tkTracks_ = *tkTracks_handle.product();
  }

  if (tkTrackAssoc_handle.isValid()) {
    tkTrackAssoc_ = *tkTrackAssoc_handle.product();
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

      // Do not decay
      //bool nodecay = (part.decayVertices().empty());

      //if (!signal)  continue;
      //if (!intime)  continue;
      if (!outoftime) continue;
      //if (!primary) continue;
      if (!secondary) continue;
      //if (!nodecay)   continue;
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
    } else if (subsystem == TriggerPrimitive::kME0) {
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

  auto get_sim_tp_matches = [&](const auto& hit) {
    int sim_tp1 = -1;
    int sim_tp2 = -1;
    if (hit.Subsystem() == TriggerPrimitive::kCSC) {
      sim_tp1 = truth_.findCSCStripSimLink(hit, trkParts_);
      sim_tp2 = truth_.findCSCWireSimLink(hit, trkParts_);
    } else if (hit.Subsystem() == TriggerPrimitive::kRPC) {
      sim_tp1 = truth_.findRPCDigiSimLink(hit, trkParts_);
      sim_tp2 = sim_tp1;
    } else if (hit.Subsystem() == TriggerPrimitive::kGEM) {
      sim_tp1 = truth_.findGEMDigiSimLink(hit, trkParts_);
      sim_tp2 = sim_tp1;
    } else if (hit.Subsystem() == TriggerPrimitive::kME0) {
      sim_tp1 = truth_.findME0DigiSimLink(hit, trkParts_);
      sim_tp2 = sim_tp1;
    }
    return std::make_pair(sim_tp1, sim_tp2);
  };

  auto get_time = [](float time) {
    return static_cast<int>(std::round(time/0.1));  // to integer unit of 0.1 ns (arbitrary)
  };

  auto get_hit_refs = [](const auto& trk, const auto& hits) {
    using namespace l1t;

    std::vector<int32_t> hit_refs = {-1, -1, -1, -1};
    EMTFHitCollection::const_iterator conv_hits_it1  = trk.Hits().begin();
    EMTFHitCollection::const_iterator conv_hits_end1 = trk.Hits().end();

    for (; conv_hits_it1 != conv_hits_end1; ++conv_hits_it1) {
      EMTFHitCollection::const_iterator conv_hits_it2  = hits.begin();
      EMTFHitCollection::const_iterator conv_hits_end2 = hits.end();

      for (; conv_hits_it2 != conv_hits_end2; ++conv_hits_it2) {
        const EMTFHit& conv_hit_i = *conv_hits_it1;
        const EMTFHit& conv_hit_j = *conv_hits_it2;

        // See L1Trigger/L1TMuonEndCap/src/PrimitiveMatching.cc
        // All these must match: [bx_history][station][chamber][segment]
        if (
          (conv_hit_i.Subsystem()  == conv_hit_j.Subsystem()) &&
          (conv_hit_i.PC_station() == conv_hit_j.PC_station()) &&
          (conv_hit_i.PC_chamber() == conv_hit_j.PC_chamber()) &&
          (conv_hit_i.Ring()       == conv_hit_j.Ring()) &&  // because of ME1/1
          (conv_hit_i.Strip()      == conv_hit_j.Strip()) &&
          (conv_hit_i.Wire()       == conv_hit_j.Wire()) &&
          (conv_hit_i.Pattern()    == conv_hit_j.Pattern()) &&
          (conv_hit_i.BX()         == conv_hit_j.BX()) &&
          (conv_hit_i.Strip_low()  == conv_hit_j.Strip_low()) && // For RPC clusters
          (conv_hit_i.Strip_hi()   == conv_hit_j.Strip_hi()) &&  // For RPC clusters
          (conv_hit_i.Roll()       == conv_hit_j.Roll()) &&
          (conv_hit_i.Endcap()     == conv_hit_j.Endcap()) && // Needed for the ntupler
          (conv_hit_i.Sector()     == conv_hit_j.Sector()) && // Needed for the ntupler
          true
        ) {
          int istation = (conv_hit_i.Station() - 1);
          auto hit_ref = std::distance(hits.begin(), conv_hits_it2);
          hit_refs.at(istation) = hit_ref;
        }  // end if
      }  // end loop over hits
    }  // end loop over trk.Hits()

    // Sanity check
    for (int istation = 0; istation < 4; ++istation) {
      bool has_hit = trk.Mode() & (1 << (3 - istation));
      bool has_hit_check = (hit_refs.at(istation) != -1);

      if (has_hit != has_hit_check) {
        std::cout << "[ERROR] mode: " << trk.Mode() << " station: " << istation << std::endl;
        EMTFHitCollection::const_iterator conv_hits_it1  = trk.Hits().begin();
        EMTFHitCollection::const_iterator conv_hits_end1 = trk.Hits().end();
        for (; conv_hits_it1 != conv_hits_end1; ++conv_hits_it1) {
          const EMTFHit& conv_hit_i = *conv_hits_it1;
          std::cout << ".. " << conv_hit_i.Subsystem() << " " << conv_hit_i.PC_station() << " " << conv_hit_i.PC_chamber() << " " << conv_hit_i.Ring() << " " << conv_hit_i.Strip() << " " << conv_hit_i.Wire() << " " << conv_hit_i.Pattern() << " " << conv_hit_i.BX() << " " << conv_hit_i.Endcap() << " " << conv_hit_i.Sector() << std::endl;
        }
        EMTFHitCollection::const_iterator conv_hits_it2  = hits.begin();
        EMTFHitCollection::const_iterator conv_hits_end2 = hits.end();
        for (; conv_hits_it2 != conv_hits_end2; ++conv_hits_it2) {
          const EMTFHit& conv_hit_j = *conv_hits_it2;
          std::cout << ".... " << conv_hit_j.Subsystem() << " " << conv_hit_j.PC_station() << " " << conv_hit_j.PC_chamber() << " " << conv_hit_j.Ring() << " " << conv_hit_j.Strip() << " " << conv_hit_j.Wire() << " " << conv_hit_j.Pattern() << " " << conv_hit_j.BX() << " " << conv_hit_j.Endcap() << " " << conv_hit_j.Sector() << std::endl;
        }
      }

      assert(has_hit == has_hit_check);
    }

    return hit_refs;
  };


  // ___________________________________________________________________________
  bool please_use_trkParts = true;

  if (verbose_ > 0) {
    std::cout << "[DEBUG] # hits: " << emuHits_.size() << " #  tracks: " << emuTracks_.size() << " # gen parts: " << genParts_.size() << " # trk parts: " << trkParts_.size() << std::endl;
  }

  // Hits
  prepare_sim_tp();  // must be called before calling get_sim_tp1() and get_sim_tp2()

  for (const auto& hit : emuHits_) {
    const std::pair<int,int>& sim_tp_pair = get_sim_tp_matches(hit);

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
    vh_sim_tp1    ->push_back(sim_tp_pair.first);
    vh_sim_tp2    ->push_back(sim_tp_pair.second);
  }
  (*vh_size) = emuHits_.size();

  // Tracks
  for (const auto& trk : emuTracks_) {
    const auto& hit_refs = get_hit_refs(trk, emuHits_);
    assert(hit_refs.size() == 4);

    const l1t::EMTFPtLUT& ptlut_data = trk.PtLUT();

    vt_pt         ->push_back(trk.Pt());
    vt_xml_pt     ->push_back(trk.Pt_XML());
    vt_phi        ->push_back(trk.Phi_glob());
    vt_eta        ->push_back(trk.Eta());
    vt_theta      ->push_back(trk.Theta());
    vt_q          ->push_back(trk.Charge());
    //
    vt_address    ->push_back(ptlut_data.address);
    vt_mode       ->push_back(trk.Mode());
    vt_endcap     ->push_back(trk.Endcap());
    vt_sector     ->push_back(trk.Sector());
    vt_bx         ->push_back(trk.BX());
    vt_hitref1    ->push_back(hit_refs.at(0));
    vt_hitref2    ->push_back(hit_refs.at(1));
    vt_hitref3    ->push_back(hit_refs.at(2));
    vt_hitref4    ->push_back(hit_refs.at(3));
  }
  (*vt_size) = emuTracks_.size();

  // L1TrackTrigger tracks
  int itkTrack = 0;
  auto tkTracks_handle = make_handle(tkTracks_);

  for (const auto& trk : tkTracks_) {
    const GlobalVector& momentum = trk.getMomentum();
    const GlobalPoint&  poca     = trk.getPOCA();
    double              rinv     = trk.getRInv();

    vu_pt    ->push_back(momentum.perp());
    vu_phi   ->push_back(momentum.phi());
    vu_eta   ->push_back(momentum.eta());
    vu_theta ->push_back(momentum.theta());
    vu_vx    ->push_back(poca.x());
    vu_vy    ->push_back(poca.y());
    vu_vz    ->push_back(poca.z());
    vu_q     ->push_back(rinv >= 0 ? 1 : -1);
    vu_rinv  ->push_back(rinv);
    vu_chi2  ->push_back(trk.getChi2());
    vu_ndof  ->push_back(trk.getStubRefs().size()*2 - 4);  // nPar=4
    vu_sector->push_back(trk.getSector());

    float sim_pt    = 0.;
    float sim_phi   = 0.;
    float sim_eta   = 0.;
    int   sim_tp    = -1;
    int   sim_pdgid = 0;
    int   sim_assoc = 0;
    {
      edm::Ptr<L1TrackTriggerTrack> trkPtr(tkTracks_handle, itkTrack);
      edm::Ptr<TrackingParticle> tpPtr = tkTrackAssoc_.findTrackingParticlePtr(trkPtr);
      if (tpPtr.isNonnull()) {
        sim_pt    = tpPtr->pt();
        sim_eta   = tpPtr->eta();
        sim_phi   = tpPtr->phi();
        sim_tp    = tpPtr.key();
        sim_pdgid = tpPtr->pdgId();
      }
      if (tkTrackAssoc_.isGenuine(trkPtr)) {
        sim_assoc = 0;
      } else if (tkTrackAssoc_.isLooselyGenuine(trkPtr)) {
        sim_assoc = 1;
      } else if (tkTrackAssoc_.isCombinatoric(trkPtr)) {
        sim_assoc = 2;
      } else if (tkTrackAssoc_.isUnknown(trkPtr)) {
        sim_assoc = 3;
      } else {
        sim_assoc = 4;
      }
    }

    vu_sim_pt   ->push_back(sim_pt);
    vu_sim_phi  ->push_back(sim_phi);
    vu_sim_eta  ->push_back(sim_eta);
    vu_sim_tp   ->push_back(sim_tp);
    vu_sim_pdgid->push_back(sim_pdgid);
    vu_sim_assoc->push_back(sim_assoc);

    ++itkTrack;
  }
  (*vu_size) = tkTracks_.size();

  // ___________________________________________________________________________
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
      vp_bx         ->push_back(0);
      vp_event      ->push_back(0);
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
      vp_bx         ->push_back(part.eventId().bunchCrossing());
      vp_event      ->push_back(part.eventId().event());
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
  vt_address    ->clear();
  vt_mode       ->clear();
  vt_endcap     ->clear();
  vt_sector     ->clear();
  vt_bx         ->clear();
  vt_hitref1    ->clear();
  vt_hitref2    ->clear();
  vt_hitref3    ->clear();
  vt_hitref4    ->clear();
  //
  (*vt_size)    = 0;

  // L1TrackTrigger tracks
  vu_pt         ->clear();
  vu_phi        ->clear();
  vu_eta        ->clear();
  vu_theta      ->clear();
  vu_q          ->clear();
  vu_vx         ->clear();
  vu_vy         ->clear();
  vu_vz         ->clear();
  vu_rinv       ->clear();
  vu_chi2       ->clear();
  vu_ndof       ->clear();
  vu_sector     ->clear();
  //
  vu_sim_pt     ->clear();
  vu_sim_phi    ->clear();
  vu_sim_eta    ->clear();
  vu_sim_tp     ->clear();
  vu_sim_pdgid  ->clear();
  vu_sim_assoc  ->clear();
  //
  (*vu_size)    = 0;

  // Gen particles
  vp_pt         ->clear();
  vp_phi        ->clear();
  vp_eta        ->clear();
  vp_theta      ->clear();
  vp_vx         ->clear();
  vp_vy         ->clear();
  vp_vz         ->clear();
  vp_q          ->clear();
  vp_bx         ->clear();
  vp_event      ->clear();
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
  vh_endcap     .reset(new std::vector<int16_t >());
  vh_station    .reset(new std::vector<int16_t >());
  vh_ring       .reset(new std::vector<int16_t >());
  vh_sector     .reset(new std::vector<int16_t >());
  vh_subsector  .reset(new std::vector<int16_t >());
  vh_chamber    .reset(new std::vector<int16_t >());
  vh_cscid      .reset(new std::vector<int16_t >());
  vh_bx         .reset(new std::vector<int16_t >());
  vh_type       .reset(new std::vector<int16_t >());
  vh_neighbor   .reset(new std::vector<int16_t >());
  //
  vh_strip      .reset(new std::vector<int16_t >());
  vh_wire       .reset(new std::vector<int16_t >());
  vh_roll       .reset(new std::vector<int16_t >());
  vh_pattern    .reset(new std::vector<int16_t >());
  vh_quality    .reset(new std::vector<int16_t >());
  vh_bend       .reset(new std::vector<int16_t >());
  vh_time       .reset(new std::vector<int16_t >());
  vh_fr         .reset(new std::vector<int16_t >());
  //
  vh_emtf_phi   .reset(new std::vector<int32_t >());
  vh_emtf_theta .reset(new std::vector<int32_t >());
  //
  vh_sim_phi    .reset(new std::vector<float   >());
  vh_sim_theta  .reset(new std::vector<float   >());
  vh_sim_eta    .reset(new std::vector<float   >());
  vh_sim_r      .reset(new std::vector<float   >());
  vh_sim_z      .reset(new std::vector<float   >());
  vh_sim_tp1    .reset(new std::vector<int32_t >());
  vh_sim_tp2    .reset(new std::vector<int32_t >());
  //
  vh_size       .reset(new int32_t(0)             );

  // Tracks
  vt_pt         .reset(new std::vector<float   >());
  vt_xml_pt     .reset(new std::vector<float   >());
  vt_phi        .reset(new std::vector<float   >());
  vt_eta        .reset(new std::vector<float   >());
  vt_theta      .reset(new std::vector<float   >());
  vt_q          .reset(new std::vector<int16_t >());
  //
  vt_address    .reset(new std::vector<uint64_t>());
  vt_mode       .reset(new std::vector<int16_t >());
  vt_endcap     .reset(new std::vector<int16_t >());
  vt_sector     .reset(new std::vector<int16_t >());
  vt_bx         .reset(new std::vector<int16_t >());
  vt_hitref1    .reset(new std::vector<int32_t >());
  vt_hitref2    .reset(new std::vector<int32_t >());
  vt_hitref3    .reset(new std::vector<int32_t >());
  vt_hitref4    .reset(new std::vector<int32_t >());
  //
  vt_size       .reset(new int32_t(0)             );

  // L1TrackTrigger trakcs
  vu_pt         .reset(new std::vector<float   >());
  vu_phi        .reset(new std::vector<float   >());
  vu_eta        .reset(new std::vector<float   >());
  vu_theta      .reset(new std::vector<float   >());
  vu_vx         .reset(new std::vector<float   >());
  vu_vy         .reset(new std::vector<float   >());
  vu_vz         .reset(new std::vector<float   >());
  vu_q          .reset(new std::vector<int16_t >());
  vu_rinv       .reset(new std::vector<float   >());
  vu_chi2       .reset(new std::vector<float   >());
  vu_ndof       .reset(new std::vector<int16_t >());
  vu_sector     .reset(new std::vector<int16_t >());
  //
  vu_sim_pt     .reset(new std::vector<float   >());
  vu_sim_phi    .reset(new std::vector<float   >());
  vu_sim_eta    .reset(new std::vector<float   >());
  vu_sim_tp     .reset(new std::vector<int     >());
  vu_sim_pdgid  .reset(new std::vector<int     >());
  vu_sim_assoc  .reset(new std::vector<int16_t >());
  //
  vu_size       .reset(new int32_t(0)             );

  // Gen particles
  vp_pt         .reset(new std::vector<float   >());
  vp_phi        .reset(new std::vector<float   >());
  vp_eta        .reset(new std::vector<float   >());
  vp_theta      .reset(new std::vector<float   >());
  vp_vx         .reset(new std::vector<float   >());
  vp_vy         .reset(new std::vector<float   >());
  vp_vz         .reset(new std::vector<float   >());
  vp_q          .reset(new std::vector<int16_t >());
  vp_bx         .reset(new std::vector<int16_t >());
  vp_event      .reset(new std::vector<int16_t >());
  vp_pdgid      .reset(new std::vector<int32_t >());
  vp_size       .reset(new int32_t(0)             );

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
  tree->Branch("vt_address"   , &(*vt_address   ));
  tree->Branch("vt_mode"      , &(*vt_mode      ));
  tree->Branch("vt_endcap"    , &(*vt_endcap    ));
  tree->Branch("vt_sector"    , &(*vt_sector    ));
  tree->Branch("vt_bx"        , &(*vt_bx        ));
  tree->Branch("vt_hitref1"   , &(*vt_hitref1   ));
  tree->Branch("vt_hitref2"   , &(*vt_hitref2   ));
  tree->Branch("vt_hitref3"   , &(*vt_hitref3   ));
  tree->Branch("vt_hitref4"   , &(*vt_hitref4   ));
  //
  tree->Branch("vt_size"      , &(*vt_size      ));

  // L1TrackTrigger tracks
  tree->Branch("vu_pt"        , &(*vu_pt        ));
  tree->Branch("vu_phi"       , &(*vu_phi       ));
  tree->Branch("vu_eta"       , &(*vu_eta       ));
  tree->Branch("vu_theta"     , &(*vu_theta     ));
  tree->Branch("vu_q"         , &(*vu_q         ));
  tree->Branch("vu_vx"        , &(*vu_vx        ));
  tree->Branch("vu_vy"        , &(*vu_vy        ));
  tree->Branch("vu_vz"        , &(*vu_vz        ));
  tree->Branch("vu_rinv"      , &(*vu_rinv      ));
  tree->Branch("vu_chi2"      , &(*vu_chi2      ));
  tree->Branch("vu_ndof"      , &(*vu_ndof      ));
  tree->Branch("vu_sector"    , &(*vu_sector    ));
  //
  tree->Branch("vu_sim_pt"    , &(*vu_sim_pt    ));
  tree->Branch("vu_sim_phi"   , &(*vu_sim_phi   ));
  tree->Branch("vu_sim_eta"   , &(*vu_sim_eta   ));
  tree->Branch("vu_sim_tp"    , &(*vu_sim_tp    ));
  tree->Branch("vu_sim_pdgid" , &(*vu_sim_pdgid ));
  tree->Branch("vu_sim_assoc" , &(*vu_sim_assoc ));
  //
  tree->Branch("vu_size"      , &(*vu_size      ));

  // Gen particles
  tree->Branch("vp_pt"        , &(*vp_pt        ));
  tree->Branch("vp_phi"       , &(*vp_phi       ));
  tree->Branch("vp_eta"       , &(*vp_eta       ));
  tree->Branch("vp_theta"     , &(*vp_theta     ));
  tree->Branch("vp_vx"        , &(*vp_vx        ));
  tree->Branch("vp_vy"        , &(*vp_vy        ));
  tree->Branch("vp_vz"        , &(*vp_vz        ));
  tree->Branch("vp_q"         , &(*vp_q         ));
  tree->Branch("vp_bx"        , &(*vp_bx        ));
  tree->Branch("vp_event"     , &(*vp_event     ));
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
