#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdint>

#include "TString.h"
#include "TFile.h"
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
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "L1Trigger/L1TMuonEndCap/interface/TrackTools.h"

#include "L1TMuonSimulations/Analyzers/interface/EMTFParticleTools.h"
#include "L1TMuonSimulations/Analyzers/interface/EMTFMCTruth.h"


// From L1Trigger/L1TMuonEndCap/interface/MuonTriggerPrimitive.h
class TriggerPrimitive {
public:
  enum subsystem_type{kDT,kCSC,kRPC,kGEM,kME0,kNSubsystems};
};


// _____________________________________________________________________________
class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
public:
  explicit NtupleMaker(const edm::ParameterSet& iConfig);
  ~NtupleMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void endRun(const edm::Run&, const edm::EventSetup&) override;

  // Main functions
  void process(const edm::Event& iEvent, const edm::EventSetup& iSetup);

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
  edm::Handle<T> make_handle(T& t)
  {
    return edm::Handle<T>();
  }

  template<typename T>
  edm::Handle<T> make_handle(T* t)
  {
    return edm::Handle<T>();
  }

  // Input parameters
  const edm::InputTag   emuHitTag_;
  const edm::InputTag   emuTrackTag_;
  const edm::InputTag   tkTrackTag_;
  const edm::InputTag   tkTrackAssocTag_;
  const edm::InputTag   genPartTag_;
  const edm::InputTag   simTrackTag_;
  const edm::InputTag   trkPartTag_;
  const edm::InputTag   pileupInfoTag_;
  const std::string     outFileName_;
  int verbose_;

  // Member data
  bool firstEvent_;
  edm::EDGetTokenT<l1t::EMTFHitCollection>          emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection>        emuTrackToken_;
  edm::EDGetTokenT<L1TrackTriggerTrackCollection>   tkTrackToken_;
  edm::EDGetTokenT<L1TrackTriggerTrackAssociator>   tkTrackAssocToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>     genPartToken_;
  edm::EDGetTokenT<edm::SimTrackContainer>          simTrackToken_;
  edm::EDGetTokenT<TrackingParticleCollection>      trkPartToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;

  l1t::EMTFHitCollection        emuHits_;
  l1t::EMTFTrackCollection      emuTracks_;
  TrackingParticleCollection    trkParts_;

  // For edm products
  const L1TrackTriggerTrackCollection*  tkTracks_;
  const L1TrackTriggerTrackAssociator*  tkTrackAssoc_;
  const reco::GenParticleCollection*    genParts_;
  const edm::SimTrackContainer*         simTracks_;
  const std::vector<PileupSummaryInfo>* pileupInfo_;

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
  std::unique_ptr<std::vector<int16_t> >  vh_quality;
  std::unique_ptr<std::vector<int16_t> >  vh_pattern;
  std::unique_ptr<std::vector<int16_t> >  vh_bend;
  std::unique_ptr<std::vector<int16_t> >  vh_time;
  std::unique_ptr<std::vector<int16_t> >  vh_fr;
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_phi;   // integer unit
  std::unique_ptr<std::vector<int32_t> >  vh_emtf_theta; // integer unit
  //
  std::unique_ptr<std::vector<float  > >  vh_sim_phi;    // in degrees
  std::unique_ptr<std::vector<float  > >  vh_sim_theta;  // in degrees
  std::unique_ptr<std::vector<float  > >  vh_sim_r;      // in cm
  std::unique_ptr<std::vector<float  > >  vh_sim_z;      // in cm
  std::unique_ptr<std::vector<int32_t> >  vh_sim_tp1;
  std::unique_ptr<std::vector<int32_t> >  vh_sim_tp2;
  std::unique_ptr<int32_t              >  vh_size;

  // EndcapSimHits
  std::unique_ptr<std::vector<int16_t> >  vc_type;
  std::unique_ptr<std::vector<int16_t> >  vc_station;
  std::unique_ptr<std::vector<int16_t> >  vc_ring;
  std::unique_ptr<std::vector<int16_t> >  vc_layer;
  std::unique_ptr<std::vector<int16_t> >  vc_chamber;
  std::unique_ptr<std::vector<float  > >  vc_phi;        // in degrees
  std::unique_ptr<std::vector<float  > >  vc_theta;      // in degrees
  std::unique_ptr<std::vector<float  > >  vc_r;          // in cm
  std::unique_ptr<std::vector<float  > >  vc_z;          // in cm
  std::unique_ptr<std::vector<int32_t> >  vc_sim_tp;
  std::unique_ptr<std::vector<int32_t> >  vc_pdgid;      // particleType
  std::unique_ptr<std::vector<int16_t> >  vc_process;    // processType
  std::unique_ptr<std::vector<float  > >  vc_mom_phi;    // momentum at entry
  std::unique_ptr<std::vector<float  > >  vc_mom_theta;  // momentum at entry
  std::unique_ptr<std::vector<float  > >  vc_tof;        // time of flight
  std::unique_ptr<int32_t              >  vc_size;

  // Tracks
  std::unique_ptr<std::vector<float  > >  vt_pt;
  std::unique_ptr<std::vector<float  > >  vt_xml_pt;
  std::unique_ptr<std::vector<float  > >  vt_pt_dxy;
  std::unique_ptr<std::vector<float  > >  vt_dxy;
  std::unique_ptr<std::vector<float  > >  vt_invpt_prompt;
  std::unique_ptr<std::vector<float  > >  vt_invpt_displ;
  std::unique_ptr<std::vector<float  > >  vt_phi;        // in degrees
  std::unique_ptr<std::vector<float  > >  vt_theta;      // in degrees
  std::unique_ptr<std::vector<float  > >  vt_eta;
  std::unique_ptr<std::vector<int16_t> >  vt_q;          // charge
  //
  std::unique_ptr<std::vector<uint64_t> > vt_address;
  std::unique_ptr<std::vector<int16_t> >  vt_mode;
  std::unique_ptr<std::vector<int16_t> >  vt_endcap;
  std::unique_ptr<std::vector<int16_t> >  vt_sector;
  std::unique_ptr<std::vector<int16_t> >  vt_bx;
  std::unique_ptr<std::vector<int16_t> >  vt_nhits;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref1;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref2;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref3;
  std::unique_ptr<std::vector<int32_t> >  vt_hitref4;
  std::unique_ptr<int32_t              >  vt_size;

  // Track trigger tracks
  std::unique_ptr<std::vector<float  > >  vu_pt;
  std::unique_ptr<std::vector<float  > >  vu_phi;        // in radians
  std::unique_ptr<std::vector<float  > >  vu_theta;      // in radians
  std::unique_ptr<std::vector<float  > >  vu_eta;
  std::unique_ptr<std::vector<float  > >  vu_vx;         // in cm
  std::unique_ptr<std::vector<float  > >  vu_vy;         // in cm
  std::unique_ptr<std::vector<float  > >  vu_vz;         // in cm
  std::unique_ptr<std::vector<int16_t> >  vu_q;          // charge
  std::unique_ptr<std::vector<float  > >  vu_rinv;
  std::unique_ptr<std::vector<float  > >  vu_chi2;
  std::unique_ptr<std::vector<int16_t> >  vu_ndof;
  std::unique_ptr<std::vector<int16_t> >  vu_sector;
  //
  std::unique_ptr<std::vector<float  > >  vu_sim_pt;
  std::unique_ptr<std::vector<float  > >  vu_sim_phi;
  std::unique_ptr<std::vector<float  > >  vu_sim_eta;
  std::unique_ptr<std::vector<int32_t> >  vu_sim_tp;
  std::unique_ptr<std::vector<int32_t> >  vu_sim_pdgid;
  std::unique_ptr<std::vector<int16_t> >  vu_sim_assoc;  // isGenuine, isLooselyGenuine, isCombinatoric, isUnknown
  std::unique_ptr<int32_t              >  vu_size;

  // Tracking particles
  std::unique_ptr<std::vector<float  > >  vp_pt;
  std::unique_ptr<std::vector<float  > >  vp_phi;        // in radians
  std::unique_ptr<std::vector<float  > >  vp_theta;      // in radians
  std::unique_ptr<std::vector<float  > >  vp_eta;
  std::unique_ptr<std::vector<float  > >  vp_vx;         // in cm
  std::unique_ptr<std::vector<float  > >  vp_vy;         // in cm
  std::unique_ptr<std::vector<float  > >  vp_vz;         // in cm
  std::unique_ptr<std::vector<float  > >  vp_invpt;
  std::unique_ptr<std::vector<float  > >  vp_d0;         // in cm
  std::unique_ptr<std::vector<float  > >  vp_beta;
  std::unique_ptr<std::vector<float  > >  vp_mass;
  std::unique_ptr<std::vector<int16_t> >  vp_q;          // charge
  std::unique_ptr<std::vector<int16_t> >  vp_bx;
  std::unique_ptr<std::vector<int16_t> >  vp_event;
  std::unique_ptr<std::vector<int32_t> >  vp_pdgid;
  std::unique_ptr<std::vector<int16_t> >  vp_status;
  std::unique_ptr<std::vector<int16_t> >  vp_decay;
  std::unique_ptr<std::vector<int32_t> >  vp_genp;
  std::unique_ptr<int32_t              >  vp_size;

  // Event info
  std::unique_ptr<std::vector<uint64_t> > ve_event;
  std::unique_ptr<std::vector<uint32_t> > ve_run;
  std::unique_ptr<std::vector<uint32_t> > ve_lumi;
  std::unique_ptr<std::vector<float  > >  ve_npv;  // getTrueNumInteractions()
  std::unique_ptr<std::vector<int32_t> >  ve_nvtx; // getPU_NumInteractions()
  std::unique_ptr<int32_t              >  ve_size;
};


// _____________________________________________________________________________
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) :
    truth_          (iConfig, consumesCollector()),
    emuHitTag_      (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_    (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    tkTrackTag_     (iConfig.getParameter<edm::InputTag>("tkTrackTag")),
    tkTrackAssocTag_(iConfig.getParameter<edm::InputTag>("tkTrackAssocTag")),
    genPartTag_     (iConfig.getParameter<edm::InputTag>("genPartTag")),
    simTrackTag_    (iConfig.getParameter<edm::InputTag>("simTrackTag")),
    trkPartTag_     (iConfig.getParameter<edm::InputTag>("trkPartTag")),
    pileupInfoTag_  (iConfig.getParameter<edm::InputTag>("pileupInfoTag")),
    outFileName_    (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_        (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");  // shared resources

  firstEvent_ = true;

  emuHitToken_       = consumes<l1t::EMTFHitCollection>         (emuHitTag_);
  emuTrackToken_     = consumes<l1t::EMTFTrackCollection>       (emuTrackTag_);
  tkTrackToken_      = consumes<L1TrackTriggerTrackCollection>  (tkTrackTag_);
  tkTrackAssocToken_ = consumes<L1TrackTriggerTrackAssociator>  (tkTrackAssocTag_);
  genPartToken_      = consumes<reco::GenParticleCollection>    (genPartTag_);
  simTrackToken_     = consumes<edm::SimTrackContainer>         (simTrackTag_);
  trkPartToken_      = consumes<TrackingParticleCollection>     (trkPartTag_);
  pileupInfoToken_   = consumes<std::vector<PileupSummaryInfo> >(pileupInfoTag_);
}

NtupleMaker::~NtupleMaker() {}

void NtupleMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void NtupleMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent, iSetup);
  process(iEvent, iSetup);
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
    if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << emuHitTag_;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << emuTrackTag_;
  }

  // Track trigger tracks
  auto tkTracks_handle = make_handle(tkTracks_);
  auto tkTrackAssoc_handle = make_handle(tkTrackAssoc_);

  if (!tkTrackToken_.isUninitialized()) {
    iEvent.getByToken(tkTrackToken_, tkTracks_handle);
  }
  if (!tkTracks_handle.isValid()) {
    if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << tkTrackTag_;
    tkTracks_ = nullptr;
  } else {
    tkTracks_ = tkTracks_handle.product();
  }

  if (!tkTrackAssocToken_.isUninitialized()) {
    iEvent.getByToken(tkTrackAssocToken_, tkTrackAssoc_handle);
  }
  if (!tkTrackAssoc_handle.isValid()) {
    if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << tkTrackAssocTag_;
    tkTrackAssoc_ = nullptr;
  } else {
    tkTrackAssoc_ = tkTrackAssoc_handle.product();
  }

  // Gen particles
  auto genParts_handle = make_handle(genParts_);

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << genPartTag_;
      genParts_ = nullptr;
    } else {
      genParts_ = genParts_handle.product();
    }
  }

  // Sim tracks
  auto simTracks_handle = make_handle(simTracks_);

  if (!iEvent.isRealData()) {
    if (!simTrackToken_.isUninitialized()) {
      iEvent.getByToken(simTrackToken_, simTracks_handle);
    }
    if (!simTracks_handle.isValid()) {
      if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << simTrackTag_;
      simTracks_ = nullptr;
    } else {
      simTracks_ = simTracks_handle.product();
    }
  }

  // Tracking particles
  auto trkParts_handle = make_handle(trkParts_);

  if (!iEvent.isRealData()) {
    if (!trkPartToken_.isUninitialized()) {
      iEvent.getByToken(trkPartToken_, trkParts_handle);
    }
    if (!trkParts_handle.isValid()) {
      if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << trkPartTag_;
    }
  }

  // Pileup info
  auto pileupInfo_handle = make_handle(pileupInfo_);

  if (!iEvent.isRealData()) {
    if (!pileupInfoToken_.isUninitialized()) {
      iEvent.getByToken(pileupInfoToken_, pileupInfo_handle);
    }
    if (!pileupInfo_handle.isValid()) {
      if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << pileupInfoTag_;
      pileupInfo_ = nullptr;
    } else {
      pileupInfo_ = pileupInfo_handle.product();
    }
  }


  // ___________________________________________________________________________
  // Object filters
  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    if (!(-2 <= hit.BX() && hit.BX() <= 2))  continue;  // only BX=[-2,+2]
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    emuTracks_.push_back(trk);
  }

  trkParts_.clear();
  for (const auto& part : (*trkParts_handle)) {

    int pdgId = std::abs(part.pdgId());
    //if (!(part.pt() >= 2.))     continue;  // only pT > 2
    if (!(pdgId == 13 || pdgId == 1000015))  continue;  // only muons (or stau)

    // Tracking particle selection
    {
      // Signal event
      //bool signal = (part.eventId().event() == 0);

      // In time bunch-crossing
      //bool intime = (part.eventId().bunchCrossing() == 0);

      // In time + out of time bunch-crossing (-2 <= BX <= +2)
      bool outoftime = (-2 <= part.eventId().bunchCrossing() && part.eventId().bunchCrossing() <= +2);

      // Primary+charged: pT > 0.2 GeV, |eta| < 3.0, |rho0| < 0.5 cm, |z0| < 30 cm
      //bool primary = (part.charge() != 0 && part.pt() > 0.2 && std::abs(part.eta()) < 3.0 && std::sqrt(part.vx() * part.vx() + part.vy() * part.vy()) < 0.5 && std::abs(part.vz()) < 30.0);

      // Primary+secondary pT > 0.5 GeV, |eta| < 3.0, |x0| < 300 cm, |y0| < 300 cm, |z0| < 500 cm
      bool secondary = (part.charge() != 0 && part.pt() > 0.5 && std::abs(part.eta()) <= 3.0 && std::abs(part.vx()) <= 300. && std::abs(part.vy()) <= 300. && std::abs(part.vz()) <= 500.);

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
void NtupleMaker::process(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
      //// 10 degree rings have even subsectors in front
      //// 20 degree rings have odd subsectors in front
      //bool is_10degree = !((station == 3 || station == 4) && (ring == 1));
      //bool isEven = (subsector % 2 == 0);
      //result = (is_10degree) ? isEven : !isEven;

      // Use the equivalent CSC chamber F/R
      bool isEven = (chamber % 2 == 0);
      result = (station < 3) ? isEven : !isEven;
    } else if (subsystem == TriggerPrimitive::kGEM) {
      //
      result = (chamber % 2 == 0);
    } else if (subsystem == TriggerPrimitive::kME0) {
      //
      result = (chamber % 2 == 0);
    } else if (subsystem == TriggerPrimitive::kDT) {
      //
      result = (chamber % 2 == 0);
    }
    return result;
  };

  auto isFront = [&](const auto& hit) {
    return isFront_detail(hit.Subsystem(), hit.Station(), hit.Ring(), hit.Chamber(), (hit.Subsystem() == TriggerPrimitive::kRPC ? hit.Subsector_RPC() : hit.Subsector()));
  };

  auto calc_invpt = EMTFInversePt();

  auto calc_d0 = EMTFDZero();

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
    } else if (hit.Subsystem() == TriggerPrimitive::kDT) {
      sim_tp1 = truth_.findDTDigiSimLink(hit, trkParts_);
      sim_tp2 = sim_tp1;
    }
    return std::make_pair(sim_tp1, sim_tp2);
  };

  auto get_pattern = [](const auto& hit) {
    int pattern = 0;
    if (hit.Subsystem() == TriggerPrimitive::kCSC) {
      pattern = hit.Pattern();
    } else if (hit.Subsystem() == TriggerPrimitive::kDT) {
      pattern = hit.Sync_err();  // syncErr was hacked to store rpc bit
    }
    return pattern;
  };

  auto get_time = [](const auto& hit) {
    float time = hit.Time();
    return static_cast<int>(std::round(time*16/25));  // integer unit is 25ns/16 (4-bit)
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
          (conv_hit_i.Endcap()     == conv_hit_j.Endcap()) && // Needed only in the ntupler
          (conv_hit_i.Sector()     == conv_hit_j.Sector()) && // Needed only in the ntupler
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
  // Verbose
  if (firstEvent_)  edm::LogInfo("NtupleMaker") << "Ready to make ntuple.";

  if (verbose_ > 0) {
    std::cout << "[DEBUG] # hits: " << emuHits_.size() << " #  tracks: " << emuTracks_.size()
              << " # trk parts: " << trkParts_.size() << std::endl;
  }

  // ___________________________________________________________________________
  // Hits
  prepare_sim_tp();  // must be called before calling get_sim_tp_matches()

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
    vh_quality    ->push_back(hit.Quality());
    vh_pattern    ->push_back(get_pattern(hit));  // modified
    vh_bend       ->push_back(hit.Bend());
    vh_time       ->push_back(get_time(hit));     // modified
    vh_fr         ->push_back(isFront(hit));      // added
    vh_emtf_phi   ->push_back(hit.Phi_fp());
    vh_emtf_theta ->push_back(hit.Theta_fp());
    //
    vh_sim_phi    ->push_back(hit.Phi_sim());
    vh_sim_theta  ->push_back(hit.Theta_sim());
    vh_sim_r      ->push_back(hit.Rho_sim());
    vh_sim_z      ->push_back(hit.Z_sim());
    vh_sim_tp1    ->push_back(sim_tp_pair.first);
    vh_sim_tp2    ->push_back(sim_tp_pair.second);
  }
  (*vh_size) = emuHits_.size();

  // ___________________________________________________________________________
  // EndcapSimHits
  const std::vector<EMTFMCTruth::EndcapSimHit>& endcapSimHits = truth_.findEndcapSimHits();
  for (const auto& hit : endcapSimHits) {
    vc_type       ->push_back(hit.type);
    vc_station    ->push_back(hit.station);
    vc_ring       ->push_back(hit.ring);
    vc_layer      ->push_back(hit.layer);
    vc_chamber    ->push_back(hit.chamber);
    vc_phi        ->push_back(emtf::rad_to_deg(hit.globalPosition.phi()));
    vc_theta      ->push_back(emtf::rad_to_deg(hit.globalPosition.theta()));
    vc_r          ->push_back(hit.globalPosition.perp());
    vc_z          ->push_back(hit.globalPosition.z());
    vc_sim_tp     ->push_back(hit.sim_tp);
    vc_pdgid      ->push_back(hit.pSimHit.particleType());
    vc_process    ->push_back(hit.pSimHit.processType());
    vc_mom_phi    ->push_back(emtf::rad_to_deg(hit.pSimHit.phiAtEntry()));
    vc_mom_theta  ->push_back(emtf::rad_to_deg(hit.pSimHit.thetaAtEntry()));
    vc_tof        ->push_back(hit.pSimHit.timeOfFlight());
  }
  (*vc_size) = endcapSimHits.size();

  // ___________________________________________________________________________
  // Tracks
  for (const auto& trk : emuTracks_) {
    const std::vector<int32_t>& hit_refs = get_hit_refs(trk, emuHits_);
    assert(hit_refs.size() == 4);

    const l1t::EMTFPtLUT& ptlut_data = trk.PtLUT();

    vt_pt         ->push_back(trk.Pt());
    vt_xml_pt     ->push_back(trk.Pt_XML());
    vt_pt_dxy     ->push_back(trk.Pt_dxy());
    vt_dxy        ->push_back(trk.Dxy());
    vt_invpt_prompt ->push_back(trk.Invpt_prompt());
    vt_invpt_displ  ->push_back(trk.Invpt_displ());
    vt_phi        ->push_back(trk.Phi_glob());
    vt_theta      ->push_back(trk.Theta());
    vt_eta        ->push_back(trk.Eta());
    vt_q          ->push_back(trk.Charge());
    //
    vt_address    ->push_back(ptlut_data.address);
    vt_mode       ->push_back(trk.Mode());
    vt_endcap     ->push_back(trk.Endcap());
    vt_sector     ->push_back(trk.Sector());
    vt_bx         ->push_back(trk.BX());
    vt_nhits      ->push_back(trk.Hits().size());
    vt_hitref1    ->push_back(hit_refs.at(0));
    vt_hitref2    ->push_back(hit_refs.at(1));
    vt_hitref3    ->push_back(hit_refs.at(2));
    vt_hitref4    ->push_back(hit_refs.at(3));
  }
  (*vt_size) = emuTracks_.size();

  // ___________________________________________________________________________
  // L1TrackTrigger tracks
  if (tkTracks_ != nullptr && tkTrackAssoc_ != nullptr) {
    int itkTrack = 0;
    auto tkTracks_handle = make_handle(tkTracks_);

    for (const auto& trk : *tkTracks_) {
      const GlobalVector& momentum = trk.getMomentum();
      const GlobalPoint&  poca     = trk.getPOCA();
      const double        rinv     = trk.getRInv();

      vu_pt    ->push_back(momentum.perp());
      vu_phi   ->push_back(momentum.phi());
      vu_theta ->push_back(momentum.theta());
      vu_eta   ->push_back(momentum.eta());
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
        edm::Ptr<TrackingParticle> tpPtr = tkTrackAssoc_->findTrackingParticlePtr(trkPtr);
        if (tpPtr.isNonnull()) {
          sim_pt    = tpPtr->pt();
          sim_eta   = tpPtr->eta();
          sim_phi   = tpPtr->phi();
          sim_tp    = tpPtr.key();
          sim_pdgid = tpPtr->pdgId();
        }
        if (tkTrackAssoc_->isGenuine(trkPtr)) {
          sim_assoc = 0;
        } else if (tkTrackAssoc_->isLooselyGenuine(trkPtr)) {
          sim_assoc = 1;
        } else if (tkTrackAssoc_->isCombinatoric(trkPtr)) {
          sim_assoc = 2;
        } else if (tkTrackAssoc_->isUnknown(trkPtr)) {
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
  }
  (*vu_size) = vu_pt->size();

  // ___________________________________________________________________________
  // Gen particles
  //int igenPart = 0;
  //for (const auto& part : *genParts_) {
  //  vp_pt         ->push_back(part.pt());
  //  vp_phi        ->push_back(part.phi());
  //  vp_theta      ->push_back(part.theta());
  //  vp_eta        ->push_back(part.eta());
  //  vp_vx         ->push_back(part.vx());
  //  vp_vy         ->push_back(part.vy());
  //  vp_vz         ->push_back(part.vz());
  //  vp_invpt      ->push_back(calc_invpt(part.charge(), part.pt()));
  //  vp_d0         ->push_back(calc_d0(calc_invpt(part.charge(), part.pt()), part.phi(), part.vx(), part.vy()));
  //  vp_beta       ->push_back(part.p()/part.energy());
  //  vp_mass       ->push_back(part.mass());
  //  vp_q          ->push_back(part.charge());
  //  vp_bx         ->push_back(0);
  //  vp_event      ->push_back(0);
  //  vp_pdgid      ->push_back(part.pdgId());
  //  vp_status     ->push_back(part.status());
  //  vp_decay      ->push_back(0);
  //  vp_genp       ->push_back(igenPart);
  //
  //  ++igenPart;
  //}
  //(*vp_size) = genParts_->size();
  //assert(static_cast<size_t>(*vp_size) == vp_pt->size());

  // ___________________________________________________________________________
  // Sim tracks
  //int isimTrack = 0;
  //for (const auto& part : *simTracks_) {
  //
  //  ++isimTrack;
  //}
  //(*vp_size) = simTracks_->size();
  //assert(static_cast<size_t>(*vp_size) == vp_pt->size());

  // ___________________________________________________________________________
  // Tracking particles
  for (const auto& part : trkParts_) {
    int igenPart = -1;
    if (!part.genParticles().empty()) {
      igenPart = (part.genParticles().begin())->key();
    }

    vp_pt         ->push_back(part.pt());
    vp_phi        ->push_back(part.phi());
    vp_theta      ->push_back(part.theta());
    vp_eta        ->push_back(part.eta());
    vp_vx         ->push_back(part.vx());
    vp_vy         ->push_back(part.vy());
    vp_vz         ->push_back(part.vz());
    vp_invpt      ->push_back(calc_invpt(part.charge(), part.pt()));
    vp_d0         ->push_back(calc_d0(calc_invpt(part.charge(), part.pt()), part.phi(), part.vx(), part.vy()));
    vp_beta       ->push_back(part.p()/part.energy());
    vp_mass       ->push_back(part.mass());
    vp_q          ->push_back(part.charge());
    vp_bx         ->push_back(part.eventId().bunchCrossing());
    vp_event      ->push_back(part.eventId().event());
    vp_pdgid      ->push_back(part.pdgId());
    vp_status     ->push_back(part.status());
    vp_decay      ->push_back(part.decayVertices().size());
    vp_genp       ->push_back(igenPart);
  }
  (*vp_size) = trkParts_.size();
  assert(static_cast<size_t>(*vp_size) == vp_pt->size());

  // ___________________________________________________________________________
  // Event info
  if (pileupInfo_ != nullptr) {
    float true_npv = 0;
    int pu_nvtx = 0;

    for (const auto& pui : *pileupInfo_) {
      if (pui.getBunchCrossing() == 0) {
        true_npv = pui.getTrueNumInteractions();
        pu_nvtx = pui.getPU_NumInteractions();
        break;
      }
    }

    ve_event      ->push_back(iEvent.id().event());
    ve_run        ->push_back(iEvent.id().run());
    ve_lumi       ->push_back(iEvent.id().luminosityBlock());
    ve_npv        ->push_back(true_npv);
    ve_nvtx       ->push_back(pu_nvtx);
  }
  (*ve_size) = 1;

  // ___________________________________________________________________________
  // Fill
  tree->Fill();

  // Hits
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
  vh_quality    ->clear();
  vh_pattern    ->clear();
  vh_bend       ->clear();
  vh_time       ->clear();
  vh_fr         ->clear();
  vh_emtf_phi   ->clear();
  vh_emtf_theta ->clear();
  //
  vh_sim_phi    ->clear();
  vh_sim_theta  ->clear();
  vh_sim_r      ->clear();
  vh_sim_z      ->clear();
  vh_sim_tp1    ->clear();
  vh_sim_tp2    ->clear();
  (*vh_size)    = 0;

  // EndcapSimHits
  vc_type       ->clear();
  vc_station    ->clear();
  vc_ring       ->clear();
  vc_layer      ->clear();
  vc_chamber    ->clear();
  vc_phi        ->clear();
  vc_theta      ->clear();
  vc_r          ->clear();
  vc_z          ->clear();
  vc_sim_tp     ->clear();
  vc_pdgid      ->clear();
  vc_process    ->clear();
  vc_mom_phi    ->clear();
  vc_mom_theta  ->clear();
  vc_tof        ->clear();
  (*vc_size)    = 0;

  // Tracks
  vt_pt         ->clear();
  vt_xml_pt     ->clear();
  vt_pt_dxy     ->clear();
  vt_dxy        ->clear();
  vt_invpt_prompt ->clear();
  vt_invpt_displ  ->clear();
  vt_phi        ->clear();
  vt_theta      ->clear();
  vt_eta        ->clear();
  vt_q          ->clear();
  //
  vt_address    ->clear();
  vt_mode       ->clear();
  vt_endcap     ->clear();
  vt_sector     ->clear();
  vt_bx         ->clear();
  vt_nhits      ->clear();
  vt_hitref1    ->clear();
  vt_hitref2    ->clear();
  vt_hitref3    ->clear();
  vt_hitref4    ->clear();
  (*vt_size)    = 0;

  // L1TrackTrigger tracks
  vu_pt         ->clear();
  vu_phi        ->clear();
  vu_theta      ->clear();
  vu_eta        ->clear();
  vu_vx         ->clear();
  vu_vy         ->clear();
  vu_vz         ->clear();
  vu_q          ->clear();
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
  (*vu_size)    = 0;

  // Tracking particles
  vp_pt         ->clear();
  vp_phi        ->clear();
  vp_theta      ->clear();
  vp_eta        ->clear();
  vp_vx         ->clear();
  vp_vy         ->clear();
  vp_vz         ->clear();
  vp_invpt      ->clear();
  vp_d0         ->clear();
  vp_beta       ->clear();
  vp_mass       ->clear();
  vp_q          ->clear();
  vp_bx         ->clear();
  vp_event      ->clear();
  vp_pdgid      ->clear();
  vp_status     ->clear();
  vp_decay      ->clear();
  vp_genp       ->clear();
  (*vp_size)    = 0;

  // Event info
  ve_event      ->clear();
  ve_run        ->clear();
  ve_lumi       ->clear();
  ve_npv        ->clear();
  ve_nvtx       ->clear();
  (*ve_size)    = 0;

  if (firstEvent_)
    firstEvent_ = false;
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

  // Create pointers
  // Hits
  vh_endcap     = std::make_unique<std::vector<int16_t> >();
  vh_station    = std::make_unique<std::vector<int16_t> >();
  vh_ring       = std::make_unique<std::vector<int16_t> >();
  vh_sector     = std::make_unique<std::vector<int16_t> >();
  vh_subsector  = std::make_unique<std::vector<int16_t> >();
  vh_chamber    = std::make_unique<std::vector<int16_t> >();
  vh_cscid      = std::make_unique<std::vector<int16_t> >();
  vh_bx         = std::make_unique<std::vector<int16_t> >();
  vh_type       = std::make_unique<std::vector<int16_t> >();
  vh_neighbor   = std::make_unique<std::vector<int16_t> >();
  //
  vh_strip      = std::make_unique<std::vector<int16_t> >();
  vh_wire       = std::make_unique<std::vector<int16_t> >();
  vh_roll       = std::make_unique<std::vector<int16_t> >();
  vh_quality    = std::make_unique<std::vector<int16_t> >();
  vh_pattern    = std::make_unique<std::vector<int16_t> >();
  vh_bend       = std::make_unique<std::vector<int16_t> >();
  vh_time       = std::make_unique<std::vector<int16_t> >();
  vh_fr         = std::make_unique<std::vector<int16_t> >();
  vh_emtf_phi   = std::make_unique<std::vector<int32_t> >();
  vh_emtf_theta = std::make_unique<std::vector<int32_t> >();
  //
  vh_sim_phi    = std::make_unique<std::vector<float  > >();
  vh_sim_theta  = std::make_unique<std::vector<float  > >();
  vh_sim_r      = std::make_unique<std::vector<float  > >();
  vh_sim_z      = std::make_unique<std::vector<float  > >();
  vh_sim_tp1    = std::make_unique<std::vector<int32_t> >();
  vh_sim_tp2    = std::make_unique<std::vector<int32_t> >();
  vh_size       = std::make_unique<int32_t>(0);

  // EndcapSimHits
  vc_type       = std::make_unique<std::vector<int16_t> >();
  vc_station    = std::make_unique<std::vector<int16_t> >();
  vc_ring       = std::make_unique<std::vector<int16_t> >();
  vc_layer      = std::make_unique<std::vector<int16_t> >();
  vc_chamber    = std::make_unique<std::vector<int16_t> >();
  vc_phi        = std::make_unique<std::vector<float  > >();
  vc_theta      = std::make_unique<std::vector<float  > >();
  vc_r          = std::make_unique<std::vector<float  > >();
  vc_z          = std::make_unique<std::vector<float  > >();
  vc_sim_tp     = std::make_unique<std::vector<int32_t> >();
  vc_pdgid      = std::make_unique<std::vector<int32_t> >();
  vc_process    = std::make_unique<std::vector<int16_t> >();
  vc_mom_phi    = std::make_unique<std::vector<float  > >();
  vc_mom_theta  = std::make_unique<std::vector<float  > >();
  vc_tof        = std::make_unique<std::vector<float  > >();
  vc_size       = std::make_unique<int32_t>(0);

  // Tracks
  vt_pt         = std::make_unique<std::vector<float  > >();
  vt_xml_pt     = std::make_unique<std::vector<float  > >();
  vt_pt_dxy     = std::make_unique<std::vector<float  > >();
  vt_dxy        = std::make_unique<std::vector<float  > >();
  vt_invpt_prompt = std::make_unique<std::vector<float  > >();
  vt_invpt_displ  = std::make_unique<std::vector<float  > >();
  vt_phi        = std::make_unique<std::vector<float  > >();
  vt_theta      = std::make_unique<std::vector<float  > >();
  vt_eta        = std::make_unique<std::vector<float  > >();
  vt_q          = std::make_unique<std::vector<int16_t> >();
  //
  vt_address    = std::make_unique<std::vector<uint64_t> >();
  vt_mode       = std::make_unique<std::vector<int16_t> >();
  vt_endcap     = std::make_unique<std::vector<int16_t> >();
  vt_sector     = std::make_unique<std::vector<int16_t> >();
  vt_bx         = std::make_unique<std::vector<int16_t> >();
  vt_nhits      = std::make_unique<std::vector<int16_t> >();
  vt_hitref1    = std::make_unique<std::vector<int32_t> >();
  vt_hitref2    = std::make_unique<std::vector<int32_t> >();
  vt_hitref3    = std::make_unique<std::vector<int32_t> >();
  vt_hitref4    = std::make_unique<std::vector<int32_t> >();
  vt_size       = std::make_unique<int32_t>(0);

  // L1TrackTrigger trakcs
  vu_pt         = std::make_unique<std::vector<float  > >();
  vu_phi        = std::make_unique<std::vector<float  > >();
  vu_theta      = std::make_unique<std::vector<float  > >();
  vu_eta        = std::make_unique<std::vector<float  > >();
  vu_vx         = std::make_unique<std::vector<float  > >();
  vu_vy         = std::make_unique<std::vector<float  > >();
  vu_vz         = std::make_unique<std::vector<float  > >();
  vu_q          = std::make_unique<std::vector<int16_t> >();
  vu_rinv       = std::make_unique<std::vector<float  > >();
  vu_chi2       = std::make_unique<std::vector<float  > >();
  vu_ndof       = std::make_unique<std::vector<int16_t> >();
  vu_sector     = std::make_unique<std::vector<int16_t> >();
  //
  vu_sim_pt     = std::make_unique<std::vector<float  > >();
  vu_sim_phi    = std::make_unique<std::vector<float  > >();
  vu_sim_eta    = std::make_unique<std::vector<float  > >();
  vu_sim_tp     = std::make_unique<std::vector<int32_t> >();
  vu_sim_pdgid  = std::make_unique<std::vector<int32_t> >();
  vu_sim_assoc  = std::make_unique<std::vector<int16_t> >();
  vu_size       = std::make_unique<int32_t>(0);

  // Tracking particles
  vp_pt         = std::make_unique<std::vector<float  > >();
  vp_phi        = std::make_unique<std::vector<float  > >();
  vp_theta      = std::make_unique<std::vector<float  > >();
  vp_eta        = std::make_unique<std::vector<float  > >();
  vp_vx         = std::make_unique<std::vector<float  > >();
  vp_vy         = std::make_unique<std::vector<float  > >();
  vp_vz         = std::make_unique<std::vector<float  > >();
  vp_invpt      = std::make_unique<std::vector<float  > >();
  vp_d0         = std::make_unique<std::vector<float  > >();
  vp_beta       = std::make_unique<std::vector<float  > >();
  vp_mass       = std::make_unique<std::vector<float  > >();
  vp_q          = std::make_unique<std::vector<int16_t> >();
  vp_bx         = std::make_unique<std::vector<int16_t> >();
  vp_event      = std::make_unique<std::vector<int16_t> >();
  vp_pdgid      = std::make_unique<std::vector<int32_t> >();
  vp_status     = std::make_unique<std::vector<int16_t> >();
  vp_decay      = std::make_unique<std::vector<int16_t> >();
  vp_genp       = std::make_unique<std::vector<int32_t> >();
  vp_size       = std::make_unique<int32_t>(0);

  // Event info
  ve_event      = std::make_unique<std::vector<uint64_t> >();
  ve_run        = std::make_unique<std::vector<uint32_t> >();
  ve_lumi       = std::make_unique<std::vector<uint32_t> >();
  ve_npv        = std::make_unique<std::vector<float  > >();
  ve_nvtx       = std::make_unique<std::vector<int32_t> >();
  ve_size       = std::make_unique<int32_t>(0);


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
  tree->Branch("vh_quality"   , &(*vh_quality   ));
  tree->Branch("vh_pattern"   , &(*vh_pattern   ));
  tree->Branch("vh_bend"      , &(*vh_bend      ));
  tree->Branch("vh_time"      , &(*vh_time      ));
  tree->Branch("vh_fr"        , &(*vh_fr        ));
  tree->Branch("vh_emtf_phi"  , &(*vh_emtf_phi  ));
  tree->Branch("vh_emtf_theta", &(*vh_emtf_theta));
  //
  tree->Branch("vh_sim_phi"   , &(*vh_sim_phi   ));
  tree->Branch("vh_sim_theta" , &(*vh_sim_theta ));
  tree->Branch("vh_sim_r"     , &(*vh_sim_r     ));
  tree->Branch("vh_sim_z"     , &(*vh_sim_z     ));
  tree->Branch("vh_sim_tp1"   , &(*vh_sim_tp1   ));
  tree->Branch("vh_sim_tp2"   , &(*vh_sim_tp2   ));
  tree->Branch("vh_size"      , &(*vh_size      ));

  // EndcapSimHits
  tree->Branch("vc_type"      , &(*vc_type      ));
  tree->Branch("vc_station"   , &(*vc_station   ));
  tree->Branch("vc_ring"      , &(*vc_ring      ));
  tree->Branch("vc_layer"     , &(*vc_layer     ));
  tree->Branch("vc_chamber"   , &(*vc_chamber   ));
  tree->Branch("vc_phi"       , &(*vc_phi       ));
  tree->Branch("vc_theta"     , &(*vc_theta     ));
  tree->Branch("vc_r"         , &(*vc_r         ));
  tree->Branch("vc_z"         , &(*vc_z         ));
  tree->Branch("vc_sim_tp"    , &(*vc_sim_tp    ));
  tree->Branch("vc_pdgid"     , &(*vc_pdgid     ));
  tree->Branch("vc_process"   , &(*vc_process   ));
  tree->Branch("vc_mom_phi"   , &(*vc_mom_phi   ));
  tree->Branch("vc_mom_theta" , &(*vc_mom_theta ));
  tree->Branch("vc_tof"       , &(*vc_tof       ));
  tree->Branch("vc_size"      , &(*vc_size      ));

  // Tracks
  tree->Branch("vt_pt"        , &(*vt_pt        ));
  tree->Branch("vt_xml_pt"    , &(*vt_xml_pt    ));
  tree->Branch("vt_pt_dxy"    , &(*vt_pt_dxy    ));
  tree->Branch("vt_dxy"       , &(*vt_dxy       ));
  tree->Branch("vt_invpt_prompt", &(*vt_invpt_prompt));
  tree->Branch("vt_invpt_displ" , &(*vt_invpt_displ ));
  tree->Branch("vt_phi"       , &(*vt_phi       ));
  tree->Branch("vt_theta"     , &(*vt_theta     ));
  tree->Branch("vt_eta"       , &(*vt_eta       ));
  tree->Branch("vt_q"         , &(*vt_q         ));
  //
  tree->Branch("vt_address"   , &(*vt_address   ));
  tree->Branch("vt_mode"      , &(*vt_mode      ));
  tree->Branch("vt_endcap"    , &(*vt_endcap    ));
  tree->Branch("vt_sector"    , &(*vt_sector    ));
  tree->Branch("vt_bx"        , &(*vt_bx        ));
  tree->Branch("vt_nhits"     , &(*vt_nhits     ));
  tree->Branch("vt_hitref1"   , &(*vt_hitref1   ));
  tree->Branch("vt_hitref2"   , &(*vt_hitref2   ));
  tree->Branch("vt_hitref3"   , &(*vt_hitref3   ));
  tree->Branch("vt_hitref4"   , &(*vt_hitref4   ));
  tree->Branch("vt_size"      , &(*vt_size      ));

  // L1TrackTrigger tracks
  tree->Branch("vu_pt"        , &(*vu_pt        ));
  tree->Branch("vu_phi"       , &(*vu_phi       ));
  tree->Branch("vu_theta"     , &(*vu_theta     ));
  tree->Branch("vu_eta"       , &(*vu_eta       ));
  tree->Branch("vu_vx"        , &(*vu_vx        ));
  tree->Branch("vu_vy"        , &(*vu_vy        ));
  tree->Branch("vu_vz"        , &(*vu_vz        ));
  tree->Branch("vu_q"         , &(*vu_q         ));
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
  tree->Branch("vu_size"      , &(*vu_size      ));

  // Tracking particles
  tree->Branch("vp_pt"        , &(*vp_pt        ));
  tree->Branch("vp_phi"       , &(*vp_phi       ));
  tree->Branch("vp_theta"     , &(*vp_theta     ));
  tree->Branch("vp_eta"       , &(*vp_eta       ));
  tree->Branch("vp_vx"        , &(*vp_vx        ));
  tree->Branch("vp_vy"        , &(*vp_vy        ));
  tree->Branch("vp_vz"        , &(*vp_vz        ));
  tree->Branch("vp_invpt"     , &(*vp_invpt     ));
  tree->Branch("vp_d0"        , &(*vp_d0        ));
  tree->Branch("vp_beta"      , &(*vp_beta      ));
  tree->Branch("vp_mass"      , &(*vp_mass      ));
  tree->Branch("vp_q"         , &(*vp_q         ));
  tree->Branch("vp_bx"        , &(*vp_bx        ));
  tree->Branch("vp_event"     , &(*vp_event     ));
  tree->Branch("vp_pdgid"     , &(*vp_pdgid     ));
  tree->Branch("vp_status"    , &(*vp_status    ));
  tree->Branch("vp_decay"     , &(*vp_decay     ));
  tree->Branch("vp_genp"      , &(*vp_genp      ));
  tree->Branch("vp_size"      , &(*vp_size      ));

  // Event info
  tree->Branch("ve_event"     , &(*ve_event     ));
  tree->Branch("ve_run"       , &(*ve_run       ));
  tree->Branch("ve_lumi"      , &(*ve_lumi      ));
  tree->Branch("ve_npv"       , &(*ve_npv       ));
  tree->Branch("ve_nvtx"      , &(*ve_nvtx      ));
  tree->Branch("ve_size"      , &(*ve_size      ));
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
