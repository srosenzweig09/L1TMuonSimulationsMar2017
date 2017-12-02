#include "L1TMuonSimulations/Analyzers/interface/EMTFMCTruth.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
//#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include <iostream>


EMTFMCTruth::EMTFMCTruth(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iConsumes) :
    cscSimHitsTag_      (iConfig.getParameter<edm::InputTag>("cscSimHitsTag")),
    cscSimHitsXFTag_    (iConfig.getParameter<edm::InputTag>("cscSimHitsXFTag")),
    cscStripSimLinksTag_(iConfig.getParameter<edm::InputTag>("cscStripSimLinksTag")),
    cscWireSimLinksTag_ (iConfig.getParameter<edm::InputTag>("cscWireSimLinksTag")),
    rpcSimHitsTag_      (iConfig.getParameter<edm::InputTag>("rpcSimHitsTag")),
    rpcSimHitsXFTag_    (iConfig.getParameter<edm::InputTag>("rpcSimHitsXFTag")),
    rpcDigiSimLinksTag_ (iConfig.getParameter<edm::InputTag>("rpcDigiSimLinksTag")),
    gemSimHitsTag_      (iConfig.getParameter<edm::InputTag>("gemSimHitsTag")),
    gemSimHitsXFTag_    (iConfig.getParameter<edm::InputTag>("gemSimHitsXFTag")),
    gemDigiSimLinksTag_ (iConfig.getParameter<edm::InputTag>("gemDigiSimLinksTag"))
{
  cscSimHitsToken_       = iConsumes.consumes<edm::PSimHitContainer>     (cscSimHitsTag_);
  cscSimHitsXFToken_     = iConsumes.consumes<CrossingFrame<PSimHit> >   (cscSimHitsXFTag_);
  cscStripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>         (cscStripSimLinksTag_);
  cscWireSimLinksToken_  = iConsumes.consumes<WireDigiSimLinks>          (cscWireSimLinksTag_);
  rpcSimHitsToken_       = iConsumes.consumes<edm::PSimHitContainer>     (rpcSimHitsTag_);
  rpcSimHitsXFToken_     = iConsumes.consumes<CrossingFrame<PSimHit> >   (rpcSimHitsXFTag_);
  rpcDigiSimLinksToken_  = iConsumes.consumes<RPCDigiSimLinks>           (rpcDigiSimLinksTag_);
  gemSimHitsToken_       = iConsumes.consumes<edm::PSimHitContainer>     (gemSimHitsTag_);
  gemSimHitsXFToken_     = iConsumes.consumes<CrossingFrame<PSimHit> >   (gemSimHitsXFTag_);
  gemDigiSimLinksToken_  = iConsumes.consumes<GEMDigiSimLinks>           (gemDigiSimLinksTag_);
}

EMTFMCTruth::~EMTFMCTruth() {}

void EMTFMCTruth::initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //// SimHits
  //edm::Handle<edm::PSimHitContainer> cscSimHits_handle;
  //
  //if (!iEvent.isRealData()) {
  //  if (!cscSimHitsToken_.isUninitialized()) {
  //    iEvent.getByToken(cscSimHitsToken_, cscSimHits_handle);
  //  }
  //  if (!cscSimHits_handle.isValid()) {
  //    edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscSimHitsTag_;
  //  }
  //}
  //
  //cscSimHitsPtr_ = cscSimHits_handle.product();


  //// SimHits (using crossing frame)
  //edm::Handle<CrossingFrame<PSimHit> > cscSimHitsXF_handle;
  //
  //if (!iEvent.isRealData()) {
  //  if (!cscSimHitsXFToken_.isUninitialized()) {
  //    iEvent.getByToken(cscSimHitsXFToken_, cscSimHitsXF_handle);
  //  }
  //  if (!cscSimHitsXF_handle.isValid()) {
  //    edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscSimHitsXFTag_;
  //  }
  //}
  //
  //cscSimHitsXFPtr_ = cscSimHitsXF_handle.product();


  // SimLinks
  edm::Handle<StripDigiSimLinks> cscStripSimLinks_handle;
  edm::Handle<WireDigiSimLinks>  cscWireSimLinks_handle;
  edm::Handle<RPCDigiSimLinks>   rpcDigiSimLinks_handle;
  edm::Handle<GEMDigiSimLinks>   gemDigiSimLinks_handle;

  if (!iEvent.isRealData()) {
    if (cscStripSimLinksTag_.encode() != "") {
      if (!cscStripSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(cscStripSimLinksToken_, cscStripSimLinks_handle);
      }
      if (!cscStripSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscStripSimLinksTag_;
      }
    }

    if (cscWireSimLinksTag_.encode() != "") {
      if (!cscWireSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(cscWireSimLinksToken_, cscWireSimLinks_handle);
      }
      if (!cscWireSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscWireSimLinksTag_;
      }
    }

    if (rpcDigiSimLinksTag_.encode() != "") {
      if (!rpcDigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(rpcDigiSimLinksToken_, rpcDigiSimLinks_handle);
      }
      if (!rpcDigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << rpcDigiSimLinksTag_;
      }
    }

    if (gemDigiSimLinksTag_.encode() != "") {
      if (!gemDigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(gemDigiSimLinksToken_, gemDigiSimLinks_handle);
      }
      if (!gemDigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << gemDigiSimLinksTag_;
      }
    }
  }

  cscStripSimLinksPtr_ = cscStripSimLinks_handle.product();
  cscWireSimLinksPtr_  = cscWireSimLinks_handle.product();
  rpcDigiSimLinksPtr_  = rpcDigiSimLinks_handle.product();
  gemDigiSimLinksPtr_  = gemDigiSimLinks_handle.product();

}

void EMTFMCTruth::makeTrackingParticleLinks(const TrackingParticleCollection& trkPartColl) {
  trackingParticleLinks_.clear();

  for (TrackingParticleCollection::const_iterator it_trkpart = trkPartColl.begin(); it_trkpart != trkPartColl.end(); ++it_trkpart) {
    for (TrackingParticle::g4t_iterator it_simtrk = it_trkpart->g4Track_begin(); it_simtrk != it_trkpart->g4Track_end(); ++it_simtrk) {
      unsigned int simTrackId = it_simtrk->trackId();
      EncodedEventId eventId = it_simtrk->eventId();
      SimHitIdpr matchId(simTrackId, eventId);
      //assert(trackingParticleLinks_.find(matchId) == trackingParticleLinks_.end());  // matchId should not already exist
      trackingParticleLinks_[matchId] = std::distance(trkPartColl.begin(), it_trkpart);
    }
  }

}


int EMTFMCTruth::findCSCStripSimLink(const l1t::EMTFHit& hit, const std::vector<TrackingParticle>& trkPartColl) const {
  std::map<SimHitIdpr, int> matches;

  // Check all 6 CSC layers. Each layer has a distinct detId
  for (unsigned csclayer = 0; csclayer < 6; ++csclayer) {
    CSCDetId cscDetId0 = hit.CSC_DetId();
    CSCDetId cscDetId1(cscDetId0.endcap(), cscDetId0.station(), cscDetId0.ring(), cscDetId0.chamber(), csclayer+1);
    int strip0 = hit.Strip();
    int strip1 = (strip0 - 1)/2 + 1;  // different convention used in CSC StripDigiSimLink

    StripDigiSimLinks::const_iterator cscStripLayerLinks = cscStripSimLinksPtr_->find(cscDetId1);
    if (cscStripLayerLinks != cscStripSimLinksPtr_->end()) {
      for (LayerLinks::const_iterator linkItr = cscStripLayerLinks->begin(); linkItr != cscStripLayerLinks->end(); ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (fraction > 0.1 && std::abs(int(strip1) - int(channel)) <= 2) {  // allow +/-2
          SimHitIdpr matchId(simTrackId, eventId);
          ++matches[matchId];
        }
      }
    }
  }

  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findCSCWireSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, int> matches;

  // Check all 6 CSC layers. Each layer has a distinct detId
  for (unsigned csclayer = 0; csclayer < 6; ++csclayer) {
    CSCDetId cscDetId0 = hit.CSC_DetId();
    CSCDetId cscDetId1(cscDetId0.endcap(), cscDetId0.station(), cscDetId0.ring(), cscDetId0.chamber(), csclayer+1);
    int wire0 = hit.Wire();
    int wire1 = (wire0 + 100) + 1;  // different convention used in CSC StripDigiSimLink

    WireDigiSimLinks::const_iterator cscWireLayerLinks = cscWireSimLinksPtr_->find(cscDetId1);
    if (cscWireLayerLinks != cscWireSimLinksPtr_->end()) {
      for (LayerLinks::const_iterator linkItr = cscWireLayerLinks->begin(); linkItr != cscWireLayerLinks->end(); ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (fraction > 0.1 && std::abs(int(wire1) - int(channel)) <= 1) {  // allow +/-1
          SimHitIdpr matchId(simTrackId, eventId);
          ++matches[matchId];
        }
      }
    }
  }

  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findRPCDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, int> matches;

  // Check all strips in the RPC cluster
  RPCDetId rpcDetId = hit.RPC_DetId();
  int stripA = hit.Strip_low();
  int stripB = hit.Strip_hi();
  int bx     = hit.BX();

#if 0
  // This is giving me an error: Assertion `std::distance(p.first, p.second) == 1' failed.
  RPCDigiSimLinks::const_iterator rpcDigiLayerLinks = rpcDigiSimLinksPtr_->find(rpcDetId);
  if (rpcDigiLayerLinks != rpcDigiSimLinksPtr_->end()) {
    for (RPCLayerLinks::const_iterator linkItr = rpcDigiLayerLinks->begin(); linkItr != rpcDigiLayerLinks->end(); ++linkItr) {
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      for (int strip1 = stripA; strip1 < stripB+1; ++strip1) {
        if (((int) simStrip == strip1) && ((int) simBX == bx)) {
          SimHitIdpr matchId(simTrackId, eventId);
          ++matches[matchId];
        }
      }
    }
  }
#else
  // My temporary fix
  for (RPCDigiSimLinks::const_iterator linkItr1 = rpcDigiSimLinksPtr_->begin(); linkItr1 != rpcDigiSimLinksPtr_->end(); ++linkItr1) {
    for (RPCLayerLinks::const_iterator linkItr = linkItr1->begin(); linkItr != linkItr1->end(); ++linkItr) {
      unsigned int detUnitId = linkItr->getDetUnitId();
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      if (detUnitId == rpcDetId.rawId()) {
        for (int strip0 = stripA; strip0 < stripB+1; ++strip0) {
          int strip1 = strip0;  // same convention
          if (((int) simStrip == strip1) && ((int) simBX == bx)) {
            SimHitIdpr matchId(simTrackId, eventId);
            ++matches[matchId];
          }
        }
      }
    }
  }
#endif

  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findGEMDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, int> matches;

  // Check all strips in the GEM cluster
  // FIXME: check both GEM layers
  GEMDetId gemDetId = hit.GEM_DetId();
  int stripA = hit.Strip_low();
  int stripB = hit.Strip_hi();
  int bx     = hit.BX();

  GEMDigiSimLinks::const_iterator gemDigiLayerLinks = gemDigiSimLinksPtr_->find(gemDetId);
  if (gemDigiLayerLinks != gemDigiSimLinksPtr_->end()) {
    for (GEMLayerLinks::const_iterator linkItr = gemDigiLayerLinks->begin(); linkItr != gemDigiLayerLinks->end(); ++linkItr) {
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      for (int strip0 = stripA; strip0 < stripB+1; ++strip0) {
        int strip1 = (strip0 - 1)/2 + 1;  // different convention used in GEMDigiSimLink
        if (((int) simStrip == strip1) && ((int) simBX == bx)) {
          SimHitIdpr matchId(simTrackId, eventId);
          ++matches[matchId];
        }
      }
    }
  }

  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findTrackingParticle(const std::map<SimHitIdpr, int>& matches, const TrackingParticleCollection& trkPartColl) const {

#if 0
  // Return highest pT
  int best_tp = -1;  // index of the tp
  double highest_pt = 0.;
  for (std::map<SimHitIdpr, int>::const_iterator it_match = matches.begin(); it_match != matches.end(); ++it_match) {
    auto found = trackingParticleLinks_.find(it_match->first);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      auto it_trkpart = trkPartColl.begin();
      std::advance(it_trkpart, tp);
      if (highest_pt < it_trkpart->pt()) {
        highest_pt = it_trkpart->pt();
        best_tp = tp;
      }
    }
  }
#else
  // Return majority
  int best_tp = -1;  // index of the tp
  int majority = 0;
  for (std::map<SimHitIdpr, int>::const_iterator it_match = matches.begin(); it_match != matches.end(); ++it_match) {
    auto found = trackingParticleLinks_.find(it_match->first);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      auto it_trkpart = trkPartColl.begin();
      std::advance(it_trkpart, tp);
      if (majority < it_match->second) {
        majority = it_match->second;
        best_tp = tp;
      }
    }
  }
#endif

  return best_tp;
}
