#include "L1TMuonSimulations/Analyzers/interface/EMTFMCTruth.h"

#include "DataFormats/MuonDetId/interface/DTWireId.h"
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
    gemDigiSimLinksTag_ (iConfig.getParameter<edm::InputTag>("gemDigiSimLinksTag")),
    gemStripSimLinksTag_(iConfig.getParameter<edm::InputTag>("gemStripSimLinksTag")),
    me0SimHitsTag_      (iConfig.getParameter<edm::InputTag>("me0SimHitsTag")),
    me0SimHitsXFTag_    (iConfig.getParameter<edm::InputTag>("me0SimHitsXFTag")),
    me0DigiSimLinksTag_ (iConfig.getParameter<edm::InputTag>("me0DigiSimLinksTag")),
    me0StripSimLinksTag_(iConfig.getParameter<edm::InputTag>("me0StripSimLinksTag")),
    dtSimHitsTag_       (iConfig.getParameter<edm::InputTag>("dtSimHitsTag")),
    dtSimHitsXFTag_     (iConfig.getParameter<edm::InputTag>("dtSimHitsXFTag")),
    dtDigiSimLinksTag_  (iConfig.getParameter<edm::InputTag>("dtDigiSimLinksTag")),
    crossingFrame_      (iConfig.getParameter<bool>         ("crossingFrame"))
{
  if (!crossingFrame_) {
    cscSimHitsToken_     = iConsumes.consumes<edm::PSimHitContainer>     (cscSimHitsTag_);
  } else {
    cscSimHitsXFToken_   = iConsumes.consumes<CrossingFrame<PSimHit> >   (cscSimHitsXFTag_);
  }
  cscStripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>         (cscStripSimLinksTag_);
  cscWireSimLinksToken_  = iConsumes.consumes<WireDigiSimLinks>          (cscWireSimLinksTag_);

  if (!crossingFrame_) {
    rpcSimHitsToken_     = iConsumes.consumes<edm::PSimHitContainer>     (rpcSimHitsTag_);
  } else {
    rpcSimHitsXFToken_   = iConsumes.consumes<CrossingFrame<PSimHit> >   (rpcSimHitsXFTag_);
  }
  rpcDigiSimLinksToken_  = iConsumes.consumes<RPCDigiSimLinks>           (rpcDigiSimLinksTag_);

  if (!crossingFrame_) {
    gemSimHitsToken_     = iConsumes.consumes<edm::PSimHitContainer>     (gemSimHitsTag_);
  } else {
    gemSimHitsXFToken_   = iConsumes.consumes<CrossingFrame<PSimHit> >   (gemSimHitsXFTag_);
  }
  gemDigiSimLinksToken_  = iConsumes.consumes<GEMDigiSimLinks>           (gemDigiSimLinksTag_);
  gemStripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>         (gemStripSimLinksTag_);

  if (!crossingFrame_) {
    me0SimHitsToken_     = iConsumes.consumes<edm::PSimHitContainer>     (me0SimHitsTag_);
  } else {
    me0SimHitsXFToken_   = iConsumes.consumes<CrossingFrame<PSimHit> >   (me0SimHitsXFTag_);
  }
  me0DigiSimLinksToken_  = iConsumes.consumes<ME0DigiSimLinks>           (me0DigiSimLinksTag_);
  me0StripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>         (me0StripSimLinksTag_);

  if (!crossingFrame_) {
    dtSimHitsToken_      = iConsumes.consumes<edm::PSimHitContainer>     (dtSimHitsTag_);
  } else {
    dtSimHitsXFToken_    = iConsumes.consumes<CrossingFrame<PSimHit> >   (dtSimHitsXFTag_);
  }
  dtDigiSimLinksToken_   = iConsumes.consumes<DTDigiSimLinks>            (dtDigiSimLinksTag_);
}

EMTFMCTruth::~EMTFMCTruth() {}

void EMTFMCTruth::initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // PSimHits
  edm::Handle<edm::PSimHitContainer> cscSimHits_handle;
  edm::Handle<edm::PSimHitContainer> rpcSimHits_handle;
  edm::Handle<edm::PSimHitContainer> gemSimHits_handle;
  edm::Handle<edm::PSimHitContainer> me0SimHits_handle;
  edm::Handle<edm::PSimHitContainer> dtSimHits_handle;

  if (!crossingFrame_) {
    if (!iEvent.isRealData()) {
      if (!cscSimHitsToken_.isUninitialized()) {
        iEvent.getByToken(cscSimHitsToken_, cscSimHits_handle);
      }
      if (!cscSimHits_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscSimHitsTag_;
      }

      if (!rpcSimHitsToken_.isUninitialized()) {
        iEvent.getByToken(rpcSimHitsToken_, rpcSimHits_handle);
      }
      if (!rpcSimHits_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << rpcSimHitsTag_;
      }

      if (!gemSimHitsToken_.isUninitialized()) {
        iEvent.getByToken(gemSimHitsToken_, gemSimHits_handle);
      }
      if (!gemSimHits_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << gemSimHitsTag_;
      }

      if (!me0SimHitsToken_.isUninitialized()) {
        iEvent.getByToken(me0SimHitsToken_, me0SimHits_handle);
      }
      if (!me0SimHits_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << me0SimHitsTag_;
      }

      if (!dtSimHitsToken_.isUninitialized()) {
        iEvent.getByToken(dtSimHitsToken_, dtSimHits_handle);
      }
      if (!dtSimHits_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << dtSimHitsTag_;
      }

      cscSimHitsPtr_ = cscSimHits_handle.product();
      rpcSimHitsPtr_ = rpcSimHits_handle.product();
      gemSimHitsPtr_ = gemSimHits_handle.product();
      me0SimHitsPtr_ = me0SimHits_handle.product();
      dtSimHitsPtr_ = dtSimHits_handle.product();
    }
  }

  // PSimHits (with crossing frame)
  edm::Handle<CrossingFrame<PSimHit> > cscSimHitsXF_handle;
  edm::Handle<CrossingFrame<PSimHit> > rpcSimHitsXF_handle;
  edm::Handle<CrossingFrame<PSimHit> > gemSimHitsXF_handle;
  edm::Handle<CrossingFrame<PSimHit> > me0SimHitsXF_handle;
  edm::Handle<CrossingFrame<PSimHit> > dtSimHitsXF_handle;

  if (crossingFrame_) {
    if (!iEvent.isRealData()) {
      if (!cscSimHitsXFToken_.isUninitialized()) {
        iEvent.getByToken(cscSimHitsXFToken_, cscSimHitsXF_handle);
      }
      if (!cscSimHitsXF_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscSimHitsXFTag_;
      }

      if (!rpcSimHitsXFToken_.isUninitialized()) {
        iEvent.getByToken(rpcSimHitsXFToken_, rpcSimHitsXF_handle);
      }
      if (!rpcSimHitsXF_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << rpcSimHitsXFTag_;
      }

      if (!gemSimHitsXFToken_.isUninitialized()) {
        iEvent.getByToken(gemSimHitsXFToken_, gemSimHitsXF_handle);
      }
      if (!gemSimHitsXF_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << gemSimHitsXFTag_;
      }

      if (!me0SimHitsXFToken_.isUninitialized()) {
        iEvent.getByToken(me0SimHitsXFToken_, me0SimHitsXF_handle);
      }
      if (!me0SimHitsXF_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << me0SimHitsXFTag_;
      }

      if (!dtSimHitsXFToken_.isUninitialized()) {
        iEvent.getByToken(dtSimHitsXFToken_, dtSimHitsXF_handle);
      }
      if (!dtSimHitsXF_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << dtSimHitsXFTag_;
      }

      cscSimHitsXFPtr_ = cscSimHitsXF_handle.product();
      rpcSimHitsXFPtr_ = rpcSimHitsXF_handle.product();
      gemSimHitsXFPtr_ = gemSimHitsXF_handle.product();
      me0SimHitsXFPtr_ = me0SimHitsXF_handle.product();
      dtSimHitsXFPtr_ = dtSimHitsXF_handle.product();
    }
  }

  // SimLinks
  edm::Handle<StripDigiSimLinks> cscStripSimLinks_handle;
  edm::Handle<WireDigiSimLinks>  cscWireSimLinks_handle;
  edm::Handle<RPCDigiSimLinks>   rpcDigiSimLinks_handle;
  edm::Handle<GEMDigiSimLinks>   gemDigiSimLinks_handle;
  edm::Handle<ME0DigiSimLinks>   me0DigiSimLinks_handle;
  edm::Handle<DTDigiSimLinks>    dtDigiSimLinks_handle;

  {
    if (!iEvent.isRealData()) {
      if (!cscStripSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(cscStripSimLinksToken_, cscStripSimLinks_handle);
      }
      if (!cscStripSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscStripSimLinksTag_;
      }

      if (!cscWireSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(cscWireSimLinksToken_, cscWireSimLinks_handle);
      }
      if (!cscWireSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << cscWireSimLinksTag_;
      }

      if (!rpcDigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(rpcDigiSimLinksToken_, rpcDigiSimLinks_handle);
      }
      if (!rpcDigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << rpcDigiSimLinksTag_;
      }

      if (!gemDigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(gemDigiSimLinksToken_, gemDigiSimLinks_handle);
      }
      if (!gemDigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << gemDigiSimLinksTag_;
      }

      if (!me0DigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(me0DigiSimLinksToken_, me0DigiSimLinks_handle);
      }
      if (!me0DigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << me0DigiSimLinksTag_;
      }

      if (!dtDigiSimLinksToken_.isUninitialized()) {
        iEvent.getByToken(dtDigiSimLinksToken_, dtDigiSimLinks_handle);
      }
      if (!dtDigiSimLinks_handle.isValid()) {
        edm::LogError("EMTFMCTruth") << "Cannot get the product: " << dtDigiSimLinksTag_;
      }

      cscStripSimLinksPtr_ = cscStripSimLinks_handle.product();
      cscWireSimLinksPtr_ = cscWireSimLinks_handle.product();
      rpcDigiSimLinksPtr_ = rpcDigiSimLinks_handle.product();
      gemDigiSimLinksPtr_ = gemDigiSimLinks_handle.product();
      me0DigiSimLinksPtr_ = me0DigiSimLinks_handle.product();
      dtDigiSimLinksPtr_ = dtDigiSimLinks_handle.product();
    }
  }

  // Geometry
  edm::ESHandle<CSCGeometry> cscGeom_handle;
  edm::ESHandle<RPCGeometry> rpcGeom_handle;
  edm::ESHandle<GEMGeometry> gemGeom_handle;
  edm::ESHandle<ME0Geometry> me0Geom_handle;
  edm::ESHandle<DTGeometry>  dtGeom_handle;

  {
    const MuonGeometryRecord& geom = iSetup.get<MuonGeometryRecord>();
    geom.get(cscGeom_handle);
    geom.get(rpcGeom_handle);
    geom.get(gemGeom_handle);
    geom.get(me0Geom_handle);
    geom.get(dtGeom_handle);

    cscGeomPtr_ = cscGeom_handle.product();
    rpcGeomPtr_ = rpcGeom_handle.product();
    gemGeomPtr_ = gemGeom_handle.product();
    me0GeomPtr_ = me0Geom_handle.product();
    dtGeomPtr_ = dtGeom_handle.product();
  }
  return;
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
  return;
}


int EMTFMCTruth::findCSCStripSimLink(const l1t::EMTFHit& hit, const std::vector<TrackingParticle>& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int strip0 = hit.Strip();
  int strip1 = (strip0 - 1)/2 + 1;  // different convention used in CSC StripDigiSimLink (full-strip)

  // Check all 6 CSC layers. Each layer has a distinct detId
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    CSCDetId detId0 = hit.CSC_DetId();
    CSCDetId detId1(detId0.endcap(), detId0.station(), detId0.ring(), detId0.chamber(), ilayer+1);

    StripDigiSimLinks::const_iterator cscStripLayerLinks = cscStripSimLinksPtr_->find(detId1);
    if (cscStripLayerLinks != cscStripSimLinksPtr_->end()) {
      for (LayerLinks::const_iterator linkItr = cscStripLayerLinks->begin(); linkItr != cscStripLayerLinks->end(); ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (std::abs(strip1 - static_cast<int>(channel)) <= 3) {  // allow +/-3
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
          matches[matchId] += fraction;
        }
      }
    }
  }
  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findCSCWireSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int wire0 = hit.Wire();
  int wire1 = (wire0 + 100) + 1;  // different convention used in CSC StripDigiSimLink

  // Check all 6 CSC layers. Each layer has a distinct detId
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    CSCDetId detId0 = hit.CSC_DetId();
    CSCDetId detId1(detId0.endcap(), detId0.station(), detId0.ring(), detId0.chamber(), ilayer+1);

    WireDigiSimLinks::const_iterator cscWireLayerLinks = cscWireSimLinksPtr_->find(detId1);
    if (cscWireLayerLinks != cscWireSimLinksPtr_->end()) {
      for (LayerLinks::const_iterator linkItr = cscWireLayerLinks->begin(); linkItr != cscWireLayerLinks->end(); ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (std::abs(wire1 - static_cast<int>(channel)) <= 3) {  // allow +/-3
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
          matches[matchId] += fraction;
        }
      }
    }
  }
  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findRPCDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  // Check all strips in the RPC cluster
  RPCDetId detId0 = hit.RPC_DetId();
  int stripA = hit.Strip_low();
  int stripB = hit.Strip_hi();
  int bx     = hit.BX();

#if 0
  // This is giving me an error: Assertion `std::distance(p.first, p.second) == 1' failed.
  RPCDigiSimLinks::const_iterator rpcDigiLayerLinks = rpcDigiSimLinksPtr_->find(detId0);
  if (rpcDigiLayerLinks != rpcDigiSimLinksPtr_->end()) {
    for (RPCLayerLinks::const_iterator linkItr = rpcDigiLayerLinks->begin(); linkItr != rpcDigiLayerLinks->end(); ++linkItr) {
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      for (int strip0 = stripA; strip0 < stripB+1; ++strip0) {
        if ((std::abs(strip0 - static_cast<int>(simStrip)) <= 1) && (bx == static_cast<int>(simBX))) {  // allow +/-1
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
          matches[matchId] += 1.0;
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

      if (detUnitId == detId0.rawId()) {
        for (int strip0 = stripA; strip0 < stripB+1; ++strip0) {
          if ((std::abs(strip0 - static_cast<int>(simStrip)) <= 1) && (bx == static_cast<int>(simBX))) {  // allow +/-1
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
            matches[matchId] += 1.0;

            //// Debug
            //std::cout << "Dump RPCDigiSimLink - strip " << linkItr->getStrip() << " bx " << linkItr->getBx()
            //          << " entry " << linkItr->getEntryPoint() << " p4 " << linkItr->getMomentumAtEntry()
            //          << " tof " << linkItr->getTimeOfFlight() << " eloss " << linkItr->getEnergyLoss()
            //          << " pdg " << linkItr->getParticleType() << " process " << linkItr->getProcessType()
            //          << " trackId " << linkItr->getTrackId() << " eventId " << linkItr->getEventId().bunchCrossing() << "," << linkItr->getEventId().event()
            //          << std::endl;
          }
        }
      }
    }
  }
#endif

  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findGEMDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  // Check all strips in the GEM cluster
  int stripA = hit.Strip_low();
  int stripB = hit.Strip_hi();
  int bx     = hit.BX();

  for (unsigned gemlayer = 0; gemlayer < 2; ++gemlayer) {  // check both GEM layers
    GEMDetId detId0 = hit.GEM_DetId();
    GEMDetId detId1(detId0.region(), detId0.ring(), detId0.station(), gemlayer+1, detId0.chamber(), detId0.roll());

    GEMDigiSimLinks::const_iterator gemDigiLayerLinks = gemDigiSimLinksPtr_->find(detId1);
    if (gemDigiLayerLinks != gemDigiSimLinksPtr_->end()) {
      for (GEMLayerLinks::const_iterator linkItr = gemDigiLayerLinks->begin(); linkItr != gemDigiLayerLinks->end(); ++linkItr) {
        unsigned int simStrip = linkItr->getStrip();
        unsigned int simBX = linkItr->getBx();
        unsigned int simTrackId = linkItr->getTrackId();
        EncodedEventId eventId = linkItr->getEventId();

        //int simPad = 1 + static_cast<int>(p->padOfStrip(simStrip));
        unsigned int simPad = (simStrip+1)/2;

        for (int strip0 = stripA; strip0 < stripB+1; ++strip0) {
          if ((detId1.station() == 1) && (std::abs(strip0 - static_cast<int>(simPad)) <= 2) && (bx == static_cast<int>(simBX))) {  // allow +/-2
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
            matches[matchId] += 1.0;

            //// Debug
            //std::cout << "Dump GEMDigiSimLink - strip " << linkItr->getStrip() << " bx " << linkItr->getBx()
            //          << " entry " << linkItr->getEntryPoint() << " p4 " << linkItr->getMomentumAtEntry()
            //          << " tof " << linkItr->getTimeOfFlight() << " eloss " << linkItr->getEnergyLoss()
            //          << " pdg " << linkItr->getParticleType() << " process " << linkItr->getProcessType()
            //          << " trackId " << linkItr->getTrackId() << " eventId " << linkItr->getEventId().bunchCrossing() << "," << linkItr->getEventId().event()
            //          << std::endl;
          } else if ((detId1.station() == 2) && (std::abs(strip0 - static_cast<int>(simPad)) <= 2) && (bx == static_cast<int>(simBX))) {  // allow +/-2
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
            matches[matchId] += 1.0;
          }
        }
      }
    }
  }
  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findME0DigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int strip0 = hit.Strip();  // 'phiposition' is in half-strip unit
  int strip1 = (strip0 >> 1);
  int bx     = hit.BX();

  // Check all 6 ME0 layers.
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    // Check neighbor rolls
    int partition = hit.Roll();  // 'partition' is in half-roll unit
    int iroll0 = (partition >> 1) + 1;
    unsigned iroll_first = (iroll0 == 1) ? iroll0 : iroll0-1;
    unsigned iroll_last = (iroll0 == 8) ? iroll0 : iroll0+1;
    for (unsigned iroll = iroll_first; iroll <= iroll_last; ++iroll) {
      ME0DetId detId0 = hit.ME0_DetId();
      ME0DetId detId1(detId0.region(), ilayer+1, detId0.chamber(), iroll);

      ME0DigiSimLinks::const_iterator me0DigiLayerLinks = me0DigiSimLinksPtr_->find(detId1);
      if (me0DigiLayerLinks != me0DigiSimLinksPtr_->end()) {
        for (ME0LayerLinks::const_iterator linkItr = me0DigiLayerLinks->begin(); linkItr != me0DigiLayerLinks->end(); ++linkItr) {
          unsigned int simStrip = linkItr->getStrip();
          unsigned int simBX = linkItr->getBx();
          unsigned int simTrackId = linkItr->getTrackId();
          EncodedEventId eventId = linkItr->getEventId();

          if ((std::abs(strip1 - static_cast<int>(simStrip)) <= 3) && (bx == static_cast<int>(simBX))) {  // allow +/-3
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
            matches[matchId] += 1.0;
          }
        }
      }
    }
  }
  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findDTDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  DTChamberId detId0 = hit.DT_DetId();
  int bx = hit.BX();

  for (DTDigiSimLinkCollection::DigiRangeIterator detUnit = dtDigiSimLinksPtr_->begin(); detUnit != dtDigiSimLinksPtr_->end(); ++detUnit) {
    const DTLayerId& layerid = (*detUnit).first;
    const DTDigiSimLinkCollection::Range& range = (*detUnit).second;

    if (detId0 == layerid.chamberId()) {

      for (DTDigiSimLinkCollection::const_iterator linkItr = range.first; linkItr != range.second; ++linkItr) {
        // Unfortunately lost all these info in the L1 data format L1MuDTChambPhDigi
        //float digitime = linkItr->time();
        //int wire = linkItr->wire();
        //int digiNumber = linkItr->number();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();

        if (bx == eventId.bunchCrossing()) {
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())  matches[matchId] = 0.;
          matches[matchId] += 1.0;
        }
      }
    }
  }
  return findTrackingParticle(matches, trkPartColl);
}

int EMTFMCTruth::findTrackingParticle(const std::map<SimHitIdpr, float>& matches, const TrackingParticleCollection& trkPartColl) const {
  // Find the matched TP with the highest weight
  int best_tp = -1;  // index of the tp
  float weight = 0.;
  for (std::map<SimHitIdpr, float>::const_iterator it_match = matches.begin(); it_match != matches.end(); ++it_match) {
    auto found = trackingParticleLinks_.find(it_match->first);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      auto it_trkpart = trkPartColl.begin();
      std::advance(it_trkpart, tp);
      if (weight < it_match->second) {
        weight = it_match->second;
        best_tp = tp;
      }
    }
  }
  return best_tp;
}

std::vector<EMTFMCTruth::EndcapSimHit> EMTFMCTruth::findEndcapSimHits() const {
  std::vector<EndcapSimHit> result;

  // CSC
  for (edm::PSimHitContainer::const_iterator it = cscSimHitsPtr_->begin(); it != cscSimHitsPtr_->end(); ++it) {
    unsigned int detUnitId = it->detUnitId();
    unsigned int simTrackId = it->trackId();
    EncodedEventId eventId = it->eventId();
    SimHitIdpr matchId(simTrackId, eventId);
    auto found = trackingParticleLinks_.find(matchId);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      const int subsystem = 1;
      const CSCLayer* layer = cscGeomPtr_->layer(detUnitId);
      const CSCDetId& id = layer->id();
      const LocalPoint& lp = it->localPosition();
      const GlobalPoint& gp = layer->toGlobal(lp);

      EndcapSimHit e_simhit = {subsystem, id.station(), id.ring(), id.layer(), id.chamber(), *it, gp, tp};
      result.emplace_back(e_simhit);
      //std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << e_simhit.type << " " << e_simhit.layer << " " << e_simhit.chamber << " " << e_simhit.phi << " " << e_simhit.theta << " " << e_simhit.eta << " " << e_simhit.r << " " << e_simhit.z << std::endl;
    }
  }

  // RPC
  for (edm::PSimHitContainer::const_iterator it = rpcSimHitsPtr_->begin(); it != rpcSimHitsPtr_->end(); ++it) {
    unsigned int detUnitId = it->detUnitId();
    unsigned int simTrackId = it->trackId();
    EncodedEventId eventId = it->eventId();
    SimHitIdpr matchId(simTrackId, eventId);
    auto found = trackingParticleLinks_.find(matchId);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      const int subsystem = 2;
      const RPCRoll* roll = rpcGeomPtr_->roll(detUnitId);
      const RPCDetId& id = roll->id();
      if (id.region() == 0)  continue;  // ignore barrel RPC
      const bool is_irpc = (id.station() == 3 || id.station() == 4) && (id.ring() == 1);
      const int id_chamber = is_irpc ? ((id.sector() - 1)*3 + id.subsector()) : ((id.sector() - 1)*6 + id.subsector());

      const LocalPoint& lp = it->localPosition();
      const GlobalPoint& gp = roll->toGlobal(lp);

      EndcapSimHit e_simhit = {subsystem, id.station(), id.ring(), id.layer(), id_chamber, *it, gp, tp};
      result.emplace_back(e_simhit);
      //std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << e_simhit.type << " " << e_simhit.layer << " " << e_simhit.chamber << " " << e_simhit.phi << " " << e_simhit.theta << " " << e_simhit.eta << " " << e_simhit.r << " " << e_simhit.z << std::endl;
    }
  }

  // GEM
  for (edm::PSimHitContainer::const_iterator it = gemSimHitsPtr_->begin(); it != gemSimHitsPtr_->end(); ++it) {
    unsigned int detUnitId = it->detUnitId();
    unsigned int simTrackId = it->trackId();
    EncodedEventId eventId = it->eventId();
    SimHitIdpr matchId(simTrackId, eventId);
    auto found = trackingParticleLinks_.find(matchId);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      const int subsystem = 3;
      const GEMEtaPartition* roll = gemGeomPtr_->etaPartition(detUnitId);
      const GEMDetId& id = roll->id();
      const LocalPoint& lp = it->localPosition();
      const GlobalPoint& gp = roll->toGlobal(lp);

      EndcapSimHit e_simhit = {subsystem, id.station(), id.ring(), id.layer(), id.chamber(), *it, gp, tp};
      result.emplace_back(e_simhit);
      //std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << e_simhit.type << " " << e_simhit.layer << " " << e_simhit.chamber << " " << e_simhit.phi << " " << e_simhit.theta << " " << e_simhit.eta << " " << e_simhit.r << " " << e_simhit.z << std::endl;
    }
  }

  // ME0
  for (edm::PSimHitContainer::const_iterator it = me0SimHitsPtr_->begin(); it != me0SimHitsPtr_->end(); ++it) {
    unsigned int detUnitId = it->detUnitId();
    unsigned int simTrackId = it->trackId();
    EncodedEventId eventId = it->eventId();
    SimHitIdpr matchId(simTrackId, eventId);
    auto found = trackingParticleLinks_.find(matchId);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      const int subsystem = 4;
      const ME0EtaPartition* roll = me0GeomPtr_->etaPartition(detUnitId);
      const ME0DetId& id = roll->id();
      const LocalPoint& lp = it->localPosition();
      const GlobalPoint& gp = roll->toGlobal(lp);

      EndcapSimHit e_simhit = {subsystem, id.station(), 1, id.layer(), id.chamber(), *it, gp, tp};
      result.emplace_back(e_simhit);
      //std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << e_simhit.type << " " << e_simhit.layer << " " << e_simhit.chamber << " " << e_simhit.phi << " " << e_simhit.theta << " " << e_simhit.eta << " " << e_simhit.r << " " << e_simhit.z << std::endl;
    }
  }

  // DT
  for (edm::PSimHitContainer::const_iterator it = dtSimHitsPtr_->begin(); it != dtSimHitsPtr_->end(); ++it) {
    unsigned int detUnitId = it->detUnitId();
    unsigned int simTrackId = it->trackId();
    EncodedEventId eventId = it->eventId();
    SimHitIdpr matchId(simTrackId, eventId);
    auto found = trackingParticleLinks_.find(matchId);
    if (found != trackingParticleLinks_.end()) {
      int tp = found->second;
      const int subsystem = 0;
      const DTWireId wireid(detUnitId);
      const DTLayer* layer = dtGeomPtr_->layer(wireid.layerId());
      const DTLayerId& id = layer->id();
      const LocalPoint& lp = it->localPosition();
      const GlobalPoint& gp = layer->toGlobal(lp);

      EndcapSimHit e_simhit = {subsystem, id.station(), 1, id.layer(), id.sector(), *it, gp, tp};
      result.emplace_back(e_simhit);
      //std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << e_simhit.type << " " << e_simhit.layer << " " << e_simhit.chamber << " " << e_simhit.phi << " " << e_simhit.theta << " " << e_simhit.eta << " " << e_simhit.r << " " << e_simhit.z << std::endl;
    }
  }

  return result;
}
