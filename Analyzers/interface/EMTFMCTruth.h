#ifndef L1TMuonSimulations_EMTFMCTruth_h
#define L1TMuonSimulations_EMTFMCTruth_h

//
// References:
//   SimMuon/MCTruth/interface/MuonTruth.h
//   SimMuon/MCTruth/interface/CSCHitAssociator.h
//   SimMuon/MCTruth/interface/RPCHitAssociator.h
//   SimMuon/MCTruth/interface/GEMHitAssociator.h
//

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"

//#include "SimDataFormats/Track/interface/SimTrack.h"
//#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

// Forwards
namespace l1t {
  class EMTFHit;
}


// Class definition
class EMTFMCTruth
{
public:

  // Typedefs
  // From SimMuon/MCTruth/interface/MuonTruth.h
  typedef edm::DetSetVector<StripDigiSimLink> StripDigiSimLinks;
  typedef edm::DetSetVector<StripDigiSimLink> WireDigiSimLinks;
  typedef edm::DetSet<StripDigiSimLink>       LayerLinks;
  typedef edm::DetSetVector<RPCDigiSimLink>   RPCDigiSimLinks;
  typedef edm::DetSet<RPCDigiSimLink>         RPCLayerLinks;
  typedef edm::DetSetVector<GEMDigiSimLink>   GEMDigiSimLinks;
  typedef edm::DetSet<GEMDigiSimLink>         GEMLayerLinks;
  typedef std::pair<uint32_t, EncodedEventId> SimHitIdpr;

  // Constructor
  EMTFMCTruth(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iConsumes);

  // Destructor
  ~EMTFMCTruth();

  // Member functions

  // Set up Handles/ESHandles
  void initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  void makeTrackingParticleLinks(const TrackingParticleCollection& trkPartColl);

  // Find the matching tracking particle (tp); returns the index of the tp
  int findCSCStripSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
  int findCSCWireSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
  int findRPCDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
  int findGEMDigiSimLink(const l1t::EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;


private:
  int findTrackingParticle(const std::map<SimHitIdpr, int>& matches, const TrackingParticleCollection& trkPartColl) const;

  const edm::InputTag   cscSimHitsTag_, cscSimHitsXFTag_, cscStripSimLinksTag_, cscWireSimLinksTag_;
  const edm::InputTag   rpcSimHitsTag_, rpcSimHitsXFTag_, rpcDigiSimLinksTag_;
  const edm::InputTag   gemSimHitsTag_, gemSimHitsXFTag_, gemDigiSimLinksTag_;

  edm::EDGetTokenT<edm::PSimHitContainer>       cscSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit> >     cscSimHitsXFToken_;
  edm::EDGetTokenT<StripDigiSimLinks>           cscStripSimLinksToken_;
  edm::EDGetTokenT<WireDigiSimLinks>            cscWireSimLinksToken_;
  edm::EDGetTokenT<edm::PSimHitContainer>       rpcSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit> >     rpcSimHitsXFToken_;
  edm::EDGetTokenT<RPCDigiSimLinks>             rpcDigiSimLinksToken_;
  edm::EDGetTokenT<edm::PSimHitContainer>       gemSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit> >     gemSimHitsXFToken_;
  edm::EDGetTokenT<GEMDigiSimLinks>             gemDigiSimLinksToken_;

  const edm::PSimHitContainer *     cscSimHitsPtr_;
  const CrossingFrame<PSimHit> *    cscSimHitsXFPtr_;
  const StripDigiSimLinks*          cscStripSimLinksPtr_;
  const WireDigiSimLinks*           cscWireSimLinksPtr_;
  const edm::PSimHitContainer *     rpcSimHitsPtr_;
  const CrossingFrame<PSimHit> *    rpcSimHitsXFPtr_;
  const RPCDigiSimLinks*            rpcDigiSimLinksPtr_;
  const edm::PSimHitContainer *     gemSimHitsPtr_;
  const CrossingFrame<PSimHit> *    gemSimHitsXFPtr_;
  const GEMDigiSimLinks*            gemDigiSimLinksPtr_;

  std::map<SimHitIdpr, int> trackingParticleLinks_;
};

#endif
