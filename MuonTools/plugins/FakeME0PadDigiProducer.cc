//
// Adapted from SimMuon/GEMDigitizer/interface/ME0PadDigiProducer.h
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"
#include "DataFormats/GEMDigi/interface/ME0PadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

class ME0Geometry;

class FakeME0PadDigiProducer : public edm::stream::EDProducer<>
{
public:

  explicit FakeME0PadDigiProducer(const edm::ParameterSet&);

  virtual ~FakeME0PadDigiProducer();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;

  virtual void produce(edm::Event&, const edm::EventSetup&) override;

private:

  void buildPads(const ME0DigiCollection &digis, ME0PadDigiCollection &out_pads) const;

  void buildPads(const ME0DigiPreRecoCollection &digis, ME0PadDigiCollection &out_pads) const;

  void buildPads(const ME0SegmentCollection &segments, ME0PadDigiCollection &out_pads) const;

  ME0DetId makeFullDetId(const ME0DetId& detId, const ME0Segment& segment) const {
    int layer = 3;
    int roll = std::round(3.757 + (-0.07915) * segment.localPosition().y());  // from my own empirical fit
    ME0DetId fullDetId(detId.region(), layer, detId.chamber(), roll);
    return fullDetId;
  }

  ME0DetId makeFullDetId(const ME0DetId& detId, int roll) const {
    int layer = 3;
    ME0DetId fullDetId(detId.region(), layer, detId.chamber(), roll);
    return fullDetId;
  }

  /// Name of input digi collection
  edm::EDGetTokenT<ME0DigiPreRecoCollection> digis_token_;
  edm::InputTag digis_;

  /// Name of input segment collection
  edm::EDGetTokenT<ME0SegmentCollection> segments_token_;
  edm::InputTag segments_;

  const ME0Geometry * geometry_;
};


// _____________________________________________________________________________
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include <set>


FakeME0PadDigiProducer::FakeME0PadDigiProducer(const edm::ParameterSet& ps)
: geometry_(nullptr)
{
  //digis_ = ps.getParameter<edm::InputTag>("InputCollection");
  digis_ = edm::InputTag("simMuonME0ReDigis");  // hardcoded for now

  digis_token_ = consumes<ME0DigiPreRecoCollection>(digis_);

  segments_ = ps.getParameter<edm::InputTag>("InputCollection");

  segments_token_ = consumes<ME0SegmentCollection>(segments_);

  produces<ME0PadDigiCollection>();
}


FakeME0PadDigiProducer::~FakeME0PadDigiProducer()
{}


void FakeME0PadDigiProducer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  eventSetup.get<MuonGeometryRecord>().get(hGeom);
  geometry_ = &*hGeom;
}


void FakeME0PadDigiProducer::produce(edm::Event& e, const edm::EventSetup& eventSetup)
{
  //edm::Handle<ME0DigiPreRecoCollection> hdigis;
  //e.getByToken(digis_token_, hdigis);

  edm::Handle<ME0SegmentCollection> hsegments;
  e.getByToken(segments_token_, hsegments);

  // Create empty output
  std::unique_ptr<ME0PadDigiCollection> pPads(new ME0PadDigiCollection());

  // build the pads
  buildPads(*(hsegments.product()), *pPads);

  // store them in the event
  e.put(std::move(pPads));
}

void FakeME0PadDigiProducer::buildPads(const ME0DigiCollection &digis, ME0PadDigiCollection &out_pads) const
{
  auto etaPartitions = geometry_->etaPartitions();
  for(const auto& p: etaPartitions)
  {
    // set of <pad, bx> pairs, sorted first by pad then by bx
    std::set<std::pair<int, int> > proto_pads;

    // walk over digis in this partition,
    // and stuff them into a set of unique pads (equivalent of OR operation)
    auto id_digis = digis.get(p->id());
    for (auto d = id_digis.first; d != id_digis.second; ++d)
    {
      int pad_num = 1 + static_cast<int>( p->padOfStrip(d->strip()) );
      proto_pads.emplace(pad_num, d->bx());
    }

    // fill the output collections
    for (const auto & d: proto_pads)
    {
      ME0PadDigi pad_digi(d.first, d.second);
      out_pads.insertDigi(p->id(), pad_digi);
    }
  }
}

void FakeME0PadDigiProducer::buildPads(const ME0DigiPreRecoCollection &digis, ME0PadDigiCollection &out_pads) const
{
  for (auto det = digis.begin(); det != digis.end(); ++det)
  {
    //auto detId = (*det).first;
    auto digi = (*det).second.first;
    auto dend = (*det).second.second;
    for (auto d = digi; d != dend; ++d)
    {
      //std::cout << "detId: " << detId << " x: " << d->x() << " y: " << d->y() << " tof: " << d->tof() << std::endl;
    }
  }
}

void FakeME0PadDigiProducer::buildPads(const ME0SegmentCollection &segments, ME0PadDigiCollection &out_pads) const
{
  for (auto id_iter = segments.id_begin(); id_iter != segments.id_end(); ++id_iter)
  {
    //// set of <pad, bx> pairs, sorted first by pad then by bx
    //std::set<std::pair<int, int> > proto_pads;

    // set of <<pad, partition>, bx> pairs, sorted first by pad, then by bx and then by roll
    std::set<std::pair<std::pair<int, int>, int> > proto_pads;

    // walk over segments in this partition,
    // and stuff them into a set of unique pads (equivalent of OR operation)
    const ME0DetId detId(*id_iter);
    auto id_segments = segments.get(detId);
    for (auto seg = id_segments.first; seg != id_segments.second; ++seg)
    {
      const ME0DetId fullDetId = makeFullDetId(detId, *seg);  // the original detId has layer 0 and roll 0
      const ME0EtaPartition* p = geometry_->etaPartition(fullDetId);

      //int strip_num = p->strip(seg->localPosition());
      //int pad_num = 1 + static_cast<int>( p->padOfStrip(strip_num) );
      LocalPoint lp(seg->localPosition().x(), 0., 0.);
      int channel = p->specificTopology().channel(lp);
      channel >>= 1;  // make 2-strip pads
      int pad_num = 1 + channel;
      const float bxWidth = 25.0;
      int bxIdx = seg->time() / bxWidth;
      proto_pads.emplace(std::make_pair(std::make_pair(pad_num, bxIdx), fullDetId.roll()) );
    }

    // fill the output collections
    for (const auto & d: proto_pads)
    {
      const ME0DetId fullDetId = makeFullDetId(detId, d.second);
      ME0PadDigi pad_digi(d.first.first, d.first.second);
      out_pads.insertDigi(fullDetId, pad_digi);
    }
  }
}


// _____________________________________________________________________________
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FakeME0PadDigiProducer);
