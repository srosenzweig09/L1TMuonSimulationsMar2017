#include "L1TMuonSimulations/TrackerTools/plugins/DummyTTTrackProducer.h"

DummyTTTrackProducer::DummyTTTrackProducer(const edm::ParameterSet& iConfig)
{
  produces<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >("Level1TTTracks");
}

DummyTTTrackProducer::~DummyTTTrackProducer()
{

}

void DummyTTTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto ttTracksForOutput = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_> >>();

  iEvent.put(std::move(ttTracksForOutput), "Level1TTTracks");
}

// Define this as a plug-in
DEFINE_FWK_MODULE(DummyTTTrackProducer);
