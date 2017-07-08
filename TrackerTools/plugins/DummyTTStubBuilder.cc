#include "L1TMuonSimulations/TrackerTools/plugins/DummyTTStubBuilder.h"

DummyTTStubBuilder::DummyTTStubBuilder(const edm::ParameterSet& iConfig)
{
  produces<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_> > >("ClusterAccepted");
  produces<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >("StubAccepted");
  produces<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >("StubRejected");
}

DummyTTStubBuilder::~DummyTTStubBuilder()
{

}

void DummyTTStubBuilder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto ttClusterDSVForOutput      = std::make_unique<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_> > >();
  auto ttStubDSVForOutputAccepted = std::make_unique<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >();
  auto ttStubDSVForOutputRejected = std::make_unique<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >();

  iEvent.put(std::move(ttClusterDSVForOutput), "ClusterAccepted");
  iEvent.put(std::move(ttStubDSVForOutputAccepted), "StubAccepted");
  iEvent.put(std::move(ttStubDSVForOutputRejected), "StubRejected");
}

// Define this as a plug-in
DEFINE_FWK_MODULE(DummyTTStubBuilder);
