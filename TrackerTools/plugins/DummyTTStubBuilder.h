#ifndef L1TMuonSimulations_DummyTTStubBuilder_h
#define L1TMuonSimulations_DummyTTStubBuilder_h

//
// Adapted from L1Trigger/TrackTrigger/plugins/TTStubBuilder.h
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"


class DummyTTStubBuilder : public edm::EDProducer
{

public:
  /// Constructor
  explicit DummyTTStubBuilder(const edm::ParameterSet& iConfig);

  /// Destructor
  ~DummyTTStubBuilder();

private:
  //virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  //virtual void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
};

#endif
