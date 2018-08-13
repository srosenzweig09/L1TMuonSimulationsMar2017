#ifndef L1TMuonSimulations_DummyTTTrackProducer_h
#define L1TMuonSimulations_DummyTTTrackProducer_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"


class DummyTTTrackProducer : public edm::EDProducer
{

public:
  /// Constructor
  explicit DummyTTTrackProducer(const edm::ParameterSet& iConfig);

  /// Destructor
  ~DummyTTTrackProducer();

private:
  //virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  //virtual void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
};

#endif
