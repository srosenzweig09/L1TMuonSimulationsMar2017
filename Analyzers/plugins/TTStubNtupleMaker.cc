#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <cassert>
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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


// _____________________________________________________________________________
class TTStubNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
public:
  explicit TTStubNtupleMaker(const edm::ParameterSet& iConfig);
  ~TTStubNtupleMaker();

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

  // Typedefs
  typedef Phase2TrackerDigi                     digi_t;
  typedef edm::DetSetVector<digi_t>             digi_dsv_t;
  typedef edm::DetSet<digi_t>                   digi_ds_t;
  typedef Ref_Phase2TrackerDigi_                digi_ref_t;

  typedef PixelDigiSimLink                      digilink_t;
  typedef edm::DetSetVector<digilink_t>         digilink_dsv_t;
  typedef edm::DetSet<digilink_t>               digilink_ds_t;

  typedef TTCluster<digi_ref_t>                 ttclus_t;
  typedef edmNew::DetSetVector<ttclus_t>        ttclus_dsv_t;
  typedef edmNew::DetSet<ttclus_t>              ttclus_ds_t;
  typedef edm::Ref<ttclus_dsv_t, ttclus_t>      ttclus_ref_t;
  typedef TTClusterAssociationMap<digi_ref_t>   ttclus_assocmap_t;
  typedef TTStub<digi_ref_t>                    ttstub_t;
  typedef edmNew::DetSetVector<ttstub_t>        ttstub_dsv_t;
  typedef edmNew::DetSet<ttstub_t>              ttstub_ds_t;
  typedef edm::Ref<ttstub_dsv_t, ttstub_t>      ttstub_ref_t;
  typedef TTStubAssociationMap<digi_ref_t>      ttstub_assocmap_t;

  template<typename T>
  edm::Handle<T> make_handle(T t)
  {
    return edm::Handle<T>();
  }

  // Input parameters
  const edm::InputTag   ttstubTag_;
  const edm::InputTag   ttstubAssocTag_;
  const edm::InputTag   genPartTag_;
  const edm::InputTag   trkPartTag_;
  const std::string     outFileName_;
  int verbose_;

  // Member data
  bool firstEvent_;
  edm::EDGetTokenT<ttstub_dsv_t>                  ttstubToken_;
  edm::EDGetTokenT<ttstub_assocmap_t>             ttstubAssocToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>   genPartToken_;
  edm::EDGetTokenT<TrackingParticleCollection>    trkPartToken_;

  reco::GenParticleCollection   genParts_;
  TrackingParticleCollection    trkParts_;

  // For edm products
  const ttstub_dsv_t *      ttstubs_;
  const ttstub_assocmap_t * ttstubAssoc_;

  // For edm event setup
  const TrackerGeometry * theGeometry;
  const TrackerTopology * theTopology;
  const MagneticField* theMagneticField;
  std::map<uint32_t, uint32_t> stackIdToGeoIdMap;

  // TTree
  TTree* tree;

  // Hits
  std::unique_ptr<std::vector<int16_t> >  vh_endcap;
  std::unique_ptr<std::vector<int16_t> >  vh_station;
  std::unique_ptr<std::vector<int16_t> >  vh_ring;
  std::unique_ptr<std::vector<int16_t> >  vh_module;
  std::unique_ptr<std::vector<int16_t> >  vh_psmodule;
  std::unique_ptr<std::vector<int16_t> >  vh_bx;
  //
  std::unique_ptr<std::vector<float> >    vh_coordx;  // float, in full-strip unit
  std::unique_ptr<std::vector<float> >    vh_coordy;  // float
  std::unique_ptr<std::vector<float> >    vh_bend;    // float
  std::unique_ptr<std::vector<int16_t> >  vh_fr;
  //
  std::unique_ptr<std::vector<float  > >  vh_sim_phi;
  std::unique_ptr<std::vector<float  > >  vh_sim_theta;
  std::unique_ptr<std::vector<float  > >  vh_sim_eta;
  std::unique_ptr<std::vector<float  > >  vh_sim_r;
  std::unique_ptr<std::vector<float  > >  vh_sim_z;
  std::unique_ptr<std::vector<int32_t> >  vh_sim_tp;
  std::unique_ptr<std::vector<int32_t> >  vh_sim_pdgid;
  std::unique_ptr<std::vector<int16_t> >  vh_sim_assoc;  // isGenuine, isLooselyGenuine, isCombinatoric, isUnknown
  //
  std::unique_ptr<int32_t              >  vh_size;
};


// _____________________________________________________________________________
TTStubNtupleMaker::TTStubNtupleMaker(const edm::ParameterSet& iConfig) :
    ttstubTag_      (iConfig.getParameter<edm::InputTag>("ttstubTag")),
    ttstubAssocTag_ (iConfig.getParameter<edm::InputTag>("ttstubAssocTag")),
    genPartTag_     (iConfig.getParameter<edm::InputTag>("genPartTag")),
    trkPartTag_     (iConfig.getParameter<edm::InputTag>("trkPartTag")),
    outFileName_    (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_        (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");  // shared resources

  firstEvent_ = true;

  ttstubToken_       = consumes<ttstub_dsv_t>                   (ttstubTag_);
  ttstubAssocToken_  = consumes<ttstub_assocmap_t>              (ttstubAssocTag_);
  genPartToken_      = consumes<reco::GenParticleCollection>    (genPartTag_);
  trkPartToken_      = consumes<TrackingParticleCollection>     (trkPartTag_);
}

TTStubNtupleMaker::~TTStubNtupleMaker() {}

void TTStubNtupleMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  /// Geometry setup
  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
  theGeometry = geometryHandle.product();

  edm::ESHandle<TrackerTopology> topologyHandle;
  iSetup.get<TrackerTopologyRcd>().get(topologyHandle);
  theTopology = topologyHandle.product();

  /// Magnetic field setup
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  theMagneticField = magneticFieldHandle.product();
}

void TTStubNtupleMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void TTStubNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent, iSetup);
  process(iEvent, iSetup);
}

// _____________________________________________________________________________
void TTStubNtupleMaker::getHandles(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // TTStub
  edm::Handle<ttstub_dsv_t> ttstubs_handle;
  edm::Handle<ttstub_assocmap_t> ttstubAssoc_handle;

  if (!ttstubToken_.isUninitialized()) {
    iEvent.getByToken(ttstubToken_, ttstubs_handle);
  }
  if (!ttstubs_handle.isValid()) {
    if (firstEvent_)  edm::LogError("TTStubNtupleMaker") << "Cannot get the product: " << ttstubTag_;
    ttstubs_ = nullptr;
  } else {
    ttstubs_ = ttstubs_handle.product();
  }

  if (!ttstubAssocToken_.isUninitialized()) {
    iEvent.getByToken(ttstubAssocToken_, ttstubAssoc_handle);
  }
  if (!ttstubAssoc_handle.isValid()) {
    if (firstEvent_)  edm::LogError("TTStubNtupleMaker") << "Cannot get the product: " << ttstubAssocTag_;
    ttstubAssoc_ = nullptr;
  } else {
    ttstubAssoc_ = ttstubAssoc_handle.product();
  }

  // Gen particles
  auto genParts_handle = make_handle(genParts_);

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      if (firstEvent_)  edm::LogError("NtupleMaker") << "Cannot get the product: " << genPartTag_;
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

  // Map stackDetId to geographicalId
  if (firstEvent_) {
    for (auto const & det_u : theGeometry->detUnits()) {
      DetId detId = det_u->geographicalId();
      if (detId.subdetId()!=StripSubdetector::TOB && detId.subdetId()!=StripSubdetector::TID)  // only run on outer tracker
        continue;
      if (!theTopology->isLower(detId))  // loop on the stacks: choose the lower arbitrarily
        continue;
      DetId stackDetId = theTopology->stack(detId);
      stackIdToGeoIdMap[stackDetId.rawId()] = detId.rawId();
      //std::cout << theTopology->print(detId) << std::endl;

      using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
      const Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u);
      assert(pixdet);
    }
    edm::LogInfo("TTStubNtupleMaker") << "Found " << stackIdToGeoIdMap.size() << " stackDetIds.";
  }
}

// _____________________________________________________________________________
void TTStubNtupleMaker::process(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // In-place functions

  auto rad_to_deg = [](double rad) {
    constexpr double factor = 180./M_PI;
    return rad * factor;
  };

  auto isBarrel = [](const auto& detid) {
    bool isBarrel = (detid.subdetId() == StripSubdetector::TOB);
    //bool isEndcap = (detid.subdetId() == StripSubdetector::TID);
    return isBarrel;
  };

  auto isPSModule = [&](const auto& detid) {
    const TrackerGeometry::ModuleType moduleType = theGeometry->getDetectorType(detid);
    bool isPSModule = (moduleType == TrackerGeometry::ModuleType::Ph2PSP) || (moduleType == TrackerGeometry::ModuleType::Ph2PSS);
    //bool isSSModule = (moduleType == TrackerGeometry::ModuleType::Ph2SS);
    return isPSModule;
  };

  auto get_region = [&](const auto& detid) {
    int region = 0;

    if (detid.subdetId() == StripSubdetector::TOB) {  // barrel
      region = 0;
    } else if (detid.subdetId() == StripSubdetector::TID) {  // endcap
      int type = theTopology->tidSide(detid);   // 1=-ve 2=+ve
      if (type == 1) {
        region = -1;
      } else if (type == 2) {
        region = +1;
      }
    }
    return region;
  };

  auto get_layer = [&](const auto& detid) {
    int layer = 0;

    if (detid.subdetId() == StripSubdetector::TOB) {  // barrel
      layer = static_cast<int>(theTopology->layer(detid));
    } else if (detid.subdetId() == StripSubdetector::TID) {  // endcap
      layer = static_cast<int>(theTopology->layer(detid));
    }
    return layer;
  };

  auto get_ring = [&](const auto& detid) {
    int ring = 0;

    if (detid.subdetId() == StripSubdetector::TOB) {  // barrel
      ring = static_cast<int>(theTopology->tobRod(detid));
    } else if (detid.subdetId() == StripSubdetector::TID) {  // endcap
      ring = static_cast<int>(theTopology->tidRing(detid));
    }
    return ring;
  };

  auto get_module = [&](const auto& detid) {
    int module = 0;

    if (detid.subdetId() == StripSubdetector::TOB) {  // barrel
      module = static_cast<int>(theTopology->module(detid));
    } else if (detid.subdetId() == StripSubdetector::TID) {  // endcap
      module = static_cast<int>(theTopology->module(detid));
    }
    return module;
  };

  // ___________________________________________________________________________
  if (firstEvent_)  edm::LogInfo("TTStubNtupleMaker") << "Ready to make ntuple.";

  // TTStub
  using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
  using Phase2TrackerTopology    = PixelTopology;

  if (ttstubs_) {
    edm::Handle<ttstub_dsv_t> ttstubs_handle;
    iEvent.getByToken(ttstubToken_, ttstubs_handle);

    for (ttstub_dsv_t::const_iterator itv = ttstubs_->begin(); itv != ttstubs_->end(); ++itv) {

      for (ttstub_ds_t::const_iterator it = itv->begin(); it != itv->end(); ++it) {
        const ttstub_t ttstub = (*it);
        const DetId stackDetId = ttstub.getDetId();
        const DetId geoId0     = stackIdToGeoIdMap.at(stackDetId.rawId());
        const DetId geoId1     = theTopology->partnerDetId(geoId0);
        if (geoId1) {}  // don't warn

        bool is_barrel   = isBarrel(geoId0);
        bool is_psmodule = isPSModule(geoId0);
        if (is_barrel) {}  // don't warn

        int region    = get_region(geoId0);     // 0 for Barrel, +/-1 for +/- Endcap
        int endcap    = (region == -1) ? 2 : region;
        int station   = get_layer(geoId0);
        int ring      = get_ring(geoId0);
        int module    = get_module(geoId0);

        // Check L1Trigger/TrackTrigger/src/TTStubAlgorithm_official.cc
        const GeomDetUnit* geoUnit = theGeometry->idToDetUnit(geoId0);
        //const GeomDetUnit* geoUnit = theGeometry->idToDetUnit(geoId1);
        const Phase2TrackerGeomDetUnit* pixUnit = dynamic_cast<const Phase2TrackerGeomDetUnit*>(geoUnit);
        const Phase2TrackerTopology* pixTopo = dynamic_cast<const Phase2TrackerTopology*>(&(pixUnit->specificTopology()));

        const MeasurementPoint& localCoord = ttstub.getClusterRef(0)->findAverageLocalCoordinatesCentered();
        //const MeasurementPoint& localCoord = aTTStub.getClusterRef(0)->findAverageLocalCoordinates();
        const GlobalPoint globalPosition = pixUnit->surface().toGlobal(pixTopo->localPosition(localCoord));

        double bend   = ttstub.getTriggerBend();  // in full-strip unit
        //double bend   = ttstub.getHardwareBend();  // in full-strip unit

        int sim_tp    = -1;
        int sim_pdgid = 0;
        int sim_assoc = 0;
        {
          const ttstub_ref_t ref = edmNew::makeRefTo(ttstubs_handle, it);
          const edm::Ptr<TrackingParticle> tpPtr = ttstubAssoc_->findTrackingParticlePtr(ref);
          if (tpPtr.isNonnull()) {
            sim_tp    = tpPtr.key();
            sim_pdgid = tpPtr->pdgId();
          }
          if (ttstubAssoc_->isGenuine(ref)) {
            sim_assoc = 0;
          //} else if (ttstubAssoc_->isLooselyGenuine(ref)) {
          //  sim_assoc = 1;
          } else if (ttstubAssoc_->isCombinatoric(ref)) {
            sim_assoc = 2;
          } else if (ttstubAssoc_->isUnknown(ref)) {
            sim_assoc = 3;
          } else {
            sim_assoc = 4;
          }
        }


        vh_endcap  ->push_back(endcap);
        vh_station ->push_back(station);
        vh_ring    ->push_back(ring);
        vh_module  ->push_back(module);
        vh_psmodule->push_back(is_psmodule);
        vh_bx      ->push_back(0);

        vh_coordx   ->push_back(localCoord.x());
        vh_coordy   ->push_back(localCoord.y());
        vh_bend     ->push_back(bend);
        vh_fr       ->push_back(0);  //FIXME

        vh_sim_phi  ->push_back(rad_to_deg(globalPosition.phi()));
        vh_sim_theta->push_back(globalPosition.theta());
        vh_sim_eta  ->push_back(globalPosition.eta());
        vh_sim_r    ->push_back(globalPosition.perp());
        vh_sim_z    ->push_back(globalPosition.z());
        vh_sim_tp   ->push_back(sim_tp);
        vh_sim_pdgid->push_back(sim_pdgid);
        vh_sim_assoc->push_back(sim_assoc);
      }
    }
  }
  (*vh_size) = vh_endcap->size();

  // Fill
  tree->Fill();

  // Clear
  vh_endcap     ->clear();
  vh_station    ->clear();
  vh_ring       ->clear();
  vh_module     ->clear();
  vh_psmodule   ->clear();
  vh_bx         ->clear();
  //
  vh_coordx     ->clear();
  vh_coordy     ->clear();
  vh_bend       ->clear();
  vh_fr         ->clear();
  //
  vh_sim_phi    ->clear();
  vh_sim_theta  ->clear();
  vh_sim_eta    ->clear();
  vh_sim_r      ->clear();
  vh_sim_z      ->clear();
  vh_sim_tp     ->clear();
  vh_sim_pdgid  ->clear();
  vh_sim_assoc  ->clear();
  //
  (*vh_size)    = 0;

  if (firstEvent_)
    firstEvent_ = false;
}

// _____________________________________________________________________________
void TTStubNtupleMaker::beginJob() {
  makeTree();
}

void TTStubNtupleMaker::endJob() {
  writeTree();
}

// _____________________________________________________________________________
void TTStubNtupleMaker::makeTree() {

  // TFileService
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  // TTStub
  vh_endcap     = std::make_unique<std::vector<int16_t > >();
  vh_station    = std::make_unique<std::vector<int16_t > >();
  vh_ring       = std::make_unique<std::vector<int16_t > >();
  vh_module     = std::make_unique<std::vector<int16_t > >();
  vh_psmodule   = std::make_unique<std::vector<int16_t > >();
  vh_bx         = std::make_unique<std::vector<int16_t > >();
  //
  vh_coordx     = std::make_unique<std::vector<float   > >();
  vh_coordy     = std::make_unique<std::vector<float   > >();
  vh_bend       = std::make_unique<std::vector<float   > >();
  vh_fr         = std::make_unique<std::vector<int16_t > >();
  //
  vh_sim_phi    = std::make_unique<std::vector<float   > >();
  vh_sim_theta  = std::make_unique<std::vector<float   > >();
  vh_sim_eta    = std::make_unique<std::vector<float   > >();
  vh_sim_r      = std::make_unique<std::vector<float   > >();
  vh_sim_z      = std::make_unique<std::vector<float   > >();
  vh_sim_tp     = std::make_unique<std::vector<int32_t > >();
  vh_sim_pdgid  = std::make_unique<std::vector<int32_t > >();
  vh_sim_assoc  = std::make_unique<std::vector<int16_t > >();
  //
  vh_size       = std::make_unique<int32_t>(0);

  // Set branches
  // TTStub
  tree->Branch("vh_endcap"    , &(*vh_endcap    ));
  tree->Branch("vh_station"   , &(*vh_station   ));
  tree->Branch("vh_ring"      , &(*vh_ring      ));
  tree->Branch("vh_module"    , &(*vh_module    ));
  tree->Branch("vh_psmodule"  , &(*vh_psmodule  ));
  tree->Branch("vh_bx"        , &(*vh_bx        ));
  //
  tree->Branch("vh_coordx"    , &(*vh_coordx    ));
  tree->Branch("vh_coordy"    , &(*vh_coordy    ));
  tree->Branch("vh_bend"      , &(*vh_bend      ));
  tree->Branch("vh_fr"        , &(*vh_fr        ));
  //
  tree->Branch("vh_sim_phi"   , &(*vh_sim_phi   ));
  tree->Branch("vh_sim_theta" , &(*vh_sim_theta ));
  tree->Branch("vh_sim_eta"   , &(*vh_sim_eta   ));
  tree->Branch("vh_sim_r"     , &(*vh_sim_r     ));
  tree->Branch("vh_sim_z"     , &(*vh_sim_z     ));
  tree->Branch("vh_sim_tp"    , &(*vh_sim_tp    ));
  tree->Branch("vh_sim_pdgid" , &(*vh_sim_pdgid ));
  tree->Branch("vh_sim_assoc" , &(*vh_sim_assoc ));
  //
  tree->Branch("vh_size"      , &(*vh_size      ));
}

void TTStubNtupleMaker::writeTree() {
  // Handled by TFileService
}

// _____________________________________________________________________________
void TTStubNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(TTStubNtupleMaker);
