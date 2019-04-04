#include "L1TMuonSimulations/ParticleGuns/interface/FlatEvtVtxGenerator2.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
//#include "CLHEP/Vector/ThreeVector.h"
#include "HepMC/SimpleVector.h"

FlatEvtVtxGenerator2::FlatEvtVtxGenerator2(const edm::ParameterSet& p )
: BaseEvtVtxGenerator2(p)
{
  sourceToken=consumes<edm::HepMCProduct>(p.getParameter<edm::InputTag>("src"));

  fMinX = p.getParameter<double>("MinX")*cm;
  fMinY = p.getParameter<double>("MinY")*cm;
  fMinZ = p.getParameter<double>("MinZ")*cm;
  fMaxX = p.getParameter<double>("MaxX")*cm;
  fMaxY = p.getParameter<double>("MaxY")*cm;
  fMaxZ = p.getParameter<double>("MaxZ")*cm;
  fMinT = p.getParameter<double>("MinT")*ns*c_light;
  fMaxT = p.getParameter<double>("MaxT")*ns*c_light;

  fVtxSpectrum = p.exists("VtxSpectrum") ? p.getParameter<std::string>("VtxSpectrum") : "flatPhi";

  if (fMinX > fMaxX) {
    throw cms::Exception("Configuration")
      << "Error in FlatEvtVtxGenerator2: "
      << "MinX is greater than MaxX";
  }
  if (fMinY > fMaxY) {
    throw cms::Exception("Configuration")
      << "Error in FlatEvtVtxGenerator2: "
      << "MinY is greater than MaxY";
  }
  if (fMinZ > fMaxZ) {
    throw cms::Exception("Configuration")
      << "Error in FlatEvtVtxGenerator2: "
      << "MinZ is greater than MaxZ";
  }
  if (fMinT > fMaxT) {
    throw cms::Exception("Configuration")
      << "Error in FlatEvtVtxGenerator2: "
      << "MinT is greater than MaxT";
  }
}

FlatEvtVtxGenerator2::~FlatEvtVtxGenerator2()
{
}

//Hep3Vector * FlatEvtVtxGenerator2::newVertex() {
HepMC::FourVector FlatEvtVtxGenerator2::newVertex(CLHEP::HepRandomEngine* engine) const {
  double aRho,aPhi,aX,aY,aZ,aT;
  //aX = CLHEP::RandFlat::shoot(engine, fMinX, fMaxX);
  //aY = CLHEP::RandFlat::shoot(engine, fMinY, fMaxY);
  aRho = CLHEP::RandFlat::shoot(engine, 0, fMaxX);
  aPhi = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);
  aX = aRho * std::cos(aPhi);
  aY = aRho * std::sin(aPhi);
  aZ = CLHEP::RandFlat::shoot(engine, fMinZ, fMaxZ);
  aT = CLHEP::RandFlat::shoot(engine, fMinT, fMaxT);

  return HepMC::FourVector(aX,aY,aZ,aT);
}

HepMC::FourVector FlatEvtVtxGenerator2::newVertexFlatD0(CLHEP::HepRandomEngine* engine, double phi) const {
  double aRho,aX,aY,aZ,aT;
  //aX = CLHEP::RandFlat::shoot(engine, fMinX, fMaxX);
  //aY = CLHEP::RandFlat::shoot(engine, fMinY, fMaxY);
  aRho = CLHEP::RandFlat::shoot(engine, -fMaxX, fMaxX);
  aX = aRho * -std::sin(phi);
  aY = aRho * std::cos(phi);
  aZ = CLHEP::RandFlat::shoot(engine, fMinZ, fMaxZ);
  aT = CLHEP::RandFlat::shoot(engine, fMinT, fMaxT);

  return HepMC::FourVector(aX,aY,aZ,aT);
}

void FlatEvtVtxGenerator2::produce( edm::Event& evt, const edm::EventSetup& )
{
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

  edm::Handle<edm::HepMCProduct> HepUnsmearedMCEvt ;

  evt.getByToken( sourceToken, HepUnsmearedMCEvt ) ;

  // Copy the HepMC::GenEvent
  HepMC::GenEvent* genevt = new HepMC::GenEvent(*HepUnsmearedMCEvt->GetEvent());
  std::unique_ptr<edm::HepMCProduct> HepMCEvt(new edm::HepMCProduct(genevt));
  // generate new vertex & apply the shift
  //
  if (fVtxSpectrum == "flatPhi") {
    HepMCEvt->applyVtxGen( newVertex(engine) ) ;

  } else if (fVtxSpectrum == "flatD0") {
    // Find the first particle momentum phi angle
    double genevt_phi = 0;
    for (auto it = genevt->particles_begin(); it != genevt->particles_end(); ++it) {
      genevt_phi = (*it)->momentum().phi();
      break;
    }
    HepMCEvt->applyVtxGen( newVertexFlatD0(engine, genevt_phi) ) ;
  }

  //HepMCEvt->LorentzBoost( 0., 142.e-6 );
  HepMCEvt->boostToLab( GetInvLorentzBoost(), "vertex" );
  HepMCEvt->boostToLab( GetInvLorentzBoost(), "momentum" );

  evt.put(std::move(HepMCEvt)) ;

  return ;
}
