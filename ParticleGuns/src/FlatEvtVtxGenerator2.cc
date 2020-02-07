#include "L1TMuonSimulations/ParticleGuns/interface/FlatEvtVtxGenerator2.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
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

  fVtxSpectrum = p.exists("VtxSpectrum") ? p.getParameter<std::string>("VtxSpectrum") : "flatVtx";

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
  double aX, aY, aZ, aT;
  aX = CLHEP::RandFlat::shoot(engine, fMinX, fMaxX);
  aY = CLHEP::RandFlat::shoot(engine, fMinY, fMaxY);
  aZ = CLHEP::RandFlat::shoot(engine, fMinZ, fMaxZ);
  aT = CLHEP::RandFlat::shoot(engine, fMinT, fMaxT);

  return HepMC::FourVector(aX,aY,aZ,aT);
}

HepMC::FourVector FlatEvtVtxGenerator2::newVertexFlatD0(CLHEP::HepRandomEngine* engine, double invpt, double phi) const {
  double aX,aY,aZ,aT;
  double fMaxRho = 150. * cm;
  aX = fMaxRho;
  aY = fMaxRho;

  int ntries = 0;  // try until rho is less than fMaxRho, rho = sqrt(aX*aX + aY*aY)

  while (!(std::hypot(aX, aY) < fMaxRho) && (ntries++ < 1000)) {
    // Throw two random numbers: d0 & rot. rot is a free rotation that allows
    // the particle to move along the circle. Assume the original phi is the
    // phi at the production vertex, the rotated phi is the phi at the point
    // of closest approach.
    double d0 = CLHEP::RandFlat::shoot(engine, fMinX, fMaxX);
    double rot = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);

    double phi_rotated = phi + rot;
    if (phi_rotated < -M_PI) {
      phi_rotated += M_PI * 2;
    }
    if (phi_rotated >= M_PI) {
      phi_rotated -= M_PI * 2;
    }

    // Find the center of the circle with phi at the point of closest approach
    double B = 3.811;  // in Tesla
    double R = -1.0 / (0.003 * B * invpt);  // R = -pT/(0.003 q B)  [cm]
    R *= cm;
    double xc = (d0 - R) * std::sin(phi_rotated);
    double yc = (d0 - R) * -std::cos(phi_rotated);

    // Find the production vertex with phi at the vertex
    aX = xc + R * std::sin(phi);  // xc = xv - R sin(phi), note that xv is aX
    aY = yc - R * std::cos(phi);  // yc = yv + R cos(phi), note that yv is aY
  }

  //aZ = CLHEP::RandFlat::shoot(engine, fMinZ, fMaxZ);
  //aT = CLHEP::RandFlat::shoot(engine, fMinT, fMaxT);
  double fSigmaZ = 5.0 * cm;  // use gaussian for vz
  aZ = CLHEP::RandGaussQ::shoot(engine, 0., fSigmaZ);
  aT = CLHEP::RandGaussQ::shoot(engine, 0., fSigmaZ);

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
  if (fVtxSpectrum == "flatVtx") {
    HepMCEvt->applyVtxGen( newVertex(engine) ) ;

  } else if (fVtxSpectrum == "flatD0") {
    // Get the q/pT and phi of the first particle
    double genevt_invpt = 0;
    double genevt_phi = 0;
    for (auto it = genevt->particles_begin(); it != genevt->particles_end(); ++it) {
      double genevt_charge = ((*it)->pdg_id() > 0) ? -1 : 1;  // Muon pdgid=13 has q=-1, pdgid=-13 has q=1. (might not work for other particles)
      double genevt_pt = (*it)->momentum().perp();
      genevt_invpt = genevt_charge / genevt_pt;
      genevt_phi = (*it)->momentum().phi();
      break;
    }
    HepMCEvt->applyVtxGen( newVertexFlatD0(engine, genevt_invpt, genevt_phi) ) ;
  }

  //HepMCEvt->LorentzBoost( 0., 142.e-6 );
  HepMCEvt->boostToLab( GetInvLorentzBoost(), "vertex" );
  HepMCEvt->boostToLab( GetInvLorentzBoost(), "momentum" );

  evt.put(std::move(HepMCEvt)) ;

  return ;
}
