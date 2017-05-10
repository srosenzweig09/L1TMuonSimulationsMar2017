/*
 *  \author Julia Yarba
 */

#include <ostream>

#include "L1TMuonSimulations/ParticleGuns/interface/FlatRandomPtGunProducer2.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

namespace {
    double get_phi0_from_phiStar(double phiStar, double invPt, double rStar) {
        constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
        double dphi = - asin(mPtFactor * rStar * invPt);
        return phiStar - dphi;
    }
    double get_eta0_from_etaStar(double etaStar, double z0, double invPt, double rStar) {
        constexpr double mPtFactor = 0.3*3.8*1e-2/2.0;
        if (std::abs(rStar) < 1e-10)  return etaStar;
        if (std::abs(invPt) < 1e-10 && invPt <  0.) invPt = -1e-10;
        if (std::abs(invPt) < 1e-10 && invPt >= 0.) invPt = +1e-10;
        double cotStar = sinh(etaStar);
        double cot0 = (rStar * cotStar - z0) / (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt));
        return asinh(cot0);
    }
}

FlatRandomPtGunProducer2::FlatRandomPtGunProducer2(const ParameterSet& pset) :
   BaseFlatGunProducer2(pset)
{


  //ParameterSet defpset ;
  ParameterSet pgun_params =
    pset.getParameter<ParameterSet>("PGunParameters") ;

  fMinPt         = pgun_params.getParameter<double>("MinPt");
  fMaxPt         = pgun_params.getParameter<double>("MaxPt");
  fMinOneOverPt  = fMaxPt != 0.0 ? 1.0/fMaxPt : 1e-9;
  fMaxOneOverPt  = fMinPt != 0.0 ? 1.0/fMinPt : 1e-9;
  fXFlatSpread   = pgun_params.exists("XFlatSpread")   ? pgun_params.getParameter<double>("XFlatSpread")     : 0.;
  fYFlatSpread   = pgun_params.exists("YFlatSpread")   ? pgun_params.getParameter<double>("YFlatSpread")     : 0.;
  fZFlatSpread   = pgun_params.exists("ZFlatSpread")   ? pgun_params.getParameter<double>("ZFlatSpread")     : 0.;
  fRStarForPhi   = pgun_params.exists("RStarForPhi")   ? pgun_params.getParameter<double>("RStarForPhi")     : 0.;
  fRStarForEta   = pgun_params.exists("RStarForEta")   ? pgun_params.getParameter<double>("RStarForEta")     : 0.;
  fRandomCharge  = pgun_params.exists("RandomCharge")  ? pgun_params.getParameter<bool>("RandomCharge")      : false;
  fPtSpectrum    = pgun_params.exists("PtSpectrum")    ? pgun_params.getParameter<std::string>("PtSpectrum") : "flatPt";

  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

FlatRandomPtGunProducer2::~FlatRandomPtGunProducer2()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

void FlatRandomPtGunProducer2::produce(Event &e, const EventSetup& es)
{
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  if ( fVerbosity > 0 )
  {
    cout << " FlatRandomPtGunProducer2 : Begin New Event Generation" << endl;
  }
  // event loop (well, another step in it...)

  // no need to clean up GenEvent memory - done in HepMCProduct
  //

  // here re-create fEvt (memory)
  //
  fEvt = new HepMC::GenEvent() ;

  // now actualy, cook up the event from PDGTable and gun parameters
  //
  // 1st, primary vertex
  //
  //HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));
  double vx = fXFlatSpread > 0. ? CLHEP::RandFlat::shoot(engine, -fXFlatSpread, fXFlatSpread) : 0.;
  double vy = fYFlatSpread > 0. ? CLHEP::RandFlat::shoot(engine, -fYFlatSpread, fYFlatSpread) : 0.;
  double vz = fZFlatSpread > 0. ? CLHEP::RandFlat::shoot(engine, -fZFlatSpread, fZFlatSpread) : 0.;
  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(vx,vy,vz));

  // loop over particles
  //
  int barcode = 1 ;
  for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
  {
    double xx     = CLHEP::RandFlat::shoot(engine, 0.0, 1.0);
    double pt     = 0.0;
    if (fPtSpectrum == "flatPt")
    {
      pt     = fMinPt + xx * (fMaxPt - fMinPt);
    }
    else if (fPtSpectrum == "flatOneOverPt")
    {
      pt     = fMinOneOverPt + xx * (fMaxOneOverPt - fMinOneOverPt);
      if (pt != 0.0)  pt = 1.0/pt;
    }
    else if (fPtSpectrum == "flatOneOverPtCMS")
    {
      pt     = std::exp((1.-xx)*std::log(fMinOneOverPt)+
                            xx*std::log(fMaxOneOverPt)) ;
      if (pt != 0.0)  pt = 1.0/pt;
    }

    int PartID = fPartIDs[ip] ;
    if (fRandomCharge && (CLHEP::RandFlat::shoot(engine, 0.0, 1.0) < 0.5))
      PartID = - PartID;
    const HepPDT::ParticleData*
      PData = fPDGTable->particle(HepPDT::ParticleID(PartID)) ;
    double mass   = PData->mass().value() ;
    double phi    = CLHEP::RandFlat::shoot(engine, fMinPhi, fMaxPhi) ;
           phi    = get_phi0_from_phiStar(phi, PData->charge()/pt, fRStarForPhi);
    double eta    = CLHEP::RandFlat::shoot(engine, fMinEta, fMaxEta) ;
           eta    = get_eta0_from_etaStar(eta, vz, PData->charge()/pt, fRStarForEta);
    double theta  = 2.*atan(exp(-eta)) ;
    double mom    = pt/sin(theta) ;
    double px     = pt*cos(phi) ;
    double py     = pt*sin(phi) ;
    double pz     = mom*cos(theta) ;
    double energy2= mom*mom + mass*mass ;
    double energy = sqrt(energy2) ;
    HepMC::FourVector p(px,py,pz,energy) ;
    HepMC::GenParticle* Part =
      new HepMC::GenParticle(p,PartID,1);
    Part->suggest_barcode( barcode ) ;
    barcode++ ;
    Vtx->add_particle_out(Part);

    if ( fAddAntiParticle )
    {
      HepMC::FourVector ap(-px,-py,-pz,energy) ;
      int APartID = -PartID ;
      if ( PartID == 22 || PartID == 23 )
      {
        APartID = PartID ;
      }
      HepMC::GenParticle* APart =
        new HepMC::GenParticle(ap,APartID,1);
      APart->suggest_barcode( barcode ) ;
      barcode++ ;
      Vtx->add_particle_out(APart) ;
    }

  }

  fEvt->add_vertex(Vtx) ;
  fEvt->set_event_number(e.id().event()) ;
  fEvt->set_signal_process_id(20) ;

  if ( fVerbosity > 0 )
  {
    fEvt->print() ;
  }

  unique_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
  BProduct->addHepMCData( fEvt );
  e.put(std::move(BProduct), "unsmeared");

  unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(std::move(genEventInfo));

  if ( fVerbosity > 0 )
  {
    // for testing purpose only
    // fEvt->print() ; // prints empty info after it's made into edm::Event
    cout << " FlatRandomPtGunProducer2 : Event Generation Done " << endl;
  }
}

//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomPtGunProducer2);
