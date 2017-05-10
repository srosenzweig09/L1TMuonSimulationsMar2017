#ifndef FlatRandomPtGunProducer2_H
#define FlatRandomPtGunProducer2_H

/** \class FlatRandomPtGunProducer2
 *
 * Generates single particle gun in HepMC format
 * Julia Yarba 12/2005
 ***************************************/

#include "L1TMuonSimulations/ParticleGuns/interface/BaseFlatGunProducer2.h"

namespace edm
{

  class FlatRandomPtGunProducer2 : public BaseFlatGunProducer2
  {

  public:
    FlatRandomPtGunProducer2(const ParameterSet & pset);
    virtual ~FlatRandomPtGunProducer2();

    virtual void produce(Event & e, const EventSetup& es) override;

  private:

    // data members

    double            fMinPt   ;
    double            fMaxPt   ;
    double            fMinOneOverPt   ;
    double            fMaxOneOverPt   ;
    double            fXFlatSpread    ;
    double            fYFlatSpread    ;
    double            fZFlatSpread    ;
    double            fRStarForPhi    ;
    double            fRStarForEta    ;
    bool              fRandomCharge   ;
    std::string       fPtSpectrum     ;

  };
}

#endif
