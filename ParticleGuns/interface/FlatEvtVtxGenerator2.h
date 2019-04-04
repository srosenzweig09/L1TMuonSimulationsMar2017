#ifndef IOMC_FlatEvtVtxGenerator2_H
#define IOMC_FlatEvtVtxGenerator2_H

/**
 * Generate event vertices according to a Flat distribution.
 * Attention: All values are assumed to be cm!
 *
 * Important note: flat independent distributions in Z and T are not correct for physics production
 * In reality, if two flat beams interact the real distribution will not be flat with independent Z and T
 * but Z and T will be correlated, as example in GaussEvtVtxGenerator.
 * Can restore correlation in configuration via MinT += (MinZ - MaxZ)/2 and MaxT += (MaxZ - MinZ)/2
 * in [ns] units (recall c_light = 29.98cm/ns)
 *
 */

#include "L1TMuonSimulations/ParticleGuns/interface/BaseEvtVtxGenerator2.h"

namespace CLHEP {
  class HepRandomEngine;
}

class FlatEvtVtxGenerator2 : public BaseEvtVtxGenerator2
{
public:
  FlatEvtVtxGenerator2(const edm::ParameterSet & p);
  ~FlatEvtVtxGenerator2() override;

  void produce( edm::Event&, const edm::EventSetup&) override;

  /// return a new event vertex
  //virtual CLHEP::Hep3Vector* newVertex();
  HepMC::FourVector newVertex(CLHEP::HepRandomEngine*) const override ;

  HepMC::FourVector newVertexFlatD0(CLHEP::HepRandomEngine* engine, double phi) const;

  const TMatrixD* GetInvLorentzBoost() const override {
    return nullptr;
  }

private:
  /** Copy constructor */
  FlatEvtVtxGenerator2(const FlatEvtVtxGenerator2 &p) = delete;
  /** Copy assignment operator */
  FlatEvtVtxGenerator2&  operator = (const FlatEvtVtxGenerator2 & rhs ) = delete;
private:
  edm::EDGetTokenT<edm::HepMCProduct> sourceToken;
  double fMinX, fMinY, fMinZ, fMinT;
  double fMaxX, fMaxY, fMaxZ, fMaxT;
  std::string fVtxSpectrum;
};

#endif
