#ifndef L1TMuonSimulations_EMTFParticleTools_h
#define L1TMuonSimulations_EMTFParticleTools_h

#include <cmath>

struct EMTFInversePt {
  double operator()(double charge, double pt) const {
    double invpt = 1.0 / (pt + 1e-15);  // protects against division by zero
    invpt = (invpt < 1e-15) ? 1e-15 : invpt;  // protects against zero
    invpt *= charge;
    return invpt;
  }
};

struct EMTFDZero {
  double operator()(double invpt, double phi, double xv, double yv) const {
    double B = 3.811;  // in Tesla
    double R = -1.0 / (0.003 * B * invpt);  // R = -pT/(0.003 q B)  [cm]
    double xc = xv - (R * std::sin(phi));   // xc = xv - R sin(phi)
    double yc = yv + (R * std::cos(phi));   // yc = yv + R cos(phi)
    double d0 = R - (R / std::abs(R)) * std::hypot(xc, yc);  // d0 = R - sign(R) * sqrt(xc^2 + yc^2)
    return d0;
  }
};

#endif
