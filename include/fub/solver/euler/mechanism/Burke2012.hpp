#ifndef FUB_EULER_SOLVER_MECHANISM_BURKE2012_HPP
#define FUB_EULER_SOLVER_MECHANISM_BURKE2012_HPP

#include "fub/core/assert.hpp"
#include "fub/core/span.hpp"
#include "fub/solver/euler/Reactor.hpp"

#include <string>
#include <vector>

namespace fub {
namespace euler {

struct Burke2012 : public FlameMasterMechanism {
  typedef enum SpeciesLabel {
    /* Computed species s.. */
    /* Steady-state species ss.. */
    sN2 = 0,
    sAR = 1,
    sH = 2,
    sO2 = 3,
    sO = 4,
    sOH = 5,
    sH2 = 6,
    sH2O = 7,
    sHE = 8,
    sHO2 = 9,
    sH2O2 = 10,
    sEnd
  } SpeciesLabel;

  typedef enum ReactionLabel {
    /* Reactions */
    r1f = 0,
    r1b = 1,
    r2f = 2,
    r2b = 3,
    r3f = 4,
    r3b = 5,
    r4f = 6,
    r4b = 7,
    r5f = 8,
    r5b = 9,
    r6f = 10,
    r6b = 11,
    r7f = 12,
    r7b = 13,
    r8f = 14,
    r8b = 15,
    r9f = 16,
    r9b = 17,
    r10f = 18,
    r10b = 19,
    r11f = 20,
    r11b = 21,
    r12f = 22,
    r12b = 23,
    r13f = 24,
    r13b = 25,
    r14f = 26,
    r14b = 27,
    r15f = 28,
    r15b = 29,
    r16f = 30,
    r16b = 31,
    r17f = 32,
    r17b = 33,
    r18f = 34,
    r18b = 35,
    r19f = 36,
    r19b = 37,
    r20f = 38,
    r20b = 39,
    r21f = 40,
    r21b = 41,
    r22f = 42,
    r22b = 43,
    r23f = 44,
    r23b = 45,
    r24f = 46,
    r24b = 47,
    r25f = 48,
    r25b = 49,
    r26f = 50,
    r26b = 51,
    r27f = 52,
    r27b = 53,
    /* PAHReactions */
    /* SootReactions */
    rEnd
  } ReactionLabel;

  typedef enum ThirdBodyLabel {
    mM1 = 0,
    mM2 = 1,
    mM3 = 2,
    mM4 = 3,
    mM5 = 4,
    mM6 = 5,
    mEnd
  } ThirdBodyLabel;

  void ComputeProductionRates(span<double> cdot, span<double> w, span<double> k,
                              span<double> c, span<double> M, double temp,
                              double pressure) const override;

  virtual void ComputeThermoData(span<double> h, span<double> cp, double T,
                                 span<double> /* s */) const override {
    ComputeThermoData(h, cp, T);
  };

  void ComputeThermoData(span<double> h, span<double> cp, double T) const;

  int GetNSpecies() const override { return sEnd; }

  int GetNSpecs() const override { return 11; }

  int GetNReactions() const override { return rEnd; }

  void GetMolarMass(span<double> W) const override {
    FUB_ASSERT(W.size() == GetNSpecies());
    W[sN2] = 2.80200000e+01;
    W[sAR] = 3.99480000e+01;
    W[sH] = 1.00800000e+00;
    W[sO2] = 3.20000000e+01;
    W[sO] = 1.60000000e+01;
    W[sOH] = 1.70080000e+01;
    W[sH2] = 2.01600000e+00;
    W[sH2O] = 1.80160000e+01;
    W[sHE] = 4.00000000e+00;
    W[sHO2] = 3.30080000e+01;
    W[sH2O2] = 3.40160000e+01;
  }

  std::vector<std::string> GetSpeciesNames() const override {
    return std::vector<std::string>{"N2", "AR",  "H",  "O2",  "O",   "OH",
                                    "H2", "H2O", "HE", "HO2", "H2O2"};
  }

  int GetNThirdBodyReactions() const override { return mEnd; }

  std::unique_ptr<FlameMasterMechanism> Clone() const override {
    return std::make_unique<Burke2012>();
  }
};

} // namespace euler
} // namespace fub

#endif
