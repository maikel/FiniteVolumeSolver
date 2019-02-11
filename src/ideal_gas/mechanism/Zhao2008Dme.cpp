#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <string>

namespace fub {
namespace ideal_gas {
static double GetLindRateCoeff(double temp, double pressure, double k0,
                               double kInf, double fc, double conc);

// static double GetPlogRateCoeff(double temp, double pressure, double lgt,
//                                double rt_inv, double* PlogP, double* PlogA,
//                                double* PlogB, double* PlogE, int np);

static double MAX_C(double X1, double X2);
// static double MIN_C(double X1, double X2);

void Zhao2008Dme::ComputeProductionRates(span<double> cdot, span<double> w,
                                         span<double> k, span<double> c,
                                         span<double> M, double temp,
                                         double pressure) const {
  /*
          This function computes rates of production cdot in [kmole/(m^3s)].
          The parameters w ( reaction rate ), k ( rate coefficient )
          and M ( third body concentrations ) are just work space for this
     function. c contains the concentrations of non steady state species in
     [kmole/m^3] and is workspace for the steady state concentrations, which are
     computed in this function. temp is the temperature in [K] and pressure is
     the pressure in [Pa]. Called functions are 'GetLindRateCoeff',
     'ComputeSteadyStates', 'CatchZero' and the functions that evaluate the
     broadening factors of the Troe formulation of pressure dependant rate
     coefficients 'Fc*'
  */

  // int	nSpec = 55;
  // int	nSpecIn = 55;
  double kTroe0, kTroeInf, fcTroe;
  double RGAS = 8314.34;
  double lgt = log(temp);
  double rt_inv = 1.0 / (RGAS * temp);
  // double PlogA[14], PlogB[14], PlogE[14] /* , PlogP[14] */;
  // int np;

  M[mM1] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H];

  M[mM2] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H];

  M[mM3] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           0x1.8p-1 * c[sAR] + 0x1.8p-1 * c[sHE];

  M[mM4] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           0x1.851eb851eb852p-2 * c[sAR] + 0x1.851eb851eb852p-2 * c[sHE];

  M[mM5] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.6p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           0x1.8f5c28f5c28f6p-1 * c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] +
           c[sHCCO] + c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] +
           0x1.e666666666666p+1 * c[sCO2] + c[sCH3HCO] + c[sOCHO] +
           c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] + c[sCH3OCH2] + c[sHCOOH] +
           c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] + c[sCH3OCO] + c[sCH3OCHO] +
           c[sCH3OCH2O] + c[sCH3OCH2OH] + c[sOCH2OCHO] + c[sHOCH2OCO] +
           c[sHOC2H4O2] + c[sCH3OCH2O2] + c[sCH2OCH2O2H] + c[sCH3OCH2O2H] +
           c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM6] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           0x1.47ae147ae147bp-1 * c[sAR] + 0x1.47ae147ae147bp-1 * c[sHE];

  M[mM7] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           0x1.bd70a3d70a3d7p-1 * c[sAR] + c[sHE];

  M[mM8] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           c[sAR] + c[sHE];

  M[mM9] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
           c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
           c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
           c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
           c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
           c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1.e666666666666p+1 * c[sCO2] +
           c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
           c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
           c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
           c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
           c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] +
           0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM10] = c[sH] + 0x1.4p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            c[sCH4] + c[sOH] + 0x1.8p+3 * c[sH2O] + c[sC2H] + c[sC2H2] +
            c[sC2H3] + 0x1.e666666666666p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
            c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
            c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] +
            0x1.e666666666666p+1 * c[sCO2] + c[sCH3HCO] + c[sOCHO] +
            c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] + c[sCH3OCH2] + c[sHCOOH] +
            c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] + c[sCH3OCO] + c[sCH3OCHO] +
            c[sCH3OCH2O] + c[sCH3OCH2OH] + c[sOCH2OCHO] + c[sHOCH2OCO] +
            c[sHOC2H4O2] + c[sCH3OCH2O2] + c[sCH2OCH2O2H] + c[sCH3OCH2O2H] +
            c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] +
            c[sHE];

  M[mM11] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + 0x1.4p+2 * c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] +
            0x1p+1 * c[sCO] + c[sN2] + c[sC2H4] + c[sHCO] + c[sC2H5] +
            c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] + c[sO2] + c[sCH3OH] +
            c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] +
            c[sCH3CO] + 0x1.8p+1 * c[sCO2] + c[sCH3HCO] + c[sOCHO] +
            c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] + c[sCH3OCH2] + c[sHCOOH] +
            c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] + c[sCH3OCO] + c[sCH3OCHO] +
            c[sCH3OCH2O] + c[sCH3OCH2OH] + c[sOCH2OCHO] + c[sHOCH2OCO] +
            c[sHOC2H4O2] + c[sCH3OCH2O2] + c[sCH2OCH2O2H] + c[sCH3OCH2O2H] +
            c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM12] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM13] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM14] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM15] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM16] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM17] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM18] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM19] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM20] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM21] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM22] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM23] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM24] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM25] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM26] = c[sH] + 0x1p+1 * c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] +
            0x1p+1 * c[sCH4] + c[sOH] + 0x1.8p+2 * c[sH2O] + c[sC2H] +
            c[sC2H2] + c[sC2H3] + 0x1.8p+0 * c[sCO] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + 0x1p+1 * c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + 0x1.6666666666666p-1 * c[sAR] + c[sHE];

  M[mM27] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sN2] + c[sC2H4] +
            c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] + c[sCH3O] +
            c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] + c[sCH2CO] +
            c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCH3HCO] + c[sOCHO] +
            c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] + c[sCH3OCH2] + c[sHCOOH] +
            c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] + c[sCH3OCO] + c[sCH3OCHO] +
            c[sCH3OCH2O] + c[sCH3OCH2OH] + c[sOCH2OCHO] + c[sHOCH2OCO] +
            c[sHOC2H4O2] + c[sCH3OCH2O2] + c[sCH2OCH2O2H] + c[sCH3OCH2O2H] +
            c[sHO2CH2OCHO] + c[sO2CH2OCH2O2H] + c[sHE];

  M[mM28] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM29] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM30] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM31] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM32] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM33] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  M[mM34] = c[sH] + c[sH2] + c[sCH2] + c[sCH2S] + c[sCH3] + c[sO] + c[sCH4] +
            c[sOH] + c[sH2O] + c[sC2H] + c[sC2H2] + c[sC2H3] + c[sCO] + c[sN2] +
            c[sC2H4] + c[sHCO] + c[sC2H5] + c[sCH2O] + c[sC2H6] + c[sCH2OH] +
            c[sCH3O] + c[sO2] + c[sCH3OH] + c[sHO2] + c[sH2O2] + c[sHCCO] +
            c[sCH2CO] + c[sHCCOH] + c[sCH2HCO] + c[sCH3CO] + c[sCO2] +
            c[sCH3HCO] + c[sOCHO] + c[sCH3CHOH] + c[sC2H4OH] + c[sCH3CH2O] +
            c[sCH3OCH2] + c[sHCOOH] + c[sCH3OCH3] + c[sC2H5OH] + c[sHOCH2O] +
            c[sCH3OCO] + c[sCH3OCHO] + c[sCH3OCH2O] + c[sCH3OCH2OH] +
            c[sOCH2OCHO] + c[sHOCH2OCO] + c[sHOC2H4O2] + c[sCH3OCH2O2] +
            c[sCH2OCH2O2H] + c[sCH3OCH2O2H] + c[sHO2CH2OCHO] +
            c[sO2CH2OCH2O2H] + c[sAR] + c[sHE];

  k[r1f] = 0x1.9cecd667p+41 *
           exp(-0x1.9fbe76c8b4396p-2 * lgt - 0x1.08ee7ap+26 * rt_inv);
  k[r1b] = 0x1.a4648c5d1ed11p+32 *
           exp(0x1.da7127a19cbap-6 * lgt + 0x1.0f77bdd2ba253p+20 * rt_inv);
  k[r2f] = 0x1.9666666666666p+5 *
           exp(0x1.55c28f5c28f5cp+1 * lgt - 0x1.91923p+24 * rt_inv);
  k[r2b] = 0x1.dcbbd77d5716bp+4 *
           exp(0x1.51468cbb3dc8p+1 * lgt - 0x1.399612bb18f4fp+24 * rt_inv);
  k[r3f] = 0x1.a5ep+17 *
           exp(0x1.828f5c28f5c29p+0 * lgt - 0x1.b5f61ffffffffp+23 * rt_inv);
  k[r3b] = 0x1.0ec8e7c618f0ap+21 *
           exp(0x1.6a979f0f7ce88p+0 * lgt - 0x1.25f639fec6edcp+26 * rt_inv);
  k[r4f] = 0x1.733ffffffffffp+11 *
           exp(0x1.028f5c28f5c29p+1 * lgt - 0x1.abbf200000001p+25 * rt_inv);
  k[r4b] = 0x1.533fab4886728p+7 *
           exp(0x1.0a0f3814c7365p+1 * lgt + 0x1.7ab76a800528bp+23 * rt_inv);
  k[r5f] = 0x1.4537251ebd4p+55 *
           exp(-0x1.6666666666666p+0 * lgt - 0x1.a07e8a0000001p+28 * rt_inv);
  k[r5b] = 0x1.64dfcd2d09c4cp+47 *
           exp(-0x1.bfb15b3d46e8ep+0 * lgt - 0x1.bd1b631bd54c3p+21 * rt_inv);
  k[r6f] = 0x1.4bf72f57dp+52 *
           exp(-0x1.199999999999ap+0 * lgt - 0x1.a07e8a0000001p+28 * rt_inv);
  k[r6b] = 0x1.6c480fafc9f0fp+44 *
           exp(-0x1.72e48e707a117p+0 * lgt - 0x1.bd1b631bd51f2p+21 * rt_inv);
  k[r7f] = 0x1.4bf72f57dp+52 *
           exp(-0x1.199999999999ap+0 * lgt - 0x1.a07e8a0000001p+28 * rt_inv);
  k[r7b] = 0x1.6c480fafca22cp+44 *
           exp(-0x1.72e48e707a159p+0 * lgt - 0x1.bd1b631bd52eap+21 * rt_inv);
  k[r8f] = 0x1.6f766f4p+32 * exp(-0x1p-1 * lgt);
  k[r8b] = 0x1.81d71449bd1ep+48 *
           exp(-0x1.3e0ce57746a8p-1 * lgt - 0x1.dacfa7a34cb29p+28 * rt_inv);
  k[r9f] = 0x1.1fc7ep+24 * exp(0x1.c89a8p+22 * rt_inv);
  k[r9b] = 0x1.2e2c59236b28ep+40 *
           exp(-0x1.f0672bba39cp-4 * lgt - 0x1.d3ad3da34cb36p+28 * rt_inv);
  k[r10f] = 0x1.1fc7ep+24 * exp(0x1.c89a8p+22 * rt_inv);
  k[r10b] = 0x1.2e2c59236b03p+40 *
            exp(-0x1.f0672bba398p-4 * lgt - 0x1.d3ad3da34cb36p+28 * rt_inv);
  k[r11f] = 0x1.126412e9p+42 * exp(-0x1p+0 * lgt);
  k[r11b] = 0x1.2552fd7938fe1p+49 *
            exp(-0x1.5f5a20d5e164p-1 * lgt - 0x1.9784916579f5bp+28 * rt_inv);
  k[r12f] = 0x1.0e0198eaeep+55 * exp(-0x1p+1 * lgt);
  k[r12b] = 0x1.3bdd3274f5b6cp+66 *
            exp(-0x1.beacc84293p+0 * lgt - 0x1.d8d230b97a21dp+28 * rt_inv);
  kTroe0 = 0x1.217dfe6e98p+49 *
           exp(-0x1.b851eb851eb85p+0 * lgt - 0x1.0c09999999999p+21 * rt_inv);
  kTroeInf = 0x1.5faadbp+30 * exp(0x1.3333333333333p-1 * lgt);
  fcTroe = 0x1.9999999999998p-3 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.999999999999ap-1 * exp(-temp / 0x1.93e5939a08ceap+99);
  k[r13f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  kTroe0 = 0x1.90fa486a5b84p+59 *
           exp(-0x1.bb9fd1d2d36cp+0 * lgt - 0x1.8a8edd9cb40c4p+27 * rt_inv);
  kTroeInf = 0x1.e718d77e6d81dp+40 *
             exp(0x1.2c976697c9cbdp-1 * lgt - 0x1.865eb7364da5ep+27 * rt_inv);
  k[r13b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  k[r14f] = 0x1.eeb7f3p+33 * exp(-0x1.a457400000001p+21 * rt_inv);
  k[r14b] = 0x1.457c1e1929d6bp+31 *
            exp(0x1.72636c926a08p-2 * lgt - 0x1.ba3b4c3d432d6p+27 * rt_inv);
  k[r15f] = 0x1.07b69ad8p+36 * exp(-0x1.2d568p+20 * rt_inv);
  k[r15b] = 0x1.9e6b84b3207f1p+23 *
            exp(0x1.85f47066ed14p-1 * lgt - 0x1.246eec1900cdbp+27 * rt_inv);
  k[r16f] = 0x1.e449a94p+34;
  k[r16b] = 0x1.75c3c10662642p+31 *
            exp(0x1.4e83578b0f8cp-2 * lgt - 0x1.a8aa6b94a647bp+27 * rt_inv);
  k[r17f] = 0x1.aea4c04p+34 * exp(0x1.fbad800000001p+20 * rt_inv);
  k[r17b] = 0x1.6bb6fdab173d6p+35 *
            exp(0x1.1284782c869p-2 * lgt - 0x1.13a7279e534f7p+28 * rt_inv);
  k[r18f] = 0x1.8727cdap+38 * exp(-0x1.7e7b68p+25 * rt_inv);
  k[r18b] = 0x1.ce90f36618f5ep+45 *
            exp(-0x1.8e9d7f16622cp-2 * lgt - 0x1.94c280edaee52p+27 * rt_inv);
  k[r19f] = 0x1.efe92p+26 * exp(0x1.a0137cccccccdp+22 * rt_inv);
  k[r19b] = 0x1.253911a921d93p+34 *
            exp(-0x1.8e9d7f1664f4p-2 * lgt - 0x1.28230b0748799p+27 * rt_inv);
  kTroe0 = 0x1.b548f9354p+46 * exp(-0x1.6b1b14p+27 * rt_inv);
  kTroeInf = 0x1.0c6452ac58p+48 * exp(-0x1.827cfap+27 * rt_inv);
  fcTroe = 0x1p-1 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1p-1 * exp(-temp / 0x1.93e5939a08ceap+99);
  k[r20f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  kTroe0 = 0x1.a38846fa5523dp+16 *
           exp(0x1.29ef7e46c6988p+0 * lgt + 0x1.72985857dd30bp+24 * rt_inv);
  kTroeInf = 0x1.017eee39e99b7p+18 *
             exp(0x1.29ef7e46c6988p+0 * lgt + 0x1.6f1250afba61cp+23 * rt_inv);
  k[r20b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  k[r21f] = 0x1.671e344p+34 * exp(-0x1.fae9600000002p+23 * rt_inv);
  k[r21b] = 0x1.930e1762cf943p+15 *
            exp(0x1.6b42b6043162p+0 * lgt - 0x1.1bf26c33fc481p+28 * rt_inv);
  k[r22f] = 0x1.671e344p+35 * exp(-0x1.fb8ccffffffffp+24 * rt_inv);
  k[r22b] = 0x1.8f972f34501a2p+25 *
            exp(0x1.808075d46812p-1 * lgt - 0x1.7befc49f289e8p+26 * rt_inv);
  k[r23f] = 0x1.2a6ffffffffffp+13 *
            exp(0x1p+1 * lgt - 0x1.fae9600000002p+23 * rt_inv);
  k[r23b] = 0x1.858a9bf2d9042p+2 *
            exp(0x1.5ba41ad42e96cp+1 * lgt - 0x1.266ab54deed3p+26 * rt_inv);
  k[r24f] = 0x1.dcd65p+29;
  k[r24b] = 0x1.548e5adeb8841p+23 *
            exp(0x1.5090fba1755ap-1 * lgt - 0x1.ec44069def81bp+26 * rt_inv);
  k[r25f] = 0x1.0e15635p+39 * exp(-0x1.31129c0000001p+25 * rt_inv);
  k[r25b] = 0x1.81c942f055162p+32 *
            exp(0x1.5090fba17662p-1 * lgt - 0x1.4266aa4ef7c3fp+27 * rt_inv);
  kTroe0 = 0x1.582b4c9a9dbp+60 *
           exp(-0x1.651eb851eb852p+1 * lgt - 0x1.0b90a80000001p+24 * rt_inv);
  kTroeInf = 0x1.12a88p+24 * exp(-0x1.3067p+23 * rt_inv);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r26f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM7]);
  kTroe0 = 0x1.ea6ef219747f3p+89 *
           exp(-0x1.e0b1f9b9f9f2p+1 * lgt - 0x1.075ddc300f9bcp+29 * rt_inv);
  kTroeInf = 0x1.8761c074bce31p+53 *
             exp(-0x1.ee4d05a039b38p-1 * lgt - 0x1.03c2f2f00f9bbp+29 * rt_inv);
  k[r26b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM7]);
  k[r27f] = 0x1.2d9979p+31 * exp(-0x1.7ca99c0000001p+27 * rt_inv);
  k[r27b] = 0x1.994d75fbaacd6p+44 *
            exp(-0x1.b0402028ee46p-1 * lgt - 0x1.c30fa879a518p+27 * rt_inv);
  k[r28f] = 0x1.c086634p+34 * exp(-0x1.6f189p+26 * rt_inv);
  k[r28b] = 0x1.d5c7e0bb0e617p+44 *
            exp(-0x1.08fe746368fcp-1 * lgt - 0x1.534e600725aaep+28 * rt_inv);
  k[r29f] = 0x1.bdcccccccccccp+7 *
            exp(0x1.e3d70a3d70a3dp+0 * lgt + 0x1.27e6033333334p+22 * rt_inv);
  k[r29b] = 0x1.2920255b95866p+30 *
            exp(0x1.38bb2fb088a3p-1 * lgt - 0x1.877a11b761ceap+26 * rt_inv);
  k[r30f] = 0x1.c4da2dp+28 *
            exp(0x1.516872b020c4ap-1 * lgt - 0x1.dacc780000001p+25 * rt_inv);
  k[r30b] = 0x1.73f3cfda0196ap+16 *
            exp(0x1.c86dacb0fe32dp-1 * lgt + 0x1.dfe65ef2f5921p+20 * rt_inv);
  k[r31f] = 0x1.c3cd9fp+32 * exp(-0x1.a2cfp+20 * rt_inv);
  k[r31b] = 0x1.01000f13a4808p+31 *
            exp(0x1.c1a5b595cab8p-3 * lgt - 0x1.0f316a7867b3cp+27 * rt_inv);
  k[r32f] = 0x1.0d56a73p+36;
  k[r32b] = 0x1.9331edb45b351p+31 *
            exp(0x1.299b23aea7c4p-1 * lgt - 0x1.5fcadddad56fap+28 * rt_inv);
  k[r33f] = 0x1.c203db8p+34;
  k[r33b] = 0x1.8b205a2313279p+29 *
            exp(0x1.17ab192afa84p-1 * lgt - 0x1.5a4b1c0686fd1p+28 * rt_inv);
  k[r34f] = 0x1.c203db8p+34;
  k[r34b] = 0x1.b0660acbfe365p+33 *
            exp(0x1.f35752f76b5p-2 * lgt - 0x1.9b98bb5a8729fp+28 * rt_inv);
  k[r35f] = 0x1.bf08ebp+34;
  k[r35b] = 0x1.059b89495555bp+52 *
            exp(-0x1.7747cb9f5dccp-1 * lgt - 0x1.c0c938812c3d4p+28 * rt_inv);
  k[r36f] = 0x1.bf08ebp+34;
  k[r36b] = 0x1.8093260caab7dp+32 *
            exp(-0x1.23f274c51efcp-2 * lgt - 0x1.749d8d50655d8p+27 * rt_inv);
  k[r37f] = 0x1.65a0bcp+31;
  k[r37b] = 0x1.b7b9061481513p+14 *
            exp(0x1.a0a05daf83e4p-1 * lgt - 0x1.2291687be2768p+28 * rt_inv);
  k[r38f] = 0x1.8ae17a4p+34;
  k[r38b] = 0x1.888c4f817eb9p+38 *
            exp(0x1.2a11943c9a4p-2 * lgt - 0x1.675caea4c7543p+28 * rt_inv);
  k[r39f] = 0x1.bf08ebp+34;
  k[r39b] = 0x1.a9704b92a26a1p+39 *
            exp(0x1.f20454f0376p-4 * lgt - 0x1.2bc40712cff68p+28 * rt_inv);
  k[r40f] = 0x1.3dc747e5d87efp+121 *
            exp(-0x1.9333333333333p+2 * lgt - 0x1.8e9e4a0000001p+28 * rt_inv);
  k[r40b] = 0x1.12421bb49bf99p+104 *
            exp(-0x1.8c1a9d46d7482p+2 * lgt - 0x1.2d066c71e916cp+25 * rt_inv);
  k[r41f] =
      0x1.1cb0a65614917p+141 * exp(-0x1p+3 * lgt - 0x1.8514f1p+28 * rt_inv);
  k[r41b] = 0x1.6fcf853e7a7a7p+119 *
            exp(-0x1.d3b4059dcf2b8p+2 * lgt - 0x1.7be2526912929p+28 * rt_inv);
  k[r42f] = 0x1.c07p+15 *
            exp(0x1.e666666666666p+0 * lgt - 0x1.5ef4dccccccccp+23 * rt_inv);
  k[r42b] = 0x1.60b0a19ffbe08p+6 *
            exp(0x1.2e09d9775b6e8p+1 * lgt - 0x1.f7f3ed735f06bp+25 * rt_inv);
  k[r43f] = 0x1.0db6054p+34 * exp(-0x1.8945800000001p+23 * rt_inv);
  k[r43b] = 0x1.f1abb2d4013e5p+23 *
            exp(0x1.b2d51d19e746p-2 * lgt - 0x1.d68a079db842cp+25 * rt_inv);
  k[r44f] =
      0x1.a2b38p+21 * exp(0x1.2e147ae147ae1p+0 * lgt + 0x1.c89a8p+20 * rt_inv);
  k[r44b] = 0x1.a6bb9103ce219p+15 *
            exp(0x1.8bca0a501f0ap+0 * lgt - 0x1.b830671edccfep+26 * rt_inv);
  k[r45f] = 0x1.338p+10 * exp(0x1.8p+1 * lgt - 0x1.9efa6p+27 * rt_inv);
  k[r45b] = 0x1.6f970d2421c08p+3 *
            exp(0x1.8c8a38b1dc3bap+1 * lgt - 0x1.4d78794b1f693p+25 * rt_inv);
  k[r46f] =
      0x1.48cccccccccccp+5 * exp(0x1.4p+1 * lgt - 0x1.45ead8p+25 * rt_inv);
  k[r46b] = 0x1.d0cf35e3ade0ep+5 *
            exp(0x1.1ab688cf0e96p+1 * lgt - 0x1.d810da03b55eep+24 * rt_inv);
  k[r47f] = 0x1.f3ba607d8040fp-29 *
            exp(0x1.5ae147ae147aep+2 * lgt - 0x1.fdb8800000001p+21 * rt_inv);
  k[r47b] = 0x1.04ff27e81995p-29 *
            exp(0x1.65ba4f9e1d3eep+2 * lgt - 0x1.fca0c48fbaf3p+25 * rt_inv);
  k[r48f] = 0x1.3a0abebp+36;
  k[r48b] = 0x1.6d3fc20435aa2p+43 *
            exp(-0x1.2daad01b6db8p-2 * lgt - 0x1.1a7070325b791p+28 * rt_inv);
  k[r49f] = 0x1.c4793ec698p+50 *
            exp(-0x1.91eb851eb851fp+0 * lgt - 0x1.d287f4p+26 * rt_inv);
  k[r49b] = 0x1.3a528b1335e65p+58 *
            exp(-0x1.06ead9f133bb1p+1 * lgt - 0x1.c50a8fca204bep+22 * rt_inv);
  k[r50f] = 0x1.64ac98p+28 * exp(-0x1.d3544p+25 * rt_inv);
  k[r50b] = 0x1.a655bd9c2d15fp+26 *
            exp(0x1.1f75724ec08p-3 * lgt - 0x1.118fe1f488bep+28 * rt_inv);
  k[r51f] = 0x1.6fbcap+24 *
            exp(0x1.851eb851eb852p-1 * lgt + 0x1.28de700000001p+23 * rt_inv);
  k[r51b] = 0x1.8a50cd55d381p+28 *
            exp(0x1.348c069011aep-1 * lgt - 0x1.7601be25ee7e4p+26 * rt_inv);
  kTroe0 = 0x1.0a7c08982a5e9p+86 *
           exp(-0x1.ep+1 * lgt - 0x1.f558333333335p+21 * rt_inv);
  kTroeInf = 0x1.0913e359p+41 *
             exp(-0x1.6147ae147ae14p-1 * lgt - 0x1.653bc7ae147afp+19 * rt_inv);
  fcTroe = 0x1p+0 * exp(-temp / 0x1.1dp+9) +
           0x1p+0 * exp(-0x1.93e5939a08ceap+99 / temp);
  k[r52f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM11]);
  kTroe0 = 0x1.2cfea836b80dep+123 *
           exp(-0x1.46ddc188fa4e8p+2 * lgt - 0x1.725a038b46b64p+28 * rt_inv);
  kTroeInf = 0x1.2b67df3fbdc91p+78 *
             exp(-0x1.060d6e9713555p+1 * lgt - 0x1.6f21f108b75a1p+28 * rt_inv);
  k[r52b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM11]);
  kTroe0 = 0x1.001da4e1fe9aep+91 *
           exp(-0x1.30a3d70a3d70ap+2 * lgt - 0x1.378d800000001p+23 * rt_inv);
  kTroeInf = 0x1.719e5fa3p+43 *
             exp(-0x1.428f5c28f5c29p-1 * lgt - 0x1.873a800000001p+20 * rt_inv);
  fcTroe = 0x1.bc6a7ef9db22cp-3 * exp(-temp / 0x1.28p+6) +
           0x1.90e5604189375p-1 * exp(-temp / 0x1.6fap+11) +
           0x1p+0 * exp(-0x1.b34p+12 / temp);
  k[r53f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM12]);
  kTroe0 = 0x1.35faa6b0f798cp+107 *
           exp(-0x1.2ce365068f44p+2 * lgt - 0x1.ae529003ba4ccp+28 * rt_inv);
  kTroeInf = 0x1.bf5a416b3244cp+59 *
             exp(-0x1.248bcc0b845d9p-1 * lgt - 0x1.a61d5e83ba4ccp+28 * rt_inv);
  k[r53b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM12]);
  k[r54f] = 0x1.ab57fffffffffp+15 *
            exp(0x1.f851eb851eb85p+0 * lgt - 0x1.65d6b80000001p+25 * rt_inv);
  k[r54b] = 0x1.41c38b03957d3p+7 *
            exp(0x1.214d8c26a646ep+1 * lgt - 0x1.294831b070f22p+25 * rt_inv);
  k[r55f] = 0x1.77825fp+31 * exp(0x1p-1 * lgt - 0x1.4878980000001p+25 * rt_inv);
  k[r55b] = 0x1.4baad8562e975p+22 *
            exp(0x1.82a24f0cae3bp-1 * lgt - 0x1.bfd8061bfabdep+24 * rt_inv);
  k[r56f] =
      0x1.658p+12 * exp(0x1.f5c28f5c28f5cp+0 * lgt - 0x1.50f65p+23 * rt_inv);
  k[r56b] = 0x1.598be1dfc2fa1p+7 *
            exp(0x1.1409ff856eb1cp+1 * lgt - 0x1.f61df9adfeb79p+25 * rt_inv);
  k[r57f] = 0x1.78b38cp+31;
  k[r57b] = 0x1.4929815c7b969p+37 *
            exp(0x1.24fae5c6d3ap-4 * lgt - 0x1.c2cd90d126f4bp+27 * rt_inv);
  k[r58f] = 0x1.593ae8p+27 * exp(-0x1.288cb8p+26 * rt_inv);
  k[r58b] = 0x1.d3380e1fe65e7p+28 *
            exp(-0x1.d7dc388815c6dp-2 * lgt - 0x1.a71c8721fcc14p+21 * rt_inv);
  k[r59f] = 0x1.74876e8p+36 * exp(-0x1.909d08p+26 * rt_inv);
  k[r59b] = 0x1.bf52411278505p+23 *
            exp(0x1.53f8ac9db47c8p-2 * lgt + 0x1.95cd88959ccecp+23 * rt_inv);
  k[r60f] = 0x1.65a0bcp+32;
  k[r60b] = 0x1.87554ff3dd129p+27 *
            exp(0x1.5c923ffca368p-1 * lgt - 0x1.2c2ea4f51b7d9p+28 * rt_inv);
  k[r61f] = 0x1.66ee8538p+36;
  k[r61b] = 0x1.8c23afaa59919p+23 *
            exp(0x1.e1779d86acd4ep-1 * lgt - 0x1.87ce5dce324fep+23 * rt_inv);
  k[r62f] = 0x1.38eca48p+35;
  k[r62b] = 0x1.91ad529758247p+29 *
            exp(0x1.4aa23578f62cp-1 * lgt - 0x1.26aee320cd0aep+28 * rt_inv);
  k[r63f] = 0x1.65a0bcp+34;
  k[r63b] = 0x1.f65ccd224c855p+32 *
            exp(0x1.2ca2c5c9b17p-1 * lgt - 0x1.67fc8274cd378p+28 * rt_inv);
  k[r64f] = 0x1.c0e5c15p+37 * exp(-0x1.404c98p+24 * rt_inv);
  k[r64b] = 0x1.754d3e7b4d59cp+35 *
            exp(0x1.46c11366dd96p-2 * lgt - 0x1.9979db59e7a1ep+26 * rt_inv);
  k[r65f] = 0x1.5f93037cp+40 * exp(-0x1p+0 * lgt);
  k[r65b] = 0x1.245e392ae9979p+38 *
            exp(-0x1.5c9f764c91abp-1 * lgt - 0x1.4966b559e79d5p+26 * rt_inv);
  k[r66f] = 0x1.65a0bcp+33;
  k[r66b] = 0x1.5fb26f91b7fcap+38 *
            exp(-0x1.1f71aebe1f6p-4 * lgt - 0x1.d9d7019aa2aebp+27 * rt_inv);
  k[r67f] = 0x1.2a05f2p+33;
  k[r67b] = 0x1.ef857c1ea2fe3p+35 *
            exp(0x1.76bdfb36b5e8p-2 * lgt - 0x1.42d11c398f44cp+28 * rt_inv);
  k[r68f] = 0x1.bf08ebp+33;
  k[r68b] = 0x1.36fb740be139dp+38 *
            exp(0x1.c4de9bb00a7p-3 * lgt - 0x1.f04f9c5a2c082p+27 * rt_inv);
  k[r69f] = 0x1.65a0bcp+31;
  k[r69b] = 0x1.b2a7139be8f3cp+33 *
            exp(0x1.dcac33d2ad68p-2 * lgt - 0x1.0f34e353d5521p+28 * rt_inv);
  k[r70f] = 0x1.65a0bcp+34;
  k[r70b] = 0x1.209677b3255edp+33 *
            exp(0x1.c6cb3d16e0f8p-2 * lgt - 0x1.2d392ce868867p+28 * rt_inv);
  k[r71f] = 0x1.7970b794fp+49 *
            exp(-0x1.3333333333333p+0 * lgt - 0x1.eec8100000001p+25 * rt_inv);
  k[r71b] = 0x1.2cea2600a420fp+33 *
            exp(-0x1.c75a8b757376cp-1 * lgt + 0x1.4f862b019b31ap+24 * rt_inv);
  k[r72f] = 0x1.dcd65p+34;
  k[r72b] = 0x1.5d6a4efd84515p+18 *
            exp(0x1.d6872228c681p-1 * lgt - 0x1.5215e418263aap+25 * rt_inv);
  k[r73f] = 0x1.65a0bcp+32;
  k[r73b] = 0x1.30caeda5e05cfp+23 *
            exp(0x1.3fb1ba1b1p-1 * lgt - 0x1.44b32cb5603f8p+28 * rt_inv);
  k[r74f] = 0x1.0c388dp+34;
  k[r74b] = 0x1.f450ccb3e209fp+28 *
            exp(0x1.21b24a6bcb58p-1 * lgt - 0x1.8600cc09606cp+28 * rt_inv);
  k[r75f] = 0x1.508166a8p+36 * exp(-0x1.7e6b100000001p+25 * rt_inv);
  k[r75b] = 0x1.739878d20861dp+30 *
            exp(0x1.30e01cab12ap-2 * lgt - 0x1.4056b1d61a3c7p+27 * rt_inv);
  k[r76f] = 0x1.4fb18p+24 * exp(-0x1.be63800000001p+22 * rt_inv);
  k[r76b] = 0x1.72b2e4361e6d7p+18 *
            exp(0x1.30e01cab1074p-2 * lgt - 0x1.dd5e13ac3471dp+26 * rt_inv);
  k[r77f] = 0x1.1e1a3p+28;
  k[r77b] = 0x1.759d6dc6c0c8p+29 *
            exp(-0x1.76f589ad518p-4 * lgt - 0x1.0aefca61e48c1p+28 * rt_inv);
  k[r78f] = 0x1.dcd65p+33 * exp(-0x1.78ac200000001p+25 * rt_inv);
  k[r78b] = 0x1.d1c5358adc045p+39 *
            exp(-0x1.70d78543217p-2 * lgt - 0x1.7facb9fb54138p+27 * rt_inv);
  k[r79f] = 0x1.4f46b04p+36;
  k[r79b] = 0x1.722076d4f5b36p+35 *
            exp(0x1.60dd047ae97p-2 * lgt - 0x1.60d565ce22795p+28 * rt_inv);
  k[r80f] = 0x1.bf08ebp+35;
  k[r80b] = 0x1.df0514e204f0ep+30 *
            exp(0x1.b0ea465b1548p-2 * lgt - 0x1.4b3d767cfbb9dp+28 * rt_inv);
  kTroe0 = 0x1.93e5939a08ceap+101 *
           exp(-0x1.7ae147ae147aep+2 * lgt - 0x1.90eecp+23 * rt_inv);
  kTroeInf = 0x1.3d2fafdd8cp+51 *
             exp(-0x1.6e147ae147ae1p+0 * lgt - 0x1.53a4cp+22 * rt_inv);
  fcTroe = 0x1.2d0e560418938p-1 * exp(-temp / 0x1.86p+7) +
           0x1.a5e353f7ced91p-2 * exp(-temp / 0x1.70cp+12) +
           0x1p+0 * exp(-0x1.8fap+12 / temp);
  k[r81f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM15]);
  kTroe0 = 0x1.72699d89e8cd4p+129 *
           exp(-0x1.ae8502eb9a308p+2 * lgt - 0x1.805394aa10abep+28 * rt_inv);
  kTroeInf = 0x1.22e3fee5a20acp+79 *
             exp(-0x1.1e51b3ebaf424p+1 * lgt - 0x1.791ab1aa10abfp+28 * rt_inv);
  k[r81b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM15]);
  kTroe0 = 0x1.208545e35d9d6p+85 *
           exp(-0x1.299999999999ap+2 * lgt - 0x1.44524p+24 * rt_inv);
  kTroeInf = 0x1.f7102ep+29 * exp(0x1p-1 * lgt - 0x1.5f64p+18 * rt_inv);
  fcTroe = 0x1.999999999999ap-2 * exp(-temp / 0x1.9p+6) +
           0x1.3333333333333p-1 * exp(-temp / 0x1.5f9p+16) +
           0x1p+0 * exp(-0x1.388p+13 / temp);
  k[r82f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM16]);
  kTroe0 = 0x1.2407d419ca666p+100 *
           exp(-0x1.210e612649aap+2 * lgt - 0x1.944fb598823e1p+28 * rt_inv);
  kTroeInf = 0x1.fd2ee5a1ecfb7p+44 *
             exp(0x1.4459c39a7f7dp-1 * lgt - 0x1.80626a98823ep+28 * rt_inv);
  k[r82b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM16]);
  kTroe0 = 0x1.66fe4ba326d8p+118 *
           exp(-0x1.dc28f5c28f5c3p+2 * lgt - 0x1.c174000000001p+25 * rt_inv);
  kTroeInf =
      0x1.21adb7p+31 * exp(0x1.07ae147ae147bp-1 * lgt - 0x1.9898p+17 * rt_inv);
  fcTroe = 0x1.3333333333334p-2 * exp(-temp / 0x1.9p+6) +
           0x1.6666666666666p-1 * exp(-temp / 0x1.5f9p+16) +
           0x1p+0 * exp(-0x1.388p+13 / temp);
  k[r83f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM17]);
  kTroe0 = 0x1.e281c11c2e502p+129 *
           exp(-0x1.d4fbccbafc138p+2 * lgt - 0x1.d63d5b2d15743p+28 * rt_inv);
  kTroeInf = 0x1.85580f28d3dcap+42 *
             exp(0x1.41175cb77b8d3p-1 * lgt - 0x1.9e41ee2d15743p+28 * rt_inv);
  k[r83b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM17]);
  k[r84f] = 0x1.dcd65p+34 * exp(-0x1.851f28p+24 * rt_inv);
  k[r84b] = 0x1.ad503fb9cb7bbp+27 *
            exp(0x1.b8f0984d353p-3 * lgt - 0x1.aa5da10a31653p+25 * rt_inv);
  k[r85f] = 0x1.dcd65p+32 * exp(-0x1.851f28p+24 * rt_inv);
  k[r85b] = 0x1.434d6a50da908p+29 *
            exp(0x1.e4b285c4ce22p-3 * lgt - 0x1.7476a8cb2f7dep+24 * rt_inv);
  k[r86f] = 0x1.83fffffffffffp+8 *
            exp(0x1.4p+1 * lgt - 0x1.8945800000001p+23 * rt_inv);
  k[r86b] = 0x1.99c984009294fp+0 *
            exp(0x1.571306e3e7e28p+1 * lgt - 0x1.1e215e67bdd2ep+25 * rt_inv);
  k[r87f] = 0x1.f4p+9 *
            exp(0x1.0cccccccccccdp+1 * lgt - 0x1.fb5f0cccccccep+20 * rt_inv);
  k[r87b] = 0x1.b330be789d618p+9 *
            exp(0x1.1f1c169c5d11cp+1 * lgt - 0x1.e5f5a4c98bef7p+25 * rt_inv);
  k[r88f] =
      0x1.bbcp+12 * exp(0x1.ccccccccccccdp+0 * lgt + 0x1.3067p+21 * rt_inv);
  k[r88b] = 0x1.007042e604224p+9 *
            exp(0x1.ebf322bcfa26p+0 * lgt - 0x1.599b4483df96cp+26 * rt_inv);
  k[r89f] = 0x1.31794b4p+34 * exp(-0x1.66514cp+27 * rt_inv);
  k[r89b] = 0x1.a20748701f2d8p+29 *
            exp(-0x1.2bd640d78f7ep-3 * lgt + 0x1.3651ffab6c444p+23 * rt_inv);
  k[r90f] = 0x1.3451eb851eb85p+3 *
            exp(0x1.7333333333333p+1 * lgt - 0x1.a27d480000001p+25 * rt_inv);
  k[r90b] = 0x1.60f390b5fae41p+5 *
            exp(0x1.53eb9673de82ap+1 * lgt - 0x1.d4293d940b48fp+24 * rt_inv);
  k[r91f] = 0x1.288879cp+35 * exp(-0x1.35a33p+26 * rt_inv);
  k[r91b] = 0x1.dfe0669eb061cp+37 *
            exp(-0x1.12444fc118a58p-1 * lgt - 0x1.58fb4bcbe0492p+25 * rt_inv);
  k[r92f] = 0x1.05532617c1bdap-5 *
            exp(0x1.95c28f5c28f5cp+1 * lgt - 0x1.c9e1600000001p+24 * rt_inv);
  k[r92b] = 0x1.387b6c4570b05p-4 *
            exp(0x1.8c2d027ce5a38p+1 * lgt - 0x1.04a6a1ace0456p+26 * rt_inv);
  k[r93f] = 0x1.1e1a3p+28 * exp(-0x1.0333ap+24 * rt_inv);
  k[r93b] = 0x1.7bea5939421b6p+24 *
            exp(-0x1.5e0f6bbcbcdp-6 * lgt - 0x1.71bc1ca499adap+25 * rt_inv);
  k[r94f] = 0x1.296d5b8p+32 *
            exp(0x1.999999999999ap-4 * lgt - 0x1.525dep+25 * rt_inv);
  k[r94b] = 0x1.466016bc3b8d2p+49 *
            exp(-0x1.1110498e64ff8p+0 * lgt - 0x1.afa907df4ff5bp+22 * rt_inv);
  k[r95f] = 0x1.338p+11 * exp(0x1p+1 * lgt - 0x1.07fd68p+25 * rt_inv);
  k[r95b] = 0x1.770908eb4ac2cp+11 *
            exp(0x1.c773d758f7ddp+0 * lgt - 0x1.b9290c381694bp+25 * rt_inv);
  k[r96f] = 0x1.dcd65p+33 * exp(0x1.231f8p+21 * rt_inv);
  k[r96b] = 0x1.a803ae4dc343fp+33 *
            exp(-0x1.21f89eea96b3p-2 * lgt - 0x1.bb941ef410373p+25 * rt_inv);
  k[r97f] =
      0x1.b58p+15 * exp(0x1.999999999999ap+0 * lgt - 0x1.5a072p+24 * rt_inv);
  k[r97b] = 0x1.5ab8d43dcffb7p+10 *
            exp(0x1.023b98f7ab5ccp+1 * lgt - 0x1.9db85175e833cp+25 * rt_inv);
  k[r98f] = 0x1.74ad942p+34;
  k[r98b] = 0x1.95178cc220d2cp+29 *
            exp(0x1.eb3e5da567c16p-2 * lgt + 0x1.5f2d8a308c1bap+22 * rt_inv);
  k[r99f] = 0x1.2a05f2p+35;
  k[r99b] = 0x1.40dd1d80b912dp+58 *
            exp(-0x1.3b2cc909b04cp+0 * lgt - 0x1.0701317a25311p+28 * rt_inv);
  k[r100f] = 0x1.65a0bcp+33 * exp(0x1.231f8p+21 * rt_inv);
  k[r100b] = 0x1.18ba4e0c9da66p+56 *
             exp(-0x1.4b1ec81d4d7cp+0 * lgt - 0x1.284e40d1a466fp+28 * rt_inv);
  k[r101f] = 0x1.dcd65p+33;
  k[r101b] = 0x1.7bce50d3e1e21p+12 *
             exp(0x1.6613287dbd494p+0 * lgt - 0x1.263032d214b0fp+25 * rt_inv);
  kTroe0 = 0x1.289c98651e77bp+107 *
           exp(-0x1.970a3d70a3d71p+2 * lgt - 0x1.41c48p+24 * rt_inv);
  kTroeInf = 0x1.b6605ec82p+48 *
             exp(-0x1.28f5c28f5c28fp+0 * lgt - 0x1.24666p+22 * rt_inv);
  fcTroe = 0x1.96d5cfaacd9e8p-2 * exp(-temp / 0x1.ap+7) +
           0x1.3495182a9930cp-1 * exp(-temp / 0x1.ea4p+11) +
           0x1p+0 * exp(-0x1.3e2p+13 / temp);
  k[r102f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM18]);
  kTroe0 = 0x1.f482d8c88e021p+139 *
           exp(-0x1.e961de887ffd8p+2 * lgt - 0x1.8d651cd2d2dc1p+28 * rt_inv);
  kTroeInf = 0x1.71dd50f3103aep+81 *
             exp(-0x1.392a237766616p+1 * lgt - 0x1.7dda6e52d2dc2p+28 * rt_inv);
  k[r102b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM18]);
  k[r103f] =
      0x1.c138p+16 * exp(0x1.e666666666666p+0 * lgt - 0x1.e0bc7p+24 * rt_inv);
  k[r103b] = 0x1.8db5106eb3991p+4 *
             exp(0x1.383f3f1c9b64p+1 * lgt - 0x1.489d79a32ac04p+25 * rt_inv);
  k[r104f] =
      0x1.5ec8p+16 * exp(0x1.eb851eb851eb8p+0 * lgt - 0x1.6b43fp+24 * rt_inv);
  k[r104b] = 0x1.6c4de3b6c1e9dp+3 *
             exp(0x1.365298a4a5c21p+1 * lgt - 0x1.c3c656016e5b1p+24 * rt_inv);
  k[r105f] = 0x1.ba8p+11 *
             exp(0x1.0f5c28f5c28f6p+1 * lgt - 0x1.bc58800000001p+21 * rt_inv);
  k[r105b] = 0x1.f6e8b4573382ap+2 *
             exp(0x1.486c56526e1ap+1 * lgt - 0x1.2939dad05c43bp+26 * rt_inv);
  k[r106f] = 0x1.2a05f2p+35 * exp(-0x1.96331c0000001p+27 * rt_inv);
  k[r106b] = 0x1.91092e87773d5p+25 *
             exp(0x1.6bf9e571c0a44p-3 * lgt + 0x1.d9c0b51e1b248p+21 * rt_inv);
  k[r107f] = 0x1.186158p+28 * exp(-0x1.dce7dp+25 * rt_inv);
  k[r107b] = 0x1.be2cbdf44ffaep+25 *
             exp(-0x1.b14118bb17a4p-3 * lgt - 0x1.d86f8326ccee5p+22 * rt_inv);
  k[r108f] =
      0x1.7fcp+12 * exp(0x1.bd70a3d70a3d7p+0 * lgt - 0x1.4d9418p+25 * rt_inv);
  k[r108b] = 0x1.c339df47504a5p+8 *
             exp(0x1.fd3f8ee1ad6fp+0 * lgt - 0x1.e261dff2b9e6ep+25 * rt_inv);
  kTroe0 = 0x1.329ba8fa506b9p+117 *
           exp(-0x1.c51eb851eb852p+2 * lgt - 0x1.aac9f8p+24 * rt_inv);
  kTroeInf = 0x1.d9d8c3ed9p+48 *
             exp(-0x1.fae147ae147aep-1 * lgt - 0x1.937c8p+22 * rt_inv);
  fcTroe = 0x1.432ca57a786c4p-3 * exp(-temp / 0x1.f4p+6) +
           0x1.af34d6a161e4fp-1 * exp(-temp / 0x1.156p+11) +
           0x1p+0 * exp(-0x1.ae2p+12 / temp);
  k[r109f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM19]);
  kTroe0 = 0x1.3b9902338262p+137 *
           exp(-0x1.d1520110e5f5p+2 * lgt - 0x1.aca90a856312cp+28 * rt_inv);
  kTroeInf = 0x1.e7bd6121d6dc5p+68 *
             exp(-0x1.2e3dc6d2f3fcfp+0 * lgt - 0x1.984a5d056312dp+28 * rt_inv);
  k[r109b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM19]);
  k[r110f] = 0x1.dcd65p+30;
  k[r110b] = 0x1.20d38a23170ccp+28 *
             exp(0x1.c949be361318p-2 * lgt - 0x1.0cd70409b33d6p+28 * rt_inv);
  k[r111f] = 0x1.ebbd028p+36;
  k[r111b] = 0x1.0498bb86513ebp+27 *
             exp(0x1.be7e5e42474cp-1 * lgt - 0x1.3dfd8812de3bbp+28 * rt_inv);
  k[r112f] = 0x1.312dp+24;
  k[r112b] = 0x1.18f5b5a2c141cp+24 *
             exp(0x1.5b99468ea27cp-4 * lgt - 0x1.981063588d3f4p+25 * rt_inv);
  k[r113f] = 0x1.4dc938p+30;
  k[r113b] = 0x1.c8bac6188f58fp+39 *
             exp(-0x1.7c5a8454b52p-4 * lgt - 0x1.01cf1bd54de78p+28 * rt_inv);
  k[r114f] = 0x1.bf08ebp+36;
  k[r114b] = 0x1.79f08d0939594p+44 *
             exp(0x1.56af40907c4p-5 * lgt - 0x1.54c2f5a67019p+28 * rt_inv);
  k[r115f] = 0x1.2ac4ae2p+36;
  k[r115b] = 0x1.70afce465a4e4p+46 *
             exp(-0x1.12ce2d13fc58p-1 * lgt - 0x1.2ea70b2d1a4d8p+28 * rt_inv);
  kTroe0 = 0x1.ea7457a7c18ap+158 *
           exp(-0x1.29eb851eb851fp+3 * lgt - 0x1.8e756e0000001p+28 * rt_inv);
  kTroeInf = 0x1.d1a94a2p+42 *
             exp(0x1.c28f5c28f5c29p-2 * lgt - 0x1.62352bp+28 * rt_inv);
  fcTroe = 0x1.0fdf3b645a1cap-2 * exp(-temp / 0x1.68p+7) +
           0x1.7810624dd2f1bp-1 * exp(-temp / 0x1.02cp+10) +
           0x1p+0 * exp(-0x1.529p+12 / temp);
  k[r116f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM20]);
  kTroe0 = 0x1.f4590ce08dd9p+139 *
           exp(-0x1.1feb1bca19a4p+3 * lgt - 0x1.c52ec75e9a032p+27 * rt_inv);
  kTroeInf = 0x1.db0df726fd4c8p+23 *
             exp(0x1.814e435e65c04p-1 * lgt - 0x1.6cae415e9a031p+27 * rt_inv);
  k[r116b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM20]);
  kTroe0 = 0x1.ce3922c2af443p+119 *
           exp(-0x1.e7ae147ae147bp+2 * lgt - 0x1.bcfbf00000001p+24 * rt_inv);
  kTroeInf =
      0x1.017df8p+30 * exp(0x1.d0e5604189375p-2 * lgt - 0x1.d0c68p+22 * rt_inv);
  fcTroe = 0x1.94af4f0d844ep-6 * exp(-temp / 0x1.a4p+7) +
           0x1.f35a858793dd9p-1 * exp(-temp / 0x1.ecp+9) +
           0x1p+0 * exp(-0x1.116p+12 / temp);
  k[r117f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM21]);
  kTroe0 = 0x1.5bb47d37e20efp+130 *
           exp(-0x1.edeff32888e5cp+2 * lgt - 0x1.57fa1c602a578p+27 * rt_inv);
  kTroeInf = 0x1.8364e9700a352p+40 *
             exp(0x1.6cc775670f565p-2 * lgt - 0x1.2ee0d2602a578p+27 * rt_inv);
  k[r117b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM21]);
  kTroe0 = 0x1.287626ee52198p+80 *
           exp(-0x1.ee147ae147ae1p+1 * lgt - 0x1.a7ea800000001p+23 * rt_inv);
  kTroeInf = 0x1.6a657p+32 *
             exp(0x1.147ae147ae148p-2 * lgt - 0x1.1e04000000001p+20 * rt_inv);
  fcTroe = 0x1.be76c8b43958p-3 * exp(-temp / 0x1.9fp+7) +
           0x1.90624dd2f1aap-1 * exp(-temp / 0x1.4cep+11) +
           0x1p+0 * exp(-0x1.7cfp+12 / temp);
  k[r118f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM22]);
  kTroe0 = 0x1.361774b5222f9p+99 *
           exp(-0x1.ff70ea53cfeap+1 * lgt - 0x1.ca2526d3e5124p+28 * rt_inv);
  kTroeInf = 0x1.7b0ec835b15fcp+51 *
             exp(0x1.132ecb66d86ap-3 * lgt - 0x1.be03d6d3e5125p+28 * rt_inv);
  k[r118b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM22]);
  k[r119f] =
      0x1.4b4p+10 * exp(0x1.43d70a3d70a3dp+1 * lgt - 0x1.86b7cp+25 * rt_inv);
  k[r119b] = 0x1.20988bdbb8b9p-1 *
             exp(0x1.81d8f41b6c426p+1 * lgt - 0x1.0f57865e35827p+24 * rt_inv);
  k[r120f] = 0x1.c1fffffffffffp+10 * exp(0x1p+1 * lgt - 0x1.3f36cp+23 * rt_inv);
  k[r120b] = 0x1.f749ebc6f5b25p+2 *
             exp(0x1.32060b513ee0ep+1 * lgt - 0x1.2f309f2ca8864p+25 * rt_inv);
  k[r121f] = 0x1.c5fffffffffffp+7 * exp(0x1p+1 * lgt - 0x1.25ad4p+25 * rt_inv);
  k[r121b] = 0x1.06a9eed34a1bap+5 *
             exp(0x1.18dd5379e4de6p+1 * lgt - 0x1.8cbf25faa7824p+23 * rt_inv);
  k[r122f] = 0x1.2cp+14 *
             exp(0x1.d47ae147ae148p+0 * lgt - 0x1.c174000000001p+19 * rt_inv);
  k[r122b] = 0x1.9cde058e7af18p-3 *
             exp(0x1.5b8a76b1ceb44p+1 * lgt - 0x1.98385344c1de5p+26 * rt_inv);
  k[r123f] = 0x1.2a05f2p+32;
  k[r123b] = 0x1.983d7432c8261p+35 *
             exp(0x1.552aea66672p-4 * lgt - 0x1.4cd5a602e3bd8p+28 * rt_inv);
  k[r124f] =
      0x1.d7ep+13 * exp(0x1.e8f5c28f5c28fp+0 * lgt - 0x1.dd8b4p+23 * rt_inv);
  k[r124b] = 0x1.e2432de95e3bbp+1 *
             exp(0x1.2e00c884be243p+1 * lgt + 0x1.674e76e6b1b5dp+24 * rt_inv);
  k[r125f] = 0x1.3a0abebp+35 * exp(-0x1.cbabp+27 * rt_inv);
  k[r125b] = 0x1.9fdcb4ff2b729p+26 *
             exp(0x1.f6af8975f058p-4 * lgt + 0x1.3e0f738be35b4p+24 * rt_inv);
  k[r126f] = 0x1.671e344p+36;
  k[r126b] = 0x1.7f34b63e64687p+36 *
             exp(0x1.6a535dfefb8p-3 * lgt - 0x1.1107c8833203p+28 * rt_inv);
  k[r127f] = 0x1.bf08ebp+35;
  k[r127b] = 0x1.d0a3200fa62b1p+47 *
             exp(-0x1.33ca466b1d54p-1 * lgt - 0x1.6990277bb881dp+28 * rt_inv);
  k[r128f] = 0x1.71434p+23 * exp(0x1.3067p+21 * rt_inv);
  k[r128b] = 0x1.d79b1e4c1b588p+24 *
             exp(0x1.10f19cb8f3ecp-2 * lgt - 0x1.730f57079b411p+26 * rt_inv);
  k[r129f] = 0x1.73eed8p+28;
  k[r129b] = 0x1.078ddcafd9d15p+37 *
             exp(-0x1.cfec1084ddcp-4 * lgt - 0x1.1899994d23e72p+28 * rt_inv);
  k[r130f] = 0x1.c9c38p+29;
  k[r130b] = 0x1.1854e67cb6059p+41 *
             exp(-0x1.3ae59ff05ccp-2 * lgt - 0x1.30e9481d4ead5p+28 * rt_inv);
  k[r131f] = 0x1.4d3d25d88p+45 *
             exp(-0x1.63d70a3d70a3dp+0 * lgt - 0x1.0333ap+22 * rt_inv);
  k[r131b] = 0x1.09abb4ad10133p+39 *
             exp(-0x1.a0b43ca39398p-1 * lgt - 0x1.65c0c4b4244e8p+28 * rt_inv);
  k[r132f] = 0x1.4e4p+10 *
             exp(0x1.9c28f5c28f5c3p+0 * lgt + 0x1.8840000000001p+20 * rt_inv);
  k[r132b] = 0x1.0f0f1d4dd8302p+13 *
             exp(0x1.6cda865dd444p+0 * lgt - 0x1.ad54872483667p+25 * rt_inv);
  k[r133f] =
      0x1.7d784p+26 * exp(0x1.28f5c28f5c28fp-2 * lgt - 0x1.679p+15 * rt_inv);
  k[r133b] = 0x1.6e3e3f70ba3d6p+32 *
             exp(-0x1.3373c083ef61p-3 * lgt - 0x1.93ab4639c1ab7p+24 * rt_inv);
  kTroe0 = 0x1.d462db8b40f69p+114 *
           exp(-0x1.d147ae147ae14p+2 * lgt - 0x1.ccf1e00000001p+24 * rt_inv);
  kTroeInf = 0x1.4dc938p+32 * exp(-0x1.3272p+23 * rt_inv);
  fcTroe = 0x1.fe90ff9724744p-3 * exp(-temp / 0x1.8ap+6) +
           0x1.805bc01a36e2fp-1 * exp(-temp / 0x1.458p+10) +
           0x1p+0 * exp(-0x1.047p+12 / temp);
  k[r134f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM23]);
  kTroe0 = 0x1.90012b739fc5p+122 *
           exp(-0x1.c6478bceb92e4p+2 * lgt - 0x1.5197516d2ccd5p+27 * rt_inv);
  kTroeInf = 0x1.1d0e24236fed8p+40 *
             exp(0x1.600448b8366p-3 * lgt - 0x1.2b20356d2ccd4p+27 * rt_inv);
  k[r134b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM23]);
  k[r135f] = 0x1.fep+13 * exp(0x1p+1 * lgt - 0x1.e5348p+22 * rt_inv);
  k[r135b] = 0x1.8a7348a9118d8p+15 *
             exp(0x1.d15282ece4908p+0 * lgt - 0x1.569a3a5da0372p+26 * rt_inv);
  k[r136f] = 0x1.46d9833756p+55 *
             exp(-0x1.68f5c28f5c28fp+0 * lgt - 0x1.ce0fe4p+26 * rt_inv);
  k[r136b] = 0x1.0fc94b2f64e36p+43 *
             exp(-0x1.b32a0a8c872e4p-1 * lgt + 0x1.ea57c32bdbd8fp+22 * rt_inv);
  k[r137f] = 0x1.fep+11 * exp(0x1p+1 * lgt - 0x1.e5348p+22 * rt_inv);
  k[r137b] = 0x1.5cdbdc0c9cc99p-7 *
             exp(0x1.a87759f35187p+1 * lgt - 0x1.7b5ef8e3f7b37p+27 * rt_inv);
  k[r138f] =
      0x1.d426c47622575p-23 * exp(0x1.2p+2 * lgt + 0x1.febep+21 * rt_inv);
  k[r138b] = 0x1.84b96e3aab8a1p-10 *
             exp(0x1.e0e43b26343ap+1 * lgt - 0x1.682a933353c5p+26 * rt_inv);
  k[r139f] =
      0x1.f8p+8 * exp(0x1.2666666666666p+1 * lgt - 0x1.aef05p+25 * rt_inv);
  k[r139b] = 0x1.5f9493f1c822ap+26 *
             exp(0x1.255d913ada3fcp+0 * lgt - 0x1.bef62b1227908p+24 * rt_inv);
  k[r140f] = 0x1.0747fffffffffp+15 * exp(0x1p+1 * lgt - 0x1.bee64p+25 * rt_inv);
  k[r140b] = 0x1.df28289f43151p+6 *
             exp(0x1.403082b8bac42p+1 * lgt + 0x1.01785c57a2e4fp+21 * rt_inv);
  k[r141f] = 0x1.034f03b80a36dp-21 * exp(0x1p+2 * lgt + 0x1.febep+22 * rt_inv);
  k[r141b] = 0x1.e9dc03a708608p-31 *
             exp(0x1.35c458ee50ebcp+2 * lgt - 0x1.a2ad732e7e004p+27 * rt_inv);
  k[r142f] = 0x1.700acp+22 * exp(-0x1.fb2abffffffffp+24 * rt_inv);
  k[r142b] = 0x1.e03e4fffb38cdp+22 *
             exp(0x1.30ab9ff7e38p-6 * lgt - 0x1.0ec4ec595558cp+28 * rt_inv);
  k[r143f] = 0x1.74876e8p+36;
  k[r143b] = 0x1.c3e4306c22af3p+16 *
             exp(0x1.8f8e300d5b988p+0 * lgt - 0x1.11d67e0c5257ep+26 * rt_inv);
  k[r144f] = 0x1.74876e8p+36;
  k[r144b] = 0x1.90f76b15ae65fp+10 *
             exp(0x1.628c110c6bcp+0 * lgt - 0x1.96fe04c6f67a3p+28 * rt_inv);
  k[r145f] = 0x1.7d784p+30 * exp(-0x1.b42c800000001p+21 * rt_inv);
  k[r145b] = 0x1.a203dac94129cp-5 *
             exp(0x1.d1e5735d1e6p+0 * lgt - 0x1.571b478923bdcp+28 * rt_inv);
  k[r146f] = 0x1.2a05f2p+33;
  k[r146b] = 0x1.6d5549a85a3e2p+12 *
             exp(0x1.ad6d08eddce8p-1 * lgt - 0x1.e59746037be98p+25 * rt_inv);
  k[r147f] = 0x1.bf08ebp+34;
  k[r147b] = 0x1.2c5be875ded0bp+35 *
             exp(0x1.ecef31e67a1p-3 * lgt - 0x1.6cdfe50b9106p+28 * rt_inv);
  k[r148f] = 0x1.2a05f2p+33;
  k[r148b] = 0x1.9ebd6fa91590ap+5 *
             exp(0x1.91398e1f8726p+0 * lgt - 0x1.48ec482f8e6d6p+28 * rt_inv);
  k[r149f] = 0x1.74876e8p+35;
  k[r149b] = 0x1.62b8d0541f4b3p+38 *
             exp(0x1.11bd9fc03998p-2 * lgt - 0x1.6f0a1f54b8fb8p+28 * rt_inv);
  kTroe0 = 0x1.83bdac6ae9bc2p+91 *
           exp(-0x1.3333333333333p+2 * lgt - 0x1.e5348p+22 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+46 * exp(-0x1p+0 * lgt);
  fcTroe = 0x1.6a161e4f765fep-2 * exp(-temp / 0x1.08p+7) +
           0x1.4af4f0d844d01p-1 * exp(-temp / 0x1.48cp+10) +
           0x1p+0 * exp(-0x1.5bep+12 / temp);
  k[r150f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM24]);
  kTroe0 = 0x1.f278acf61d604p+110 *
           exp(-0x1.42f6a6a0352dp+2 * lgt - 0x1.0d235db914b5fp+29 * rt_inv);
  kTroeInf = 0x1.d3b0d47b217dap+65 *
             exp(-0x1.3f0dcdb407e74p+0 * lgt - 0x1.0958f4b914b5fp+29 * rt_inv);
  k[r150b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM24]);
  k[r151f] = 0x1.2a05f2p+34;
  k[r151b] = 0x1.15331e1c898b1p+48 *
             exp(-0x1.7c1c74b8653cp-1 * lgt - 0x1.927e29482f044p+27 * rt_inv);
  k[r152f] = 0x1.74876e8p+35 * exp(-0x1.7f0e800000001p+22 * rt_inv);
  k[r152b] = 0x1.ce47501a2d189p+26 *
             exp(0x1.b0a93800fbfp-1 * lgt - 0x1.3013d946171bap+29 * rt_inv);
  k[r153f] = 0x1.eap+8 * exp(0x1.4p+1 * lgt - 0x1.1e04000000001p+21 * rt_inv);
  k[r153b] = 0x1.59a0bc45b10c9p+20 *
             exp(0x1.e7a73d751238p+0 * lgt - 0x1.dfa678e1840ccp+26 * rt_inv);
  kTroe0 = 0x1.5af1d78b58c4p+71 *
           exp(-0x1.91eb851eb851fp+1 * lgt - 0x1.3a1b400000001p+22 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+44 * exp(-0x1.999999999999ap-1 * lgt);
  fcTroe = 0x1.47ae147ae147ap-2 * exp(-temp / 0x1.38p+6) +
           0x1.5c28f5c28f5c3p-1 * exp(-temp / 0x1.f2cp+10) +
           0x1p+0 * exp(-0x1.5d6p+12 / temp);
  k[r154f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM25]);
  kTroe0 = 0x1.00111bed5dd25p+88 *
           exp(-0x1.a6b0b56ae04fp+1 * lgt - 0x1.bfa4058abd1e1p+28 * rt_inv);
  kTroeInf = 0x1.0c816abb210e5p+61 *
             exp(-0x1.ecae5aca398dep-1 * lgt - 0x1.babb988abd1e1p+28 * rt_inv);
  k[r154b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM25]);
  k[r155f] = 0x1.2a05f2p+36;
  k[r155b] = 0x1.b994980f0eb34p+42 *
             exp(-0x1.624af3b703ep-2 * lgt - 0x1.6c2e8c4b55a85p+28 * rt_inv);
  k[r156f] = 0x1.2a05f2p+34;
  k[r156b] = 0x1.dea04d0647307p+50 *
             exp(-0x1.8a900868753p-1 * lgt - 0x1.3da777579ea04p+28 * rt_inv);
  k[r157f] = 0x1.f4p+8 * exp(0x1p+1 * lgt - 0x1.cd95500000001p+24 * rt_inv);
  k[r157b] = 0x1.94f492398fb1bp+17 *
             exp(0x1.7d2aaa90ca7cp+0 * lgt - 0x1.d484d287a5b5bp+25 * rt_inv);
  k[r158f] = 0x1.896402p+33 * exp(-0x1.7f0e800000001p+22 * rt_inv);
  k[r158b] = 0x1.28b6b84d73fffp+31 *
             exp(0x1.6c6a562f1fp-4 * lgt - 0x1.2edfb00d82ea1p+28 * rt_inv);
  k[r159f] = 0x1.2a05f2p+34;
  k[r159b] = 0x1.78135831ed8f8p+38 *
             exp(-0x1.26e60061fap-7 * lgt - 0x1.ceb196e41f064p+28 * rt_inv);
  kTroe0 = 0x1.1623b50660ab2p+91 *
           exp(-0x1.470a3d70a3d71p+2 * lgt - 0x1.c4f6e80000001p+24 * rt_inv);
  kTroeInf = 0x1.823cf4p+29 * exp(0x1p-1 * lgt - 0x1.1fee5p+24 * rt_inv);
  fcTroe = 0x1.a31f8a0902dep-2 * exp(-temp / 0x1.13p+8) +
           0x1.2e703afb7e91p-1 * exp(-temp / 0x1.328p+10) +
           0x1p+0 * exp(-0x1.441p+12 / temp);
  k[r160f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM26]);
  kTroe0 = 0x1.68ecb9ae7f7fcp+129 *
           exp(-0x1.b6bf10f1ee8c8p+2 * lgt - 0x1.5bc17640530efp+28 * rt_inv);
  kTroeInf = 0x1.f5329d96818c1p+67 *
             exp(-0x1.3ed34e052ad5cp+0 * lgt - 0x1.5170ecc0530fp+28 * rt_inv);
  k[r160b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM26]);
  k[r161f] = 0x1.dcd65p+34;
  k[r161b] = 0x1.828d1c44d8489p+55 *
             exp(-0x1.14b3defd1014p+0 * lgt - 0x1.0aef5fda179f7p+29 * rt_inv);
  k[r162f] = 0x1.0c388dp+33 * exp(-0x1.3272p+21 * rt_inv);
  k[r162b] = 0x1.871dcffe514dep+32 *
             exp(-0x1.fe3fe273a31p-5 * lgt - 0x1.2fc192bbf9b2ap+25 * rt_inv);
  k[r163f] = 0x1.bf08ebp+34;
  k[r163b] = 0x1.45ee2d53ee1f1p+34 *
             exp(-0x1.fe3fe273a478p-5 * lgt - 0x1.1c9a72bbf9b36p+25 * rt_inv);
  k[r164f] = 0x1.0c388dp+33;
  k[r164b] = 0x1.871dcffe50f85p+32 *
             exp(-0x1.fe3fe273a49p-5 * lgt - 0x1.1c9a72bbf9b33p+25 * rt_inv);
  k[r165f] = 0x1.a13b86p+32;
  k[r165b] = 0x1.3033a1c5cd53fp+32 *
             exp(-0x1.fe3fe273a508p-5 * lgt - 0x1.1c9a72bbf9b3bp+25 * rt_inv);
  k[r166f] = 0x1.0c388dp+33 * exp(-0x1.3272p+21 * rt_inv);
  k[r166b] = 0x1.871dcffe519bbp+32 *
             exp(-0x1.fe3fe273a3ep-5 * lgt - 0x1.2fc192bbf9b42p+25 * rt_inv);
  k[r167f] = 0x1.bf08ebp+33;
  k[r167b] = 0x1.69781f033b037p+35 *
             exp(0x1.6246aeafb2ap-3 * lgt - 0x1.77c65c3ed5255p+29 * rt_inv);
  k[r168f] = 0x1.bf08ebp+33;
  k[r168b] = 0x1.e2ee304863881p+39 *
             exp(-0x1.a212f0057818p-2 * lgt - 0x1.8fc1daa2d4de5p+28 * rt_inv);
  k[r169f] = 0x1.bf08ebp+34;
  k[r169b] = 0x1.05b8faa349a8ep+51 *
             exp(-0x1.aa74068faf5p-1 * lgt - 0x1.613ac5af1dd61p+28 * rt_inv);
  k[r170f] = 0x1.04c533cp+36;
  k[r170b] = 0x1.33f862bafae8ep+44 *
             exp(-0x1.258ea905a5f1p-1 * lgt - 0x1.052a4ea1cfad8p+26 * rt_inv);
  k[r171f] = 0x1.a13b86p+34;
  k[r171b] = 0x1.78e89b67a31ep+19 *
             exp(0x1.095d0d3f0af8p-2 * lgt - 0x1.0f3d4f060f277p+28 * rt_inv);
  k[r172f] = 0x1.65a0bcp+33;
  k[r172b] = 0x1.79ef0a93443cbp+29 *
             exp(0x1.0754f61a5fbp-1 * lgt - 0x1.7407bfdfc4a45p+29 * rt_inv);
  k[r173f] = 0x1.a13b86p+33;
  k[r173b] = 0x1.6e80de3029a8dp+27 *
             exp(0x1.c8fdbc755228p-2 * lgt - 0x1.f5795268f12c3p+27 * rt_inv);
  k[r174f] = 0x1.8de76816d8p+52 * exp(-0x1.45e4b7p+28 * rt_inv);
  k[r174b] = 0x1.26f45fd9121f2p+15 *
             exp(0x1.85089d5c909dp+0 * lgt + 0x1.3c2488bfef25ep+24 * rt_inv);
  kTroe0 = 0x1.176592ep+40 * exp(-0x1.8f9788p+25 * rt_inv);
  kTroeInf = 0x1.5d3ef798p+41 * exp(-0x1.0ae50cp+26 * rt_inv);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r175f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM28]);
  kTroe0 = 0x1.098f9af963235p+12 *
           exp(0x1.271ce153d850dp+0 * lgt - 0x1.28b43d6421404p+21 * rt_inv);
  kTroeInf = 0x1.4bf381b7bbec2p+13 *
             exp(0x1.271ce153d850dp+0 * lgt - 0x1.317ba7ac8428p+24 * rt_inv);
  k[r175b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM28]);
  k[r176f] = 0x1.91bc3dp+31 * exp(0x1.3ca9p+21 * rt_inv);
  k[r176b] = 0x1.2d0b801b4c3cep+21 *
             exp(0x1.b583218d1f28p-1 * lgt - 0x1.bc67c94effa2fp+26 * rt_inv);
  k[r177f] = 0x1.416364p+28 * exp(0x1.3ca9p+21 * rt_inv);
  k[r177b] = 0x1.865a0a084babp+21 *
             exp(0x1.db2131cc8904p-2 * lgt - 0x1.3a77443d39974p+26 * rt_inv);
  k[r178f] = 0x1.92738f5028p+50 *
             exp(-0x1.e666666666666p+0 * lgt - 0x1.7bdd5p+23 * rt_inv);
  k[r178b] = 0x1.1395e8cfab88ep+36 *
             exp(-0x1.f94a3b9068b5p-1 * lgt - 0x1.e1247bfdfde85p+25 * rt_inv);
  k[r179f] = 0x1.15295e8p+35 *
             exp(-0x1.999999999999ap-3 * lgt - 0x1.c60ccp+23 * rt_inv);
  k[r179b] = 0x1.339e0e1a31a11p+24 *
             exp(0x1.4a53445e4634p-2 * lgt - 0x1.df9e9bb4e3a47p+24 * rt_inv);
  k[r180f] = 0x1.5b32724p+35 *
             exp(-0x1.6666666666666p-2 * lgt - 0x1.7d864p+23 * rt_inv);
  k[r180b] = 0x1.9558cf550ce4ap+21 *
             exp(0x1.323f688cde3ep-1 * lgt - 0x1.06c6635038bc7p+26 * rt_inv);
  k[r181f] = 0x1.b9130ap+30 *
             exp(0x1.999999999999ap-2 * lgt - 0x1.562228p+24 * rt_inv);
  k[r181b] = 0x1.a1517e2afb8a5p+20 *
             exp(0x1.ea4cdfe603bdp-1 * lgt - 0x1.555b407ce5683p+25 * rt_inv);
  k[r182f] = 0x1.accf3dacbf296p-32 *
             exp(0x1.7333333333333p+2 * lgt - 0x1.18e8800000001p+23 * rt_inv);
  k[r182b] = 0x1.4c7303f124aefp-37 *
             exp(0x1.9d4f3b792a068p+2 * lgt - 0x1.1879ee78004eep+26 * rt_inv);
  k[r183f] = 0x1.916872b020c4ap-6 *
             exp(0x1.9333333333333p+1 * lgt - 0x1.6da0a80000001p+24 * rt_inv);
  k[r183b] = 0x1.f86864bc7d514p-8 *
             exp(0x1.b56ea1956a4b2p+1 * lgt - 0x1.9da906cc748a4p+25 * rt_inv);
  k[r184f] = 0x1.ff973cafa8p+54 *
             exp(-0x1.199999999999ap+1 * lgt - 0x1.bfdb680000001p+25 * rt_inv);
  k[r184b] = 0x1.0c63b44c2d746p+51 *
             exp(-0x1.005d101eaec8ep+1 * lgt - 0x1.73ed7d6220591p+25 * rt_inv);
  k[r185f] = 0x1.ba814p+27 *
             exp(0x1.999999999999ap-2 * lgt - 0x1.da7ac00000001p+25 * rt_inv);
  k[r185b] = 0x1.78443e1147b3cp+27 *
             exp(0x1.a731a8467613p-3 * lgt - 0x1.1557967d2885cp+24 * rt_inv);
  k[r186f] = 0x1.74876e8p+36 * exp(-0x1.50c5480000001p+27 * rt_inv);
  k[r186b] = 0x1.4a8732448af0cp+25 *
             exp(0x1.2c40e576e0639p-1 * lgt - 0x1.14c4cd5b28168p+22 * rt_inv);
  k[r187f] = 0x1.a13b86p+30 * exp(-0x1.58c0400000001p+22 * rt_inv);
  k[r187b] = 0x1.ca2d09b143171p+21 *
             exp(0x1.8f59966a1c9ap-1 * lgt - 0x1.87e74e3f98516p+27 * rt_inv);
  k[r188f] =
      0x1.b199999999999p+4 * exp(0x1.6p+1 * lgt - 0x1.6cab8p+21 * rt_inv);
  k[r188b] = 0x1.ed3d76dc2e00ap-18 *
             exp(0x1.15523b5b36dp+2 * lgt - 0x1.f88b9f29a83efp+26 * rt_inv);
  k[r189f] = 0x1.74876e8p+37 * exp(-0x1.febep+24 * rt_inv);
  k[r189b] = 0x1.27cdefb7ae12ap+27 *
             exp(0x1.310423c4a5ddp-1 * lgt - 0x1.5753d9ee1904fp+24 * rt_inv);
  k[r190f] = 0x1.2a05f2p+33 * exp(-0x1.febep+24 * rt_inv);
  k[r190b] = 0x1.15992d6e0d584p+22 *
             exp(0x1.1f141940f89acp-1 * lgt - 0x1.feaf795263bbbp+23 * rt_inv);
  k[r191f] = 0x1.2a05f2p+33 * exp(-0x1.febep+22 * rt_inv);
  k[r191b] = 0x1.2fc8bf6a30a6ep+26 *
             exp(0x1.0114a991b347p-1 * lgt - 0x1.ca9198f49a4d3p+25 * rt_inv);
  k[r192f] = 0x1.bca691p+31 * exp(0x1.02b0ep+22 * rt_inv);
  k[r192b] = 0x1.ca4d5567b449dp+22 *
             exp(0x1.491a3d53098fp-1 * lgt - 0x1.ac01696fe1f35p+26 * rt_inv);
  k[r193f] = 0x1.74876e8p+35;
  k[r193b] = 0x1.09f90e18e8ab4p+16 *
             exp(0x1.4f9388a6d7004p+0 * lgt - 0x1.70a070f165df4p+25 * rt_inv);
  k[r194f] = 0x1.2a05f2p+34;
  k[r194b] = 0x1.1802cc5674e39p+32 *
             exp(0x1.3c60b4453f28p-2 * lgt - 0x1.12955b2e981afp+28 * rt_inv);
  k[r195f] = 0x1.74876e8p+36;
  k[r195b] = 0x1.3557869ff4a1cp+24 *
             exp(0x1.0428d49ffbb8p+0 * lgt - 0x1.48847e508833dp+28 * rt_inv);
  k[r196f] = 0x1.bf08ebp+34;
  k[r196b] = 0x1.0d9781abd23a9p+36 *
             exp(0x1.b9037fbeb5bp-3 * lgt - 0x1.4e6338ae49d61p+28 * rt_inv);
  k[r197f] = 0x1.c9c38p+24;
  k[r197b] = 0x1.3ddbef4bf7648p-9 *
             exp(0x1.af04d3f11c5ap+0 * lgt - 0x1.8fffe56784fc6p+27 * rt_inv);
  k[r198f] = 0x1.c8591a9p+38 * exp(-0x1p-1 * lgt);
  k[r198b] = 0x1.25a84b7f04427p+24 *
             exp(-0x1.f93bd80ba268p-4 * lgt + 0x1.9593f90a47f77p+25 * rt_inv);
  k[r199f] = 0x1.a13b86p+32;
  k[r199b] = 0x1.fd50e3718f111p+0 *
             exp(0x1.76d01d3e62ff8p+0 * lgt - 0x1.082831de3a93ap+26 * rt_inv);
  k[r200f] = 0x1.65a0bcp+31;
  k[r200b] = 0x1.f15e55bfe25e3p+38 *
             exp(-0x1.9139734005fp-3 * lgt - 0x1.8917641e09d56p+27 * rt_inv);
  k[r201f] = 0x1.0c9e6b65ddbadp+143 *
             exp(-0x1.3a8f5c28f5c29p+3 * lgt - 0x1.5d30240000001p+27 * rt_inv);
  k[r201b] = 0x1.3b0bdf4a919bap+111 *
             exp(-0x1.092c97740cd82p+3 * lgt - 0x1.3ee5557e7394dp+27 * rt_inv);
  k[r202f] = 0x1.9f8e3d1eae8a1p+143 *
             exp(-0x1.33851eb851eb8p+3 * lgt - 0x1.6e0ae40000001p+27 * rt_inv);
  k[r202b] = 0x1.ac72aefb0be67p+133 *
             exp(-0x1.34cb77b104671p+3 * lgt - 0x1.64b3cfa67ddefp+25 * rt_inv);
  k[r203f] = 0x1.af10410a711a5p+169 *
             exp(-0x1.52e147ae147aep+3 * lgt - 0x1.927c1b8p+28 * rt_inv);
  k[r203b] = 0x1.3f995851b060ep+130 *
             exp(-0x1.1ea4edb24d5bep+3 * lgt - 0x1.847c6c5c2353p+25 * rt_inv);
  k[r204f] = 0x1.232ae63c59c6cp+86 *
             exp(-0x1.d70a3d70a3d71p+1 * lgt - 0x1.1a80128p+28 * rt_inv);
  k[r204b] = 0x1.9694c570f9ab2p+51 *
             exp(-0x1.ebc50cf3ad22p+0 * lgt - 0x1.d3d38245b0593p+27 * rt_inv);
  k[r205f] =
      0x1.593ae8p+27 * exp(0x1.8f5c28f5c28f6p-2 * lgt - 0x1.6df26p+21 * rt_inv);
  k[r205b] = 0x1.2836a690844dep+11 *
             exp(0x1.e6d1bc27809dp+0 * lgt - 0x1.d1c1286279356p+25 * rt_inv);
  k[r206f] = 0x1.d77f2p+24 *
             exp(0x1.f5c28f5c28f5cp-2 * lgt + 0x1.83f5b33333333p+20 * rt_inv);
  k[r206b] = 0x1.48e7126367b0bp+16 *
             exp(0x1.37909f166b9a8p+0 * lgt - 0x1.cc62cf9f862aap+25 * rt_inv);
  k[r207f] = 0x1.406f4p+23 *
             exp(0x1.947ae147ae148p-1 * lgt - 0x1.6e26acccccccdp+21 * rt_inv);
  k[r207b] = 0x1.8a94074f4ac7fp+12 *
             exp(0x1.fbfe5eb31dec8p+0 * lgt - 0x1.ca74815d0b93fp+25 * rt_inv);
  k[r208f] = 0x1.28ep+14 *
             exp(0x1.ccccccccccccdp+0 * lgt - 0x1.4578700000001p+24 * rt_inv);
  k[r208b] = 0x1.8cda4926357efp-6 *
             exp(0x1.b3df9de82af6bp+1 * lgt - 0x1.fcbd3993adaf6p+23 * rt_inv);
  k[r209f] = 0x1.931ffffffffffp+14 *
             exp(0x1.a666666666666p+0 * lgt - 0x1.68f79p+23 * rt_inv);
  k[r209b] = 0x1.b61d6beb93f0dp+2 *
             exp(0x1.3c3f0f5fa072p+1 * lgt - 0x1.5145d4ee47da6p+23 * rt_inv);
  k[r210f] =
      0x1.d4cp+13 * exp(0x1.999999999999ap+0 * lgt - 0x1.83e8ap+23 * rt_inv);
  k[r210b] = 0x1.c1a3be8558a1ep-1 *
             exp(0x1.71a922612cc66p+1 * lgt - 0x1.b0ea949587cc2p+22 * rt_inv);
  k[r211f] =
      0x1.6f94p+16 * exp(0x1.b333333333333p+0 * lgt - 0x1.5c8488p+24 * rt_inv);
  k[r211b] = 0x1.2033ceeb922bap-4 *
             exp(0x1.a296ce7a72d9ep+1 * lgt - 0x1.7add2f09df678p+23 * rt_inv);
  k[r212f] =
      0x1.25cp+14 * exp(0x1.d99999999999ap+0 * lgt - 0x1.d1ccp+22 * rt_inv);
  k[r212b] = 0x1.767f3e9d99fc3p+1 *
             exp(0x1.515ca6584eb3ep+1 * lgt - 0x1.09e05323cc5bfp+20 * rt_inv);
  k[r213f] = 0x1.edbffffffffffp+13 * exp(0x1p+1 * lgt - 0x1.1bf9p+24 * rt_inv);
  k[r213b] = 0x1.15cb39484e612p-1 *
             exp(0x1.a06052f374b5ep+1 * lgt - 0x1.b90cdf81eb35fp+22 * rt_inv);
  k[r214f] = 0x1.c083126e978d4p-3 *
             exp(0x1.970a3d70a3d71p+1 * lgt - 0x1.3325c80000001p+25 * rt_inv);
  k[r214b] = 0x1.8e24cb8171c42p-14 *
             exp(0x1.1faf6f4728edp+2 * lgt - 0x1.4c2764b47a953p+25 * rt_inv);
  k[r215f] = 0x1.74bc6a7ef9db1p-1 *
             exp(0x1.7eb851eb851ecp+1 * lgt - 0x1.fb6c2p+24 * rt_inv);
  k[r215b] = 0x1.0d0130a0767d7p-4 *
             exp(0x1.c29f97b3dbd58p+1 * lgt - 0x1.3458278b2120ep+25 * rt_inv);
  k[r216f] = 0x1.28f5c28f5c28fp-3 *
             exp(0x1.7eb851eb851ecp+1 * lgt - 0x1.e85558p+24 * rt_inv);
  k[r216b] = 0x1.7a5261c20fd16p-9 *
             exp(0x1.fe70111bce91ap+1 * lgt - 0x1.05dc5ce24020ep+25 * rt_inv);
  k[r217f] = 0x1.0666666666666p+3 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.5727a8p+25 * rt_inv);
  k[r217b] = 0x1.004a9b085efc8p+1 *
             exp(0x1.4f52251dba044p+1 * lgt + 0x1.51bbd0057e5cap+24 * rt_inv);
  k[r218f] = 0x1.84cccccccccccp+4 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.f6c3080000002p+25 * rt_inv);
  k[r218b] = 0x1.d317f3545bfa9p-6 *
             exp(0x1.b3bf8073116cp+1 * lgt + 0x1.371816cb2db3bp+22 * rt_inv);
  k[r219f] = 0x1.c4feccp+31 * exp(-0x1.7f0e800000001p+26 * rt_inv);
  k[r219b] = 0x1.86841e03394dep+27 *
             exp(0x1.12f0e07d1ae2p-1 * lgt - 0x1.b24e12a8bfc64p+24 * rt_inv);
  k[r220f] = 0x1.6168e126c7b4dp+105 *
             exp(-0x1.78f5c28f5c28fp+2 * lgt - 0x1.9363fc0000001p+26 * rt_inv);
  k[r220b] = 0x1.5c02d4b2e66a9p+95 *
             exp(-0x1.918282cd32ecp+2 * lgt - 0x1.29019e442e75ap+25 * rt_inv);
  k[r221f] = 0x1.49b81dea72c89p+115 *
             exp(-0x1.bd70a3d70a3d7p+2 * lgt - 0x1.7bdd5p+26 * rt_inv);
  k[r221b] = 0x1.16de8291a824dp+85 *
             exp(-0x1.7bd3d2aa1884dp+2 * lgt - 0x1.74a82d724ddc6p+25 * rt_inv);
  k[r222f] = 0x1.312dp+25 * exp(-0x1.18e8800000001p+22 * rt_inv);
  k[r222b] = 0x1.a03e79ccf2e43p+25 *
             exp(-0x1.96039d144394p-2 * lgt - 0x1.0fb464c759389p+27 * rt_inv);
  k[r223f] = 0x1.df3b645a1cacp-2 *
             exp(0x1.947ae147ae148p+1 * lgt - 0x1.5779600000001p+24 * rt_inv);
  k[r223b] = 0x1.1077cdd9856e6p+8 *
             exp(0x1.2c81aaa8f0a8p+1 * lgt - 0x1.4a34daf715677p+27 * rt_inv);
  k[r224f] = 0x1.bf08ebp+34;
  k[r224b] = 0x1.3adfffe954b9ep+17 *
             exp(0x1.62ea3318b392p-1 * lgt - 0x1.01cd7fcbda8a3p+26 * rt_inv);
  k[r225f] = 0x1.bf08ebp+34;
  k[r225b] = 0x1.28841cb98791p+22 *
             exp(0x1.a3c36f3779f4p-1 * lgt - 0x1.71c99c8a4a69cp+28 * rt_inv);
  k[r226f] = 0x1.2a05f2p+33;
  k[r226b] = 0x1.57504d246ec89p+34 *
             exp(-0x1.06fe49cf7bap-3 * lgt - 0x1.99196581ffec9p+28 * rt_inv);
  k[r227f] = 0x1.671e344p+35 * exp(-0x1.404c98p+24 * rt_inv);
  k[r227b] = 0x1.b053230a5a04cp+33 *
             exp(0x1.2200b8ad4bdp-4 * lgt - 0x1.25ba9a1d20faep+27 * rt_inv);
  k[r228f] = 0x1.888d6e5cp+39 * exp(-0x1.3333333333333p+0 * lgt);
  k[r228b] = 0x1.d8932d0210783p+37 *
             exp(-0x1.211327a85eb2p+0 * lgt - 0x1.fb620e3a41f12p+26 * rt_inv);
  k[r229f] = 0x1.74876e8p+36;
  k[r229b] = 0x1.5a1ef3e12b1c8p+31 *
             exp(0x1.970385b6613p-2 * lgt - 0x1.532db958e3a0cp+28 * rt_inv);
  k[r230f] = 0x1.bf08ebp+34;
  k[r230b] = 0x1.05b63488564ffp+20 *
             exp(0x1.4982aa6ba26cp+0 * lgt - 0x1.6d2b8fb52e4a3p+28 * rt_inv);
  k[r231f] = 0x1.bf08ebp+34;
  k[r231b] = 0x1.15ea6e73e6b74p+15 *
             exp(0x1.29160c5c3f428p+0 * lgt - 0x1.deaa98eed4158p+25 * rt_inv);
  k[r232f] = 0x1.2a05f2p+35;
  k[r232b] = 0x1.970cf0fef94c9p+10 *
             exp(0x1.b0d05417ff68p-1 * lgt - 0x1.32cd1dffa8428p+26 * rt_inv);
  k[r233f] = 0x1.2a05f2p+32;
  k[r233b] = 0x1.2f043074d53b8p+31 *
             exp(0x1.5b04a657d828p-2 * lgt - 0x1.947b58ace3cccp+28 * rt_inv);
  k[r234f] = 0x1.74876e8p+36 * exp(-0x1.8f04700000001p+26 * rt_inv);
  k[r234b] = 0x1.43c79e48bb71p+24 *
             exp(0x1.56df1d88a6e8p-4 * lgt - 0x1.f6a43f369aeb6p+24 * rt_inv);
  k[r235f] = 0x1.dcd65p+29 * exp(0x1.18e8800000001p+22 * rt_inv);
  k[r235b] = 0x1.efec1e0a3446p+61 *
             exp(-0x1.8f876c36554dp+0 * lgt - 0x1.0a4edaf163864p+27 * rt_inv);
  k[r236f] = 0x1.4f46b04p+37 * exp(-0x1.8709780000001p+26 * rt_inv);
  k[r236b] = 0x1.cf569d8976727p-22 *
             exp(0x1.3457cdac952ap+1 * lgt - 0x1.f84172bf65994p+26 * rt_inv);
  k[r237f] = 0x1.cbabc8p+27 * exp(0x1.3087bp+23 * rt_inv);
  k[r237b] = 0x1.1a72c82f94c87p+46 *
             exp(-0x1.faa5e01c3828p-3 * lgt - 0x1.79b950a5dc179p+26 * rt_inv);
  k[r238f] = 0x1.2a05f2p+35;
  k[r238b] = 0x1.1282b4c18eb29p+36 *
             exp(0x1.2dcccc2f15ap-2 * lgt - 0x1.9f95922e6bd1bp+26 * rt_inv);
  k[r239f] = 0x1.380d2f7c2f57bp+140 *
             exp(-0x1.fd079e59f2baap+2 * lgt - 0x1.6e5301e666667p+28 * rt_inv);
  k[r239b] = 0x1.f0b924de84c57p+103 *
             exp(-0x1.9d7f8f72130a4p+2 * lgt - 0x1.f6ff838700211p+24 * rt_inv);
  k[r240f] = 0x1.a36p+12 * exp(0x1p+1 * lgt + 0x1.41b4cf5c28f5cp+21 * rt_inv);
  k[r240b] = 0x1.c846988b8844ap+7 *
             exp(0x1.0f31219598ebcp+1 * lgt - 0x1.4a50db3eca4f9p+26 * rt_inv);
  k[r241f] = 0x1.d01p+14 * exp(0x1p+1 * lgt - 0x1.0184503d70a3ep+24 * rt_inv);
  k[r241b] = 0x1.894e67ed1ca05p+6 *
             exp(0x1.1b2d002255812p+1 * lgt - 0x1.4b103f9481b7p+25 * rt_inv);
  k[r242f] = 0x1.b71758e219653p-6 *
             exp(0x1.e3923a29c779ap+1 * lgt - 0x1.3371c79999999p+25 * rt_inv);
  k[r242b] = 0x1.ee40bffb7afe4p-6 *
             exp(0x1.d99aa3e806748p+1 * lgt - 0x1.1d2732af7915p+26 * rt_inv);
  k[r243f] = 0x1.f1f2a11cc5b38p-20 *
             exp(0x1.528f5c28f5c29p+2 * lgt + 0x1.bd5e000000001p+18 * rt_inv);
  k[r243b] = 0x1.ef101b508969ep-29 *
             exp(0x1.5de7dae9aac7fp+2 * lgt - 0x1.35aa99a6aba16p+24 * rt_inv);
  k[r244f] = 0x1.2a05f2p+34 * exp(-0x1.0759f8p+26 * rt_inv);
  k[r244b] = 0x1.c60027c531581p+35 *
             exp(-0x1.13cc754b10124p-1 * lgt - 0x1.bdd1cc6ef090cp+24 * rt_inv);
  k[r245f] = 0x1.31794b4p+35 * exp(-0x1.6665ba0000001p+27 * rt_inv);
  k[r245b] = 0x1.8982919541065p+29 *
             exp(-0x1.31f6d6ff6b82p-3 * lgt + 0x1.ab0af5fd0c407p+23 * rt_inv);
  k[r246f] = 0x1.1f0e54p+29 * exp(-0x1.04187p+24 * rt_inv);
  k[r246b] = 0x1.66d31522c2cd5p+24 *
             exp(-0x1.8f141cfb9f4p-6 * lgt - 0x1.54ae8f1031b0fp+25 * rt_inv);
  k[r247f] = 0x1.2a05f2p+32 * exp(-0x1.1a583cp+26 * rt_inv);
  k[r247b] = 0x1.045875bb6d68cp+39 *
             exp(-0x1.103288c15e7eep+0 * lgt - 0x1.b929ddb908ad8p+24 * rt_inv);
  k[r248f] = 0x1.5d3ef798p+43 * exp(-0x1.9afce4p+26 * rt_inv);
  k[r248b] = 0x1.dc8ef5c315237p+6 *
             exp(0x1.f0971daa2f0ap+0 * lgt - 0x1.23156bb3fc5b1p+26 * rt_inv);
  k[r249f] = 0x1.671e344p+34;
  k[r249b] = 0x1.33d89e034a191p+34 *
             exp(0x1.c9db882acf1p-2 * lgt - 0x1.30e92b9af5868p+28 * rt_inv);
  k[r250f] = 0x1.5f5c28f5c28f6p+2 *
             exp(0x1.6666666666666p+1 * lgt - 0x1.763f1p+24 * rt_inv);
  k[r250b] = 0x1.460dbdc304a97p+1 *
             exp(0x1.86100c883948cp+1 * lgt - 0x1.910826ca627c7p+25 * rt_inv);
  k[r251f] = 0x1.0c388dp+33;
  k[r251b] = 0x1.ebc5d6d97ce22p+38 *
             exp(-0x1.2d5b520f3bfp-2 * lgt - 0x1.dad4c20e6ed8ep+26 * rt_inv);
  k[r252f] = 0x1.eff575daa5p+53 *
             exp(-0x1.51eb851eb851fp-1 * lgt - 0x1.761e6p+25 * rt_inv);
  k[r252b] = 0x1.fb6b56b37fc6bp+41 *
             exp(-0x1.751315bb83a28p-1 * lgt - 0x1.34914aacb4529p+25 * rt_inv);
  k[r253f] = 0x1.351609ff758p+60 *
             exp(-0x1.fae147ae147aep-1 * lgt - 0x1.3bc846p+28 * rt_inv);
  k[r253b] = 0x1.2d8550a146417p+25 *
             exp(0x1.4609cf789044ap-1 * lgt + 0x1.97b9c9010f7cap+23 * rt_inv);
  k[r254f] = 0x1.2a05f2p+33 * exp(-0x1.8c9f8c0000001p+27 * rt_inv);
  k[r254b] = 0x1.61f41fa184d42p+24 *
             exp(-0x1.1b81c56fb3614p+0 * lgt - 0x1.828980df3f295p+21 * rt_inv);
  k[r255f] = 0x1.6dap+14 *
             exp(0x1.9c28f5c28f5c3p+0 * lgt + 0x1.1e04000000001p+17 * rt_inv);
  k[r255b] = 0x1.6ec11602e154ap+6 *
             exp(0x1.8a909cbbf5f8p-1 * lgt - 0x1.48d1868047173p+26 * rt_inv);
  k[r256f] = 0x1.22dee4p+30 * exp(-0x1.0f54fp+26 * rt_inv);
  k[r256b] = 0x1.9886e253decf2p+28 *
             exp(-0x1.7f2925354e346p+0 * lgt - 0x1.b1c5c7895e8aep+24 * rt_inv);
  k[r257f] = 0x1.d6p+7 * exp(0x1.4p+1 * lgt - 0x1.1cbd2p+23 * rt_inv);
  k[r257b] = 0x1.aed024a5bcf64p-5 *
             exp(0x1.b81f10730dffep+0 * lgt - 0x1.9f06bcc119a12p+24 * rt_inv);
  k[r258f] = 0x1.1c6p+12 * exp(0x1p+1 * lgt - 0x1.3f36cp+24 * rt_inv);
  k[r258b] = 0x1.bc6aa1d26ef3fp+0 *
             exp(0x1.411715b4e4d74p+0 * lgt - 0x1.53ed85030063cp+25 * rt_inv);
  k[r259f] = 0x1.8bd66277c45cbp-11 *
             exp(0x1.bae147ae147aep+1 * lgt - 0x1.5dec18p+24 * rt_inv);
  k[r259b] = 0x1.9acbac9a55618p-14 *
             exp(0x1.36483c24703e8p+1 * lgt - 0x1.9fd6b7528f87ep+25 * rt_inv);
  k[r260f] = 0x1.054e88p+29 * exp(-0x1.3f36cp+24 * rt_inv);
  k[r260b] = 0x1.2d25ff82ba5e8p+21 *
             exp(-0x1.f6fe76076956p-1 * lgt - 0x1.5c41c49d68b01p+25 * rt_inv);
  k[r261f] = 0x1.b1b48d538p+42 *
             exp(-0x1.c28f5c28f5c29p+0 * lgt - 0x1.11b9d4p+26 * rt_inv);
  k[r261b] = 0x1.17ab97e49243p+6 *
             exp(0x1.d2c1d9465b008p+0 * lgt + 0x1.04771a051386cp+20 * rt_inv);
  k[r262f] = 0x1.60816ea4p+40 *
             exp(-0x1.c7ae147ae147bp+0 * lgt - 0x1.b927500000001p+25 * rt_inv);
  k[r262b] = 0x1.bc11614f45efp+9 *
             exp(0x1.716d3fa3a65ap+0 * lgt - 0x1.04e5adc749e6cp+27 * rt_inv);
  k[r263f] = 0x1.1c6712d7p+41 * exp(-0x1p-1 * lgt - 0x1.a6f558p+26 * rt_inv);
  k[r263b] = 0x1.3aea00b3e70ccp+26 *
             exp(0x1.c7a5a4e3p-10 * lgt - 0x1.dc5d2249148aap+26 * rt_inv);
  k[r264f] = 0x1.dcd65p+30;
  k[r264b] = 0x1.08b55543c55a2p+66 *
             exp(-0x1.a876aa27b1e2p+0 * lgt - 0x1.229b87cd14533p+27 * rt_inv);
  k[r265f] = 0x1.dcd65p+29 * exp(-0x1.7485c8p+25 * rt_inv);
  k[r265b] = 0x1.828d15375209cp+35 *
             exp(-0x1.a1be78fb72884p-1 * lgt - 0x1.e4a5bb4dcd7b8p+24 * rt_inv);
  k[r266f] = 0x1.15090bfbf0bc3p+67 * exp(-0x1.2p+2 * lgt);
  k[r266b] = 0x1.4cf129bad45b8p+31 *
             exp(-0x1.4604399cb2b83p+1 * lgt + 0x1.85f55b9b3d6bp+21 * rt_inv);
  k[r267f] = 0x1.dae5dcb1d5dep+65 * exp(-0x1.2p+2 * lgt);
  k[r267b] = 0x1.5c3e036983fcbp+30 *
             exp(-0x1.3df5eddd659fp+1 * lgt - 0x1.93197f7b12091p+28 * rt_inv);
  k[r268f] = 0x1.1d6a8d84d59fap+74 *
             exp(-0x1.0f5c28f5c28f6p+1 * lgt - 0x1.5dc7520000001p+27 * rt_inv);
  k[r268b] = 0x1.5b8419ccae71dp+32 *
             exp(0x1.42d98cce7d4bap-3 * lgt + 0x1.3cadee675d12bp+20 * rt_inv);
  k[r269f] = 0x1.1450dc240dp+53 *
             exp(-0x1.199999999999ap+0 * lgt - 0x1.496dcp+26 * rt_inv);
  k[r269b] = 0x1.b904697d528cap+14 *
             exp(0x1.f415ff66eeb68p-1 * lgt - 0x1.239e2396f7fcap+25 * rt_inv);
  k[r270f] = 0x1.7d784p+25 * exp(-0x1.febep+20 * rt_inv);
  k[r270b] = 0x1.0e4b0243598dfp+24 *
             exp(-0x1.4e1ae9c1ad1p-4 * lgt - 0x1.79f8ede17ab4p+27 * rt_inv);
  k[r271f] = 0x1.bf08ebp+35 * exp(-0x1.5727a8p+26 * rt_inv);
  k[r271b] = 0x1.3a7a52a5243d6p+41 *
             exp(-0x1.c03730b590058p-1 * lgt - 0x1.609743cc650f2p+25 * rt_inv);
  k[r272f] = 0x1.b48eb57ep+43 * exp(-0x1.4731b8p+26 * rt_inv);
  k[r272b] = 0x1.c38ca87dcf15ep-36 *
             exp(0x1.2746039da03c8p+2 * lgt - 0x1.46c24f02e2152p+27 * rt_inv);
  k[r273f] = 0x1.4dc938p+29;
  k[r273b] = 0x1.a5a6291e3c3cfp+64 *
             exp(-0x1.ace890a0ef64p+0 * lgt - 0x1.22cd7c96faae4p+27 * rt_inv);
  k[r274f] = 0x1.2a05f2p+35 * exp(-0x1.2745d80000001p+26 * rt_inv);
  k[r274b] = 0x1.b59bd0927a53ap+1 *
             exp(0x1.749b1afca846p+0 * lgt - 0x1.d9e08e36510f4p+27 * rt_inv);
  k[r275f] = 0x1.aa535d3d0cp+54 * exp(-0x1.3f36cp+27 * rt_inv);
  k[r275b] = 0x1.2f825d9a646f5p+16 *
             exp(0x1.04dbc068e2902p+1 * lgt + 0x1.295817e4901d8p+24 * rt_inv);
  k[r276f] = 0x1.74876e8p+36 * exp(-0x1.bee64p+25 * rt_inv);
  k[r276b] = 0x1.0959a3065fadfp+30 *
             exp(0x1.0df535d35735p-1 * lgt - 0x1.47e7488e4bcd1p+26 * rt_inv);
  k[r277f] = 0x1.355ecc730a8p+54 *
             exp(-0x1.5851eb851eb85p+1 * lgt - 0x1.12862p+26 * rt_inv);
  k[r277b] = 0x1.d69bf53d20066p+16 *
             exp(-0x1.214cf475307e4p-1 * lgt - 0x1.a845ecd2f0831p+23 * rt_inv);
  k[r278f] = 0x1.2de538c66fp+52 *
             exp(-0x1.4e147ae147ae1p+1 * lgt - 0x1.4c245c0000001p+26 * rt_inv);
  k[r278b] = 0x1.1c5271a615e65p+20 *
             exp(-0x1.b8b7c8b063f8p-1 * lgt - 0x1.1b04d9178989bp+27 * rt_inv);
  k[r279f] = 0x1.6bcc41e9p+46 * exp(-0x1.dba0f00000001p+25 * rt_inv);
  k[r279b] = 0x1.0ae1393bcf37cp+37 *
             exp(-0x1.a875778dcd94p-4 * lgt - 0x1.3b9280dac00fp+25 * rt_inv);
  k[r280f] = 0x1.05ef39b2p+42 * exp(-0x1.1c28f5c28f5c3p+0 * lgt);
  k[r280b] = 0x1.d5aafd97a2d7bp+71 *
             exp(-0x1.2c3823f4efb9cp+1 * lgt - 0x1.95ec726893ae7p+26 * rt_inv);
  k[r281f] = 0x1.56ba098p+34 * exp(-0x1.8f04700000001p+27 * rt_inv);
  k[r281b] = 0x1.b021da364f0aep-5 *
             exp(0x1.f006d60f75a2p+0 * lgt - 0x1.514850668eb38p+27 * rt_inv);
  k[r282f] = 0x1.b48eb57ep+43 * exp(-0x1.c6e138p+27 * rt_inv);
  k[r282b] = 0x1.1dc8fe2009a17p+24 *
             exp(0x1.810a4187857p-1 * lgt - 0x1.de85965c75c24p+27 * rt_inv);
  k[r283f] = 0x1.fdece9526af4p+61 *
             exp(-0x1.d70a3d70a3d71p-2 * lgt - 0x1.b022c2p+28 * rt_inv);
  k[r283b] = 0x1.4e9226c43c817p+24 *
             exp(0x1.fadce3eae9698p-1 * lgt + 0x1.4a8124e7f8a98p+23 * rt_inv);
  k[r284f] = 0x1.477ffffffffffp+11 *
             exp(0x1.07ae147ae147bp+1 * lgt - 0x1.d3d7000000001p+21 * rt_inv);
  k[r284b] = 0x1.133851a6f83fdp-5 *
             exp(0x1.5bf4c65004aep+1 * lgt - 0x1.2d1eeab7b247ap+26 * rt_inv);
  k[r285f] = 0x1.211p+14 *
             exp(0x1.828f5c28f5c29p+0 * lgt + 0x1.eb55800000002p+21 * rt_inv);
  k[r285b] = 0x1.6c780984c30cep-25 *
             exp(0x1.b94b191c34588p+1 * lgt + 0x1.15a5d665c5763p+25 * rt_inv);
  k[r286f] =
      0x1.09p+12 * exp(0x1.0cccccccccccdp+1 * lgt - 0x1.36c96p+24 * rt_inv);
  k[r286b] = 0x1.5af4abb31eb68p-8 *
             exp(0x1.6d0f5d2eacc7p+1 * lgt - 0x1.f3ec52e3ad868p+24 * rt_inv);
  k[r287f] = 0x1.c1451f6p+35 *
             exp(-0x1.6666666666666p-2 * lgt - 0x1.7d864p+23 * rt_inv);
  k[r287b] = 0x1.b94573b1f7dc9p-7 *
             exp(0x1.ae64f98f5289p+0 * lgt + 0x1.3afeed31a9a31p+26 * rt_inv);
  k[r288f] = 0x1.accf3dacbf296p-32 *
             exp(0x1.7333333333333p+2 * lgt - 0x1.18e8800000001p+23 * rt_inv);
  k[r288b] = 0x1.17aff57e25194p-65 *
             exp(0x1.e2a08ccb62e19p+2 * lgt + 0x1.294b6209e2107p+26 * rt_inv);
  k[r289f] = 0x1.dcd65p+29 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r289b] = 0x1.a4e9370c66111p-3 *
             exp(0x1.47be583eb8f48p+0 * lgt + 0x1.a97be5d0d234bp+26 * rt_inv);
  k[r290f] = 0x1.92738f5028p+50 *
             exp(-0x1.e666666666666p+0 * lgt - 0x1.7bdd5p+23 * rt_inv);
  k[r290b] = 0x1.cfb253d33cc36p+7 *
             exp(0x1.8a02780af1p-4 * lgt + 0x1.51331282e36aep+26 * rt_inv);

  w[r1f] = k[r1f] * c[sH] * c[sO2];
  w[r1b] = k[r1b] * c[sOH] * c[sO];
  w[r2f] = k[r2f] * c[sO] * c[sH2];
  w[r2b] = k[r2b] * c[sOH] * c[sH];
  w[r3f] = k[r3f] * c[sH2] * c[sOH];
  w[r3b] = k[r3b] * c[sH] * c[sH2O];
  w[r4f] = k[r4f] * c[sO] * c[sH2O];
  w[r4b] = k[r4b] * c[sOH] * c[sOH];
  w[r5f] = k[r5f] * c[sH2] * M[mM1];
  w[r5b] = k[r5b] * c[sH] * c[sH] * M[mM1];
  w[r6f] = k[r6f] * c[sH2] * c[sAR];
  w[r6b] = k[r6b] * c[sAR] * c[sH] * c[sH];
  w[r7f] = k[r7f] * c[sH2] * c[sHE];
  w[r7b] = k[r7b] * c[sHE] * c[sH] * c[sH];
  w[r8f] = k[r8f] * c[sO] * c[sO] * M[mM2];
  w[r8b] = k[r8b] * c[sO2] * M[mM2];
  w[r9f] = k[r9f] * c[sO] * c[sO] * c[sAR];
  w[r9b] = k[r9b] * c[sAR] * c[sO2];
  w[r10f] = k[r10f] * c[sO] * c[sO] * c[sHE];
  w[r10b] = k[r10b] * c[sHE] * c[sO2];
  w[r11f] = k[r11f] * c[sO] * c[sH] * M[mM3];
  w[r11b] = k[r11b] * c[sOH] * M[mM3];
  w[r12f] = k[r12f] * c[sH] * c[sOH] * M[mM4];
  w[r12b] = k[r12b] * c[sH2O] * M[mM4];
  w[r13f] = k[r13f] * c[sH] * c[sO2];
  w[r13b] = k[r13b] * c[sHO2];
  w[r14f] = k[r14f] * c[sHO2] * c[sH];
  w[r14b] = k[r14b] * c[sO2] * c[sH2];
  w[r15f] = k[r15f] * c[sHO2] * c[sH];
  w[r15b] = k[r15b] * c[sOH] * c[sOH];
  w[r16f] = k[r16f] * c[sHO2] * c[sO];
  w[r16b] = k[r16b] * c[sOH] * c[sO2];
  w[r17f] = k[r17f] * c[sHO2] * c[sOH];
  w[r17b] = k[r17b] * c[sO2] * c[sH2O];
  w[r18f] = k[r18f] * c[sHO2] * c[sHO2];
  w[r18b] = k[r18b] * c[sO2] * c[sH2O2];
  w[r19f] = k[r19f] * c[sHO2] * c[sHO2];
  w[r19b] = k[r19b] * c[sO2] * c[sH2O2];
  w[r20f] = k[r20f] * c[sH2O2];
  w[r20b] = k[r20b] * c[sOH] * c[sOH];
  w[r21f] = k[r21f] * c[sH2O2] * c[sH];
  w[r21b] = k[r21b] * c[sOH] * c[sH2O];
  w[r22f] = k[r22f] * c[sH2O2] * c[sH];
  w[r22b] = k[r22b] * c[sH2] * c[sHO2];
  w[r23f] = k[r23f] * c[sH2O2] * c[sO];
  w[r23b] = k[r23b] * c[sHO2] * c[sOH];
  w[r24f] = k[r24f] * c[sH2O2] * c[sOH];
  w[r24b] = k[r24b] * c[sH2O] * c[sHO2];
  w[r25f] = k[r25f] * c[sH2O2] * c[sOH];
  w[r25b] = k[r25b] * c[sH2O] * c[sHO2];
  w[r26f] = k[r26f] * c[sCO] * c[sO];
  w[r26b] = k[r26b] * c[sCO2];
  w[r27f] = k[r27f] * c[sCO] * c[sO2];
  w[r27b] = k[r27b] * c[sO] * c[sCO2];
  w[r28f] = k[r28f] * c[sCO] * c[sHO2];
  w[r28b] = k[r28b] * c[sOH] * c[sCO2];
  w[r29f] = k[r29f] * c[sCO] * c[sOH];
  w[r29b] = k[r29b] * c[sH] * c[sCO2];
  w[r30f] = k[r30f] * c[sHCO] * M[mM8];
  w[r30b] = k[r30b] * c[sCO] * c[sH] * M[mM8];
  w[r31f] = k[r31f] * c[sHCO] * c[sO2];
  w[r31b] = k[r31b] * c[sHO2] * c[sCO];
  w[r32f] = k[r32f] * c[sHCO] * c[sH];
  w[r32b] = k[r32b] * c[sH2] * c[sCO];
  w[r33f] = k[r33f] * c[sHCO] * c[sO];
  w[r33b] = k[r33b] * c[sOH] * c[sCO];
  w[r34f] = k[r34f] * c[sHCO] * c[sOH];
  w[r34b] = k[r34b] * c[sH2O] * c[sCO];
  w[r35f] = k[r35f] * c[sHCO] * c[sO];
  w[r35b] = k[r35b] * c[sH] * c[sCO2];
  w[r36f] = k[r36f] * c[sHCO] * c[sHO2];
  w[r36b] = k[r36b] * c[sH] * c[sOH] * c[sCO2];
  w[r37f] = k[r37f] * c[sHCO] * c[sHCO];
  w[r37b] = k[r37b] * c[sCO] * c[sCO] * c[sH2];
  w[r38f] = k[r38f] * c[sHCO] * c[sCH3];
  w[r38b] = k[r38b] * c[sCH4] * c[sCO];
  w[r39f] = k[r39f] * c[sHCO] * c[sHCO];
  w[r39b] = k[r39b] * c[sCO] * c[sCH2O];
  w[r40f] = k[r40f] * c[sCH2O] * M[mM9];
  w[r40b] = k[r40b] * c[sH] * c[sHCO] * M[mM9];
  w[r41f] = k[r41f] * c[sCH2O] * M[mM10];
  w[r41b] = k[r41b] * c[sH2] * c[sCO] * M[mM10];
  w[r42f] = k[r42f] * c[sCH2O] * c[sH];
  w[r42b] = k[r42b] * c[sH2] * c[sHCO];
  w[r43f] = k[r43f] * c[sCH2O] * c[sO];
  w[r43b] = k[r43b] * c[sOH] * c[sHCO];
  w[r44f] = k[r44f] * c[sCH2O] * c[sOH];
  w[r44b] = k[r44b] * c[sH2O] * c[sHCO];
  w[r45f] = k[r45f] * c[sCH2O] * c[sO2];
  w[r45b] = k[r45b] * c[sHO2] * c[sHCO];
  w[r46f] = k[r46f] * c[sCH2O] * c[sHO2];
  w[r46b] = k[r46b] * c[sH2O2] * c[sHCO];
  w[r47f] = k[r47f] * c[sCH2O] * c[sCH3];
  w[r47b] = k[r47b] * c[sCH4] * c[sHCO];
  w[r48f] = k[r48f] * c[sCH3] * c[sO];
  w[r48b] = k[r48b] * c[sH] * c[sCH2O];
  w[r49f] = k[r49f] * c[sCH3] * c[sO2];
  w[r49b] = k[r49b] * c[sO] * c[sCH3O];
  w[r50f] = k[r50f] * c[sCH3] * c[sO2];
  w[r50b] = k[r50b] * c[sOH] * c[sCH2O];
  w[r51f] = k[r51f] * c[sCH3] * c[sHO2];
  w[r51b] = k[r51b] * c[sOH] * c[sCH3O];
  w[r52f] = k[r52f] * c[sCH3] * c[sCH3];
  w[r52b] = k[r52b] * c[sC2H6];
  w[r53f] = k[r53f] * c[sCH3] * c[sH];
  w[r53b] = k[r53b] * c[sCH4];
  w[r54f] = k[r54f] * c[sCH4] * c[sH];
  w[r54b] = k[r54b] * c[sH2] * c[sCH3];
  w[r55f] = k[r55f] * c[sCH4] * c[sO];
  w[r55b] = k[r55b] * c[sOH] * c[sCH3];
  w[r56f] = k[r56f] * c[sCH4] * c[sOH];
  w[r56b] = k[r56b] * c[sH2O] * c[sCH3];
  w[r57f] = k[r57f] * c[sCH3] * c[sHO2];
  w[r57b] = k[r57b] * c[sO2] * c[sCH4];
  w[r58f] = k[r58f] * c[sCH4] * c[sHO2];
  w[r58b] = k[r58b] * c[sH2O2] * c[sCH3];
  w[r59f] = k[r59f] * c[sCH2OH] * M[mM13];
  w[r59b] = k[r59b] * c[sH] * c[sCH2O] * M[mM13];
  w[r60f] = k[r60f] * c[sCH2OH] * c[sH];
  w[r60b] = k[r60b] * c[sH2] * c[sCH2O];
  w[r61f] = k[r61f] * c[sCH2OH] * c[sH];
  w[r61b] = k[r61b] * c[sOH] * c[sCH3];
  w[r62f] = k[r62f] * c[sCH2OH] * c[sO];
  w[r62b] = k[r62b] * c[sOH] * c[sCH2O];
  w[r63f] = k[r63f] * c[sCH2OH] * c[sOH];
  w[r63b] = k[r63b] * c[sH2O] * c[sCH2O];
  w[r64f] = k[r64f] * c[sCH2OH] * c[sO2];
  w[r64b] = k[r64b] * c[sHO2] * c[sCH2O];
  w[r65f] = k[r65f] * c[sCH2OH] * c[sO2];
  w[r65b] = k[r65b] * c[sHO2] * c[sCH2O];
  w[r66f] = k[r66f] * c[sCH2OH] * c[sHO2];
  w[r66b] = k[r66b] * c[sH2O2] * c[sCH2O];
  w[r67f] = k[r67f] * c[sCH2OH] * c[sHCO];
  w[r67b] = k[r67b] * c[sCO] * c[sCH3OH];
  w[r68f] = k[r68f] * c[sCH2OH] * c[sHCO];
  w[r68b] = k[r68b] * c[sCH2O] * c[sCH2O];
  w[r69f] = k[r69f] * c[sCH2OH] * c[sCH2OH];
  w[r69b] = k[r69b] * c[sCH2O] * c[sCH3OH];
  w[r70f] = k[r70f] * c[sCH2OH] * c[sCH3O];
  w[r70b] = k[r70b] * c[sCH2O] * c[sCH3OH];
  w[r71f] = k[r71f] * c[sCH3O] * M[mM14];
  w[r71b] = k[r71b] * c[sH] * c[sCH2O] * M[mM14];
  w[r72f] = k[r72f] * c[sCH3O] * c[sH];
  w[r72b] = k[r72b] * c[sOH] * c[sCH3];
  w[r73f] = k[r73f] * c[sCH3O] * c[sO];
  w[r73b] = k[r73b] * c[sOH] * c[sCH2O];
  w[r74f] = k[r74f] * c[sCH3O] * c[sOH];
  w[r74b] = k[r74b] * c[sH2O] * c[sCH2O];
  w[r75f] = k[r75f] * c[sCH3O] * c[sO2];
  w[r75b] = k[r75b] * c[sHO2] * c[sCH2O];
  w[r76f] = k[r76f] * c[sCH3O] * c[sO2];
  w[r76b] = k[r76b] * c[sHO2] * c[sCH2O];
  w[r77f] = k[r77f] * c[sCH3O] * c[sHO2];
  w[r77b] = k[r77b] * c[sH2O2] * c[sCH2O];
  w[r78f] = k[r78f] * c[sCH3O] * c[sCO];
  w[r78b] = k[r78b] * c[sCO2] * c[sCH3];
  w[r79f] = k[r79f] * c[sCH3O] * c[sHCO];
  w[r79b] = k[r79b] * c[sCO] * c[sCH3OH];
  w[r80f] = k[r80f] * c[sCH3O] * c[sCH3O];
  w[r80b] = k[r80b] * c[sCH2O] * c[sCH3OH];
  w[r81f] = k[r81f] * c[sOH] * c[sCH3];
  w[r81b] = k[r81b] * c[sCH3OH];
  w[r82f] = k[r82f] * c[sH] * c[sCH2OH];
  w[r82b] = k[r82b] * c[sCH3OH];
  w[r83f] = k[r83f] * c[sH] * c[sCH3O];
  w[r83b] = k[r83b] * c[sCH3OH];
  w[r84f] = k[r84f] * c[sCH3OH] * c[sH];
  w[r84b] = k[r84b] * c[sH2] * c[sCH2OH];
  w[r85f] = k[r85f] * c[sCH3OH] * c[sH];
  w[r85b] = k[r85b] * c[sH2] * c[sCH3O];
  w[r86f] = k[r86f] * c[sCH3OH] * c[sO];
  w[r86b] = k[r86b] * c[sOH] * c[sCH2OH];
  w[r87f] = k[r87f] * c[sCH3OH] * c[sOH];
  w[r87b] = k[r87b] * c[sH2O] * c[sCH3O];
  w[r88f] = k[r88f] * c[sCH3OH] * c[sOH];
  w[r88b] = k[r88b] * c[sH2O] * c[sCH2OH];
  w[r89f] = k[r89f] * c[sCH3OH] * c[sO2];
  w[r89b] = k[r89b] * c[sHO2] * c[sCH2OH];
  w[r90f] = k[r90f] * c[sCH3OH] * c[sHCO];
  w[r90b] = k[r90b] * c[sCH2O] * c[sCH2OH];
  w[r91f] = k[r91f] * c[sCH3OH] * c[sHO2];
  w[r91b] = k[r91b] * c[sH2O2] * c[sCH2OH];
  w[r92f] = k[r92f] * c[sCH3OH] * c[sCH3];
  w[r92b] = k[r92b] * c[sCH4] * c[sCH2OH];
  w[r93f] = k[r93f] * c[sCH3O] * c[sCH3OH];
  w[r93b] = k[r93b] * c[sCH2OH] * c[sCH3OH];
  w[r94f] = k[r94f] * c[sCH3] * c[sCH3];
  w[r94b] = k[r94b] * c[sC2H5] * c[sH];
  w[r95f] = k[r95f] * c[sCH4] * c[sCH2];
  w[r95b] = k[r95b] * c[sCH3] * c[sCH3];
  w[r96f] = k[r96f] * c[sCH4] * c[sCH2S];
  w[r96b] = k[r96b] * c[sCH3] * c[sCH3];
  w[r97f] = k[r97f] * c[sCH3] * c[sOH];
  w[r97b] = k[r97b] * c[sH2O] * c[sCH2];
  w[r98f] = k[r98f] * c[sCH3] * c[sOH];
  w[r98b] = k[r98b] * c[sH2O] * c[sCH2S];
  w[r99f] = k[r99f] * c[sCH3] * c[sCH2];
  w[r99b] = k[r99b] * c[sH] * c[sC2H4];
  w[r100f] = k[r100f] * c[sCH3] * c[sCH2S];
  w[r100b] = k[r100b] * c[sH] * c[sC2H4];
  w[r101f] = k[r101f] * c[sCH3O] * c[sH];
  w[r101b] = k[r101b] * c[sH2O] * c[sCH2S];
  w[r102f] = k[r102f] * c[sCH2S] * c[sH2O];
  w[r102b] = k[r102b] * c[sCH3OH];
  w[r103f] = k[r103f] * c[sC2H6] * c[sH];
  w[r103b] = k[r103b] * c[sH2] * c[sC2H5];
  w[r104f] = k[r104f] * c[sC2H6] * c[sO];
  w[r104b] = k[r104b] * c[sOH] * c[sC2H5];
  w[r105f] = k[r105f] * c[sC2H6] * c[sOH];
  w[r105b] = k[r105b] * c[sH2O] * c[sC2H5];
  w[r106f] = k[r106f] * c[sC2H6] * c[sO2];
  w[r106b] = k[r106b] * c[sHO2] * c[sC2H5];
  w[r107f] = k[r107f] * c[sC2H6] * c[sHO2];
  w[r107b] = k[r107b] * c[sH2O2] * c[sC2H5];
  w[r108f] = k[r108f] * c[sC2H6] * c[sCH3];
  w[r108b] = k[r108b] * c[sCH4] * c[sC2H5];
  w[r109f] = k[r109f] * c[sC2H5] * c[sH];
  w[r109b] = k[r109b] * c[sC2H6];
  w[r110f] = k[r110f] * c[sC2H5] * c[sH];
  w[r110b] = k[r110b] * c[sH2] * c[sC2H4];
  w[r111f] = k[r111f] * c[sC2H5] * c[sO];
  w[r111b] = k[r111b] * c[sCH2O] * c[sCH3];
  w[r112f] = k[r112f] * c[sC2H5] * c[sO2];
  w[r112b] = k[r112b] * c[sHO2] * c[sC2H4];
  w[r113f] = k[r113f] * c[sC2H5] * c[sC2H5];
  w[r113b] = k[r113b] * c[sC2H6] * c[sC2H4];
  w[r114f] = k[r114f] * c[sC2H5] * c[sHCO];
  w[r114b] = k[r114b] * c[sCO] * c[sC2H6];
  w[r115f] = k[r115f] * c[sC2H5] * c[sO];
  w[r115b] = k[r115b] * c[sH] * c[sCH3HCO];
  w[r116f] = k[r116f] * c[sC2H4];
  w[r116b] = k[r116b] * c[sC2H2] * c[sH2];
  w[r117f] = k[r117f] * c[sC2H4] * c[sH];
  w[r117b] = k[r117b] * c[sC2H5];
  w[r118f] = k[r118f] * c[sC2H3] * c[sH];
  w[r118b] = k[r118b] * c[sC2H4];
  w[r119f] = k[r119f] * c[sC2H4] * c[sH];
  w[r119b] = k[r119b] * c[sH2] * c[sC2H3];
  w[r120f] = k[r120f] * c[sC2H4] * c[sOH];
  w[r120b] = k[r120b] * c[sH2O] * c[sC2H3];
  w[r121f] = k[r121f] * c[sC2H4] * c[sCH3];
  w[r121b] = k[r121b] * c[sCH4] * c[sC2H3];
  w[r122f] = k[r122f] * c[sC2H4] * c[sO];
  w[r122b] = k[r122b] * c[sHCO] * c[sCH3];
  w[r123f] = k[r123f] * c[sC2H3] * c[sOH];
  w[r123b] = k[r123b] * c[sH2O] * c[sC2H2];
  w[r124f] = k[r124f] * c[sC2H4] * c[sO];
  w[r124b] = k[r124b] * c[sC2H3] * c[sOH];
  w[r125f] = k[r125f] * c[sC2H4] * c[sO2];
  w[r125b] = k[r125b] * c[sHO2] * c[sC2H3];
  w[r126f] = k[r126f] * c[sC2H3] * c[sH];
  w[r126b] = k[r126b] * c[sH2] * c[sC2H2];
  w[r127f] = k[r127f] * c[sC2H3] * c[sO];
  w[r127b] = k[r127b] * c[sH] * c[sCH2CO];
  w[r128f] = k[r128f] * c[sC2H3] * c[sH2O2];
  w[r128b] = k[r128b] * c[sHO2] * c[sC2H4];
  w[r129f] = k[r129f] * c[sC2H3] * c[sCH3];
  w[r129b] = k[r129b] * c[sCH4] * c[sC2H2];
  w[r130f] = k[r130f] * c[sC2H3] * c[sC2H3];
  w[r130b] = k[r130b] * c[sC2H2] * c[sC2H4];
  w[r131f] = k[r131f] * c[sC2H3] * c[sO2];
  w[r131b] = k[r131b] * c[sCH2O] * c[sHCO];
  w[r132f] = k[r132f] * c[sC2H3] * c[sO2];
  w[r132b] = k[r132b] * c[sC2H2] * c[sHO2];
  w[r133f] = k[r133f] * c[sC2H3] * c[sO2];
  w[r133b] = k[r133b] * c[sCH2HCO] * c[sO];
  w[r134f] = k[r134f] * c[sC2H2] * c[sH];
  w[r134b] = k[r134b] * c[sC2H3];
  w[r135f] = k[r135f] * c[sC2H2] * c[sO];
  w[r135b] = k[r135b] * c[sH] * c[sHCCO];
  w[r136f] = k[r136f] * c[sC2H2] * c[sO];
  w[r136b] = k[r136b] * c[sOH] * c[sC2H];
  w[r137f] = k[r137f] * c[sC2H2] * c[sO];
  w[r137b] = k[r137b] * c[sCO] * c[sCH2];
  w[r138f] = k[r138f] * c[sC2H2] * c[sOH];
  w[r138b] = k[r138b] * c[sH] * c[sCH2CO];
  w[r139f] = k[r139f] * c[sC2H2] * c[sOH];
  w[r139b] = k[r139b] * c[sH] * c[sHCCOH];
  w[r140f] = k[r140f] * c[sC2H2] * c[sOH];
  w[r140b] = k[r140b] * c[sH2O] * c[sC2H];
  w[r141f] = k[r141f] * c[sC2H2] * c[sOH];
  w[r141b] = k[r141b] * c[sCO] * c[sCH3];
  w[r142f] = k[r142f] * c[sC2H2] * c[sHO2];
  w[r142b] = k[r142b] * c[sOH] * c[sCH2CO];
  w[r143f] = k[r143f] * c[sHCCO] * c[sH];
  w[r143b] = k[r143b] * c[sCO] * c[sCH2S];
  w[r144f] = k[r144f] * c[sHCCO] * c[sO];
  w[r144b] = k[r144b] * c[sCO] * c[sCO] * c[sH];
  w[r145f] = k[r145f] * c[sHCCO] * c[sO2];
  w[r145b] = k[r145b] * c[sCO] * c[sCO] * c[sOH];
  w[r146f] = k[r146f] * c[sHCCO] * c[sOH];
  w[r146b] = k[r146b] * c[sCO] * c[sH] * c[sHCO];
  w[r147f] = k[r147f] * c[sHCCO] * c[sCH2];
  w[r147b] = k[r147b] * c[sCO] * c[sC2H3];
  w[r148f] = k[r148f] * c[sHCCO] * c[sHCCO];
  w[r148b] = k[r148b] * c[sCO] * c[sCO] * c[sC2H2];
  w[r149f] = k[r149f] * c[sCH3] * c[sHCCO];
  w[r149b] = k[r149b] * c[sCO] * c[sC2H4];
  w[r150f] = k[r150f] * c[sC2H] * c[sH];
  w[r150b] = k[r150b] * c[sC2H2];
  w[r151f] = k[r151f] * c[sC2H] * c[sOH];
  w[r151b] = k[r151b] * c[sHCCO] * c[sH];
  w[r152f] = k[r152f] * c[sC2H] * c[sO2];
  w[r152b] = k[r152b] * c[sCO] * c[sHCO];
  w[r153f] = k[r153f] * c[sC2H] * c[sH2];
  w[r153b] = k[r153b] * c[sC2H2] * c[sH];
  w[r154f] = k[r154f] * c[sCH2] * c[sH];
  w[r154b] = k[r154b] * c[sCH3];
  w[r155f] = k[r155f] * c[sCH2] * c[sO];
  w[r155b] = k[r155b] * c[sH] * c[sHCO];
  w[r156f] = k[r156f] * c[sCH2] * c[sOH];
  w[r156b] = k[r156b] * c[sH] * c[sCH2O];
  w[r157f] = k[r157f] * c[sCH2] * c[sH2];
  w[r157b] = k[r157b] * c[sCH3] * c[sH];
  w[r158f] = k[r158f] * c[sCH2] * c[sO2];
  w[r158b] = k[r158b] * c[sOH] * c[sHCO];
  w[r159f] = k[r159f] * c[sCH2] * c[sHO2];
  w[r159b] = k[r159b] * c[sOH] * c[sCH2O];
  w[r160f] = k[r160f] * c[sCH2] * c[sCO];
  w[r160b] = k[r160b] * c[sCH2CO];
  w[r161f] = k[r161f] * c[sCH2] * c[sCH2];
  w[r161b] = k[r161b] * c[sH2] * c[sC2H2];
  w[r162f] = k[r162f] * c[sCH2S] * M[mM27];
  w[r162b] = k[r162b] * c[sCH2] * M[mM27];
  w[r163f] = k[r163f] * c[sCH2S] * c[sH2O];
  w[r163b] = k[r163b] * c[sH2O] * c[sCH2];
  w[r164f] = k[r164f] * c[sCH2S] * c[sCO];
  w[r164b] = k[r164b] * c[sCO] * c[sCH2];
  w[r165f] = k[r165f] * c[sCH2S] * c[sCO2];
  w[r165b] = k[r165b] * c[sCO2] * c[sCH2];
  w[r166f] = k[r166f] * c[sCH2S] * c[sAR];
  w[r166b] = k[r166b] * c[sAR] * c[sCH2];
  w[r167f] = k[r167f] * c[sCH2S] * c[sO];
  w[r167b] = k[r167b] * c[sH2] * c[sCO];
  w[r168f] = k[r168f] * c[sCH2S] * c[sO];
  w[r168b] = k[r168b] * c[sH] * c[sHCO];
  w[r169f] = k[r169f] * c[sCH2S] * c[sOH];
  w[r169b] = k[r169b] * c[sH] * c[sCH2O];
  w[r170f] = k[r170f] * c[sCH2S] * c[sH2];
  w[r170b] = k[r170b] * c[sH] * c[sCH3];
  w[r171f] = k[r171f] * c[sCH2S] * c[sO2];
  w[r171b] = k[r171b] * c[sCO] * c[sOH] * c[sH];
  w[r172f] = k[r172f] * c[sCH2S] * c[sO2];
  w[r172b] = k[r172b] * c[sH2O] * c[sCO];
  w[r173f] = k[r173f] * c[sCH2S] * c[sCO2];
  w[r173b] = k[r173b] * c[sCO] * c[sCH2O];
  w[r174f] = k[r174f] * c[sCH3HCO];
  w[r174b] = k[r174b] * c[sHCO] * c[sCH3];
  w[r175f] = k[r175f] * c[sCH3CO];
  w[r175b] = k[r175b] * c[sCO] * c[sCH3];
  w[r176f] = k[r176f] * c[sCH3HCO] * c[sOH];
  w[r176b] = k[r176b] * c[sH2O] * c[sCH3CO];
  w[r177f] = k[r177f] * c[sCH3HCO] * c[sOH];
  w[r177b] = k[r177b] * c[sH2O] * c[sCH2HCO];
  w[r178f] = k[r178f] * c[sCH3HCO] * c[sO];
  w[r178b] = k[r178b] * c[sOH] * c[sCH3CO];
  w[r179f] = k[r179f] * c[sCH3HCO] * c[sO];
  w[r179b] = k[r179b] * c[sOH] * c[sCH2HCO];
  w[r180f] = k[r180f] * c[sCH3HCO] * c[sH];
  w[r180b] = k[r180b] * c[sH2] * c[sCH3CO];
  w[r181f] = k[r181f] * c[sCH3HCO] * c[sH];
  w[r181b] = k[r181b] * c[sH2] * c[sCH2HCO];
  w[r182f] = k[r182f] * c[sCH3HCO] * c[sCH3];
  w[r182b] = k[r182b] * c[sCH4] * c[sCH3CO];
  w[r183f] = k[r183f] * c[sCH3HCO] * c[sCH3];
  w[r183b] = k[r183b] * c[sCH4] * c[sCH2HCO];
  w[r184f] = k[r184f] * c[sCH3HCO] * c[sHO2];
  w[r184b] = k[r184b] * c[sH2O2] * c[sCH3CO];
  w[r185f] = k[r185f] * c[sCH3HCO] * c[sHO2];
  w[r185b] = k[r185b] * c[sH2O2] * c[sCH2HCO];
  w[r186f] = k[r186f] * c[sCH3HCO] * c[sO2];
  w[r186b] = k[r186b] * c[sHO2] * c[sCH3CO];
  w[r187f] = k[r187f] * c[sCH2CO] * c[sO];
  w[r187b] = k[r187b] * c[sCH2] * c[sCO2];
  w[r188f] = k[r188f] * c[sCH2CO] * c[sH];
  w[r188b] = k[r188b] * c[sCO] * c[sCH3];
  w[r189f] = k[r189f] * c[sCH2CO] * c[sH];
  w[r189b] = k[r189b] * c[sH2] * c[sHCCO];
  w[r190f] = k[r190f] * c[sCH2CO] * c[sO];
  w[r190b] = k[r190b] * c[sOH] * c[sHCCO];
  w[r191f] = k[r191f] * c[sCH2CO] * c[sOH];
  w[r191b] = k[r191b] * c[sH2O] * c[sHCCO];
  w[r192f] = k[r192f] * c[sCH2CO] * c[sOH];
  w[r192b] = k[r192b] * c[sCO] * c[sCH2OH];
  w[r193f] = k[r193f] * c[sCH2HCO] * c[sH];
  w[r193b] = k[r193b] * c[sHCO] * c[sCH3];
  w[r194f] = k[r194f] * c[sCH2HCO] * c[sH];
  w[r194b] = k[r194b] * c[sH2] * c[sCH2CO];
  w[r195f] = k[r195f] * c[sCH2HCO] * c[sO];
  w[r195b] = k[r195b] * c[sHCO] * c[sCH2O];
  w[r196f] = k[r196f] * c[sCH2HCO] * c[sOH];
  w[r196b] = k[r196b] * c[sH2O] * c[sCH2CO];
  w[r197f] = k[r197f] * c[sCH2HCO] * c[sO2];
  w[r197b] = k[r197b] * c[sOH] * c[sCO] * c[sCH2O];
  w[r198f] = k[r198f] * c[sCH2HCO] * c[sCH3];
  w[r198b] = k[r198b] * c[sH] * c[sCO] * c[sC2H5];
  w[r199f] = k[r199f] * c[sCH2HCO] * c[sHO2];
  w[r199b] = k[r199b] * c[sOH] * c[sHCO] * c[sCH2O];
  w[r200f] = k[r200f] * c[sCH2HCO] * c[sHO2];
  w[r200b] = k[r200b] * c[sO2] * c[sCH3HCO];
  w[r201f] = k[r201f] * c[sCH2HCO];
  w[r201b] = k[r201b] * c[sCO] * c[sCH3];
  w[r202f] = k[r202f] * c[sCH2HCO];
  w[r202b] = k[r202b] * c[sH] * c[sCH2CO];
  w[r203f] = k[r203f] * c[sC2H5OH];
  w[r203b] = k[r203b] * c[sCH2OH] * c[sCH3];
  w[r204f] = k[r204f] * c[sC2H5OH];
  w[r204b] = k[r204b] * c[sH2O] * c[sC2H4];
  w[r205f] = k[r205f] * c[sC2H5OH] * c[sOH];
  w[r205b] = k[r205b] * c[sH2O] * c[sC2H4OH];
  w[r206f] = k[r206f] * c[sC2H5OH] * c[sOH];
  w[r206b] = k[r206b] * c[sH2O] * c[sCH3CHOH];
  w[r207f] = k[r207f] * c[sC2H5OH] * c[sOH];
  w[r207b] = k[r207b] * c[sH2O] * c[sCH3CH2O];
  w[r208f] = k[r208f] * c[sC2H5OH] * c[sH];
  w[r208b] = k[r208b] * c[sH2] * c[sC2H4OH];
  w[r209f] = k[r209f] * c[sC2H5OH] * c[sH];
  w[r209b] = k[r209b] * c[sH2] * c[sCH3CHOH];
  w[r210f] = k[r210f] * c[sC2H5OH] * c[sH];
  w[r210b] = k[r210b] * c[sH2] * c[sCH3CH2O];
  w[r211f] = k[r211f] * c[sC2H5OH] * c[sO];
  w[r211b] = k[r211b] * c[sOH] * c[sC2H4OH];
  w[r212f] = k[r212f] * c[sC2H5OH] * c[sO];
  w[r212b] = k[r212b] * c[sOH] * c[sCH3CHOH];
  w[r213f] = k[r213f] * c[sC2H5OH] * c[sO];
  w[r213b] = k[r213b] * c[sOH] * c[sCH3CH2O];
  w[r214f] = k[r214f] * c[sC2H5OH] * c[sCH3];
  w[r214b] = k[r214b] * c[sCH4] * c[sC2H4OH];
  w[r215f] = k[r215f] * c[sC2H5OH] * c[sCH3];
  w[r215b] = k[r215b] * c[sCH4] * c[sCH3CHOH];
  w[r216f] = k[r216f] * c[sC2H5OH] * c[sCH3];
  w[r216b] = k[r216b] * c[sCH4] * c[sCH3CH2O];
  w[r217f] = k[r217f] * c[sC2H5OH] * c[sHO2];
  w[r217b] = k[r217b] * c[sH2O2] * c[sCH3CHOH];
  w[r218f] = k[r218f] * c[sC2H5OH] * c[sHO2];
  w[r218b] = k[r218b] * c[sH2O2] * c[sC2H4OH];
  w[r219f] = k[r219f] * c[sC2H5OH] * c[sHO2];
  w[r219b] = k[r219b] * c[sH2O2] * c[sCH3CH2O];
  w[r220f] = k[r220f] * c[sCH3CH2O] * M[mM29];
  w[r220b] = k[r220b] * c[sH] * c[sCH3HCO] * M[mM29];
  w[r221f] = k[r221f] * c[sCH3CH2O] * M[mM30];
  w[r221b] = k[r221b] * c[sCH2O] * c[sCH3] * M[mM30];
  w[r222f] = k[r222f] * c[sCH3CH2O] * c[sO2];
  w[r222b] = k[r222b] * c[sHO2] * c[sCH3HCO];
  w[r223f] = k[r223f] * c[sCH3CH2O] * c[sCO];
  w[r223b] = k[r223b] * c[sCO2] * c[sC2H5];
  w[r224f] = k[r224f] * c[sCH3CH2O] * c[sH];
  w[r224b] = k[r224b] * c[sCH2OH] * c[sCH3];
  w[r225f] = k[r225f] * c[sCH3CH2O] * c[sH];
  w[r225b] = k[r225b] * c[sH2O] * c[sC2H4];
  w[r226f] = k[r226f] * c[sCH3CH2O] * c[sOH];
  w[r226b] = k[r226b] * c[sH2O] * c[sCH3HCO];
  w[r227f] = k[r227f] * c[sCH3CHOH] * c[sO2];
  w[r227b] = k[r227b] * c[sHO2] * c[sCH3HCO];
  w[r228f] = k[r228f] * c[sCH3CHOH] * c[sO2];
  w[r228b] = k[r228b] * c[sHO2] * c[sCH3HCO];
  w[r229f] = k[r229f] * c[sCH3CHOH] * c[sO];
  w[r229b] = k[r229b] * c[sOH] * c[sCH3HCO];
  w[r230f] = k[r230f] * c[sCH3CHOH] * c[sH];
  w[r230b] = k[r230b] * c[sH2O] * c[sC2H4];
  w[r231f] = k[r231f] * c[sCH3CHOH] * c[sH];
  w[r231b] = k[r231b] * c[sCH2OH] * c[sCH3];
  w[r232f] = k[r232f] * c[sCH3CHOH] * c[sHO2];
  w[r232b] = k[r232b] * c[sOH] * c[sOH] * c[sCH3HCO];
  w[r233f] = k[r233f] * c[sCH3CHOH] * c[sOH];
  w[r233b] = k[r233b] * c[sH2O] * c[sCH3HCO];
  w[r234f] = k[r234f] * c[sCH3CHOH] * M[mM31];
  w[r234b] = k[r234b] * c[sH] * c[sCH3HCO] * M[mM31];
  w[r235f] = k[r235f] * c[sC2H4OH] * c[sO2];
  w[r235b] = k[r235b] * c[sHOC2H4O2];
  w[r236f] = k[r236f] * c[sHOC2H4O2];
  w[r236b] = k[r236b] * c[sOH] * c[sCH2O] * c[sCH2O];
  w[r237f] = k[r237f] * c[sC2H4] * c[sOH];
  w[r237b] = k[r237b] * c[sC2H4OH];
  w[r238f] = k[r238f] * c[sC2H5] * c[sHO2];
  w[r238b] = k[r238b] * c[sOH] * c[sCH3CH2O];
  w[r239f] = k[r239f] * c[sCH3OCH3];
  w[r239b] = k[r239b] * c[sCH3O] * c[sCH3];
  w[r240f] = k[r240f] * c[sCH3OCH3] * c[sOH];
  w[r240b] = k[r240b] * c[sH2O] * c[sCH3OCH2];
  w[r241f] = k[r241f] * c[sCH3OCH3] * c[sH];
  w[r241b] = k[r241b] * c[sH2] * c[sCH3OCH2];
  w[r242f] = k[r242f] * c[sCH3OCH3] * c[sCH3];
  w[r242b] = k[r242b] * c[sCH4] * c[sCH3OCH2];
  w[r243f] = k[r243f] * c[sCH3OCH3] * c[sO];
  w[r243b] = k[r243b] * c[sOH] * c[sCH3OCH2];
  w[r244f] = k[r244f] * c[sCH3OCH3] * c[sHO2];
  w[r244b] = k[r244b] * c[sH2O2] * c[sCH3OCH2];
  w[r245f] = k[r245f] * c[sCH3OCH3] * c[sO2];
  w[r245b] = k[r245b] * c[sHO2] * c[sCH3OCH2];
  w[r246f] = k[r246f] * c[sCH3OCH3] * c[sCH3O];
  w[r246b] = k[r246b] * c[sCH3OH] * c[sCH3OCH2];
  w[r247f] = k[r247f] * c[sCH3OCH3] * c[sCH3OCH2O2];
  w[r247b] = k[r247b] * c[sCH3OCH2O2H] * c[sCH3OCH2];
  w[r248f] = k[r248f] * c[sCH3OCH2];
  w[r248b] = k[r248b] * c[sCH3] * c[sCH2O];
  w[r249f] = k[r249f] * c[sCH3OCH2] * c[sCH3O];
  w[r249b] = k[r249b] * c[sCH2O] * c[sCH3OCH3];
  w[r250f] = k[r250f] * c[sCH3OCH2] * c[sCH2O];
  w[r250b] = k[r250b] * c[sHCO] * c[sCH3OCH3];
  w[r251f] = k[r251f] * c[sCH3OCH2] * c[sHO2];
  w[r251b] = k[r251b] * c[sOH] * c[sCH3OCH2O];
  w[r252f] = k[r252f] * c[sCH3OCH2O];
  w[r252b] = k[r252b] * c[sH] * c[sCH3OCHO];
  w[r253f] = k[r253f] * c[sCH3OCHO];
  w[r253b] = k[r253b] * c[sOCHO] * c[sCH3];
  w[r254f] = k[r254f] * c[sCH3OCHO] * c[sO2];
  w[r254b] = k[r254b] * c[sHO2] * c[sCH3OCO];
  w[r255f] = k[r255f] * c[sCH3OCHO] * c[sOH];
  w[r255b] = k[r255b] * c[sH2O] * c[sCH3OCO];
  w[r256f] = k[r256f] * c[sCH3OCHO] * c[sHO2];
  w[r256b] = k[r256b] * c[sH2O2] * c[sCH3OCO];
  w[r257f] = k[r257f] * c[sCH3OCHO] * c[sO];
  w[r257b] = k[r257b] * c[sOH] * c[sCH3OCO];
  w[r258f] = k[r258f] * c[sCH3OCHO] * c[sH];
  w[r258b] = k[r258b] * c[sH2] * c[sCH3OCO];
  w[r259f] = k[r259f] * c[sCH3OCHO] * c[sCH3];
  w[r259b] = k[r259b] * c[sCH4] * c[sCH3OCO];
  w[r260f] = k[r260f] * c[sCH3OCHO] * c[sCH3O];
  w[r260b] = k[r260b] * c[sCH3OH] * c[sCH3OCO];
  w[r261f] = k[r261f] * c[sCH3OCO];
  w[r261b] = k[r261b] * c[sCO] * c[sCH3O];
  w[r262f] = k[r262f] * c[sCH3OCO];
  w[r262b] = k[r262b] * c[sCO2] * c[sCH3];
  w[r263f] = k[r263f] * c[sOCHO] * M[mM32];
  w[r263b] = k[r263b] * c[sCO2] * c[sH] * M[mM32];
  w[r264f] = k[r264f] * c[sCH3OCH2] * c[sO2];
  w[r264b] = k[r264b] * c[sCH3OCH2O2];
  w[r265f] = k[r265f] * c[sCH3OCH2O2] * c[sCH2O];
  w[r265b] = k[r265b] * c[sHCO] * c[sCH3OCH2O2H];
  w[r266f] = k[r266f] * c[sCH3OCH2O2] * c[sCH3OCH2O2];
  w[r266b] = k[r266b] * c[sCH3OCH2O] * c[sCH3OCH2O] * c[sO2];
  w[r267f] = k[r267f] * c[sCH3OCH2O2] * c[sCH3OCH2O2];
  w[r267b] = k[r267b] * c[sCH3OCH2OH] * c[sCH3OCHO] * c[sO2];
  w[r268f] = k[r268f] * c[sCH3OCH2O2H];
  w[r268b] = k[r268b] * c[sOH] * c[sCH3OCH2O];
  w[r269f] = k[r269f] * c[sCH3OCH2O];
  w[r269b] = k[r269b] * c[sCH2O] * c[sCH3O];
  w[r270f] = k[r270f] * c[sCH3OCH2O] * c[sO2];
  w[r270b] = k[r270b] * c[sHO2] * c[sCH3OCHO];
  w[r271f] = k[r271f] * c[sCH3OCH2O2];
  w[r271b] = k[r271b] * c[sCH2OCH2O2H];
  w[r272f] = k[r272f] * c[sCH2OCH2O2H];
  w[r272b] = k[r272b] * c[sCH2O] * c[sCH2O] * c[sOH];
  w[r273f] = k[r273f] * c[sCH2OCH2O2H] * c[sO2];
  w[r273b] = k[r273b] * c[sO2CH2OCH2O2H];
  w[r274f] = k[r274f] * c[sO2CH2OCH2O2H];
  w[r274b] = k[r274b] * c[sOH] * c[sHO2CH2OCHO];
  w[r275f] = k[r275f] * c[sHO2CH2OCHO];
  w[r275b] = k[r275b] * c[sOH] * c[sOCH2OCHO];
  w[r276f] = k[r276f] * c[sOCH2OCHO];
  w[r276b] = k[r276b] * c[sHOCH2OCO];
  w[r277f] = k[r277f] * c[sHOCH2OCO];
  w[r277b] = k[r277b] * c[sCO] * c[sHOCH2O];
  w[r278f] = k[r278f] * c[sHOCH2OCO];
  w[r278b] = k[r278b] * c[sCO2] * c[sCH2OH];
  w[r279f] = k[r279f] * c[sHOCH2O];
  w[r279b] = k[r279b] * c[sH] * c[sHCOOH];
  w[r280f] = k[r280f] * c[sCH2O] * c[sOH];
  w[r280b] = k[r280b] * c[sHOCH2O];
  w[r281f] = k[r281f] * c[sHCOOH] * M[mM33];
  w[r281b] = k[r281b] * c[sH2O] * c[sCO] * M[mM33];
  w[r282f] = k[r282f] * c[sHCOOH] * M[mM34];
  w[r282b] = k[r282b] * c[sH2] * c[sCO2] * M[mM34];
  w[r283f] = k[r283f] * c[sHCOOH];
  w[r283b] = k[r283b] * c[sOH] * c[sHCO];
  w[r284f] = k[r284f] * c[sHCOOH] * c[sOH];
  w[r284b] = k[r284b] * c[sH] * c[sCO2] * c[sH2O];
  w[r285f] = k[r285f] * c[sHCOOH] * c[sOH];
  w[r285b] = k[r285b] * c[sOH] * c[sCO] * c[sH2O];
  w[r286f] = k[r286f] * c[sHCOOH] * c[sH];
  w[r286b] = k[r286b] * c[sH] * c[sCO2] * c[sH2];
  w[r287f] = k[r287f] * c[sHCOOH] * c[sH];
  w[r287b] = k[r287b] * c[sOH] * c[sCO] * c[sH2];
  w[r288f] = k[r288f] * c[sHCOOH] * c[sCH3];
  w[r288b] = k[r288b] * c[sOH] * c[sCO] * c[sCH4];
  w[r289f] = k[r289f] * c[sHCOOH] * c[sHO2];
  w[r289b] = k[r289b] * c[sOH] * c[sCO] * c[sH2O2];
  w[r290f] = k[r290f] * c[sHCOOH] * c[sO];
  w[r290b] = k[r290b] * c[sOH] * c[sOH] * c[sCO];

  cdot[sH] =
      -w[r1f] + w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] + 0x1p+1 * w[r5f] -
      0x1p+1 * w[r5b] + 0x1p+1 * w[r6f] - 0x1p+1 * w[r6b] + 0x1p+1 * w[r7f] -
      0x1p+1 * w[r7b] - w[r11f] + w[r11b] - w[r12f] + w[r12b] - w[r13f] +
      w[r13b] - w[r14f] + w[r14b] - w[r15f] + w[r15b] - w[r21f] + w[r21b] -
      w[r22f] + w[r22b] + w[r29f] - w[r29b] + w[r30f] - w[r30b] - w[r32f] +
      w[r32b] + w[r35f] - w[r35b] + w[r36f] - w[r36b] + w[r40f] - w[r40b] -
      w[r42f] + w[r42b] + w[r48f] - w[r48b] - w[r53f] + w[r53b] - w[r54f] +
      w[r54b] + w[r59f] - w[r59b] - w[r60f] + w[r60b] - w[r61f] + w[r61b] +
      w[r71f] - w[r71b] - w[r72f] + w[r72b] - w[r82f] + w[r82b] - w[r83f] +
      w[r83b] - w[r84f] + w[r84b] - w[r85f] + w[r85b] + w[r94f] - w[r94b] +
      w[r99f] - w[r99b] + w[r100f] - w[r100b] - w[r101f] + w[r101b] - w[r103f] +
      w[r103b] - w[r109f] + w[r109b] - w[r110f] + w[r110b] + w[r115f] -
      w[r115b] - w[r117f] + w[r117b] - w[r118f] + w[r118b] - w[r119f] +
      w[r119b] - w[r126f] + w[r126b] + w[r127f] - w[r127b] - w[r134f] +
      w[r134b] + w[r135f] - w[r135b] + w[r138f] - w[r138b] + w[r139f] -
      w[r139b] - w[r143f] + w[r143b] + w[r144f] - w[r144b] + w[r146f] -
      w[r146b] - w[r150f] + w[r150b] + w[r151f] - w[r151b] + w[r153f] -
      w[r153b] - w[r154f] + w[r154b] + w[r155f] - w[r155b] + w[r156f] -
      w[r156b] + w[r157f] - w[r157b] + w[r168f] - w[r168b] + w[r169f] -
      w[r169b] + w[r170f] - w[r170b] + w[r171f] - w[r171b] - w[r180f] +
      w[r180b] - w[r181f] + w[r181b] - w[r188f] + w[r188b] - w[r189f] +
      w[r189b] - w[r193f] + w[r193b] - w[r194f] + w[r194b] + w[r198f] -
      w[r198b] + w[r202f] - w[r202b] - w[r208f] + w[r208b] - w[r209f] +
      w[r209b] - w[r210f] + w[r210b] + w[r220f] - w[r220b] - w[r224f] +
      w[r224b] - w[r225f] + w[r225b] - w[r230f] + w[r230b] - w[r231f] +
      w[r231b] + w[r234f] - w[r234b] - w[r241f] + w[r241b] + w[r252f] -
      w[r252b] - w[r258f] + w[r258b] + w[r263f] - w[r263b] + w[r279f] -
      w[r279b] + w[r284f] - w[r284b] - w[r286f] + w[r286f] - w[r286b] +
      w[r286b] - w[r287f] + w[r287b];

  cdot[sH2] = -w[r2f] + w[r2b] - w[r3f] + w[r3b] - w[r5f] + w[r5b] - w[r6f] +
              w[r6b] - w[r7f] + w[r7b] + w[r14f] - w[r14b] + w[r22f] - w[r22b] +
              w[r32f] - w[r32b] + w[r37f] - w[r37b] + w[r41f] - w[r41b] +
              w[r42f] - w[r42b] + w[r54f] - w[r54b] + w[r60f] - w[r60b] +
              w[r84f] - w[r84b] + w[r85f] - w[r85b] + w[r103f] - w[r103b] +
              w[r110f] - w[r110b] + w[r116f] - w[r116b] + w[r119f] - w[r119b] +
              w[r126f] - w[r126b] - w[r153f] + w[r153b] - w[r157f] + w[r157b] +
              w[r161f] - w[r161b] + w[r167f] - w[r167b] - w[r170f] + w[r170b] +
              w[r180f] - w[r180b] + w[r181f] - w[r181b] + w[r189f] - w[r189b] +
              w[r194f] - w[r194b] + w[r208f] - w[r208b] + w[r209f] - w[r209b] +
              w[r210f] - w[r210b] + w[r241f] - w[r241b] + w[r258f] - w[r258b] +
              w[r282f] - w[r282b] + w[r286f] - w[r286b] + w[r287f] - w[r287b];

  cdot[sCH2] = -w[r95f] + w[r95b] + w[r97f] - w[r97b] - w[r99f] + w[r99b] +
               w[r137f] - w[r137b] - w[r147f] + w[r147b] - w[r154f] + w[r154b] -
               w[r155f] + w[r155b] - w[r156f] + w[r156b] - w[r157f] + w[r157b] -
               w[r158f] + w[r158b] - w[r159f] + w[r159b] - w[r160f] + w[r160b] -
               0x1p+1 * w[r161f] + 0x1p+1 * w[r161b] + w[r162f] - w[r162b] +
               w[r163f] - w[r163b] + w[r164f] - w[r164b] + w[r165f] - w[r165b] +
               w[r166f] - w[r166b] + w[r187f] - w[r187b];

  cdot[sCH2S] = -w[r96f] + w[r96b] + w[r98f] - w[r98b] - w[r100f] + w[r100b] +
                w[r101f] - w[r101b] - w[r102f] + w[r102b] + w[r143f] -
                w[r143b] - w[r162f] + w[r162b] - w[r163f] + w[r163b] -
                w[r164f] + w[r164b] - w[r165f] + w[r165b] - w[r166f] +
                w[r166b] - w[r167f] + w[r167b] - w[r168f] + w[r168b] -
                w[r169f] + w[r169b] - w[r170f] + w[r170b] - w[r171f] +
                w[r171b] - w[r172f] + w[r172b] - w[r173f] + w[r173b];

  cdot[sCH3] = -w[r38f] + w[r38b] - w[r47f] + w[r47b] - w[r48f] + w[r48b] -
               w[r49f] + w[r49b] - w[r50f] + w[r50b] - w[r51f] + w[r51b] -
               0x1p+1 * w[r52f] + 0x1p+1 * w[r52b] - w[r53f] + w[r53b] +
               w[r54f] - w[r54b] + w[r55f] - w[r55b] + w[r56f] - w[r56b] -
               w[r57f] + w[r57b] + w[r58f] - w[r58b] + w[r61f] - w[r61b] +
               w[r72f] - w[r72b] + w[r78f] - w[r78b] - w[r81f] + w[r81b] -
               w[r92f] + w[r92b] - 0x1p+1 * w[r94f] + 0x1p+1 * w[r94b] +
               0x1p+1 * w[r95f] - 0x1p+1 * w[r95b] + 0x1p+1 * w[r96f] -
               0x1p+1 * w[r96b] - w[r97f] + w[r97b] - w[r98f] + w[r98b] -
               w[r99f] + w[r99b] - w[r100f] + w[r100b] - w[r108f] + w[r108b] +
               w[r111f] - w[r111b] - w[r121f] + w[r121b] + w[r122f] - w[r122b] -
               w[r129f] + w[r129b] + w[r141f] - w[r141b] - w[r149f] + w[r149b] +
               w[r154f] - w[r154b] + w[r157f] - w[r157b] + w[r170f] - w[r170b] +
               w[r174f] - w[r174b] + w[r175f] - w[r175b] - w[r182f] + w[r182b] -
               w[r183f] + w[r183b] + w[r188f] - w[r188b] + w[r193f] - w[r193b] -
               w[r198f] + w[r198b] + w[r201f] - w[r201b] + w[r203f] - w[r203b] -
               w[r214f] + w[r214b] - w[r215f] + w[r215b] - w[r216f] + w[r216b] +
               w[r221f] - w[r221b] + w[r224f] - w[r224b] + w[r231f] - w[r231b] +
               w[r239f] - w[r239b] - w[r242f] + w[r242b] + w[r248f] - w[r248b] +
               w[r253f] - w[r253b] - w[r259f] + w[r259b] + w[r262f] - w[r262b] -
               w[r288f] + w[r288b];

  cdot[sO] =
      w[r1f] - w[r1b] - w[r2f] + w[r2b] - w[r4f] + w[r4b] - 0x1p+1 * w[r8f] +
      0x1p+1 * w[r8b] - 0x1p+1 * w[r9f] + 0x1p+1 * w[r9b] - 0x1p+1 * w[r10f] +
      0x1p+1 * w[r10b] - w[r11f] + w[r11b] - w[r16f] + w[r16b] - w[r23f] +
      w[r23b] - w[r26f] + w[r26b] + w[r27f] - w[r27b] - w[r33f] + w[r33b] -
      w[r35f] + w[r35b] - w[r43f] + w[r43b] - w[r48f] + w[r48b] + w[r49f] -
      w[r49b] - w[r55f] + w[r55b] - w[r62f] + w[r62b] - w[r73f] + w[r73b] -
      w[r86f] + w[r86b] - w[r104f] + w[r104b] - w[r111f] + w[r111b] - w[r115f] +
      w[r115b] - w[r122f] + w[r122b] - w[r124f] + w[r124b] - w[r127f] +
      w[r127b] + w[r133f] - w[r133b] - w[r135f] + w[r135b] - w[r136f] +
      w[r136b] - w[r137f] + w[r137b] - w[r144f] + w[r144b] - w[r155f] +
      w[r155b] - w[r167f] + w[r167b] - w[r168f] + w[r168b] - w[r178f] +
      w[r178b] - w[r179f] + w[r179b] - w[r187f] + w[r187b] - w[r190f] +
      w[r190b] - w[r195f] + w[r195b] - w[r211f] + w[r211b] - w[r212f] +
      w[r212b] - w[r213f] + w[r213b] - w[r229f] + w[r229b] - w[r243f] +
      w[r243b] - w[r257f] + w[r257b] - w[r290f] + w[r290b];

  cdot[sCH4] = w[r38f] - w[r38b] + w[r47f] - w[r47b] + w[r53f] - w[r53b] -
               w[r54f] + w[r54b] - w[r55f] + w[r55b] - w[r56f] + w[r56b] +
               w[r57f] - w[r57b] - w[r58f] + w[r58b] + w[r92f] - w[r92b] -
               w[r95f] + w[r95b] - w[r96f] + w[r96b] + w[r108f] - w[r108b] +
               w[r121f] - w[r121b] + w[r129f] - w[r129b] + w[r182f] - w[r182b] +
               w[r183f] - w[r183b] + w[r214f] - w[r214b] + w[r215f] - w[r215b] +
               w[r216f] - w[r216b] + w[r242f] - w[r242b] + w[r259f] - w[r259b] +
               w[r288f] - w[r288b];

  cdot[sOH] =
      w[r1f] - w[r1b] + w[r2f] - w[r2b] - w[r3f] + w[r3b] + 0x1p+1 * w[r4f] -
      0x1p+1 * w[r4b] + w[r11f] - w[r11b] - w[r12f] + w[r12b] +
      0x1p+1 * w[r15f] - 0x1p+1 * w[r15b] + w[r16f] - w[r16b] - w[r17f] +
      w[r17b] + 0x1p+1 * w[r20f] - 0x1p+1 * w[r20b] + w[r21f] - w[r21b] +
      w[r23f] - w[r23b] - w[r24f] + w[r24b] - w[r25f] + w[r25b] + w[r28f] -
      w[r28b] - w[r29f] + w[r29b] + w[r33f] - w[r33b] - w[r34f] + w[r34b] +
      w[r36f] - w[r36b] + w[r43f] - w[r43b] - w[r44f] + w[r44b] + w[r50f] -
      w[r50b] + w[r51f] - w[r51b] + w[r55f] - w[r55b] - w[r56f] + w[r56b] +
      w[r61f] - w[r61b] + w[r62f] - w[r62b] - w[r63f] + w[r63b] + w[r72f] -
      w[r72b] + w[r73f] - w[r73b] - w[r74f] + w[r74b] - w[r81f] + w[r81b] +
      w[r86f] - w[r86b] - w[r87f] + w[r87b] - w[r88f] + w[r88b] - w[r97f] +
      w[r97b] - w[r98f] + w[r98b] + w[r104f] - w[r104b] - w[r105f] + w[r105b] -
      w[r120f] + w[r120b] - w[r123f] + w[r123b] + w[r124f] - w[r124b] +
      w[r136f] - w[r136b] - w[r138f] + w[r138b] - w[r139f] + w[r139b] -
      w[r140f] + w[r140b] - w[r141f] + w[r141b] + w[r142f] - w[r142b] +
      w[r145f] - w[r145b] - w[r146f] + w[r146b] - w[r151f] + w[r151b] -
      w[r156f] + w[r156b] + w[r158f] - w[r158b] + w[r159f] - w[r159b] -
      w[r169f] + w[r169b] + w[r171f] - w[r171b] - w[r176f] + w[r176b] -
      w[r177f] + w[r177b] + w[r178f] - w[r178b] + w[r179f] - w[r179b] +
      w[r190f] - w[r190b] - w[r191f] + w[r191b] - w[r192f] + w[r192b] -
      w[r196f] + w[r196b] + w[r197f] - w[r197b] + w[r199f] - w[r199b] -
      w[r205f] + w[r205b] - w[r206f] + w[r206b] - w[r207f] + w[r207b] +
      w[r211f] - w[r211b] + w[r212f] - w[r212b] + w[r213f] - w[r213b] -
      w[r226f] + w[r226b] + w[r229f] - w[r229b] + 0x1p+1 * w[r232f] -
      0x1p+1 * w[r232b] - w[r233f] + w[r233b] + w[r236f] - w[r236b] - w[r237f] +
      w[r237b] + w[r238f] - w[r238b] - w[r240f] + w[r240b] + w[r243f] -
      w[r243b] + w[r251f] - w[r251b] - w[r255f] + w[r255b] + w[r257f] -
      w[r257b] + w[r268f] - w[r268b] + w[r272f] - w[r272b] + w[r274f] -
      w[r274b] + w[r275f] - w[r275b] - w[r280f] + w[r280b] + w[r283f] -
      w[r283b] - w[r284f] + w[r284b] - w[r285f] + w[r285f] - w[r285b] +
      w[r285b] + w[r287f] - w[r287b] + w[r288f] - w[r288b] + w[r289f] -
      w[r289b] + 0x1p+1 * w[r290f] - 0x1p+1 * w[r290b];

  cdot[sH2O] = w[r3f] - w[r3b] - w[r4f] + w[r4b] + w[r12f] - w[r12b] + w[r17f] -
               w[r17b] + w[r21f] - w[r21b] + w[r24f] - w[r24b] + w[r25f] -
               w[r25b] + w[r34f] - w[r34b] + w[r44f] - w[r44b] + w[r56f] -
               w[r56b] + w[r63f] - w[r63b] + w[r74f] - w[r74b] + w[r87f] -
               w[r87b] + w[r88f] - w[r88b] + w[r97f] - w[r97b] + w[r98f] -
               w[r98b] + w[r101f] - w[r101b] - w[r102f] + w[r102b] + w[r105f] -
               w[r105b] + w[r120f] - w[r120b] + w[r123f] - w[r123b] + w[r140f] -
               w[r140b] - w[r163f] + w[r163f] - w[r163b] + w[r163b] + w[r172f] -
               w[r172b] + w[r176f] - w[r176b] + w[r177f] - w[r177b] + w[r191f] -
               w[r191b] + w[r196f] - w[r196b] + w[r204f] - w[r204b] + w[r205f] -
               w[r205b] + w[r206f] - w[r206b] + w[r207f] - w[r207b] + w[r225f] -
               w[r225b] + w[r226f] - w[r226b] + w[r230f] - w[r230b] + w[r233f] -
               w[r233b] + w[r240f] - w[r240b] + w[r255f] - w[r255b] + w[r281f] -
               w[r281b] + w[r284f] - w[r284b] + w[r285f] - w[r285b];

  cdot[sC2H] = w[r136f] - w[r136b] + w[r140f] - w[r140b] - w[r150f] + w[r150b] -
               w[r151f] + w[r151b] - w[r152f] + w[r152b] - w[r153f] + w[r153b];

  cdot[sC2H2] = w[r116f] - w[r116b] + w[r123f] - w[r123b] + w[r126f] -
                w[r126b] + w[r129f] - w[r129b] + w[r130f] - w[r130b] +
                w[r132f] - w[r132b] - w[r134f] + w[r134b] - w[r135f] +
                w[r135b] - w[r136f] + w[r136b] - w[r137f] + w[r137b] -
                w[r138f] + w[r138b] - w[r139f] + w[r139b] - w[r140f] +
                w[r140b] - w[r141f] + w[r141b] - w[r142f] + w[r142b] +
                w[r148f] - w[r148b] + w[r150f] - w[r150b] + w[r153f] -
                w[r153b] + w[r161f] - w[r161b];

  cdot[sC2H3] = -w[r118f] + w[r118b] + w[r119f] - w[r119b] + w[r120f] -
                w[r120b] + w[r121f] - w[r121b] - w[r123f] + w[r123b] +
                w[r124f] - w[r124b] + w[r125f] - w[r125b] - w[r126f] +
                w[r126b] - w[r127f] + w[r127b] - w[r128f] + w[r128b] -
                w[r129f] + w[r129b] - 0x1p+1 * w[r130f] + 0x1p+1 * w[r130b] -
                w[r131f] + w[r131b] - w[r132f] + w[r132b] - w[r133f] +
                w[r133b] + w[r134f] - w[r134b] + w[r147f] - w[r147b];

  cdot[sCO] = -w[r26f] + w[r26b] - w[r27f] + w[r27b] - w[r28f] + w[r28b] -
              w[r29f] + w[r29b] + w[r30f] - w[r30b] + w[r31f] - w[r31b] +
              w[r32f] - w[r32b] + w[r33f] - w[r33b] + w[r34f] - w[r34b] +
              0x1p+1 * w[r37f] - 0x1p+1 * w[r37b] + w[r38f] - w[r38b] +
              w[r39f] - w[r39b] + w[r41f] - w[r41b] + w[r67f] - w[r67b] -
              w[r78f] + w[r78b] + w[r79f] - w[r79b] + w[r114f] - w[r114b] +
              w[r137f] - w[r137b] + w[r141f] - w[r141b] + w[r143f] - w[r143b] +
              0x1p+1 * w[r144f] - 0x1p+1 * w[r144b] + 0x1p+1 * w[r145f] -
              0x1p+1 * w[r145b] + w[r146f] - w[r146b] + w[r147f] - w[r147b] +
              0x1p+1 * w[r148f] - 0x1p+1 * w[r148b] + w[r149f] - w[r149b] +
              w[r152f] - w[r152b] - w[r160f] + w[r160b] - w[r164f] + w[r164f] -
              w[r164b] + w[r164b] + w[r167f] - w[r167b] + w[r171f] - w[r171b] +
              w[r172f] - w[r172b] + w[r173f] - w[r173b] + w[r175f] - w[r175b] +
              w[r188f] - w[r188b] + w[r192f] - w[r192b] + w[r197f] - w[r197b] +
              w[r198f] - w[r198b] + w[r201f] - w[r201b] - w[r223f] + w[r223b] +
              w[r261f] - w[r261b] + w[r277f] - w[r277b] + w[r281f] - w[r281b] +
              w[r285f] - w[r285b] + w[r287f] - w[r287b] + w[r288f] - w[r288b] +
              w[r289f] - w[r289b] + w[r290f] - w[r290b];

  cdot[sN2] = 0.0;

  cdot[sC2H4] = w[r99f] - w[r99b] + w[r100f] - w[r100b] + w[r110f] - w[r110b] +
                w[r112f] - w[r112b] + w[r113f] - w[r113b] - w[r116f] +
                w[r116b] - w[r117f] + w[r117b] + w[r118f] - w[r118b] -
                w[r119f] + w[r119b] - w[r120f] + w[r120b] - w[r121f] +
                w[r121b] - w[r122f] + w[r122b] - w[r124f] + w[r124b] -
                w[r125f] + w[r125b] + w[r128f] - w[r128b] + w[r130f] -
                w[r130b] + w[r149f] - w[r149b] + w[r204f] - w[r204b] +
                w[r225f] - w[r225b] + w[r230f] - w[r230b] - w[r237f] + w[r237b];

  cdot[sHCO] = -w[r30f] + w[r30b] - w[r31f] + w[r31b] - w[r32f] + w[r32b] -
               w[r33f] + w[r33b] - w[r34f] + w[r34b] - w[r35f] + w[r35b] -
               w[r36f] + w[r36b] - 0x1p+1 * w[r37f] + 0x1p+1 * w[r37b] -
               w[r38f] + w[r38b] - 0x1p+1 * w[r39f] + 0x1p+1 * w[r39b] +
               w[r40f] - w[r40b] + w[r42f] - w[r42b] + w[r43f] - w[r43b] +
               w[r44f] - w[r44b] + w[r45f] - w[r45b] + w[r46f] - w[r46b] +
               w[r47f] - w[r47b] - w[r67f] + w[r67b] - w[r68f] + w[r68b] -
               w[r79f] + w[r79b] - w[r90f] + w[r90b] - w[r114f] + w[r114b] +
               w[r122f] - w[r122b] + w[r131f] - w[r131b] + w[r146f] - w[r146b] +
               w[r152f] - w[r152b] + w[r155f] - w[r155b] + w[r158f] - w[r158b] +
               w[r168f] - w[r168b] + w[r174f] - w[r174b] + w[r193f] - w[r193b] +
               w[r195f] - w[r195b] + w[r199f] - w[r199b] + w[r250f] - w[r250b] +
               w[r265f] - w[r265b] + w[r283f] - w[r283b];

  cdot[sC2H5] = w[r94f] - w[r94b] + w[r103f] - w[r103b] + w[r104f] - w[r104b] +
                w[r105f] - w[r105b] + w[r106f] - w[r106b] + w[r107f] -
                w[r107b] + w[r108f] - w[r108b] - w[r109f] + w[r109b] -
                w[r110f] + w[r110b] - w[r111f] + w[r111b] - w[r112f] +
                w[r112b] - 0x1p+1 * w[r113f] + 0x1p+1 * w[r113b] - w[r114f] +
                w[r114b] - w[r115f] + w[r115b] + w[r117f] - w[r117b] +
                w[r198f] - w[r198b] + w[r223f] - w[r223b] - w[r238f] + w[r238b];

  cdot[sCH2O] =
      w[r39f] - w[r39b] - w[r40f] + w[r40b] - w[r41f] + w[r41b] - w[r42f] +
      w[r42b] - w[r43f] + w[r43b] - w[r44f] + w[r44b] - w[r45f] + w[r45b] -
      w[r46f] + w[r46b] - w[r47f] + w[r47b] + w[r48f] - w[r48b] + w[r50f] -
      w[r50b] + w[r59f] - w[r59b] + w[r60f] - w[r60b] + w[r62f] - w[r62b] +
      w[r63f] - w[r63b] + w[r64f] - w[r64b] + w[r65f] - w[r65b] + w[r66f] -
      w[r66b] + 0x1p+1 * w[r68f] - 0x1p+1 * w[r68b] + w[r69f] - w[r69b] +
      w[r70f] - w[r70b] + w[r71f] - w[r71b] + w[r73f] - w[r73b] + w[r74f] -
      w[r74b] + w[r75f] - w[r75b] + w[r76f] - w[r76b] + w[r77f] - w[r77b] +
      w[r80f] - w[r80b] + w[r90f] - w[r90b] + w[r111f] - w[r111b] + w[r131f] -
      w[r131b] + w[r156f] - w[r156b] + w[r159f] - w[r159b] + w[r169f] -
      w[r169b] + w[r173f] - w[r173b] + w[r195f] - w[r195b] + w[r197f] -
      w[r197b] + w[r199f] - w[r199b] + w[r221f] - w[r221b] + 0x1p+1 * w[r236f] -
      0x1p+1 * w[r236b] + w[r248f] - w[r248b] + w[r249f] - w[r249b] - w[r250f] +
      w[r250b] - w[r265f] + w[r265b] + w[r269f] - w[r269b] + 0x1p+1 * w[r272f] -
      0x1p+1 * w[r272b] - w[r280f] + w[r280b];

  cdot[sC2H6] = w[r52f] - w[r52b] - w[r103f] + w[r103b] - w[r104f] + w[r104b] -
                w[r105f] + w[r105b] - w[r106f] + w[r106b] - w[r107f] +
                w[r107b] - w[r108f] + w[r108b] + w[r109f] - w[r109b] +
                w[r113f] - w[r113b] + w[r114f] - w[r114b];

  cdot[sCH2OH] =
      -w[r59f] + w[r59b] - w[r60f] + w[r60b] - w[r61f] + w[r61b] - w[r62f] +
      w[r62b] - w[r63f] + w[r63b] - w[r64f] + w[r64b] - w[r65f] + w[r65b] -
      w[r66f] + w[r66b] - w[r67f] + w[r67b] - w[r68f] + w[r68b] -
      0x1p+1 * w[r69f] + 0x1p+1 * w[r69b] - w[r70f] + w[r70b] - w[r82f] +
      w[r82b] + w[r84f] - w[r84b] + w[r86f] - w[r86b] + w[r88f] - w[r88b] +
      w[r89f] - w[r89b] + w[r90f] - w[r90b] + w[r91f] - w[r91b] + w[r92f] -
      w[r92b] + w[r93f] - w[r93b] + w[r192f] - w[r192b] + w[r203f] - w[r203b] +
      w[r224f] - w[r224b] + w[r231f] - w[r231b] + w[r278f] - w[r278b];

  cdot[sCH3O] = w[r49f] - w[r49b] + w[r51f] - w[r51b] - w[r70f] + w[r70b] -
                w[r71f] + w[r71b] - w[r72f] + w[r72b] - w[r73f] + w[r73b] -
                w[r74f] + w[r74b] - w[r75f] + w[r75b] - w[r76f] + w[r76b] -
                w[r77f] + w[r77b] - w[r78f] + w[r78b] - w[r79f] + w[r79b] -
                0x1p+1 * w[r80f] + 0x1p+1 * w[r80b] - w[r83f] + w[r83b] +
                w[r85f] - w[r85b] + w[r87f] - w[r87b] - w[r93f] + w[r93b] -
                w[r101f] + w[r101b] + w[r239f] - w[r239b] - w[r246f] +
                w[r246b] - w[r249f] + w[r249b] - w[r260f] + w[r260b] +
                w[r261f] - w[r261b] + w[r269f] - w[r269b];

  cdot[sO2] = -w[r1f] + w[r1b] + w[r8f] - w[r8b] + w[r9f] - w[r9b] + w[r10f] -
              w[r10b] - w[r13f] + w[r13b] + w[r14f] - w[r14b] + w[r16f] -
              w[r16b] + w[r17f] - w[r17b] + w[r18f] - w[r18b] + w[r19f] -
              w[r19b] - w[r27f] + w[r27b] - w[r31f] + w[r31b] - w[r45f] +
              w[r45b] - w[r49f] + w[r49b] - w[r50f] + w[r50b] + w[r57f] -
              w[r57b] - w[r64f] + w[r64b] - w[r65f] + w[r65b] - w[r75f] +
              w[r75b] - w[r76f] + w[r76b] - w[r89f] + w[r89b] - w[r106f] +
              w[r106b] - w[r112f] + w[r112b] - w[r125f] + w[r125b] - w[r131f] +
              w[r131b] - w[r132f] + w[r132b] - w[r133f] + w[r133b] - w[r145f] +
              w[r145b] - w[r152f] + w[r152b] - w[r158f] + w[r158b] - w[r171f] +
              w[r171b] - w[r172f] + w[r172b] - w[r186f] + w[r186b] - w[r197f] +
              w[r197b] + w[r200f] - w[r200b] - w[r222f] + w[r222b] - w[r227f] +
              w[r227b] - w[r228f] + w[r228b] - w[r235f] + w[r235b] - w[r245f] +
              w[r245b] - w[r254f] + w[r254b] - w[r264f] + w[r264b] + w[r266f] -
              w[r266b] + w[r267f] - w[r267b] - w[r270f] + w[r270b] - w[r273f] +
              w[r273b];

  cdot[sCH3OH] = w[r67f] - w[r67b] + w[r69f] - w[r69b] + w[r70f] - w[r70b] +
                 w[r79f] - w[r79b] + w[r80f] - w[r80b] + w[r81f] - w[r81b] +
                 w[r82f] - w[r82b] + w[r83f] - w[r83b] - w[r84f] + w[r84b] -
                 w[r85f] + w[r85b] - w[r86f] + w[r86b] - w[r87f] + w[r87b] -
                 w[r88f] + w[r88b] - w[r89f] + w[r89b] - w[r90f] + w[r90b] -
                 w[r91f] + w[r91b] - w[r92f] + w[r92b] - w[r93f] + w[r93f] -
                 w[r93b] + w[r93b] + w[r102f] - w[r102b] + w[r246f] - w[r246b] +
                 w[r260f] - w[r260b];

  cdot[sHO2] = w[r13f] - w[r13b] - w[r14f] + w[r14b] - w[r15f] + w[r15b] -
               w[r16f] + w[r16b] - w[r17f] + w[r17b] - 0x1p+1 * w[r18f] +
               0x1p+1 * w[r18b] - 0x1p+1 * w[r19f] + 0x1p+1 * w[r19b] +
               w[r22f] - w[r22b] + w[r23f] - w[r23b] + w[r24f] - w[r24b] +
               w[r25f] - w[r25b] - w[r28f] + w[r28b] + w[r31f] - w[r31b] -
               w[r36f] + w[r36b] + w[r45f] - w[r45b] - w[r46f] + w[r46b] -
               w[r51f] + w[r51b] - w[r57f] + w[r57b] - w[r58f] + w[r58b] +
               w[r64f] - w[r64b] + w[r65f] - w[r65b] - w[r66f] + w[r66b] +
               w[r75f] - w[r75b] + w[r76f] - w[r76b] - w[r77f] + w[r77b] +
               w[r89f] - w[r89b] - w[r91f] + w[r91b] + w[r106f] - w[r106b] -
               w[r107f] + w[r107b] + w[r112f] - w[r112b] + w[r125f] - w[r125b] +
               w[r128f] - w[r128b] + w[r132f] - w[r132b] - w[r142f] + w[r142b] -
               w[r159f] + w[r159b] - w[r184f] + w[r184b] - w[r185f] + w[r185b] +
               w[r186f] - w[r186b] - w[r199f] + w[r199b] - w[r200f] + w[r200b] -
               w[r217f] + w[r217b] - w[r218f] + w[r218b] - w[r219f] + w[r219b] +
               w[r222f] - w[r222b] + w[r227f] - w[r227b] + w[r228f] - w[r228b] -
               w[r232f] + w[r232b] - w[r238f] + w[r238b] - w[r244f] + w[r244b] +
               w[r245f] - w[r245b] - w[r251f] + w[r251b] + w[r254f] - w[r254b] -
               w[r256f] + w[r256b] + w[r270f] - w[r270b] - w[r289f] + w[r289b];

  cdot[sH2O2] = w[r18f] - w[r18b] + w[r19f] - w[r19b] - w[r20f] + w[r20b] -
                w[r21f] + w[r21b] - w[r22f] + w[r22b] - w[r23f] + w[r23b] -
                w[r24f] + w[r24b] - w[r25f] + w[r25b] + w[r46f] - w[r46b] +
                w[r58f] - w[r58b] + w[r66f] - w[r66b] + w[r77f] - w[r77b] +
                w[r91f] - w[r91b] + w[r107f] - w[r107b] - w[r128f] + w[r128b] +
                w[r184f] - w[r184b] + w[r185f] - w[r185b] + w[r217f] -
                w[r217b] + w[r218f] - w[r218b] + w[r219f] - w[r219b] +
                w[r244f] - w[r244b] + w[r256f] - w[r256b] + w[r289f] - w[r289b];

  cdot[sHCCO] = w[r135f] - w[r135b] - w[r143f] + w[r143b] - w[r144f] +
                w[r144b] - w[r145f] + w[r145b] - w[r146f] + w[r146b] -
                w[r147f] + w[r147b] - 0x1p+1 * w[r148f] + 0x1p+1 * w[r148b] -
                w[r149f] + w[r149b] + w[r151f] - w[r151b] + w[r189f] -
                w[r189b] + w[r190f] - w[r190b] + w[r191f] - w[r191b];

  cdot[sCH2CO] = w[r127f] - w[r127b] + w[r138f] - w[r138b] + w[r142f] -
                 w[r142b] + w[r160f] - w[r160b] - w[r187f] + w[r187b] -
                 w[r188f] + w[r188b] - w[r189f] + w[r189b] - w[r190f] +
                 w[r190b] - w[r191f] + w[r191b] - w[r192f] + w[r192b] +
                 w[r194f] - w[r194b] + w[r196f] - w[r196b] + w[r202f] -
                 w[r202b];

  cdot[sHCCOH] = w[r139f] - w[r139b];

  cdot[sCH2HCO] = w[r133f] - w[r133b] + w[r177f] - w[r177b] + w[r179f] -
                  w[r179b] + w[r181f] - w[r181b] + w[r183f] - w[r183b] +
                  w[r185f] - w[r185b] - w[r193f] + w[r193b] - w[r194f] +
                  w[r194b] - w[r195f] + w[r195b] - w[r196f] + w[r196b] -
                  w[r197f] + w[r197b] - w[r198f] + w[r198b] - w[r199f] +
                  w[r199b] - w[r200f] + w[r200b] - w[r201f] + w[r201b] -
                  w[r202f] + w[r202b];

  cdot[sCH3CO] = -w[r175f] + w[r175b] + w[r176f] - w[r176b] + w[r178f] -
                 w[r178b] + w[r180f] - w[r180b] + w[r182f] - w[r182b] +
                 w[r184f] - w[r184b] + w[r186f] - w[r186b];

  cdot[sCO2] = w[r26f] - w[r26b] + w[r27f] - w[r27b] + w[r28f] - w[r28b] +
               w[r29f] - w[r29b] + w[r35f] - w[r35b] + w[r36f] - w[r36b] +
               w[r78f] - w[r78b] - w[r165f] + w[r165f] - w[r165b] + w[r165b] -
               w[r173f] + w[r173b] + w[r187f] - w[r187b] + w[r223f] - w[r223b] +
               w[r262f] - w[r262b] + w[r263f] - w[r263b] + w[r278f] - w[r278b] +
               w[r282f] - w[r282b] + w[r284f] - w[r284b] + w[r286f] - w[r286b];

  cdot[sCH3HCO] =
      w[r115f] - w[r115b] - w[r174f] + w[r174b] - w[r176f] + w[r176b] -
      w[r177f] + w[r177b] - w[r178f] + w[r178b] - w[r179f] + w[r179b] -
      w[r180f] + w[r180b] - w[r181f] + w[r181b] - w[r182f] + w[r182b] -
      w[r183f] + w[r183b] - w[r184f] + w[r184b] - w[r185f] + w[r185b] -
      w[r186f] + w[r186b] + w[r200f] - w[r200b] + w[r220f] - w[r220b] +
      w[r222f] - w[r222b] + w[r226f] - w[r226b] + w[r227f] - w[r227b] +
      w[r228f] - w[r228b] + w[r229f] - w[r229b] + w[r232f] - w[r232b] +
      w[r233f] - w[r233b] + w[r234f] - w[r234b];

  cdot[sOCHO] = w[r253f] - w[r253b] - w[r263f] + w[r263b];

  cdot[sCH3CHOH] = w[r206f] - w[r206b] + w[r209f] - w[r209b] + w[r212f] -
                   w[r212b] + w[r215f] - w[r215b] + w[r217f] - w[r217b] -
                   w[r227f] + w[r227b] - w[r228f] + w[r228b] - w[r229f] +
                   w[r229b] - w[r230f] + w[r230b] - w[r231f] + w[r231b] -
                   w[r232f] + w[r232b] - w[r233f] + w[r233b] - w[r234f] +
                   w[r234b];

  cdot[sC2H4OH] = w[r205f] - w[r205b] + w[r208f] - w[r208b] + w[r211f] -
                  w[r211b] + w[r214f] - w[r214b] + w[r218f] - w[r218b] -
                  w[r235f] + w[r235b] + w[r237f] - w[r237b];

  cdot[sCH3CH2O] = w[r207f] - w[r207b] + w[r210f] - w[r210b] + w[r213f] -
                   w[r213b] + w[r216f] - w[r216b] + w[r219f] - w[r219b] -
                   w[r220f] + w[r220b] - w[r221f] + w[r221b] - w[r222f] +
                   w[r222b] - w[r223f] + w[r223b] - w[r224f] + w[r224b] -
                   w[r225f] + w[r225b] - w[r226f] + w[r226b] + w[r238f] -
                   w[r238b];

  cdot[sCH3OCH2] = w[r240f] - w[r240b] + w[r241f] - w[r241b] + w[r242f] -
                   w[r242b] + w[r243f] - w[r243b] + w[r244f] - w[r244b] +
                   w[r245f] - w[r245b] + w[r246f] - w[r246b] + w[r247f] -
                   w[r247b] - w[r248f] + w[r248b] - w[r249f] + w[r249b] -
                   w[r250f] + w[r250b] - w[r251f] + w[r251b] - w[r264f] +
                   w[r264b];

  cdot[sHCOOH] = w[r279f] - w[r279b] - w[r281f] + w[r281b] - w[r282f] +
                 w[r282b] - w[r283f] + w[r283b] - w[r284f] + w[r284b] -
                 w[r285f] + w[r285b] - w[r286f] + w[r286b] - w[r287f] +
                 w[r287b] - w[r288f] + w[r288b] - w[r289f] + w[r289b] -
                 w[r290f] + w[r290b];

  cdot[sCH3OCH3] = -w[r239f] + w[r239b] - w[r240f] + w[r240b] - w[r241f] +
                   w[r241b] - w[r242f] + w[r242b] - w[r243f] + w[r243b] -
                   w[r244f] + w[r244b] - w[r245f] + w[r245b] - w[r246f] +
                   w[r246b] - w[r247f] + w[r247b] + w[r249f] - w[r249b] +
                   w[r250f] - w[r250b];

  cdot[sC2H5OH] = -w[r203f] + w[r203b] - w[r204f] + w[r204b] - w[r205f] +
                  w[r205b] - w[r206f] + w[r206b] - w[r207f] + w[r207b] -
                  w[r208f] + w[r208b] - w[r209f] + w[r209b] - w[r210f] +
                  w[r210b] - w[r211f] + w[r211b] - w[r212f] + w[r212b] -
                  w[r213f] + w[r213b] - w[r214f] + w[r214b] - w[r215f] +
                  w[r215b] - w[r216f] + w[r216b] - w[r217f] + w[r217b] -
                  w[r218f] + w[r218b] - w[r219f] + w[r219b];

  cdot[sHOCH2O] =
      w[r277f] - w[r277b] - w[r279f] + w[r279b] + w[r280f] - w[r280b];

  cdot[sCH3OCO] = w[r254f] - w[r254b] + w[r255f] - w[r255b] + w[r256f] -
                  w[r256b] + w[r257f] - w[r257b] + w[r258f] - w[r258b] +
                  w[r259f] - w[r259b] + w[r260f] - w[r260b] - w[r261f] +
                  w[r261b] - w[r262f] + w[r262b];

  cdot[sCH3OCHO] = w[r252f] - w[r252b] - w[r253f] + w[r253b] - w[r254f] +
                   w[r254b] - w[r255f] + w[r255b] - w[r256f] + w[r256b] -
                   w[r257f] + w[r257b] - w[r258f] + w[r258b] - w[r259f] +
                   w[r259b] - w[r260f] + w[r260b] + w[r267f] - w[r267b] +
                   w[r270f] - w[r270b];

  cdot[sCH3OCH2O] = w[r251f] - w[r251b] - w[r252f] + w[r252b] +
                    0x1p+1 * w[r266f] - 0x1p+1 * w[r266b] + w[r268f] -
                    w[r268b] - w[r269f] + w[r269b] - w[r270f] + w[r270b];

  cdot[sCH3OCH2OH] = w[r267f] - w[r267b];

  cdot[sOCH2OCHO] = w[r275f] - w[r275b] - w[r276f] + w[r276b];

  cdot[sHOCH2OCO] =
      w[r276f] - w[r276b] - w[r277f] + w[r277b] - w[r278f] + w[r278b];

  cdot[sHOC2H4O2] = w[r235f] - w[r235b] - w[r236f] + w[r236b];

  cdot[sCH3OCH2O2] = -w[r247f] + w[r247b] + w[r264f] - w[r264b] - w[r265f] +
                     w[r265b] - 0x1p+1 * w[r266f] + 0x1p+1 * w[r266b] -
                     0x1p+1 * w[r267f] + 0x1p+1 * w[r267b] - w[r271f] +
                     w[r271b];

  cdot[sCH2OCH2O2H] =
      w[r271f] - w[r271b] - w[r272f] + w[r272b] - w[r273f] + w[r273b];

  cdot[sCH3OCH2O2H] =
      w[r247f] - w[r247b] + w[r265f] - w[r265b] - w[r268f] + w[r268b];

  cdot[sHO2CH2OCHO] = w[r274f] - w[r274b] - w[r275f] + w[r275b];

  cdot[sO2CH2OCH2O2H] = w[r273f] - w[r273b] - w[r274f] + w[r274b];

  cdot[sAR] = -w[r6f] + w[r6f] - w[r6b] + w[r6b] - w[r9f] + w[r9f] - w[r9b] +
              w[r9b] - w[r166f] + w[r166f] - w[r166b] + w[r166b];

  cdot[sHE] = -w[r7f] + w[r7f] - w[r7b] + w[r7b] - w[r10f] + w[r10f] - w[r10b] +
              w[r10b];
}

double GetLindRateCoeff(double /* temp */, double /* pressure */, double k0,
                        double kInf, double fc, double conc) {
  // const double R = 8314.34; /* [J / kmole K] */
  double Ntmp;
  double kl;
  double f;
  double cCoeff, DCoeff, log10kNull;

  int iTroe = 1;

  if (std::isinf(conc)) {
    /*conc = pressure / ( R * temp );*/
    return 0.0;
  }
  Ntmp = 0.75 - 1.27 * (fc < 1e-200 ? -200 : log10(fc));
  if (iTroe) {
    cCoeff = -0.4 - 0.67 * (fc < 1e-200 ? -200 : log10(fc));
    DCoeff = 0.14;
    k0 *= conc / MAX_C(kInf, 1.0e-60);
    log10kNull = log10(MAX_C(k0, 1.0e-60));
    f = (log10kNull + cCoeff) / (Ntmp - DCoeff * (log10kNull + cCoeff));
    f = pow(fc, 1.0 / (f * f + 1.0));
    kInf *= f * k0 / (1.0 + k0);
  } else {
    k0 = k0 * conc / kInf;
    kl = k0 / (1.0 + k0);
    f = (k0 < 1e-200 ? -200 : log10(k0)) / Ntmp;
    f = pow(fc, 1.0 / (f * f + 1.0));
    kInf = kInf * f * kl;
  }
  return kInf;
}

/*
double CatchZero( double a )
{
        return ( a == 0.0 ) ? 1.0e-20 : a;
}

*/
void Zhao2008Dme::getMolarMass(span<double> W) const {
  W[sH] = 0x1.020c49ba5e354p+0;
  W[sH2] = 0x1.020c49ba5e354p+1;
  W[sCH2] = 0x1.c0d4fdf3b645ap+3;
  W[sCH2S] = 0x1.c0d4fdf3b645ap+3;
  W[sCH3] = 0x1.e116872b020c4p+3;
  W[sO] = 0x1p+4;
  W[sCH4] = 0x1.00ac083126e98p+4;
  W[sOH] = 0x1.1020c49ba5e35p+4;
  W[sH2O] = 0x1.204189374bc6ap+4;
  W[sC2H] = 0x1.9072b020c49bap+4;
  W[sC2H2] = 0x1.a09374bc6a7fp+4;
  W[sC2H3] = 0x1.b0b4395810625p+4;
  W[sCO] = 0x1.c028f5c28f5c2p+4;
  W[sN2] = 0x1.c051eb851eb85p+4;
  W[sC2H4] = 0x1.c0d4fdf3b645ap+4;
  W[sHCO] = 0x1.d049ba5e353f8p+4;
  W[sC2H5] = 0x1.d0f5c28f5c28fp+4;
  W[sCH2O] = 0x1.e06a7ef9db22cp+4;
  W[sC2H6] = 0x1.e116872b020c4p+4;
  W[sCH2OH] = 0x1.f08b439581062p+4;
  W[sCH3O] = 0x1.f08b439581062p+4;
  W[sO2] = 0x1p+5;
  W[sCH3OH] = 0x1.005604189374cp+5;
  W[sHO2] = 0x1.0810624dd2f1bp+5;
  W[sH2O2] = 0x1.1020c49ba5e35p+5;
  W[sHCCO] = 0x1.48395810624ddp+5;
  W[sCH2CO] = 0x1.5049ba5e353f8p+5;
  W[sHCCOH] = 0x1.e170a3d70a3d6p+5;
  W[sCH2HCO] = 0x1.585a1cac08312p+5;
  W[sCH3CO] = 0x1.585a1cac08312p+5;
  W[sCO2] = 0x1.60147ae147ae1p+5;
  W[sCH3HCO] = 0x1.606a7ef9db22dp+5;
  W[sOCHO] = 0x1.6824dd2f1a9fcp+5;
  W[sCH3CHOH] = 0x1.687ae147ae148p+5;
  W[sC2H4OH] = 0x1.687ae147ae148p+5;
  W[sCH3CH2O] = 0x1.687ae147ae148p+5;
  W[sCH3OCH2] = 0x1.687ae147ae148p+5;
  W[sHCOOH] = 0x1.70353f7ced916p+5;
  W[sCH3OCH3] = 0x1.708b439581062p+5;
  W[sC2H5OH] = 0x1.708b439581062p+5;
  W[sHOCH2O] = 0x1.7845a1cac0831p+5;
  W[sCH3OCO] = 0x1.d85a1cac08312p+5;
  W[sCH3OCHO] = 0x1.e06a7ef9db22cp+5;
  W[sCH3OCH2O] = 0x1.e87ae147ae148p+5;
  W[sCH3OCH2OH] = 0x1.f08b439581062p+5;
  W[sOCH2OCHO] = 0x1.2c2d0e5604189p+6;
  W[sHOCH2OCO] = 0x1.2c2d0e5604189p+6;
  W[sHOC2H4O2] = 0x1.343d70a3d70a4p+6;
  W[sCH3OCH2O2] = 0x1.343d70a3d70a4p+6;
  W[sCH2OCH2O2H] = 0x1.343d70a3d70a4p+6;
  W[sCH3OCH2O2H] = 0x1.3845a1cac0831p+6;
  W[sHO2CH2OCHO] = 0x1.70353f7ced916p+6;
  W[sO2CH2OCH2O2H] = 0x1.b43d70a3d70a4p+6;
  W[sAR] = 0x1.3f95810624dd3p+5;
  W[sHE] = 0x1.42e00d1b71759p+4;
}

std::vector<std::string> Zhao2008Dme::getSpeciesNames() const {
  return std::vector<std::string>{
      "H",          "H2",         "CH2",          "CH2S",      "CH3",
      "O",          "CH4",        "OH ",          "H2O",       "C2H",
      "C2H2",       "C2H3",       "CO",           "N2",        "C2H4",
      "HCO",        "C2H5",       "CH2O",         "C2H6",      "CH2OH",
      "CH3O",       "O2",         "CH3OH",        "HO2",       "H2O2",
      "HCCO",       "CH2CO",      "HCCOH",        "CH2HCO ",   "CH3CO",
      "CO2",        "CH3HCO ",    "OCHO",         "CH3CHOH",   "C2H4OH ",
      "CH3CH2O",    "CH3OCH2",    "HCOOH",        "CH3OCH3",   "C2H5OH",
      "HOCH2O",     "CH3OCO",     "CH3OCHO",      "CH3OCH2O",  "CH3OCH2OH",
      "OCH2OCHO",   "HOCH2OCO",   "HOC2H4O2",     "CH3OCH2O2", "CH2OCH2O2H",
      "CH3OCH2O2H", "HO2CH2OCHO", "O2CH2OCH2O2H", "AR",        "HE"};
}

void Zhao2008Dme::ComputeThermoData(span<double> h, span<double> cp, double T,
                                    span<double>) const {
  /*
 This function computes enthalpy 'h' and heat capacity 'cp' as
  function of temperature 'T' for all non steady state species
      in units [J/kg] and [J/(kg K)], respectively.
          The parameter h and cp should provide workspace of length 55 */

  int i;
  if (T > 1000.0) {
    h[sH] = 0x1.01c3c6b46b46cp+13 *
            (T * (0x1.4p+1 +
                  T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) +
             0x1.8dfe851eb851fp+14);
    cp[sH] =
        0x1.01c3c6b46b46cp+13 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    h[sH2] = 0x1.01c3c6b46b46cp+12 *
             (T * (0x1.7ee6f2e8c0486p+1 +
                   T * (0x1.6f090d9fefacdp-12 +
                        T * (-0x1.42a0ce7cd092dp-26 +
                             T * (-0x1.44cea998664c2p-39 +
                                  T * 0x1.6cf52ffc4f362p-52)))) -
              0x1.a1845a1cac083p+9);
    cp[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.7ee6f2e8c0486p+1 +
         T * (0x1.6f090d9fefacdp-11 +
              T * (-0x1.e3f135bb38dc3p-25 +
                   T * (-0x1.44cea998664c2p-37 + T * 0x1.c8327bfb6303ap-50))));
    h[sCH2] = 0x1.286500f4d652p+9 *
              (T * (0x1.6fe28bbb5f921p+1 +
                    T * (0x1.df4030068403dp-10 +
                         T * (-0x1.f8480a2ca2df1p-22 +
                              T * (0x1.1e1208519dbd9p-34 +
                                   T * -0x1.0e8b3f67fea4fp-48)))) +
               0x1.696f353f7ced9p+15);
    cp[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.6fe28bbb5f921p+1 +
         T * (0x1.df4030068403dp-9 +
              T * (-0x1.7a3607a17a275p-20 +
                   T * (0x1.1e1208519dbd9p-32 + T * -0x1.522e0f41fe4e3p-46))));
    h[sCH2S] = 0x1.286500f4d652p+9 *
               (T * (0x1.256183d389aa6p+1 +
                     T * (0x1.3120cfb16b35fp-9 +
                          T * (-0x1.680c0914ed793p-21 +
                               T * (0x1.cb7e14e4de1f4p-34 +
                                    T * -0x1.e9953783e1eedp-48)))) +
                0x1.8ddbffd8adabap+15);
    cp[sCH2S] =
        0x1.286500f4d652p+9 *
        (0x1.256183d389aa6p+1 +
         T * (0x1.3120cfb16b35fp-8 +
              T * (-0x1.0e0906cfb21aep-19 +
                   T * (0x1.cb7e14e4de1f4p-32 + T * -0x1.31fd42b26d354p-45))));
    h[sCH3] = 0x1.148599bc79368p+9 *
              (T * (0x1.7d330e4a459e7p+1 +
                    T * (0x1.7bf7d0ba1fd49p-9 +
                         T * (-0x1.618b37453a214p-21 +
                              T * (0x1.51e0ab53b03a6p-34 +
                                   T * -0x1.023798f4a2ap-48)))) +
               0x1.01f60d4fdf3b6p+14);
    cp[sCH3] =
        0x1.148599bc79368p+9 *
        (0x1.7d330e4a459e7p+1 +
         T * (0x1.7bf7d0ba1fd49p-8 +
              T * (-0x1.09286973eb98fp-19 +
                   T * (0x1.51e0ab53b03a6p-32 + T * -0x1.42c57f31cb48p-46))));
    h[sO] = 0x1.03d3adab9f55ap+9 *
            (T * (0x1.456238da3c212p+1 +
                  T * (-0x1.ce39030add37ep-17 +
                       T * (-0x1.1c4c1de9ad9d8p-30 +
                            T * (0x1.4040bc0d8605ep-40 +
                                 T * -0x1.92e1b63899c98p-54)))) +
             0x1.c8bb333333333p+14);
    cp[sO] =
        0x1.03d3adab9f55ap+9 *
        (0x1.456238da3c212p+1 +
         T * (-0x1.ce39030add37ep-16 +
              T * (-0x1.aa722cde846c4p-29 +
                   T * (0x1.4040bc0d8605ep-38 + T * -0x1.f79a23c6c03bdp-52))));
    h[sCH4] = 0x1.0325882935d13p+9 *
              (T * (0x1.aef87ad080b67p+0 +
                    T * (0x1.4f7431802ab2bp-8 +
                         T * (-0x1.5abd9cc5788acp-20 +
                              T * (0x1.750a9e764526bp-33 +
                                   T * -0x1.44817dd558abep-47)))) -
               0x1.3b0651eb851ecp+13);
    cp[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.aef87ad080b67p+0 +
         T * (0x1.4f7431802ab2bp-7 +
              T * (-0x1.040e35941a681p-18 +
                   T * (0x1.750a9e764526bp-31 + T * -0x1.95a1dd4aaed6dp-45))));
    h[sOH] = 0x1.e8db171241376p+8 *
             (T * (0x1.6eaf6f6ecdbe1p+1 +
                   T * (0x1.14f4d0c23ad0ap-11 +
                        T * (-0x1.72ead5c210468p-24 +
                             T * (0x1.0c7922a4a2d04p-37 +
                                  T * -0x1.3320f7b200b9bp-52)))) +
              0x1.cc741eb851eb8p+11);
    cp[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.6eaf6f6ecdbe1p+1 +
         T * (0x1.14f4d0c23ad0ap-10 +
              T * (-0x1.163020518c34ep-22 +
                   T * (0x1.0c7922a4a2d04p-35 + T * -0x1.7fe9359e80e81p-50))));
    h[sH2O] = 0x1.cd8113bbf9d1p+8 *
              (T * (0x1.5608e15011905p+1 +
                    T * (0x1.90982cf6c7e1cp-10 +
                         T * (-0x1.3877da6b0bdcfp-22 +
                              T * (0x1.081a10711c1a4p-35 +
                                   T * -0x1.7073a217736dcp-50)))) -
               0x1.d32cd70a3d70ap+14);
    cp[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.5608e15011905p+1 +
         T * (0x1.90982cf6c7e1cp-9 +
              T * (-0x1.d4b3c7a091cb6p-21 +
                   T * (0x1.081a10711c1a4p-33 + T * -0x1.cc908a9d50493p-48))));
    h[sC2H] = 0x1.4c34d172a6531p+8 *
              (T * (0x1.957aaf1dba502p+1 +
                    T * (0x1.377101463a6fcp-9 +
                         T * (-0x1.48e6585748e57p-21 +
                              T * (0x1.4e75f1b05c789p-34 +
                                   T * -0x1.fed6b3b564915p-49)))) +
               0x1.063110a3d70a4p+16);
    cp[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.957aaf1dba502p+1 +
         T * (0x1.377101463a6fcp-8 +
              T * (-0x1.ed598482ed583p-20 +
                   T * (0x1.4e75f1b05c789p-32 + T * -0x1.3f4630515edadp-46))));
    h[sC2H2] = 0x1.3f584141b0716p+8 *
               (T * (0x1.0971c7ee6baep+2 +
                     T * (0x1.86b42b3f9ab1ep-9 +
                          T * (-0x1.a8a7da8ef5fe7p-21 +
                               T * (0x1.00f66a3bafde2p-33 +
                                    T * -0x1.044c229f7f34fp-47)))) +
                0x1.953fff2e48e8ap+14);
    cp[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.0971c7ee6baep+2 +
         T * (0x1.86b42b3f9ab1ep-8 +
              T * (-0x1.3e7de3eb387edp-19 +
                   T * (0x1.00f66a3bafde2p-31 + T * -0x1.455f2b475f023p-45))));
    h[sC2H3] = 0x1.337122eb1b451p+8 *
               (T * (0x1.8224031487768p+1 +
                     T * (0x1.52803e497ede9p-8 +
                          T * (-0x1.a2d53f39c89fcp-20 +
                               T * (0x1.17b98c3c9686ep-32 +
                                    T * -0x1.36c974e43b3b7p-46)))) +
                0x1.0e69bf6fd21ffp+15);
    cp[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.8224031487768p+1 +
         T * (0x1.52803e497ede9p-7 +
              T * (-0x1.3a1fef6b5677dp-18 +
                   T * (0x1.17b98c3c9686ep-30 + T * -0x1.847bd21d4a0a4p-44))));
    h[sCO] = 0x1.28d6c752d48edp+8 *
             (T * (0x1.8335c182ecaefp+1 +
                   T * (0x1.7a31384b0ee0ep-11 +
                        T * (-0x1.931203ab5b4cp-23 +
                             T * (0x1.bffa067a23237p-36 +
                                  T * -0x1.8e63a678bbe02p-50)))) -
              0x1.bde2ccccccccdp+13);
    cp[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.8335c182ecaefp+1 +
         T * (0x1.7a31384b0ee0ep-10 +
              T * (-0x1.2e4d82c08479p-21 +
                   T * (0x1.bffa067a23237p-34 + T * -0x1.f1fc9016ead82p-48))));
    h[sN2] = 0x1.28bba88de9adap+8 *
             (T * (0x1.769c23b7952d2p+1 +
                   T * (0x1.8610723573f79p-11 +
                        T * (-0x1.96ee58d5ad6b7p-23 +
                             T * (0x1.bc12905f47c3dp-36 +
                                  T * -0x1.854ddeba4197bp-50)))) -
              0x1.cd661b089a027p+9);
    cp[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.769c23b7952d2p+1 +
         T * (0x1.8610723573f79p-10 +
              T * (-0x1.3132c2a042109p-21 +
                   T * (0x1.bc12905f47c3dp-34 + T * -0x1.e6a15668d1fd9p-48))));
    h[sC2H4] = 0x1.286500f4d652p+8 *
               (T * (0x1.049f4a5d9c3d6p+1 +
                     T * (0x1.dfe6a5720731ep-8 +
                          T * (-0x1.2c3c348d37b5bp-19 +
                               T * (0x1.94aeec0be93f2p-32 +
                                    T * -0x1.c4e760753f615p-46)))) +
                0x1.34be2da122fadp+12);
    cp[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.049f4a5d9c3d6p+1 +
         T * (0x1.dfe6a5720731ep-7 +
              T * (-0x1.c25a4ed3d3908p-18 +
                   T * (0x1.94aeec0be93f2p-30 + T * -0x1.1b109c49479cdp-43))));
    h[sHCO] = 0x1.1e87151786dcap+8 *
              (T * (0x1.c754a7f8012ep+1 +
                    T * (0x1.b682cd3e25996p-10 +
                         T * (-0x1.ddd137e584398p-22 +
                              T * (0x1.0fa472a0a4455p-34 +
                                   T * -0x1.edfbe2137d49dp-49)))) +
               0x1.e98a5e353f7cfp+11);
    cp[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.c754a7f8012ep+1 +
         T * (0x1.b682cd3e25996p-9 +
              T * (-0x1.665ce9ec232b2p-20 +
                   T * (0x1.0fa472a0a4455p-32 + T * -0x1.34bd6d4c2e4e2p-46))));
    h[sC2H5] = 0x1.1e1d11b23ceecp+8 *
               (T * (0x1.126ca61b881bcp+2 +
                     T * (0x1.976f0de602e4fp-8 +
                          T * (-0x1.8af33d1c075bp-20 +
                               T * (0x1.846ccefec2441p-33 +
                                    T * -0x1.2ee523bc24286p-47)))) +
                0x1.78c3a3d70a3d7p+13);
    cp[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.126ca61b881bcp+2 +
         T * (0x1.976f0de602e4fp-7 +
              T * (-0x1.28366dd505844p-18 +
                   T * (0x1.846ccefec2441p-31 + T * -0x1.7a9e6cab2d327p-45))));
    h[sCH2O] = 0x1.14e89ec288e65p+8 *
               (T * (0x1.497bf401c4fc2p+2 +
                     T * (0x1.77e3742ab5906p-10 +
                          T * (-0x1.547c4590e3dd3p-24 +
                               T * (-0x1.624a929eb324bp-35 +
                                    T * 0x1.9bb0a660b7c0ap-48)))) -
                0x1.fb31624dd2f1bp+13);
    cp[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.497bf401c4fc2p+2 +
         T * (0x1.77e3742ab5906p-9 +
              T * (-0x1.feba685955cbdp-23 +
                   T * (-0x1.624a929eb324bp-33 + T * 0x1.014e67fc72d86p-45))));
    h[sC2H6] = 0x1.148599bc79368p+8 *
               (T * (0x1.34dc2b0ea1837p+2 +
                     T * (0x1.c585f11b35bb8p-8 +
                          T * (-0x1.97c6d164be5d7p-20 +
                               T * (0x1.71b57f19a838bp-33 +
                                    T * -0x1.03465ae807c7fp-47)))) -
                0x1.8d6e51eb851ecp+13);
    cp[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.34dc2b0ea1837p+2 +
         T * (0x1.c585f11b35bb8p-7 +
              T * (-0x1.31d51d0b8ec61p-18 +
                   T * (0x1.71b57f19a838bp-31 + T * -0x1.4417f1a209b9ep-45))));
    h[sCH2OH] = 0x1.0bea1f388b9cfp+8 *
                (T * (0x1.df9ac1b7c5d49p+1 +
                      T * (0x1.2279c187d63ffp-8 +
                           T * (-0x1.7d017f71bedccp-20 +
                                T * (0x1.154c454622351p-32 +
                                     T * -0x1.547a4e3656df4p-46)))) -
                 0x1.ca4f6fd21ff2ep+11);
    cp[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.df9ac1b7c5d49p+1 +
         T * (0x1.2279c187d63ffp-7 +
              T * (-0x1.1dc11f954f259p-18 +
                   T * (0x1.154c454622351p-30 + T * -0x1.a998e1c3ec971p-44))));
    h[sCH3O] = 0x1.0bea1f388b9cfp+8 *
               (T * (0x1.e2a9930be0dedp+1 +
                     T * (0x1.01eee717c07fdp-8 +
                          T * (-0x1.db60e105b3b4p-21 +
                               T * (0x1.b1b1dcc557332p-34 +
                                    T * -0x1.3075c5faa9763p-48)))) +
                0x1.ff547ae147ae1p+6);
    cp[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.e2a9930be0dedp+1 +
         T * (0x1.01eee717c07fdp-7 +
              T * (-0x1.6488a8c446c7p-19 +
                   T * (0x1.b1b1dcc557332p-32 + T * -0x1.7c93377953d3bp-46))));
    h[sO2] = 0x1.03d3adab9f55ap+8 *
             (T * (0x1.d94a3c64345dp+1 +
                   T * (0x1.41a93860283abp-12 +
                        T * (-0x1.6872182ced1cbp-25 +
                             T * (0x1.384f8c6a8c0c9p-38 +
                                  T * -0x1.060b482c3439ap-52)))) -
              0x1.347b851eb851fp+10);
    cp[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.d94a3c64345dp+1 +
         T * (0x1.41a93860283abp-11 +
              T * (-0x1.0e559221b1d58p-23 +
                   T * (0x1.384f8c6a8c0c9p-36 + T * -0x1.478e1a374148p-50))));
    h[sCH3OH] = 0x1.037c7db28a41cp+8 *
                (T * (0x1.01dc22ab25b31p+2 +
                      T * (0x1.3340902436267p-8 +
                           T * (-0x1.10ee9e9dedaf3p-20 +
                                T * (0x1.df411da620aedp-34 +
                                     T * -0x1.409dcaa8618b7p-48)))) -
                 0x1.98b7a3d70a3d7p+14);
    cp[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.01dc22ab25b31p+2 +
         T * (0x1.3340902436267p-7 +
              T * (-0x1.9965edece486dp-19 +
                   T * (0x1.df411da620aedp-32 + T * -0x1.90c53d5279ee4p-46))));
    h[sHO2] = 0x1.f7c8d6a0c5c07p+7 *
              (T * (0x1.0119fbbf289f6p+2 +
                    T * (0x1.2593e46a1f9dfp-10 +
                         T * (-0x1.c597158096605p-23 +
                              T * (0x1.f675fa32f618cp-36 +
                                   T * -0x1.370671f4bb474p-49)))) +
               0x1.bf6d462c343b7p+6);
    cp[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.0119fbbf289f6p+2 +
         T * (0x1.2593e46a1f9dfp-9 +
              T * (-0x1.5431502070c84p-21 +
                   T * (0x1.f675fa32f618cp-34 + T * -0x1.84c80e71ea191p-47))));
    h[sH2O2] = 0x1.e8db171241376p+7 *
               (T * (0x1.24aec4a4095f2p+2 +
                     T * (0x1.1c2c4a4f9e3cbp-9 +
                          T * (-0x1.07e7e77f421c7p-21 +
                               T * (0x1.0243c5162bd81p-34 +
                                    T * -0x1.9ca56b756b127p-49)))) -
                0x1.195bd70a3d70ap+14);
    cp[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.24aec4a4095f2p+2 +
         T * (0x1.1c2c4a4f9e3cbp-8 +
              T * (-0x1.8bdbdb3ee32aap-20 +
                   T * (0x1.0243c5162bd81p-32 + T * -0x1.01e7632962eb8p-46))));
    h[sHCCO] = 0x1.954e7e0035abdp+7 *
               (T * (0x1.683486198a14cp+2 +
                     T * (0x1.0bbca21f5e9bfp-9 +
                          T * (-0x1.1d28ea5b6f3c1p-21 +
                               T * (0x1.3abf2c56d8d5fp-34 +
                                    T * -0x1.17b243118233fp-48)))) +
                0x1.2dfcdc28f5c29p+14);
    cp[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.683486198a14cp+2 +
         T * (0x1.0bbca21f5e9bfp-8 +
              T * (-0x1.abbd5f8926da2p-20 +
                   T * (0x1.3abf2c56d8d5fp-32 + T * -0x1.5d9ed3d5e2c0fp-46))));
    h[sCH2CO] = 0x1.8b966bc75924bp+7 *
                (T * (0x1.20b91864fbad3p+2 +
                      T * (0x1.2707a64c0b6f9p-8 +
                           T * (-0x1.75123ec2bd6e9p-20 +
                                T * (0x1.fb9d615c62415p-33 +
                                     T * -0x1.1e5ee26606538p-46)))) -
                 0x1.d7f0d989df117p+12);
    cp[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.20b91864fbad3p+2 +
         T * (0x1.2707a64c0b6f9p-7 +
              T * (-0x1.17cdaf120e12fp-18 +
                   T * (0x1.fb9d615c62415p-31 + T * -0x1.65f69aff87e86p-44))));
    h[sHCCOH] = 0x1.1451d7eb0270bp+7 *
                (T * (0x1.7b200416e5f59p+2 +
                      T * (0x1.bd24e4100a643p-9 +
                           T * (-0x1.cb2d8a19f367cp-21 +
                                T * (0x1.eea583d58eb03p-34 +
                                     T * -0x1.af7b79e42769ap-48)))) +
                 0x1.c60a04189374cp+12);
    cp[sHCCOH] =
        0x1.1451d7eb0270bp+7 *
        (0x1.7b200416e5f59p+2 +
         T * (0x1.bd24e4100a643p-8 +
              T * (-0x1.58622793768ddp-19 +
                   T * (0x1.eea583d58eb03p-32 + T * -0x1.0dad2c2e98a2p-45))));
    h[sCH2HCO] = 0x1.8252e1709b12p+7 *
                 (T * (0x1.7e71602a0c43p+2 +
                       T * (0x1.0a6c58147f6e1p-8 +
                            T * (-0x1.eafda06d3ace3p-21 +
                                 T * (0x1.bf88e00b7cb09p-34 +
                                      T * -0x1.3998dc7039d47p-48)))) +
                  0x1.ea52602c9081cp+8);
    cp[sCH2HCO] =
        0x1.8252e1709b12p+7 *
        (0x1.7e71602a0c43p+2 +
         T * (0x1.0a6c58147f6e1p-7 +
              T * (-0x1.703e3851ec1aap-19 +
                   T * (0x1.bf88e00b7cb09p-32 + T * -0x1.87ff138c48498p-46))));
    h[sCH3CO] = 0x1.8252e1709b12p+7 *
                (T * (0x1.7c772997a8fe9p+2 +
                      T * (0x1.01c6d5a31b5ecp-8 +
                           T * (-0x1.02499c374c1f9p-20 +
                                T * (0x1.03dfd79a539a7p-33 +
                                     T * -0x1.9c2add000ed22p-48)))) -
                 0x1.d969d70a3d70ap+11);
    cp[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.7c772997a8fe9p+2 +
         T * (0x1.01c6d5a31b5ecp-7 +
              T * (-0x1.836e6a52f22f6p-19 +
                   T * (0x1.03dfd79a539a7p-31 + T * -0x1.019aca2009435p-45))));
    h[sCO2] = 0x1.79d818115de7cp+7 *
              (T * (0x1.1d0828c36da88p+2 +
                    T * (0x1.9b9696515d0c4p-10 +
                         T * (-0x1.c98fa589a2978p-22 +
                              T * (0x1.07390665390c5p-34 +
                                   T * -0x1.e110e5fdf498ap-49)))) -
               0x1.7e8deb851eb85p+15);
    cp[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.1d0828c36da88p+2 +
         T * (0x1.9b9696515d0c4p-9 +
              T * (-0x1.572bbc2739f1ap-20 +
                   T * (0x1.07390665390c5p-32 + T * -0x1.2caa8fbeb8df6p-46))));
    h[sCH3HCO] = 0x1.797bdf23f83eep+7 *
                 (T * (0x1.59dcf38b7d772p+2 +
                       T * (0x1.80242581cd52bp-8 +
                            T * (-0x1.7a2a05a17d829p-20 +
                                 T * (0x1.77e1ab96723c6p-33 +
                                      T * -0x1.2753ba556825fp-47)))) -
                  0x1.61047ced91687p+14);
    cp[sCH3HCO] =
        0x1.797bdf23f83eep+7 *
        (0x1.59dcf38b7d772p+2 +
         T * (0x1.80242581cd52bp-7 +
              T * (-0x1.1b9f84391e21fp-18 +
                   T * (0x1.77e1ab96723c6p-31 + T * -0x1.7128a8eac22f7p-45))));
    h[sOCHO] = 0x1.716240449121bp+7 *
               (T * (0x1.88151982321eep+2 +
                     T * (0x1.ec4f6e2f37151p-10 +
                          T * (-0x1.fc465c6f552c4p-22 +
                               T * (0x1.03f4e74ff3af1p-34 +
                                    T * -0x1.9f88f351e235cp-49)))) -
                0x1.542762eb1c433p+14);
    cp[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.88151982321eep+2 +
         T * (0x1.ec4f6e2f37151p-9 +
              T * (-0x1.7d34c5537fe13p-20 +
                   T * (0x1.03f4e74ff3af1p-32 + T * -0x1.03b598132d619p-46))));
    h[sCH3CHOH] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.b10f07dabfd49p+2 +
                        T * (0x1.7d3cb9104ecf1p-8 +
                             T * (-0x1.5225201ab9133p-20 +
                                  T * (0x1.27ed681c80e62p-33 +
                                       T * -0x1.89a7d5f06faf9p-48)))) -
                   0x1.5e94c01a36e2fp+12);
    cp[sCH3CHOH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.b10f07dabfd49p+2 +
         T * (0x1.7d3cb9104ecf1p-7 +
              T * (-0x1.fb37b028159ccp-19 +
                   T * (0x1.27ed681c80e62p-31 + T * -0x1.ec11cb6c8b9b7p-46))));
    h[sC2H4OH] = 0x1.710a1c4780b3cp+7 *
                 (T * (0x1.e60aac2b6c955p+2 +
                       T * (0x1.317e701a04be2p-8 +
                            T * (-0x1.0f277fea138a6p-20 +
                                 T * (0x1.db2b1fb40206bp-34 +
                                      T * -0x1.3c9f232902112p-48)))) -
                  0x1.68cc902de00d2p+12);
    cp[sC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.e60aac2b6c955p+2 +
         T * (0x1.317e701a04be2p-7 +
              T * (-0x1.96bb3fdf1d4f9p-19 +
                   T * (0x1.db2b1fb40206bp-32 + T * -0x1.8bc6ebf342956p-46))));
    h[sCH3CH2O] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.80bb583042bdfp+2 +
                        T * (0x1.8ea140e07883dp-8 +
                             T * (-0x1.69efdaa20e5fp-20 +
                                  T * (0x1.44c6e7f9eb7edp-33 +
                                       T * -0x1.be51a37bb39dcp-48)))) -
                   0x1.348b2fec56d5dp+12);
    cp[sCH3CH2O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.80bb583042bdfp+2 +
         T * (0x1.8ea140e07883dp-7 +
              T * (-0x1.0f73e3f98ac74p-18 +
                   T * (0x1.44c6e7f9eb7edp-31 + T * -0x1.16f3062d50429p-45))));
    h[sCH3OCH2] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.057bee98a47c2p+3 +
                        T * (0x1.68bafb5b5900dp-8 +
                             T * (-0x1.561f7ea367a58p-20 +
                                  T * (0x1.49a76f092f448p-33 +
                                       T * -0x1.f8dc59a293e24p-48)))) -
                   0x1.ab6d504816fp+11);
    cp[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.057bee98a47c2p+3 +
         T * (0x1.68bafb5b5900dp-7 +
              T * (-0x1.00979efa8dbc2p-18 +
                   T * (0x1.49a76f092f448p-31 + T * -0x1.3b89b8059c6d6p-45))));
    h[sHCOOH] = 0x1.694b465e7f561p+7 *
                (T * (0x1.abfd378379f29p+2 +
                      T * (0x1.510b702993614p-9 +
                           T * (-0x1.4620e6d5987f1p-21 +
                                T * (0x1.3e8cb1d51c3f6p-34 +
                                     T * -0x1.ec902fb5ac402p-49)))) -
                 0x1.7a1f147ae147bp+15);
    cp[sHCOOH] =
        0x1.694b465e7f561p+7 *
        (0x1.abfd378379f29p+2 +
         T * (0x1.510b702993614p-8 +
              T * (-0x1.e9315a4064be9p-20 +
                   T * (0x1.3e8cb1d51c3f6p-32 + T * -0x1.33da1dd18ba81p-46))));
    h[sCH3OCH3] = 0x1.68f6f3733777dp+7 *
                  (T * (0x1.a960a7be2821ep-1 +
                        T * (0x1.b9037309df683p-7 +
                             T * (-0x1.36a8599d91d61p-18 +
                                  T * (0x1.dd9eff9b15bfp-31 +
                                       T * -0x1.33c8384bf040cp-44)))) -
                   0x1.6dd063d70a3d7p+14);
    cp[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.a960a7be2821ep-1 +
         T * (0x1.b9037309df683p-6 +
              T * (-0x1.d1fc866c5ac12p-17 +
                   T * (0x1.dd9eff9b15bfp-29 + T * -0x1.80ba465eec50fp-42))));
    h[sC2H5OH] = 0x1.68f6f3733777dp+7 *
                 (T * (0x1.a3fef5a964e8bp+2 +
                       T * (0x1.f236422024d26p-8 +
                            T * (-0x1.e242a6d3f2e49p-20 +
                                 T * (0x1.da06ea3159551p-33 +
                                      T * -0x1.7194f51aaef84p-47)))) -
                  0x1.ec967be76c8b4p+14);
    cp[sC2H5OH] =
        0x1.68f6f3733777dp+7 *
        (0x1.a3fef5a964e8bp+2 +
         T * (0x1.f236422024d26p-7 +
              T * (-0x1.69b1fd1ef62b7p-18 +
                   T * (0x1.da06ea3159551p-31 + T * -0x1.cdfa32615ab65p-45))));
    h[sHOCH2O] = 0x1.618d10efa02bbp+7 *
                 (T * (0x1.994b347c088f2p+2 +
                       T * (0x1.e75fa1fc6bad4p-9 +
                            T * (-0x1.c025f44665ca3p-21 +
                                 T * (0x1.a72e020a6d527p-34 +
                                      T * -0x1.3f9de4f37299cp-48)))) -
                  0x1.82b8276c8b439p+14);
    cp[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.994b347c088f2p+2 +
         T * (0x1.e75fa1fc6bad4p-8 +
              T * (-0x1.501c7734cc57ap-19 +
                   T * (0x1.a72e020a6d527p-32 + T * -0x1.8f855e304f403p-46))));
    h[sCH3OCO] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.a2ceee0f3cb3ep+3 +
                       T * (0x1.293c374670db4p-9 +
                            T * (-0x1.277392569aaa2p-21 +
                                 T * (0x1.25c958ea02d62p-34 +
                                      T * -0x1.cbf2f12bc093dp-49)))) -
                  0x1.81568f5c28f5cp+14);
    cp[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.a2ceee0f3cb3ep+3 +
         T * (0x1.293c374670db4p-8 +
              T * (-0x1.bb2d5b81e7ff3p-20 +
                   T * (0x1.25c958ea02d62p-32 + T * -0x1.1f77d6bb585c6p-46))));
    h[sCH3OCHO] = 0x1.14e89ec288e65p+7 *
                  (T * (0x1.161e993d7e3acp+3 +
                        T * (0x1.7a7b0a942be72p-8 +
                             T * (-0x1.7ec5f6447beccp-20 +
                                  T * (0x1.8238bcfb58bddp-33 +
                                       T * -0x1.31c3b96f6a766p-47)))) -
                   0x1.6ac8f42c3c9efp+15);
    cp[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.161e993d7e3acp+3 +
         T * (0x1.7a7b0a942be72p-7 +
              T * (-0x1.1f1478b35cf19p-18 +
                   T * (0x1.8238bcfb58bddp-31 + T * -0x1.7e34a7cb4513fp-45))));
    h[sCH3OCH2O] = 0x1.10565da70e7acp+7 *
                   (T * (0x1.1348a67cd6eb3p+3 +
                         T * (0x1.bce5f8dc8efb7p-8 +
                              T * (-0x1.b1aaf8e73a043p-20 +
                                   T * (0x1.ab94dc4bba297p-33 +
                                        T * -0x1.4d5cd732c13eap-47)))) -
                    0x1.4e00fa43fe5c9p+14);
    cp[sCH3OCH2O] =
        0x1.10565da70e7acp+7 *
        (0x1.1348a67cd6eb3p+3 +
         T * (0x1.bce5f8dc8efb7p-7 +
              T * (-0x1.45403aad6b832p-18 +
                   T * (0x1.ab94dc4bba297p-31 + T * -0x1.a0b40cff718e4p-45))));
    h[sCH3OCH2OH] = 0x1.0bea1f388b9cfp+7 *
                    (T * (0x1.16b6cf6a35fecp+3 +
                          T * (0x1.f753023766137p-8 +
                               T * (-0x1.e414fda179e1p-20 +
                                    T * (0x1.d91af2017f03fp-33 +
                                         T * -0x1.6ea4ad6d3a24cp-47)))) -
                     0x1.74596c49ba5e3p+15);
    cp[sCH3OCH2OH] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.16b6cf6a35fecp+3 +
         T * (0x1.f753023766137p-7 +
              T * (-0x1.6b0fbe391b68cp-18 +
                   T * (0x1.d91af2017f03fp-31 + T * -0x1.ca4dd8c888adfp-45))));
    h[sOCH2OCHO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.80bf9fbda0092p+3 +
                         T * (0x1.09d5a4f1c13fdp-8 +
                              T * (-0x1.04b38b9da68f1p-20 +
                                   T * (0x1.00ec4fd57a83cp-33 +
                                        T * -0x1.8fbd81a60883dp-48)))) -
                    0x1.52c9723a29c78p+15);
    cp[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.80bf9fbda0092p+3 +
         T * (0x1.09d5a4f1c13fdp-7 +
              T * (-0x1.870d516c79d69p-19 +
                   T * (0x1.00ec4fd57a83cp-31 + T * -0x1.f3ace20f8aa4cp-46))));
    h[sHOCH2OCO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.6bf5abb377913p+3 +
                         T * (0x1.0bee9e8151789p-8 +
                              T * (-0x1.054ec018845bdp-20 +
                                   T * (0x1.0091919e564d2p-33 +
                                        T * -0x1.8e28282db74cep-48)))) -
                    0x1.6bbb260aa64c3p+15);
    cp[sHOCH2OCO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.6bf5abb377913p+3 +
         T * (0x1.0bee9e8151789p-7 +
              T * (-0x1.87f62024c689bp-19 +
                   T * (0x1.0091919e564d2p-31 + T * -0x1.f1b2323925201p-46))));
    h[sHOC2H4O2] = 0x1.af956cbfc9bdbp+6 *
                   (T * (0x1.4303562b85bfdp+3 +
                         T * (0x1.95ed3fe77f334p-8 +
                              T * (-0x1.4e7b24e8e297bp-20 +
                                   T * (0x1.2ca5c0315c64bp-33 +
                                        T * -0x1.bead2283f207cp-48)))) -
                    0x1.736c3573eab36p+14);
    cp[sHOC2H4O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.4303562b85bfdp+3 +
         T * (0x1.95ed3fe77f334p-7 +
              T * (-0x1.f5b8b75d53e38p-19 +
                   T * (0x1.2ca5c0315c64bp-31 + T * -0x1.172c35927744dp-45))));
    h[sCH3OCH2O2] = 0x1.af956cbfc9bdbp+6 *
                    (T * (0x1.8d9960c465f6p+3 +
                          T * (0x1.84f9cc62ae45ep-8 +
                               T * (-0x1.6cfd13ffb9563p-20 +
                                    T * (0x1.5d440c2e71331p-33 +
                                         T * -0x1.0a3369c503bc5p-47)))) -
                     0x1.66dfb1f8a0903p+14);
    cp[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.8d9960c465f6p+3 +
         T * (0x1.84f9cc62ae45ep-7 +
              T * (-0x1.11bdceffcb00ap-18 +
                   T * (0x1.5d440c2e71331p-31 + T * -0x1.4cc0443644ab6p-45))));
    h[sCH2OCH2O2H] = 0x1.af956cbfc9bdbp+6 *
                     (T * (0x1.e3d04f029c927p+3 +
                           T * (0x1.2eaf27f746174p-8 +
                                T * (-0x1.1d8ceaccac18dp-20 +
                                     T * (0x1.1264256c3fbcep-33 +
                                          T * -0x1.a39bff7cedbedp-48)))) -
                      0x1.1fadf2617c1bep+14);
    cp[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.e3d04f029c927p+3 +
         T * (0x1.2eaf27f746174p-7 +
              T * (-0x1.ac53603302254p-19 +
                   T * (0x1.1264256c3fbcep-31 + T * -0x1.06417fae14974p-45))));
    h[sCH3OCH2O2H] = 0x1.aa02db5960032p+6 *
                     (T * (0x1.ddfcb196e660fp+3 +
                           T * (0x1.8777336d5eb5fp-8 +
                                T * (-0x1.7151b674a2ed1p-20 +
                                     T * (0x1.62d3277c8f756p-33 +
                                          T * -0x1.0f3eb2d2d06fep-47)))) -
                      0x1.411836ae7d567p+15);
    cp[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.ddfcb196e660fp+3 +
         T * (0x1.8777336d5eb5fp-7 +
              T * (-0x1.14fd48d77a31dp-18 +
                   T * (0x1.62d3277c8f756p-31 + T * -0x1.530e5f87848bdp-45))));
    h[sHO2CH2OCHO] = 0x1.694b465e7f561p+6 *
                     (T * (0x1.0755ba7c68307p+4 +
                           T * (0x1.176846f7c96d5p-8 +
                                T * (-0x1.101dbc3c42d91p-20 +
                                     T * (0x1.0af5b06e6e2adp-33 +
                                          T * -0x1.9e1104836475p-48)))) -
                      0x1.e777ebedfa44p+15);
    cp[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.0755ba7c68307p+4 +
         T * (0x1.176846f7c96d5p-7 +
              T * (-0x1.982c9a5a64459p-19 +
                   T * (0x1.0af5b06e6e2adp-31 + T * -0x1.02caa2d21ec92p-45))));
    h[sO2CH2OCH2O2H] = 0x1.30f32e834a014p+6 *
                       (T * (0x1.3342c89cbc63dp+4 +
                             T * (0x1.5614bd65f761cp-8 +
                                  T * (-0x1.42a4e72863769p-20 +
                                       T * (0x1.35f2c9e6b47e8p-33 +
                                            T * -0x1.d9dc73559a61dp-48)))) -
                        0x1.284169374bc6ap+15);
    cp[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (0x1.3342c89cbc63dp+4 +
         T * (0x1.5614bd65f761cp-7 +
              T * (-0x1.e3f75abc9531dp-19 +
                   T * (0x1.35f2c9e6b47e8p-31 + T * -0x1.2829c815807d2p-45))));
    h[sAR] = 0x1.a0439e3ea879p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a0439e3ea879p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    h[sHE] = 0x1.9c055f0ad864p+8 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sHE] =
        0x1.9c055f0ad864p+8 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
  } else if (T >= 299.999999) {
    h[sH] = 0x1.01c3c6b46b46cp+13 *
            (T * (0x1.4p+1 +
                  T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) +
             0x1.8dfe851eb851fp+14);
    cp[sH] =
        0x1.01c3c6b46b46cp+13 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    h[sH2] = 0x1.01c3c6b46b46cp+12 *
             (T * (0x1.a628ed5f138bdp+1 +
                   T * (0x1.b08222e154b9ep-12 +
                        T * (-0x1.237329fa9f7bbp-22 +
                             T * (-0x1.a0bbe70f69bcfp-36 +
                                  T * 0x1.746fa82e00818p-44)))) -
              0x1.fa42b020c49bap+9);
    cp[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.a628ed5f138bdp+1 +
         T * (0x1.b08222e154b9ep-11 +
              T * (-0x1.b52cbef7ef399p-21 +
                   T * (-0x1.a0bbe70f69bcfp-34 + T * 0x1.d18b923980a1ep-42))));
    h[sCH2] = 0x1.286500f4d652p+9 *
              (T * (0x1.e19f746480de1p+1 +
                    T * (0x1.fbf7d158743a2p-12 +
                         T * (0x1.f42aa333e052bp-21 +
                              T * (-0x1.08a1f3bb46844p-30 +
                                   T * 0x1.7bf8fa708a0d4p-42)))) +
               0x1.67681487fcb92p+15);
    cp[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.e19f746480de1p+1 +
         T * (0x1.fbf7d158743a2p-11 +
              T * (0x1.771ffa66e83ep-19 +
                   T * (-0x1.08a1f3bb46844p-28 + T * 0x1.daf7390cac908p-40))));
    h[sCH2S] = 0x1.286500f4d652p+9 *
               (T * (0x1.0cb5ee035346ap+2 +
                     T * (-0x1.36326518bab75p-10 +
                          T * (0x1.7056247403375p-19 +
                               T * (-0x1.cb9b5a07323d1p-30 +
                                    T * 0x1.b58ed1c91768cp-42)))) +
                0x1.8a81a1f212d77p+15);
    cp[sCH2S] =
        0x1.286500f4d652p+9 *
        (0x1.0cb5ee035346ap+2 +
         T * (-0x1.36326518bab75p-9 +
              T * (0x1.14409b5702698p-17 +
                   T * (-0x1.cb9b5a07323d1p-28 + T * 0x1.1179431daea17p-39))));
    h[sCH3] = 0x1.148599bc79368p+9 *
              (T * (0x1.d41e76e38c2bfp+1 +
                    T * (0x1.16bcc8dd4edc2p-10 +
                         T * (0x1.e868877df9784p-20 +
                              T * (-0x1.c6cada0b83d8bp-30 +
                                   T * 0x1.159d2d58edf15p-41)))) +
               0x1.009add2f1a9fcp+14);
    cp[sCH3] =
        0x1.148599bc79368p+9 *
        (0x1.d41e76e38c2bfp+1 +
         T * (0x1.16bcc8dd4edc2p-9 +
              T * (0x1.6e4e659e7b1a3p-18 +
                   T * (-0x1.c6cada0b83d8bp-28 + T * 0x1.5b0478af296dap-39))));
    h[sO] = 0x1.03d3adab9f55ap+9 *
            (T * (0x1.792495e17e34cp+1 +
                  T * (-0x1.ad6f7594e8c1p-11 +
                       T * (0x1.b142b3935e374p-21 +
                            T * (-0x1.b8960c8cd58ecp-32 +
                                 T * 0x1.5e71577b56d5cp-44)))) +
             0x1.c76e8f5c28f5cp+14);
    cp[sO] =
        0x1.03d3adab9f55ap+9 *
        (0x1.792495e17e34cp+1 +
         T * (-0x1.ad6f7594e8c1p-10 +
              T * (0x1.44f206ae86a97p-19 +
                   T * (-0x1.b8960c8cd58ecp-30 + T * 0x1.b60dad5a2c8b3p-42))));
    h[sCH4] = 0x1.0325882935d13p+9 *
              (T * (0x1.8eb734b51372ap-1 +
                    T * (0x1.1e568242bae5ap-7 +
                         T * (-0x1.3751abbc576eep-17 +
                              T * (0x1.05f7c9c3f6e46p-27 +
                                   T * -0x1.588185ebbec0cp-39)))) -
               0x1.3309d4fdf3b64p+13);
    cp[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.8eb734b51372ap-1 +
         T * (0x1.1e568242bae5ap-6 +
              T * (-0x1.d2fa819a83265p-16 +
                   T * (0x1.05f7c9c3f6e46p-25 + T * -0x1.aea1e766ae70ep-37))));
    h[sOH] = 0x1.e8db171241376p+8 *
             (T * (0x1.080501d23d242p+2 +
                   T * (-0x1.a6c41f4a374cap-10 +
                        T * (0x1.240abf253f3d4p-19 +
                             T * (-0x1.8e78ee8164a9cp-30 +
                                  T * 0x1.d067c130d5fbcp-42)))) +
              0x1.a249e4649906dp+11);
    cp[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.080501d23d242p+2 +
         T * (-0x1.a6c41f4a374cap-9 +
              T * (0x1.b6101eb7dedbep-18 +
                   T * (-0x1.8e78ee8164a9cp-28 + T * 0x1.2240d8be85bd5p-39))));
    h[sH2O] = 0x1.cd8113bbf9d1p+8 *
              (T * (0x1.b18409e55c0fdp+1 +
                    T * (0x1.c7790c169fe5dp-10 +
                         T * (-0x1.1c4de5b6f12f7p-19 +
                              T * (0x1.dee092cb7ea3cp-30 +
                                   T * -0x1.1a377aef067a5p-41)))) -
               0x1.d80070a3d70a4p+14);
    cp[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.b18409e55c0fdp+1 +
         T * (0x1.c7790c169fe5dp-9 +
              T * (-0x1.aa74d89269c73p-18 +
                   T * (0x1.dee092cb7ea3cp-28 + T * -0x1.60c559aac818ep-39))));
    h[sC2H] = 0x1.4c34d172a6531p+8 *
              (T * (0x1.71e04a987f933p+1 +
                    T * (0x1.b76ae82ebca6cp-8 +
                         T * (-0x1.3e82612c20167p-17 +
                              T * (0x1.fa727902a68efp-28 +
                                   T * -0x1.33bda806adcfp-39)))) +
               0x1.051764a8c154dp+16);
    cp[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.71e04a987f933p+1 +
         T * (0x1.b76ae82ebca6cp-7 +
              T * (-0x1.ddc391c23021ap-16 +
                   T * (0x1.fa727902a68efp-26 + T * -0x1.80ad12085942bp-37))));
    h[sC2H2] = 0x1.3f584141b0716p+8 *
               (T * (0x1.9e0b72c73f3bap-1 +
                     T * (0x1.7ec17f28e481ap-7 +
                          T * (-0x1.8d40c15d06efbp-17 +
                               T * (0x1.e14c5845a203fp-28 +
                                    T * -0x1.de8c6d30d264p-40)))) +
                0x1.9cf3ec3c9eeccp+14);
    cp[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.9e0b72c73f3bap-1 +
         T * (0x1.7ec17f28e481ap-6 +
              T * (-0x1.29f09105c533cp-15 +
                   T * (0x1.e14c5845a203fp-26 + T * -0x1.2b17c43e837e8p-37))));
    h[sC2H3] = 0x1.337122eb1b451p+8 *
               (T * (0x1.9b3219c31fa4ep+1 +
                     T * (0x1.8d17f1df63fcdp-11 +
                          T * (0x1.21ebbad5b9bbdp-17 +
                               T * (-0x1.3339cad4b897fp-27 +
                                    T * 0x1.9e3160f1cd55p-39)))) +
                0x1.1057b18fc5048p+15);
    cp[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.9b3219c31fa4ep+1 +
         T * (0x1.8d17f1df63fcdp-10 +
              T * (0x1.b2e198409699cp-16 +
                   T * (-0x1.3339cad4b897fp-25 + T * 0x1.02dedc9720552p-36))));
    h[sCO] = 0x1.28d6c752d48edp+8 *
             (T * (0x1.a19806f262889p+1 +
                   T * (0x1.8c58a4980b8b3p-11 +
                        T * (-0x1.5b55640bd90cfp-20 +
                             T * (0x1.7f9698eb1c752p-30 +
                                  T * -0x1.16a79b13a86b2p-41)))) -
              0x1.bf3451eb851ecp+13);
    cp[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.a19806f262889p+1 +
         T * (0x1.8c58a4980b8b3p-10 +
              T * (-0x1.04800b08e2c9bp-18 +
                   T * (0x1.7f9698eb1c752p-28 + T * -0x1.5c5181d89285ep-39))));
    h[sN2] = 0x1.28bba88de9adap+8 *
             (T * (0x1.a63b0c4588a05p+1 +
                   T * (0x1.712962facc0e9p-11 +
                        T * (-0x1.629f839620fc5p-20 +
                             T * (0x1.83ae94da0fe59p-30 +
                                  T * -0x1.134425cafd91ep-41)))) -
              0x1.fe73333333333p+9);
    cp[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.a63b0c4588a05p+1 +
         T * (0x1.712962facc0e9p-10 +
              T * (-0x1.09f7a2b098bd4p-18 +
                   T * (0x1.83ae94da0fe59p-28 + T * -0x1.58152f3dbcf65p-39))));
    h[sC2H4] = 0x1.286500f4d652p+8 *
               (T * (0x1.fac71d356ff96p+1 +
                     T * (-0x1.f0244a6c1abf8p-9 +
                          T * (0x1.3f5227836c76ep-16 +
                               T * (-0x1.2908fcd07cb16p-26 +
                                    T * 0x1.7bd417ca8706dp-38)))) +
                0x1.3e1c6a35935fcp+12);
    cp[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.fac71d356ff96p+1 +
         T * (-0x1.f0244a6c1abf8p-8 +
              T * (0x1.defb3b4522b25p-15 +
                   T * (-0x1.2908fcd07cb16p-24 + T * 0x1.dac91dbd28c88p-36))));
    h[sHCO] = 0x1.1e87151786dcap+8 *
              (T * (0x1.72fc7a398201dp+1 +
                    T * (0x1.96446da0caeeep-9 +
                         T * (-0x1.ae878cacc2589p-19 +
                              T * (0x1.76760551ca70ep-29 +
                                   T * -0x1.018b0a8d6588ep-40)))) +
               0x1.03fec083126e9p+12);
    cp[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.72fc7a398201dp+1 +
         T * (0x1.96446da0caeeep-8 +
              T * (-0x1.42e5a98191c27p-17 +
                   T * (0x1.76760551ca70ep-27 + T * -0x1.41edcd30beeb1p-38))));
    h[sC2H5] = 0x1.1e1d11b23ceecp+8 *
               (T * (0x1.13932d6ece13fp+2 +
                     T * (-0x1.122932b551339p-9 +
                          T * (0x1.15fbb31ddde43p-16 +
                               T * (-0x1.014b347d2de42p-26 +
                                    T * 0x1.4460e4bf6afcdp-38)))) +
                0x1.914db645a1cacp+13);
    cp[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.13932d6ece13fp+2 +
         T * (-0x1.122932b551339p-8 +
              T * (0x1.a0f98cacccd65p-15 +
                   T * (-0x1.014b347d2de42p-24 + T * 0x1.95791def45bcp-36))));
    h[sCH2O] = 0x1.14e89ec288e65p+8 *
               (T * (0x1.591f1645bca2bp+1 +
                     T * (0x1.42d6f412ede0ap-9 +
                          T * (0x1.287293e432887p-22 +
                               T * (-0x1.2e935feb740dp-33 +
                                    T * -0x1.64c7317c1815cp-44)))) -
                0x1.d3d65810624ddp+13);
    cp[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.591f1645bca2bp+1 +
         T * (0x1.42d6f412ede0ap-8 +
              T * (0x1.bcabddd64bccbp-21 +
                   T * (-0x1.2e935feb740dp-31 + T * -0x1.bdf8fddb1e1b3p-42))));
    h[sC2H6] = 0x1.148599bc79368p+8 *
               (T * (0x1.7668f4b61fe22p+0 +
                     T * (0x1.fbbab674c6a67p-8 +
                          T * (0x1.029d91c333bc3p-19 +
                               T * (-0x1.b03012a1e9c75p-29 +
                                    T * 0x1.022f12b8754e3p-40)))) -
                0x1.5f3970a3d70a4p+13);
    cp[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.7668f4b61fe22p+0 +
         T * (0x1.fbbab674c6a67p-7 +
              T * (0x1.83ec5aa4cd9a5p-18 +
                   T * (-0x1.b03012a1e9c75p-27 + T * 0x1.42bad76692a1bp-38))));
    h[sCH2OH] = 0x1.0bea1f388b9cfp+8 *
                (T * (0x1.272aaace75438p+2 +
                      T * (-0x1.98fe71be10fadp-10 +
                           T * (0x1.8d6a44e053b58p-17 +
                                T * (-0x1.a82a71cdd7b0dp-27 +
                                     T * 0x1.360185672271cp-38)))) -
                 0x1.c282594af4f0ep+11);
    cp[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.272aaace75438p+2 +
         T * (-0x1.98fe71be10fadp-9 +
              T * (0x1.2a0fb3a83ec82p-15 +
                   T * (-0x1.a82a71cdd7b0dp-25 + T * 0x1.8381e6c0eb0e2p-36))));
    h[sCH3O] = 0x1.0bea1f388b9cfp+8 *
               (T * (0x1.0d9817b95a294p+1 +
                     T * (0x1.d8f25f83733c9p-9 +
                          T * (0x1.ddadaadf4f2cbp-20 +
                               T * (-0x1.fafcbebd8fc7ap-30 +
                                    T * 0x1.d362d3ee4ce98p-42)))) +
                0x1.e94cf0d844d01p+9);
    cp[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.0d9817b95a294p+1 +
         T * (0x1.d8f25f83733c9p-8 +
              T * (0x1.664240277b618p-18 +
                   T * (-0x1.fafcbebd8fc7ap-28 + T * 0x1.241dc474f011fp-39))));
    h[sO2] = 0x1.03d3adab9f55ap+8 *
             (T * (0x1.9b417ca2120e2p+1 +
                   T * (0x1.27904dfc9e5cep-11 +
                        T * (-0x1.9c0a9074060ddp-23 +
                             T * (0x1.6927dfded5549p-32 +
                                  T * -0x1.8ae68b720d9f3p-43)))) -
              0x1.f69fdf3b645a2p+9);
    cp[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.9b417ca2120e2p+1 +
         T * (0x1.27904dfc9e5cep-10 +
              T * (-0x1.3507ec57048a6p-21 +
                   T * (0x1.6927dfded5549p-30 + T * -0x1.eda02e4e9106fp-41))));
    h[sCH3OH] = 0x1.037c7db28a41cp+8 *
                (T * (0x1.547ea5f84cad5p+1 +
                      T * (0x1.e12210c369184p-9 +
                           T * (0x1.40c85bdb4df8bp-19 +
                                T * (-0x1.2e21c09c63d3cp-29 +
                                     T * 0x1.0d277d312e753p-41)))) -
                 0x1.8c25eb851eb85p+14);
    cp[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.547ea5f84cad5p+1 +
         T * (0x1.e12210c369184p-8 +
              T * (0x1.e12c89c8f4f5p-18 +
                   T * (-0x1.2e21c09c63d3cp-27 + T * 0x1.50715c7d7a127p-39))));
    h[sHO2] = 0x1.f7c8d6a0c5c07p+7 *
              (T * (0x1.1350a899bcaa1p+2 +
                    T * (-0x1.373d054674594p-9 +
                         T * (0x1.d94d8bda368ddp-18 +
                              T * (-0x1.a110b0905e0fcp-28 +
                                   T * 0x1.058dba0c9e34ap-39)))) +
               0x1.26cedbb59ddc2p+8);
    cp[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.1350a899bcaa1p+2 +
         T * (-0x1.373d054674594p-8 +
              T * (0x1.62fa28e3a8ea6p-16 +
                   T * (-0x1.a110b0905e0fcp-26 + T * 0x1.46f1288fc5c1cp-37))));
    h[sH2O2] = 0x1.e8db171241376p+7 *
               (T * (0x1.b1c2b0ea18373p+1 +
                     T * (0x1.ae8552d47d093p-9 +
                          T * (-0x1.a9349aa3ab715p-25 +
                               T * (-0x1.3de20a2d8a14dp-30 +
                                    T * 0x1.164491df29f0bp-41)))) -
                0x1.13fc99999999ap+14);
    cp[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.b1c2b0ea18373p+1 +
         T * (0x1.ae8552d47d093p-8 +
              T * (-0x1.3ee773fac095p-23 +
                   T * (-0x1.3de20a2d8a14dp-28 + T * 0x1.5bd5b656f46cdp-39))));
    h[sHCCO] = 0x1.954e7e0035abdp+7 *
               (T * (0x1.203868265a06ep+1 +
                     T * (0x1.2142867388492p-7 +
                          T * (-0x1.0967cefa39fe4p-17 +
                               T * (0x1.28cb9772e0568p-28 +
                                    T * -0x1.1d37b00a8be7ep-40)))) +
                0x1.396dcbc6a7efap+14);
    cp[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.203868265a06ep+1 +
         T * (0x1.2142867388492p-6 +
              T * (-0x1.8e1bb67756fd6p-16 +
                   T * (0x1.28cb9772e0568p-26 + T * -0x1.64859c0d2ee1dp-38))));
    h[sCH2CO] = 0x1.8b966bc75924bp+7 *
                (T * (0x1.116315790e08dp+1 +
                      T * (0x1.28dc0ec708b6bp-7 +
                           T * (-0x1.851d295f1200dp-18 +
                                T * (0x1.410e7ab1ff68dp-29 +
                                     T * -0x1.c5a468865bdep-42)))) -
                 0x1.b82eb04ab606bp+12);
    cp[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.116315790e08dp+1 +
         T * (0x1.28dc0ec708b6bp-6 +
              T * (-0x1.23d5df074d80ap-16 +
                   T * (0x1.410e7ab1ff68dp-27 + T * -0x1.1b86c153f96acp-39))));
    h[sHCCOH] = 0x1.1451d7eb0270bp+7 *
                (T * (0x1.3e0c2d34ec70dp+0 +
                      T * (0x1.fd1641c705f4ap-7 +
                           T * (-0x1.1c77d6cfea29fp-16 +
                                T * (0x1.728b8de31e2a5p-27 +
                                     T * -0x1.8a79cae1a9a06p-39)))) +
                 0x1.f5f9d42c3c9efp+12);
    cp[sHCCOH] =
        0x1.1451d7eb0270bp+7 *
        (0x1.3e0c2d34ec70dp+0 +
         T * (0x1.fd1641c705f4ap-6 +
              T * (-0x1.aab3c237df3eep-15 +
                   T * (0x1.728b8de31e2a5p-25 + T * -0x1.ed183d9a14087p-37))));
    h[sCH2HCO] = 0x1.8252e1709b12p+7 *
                 (T * (0x1.b45c281f02fa8p+1 +
                       T * (0x1.5fe1b0115dd4p-8 +
                            T * (0x1.527eeaa41db6p-21 +
                                 T * (-0x1.ebef202e81515p-30 +
                                      T * 0x1.42d6bfa3e34f7p-41)))) +
                  0x1.7c5e809d49518p+10);
    cp[sCH2HCO] =
        0x1.8252e1709b12p+7 *
        (0x1.b45c281f02fa8p+1 +
         T * (0x1.5fe1b0115dd4p-7 +
              T * (0x1.fbbe5ff62c91p-20 +
                   T * (-0x1.ebef202e81515p-28 + T * 0x1.938c6f8cdc234p-39))));
    h[sCH3CO] = 0x1.8252e1709b12p+7 *
                (T * (0x1.0a75911134dbbp+2 +
                      T * (-0x1.e7d4d15840c19p-14 +
                           T * (0x1.7f4771b2f3311p-17 +
                                T * (-0x1.7adc6b7c3bbaap-27 +
                                     T * 0x1.e643e753fc6ecp-39)))) -
                 0x1.4c2e7e28240b8p+11);
    cp[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.0a75911134dbbp+2 +
         T * (-0x1.e7d4d15840c19p-13 +
              T * (0x1.1f7595463664dp-15 +
                   T * (-0x1.7adc6b7c3bbaap-25 + T * 0x1.2fea70947dc53p-36))));
    h[sCO2] = 0x1.79d818115de7cp+7 *
              (T * (0x1.234af4f0d844dp+1 +
                    T * (0x1.45205f5fd0b9ep-8 +
                         T * (-0x1.d1b219478b079p-19 +
                              T * (0x1.d7e00903b44b7p-30 +
                                   T * -0x1.dcc4e1d57aeep-42)))) -
               0x1.79ea47ae147aep+15);
    cp[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.234af4f0d844dp+1 +
         T * (0x1.45205f5fd0b9ep-7 +
              T * (-0x1.5d4592f5a845bp-17 +
                   T * (0x1.d7e00903b44b7p-28 + T * -0x1.29fb0d256cd4cp-39))));
    h[sCH3HCO] = 0x1.797bdf23f83eep+7 *
                 (T * (0x1.2eaf76e6106abp+2 +
                       T * (-0x1.a28ce427d2efep-10 +
                            T * (0x1.09d5a4c9ce534p-16 +
                                 T * (-0x1.ed90d262d7aa3p-27 +
                                      T * 0x1.34a728840b03p-38)))) -
                  0x1.511383126e979p+14);
    cp[sCH3HCO] =
        0x1.797bdf23f83eep+7 *
        (0x1.2eaf76e6106abp+2 +
         T * (-0x1.a28ce427d2efep-9 +
              T * (0x1.8ec0772eb57cep-15 +
                   T * (-0x1.ed90d262d7aa3p-25 + T * 0x1.81d0f2a50dc3bp-36))));
    h[sOCHO] = 0x1.716240449121bp+7 *
               (T * (0x1.5a257ce73152p+0 +
                     T * (0x1.ebc9e8f1e56ddp-8 +
                          T * (-0x1.ebaac14cf36fdp-19 +
                               T * (0x1.00ca73e1bd8ebp-30 +
                                    T * -0x1.b142644a03bep-44)))) -
                0x1.3c057573eab36p+14);
    cp[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.5a257ce73152p+0 +
         T * (0x1.ebc9e8f1e56ddp-7 +
              T * (-0x1.70c010f9b693ep-17 +
                   T * (0x1.00ca73e1bd8ebp-28 + T * -0x1.0ec97eae4256cp-41))));
    h[sCH3CHOH] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.3d9c502d3efd7p+1 +
                        T * (0x1.131683e131c4ep-7 +
                             T * (0x1.51d49ac0eb69p-20 +
                                  T * (-0x1.de685ff350807p-29 +
                                       T * 0x1.524e4ce8187c9p-40)))) -
                   0x1.f5802c3c9eeccp+11);
    cp[sCH3CHOH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.3d9c502d3efd7p+1 +
         T * (0x1.131683e131c4ep-6 +
              T * (0x1.fabee821611d8p-19 +
                   T * (-0x1.de685ff350807p-27 + T * 0x1.a6e1e0221e9bbp-38))));
    h[sC2H4OH] = 0x1.710a1c4780b3cp+7 *
                 (T * (0x1.66e63f6499dd9p+0 +
                       T * (0x1.60f6a004eda59p-7 +
                            T * (-0x1.8f8c6a303edadp-21 +
                                 T * (-0x1.f0fb7fabb4b7cp-29 +
                                      T * 0x1.c51c09b157bafp-40)))) -
                  0x1.e0ce75f6fd22p+11);
    cp[sC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.66e63f6499dd9p+0 +
         T * (0x1.60f6a004eda59p-6 +
              T * (-0x1.2ba94fa42f242p-19 +
                   T * (-0x1.f0fb7fabb4b7cp-27 + T * 0x1.1b31860ed6d4dp-37))));
    h[sCH3CH2O] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.baf1b0b1e4134p+0 +
                        T * (0x1.150757d1e7d88p-7 +
                             T * (0x1.65e15007ba44p-20 +
                                  T * (-0x1.d71cd091c485p-29 +
                                       T * 0x1.44812fd72e4e5p-40)))) -
                   0x1.9b87f212d7732p+11);
    cp[sCH3CH2O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.baf1b0b1e4134p+0 +
         T * (0x1.150757d1e7d88p-6 +
              T * (0x1.0c68fc05cbb3p-18 +
                   T * (-0x1.d71cd091c485p-27 + T * 0x1.95a17bccf9e1ep-38))));
    h[sCH3OCH2] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.74e62ad7441f1p+1 +
                        T * (0x1.4d3151fd3ae33p-7 +
                             T * (-0x1.ad5e377306007p-19 +
                                  T * (0x1.1d280428550a3p-31 +
                                       T * -0x1.34aa2dbe88ceap-45)))) -
                   0x1.291c504816fp+10);
    cp[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.74e62ad7441f1p+1 +
         T * (0x1.4d3151fd3ae33p-6 +
              T * (-0x1.4206a99644805p-17 +
                   T * (0x1.1d280428550a3p-29 + T * -0x1.81d4b92e2b024p-43))));
    h[sHCOOH] = 0x1.694b465e7f561p+7 *
                (T * (0x1.6f7bbd0fc0676p+0 +
                      T * (0x1.0ba76a4703902p-7 +
                           T * (-0x1.db633d505962fp-19 +
                                T * (0x1.c87ae70772a4cp-31 +
                                     T * -0x1.6a3f7edcdbd0fp-44)))) -
                 0x1.6afb4d013a92ap+15);
    cp[sHCOOH] =
        0x1.694b465e7f561p+7 *
        (0x1.6f7bbd0fc0676p+0 +
         T * (0x1.0ba76a4703902p-6 +
              T * (-0x1.648a6dfc430a3p-17 +
                   T * (0x1.c87ae70772a4cp-29 + T * -0x1.c4cf5e9412c52p-42))));
    h[sCH3OCH3] = 0x1.68f6f3733777dp+7 *
                  (T * (0x1.6b9515f183e65p+2 +
                        T * (-0x1.61862223826f7p-9 +
                             T * (0x1.6b362404ae8e3p-16 +
                                  T * (-0x1.59c5de29d88f4p-26 +
                                       T * 0x1.cce0f30111bacp-38)))) -
                   0x1.769e2e978d4fep+14);
    cp[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.6b9515f183e65p+2 +
         T * (-0x1.61862223826f7p-8 +
              T * (0x1.10689b0382eaap-14 +
                   T * (-0x1.59c5de29d88f4p-24 + T * 0x1.200c97e0ab14bp-35))));
    h[sC2H5OH] = 0x1.68f6f3733777dp+7 *
                 (T * (0x1.36f4decf2dd02p+2 +
                       T * (-0x1.ea3b5dff2daaep-10 +
                            T * (0x1.84fb5ba0e2354p-16 +
                                 T * (-0x1.7cc4faa3aeedep-26 +
                                      T * 0x1.eef514ee47554p-38)))) -
                  0x1.d4b0872b020c5p+14);
    cp[sC2H5OH] =
        0x1.68f6f3733777dp+7 *
        (0x1.36f4decf2dd02p+2 +
         T * (-0x1.ea3b5dff2daaep-9 +
              T * (0x1.23bc84b8a9a7fp-14 +
                   T * (-0x1.7cc4faa3aeedep-24 + T * 0x1.35592d14ec954p-35))));
    h[sHOCH2O] = 0x1.618d10efa02bbp+7 *
                 (T * (0x1.07283f191a834p+2 +
                       T * (0x1.ee0b28e595de3p-9 +
                            T * (0x1.51a2c161bff7fp-20 +
                                 T * (-0x1.72393312daf39p-30 +
                                      T * 0x1.47e5d8012317p-42)))) -
                  0x1.6e45d182a9931p+14);
    cp[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.07283f191a834p+2 +
         T * (0x1.ee0b28e595de3p-8 +
              T * (0x1.fa7422129ff3ep-19 +
                   T * (-0x1.72393312daf39p-28 + T * 0x1.99df4e016bdccp-40))));
    h[sCH3OCO] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.f8932e301419fp+1 +
                       T * (0x1.8ed7fda31b87fp-7 +
                            T * (-0x1.726e55d495c74p-18 +
                                 T * (0x1.3b1ac1221623bp-30 +
                                      T * -0x1.2adae1a1fb407p-44)))) -
                  0x1.4f01ee7d566cfp+14);
    cp[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.f8932e301419fp+1 +
         T * (0x1.8ed7fda31b87fp-6 +
              T * (-0x1.15d2c05f70557p-16 +
                   T * (0x1.3b1ac1221623bp-28 + T * -0x1.75919a0a7a109p-42))));
    h[sCH3OCHO] = 0x1.14e89ec288e65p+7 *
                  (T * (0x1.8b509ebe71954p+1 +
                        T * (0x1.4dd7288f4e0dcp-7 +
                             T * (-0x1.325d3324f5b48p-19 +
                                  T * (-0x1.905318e65fa7ep-33 +
                                       T * 0x1.fa526723596f7p-44)))) -
                   0x1.5933088ce703bp+15);
    cp[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.8b509ebe71954p+1 +
         T * (0x1.4dd7288f4e0dcp-6 +
              T * (-0x1.cb8bccb7708ecp-18 +
                   T * (-0x1.905318e65fa7ep-31 + T * 0x1.3c73807617e5ap-41))));
    h[sCH3OCH2O] = 0x1.10565da70e7acp+7 *
                   (T * (0x1.a1236b2999ac6p+1 +
                         T * (0x1.6bf6efabb70dp-7 +
                              T * (-0x1.5c51f94c84aabp-19 +
                                   T * (-0x1.0983bf6c7b32ap-34 +
                                        T * 0x1.970c635ce81a4p-44)))) -
                    0x1.2c96e28240b78p+14);
    cp[sCH3OCH2O] =
        0x1.10565da70e7acp+7 *
        (0x1.a1236b2999ac6p+1 +
         T * (0x1.6bf6efabb70dp-6 +
              T * (-0x1.053d7af9638p-17 +
                   T * (-0x1.0983bf6c7b32ap-32 + T * 0x1.fccf7c342220dp-42))));
    h[sCH3OCH2OH] = 0x1.0bea1f388b9cfp+7 *
                    (T * (0x1.944a57bd00511p+1 +
                          T * (0x1.904da5c0bd76ap-7 +
                               T * (-0x1.83e1e6b243898p-19 +
                                    T * (-0x1.04f1d14f0a606p-36 +
                                         T * 0x1.8912fd5755201p-44)))) -
                     0x1.6311c7a0f9097p+15);
    cp[sCH3OCH2OH] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.944a57bd00511p+1 +
         T * (0x1.904da5c0bd76ap-6 +
              T * (-0x1.22e96d05b2a72p-17 +
                   T * (-0x1.04f1d14f0a606p-34 + T * 0x1.eb57bcad2a681p-42))));
    h[sOCH2OCHO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.4c9a25905d5adp+2 +
                         T * (0x1.043e3563c7cd2p-7 +
                              T * (0x1.fa25fb13e4183p-24 +
                                   T * (-0x1.a380b7b5a0d5cp-30 +
                                        T * 0x1.b656e3644caddp-42)))) -
                    0x1.3a408ef34d6a1p+15);
    cp[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.4c9a25905d5adp+2 +
         T * (0x1.043e3563c7cd2p-6 +
              T * (0x1.7b9c7c4eeb122p-22 +
                   T * (-0x1.a380b7b5a0d5cp-28 + T * 0x1.11f64e1eafecap-39))));
    h[sHOCH2OCO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.853c57a9e00dcp+2 +
                         T * (0x1.a5f2ba8881b83p-8 +
                              T * (0x1.6dd2a2e236c5fp-21 +
                                   T * (-0x1.a34b96bfa03c3p-30 +
                                        T * 0x1.94eb7ea605162p-42)))) -
                    0x1.57613c91d14e4p+15);
    cp[sHOCH2OCO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.853c57a9e00dcp+2 +
         T * (0x1.a5f2ba8881b83p-7 +
              T * (0x1.125dfa29a9147p-19 +
                   T * (-0x1.a34b96bfa03c3p-28 + T * 0x1.fa265e4f865bap-40))));
    h[sHOC2H4O2] = 0x1.af956cbfc9bdbp+6 *
                   (T * (0x1.1c4b4aa163e86p+2 +
                         T * (0x1.9e51b85e9fab5p-7 +
                              T * (-0x1.5322a05bfa74fp-18 +
                                   T * (0x1.68b91d1dfb5cep-30 +
                                        T * -0x1.64726065d84d6p-43)))) -
                    0x1.53dad9b3d07c8p+14);
    cp[sHOC2H4O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.1c4b4aa163e86p+2 +
         T * (0x1.9e51b85e9fab5p-6 +
              T * (-0x1.fcb3f089f7af6p-17 +
                   T * (0x1.68b91d1dfb5cep-28 + T * -0x1.bd8ef87f4e60bp-41))));
    h[sCH3OCH2O2] = 0x1.af956cbfc9bdbp+6 *
                    (T * (0x1.1aeafbb6f016bp+1 +
                          T * (0x1.2e2f3583b7b37p-6 +
                               T * (-0x1.3c0a2d04dd0b1p-17 +
                                    T * (0x1.8da5a743c3ff8p-29 +
                                         T * -0x1.bbe5faa23e23p-42)))) -
                     0x1.3098604189375p+14);
    cp[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.1aeafbb6f016bp+1 +
         T * (0x1.2e2f3583b7b37p-5 +
              T * (-0x1.da0f43874b909p-16 +
                   T * (0x1.8da5a743c3ff8p-27 + T * -0x1.156fbca566d5ep-39))));
    h[sCH2OCH2O2H] = 0x1.af956cbfc9bdbp+6 *
                     (T * (0x1.43b4ccbb5a08bp+1 +
                           T * (0x1.5b72262f8c012p-6 +
                                T * (-0x1.a1a5e18d87dacp-17 +
                                     T * (0x1.1e48c0f339dd7p-28 +
                                          T * -0x1.4dc3f7a2a3804p-41)))) -
                      0x1.c2eaa5119ce07p+13);
    cp[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.43b4ccbb5a08bp+1 +
         T * (0x1.5b72262f8c012p-5 +
              T * (-0x1.393c692a25e41p-15 +
                   T * (0x1.1e48c0f339dd7p-26 + T * -0x1.a134f58b4c605p-39))));
    h[sCH3OCH2O2H] = 0x1.aa02db5960032p+6 *
                     (T * (0x1.32d4abe952698p+0 +
                           T * (0x1.7810051a6bf31p-6 +
                                T * (-0x1.99a57a4ca5e97p-17 +
                                     T * (0x1.00872d70f7f74p-28 +
                                          T * -0x1.15093adc83a05p-41)))) -
                      0x1.1d70a1d7dbf48p+15);
    cp[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.32d4abe952698p+0 +
         T * (0x1.7810051a6bf31p-5 +
              T * (-0x1.333c1bb97c6f1p-15 +
                   T * (0x1.00872d70f7f74p-26 + T * -0x1.5a4b8993a4886p-39))));
    h[sHO2CH2OCHO] = 0x1.694b465e7f561p+6 *
                     (T * (0x1.bd5b92377a95ap+1 +
                           T * (0x1.4a193dd18d4ap-6 +
                                T * (-0x1.7138934ae79cdp-17 +
                                     T * (0x1.cda868a48a234p-29 +
                                          T * -0x1.ec3f3b33f46c7p-42)))) -
                      0x1.c59dfc9eecbfbp+15);
    cp[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.bd5b92377a95ap+1 +
         T * (0x1.4a193dd18d4ap-5 +
              T * (-0x1.14ea6e782db5ap-15 +
                   T * (0x1.cda868a48a234p-27 + T * -0x1.33a7850078c3cp-39))));
    h[sO2CH2OCH2O2H] = 0x1.30f32e834a014p+6 *
                       (T * (0x1.ff146e7701135p+0 +
                             T * (0x1.ddc767e8598bbp-6 +
                                  T * (-0x1.3567c14d95457p-16 +
                                       T * (0x1.be59e2473ce95p-28 +
                                            T * -0x1.0c9b44033abefp-40)))) -
                        0x1.ffeb7f2e48e8ap+14);
    cp[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (0x1.ff146e7701135p+0 +
         T * (0x1.ddc767e8598bbp-5 +
              T * (-0x1.d01ba1f45fe82p-15 +
                   T * (0x1.be59e2473ce95p-26 + T * -0x1.4fc21504096eap-38))));
    h[sAR] = 0x1.a0439e3ea879p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a0439e3ea879p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    h[sHE] = 0x1.9c055f0ad864p+8 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sHE] =
        0x1.9c055f0ad864p+8 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
  } else {
    ComputeThermoData(h, cp, 300.0);
    for (i = 0; i < sEnd; i++) {
      h[i] = (T - 300.) * cp[i] + h[i];
    }
  }
}

double MAX_C(double X1, double X2) { return ((X1 > X2) ? X1 : X2); }
// double MIN_C(double X1, double X2) { return ((X2 > X1) ? X1 : X2); }

// double GetPlogRateCoeff(double /* temp */, double pressure, double lgt, double rt_inv,
//                         double* PlogP, double* PlogA, double* PlogB,
//                         double* PlogE, int np) {
//   double kR, kR_l, kR_r;
//   int i;
//   if (pressure <= PlogP[0])
//     kR = PlogA[0] * exp(PlogB[0] * lgt - PlogE[0] * rt_inv);
//   else if (pressure >= PlogP[np - 1])
//     kR = PlogA[np - 1] * exp(PlogB[np - 1] * lgt - PlogE[np - 1] * rt_inv);
//   else {
//     /* interpolate */
//     for (i = 0; i < np; i++) {
//       if (pressure <= PlogP[i])
//         break;
//     }

//     kR_l = PlogA[i - 1] * exp(PlogB[i - 1] * lgt - PlogE[i - 1] * rt_inv);
//     kR_r = PlogA[i] * exp(PlogB[i] * lgt - PlogE[i] * rt_inv);
//     kR = exp(log(MAX_C(kR_l, 1e-60)) +
//              (log(MAX_C(kR_r, 1e-60)) - log(MAX_C(kR_l, 1e-60))) /
//                  (log(PlogP[i]) - log(PlogP[i - 1])) *
//                  (log(pressure) - log(PlogP[i - 1])));
//   }
//   return MIN_C(kR, DBL_MAX);
// }

const char* Zhao2008Dme::GetMechanismResource() const {
  return R"RESOURCE(ELEMENTS
H  O  C  AR N  HE 
END
SPECIES
H                 H2                CH2               CH2S              CH3               O                 CH4               OH                
H2O               C2H               C2H2              C2H3              CO                N2                C2H4              HCO               
C2H5              CH2O              C2H6              CH2OH             CH3O              O2                CH3OH             HO2               
H2O2              HCCO              CH2CO             HCCOH             CH2HCO            CH3CO             CO2               CH3HCO            
OCHO              CH3CHOH           C2H4OH            CH3CH2O           CH3OCH2           HCOOH             CH3OCH3           C2H5OH            
HOCH2O            CH3OCO            CH3OCHO           CH3OCH2O          CH3OCH2OH         OCH2OCHO          HOCH2OCO          HOC2H4O2          
CH3OCH2O2         CH2OCH2O2H        CH3OCH2O2H        HO2CH2OCHO        O2CH2OCH2O2H      AR                HE                
END
REACTIONS
H+O2<=>O+OH                              3.547e+15   -0.406  16599.00
O+H2<=>H+OH                              5.080e+04    2.670   6290.00
H2+OH<=>H2O+H                            2.160e+08    1.510   3430.00
O+H2O<=>2OH                              2.970e+06    2.020  13400.00
H2+M<=>2H+M                              4.577e+19   -1.400 104380.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.00/ HE/0.00/ 
H2+AR<=>2H+AR                            5.840e+18   -1.100 104380.00
H2+HE<=>2H+HE                            5.840e+18   -1.100 104380.00
2O+M<=>O2+M                              6.165e+15   -0.500      0.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.00/ HE/0.00/ 
2O+AR<=>O2+AR                            1.886e+13    0.000  -1788.00
2O+HE<=>O2+HE                            1.886e+13    0.000  -1788.00
O+H+M<=>OH+M                             4.714e+18   -1.000      0.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.75/ HE/0.75/ 
H+OH+M<=>H2O+M                           3.800e+22   -2.000      0.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.38/ HE/0.38/ 
H+O2(+M)<=>HO2(+M)                       1.475e+12    0.600      0.00
H2/2.00/ H2O/11.00/ CO/1.90/ O2/0.78/ CO2/3.80/ 
     LOW  /  6.366e+20   -1.720    524.80 /
!     FCCHECK/     0.2     1e-30    0.8     1e+30      0         0
     TROE/     0.8    1e-30     1e+30        /
HO2+H<=>H2+O2                            1.660e+13    0.000    823.00
HO2+H<=>2OH                              7.079e+13    0.000    295.00
HO2+O<=>O2+OH                            3.250e+13    0.000      0.00
HO2+OH<=>H2O+O2                          2.890e+13    0.000   -497.00
2HO2<=>H2O2+O2                           4.200e+14    0.000  11982.00
     DUPLICATE
2HO2<=>H2O2+O2                           1.300e+11    0.000  -1629.30
     DUPLICATE
H2O2(+M)<=>2OH(+M)                       2.951e+14    0.000  48430.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.64/ HE/0.64/ 
     LOW  /  1.202e+17    0.000  45500.00 /
!     FCCHECK/     0.5     1e-30    0.5     1e+30      0         0
     TROE/     0.5    1e-30     1e+30        /
H2O2+H<=>H2O+OH                          2.410e+13    0.000   3970.00
H2O2+H<=>HO2+H2                          4.820e+13    0.000   7950.00
H2O2+O<=>OH+HO2                          9.550e+06    2.000   3970.00
H2O2+OH<=>HO2+H2O                        1.000e+12    0.000      0.00
     DUPLICATE
H2O2+OH<=>HO2+H2O                        5.800e+14    0.000   9557.00
     DUPLICATE
CO+O(+M)<=>CO2(+M)                       1.800e+10    0.000   2384.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.87/ 
     LOW  /  1.550e+24   -2.790   4191.00 /
!     FCCHECK/       0         0      0         0      1         0
     TROE/       1     1.00 10000000.00 10000000.00 /
CO+O2<=>CO2+O                            2.530e+12    0.000  47700.00
CO+HO2<=>CO2+OH                          3.010e+13    0.000  23000.00
CO+OH<=>CO2+H                            2.229e+05    1.890  -1158.70
HCO+M<=>H+CO+M                           4.748e+11    0.659  14874.00
H2/2.50/ H2O/6.00/ CO/1.90/ CO2/3.80/ 
HCO+O2<=>CO+HO2                          7.580e+12    0.000    410.00
HCO+H<=>CO+H2                            7.230e+13    0.000      0.00
HCO+O<=>CO+OH                            3.020e+13    0.000      0.00
HCO+OH<=>CO+H2O                          3.020e+13    0.000      0.00
HCO+O<=>CO2+H                            3.000e+13    0.000      0.00
HCO+HO2<=>CO2+OH+H                       3.000e+13    0.000      0.00
2HCO<=>H2+2CO                            3.000e+12    0.000      0.00
HCO+CH3<=>CO+CH4                         2.650e+13    0.000      0.00
2HCO<=>CH2O+CO                           3.000e+13    0.000      0.00
CH2O+M<=>HCO+H+M                         3.300e+39   -6.300  99900.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.70/ 
CH2O+M<=>CO+H2+M                         3.100e+45   -8.000  97510.00
H2/2.50/ H2O/12.00/ CO/1.90/ CO2/3.80/ AR/0.70/ 
CH2O+H<=>HCO+H2                          5.740e+07    1.900   2748.60
CH2O+O<=>HCO+OH                          1.810e+13    0.000   3080.00
CH2O+OH<=>HCO+H2O                        3.430e+09    1.180   -447.00
CH2O+O2<=>HCO+HO2                        1.230e+06    3.000  52000.00
CH2O+HO2<=>HCO+H2O2                      4.110e+04    2.500  10210.00
CH2O+CH3<=>HCO+CH4                       3.636e-06    5.420    998.00
CH3+O<=>CH2O+H                           8.430e+13    0.000      0.00
CH3+O2<=>CH3O+O                          1.990e+18   -1.570  29230.00
CH3+O2<=>CH2O+OH                         3.740e+11    0.000  14640.00
CH3+HO2<=>CH3O+OH                        2.410e+10    0.760  -2325.00
2CH3(+M)<=>C2H6(+M)                      2.277e+15   -0.690    174.86
H2O/5.00/ CO/2.00/ CO2/3.00/ 
     LOW  /  8.054e+31   -3.750    981.60 /
     TROE/       0      570         0     1e+30 /
CH3+H(+M)<=>CH4(+M)                      1.270e+16   -0.630    383.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  2.477e+33   -4.760   2440.00 /
     TROE/   0.783       74      2941      6964 /
CH4+H<=>CH3+H2                           5.470e+07    1.970  11210.00
CH4+O<=>CH3+OH                           3.150e+12    0.500  10290.00
CH4+OH<=>CH3+H2O                         5.720e+06    1.960   2639.00
CH3+HO2<=>CH4+O2                         3.160e+12    0.000      0.00
CH4+HO2<=>CH3+H2O2                       1.810e+11    0.000  18580.00
CH2OH+M<=>CH2O+H+M                       1.000e+14    0.000  25100.00

CH2OH+H<=>CH2O+H2                        6.000e+12    0.000      0.00
CH2OH+H<=>CH3+OH                         9.635e+13    0.000      0.00
CH2OH+O<=>CH2O+OH                        4.200e+13    0.000      0.00
CH2OH+OH<=>CH2O+H2O                      2.400e+13    0.000      0.00
CH2OH+O2<=>CH2O+HO2                      2.410e+14    0.000   5017.00
     DUPLICATE
CH2OH+O2<=>CH2O+HO2                      1.510e+15   -1.000      0.00
     DUPLICATE
CH2OH+HO2<=>CH2O+H2O2                    1.200e+13    0.000      0.00
CH2OH+HCO<=>CH3OH+CO                     1.000e+13    0.000      0.00
CH2OH+HCO<=>2CH2O                        1.500e+13    0.000      0.00
2CH2OH<=>CH3OH+CH2O                      3.000e+12    0.000      0.00
CH2OH+CH3O<=>CH3OH+CH2O                  2.400e+13    0.000      0.00
CH3O+M<=>CH2O+H+M                        8.300e+17   -1.200  15500.00

CH3O+H<=>CH3+OH                          3.200e+13    0.000      0.00
CH3O+O<=>CH2O+OH                         6.000e+12    0.000      0.00
CH3O+OH<=>CH2O+H2O                       1.800e+13    0.000      0.00
CH3O+O2<=>CH2O+HO2                       9.033e+13    0.000  11980.00
     DUPLICATE
CH3O+O2<=>CH2O+HO2                       2.200e+10    0.000   1748.00
     DUPLICATE
CH3O+HO2<=>CH2O+H2O2                     3.000e+11    0.000      0.00
CH3O+CO<=>CH3+CO2                        1.600e+13    0.000  11800.00
CH3O+HCO<=>CH3OH+CO                      9.000e+13    0.000      0.00
2CH3O<=>CH3OH+CH2O                       6.000e+13    0.000      0.00
OH+CH3(+M)<=>CH3OH(+M)                   2.790e+18   -1.430   1330.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ 
     LOW  /  4.000e+36   -5.920   3140.00 /
     TROE/   0.412      195      5900      6394 /
H+CH2OH(+M)<=>CH3OH(+M)                  1.055e+12    0.500     86.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ 
     LOW  /  4.360e+31   -4.650   5080.00 /
     TROE/     0.6      100     90000     10000 /
H+CH3O(+M)<=>CH3OH(+M)                   2.430e+12    0.515     50.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ 
     LOW  /  4.660e+41   -7.440  14080.00 /
     TROE/     0.7      100     90000     10000 /
CH3OH+H<=>CH2OH+H2                       3.200e+13    0.000   6095.00
CH3OH+H<=>CH3O+H2                        8.000e+12    0.000   6095.00
CH3OH+O<=>CH2OH+OH                       3.880e+05    2.500   3080.00
CH3OH+OH<=>CH3O+H2O                      1.000e+06    2.100    496.70
CH3OH+OH<=>CH2OH+H2O                     7.100e+06    1.800   -596.00
CH3OH+O2<=>CH2OH+HO2                     2.050e+13    0.000  44900.00
CH3OH+HCO<=>CH2OH+CH2O                   9.635e+03    2.900  13110.00
CH3OH+HO2<=>CH2OH+H2O2                   3.980e+13    0.000  19400.00
CH3OH+CH3<=>CH2OH+CH4                    3.190e+01    3.170   7172.00
CH3O+CH3OH<=>CH3OH+CH2OH                 3.000e+11    0.000   4060.00
2CH3<=>H+C2H5                            4.990e+12    0.100  10600.00
CH4+CH2<=>2CH3                           2.460e+06    2.000   8270.00
CH4+CH2S<=>2CH3                          1.600e+13    0.000   -570.00
CH3+OH<=>CH2+H2O                         5.600e+07    1.600   5420.00
CH3+OH<=>CH2S+H2O                        2.501e+13    0.000      0.00
CH3+CH2<=>C2H4+H                         4.000e+13    0.000      0.00
CH3+CH2S<=>C2H4+H                        1.200e+13    0.000   -570.00
CH3O+H<=>CH2S+H2O                        1.600e+13    0.000      0.00
CH2S+H2O(+M)<=>CH3OH(+M)                 4.820e+17   -1.160   1145.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ 
     LOW  /  1.880e+38   -6.360   5040.00 /
     TROE/  0.6027      208      3922     10180 /
C2H6+H<=>C2H5+H2                         1.150e+08    1.900   7530.00
C2H6+O<=>C2H5+OH                         8.980e+07    1.920   5690.00
C2H6+OH<=>C2H5+H2O                       3.540e+06    2.120    870.00
C2H6+O2<=>C2H5+HO2                       4.000e+13    0.000  50900.00
C2H6+HO2<=>C2H5+H2O2                     2.940e+11    0.000  14940.00
C2H6+CH3<=>C2H5+CH4                      6.140e+06    1.740  10450.00
C2H5+H(+M)<=>C2H6(+M)                    5.210e+17   -0.990   1580.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  1.990e+41   -7.080   6685.00 /
     TROE/  0.8422      125      2219      6882 /
C2H5+H<=>C2H4+H2                         2.000e+12    0.000      0.00
C2H5+O<=>CH3+CH2O                        1.320e+14    0.000      0.00
C2H5+O2<=>C2H4+HO2                       2.000e+10    0.000      0.00
2C2H5<=>C2H4+C2H6                        1.400e+12    0.000      0.00
C2H5+HCO<=>C2H6+CO                       1.200e+14    0.000      0.00
C2H5+O<=>CH3HCO+H                        8.020e+13    0.000      0.00
C2H4(+M)<=>H2+C2H2(+M)                   8.000e+12    0.440  88770.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  7.000e+50   -9.310  99860.00 /
     TROE/  0.7345      180      1035      5417 /
C2H4+H(+M)<=>C2H5(+M)                    1.080e+12    0.454   1820.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  1.200e+42   -7.620   6970.00 /
     TROE/  0.9753      210       984      4374 /
C2H3+H(+M)<=>C2H4(+M)                    6.080e+12    0.270    280.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  1.400e+30   -3.860   3320.00 /
     TROE/   0.782    207.5      2663      6095 /
C2H4+H<=>C2H3+H2                         1.325e+06    2.530  12240.00
C2H4+OH<=>C2H3+H2O                       1.800e+06    2.000   2500.00
C2H4+CH3<=>C2H3+CH4                      2.270e+05    2.000   9200.00
C2H4+O<=>CH3+HCO                         1.920e+07    1.830    220.00
C2H3+OH<=>C2H2+H2O                       5.000e+12    0.000      0.00
C2H4+O<=>OH+C2H3                         1.510e+07    1.910   3740.00
C2H4+O2<=>C2H3+HO2                       4.215e+13    0.000  57600.00
C2H3+H<=>C2H2+H2                         9.640e+13    0.000      0.00
C2H3+O<=>CH2CO+H                         6.000e+13    0.000      0.00
C2H3+H2O2<=>C2H4+HO2                     1.210e+10    0.000   -596.00
C2H3+CH3<=>C2H2+CH4                      3.900e+11    0.000      0.00
2C2H3<=>C2H4+C2H2                        9.600e+11    0.000      0.00
C2H3+O2<=>HCO+CH2O                       4.580e+16   -1.390   1015.00
C2H3+O2<=>HO2+C2H2                       1.337e+06    1.610   -384.00
C2H3+O2<=>O+CH2HCO                       1.000e+11    0.290     11.00
C2H2+H(+M)<=>C2H3(+M)                    5.600e+12    0.000   2400.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  3.800e+40   -7.270   7220.00 /
     TROE/  0.7507     98.5      1302      4167 /
C2H2+O<=>HCCO+H                          1.632e+07    2.000   1900.00
C2H2+O<=>C2H+OH                          4.600e+19   -1.410  28950.00
C2H2+O<=>CH2+CO                          4.080e+06    2.000   1900.00
C2H2+OH<=>CH2CO+H                        2.180e-04    4.500  -1000.00
C2H2+OH<=>HCCOH+H                        5.040e+05    2.300  13500.00
C2H2+OH<=>C2H+H2O                        3.370e+07    2.000  14000.00
C2H2+OH<=>CH3+CO                         4.830e-04    4.000  -2000.00
C2H2+HO2<=>CH2CO+OH                      6.030e+09    0.000   7944.00
HCCO+H<=>CH2S+CO                         1.000e+14    0.000      0.00
HCCO+O<=>H+2CO                           1.000e+14    0.000      0.00
HCCO+O2<=>OH+2CO                         1.600e+12    0.000    854.00
HCCO+OH<=>HCO+H+CO                       1.000e+13    0.000      0.00
HCCO+CH2<=>C2H3+CO                       3.000e+13    0.000      0.00
2HCCO<=>C2H2+2CO                         1.000e+13    0.000      0.00
CH3+HCCO<=>C2H4+CO                       5.000e+13    0.000      0.00
C2H+H(+M)<=>C2H2(+M)                     1.000e+17   -1.000      0.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  3.750e+33   -4.800   1900.00 /
     TROE/  0.6464      132      1315      5566 /
C2H+OH<=>H+HCCO                          2.000e+13    0.000      0.00
C2H+O2<=>HCO+CO                          5.000e+13    0.000   1500.00
C2H+H2<=>H+C2H2                          4.900e+05    2.500    560.00
CH2+H(+M)<=>CH3(+M)                      2.500e+16   -0.800      0.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  3.200e+27   -3.140   1230.00 /
     TROE/    0.68       78      1995      5590 /
CH2+O<=>HCO+H                            8.000e+13    0.000      0.00
CH2+OH<=>CH2O+H                          2.000e+13    0.000      0.00
CH2+H2<=>H+CH3                           5.000e+05    2.000   7230.00
CH2+O2<=>HCO+OH                          1.320e+13    0.000   1500.00
CH2+HO2<=>CH2O+OH                        2.000e+13    0.000      0.00
CH2+CO(+M)<=>CH2CO(+M)                   8.100e+11    0.500   4510.00
H2/2.00/ CH4/2.00/ H2O/6.00/ CO/1.50/ C2H6/3.00/ CO2/2.00/ AR/0.70/ 
     LOW  /  2.690e+33   -5.110   7095.00 /
     TROE/  0.5907      275      1226      5185 /
2CH2<=>C2H2+H2                           3.200e+13    0.000      0.00
CH2S+M<=>CH2+M                           9.000e+12    0.000    600.00
H2O/0.00/ CO/0.00/ CO2/0.00/ AR/0.00/ 
CH2S+H2O<=>CH2+H2O                       3.000e+13    0.000      0.00
CH2S+CO<=>CH2+CO                         9.000e+12    0.000      0.00
CH2S+CO2<=>CH2+CO2                       7.000e+12    0.000      0.00
CH2S+AR<=>CH2+AR                         9.000e+12    0.000    600.00
CH2S+O<=>CO+H2                           1.500e+13    0.000      0.00
CH2S+O<=>HCO+H                           1.500e+13    0.000      0.00
CH2S+OH<=>CH2O+H                         3.000e+13    0.000      0.00
CH2S+H2<=>CH3+H                          7.000e+13    0.000      0.00
CH2S+O2<=>H+OH+CO                        2.800e+13    0.000      0.00
CH2S+O2<=>CO+H2O                         1.200e+13    0.000      0.00
CH2S+CO2<=>CH2O+CO                       1.400e+13    0.000      0.00
CH3HCO<=>CH3+HCO                         7.000e+15    0.000  81674.00
CH3CO(+M)<=>CH3+CO(+M)                   3.000e+12    0.000  16722.00

     LOW  /  1.200e+15    0.000  12518.00 /
!     FCCHECK/       0         0      0         0      1         0
     TROE/       1     1.00 10000000.00 10000000.00 /
CH3HCO+OH<=>CH3CO+H2O                    3.370e+12    0.000   -620.00
CH3HCO+OH<=>CH2HCO+H2O                   3.370e+11    0.000   -620.00
CH3HCO+O<=>CH3CO+OH                      1.770e+18   -1.900   2975.00
CH3HCO+O<=>CH2HCO+OH                     3.720e+13   -0.200   3556.00
CH3HCO+H<=>CH3CO+H2                      4.660e+13   -0.350   2988.00
CH3HCO+H<=>CH2HCO+H2                     1.850e+12    0.400   5359.00
CH3HCO+CH3<=>CH3CO+CH4                   3.900e-07    5.800   2200.00
CH3HCO+CH3<=>CH2HCO+CH4                  2.450e+01    3.150   5727.00
CH3HCO+HO2<=>CH3CO+H2O2                  3.600e+19   -2.200  14030.00
CH3HCO+HO2<=>CH2HCO+H2O2                 2.320e+11    0.400  14864.00
CH3HCO+O2<=>CH3CO+HO2                    1.000e+14    0.000  42200.00
CH2CO+O<=>CO2+CH2                        1.750e+12    0.000   1350.00
CH2CO+H<=>CH3+CO                         2.710e+04    2.750    714.00
CH2CO+H<=>HCCO+H2                        2.000e+14    0.000   8000.00
CH2CO+O<=>HCCO+OH                        1.000e+13    0.000   8000.00
CH2CO+OH<=>HCCO+H2O                      1.000e+13    0.000   2000.00
CH2CO+OH<=>CH2OH+CO                      3.730e+12    0.000  -1013.00
CH2HCO+H<=>CH3+HCO                       5.000e+13    0.000      0.00
CH2HCO+H<=>CH2CO+H2                      2.000e+13    0.000      0.00
CH2HCO+O<=>CH2O+HCO                      1.000e+14    0.000      0.00
CH2HCO+OH<=>CH2CO+H2O                    3.000e+13    0.000      0.00
CH2HCO+O2<=>CH2O+CO+OH                   3.000e+10    0.000      0.00
CH2HCO+CH3<=>C2H5+CO+H                   4.900e+14   -0.500      0.00
CH2HCO+HO2<=>CH2O+HCO+OH                 7.000e+12    0.000      0.00
CH2HCO+HO2<=>CH3HCO+O2                   3.000e+12    0.000      0.00
CH2HCO<=>CH3+CO                          1.170e+43   -9.830  43756.00
CH2HCO<=>CH2CO+H                         1.810e+43   -9.610  45868.00
C2H5OH<=>CH3+CH2OH                       1.260e+51  -10.590 100869.00
C2H5OH<=>C2H4+H2O                        8.800e+25   -3.680  70799.00
C2H5OH+OH<=>C2H4OH+H2O                   1.810e+11    0.390    716.50
C2H5OH+OH<=>CH3CHOH+H2O                  3.090e+10    0.490   -379.80
C2H5OH+OH<=>CH3CH2O+H2O                  1.050e+10    0.790    716.90
C2H5OH+H<=>C2H4OH+H2                     1.900e+07    1.800   5098.00
C2H5OH+H<=>CH3CHOH+H2                    2.580e+07    1.650   2827.00
C2H5OH+H<=>CH3CH2O+H2                    1.500e+07    1.600   3038.00
C2H5OH+O<=>C2H4OH+OH                     9.410e+07    1.700   5459.00
C2H5OH+O<=>CH3CHOH+OH                    1.880e+07    1.850   1824.00
C2H5OH+O<=>CH3CH2O+OH                    1.580e+07    2.000   4448.00
C2H5OH+CH3<=>C2H4OH+CH4                  2.190e+02    3.180   9622.00
C2H5OH+CH3<=>CH3CHOH+CH4                 7.280e+02    2.990   7948.00
C2H5OH+CH3<=>CH3CH2O+CH4                 1.450e+02    2.990   7649.00
C2H5OH+HO2<=>CH3CHOH+H2O2                8.200e+03    2.550  10750.00
C2H5OH+HO2<=>C2H4OH+H2O2                 2.430e+04    2.550  15750.00
C2H5OH+HO2<=>CH3CH2O+H2O2                3.800e+12    0.000  24000.00
CH3CH2O+M<=>CH3HCO+H+M                   5.600e+34   -5.890  25274.00

CH3CH2O+M<=>CH3+CH2O+M                   5.350e+37   -6.960  23800.00

CH3CH2O+O2<=>CH3HCO+HO2                  4.000e+10    0.000   1100.00
CH3CH2O+CO<=>C2H5+CO2                    4.680e+02    3.160   5380.00
CH3CH2O+H<=>CH3+CH2OH                    3.000e+13    0.000      0.00
CH3CH2O+H<=>C2H4+H2O                     3.000e+13    0.000      0.00
CH3CH2O+OH<=>CH3HCO+H2O                  1.000e+13    0.000      0.00
CH3CHOH+O2<=>CH3HCO+HO2                  4.820e+13    0.000   5017.00
     DUPLICATE
CH3CHOH+O2<=>CH3HCO+HO2                  8.430e+14   -1.200      0.00
     DUPLICATE
CH3CHOH+O<=>CH3HCO+OH                    1.000e+14    0.000      0.00
CH3CHOH+H<=>C2H4+H2O                     3.000e+13    0.000      0.00
CH3CHOH+H<=>CH3+CH2OH                    3.000e+13    0.000      0.00
CH3CHOH+HO2<=>CH3HCO+2OH                 4.000e+13    0.000      0.00
CH3CHOH+OH<=>CH3HCO+H2O                  5.000e+12    0.000      0.00
CH3CHOH+M<=>CH3HCO+H+M                   1.000e+14    0.000  25000.00

C2H4OH+O2<=>HOC2H4O2                     1.000e+12    0.000  -1100.00
HOC2H4O2<=>2CH2O+OH                      1.800e+11    0.000  24500.00
C2H4+OH<=>C2H4OH                         2.410e+11    0.000  -2385.00
C2H5+HO2<=>CH3CH2O+OH                    4.000e+13    0.000      0.00
CH3OCH3<=>CH3+CH3O                       1.699e+42   -7.954  91806.60
CH3OCH3+OH<=>CH3OCH2+H2O                 6.710e+06    2.000   -629.88
CH3OCH3+H<=>CH3OCH2+H2                   2.970e+07    2.000   4033.61
CH3OCH3+CH3<=>CH3OCH2+CH4                2.680e+01    3.778   9631.30
CH3OCH3+O<=>CH3OCH2+OH                   1.855e-03    5.290   -109.00
CH3OCH3+HO2<=>CH3OCH2+H2O2               2.000e+13    0.000  16500.00
CH3OCH3+O2<=>CH3OCH2+HO2                 4.100e+13    0.000  44910.00
CH3OCH3+CH3O<=>CH3OCH2+CH3OH             6.020e+11    0.000   4074.00
CH3OCH3+CH3OCH2O2<=>CH3OCH2+CH3OCH2O2H   5.000e+12    0.000  17690.00
CH3OCH2<=>CH2O+CH3                       1.200e+13    0.000  25750.00
CH3OCH2+CH3O<=>CH3OCH3+CH2O              2.410e+13    0.000      0.00
CH3OCH2+CH2O<=>CH3OCH3+HCO               5.490e+03    2.800   5862.00
CH3OCH2+HO2<=>CH3OCH2O+OH                9.000e+12    0.000      0.00
CH3OCH2O<=>CH3OCHO+H                     1.745e+16   -0.660  11720.00
CH3OCHO<=>CH3+OCHO                       1.392e+18   -0.990  79140.00
CH3OCHO+O2<=>CH3OCO+HO2                  1.000e+13    0.000  49700.00
CH3OCHO+OH<=>CH3OCO+H2O                  2.340e+07    1.610    -35.00
CH3OCHO+HO2<=>CH3OCO+H2O2                1.220e+12    0.000  17000.00
CH3OCHO+O<=>CH3OCO+OH                    2.350e+05    2.500   2230.00
CH3OCHO+H<=>CH3OCO+H2                    4.550e+06    2.000   5000.00
CH3OCHO+CH3<=>CH3OCO+CH4                 7.550e-01    3.460   5481.00
CH3OCHO+CH3O<=>CH3OCO+CH3OH              5.480e+11    0.000   5000.00
CH3OCO<=>CH3O+CO                         7.451e+12   -1.760  17150.00
CH3OCO<=>CH3+CO2                         1.514e+12   -1.780  13820.00
OCHO+M<=>H+CO2+M                         2.443e+15   -0.500  26500.00

CH3OCH2+O2<=>CH3OCH2O2                   2.000e+12    0.000      0.00
CH3OCH2O2+CH2O<=>CH3OCH2O2H+HCO          1.000e+12    0.000  11670.00
2CH3OCH2O2<=>O2+2CH3OCH2O                1.597e+23   -4.500      0.00
2CH3OCH2O2<=>O2+CH3OCHO+CH3OCH2OH        6.844e+22   -4.500      0.00
CH3OCH2O2H<=>CH3OCH2O+OH                 2.106e+22   -2.120  43830.00
CH3OCH2O<=>CH3O+CH2O                     9.722e+15   -1.100  20640.00
CH3OCH2O+O2<=>CH3OCHO+HO2                5.000e+10    0.000    500.00
CH3OCH2O2<=>CH2OCH2O2H                   6.000e+10    0.000  21500.00
CH2OCH2O2H<=>OH+2CH2O                    1.500e+13    0.000  20500.00
CH2OCH2O2H+O2<=>O2CH2OCH2O2H             7.000e+11    0.000      0.00
O2CH2OCH2O2H<=>HO2CH2OCHO+OH             4.000e+10    0.000  18500.00
HO2CH2OCHO<=>OCH2OCHO+OH                 3.000e+16    0.000  40000.00
OCH2OCHO<=>HOCH2OCO                      1.000e+11    0.000  14000.00
HOCH2OCO<=>HOCH2O+CO                     2.177e+16   -2.690  17200.00
HOCH2OCO<=>CH2OH+CO2                     5.311e+15   -2.610  20810.00
HOCH2O<=>HCOOH+H                         1.000e+14    0.000  14900.00
CH2O+OH<=>HOCH2O                         4.500e+15   -1.110      0.00
HCOOH+M<=>CO+H2O+M                       2.300e+13    0.000  50000.00

HCOOH+M<=>CO2+H2+M                       1.500e+16    0.000  57000.00

HCOOH<=>HCO+OH                           4.593e+18   -0.460 108300.00
HCOOH+OH<=>H2O+CO2+H                     2.620e+06    2.060    916.00
HCOOH+OH<=>H2O+CO+OH                     1.850e+07    1.510   -962.00
HCOOH+H<=>H2+CO2+H                       4.240e+06    2.100   4868.00
HCOOH+H<=>H2+CO+OH                       6.030e+13   -0.350   2988.00
HCOOH+CH3<=>CH4+CO+OH                    3.900e-07    5.800   2200.00
HCOOH+HO2<=>H2O2+CO+OH                   1.000e+12    0.000  11920.00
HCOOH+O<=>CO+2OH                         1.770e+18   -1.900   2975.00
END
)RESOURCE";
}

const char* Zhao2008Dme::GetThermoResource() const {
  return R"RESOURCE(THERMO
   300.000  1000.000  5000.000
H                 000000H   1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54716300E+04-4.60117600E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54716300E+04-4.60117600E-01                   4
H2                000000H   2               G    300.00   5000.00 1000.00      1
 2.99142300E+00 7.00064400E-04-5.63382900E-08-9.23157800E-12 1.58275200E-15    2
-8.35034000E+02-1.35511000E+00 3.29812400E+00 8.24944200E-04-8.14301500E-07    3
-9.47543400E-11 4.13487200E-13-1.01252100E+03-3.29409400E+00                   4
CH2               000000H   2C   1          G    300.00   5000.00 1000.00      1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2S              000000H   2C   1          G    300.00   5000.00 1000.00      1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               000000H   3C   1          G    300.00   5000.00 1000.00      1
 2.97812060E+00 5.79785200E-03-1.97558000E-06 3.07297900E-10-1.79174160E-14    2
 1.65095130E+04 4.72247990E+00 3.65717970E+00 2.12659790E-03 5.45838830E-06    3
-6.61810030E-09 2.46570740E-12 1.64227160E+04 1.67353540E+00                   4
O                 000000O   1               G    300.00   5000.00 1000.00      1
 2.54206000E+00-2.75506200E-05-3.10280300E-09 4.55106700E-12-4.36805200E-16    2
 2.92308000E+04 4.92030800E+00 2.94642900E+00-1.63816600E-03 2.42103200E-06    3
-1.60284300E-09 3.89069600E-13 2.91476400E+04 2.96399500E+00                   4
CH4               000000H   4C   1          G    300.00   5000.00 1000.00      1
 1.68347900E+00 1.02372400E-02-3.87512900E-06 6.78558500E-10-4.50342300E-14    2
-1.00807900E+04 9.62339500E+00 7.78741500E-01 1.74766800E-02-2.78340900E-05    3
 3.04970800E-08-1.22393100E-11-9.82522900E+03 1.37221900E+01                   4
OH                000000H   1O   1          G    300.00   5000.00 1000.00      1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01                   4
H2O               000000H   2O   1          G    300.00   5000.00 1000.00      1
 2.67214600E+00 3.05629300E-03-8.73026000E-07 1.20099600E-10-6.39161800E-15    2
-2.98992100E+04 6.86281700E+00 3.38684200E+00 3.47498200E-03-6.35469600E-06    3
 6.96858100E-09-2.50658800E-12-3.02081100E+04 2.59023300E+00                   4
C2H               000000H   1C   2          G    300.00   5000.00 1000.00      1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00                   4
C2H2              000000H   2C   2          G    300.00   5000.00 1000.00      1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              000000H   3C   2          G    300.00   5000.00 1000.00      1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
CO                000000O   1C   1          G    300.00   5000.00 1000.00      1
 3.02507800E+00 1.44268900E-03-5.63082800E-07 1.01858100E-10-6.91095200E-15    2
-1.42683500E+04 6.10821800E+00 3.26245200E+00 1.51194100E-03-3.88175500E-06    3
 5.58194400E-09-2.47495100E-12-1.43105400E+04 4.84889700E+00                   4
N2                000000N   2               G    300.00   5000.00 1000.00      1
 2.92664000E+00 1.48797700E-03-5.68476100E-07 1.00970400E-10-6.75335100E-15    2
-9.22797700E+02 5.98052800E+00 3.29867700E+00 1.40824000E-03-3.96322200E-06    3
 5.64151500E-09-2.44485500E-12-1.02090000E+03 3.95037200E+00                   4
C2H4              000000H   4C   2          G    300.00   5000.00 1000.00      1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
HCO               000000H   1O   1C   1     G    300.00   5000.00 1000.00      1
 3.55727100E+00 3.34557300E-03-1.33500600E-06 2.47057300E-10-1.71385100E-14    2
 3.91632400E+03 5.55229900E+00 2.89833000E+00 6.19914700E-03-9.62308400E-06    3
 1.08982500E-08-4.57488500E-12 4.15992200E+03 8.98361400E+00                   4
C2H5              000000H   5C   2          G    300.00   5000.00 1000.00      1
 4.28788140E+00 1.24338930E-02-4.41391190E-06 7.06541020E-10-4.20351360E-14    2
 1.20564550E+04 8.46025830E-01 4.30585800E+00-4.18336380E-03 4.97072700E-05    3
-5.99058740E-08 2.30484780E-11 1.28417140E+04 4.71002360E+00                   4
CH2O              000000H   2O   1C   1     G    300.00   5000.00 1000.00      1
 5.14819050E+00 2.86780160E-03-2.37826330E-07-1.61113030E-10 2.85667350E-14    2
-1.62301730E+04-5.12138130E+00 2.69626120E+00 4.92614230E-03 8.28264940E-07    3
-5.50381960E-10-3.96103260E-13-1.49707930E+04 9.46975990E+00                   4
C2H6              000000H   6C   2          G    300.00   5000.00 1000.00      1
 4.82593800E+00 1.38404300E-02-4.55725900E-06 6.72496700E-10-3.59816100E-14    2
-1.27177900E+04-5.23950700E+00 1.46253900E+00 1.54946700E-02 5.78050700E-06    3
-1.25783200E-08 4.58626700E-12-1.12391800E+04 1.44322900E+01                   4
CH2OH             000000H   3O   1C   1     G    300.00   5000.00 1000.00      1
 3.74691030E+00 8.86461210E-03-4.25807220E-06 1.00880400E-09-9.45015610E-14    2
-3.66648240E+03 5.42810950E+00 4.61197920E+00-3.12037600E-03 3.55316800E-05    3
-4.93793980E-08 2.20272470E-11-3.60407340E+03 2.83513990E+00                   4
CH3O              000000H   3O   1C   1     G    300.00   5000.00 1000.00      1
 3.77080000E+00 7.87149700E-03-2.65638400E-06 3.94443100E-10-2.11261600E-14    2
 1.27832500E+02 2.92957500E+00 2.10620400E+00 7.21659500E-03 5.33847200E-06    3
-7.37763600E-09 2.07561100E-12 9.78601100E+02 1.31521800E+01                   4
O2                000000O   2               G    300.00   5000.00 1000.00      1
 3.69757800E+00 6.13519700E-04-1.25884200E-07 1.77528100E-11-1.13643500E-15    2
-1.23393000E+03 3.18916600E+00 3.21293600E+00 1.12748600E-03-5.75615000E-07    3
 1.31387700E-09-8.76855400E-13-1.00524900E+03 6.03473800E+00                   4
CH3OH             000000H   4O   1C   1     G    300.00   5000.00 1000.00      1
 4.02906100E+00 9.37659300E-03-3.05025400E-06 4.35879300E-10-2.22472300E-14    2
-2.61579100E+04 2.37819600E+00 2.66011500E+00 7.34150800E-03 7.17005100E-06    3
-8.79319400E-09 2.39057000E-12-2.53534800E+04 1.12326300E+01                   4
HO2               000000H   1O   2          G    300.00   5000.00 1000.00      1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              000000H   2O   2          G    300.00   5000.00 1000.00      1
 4.57316700E+00 4.33613600E-03-1.47468900E-06 2.34890400E-10-1.43165400E-14    2
-1.80069600E+04 5.01137000E-01 3.38875400E+00 6.56922600E-03-1.48501300E-07    3
-4.62580600E-09 2.47151500E-12-1.76631500E+04 6.78536300E+00                   4
HCCO              000000H   1O   1C   2     G    300.00   5000.00 1000.00      1
 5.62820580E+00 4.08534010E-03-1.59345470E-06 2.86260520E-10-1.94078320E-14    2
 1.93272150E+04-3.93025950E+00 2.25172140E+00 1.76550210E-02-2.37291010E-05    3
 1.72757590E-08-5.06648110E-12 2.00594490E+04 1.24904170E+01                   4
CH2CO             000000H   2O   1C   2     G    300.00   5000.00 1000.00      1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01                   4
HCCOH             000000H  20O   1C   2     G    300.00   5000.00 1000.00      1
 5.92382910E+00 6.79236000E-03-2.56585640E-06 4.49878410E-10-2.99401010E-14    2
 7.26462600E+03-7.60177420E+00 1.24237330E+00 3.10722010E-02-5.08668640E-05    3
 4.31371310E-08-1.40145940E-11 8.03161430E+03 1.38743190E+01                   4
CH2HCO            000000H   3O   1C   2     G    300.00   5000.00 1000.00      1
 5.97566990E+00 8.13059140E-03-2.74362450E-06 4.07030410E-10-2.17601710E-14    2
 4.90321780E+02-5.03208790E+00 3.40906240E+00 1.07385740E-02 1.89149250E-06    3
-7.15858310E-09 2.86738510E-12 1.52147660E+03 9.57145350E+00                   4
CH3CO             000000H   3O   1C   2     G    300.00   5000.00 1000.00      1
 5.94477310E+00 7.86672050E-03-2.88658820E-06 4.72708750E-10-2.85998610E-14    2
-3.78730750E+03-5.01367510E+00 4.16342570E+00-2.32616100E-04 3.42678200E-05    3
-4.41052270E-08 1.72756120E-11-2.65745290E+03 7.34682800E+00                   4
CO2               000000O   2C   1          G    300.00   5000.00 1000.00      1
 4.45362300E+00 3.14016900E-03-1.27841100E-06 2.39399700E-10-1.66903300E-14    2
-4.89669600E+04-9.55395900E-01 2.27572500E+00 9.92207200E-03-1.04091100E-05    3
 6.86668700E-09-2.11728000E-12-4.83731400E+04 1.01884900E+01                   4
CH3HCO            000000H   4O   1C   2     G    300.00   5000.00 1000.00      1
 5.40411080E+00 1.17230590E-02-4.22631370E-06 6.83724510E-10-4.09848630E-14    2
-2.25931220E+04-3.48079170E+00 4.72945950E+00-3.19328580E-03 4.75349210E-05    3
-5.74586110E-08 2.19311120E-11-2.15728780E+04 4.10301590E+00                   4
OCHO              000000H   1O   2C   1     G    300.00   5000.00 1000.00      1
 6.12628782E+00 3.75602932E-03-1.42010352E-06 2.36429200E-10-1.44167651E-14    2
-2.17698466E+04-8.01574694E+00 1.35213452E+00 1.50082004E-02-1.09896141E-05    3
 3.73679840E-09-4.81014498E-13-2.02253647E+04 1.74373147E+01                   4
CH3CHOH           000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 6.76654240E+00 1.16344360E-02-3.77906510E-06 5.38288750E-10-2.73153450E-14    2
-5.60929690E+03-9.39804420E+00 2.48133280E+00 1.67900360E-02 3.77554990E-06    3
-1.39234970E-08 6.00951930E-12-4.01200540E+03 1.45816220E+01                   4
C2H4OH            000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 7.59440140E+00 9.32293390E-03-3.03038540E-06 4.32163190E-10-2.19700390E-14    2
-5.77278520E+03-1.39555720E+01 1.40195080E+00 2.15431750E-02-2.23265120E-06    3
-1.44640920E-08 8.04884200E-12-3.84645190E+03 1.91489810E+01                   4
CH3CH2O           000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 6.01143460E+00 1.21652190E-02-4.04496040E-06 5.90765880E-10-3.09695950E-14    2
-4.93669920E+03-6.79017980E+00 1.73025040E+00 1.69084890E-02 3.99962210E-06    3
-1.37111800E-08 5.76436030E-12-3.29224830E+03 1.73361150E+01                   4
CH3OCH2           000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 8.17137842E+00 1.10086181E-02-3.82352277E-06 5.99637202E-10-3.50317513E-14    2
-3.41941605E+03-1.78650856E+01 2.91327415E+00 2.03364659E-02-9.59712342E-06    3
 2.07478525E-09-1.71343362E-13-1.18844240E+03 1.16066817E+01                   4
HCOOH             000000H   2O   2C   1     G    300.00   5000.00 1000.00      1
 6.68733013E+00 5.14289368E-03-1.82238513E-06 2.89719163E-10-1.70892199E-14    2
-4.83995400E+04-1.13104798E+01 1.43548185E+00 1.63363016E-02-1.06257421E-05    3
 3.32132977E-09-4.02176103E-13-4.64616504E+04 1.72885798E+01                   4
CH3OCH3           000000H   6O   1C   2     G    300.00   5000.00 1000.00      1
 8.30815546E-01 2.69173263E-02-1.38874777E-05 3.47515079E-09-3.41706784E-13    2
-2.34120975E+04 2.02174360E+01 5.68097447E+00-5.39434751E-03 6.49472750E-05    3
-8.05065318E-08 3.27474018E-11-2.39755455E+04-6.36955496E-01                   4
C2H5OH            000000H   6O   1C   2     G    300.00   5000.00 1000.00      1
 6.56243650E+00 1.52042220E-02-5.38967950E-06 8.62250110E-10-5.12897870E-14    2
-3.15256210E+04-9.47302020E+00 4.85869570E+00-3.74017260E-03 6.95553780E-05    3
-8.86547960E-08 3.51688350E-11-2.99961320E+04 4.80185450E+00                   4
HOCH2O            000000H   3O   2C   1     G    300.00   5000.00 1000.00      1
 6.39521515E+00 7.43673043E-03-2.50422354E-06 3.84879712E-10-2.21778689E-14    2
-2.47500385E+04-7.29290847E+00 4.11183145E+00 7.53850697E-03 3.77337370E-06    3
-5.38746005E-09 1.45615887E-12-2.34414546E+04 6.81381989E+00                   4
CH3OCO            000000H   3O   2C   2     G    300.00   5000.00 1000.00      1
 1.30877600E+01 4.53544950E-03-1.65096364E-06 2.67197277E-10-1.59576863E-14    2
-2.46616400E+04-3.27914051E+01 3.94199159E+00 2.43434884E-02-1.65595560E-05    3
 4.58537411E-09-3.31795708E-13-2.14404829E+04 1.66954362E+01                   4
CH3OCHO           000000H   4O   2C   2     G    300.00   5000.00 1000.00      1
 8.69123518E+00 1.15503122E-02-4.27782486E-06 7.02533059E-10-4.24333552E-14    2
-4.64364769E+04-1.89301478E+01 3.08839783E+00 2.03760048E-02-6.84777040E-06    3
-7.28186203E-10 5.62130216E-13-4.41855167E+04 1.25364719E+01                   4
CH3OCH2O          000000H   5O   2C   2     G    300.00   5000.00 1000.00      1
 8.60261845E+00 1.35772195E-02-4.84661602E-06 7.77766193E-10-4.62633624E-14    2
-2.13762444E+04-1.75775023E+01 3.25889339E+00 2.22146359E-02-7.78556340E-06    3
-2.41484158E-10 4.51914496E-13-1.92377212E+04 1.23680069E+01                   4
CH3OCH2OH         000000H   6O   2C   2     G    300.00   5000.00 1000.00      1
 8.70981570E+00 1.53602372E-02-5.41003788E-06 8.60573446E-10-5.08819752E-14    2
-4.76607115E+04-1.80226702E+01 3.15851876E+00 2.44325751E-02-8.66984784E-06    3
-5.93319328E-11 4.36400003E-13-4.54488899E+04 1.30511235E+01                   4
OCH2OCHO          000000H   3O   3C   2     G    300.00   5000.00 1000.00      1
 1.20233916E+01 8.11262659E-03-2.91356462E-06 4.67340384E-10-2.77375525E-14    2
-4.33647231E+04-3.33691809E+01 5.19690837E+00 1.58839723E-02 3.53540547E-07    3
-6.10456923E-09 1.94661801E-12-4.02242792E+04 6.11645828E+00                   4
HOCH2OCO          000000H   3O   3C   2     G    300.00   5000.00 1000.00      1
 1.13737391E+01 8.17663898E-03-2.92034021E-06 4.66695616E-10-2.76276823E-14    2
-4.65575743E+04-2.86035265E+01 6.08180801E+00 1.28768359E-02 2.04419418E-06    3
-6.10154921E-09 1.79820559E-12-4.39526183E+04 2.54054449E+00                   4
HOC2H4O2          000000H   5O   3C   2     G    300.00   5000.00 1000.00      1
 1.00941573E+01 1.23879015E-02-3.73811683E-06 5.46874551E-10-3.09943951E-14    2
-2.37710522E+04-2.00956526E+01 4.44209543E+00 2.52880383E-02-1.51605275E-05    3
 5.24921198E-09-7.91470852E-13-2.17507126E+04 1.04122371E+01                   4
CH3OCH2O2         000000H   5O   3C   2     G    300.00   5000.00 1000.00      1
 1.24249729E+01 1.18705986E-02-4.07906532E-06 6.35310809E-10-3.69427867E-14    2
-2.29679238E+04-3.53740145E+01 2.21029612E+00 3.68877454E-02-2.82561555E-05    3
 1.15730533E-08-1.97130470E-12-1.94940940E+04 1.91463601E+01                   4
CH2OCH2O2H        000000H   5O   3C   2     G    300.00   5000.00 1000.00      1
 1.51191783E+01 9.23718883E-03-3.19127505E-06 4.99114678E-10-2.91162488E-14    2
-1.84114867E+04-4.85706618E+01 2.52895507E+00 4.24128290E-02-3.73406386E-05    3
 1.66639333E-08-2.96443312E-12-1.44293306E+04 1.76899251E+01                   4
CH3OCH2O2H        000000H   6O   3C   2     G    300.00   5000.00 1000.00      1
 1.49370964E+01 1.19465829E-02-4.12746359E-06 6.45422590E-10-3.76427939E-14    2
-4.11001068E+04-4.99552737E+01 1.19855761E+00 4.59060764E-02-3.66252420E-05    3
 1.49318970E-08-2.46057445E-12-3.65363161E+04 2.31339904E+01                   4
HO2CH2OCHO        000000H   4O   4C   2     G    300.00   5000.00 1000.00      1
 1.64584298E+01 8.52683511E-03-3.04113500E-06 4.85596908E-10-2.87316334E-14    2
-6.23959608E+04-5.38924139E+01 3.47935703E+00 4.02952392E-02-3.30109296E-05    3
 1.34360117E-08-2.18601580E-12-5.80629934E+04 1.52521392E+01                   4
O2CH2OCH2O2H      000000H   5O   5C   2     G    300.00   5000.00 1000.00      1
 1.92038046E+01 1.04394841E-02-3.60582939E-06 5.63792843E-10-3.28807214E-14    2
-3.79207055E+04-6.51847273E+01 1.99640551E+00 5.83226232E-02-5.53259778E-05    3
 2.59810540E-08-4.77141005E-12-3.27628742E+04 2.44215005E+01                   4
AR                000000AR  1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.36600100E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.36600100E+00                   4
HE                000000HE  1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 9.15348900E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.15348800E-01                   4
)RESOURCE";
}

} // namespace ideal_gas
} // namespace fub
