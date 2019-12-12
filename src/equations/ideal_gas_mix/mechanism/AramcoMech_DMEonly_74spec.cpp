#include "fub/ideal_gas/mechanism/AramcoMech_DMEonly_74spec.hpp"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace fub {
namespace ideal_gas {

static double GetLindRateCoeff(double temp, double pressure, double k0,
                               double kInf, double fc, double conc);

static double GetPlogRateCoeff(double temp, double pressure, double lgt,
                               double rt_inv, double* PlogP, double* PlogA,
                               double* PlogB, double* PlogE, int np);

static double MAX_C(double X1, double X2);
static double MIN_C(double X1, double X2);

void AramcoMech_DMEonly_74spec::ComputeProductionRates(
    span<double> cdot, span<double> w, span<double> k, span<double> c,
    span<double> M, double temp, double pressure) const {
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

  int nSpec = 76;
  int nSpecIn = 76;
  double kTroe0, kTroeInf, fcTroe;
  double RGAS = 8314.34;
  double lgt = log(temp);
  double rt_inv = 1.0 / (RGAS * temp);
  double PlogA[14], PlogB[14], PlogE[14], PlogP[14];
  int np;

  M[mM1] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1.4p+1 * c[sH2] +
           c[sO2] + 0x1.e666666666666p+0 * c[sCO] +
           0x1.e666666666666p+1 * c[sCO2] + 0x1.8p+3 * c[sH2O] + c[sH] + c[sO] +
           c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
           c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] +
           c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] +
           c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] +
           c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] +
           c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] + c[sCH2OCHO] +
           c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] + c[sOHY] +
           c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] + c[sC2H3CHO] +
           c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] + c[sCH2CHO] +
           c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] + c[sOCHO] +
           c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
           c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] +
           c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
           0x1.a8f5c28f5c28fp-1 * c[sHE] + c[sAR];

  M[mM2] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1.4p+1 * c[sH2] +
           c[sO2] + 0x1.e666666666666p+0 * c[sCO] +
           0x1.e666666666666p+1 * c[sCO2] + 0x1.8p+3 * c[sH2O] + c[sH] + c[sO] +
           c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
           c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] +
           c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] +
           c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] +
           c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] +
           c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] + c[sCH2OCHO] +
           c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] + c[sOHY] +
           c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] + c[sC2H3CHO] +
           c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] + c[sCH2CHO] +
           c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] + c[sOCHO] +
           c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
           c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] +
           c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
           0x1.a8f5c28f5c28fp-1 * c[sHE] + 0x1.a8f5c28f5c28fp-1 * c[sAR];

  M[mM3] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1.4p+1 * c[sH2] +
           c[sO2] + 0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+3 * c[sH2O] +
           c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] +
           c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] +
           c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] +
           c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] +
           c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] +
           c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] +
           c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] +
           c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] +
           c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] +
           c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
           c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] +
           c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] +
           c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] +
           c[sOCH2O2H] + 0x1.8p-1 * c[sHE] + 0x1.8p-1 * c[sAR];

  M[mM4] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1.75c28f5c28f5cp-1 * c[sH2] +
      c[sO2] + c[sCO] + c[sCO2] + 0x1.d333333333333p+1 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] +
      0x1.851eb851eb852p-2 * c[sAR];

  M[mM5] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1.4cccccccccccdp+0 * c[sH2] +
      c[sO2] + 0x1.e666666666666p+0 * c[sCO] + 0x1.e666666666666p+1 * c[sCO2] +
      0x1.4p+3 * c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
      0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
      c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] +
      c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] +
      c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] +
      c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] +
      c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] + c[sOHY] +
      c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] + c[sC2H3CHO] +
      c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] + c[sCH2CHO] +
      c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] + c[sOCHO] + c[sSC2H4OH] +
      c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] +
      c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] +
      c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H];

  M[mM6] = +c[sAR];

  M[mM7] = +c[sHE];

  M[mM8] = +c[sH2O];

  M[mM9] = 0x1.8p+0 * c[sN2] + c[sCH3OCH3] + c[sCH4] +
           0x1.d99999999999ap+1 * c[sH2] + 0x1.3333333333333p+0 * c[sO2] +
           0x1.6666666666666p+1 * c[sCO] + 0x1.999999999999ap+0 * c[sCO2] +
           c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + c[sC2H6] +
           c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] +
           0x1.ecccccccccccdp+2 * c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] +
           c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] +
           c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
           c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] +
           c[sCH] + c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] +
           c[sC2H3] + c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] +
           c[sHOCH2O2] + c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] +
           c[sC2H5O2] + c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] +
           c[sH2CC] + c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] +
           c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] +
           c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] +
           c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
           0x1.4cccccccccccdp-1 * c[sHE] + c[sAR];

  M[mM10] =
      c[sN2] + c[sCH3OCH3] + c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.cp+0 * c[sCO] + 0x1.ccccccccccccdp+1 * c[sCO2] + 0x1.8p+3 * c[sH2O] +
      c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM11] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+3 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM12] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM13] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM14] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM15] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM16] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM17] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM18] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM19] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM20] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM21] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM22] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + 0x1p+1 * c[sCO] +
            0x1.8p+1 * c[sCO2] + 0x1.4p+2 * c[sH2O] + c[sH] + c[sO] + c[sOH] +
            c[sCH3] + c[sCH2O] + c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] +
            c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] +
            c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] +
            c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
            c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] +
            c[sCH] + c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] +
            c[sC2H3] + c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] +
            c[sHOCH2O2] + c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] +
            c[sC2H5O2] + c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] +
            c[sH2CC] + c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] +
            c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] +
            c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] +
            c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM23] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM24] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
            0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
            c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] +
            c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] +
            c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] +
            c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] +
            c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] +
            c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] +
            c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] +
            c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] +
            c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] +
            c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
            c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] +
            c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] +
            c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] +
            c[sOCH2O2H] + c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM25] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM26] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM27] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM28] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM29] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM30] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM31] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM32] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM33] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM34] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM35] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
            0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
            c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] +
            c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] +
            c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] +
            c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] +
            c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] +
            c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] +
            c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] + c[sC2H5] +
            c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] + c[sC2H3CO] +
            c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] + c[sCH2CO] +
            c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
            c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] +
            c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] +
            c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] +
            c[sOCH2O2H] + c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM36] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM37] =
      c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
      0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
      c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] + c[sCH3OCH2] +
      c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] +
      c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] +
      c[sHCO] + c[sO2CHO] + c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] +
      c[sCH3O2H] + c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
      c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
      c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
      c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
      c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] + c[sC2H2] +
      c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] +
      c[sHCOH] + c[sC2H3OH] + c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] +
      c[sCH3COCH3] + c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] +
      0x1.6666666666666p-1 * c[sHE] + 0x1.6666666666666p-1 * c[sAR];

  M[mM38] = c[sN2] + c[sCH3OCH3] + 0x1p+1 * c[sCH4] + 0x1p+1 * c[sH2] + c[sO2] +
            0x1.8p+0 * c[sCO] + 0x1p+1 * c[sCO2] + 0x1.8p+2 * c[sH2O] + c[sH] +
            c[sO] + c[sOH] + c[sCH3] + c[sCH2O] + 0x1.8p+1 * c[sC2H6] +
            c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] + c[sHO2] +
            c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] + c[sCH3OCH2O] +
            c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] + c[sHOCHO] +
            c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] + c[sCH3O2] +
            c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] + c[sCH3OCHO] +
            c[sCH2OCHO] + c[sC2H4O1X2] + 0x1.4p+1 * c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            0x1.4p+1 * c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] +
            c[sCH3CO2] + c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] +
            c[sC2H5O] + c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] +
            c[sCH3COCH2] + c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM39] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  M[mM40] = c[sN2] + c[sCH3OCH3] + c[sCH4] + c[sH2] + c[sO2] + c[sCO] +
            c[sCO2] + c[sH2O] + c[sH] + c[sO] + c[sOH] + c[sCH3] + c[sCH2O] +
            c[sC2H6] + c[sCH3OCH2] + c[sO2CH2OCH2O2H] + c[sCH2OCH2O2H] +
            c[sHO2] + c[sH2O2] + c[sCH3OCH2O2] + c[sCH3O] + c[sCH3OH] +
            c[sCH3OCH2O] + c[sOCH2OCHO] + c[sHO2CH2OCHO] + c[sHCO] + c[sO2CHO] +
            c[sHOCHO] + c[sHOCH2O] + c[sHOCH2OCO] + c[sCH2OH] + c[sCH3O2H] +
            c[sCH3O2] + c[sHO2CHO] + c[sCH3OCH2O2H] + c[sCH2] + c[sCH] +
            c[sCH3OCHO] + c[sCH2OCHO] + c[sC2H4O1X2] + c[sC2H4] + c[sC2H3] +
            c[sC2H5] + c[sOHY] + c[sCH3OCO] + c[sHOCH2O2H] + c[sHOCH2O2] +
            c[sC2H3CO] + c[sC2H3CHO] + c[sCH3CO] + c[sC2H5O2H] + c[sC2H5O2] +
            c[sCH2CO] + c[sCH2CHO] + c[sCH3CO3H] + c[sCH3CO3] + c[sH2CC] +
            c[sC2H2] + c[sOCHO] + c[sSC2H4OH] + c[sC2H5OH] + c[sCH3CO2] +
            c[sCH3CHO] + c[sHCCO] + c[sHCOH] + c[sC2H3OH] + c[sC2H5O] +
            c[sC2H5CO] + c[sC2H5CHO] + c[sC2H] + c[sCH3COCH3] + c[sCH3COCH2] +
            c[sCH2Y] + c[sOCH2O2H] + c[sHE] + c[sAR];

  k[r1f] = 0x1.836e21p+36 * exp(-0x1.e7f348p+25 * rt_inv);
  k[r1b] = 0x1.998b126e2247ap+27 *
           exp(0x1.b9e296fd9458cp-2 * lgt + 0x1.91de0bc6a170cp+22 * rt_inv);
  k[r2f] = 0x1.9666666666666p+5 *
           exp(0x1.55c28f5c28f5cp+1 * lgt - 0x1.91b2e00000001p+24 * rt_inv);
  k[r2b] = 0x1.dd3457cf8f588p+4 *
           exp(0x1.517e317544529p+1 * lgt - 0x1.386aaa7262f73p+24 * rt_inv);
  k[r3f] = 0x1.4655decp+35 * exp(-0x1.be42d00000001p+24 * rt_inv);
  k[r3b] = 0x1.b95b6cbde50f2p+38 *
           exp(-0x1.a04cae3eb9b8p-4 * lgt - 0x1.5f0aa069bfe67p+26 * rt_inv);
  k[r4f] = 0x1.733ffffffffffp+11 *
           exp(0x1.028f5c28f5c29p+1 * lgt - 0x1.abbf200000001p+25 * rt_inv);
  k[r4b] = 0x1.4252b270def85p+7 *
           exp(0x1.0b4d63b407354p+1 * lgt + 0x1.7f634e6938d44p+23 * rt_inv);
  k[r5f] = 0x1.4537251ebd4p+55 *
           exp(-0x1.6666666666666p+0 * lgt - 0x1.a092f80000001p+28 * rt_inv);
  k[r5b] = 0x1.5aa81392ce7ebp+47 *
           exp(-0x1.bf25c679c4a82p+0 * lgt - 0x1.c27a52d2ccd8p+21 * rt_inv);
  k[r6f] = 0x1.6f766f4p+32 * exp(-0x1p-1 * lgt);
  k[r6b] = 0x1.7ef055a6f4cdep+48 *
           exp(-0x1.3c8402f3951p-1 * lgt - 0x1.dabf61309b311p+28 * rt_inv);
  k[r7f] = 0x1.126412e9p+42 * exp(-0x1p+0 * lgt);
  k[r7b] = 0x1.2e449e5f13778p+49 *
           exp(-0x1.5f92b774cbep-1 * lgt - 0x1.9779800180a71p+28 * rt_inv);
  k[r8f] = 0x1.f161421c8ep+54 * exp(-0x1p+1 * lgt);
  k[r8b] = 0x1.3b8a82906173bp+66 *
           exp(-0x1.c1456ad0883cp+0 * lgt - 0x1.d8ec7e74ca706p+28 * rt_inv);
  kTroe0 = 0x1.f98895c08p+43 * exp(-0x1.3ae147ae147aep+0 * lgt);
  kTroeInf = 0x1.15295e8p+32 * exp(0x1.c28f5c28f5c29p-2 * lgt);
  fcTroe = 0x1.51eb851eb851ep-2 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.570a3d70a3d71p-1 * exp(-temp / 0x1.93e5939a08ceap+99) +
           0x1p+0 * exp(-0x1.93e5939a08ceap+99 / temp);
  k[r9f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  kTroe0 = 0x1.6a6a1ff809bf5p+54 *
           exp(-0x1.3e8ceb4977bbp+0 * lgt - 0x1.86e54c1a783ecp+27 * rt_inv);
  kTroeInf = 0x1.8d64416721a96p+42 *
             exp(0x1.b3e0cdbb68c21p-2 * lgt - 0x1.86e54c1a783ecp+27 * rt_inv);
  k[r9b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  kTroe0 = 0x1.8c64e861p+42 * exp(-0x1.3333333333333p+0 * lgt);
  kTroeInf = 0x1.15295e8p+32 * exp(0x1.c28f5c28f5c29p-2 * lgt);
  fcTroe = 0x1.3333333333334p-2 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.6666666666666p-1 * exp(-temp / 0x1.93e5939a08ceap+99) +
           0x1p+0 * exp(-0x1.93e5939a08ceap+99 / temp);
  k[r10f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  kTroe0 = 0x1.1c2c4e1fc7f7dp+53 *
           exp(-0x1.36ded6ce9674p+0 * lgt - 0x1.86e54c1a783e8p+27 * rt_inv);
  kTroeInf = 0x1.8d64416721a71p+42 *
             exp(0x1.b3e0cdbb68bf5p-2 * lgt - 0x1.86e54c1a783e8p+27 * rt_inv);
  k[r10b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  kTroe0 = 0x1.0b85bca2p+43 * exp(-0x1.3333333333333p+0 * lgt);
  kTroeInf = 0x1.15295e8p+32 * exp(0x1.c28f5c28f5c29p-2 * lgt);
  fcTroe = 0x1.a3d70a3d70a3ep-2 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.2e147ae147ae1p-1 * exp(-temp / 0x1.93e5939a08ceap+99) +
           0x1p+0 * exp(-0x1.93e5939a08ceap+99 / temp);
  k[r11f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM7]);
  kTroe0 = 0x1.7f92302918746p+53 *
           exp(-0x1.36ded6ce9681p+0 * lgt - 0x1.86e54c1a783f4p+27 * rt_inv);
  kTroeInf = 0x1.8d644167224aep+42 *
             exp(0x1.b3e0cdbb688b5p-2 * lgt - 0x1.86e54c1a783f4p+27 * rt_inv);
  k[r11b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM7]);
  k[r12f] = 0x1.07b69ad8p+36 * exp(-0x1.2d568p+20 * rt_inv);
  k[r12b] = 0x1.ac5b38eeee987p+23 *
            exp(0x1.84b5db40c32cp-1 * lgt - 0x1.23dc9e8a53ffap+27 * rt_inv);
  k[r13f] = 0x1.02cccccccccccp+9 *
            exp(0x1.376c8b439581p+1 * lgt - 0x1.aaf6ea0000001p+27 * rt_inv);
  k[r13b] = 0x1.8b8765a38e658p+11 *
            exp(0x1.0937096c336e8p+1 * lgt + 0x1.07fa13479404cp+22 * rt_inv);
  k[r14f] = 0x1.e449a94p+34;
  k[r14b] = 0x1.74155ba0e0b57p+31 *
            exp(0x1.4f891f83f4f8p-2 * lgt - 0x1.a80db3e88910dp+27 * rt_inv);
  k[r15f] = 0x1.6df8f7p+34 * exp(0x1.fbad800000001p+20 * rt_inv);
  k[r15b] = 0x1.43dcc24aa03ecp+35 *
            exp(0x1.0998e32b6cp-2 * lgt - 0x1.137e2ae78e519p+28 * rt_inv);
  k[r16f] = 0x1.efe92p+26 * exp(0x1.a041400000001p+22 * rt_inv);
  k[r16b] = 0x1.6dd8ca4aa29cep+32 *
            exp(-0x1.f38a6281b8e8p-3 * lgt - 0x1.262cff24238a7p+27 * rt_inv);
  k[r17f] = 0x1.54ad8428p+38 * exp(-0x1.7f0e800000001p+25 * rt_inv);
  k[r17b] = 0x1.f6a796dc87bf9p+43 *
            exp(-0x1.f38a6281b38p-3 * lgt - 0x1.92f2a924238fdp+27 * rt_inv);
  kTroe0 = 0x1.f9825f0e08559p+73 *
           exp(-0x1.2666666666666p+1 * lgt - 0x1.8508af0000001p+27 * rt_inv);
  kTroeInf = 0x1.d1a94a2p+40 *
             exp(0x1.ccccccccccccdp-1 * lgt - 0x1.8508af0000001p+27 * rt_inv);
  fcTroe = 0x1.f5c28f5c28f5cp-2 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.051eb851eb852p-1 * exp(-temp / 0x1.93e5939a08ceap+99);
  k[r18f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM8]);
  kTroe0 = 0x1.84246132aca7dp+45 *
           exp(-0x1.4854ef40ce735p+0 * lgt + 0x1.389b4b447b80dp+23 * rt_inv);
  kTroeInf = 0x1.658bbeb8ff583p+12 *
             exp(0x1.eade43f264bfep+0 * lgt + 0x1.389b4b447b80cp+23 * rt_inv);
  k[r18b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM8]);
  kTroe0 = 0x1.0df7621ed4455p+71 *
           exp(-0x1.2666666666666p+1 * lgt - 0x1.8508af0000001p+27 * rt_inv);
  kTroeInf = 0x1.d1a94a2p+40 *
             exp(0x1.ccccccccccccdp-1 * lgt - 0x1.8508af0000001p+27 * rt_inv);
  fcTroe = 0x1.23d70a3d70a3ep-1 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.b851eb851eb85p-2 * exp(-temp / 0x1.93e5939a08ceap+99);
  k[r19f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM9]);
  kTroe0 = 0x1.9e92bf6461539p+42 *
           exp(-0x1.4854ef40ce621p+0 * lgt + 0x1.389b4b447b917p+23 * rt_inv);
  kTroeInf = 0x1.658bbeb8fe951p+12 *
             exp(0x1.eade43f264d12p+0 * lgt + 0x1.389b4b447b90fp+23 * rt_inv);
  k[r19b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM9]);
  k[r20f] = 0x1.671e344p+34 * exp(-0x1.fae9600000002p+23 * rt_inv);
  k[r20b] = 0x1.5ddc9971ac6b5p+17 *
            exp(0x1.433272bb73b8p+0 * lgt - 0x1.1c7a979aa68d3p+28 * rt_inv);
  k[r21f] = 0x1.48106p+24 * exp(0x1p+0 * lgt - 0x1.7f0e800000001p+24 * rt_inv);
  k[r21b] = 0x1.22f867c67033ap+16 *
            exp(0x1.9adc4ffefe21p+0 * lgt - 0x1.5fd302ec32592p+26 * rt_inv);
  k[r22f] = 0x1.2a6ffffffffffp+13 *
            exp(0x1p+1 * lgt - 0x1.fae9600000002p+23 * rt_inv);
  k[r22b] = 0x1.36cf3f9fa7108p+4 *
            exp(0x1.4929ca189a39cp+1 * lgt - 0x1.291a8188cb0fcp+26 * rt_inv);
  k[r23f] = 0x1.9ed92cp+30 * exp(-0x1.44d5000000001p+20 * rt_inv);
  k[r23b] = 0x1.f1a0f445fb656p+25 *
            exp(0x1.01af0a3623ecp-1 * lgt - 0x1.f49ca355f2338p+26 * rt_inv);
  k[r24f] = 0x1.1abfe17p+36 * exp(-0x1.d012b80000001p+24 * rt_inv);
  k[r24b] = 0x1.532ba2139441dp+31 *
            exp(0x1.01af0a3624bap-1 * lgt - 0x1.31c6feaaf91bcp+27 * rt_inv);
  kTroe0 = 0x1.047554901888p+60 *
           exp(-0x1.651eb851eb852p+1 * lgt - 0x1.0b90a80000001p+24 * rt_inv);
  kTroeInf = 0x1.9fa64p+23 * exp(-0x1.3067p+23 * rt_inv);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r25f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM10]);
  kTroe0 = 0x1.6879f052d1b9p+89 *
           exp(-0x1.dfe4427990bap+1 * lgt - 0x1.0751d6e9795a2p+29 * rt_inv);
  kTroeInf = 0x1.1fa188e92b4a7p+53 *
             exp(-0x1.eb16289e94d38p-1 * lgt - 0x1.03b6eda9795a2p+29 * rt_inv);
  k[r25b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM10]);
  k[r26f] = 0x1.0aca57p+30 * exp(-0x1.7ca99c0000001p+27 * rt_inv);
  k[r26b] = 0x1.625141dd64adcp+43 *
            exp(-0x1.ae9225aafbf2p-1 * lgt - 0x1.c3002044af1a1p+27 * rt_inv);
  k[r27f] = 0x1.1899999999999p+6 *
            exp(0x1.06c8b43958106p+1 * lgt + 0x1.6b578cccccccdp+20 * rt_inv);
  k[r27b] = 0x1.608941a8b2432p+28 *
            exp(0x1.8f9f5fbb96b8p-1 * lgt - 0x1.94172f1295013p+26 * rt_inv);
  k[r28f] = 0x1.5724d94p+32 *
            exp(-0x1.53f7ced916873p-1 * lgt - 0x1.52edb33333333p+20 * rt_inv);
  k[r28b] = 0x1.af1d3c9f6c051p+54 *
            exp(-0x1.efbda0016ff78p+0 * lgt - 0x1.9f1044129500bp+26 * rt_inv);
  k[r29f] =
      0x1.3ap+7 * exp(0x1.170a3d70a3d71p+1 * lgt - 0x1.1e55b8p+26 * rt_inv);
  k[r29b] = 0x1.4065b3ebf9299p+17 *
            exp(0x1.aaadafecc594p+0 * lgt - 0x1.3ec78a169c0f2p+28 * rt_inv);
  k[r30f] = 0x1.0fcc14p+29 *
            exp(0x1.51eb851eb851fp-1 * lgt - 0x1.daabc8p+25 * rt_inv);
  k[r30b] = 0x1.b22591d7b4ad8p+15 *
            exp(0x1.f693c540402dp-1 * lgt + 0x1.c35a7097e27a7p+21 * rt_inv);
  k[r31f] = 0x1.c3cd9fp+32 * exp(-0x1.a2cfp+20 * rt_inv);
  k[r31b] = 0x1.02aea62171da8p+30 *
            exp(0x1.3aa1f1d57fc4p-2 * lgt - 0x1.0c728e5818af1p+27 * rt_inv);
  k[r32f] = 0x1.116fb1ep+36;
  k[r32b] = 0x1.99c0216dfaea9p+30 *
            exp(0x1.562700484d88p-1 * lgt - 0x1.5e31d5792ab1cp+28 * rt_inv);
  k[r33f] = 0x1.c203db8p+34;
  k[r33b] = 0x1.8bec6490825fep+28 *
            exp(0x1.451588acba4cp-1 * lgt - 0x1.589d522050dfep+28 * rt_inv);
  k[r34f] = 0x1.bf08ebp+34;
  k[r34b] = 0x1.ee215be7c88f1p+50 *
            exp(-0x1.466de87d0f2cp-1 * lgt - 0x1.bf0e7571c2ed3p+28 * rt_inv);
  k[r35f] = 0x1.7bfac7cp+36;
  k[r35b] = 0x1.810d58c202794p+34 *
            exp(0x1.221d6a8075fp-1 * lgt - 0x1.9a1050939aa8ep+28 * rt_inv);
  k[r36] = 0x1.bf08ebp+34;
  k[r37] = 0x1.65a0bcp+31;
  k[r38f] = 0x1.8ae17a4p+34;
  k[r38b] = 0x1.055dc1b3d66b1p+40 *
            exp(0x1.3685dcde18fp-3 * lgt - 0x1.66ac5c7232a9p+28 * rt_inv);
  k[r39f] = 0x1.d5bc5eefp+42 * exp(-0x1.aa4f640000001p+27 * rt_inv);
  k[r39b] = 0x1.4a22b92668f39p+37 *
            exp(-0x1.de4fe3166ep-7 * lgt - 0x1.be0309c13e81p+25 * rt_inv);
  k[r40f] = 0x1.c9c38p+26 * exp(0x1.18e8800000001p+22 * rt_inv);
  k[r40b] = 0x1.b8f5192cd208ep+53 *
            exp(-0x1.2fd49aa2012fp+0 * lgt - 0x1.464540400c1cap+27 * rt_inv);
  k[r41f] = 0x1.da73f6p+30 * exp(-0x1.74341p+25 * rt_inv);
  k[r41b] = 0x1.6a876b3cd7529p+39 *
            exp(-0x1.fe9ce05dad0c8p-1 * lgt - 0x1.7ed0bad9b865p+25 * rt_inv);
  k[r42f] = 0x1.2a05f2p+34;
  k[r42b] = 0x1.43a680e2ceb38p+71 *
            exp(-0x1.024eaa2f85118p+1 * lgt - 0x1.73b3ae1cf95a4p+27 * rt_inv);
  k[r43f] = 0x1.176592ep+36 * exp(-0x1.cedc300000001p+26 * rt_inv);
  k[r43b] = 0x1.4aa0fa85e8ce9p+43 *
            exp(0x1.493287045abp-5 * lgt - 0x1.1f21bbfcb19bp+26 * rt_inv);
  k[r44f] = 0x1.0c388dp+34;
  k[r44b] = 0x1.b50514b9dde3bp+36 *
            exp(0x1.499470ee3e1p-2 * lgt - 0x1.21fdc8f3e4914p+28 * rt_inv);
  k[r45f] =
      0x1.61749e8p+32 * exp(0x1p-1 * lgt + 0x1.b73d000000001p+21 * rt_inv);
  k[r45b] = 0x1.2b1ab81701da2p+32 *
            exp(0x1.0ca7ba6ea978p-1 * lgt - 0x1.6e3b7c8f92db3p+28 * rt_inv);
  k[r46f] = 0x1.5faadbp+31 * exp(0x1p-1 * lgt + 0x1.c58ap+20 * rt_inv);
  k[r46b] = 0x1.2997589d94771p+31 *
            exp(0x1.0ca7ba6ea9bcp-1 * lgt - 0x1.6fe46c8f92dadp+28 * rt_inv);
  k[r47f] = 0x1.9bfccp+26 * exp(0x1p-1 * lgt + 0x1.3d2bcp+22 * rt_inv);
  k[r47b] = 0x1.5ca2ba151d611p+26 *
            exp(0x1.0ca7ba6ea99p-1 * lgt - 0x1.6cb5478f92dabp+28 * rt_inv);
  k[r48f] =
      0x1.6639528p+32 * exp(0x1p-1 * lgt + 0x1.8635000000001p+21 * rt_inv);
  k[r48b] = 0x1.2f23b6b0260ap+32 *
            exp(0x1.0ca7ba6ea91p-1 * lgt - 0x1.6e9d8c8f92dc8p+28 * rt_inv);
  k[r49f] = 0x1.38540ep+30 * exp(0x1p-1 * lgt + 0x1.552d000000001p+19 * rt_inv);
  k[r49b] = 0x1.084d242e3bf73p+30 *
            exp(0x1.0ca7ba6ea988p-1 * lgt - 0x1.70ff600f92db8p+28 * rt_inv);
  k[r50f] = 0x1.92ed6ap+30 * exp(-0x1.07fd68p+24 * rt_inv);
  k[r50b] = 0x1.54f8098bc8009p+30 *
            exp(0x1.94f74dd53fp-6 * lgt - 0x1.8229cd0f92dc6p+28 * rt_inv);
  k[r51f] = 0x1.6201p+20;
  k[r51b] = 0x1.2b918372dbaffp+20 *
            exp(0x1.94f74dd535p-6 * lgt - 0x1.71a9f68f92dbp+28 * rt_inv);
  k[r52f] = 0x1.f4add4p+30 * exp(0x1p-1 * lgt + 0x1.e845000000001p+20 * rt_inv);
  k[r52b] = 0x1.a7b06ccefe827p+30 *
            exp(0x1.0ca7ba6ea9cp-1 * lgt - 0x1.6fc1b18f92da4p+28 * rt_inv);
  k[r53f] = 0x1.47d357p+31 * exp(0x1p-1 * lgt + 0x1.ee66000000001p+21 * rt_inv);
  k[r53b] = 0x1.156a5f9feb9ffp+31 *
            exp(0x1.0ca7ba6ea928p-1 * lgt - 0x1.6dcd2a8f92dc1p+28 * rt_inv);
  k[r54f] = 0x1.810bc7p+31 * exp(0x1p-1 * lgt + 0x1.91f44p+21 * rt_inv);
  k[r54b] = 0x1.45d64ecd83c6dp+31 *
            exp(0x1.0ca7ba6ea98cp-1 * lgt - 0x1.6e860e0f92dafp+28 * rt_inv);
  k[r55f] = 0x1.908b1p+31 * exp(0x1p-1 * lgt + 0x1.44524p+21 * rt_inv);
  k[r55b] = 0x1.52f38a3f335b4p+31 *
            exp(0x1.0ca7ba6ea974p-1 * lgt - 0x1.6f21520f92db9p+28 * rt_inv);
  k[r56f] = 0x1.2d00e28p+35;
  k[r56b] = 0x1.8a3b65a5d394ep+32 *
            exp(0x1.6e3881371d58p-2 * lgt - 0x1.0ab1c0eb92c43p+28 * rt_inv);
  kTroe0 = 0x1.2bc29d8eec7p+60 *
           exp(-0x1.48f5c28f5c28fp+1 * lgt - 0x1.6be76p+22 * rt_inv);
  kTroeInf = 0x1.03e052p+30 *
             exp(0x1.eb851eb851eb8p-2 * lgt + 0x1.0996000000001p+20 * rt_inv);
  fcTroe = 0x1.bda5119ce076p-3 * exp(-temp / 0x1.0fp+8) +
           0x1.9096bb98c7e28p-1 * exp(-temp / 0x1.586p+11) +
           0x1p+0 * exp(-0x1.9aap+12 / temp);
  k[r57f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM12]);
  kTroe0 = 0x1.31c46874bf09fp+76 *
           exp(-0x1.48ed4479f5f7p+1 * lgt - 0x1.6689945514598p+28 * rt_inv);
  kTroeInf = 0x1.0915814fd653bp+46 *
             exp(0x1.ebc90f63837bp-2 * lgt - 0x1.5fd060d514598p+28 * rt_inv);
  k[r57b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM12]);
  kTroe0 = 0x1.12d86259fd1e8p+72 *
           exp(-0x1.b5c28f5c28f5cp+1 * lgt - 0x1.50902ap+28 * rt_inv);
  kTroeInf = 0x1.4ffp+15 * exp(0x1.8p+0 * lgt - 0x1.3d9e28p+28 * rt_inv);
  fcTroe = 0x1.16872b020c498p-4 * exp(-temp / 0x1.8ap+7) +
           0x1.dd2f1a9fbe76dp-1 * exp(-temp / 0x1.81p+10) +
           0x1p+0 * exp(-0x1.41ep+13 / temp);
  k[r58f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM13]);
  kTroe0 = 0x1.762c278483cf9p+93 *
           exp(-0x1.05a1e8ac69fa8p+2 * lgt - 0x1.53384b5be9b61p+28 * rt_inv);
  kTroeInf = 0x1.c957d2806e4a6p+36 *
             exp(0x1.a9faf80d5403p-1 * lgt - 0x1.4046495be9b61p+28 * rt_inv);
  k[r58b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM13]);
  k[r59f] = 0x1.3178p+16 *
            exp(0x1.a147ae147ae14p+0 * lgt + 0x1.0d6aa00000001p+22 * rt_inv);
  k[r59b] = 0x1.7bf7e93fdf104p+11 *
            exp(0x1.dff147192625p+0 * lgt - 0x1.cf73747ed861dp+26 * rt_inv);
  k[r60f] = 0x1.c07p+15 *
            exp(0x1.e666666666666p+0 * lgt - 0x1.5ddbc00000001p+23 * rt_inv);
  k[r60b] = 0x1.9c6fa615e739ep+7 *
            exp(0x1.1f8a65277edbcp+1 * lgt - 0x1.1c8baa1518865p+26 * rt_inv);
  k[r61f] =
      0x1.7e148p+22 * exp(0x1.2666666666666p+0 * lgt - 0x1.2091cp+23 * rt_inv);
  k[r61b] = 0x1.9ca114167f09ap+13 *
            exp(0x1.768c0e81341bp+0 * lgt - 0x1.fd20b963627ccp+25 * rt_inv);
  k[r62f] = 0x1.39c0ebedfa43fp-5 *
            exp(0x1.ae147ae147ae1p+1 * lgt - 0x1.134a4p+24 * rt_inv);
  k[r62b] = 0x1.fdd3f64faf5ffp-3 *
            exp(0x1.984a4a91618ep+1 * lgt - 0x1.578cddf938641p+26 * rt_inv);
  k[r63f] = 0x1.2cccccccccccdp+4 *
            exp(0x1.599999999999ap+1 * lgt - 0x1.6fbcp+25 * rt_inv);
  k[r63b] = 0x1.37eb5fb92477cp+4 *
            exp(0x1.3882a38e6689cp+1 * lgt - 0x1.513d9e51cc74p+25 * rt_inv);
  k[r64f] = 0x1.05ef39b2p+42 * exp(-0x1.199999999999ap+0 * lgt);
  k[r64b] = 0x1.1dbccd059d803p+71 *
            exp(-0x1.2e196f1a25e3cp+1 * lgt - 0x1.9c04f91249508p+26 * rt_inv);
  k[r65f] = 0x1.6bcc41e9p+46 * exp(-0x1.dba0f00000001p+25 * rt_inv);
  k[r65b] = 0x1.016c2370c2364p+29 *
            exp(0x1.a52ba59e3e6a8p-1 * lgt - 0x1.562b87dc65f87p+25 * rt_inv);
  k[r66f] = 0x1.1d37b09ap+41 * exp(-0x1.e29252p+27 * rt_inv);
  k[r66b] = 0x1.6f13d6c1b4588p+11 *
            exp(0x1.01012ae10557p+0 * lgt - 0x1.a859e2fcce96dp+27 * rt_inv);
  k[r67f] = 0x1.5faadbp+31 * exp(-0x1.8334d80000001p+27 * rt_inv);
  k[r67b] = 0x1.a4709113c814cp+20 *
            exp(-0x1.55de167fa07p-3 * lgt - 0x1.9e21b96ad2bbdp+27 * rt_inv);
  k[r68f] = 0x1.0b076p+25 * exp(0x1.a22b9p+23 * rt_inv);
  k[r68b] = 0x1.edfc00e1fd5bap+25 *
            exp(0x1.f374c91a6fa8p-2 * lgt - 0x1.d6045a3111457p+27 * rt_inv);
  k[r69] = 0x1.477ffffffffffp+11 *
           exp(0x1.07ae147ae147bp+1 * lgt - 0x1.d3d7000000001p+21 * rt_inv);
  k[r70] = 0x1.211p+14 *
           exp(0x1.828f5c28f5c29p+0 * lgt + 0x1.eb55800000002p+21 * rt_inv);
  k[r71] =
      0x1.09p+12 * exp(0x1.0cccccccccccdp+1 * lgt - 0x1.36c96p+24 * rt_inv);
  k[r72] = 0x1.c1451f6p+35 *
           exp(-0x1.6666666666666p-2 * lgt - 0x1.7d864p+23 * rt_inv);
  k[r73] = 0x1.accf3dacbf296p-32 *
           exp(0x1.7333333333333p+2 * lgt - 0x1.18e8800000001p+23 * rt_inv);
  k[r74f] = 0x1.1e1a3p+31 * exp(-0x1.3f36cp+25 * rt_inv);
  k[r74b] = 0x1.66b7428da7c2bp+26 *
            exp(0x1.769cfd2da738p-1 * lgt - 0x1.0cc5ba0cedbf6p+27 * rt_inv);
  k[r75] = 0x1.dcd65p+29 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r76] = 0x1.92738f5028p+50 *
           exp(-0x1.e666666666666p+0 * lgt - 0x1.7bdd5p+23 * rt_inv);
  k[r77f] = 0x1.4dc938p+32 * exp(-0x1.b2218p+25 * rt_inv);
  k[r77b] = 0x1.b1f8e5991003ap+27 *
            exp(0x1.e4824a01b47cp-2 * lgt - 0x1.21e0d1a160d8ep+27 * rt_inv);
  kTroe0 = 0x1.fa0d26310cdf8p+73 *
           exp(-0x1.8p+1 * lgt - 0x1.83f4e20000001p+26 * rt_inv);
  kTroeInf = 0x1.eec3dec2p+45 * exp(-0x1.a1b0fc0000001p+26 * rt_inv);
  fcTroe = 0x1.9999999999998p-4 * exp(-temp / 0x1.388p+11) +
           0x1.ccccccccccccdp-1 * exp(-temp / 0x1.45p+10) +
           0x1p+0 * exp(-0x1.d42aea2879f2ep+328 / temp);
  k[r78f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM14]);
  kTroe0 = 0x1.2944901732db5p+53 *
           exp(-0x1.0e4441ce88ed8p+1 * lgt - 0x1.6599d0b401182p+23 * rt_inv);
  kTroeInf = 0x1.22a3531e6833bp+25 *
             exp(0x1.c6eef8c5dc4ap-1 * lgt - 0x1.29bd505a008c3p+24 * rt_inv);
  k[r78b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM14]);
  k[r79f] =
      0x1.08c1379733209p-71 * exp(0x1.3p+3 * lgt + 0x1.5f32f8p+24 * rt_inv);
  k[r79b] = 0x1.bdfa2df1007b9p-82 *
            exp(0x1.4bf97b18f128fp+3 * lgt - 0x1.5ebc324b708d4p+26 * rt_inv);
  k[r80f] = 0x1.3baa8cp+29 * exp(-0x1.24e92p+23 * rt_inv);
  k[r80b] = 0x1.3f24e046c7cccp+20 *
            exp(0x1.5c20631bd221p-1 * lgt - 0x1.283d8e06971b4p+26 * rt_inv);
  k[r81f] = 0x1.d7dbf487fcb92p-7 *
            exp(0x1.8cccccccccccdp+1 * lgt - 0x1.babfe8p+24 * rt_inv);
  k[r81b] = 0x1.7b310ccf9243p+5 *
            exp(0x1.1ffa83b5f23fap+1 * lgt - 0x1.f72777ca85367p+24 * rt_inv);
  k[r82f] = 0x1.65a0bcp+33;
  k[r82b] = 0x1.5c34cb506d26fp+31 *
            exp(0x1.6fe82fdbd9b4p-1 * lgt - 0x1.4fb820590274dp+28 * rt_inv);
  k[r83f] = 0x1.2a05f2p+34;
  k[r83b] = 0x1.487a241715b09p+21 *
            exp(0x1.3c36dc76508ep+0 * lgt - 0x1.473d995ffa7ddp+28 * rt_inv);
  k[r84f] = 0x1.1f0e54p+28;
  k[r84b] = 0x1.64b92667b705dp+23 *
            exp(0x1.42b518eea55p-1 * lgt - 0x1.0739c0a4ede8ep+28 * rt_inv);
  kTroe0 = 0x1.a4352f2e5ea79p+86 *
           exp(-0x1.347ae147ae148p+2 * lgt - 0x1.a0e4b00000001p+24 * rt_inv);
  kTroeInf =
      0x1.017df8p+29 * exp(0x1.d0e5604189375p-2 * lgt - 0x1.cbabp+23 * rt_inv);
  fcTroe = 0x1.200d1b71758e2p-2 * exp(-temp / 0x1.9cp+6) +
           0x1.6ff972474538fp-1 * exp(-temp / 0x1.42cp+10) +
           0x1p+0 * exp(-0x1.04p+12 / temp);
  k[r85f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM15]);
  kTroe0 = 0x1.9235e895669cap+99 *
           exp(-0x1.4afd62a404f62p+2 * lgt - 0x1.215c4ecd2b1ecp+27 * rt_inv);
  kTroeInf = 0x1.eced81c3881dep+41 *
             exp(0x1.a2f529f06c754p-4 * lgt - 0x1.09fa68cd2b1ecp+27 * rt_inv);
  k[r85b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM15]);
  k[r86f] = 0x1.5f93037cp+40 * exp(-0x1p+0 * lgt);
  k[r86b] = 0x1.0751d633ec1bdp+38 *
            exp(-0x1.53433c540ee2p-1 * lgt - 0x1.334b269a9a45cp+26 * rt_inv);
  k[r87f] = 0x1.c0e5c15p+37 * exp(-0x1.404c98p+24 * rt_inv);
  k[r87b] = 0x1.50364edef6c77p+35 *
            exp(0x1.59798757e31ap-2 * lgt - 0x1.835e4c9a9a4acp+26 * rt_inv);
  k[r88f] = 0x1.65a0bcp+32;
  k[r88b] = 0x1.5e854a7c5def4p+27 *
            exp(0x1.6592cb097edcp-1 * lgt - 0x1.266e26f3c4eb7p+28 * rt_inv);
  k[r89f] = 0x1.65a0bcp+33;
  k[r89b] = 0x1.8b34b78eaf223p+36 *
            exp(0x1.7ed1585c198p-4 * lgt - 0x1.ccd49c7170acfp+27 * rt_inv);
  k[r90f] = 0x1.4f46b04p+37;
  k[r90b] = 0x1.654c311b8dafp+40 *
            exp(0x1.686c0670a088p-2 * lgt - 0x1.d47434dcfd968p+27 * rt_inv);
  k[r91f] = 0x1.65a0bcp+34;
  k[r91b] = 0x1.81510117dea57p+28 *
            exp(0x1.082b332a1122p+0 * lgt - 0x1.2b2234f024912p+28 * rt_inv);
  k[r92f] = 0x1.3428f5c28f5c2p+3 *
            exp(0x1.7333333333333p+1 * lgt - 0x1.a27d480000001p+25 * rt_inv);
  k[r92b] = 0x1.56c2996b11cfap+4 *
            exp(0x1.60e1d5e5078e6p+1 * lgt - 0x1.43700ea8fcf24p+24 * rt_inv);
  k[r93f] = 0x1.65a0bcp+34;
  k[r93b] = 0x1.da111a3d1ba9ap+32 *
            exp(0x1.31893541a734p-1 * lgt - 0x1.624ca20e34e2cp+28 * rt_inv);
  k[r94f] = 0x1.38eca48p+35;
  k[r94b] = 0x1.682414095f80fp+29 *
            exp(0x1.5481536deb9cp-1 * lgt - 0x1.20d9a39aeb19ep+28 * rt_inv);
  k[r95f] = 0x1.65a0bcp+31;
  k[r95b] = 0x1.56a532e52a9cep+33 *
            exp(0x1.faf6f0e2005p-2 * lgt - 0x1.0a52c283eefe9p+28 * rt_inv);
  k[r96f] = 0x1.2a05f2p+33;
  k[r96b] = 0x1.13da7c11f18ecp+37 *
            exp(-0x1.31a28c47a63p-3 * lgt - 0x1.0244b5464d8bbp+27 * rt_inv);
  k[r97f] = 0x1.1e1a3p+27 * exp(-0x1.7bdd5p+25 * rt_inv);
  k[r97b] = 0x1.662406b7a5be1p+73 *
            exp(-0x1.7546d0f260968p+1 * lgt - 0x1.c595f06b572ccp+25 * rt_inv);
  k[r98f] = 0x1.176592ep+38 * exp(-0x1.12862p+25 * rt_inv);
  k[r98b] = 0x1.3655379b3bd43p+27 *
            exp(0x1.13c7459dbb918p+0 * lgt - 0x1.9f7d754927f53p+26 * rt_inv);
  k[r99f] = 0x1.0b076p+25 * exp(0x1.a22b9p+23 * rt_inv);
  k[r99b] = 0x1.869012fb2f329p+37 *
            exp(-0x1.da47d14de17ap-1 * lgt - 0x1.108513e26608bp+27 * rt_inv);
  k[r100f] = 0x1.2a05f2p+33;
  k[r100b] = 0x1.ea6905d698fefp+74 *
             exp(-0x1.23abc8c784f18p+1 * lgt - 0x1.5f9405a8cf5f5p+27 * rt_inv);
  kTroe0 = 0x1.60a5f7552857cp+133 *
           exp(-0x1.bfae147ae147bp+2 * lgt - 0x1.8701804cccccdp+28 * rt_inv);
  kTroeInf = 0x1.cebdaf55feap+60 *
             exp(-0x1.3ae147ae147aep-1 * lgt - 0x1.7140c6e666667p+28 * rt_inv);
  fcTroe = 0x1.798c7e28240b8p+0 * exp(-temp / 0x1.15f8p+15) +
           -0x1.e631f8a0902dep-2 * exp(-temp / 0x1.17p+10) +
           0x1p+0 * exp(-0x1.19f8p+13 / temp);
  k[r101f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM16]);
  kTroe0 = 0x1.c2879a22a7f5ep+105 *
           exp(-0x1.8d1ae6b6bc7f4p+2 * lgt - 0x1.3133c3d052587p+24 * rt_inv);
  kTroeInf = 0x1.2797028dd4008p+33 *
             exp(0x1.66e099cc47228p-3 * lgt + 0x1.56be94b0a0742p+21 * rt_inv);
  k[r101b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM16]);
  kTroe0 = 0x1.9a63cebf2f32bp+146 *
           exp(-0x1.074395810624ep+3 * lgt - 0x1.8cb103a666667p+28 * rt_inv);
  kTroeInf = 0x1.5a8027b6bff4p+61 *
             exp(-0x1.045a1cac08312p+0 * lgt - 0x1.6df26p+28 * rt_inv);
  fcTroe = -0x1.8b851eb851eb8p+0 * exp(-temp / 0x1.9b4p+11) +
           0x1.45c28f5c28f5cp+1 * exp(-temp / 0x1.71bp+15) +
           0x1p+0 * exp(-0x1.700cp+15 / temp);
  k[r102f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM17]);
  kTroe0 = 0x1.2b486a68ea1f5p+114 *
           exp(-0x1.bdbdab314410cp+2 * lgt - 0x1.462f788f0343ap+24 * rt_inv);
  kTroeInf = 0x1.f961a11906a65p+28 *
             exp(0x1.f65f14b8c597p-3 * lgt + 0x1.4b7583aec646ep+23 * rt_inv);
  k[r102b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM17]);
  kTroe0 = 0x1.3ecb480966017p+131 *
           exp(-0x1.cf9db22d0e56p+2 * lgt - 0x1.a3e31b7333334p+28 * rt_inv);
  kTroeInf = 0x1.02bc72e275ab8p-7 *
             exp(0x1.426e978d4fdf4p+2 * lgt - 0x1.510a21199999ap+28 * rt_inv);
  fcTroe = 0x1.2ba3d70a3d70ap+6 * exp(-temp / 0x1.2174p+15) +
           -0x1.27a3d70a3d70ap+6 * exp(-temp / 0x1.4438p+15) +
           0x1p+0 * exp(-0x1.464p+12 / temp);
  k[r103f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM18]);
  kTroe0 = 0x1.5b9f414039d28p+116 *
           exp(-0x1.d8ca9fded6082p+2 * lgt - 0x1.1783e44575c95p+25 * rt_inv);
  kTroeInf = 0x1.1a2219fd320c8p-22 *
             exp(0x1.3941a9db882d2p+2 * lgt + 0x1.7f43ee875703dp+25 * rt_inv);
  k[r103b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM18]);
  k[r104f] =
      0x1.33p+8 * exp(0x1.4666666666666p+1 * lgt - 0x1.5b4ep+24 * rt_inv);
  k[r104b] = 0x1.3a0e5bcba5a13p+1 *
             exp(0x1.606c3b0c862fap+1 * lgt - 0x1.8e82237eaf713p+25 * rt_inv);
  k[r105f] = 0x1.8dfffffffffffp+7 *
             exp(0x1.47ae147ae147bp+1 * lgt - 0x1.48ca500000001p+25 * rt_inv);
  k[r105b] = 0x1.6a0f119cefb66p+8 *
             exp(0x1.1cfd2da838a4cp+1 * lgt - 0x1.2329e01d02e66p+25 * rt_inv);
  k[r106f] = 0x1.83fffffffffffp+8 *
             exp(0x1.4p+1 * lgt - 0x1.8945800000001p+23 * rt_inv);
  k[r106b] = 0x1.d211fac86a985p+0 *
             exp(0x1.55c176bf3aefap+1 * lgt - 0x1.168868b7e0de5p+25 * rt_inv);
  k[r107f] = 0x1.eccccccccccccp+4 *
             exp(0x1.5333333333333p+1 * lgt + 0x1.9c04066666668p+21 * rt_inv);
  k[r107b] = 0x1.54e819a999a1bp+1 *
             exp(0x1.6036a2675ce08p+1 * lgt - 0x1.53075df5e4606p+26 * rt_inv);
  k[r108f] = 0x1.3333333333333p-3 *
             exp(0x1.83d70a3d70a3dp+1 * lgt + 0x1.85b2400000001p+21 * rt_inv);
  k[r108b] = 0x1.79f52e35ce4ffp+1 *
             exp(0x1.4c23bdf8d1c34p+1 * lgt - 0x1.a0f844f082965p+25 * rt_inv);
  k[r109f] = 0x1.31794b4p+34 * exp(-0x1.66514cp+27 * rt_inv);
  k[r109b] = 0x1.dd9778cc48dc2p+29 *
             exp(-0x1.42fad31429e7p-3 * lgt + 0x1.4aea5ba90c9f1p+23 * rt_inv);
  k[r110f] = 0x1.5999999999999p+3 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.5021d80000001p+25 * rt_inv);
  k[r110b] = 0x1.8e9d25c6f6de5p+4 *
             exp(0x1.12fe130d0765cp+1 * lgt - 0x1.86f1ad325685ep+22 * rt_inv);
  k[r111f] = 0x1.05532617c1bdap-5 *
             exp(0x1.95c28f5c28f5cp+1 * lgt - 0x1.c9e1600000001p+24 * rt_inv);
  k[r111b] = 0x1.d8500e04eb2bep-2 *
             exp(0x1.6da701be16f64p+1 * lgt - 0x1.04d005a37796bp+26 * rt_inv);
  k[r112f] = 0x1.1e1a3p+28 * exp(-0x1.04187p+24 * rt_inv);
  k[r112b] = 0x1.41bb5b4d7fa6fp+20 *
             exp(0x1.12daede322838p-1 * lgt - 0x1.8887cb61ac9acp+25 * rt_inv);
  k[r113f] = 0x1.c1451f6p+35;
  k[r113b] = 0x1.102ae3fb06475p+22 *
             exp(0x1.9198aa1ba25p+0 * lgt - 0x1.4bf1a75c5a22fp+28 * rt_inv);
  kTroe0 = 0x1.001da4e1fe9aep+91 *
           exp(-0x1.30a3d70a3d70ap+2 * lgt - 0x1.378d800000001p+23 * rt_inv);
  kTroeInf = 0x1.719e5fa3p+43 *
             exp(-0x1.428f5c28f5c29p-1 * lgt - 0x1.873a800000001p+20 * rt_inv);
  fcTroe = 0x1.bc6a7ef9db22cp-3 * exp(-temp / 0x1.28p+6) +
           0x1.90e5604189375p-1 * exp(-temp / 0x1.6fap+11) +
           0x1p+0 * exp(-0x1.b34p+12 / temp);
  k[r114f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM19]);
  kTroe0 = 0x1.a882a3a7d5044p+109 *
           exp(-0x1.3b84b0277d718p+2 * lgt - 0x1.af44f6536271ap+28 * rt_inv);
  kTroeInf = 0x1.3251fdf9a212dp+62 *
             exp(-0x1.99962512f5c99p-1 * lgt - 0x1.a70fc4d36271ap+28 * rt_inv);
  k[r114b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM19]);
  k[r115f] = 0x1.33p+9 * exp(0x1.4p+1 * lgt - 0x1.3207c4p+25 * rt_inv);
  k[r115b] = 0x1.5b86d7341dafcp-2 *
             exp(0x1.82216244320f8p+1 * lgt - 0x1.dc67186f80acfp+24 * rt_inv);
  k[r116f] = 0x1.d266666666666p+5 *
             exp(0x1.4cccccccccccdp+1 * lgt - 0x1.17a1ap+23 * rt_inv);
  k[r116b] = 0x1.65079782467p-2 *
             exp(0x1.81ebc99f08b38p+1 * lgt - 0x1.e108090b400c7p+25 * rt_inv);
  k[r117f] = 0x1.f20cp+19 * exp(0x1.8p+0 * lgt - 0x1.12862p+25 * rt_inv);
  k[r117b] = 0x1.4b02822f115e7p+8 *
             exp(0x1.fbba08ba9a7bcp+0 * lgt - 0x1.441b9ae1e38bcp+24 * rt_inv);
  k[r118f] = 0x1.15b573eab367ap-6 *
             exp(0x1.deb851eb851ecp+1 * lgt - 0x1.4f558cp+26 * rt_inv);
  k[r118b] = 0x1.627202a79c47fp-9 *
             exp(0x1.d36b8c30387ebp+1 * lgt - 0x1.6ae0697d6f172p+23 * rt_inv);
  k[r119f] = 0x1.338p+11 * exp(0x1p+1 * lgt - 0x1.07fd68p+25 * rt_inv);
  k[r119b] = 0x1.6576defc44bf8p+9 *
             exp(0x1.f93d2e49ecc1cp+0 * lgt - 0x1.aa63c20e16477p+25 * rt_inv);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 4.936e+11;
  PlogA[1] = 1.207e+12;
  PlogA[2] = 5.282e+14;
  PlogA[3] = 4.788e+20;
  PlogA[4] = 8.433e+16;
  PlogB[0] = -0.669;
  PlogB[1] = -0.778;
  PlogB[2] = -1.518;
  PlogB[3] = -3.155;
  PlogB[4] = -1.962;
  PlogE[0] = -1.86523e+06;
  PlogE[1] = -734710;
  PlogE[2] = 7.41405e+06;
  PlogE[3] = 2.93006e+07;
  PlogE[4] = 3.44929e+07;
  np = 5;
  k[r120f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.761e+10;
  PlogA[1] = 4.30617e+10;
  PlogA[2] = 1.88444e+13;
  PlogA[3] = 1.7082e+19;
  PlogA[4] = 3.00861e+15;
  PlogB[0] = -0.196935;
  PlogB[1] = -0.305935;
  PlogB[2] = -1.04593;
  PlogB[3] = -2.68293;
  PlogB[4] = -1.48993;
  PlogE[0] = -6.45185e+06;
  PlogE[1] = -5.32134e+06;
  PlogE[2] = 2.82742e+06;
  PlogE[3] = 2.47139e+07;
  PlogE[4] = 2.99063e+07;
  np = 5;
  k[r120b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 350.2;
  PlogA[1] = 885.4;
  PlogA[2] = 16500;
  PlogA[3] = 5.374e+06;
  PlogA[4] = 9.494e+15;
  PlogB[0] = 1.441;
  PlogB[1] = 1.327;
  PlogB[2] = 0.973;
  PlogB[3] = 0.287;
  PlogB[4] = -2.199;
  PlogE[0] = -1.35729e+07;
  PlogE[1] = -1.24474e+07;
  PlogE[2] = -8.40984e+06;
  PlogE[3] = 1.17152e+06;
  PlogE[4] = 4.08735e+07;
  np = 5;
  k[r121f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 74998.7;
  PlogA[1] = 189617;
  PlogA[2] = 3.53364e+06;
  PlogA[3] = 1.15089e+09;
  PlogA[4] = 2.03323e+18;
  PlogB[0] = 1.20579;
  PlogB[1] = 1.09179;
  PlogB[2] = 0.737785;
  PlogB[3] = 0.0517853;
  PlogB[4] = -2.43421;
  PlogE[0] = 2.8151e+08;
  PlogE[1] = 2.82636e+08;
  PlogE[2] = 2.86673e+08;
  PlogE[3] = 2.96255e+08;
  PlogE[4] = 3.35957e+08;
  np = 5;
  k[r121b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.621e+07;
  PlogA[1] = 1.807e+07;
  PlogA[2] = 4.686e+07;
  PlogA[3] = 1.525e+10;
  PlogA[4] = 3.59e+11;
  PlogB[0] = 0.965;
  PlogB[1] = 0.95;
  PlogB[2] = 0.833;
  PlogB[3] = 0.134;
  PlogB[4] = -0.186;
  PlogE[0] = 1.34474e+07;
  PlogE[1] = 1.35854e+07;
  PlogE[2] = 1.49201e+07;
  PlogE[3] = 2.36019e+07;
  PlogE[4] = 3.59866e+07;
  np = 5;
  k[r122f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.13341e+11;
  PlogA[1] = 1.26347e+11;
  PlogA[2] = 3.27648e+11;
  PlogA[3] = 1.06629e+14;
  PlogA[4] = 2.51015e+15;
  PlogB[0] = 0.0313997;
  PlogB[1] = 0.0163997;
  PlogB[2] = -0.1006;
  PlogB[3] = -0.7996;
  PlogB[4] = -1.1196;
  PlogE[0] = -201950;
  PlogE[1] = -63877.7;
  PlogE[2] = 1.27082e+06;
  PlogE[3] = 9.95262e+06;
  PlogE[4] = 2.23373e+07;
  np = 5;
  k[r122b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.186e+06;
  PlogA[1] = 1.188e+06;
  PlogA[2] = 1.23e+06;
  PlogA[3] = 1.798e+06;
  PlogA[4] = 5.242e+07;
  PlogB[0] = 1.016;
  PlogB[1] = 1.016;
  PlogB[2] = 1.011;
  PlogB[3] = 0.965;
  PlogB[4] = 0.551;
  PlogE[0] = 4.9957e+07;
  PlogE[1] = 4.9957e+07;
  PlogE[2] = 4.99988e+07;
  PlogE[3] = 5.0459e+07;
  PlogE[4] = 5.46849e+07;
  np = 5;
  k[r123f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.88781e+12;
  PlogA[1] = 1.89099e+12;
  PlogA[2] = 1.95784e+12;
  PlogA[3] = 2.86195e+12;
  PlogA[4] = 8.34391e+13;
  PlogB[0] = -0.454427;
  PlogB[1] = -0.454427;
  PlogB[2] = -0.459427;
  PlogB[3] = -0.505427;
  PlogB[4] = -0.919427;
  PlogE[0] = 1.9035e+06;
  PlogE[1] = 1.9035e+06;
  PlogE[2] = 1.94534e+06;
  PlogE[3] = 2.40558e+06;
  PlogE[4] = 6.63142e+06;
  np = 5;
  k[r123b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 867400;
  PlogA[1] = 3.115e+06;
  PlogA[2] = 1.557e+08;
  PlogA[3] = 1.704e+18;
  PlogA[4] = 7.25e+17;
  PlogB[0] = 0.787;
  PlogB[1] = 0.63;
  PlogB[2] = 0.156;
  PlogB[3] = -2.641;
  PlogB[4] = -2.402;
  PlogE[0] = -1.27445e+07;
  PlogE[1] = -1.11671e+07;
  PlogE[2] = -5.72371e+06;
  PlogE[3] = 2.68278e+07;
  PlogE[4] = 4.03296e+07;
  np = 5;
  k[r124f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.00073e+15;
  PlogA[1] = 3.5938e+15;
  PlogA[2] = 1.79632e+17;
  PlogA[3] = 1.96592e+27;
  PlogA[4] = 8.36438e+26;
  PlogB[0] = -1.94061;
  PlogB[1] = -2.09761;
  PlogB[2] = -2.57161;
  PlogB[3] = -5.36861;
  PlogB[4] = -5.12961;
  PlogE[0] = 8.27503e+07;
  PlogE[1] = 8.43277e+07;
  PlogE[2] = 8.9771e+07;
  PlogE[3] = 1.22323e+08;
  PlogE[4] = 1.35824e+08;
  np = 5;
  k[r124b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r125f] = 0x1.2a05f2p+34;
  k[r125b] = 0x1.209fc768b4bdfp+7 *
             exp(0x1.5e5b965cdd1cp+1 * lgt - 0x1.366a2df07f0e9p+28 * rt_inv);
  k[r126f] = 0x1.74876e8p+37;
  k[r126b] = 0x1.220aef77aa188p+15 *
             exp(0x1.3f06c9da87698p+1 * lgt - 0x1.7caf4ca191ed4p+27 * rt_inv);
  k[r127] = 0x1.74876e8p+35;
  k[r128] = 0x1.bf08ebp+34;
  k[r129] = 0x1.2a05f2p+32;
  k[r130f] = 0x1.bf08ebp+34;
  k[r130b] = 0x1.f9dc29366f2bcp+14 *
             exp(0x1.43fc6f1d4b9ep+1 * lgt - 0x1.5919611993b88p+29 * rt_inv);
  k[r131f] =
      0x1.dcd65p+29 * exp(0x1.1374bc6a7ef9ep-2 * lgt + 0x1.5f22ap+21 * rt_inv);
  k[r131b] = 0x1.25f02e55a446p+38 *
             exp(-0x1.c4d72039088ap-2 * lgt - 0x1.80bb98f8de65bp+26 * rt_inv);
  k[r132f] = 0x1.cffffffffffffp+6 *
             exp(0x1.1d70a3d70a3d7p+1 * lgt + 0x1.81ddap+23 * rt_inv);
  k[r132b] = 0x1.0c327a7c0885ep+15 *
             exp(0x1.0984c36a3bbep+1 * lgt - 0x1.ac0dee8c4ca24p+27 * rt_inv);
  k[r133f] =
      0x1.9cc31d4p+35 * exp(0x1.999999999999ap-5 * lgt + 0x1.15d8p+19 * rt_inv);
  k[r133b] = 0x1.9575de5980971p+42 *
             exp(-0x1.bf97ba9f7ddp-3 * lgt - 0x1.134a5cc02e45dp+28 * rt_inv);
  k[r134f] = 0x1.c1c6d28p+32 * exp(-0x1.c401cp+26 * rt_inv);
  k[r134b] = 0x1.68ddb3fd0f47cp+44 *
             exp(-0x1.09f53f09dda2p+0 * lgt + 0x1.93e760ce9b7dfp+16 * rt_inv);
  k[r135f] = 0x1.5a294141e9af5p-9 *
             exp(0x1.a4395810624ddp+1 * lgt - 0x1.02b90c0000001p+25 * rt_inv);
  k[r135b] = 0x1.677215de942efp-11 *
             exp(0x1.b915c8dfb67p+1 * lgt - 0x1.e1cd1222277cbp+27 * rt_inv);
  kTroe0 = 0x1.7c405ad41db4p+62 * exp(-0x1.8p+1 * lgt);
  kTroeInf = 0x1.dcce8p+22 * exp(0x1.ccccccccccccdp-1 * lgt);
  fcTroe = 0x1.999999999999ap-2 * exp(-temp / 0x1.f4p+9) +
           0x1.3333333333333p-1 * exp(-temp / 0x1.18p+6) +
           0x1p+0 * exp(-0x1.a9p+10 / temp);
  k[r136f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM20]);
  kTroe0 = 0x1.afeda0386bdb3p+91 *
           exp(-0x1.09a049e25fb94p+2 * lgt - 0x1.08fed9cb0971ep+27 * rt_inv);
  kTroeInf = 0x1.0ecd8167fa91p+52 *
             exp(-0x1.006b048c61fa6p-2 * lgt - 0x1.08fed9cb0971ep+27 * rt_inv);
  k[r136b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM20]);
  k[r137f] = 0x1.da73f6p+30 * exp(-0x1.74341p+25 * rt_inv);
  k[r137b] = 0x1.7720b52e5bb17p+35 *
             exp(-0x1.a6e575210a8ep-1 * lgt - 0x1.0aae77b3351c9p+25 * rt_inv);
  k[r138f] = 0x1.f75104d551d68p-11 *
             exp(0x1.e28f5c28f5c29p+1 * lgt - 0x1.1c428c0000001p+26 * rt_inv);
  k[r138b] = 0x1.e9cdcc45b1f69p-9 *
             exp(0x1.8ea02f30997bcp+1 * lgt + 0x1.59d470fcee60dp+23 * rt_inv);
  k[r139f] = 0x1.af89a2p+30 * exp(-0x1.b5a4680000001p+25 * rt_inv);
  k[r139b] = 0x1.7b819054a7bb4p+36 *
             exp(-0x1.f02aea59ba39ep-1 * lgt - 0x1.2d663c1ece0dp+23 * rt_inv);
  k[r140f] = 0x1.2ecaa6p+32 * exp(0x1.6854200000001p+22 * rt_inv);
  k[r140b] = 0x1.65a41ec1a2cf1p+42 *
             exp(-0x1.0bab580409b4p+0 * lgt - 0x1.b9ad1b7bf265ap+26 * rt_inv);
  k[r141f] = 0x1.d71d78p+27 * exp(0x1.90eecp+22 * rt_inv);
  k[r141b] = 0x1.08ffc31a7a2a8p+38 *
             exp(-0x1.9f6c3594acb6p-1 * lgt - 0x1.13e5c57c7db58p+27 * rt_inv);
  k[r142] = 0x1.21a42d98p+38 *
            exp(-0x1.9c28f5c28f5c3p+0 * lgt + 0x1.0c65200000001p+22 * rt_inv);
  k[r143] = 0x1.977420dcp+43 *
            exp(-0x1.9c28f5c28f5c3p+0 * lgt - 0x1.dafd800000001p+22 * rt_inv);
  k[r144f] = 0x1.65a0bcp+36;
  k[r144b] = 0x1.16446dd007b74p+26 *
             exp(0x1.b30a3314e72cp-2 * lgt - 0x1.43c0c94bddffcp+27 * rt_inv);
  k[r145f] = 0x1.0c388dp+35;
  k[r145b] = 0x1.8add16301b4p+33 *
             exp(-0x1.b618fa2a92p-8 * lgt - 0x1.ca4c8baa1311bp+27 * rt_inv);
  k[r146f] = 0x1.bf08ebp+35;
  k[r146b] = 0x1.340c6460186f7p+34 *
             exp(0x1.70d1c1e3c2ep-2 * lgt - 0x1.deddae548600cp+27 * rt_inv);
  k[r147f] = 0x1.1ef2116d38p+49 * exp(-0x1.519194p+27 * rt_inv);
  k[r147b] = 0x1.14d6cdb14f25ap+18 *
             exp(0x1.40244b2af5737p+0 * lgt + 0x1.2002a4b17e509p+23 * rt_inv);
  k[r148f] = 0x1.bf08ebp+33 * exp(-0x1.3272p+21 * rt_inv);
  k[r148b] = 0x1.01d969da34351p+33 *
             exp(-0x1.f61f36009adp-6 * lgt - 0x1.2edea76a9e2ebp+25 * rt_inv);
  k[r149f] = 0x1.0c388dp+33 * exp(-0x1.3272p+21 * rt_inv);
  k[r149b] = 0x1.356b4bd2a48bdp+32 *
             exp(-0x1.f61f3600975p-6 * lgt - 0x1.2edea76a9e2cp+25 * rt_inv);
  k[r150f] = 0x1.bf08ebp+34;
  k[r150b] = 0x1.55e8417b68d9cp+30 *
             exp(0x1.c5d5a87e8122p-3 * lgt - 0x1.76bea9e3eb96fp+25 * rt_inv);
  k[r151f] = 0x1.bf08ebp+33;
  k[r151b] = 0x1.667c07adb321ap+35 *
             exp(0x1.612066ac3d8p-3 * lgt - 0x1.773cb6f35ece1p+29 * rt_inv);
  k[r152f] = 0x1.bf08ebp+33;
  k[r152b] = 0x1.de7356739c58cp+40 *
             exp(-0x1.fbbdcd3a7b18p-2 * lgt - 0x1.9047986d92e87p+28 * rt_inv);
  k[r153f] = 0x1.bf08ebp+34;
  k[r153b] = 0x1.bb073eb4d7183p+50 *
             exp(-0x1.9e2a36d2d93cp-1 * lgt - 0x1.59a80f41269b6p+28 * rt_inv);
  k[r154f] = 0x1.04c533cp+36;
  k[r154b] = 0x1.34ec547c41e4cp+44 *
             exp(-0x1.25bc262cf35fp-1 * lgt - 0x1.00f90ca07a0f5p+26 * rt_inv);
  k[r155] = 0x1.a13b86p+34;
  k[r156f] = 0x1.65a0bcp+33;
  k[r156b] = 0x1.9a021c46692e1p+29 *
             exp(0x1.012fcf620108p-1 * lgt - 0x1.738903e909849p+29 * rt_inv);
  k[r157f] = 0x1.bf08ebp+34;
  k[r157b] = 0x1.01d969da338f6p+34 *
             exp(-0x1.f61f36009a5p-6 * lgt - 0x1.1bb7876a9e2bdp+25 * rt_inv);
  k[r158f] = 0x1.0c388dp+33;
  k[r158b] = 0x1.356b4bd2a41c9p+32 *
             exp(-0x1.f61f3600992p-6 * lgt - 0x1.1bb7876a9e2b8p+25 * rt_inv);
  k[r159f] = 0x1.a13b86p+32;
  k[r159b] = 0x1.e15192641bc1fp+31 *
             exp(-0x1.f61f3600996p-6 * lgt - 0x1.1bb7876a9e2b7p+25 * rt_inv);
  k[r160f] = 0x1.a13b86p+33;
  k[r160b] = 0x1.491e5391cdd05p+27 *
             exp(0x1.dab274ade1p-2 * lgt - 0x1.e66dd7df691b7p+27 * rt_inv);
  kTroe0 = 0x1.5af1d78b58c4p+71 *
           exp(-0x1.91eb851eb851fp+1 * lgt - 0x1.3a1b400000001p+22 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+44 * exp(-0x1.999999999999ap-1 * lgt);
  fcTroe = 0x1.47ae147ae147ap-2 * exp(-temp / 0x1.38p+6) +
           0x1.5c28f5c28f5c3p-1 * exp(-temp / 0x1.f2cp+10) +
           0x1p+0 * exp(-0x1.5d6p+12 / temp);
  k[r161f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM21]);
  kTroe0 = 0x1.4e3fab3860958p+88 *
           exp(-0x1.ab0ea034422ap+1 * lgt - 0x1.bebdc29525385p+28 * rt_inv);
  kTroeInf = 0x1.5e7c343372ba3p+61 *
             exp(-0x1.fe2605efc0f9ep-1 * lgt - 0x1.b9d5559525385p+28 * rt_inv);
  k[r161b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM21]);
  k[r162f] = 0x1.3be795p+33 * exp(-0x1.7f0e800000001p+22 * rt_inv);
  k[r162b] = 0x1.35d0d20197661p+32 *
             exp(-0x1.13ca16e7038p-5 * lgt - 0x1.2f870051249bdp+28 * rt_inv);
  k[r163] = 0x1.3ab668p+31 * exp(-0x1.7f0e800000001p+22 * rt_inv);
  k[r164] = 0x1.74876e8p+35;
  k[r165f] = 0x1.c6bf52634p+49 * exp(-0x1.8f5c28f5c28f6p+0 * lgt);
  k[r165b] = 0x1.2d7f289d2867fp+46 *
             exp(-0x1.4ec8f70df00ap+0 * lgt - 0x1.6c1c89e535b17p+23 * rt_inv);
  k[r166f] = 0x1.017df8p+28 *
             exp(0x1.570a3d70a3d71p-1 * lgt - 0x1.9a3098p+26 * rt_inv);
  k[r166b] = 0x1.556f001936368p+24 *
             exp(0x1.d830a1404b9p-1 * lgt - 0x1.c7b4293ca6c81p+26 * rt_inv);
  k[r167f] = 0x1.611ffffffffffp+13 *
             exp(0x1p+1 * lgt - 0x1.7f0e800000001p+23 * rt_inv);
  k[r167b] = 0x1.3ca4074fe040ep+13 *
             exp(0x1.13473381f361p+1 * lgt - 0x1.4cdf4da666986p+26 * rt_inv);
  k[r168f] = 0x1.ebbd028p+34;
  k[r168b] = 0x1.35bbcc9f2cb2fp+38 *
             exp(-0x1.02a31b4504cp-2 * lgt - 0x1.23be655ad4c06p+28 * rt_inv);
  k[r169f] = 0x1.a8aedf4p+35;
  k[r169b] = 0x1.bd4610b428d7dp+41 *
             exp(-0x1.92d50749078p-5 * lgt - 0x1.5fd0cc552013fp+29 * rt_inv);
  k[r170f] = 0x1.bf08ebp+34;
  k[r170b] = 0x1.0a5f591d011a1p+47 *
             exp(-0x1.5e42d9214b2cp-1 * lgt - 0x1.67044689ef498p+28 * rt_inv);
  k[r171f] = 0x1.fe83874p+33 * exp(0x1.819c4p+21 * rt_inv);
  k[r171b] = 0x1.e91ba217732p+50 *
             exp(-0x1.db960b2aa1a8p-1 * lgt - 0x1.d7dd0cd472606p+27 * rt_inv);
  k[r172f] = 0x1.954fc4p+30 * exp(-0x1.5ddbc00000001p+21 * rt_inv);
  k[r172b] = 0x1.8075b696d1299p+20 *
             exp(0x1.2d4098087ea8p-1 * lgt - 0x1.034edab87d3c1p+28 * rt_inv);
  kTroe0 = 0x1.0a7c08982a5e9p+86 *
           exp(-0x1.ep+1 * lgt - 0x1.f558333333335p+21 * rt_inv);
  kTroeInf = 0x1.0913e359p+41 *
             exp(-0x1.6147ae147ae14p-1 * lgt - 0x1.6550b33333333p+19 * rt_inv);
  fcTroe = 0x1p+0 * exp(-temp / 0x1.1dp+9) +
           0x1p+0 * exp(-0x1.93e5939a08ceap+99 / temp);
  k[r173f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM22]);
  kTroe0 = 0x1.45b44519de15p+123 *
           exp(-0x1.47a1c5059e178p+2 * lgt - 0x1.72757f73f0bedp+28 * rt_inv);
  kTroeInf = 0x1.43fc174495a92p+78 *
             exp(-0x1.079575905ae75p+1 * lgt - 0x1.6f3d776723f2p+28 * rt_inv);
  k[r173b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM22]);
  kTroe0 = 0x1.329ba8fa506b9p+117 *
           exp(-0x1.c51eb851eb852p+2 * lgt - 0x1.aac9f8p+24 * rt_inv);
  kTroeInf = 0x1.d9d8c3ed9p+48 *
             exp(-0x1.fae147ae147aep-1 * lgt - 0x1.937c8p+22 * rt_inv);
  fcTroe = 0x1.4395810624dd4p-3 * exp(-temp / 0x1.f4p+6) +
           0x1.af1a9fbe76c8bp-1 * exp(-temp / 0x1.156p+11) +
           0x1p+0 * exp(-0x1.ae2p+12 / temp);
  k[r174f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM23]);
  kTroe0 = 0x1.2976328b8db13p+137 *
           exp(-0x1.d095f9638ad4p+2 * lgt - 0x1.aef0ca5432e0cp+28 * rt_inv);
  kTroeInf = 0x1.cbb628325a0aap+68 *
             exp(-0x1.2b4da81d8778fp+0 * lgt - 0x1.9a921cd432e0cp+28 * rt_inv);
  k[r174b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM23]);
  k[r175f] =
      0x1.c138p+16 * exp(0x1.e666666666666p+0 * lgt - 0x1.e0bc7p+24 * rt_inv);
  k[r175b] = 0x1.b2646a25dedep+4 *
             exp(0x1.3681656023f0ap+1 * lgt - 0x1.36acfc313cdcfp+25 * rt_inv);
  k[r176f] =
      0x1.bbcp+11 * exp(0x1.3333333333333p+1 * lgt - 0x1.74341p+24 * rt_inv);
  k[r176b] = 0x1.f7dd30de186fp-2 *
             exp(0x1.723d07793f1ccp+1 * lgt - 0x1.a78962d4dc952p+24 * rt_inv);
  k[r177f] = 0x1.ce7ffffffffffp+13 *
             exp(0x1.e666666666666p+0 * lgt - 0x1.e5348p+21 * rt_inv);
  k[r177b] = 0x1.2e6f52e901bacp+5 *
             exp(0x1.297effee2dd74p+1 * lgt - 0x1.21caf2825e48p+26 * rt_inv);
  k[r178f] = 0x1.c1451f6p+35 * exp(-0x1.9df0cap+27 * rt_inv);
  k[r178b] = 0x1.4bfb97dc72a3ap+26 *
             exp(0x1.518b0558e848p-3 * lgt + 0x1.d91fc6f6b78a5p+20 * rt_inv);
  k[r179f] = 0x1.1f4f50a02b841p-11 *
             exp(0x1p+2 * lgt - 0x1.084f200000001p+25 * rt_inv);
  k[r179b] = 0x1.eadb58c30d3aap-13 *
             exp(0x1.009667f45f79bp+2 * lgt - 0x1.92721bf97c95ap+25 * rt_inv);
  k[r180f] = 0x1.1b71758e21965p-5 *
             exp(0x1.ce147ae147ae1p+1 * lgt - 0x1.0e0e100000001p+26 * rt_inv);
  k[r180b] = 0x1.350792621a99bp-9 *
             exp(0x1.c3f4850eb9b9ep+1 * lgt - 0x1.893079636111cp+23 * rt_inv);
  k[r181f] = 0x1.3dd97f62b6ae7p-6 *
             exp(0x1.d1eb851eb851fp+1 * lgt - 0x1.10ed880000001p+26 * rt_inv);
  k[r181b] = 0x1.083a19944352ep-5 *
             exp(0x1.7f29280f1ac05p+1 * lgt - 0x1.d03d7ba40f149p+21 * rt_inv);
  k[r182f] = 0x1.cbabc8p+27 * exp(-0x1.c4a53p+24 * rt_inv);
  k[r182b] = 0x1.e89fdc889778p+14 *
             exp(0x1.b7fc63fe667cp-1 * lgt - 0x1.4e41cc143a062p+25 * rt_inv);
  k[r183f] = 0x1.99c82ccp+36 * exp(0x1.0996000000001p+20 * rt_inv);
  k[r183b] = 0x1.2ad64fd7c2863p+28 *
             exp(0x1.1824c9c83a614p-2 * lgt + 0x1.d050e48109c9ep+21 * rt_inv);
  k[r184f] = 0x1.bf08ebp+36;
  k[r184b] = 0x1.000d444c1c5adp+33 *
             exp(-0x1.8835d7930f4p-5 * lgt - 0x1.24206eb91874ap+26 * rt_inv);
  kTroe0 = 0x1.17d93124d9da2p+110 *
           exp(-0x1.a916872b020c5p+2 * lgt - 0x1.704f180000001p+24 * rt_inv);
  kTroeInf =
      0x1.d33c8p+19 * exp(0x1.76872b020c49cp+0 * lgt - 0x1.5a072p+22 * rt_inv);
  fcTroe = 0x1.91a9fbe76c8b4p+0 * exp(-temp / 0x1.2bp+8) +
           -0x1.2353f7ced9168p-1 * exp(-temp / -0x1.1dd8p+13) +
           0x1p+0 * exp(-0x1.30ccccccccccdp+7 / temp);
  k[r185f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM24]);
  kTroe0 = 0x1.204b6e489e391p+121 *
           exp(-0x1.b219c9a078f18p+2 * lgt - 0x1.4a5675c22b207p+27 * rt_inv);
  kTroeInf = 0x1.e15695c1d140bp+30 *
             exp(0x1.527a212c30b5p+0 * lgt - 0x1.271ccbc22b207p+27 * rt_inv);
  k[r185b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM24]);
  k[r186f] = 0x1.176592ep+37 * exp(-0x1.9f74f4p+26 * rt_inv);
  k[r186b] = 0x1.e060a28630d36p+49 *
             exp(-0x1.2c211e791bc58p+0 * lgt - 0x1.e787d712085adp+24 * rt_inv);
  k[r187f] = 0x1.176592ep+37 * exp(-0x1.9f74f4p+26 * rt_inv);
  k[r187b] = 0x1.07771e887a7efp+41 *
             exp(-0x1.2d3d1a3eacep-3 * lgt - 0x1.98e1cc3fd615cp+24 * rt_inv);
  k[r188f] = 0x1.c0e5c15p+38 * exp(-0x1.1d6ac7p+28 * rt_inv);
  k[r188b] = 0x1.6e56a77ca8197p+29 *
             exp(0x1.03857eeb8729p-3 * lgt + 0x1.d5970067ee317p+23 * rt_inv);
  k[r189f] = 0x1.7999999999999p+3 *
             exp(0x1.399999999999ap+1 * lgt + 0x1.74f8300000001p+23 * rt_inv);
  k[r189b] = 0x1.2fc4445ed63c8p+11 *
             exp(0x1.35de6c4a0778p+1 * lgt - 0x1.0bba7f724ce0ep+28 * rt_inv);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 4.74e+09;
  PlogA[1] = 2.57e+10;
  PlogA[2] = 3.1e+11;
  PlogA[3] = 2.15e+07;
  PlogA[4] = 0.1032;
  PlogB[0] = 0.105;
  PlogB[1] = -0.096;
  PlogB[2] = -0.362;
  PlogB[3] = 0.885;
  PlogB[4] = 3.23;
  PlogE[0] = 4.46194e+07;
  PlogE[1] = 4.77231e+07;
  PlogE[2] = 5.59505e+07;
  PlogE[3] = 5.662e+07;
  PlogE[4] = 4.70118e+07;
  np = 5;
  k[r190f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 7.82694e+14;
  PlogA[1] = 4.24372e+15;
  PlogA[2] = 5.11889e+16;
  PlogA[3] = 3.5502e+12;
  PlogA[4] = 17040.9;
  PlogB[0] = -1.08509;
  PlogB[1] = -1.28609;
  PlogB[2] = -1.55209;
  PlogB[3] = -0.305095;
  PlogB[4] = 2.03991;
  PlogE[0] = 5.06289e+06;
  PlogE[1] = 8.16658e+06;
  PlogE[2] = 1.6394e+07;
  PlogE[3] = 1.70634e+07;
  PlogE[4] = 7.4553e+06;
  np = 5;
  k[r190b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r191f] = 0x1.dcd65p+30;
  k[r191b] = 0x1.b23cfc50b2241p+27 *
             exp(0x1.f331a7a4fdfp-2 * lgt - 0x1.0ee7b9f944ea7p+28 * rt_inv);
  k[r192f] = 0x1.99c82ccp+36;
  k[r192b] = 0x1.b0c1650fb8438p+46 *
             exp(-0x1.069121fbd08p-1 * lgt - 0x1.30c3bb35e53aap+28 * rt_inv);
  k[r193f] = 0x1.47d357p+33;
  k[r193b] = 0x1.406679f07806ap+43 *
             exp(-0x1.9dd91531aa37p-1 * lgt - 0x1.a5a6dea824b41p+26 * rt_inv);
  k[r194f] = 0x1.dcd65p+32 * exp(0x1.febep+21 * rt_inv);
  k[r194b] = 0x1.be7c6b9f852e4p+44 *
             exp(-0x1.2484eb73fd01p+0 * lgt - 0x1.da2e9e2b38b59p+26 * rt_inv);
  k[r195f] = 0x1.4689cp+25 * exp(-0x1.18246p+22 * rt_inv);
  k[r195b] = 0x1.74f99c7f4f098p+16 *
             exp(0x1.74393eb4a27ap-1 * lgt - 0x1.10e967b9830a3p+27 * rt_inv);
  k[r196f] = 0x1.1e1a3p+28 * exp(-0x1.9482p+24 * rt_inv);
  k[r196b] = 0x1.e669fb5b96f15p+67 *
             exp(-0x1.167d8fd16d5d4p+1 * lgt - 0x1.3f6ee97e24245p+26 * rt_inv);
  k[r197f] =
      0x1.6828p+15 * exp(0x1.b5c28f5c28f5cp+0 * lgt - 0x1.c4a53p+24 * rt_inv);
  k[r197b] = 0x1.c4189431186afp+34 *
             exp(0x1.eff498cce99fp-1 * lgt - 0x1.6ea35ac1ea718p+26 * rt_inv);
  k[r198f] = 0x1.da73f6p+30 * exp(-0x1.74341p+25 * rt_inv);
  k[r198b] = 0x1.9b7b306406e39p+26 *
             exp(0x1.983605060655p-3 * lgt - 0x1.c6b6e49437ff7p+24 * rt_inv);
  k[r199f] = 0x1.593ae8p+27 * exp(-0x1.26f42p+26 * rt_inv);
  k[r199b] = 0x1.708501a93b8b1p+20 *
             exp(0x1.7a6c85023671cp-2 * lgt + 0x1.a193e6a152ddp+23 * rt_inv);
  k[r200f] = 0x1.af89a2p+30 * exp(-0x1.b5a4680000001p+25 * rt_inv);
  k[r200b] = 0x1.a048aaab2aa44p+27 *
             exp(0x1.cc80c08d2229p-5 * lgt - 0x1.20344cf4d2ff8p+22 * rt_inv);
  k[r201f] = 0x1.0b076p+24 * exp(0x1.a22b9p+23 * rt_inv);
  k[r201b] = 0x1.49840a90d71dp+25 *
             exp(0x1.b61b03377d8p-3 * lgt - 0x1.f8eb82446ed9fp+26 * rt_inv);
  k[r202f] = 0x1.19ce075f6fd21p-7 *
             exp(0x1.e147ae147ae14p+1 * lgt - 0x1.12862p+26 * rt_inv);
  k[r202b] = 0x1.00f723a62739ep-15 *
             exp(0x1.08e1074ec031cp+2 * lgt + 0x1.c77f6bb60bf5dp+19 * rt_inv);
  k[r203f] = 0x1.1ef2116d38p+49 * exp(-0x1.519194p+27 * rt_inv);
  k[r203b] = 0x1.42119564acb64p+18 *
             exp(0x1.0833b833932f6p+0 * lgt + 0x1.847a276d47ed5p+23 * rt_inv);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 3.398e+50;
  PlogA[1] = 9.362e+56;
  PlogA[2] = 1.262e+57;
  PlogB[0] = -13.9;
  PlogB[1] = -15.28;
  PlogB[2] = -14.91;
  PlogE[0] = 3.88233e+07;
  PlogE[1] = 5.95802e+07;
  PlogE[2] = 6.79482e+07;
  np = 3;
  k[r204f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 2.63637e+62;
  PlogA[1] = 7.26359e+68;
  PlogA[2] = 9.79134e+68;
  PlogB[0] = -15.9543;
  PlogB[1] = -17.3343;
  PlogB[2] = -16.9643;
  PlogE[0] = 1.93005e+08;
  PlogE[1] = 2.13762e+08;
  PlogE[2] = 2.2213e+08;
  np = 3;
  k[r204b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 2.094e+06;
  PlogA[1] = 18430;
  PlogA[2] = 7.561e+11;
  PlogB[0] = 0.49;
  PlogB[1] = 1.13;
  PlogB[2] = -1.01;
  PlogE[0] = -1.63762e+06;
  PlogE[1] = -3.01499e+06;
  PlogE[2] = 1.98698e+07;
  np = 3;
  k[r205f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.4572e+06;
  PlogA[1] = 12825.3;
  PlogA[2] = 5.26164e+11;
  PlogB[0] = 0.616486;
  PlogB[1] = 1.25649;
  PlogB[2] = -0.883514;
  PlogE[0] = 5.42497e+07;
  PlogE[1] = 5.28723e+07;
  PlogE[2] = 7.57571e+07;
  np = 3;
  k[r205b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r206f] = 0x1.b1209edbf8b9bp-8 *
             exp(0x1.c147ae147ae14p+1 * lgt - 0x1.c401cp+25 * rt_inv);
  k[r206b] = 0x1.2d68f2b81cbd2p-8 *
             exp(0x1.d1786131b7888p+1 * lgt - 0x1.b73252b09a4bcp+26 * rt_inv);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.303;
  PlogA[1] = 0.2438;
  PlogA[2] = 4.621e+06;
  PlogB[0] = 1.93;
  PlogB[1] = 2.18;
  PlogB[2] = 0.15;
  PlogE[0] = -2.1033e+06;
  PlogE[1] = -261500;
  PlogE[2] = 2.26313e+07;
  np = 3;
  k[r207f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 376.971;
  PlogA[1] = 70.5339;
  PlogA[2] = 1.3369e+09;
  PlogB[0] = 1.57804;
  PlogB[1] = 1.82804;
  PlogB[2] = -0.201956;
  PlogE[0] = 1.35758e+08;
  PlogE[1] = 1.37599e+08;
  PlogE[2] = 1.60492e+08;
  np = 3;
  k[r207b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 4.908e-09;
  PlogA[1] = 6.803e-05;
  PlogA[2] = 0.8265;
  PlogB[0] = 4.76;
  PlogB[1] = 3.57;
  PlogB[2] = 2.41;
  PlogE[0] = 1.06399e+06;
  PlogE[1] = 1.10583e+07;
  PlogE[2] = 2.21124e+07;
  np = 3;
  k[r208f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.0958e-08;
  PlogA[1] = 0.000151889;
  PlogA[2] = 1.84531;
  PlogB[0] = 4.6787;
  PlogB[1] = 3.4887;
  PlogB[2] = 2.3287;
  PlogE[0] = 2.50092e+08;
  PlogE[1] = 2.60086e+08;
  PlogE[2] = 2.7114e+08;
  np = 3;
  k[r208b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.237e+35;
  PlogA[1] = 1.687e+36;
  PlogA[2] = 2.52e+41;
  PlogB[0] = -9.42;
  PlogB[1] = -9.22;
  PlogB[2] = -10.2;
  PlogE[0] = 1.5213e+08;
  PlogE[1] = 1.61921e+08;
  PlogE[2] = 1.82883e+08;
  np = 3;
  k[r209f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 3.5597e+23;
  PlogA[1] = 4.85465e+24;
  PlogA[2] = 7.25176e+29;
  PlogB[0] = -7.44704;
  PlogB[1] = -7.24704;
  PlogB[2] = -8.22704;
  PlogE[0] = 2.46976e+08;
  PlogE[1] = 2.56767e+08;
  PlogE[2] = 2.77729e+08;
  np = 3;
  k[r209b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.782e+32;
  PlogA[1] = 2.701e+37;
  PlogA[2] = 1.98e+38;
  PlogB[0] = -7.1;
  PlogB[1] = -8.47;
  PlogB[2] = -8.46;
  PlogE[0] = 1.37403e+08;
  PlogE[1] = 1.49955e+08;
  PlogE[2] = 1.58574e+08;
  np = 3;
  k[r210f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 1.59833e+20;
  PlogA[1] = 2.42261e+25;
  PlogA[2] = 1.77592e+26;
  PlogB[0] = -4.91926;
  PlogB[1] = -6.28926;
  PlogB[2] = -6.27926;
  PlogE[0] = 3.91079e+07;
  PlogE[1] = 5.16599e+07;
  PlogE[2] = 6.02789e+07;
  np = 3;
  k[r210b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 5.778e+45;
  PlogA[1] = 1.916e+43;
  PlogA[2] = 3.965e+43;
  PlogB[0] = -11.9;
  PlogB[1] = -10.75;
  PlogB[2] = -10.46;
  PlogE[0] = 1.72046e+07;
  PlogE[1] = 1.77402e+08;
  PlogE[2] = 1.90707e+08;
  np = 3;
  k[r211f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 4000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+06;
  PlogA[0] = 2.15456e+36;
  PlogA[1] = 7.14458e+33;
  PlogA[2] = 1.47851e+34;
  PlogB[0] = -10.1977;
  PlogB[1] = -9.0477;
  PlogB[2] = -8.7577;
  PlogE[0] = 883563;
  PlogE[1] = 1.61081e+08;
  PlogE[2] = 1.74386e+08;
  np = 3;
  k[r211b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r212f] = 0x1.081e04d7cp+45 * exp(-0x1.c879d00000001p+27 * rt_inv);
  k[r212b] = 0x1.79b0f3591a5eep+1 *
             exp(0x1.b46bbdfd977eap+0 * lgt + 0x1.39b8007c080eap+23 * rt_inv);
  k[r213f] = 0x1.af24e6a58p+42 * exp(-0x1.ad57b8p+27 * rt_inv);
  k[r213b] = 0x1.a9e31b257c34ap+35 *
             exp(0x1.152775ecf458p-2 * lgt - 0x1.40b041804284cp+28 * rt_inv);
  kTroe0 = 0x1.0cd7b0f80ba63p+186 *
           exp(-0x1.699999999999ap+3 * lgt - 0x1.7eb51ecp+28 * rt_inv);
  kTroeInf = 0x1.4c0973485bf39p+74 *
             exp(-0x1.bd70a3d70a3d7p+0 * lgt - 0x1.5892488p+28 * rt_inv);
  fcTroe = 0x1.feb9a176ddacfp-1 * exp(-temp / 0x1.670cccccccccdp+9) +
           0x1.465e892253112p-9 * exp(-temp / 0x1.85b22d0e56042p+2) +
           0x1p+0 * exp(-0x1.d88p+11 / temp);
  k[r214f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM25]);
  kTroe0 = 0x1.8531a5ce11c3ep+149 *
           exp(-0x1.3bb55d894db2p+3 * lgt - 0x1.353089deea582p+25 * rt_inv);
  kTroeInf = 0x1.e0ad9f75d4e26p+37 *
             exp(-0x1.393b0d52ac01cp-2 * lgt - 0x1.0675f7ba95f36p+19 * rt_inv);
  k[r214b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM25]);
  kTroe0 = 0x1.ddc1b68ff01dbp+182 *
           exp(-0x1.699999999999ap+3 * lgt - 0x1.7eb51ecp+28 * rt_inv);
  kTroeInf = 0x1.26e72a69a50dp+71 *
             exp(-0x1.bd70a3d70a3d7p+0 * lgt - 0x1.5892488p+28 * rt_inv);
  fcTroe = 0x1.feb9a176ddacfp-1 * exp(-temp / 0x1.670cccccccccdp+9) +
           0x1.465e892253112p-9 * exp(-temp / 0x1.85b22d0e56042p+2) +
           0x1p+0 * exp(-0x1.d88p+11 / temp);
  k[r215f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM26]);
  kTroe0 = 0x1.c9c817da78a17p+151 *
           exp(-0x1.36db4615d54cp+3 * lgt - 0x1.8d526dae0ff35p+28 * rt_inv);
  kTroeInf = 0x1.1a92b39e18a6p+40 *
             exp(-0x1.3bf03dc73e838p-3 * lgt - 0x1.672f976e0ff34p+28 * rt_inv);
  k[r215b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM26]);
  k[r216f] = 0x1.06p+7 *
             exp(0x1.4a3d70a3d70a4p+1 * lgt - 0x1.378d800000001p+22 * rt_inv);
  k[r216b] = 0x1.48b9ec4a246ccp-4 *
             exp(0x1.9ed31a89fece4p+1 * lgt - 0x1.e3648c003fcb6p+25 * rt_inv);
  k[r217f] = 0x1.5c28f5c28f5c2p+1 *
             exp(0x1.8cccccccccccdp+1 * lgt - 0x1.4c9ef00000001p+24 * rt_inv);
  k[r217b] = 0x1.89bc1a01ae554p-6 *
             exp(0x1.b1dd0a1341bccp+1 * lgt - 0x1.b39d32dc56909p+25 * rt_inv);
  k[r218f] = 0x1.620d35p+32 * exp(-0x1.dd088p+22 * rt_inv);
  k[r218b] = 0x1.04cec7833402ap+21 *
             exp(0x1.41452ffd0bf4p-1 * lgt - 0x1.cb6fd139713cfp+25 * rt_inv);
  k[r219f] = 0x1.91bc3dp+31 * exp(0x1.3c264p+21 * rt_inv);
  k[r219b] = 0x1.54dac082e62e4p+24 *
             exp(0x1.1e4d11d0c6fcp-1 * lgt - 0x1.c3d22869dfc4ep+26 * rt_inv);
  k[r220f] = 0x1.c086634p+34 * exp(-0x1.386e3ap+27 * rt_inv);
  k[r220b] = 0x1.ae08ecd1bed5ep+26 *
             exp(0x1.330140762a948p-2 * lgt + 0x1.7579334596821p+22 * rt_inv);
  k[r221f] = 0x1.7c1ac7705b4b9p-21 *
             exp(0x1.251eb851eb852p+2 * lgt - 0x1.f60f400000001p+22 * rt_inv);
  k[r221b] = 0x1.a54b8c1027283p-21 *
             exp(0x1.2e58dc22e67f2p+2 * lgt - 0x1.1f847de43fc33p+26 * rt_inv);
  k[r222f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r222b] = 0x1.fb98f97079b4p+28 *
             exp(0x1.c9e079aa465p-5 * lgt - 0x1.38d4d627db34ap+25 * rt_inv);
  k[r223f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r223b] = 0x1.8306d7fa30df9p+33 *
             exp(-0x1.05eb955999fep-1 * lgt - 0x1.db9b3f1287befp+24 * rt_inv);
  k[r224f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r224b] = 0x1.cbac035d1cc95p+31 *
             exp(-0x1.3dab628e2068p-4 * lgt - 0x1.9bbe0000eb2f2p+25 * rt_inv);
  k[r225f] = 0x1.5d3ef798p+41 * exp(-0x1.1374bc6a7ef9ep+0 * lgt);
  k[r225b] = 0x1.8e161478a26fcp+32 *
             exp(-0x1.4455188ec1dp-4 * lgt - 0x1.7c75ea44428c9p+26 * rt_inv);
  k[r226f] = 0x1.58p+7 *
             exp(0x1.3333333333333p+1 * lgt - 0x1.a041400000001p+21 * rt_inv);
  k[r226b] = 0x1.0713269107c3dp+4 *
             exp(0x1.4b410b07b214cp+1 * lgt - 0x1.8322d3d7eb241p+26 * rt_inv);
  kTroe0 = 0x1.412a522fb2p+52 *
           exp(-0x1.f0a3d70a3d70ap-1 * lgt - 0x1.d20d6p+25 * rt_inv);
  kTroeInf = 0x1.f241f098p+39 *
             exp(0x1.428f5c28f5c29p-1 * lgt - 0x1.0dbc58p+26 * rt_inv);
  fcTroe = 0x1.7be76c8b43958p-2 * exp(-temp / 0x1.042c9d4p+33) +
           0x1.420c49ba5e354p-1 * exp(-temp / 0x1.6147ae147ae14p+2) +
           0x1p+0 * exp(-0x1.21eacp+26 / temp);
  k[r227f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM27]);
  kTroe0 = 0x1.15a646cbcc53ap+21 *
           exp(0x1.e2e085544dfccp-2 * lgt - 0x1.1ac30e9ffeb8p+23 * rt_inv);
  kTroeInf = 0x1.aebf3734527fcp+8 *
             exp(0x1.0928dd77568c6p+1 * lgt - 0x1.2038274fff5c1p+24 * rt_inv);
  k[r227b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM27]);
  k[r228f] = 0x1.2a05f2p+34;
  k[r228b] = 0x1.ed29b764b183cp+31 *
             exp(0x1.31000d382428p-2 * lgt - 0x1.e134fc9777179p+27 * rt_inv);
  k[r229f] = 0x1.2a05f2p+34;
  k[r229b] = 0x1.218ab280c49b5p+31 *
             exp(0x1.0edd1e00fe08p-2 * lgt - 0x1.d60bf5e5c373p+27 * rt_inv);
  k[r230f] = 0x1.74876e8p+35;
  k[r230b] = 0x1.1048772713a8fp+44 *
             exp(-0x1.c01609d2d49p-3 * lgt - 0x1.f22a0a8987052p+27 * rt_inv);
  k[r231f] = 0x1.c9c38p+26 * exp(0x1.18e8800000001p+22 * rt_inv);
  k[r231b] = 0x1.44bb099fd7de1p+66 *
             exp(-0x1.0214d6412fa38p+1 * lgt - 0x1.3b6060d4c5edfp+27 * rt_inv);
  k[r232f] = 0x1.0b076p+24 * exp(0x1.a22b9p+23 * rt_inv);
  k[r232b] = 0x1.64ca312e97913p+32 *
             exp(-0x1.826c1919addcp-2 * lgt - 0x1.31c69a9a6788ep+27 * rt_inv);
  k[r233f] = 0x1.1f4b5dp+31 * exp(-0x1.3d2bcp+25 * rt_inv);
  k[r233b] = 0x1.042b45a93a8b5p+34 *
             exp(-0x1.114dcfb19ddcp-3 * lgt - 0x1.a014e9d9100abp+25 * rt_inv);
  k[r234f] = 0x1.593ae8p+27 * exp(-0x1.26f42p+26 * rt_inv);
  k[r234b] = 0x1.8f059cf16e622p+27 *
             exp(-0x1.c61a2b666c63p-3 * lgt - 0x1.b379b0e1aee4p+23 * rt_inv);
  k[r235f] = 0x1.da73f6p+30 * exp(-0x1.74341p+25 * rt_inv);
  k[r235b] = 0x1.bd8a1d86a6c32p+33 *
             exp(-0x1.915e9832690fp-2 * lgt - 0x1.b89ed82adc6aep+25 * rt_inv);
  k[r236f] = 0x1.faa3b5p+33 * exp(-0x1.468e48p+26 * rt_inv);
  k[r236b] = 0x1.f43873433cc1p+32 *
             exp(-0x1.b34d2cda7f1ap-3 * lgt - 0x1.3635b831e8404p+25 * rt_inv);
  k[r237f] = 0x1.c7a827085p+48 * exp(-0x1.406932p+27 * rt_inv);
  k[r237b] = 0x1.087e14e7e96f8p+13 *
             exp(0x1.d958c0bbea296p+0 * lgt - 0x1.0d3884bd64023p+22 * rt_inv);
  k[r238f] = 0x1.001d1bf8p+42 * exp(-0x1.4f2cbp+25 * rt_inv);
  k[r238b] = 0x1.9aa0c4df7b6d6p+8 *
             exp(0x1.790f33f2913p+0 * lgt - 0x1.0e8d65043ef53p+26 * rt_inv);
  kTroe0 = 0x1.f04ef12cb04cep+88 *
           exp(-0x1.e666666666666p+1 * lgt - 0x1.5a89abb333333p+27 * rt_inv);
  kTroeInf = 0x1.4524f481d8p+50 *
             exp(-0x1.3333333333333p-3 * lgt - 0x1.6be76p+27 * rt_inv);
  fcTroe = 0x1.eb851eb851ecp-7 * exp(-temp / 0x1.89p+8) +
           0x1.f851eb851eb85p-1 * exp(-temp / 0x1.241011p+33) +
           0x1p+0 * exp(-0x1.2a05f2p+32 / temp);
  k[r239f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM29]);
  kTroe0 = 0x1.e59f42191a614p+74 *
           exp(-0x1.bd20a8295f93dp+1 * lgt - 0x1.6b5f4ef77db0bp+24 * rt_inv);
  kTroeInf = 0x1.3e24bd03e5ec3p+36 *
             exp(0x1.6128b09d39f5dp-3 * lgt - 0x1.f64cf15de416cp+24 * rt_inv);
  k[r239b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM29]);
  kTroe0 = 0x1.e0a31db9dc668p+102 *
           exp(-0x1.447ae147ae148p+2 * lgt - 0x1.49969cp+27 * rt_inv);
  kTroeInf = 0x1.5518cdfap+41 *
             exp(0x1.28f5c28f5c28fp-2 * lgt - 0x1.419ba4p+27 * rt_inv);
  fcTroe = 0x1.fffffffffffffp-1 * exp(-temp / 0x1.1f8p+10) +
           0x1.48d02ebc10c3p-54 * exp(-temp / 0x1.296d5b8p+32) +
           0x1p+0 * exp(-0x1.aac4eep+30 / temp);
  k[r240f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM30]);
  kTroe0 = 0x1.ccfe141313b5bp+67 *
           exp(-0x1.a0eb4f828f6e8p+1 * lgt - 0x1.1288bd32fa45p+27 * rt_inv);
  kTroeInf = 0x1.4727d09ab6c18p+6 *
             exp(0x1.0d292b5eb83fap+1 * lgt - 0x1.0a8dc532fa45p+27 * rt_inv);
  k[r240b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM30]);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogA[0] = 188;
  PlogA[1] = 188;
  PlogA[2] = 251;
  PlogA[3] = 70500;
  PlogB[0] = 2.37;
  PlogB[1] = 2.37;
  PlogB[2] = 2.33;
  PlogB[3] = 1.63;
  PlogE[0] = 9.92863e+07;
  PlogE[1] = 1.14516e+08;
  PlogE[2] = 9.95792e+07;
  PlogE[3] = 1.05813e+08;
  np = 4;
  k[r241f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogA[0] = 16.4843;
  PlogA[1] = 16.4843;
  PlogA[2] = 22.0083;
  PlogA[3] = 6181.61;
  PlogB[0] = 2.6781;
  PlogB[1] = 2.6781;
  PlogB[2] = 2.6381;
  PlogB[3] = 1.9381;
  PlogE[0] = 1.46357e+08;
  PlogE[1] = 1.61586e+08;
  PlogE[2] = 1.46649e+08;
  PlogE[3] = 1.52884e+08;
  np = 4;
  k[r241b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogA[0] = 2.68e+14;
  PlogA[1] = 1.52e+17;
  PlogA[2] = 1.65e+16;
  PlogA[3] = 8.953e+10;
  PlogB[0] = -1.84;
  PlogB[1] = -2.58;
  PlogB[2] = -2.22;
  PlogB[3] = -0.6;
  PlogE[0] = 2.73215e+07;
  PlogE[1] = 3.75723e+07;
  PlogE[2] = 4.32626e+07;
  PlogE[3] = 4.23421e+07;
  np = 4;
  k[r242] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                             PlogE, np);
  kTroe0 = 0x1.1623b50660ab2p+91 *
           exp(-0x1.470a3d70a3d71p+2 * lgt - 0x1.c4f6e80000001p+24 * rt_inv);
  kTroeInf = 0x1.823cf4p+29;
  fcTroe = 0x1.a31f8a0902dep-2 * exp(-temp / 0x1.13p+8) +
           0x1.2e703afb7e91p-1 * exp(-temp / 0x1.328p+10) +
           0x1p+0 * exp(-0x1.441p+12 / temp);
  k[r243f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM31]);
  kTroe0 = 0x1.115d5ea388d24p+129 *
           exp(-0x1.b2fe25634ba18p+2 * lgt - 0x1.5b1cd29186575p+28 * rt_inv);
  kTroeInf = 0x1.7b9b8e363639p+67 *
             exp(-0x1.afcf9fca9f29cp+0 * lgt - 0x1.3ecd641186575p+28 * rt_inv);
  k[r243b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM31]);
  kTroe0 = 0x1.098bca859dd0fp+160 *
           exp(-0x1.48a3d70a3d70ap+3 * lgt - 0x1.ba080ap+27 * rt_inv);
  kTroeInf = 0x1.6713d4p+26 *
             exp(0x1.eac083126e979p+0 * lgt - 0x1.670371999999ap+27 * rt_inv);
  fcTroe = 0x1.98adab9f559b4p-2 * exp(-temp / 0x1.e2f9f7cp+32) +
           0x1.33a92a3055326p-1 * exp(-temp / 0x1.4dd999999999ap+9) +
           0x1p+0 * exp(-0x1.2a05f2p+32 / temp);
  k[r244f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM32]);
  kTroe0 = 0x1.d46414f4b9796p+149 *
           exp(-0x1.4a33c2a2e862ap+3 * lgt - 0x1.8483ff8b08dap+25 * rt_inv);
  kTroeInf = 0x1.3caf39a7149c6p+16 *
             exp(0x1.de41264d17079p+0 * lgt - 0x1.c38cef8b7a012p+22 * rt_inv);
  k[r244b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM32]);
  k[r245f] = 0x1.46321b7ap+40 *
             exp(-0x1.5e353f7ced917p-3 * lgt - 0x1.185f366666667p+25 * rt_inv);
  k[r245b] = 0x1.2b08b6ee23127p+32 *
             exp(0x1.7f523e79f51cp-3 * lgt - 0x1.7c655babdb6bbp+24 * rt_inv);
  k[r246f] = 0x1.1eff11cp+36 *
             exp(-0x1.5e353f7ced917p-3 * lgt - 0x1.0b112cccccccdp+24 * rt_inv);
  k[r246b] = 0x1.1952d25851aa4p+15 *
             exp(0x1.51c2c1afee72p+0 * lgt - 0x1.177208a0d75ep+27 * rt_inv);
  k[r247f] = 0x1.a13b86p+30 * exp(-0x1.58c0400000001p+22 * rt_inv);
  k[r247b] = 0x1.25c50333e0558p+22 *
             exp(0x1.748916f6aa0ap-1 * lgt - 0x1.89008082d8bc7p+27 * rt_inv);
  k[r248f] = 0x1.2a05f2p+33 * exp(-0x1.febep+24 * rt_inv);
  k[r248b] = 0x1.40ce7bd80ef9p+24 *
             exp(0x1.4ca0cfc44af1p-2 * lgt - 0x1.e23972a2e2f79p+23 * rt_inv);
  k[r249f] = 0x1.2a05f2p+33 * exp(-0x1.febep+22 * rt_inv);
  k[r249b] = 0x1.7180dba3f948ap+28 *
             exp(0x1.06b0936bc076p-2 * lgt - 0x1.c49f104307027p+25 * rt_inv);
  k[r250f] = 0x1.dcd65p+30 * exp(0x1.01eccp+22 * rt_inv);
  k[r250b] = 0x1.8ef271ad7d6c2p+22 *
             exp(0x1.1d11f78fb32dp-1 * lgt - 0x1.a7ef8ea388347p+26 * rt_inv);
  k[r251f] = 0x1.7d851eb851eb8p+5 *
             exp(0x1.27ef9db22d0e5p+1 * lgt - 0x1.2e3b5p+25 * rt_inv);
  k[r251b] = 0x1.d724121d8baeap+1 *
             exp(0x1.4e5f4a99f5df8p+1 * lgt - 0x1.ec57fef3d96ccp+26 * rt_inv);
  k[r252f] = 0x1.2a05f2p+37;
  k[r252b] = 0x1.8d0a268e603c9p+38 *
             exp(0x1.b66d8512b78p-3 * lgt - 0x1.a47492de59c6fp+28 * rt_inv);
  k[r253] = 0x1.74876e8p+36;
  k[r254] = 0x1.2a05f2p+36;
  k[r255f] = 0x1.74876e8p+36;
  k[r255b] = 0x1.503e932a5e9cfp+15 *
             exp(0x1.b4b68cf72999p+0 * lgt - 0x1.183cfdb63dd9bp+26 * rt_inv);
  k[r256] = 0x1.6c4db8p+27 *
            exp(-0x1.47ae147ae147bp-6 * lgt - 0x1.047a8p+22 * rt_inv);
  k[r257] = 0x1.1ce903p+32 *
            exp(-0x1.22d0e56041893p-3 * lgt - 0x1.25ad4p+22 * rt_inv);
  k[r258f] = 0x1.0cf0bfb97f4p+56 * exp(-0x1.e666666666666p+0 * lgt);
  k[r258b] = 0x1.6d7b5db383b1dp+89 *
             exp(-0x1.bd8c242cfe51p+1 * lgt - 0x1.2826eeb04d91p+28 * rt_inv);
  k[r259f] = 0x1.6069972p+36 * exp(0x1.07084p+21 * rt_inv);
  k[r259b] = 0x1.7fbcfd82d95cdp+56 *
             exp(-0x1.45604db917cep+0 * lgt - 0x1.2eb64de67301ep+28 * rt_inv);
  k[r260f] = 0x1.74876e8p+35;
  k[r260b] = 0x1.7a244295da0e1p+38 *
             exp(0x1.1be60c0f0cap-2 * lgt - 0x1.3914533b89b7cp+29 * rt_inv);
  kTroe0 = 0x1.287626ee52198p+80 *
           exp(-0x1.ee147ae147ae1p+1 * lgt - 0x1.a7ea800000001p+23 * rt_inv);
  kTroeInf = 0x1.6a657p+32 *
             exp(0x1.147ae147ae148p-2 * lgt - 0x1.1e04000000001p+20 * rt_inv);
  fcTroe = 0x1.be76c8b43958p-3 * exp(-temp / 0x1.9fp+7) +
           0x1.90624dd2f1aap-1 * exp(-temp / 0x1.4cep+11) +
           0x1p+0 * exp(-0x1.7cfp+12 / temp);
  k[r261f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM34]);
  kTroe0 = 0x1.763cab3fcc0b6p+100 *
           exp(-0x1.0829abdd760c8p+2 * lgt - 0x1.c77d1c64550d5p+28 * rt_inv);
  kTroeInf = 0x1.c9782a42b40b8p+52 *
             exp(0x1.41fd3d455e8p-9 * lgt - 0x1.bb5bcc64550d6p+28 * rt_inv);
  k[r261b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM34]);
  kTroe0 = 0x1.ea7457a7c18ap+158 *
           exp(-0x1.29eb851eb851fp+3 * lgt - 0x1.8e756e0000001p+28 * rt_inv);
  kTroeInf = 0x1.d1a94a2p+42 *
             exp(0x1.c28f5c28f5c29p-2 * lgt - 0x1.62352bp+28 * rt_inv);
  fcTroe = 0x1.0fdf3b645a1cap-2 * exp(-temp / 0x1.68p+7) +
           0x1.7810624dd2f1bp-1 * exp(-temp / 0x1.02cp+10) +
           0x1p+0 * exp(-0x1.529p+12 / temp);
  k[r262f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM35]);
  kTroe0 = 0x1.3fa8df89b4b2bp+131 *
           exp(-0x1.0f4bd1a85fc3ep+3 * lgt - 0x1.796ad42ccce26p+25 * rt_inv);
  kTroeInf = 0x1.2f801ab3b2da5p+15 *
             exp(0x1.45a172bd01e12p+0 * lgt - 0x1.768bc2ccce1dbp+21 * rt_inv);
  k[r262b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM35]);
  k[r263f] = 0x1.8c17fffffffffp+15 *
             exp(0x1.ee147ae147ae1p+0 * lgt - 0x1.9d61c8p+25 * rt_inv);
  k[r263b] = 0x1.265e560eb8c2fp+3 *
             exp(0x1.45a8ca53fa7dep+1 * lgt - 0x1.67c73f6056ebfp+24 * rt_inv);
  k[r264f] = 0x1.d1cffffffffffp+12 *
             exp(0x1.e147ae147ae14p+0 * lgt - 0x1.75ddp+19 * rt_inv);
  k[r264b] = 0x1.6ed1da626f1f2p-3 *
             exp(0x1.5489f9e18881cp+1 * lgt - 0x1.9e57964b600e8p+26 * rt_inv);
  k[r265f] =
      0x1.7d2p+12 * exp(0x1.e147ae147ae14p+0 * lgt - 0x1.75ddp+19 * rt_inv);
  k[r265b] = 0x1.f3d5021b0e289p+18 *
             exp(0x1.2b532dba3cd58p+0 * lgt - 0x1.11029c60ac839p+26 * rt_inv);
  k[r266f] = 0x1.64ccccccccccdp+4 *
             exp(0x1.5f5c28f5c28f6p+1 * lgt - 0x1.1ae328p+23 * rt_inv);
  k[r266b] = 0x1.66a0c000f30cap-5 *
             exp(0x1.a0f8506722fdap+1 * lgt - 0x1.3c2e7a83ab24bp+25 * rt_inv);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 0.00535;
  PlogA[1] = 0.0319;
  PlogA[2] = 0.555;
  PlogA[3] = 178;
  PlogA[4] = 2.37e+06;
  PlogA[5] = 2.76e+10;
  PlogB[0] = 2.92;
  PlogB[1] = 2.71;
  PlogB[2] = 2.36;
  PlogB[3] = 1.68;
  PlogB[4] = 0.56;
  PlogB[5] = -0.5;
  PlogE[0] = -7.24962e+06;
  PlogE[1] = -4.9049e+06;
  PlogE[2] = -756467;
  PlogE[3] = 8.62113e+06;
  PlogE[4] = 2.5132e+07;
  PlogE[5] = 4.79281e+07;
  np = 6;
  k[r267f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 6.09551e-05;
  PlogA[1] = 0.000363452;
  PlogA[2] = 0.00632338;
  PlogA[3] = 2.02804;
  PlogA[4] = 27002.5;
  PlogA[5] = 3.1446e+08;
  PlogB[0] = 3.38739;
  PlogB[1] = 3.17739;
  PlogB[2] = 2.82739;
  PlogB[3] = 2.14739;
  PlogB[4] = 1.02739;
  PlogB[5] = -0.0326139;
  PlogE[0] = 4.33255e+07;
  PlogE[1] = 4.56702e+07;
  PlogE[2] = 4.98186e+07;
  PlogE[3] = 5.91962e+07;
  PlogE[4] = 7.57071e+07;
  PlogE[5] = 9.85032e+07;
  np = 6;
  k[r267b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 2.37e-10;
  PlogA[1] = 8.73e-08;
  PlogA[2] = 0.000403;
  PlogA[3] = 2.38e-05;
  PlogA[4] = 825000;
  PlogA[5] = 6.8e+06;
  PlogB[0] = 5.3;
  PlogB[1] = 4.57;
  PlogB[2] = 3.54;
  PlogB[3] = 3.91;
  PlogB[4] = 1.01;
  PlogB[5] = 0.81;
  PlogE[0] = -8.57971e+06;
  PlogE[1] = -2.58571e+06;
  PlogE[2] = 7.87303e+06;
  PlogE[3] = 7.20778e+06;
  PlogE[4] = 4.39625e+07;
  PlogE[5] = 5.80208e+07;
  np = 6;
  k[r268f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 3.83486e-06;
  PlogA[1] = 0.00141259;
  PlogA[2] = 6.52088;
  PlogA[3] = 0.385104;
  PlogA[4] = 1.33492e+10;
  PlogA[5] = 1.1003e+11;
  PlogB[0] = 4.33302;
  PlogB[1] = 3.60302;
  PlogB[2] = 2.57302;
  PlogB[3] = 2.94302;
  PlogB[4] = 0.0430163;
  PlogB[5] = -0.156984;
  PlogE[0] = 3.27756e+07;
  PlogE[1] = 3.87696e+07;
  PlogE[2] = 4.92284e+07;
  PlogE[3] = 4.85631e+07;
  PlogE[4] = 8.53179e+07;
  PlogE[5] = 9.93761e+07;
  np = 6;
  k[r268b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 10.4;
  PlogA[1] = 10.7;
  PlogA[2] = 15.2;
  PlogA[3] = 319;
  PlogA[4] = 194000;
  PlogA[5] = 8.55e+07;
  PlogB[0] = 2.6;
  PlogB[1] = 2.6;
  PlogB[2] = 2.56;
  PlogB[3] = 2.19;
  PlogB[4] = 1.43;
  PlogB[5] = 0.75;
  PlogE[0] = 1.72423e+07;
  PlogE[1] = 1.72757e+07;
  PlogE[2] = 1.7733e+07;
  PlogE[3] = 2.19894e+07;
  PlogE[4] = 3.27557e+07;
  PlogE[5] = 4.80775e+07;
  np = 6;
  k[r269f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 6.83969e+06;
  PlogA[1] = 7.03699e+06;
  PlogA[2] = 9.99647e+06;
  PlogA[3] = 2.09794e+08;
  PlogA[4] = 1.27587e+11;
  PlogA[5] = 5.62302e+13;
  PlogB[0] = 1.10902;
  PlogB[1] = 1.10902;
  PlogB[2] = 1.06902;
  PlogB[3] = 0.699022;
  PlogB[4] = -0.0609783;
  PlogB[5] = -0.740978;
  PlogE[0] = 1.48337e+07;
  PlogE[1] = 1.48672e+07;
  PlogE[2] = 1.53245e+07;
  PlogE[3] = 1.95808e+07;
  PlogE[4] = 3.03471e+07;
  PlogE[5] = 4.56689e+07;
  np = 6;
  k[r269b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r270f] = 0x1.fa66acp+28 *
             exp(0x1.ae147ae147ae1p-3 * lgt - 0x1.3ddb72p+27 * rt_inv);
  k[r270b] = 0x1.588c9390f1d35p+19 *
             exp(0x1.5339102bd390cp-1 * lgt - 0x1.0b89fecca4a7cp+24 * rt_inv);
  k[r271f] = 0x1.d4cp+10 *
             exp(0x1.e666666666666p+0 * lgt + 0x1.b73d000000001p+21 * rt_inv);
  k[r271b] = 0x1.ea13d0a663371p-3 *
             exp(0x1.57115380912ecp+1 * lgt - 0x1.098d1b843b2e9p+26 * rt_inv);
  k[r272f] =
      0x1.967e8p+21 * exp(0x1.199999999999ap+0 * lgt - 0x1.140e6p+21 * rt_inv);
  k[r272b] = 0x1.e980080a0bb6ep+12 *
             exp(0x1.cfd9cb1e3356p+0 * lgt - 0x1.12d9b828b12b2p+27 * rt_inv);
  k[r273f] = 0x1.651f128b240f9p-36 *
             exp(0x1.799999999999ap+2 * lgt - 0x1.0ca68p+22 * rt_inv);
  k[r273b] = 0x1.18e43a318151bp-37 *
             exp(0x1.8c9a2791a2182p+2 * lgt - 0x1.604d94cbc254dp+26 * rt_inv);
  k[r274f] =
      0x1.b333333333333p+1 * exp(0x1.4p+1 * lgt - 0x1.1ccd78p+25 * rt_inv);
  k[r274b] = 0x1.4d1d27f770d29p+2 *
             exp(0x1.1211eef7b47aep+1 * lgt - 0x1.2cd99d5849003p+25 * rt_inv);
  k[r275f] = 0x1.7ae147ae147aep+0 *
             exp(0x1.89db22d0e5604p+1 * lgt - 0x1.ccf1e00000001p+24 * rt_inv);
  k[r275b] = 0x1.515851f6bd5efp-12 *
             exp(0x1.f1fda105285d8p+1 * lgt - 0x1.a0d588e7a27bbp+26 * rt_inv);
  k[r276f] =
      0x1.2ap+7 * exp(0x1.ab851eb851eb8p+0 * lgt - 0x1.b2c4fp+24 * rt_inv);
  k[r276b] = 0x1.d53d2536b071bp+1 *
             exp(0x1.18d4d049f71acp+1 * lgt - 0x1.13a36f797738ap+26 * rt_inv);
  PlogP[0] = 10000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+07;
  PlogA[0] = 7.42e+46;
  PlogA[1] = 4.42e+42;
  PlogA[2] = 2.9e+27;
  PlogB[0] = -10.56;
  PlogB[1] = -9.09;
  PlogB[2] = -4.35;
  PlogE[0] = 2.82085e+08;
  PlogE[1] = 2.80618e+08;
  PlogE[2] = 2.57788e+08;
  np = 3;
  k[r277f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 10000;
  PlogP[1] = 100000;
  PlogP[2] = 1e+07;
  PlogA[0] = 1.82558e+45;
  PlogA[1] = 1.08748e+41;
  PlogA[2] = 7.13503e+25;
  PlogB[0] = -10.036;
  PlogB[1] = -8.56601;
  PlogB[2] = -3.82601;
  PlogE[0] = 3.25849e+08;
  PlogE[1] = 3.24381e+08;
  PlogE[2] = 3.01552e+08;
  np = 3;
  k[r277b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r278f] = 0x1.b1d92b7fe08afp-8 *
             exp(0x1.d99999999999ap+1 * lgt - 0x1.2f40d00000001p+25 * rt_inv);
  k[r278b] = 0x1.1cd3c26ff2cb5p-9 *
             exp(0x1.e616c438be4fbp+1 * lgt - 0x1.132dbef0d65a2p+24 * rt_inv);
  k[r279f] = 0x1.3a6a1ccp+35 * exp(-0x1.cbda314cccccdp+27 * rt_inv);
  k[r279b] = 0x1.651dcfc60db85p+25 *
             exp(0x1.0348585fa31p-2 * lgt + 0x1.0de09b0b27d29p+24 * rt_inv);
  k[r280f] = 0x1.c9c38p+26 * exp(-0x1.aef05p+24 * rt_inv);
  k[r280b] = 0x1.75f8b73bb5ebep+13 *
             exp(0x1.e53dced7fd6d8p-1 * lgt - 0x1.39a6f9328992bp+21 * rt_inv);
  k[r281f] = 0x1.197a24894c447p-7 *
             exp(0x1.e083126e978d5p+1 * lgt - 0x1.b10ba80000001p+26 * rt_inv);
  k[r281b] = 0x1.67aa2380ca222p-7 *
             exp(0x1.99111015600c4p+1 * lgt - 0x1.6b9959c97c6a5p+22 * rt_inv);
  k[r282f] = 0x1.197a24894c447p-7 *
             exp(0x1.e083126e978d5p+1 * lgt - 0x1.b10ba80000001p+26 * rt_inv);
  k[r282b] = 0x1.8a84ff3e09c7fp-16 *
             exp(0x1.0e26e6d7018c1p+2 * lgt - 0x1.880974059bb62p+19 * rt_inv);
  k[r283f] = 0x1.50c4288p+33 * exp(-0x1.e5af140000001p+26 * rt_inv);
  k[r283b] = 0x1.ff14ddda35339p+31 *
             exp(-0x1.fc8f02e83dfp-4 * lgt - 0x1.44aa63b0d6e25p+25 * rt_inv);
  k[r284f] = 0x1.502b92p+31 * exp(-0x1.111664p+26 * rt_inv);
  k[r284b] = 0x1.058434cf9d73ap+42 *
             exp(-0x1.a0271f9ebb94p-1 * lgt - 0x1.472439764d53bp+27 * rt_inv);
  k[r285f] = 0x1.502b92p+31 * exp(-0x1.111664p+26 * rt_inv);
  k[r285b] = 0x1.4dba059b81378p+33 *
             exp(-0x1.8aa79579d5p-8 * lgt - 0x1.3707fff04a709p+27 * rt_inv);
  k[r286f] = 0x1.09d633p+29 * exp(-0x1.125d44p+26 * rt_inv);
  k[r286b] = 0x1.afb706411febp+37 *
             exp(-0x1.e9ecbbd0d7bcp-2 * lgt - 0x1.2588d1b4c352ap+27 * rt_inv);
  k[r287f] = 0x1.bf08ebp+35;
  k[r287b] = 0x1.c19f737164761p+51 *
             exp(-0x1.f66753094334p-1 * lgt - 0x1.e4348a4a6af9ap+27 * rt_inv);
  k[r288f] = 0x1.2a05f2p+34;
  k[r288b] = 0x1.950bbf7219093p+56 *
             exp(-0x1.46bbb8fd3506p+0 * lgt - 0x1.296ca15abae63p+28 * rt_inv);
  kTroe0 = 0x1.a3f16bf245702p+85 *
           exp(-0x1.2a7ef9db22d0ep+2 * lgt - 0x1.e2a6cp+23 * rt_inv);
  kTroeInf = 0x1.04ecep+24 *
             exp(0x1.44189374bc6a8p+0 * lgt - 0x1.59e6700000001p+23 * rt_inv);
  fcTroe = 0x1.b22d0e5604188p-3 * exp(-temp / -0x1.3ecp+13) +
           0x1.9374bc6a7ef9ep-1 * exp(-temp / 0x1.4484bfeebc2ap-100);
  k[r289f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM36]);
  kTroe0 = 0x1.1208e1bc7d9adp+92 *
           exp(-0x1.1625b930a4d34p+2 * lgt - 0x1.3b5a391a6fc9bp+27 * rt_inv);
  kTroeInf = 0x1.5488dfa10f2f8p+30 *
             exp(0x1.957d961eb461p+0 * lgt - 0x1.32ce341a6fc9ap+27 * rt_inv);
  k[r289b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM36]);
  k[r290f] = 0x1.193dcceea82b9p+87 *
             exp(-0x1.53f7ced916873p+2 * lgt - 0x1.9f2d0a6666667p+24 * rt_inv);
  k[r290b] = 0x1.0787ea5523138p+83 *
             exp(-0x1.3cc39c83ffb98p+2 * lgt - 0x1.742177ad268dcp+28 * rt_inv);
  k[r291f] = 0x1.45f680bp+39 *
             exp(-0x1.38d4fdf3b645ap-1 * lgt - 0x1.4ff759999999ap+24 * rt_inv);
  k[r291b] = 0x1.02ea35f17a6ap+50 *
             exp(-0x1.789ab78d97a7p+0 * lgt - 0x1.c01cab2bf4c2bp+25 * rt_inv);
  k[r292] = 0x1.2e1906cfp+42 *
            exp(-0x1.428f5c28f5c29p+0 * lgt - 0x1.a6f89cccccccdp+23 * rt_inv);
  k[r293f] = 0x1.75d72p+28;
  k[r293b] = 0x1.dac7cfa58ece5p+40 *
             exp(-0x1.f3a19c7bdf28p-2 * lgt - 0x1.16f0a3c62a8c2p+28 * rt_inv);
  k[r294f] = 0x1.4f46b04p+36;
  k[r294b] = 0x1.e20373d0b2088p+37 *
             exp(0x1.d6975a5af08p-6 * lgt - 0x1.0e761ccd2295ap+28 * rt_inv);
  k[r295f] = 0x1.bf08ebp+35;
  k[r295b] = 0x1.6fcc2e805e8e2p+28 *
             exp(0x1.20ffc3feec5ep-1 * lgt - 0x1.6bd6d3a7ba602p+26 * rt_inv);
  k[r296f] = 0x1.c0ac88ep+34;
  k[r296b] = 0x1.b4329af531687p+39 *
             exp(-0x1.2aa6d7a8036p-4 * lgt - 0x1.4a5497e7928ep+28 * rt_inv);
  k[r297f] = 0x1.c9c38p+29;
  k[r297b] = 0x1.bac42b10d311fp+43 *
             exp(-0x1.2bc578ba81ep-1 * lgt - 0x1.2ba5e1d71d2ap+28 * rt_inv);
  kTroe0 = 0x1.83bdac6ae9bc2p+91 *
           exp(-0x1.3333333333333p+2 * lgt - 0x1.e5348p+22 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+46;
  fcTroe = 0x1.6a7ef9db22d0ep-2 * exp(-temp / 0x1.08p+7) +
           0x1.4ac083126e979p-1 * exp(-temp / 0x1.48cp+10) +
           0x1p+0 * exp(-0x1.5bep+12 / temp);
  k[r298f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM37]);
  kTroe0 = 0x1.5cc307f4e24a5p+111 *
           exp(-0x1.45bf769afe5p+2 * lgt - 0x1.0e25716ac6303p+29 * rt_inv);
  kTroeInf = 0x1.4739cc06e2061p+66 *
             exp(-0x1.28c4367cb1cdp-2 * lgt - 0x1.0a5b086ac6302p+29 * rt_inv);
  k[r298b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM37]);
  k[r299f] = 0x1.74876e8p+35;
  k[r299b] = 0x1.fd5f50c46aeedp+24 *
             exp(0x1.ede1522cbc7p-1 * lgt - 0x1.3a0814586b893p+28 * rt_inv);
  k[r300f] = 0x1.2a05f2p+34;
  k[r300b] = 0x1.f6b44ddcad939p+49 *
             exp(-0x1.dbefba45a416p-1 * lgt - 0x1.956b060e70e53p+27 * rt_inv);
  k[r301f] = 0x1.74876e8p+35 * exp(-0x1.7f0e800000001p+22 * rt_inv);
  k[r301b] = 0x1.40d73f329024ep+28 *
             exp(0x1.6c8fc48a3acp-1 * lgt - 0x1.31e159d9a023ep+29 * rt_inv);
  k[r302f] = 0x1.eap+8 * exp(0x1.4p+1 * lgt - 0x1.1e04000000001p+21 * rt_inv);
  k[r302b] = 0x1.d5cc702497b8ep+20 *
             exp(0x1.dd0f924d6fddp+0 * lgt - 0x1.e79055ecc79ep+26 * rt_inv);
  kTroe0 = 0x1.1d37b09ap+41 *
           exp(-0x1.47ae147ae147bp-1 * lgt - 0x1.8c9f8c0000001p+27 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+49 *
             exp(-0x1.0a3d70a3d70a4p-1 * lgt - 0x1.9500aa0000001p+27 * rt_inv);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r303f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM38]);
  kTroe0 = 0x1.46734baff003ep+32 *
           exp(-0x1.ab185a763a6ep-4 * lgt - 0x1.2cf5e1ccc0b3ep+24 * rt_inv);
  kTroeInf = 0x1.a0640a55f60f7p+40 *
             exp(0x1.01b311085df6p-6 * lgt - 0x1.6ffed1ccc0b42p+24 * rt_inv);
  k[r303b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM38]);
  k[r304f] =
      0x1.69158p+19 * exp(0x1.47ae147ae147bp+0 * lgt - 0x1.3ba38p+23 * rt_inv);
  k[r304b] = 0x1.845a15977e264p+0 *
             exp(0x1.507dfcadfe8fp+1 * lgt - 0x1.81b87d1c1feb5p+27 * rt_inv);
  k[r305f] =
      0x1.69158p+21 * exp(0x1.47ae147ae147bp+0 * lgt - 0x1.3ba38p+23 * rt_inv);
  k[r305b] = 0x1.74f8c2d9d75dbp+24 *
             exp(0x1.e83bd279aba6p-1 * lgt - 0x1.5d5838ccb2eabp+26 * rt_inv);
  k[r306f] = 0x1.48fffffffffffp+11 *
             exp(0x1.11eb851eb851fp+1 * lgt - 0x1.104a180000001p+26 * rt_inv);
  k[r306b] = 0x1.d01820ec633bcp+2 *
             exp(0x1.566156860ad9p+1 * lgt - 0x1.091e73e7c287bp+23 * rt_inv);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 1.578;
  PlogA[1] = 15.18;
  PlogA[2] = 301.7;
  PlogA[3] = 7528;
  PlogA[4] = 5101;
  PlogA[5] = 14.57;
  PlogB[0] = 2.56;
  PlogB[1] = 2.28;
  PlogB[2] = 1.92;
  PlogB[3] = 1.55;
  PlogB[4] = 1.65;
  PlogB[5] = 2.45;
  PlogE[0] = -3.53339e+06;
  PlogE[1] = -1.22215e+06;
  PlogE[2] = 2.50245e+06;
  PlogE[3] = 8.8115e+06;
  PlogE[4] = 1.42256e+07;
  PlogE[5] = 1.87318e+07;
  np = 6;
  k[r307f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 6202.14;
  PlogA[1] = 59663.2;
  PlogA[2] = 1.1858e+06;
  PlogA[3] = 2.95879e+07;
  PlogA[4] = 2.00489e+07;
  PlogA[5] = 57265.6;
  PlogB[0] = 1.90875;
  PlogB[1] = 1.62875;
  PlogB[2] = 1.26875;
  PlogB[3] = 0.898749;
  PlogB[4] = 0.998749;
  PlogB[5] = 1.79875;
  PlogE[0] = 9.53728e+07;
  PlogE[1] = 9.76841e+07;
  PlogE[2] = 1.01409e+08;
  PlogE[3] = 1.07718e+08;
  PlogE[4] = 1.13132e+08;
  PlogE[5] = 1.17638e+08;
  np = 6;
  k[r307b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 475.7;
  PlogA[1] = 4372;
  PlogA[2] = 76480;
  PlogA[3] = 1.277e+06;
  PlogA[4] = 431200;
  PlogA[5] = 825;
  PlogB[0] = 1.68;
  PlogB[1] = 1.4;
  PlogB[2] = 1.05;
  PlogB[3] = 0.73;
  PlogB[4] = 0.92;
  PlogB[5] = 1.77;
  PlogE[0] = -1.37988e+06;
  PlogE[1] = 947676;
  PlogE[2] = 4.66516e+06;
  PlogE[3] = 1.07905e+07;
  PlogE[4] = 1.56314e+07;
  PlogE[5] = 1.96522e+07;
  np = 6;
  k[r308f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 2500;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 0.873911;
  PlogA[1] = 8.03183;
  PlogA[2] = 140.502;
  PlogA[3] = 2345.98;
  PlogA[4] = 792.16;
  PlogA[5] = 1.51561;
  PlogB[0] = 2.51913;
  PlogB[1] = 2.23913;
  PlogB[2] = 1.88913;
  PlogB[3] = 1.56913;
  PlogB[4] = 1.75913;
  PlogB[5] = 2.60913;
  PlogE[0] = 2.26534e+08;
  PlogE[1] = 2.28861e+08;
  PlogE[2] = 2.32579e+08;
  PlogE[3] = 2.38704e+08;
  PlogE[4] = 2.43545e+08;
  PlogE[5] = 2.47566e+08;
  np = 6;
  k[r308b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r309f] = 0x1.388p+13 * exp(0x1p+1 * lgt - 0x1.7f0e800000001p+24 * rt_inv);
  k[r309b] = 0x1.45bab161fa928p+6 *
             exp(0x1.51dc915d5da9cp+1 * lgt - 0x1.9eb282b02073ep+26 * rt_inv);
  k[r310f] = 0x1.74876e8p+36;
  k[r310b] = 0x1.4579e564aa681p+45 *
             exp(-0x1.124b092c1506p-1 * lgt - 0x1.6700cfc667fcep+27 * rt_inv);
  k[r311f] = 0x1.2a05f2p+34;
  k[r311b] = 0x1.f3b47a4d0e6f8p+54 *
             exp(-0x1.2fdde3c52864p+0 * lgt - 0x1.11d36e8149a1ep+28 * rt_inv);
  k[r312f] = 0x1.2a05f2p+33;
  k[r312b] = 0x1.382b3fe1bed27p+28 *
             exp(0x1.27fa59ebdffp-3 * lgt - 0x1.3b6cfea217ae5p+28 * rt_inv);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 3.41e+59;
  PlogA[1] = 2.62e+57;
  PlogA[2] = 1.65e+52;
  PlogA[3] = 5.23e+43;
  PlogA[4] = 4.59e+32;
  PlogA[5] = 3.84e+20;
  PlogB[0] = -14.2;
  PlogB[1] = -13.3;
  PlogB[2] = -11.5;
  PlogB[3] = -8.9;
  PlogB[4] = -5.6;
  PlogB[5] = -2.06;
  PlogE[0] = 3.50086e+08;
  PlogE[1] = 3.56737e+08;
  PlogE[2] = 3.54576e+08;
  PlogE[3] = 3.41024e+08;
  PlogE[4] = 3.18245e+08;
  PlogE[5] = 2.90644e+08;
  np = 6;
  k[r313f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 1.28657e+49;
  PlogA[1] = 9.88509e+46;
  PlogA[2] = 6.22535e+41;
  PlogA[3] = 1.97325e+33;
  PlogA[4] = 1.73178e+22;
  PlogA[5] = 1.44881e+10;
  PlogB[0] = -12.4324;
  PlogB[1] = -11.5324;
  PlogB[2] = -9.7324;
  PlogB[3] = -7.1324;
  PlogB[4] = -3.8324;
  PlogB[5] = -0.292398;
  PlogE[0] = 2.99079e+08;
  PlogE[1] = 3.0573e+08;
  PlogE[2] = 3.03568e+08;
  PlogE[3] = 2.90017e+08;
  PlogE[4] = 2.67238e+08;
  PlogE[5] = 2.39636e+08;
  np = 6;
  k[r313b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 1.2e+54;
  PlogA[1] = 5.18e+59;
  PlogA[2] = 1.62e+66;
  PlogA[3] = 5.55e+64;
  PlogA[4] = 1.55e+58;
  PlogA[5] = 1.78e+47;
  PlogB[0] = -12.9;
  PlogB[1] = -14;
  PlogB[2] = -15.3;
  PlogB[3] = -14.5;
  PlogB[4] = -12.3;
  PlogB[5] = -8.96;
  PlogE[0] = 4.18424e+08;
  PlogE[1] = 4.18008e+08;
  PlogE[2] = 4.40954e+08;
  PlogE[3] = 4.4427e+08;
  PlogE[4] = 4.42533e+08;
  PlogE[5] = 4.2283e+08;
  np = 6;
  k[r314f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 1.55657e+42;
  PlogA[1] = 6.71919e+47;
  PlogA[2] = 2.10137e+54;
  PlogA[3] = 7.19913e+52;
  PlogA[4] = 2.01057e+46;
  PlogA[5] = 2.30891e+35;
  PlogB[0] = -11.2618;
  PlogB[1] = -12.3618;
  PlogB[2] = -13.6618;
  PlogB[3] = -12.8618;
  PlogB[4] = -10.6618;
  PlogB[5] = -7.32176;
  PlogE[0] = 4.64817e+07;
  PlogE[1] = 4.60662e+07;
  PlogE[2] = 6.90117e+07;
  PlogE[3] = 7.23275e+07;
  PlogE[4] = 7.05911e+07;
  PlogE[5] = 5.08879e+07;
  np = 6;
  k[r314b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 8.1e+46;
  PlogA[1] = 1.86e+56;
  PlogA[2] = 4.65e+63;
  PlogA[3] = 4.46e+65;
  PlogA[4] = 2.79e+61;
  PlogA[5] = 6.17e+51;
  PlogB[0] = -11.3;
  PlogB[1] = -13.5;
  PlogB[2] = -15;
  PlogB[3] = -14.9;
  PlogB[4] = -13.4;
  PlogB[5] = -10.3;
  PlogE[0] = 4.64647e+08;
  PlogE[1] = 4.48685e+08;
  PlogE[2] = 4.58662e+08;
  PlogE[3] = 4.70051e+08;
  PlogE[4] = 4.73128e+08;
  PlogE[5] = 4.59992e+08;
  np = 6;
  k[r315f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 1e+07;
  PlogA[0] = 2.48131e+36;
  PlogA[1] = 5.69782e+45;
  PlogA[2] = 1.42445e+53;
  PlogA[3] = 1.36625e+55;
  PlogA[4] = 8.54673e+50;
  PlogA[5] = 1.89008e+41;
  PlogB[0] = -9.91826;
  PlogB[1] = -12.1183;
  PlogB[2] = -13.6183;
  PlogB[3] = -13.5183;
  PlogB[4] = -12.0183;
  PlogB[5] = -8.91826;
  PlogE[0] = 6.6798e+07;
  PlogE[1] = 5.08361e+07;
  PlogE[2] = 6.08124e+07;
  PlogE[3] = 7.22021e+07;
  PlogE[4] = 7.52782e+07;
  PlogE[5] = 6.21425e+07;
  np = 6;
  k[r315b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r316f] = 0x1.bf08ebp+33 * exp(-0x1.9036e2p+27 * rt_inv);
  k[r316b] = 0x1.fc358558b5834p+22 *
             exp(0x1.d37011cab9a1cp-2 * lgt - 0x1.9e045aecbb0adp+23 * rt_inv);
  k[r317f] = 0x1.5f99999999999p+6 *
             exp(0x1.570a3d70a3d71p+1 * lgt - 0x1.7390ap+23 * rt_inv);
  k[r317b] = 0x1.0589e85adf731p-7 *
             exp(0x1.bfadc1815d3b6p+1 * lgt - 0x1.5064a12421523p+25 * rt_inv);
  k[r318f] = 0x1.e3d70a3d70a3dp-1 *
             exp(0x1.91eb851eb851fp+1 * lgt - 0x1.15c04d3333334p+25 * rt_inv);
  k[r318b] = 0x1.1f78e3773ba6p-5 *
             exp(0x1.a6846e5aceebcp+1 * lgt - 0x1.d0399c6fffb32p+24 * rt_inv);
  k[r319f] = 0x1.1eccccccccccdp+6 *
             exp(0x1.451eb851eb852p+1 * lgt + 0x1.87bd4p+22 * rt_inv);
  k[r319b] = 0x1.2087e1b502932p-4 *
             exp(0x1.a0bfd6f0aed54p+1 * lgt - 0x1.50be54fbd0861p+26 * rt_inv);
  k[r320f] = 0x1.85e70a3ac171fp-18 *
             exp(0x1.11eb851eb851fp+2 * lgt + 0x1.c68f8p+23 * rt_inv);
  k[r320b] = 0x1.394fb10722d8ap-19 *
             exp(0x1.15b6c703c87cp+2 * lgt - 0x1.3fac79d84c543p+25 * rt_inv);
  k[r321f] = 0x1.2ca5d05ea7ab2p-25 *
             exp(0x1.50a3d70a3d70ap+2 * lgt - 0x1.dd3b2a6666667p+24 * rt_inv);
  k[r321b] = 0x1.f84b907c04ecfp-31 *
             exp(0x1.5e3e8512dab8bp+2 * lgt + 0x1.e00b7811018cdp+21 * rt_inv);
  k[r322f] = 0x1.63b127abd8548p-31 *
             exp(0x1.5333333333333p+2 * lgt - 0x1.503b2d3333334p+25 * rt_inv);
  k[r322b] = 0x1.dc8cad8d2b08cp-28 *
             exp(0x1.36c893d17f0edp+2 * lgt + 0x1.bb0e2f40c994ap+24 * rt_inv);
  k[r323f] = 0x1.0666666666666p+3 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.5727a8p+25 * rt_inv);
  k[r323b] = 0x1.4f97aad84bd5dp+2 *
             exp(0x1.18f95b3a9172p+1 * lgt + 0x1.f76d4b6a736edp+15 * rt_inv);
  k[r324f] = 0x1.2a05f2p+31 * exp(-0x1.7f0e800000001p+26 * rt_inv);
  k[r324b] = 0x1.307150c83ca4ep+39 *
             exp(-0x1.02ef4c00ef13dp+0 * lgt - 0x1.0aa7091ba181ep+24 * rt_inv);
  k[r325f] =
      0x1.22p+7 * exp(0x1.3c28f5c28f5c3p+1 * lgt - 0x1.bf69p+21 * rt_inv);
  k[r325b] = 0x1.fa99cf72e9e5dp-8 *
             exp(0x1.a0881bec63e4bp+1 * lgt - 0x1.c5a5dcbaa583dp+24 * rt_inv);
  k[r326f] = 0x1.87ea6f9ff5ff2p-20 *
             exp(0x1.2eb851eb851ecp+2 * lgt - 0x1.b906ap+22 * rt_inv);
  k[r326b] = 0x1.116cb52c2b12ap-25 *
             exp(0x1.36e297961dec6p+2 * lgt + 0x1.19362e100f835p+22 * rt_inv);
  k[r327f] = 0x1.46887a8d64d7ep-6 *
             exp(0x1.af5c28f5c28f6p+1 * lgt - 0x1.e7603p+24 * rt_inv);
  k[r327b] = 0x1.ad22b40b6593ap-9 *
             exp(0x1.d5de4ac24a3e8p+1 * lgt - 0x1.1582647630896p+26 * rt_inv);
  k[r328f] = 0x1.0abb44e50c5ebp-9 *
             exp(0x1.c8f5c28f5c28fp+1 * lgt - 0x1.ecee180000001p+24 * rt_inv);
  k[r328b] = 0x1.17fe7275eae82p-3 *
             exp(0x1.9b6d498740e6cp+1 * lgt - 0x1.0ca7c4cd0c5fcp+25 * rt_inv);
  k[r329f] = 0x1.7d784p+25 * exp(-0x1.4bfb800000001p+25 * rt_inv);
  k[r329b] = 0x1.2571467daddbp+24 *
             exp(0x1.2aaa8f1e47f6p-2 * lgt - 0x1.f92d34f2e48fcp+25 * rt_inv);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 5.69e+52;
  PlogA[1] = 3.29e+56;
  PlogA[2] = 8.58e+57;
  PlogA[3] = 5.36e+55;
  PlogA[4] = 1.66e+48;
  PlogA[5] = 8.26e+44;
  PlogA[6] = 1.01e+40;
  PlogA[7] = 1.1e+36;
  PlogB[0] = -13.38;
  PlogB[1] = -14.12;
  PlogB[2] = -14.16;
  PlogB[3] = -13.15;
  PlogB[4] = -10.64;
  PlogB[5] = -9.59;
  PlogB[6] = -8.06;
  PlogB[7] = -6.84;
  PlogE[0] = 1.88485e+08;
  PlogE[1] = 2.01372e+08;
  PlogE[2] = 2.12309e+08;
  PlogE[3] = 2.17091e+08;
  PlogE[4] = 2.10443e+08;
  PlogE[5] = 2.05928e+08;
  PlogE[6] = 1.98485e+08;
  PlogE[7] = 1.92041e+08;
  np = 8;
  k[r330f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 3.53573e+49;
  PlogA[1] = 2.04439e+53;
  PlogA[2] = 5.33156e+54;
  PlogA[3] = 3.33067e+52;
  PlogA[4] = 1.03151e+45;
  PlogA[5] = 5.13272e+41;
  PlogA[6] = 6.27608e+36;
  PlogA[7] = 6.83534e+32;
  PlogB[0] = -13.2952;
  PlogB[1] = -14.0352;
  PlogB[2] = -14.0752;
  PlogB[3] = -13.0652;
  PlogB[4] = -10.5552;
  PlogB[5] = -9.50524;
  PlogB[6] = -7.97524;
  PlogB[7] = -6.75524;
  PlogE[0] = 8.41394e+07;
  PlogE[1] = 9.70261e+07;
  PlogE[2] = 1.07963e+08;
  PlogE[3] = 1.12745e+08;
  PlogE[4] = 1.06097e+08;
  PlogE[5] = 1.01582e+08;
  PlogE[6] = 9.41391e+07;
  PlogE[7] = 8.76958e+07;
  np = 8;
  k[r330b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 5.4e+46;
  PlogA[1] = 1.21e+51;
  PlogA[2] = 2.87e+54;
  PlogA[3] = 3.79e+53;
  PlogA[4] = 6.33e+46;
  PlogA[5] = 3.87e+43;
  PlogA[6] = 5.08e+38;
  PlogA[7] = 5.12e+34;
  PlogB[0] = -11.63;
  PlogB[1] = -12.55;
  PlogB[2] = -13.15;
  PlogB[3] = -12.51;
  PlogB[4] = -10.2;
  PlogB[5] = -9.17;
  PlogB[6] = -7.65;
  PlogB[7] = -6.41;
  PlogE[0] = 1.85447e+08;
  PlogE[1] = 1.97652e+08;
  PlogE[2] = 2.12137e+08;
  PlogE[3] = 2.19911e+08;
  PlogE[4] = 2.15229e+08;
  PlogE[5] = 2.11041e+08;
  PlogE[6] = 2.03815e+08;
  PlogE[7] = 1.97409e+08;
  np = 8;
  k[r331f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 1.36384e+45;
  PlogA[1] = 3.05601e+49;
  PlogA[2] = 7.24856e+52;
  PlogA[3] = 9.57214e+51;
  PlogA[4] = 1.59872e+45;
  PlogA[5] = 9.77419e+41;
  PlogA[6] = 1.28302e+37;
  PlogA[7] = 1.29312e+33;
  PlogB[0] = -12.0692;
  PlogB[1] = -12.9892;
  PlogB[2] = -13.5892;
  PlogB[3] = -12.9492;
  PlogB[4] = -10.6392;
  PlogB[5] = -9.60923;
  PlogB[6] = -8.08923;
  PlogB[7] = -6.84923;
  PlogE[0] = 3.73379e+07;
  PlogE[1] = 4.95426e+07;
  PlogE[2] = 6.40276e+07;
  PlogE[3] = 7.18015e+07;
  PlogE[4] = 6.71196e+07;
  PlogE[5] = 6.29314e+07;
  PlogE[6] = 5.57056e+07;
  PlogE[7] = 4.92999e+07;
  np = 8;
  k[r331b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 5.48e+45;
  PlogA[1] = 2.54e+49;
  PlogA[2] = 1.65e+54;
  PlogA[3] = 1.81e+55;
  PlogA[4] = 4.58e+49;
  PlogA[5] = 4.11e+46;
  PlogA[6] = 6.68e+41;
  PlogA[7] = 6.54e+37;
  PlogB[0] = -11.63;
  PlogB[1] = -12.37;
  PlogB[2] = -13.4;
  PlogB[3] = -13.31;
  PlogB[4] = -11.32;
  PlogB[5] = -10.33;
  PlogB[6] = -8.83;
  PlogB[7] = -7.58;
  PlogE[0] = 1.85468e+08;
  PlogE[1] = 1.94326e+08;
  PlogE[2] = 2.10581e+08;
  PlogE[3] = 2.22304e+08;
  PlogE[4] = 2.20555e+08;
  PlogE[5] = 2.16873e+08;
  PlogE[6] = 2.10045e+08;
  PlogE[7] = 2.03748e+08;
  np = 8;
  k[r332f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 100;
  PlogP[1] = 1000;
  PlogP[2] = 10000;
  PlogP[3] = 100000;
  PlogP[4] = 1e+06;
  PlogP[5] = 2e+06;
  PlogP[6] = 5e+06;
  PlogP[7] = 1e+07;
  PlogA[0] = 2.24108e+48;
  PlogA[1] = 1.03875e+52;
  PlogA[2] = 6.74778e+56;
  PlogA[3] = 7.40211e+57;
  PlogA[4] = 1.87302e+52;
  PlogA[5] = 1.68081e+49;
  PlogA[6] = 2.73183e+44;
  PlogA[7] = 2.67457e+40;
  PlogB[0] = -12.2866;
  PlogB[1] = -13.0266;
  PlogB[2] = -14.0566;
  PlogB[3] = -13.9666;
  PlogB[4] = -11.9766;
  PlogB[5] = -10.9866;
  PlogB[6] = -9.48657;
  PlogB[7] = -8.23657;
  PlogE[0] = 1.4757e+08;
  PlogE[1] = 1.56428e+08;
  PlogE[2] = 1.72682e+08;
  PlogE[3] = 1.84406e+08;
  PlogE[4] = 1.82657e+08;
  PlogE[5] = 1.78975e+08;
  PlogE[6] = 1.72147e+08;
  PlogE[7] = 1.6585e+08;
  np = 8;
  k[r332b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 5.26e+14;
  PlogA[1] = 5.26e+14;
  PlogA[2] = 5.28e+14;
  PlogA[3] = 1.54e+15;
  PlogA[4] = 3.78e+17;
  PlogB[0] = -1.637;
  PlogB[1] = -1.637;
  PlogB[2] = -1.638;
  PlogB[3] = -1.771;
  PlogB[4] = -2.429;
  PlogE[0] = 3.50619e+06;
  PlogE[1] = 3.50619e+06;
  PlogE[2] = 3.51038e+06;
  PlogE[3] = 4.68608e+06;
  PlogE[4] = 1.29286e+07;
  np = 5;
  k[r333f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 4.79887e+14;
  PlogA[1] = 4.79887e+14;
  PlogA[2] = 4.81711e+14;
  PlogA[3] = 1.40499e+15;
  PlogA[4] = 3.44862e+17;
  PlogB[0] = -1.56657;
  PlogB[1] = -1.56657;
  PlogB[2] = -1.56757;
  PlogB[3] = -1.70057;
  PlogB[4] = -2.35857;
  PlogE[0] = 1.04102e+08;
  PlogE[1] = 1.04102e+08;
  PlogE[2] = 1.04107e+08;
  PlogE[3] = 1.05282e+08;
  PlogE[4] = 1.13525e+08;
  np = 5;
  k[r333b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 0.512;
  PlogA[1] = 0.533;
  PlogA[2] = 0.762;
  PlogA[3] = 8.92;
  PlogA[4] = 438;
  PlogB[0] = 2.496;
  PlogB[1] = 2.49;
  PlogB[2] = 2.446;
  PlogB[3] = 2.146;
  PlogB[4] = 1.699;
  PlogE[0] = -1.73218e+06;
  PlogE[1] = -1.68197e+06;
  PlogE[2] = -1.23846e+06;
  PlogE[3] = 1.96648e+06;
  PlogE[4] = 9.74872e+06;
  np = 5;
  k[r334f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 18.9856;
  PlogA[1] = 19.7643;
  PlogA[2] = 28.256;
  PlogA[3] = 330.765;
  PlogA[4] = 16241.6;
  PlogB[0] = 2.04243;
  PlogB[1] = 2.03643;
  PlogB[2] = 1.99243;
  PlogB[3] = 1.69243;
  PlogB[4] = 1.24543;
  PlogE[0] = 5.51002e+07;
  PlogE[1] = 5.51504e+07;
  PlogE[2] = 5.55939e+07;
  PlogE[3] = 5.87988e+07;
  PlogE[4] = 6.65811e+07;
  np = 5;
  k[r334b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 2.05e+58;
  PlogA[1] = 3.3e+51;
  PlogA[2] = 1.31e+42;
  PlogA[3] = 2.16e+33;
  PlogA[4] = 9.4e+28;
  PlogB[0] = -12.796;
  PlogB[1] = -10.574;
  PlogB[2] = -7.657;
  PlogB[3] = -4.989;
  PlogB[4] = -3.669;
  PlogE[0] = 4.18526e+08;
  PlogE[1] = 4.10958e+08;
  PlogE[2] = 3.9606e+08;
  PlogE[3] = 3.80395e+08;
  PlogE[4] = 3.72471e+08;
  np = 5;
  k[r335f] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  PlogP[0] = 1000;
  PlogP[1] = 10000;
  PlogP[2] = 100000;
  PlogP[3] = 1e+06;
  PlogP[4] = 1e+07;
  PlogA[0] = 1.69752e+46;
  PlogA[1] = 2.73259e+39;
  PlogA[2] = 1.08476e+30;
  PlogA[3] = 1.78861e+21;
  PlogA[4] = 7.78375e+16;
  PlogB[0] = -11.1823;
  PlogB[1] = -8.9603;
  PlogB[2] = -6.0433;
  PlogB[3] = -3.3753;
  PlogB[4] = -2.0553;
  PlogE[0] = 6.13946e+07;
  PlogE[1] = 5.38262e+07;
  PlogE[2] = 3.89286e+07;
  PlogE[3] = 2.32633e+07;
  PlogE[4] = 1.53401e+07;
  np = 5;
  k[r335b] = GetPlogRateCoeff(temp, pressure, lgt, rt_inv, PlogP, PlogA, PlogB,
                              PlogE, np);
  k[r336f] =
      0x1.f4p+6 * exp(0x1.3dd2f1a9fbe77p+1 * lgt - 0x1.c68f8p+20 * rt_inv);
  k[r336b] = 0x1.c12bd7dc6c438p-3 *
             exp(0x1.8699fd711facp+1 * lgt - 0x1.6f8f937e4eac7p+26 * rt_inv);
  k[r337f] =
      0x1.eap+9 * exp(0x1.370a3d70a3d71p+1 * lgt - 0x1.496dcp+24 * rt_inv);
  k[r337b] = 0x1.45788ac195589p-3 *
             exp(0x1.8cd3aea9bdb2ap+1 * lgt - 0x1.96adb2291d9fp+25 * rt_inv);
  k[r338f] = 0x1.e93c24p+28 *
             exp(0x1.b020c49ba5e35p-3 * lgt - 0x1.3830f00000001p+24 * rt_inv);
  k[r338b] = 0x1.7d942f5b99e7bp+15 *
             exp(0x1.b21c7e6fbda9p-1 * lgt - 0x1.616b2f624f108p+25 * rt_inv);
  k[r339f] = 0x1.79a7bp+28 * exp(-0x1.3851ap+25 * rt_inv);
  k[r339b] = 0x1.bb310a1edffa6p+26 *
             exp(0x1.3a80ef4e8284p-3 * lgt - 0x1.370e54f8aeaf1p+26 * rt_inv);
  k[r340f] = 0x1.9de508p+28 * exp(-0x1.9c6ca00000001p+24 * rt_inv);
  k[r340b] = 0x1.2e35ec547d9f9p+15 *
             exp(0x1.00f4b01785808p+0 * lgt - 0x1.e5cd920c1ace8p+25 * rt_inv);
  k[r341f] = 0x1.c1451f6p+35 * exp(-0x1.830bfcp+27 * rt_inv);
  k[r341b] = 0x1.c8138225923e2p+25 *
             exp(0x1.3c9f7b0dbc654p-2 * lgt - 0x1.8a5ebe01585aap+22 * rt_inv);
  k[r342f] = 0x1.faa3b5p+33 * exp(-0x1.468e48p+26 * rt_inv);
  k[r342b] = 0x1.7b6c906dc27a2p+29 *
             exp(0x1.0b6927336ac8p-4 * lgt - 0x1.7ef49c50b90bep+25 * rt_inv);
  k[r343f] = 0x1.faa3b5p+33 * exp(-0x1.468e48p+26 * rt_inv);
  k[r343b] = 0x1.214c70f5ca165p+34 *
             exp(-0x1.011c780dd11ep-1 * lgt - 0x1.33ed65b221b8cp+25 * rt_inv);
  k[r344f] = 0x1.1999999999999p+4 *
             exp(0x1.3d70a3d70a3d7p+1 * lgt - 0x1.875b300000001p+24 * rt_inv);
  k[r344b] = 0x1.f099cbc8441eep+33 *
             exp(0x1.3d18663f04f9p+0 * lgt - 0x1.3560a77e4ddacp+27 * rt_inv);
  k[r345f] = 0x1.0db6054p+34;
  k[r345b] = 0x1.f0d2cc7476f72p+82 *
             exp(-0x1.331f68e08a86p+1 * lgt - 0x1.9a47a06be8c43p+28 * rt_inv);
  k[r346f] = 0x1.8f59e3p+33 * exp(-0x1.a55ccp+23 * rt_inv);
  k[r346b] = 0x1.c6ff8a4e2916p+23 *
             exp(0x1.68155e9c421fp-1 * lgt - 0x1.64d3226a4867ap+26 * rt_inv);
  k[r347f] = 0x1.620d35p+32 * exp(-0x1.dd088p+22 * rt_inv);
  k[r347b] = 0x1.d9aa78af11f09p+21 *
             exp(0x1.5703e700aed1p-1 * lgt - 0x1.37a60506e11f1p+26 * rt_inv);
  k[r348f] = 0x1.20cp+13 * exp(0x1.8p+0 * lgt + 0x1.eb55800000002p+21 * rt_inv);
  k[r348b] = 0x1.bcf0f6cbb11d9p+6 *
             exp(0x1.0d02f2351a7cp+1 * lgt - 0x1.0823656a04219p+27 * rt_inv);
  k[r349f] = 0x1.2b836a4p+33 * exp(-0x1.44ccd40000001p+27 * rt_inv);
  k[r349b] = 0x1.04c4b32ac25f2p+26 *
             exp(0x1.5e7eae7d70b48p-2 * lgt - 0x1.4d4ef4d73c631p+24 * rt_inv);
  k[r350f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r350b] = 0x1.ccefc4a33ed7fp+29 *
             exp(0x1.92e5f4f23a9p-4 * lgt - 0x1.dcb10efc2c343p+25 * rt_inv);
  k[r351f] =
      0x1.46p+11 * exp(0x1.c7ae147ae147bp+0 * lgt - 0x1.795fe8p+24 * rt_inv);
  k[r351b] = 0x1.481cc187c9c09p+12 *
             exp(0x1.f775ff409ed48p+0 * lgt - 0x1.b069a04e68464p+26 * rt_inv);
  k[r352f] = 0x1.9ed92cp+30 * exp(-0x1.0d6aa00000001p+25 * rt_inv);
  k[r352b] = 0x1.3dfecfa73c37cp+33 *
             exp(0x1.6cd958774b2p-4 * lgt - 0x1.15cdf749195cfp+27 * rt_inv);
  k[r353f] = 0x1.dcd65p+29 * exp(-0x1.a55ccp+23 * rt_inv);
  k[r353b] = 0x1.2a9b2d35885acp+19 *
             exp(0x1.096c7cf372cep+0 * lgt - 0x1.77a35a5bc6fddp+26 * rt_inv);
  k[r354f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r354b] = 0x1.5f730289920ecp+34 *
             exp(-0x1.e059bcabee43p-2 * lgt - 0x1.91a9d85d94e07p+25 * rt_inv);
  k[r355f] = 0x1.200278p+27 * exp(-0x1.33157p+24 * rt_inv);
  k[r355b] = 0x1.935e4f2178e53p+71 *
             exp(-0x1.2ea3d14b8d51p+1 * lgt - 0x1.36a2091aa058cp+27 * rt_inv);
  k[r356f] = 0x1.0db6054p+34;
  k[r356b] = 0x1.cb04ae5383a1ep+78 *
             exp(-0x1.d5f1d70523dp+0 * lgt - 0x1.51708ba851515p+28 * rt_inv);
  k[r357f] = 0x1.2a05f2p+35 * exp(-0x1.0c23cp+24 * rt_inv);
  k[r357b] = 0x1.5cc0d8caaa287p+14 *
             exp(0x1.4ddf74587bfc8p+0 * lgt - 0x1.283a704fee634p+26 * rt_inv);
  k[r358f] = 0x1.2a05f2p+32 * exp(-0x1.c91d400000001p+22 * rt_inv);
  k[r358b] = 0x1.9983ad10d26eap+10 *
             exp(0x1.4556b88ab24ap+0 * lgt - 0x1.d6e28dd90e333p+25 * rt_inv);
  k[r359f] =
      0x1.9a762p+24 * exp(0x1.851eb851eb852p-1 * lgt + 0x1.5b4ep+20 * rt_inv);
  k[r359b] = 0x1.44d0b8801ec2bp+7 *
             exp(0x1.f66a059d859ep+0 * lgt - 0x1.cf3e34b9ae3ecp+26 * rt_inv);
  k[r360f] =
      0x1.46p+11 * exp(0x1.c7ae147ae147bp+0 * lgt - 0x1.795fe8p+24 * rt_inv);
  k[r360b] = 0x1.5101362827852p+1 *
             exp(0x1.48a562257cd84p+1 * lgt - 0x1.657396340e41ap+26 * rt_inv);
  k[r361f] = 0x1.4dc938p+31 * exp(-0x1.b2218p+25 * rt_inv);
  k[r361b] = 0x1.b8661b01623bep+18 *
             exp(0x1.660648b2fd3ap-1 * lgt - 0x1.7c65bac7782cbp+25 * rt_inv);
  k[r362f] = 0x1.dcd65p+29 * exp(-0x1.a55ccp+23 * rt_inv);
  k[r362b] = 0x1.32b2ead1ec5cep+8 *
             exp(0x1.a34141fdcdb98p+0 * lgt - 0x1.2cad50416cf75p+26 * rt_inv);
  k[r363f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r363b] = 0x1.68f960c2adeecp+23 *
             exp(0x1.0df2aefafaa2p-3 * lgt - 0x1.f77b8851c1aa5p+24 * rt_inv);
  k[r364f] = 0x1.dcd65p+29 * exp(-0x1.febep+24 * rt_inv);
  k[r364b] = 0x1.2086843858ed8p+21 *
             exp(0x1.8e861ffd360fp-1 * lgt - 0x1.41b99e374ffecp+26 * rt_inv);
  k[r365f] = 0x1.1f5792p+29 * exp(-0x1.a55ccp+23 * rt_inv);
  k[r365b] = 0x1.1af8cf7c0fbcap+13 *
             exp(0x1.24ada1e04f29p+0 * lgt - 0x1.30aed7cd88167p+26 * rt_inv);
  k[r366f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r366b] = 0x1.8bf4b964fe8b9p+14 *
             exp(0x1.2837d110a560ap+0 * lgt - 0x1.a8d57d7f8f722p+24 * rt_inv);
  k[r367f] = 0x1.2b836a4p+33 * exp(-0x1.44ccd40000001p+27 * rt_inv);
  k[r367b] = 0x1.0bd5eba41c489p+15 *
             exp(0x1.e2e8e1536e18ap-1 * lgt - 0x1.0bb6636ea2812p+21 * rt_inv);
  k[r368f] = 0x1.66d1e9p+31 * exp(-0x1.7c80cp+25 * rt_inv);
  k[r368b] = 0x1.acbaaf1a0ea9fp+21 *
             exp(0x1.21b2d4c6949bp-1 * lgt - 0x1.a9ae24a088253p+25 * rt_inv);
  k[r369f] = 0x1.954fc4p+30 * exp(-0x1.0d6aa00000001p+25 * rt_inv);
  k[r369b] = 0x1.3f1aef6535bap+22 *
             exp(0x1.6144b5239f18p-1 * lgt - 0x1.e0a5e477d8b59p+26 * rt_inv);
  k[r370f] = 0x1.200278p+27 * exp(-0x1.33157p+24 * rt_inv);
  k[r370b] = 0x1.7ec6534214335p+56 *
             exp(-0x1.3325e2d0ce85cp+0 * lgt - 0x1.fde36a1911c1fp+25 * rt_inv);
  kTroe0 = 0x1.b5b8bb7p+42 * exp(-0x1.557aa20000001p+27 * rt_inv);
  kTroeInf = 0x1.dae16bd387e5ep+71 *
             exp(-0x1.91eb851eb851fp+0 * lgt - 0x1.4ebc53p+28 * rt_inv);
  fcTroe = 0x1.178d4fdf3b646p-1 * exp(-temp / 0x1.4484bfeebc2ap-100) +
           0x1.d0e5604189375p-2 * exp(-temp / 0x1.39cp+11);
  k[r371f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM39]);
  kTroe0 = 0x1.2dea127d63e5cp+12 *
           exp(0x1.7b7bd87ff382p-1 * lgt + 0x1.4f79fefde4c2p+27 * rt_inv);
  kTroeInf = 0x1.478b5c0e541d7p+41 *
             exp(-0x1.a85b31bd7d21ep-1 * lgt + 0x1.defebf7930848p+21 * rt_inv);
  k[r371b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM39]);
  k[r372f] = 0x1.8b4p+12 * exp(0x1p+1 * lgt + 0x1.4cd9c66666668p+21 * rt_inv);
  k[r372b] = 0x1.20f110e45fddcp+2 *
             exp(0x1.52aeadaa8f1a4p+1 * lgt - 0x1.538b5d96b26b2p+26 * rt_inv);
  k[r373f] =
      0x1.e29p+12 * exp(0x1.0b851eb851eb8p+1 * lgt - 0x1.b0168p+23 * rt_inv);
  k[r373b] = 0x1.04d5c595b07d3p-1 *
             exp(0x1.6b3631d4d713cp+1 * lgt - 0x1.48f61ec04b83ap+25 * rt_inv);
  k[r374f] =
      0x1.7a6bp+19 * exp(0x1.5c28f5c28f5c3p+0 * lgt - 0x1.1f4aep+23 * rt_inv);
  k[r374b] = 0x1.e05b9df626688p+4 *
             exp(0x1.09813016e8074p+1 * lgt - 0x1.f03e37f2f9e34p+24 * rt_inv);
  k[r375f] = 0x1.bbe76c8b43957p-1 *
             exp(0x1.8147ae147ae14p+1 * lgt - 0x1.81edf80000001p+25 * rt_inv);
  k[r375b] = 0x1.0e870350f97eep-6 *
             exp(0x1.938a993181588p+1 * lgt - 0x1.7afec39f9b9f2p+23 * rt_inv);
  k[r376f] = 0x1.3f7ced916872bp-2 *
             exp(0x1.8f5c28f5c28f6p+1 * lgt - 0x1.a50b080000001p+25 * rt_inv);
  k[r376b] = 0x1.28e98c084484bp-3 *
             exp(0x1.58fcacd5b96a2p+1 * lgt - 0x1.b6ac524a7cc17p+22 * rt_inv);
  k[r377f] = 0x1.8d32d5d8a0a8cp-30 *
             exp(0x1.6eb851eb851ecp+2 * lgt - 0x1.6be76p+24 * rt_inv);
  k[r377b] = 0x1.7b50afade657ap-33 *
             exp(0x1.7d802a57aed06p+2 * lgt - 0x1.d6b866888b40cp+25 * rt_inv);
  k[r378f] = 0x1.31794b4p+35 * exp(-0x1.6665ba0000001p+27 * rt_inv);
  k[r378b] = 0x1.f8b24d5706e62p+23 *
             exp(0x1.8bdc8a2916dfp-2 * lgt + 0x1.594e0ea29c649p+23 * rt_inv);
  k[r379f] = 0x1.1f0e54p+29 * exp(-0x1.04187p+24 * rt_inv);
  k[r379b] = 0x1.551fdb9235989p+14 *
             exp(0x1.14c3f3de5c1f8p+0 * lgt - 0x1.849d26a348af4p+25 * rt_inv);
  k[r380f] = 0x1.2a05f2p+32 * exp(-0x1.1a583cp+26 * rt_inv);
  k[r380b] = 0x1.26c8e7ef954f9p+33 *
             exp(-0x1.0695ab29f8a44p-1 * lgt - 0x1.bd5bdfa8f36fbp+24 * rt_inv);
  k[r381f] =
      0x1.62p+5 * exp(0x1.4cccccccccccdp+1 * lgt - 0x1.bc06c8p+25 * rt_inv);
  k[r381b] = 0x1.3defbf7db9642p+8 *
             exp(0x1.007f75dd9b11bp+1 * lgt - 0x1.83e71adfa5bdap+24 * rt_inv);
  k[r382f] = 0x1.2a05f2p+33 * exp(-0x1.1a583cp+26 * rt_inv);
  k[r382b] = 0x1.c770bca12efb1p+22 *
             exp(0x1.bfa8a9a1c09p-1 * lgt - 0x1.01589646e7782p+27 * rt_inv);
  k[r383f] = 0x1.d1a94a2p+43 * exp(-0x1.96ff68p+26 * rt_inv);
  k[r383b] = 0x1.4777be446f68ep+14 *
             exp(0x1.3a92a27d3f348p+0 * lgt - 0x1.f7289847f54fp+25 * rt_inv);
  k[r384f] = 0x1.671e344p+34;
  k[r384b] = 0x1.6e246790168c5p+35 *
             exp(0x1.f352d8f5194p-2 * lgt - 0x1.2b9f8987f10e8p+28 * rt_inv);
  k[r385f] = 0x1.5f5c28f5c28f6p+2 *
             exp(0x1.6666666666666p+1 * lgt - 0x1.763f1p+24 * rt_inv);
  k[r385b] = 0x1.2aed3d196a59cp+8 *
             exp(0x1.330c853e2d096p+1 * lgt - 0x1.bfcf6d69e596ap+25 * rt_inv);
  k[r386f] = 0x1.2c684cp+30 * exp(-0x1.0f4cc40000001p+25 * rt_inv);
  k[r386b] = 0x1.5ca8ce91ae67dp+33 *
             exp(-0x1.636d26cba56p-4 * lgt - 0x1.eecf213ff45cfp+25 * rt_inv);
  k[r387f] = 0x1.dcd65p+30;
  k[r387b] = 0x1.c15284fcfa9f6p+71 *
             exp(-0x1.21abe178f0d7p+1 * lgt - 0x1.2fec30878771p+27 * rt_inv);
  k[r388f] = 0x1.dcd65p+29 * exp(-0x1.74341p+25 * rt_inv);
  k[r388b] = 0x1.91454d55b16b2p+36 *
             exp(-0x1.d3fd2fcadf138p-1 * lgt - 0x1.22e16d3e5f36dp+25 * rt_inv);
  k[r389f] = 0x1.4dc938p+31 * exp(-0x1.b2218p+25 * rt_inv);
  k[r389b] = 0x1.7f30908b5ca1ap+35 *
             exp(-0x1.330350036e77p-1 * lgt - 0x1.3ba155146dfedp+25 * rt_inv);
  k[r390] = 0x1.0c5c99cba81dcp+67 * exp(-0x1.2p+2 * lgt);
  k[r391f] = 0x1.2a05f2p+34;
  k[r391b] = 0x1.2a26b570879b8p+74 *
             exp(-0x1.10a696ff6cadp+1 * lgt - 0x1.6e290c1d77464p+27 * rt_inv);
  k[r392f] = 0x1.7d784p+26 * exp(-0x1.fc304p+24 * rt_inv);
  k[r392b] = 0x1.f038c7eee1ddcp+60 *
             exp(-0x1.8ccc3866a49e8p+0 * lgt - 0x1.5f3cb875d7dcep+26 * rt_inv);
  k[r393f] = 0x1.7d784p+25 * exp(-0x1.febep+20 * rt_inv);
  k[r393b] = 0x1.1b6a761c48acap+19 *
             exp(0x1.5a21c1014de8p-1 * lgt - 0x1.8ccbc9e59ec96p+27 * rt_inv);
  k[r394f] = 0x1.2a05f2p+33 * exp(-0x1.f466500000001p+24 * rt_inv);
  k[r394b] = 0x1.1f91996320d9ep+50 *
             exp(-0x1.6179083813cp-1 * lgt - 0x1.e51e41a6cbcc2p+24 * rt_inv);
  k[r395f] = 0x1.bf08ebp+35 * exp(-0x1.586e88p+26 * rt_inv);
  k[r395b] = 0x1.e637667f0cb24p+41 *
             exp(-0x1.d4fe1b8754cb8p-1 * lgt - 0x1.66e5d971c110fp+25 * rt_inv);
  k[r396] = 0x1.b48eb57ep+43 * exp(-0x1.4b58100000001p+26 * rt_inv);
  k[r397f] = 0x1.4dc938p+29;
  k[r397b] = 0x1.5e8b149413b79p+66 *
             exp(-0x1.e30c198918a2p+0 * lgt - 0x1.2d2afeae86384p+27 * rt_inv);
  k[r398f] = 0x1.2a05f2p+35 * exp(-0x1.288cb8p+26 * rt_inv);
  k[r398b] = 0x1.967979ed1563ep+9 *
             exp(0x1.12b0627529aep+0 * lgt - 0x1.d1143f04b6213p+27 * rt_inv);
  k[r399f] = 0x1.1c37937e08p+54 * exp(-0x1.43343cp+27 * rt_inv);
  k[r399b] = 0x1.28635c156a192p+12 *
             exp(0x1.337d6fae2e6e1p+1 * lgt + 0x1.3371ef1c9b703p+24 * rt_inv);
  k[r400f] = 0x1.dcd65p+26 * exp(-0x1.7bdd5p+25 * rt_inv);
  k[r400b] = 0x1.70bdb42d544ccp+76 *
             exp(-0x1.687ebb0b9c954p+1 * lgt - 0x1.0b52dca5fae9ap+27 * rt_inv);
  k[r401f] = 0x1.74876e8p+36 * exp(-0x1.bee64p+25 * rt_inv);
  k[r401b] = 0x1.d07de3e4f1cebp+24 *
             exp(0x1.51e864f0e25p-1 * lgt - 0x1.2a949b5ae5cb3p+26 * rt_inv);
  k[r402f] = 0x1.1e1a3p+27 * exp(-0x1.3272p+24 * rt_inv);
  k[r402b] = 0x1.77f7e952e77b2p+65 *
             exp(-0x1.107411b06d14p+1 * lgt - 0x1.3e7a2cd70bf49p+26 * rt_inv);
  k[r403f] = 0x1.1e1a3p+27 * exp(-0x1.1d0ed8p+27 * rt_inv);
  k[r403b] = 0x1.550d4ac262c89p+59 *
             exp(-0x1.c3b5a9f54a42p+0 * lgt - 0x1.53bc570936f0dp+26 * rt_inv);
  k[r404f] = 0x1.74876e8p+36;
  k[r404b] = 0x1.dd242504ad1f7p+47 *
             exp(0x1.04b4dfcbdbbcp-1 * lgt - 0x1.93bacb2f43f3ep+28 * rt_inv);
  k[r405f] = 0x1.74876e8p+36;
  k[r405b] = 0x1.67ddb2832f16ap+52 *
             exp(0x1.798c5f208p-14 * lgt - 0x1.905f827455793p+28 * rt_inv);
  kTroe0 = 0x1.f50fd33a6ad68p+191 *
           exp(-0x1.823d70a3d70a4p+3 * lgt - 0x1.2cdbec0000001p+28 * rt_inv);
  kTroeInf = 0x1.6bcc41e9p+46 * exp(-0x1.f2c58cp+27 * rt_inv);
  fcTroe = 0x1.c28f5c28f5c28p-3 * exp(-temp / 0x1.ed86c6p+32) +
           0x1.8f5c28f5c28f6p-1 * exp(-temp / 0x1.b6e6666666666p+8) +
           0x1p+0 * exp(-0x1.3f7b1cp+29 / temp);
  k[r406f] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM40]);
  kTroe0 = 0x1.2c1d635a90165p+152 *
           exp(-0x1.46aeefe6adae4p+3 * lgt - 0x1.ebce5b3ad9a18p+27 * rt_inv);
  kTroeInf = 0x1.b3cc8d51cc17dp+6 *
             exp(0x1.dc7405e94aep+0 * lgt - 0x1.84dc0f3ad9a17p+27 * rt_inv);
  k[r406b] =
      GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM40]);
  k[r407f] = 0x1.bf08ebp+34;
  k[r407b] = 0x1.335d1755ea921p+68 *
             exp(-0x1.b7fd703f9ce4p-1 * lgt - 0x1.99daa1d81d938p+28 * rt_inv);
  k[r408f] = 0x1.2a05f2p+33;
  k[r408b] = 0x1.05f614c0d799cp+71 *
             exp(-0x1.6b48989d11aap+0 * lgt - 0x1.70607c21bc872p+28 * rt_inv);
  k[r409f] = 0x1.2a05f2p+33 * exp(-0x1.8c9f8c0000001p+27 * rt_inv);
  k[r409b] = 0x1.ba56bd1a1e6e6p+27 *
             exp(-0x1.d8c4e66ea1d4p-7 * lgt + 0x1.a74599c654093p+22 * rt_inv);
  k[r410f] = 0x1.31794b4p+34 * exp(-0x1.9efa6p+27 * rt_inv);
  k[r410b] = 0x1.55f54306e4981p+33 *
             exp(-0x1.0c0c27029cf23p-1 * lgt + 0x1.95ea440f96194p+19 * rt_inv);
  k[r411f] = 0x1.edbffffffffffp+13 *
             exp(0x1.ccccccccccccdp+0 * lgt - 0x1.dd088p+21 * rt_inv);
  k[r411b] = 0x1.4442af8e1459ap+9 *
             exp(0x1.05c0bde564024p+1 * lgt - 0x1.311c3401d3e34p+26 * rt_inv);
  k[r412f] = 0x1.41a7cp+22 *
             exp(0x1.f0a3d70a3d70ap-1 * lgt - 0x1.9504c00000001p+22 * rt_inv);
  k[r412b] = 0x1.3ea411acefe58p+22 *
             exp(0x1.6964219d511bp-1 * lgt - 0x1.2e17191619f91p+26 * rt_inv);
  k[r413f] = 0x1.347ae147ae147p+2 *
             exp(0x1.4cccccccccccdp+1 * lgt - 0x1.bc06c8p+25 * rt_inv);
  k[r413b] = 0x1.51c677606acc5p+2 *
             exp(0x1.2bbb61be41bf5p+1 * lgt - 0x1.0ae04abe1bac8p+22 * rt_inv);
  k[r414f] = 0x1.7ccccccccccccp+4 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.07311cp+26 * rt_inv);
  k[r414b] = 0x1.3a7aa10623d5p+9 *
             exp(0x1.c8556cfb456a2p+0 * lgt - 0x1.6374ce013e8cbp+23 * rt_inv);
  k[r415f] = 0x1.137ffffffffffp+8 *
             exp(0x1.399999999999ap+1 * lgt - 0x1.6959a00000001p+23 * rt_inv);
  k[r415b] = 0x1.3a2b47502393p-1 *
             exp(0x1.61b1f8a3a85dfp+1 * lgt - 0x1.264ca8d2b300bp+24 * rt_inv);
  k[r416f] = 0x1.eap+9 *
             exp(0x1.370a3d70a3d71p+1 * lgt - 0x1.2f40d00000001p+24 * rt_inv);
  k[r416b] = 0x1.a56fa8591f621p+5 *
             exp(0x1.1df857a079f9p+1 * lgt - 0x1.6b2c1d23cb564p+24 * rt_inv);
  k[r417f] =
      0x1.45p+9 * exp(0x1.3333333333333p+1 * lgt - 0x1.1d70e8p+24 * rt_inv);
  k[r417b] = 0x1.3ba08d8782eb2p+1 *
             exp(0x1.5f8ff02426d38p+1 * lgt - 0x1.e858f660501e9p+24 * rt_inv);
  k[r418f] = 0x1.4c7ffffffffffp+9 *
             exp(0x1.451eb851eb852p+1 * lgt - 0x1.af526p+24 * rt_inv);
  k[r418b] = 0x1.e716478e17d1ap+5 *
             exp(0x1.30513068a6886p+1 * lgt - 0x1.2242f158b4393p+25 * rt_inv);
  k[r419f] = 0x1.8bd66277c45cbp-11 *
             exp(0x1.bae147ae147aep+1 * lgt - 0x1.5dec18p+24 * rt_inv);
  k[r419b] = 0x1.5397abe7b433ap-8 *
             exp(0x1.a51ca25ad64b6p+1 * lgt - 0x1.583e4af867cb6p+25 * rt_inv);
  k[r420f] = 0x1.d9f4d37c1376dp-12 *
             exp(0x1.d333333333333p+1 * lgt - 0x1.c8bb300000001p+24 * rt_inv);
  k[r420b] = 0x1.32abe8799e763p-4 *
             exp(0x1.7c444905bc5d8p+1 * lgt - 0x1.72cb9120f3f66p+25 * rt_inv);
  k[r421f] = 0x1.054e88p+29 * exp(-0x1.3f36cp+24 * rt_inv);
  k[r421b] = 0x1.16f65a8d3c1b9p+20 *
             exp(0x1.5c368f0e720cp-1 * lgt - 0x1.2aafd713253cp+25 * rt_inv);
  k[r422f] = 0x1.9de508p+27 * exp(-0x1.9c4bfp+24 * rt_inv);
  k[r422b] = 0x1.4d41a9cd6521p+23 *
             exp(0x1.5e35ee963d78p-3 * lgt - 0x1.3e60293bb1687p+25 * rt_inv);
  k[r423f] = 0x1.347ae147ae147p+2 *
             exp(0x1.4cccccccccccdp+1 * lgt - 0x1.bc06c8p+25 * rt_inv);
  k[r423b] = 0x1.018aefaa09e7ep+7 *
             exp(0x1.c631f5026441dp+0 * lgt + 0x1.4d596a369ef7p+22 * rt_inv);
  k[r424f] = 0x1.7ccccccccccccp+4 *
             exp(0x1.4666666666666p+1 * lgt - 0x1.07311cp+26 * rt_inv);
  k[r424b] = 0x1.df8f67eabe1d7p+13 *
             exp(0x1.37109e8126266p+0 * lgt - 0x1.babf9c370a1dep+20 * rt_inv);
  k[r425f] = 0x1.517ffffffffffp+12 *
             exp(0x1.e666666666666p+0 * lgt - 0x1.0f7dcc0000001p+26 * rt_inv);
  k[r425b] = 0x1.646095363db77p+12 *
             exp(0x1.e6717c5fb70dp+0 * lgt - 0x1.459e760bee39ap+24 * rt_inv);
  k[r426f] = 0x1.9ap+6 * exp(0x1.4p+1 * lgt - 0x1.2627d4p+26 * rt_inv);
  k[r426b] = 0x1.4685ee2b55b6ep+11 *
             exp(0x1.fdb68c44df506p+0 * lgt - 0x1.6a920a5d06923p+24 * rt_inv);
  k[r427f] = 0x1.7b47f762p+40 *
             exp(-0x1.70a3d70a3d70ap-3 * lgt - 0x1.448f8a0000001p+27 * rt_inv);
  k[r427b] = 0x1.1e0f05ed67b38p+45 *
             exp(-0x1.60d2092b6dfep-1 * lgt - 0x1.3dd8f88a23194p+27 * rt_inv);
  k[r428f] = 0x1.73dffffffffffp+15 *
             exp(0x1.8a3d70a3d70a4p+0 * lgt - 0x1.14eb04p+27 * rt_inv);
  k[r428b] = 0x1.906d1fefe8086p+44 *
             exp(0x1.49c4306b5f88p-3 * lgt - 0x1.f43ef5649bb1p+25 * rt_inv);
  k[r429f] = 0x1.837ffffffffffp+10 *
             exp(0x1.028f5c28f5c29p+1 * lgt - 0x1.6dd1b00000001p+24 * rt_inv);
  k[r429b] = 0x1.5955bddb4b3bfp+41 *
             exp(0x1.ad8bf3dfbb5p-1 * lgt - 0x1.7cd1a113df8fep+26 * rt_inv);
  k[r430f] = 0x1.1e1a3p+27 * exp(-0x1.7bdd5p+25 * rt_inv);
  k[r430b] = 0x1.0573c6178f796p+70 *
             exp(-0x1.20e852345452p+1 * lgt - 0x1.16d7d546733bap+27 * rt_inv);

  w[r1f] = k[r1f] * c[sH] * c[sO2];
  w[r1b] = k[r1b] * c[sOH] * c[sO];
  w[r2f] = k[r2f] * c[sO] * c[sH2];
  w[r2b] = k[r2b] * c[sOH] * c[sH];
  w[r3f] = k[r3f] * c[sOH] * c[sH2];
  w[r3b] = k[r3b] * c[sH2O] * c[sH];
  w[r4f] = k[r4f] * c[sO] * c[sH2O];
  w[r4b] = k[r4b] * c[sOH] * c[sOH];
  w[r5f] = k[r5f] * c[sH2] * M[mM1];
  w[r5b] = k[r5b] * c[sH] * c[sH] * M[mM1];
  w[r6f] = k[r6f] * c[sO] * c[sO] * M[mM2];
  w[r6b] = k[r6b] * c[sO2] * M[mM2];
  w[r7f] = k[r7f] * c[sO] * c[sH] * M[mM3];
  w[r7b] = k[r7b] * c[sOH] * M[mM3];
  w[r8f] = k[r8f] * c[sH] * c[sOH] * M[mM4];
  w[r8b] = k[r8b] * c[sH2O] * M[mM4];
  w[r9f] = k[r9f] * c[sH] * c[sO2];
  w[r9b] = k[r9b] * c[sHO2];
  w[r10f] = k[r10f] * c[sH] * c[sO2];
  w[r10b] = k[r10b] * c[sHO2];
  w[r11f] = k[r11f] * c[sH] * c[sO2];
  w[r11b] = k[r11b] * c[sHO2];
  w[r12f] = k[r12f] * c[sHO2] * c[sH];
  w[r12b] = k[r12b] * c[sOH] * c[sOH];
  w[r13f] = k[r13f] * c[sH2] * c[sO2];
  w[r13b] = k[r13b] * c[sHO2] * c[sH];
  w[r14f] = k[r14f] * c[sHO2] * c[sO];
  w[r14b] = k[r14b] * c[sO2] * c[sOH];
  w[r15f] = k[r15f] * c[sHO2] * c[sOH];
  w[r15b] = k[r15b] * c[sO2] * c[sH2O];
  w[r16f] = k[r16f] * c[sHO2] * c[sHO2];
  w[r16b] = k[r16b] * c[sO2] * c[sH2O2];
  w[r17f] = k[r17f] * c[sHO2] * c[sHO2];
  w[r17b] = k[r17b] * c[sO2] * c[sH2O2];
  w[r18f] = k[r18f] * c[sH2O2];
  w[r18b] = k[r18b] * c[sOH] * c[sOH];
  w[r19f] = k[r19f] * c[sH2O2];
  w[r19b] = k[r19b] * c[sOH] * c[sOH];
  w[r20f] = k[r20f] * c[sH2O2] * c[sH];
  w[r20b] = k[r20b] * c[sOH] * c[sH2O];
  w[r21f] = k[r21f] * c[sH2O2] * c[sH];
  w[r21b] = k[r21b] * c[sHO2] * c[sH2];
  w[r22f] = k[r22f] * c[sH2O2] * c[sO];
  w[r22b] = k[r22b] * c[sHO2] * c[sOH];
  w[r23f] = k[r23f] * c[sH2O2] * c[sOH];
  w[r23b] = k[r23b] * c[sHO2] * c[sH2O];
  w[r24f] = k[r24f] * c[sH2O2] * c[sOH];
  w[r24b] = k[r24b] * c[sHO2] * c[sH2O];
  w[r25f] = k[r25f] * c[sCO] * c[sO];
  w[r25b] = k[r25b] * c[sCO2];
  w[r26f] = k[r26f] * c[sCO] * c[sO2];
  w[r26b] = k[r26b] * c[sO] * c[sCO2];
  w[r27f] = k[r27f] * c[sCO] * c[sOH];
  w[r27b] = k[r27b] * c[sH] * c[sCO2];
  w[r28f] = k[r28f] * c[sCO] * c[sOH];
  w[r28b] = k[r28b] * c[sH] * c[sCO2];
  w[r29f] = k[r29f] * c[sCO] * c[sHO2];
  w[r29b] = k[r29b] * c[sOH] * c[sCO2];
  w[r30f] = k[r30f] * c[sHCO] * M[mM11];
  w[r30b] = k[r30b] * c[sCO] * c[sH] * M[mM11];
  w[r31f] = k[r31f] * c[sHCO] * c[sO2];
  w[r31b] = k[r31b] * c[sHO2] * c[sCO];
  w[r32f] = k[r32f] * c[sHCO] * c[sH];
  w[r32b] = k[r32b] * c[sH2] * c[sCO];
  w[r33f] = k[r33f] * c[sHCO] * c[sO];
  w[r33b] = k[r33b] * c[sOH] * c[sCO];
  w[r34f] = k[r34f] * c[sHCO] * c[sO];
  w[r34b] = k[r34b] * c[sH] * c[sCO2];
  w[r35f] = k[r35f] * c[sHCO] * c[sOH];
  w[r35b] = k[r35b] * c[sH2O] * c[sCO];
  w[r36] = k[r36] * c[sHCO] * c[sHO2];
  w[r37] = k[r37] * c[sHCO] * c[sHCO];
  w[r38f] = k[r38f] * c[sHCO] * c[sCH3];
  w[r38b] = k[r38b] * c[sCO] * c[sCH4];
  w[r39f] = k[r39f] * c[sCH2O] * c[sO2];
  w[r39b] = k[r39b] * c[sHO2] * c[sHCO];
  w[r40f] = k[r40f] * c[sHCO] * c[sO2];
  w[r40b] = k[r40b] * c[sO2CHO];
  w[r41f] = k[r41f] * c[sCH2O] * c[sO2CHO];
  w[r41b] = k[r41b] * c[sHO2CHO] * c[sHCO];
  w[r42f] = k[r42f] * c[sOCHO] * c[sOH];
  w[r42b] = k[r42b] * c[sHO2CHO];
  w[r43f] = k[r43f] * c[sH] * c[sCO2];
  w[r43b] = k[r43b] * c[sOCHO];
  w[r44f] = k[r44f] * c[sHCO] * c[sHCO];
  w[r44b] = k[r44b] * c[sCO] * c[sCH2O];
  w[r45f] = k[r45f] * c[sOHY] * c[sH2O];
  w[r45b] = k[r45b] * c[sH2O] * c[sOH];
  w[r46f] = k[r46f] * c[sOHY] * c[sH2];
  w[r46b] = k[r46b] * c[sH2] * c[sOH];
  w[r47f] = k[r47f] * c[sOHY] * c[sN2];
  w[r47b] = k[r47b] * c[sN2] * c[sOH];
  w[r48f] = k[r48f] * c[sOHY] * c[sOH];
  w[r48b] = k[r48b] * c[sOH] * c[sOH];
  w[r49f] = k[r49f] * c[sOHY] * c[sH];
  w[r49b] = k[r49b] * c[sH] * c[sOH];
  w[r50f] = k[r50f] * c[sOHY] * c[sAR];
  w[r50b] = k[r50b] * c[sAR] * c[sOH];
  w[r51f] = k[r51f] * c[sOHY];
  w[r51b] = k[r51b] * c[sOH];
  w[r52f] = k[r52f] * c[sOHY] * c[sO2];
  w[r52b] = k[r52b] * c[sO2] * c[sOH];
  w[r53f] = k[r53f] * c[sOHY] * c[sCO2];
  w[r53b] = k[r53b] * c[sCO2] * c[sOH];
  w[r54f] = k[r54f] * c[sOHY] * c[sCO];
  w[r54b] = k[r54b] * c[sCO] * c[sOH];
  w[r55f] = k[r55f] * c[sOHY] * c[sCH4];
  w[r55b] = k[r55b] * c[sCH4] * c[sOH];
  w[r56f] = k[r56f] * c[sCH] * c[sO2];
  w[r56b] = k[r56b] * c[sOHY] * c[sCO];
  w[r57f] = k[r57f] * c[sHCO] * c[sH];
  w[r57b] = k[r57b] * c[sCH2O];
  w[r58f] = k[r58f] * c[sCO] * c[sH2];
  w[r58b] = k[r58b] * c[sCH2O];
  w[r59f] = k[r59f] * c[sCH2O] * c[sOH];
  w[r59b] = k[r59b] * c[sH2O] * c[sHCO];
  w[r60f] = k[r60f] * c[sCH2O] * c[sH];
  w[r60b] = k[r60b] * c[sH2] * c[sHCO];
  w[r61f] = k[r61f] * c[sCH2O] * c[sO];
  w[r61b] = k[r61b] * c[sOH] * c[sHCO];
  w[r62f] = k[r62f] * c[sCH2O] * c[sCH3];
  w[r62b] = k[r62b] * c[sCH4] * c[sHCO];
  w[r63f] = k[r63f] * c[sCH2O] * c[sHO2];
  w[r63b] = k[r63b] * c[sH2O2] * c[sHCO];
  w[r64f] = k[r64f] * c[sCH2O] * c[sOH];
  w[r64b] = k[r64b] * c[sHOCH2O];
  w[r65f] = k[r65f] * c[sHOCH2O];
  w[r65b] = k[r65b] * c[sH] * c[sHOCHO];
  w[r66f] = k[r66f] * c[sHOCHO];
  w[r66b] = k[r66b] * c[sH2O] * c[sCO];
  w[r67f] = k[r67f] * c[sHOCHO];
  w[r67b] = k[r67b] * c[sH2] * c[sCO2];
  w[r68f] = k[r68f] * c[sOCHO] * c[sHO2];
  w[r68b] = k[r68b] * c[sO2] * c[sHOCHO];
  w[r69] = k[r69] * c[sHOCHO] * c[sOH];
  w[r70] = k[r70] * c[sHOCHO] * c[sOH];
  w[r71] = k[r71] * c[sHOCHO] * c[sH];
  w[r72] = k[r72] * c[sHOCHO] * c[sH];
  w[r73] = k[r73] * c[sHOCHO] * c[sCH3];
  w[r74f] = k[r74f] * c[sOCHO] * c[sH2O2];
  w[r74b] = k[r74b] * c[sHO2] * c[sHOCHO];
  w[r75] = k[r75] * c[sHOCHO] * c[sHO2];
  w[r76] = k[r76] * c[sHOCHO] * c[sO];
  w[r77f] = k[r77f] * c[sCH2O] * c[sOCHO];
  w[r77b] = k[r77b] * c[sHCO] * c[sHOCHO];
  w[r78f] = k[r78f] * c[sCH3O];
  w[r78b] = k[r78b] * c[sH] * c[sCH2O];
  w[r79f] = k[r79f] * c[sCH3O] * c[sO2];
  w[r79b] = k[r79b] * c[sHO2] * c[sCH2O];
  w[r80f] = k[r80f] * c[sCH2O] * c[sCH3O];
  w[r80b] = k[r80b] * c[sHCO] * c[sCH3OH];
  w[r81f] = k[r81f] * c[sCH3] * c[sCH3OH];
  w[r81b] = k[r81b] * c[sCH3O] * c[sCH4];
  w[r82f] = k[r82f] * c[sCH3O] * c[sCH3];
  w[r82b] = k[r82b] * c[sCH4] * c[sCH2O];
  w[r83f] = k[r83f] * c[sCH3O] * c[sH];
  w[r83b] = k[r83b] * c[sH2] * c[sCH2O];
  w[r84f] = k[r84f] * c[sCH3O] * c[sHO2];
  w[r84b] = k[r84b] * c[sH2O2] * c[sCH2O];
  w[r85f] = k[r85f] * c[sCH2O] * c[sH];
  w[r85b] = k[r85b] * c[sCH2OH];
  w[r86f] = k[r86f] * c[sCH2OH] * c[sO2];
  w[r86b] = k[r86b] * c[sHO2] * c[sCH2O];
  w[r87f] = k[r87f] * c[sCH2OH] * c[sO2];
  w[r87b] = k[r87b] * c[sHO2] * c[sCH2O];
  w[r88f] = k[r88f] * c[sCH2OH] * c[sH];
  w[r88b] = k[r88b] * c[sH2] * c[sCH2O];
  w[r89f] = k[r89f] * c[sCH2OH] * c[sHO2];
  w[r89b] = k[r89b] * c[sH2O2] * c[sCH2O];
  w[r90f] = k[r90f] * c[sCH2OH] * c[sHCO];
  w[r90b] = k[r90b] * c[sCH2O] * c[sCH2O];
  w[r91f] = k[r91f] * c[sCH2OH] * c[sCH3O];
  w[r91b] = k[r91b] * c[sCH3OH] * c[sCH2O];
  w[r92f] = k[r92f] * c[sCH3OH] * c[sHCO];
  w[r92b] = k[r92b] * c[sCH2O] * c[sCH2OH];
  w[r93f] = k[r93f] * c[sOH] * c[sCH2OH];
  w[r93b] = k[r93b] * c[sCH2O] * c[sH2O];
  w[r94f] = k[r94f] * c[sO] * c[sCH2OH];
  w[r94b] = k[r94b] * c[sCH2O] * c[sOH];
  w[r95f] = k[r95f] * c[sCH2OH] * c[sCH2OH];
  w[r95b] = k[r95b] * c[sCH3OH] * c[sCH2O];
  w[r96f] = k[r96f] * c[sCH2OH] * c[sHO2];
  w[r96b] = k[r96b] * c[sOH] * c[sHOCH2O];
  w[r97f] = k[r97f] * c[sCH2O] * c[sHO2];
  w[r97b] = k[r97b] * c[sOCH2O2H];
  w[r98f] = k[r98f] * c[sOCH2O2H];
  w[r98b] = k[r98b] * c[sHOCH2O2];
  w[r99f] = k[r99f] * c[sHOCH2O2] * c[sHO2];
  w[r99b] = k[r99b] * c[sO2] * c[sHOCH2O2H];
  w[r100f] = k[r100f] * c[sHOCH2O] * c[sOH];
  w[r100b] = k[r100b] * c[sHOCH2O2H];
  w[r101f] = k[r101f] * c[sCH3OH];
  w[r101b] = k[r101b] * c[sOH] * c[sCH3];
  w[r102f] = k[r102f] * c[sCH3OH];
  w[r102b] = k[r102b] * c[sH2O] * c[sCH2Y];
  w[r103f] = k[r103f] * c[sCH3OH];
  w[r103b] = k[r103b] * c[sH] * c[sCH2OH];
  w[r104f] = k[r104f] * c[sCH3OH] * c[sH];
  w[r104b] = k[r104b] * c[sH2] * c[sCH2OH];
  w[r105f] = k[r105f] * c[sCH3OH] * c[sH];
  w[r105b] = k[r105b] * c[sH2] * c[sCH3O];
  w[r106f] = k[r106f] * c[sCH3OH] * c[sO];
  w[r106b] = k[r106b] * c[sOH] * c[sCH2OH];
  w[r107f] = k[r107f] * c[sCH3OH] * c[sOH];
  w[r107b] = k[r107b] * c[sH2O] * c[sCH2OH];
  w[r108f] = k[r108f] * c[sCH3OH] * c[sOH];
  w[r108b] = k[r108b] * c[sH2O] * c[sCH3O];
  w[r109f] = k[r109f] * c[sCH3OH] * c[sO2];
  w[r109b] = k[r109b] * c[sHO2] * c[sCH2OH];
  w[r110f] = k[r110f] * c[sCH3OH] * c[sHO2];
  w[r110b] = k[r110b] * c[sH2O2] * c[sCH2OH];
  w[r111f] = k[r111f] * c[sCH3OH] * c[sCH3];
  w[r111b] = k[r111b] * c[sCH4] * c[sCH2OH];
  w[r112f] = k[r112f] * c[sCH3O] * c[sCH3OH];
  w[r112b] = k[r112b] * c[sCH3OH] * c[sCH2OH];
  w[r113f] = k[r113f] * c[sCH3O] * c[sCH3O];
  w[r113b] = k[r113b] * c[sCH2O] * c[sCH3OH];
  w[r114f] = k[r114f] * c[sCH3] * c[sH];
  w[r114b] = k[r114b] * c[sCH4];
  w[r115f] = k[r115f] * c[sCH4] * c[sH];
  w[r115b] = k[r115b] * c[sH2] * c[sCH3];
  w[r116f] = k[r116f] * c[sCH4] * c[sOH];
  w[r116b] = k[r116b] * c[sH2O] * c[sCH3];
  w[r117f] = k[r117f] * c[sCH4] * c[sO];
  w[r117b] = k[r117b] * c[sOH] * c[sCH3];
  w[r118f] = k[r118f] * c[sCH4] * c[sHO2];
  w[r118b] = k[r118b] * c[sH2O2] * c[sCH3];
  w[r119f] = k[r119f] * c[sCH4] * c[sCH2];
  w[r119b] = k[r119b] * c[sCH3] * c[sCH3];
  w[r120f] = k[r120f] * c[sCH3] * c[sOH];
  w[r120b] = k[r120b] * c[sH2O] * c[sCH2Y];
  w[r121f] = k[r121f] * c[sCH3] * c[sOH];
  w[r121b] = k[r121b] * c[sH2] * c[sCH2O];
  w[r122f] = k[r122f] * c[sCH3] * c[sOH];
  w[r122b] = k[r122b] * c[sH] * c[sCH2OH];
  w[r123f] = k[r123f] * c[sCH3] * c[sOH];
  w[r123b] = k[r123b] * c[sCH3O] * c[sH];
  w[r124f] = k[r124f] * c[sCH3] * c[sOH];
  w[r124b] = k[r124b] * c[sH2] * c[sHCOH];
  w[r125f] = k[r125f] * c[sHCOH] * c[sOH];
  w[r125b] = k[r125b] * c[sH2O] * c[sHCO];
  w[r126f] = k[r126f] * c[sHCOH] * c[sH];
  w[r126b] = k[r126b] * c[sH] * c[sCH2O];
  w[r127] = k[r127] * c[sHCOH] * c[sO];
  w[r128] = k[r128] * c[sHCOH] * c[sO];
  w[r129] = k[r129] * c[sHCOH] * c[sO2];
  w[r130f] = k[r130f] * c[sHCOH] * c[sO2];
  w[r130b] = k[r130b] * c[sH2O] * c[sCO2];
  w[r131f] = k[r131f] * c[sCH3] * c[sHO2];
  w[r131b] = k[r131b] * c[sOH] * c[sCH3O];
  w[r132f] = k[r132f] * c[sCH3] * c[sHO2];
  w[r132b] = k[r132b] * c[sO2] * c[sCH4];
  w[r133f] = k[r133f] * c[sCH3] * c[sO];
  w[r133b] = k[r133b] * c[sH] * c[sCH2O];
  w[r134f] = k[r134f] * c[sCH3] * c[sO2];
  w[r134b] = k[r134b] * c[sO] * c[sCH3O];
  w[r135f] = k[r135f] * c[sCH3] * c[sO2];
  w[r135b] = k[r135b] * c[sOH] * c[sCH2O];
  w[r136f] = k[r136f] * c[sCH3] * c[sO2];
  w[r136b] = k[r136b] * c[sCH3O2];
  w[r137f] = k[r137f] * c[sCH3O2] * c[sCH2O];
  w[r137b] = k[r137b] * c[sHCO] * c[sCH3O2H];
  w[r138f] = k[r138f] * c[sCH4] * c[sCH3O2];
  w[r138b] = k[r138b] * c[sCH3O2H] * c[sCH3];
  w[r139f] = k[r139f] * c[sCH3OH] * c[sCH3O2];
  w[r139b] = k[r139b] * c[sCH3O2H] * c[sCH2OH];
  w[r140f] = k[r140f] * c[sCH3O2] * c[sCH3];
  w[r140b] = k[r140b] * c[sCH3O] * c[sCH3O];
  w[r141f] = k[r141f] * c[sCH3O2] * c[sHO2];
  w[r141b] = k[r141b] * c[sO2] * c[sCH3O2H];
  w[r142] = k[r142] * c[sCH3O2] * c[sCH3O2];
  w[r143] = k[r143] * c[sCH3O2] * c[sCH3O2];
  w[r144f] = k[r144f] * c[sCH3O2] * c[sH];
  w[r144b] = k[r144b] * c[sOH] * c[sCH3O];
  w[r145f] = k[r145f] * c[sCH3O2] * c[sO];
  w[r145b] = k[r145b] * c[sO2] * c[sCH3O];
  w[r146f] = k[r146f] * c[sCH3O2] * c[sOH];
  w[r146b] = k[r146b] * c[sO2] * c[sCH3OH];
  w[r147f] = k[r147f] * c[sCH3O2H];
  w[r147b] = k[r147b] * c[sOH] * c[sCH3O];
  w[r148f] = k[r148f] * c[sCH2Y] * c[sN2];
  w[r148b] = k[r148b] * c[sN2] * c[sCH2];
  w[r149f] = k[r149f] * c[sCH2Y] * c[sAR];
  w[r149b] = k[r149b] * c[sAR] * c[sCH2];
  w[r150f] = k[r150f] * c[sCH2Y] * c[sH];
  w[r150b] = k[r150b] * c[sH2] * c[sCH];
  w[r151f] = k[r151f] * c[sCH2Y] * c[sO];
  w[r151b] = k[r151b] * c[sH2] * c[sCO];
  w[r152f] = k[r152f] * c[sCH2Y] * c[sO];
  w[r152b] = k[r152b] * c[sH] * c[sHCO];
  w[r153f] = k[r153f] * c[sCH2Y] * c[sOH];
  w[r153b] = k[r153b] * c[sH] * c[sCH2O];
  w[r154f] = k[r154f] * c[sCH2Y] * c[sH2];
  w[r154b] = k[r154b] * c[sH] * c[sCH3];
  w[r155] = k[r155] * c[sCH2Y] * c[sO2];
  w[r156f] = k[r156f] * c[sCH2Y] * c[sO2];
  w[r156b] = k[r156b] * c[sH2O] * c[sCO];
  w[r157f] = k[r157f] * c[sCH2Y] * c[sH2O];
  w[r157b] = k[r157b] * c[sH2O] * c[sCH2];
  w[r158f] = k[r158f] * c[sCH2Y] * c[sCO];
  w[r158b] = k[r158b] * c[sCO] * c[sCH2];
  w[r159f] = k[r159f] * c[sCH2Y] * c[sCO2];
  w[r159b] = k[r159b] * c[sCO2] * c[sCH2];
  w[r160f] = k[r160f] * c[sCH2Y] * c[sCO2];
  w[r160b] = k[r160b] * c[sCO] * c[sCH2O];
  w[r161f] = k[r161f] * c[sCH2] * c[sH];
  w[r161b] = k[r161b] * c[sCH3];
  w[r162f] = k[r162f] * c[sCH2] * c[sO2];
  w[r162b] = k[r162b] * c[sOH] * c[sHCO];
  w[r163] = k[r163] * c[sCH2] * c[sO2];
  w[r164] = k[r164] * c[sCH2] * c[sO];
  w[r165f] = k[r165f] * c[sCH2] * c[sH];
  w[r165b] = k[r165b] * c[sH2] * c[sCH];
  w[r166f] = k[r166f] * c[sCH2] * c[sH];
  w[r166b] = k[r166b] * c[sH2] * c[sCH];
  w[r167f] = k[r167f] * c[sCH2] * c[sOH];
  w[r167b] = k[r167b] * c[sH2O] * c[sCH];
  w[r168f] = k[r168f] * c[sCH] * c[sO2];
  w[r168b] = k[r168b] * c[sO] * c[sHCO];
  w[r169f] = k[r169f] * c[sCH] * c[sO];
  w[r169b] = k[r169b] * c[sH] * c[sCO];
  w[r170f] = k[r170f] * c[sCH] * c[sOH];
  w[r170b] = k[r170b] * c[sH] * c[sHCO];
  w[r171f] = k[r171f] * c[sCH] * c[sH2O];
  w[r171b] = k[r171b] * c[sCH2O] * c[sH];
  w[r172f] = k[r172f] * c[sCH] * c[sCO2];
  w[r172b] = k[r172b] * c[sCO] * c[sHCO];
  w[r173f] = k[r173f] * c[sCH3] * c[sCH3];
  w[r173b] = k[r173b] * c[sC2H6];
  w[r174f] = k[r174f] * c[sC2H5] * c[sH];
  w[r174b] = k[r174b] * c[sC2H6];
  w[r175f] = k[r175f] * c[sC2H6] * c[sH];
  w[r175b] = k[r175b] * c[sH2] * c[sC2H5];
  w[r176f] = k[r176f] * c[sC2H6] * c[sO];
  w[r176b] = k[r176b] * c[sOH] * c[sC2H5];
  w[r177f] = k[r177f] * c[sC2H6] * c[sOH];
  w[r177b] = k[r177b] * c[sH2O] * c[sC2H5];
  w[r178f] = k[r178f] * c[sC2H6] * c[sO2];
  w[r178b] = k[r178b] * c[sHO2] * c[sC2H5];
  w[r179f] = k[r179f] * c[sC2H6] * c[sCH3];
  w[r179b] = k[r179b] * c[sCH4] * c[sC2H5];
  w[r180f] = k[r180f] * c[sC2H6] * c[sHO2];
  w[r180b] = k[r180b] * c[sH2O2] * c[sC2H5];
  w[r181f] = k[r181f] * c[sC2H6] * c[sCH3O2];
  w[r181b] = k[r181b] * c[sCH3O2H] * c[sC2H5];
  w[r182f] = k[r182f] * c[sC2H6] * c[sCH3O];
  w[r182b] = k[r182b] * c[sCH3OH] * c[sC2H5];
  w[r183f] = k[r183f] * c[sC2H6] * c[sCH];
  w[r183b] = k[r183b] * c[sCH2] * c[sC2H5];
  w[r184f] = k[r184f] * c[sCH2Y] * c[sC2H6];
  w[r184b] = k[r184b] * c[sC2H5] * c[sCH3];
  w[r185f] = k[r185f] * c[sC2H4] * c[sH];
  w[r185b] = k[r185b] * c[sC2H5];
  w[r186f] = k[r186f] * c[sH2] * c[sCH3O2];
  w[r186b] = k[r186b] * c[sCH3O2H] * c[sH];
  w[r187f] = k[r187f] * c[sH2] * c[sC2H5O2];
  w[r187b] = k[r187b] * c[sC2H5O2H] * c[sH];
  w[r188f] = k[r188f] * c[sC2H4] * c[sC2H4];
  w[r188b] = k[r188b] * c[sC2H3] * c[sC2H5];
  w[r189f] = k[r189f] * c[sCH3] * c[sC2H5];
  w[r189b] = k[r189b] * c[sC2H4] * c[sCH4];
  w[r190f] = k[r190f] * c[sCH3] * c[sCH3];
  w[r190b] = k[r190b] * c[sC2H5] * c[sH];
  w[r191f] = k[r191f] * c[sC2H5] * c[sH];
  w[r191b] = k[r191b] * c[sH2] * c[sC2H4];
  w[r192f] = k[r192f] * c[sC2H5] * c[sO];
  w[r192b] = k[r192b] * c[sH] * c[sCH3CHO];
  w[r193f] = k[r193f] * c[sC2H5] * c[sHO2];
  w[r193b] = k[r193b] * c[sOH] * c[sC2H5O];
  w[r194f] = k[r194f] * c[sCH3O2] * c[sC2H5];
  w[r194b] = k[r194b] * c[sC2H5O] * c[sCH3O];
  w[r195f] = k[r195f] * c[sC2H5O] * c[sO2];
  w[r195b] = k[r195b] * c[sHO2] * c[sCH3CHO];
  w[r196f] = k[r196f] * c[sCH3] * c[sCH2O];
  w[r196b] = k[r196b] * c[sC2H5O];
  w[r197f] = k[r197f] * c[sCH3CHO] * c[sH];
  w[r197b] = k[r197b] * c[sC2H5O];
  w[r198f] = k[r198f] * c[sC2H5O2] * c[sCH2O];
  w[r198b] = k[r198b] * c[sHCO] * c[sC2H5O2H];
  w[r199f] = k[r199f] * c[sCH4] * c[sC2H5O2];
  w[r199b] = k[r199b] * c[sC2H5O2H] * c[sCH3];
  w[r200f] = k[r200f] * c[sCH3OH] * c[sC2H5O2];
  w[r200b] = k[r200b] * c[sC2H5O2H] * c[sCH2OH];
  w[r201f] = k[r201f] * c[sC2H5O2] * c[sHO2];
  w[r201b] = k[r201b] * c[sO2] * c[sC2H5O2H];
  w[r202f] = k[r202f] * c[sC2H6] * c[sC2H5O2];
  w[r202b] = k[r202b] * c[sC2H5O2H] * c[sC2H5];
  w[r203f] = k[r203f] * c[sC2H5O2H];
  w[r203b] = k[r203b] * c[sOH] * c[sC2H5O];
  w[r204f] = k[r204f] * c[sC2H5] * c[sO2];
  w[r204b] = k[r204b] * c[sC2H5O2];
  w[r205f] = k[r205f] * c[sC2H5] * c[sO2];
  w[r205b] = k[r205b] * c[sHO2] * c[sC2H4];
  w[r206f] = k[r206f] * c[sC2H5] * c[sO2];
  w[r206b] = k[r206b] * c[sHO2] * c[sC2H4];
  w[r207f] = k[r207f] * c[sC2H5] * c[sO2];
  w[r207b] = k[r207b] * c[sOH] * c[sC2H4O1X2];
  w[r208f] = k[r208f] * c[sC2H5] * c[sO2];
  w[r208b] = k[r208b] * c[sOH] * c[sCH3CHO];
  w[r209f] = k[r209f] * c[sC2H5O2];
  w[r209b] = k[r209b] * c[sOH] * c[sCH3CHO];
  w[r210f] = k[r210f] * c[sC2H5O2];
  w[r210b] = k[r210b] * c[sHO2] * c[sC2H4];
  w[r211f] = k[r211f] * c[sC2H5O2];
  w[r211b] = k[r211b] * c[sOH] * c[sC2H4O1X2];
  w[r212f] = k[r212f] * c[sC2H4O1X2];
  w[r212b] = k[r212b] * c[sHCO] * c[sCH3];
  w[r213f] = k[r213f] * c[sC2H4O1X2];
  w[r213b] = k[r213b] * c[sCH3CHO];
  w[r214f] = k[r214f] * c[sCH3CHO];
  w[r214b] = k[r214b] * c[sHCO] * c[sCH3];
  w[r215f] = k[r215f] * c[sCH3CHO];
  w[r215b] = k[r215b] * c[sCO] * c[sCH4];
  w[r216f] = k[r216f] * c[sCH3CHO] * c[sH];
  w[r216b] = k[r216b] * c[sH2] * c[sCH3CO];
  w[r217f] = k[r217f] * c[sCH3CHO] * c[sH];
  w[r217b] = k[r217b] * c[sH2] * c[sCH2CHO];
  w[r218f] = k[r218f] * c[sCH3CHO] * c[sO];
  w[r218b] = k[r218b] * c[sOH] * c[sCH3CO];
  w[r219f] = k[r219f] * c[sCH3CHO] * c[sOH];
  w[r219b] = k[r219b] * c[sH2O] * c[sCH3CO];
  w[r220f] = k[r220f] * c[sCH3CHO] * c[sO2];
  w[r220b] = k[r220b] * c[sHO2] * c[sCH3CO];
  w[r221f] = k[r221f] * c[sCH3CHO] * c[sCH3];
  w[r221b] = k[r221b] * c[sCH4] * c[sCH3CO];
  w[r222f] = k[r222f] * c[sCH3CHO] * c[sHO2];
  w[r222b] = k[r222b] * c[sH2O2] * c[sCH3CO];
  w[r223f] = k[r223f] * c[sCH3O2] * c[sCH3CHO];
  w[r223b] = k[r223b] * c[sCH3CO] * c[sCH3O2H];
  w[r224f] = k[r224f] * c[sCH3CHO] * c[sCH3CO3];
  w[r224b] = k[r224b] * c[sCH3CO3H] * c[sCH3CO];
  w[r225f] = k[r225f] * c[sCH3CHO] * c[sOH];
  w[r225b] = k[r225b] * c[sHOCHO] * c[sCH3];
  w[r226f] = k[r226f] * c[sCH3CHO] * c[sOH];
  w[r226b] = k[r226b] * c[sH2O] * c[sCH2CHO];
  w[r227f] = k[r227f] * c[sCH3CO];
  w[r227b] = k[r227b] * c[sCO] * c[sCH3];
  w[r228f] = k[r228f] * c[sCH3CO] * c[sH];
  w[r228b] = k[r228b] * c[sH2] * c[sCH2CO];
  w[r229f] = k[r229f] * c[sCH3CO] * c[sO];
  w[r229b] = k[r229b] * c[sOH] * c[sCH2CO];
  w[r230f] = k[r230f] * c[sCH3CO] * c[sCH3];
  w[r230b] = k[r230b] * c[sCH4] * c[sCH2CO];
  w[r231f] = k[r231f] * c[sCH3CO] * c[sO2];
  w[r231b] = k[r231b] * c[sCH3CO3];
  w[r232f] = k[r232f] * c[sCH3CO3] * c[sHO2];
  w[r232b] = k[r232b] * c[sO2] * c[sCH3CO3H];
  w[r233f] = k[r233f] * c[sH2O2] * c[sCH3CO3];
  w[r233b] = k[r233b] * c[sCH3CO3H] * c[sHO2];
  w[r234f] = k[r234f] * c[sCH4] * c[sCH3CO3];
  w[r234b] = k[r234b] * c[sCH3CO3H] * c[sCH3];
  w[r235f] = k[r235f] * c[sCH2O] * c[sCH3CO3];
  w[r235b] = k[r235b] * c[sCH3CO3H] * c[sHCO];
  w[r236f] = k[r236f] * c[sC2H6] * c[sCH3CO3];
  w[r236b] = k[r236b] * c[sCH3CO3H] * c[sC2H5];
  w[r237f] = k[r237f] * c[sCH3CO3H];
  w[r237b] = k[r237b] * c[sOH] * c[sCH3CO2];
  w[r238f] = k[r238f] * c[sCH3CO2] * M[mM28];
  w[r238b] = k[r238b] * c[sCO2] * c[sCH3] * M[mM28];
  w[r239f] = k[r239f] * c[sCH2CHO];
  w[r239b] = k[r239b] * c[sH] * c[sCH2CO];
  w[r240f] = k[r240f] * c[sCH2CHO];
  w[r240b] = k[r240b] * c[sCO] * c[sCH3];
  w[r241f] = k[r241f] * c[sCH2CHO] * c[sO2];
  w[r241b] = k[r241b] * c[sHO2] * c[sCH2CO];
  w[r242] = k[r242] * c[sCH2CHO] * c[sO2];
  w[r243f] = k[r243f] * c[sCH2] * c[sCO];
  w[r243b] = k[r243b] * c[sCH2CO];
  w[r244f] = k[r244f] * c[sCH3CO];
  w[r244b] = k[r244b] * c[sH] * c[sCH2CO];
  w[r245f] = k[r245f] * c[sCH2CO] * c[sH];
  w[r245b] = k[r245b] * c[sH2] * c[sHCCO];
  w[r246f] = k[r246f] * c[sCH2CO] * c[sH];
  w[r246b] = k[r246b] * c[sCO] * c[sCH3];
  w[r247f] = k[r247f] * c[sCH2CO] * c[sO];
  w[r247b] = k[r247b] * c[sCO2] * c[sCH2];
  w[r248f] = k[r248f] * c[sCH2CO] * c[sO];
  w[r248b] = k[r248b] * c[sOH] * c[sHCCO];
  w[r249f] = k[r249f] * c[sCH2CO] * c[sOH];
  w[r249b] = k[r249b] * c[sH2O] * c[sHCCO];
  w[r250f] = k[r250f] * c[sCH2CO] * c[sOH];
  w[r250b] = k[r250b] * c[sCO] * c[sCH2OH];
  w[r251f] = k[r251f] * c[sCH2CO] * c[sCH3];
  w[r251b] = k[r251b] * c[sCO] * c[sC2H5];
  w[r252f] = k[r252f] * c[sCH2Y] * c[sCH2CO];
  w[r252b] = k[r252b] * c[sCO] * c[sC2H4];
  w[r253] = k[r253] * c[sHCCO] * c[sOH];
  w[r254] = k[r254] * c[sHCCO] * c[sO];
  w[r255f] = k[r255f] * c[sHCCO] * c[sH];
  w[r255b] = k[r255b] * c[sCO] * c[sCH2Y];
  w[r256] = k[r256] * c[sHCCO] * c[sO2];
  w[r257] = k[r257] * c[sHCCO] * c[sO2];
  w[r258f] = k[r258f] * c[sCH] * c[sCO] * M[mM33];
  w[r258b] = k[r258b] * c[sHCCO] * M[mM33];
  w[r259f] = k[r259f] * c[sCH] * c[sCH2O];
  w[r259b] = k[r259b] * c[sCH2CO] * c[sH];
  w[r260f] = k[r260f] * c[sCH] * c[sHCCO];
  w[r260b] = k[r260b] * c[sC2H2] * c[sCO];
  w[r261f] = k[r261f] * c[sC2H3] * c[sH];
  w[r261b] = k[r261b] * c[sC2H4];
  w[r262f] = k[r262f] * c[sC2H4];
  w[r262b] = k[r262b] * c[sH2CC] * c[sH2];
  w[r263f] = k[r263f] * c[sC2H4] * c[sH];
  w[r263b] = k[r263b] * c[sH2] * c[sC2H3];
  w[r264f] = k[r264f] * c[sC2H4] * c[sO];
  w[r264b] = k[r264b] * c[sHCO] * c[sCH3];
  w[r265f] = k[r265f] * c[sC2H4] * c[sO];
  w[r265b] = k[r265b] * c[sH] * c[sCH2CHO];
  w[r266f] = k[r266f] * c[sC2H4] * c[sOH];
  w[r266b] = k[r266b] * c[sH2O] * c[sC2H3];
  w[r267f] = k[r267f] * c[sC2H4] * c[sOH];
  w[r267b] = k[r267b] * c[sCH2O] * c[sCH3];
  w[r268f] = k[r268f] * c[sC2H4] * c[sOH];
  w[r268b] = k[r268b] * c[sH] * c[sCH3CHO];
  w[r269f] = k[r269f] * c[sC2H4] * c[sOH];
  w[r269b] = k[r269b] * c[sH] * c[sC2H3OH];
  w[r270f] = k[r270f] * c[sC2H3OH] * c[sO2];
  w[r270b] = k[r270b] * c[sHO2] * c[sCH2CHO];
  w[r271f] = k[r271f] * c[sC2H3OH] * c[sO];
  w[r271b] = k[r271b] * c[sOH] * c[sCH2CHO];
  w[r272f] = k[r272f] * c[sC2H3OH] * c[sOH];
  w[r272b] = k[r272b] * c[sH2O] * c[sCH2CHO];
  w[r273f] = k[r273f] * c[sC2H3OH] * c[sCH3];
  w[r273b] = k[r273b] * c[sCH4] * c[sCH2CHO];
  w[r274f] = k[r274f] * c[sC2H3OH] * c[sCH3O2];
  w[r274b] = k[r274b] * c[sCH3O2H] * c[sCH2CHO];
  w[r275f] = k[r275f] * c[sC2H3OH] * c[sH];
  w[r275b] = k[r275b] * c[sH2] * c[sCH2CHO];
  w[r276f] = k[r276f] * c[sC2H3OH] * c[sHO2];
  w[r276b] = k[r276b] * c[sHO2] * c[sCH3CHO];
  w[r277f] = k[r277f] * c[sC2H3OH];
  w[r277b] = k[r277b] * c[sCH3CHO];
  w[r278f] = k[r278f] * c[sC2H4] * c[sCH3];
  w[r278b] = k[r278b] * c[sCH4] * c[sC2H3];
  w[r279f] = k[r279f] * c[sC2H4] * c[sO2];
  w[r279b] = k[r279b] * c[sHO2] * c[sC2H3];
  w[r280f] = k[r280f] * c[sC2H4] * c[sCH3O];
  w[r280b] = k[r280b] * c[sCH3OH] * c[sC2H3];
  w[r281f] = k[r281f] * c[sC2H4] * c[sCH3O2];
  w[r281b] = k[r281b] * c[sCH3O2H] * c[sC2H3];
  w[r282f] = k[r282f] * c[sC2H4] * c[sC2H5O2];
  w[r282b] = k[r282b] * c[sC2H5O2H] * c[sC2H3];
  w[r283f] = k[r283f] * c[sC2H4] * c[sCH3CO3];
  w[r283b] = k[r283b] * c[sCH3CO3H] * c[sC2H3];
  w[r284f] = k[r284f] * c[sC2H4] * c[sCH3O2];
  w[r284b] = k[r284b] * c[sCH3O] * c[sC2H4O1X2];
  w[r285f] = k[r285f] * c[sC2H4] * c[sC2H5O2];
  w[r285b] = k[r285b] * c[sC2H5O] * c[sC2H4O1X2];
  w[r286f] = k[r286f] * c[sC2H4] * c[sHO2];
  w[r286b] = k[r286b] * c[sOH] * c[sC2H4O1X2];
  w[r287f] = k[r287f] * c[sCH] * c[sCH4];
  w[r287b] = k[r287b] * c[sH] * c[sC2H4];
  w[r288f] = k[r288f] * c[sCH2Y] * c[sCH3];
  w[r288b] = k[r288b] * c[sH] * c[sC2H4];
  w[r289f] = k[r289f] * c[sC2H2] * c[sH];
  w[r289b] = k[r289b] * c[sC2H3];
  w[r290f] = k[r290f] * c[sC2H3] * c[sO2];
  w[r290b] = k[r290b] * c[sHCO] * c[sCH2O];
  w[r291f] = k[r291f] * c[sC2H3] * c[sO2];
  w[r291b] = k[r291b] * c[sO] * c[sCH2CHO];
  w[r292] = k[r292] * c[sC2H3] * c[sO2];
  w[r293f] = k[r293f] * c[sCH3] * c[sC2H3];
  w[r293b] = k[r293b] * c[sC2H2] * c[sCH4];
  w[r294f] = k[r294f] * c[sC2H3] * c[sH];
  w[r294b] = k[r294b] * c[sH2] * c[sC2H2];
  w[r295f] = k[r295f] * c[sC2H3] * c[sH];
  w[r295b] = k[r295b] * c[sH2] * c[sH2CC];
  w[r296f] = k[r296f] * c[sC2H3] * c[sOH];
  w[r296b] = k[r296b] * c[sH2O] * c[sC2H2];
  w[r297f] = k[r297f] * c[sC2H3] * c[sC2H3];
  w[r297b] = k[r297b] * c[sC2H4] * c[sC2H2];
  w[r298f] = k[r298f] * c[sC2H] * c[sH];
  w[r298b] = k[r298b] * c[sC2H2];
  w[r299f] = k[r299f] * c[sC2H] * c[sO];
  w[r299b] = k[r299b] * c[sCO] * c[sCH];
  w[r300f] = k[r300f] * c[sC2H] * c[sOH];
  w[r300b] = k[r300b] * c[sHCCO] * c[sH];
  w[r301f] = k[r301f] * c[sC2H] * c[sO2];
  w[r301b] = k[r301b] * c[sCO] * c[sHCO];
  w[r302f] = k[r302f] * c[sC2H] * c[sH2];
  w[r302b] = k[r302b] * c[sC2H2] * c[sH];
  w[r303f] = k[r303f] * c[sC2H2];
  w[r303b] = k[r303b] * c[sH2CC];
  w[r304f] = k[r304f] * c[sC2H2] * c[sO];
  w[r304b] = k[r304b] * c[sCO] * c[sCH2];
  w[r305f] = k[r305f] * c[sC2H2] * c[sO];
  w[r305b] = k[r305b] * c[sH] * c[sHCCO];
  w[r306f] = k[r306f] * c[sC2H2] * c[sOH];
  w[r306b] = k[r306b] * c[sH2O] * c[sC2H];
  w[r307f] = k[r307f] * c[sC2H2] * c[sOH];
  w[r307b] = k[r307b] * c[sH] * c[sCH2CO];
  w[r308f] = k[r308f] * c[sC2H2] * c[sOH];
  w[r308b] = k[r308b] * c[sCO] * c[sCH3];
  w[r309f] = k[r309f] * c[sC2H2] * c[sHCO];
  w[r309b] = k[r309b] * c[sCO] * c[sC2H3];
  w[r310f] = k[r310f] * c[sH2CC] * c[sH];
  w[r310b] = k[r310b] * c[sH] * c[sC2H2];
  w[r311f] = k[r311f] * c[sH2CC] * c[sOH];
  w[r311b] = k[r311b] * c[sH] * c[sCH2CO];
  w[r312f] = k[r312f] * c[sH2CC] * c[sO2];
  w[r312b] = k[r312b] * c[sHCO] * c[sHCO];
  w[r313f] = k[r313f] * c[sC2H5OH];
  w[r313b] = k[r313b] * c[sH2O] * c[sC2H4];
  w[r314f] = k[r314f] * c[sC2H5OH];
  w[r314b] = k[r314b] * c[sCH2OH] * c[sCH3];
  w[r315f] = k[r315f] * c[sC2H5OH];
  w[r315b] = k[r315b] * c[sOH] * c[sC2H5];
  w[r316f] = k[r316f] * c[sC2H5OH] * c[sO2];
  w[r316b] = k[r316b] * c[sHO2] * c[sSC2H4OH];
  w[r317f] = k[r317f] * c[sC2H5OH] * c[sH];
  w[r317b] = k[r317b] * c[sH2] * c[sSC2H4OH];
  w[r318f] = k[r318f] * c[sC2H5OH] * c[sH];
  w[r318b] = k[r318b] * c[sH2] * c[sC2H5O];
  w[r319f] = k[r319f] * c[sC2H5OH] * c[sOH];
  w[r319b] = k[r319b] * c[sH2O] * c[sSC2H4OH];
  w[r320f] = k[r320f] * c[sC2H5OH] * c[sOH];
  w[r320b] = k[r320b] * c[sH2O] * c[sC2H5O];
  w[r321f] = k[r321f] * c[sC2H5OH] * c[sHO2];
  w[r321b] = k[r321b] * c[sH2O2] * c[sSC2H4OH];
  w[r322f] = k[r322f] * c[sC2H5OH] * c[sHO2];
  w[r322b] = k[r322b] * c[sH2O2] * c[sC2H5O];
  w[r323f] = k[r323f] * c[sC2H5OH] * c[sCH3O2];
  w[r323b] = k[r323b] * c[sCH3O2H] * c[sSC2H4OH];
  w[r324f] = k[r324f] * c[sC2H5OH] * c[sCH3O2];
  w[r324b] = k[r324b] * c[sCH3O2H] * c[sC2H5O];
  w[r325f] = k[r325f] * c[sC2H5OH] * c[sO];
  w[r325b] = k[r325b] * c[sOH] * c[sSC2H4OH];
  w[r326f] = k[r326f] * c[sC2H5OH] * c[sO];
  w[r326b] = k[r326b] * c[sOH] * c[sC2H5O];
  w[r327f] = k[r327f] * c[sC2H5OH] * c[sCH3];
  w[r327b] = k[r327b] * c[sCH4] * c[sSC2H4OH];
  w[r328f] = k[r328f] * c[sC2H5OH] * c[sCH3];
  w[r328b] = k[r328b] * c[sCH4] * c[sC2H5O];
  w[r329f] = k[r329f] * c[sC2H5OH] * c[sC2H5];
  w[r329b] = k[r329b] * c[sC2H6] * c[sSC2H4OH];
  w[r330f] = k[r330f] * c[sSC2H4OH];
  w[r330b] = k[r330b] * c[sH] * c[sCH3CHO];
  w[r331f] = k[r331f] * c[sSC2H4OH];
  w[r331b] = k[r331b] * c[sH] * c[sC2H3OH];
  w[r332f] = k[r332f] * c[sSC2H4OH];
  w[r332b] = k[r332b] * c[sC2H5O];
  w[r333f] = k[r333f] * c[sSC2H4OH] * c[sO2];
  w[r333b] = k[r333b] * c[sHO2] * c[sCH3CHO];
  w[r334f] = k[r334f] * c[sSC2H4OH] * c[sO2];
  w[r334b] = k[r334b] * c[sHO2] * c[sC2H3OH];
  w[r335f] = k[r335f] * c[sCH3COCH3];
  w[r335b] = k[r335b] * c[sCH3] * c[sCH3CO];
  w[r336f] = k[r336f] * c[sCH3COCH3] * c[sOH];
  w[r336b] = k[r336b] * c[sH2O] * c[sCH3COCH2];
  w[r337f] = k[r337f] * c[sCH3COCH3] * c[sH];
  w[r337b] = k[r337b] * c[sH2] * c[sCH3COCH2];
  w[r338f] = k[r338f] * c[sCH3COCH3] * c[sO];
  w[r338b] = k[r338b] * c[sOH] * c[sCH3COCH2];
  w[r339f] = k[r339f] * c[sCH3COCH3] * c[sCH3];
  w[r339b] = k[r339b] * c[sCH4] * c[sCH3COCH2];
  w[r340f] = k[r340f] * c[sCH3COCH3] * c[sCH3O];
  w[r340b] = k[r340b] * c[sCH3OH] * c[sCH3COCH2];
  w[r341f] = k[r341f] * c[sCH3COCH3] * c[sO2];
  w[r341b] = k[r341b] * c[sHO2] * c[sCH3COCH2];
  w[r342f] = k[r342f] * c[sCH3COCH3] * c[sHO2];
  w[r342b] = k[r342b] * c[sH2O2] * c[sCH3COCH2];
  w[r343f] = k[r343f] * c[sCH3COCH3] * c[sCH3O2];
  w[r343b] = k[r343b] * c[sCH3O2H] * c[sCH3COCH2];
  w[r344f] = k[r344f] * c[sCH2CO] * c[sCH3];
  w[r344b] = k[r344b] * c[sCH3COCH2];
  w[r345f] = k[r345f] * c[sC2H3] * c[sHCO];
  w[r345b] = k[r345b] * c[sC2H3CHO];
  w[r346f] = k[r346f] * c[sC2H3CHO] * c[sH];
  w[r346b] = k[r346b] * c[sH2] * c[sC2H3CO];
  w[r347f] = k[r347f] * c[sC2H3CHO] * c[sO];
  w[r347b] = k[r347b] * c[sOH] * c[sC2H3CO];
  w[r348f] = k[r348f] * c[sC2H3CHO] * c[sOH];
  w[r348b] = k[r348b] * c[sH2O] * c[sC2H3CO];
  w[r349f] = k[r349f] * c[sC2H3CHO] * c[sO2];
  w[r349b] = k[r349b] * c[sHO2] * c[sC2H3CO];
  w[r350f] = k[r350f] * c[sC2H3CHO] * c[sHO2];
  w[r350b] = k[r350b] * c[sH2O2] * c[sC2H3CO];
  w[r351f] = k[r351f] * c[sC2H3CHO] * c[sCH3];
  w[r351b] = k[r351b] * c[sCH4] * c[sC2H3CO];
  w[r352f] = k[r352f] * c[sC2H3CHO] * c[sC2H3];
  w[r352b] = k[r352b] * c[sC2H4] * c[sC2H3CO];
  w[r353f] = k[r353f] * c[sC2H3CHO] * c[sCH3O];
  w[r353b] = k[r353b] * c[sCH3OH] * c[sC2H3CO];
  w[r354f] = k[r354f] * c[sC2H3CHO] * c[sCH3O2];
  w[r354b] = k[r354b] * c[sCH3O2H] * c[sC2H3CO];
  w[r355f] = k[r355f] * c[sC2H3] * c[sCO];
  w[r355b] = k[r355b] * c[sC2H3CO];
  w[r356f] = k[r356f] * c[sC2H5] * c[sHCO];
  w[r356b] = k[r356b] * c[sC2H5CHO];
  w[r357f] = k[r357f] * c[sC2H5CHO] * c[sH];
  w[r357b] = k[r357b] * c[sH2] * c[sC2H5CO];
  w[r358f] = k[r358f] * c[sC2H5CHO] * c[sO];
  w[r358b] = k[r358b] * c[sOH] * c[sC2H5CO];
  w[r359f] = k[r359f] * c[sC2H5CHO] * c[sOH];
  w[r359b] = k[r359b] * c[sH2O] * c[sC2H5CO];
  w[r360f] = k[r360f] * c[sC2H5CHO] * c[sCH3];
  w[r360b] = k[r360b] * c[sCH4] * c[sC2H5CO];
  w[r361f] = k[r361f] * c[sC2H5CHO] * c[sHO2];
  w[r361b] = k[r361b] * c[sH2O2] * c[sC2H5CO];
  w[r362f] = k[r362f] * c[sC2H5CHO] * c[sCH3O];
  w[r362b] = k[r362b] * c[sCH3OH] * c[sC2H5CO];
  w[r363f] = k[r363f] * c[sC2H5CHO] * c[sCH3O2];
  w[r363b] = k[r363b] * c[sCH3O2H] * c[sC2H5CO];
  w[r364f] = k[r364f] * c[sC2H5CHO] * c[sC2H5];
  w[r364b] = k[r364b] * c[sC2H6] * c[sC2H5CO];
  w[r365f] = k[r365f] * c[sC2H5CHO] * c[sC2H5O];
  w[r365b] = k[r365b] * c[sC2H5OH] * c[sC2H5CO];
  w[r366f] = k[r366f] * c[sC2H5CHO] * c[sC2H5O2];
  w[r366b] = k[r366b] * c[sC2H5O2H] * c[sC2H5CO];
  w[r367f] = k[r367f] * c[sC2H5CHO] * c[sO2];
  w[r367b] = k[r367b] * c[sHO2] * c[sC2H5CO];
  w[r368f] = k[r368f] * c[sC2H5CHO] * c[sCH3CO3];
  w[r368b] = k[r368b] * c[sCH3CO3H] * c[sC2H5CO];
  w[r369f] = k[r369f] * c[sC2H5CHO] * c[sC2H3];
  w[r369b] = k[r369b] * c[sC2H4] * c[sC2H5CO];
  w[r370f] = k[r370f] * c[sC2H5] * c[sCO];
  w[r370b] = k[r370b] * c[sC2H5CO];
  w[r371f] = k[r371f] * c[sCH3OCH3];
  w[r371b] = k[r371b] * c[sCH3O] * c[sCH3];
  w[r372f] = k[r372f] * c[sCH3OCH3] * c[sOH];
  w[r372b] = k[r372b] * c[sH2O] * c[sCH3OCH2];
  w[r373f] = k[r373f] * c[sCH3OCH3] * c[sH];
  w[r373b] = k[r373b] * c[sH2] * c[sCH3OCH2];
  w[r374f] = k[r374f] * c[sCH3OCH3] * c[sO];
  w[r374b] = k[r374b] * c[sOH] * c[sCH3OCH2];
  w[r375f] = k[r375f] * c[sCH3OCH3] * c[sHO2];
  w[r375b] = k[r375b] * c[sH2O2] * c[sCH3OCH2];
  w[r376f] = k[r376f] * c[sCH3OCH3] * c[sCH3O2];
  w[r376b] = k[r376b] * c[sCH3O2H] * c[sCH3OCH2];
  w[r377f] = k[r377f] * c[sCH3OCH3] * c[sCH3];
  w[r377b] = k[r377b] * c[sCH4] * c[sCH3OCH2];
  w[r378f] = k[r378f] * c[sCH3OCH3] * c[sO2];
  w[r378b] = k[r378b] * c[sHO2] * c[sCH3OCH2];
  w[r379f] = k[r379f] * c[sCH3OCH3] * c[sCH3O];
  w[r379b] = k[r379b] * c[sCH3OH] * c[sCH3OCH2];
  w[r380f] = k[r380f] * c[sCH3OCH3] * c[sCH3OCH2O2];
  w[r380b] = k[r380b] * c[sCH3OCH2O2H] * c[sCH3OCH2];
  w[r381f] = k[r381f] * c[sCH3OCH3] * c[sO2CHO];
  w[r381b] = k[r381b] * c[sHO2CHO] * c[sCH3OCH2];
  w[r382f] = k[r382f] * c[sCH3OCH3] * c[sOCHO];
  w[r382b] = k[r382b] * c[sHOCHO] * c[sCH3OCH2];
  w[r383f] = k[r383f] * c[sCH3OCH2];
  w[r383b] = k[r383b] * c[sCH3] * c[sCH2O];
  w[r384f] = k[r384f] * c[sCH3OCH2] * c[sCH3O];
  w[r384b] = k[r384b] * c[sCH2O] * c[sCH3OCH3];
  w[r385f] = k[r385f] * c[sCH3OCH2] * c[sCH2O];
  w[r385b] = k[r385b] * c[sHCO] * c[sCH3OCH3];
  w[r386f] = k[r386f] * c[sCH3OCH2] * c[sCH3CHO];
  w[r386b] = k[r386b] * c[sCH3CO] * c[sCH3OCH3];
  w[r387f] = k[r387f] * c[sCH3OCH2] * c[sO2];
  w[r387b] = k[r387b] * c[sCH3OCH2O2];
  w[r388f] = k[r388f] * c[sCH3OCH2O2] * c[sCH2O];
  w[r388b] = k[r388b] * c[sHCO] * c[sCH3OCH2O2H];
  w[r389f] = k[r389f] * c[sCH3OCH2O2] * c[sCH3CHO];
  w[r389b] = k[r389b] * c[sCH3CO] * c[sCH3OCH2O2H];
  w[r390] = k[r390] * c[sCH3OCH2O2] * c[sCH3OCH2O2];
  w[r391f] = k[r391f] * c[sCH3OCH2O] * c[sOH];
  w[r391b] = k[r391b] * c[sCH3OCH2O2H];
  w[r392f] = k[r392f] * c[sCH3O] * c[sCH2O];
  w[r392b] = k[r392b] * c[sCH3OCH2O];
  w[r393f] = k[r393f] * c[sCH3OCH2O] * c[sO2];
  w[r393b] = k[r393b] * c[sHO2] * c[sCH3OCHO];
  w[r394f] = k[r394f] * c[sCH3OCHO] * c[sH];
  w[r394b] = k[r394b] * c[sCH3OCH2O];
  w[r395f] = k[r395f] * c[sCH3OCH2O2];
  w[r395b] = k[r395b] * c[sCH2OCH2O2H];
  w[r396] = k[r396] * c[sCH2OCH2O2H];
  w[r397f] = k[r397f] * c[sCH2OCH2O2H] * c[sO2];
  w[r397b] = k[r397b] * c[sO2CH2OCH2O2H];
  w[r398f] = k[r398f] * c[sO2CH2OCH2O2H];
  w[r398b] = k[r398b] * c[sOH] * c[sHO2CH2OCHO];
  w[r399f] = k[r399f] * c[sHO2CH2OCHO];
  w[r399b] = k[r399b] * c[sOH] * c[sOCH2OCHO];
  w[r400f] = k[r400f] * c[sCH2O] * c[sOCHO];
  w[r400b] = k[r400b] * c[sOCH2OCHO];
  w[r401f] = k[r401f] * c[sOCH2OCHO];
  w[r401b] = k[r401b] * c[sHOCH2OCO];
  w[r402f] = k[r402f] * c[sHOCH2O] * c[sCO];
  w[r402b] = k[r402b] * c[sHOCH2OCO];
  w[r403f] = k[r403f] * c[sCH2OH] * c[sCO2];
  w[r403b] = k[r403b] * c[sHOCH2OCO];
  w[r404f] = k[r404f] * c[sCH2OCHO] * c[sH];
  w[r404b] = k[r404b] * c[sCH3OCHO];
  w[r405f] = k[r405f] * c[sCH3OCO] * c[sH];
  w[r405b] = k[r405b] * c[sCH3OCHO];
  w[r406f] = k[r406f] * c[sCH3OCHO];
  w[r406b] = k[r406b] * c[sCO] * c[sCH3OH];
  w[r407f] = k[r407f] * c[sCH3O] * c[sHCO];
  w[r407b] = k[r407b] * c[sCH3OCHO];
  w[r408f] = k[r408f] * c[sCH3] * c[sOCHO];
  w[r408b] = k[r408b] * c[sCH3OCHO];
  w[r409f] = k[r409f] * c[sCH3OCHO] * c[sO2];
  w[r409b] = k[r409b] * c[sHO2] * c[sCH3OCO];
  w[r410f] = k[r410f] * c[sCH3OCHO] * c[sO2];
  w[r410b] = k[r410b] * c[sHO2] * c[sCH2OCHO];
  w[r411f] = k[r411f] * c[sCH3OCHO] * c[sOH];
  w[r411b] = k[r411b] * c[sH2O] * c[sCH3OCO];
  w[r412f] = k[r412f] * c[sCH3OCHO] * c[sOH];
  w[r412b] = k[r412b] * c[sH2O] * c[sCH2OCHO];
  w[r413f] = k[r413f] * c[sCH3OCHO] * c[sHO2];
  w[r413b] = k[r413b] * c[sH2O2] * c[sCH3OCO];
  w[r414f] = k[r414f] * c[sCH3OCHO] * c[sHO2];
  w[r414b] = k[r414b] * c[sH2O2] * c[sCH2OCHO];
  w[r415f] = k[r415f] * c[sCH3OCHO] * c[sO];
  w[r415b] = k[r415b] * c[sOH] * c[sCH3OCO];
  w[r416f] = k[r416f] * c[sCH3OCHO] * c[sO];
  w[r416b] = k[r416b] * c[sOH] * c[sCH2OCHO];
  w[r417f] = k[r417f] * c[sCH3OCHO] * c[sH];
  w[r417b] = k[r417b] * c[sH2] * c[sCH3OCO];
  w[r418f] = k[r418f] * c[sCH3OCHO] * c[sH];
  w[r418b] = k[r418b] * c[sH2] * c[sCH2OCHO];
  w[r419f] = k[r419f] * c[sCH3OCHO] * c[sCH3];
  w[r419b] = k[r419b] * c[sCH4] * c[sCH3OCO];
  w[r420f] = k[r420f] * c[sCH3OCHO] * c[sCH3];
  w[r420b] = k[r420b] * c[sCH4] * c[sCH2OCHO];
  w[r421f] = k[r421f] * c[sCH3OCHO] * c[sCH3O];
  w[r421b] = k[r421b] * c[sCH3OH] * c[sCH3OCO];
  w[r422f] = k[r422f] * c[sCH3OCHO] * c[sCH3O];
  w[r422b] = k[r422b] * c[sCH3OH] * c[sCH2OCHO];
  w[r423f] = k[r423f] * c[sCH3OCHO] * c[sCH3O2];
  w[r423b] = k[r423b] * c[sCH3O2H] * c[sCH3OCO];
  w[r424f] = k[r424f] * c[sCH3OCHO] * c[sCH3O2];
  w[r424b] = k[r424b] * c[sCH3O2H] * c[sCH2OCHO];
  w[r425f] = k[r425f] * c[sCH3OCHO] * c[sHCO];
  w[r425b] = k[r425b] * c[sCH2O] * c[sCH3OCO];
  w[r426f] = k[r426f] * c[sCH3OCHO] * c[sHCO];
  w[r426b] = k[r426b] * c[sCH2O] * c[sCH2OCHO];
  w[r427f] = k[r427f] * c[sCH3OCO];
  w[r427b] = k[r427b] * c[sCH2OCHO];
  w[r428f] = k[r428f] * c[sCH3] * c[sCO2];
  w[r428b] = k[r428b] * c[sCH3OCO];
  w[r429f] = k[r429f] * c[sCH3O] * c[sCO];
  w[r429b] = k[r429b] * c[sCH3OCO];
  w[r430f] = k[r430f] * c[sCH2O] * c[sHCO];
  w[r430b] = k[r430b] * c[sCH2OCHO];

  cdot[sN2] = -w[r47f] + w[r47f] - w[r47b] + w[r47b] - w[r148f] + w[r148f] -
              w[r148b] + w[r148b];

  cdot[sCH3OCH3] = -w[r371f] + w[r371b] - w[r372f] + w[r372b] - w[r373f] +
                   w[r373b] - w[r374f] + w[r374b] - w[r375f] + w[r375b] -
                   w[r376f] + w[r376b] - w[r377f] + w[r377b] - w[r378f] +
                   w[r378b] - w[r379f] + w[r379b] - w[r380f] + w[r380b] -
                   w[r381f] + w[r381b] - w[r382f] + w[r382b] + w[r384f] -
                   w[r384b] + w[r385f] - w[r385b] + w[r386f] - w[r386b];

  cdot[sCH4] = w[r38f] - w[r38b] - w[r55f] + w[r55f] - w[r55b] + w[r55b] +
               w[r62f] - w[r62b] + w[r73] + w[r81f] - w[r81b] + w[r82f] -
               w[r82b] + w[r111f] - w[r111b] + w[r114f] - w[r114b] - w[r115f] +
               w[r115b] - w[r116f] + w[r116b] - w[r117f] + w[r117b] - w[r118f] +
               w[r118b] - w[r119f] + w[r119b] + w[r132f] - w[r132b] - w[r138f] +
               w[r138b] + w[r179f] - w[r179b] + w[r189f] - w[r189b] - w[r199f] +
               w[r199b] + w[r215f] - w[r215b] + w[r221f] - w[r221b] + w[r230f] -
               w[r230b] - w[r234f] + w[r234b] + w[r273f] - w[r273b] + w[r278f] -
               w[r278b] - w[r287f] + w[r287b] + w[r293f] - w[r293b] + w[r327f] -
               w[r327b] + w[r328f] - w[r328b] + w[r339f] - w[r339b] + w[r351f] -
               w[r351b] + w[r360f] - w[r360b] + w[r377f] - w[r377b] + w[r419f] -
               w[r419b] + w[r420f] - w[r420b];

  cdot[sH2] = -w[r2f] + w[r2b] - w[r3f] + w[r3b] - w[r5f] + w[r5b] - w[r13f] +
              w[r13b] + w[r21f] - w[r21b] + w[r32f] - w[r32b] + w[r37] -
              w[r46f] + w[r46f] - w[r46b] + w[r46b] - w[r58f] + w[r58b] +
              w[r60f] - w[r60b] + w[r67f] - w[r67b] + w[r71] + w[r72] +
              w[r83f] - w[r83b] + w[r88f] - w[r88b] + w[r104f] - w[r104b] +
              w[r105f] - w[r105b] + w[r115f] - w[r115b] + w[r121f] - w[r121b] +
              w[r124f] - w[r124b] + w[r150f] - w[r150b] + w[r151f] - w[r151b] -
              w[r154f] + w[r154b] + w[r165f] - w[r165b] + w[r166f] - w[r166b] +
              w[r175f] - w[r175b] - w[r186f] + w[r186b] - w[r187f] + w[r187b] +
              w[r191f] - w[r191b] + w[r216f] - w[r216b] + w[r217f] - w[r217b] +
              w[r228f] - w[r228b] + w[r245f] - w[r245b] + w[r253] + w[r262f] -
              w[r262b] + w[r263f] - w[r263b] + w[r275f] - w[r275b] + w[r294f] -
              w[r294b] + w[r295f] - w[r295b] - w[r302f] + w[r302b] + w[r317f] -
              w[r317b] + w[r318f] - w[r318b] + w[r337f] - w[r337b] + w[r346f] -
              w[r346b] + w[r357f] - w[r357b] + w[r373f] - w[r373b] + w[r417f] -
              w[r417b] + w[r418f] - w[r418b];

  cdot[sO2] = -w[r1f] + w[r1b] + w[r6f] - w[r6b] - w[r9f] + w[r9b] - w[r10f] +
              w[r10b] - w[r11f] + w[r11b] - w[r13f] + w[r13b] + w[r14f] -
              w[r14b] + w[r15f] - w[r15b] + w[r16f] - w[r16b] + w[r17f] -
              w[r17b] - w[r26f] + w[r26b] - w[r31f] + w[r31b] - w[r39f] +
              w[r39b] - w[r40f] + w[r40b] - w[r52f] + w[r52f] - w[r52b] +
              w[r52b] - w[r56f] + w[r56b] + w[r68f] - w[r68b] - w[r79f] +
              w[r79b] - w[r86f] + w[r86b] - w[r87f] + w[r87b] + w[r99f] -
              w[r99b] - w[r109f] + w[r109b] - w[r129] - w[r130f] + w[r130b] +
              w[r132f] - w[r132b] - w[r134f] + w[r134b] - w[r135f] + w[r135b] -
              w[r136f] + w[r136b] + w[r141f] - w[r141b] + w[r142] + w[r143] +
              w[r145f] - w[r145b] + w[r146f] - w[r146b] - w[r155] - w[r156f] +
              w[r156b] - w[r162f] + w[r162b] - w[r163] - w[r168f] + w[r168b] -
              w[r178f] + w[r178b] - w[r195f] + w[r195b] + w[r201f] - w[r201b] -
              w[r204f] + w[r204b] - w[r205f] + w[r205b] - w[r206f] + w[r206b] -
              w[r207f] + w[r207b] - w[r208f] + w[r208b] - w[r220f] + w[r220b] -
              w[r231f] + w[r231b] + w[r232f] - w[r232b] - w[r241f] + w[r241b] -
              w[r242] - w[r256] - w[r257] - w[r270f] + w[r270b] - w[r279f] +
              w[r279b] - w[r290f] + w[r290b] - w[r291f] + w[r291b] - w[r292] -
              w[r301f] + w[r301b] - w[r312f] + w[r312b] - w[r316f] + w[r316b] -
              w[r333f] + w[r333b] - w[r334f] + w[r334b] - w[r341f] + w[r341b] -
              w[r349f] + w[r349b] - w[r367f] + w[r367b] - w[r378f] + w[r378b] -
              w[r387f] + w[r387b] + w[r390] - w[r393f] + w[r393b] - w[r397f] +
              w[r397b] - w[r409f] + w[r409b] - w[r410f] + w[r410b];

  cdot[sCO] = -w[r25f] + w[r25b] - w[r26f] + w[r26b] - w[r27f] + w[r27b] -
              w[r28f] + w[r28b] - w[r29f] + w[r29b] + w[r30f] - w[r30b] +
              w[r31f] - w[r31b] + w[r32f] - w[r32b] + w[r33f] - w[r33b] +
              w[r35f] - w[r35b] + 0x1p+1 * w[r37] + w[r38f] - w[r38b] +
              w[r44f] - w[r44b] - w[r54f] + w[r54f] - w[r54b] + w[r54b] +
              w[r56f] - w[r56b] - w[r58f] + w[r58b] + w[r66f] - w[r66b] +
              w[r70] + w[r72] + w[r73] + w[r75] + w[r76] + w[r128] + w[r151f] -
              w[r151b] + w[r155] + w[r156f] - w[r156b] - w[r158f] + w[r158f] -
              w[r158b] + w[r158b] + w[r160f] - w[r160b] + w[r164] + w[r169f] -
              w[r169b] + w[r172f] - w[r172b] + w[r215f] - w[r215b] + w[r227f] -
              w[r227b] + w[r240f] - w[r240b] + w[r242] - w[r243f] + w[r243b] +
              w[r246f] - w[r246b] + w[r250f] - w[r250b] + w[r251f] - w[r251b] +
              w[r252f] - w[r252b] + 0x1p+1 * w[r253] + 0x1p+1 * w[r254] +
              w[r255f] - w[r255b] + 0x1p+1 * w[r256] + w[r257] - w[r258f] +
              w[r258b] + w[r260f] - w[r260b] + w[r292] + w[r299f] - w[r299b] +
              w[r301f] - w[r301b] + w[r304f] - w[r304b] + w[r308f] - w[r308b] +
              w[r309f] - w[r309b] - w[r355f] + w[r355b] - w[r370f] + w[r370b] -
              w[r402f] + w[r402b] + w[r406f] - w[r406b] - w[r429f] + w[r429b];

  cdot[sCO2] = w[r25f] - w[r25b] + w[r26f] - w[r26b] + w[r27f] - w[r27b] +
               w[r28f] - w[r28b] + w[r29f] - w[r29b] + w[r34f] - w[r34b] +
               w[r36] - w[r43f] + w[r43b] - w[r53f] + w[r53f] - w[r53b] +
               w[r53b] + w[r67f] - w[r67b] + w[r69] + w[r71] + w[r127] +
               w[r129] + w[r130f] - w[r130b] - w[r159f] + w[r159f] - w[r159b] +
               w[r159b] - w[r160f] + w[r160b] + w[r163] - w[r172f] + w[r172b] +
               w[r238f] - w[r238b] + w[r247f] - w[r247b] + w[r257] - w[r403f] +
               w[r403b] - w[r428f] + w[r428b];

  cdot[sH2O] = w[r3f] - w[r3b] - w[r4f] + w[r4b] + w[r8f] - w[r8b] + w[r15f] -
               w[r15b] + w[r20f] - w[r20b] + w[r23f] - w[r23b] + w[r24f] -
               w[r24b] + w[r35f] - w[r35b] - w[r45f] + w[r45f] - w[r45b] +
               w[r45b] + w[r59f] - w[r59b] + w[r66f] - w[r66b] + w[r69] +
               w[r70] + w[r93f] - w[r93b] + w[r102f] - w[r102b] + w[r107f] -
               w[r107b] + w[r108f] - w[r108b] + w[r116f] - w[r116b] + w[r120f] -
               w[r120b] + w[r125f] - w[r125b] + w[r130f] - w[r130b] + w[r156f] -
               w[r156b] - w[r157f] + w[r157f] - w[r157b] + w[r157b] + w[r167f] -
               w[r167b] - w[r171f] + w[r171b] + w[r177f] - w[r177b] + w[r219f] -
               w[r219b] + w[r226f] - w[r226b] + w[r249f] - w[r249b] + w[r266f] -
               w[r266b] + w[r272f] - w[r272b] + w[r296f] - w[r296b] + w[r306f] -
               w[r306b] + w[r313f] - w[r313b] + w[r319f] - w[r319b] + w[r320f] -
               w[r320b] + w[r336f] - w[r336b] + w[r348f] - w[r348b] + w[r359f] -
               w[r359b] + w[r372f] - w[r372b] + w[r411f] - w[r411b] + w[r412f] -
               w[r412b];

  cdot[sH] =
      -w[r1f] + w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] + 0x1p+1 * w[r5f] -
      0x1p+1 * w[r5b] - w[r7f] + w[r7b] - w[r8f] + w[r8b] - w[r9f] + w[r9b] -
      w[r10f] + w[r10b] - w[r11f] + w[r11b] - w[r12f] + w[r12b] + w[r13f] -
      w[r13b] - w[r20f] + w[r20b] - w[r21f] + w[r21b] + w[r27f] - w[r27b] +
      w[r28f] - w[r28b] + w[r30f] - w[r30b] - w[r32f] + w[r32b] + w[r34f] -
      w[r34b] + w[r36] - w[r43f] + w[r43b] - w[r49f] + w[r49f] - w[r49b] +
      w[r49b] - w[r57f] + w[r57b] - w[r60f] + w[r60b] + w[r65f] - w[r65b] +
      w[r69] - w[r71] + w[r71] - w[r72] + w[r78f] - w[r78b] - w[r83f] +
      w[r83b] - w[r85f] + w[r85b] - w[r88f] + w[r88b] + w[r103f] - w[r103b] -
      w[r104f] + w[r104b] - w[r105f] + w[r105b] - w[r114f] + w[r114b] -
      w[r115f] + w[r115b] + w[r122f] - w[r122b] + w[r123f] - w[r123b] -
      w[r126f] + w[r126f] - w[r126b] + w[r126b] + 0x1p+1 * w[r127] + w[r128] +
      w[r129] + w[r133f] - w[r133b] - w[r144f] + w[r144b] - w[r150f] +
      w[r150b] + w[r152f] - w[r152b] + w[r153f] - w[r153b] + w[r154f] -
      w[r154b] + w[r155] - w[r161f] + w[r161b] + 0x1p+1 * w[r163] +
      0x1p+1 * w[r164] - w[r165f] + w[r165b] - w[r166f] + w[r166b] + w[r169f] -
      w[r169b] + w[r170f] - w[r170b] + w[r171f] - w[r171b] - w[r174f] +
      w[r174b] - w[r175f] + w[r175b] - w[r185f] + w[r185b] + w[r186f] -
      w[r186b] + w[r187f] - w[r187b] + w[r190f] - w[r190b] - w[r191f] +
      w[r191b] + w[r192f] - w[r192b] - w[r197f] + w[r197b] - w[r216f] +
      w[r216b] - w[r217f] + w[r217b] - w[r228f] + w[r228b] + w[r239f] -
      w[r239b] + w[r244f] - w[r244b] - w[r245f] + w[r245b] - w[r246f] +
      w[r246b] + w[r254] - w[r255f] + w[r255b] + w[r257] + w[r259f] - w[r259b] -
      w[r261f] + w[r261b] - w[r263f] + w[r263b] + w[r265f] - w[r265b] +
      w[r268f] - w[r268b] + w[r269f] - w[r269b] - w[r275f] + w[r275b] +
      w[r287f] - w[r287b] + w[r288f] - w[r288b] - w[r289f] + w[r289b] +
      w[r292] - w[r294f] + w[r294b] - w[r295f] + w[r295b] - w[r298f] +
      w[r298b] + w[r300f] - w[r300b] + w[r302f] - w[r302b] + w[r305f] -
      w[r305b] + w[r307f] - w[r307b] - w[r310f] + w[r310f] - w[r310b] +
      w[r310b] + w[r311f] - w[r311b] - w[r317f] + w[r317b] - w[r318f] +
      w[r318b] + w[r330f] - w[r330b] + w[r331f] - w[r331b] - w[r337f] +
      w[r337b] - w[r346f] + w[r346b] - w[r357f] + w[r357b] - w[r373f] +
      w[r373b] - w[r394f] + w[r394b] - w[r404f] + w[r404b] - w[r405f] +
      w[r405b] - w[r417f] + w[r417b] - w[r418f] + w[r418b];

  cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b] - w[r4f] + w[r4b] -
             0x1p+1 * w[r6f] + 0x1p+1 * w[r6b] - w[r7f] + w[r7b] - w[r14f] +
             w[r14b] - w[r22f] + w[r22b] - w[r25f] + w[r25b] + w[r26f] -
             w[r26b] - w[r33f] + w[r33b] - w[r34f] + w[r34b] - w[r61f] +
             w[r61b] - w[r76] - w[r94f] + w[r94b] - w[r106f] + w[r106b] -
             w[r117f] + w[r117b] - w[r127] - w[r128] - w[r133f] + w[r133b] +
             w[r134f] - w[r134b] - w[r145f] + w[r145b] - w[r151f] + w[r151b] -
             w[r152f] + w[r152b] - w[r164] + w[r168f] - w[r168b] - w[r169f] +
             w[r169b] - w[r176f] + w[r176b] - w[r192f] + w[r192b] - w[r218f] +
             w[r218b] - w[r229f] + w[r229b] - w[r247f] + w[r247b] - w[r248f] +
             w[r248b] - w[r254] - w[r264f] + w[r264b] - w[r265f] + w[r265b] -
             w[r271f] + w[r271b] + w[r291f] - w[r291b] - w[r299f] + w[r299b] -
             w[r304f] + w[r304b] - w[r305f] + w[r305b] - w[r325f] + w[r325b] -
             w[r326f] + w[r326b] - w[r338f] + w[r338b] - w[r347f] + w[r347b] -
             w[r358f] + w[r358b] - w[r374f] + w[r374b] - w[r415f] + w[r415b] -
             w[r416f] + w[r416b];

  cdot[sOH] =
      w[r1f] - w[r1b] + w[r2f] - w[r2b] - w[r3f] + w[r3b] + 0x1p+1 * w[r4f] -
      0x1p+1 * w[r4b] + w[r7f] - w[r7b] - w[r8f] + w[r8b] + 0x1p+1 * w[r12f] -
      0x1p+1 * w[r12b] + w[r14f] - w[r14b] - w[r15f] + w[r15b] +
      0x1p+1 * w[r18f] - 0x1p+1 * w[r18b] + 0x1p+1 * w[r19f] -
      0x1p+1 * w[r19b] + w[r20f] - w[r20b] + w[r22f] - w[r22b] - w[r23f] +
      w[r23b] - w[r24f] + w[r24b] - w[r27f] + w[r27b] - w[r28f] + w[r28b] +
      w[r29f] - w[r29b] + w[r33f] - w[r33b] - w[r35f] + w[r35b] + w[r36] -
      w[r42f] + w[r42b] + w[r45f] - w[r45b] + w[r46f] - w[r46b] + w[r47f] -
      w[r47b] - w[r48f] + 0x1p+1 * w[r48f] - 0x1p+1 * w[r48b] + w[r48b] +
      w[r49f] - w[r49b] + w[r50f] - w[r50b] + w[r51f] - w[r51b] + w[r52f] -
      w[r52b] + w[r53f] - w[r53b] + w[r54f] - w[r54b] + w[r55f] - w[r55b] -
      w[r59f] + w[r59b] + w[r61f] - w[r61b] - w[r64f] + w[r64b] - w[r69] -
      w[r70] + w[r70] + w[r72] + w[r73] + w[r75] + 0x1p+1 * w[r76] - w[r93f] +
      w[r93b] + w[r94f] - w[r94b] + w[r96f] - w[r96b] - w[r100f] + w[r100b] +
      w[r101f] - w[r101b] + w[r106f] - w[r106b] - w[r107f] + w[r107b] -
      w[r108f] + w[r108b] - w[r116f] + w[r116b] + w[r117f] - w[r117b] -
      w[r120f] + w[r120b] - w[r121f] + w[r121b] - w[r122f] + w[r122b] -
      w[r123f] + w[r123b] - w[r124f] + w[r124b] - w[r125f] + w[r125b] +
      w[r128] + w[r129] + w[r131f] - w[r131b] + w[r135f] - w[r135b] + w[r144f] -
      w[r144b] - w[r146f] + w[r146b] + w[r147f] - w[r147b] - w[r153f] +
      w[r153b] + w[r155] + w[r162f] - w[r162b] - w[r167f] + w[r167b] -
      w[r170f] + w[r170b] + w[r176f] - w[r176b] - w[r177f] + w[r177b] +
      w[r193f] - w[r193b] + w[r203f] - w[r203b] + w[r207f] - w[r207b] +
      w[r208f] - w[r208b] + w[r209f] - w[r209b] + w[r211f] - w[r211b] +
      w[r218f] - w[r218b] - w[r219f] + w[r219b] - w[r225f] + w[r225b] -
      w[r226f] + w[r226b] + w[r229f] - w[r229b] + w[r237f] - w[r237b] +
      w[r242] + w[r248f] - w[r248b] - w[r249f] + w[r249b] - w[r250f] +
      w[r250b] - w[r253] + w[r256] - w[r266f] + w[r266b] - w[r267f] + w[r267b] -
      w[r268f] + w[r268b] - w[r269f] + w[r269b] + w[r271f] - w[r271b] -
      w[r272f] + w[r272b] + w[r286f] - w[r286b] - w[r296f] + w[r296b] -
      w[r300f] + w[r300b] - w[r306f] + w[r306b] - w[r307f] + w[r307b] -
      w[r308f] + w[r308b] - w[r311f] + w[r311b] + w[r315f] - w[r315b] -
      w[r319f] + w[r319b] - w[r320f] + w[r320b] + w[r325f] - w[r325b] +
      w[r326f] - w[r326b] - w[r336f] + w[r336b] + w[r338f] - w[r338b] +
      w[r347f] - w[r347b] - w[r348f] + w[r348b] + w[r358f] - w[r358b] -
      w[r359f] + w[r359b] - w[r372f] + w[r372b] + w[r374f] - w[r374b] -
      w[r391f] + w[r391b] + w[r396] + w[r398f] - w[r398b] + w[r399f] -
      w[r399b] - w[r411f] + w[r411b] - w[r412f] + w[r412b] + w[r415f] -
      w[r415b] + w[r416f] - w[r416b];

  cdot[sCH3] = -w[r38f] + w[r38b] - w[r62f] + w[r62b] - w[r73] - w[r81f] +
               w[r81b] - w[r82f] + w[r82b] + w[r101f] - w[r101b] - w[r111f] +
               w[r111b] - w[r114f] + w[r114b] + w[r115f] - w[r115b] + w[r116f] -
               w[r116b] + w[r117f] - w[r117b] + w[r118f] - w[r118b] +
               0x1p+1 * w[r119f] - 0x1p+1 * w[r119b] - w[r120f] + w[r120b] -
               w[r121f] + w[r121b] - w[r122f] + w[r122b] - w[r123f] + w[r123b] -
               w[r124f] + w[r124b] - w[r131f] + w[r131b] - w[r132f] + w[r132b] -
               w[r133f] + w[r133b] - w[r134f] + w[r134b] - w[r135f] + w[r135b] -
               w[r136f] + w[r136b] + w[r138f] - w[r138b] - w[r140f] + w[r140b] +
               w[r154f] - w[r154b] + w[r161f] - w[r161b] - 0x1p+1 * w[r173f] +
               0x1p+1 * w[r173b] - w[r179f] + w[r179b] + w[r184f] - w[r184b] -
               w[r189f] + w[r189b] - 0x1p+1 * w[r190f] + 0x1p+1 * w[r190b] -
               w[r196f] + w[r196b] + w[r199f] - w[r199b] + w[r212f] - w[r212b] +
               w[r214f] - w[r214b] - w[r221f] + w[r221b] + w[r225f] - w[r225b] +
               w[r227f] - w[r227b] - w[r230f] + w[r230b] + w[r234f] - w[r234b] +
               w[r238f] - w[r238b] + w[r240f] - w[r240b] + w[r246f] - w[r246b] -
               w[r251f] + w[r251b] + w[r264f] - w[r264b] + w[r267f] - w[r267b] -
               w[r273f] + w[r273b] - w[r278f] + w[r278b] - w[r288f] + w[r288b] -
               w[r293f] + w[r293b] + w[r308f] - w[r308b] + w[r314f] - w[r314b] -
               w[r327f] + w[r327b] - w[r328f] + w[r328b] + w[r335f] - w[r335b] -
               w[r339f] + w[r339b] - w[r344f] + w[r344b] - w[r351f] + w[r351b] -
               w[r360f] + w[r360b] + w[r371f] - w[r371b] - w[r377f] + w[r377b] +
               w[r383f] - w[r383b] - w[r408f] + w[r408b] - w[r419f] + w[r419b] -
               w[r420f] + w[r420b] - w[r428f] + w[r428b];

  cdot[sCH2O] =
      -w[r39f] + w[r39b] - w[r41f] + w[r41b] + w[r44f] - w[r44b] + w[r57f] -
      w[r57b] + w[r58f] - w[r58b] - w[r59f] + w[r59b] - w[r60f] + w[r60b] -
      w[r61f] + w[r61b] - w[r62f] + w[r62b] - w[r63f] + w[r63b] - w[r64f] +
      w[r64b] - w[r77f] + w[r77b] + w[r78f] - w[r78b] + w[r79f] - w[r79b] -
      w[r80f] + w[r80b] + w[r82f] - w[r82b] + w[r83f] - w[r83b] + w[r84f] -
      w[r84b] - w[r85f] + w[r85b] + w[r86f] - w[r86b] + w[r87f] - w[r87b] +
      w[r88f] - w[r88b] + w[r89f] - w[r89b] + 0x1p+1 * w[r90f] -
      0x1p+1 * w[r90b] + w[r91f] - w[r91b] + w[r92f] - w[r92b] + w[r93f] -
      w[r93b] + w[r94f] - w[r94b] + w[r95f] - w[r95b] - w[r97f] + w[r97b] +
      w[r113f] - w[r113b] + w[r121f] - w[r121b] + w[r126f] - w[r126b] +
      w[r133f] - w[r133b] + w[r135f] - w[r135b] - w[r137f] + w[r137b] +
      w[r142] + w[r153f] - w[r153b] + w[r160f] - w[r160b] + w[r171f] -
      w[r171b] - w[r196f] + w[r196b] - w[r198f] + w[r198b] - w[r235f] +
      w[r235b] + w[r242] - w[r259f] + w[r259b] + w[r267f] - w[r267b] +
      w[r290f] - w[r290b] + w[r292] + w[r383f] - w[r383b] + w[r384f] -
      w[r384b] - w[r385f] + w[r385b] - w[r388f] + w[r388b] - w[r392f] +
      w[r392b] + 0x1p+1 * w[r396] - w[r400f] + w[r400b] + w[r425f] - w[r425b] +
      w[r426f] - w[r426b] - w[r430f] + w[r430b];

  cdot[sC2H6] = w[r173f] - w[r173b] + w[r174f] - w[r174b] - w[r175f] +
                w[r175b] - w[r176f] + w[r176b] - w[r177f] + w[r177b] -
                w[r178f] + w[r178b] - w[r179f] + w[r179b] - w[r180f] +
                w[r180b] - w[r181f] + w[r181b] - w[r182f] + w[r182b] -
                w[r183f] + w[r183b] - w[r184f] + w[r184b] - w[r202f] +
                w[r202b] - w[r236f] + w[r236b] + w[r329f] - w[r329b] +
                w[r364f] - w[r364b];

  cdot[sCH3OCH2] = w[r372f] - w[r372b] + w[r373f] - w[r373b] + w[r374f] -
                   w[r374b] + w[r375f] - w[r375b] + w[r376f] - w[r376b] +
                   w[r377f] - w[r377b] + w[r378f] - w[r378b] + w[r379f] -
                   w[r379b] + w[r380f] - w[r380b] + w[r381f] - w[r381b] +
                   w[r382f] - w[r382b] - w[r383f] + w[r383b] - w[r384f] +
                   w[r384b] - w[r385f] + w[r385b] - w[r386f] + w[r386b] -
                   w[r387f] + w[r387b];

  cdot[sO2CH2OCH2O2H] = w[r397f] - w[r397b] - w[r398f] + w[r398b];

  cdot[sCH2OCH2O2H] = w[r395f] - w[r395b] - w[r396] - w[r397f] + w[r397b];

  cdot[sHO2] = w[r9f] - w[r9b] + w[r10f] - w[r10b] + w[r11f] - w[r11b] -
               w[r12f] + w[r12b] + w[r13f] - w[r13b] - w[r14f] + w[r14b] -
               w[r15f] + w[r15b] - 0x1p+1 * w[r16f] + 0x1p+1 * w[r16b] -
               0x1p+1 * w[r17f] + 0x1p+1 * w[r17b] + w[r21f] - w[r21b] +
               w[r22f] - w[r22b] + w[r23f] - w[r23b] + w[r24f] - w[r24b] -
               w[r29f] + w[r29b] + w[r31f] - w[r31b] - w[r36] + w[r39f] -
               w[r39b] - w[r63f] + w[r63b] - w[r68f] + w[r68b] + w[r74f] -
               w[r74b] - w[r75] + w[r79f] - w[r79b] - w[r84f] + w[r84b] +
               w[r86f] - w[r86b] + w[r87f] - w[r87b] - w[r89f] + w[r89b] -
               w[r96f] + w[r96b] - w[r97f] + w[r97b] - w[r99f] + w[r99b] +
               w[r109f] - w[r109b] - w[r110f] + w[r110b] - w[r118f] + w[r118b] -
               w[r131f] + w[r131b] - w[r132f] + w[r132b] - w[r141f] + w[r141b] +
               w[r178f] - w[r178b] - w[r180f] + w[r180b] - w[r193f] + w[r193b] +
               w[r195f] - w[r195b] - w[r201f] + w[r201b] + w[r205f] - w[r205b] +
               w[r206f] - w[r206b] + w[r210f] - w[r210b] + w[r220f] - w[r220b] -
               w[r222f] + w[r222b] - w[r232f] + w[r232b] + w[r233f] - w[r233b] +
               w[r241f] - w[r241b] + w[r270f] - w[r270b] - w[r276f] + w[r276f] -
               w[r276b] + w[r276b] + w[r279f] - w[r279b] - w[r286f] + w[r286b] +
               w[r316f] - w[r316b] - w[r321f] + w[r321b] - w[r322f] + w[r322b] +
               w[r333f] - w[r333b] + w[r334f] - w[r334b] + w[r341f] - w[r341b] -
               w[r342f] + w[r342b] + w[r349f] - w[r349b] - w[r350f] + w[r350b] -
               w[r361f] + w[r361b] + w[r367f] - w[r367b] - w[r375f] + w[r375b] +
               w[r378f] - w[r378b] + w[r393f] - w[r393b] + w[r409f] - w[r409b] +
               w[r410f] - w[r410b] - w[r413f] + w[r413b] - w[r414f] + w[r414b];

  cdot[sH2O2] = w[r16f] - w[r16b] + w[r17f] - w[r17b] - w[r18f] + w[r18b] -
                w[r19f] + w[r19b] - w[r20f] + w[r20b] - w[r21f] + w[r21b] -
                w[r22f] + w[r22b] - w[r23f] + w[r23b] - w[r24f] + w[r24b] +
                w[r63f] - w[r63b] - w[r74f] + w[r74b] + w[r75] + w[r84f] -
                w[r84b] + w[r89f] - w[r89b] + w[r110f] - w[r110b] + w[r118f] -
                w[r118b] + w[r180f] - w[r180b] + w[r222f] - w[r222b] -
                w[r233f] + w[r233b] + w[r321f] - w[r321b] + w[r322f] -
                w[r322b] + w[r342f] - w[r342b] + w[r350f] - w[r350b] +
                w[r361f] - w[r361b] + w[r375f] - w[r375b] + w[r413f] -
                w[r413b] + w[r414f] - w[r414b];

  cdot[sCH3OCH2O2] = -w[r380f] + w[r380b] + w[r387f] - w[r387b] - w[r388f] +
                     w[r388b] - w[r389f] + w[r389b] - 0x1p+1 * w[r390] -
                     w[r395f] + w[r395b];

  cdot[sCH3O] =
      -w[r78f] + w[r78b] - w[r79f] + w[r79b] - w[r80f] + w[r80b] + w[r81f] -
      w[r81b] - w[r82f] + w[r82b] - w[r83f] + w[r83b] - w[r84f] + w[r84b] -
      w[r91f] + w[r91b] + w[r105f] - w[r105b] + w[r108f] - w[r108b] - w[r112f] +
      w[r112b] - 0x1p+1 * w[r113f] + 0x1p+1 * w[r113b] + w[r123f] - w[r123b] +
      w[r131f] - w[r131b] + w[r134f] - w[r134b] + 0x1p+1 * w[r140f] -
      0x1p+1 * w[r140b] + 0x1p+1 * w[r143] + w[r144f] - w[r144b] + w[r145f] -
      w[r145b] + w[r147f] - w[r147b] - w[r182f] + w[r182b] + w[r194f] -
      w[r194b] - w[r280f] + w[r280b] + w[r284f] - w[r284b] - w[r340f] +
      w[r340b] - w[r353f] + w[r353b] - w[r362f] + w[r362b] + w[r371f] -
      w[r371b] - w[r379f] + w[r379b] - w[r384f] + w[r384b] - w[r392f] +
      w[r392b] - w[r407f] + w[r407b] - w[r421f] + w[r421b] - w[r422f] +
      w[r422b] - w[r429f] + w[r429b];

  cdot[sCH3OH] =
      w[r80f] - w[r80b] - w[r81f] + w[r81b] + w[r91f] - w[r91b] - w[r92f] +
      w[r92b] + w[r95f] - w[r95b] - w[r101f] + w[r101b] - w[r102f] + w[r102b] -
      w[r103f] + w[r103b] - w[r104f] + w[r104b] - w[r105f] + w[r105b] -
      w[r106f] + w[r106b] - w[r107f] + w[r107b] - w[r108f] + w[r108b] -
      w[r109f] + w[r109b] - w[r110f] + w[r110b] - w[r111f] + w[r111b] -
      w[r112f] + w[r112f] - w[r112b] + w[r112b] + w[r113f] - w[r113b] -
      w[r139f] + w[r139b] + w[r142] + w[r146f] - w[r146b] + w[r182f] -
      w[r182b] - w[r200f] + w[r200b] + w[r280f] - w[r280b] + w[r340f] -
      w[r340b] + w[r353f] - w[r353b] + w[r362f] - w[r362b] + w[r379f] -
      w[r379b] + w[r406f] - w[r406b] + w[r421f] - w[r421b] + w[r422f] -
      w[r422b];

  cdot[sCH3OCH2O] = 0x1p+1 * w[r390] - w[r391f] + w[r391b] + w[r392f] -
                    w[r392b] - w[r393f] + w[r393b] + w[r394f] - w[r394b];

  cdot[sOCH2OCHO] =
      w[r399f] - w[r399b] + w[r400f] - w[r400b] - w[r401f] + w[r401b];

  cdot[sHO2CH2OCHO] = w[r398f] - w[r398b] - w[r399f] + w[r399b];

  cdot[sHCO] =
      -w[r30f] + w[r30b] - w[r31f] + w[r31b] - w[r32f] + w[r32b] - w[r33f] +
      w[r33b] - w[r34f] + w[r34b] - w[r35f] + w[r35b] - w[r36] -
      0x1p+1 * w[r37] - w[r38f] + w[r38b] + w[r39f] - w[r39b] - w[r40f] +
      w[r40b] + w[r41f] - w[r41b] - 0x1p+1 * w[r44f] + 0x1p+1 * w[r44b] -
      w[r57f] + w[r57b] + w[r59f] - w[r59b] + w[r60f] - w[r60b] + w[r61f] -
      w[r61b] + w[r62f] - w[r62b] + w[r63f] - w[r63b] + w[r77f] - w[r77b] +
      w[r80f] - w[r80b] - w[r90f] + w[r90b] - w[r92f] + w[r92b] + w[r125f] -
      w[r125b] + w[r137f] - w[r137b] + w[r152f] - w[r152b] + w[r162f] -
      w[r162b] + w[r168f] - w[r168b] + w[r170f] - w[r170b] + w[r172f] -
      w[r172b] + w[r198f] - w[r198b] + w[r212f] - w[r212b] + w[r214f] -
      w[r214b] + w[r235f] - w[r235b] + w[r264f] - w[r264b] + w[r290f] -
      w[r290b] + w[r301f] - w[r301b] - w[r309f] + w[r309b] + 0x1p+1 * w[r312f] -
      0x1p+1 * w[r312b] - w[r345f] + w[r345b] - w[r356f] + w[r356b] + w[r385f] -
      w[r385b] + w[r388f] - w[r388b] - w[r407f] + w[r407b] - w[r425f] +
      w[r425b] - w[r426f] + w[r426b] - w[r430f] + w[r430b];

  cdot[sO2CHO] = w[r40f] - w[r40b] - w[r41f] + w[r41b] - w[r381f] + w[r381b];

  cdot[sHOCHO] = w[r65f] - w[r65b] - w[r66f] + w[r66b] - w[r67f] + w[r67b] +
                 w[r68f] - w[r68b] - w[r69] - w[r70] - w[r71] - w[r72] -
                 w[r73] + w[r74f] - w[r74b] - w[r75] - w[r76] + w[r77f] -
                 w[r77b] + w[r225f] - w[r225b] + w[r382f] - w[r382b];

  cdot[sHOCH2O] = w[r64f] - w[r64b] - w[r65f] + w[r65b] + w[r96f] - w[r96b] -
                  w[r100f] + w[r100b] - w[r402f] + w[r402b];

  cdot[sHOCH2OCO] =
      w[r401f] - w[r401b] + w[r402f] - w[r402b] + w[r403f] - w[r403b];

  cdot[sCH2OH] = w[r85f] - w[r85b] - w[r86f] + w[r86b] - w[r87f] + w[r87b] -
                 w[r88f] + w[r88b] - w[r89f] + w[r89b] - w[r90f] + w[r90b] -
                 w[r91f] + w[r91b] + w[r92f] - w[r92b] - w[r93f] + w[r93b] -
                 w[r94f] + w[r94b] - 0x1p+1 * w[r95f] + 0x1p+1 * w[r95b] -
                 w[r96f] + w[r96b] + w[r103f] - w[r103b] + w[r104f] - w[r104b] +
                 w[r106f] - w[r106b] + w[r107f] - w[r107b] + w[r109f] -
                 w[r109b] + w[r110f] - w[r110b] + w[r111f] - w[r111b] +
                 w[r112f] - w[r112b] + w[r122f] - w[r122b] + w[r139f] -
                 w[r139b] + w[r200f] - w[r200b] + w[r250f] - w[r250b] +
                 w[r314f] - w[r314b] - w[r403f] + w[r403b];

  cdot[sCH3O2H] =
      w[r137f] - w[r137b] + w[r138f] - w[r138b] + w[r139f] - w[r139b] +
      w[r141f] - w[r141b] - w[r147f] + w[r147b] + w[r181f] - w[r181b] +
      w[r186f] - w[r186b] + w[r223f] - w[r223b] + w[r274f] - w[r274b] +
      w[r281f] - w[r281b] + w[r323f] - w[r323b] + w[r324f] - w[r324b] +
      w[r343f] - w[r343b] + w[r354f] - w[r354b] + w[r363f] - w[r363b] +
      w[r376f] - w[r376b] + w[r423f] - w[r423b] + w[r424f] - w[r424b];

  cdot[sCH3O2] =
      w[r136f] - w[r136b] - w[r137f] + w[r137b] - w[r138f] + w[r138b] -
      w[r139f] + w[r139b] - w[r140f] + w[r140b] - w[r141f] + w[r141b] -
      0x1p+1 * w[r142] - 0x1p+1 * w[r143] - w[r144f] + w[r144b] - w[r145f] +
      w[r145b] - w[r146f] + w[r146b] - w[r181f] + w[r181b] - w[r186f] +
      w[r186b] - w[r194f] + w[r194b] - w[r223f] + w[r223b] - w[r274f] +
      w[r274b] - w[r281f] + w[r281b] - w[r284f] + w[r284b] - w[r323f] +
      w[r323b] - w[r324f] + w[r324b] - w[r343f] + w[r343b] - w[r354f] +
      w[r354b] - w[r363f] + w[r363b] - w[r376f] + w[r376b] - w[r423f] +
      w[r423b] - w[r424f] + w[r424b];

  cdot[sHO2CHO] = w[r41f] - w[r41b] + w[r42f] - w[r42b] + w[r381f] - w[r381b];

  cdot[sCH3OCH2O2H] = w[r380f] - w[r380b] + w[r388f] - w[r388b] + w[r389f] -
                      w[r389b] + w[r391f] - w[r391b];

  cdot[sCH2] = -w[r119f] + w[r119b] + w[r148f] - w[r148b] + w[r149f] -
               w[r149b] + w[r157f] - w[r157b] + w[r158f] - w[r158b] + w[r159f] -
               w[r159b] - w[r161f] + w[r161b] - w[r162f] + w[r162b] - w[r163] -
               w[r164] - w[r165f] + w[r165b] - w[r166f] + w[r166b] - w[r167f] +
               w[r167b] + w[r183f] - w[r183b] - w[r243f] + w[r243b] + w[r247f] -
               w[r247b] + w[r304f] - w[r304b];

  cdot[sCH] = -w[r56f] + w[r56b] + w[r150f] - w[r150b] + w[r165f] - w[r165b] +
              w[r166f] - w[r166b] + w[r167f] - w[r167b] - w[r168f] + w[r168b] -
              w[r169f] + w[r169b] - w[r170f] + w[r170b] - w[r171f] + w[r171b] -
              w[r172f] + w[r172b] - w[r183f] + w[r183b] - w[r258f] + w[r258b] -
              w[r259f] + w[r259b] - w[r260f] + w[r260b] - w[r287f] + w[r287b] +
              w[r299f] - w[r299b];

  cdot[sCH3OCHO] = w[r393f] - w[r393b] - w[r394f] + w[r394b] + w[r404f] -
                   w[r404b] + w[r405f] - w[r405b] - w[r406f] + w[r406b] +
                   w[r407f] - w[r407b] + w[r408f] - w[r408b] - w[r409f] +
                   w[r409b] - w[r410f] + w[r410b] - w[r411f] + w[r411b] -
                   w[r412f] + w[r412b] - w[r413f] + w[r413b] - w[r414f] +
                   w[r414b] - w[r415f] + w[r415b] - w[r416f] + w[r416b] -
                   w[r417f] + w[r417b] - w[r418f] + w[r418b] - w[r419f] +
                   w[r419b] - w[r420f] + w[r420b] - w[r421f] + w[r421b] -
                   w[r422f] + w[r422b] - w[r423f] + w[r423b] - w[r424f] +
                   w[r424b] - w[r425f] + w[r425b] - w[r426f] + w[r426b];

  cdot[sCH2OCHO] = -w[r404f] + w[r404b] + w[r410f] - w[r410b] + w[r412f] -
                   w[r412b] + w[r414f] - w[r414b] + w[r416f] - w[r416b] +
                   w[r418f] - w[r418b] + w[r420f] - w[r420b] + w[r422f] -
                   w[r422b] + w[r424f] - w[r424b] + w[r426f] - w[r426b] +
                   w[r427f] - w[r427b] + w[r430f] - w[r430b];

  cdot[sC2H4O1X2] = w[r207f] - w[r207b] + w[r211f] - w[r211b] - w[r212f] +
                    w[r212b] - w[r213f] + w[r213b] + w[r284f] - w[r284b] +
                    w[r285f] - w[r285b] + w[r286f] - w[r286b];

  cdot[sC2H4] =
      -w[r185f] + w[r185b] - 0x1p+1 * w[r188f] + 0x1p+1 * w[r188b] + w[r189f] -
      w[r189b] + w[r191f] - w[r191b] + w[r205f] - w[r205b] + w[r206f] -
      w[r206b] + w[r210f] - w[r210b] + w[r252f] - w[r252b] + w[r261f] -
      w[r261b] - w[r262f] + w[r262b] - w[r263f] + w[r263b] - w[r264f] +
      w[r264b] - w[r265f] + w[r265b] - w[r266f] + w[r266b] - w[r267f] +
      w[r267b] - w[r268f] + w[r268b] - w[r269f] + w[r269b] - w[r278f] +
      w[r278b] - w[r279f] + w[r279b] - w[r280f] + w[r280b] - w[r281f] +
      w[r281b] - w[r282f] + w[r282b] - w[r283f] + w[r283b] - w[r284f] +
      w[r284b] - w[r285f] + w[r285b] - w[r286f] + w[r286b] + w[r287f] -
      w[r287b] + w[r288f] - w[r288b] + w[r297f] - w[r297b] + w[r313f] -
      w[r313b] + w[r352f] - w[r352b] + w[r369f] - w[r369b];

  cdot[sC2H3] = w[r188f] - w[r188b] - w[r261f] + w[r261b] + w[r263f] -
                w[r263b] + w[r266f] - w[r266b] + w[r278f] - w[r278b] +
                w[r279f] - w[r279b] + w[r280f] - w[r280b] + w[r281f] -
                w[r281b] + w[r282f] - w[r282b] + w[r283f] - w[r283b] +
                w[r289f] - w[r289b] - w[r290f] + w[r290b] - w[r291f] +
                w[r291b] - w[r292] - w[r293f] + w[r293b] - w[r294f] + w[r294b] -
                w[r295f] + w[r295b] - w[r296f] + w[r296b] - 0x1p+1 * w[r297f] +
                0x1p+1 * w[r297b] + w[r309f] - w[r309b] - w[r345f] + w[r345b] -
                w[r352f] + w[r352b] - w[r355f] + w[r355b] - w[r369f] + w[r369b];

  cdot[sC2H5] =
      -w[r174f] + w[r174b] + w[r175f] - w[r175b] + w[r176f] - w[r176b] +
      w[r177f] - w[r177b] + w[r178f] - w[r178b] + w[r179f] - w[r179b] +
      w[r180f] - w[r180b] + w[r181f] - w[r181b] + w[r182f] - w[r182b] +
      w[r183f] - w[r183b] + w[r184f] - w[r184b] + w[r185f] - w[r185b] +
      w[r188f] - w[r188b] - w[r189f] + w[r189b] + w[r190f] - w[r190b] -
      w[r191f] + w[r191b] - w[r192f] + w[r192b] - w[r193f] + w[r193b] -
      w[r194f] + w[r194b] + w[r202f] - w[r202b] - w[r204f] + w[r204b] -
      w[r205f] + w[r205b] - w[r206f] + w[r206b] - w[r207f] + w[r207b] -
      w[r208f] + w[r208b] + w[r236f] - w[r236b] + w[r251f] - w[r251b] +
      w[r315f] - w[r315b] - w[r329f] + w[r329b] - w[r356f] + w[r356b] -
      w[r364f] + w[r364b] - w[r370f] + w[r370b];

  cdot[sOHY] = -w[r45f] + w[r45b] - w[r46f] + w[r46b] - w[r47f] + w[r47b] -
               w[r48f] + w[r48b] - w[r49f] + w[r49b] - w[r50f] + w[r50b] -
               w[r51f] + w[r51b] - w[r52f] + w[r52b] - w[r53f] + w[r53b] -
               w[r54f] + w[r54b] - w[r55f] + w[r55b] + w[r56f] - w[r56b];

  cdot[sCH3OCO] = -w[r405f] + w[r405b] + w[r409f] - w[r409b] + w[r411f] -
                  w[r411b] + w[r413f] - w[r413b] + w[r415f] - w[r415b] +
                  w[r417f] - w[r417b] + w[r419f] - w[r419b] + w[r421f] -
                  w[r421b] + w[r423f] - w[r423b] + w[r425f] - w[r425b] -
                  w[r427f] + w[r427b] + w[r428f] - w[r428b] + w[r429f] -
                  w[r429b];

  cdot[sHOCH2O2H] = w[r99f] - w[r99b] + w[r100f] - w[r100b];

  cdot[sHOCH2O2] = w[r98f] - w[r98b] - w[r99f] + w[r99b];

  cdot[sC2H3CO] = w[r346f] - w[r346b] + w[r347f] - w[r347b] + w[r348f] -
                  w[r348b] + w[r349f] - w[r349b] + w[r350f] - w[r350b] +
                  w[r351f] - w[r351b] + w[r352f] - w[r352b] + w[r353f] -
                  w[r353b] + w[r354f] - w[r354b] + w[r355f] - w[r355b];

  cdot[sC2H3CHO] = w[r345f] - w[r345b] - w[r346f] + w[r346b] - w[r347f] +
                   w[r347b] - w[r348f] + w[r348b] - w[r349f] + w[r349b] -
                   w[r350f] + w[r350b] - w[r351f] + w[r351b] - w[r352f] +
                   w[r352b] - w[r353f] + w[r353b] - w[r354f] + w[r354b];

  cdot[sCH3CO] = w[r216f] - w[r216b] + w[r218f] - w[r218b] + w[r219f] -
                 w[r219b] + w[r220f] - w[r220b] + w[r221f] - w[r221b] +
                 w[r222f] - w[r222b] + w[r223f] - w[r223b] + w[r224f] -
                 w[r224b] - w[r227f] + w[r227b] - w[r228f] + w[r228b] -
                 w[r229f] + w[r229b] - w[r230f] + w[r230b] - w[r231f] +
                 w[r231b] - w[r244f] + w[r244b] + w[r335f] - w[r335b] +
                 w[r386f] - w[r386b] + w[r389f] - w[r389b];

  cdot[sC2H5O2H] = w[r187f] - w[r187b] + w[r198f] - w[r198b] + w[r199f] -
                   w[r199b] + w[r200f] - w[r200b] + w[r201f] - w[r201b] +
                   w[r202f] - w[r202b] - w[r203f] + w[r203b] + w[r282f] -
                   w[r282b] + w[r366f] - w[r366b];

  cdot[sC2H5O2] = -w[r187f] + w[r187b] - w[r198f] + w[r198b] - w[r199f] +
                  w[r199b] - w[r200f] + w[r200b] - w[r201f] + w[r201b] -
                  w[r202f] + w[r202b] + w[r204f] - w[r204b] - w[r209f] +
                  w[r209b] - w[r210f] + w[r210b] - w[r211f] + w[r211b] -
                  w[r282f] + w[r282b] - w[r285f] + w[r285b] - w[r366f] +
                  w[r366b];

  cdot[sCH2CO] = w[r228f] - w[r228b] + w[r229f] - w[r229b] + w[r230f] -
                 w[r230b] + w[r239f] - w[r239b] + w[r241f] - w[r241b] +
                 w[r243f] - w[r243b] + w[r244f] - w[r244b] - w[r245f] +
                 w[r245b] - w[r246f] + w[r246b] - w[r247f] + w[r247b] -
                 w[r248f] + w[r248b] - w[r249f] + w[r249b] - w[r250f] +
                 w[r250b] - w[r251f] + w[r251b] - w[r252f] + w[r252b] +
                 w[r259f] - w[r259b] + w[r307f] - w[r307b] + w[r311f] -
                 w[r311b] - w[r344f] + w[r344b];

  cdot[sCH2CHO] = w[r217f] - w[r217b] + w[r226f] - w[r226b] - w[r239f] +
                  w[r239b] - w[r240f] + w[r240b] - w[r241f] + w[r241b] -
                  w[r242] + w[r265f] - w[r265b] + w[r270f] - w[r270b] +
                  w[r271f] - w[r271b] + w[r272f] - w[r272b] + w[r273f] -
                  w[r273b] + w[r274f] - w[r274b] + w[r275f] - w[r275b] +
                  w[r291f] - w[r291b];

  cdot[sCH3CO3H] = w[r224f] - w[r224b] + w[r232f] - w[r232b] + w[r233f] -
                   w[r233b] + w[r234f] - w[r234b] + w[r235f] - w[r235b] +
                   w[r236f] - w[r236b] - w[r237f] + w[r237b] + w[r283f] -
                   w[r283b] + w[r368f] - w[r368b];

  cdot[sCH3CO3] = -w[r224f] + w[r224b] + w[r231f] - w[r231b] - w[r232f] +
                  w[r232b] - w[r233f] + w[r233b] - w[r234f] + w[r234b] -
                  w[r235f] + w[r235b] - w[r236f] + w[r236b] - w[r283f] +
                  w[r283b] - w[r368f] + w[r368b];

  cdot[sH2CC] = w[r262f] - w[r262b] + w[r295f] - w[r295b] + w[r303f] -
                w[r303b] - w[r310f] + w[r310b] - w[r311f] + w[r311b] -
                w[r312f] + w[r312b];

  cdot[sC2H2] = w[r260f] - w[r260b] - w[r289f] + w[r289b] + w[r293f] -
                w[r293b] + w[r294f] - w[r294b] + w[r296f] - w[r296b] +
                w[r297f] - w[r297b] + w[r298f] - w[r298b] + w[r302f] -
                w[r302b] - w[r303f] + w[r303b] - w[r304f] + w[r304b] -
                w[r305f] + w[r305b] - w[r306f] + w[r306b] - w[r307f] +
                w[r307b] - w[r308f] + w[r308b] - w[r309f] + w[r309b] +
                w[r310f] - w[r310b];

  cdot[sOCHO] = -w[r42f] + w[r42b] + w[r43f] - w[r43b] - w[r68f] + w[r68b] -
                w[r74f] + w[r74b] - w[r77f] + w[r77b] - w[r382f] + w[r382b] -
                w[r400f] + w[r400b] - w[r408f] + w[r408b];

  cdot[sSC2H4OH] = w[r316f] - w[r316b] + w[r317f] - w[r317b] + w[r319f] -
                   w[r319b] + w[r321f] - w[r321b] + w[r323f] - w[r323b] +
                   w[r325f] - w[r325b] + w[r327f] - w[r327b] + w[r329f] -
                   w[r329b] - w[r330f] + w[r330b] - w[r331f] + w[r331b] -
                   w[r332f] + w[r332b] - w[r333f] + w[r333b] - w[r334f] +
                   w[r334b];

  cdot[sC2H5OH] =
      -w[r313f] + w[r313b] - w[r314f] + w[r314b] - w[r315f] + w[r315b] -
      w[r316f] + w[r316b] - w[r317f] + w[r317b] - w[r318f] + w[r318b] -
      w[r319f] + w[r319b] - w[r320f] + w[r320b] - w[r321f] + w[r321b] -
      w[r322f] + w[r322b] - w[r323f] + w[r323b] - w[r324f] + w[r324b] -
      w[r325f] + w[r325b] - w[r326f] + w[r326b] - w[r327f] + w[r327b] -
      w[r328f] + w[r328b] - w[r329f] + w[r329b] + w[r365f] - w[r365b];

  cdot[sCH3CO2] = w[r237f] - w[r237b] - w[r238f] + w[r238b];

  cdot[sCH3CHO] =
      w[r192f] - w[r192b] + w[r195f] - w[r195b] - w[r197f] + w[r197b] +
      w[r208f] - w[r208b] + w[r209f] - w[r209b] + w[r213f] - w[r213b] -
      w[r214f] + w[r214b] - w[r215f] + w[r215b] - w[r216f] + w[r216b] -
      w[r217f] + w[r217b] - w[r218f] + w[r218b] - w[r219f] + w[r219b] -
      w[r220f] + w[r220b] - w[r221f] + w[r221b] - w[r222f] + w[r222b] -
      w[r223f] + w[r223b] - w[r224f] + w[r224b] - w[r225f] + w[r225b] -
      w[r226f] + w[r226b] + w[r268f] - w[r268b] + w[r276f] - w[r276b] +
      w[r277f] - w[r277b] + w[r330f] - w[r330b] + w[r333f] - w[r333b] -
      w[r386f] + w[r386b] - w[r389f] + w[r389b];

  cdot[sHCCO] = w[r245f] - w[r245b] + w[r248f] - w[r248b] + w[r249f] -
                w[r249b] - w[r253] - w[r254] - w[r255f] + w[r255b] - w[r256] -
                w[r257] + w[r258f] - w[r258b] - w[r260f] + w[r260b] + w[r300f] -
                w[r300b] + w[r305f] - w[r305b];

  cdot[sHCOH] = w[r124f] - w[r124b] - w[r125f] + w[r125b] - w[r126f] +
                w[r126b] - w[r127] - w[r128] - w[r129] - w[r130f] + w[r130b];

  cdot[sC2H3OH] = w[r269f] - w[r269b] - w[r270f] + w[r270b] - w[r271f] +
                  w[r271b] - w[r272f] + w[r272b] - w[r273f] + w[r273b] -
                  w[r274f] + w[r274b] - w[r275f] + w[r275b] - w[r276f] +
                  w[r276b] - w[r277f] + w[r277b] + w[r331f] - w[r331b] +
                  w[r334f] - w[r334b];

  cdot[sC2H5O] = w[r193f] - w[r193b] + w[r194f] - w[r194b] - w[r195f] +
                 w[r195b] + w[r196f] - w[r196b] + w[r197f] - w[r197b] +
                 w[r203f] - w[r203b] + w[r285f] - w[r285b] + w[r318f] -
                 w[r318b] + w[r320f] - w[r320b] + w[r322f] - w[r322b] +
                 w[r324f] - w[r324b] + w[r326f] - w[r326b] + w[r328f] -
                 w[r328b] + w[r332f] - w[r332b] - w[r365f] + w[r365b];

  cdot[sC2H5CO] = w[r357f] - w[r357b] + w[r358f] - w[r358b] + w[r359f] -
                  w[r359b] + w[r360f] - w[r360b] + w[r361f] - w[r361b] +
                  w[r362f] - w[r362b] + w[r363f] - w[r363b] + w[r364f] -
                  w[r364b] + w[r365f] - w[r365b] + w[r366f] - w[r366b] +
                  w[r367f] - w[r367b] + w[r368f] - w[r368b] + w[r369f] -
                  w[r369b] + w[r370f] - w[r370b];

  cdot[sC2H5CHO] = w[r356f] - w[r356b] - w[r357f] + w[r357b] - w[r358f] +
                   w[r358b] - w[r359f] + w[r359b] - w[r360f] + w[r360b] -
                   w[r361f] + w[r361b] - w[r362f] + w[r362b] - w[r363f] +
                   w[r363b] - w[r364f] + w[r364b] - w[r365f] + w[r365b] -
                   w[r366f] + w[r366b] - w[r367f] + w[r367b] - w[r368f] +
                   w[r368b] - w[r369f] + w[r369b];

  cdot[sC2H] = -w[r298f] + w[r298b] - w[r299f] + w[r299b] - w[r300f] +
               w[r300b] - w[r301f] + w[r301b] - w[r302f] + w[r302b] + w[r306f] -
               w[r306b];

  cdot[sCH3COCH3] = -w[r335f] + w[r335b] - w[r336f] + w[r336b] - w[r337f] +
                    w[r337b] - w[r338f] + w[r338b] - w[r339f] + w[r339b] -
                    w[r340f] + w[r340b] - w[r341f] + w[r341b] - w[r342f] +
                    w[r342b] - w[r343f] + w[r343b];

  cdot[sCH3COCH2] = w[r336f] - w[r336b] + w[r337f] - w[r337b] + w[r338f] -
                    w[r338b] + w[r339f] - w[r339b] + w[r340f] - w[r340b] +
                    w[r341f] - w[r341b] + w[r342f] - w[r342b] + w[r343f] -
                    w[r343b] + w[r344f] - w[r344b];

  cdot[sCH2Y] = w[r102f] - w[r102b] + w[r120f] - w[r120b] - w[r148f] +
                w[r148b] - w[r149f] + w[r149b] - w[r150f] + w[r150b] -
                w[r151f] + w[r151b] - w[r152f] + w[r152b] - w[r153f] +
                w[r153b] - w[r154f] + w[r154b] - w[r155] - w[r156f] + w[r156b] -
                w[r157f] + w[r157b] - w[r158f] + w[r158b] - w[r159f] +
                w[r159b] - w[r160f] + w[r160b] - w[r184f] + w[r184b] -
                w[r252f] + w[r252b] + w[r255f] - w[r255b] - w[r288f] + w[r288b];

  cdot[sOCH2O2H] = w[r97f] - w[r97b] - w[r98f] + w[r98b];

  cdot[sHE] = 0.0;

  cdot[sAR] = -w[r50f] + w[r50f] - w[r50b] + w[r50b] - w[r149f] + w[r149f] -
              w[r149b] + w[r149b];
}

static double GetLindRateCoeff(double temp, double pressure, double k0,
                               double kInf, double fc, double conc) {
  const double R = 8314.34; /* [J / kmole K] */
  double Ntmp;
  double kl;
  double f;
  double cCoeff, DCoeff, log10kNull;

  int iTroe = 1;

  if (isinf(conc)) {
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
void AramcoMech_DMEonly_74spec::getMolarMass(span<double> W) const {
  W[sN2] = 0x1.c051eb851eb85p+4;
  W[sCH3OCH3] = 0x1.708b439581062p+5;
  W[sCH4] = 0x1.00ac083126e98p+4;
  W[sH2] = 0x1.020c49ba5e354p+1;
  W[sO2] = 0x1p+5;
  W[sCO] = 0x1.c028f5c28f5c2p+4;
  W[sCO2] = 0x1.60147ae147ae1p+5;
  W[sH2O] = 0x1.204189374bc6ap+4;
  W[sH] = 0x1.020c49ba5e354p+0;
  W[sO] = 0x1p+4;
  W[sOH] = 0x1.1020c49ba5e35p+4;
  W[sCH3] = 0x1.e116872b020c4p+3;
  W[sCH2O] = 0x1.e06a7ef9db22cp+4;
  W[sC2H6] = 0x1.e116872b020c4p+4;
  W[sCH3OCH2] = 0x1.687ae147ae148p+5;
  W[sO2CH2OCH2O2H] = 0x1.b43d70a3d70a4p+6;
  W[sCH2OCH2O2H] = 0x1.343d70a3d70a4p+6;
  W[sHO2] = 0x1.0810624dd2f1bp+5;
  W[sH2O2] = 0x1.1020c49ba5e35p+5;
  W[sCH3OCH2O2] = 0x1.343d70a3d70a4p+6;
  W[sCH3O] = 0x1.f08b439581062p+4;
  W[sCH3OH] = 0x1.005604189374cp+5;
  W[sCH3OCH2O] = 0x1.e87ae147ae148p+5;
  W[sOCH2OCHO] = 0x1.2c2d0e5604189p+6;
  W[sHO2CH2OCHO] = 0x1.70353f7ced916p+6;
  W[sHCO] = 0x1.d049ba5e353f8p+4;
  W[sO2CHO] = 0x1.e824dd2f1a9fcp+5;
  W[sHOCHO] = 0x1.70353f7ced916p+5;
  W[sHOCH2O] = 0x1.7845a1cac0831p+5;
  W[sHOCH2OCO] = 0x1.2c2d0e5604189p+6;
  W[sCH2OH] = 0x1.f08b439581062p+4;
  W[sCH3O2H] = 0x1.805604189374bp+5;
  W[sCH3O2] = 0x1.7845a1cac0831p+5;
  W[sHO2CHO] = 0x1.f0353f7ced916p+5;
  W[sCH3OCH2O2H] = 0x1.3845a1cac0831p+6;
  W[sCH2] = 0x1.c0d4fdf3b645ap+3;
  W[sCH] = 0x1.a09374bc6a7fp+3;
  W[sCH3OCHO] = 0x1.e06a7ef9db22cp+5;
  W[sCH2OCHO] = 0x1.d85a1cac08312p+5;
  W[sC2H4O1X2] = 0x1.606a7ef9db22dp+5;
  W[sC2H4] = 0x1.c0d4fdf3b645ap+4;
  W[sC2H3] = 0x1.b0b4395810625p+4;
  W[sC2H5] = 0x1.d0f5c28f5c28fp+4;
  W[sOHY] = 0x1.1020c49ba5e35p+4;
  W[sCH3OCO] = 0x1.d85a1cac08312p+5;
  W[sHOCH2O2H] = 0x1.002b020c49ba6p+6;
  W[sHOCH2O2] = 0x1.f845a1cac0831p+5;
  W[sC2H3CO] = 0x1.b86e978d4fdf4p+5;
  W[sC2H3CHO] = 0x1.c07ef9db22d0ep+5;
  W[sCH3CO] = 0x1.585a1cac08312p+5;
  W[sC2H5O2H] = 0x1.f08b439581062p+5;
  W[sC2H5O2] = 0x1.e87ae147ae148p+5;
  W[sCH2CO] = 0x1.5049ba5e353f8p+5;
  W[sCH2CHO] = 0x1.585a1cac08312p+5;
  W[sCH3CO3H] = 0x1.30353f7ced916p+6;
  W[sCH3CO3] = 0x1.2c2d0e5604189p+6;
  W[sH2CC] = 0x1.a09374bc6a7fp+4;
  W[sC2H2] = 0x1.a09374bc6a7fp+4;
  W[sOCHO] = 0x1.6824dd2f1a9fcp+5;
  W[sSC2H4OH] = 0x1.687ae147ae148p+5;
  W[sC2H5OH] = 0x1.708b439581062p+5;
  W[sCH3CO2] = 0x1.d85a1cac08312p+5;
  W[sCH3CHO] = 0x1.606a7ef9db22dp+5;
  W[sHCCO] = 0x1.48395810624ddp+5;
  W[sHCOH] = 0x1.e06a7ef9db22cp+4;
  W[sC2H3OH] = 0x1.606a7ef9db22dp+5;
  W[sC2H5O] = 0x1.687ae147ae148p+5;
  W[sC2H5CO] = 0x1.c88f5c28f5c29p+5;
  W[sC2H5CHO] = 0x1.d09fbe76c8b44p+5;
  W[sC2H] = 0x1.9072b020c49bap+4;
  W[sCH3COCH3] = 0x1.d09fbe76c8b44p+5;
  W[sCH3COCH2] = 0x1.c88f5c28f5c29p+5;
  W[sCH2Y] = 0x1.c0d4fdf3b645ap+3;
  W[sOCH2O2H] = 0x1.f845a1cac0831p+5;
  W[sHE] = 0x1.42e00d1b71759p+4;
  W[sAR] = 0x1.3f95810624dd3p+5;
}

std::vector<std::string> AramcoMech_DMEonly_74spec::getSpeciesNames() const {
  return {"N2",         "CH3OCH3",    "CH4",      "H2",       "O2",
          "CO",         "CO2",        "H2O",      "H",        "O",
          "OH",         "CH3",        "CH2O",     "C2H6",     "CH3OCH2",
          "O2CH2OCH2O", "CH2OCH2O2H", "HO2",      "H2O2",     "CH3OCH2O2",
          "CH3O",       "CH3OH",      "CH3OCH2O", "OCH2OCHO", "HO2CH2OCHO",
          "HCO",        "O2CHO",      "HOCHO",    "HOCH2O",   "HOCH2OCO",
          "CH2OH",      "CH3O2H",     "CH3O2",    "HO2CHO",   "CH3OCH2O2H",
          "CH2",        "CH",         "CH3OCHO",  "CH2OCHO",  "C2H4O1-2",
          "C2H4",       "C2H3",       "C2H5",     "OH*",      "CH3OCO",
          "HOCH2O2H",   "HOCH2O2",    "C2H3CO",   "C2H3CHO",  "CH3CO",
          "C2H5O2H",    "C2H5O2",     "CH2CO",    "CH2CHO",   "CH3CO3H",
          "CH3CO3",     "H2CC",       "C2H2",     "OCHO",     "SC2H4OH",
          "C2H5OH",     "CH3CO2",     "CH3CHO",   "HCCO",     "HCOH",
          "C2H3OH",     "C2H5O",      "C2H5CO",   "C2H5CHO",  "C2H",
          "CH3COCH3",   "CH3COCH2",   "CH2*",     "OCH2O2H",  "HE",
          "AR"};
}

void AramcoMech_DMEonly_74spec::ComputeThermoData(span<double> h,
                                                  span<double> cp, double T,
                                                  span<double> s) const {
  /*
          This function computes entropy 's', enthalpy 'h' and heat
          capacity 'cp' as function of temperature 'T' for all non steady
          state species in units [J/(kg K)], [J/kg] and [J/(kg K)],
     respectively. The parameter s, h and cp should provide workspace of length
     76 If you compile using cdecl, 's' is optional for backwards compatibility:
          This function will only store entropies if *s == 0, to detect omitted
          arguments (in most cases). On the other hand, if *s == 0 after
     invocation,
          you were using an old mechanism without support for entropies. */

  int i;
  if (s.empty() || (s[0] != 0))
    s = span<double>((double*)alloca(sizeof(double) * 76), 76);
  if (T > 1000.0) {
    h[sN2] = 0x1.28bba88de9adap+8 *
             (T * (0x1.79ee05c20bd3ap+1 +
                   T * (0x1.6e306622946cap-11 +
                        T * (-0x1.60a3b229df0c4p-23 +
                             T * (0x1.59b0e16a98bbdp-36 +
                                  T * -0x1.099b777519683p-50)))) -
              0x1.cdf96e9bbf0dcp+9);
    cp[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.79ee05c20bd3ap+1 +
         T * (0x1.6e306622946cap-10 +
              T * (-0x1.087ac59f67493p-21 +
                   T * (0x1.59b0e16a98bbdp-34 + T * -0x1.4c0255525fc23p-48))));
    s[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.79ee05c20bd3ap+1 * log(T) +
         T * (0x1.6e306622946cap-10 +
              T * (-0x1.087ac59f67493p-22 +
                   T * (0x1.ccebd738cba51p-36 + T * -0x1.4c0255525fc23p-50))) +
         0x1.77cd01bb6bfc5p+2);
    h[sCH3OCH3] = 0x1.68f6f3733777dp+7 *
                  (T * (0x1.d745018713b01p+2 +
                        T * (0x1.c72e4835fb6eep-8 +
                             T * (-0x1.a833d83139d93p-20 +
                                  T * (0x1.9400641f736d9p-33 +
                                       T * -0x1.32deb3c3cf397p-47)))) -
                   0x1.980b3ac710cb3p+14);
    cp[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.d745018713b01p+2 +
         T * (0x1.c72e4835fb6eep-7 +
              T * (-0x1.3e26e224eb62ep-18 +
                   T * (0x1.9400641f736d9p-31 + T * -0x1.7f9660b4c307dp-45))));
    s[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.d745018713b01p+2 * log(T) +
         T * (0x1.c72e4835fb6eep-7 +
              T * (-0x1.3e26e224eb62ep-19 +
                   T * (0x1.0d559814f79e6p-32 + T * -0x1.7f9660b4c307dp-47))) -
         0x1.0221f22daf65bp+4);
    h[sCH4] = 0x1.0325882935d13p+9 *
              (T * (0x1.a73c320a693a8p+0 +
                    T * (0x1.488ac88f5e28fp-8 +
                         T * (-0x1.28c3f4ace66e9p-20 +
                              T * (0x1.26ef4a1545f72p-33 +
                                   T * -0x1.c5869c71e7825p-48)))) -
               0x1.38ccbfb15b574p+13);
    cp[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.a73c320a693a8p+0 +
         T * (0x1.488ac88f5e28fp-7 +
              T * (-0x1.bd25ef0359a5ep-19 +
                   T * (0x1.26ef4a1545f72p-31 + T * -0x1.1b7421c730b17p-45))));
    s[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.a73c320a693a8p+0 * log(T) +
         T * (0x1.488ac88f5e28fp-7 +
              T * (-0x1.bd25ef0359a5ep-20 +
                   T * (0x1.893f0d71b29edp-33 + T * -0x1.1b7421c730b17p-47))) +
         0x1.3cf64652f59b5p+3);
    h[sH2] = 0x1.01c3c6b46b46cp+12 *
             (T * (0x1.77682517e77d5p+1 +
                   T * (0x1.b16173819987dp-12 +
                        T * (-0x1.a3321154deb67p-25 +
                             T * (0x1.0f18a81defe5ep-38 +
                                  T * -0x1.3da7b7af53224p-53)))) -
              0x1.968864f54d1e9p+9);
    cp[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.77682517e77d5p+1 +
         T * (0x1.b16173819987dp-11 +
              T * (-0x1.3a658cffa708dp-23 +
                   T * (0x1.0f18a81defe5ep-36 + T * -0x1.8d11a59b27eacp-51))));
    s[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.77682517e77d5p+1 * log(T) +
         T * (0x1.b16173819987dp-11 +
              T * (-0x1.3a658cffa708dp-24 +
                   T * (0x1.6976357d3fdd3p-38 + T * -0x1.8d11a59b27eacp-53))) -
         0x1.063a67041b17bp+0);
    h[sO2] = 0x1.03d3adab9f55ap+8 *
             (T * (0x1.d49a5bcbd8b1p+1 +
                   T * (0x1.581fed8a272abp-12 +
                        T * (-0x1.9427c550a5b58p-25 +
                             T * (0x1.6a0b267c8c751p-38 +
                                  T * -0x1.2b8f5b630c4f6p-52)))) -
              0x1.2ffe8a1dfb939p+10);
    cp[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.d49a5bcbd8b1p+1 +
         T * (0x1.581fed8a272abp-11 +
              T * (-0x1.2f1dd3fc7c482p-23 +
                   T * (0x1.6a0b267c8c751p-36 + T * -0x1.7673323bcf633p-50))));
    s[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.d49a5bcbd8b1p+1 * log(T) +
         T * (0x1.581fed8a272abp-11 +
              T * (-0x1.2f1dd3fc7c482p-24 +
                   T * (0x1.e2b988a6109c1p-38 + T * -0x1.7673323bcf633p-52))) +
         0x1.b52a9b9f833d9p+1);
    h[sCO] = 0x1.28d6c752d48edp+8 *
             (T * (0x1.8634c9356897ap+1 +
                   T * (0x1.6258efee2a8eap-11 +
                        T * (-0x1.5bbeb37c2824cp-23 +
                             T * (0x1.5acd4f7896f0cp-36 +
                                  T * -0x1.0ed3580521edfp-50)))) -
              0x1.bdd0ef9db22d1p+13);
    cp[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.8634c9356897ap+1 +
         T * (0x1.6258efee2a8eap-10 +
              T * (-0x1.04cf069d1e1b9p-21 +
                   T * (0x1.5acd4f7896f0cp-34 + T * -0x1.52882e066a697p-48))));
    s[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.8634c9356897ap+1 * log(T) +
         T * (0x1.6258efee2a8eap-10 +
              T * (-0x1.04cf069d1e1b9p-22 +
                   T * (0x1.ce6714a0c941p-36 + T * -0x1.52882e066a697p-50))) +
         0x1.811820f3958e7p+2);
    h[sCO2] = 0x1.79d818115de7cp+7 *
              (T * (0x1.28bc990d829f8p+2 +
                    T * (0x1.675407753188dp-10 +
                         T * (-0x1.647214e9a182bp-22 +
                              T * (0x1.60b1a9abb6d1cp-35 +
                                   T * -0x1.081385acab026p-49)))) -
               0x1.7f01ced916873p+15);
    cp[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.28bc990d829f8p+2 +
         T * (0x1.675407753188dp-9 +
              T * (-0x1.0b558faf3922p-20 +
                   T * (0x1.60b1a9abb6d1cp-33 + T * -0x1.4a186717d5c2fp-47))));
    s[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.28bc990d829f8p+2 * log(T) +
         T * (0x1.675407753188dp-9 +
              T * (-0x1.0b558faf3922p-21 +
                   T * (0x1.d642378f9e6dp-35 + T * -0x1.4a186717d5c2fp-49))) -
         0x1.ef554fbdad752p+0);
    h[sH2O] = 0x1.cd8113bbf9d1p+8 *
              (T * (0x1.56a935eecf561p+1 +
                    T * (0x1.85b36b79611a1p-10 +
                         T * (-0x1.14f153063bb17p-22 +
                              T * (0x1.9f52af42f4ce4p-36 +
                                   T * -0x1.ec2e9d39714d2p-51)))) -
               0x1.d2f79374bc6a8p+14);
    cp[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.56a935eecf561p+1 +
         T * (0x1.85b36b79611a1p-9 +
              T * (-0x1.9f69fc89598a3p-21 +
                   T * (0x1.9f52af42f4ce4p-34 + T * -0x1.339d2243e6d03p-48))));
    s[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.56a935eecf561p+1 * log(T) +
         T * (0x1.85b36b79611a1p-9 +
              T * (-0x1.9f69fc89598a3p-22 +
                   T * (0x1.14e1ca2ca3343p-35 + T * -0x1.339d2243e6d03p-50))) +
         0x1.b87bb2fec56d6p+2);
    h[sH] = 0x1.01c3c6b46b46cp+13 *
            (T * (0x1.4p+1 +
                  T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) +
             0x1.8e06a3d70a3d7p+14);
    cp[sH] =
        0x1.01c3c6b46b46cp+13 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sH] = 0x1.01c3c6b46b46cp+13 *
            (0x1.4p+1 * log(T) +
             T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) -
             0x1.c9673aa1bc7ddp-2);
    h[sO] = 0x1.03d3adab9f55ap+9 *
            (T * (0x1.4595e56f85f4dp+1 +
                  T * (-0x1.ca4a64f681eadp-17 +
                       T * (-0x1.7ff09a88f9768p-30 +
                            T * (0x1.5caa12f83efbfp-40 +
                                 T * -0x1.ba4f6a20df63cp-54)))) +
             0x1.c8a80c49ba5e3p+14);
    cp[sO] =
        0x1.03d3adab9f55ap+9 *
        (0x1.4595e56f85f4dp+1 +
         T * (-0x1.ca4a64f681eadp-16 +
              T * (-0x1.1ff473e6bb18ep-28 +
                   T * (0x1.5caa12f83efbfp-38 + T * -0x1.1471a2548b9e5p-51))));
    s[sO] =
        0x1.03d3adab9f55ap+9 *
        (0x1.4595e56f85f4dp+1 * log(T) +
         T * (-0x1.ca4a64f681eadp-16 +
              T * (-0x1.1ff473e6bb18ep-29 +
                   T * (0x1.d0e2c3f5a94ffp-40 + T * -0x1.1471a2548b9e5p-53))) +
         0x1.3b06dfcddb6aap+2);
    h[sOH] = 0x1.e8db171241376p+8 *
             (T * (0x1.6b54f63c06ec6p+1 +
                   T * (0x1.224d389537523p-11 +
                        T * (-0x1.a4e83753bab4dp-24 +
                             T * (0x1.720cee94cc5e3p-37 +
                                  T * -0x1.17575a7da433bp-51)))) +
              0x1.ce39dbca9691ap+11);
    cp[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.6b54f63c06ec6p+1 +
         T * (0x1.224d389537523p-10 +
              T * (-0x1.3bae297ecc07ap-22 +
                   T * (0x1.720cee94cc5e3p-35 + T * -0x1.5d2d311d0d40ap-49))));
    s[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.6b54f63c06ec6p+1 * log(T) +
         T * (0x1.224d389537523p-10 +
              T * (-0x1.3bae297ecc07ap-23 +
                   T * (0x1.ed669371107d9p-37 + T * -0x1.5d2d311d0d40ap-51))) +
         0x1.76139a9191377p+2);
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
    s[sCH3] =
        0x1.148599bc79368p+9 *
        (0x1.7d330e4a459e7p+1 * log(T) +
         T * (0x1.7bf7d0ba1fd49p-8 +
              T * (-0x1.09286973eb98fp-20 +
                   T * (0x1.c280e46f95a33p-34 + T * -0x1.42c57f31cb48p-48))) +
         0x1.2e3d1c55a11c4p+2);
    h[sCH2O] = 0x1.14e89ec288e65p+8 *
               (T * (0x1.95b30c9cc674cp+1 +
                     T * (0x1.95e0bf82684dep-9 +
                          T * (-0x1.92c109bdb320dp-21 +
                               T * (0x1.9264fd893a18ep-34 +
                                    T * -0x1.3d44d06f01e04p-48)))) -
                0x1.c6a576fd21ff3p+13);
    cp[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.95b30c9cc674cp+1 +
         T * (0x1.95e0bf82684dep-8 +
              T * (-0x1.2e10c74e4658ap-19 +
                   T * (0x1.9264fd893a18ep-32 + T * -0x1.8c96048ac2584p-46))));
    s[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.95b30c9cc674cp+1 * log(T) +
         T * (0x1.95e0bf82684dep-8 +
              T * (-0x1.2e10c74e4658ap-20 +
                   T * (0x1.0c4353b0d165fp-33 + T * -0x1.8c96048ac2584p-48))) +
         0x1.82b16c08bcbd1p+2);
    h[sC2H6] = 0x1.148599bc79368p+8 *
               (T * (0x1.02fc8b7696346p+2 +
                     T * (0x1.f71daea9b099ep-8 +
                          T * (-0x1.e97b8ea557d95p-20 +
                               T * (0x1.e2971a5f2372bp-33 +
                                    T * -0x1.78fb602219704p-47)))) -
                0x1.84facc985f06fp+13);
    cp[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.02fc8b7696346p+2 +
         T * (0x1.f71daea9b099ep-7 +
              T * (-0x1.6f1caafc01e3p-18 +
                   T * (0x1.e2971a5f2372bp-31 + T * -0x1.d73a382a9fcc5p-45))));
    s[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.02fc8b7696346p+2 * log(T) +
         T * (0x1.f71daea9b099ep-7 +
              T * (-0x1.6f1caafc01e3p-19 +
                   T * (0x1.41ba1194c24c7p-32 + T * -0x1.d73a382a9fcc5p-47))) -
         0x1.eff939ac0cfd8p-1);
    h[sCH3OCH2] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.aece6263a266cp+2 +
                        T * (0x1.82bc7953ba229p-8 +
                             T * (-0x1.6601e4223d2b7p-20 +
                                  T * (0x1.53530601c9a51p-33 +
                                       T * -0x1.00d6a80f32437p-47)))) -
                   0x1.a2aab11c6d1e1p+11);
    cp[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.aece6263a266cp+2 +
         T * (0x1.82bc7953ba229p-7 +
              T * (-0x1.0c816b19ade09p-18 +
                   T * (0x1.53530601c9a51p-31 + T * -0x1.410c5212fed45p-45))));
    s[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.aece6263a266cp+2 * log(T) +
         T * (0x1.82bc7953ba229p-7 +
              T * (-0x1.0c816b19ade09p-19 +
                   T * (0x1.c46eb2ad0cdc1p-33 + T * -0x1.410c5212fed45p-47))) -
         0x1.2d0ce7497ac85p+3);
    h[sO2CH2OCH2O2H] = 0x1.30f32e834a014p+6 *
                       (T * (0x1.39eeb702602c9p+4 +
                             T * (0x1.46699569b9fb1p-8 +
                                  T * (-0x1.318d32d28c33ep-20 +
                                       T * (0x1.241ef4654de51p-33 +
                                            T * -0x1.bd2cc27876aecp-48)))) -
                        0x1.3777ed9e83e42p+15);
    cp[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (0x1.39eeb702602c9p+4 +
         T * (0x1.46699569b9fb1p-7 +
              T * (-0x1.ca53cc3bd24ddp-19 +
                   T * (0x1.241ef4654de51p-31 + T * -0x1.163bf98b4a2d3p-45))));
    s[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (0x1.39eeb702602c9p+4 * log(T) +
         T * (0x1.46699569b9fb1p-7 +
              T * (-0x1.ca53cc3bd24ddp-20 +
                   T * (0x1.857e9b31bd317p-33 + T * -0x1.163bf98b4a2d3p-47))) -
         0x1.0d16abffcdab2p+6);
    h[sCH2OCH2O2H] = 0x1.af956cbfc9bdbp+6 *
                     (T * (0x1.e34cc93c1e944p+3 +
                           T * (0x1.2ce582ba7dbf7p-8 +
                                T * (-0x1.1aa88bdef57a1p-20 +
                                     T * (0x1.0ee9f1a979bcbp-33 +
                                          T * -0x1.9d9cca344d604p-48)))) -
                      0x1.3334bdd97f62bp+14);
    cp[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.e34cc93c1e944p+3 +
         T * (0x1.2ce582ba7dbf7p-7 +
              T * (-0x1.a7fcd1ce70372p-19 +
                   T * (0x1.0ee9f1a979bcbp-31 + T * -0x1.0281fe60b05c2p-45))));
    s[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.e34cc93c1e944p+3 * log(T) +
         T * (0x1.2ce582ba7dbf7p-7 +
              T * (-0x1.a7fcd1ce70372p-20 +
                   T * (0x1.6937ece1f7a64p-33 + T * -0x1.0281fe60b05c2p-47))) -
         0x1.84235a2793c2cp+5);
    h[sHO2] = 0x1.f7c8d6a0c5c07p+7 *
              (T * (0x1.0b06c1c5dd9a5p+2 +
                    T * (0x1.ed239a3b0b8e1p-11 +
                         T * (-0x1.efbff7abb6e8cp-24 +
                              T * (0x1.5671f76d783fcp-38 +
                                   T * 0x1.4522f855802p-55)))) +
               0x1.f054b8a420dc2p+4);
    cp[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.0b06c1c5dd9a5p+2 +
         T * (0x1.ed239a3b0b8e1p-10 +
              T * (-0x1.73cff9c0c92e9p-22 +
                   T * (0x1.5671f76d783fcp-36 + T * 0x1.966bb66ae028p-53))));
    s[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.0b06c1c5dd9a5p+2 * log(T) +
         T * (0x1.ed239a3b0b8e1p-10 +
              T * (-0x1.73cff9c0c92e9p-23 +
                   T * (0x1.c897f491f5aa5p-38 + T * 0x1.966bb66ae028p-55))) +
         0x1.7a9526984530bp+1);
    h[sH2O2] = 0x1.e8db171241376p+7 *
               (T * (0x1.251b006c368ap+2 +
                     T * (0x1.09a26b45627bap-9 +
                          T * (-0x1.d0bb7df5eb98fp-22 +
                               T * (0x1.b3df193095563p-35 +
                                    T * -0x1.487e202329df9p-49)))) -
                0x1.195cb5c28f5c3p+14);
    cp[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.251b006c368ap+2 +
         T * (0x1.09a26b45627bap-8 +
              T * (-0x1.5c8c9e7870b2bp-20 +
                   T * (0x1.b3df193095563p-33 + T * -0x1.9a9da82bf4577p-47))));
    s[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.251b006c368ap+2 * log(T) +
         T * (0x1.09a26b45627bap-8 +
              T * (-0x1.5c8c9e7870b2bp-21 +
                   T * (0x1.2294bb75b8e42p-34 + T * -0x1.9a9da82bf4577p-49))) +
         0x1.547709ef0e8d8p-1);
    h[sCH3OCH2O2] = 0x1.af956cbfc9bdbp+6 *
                    (T * (0x1.867a2ecff723cp+3 +
                          T * (0x1.86ad44208a46ep-8 +
                               T * (-0x1.6bd775567cf53p-20 +
                                    T * (0x1.5a63d126458fp-33 +
                                         T * -0x1.070ea9ef39863p-47)))) -
                     0x1.7782deb851eb8p+14);
    cp[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.867a2ecff723cp+3 +
         T * (0x1.86ad44208a46ep-7 +
              T * (-0x1.10e19800ddb7ep-18 +
                   T * (0x1.5a63d126458fp-31 + T * -0x1.48d2546b07e7bp-45))));
    s[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.867a2ecff723cp+3 * log(T) +
         T * (0x1.86ad44208a46ep-7 +
              T * (-0x1.10e19800ddb7ep-19 +
                   T * (0x1.cdda6c3307695p-33 + T * -0x1.48d2546b07e7bp-47))) -
         0x1.0eeeafc62bc8ep+5);
    h[sCH3O] = 0x1.0bea1f388b9cfp+8 *
               (T * (0x1.307fab9c50832p+2 +
                     T * (0x1.e7ae63e4d9997p-9 +
                          T * (-0x1.e2a7fd8e3929cp-21 +
                               T * (0x1.e1af83b06912p-34 +
                                    T * -0x1.7bcc072bf8603p-48)))) +
                0x1.7a1ca8198f1d4p+8);
    cp[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.307fab9c50832p+2 +
         T * (0x1.e7ae63e4d9997p-8 +
              T * (-0x1.69fdfe2aaadf5p-19 +
                   T * (0x1.e1af83b06912p-32 + T * -0x1.dabf08f6f6783p-46))));
    s[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.307fab9c50832p+2 * log(T) +
         T * (0x1.e7ae63e4d9997p-8 +
              T * (-0x1.69fdfe2aaadf5p-20 +
                   T * (0x1.411fad20460cp-33 + T * -0x1.dabf08f6f6783p-48))) -
         0x1.f78039205d17bp+0);
    h[sCH3OH] = 0x1.037c7db28a41cp+8 *
                (T * (0x1.c37d8424bd4ebp+1 +
                      T * (0x1.5218a2ee69f1p-8 +
                           T * (-0x1.44b60b48944b5p-20 +
                                T * (0x1.3d74955b82fccp-33 +
                                     T * -0x1.ed231bd8f8464p-48)))) -
                 0x1.964b889a02752p+14);
    cp[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.c37d8424bd4ebp+1 +
         T * (0x1.5218a2ee69f1p-7 +
              T * (-0x1.e71110ecde71p-19 +
                   T * (0x1.3d74955b82fccp-31 + T * -0x1.3435f1679b2bep-45))));
    s[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.c37d8424bd4ebp+1 * log(T) +
         T * (0x1.5218a2ee69f1p-7 +
              T * (-0x1.e71110ecde71p-20 +
                   T * (0x1.a7461c7a03fbbp-33 + T * -0x1.3435f1679b2bep-47))) +
         0x1.4ab9be87e5921p+2);
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
    s[sCH3OCH2O] =
        0x1.10565da70e7acp+7 *
        (0x1.1348a67cd6eb3p+3 * log(T) +
         T * (0x1.bce5f8dc8efb7p-7 +
              T * (-0x1.45403aad6b832p-19 +
                   T * (0x1.1d0de8327c1bap-32 + T * -0x1.a0b40cff718e4p-47))) -
         0x1.193d730d3dd62p+4);
    h[sOCH2OCHO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.8e027df421a57p+3 +
                         T * (0x1.0054a662e73e3p-8 +
                              T * (-0x1.fa3929569ac63p-21 +
                                   T * (0x1.f52acf3274824p-34 +
                                        T * -0x1.871dd401824fcp-48)))) -
                    0x1.5b4945aee632p+15);
    cp[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.8e027df421a57p+3 +
         T * (0x1.0054a662e73e3p-7 +
              T * (-0x1.7baadf00f414ap-19 +
                   T * (0x1.f52acf3274824p-32 + T * -0x1.e8e54901e2e3ap-46))));
    s[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.8e027df421a57p+3 * log(T) +
         T * (0x1.0054a662e73e3p-7 +
              T * (-0x1.7baadf00f414ap-20 +
                   T * (0x1.4e1c8a21a3018p-33 + T * -0x1.e8e54901e2e3ap-48))) -
         0x1.3445170351c41p+5);
    h[sHO2CH2OCHO] = 0x1.694b465e7f561p+6 *
                     (T * (0x1.1da534f326d3bp+4 +
                           T * (0x1.ee2b51f66b797p-9 +
                                T * (-0x1.ebf1e3bc14c98p-21 +
                                     T * (0x1.e9a9038bf91dap-34 +
                                          T * -0x1.7f917f3d1fb5p-48)))) -
                      0x1.f5140ce703afbp+15);
    cp[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.1da534f326d3bp+4 +
         T * (0x1.ee2b51f66b797p-8 +
              T * (-0x1.70f56acd0f972p-19 +
                   T * (0x1.e9a9038bf91dap-32 + T * -0x1.df75df0c67a23p-46))));
    s[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.1da534f326d3bp+4 * log(T) +
         T * (0x1.ee2b51f66b797p-8 +
              T * (-0x1.70f56acd0f972p-20 +
                   T * (0x1.4670ad07fb691p-33 + T * -0x1.df75df0c67a23p-48))) -
         0x1.02f1651e67ca1p+6);
    h[sHCO] = 0x1.1e87151786dcap+8 *
              (T * (0x1.f5c310b65b9c3p+1 +
                    T * (0x1.4aaae4ebb599dp-10 +
                         T * (-0x1.e052d03d363d1p-23 +
                              T * (0x1.d080fc6b6ed04p-36 +
                                   T * -0x1.acc53fd27579ap-50)))) +
               0x1.c8adbca9691a7p+11);
    cp[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.f5c310b65b9c3p+1 +
         T * (0x1.4aaae4ebb599dp-9 +
              T * (-0x1.683e1c2de8addp-21 +
                   T * (0x1.d080fc6b6ed04p-34 + T * -0x1.0bfb47e3896cp-47))));
    s[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.f5c310b65b9c3p+1 * log(T) +
         T * (0x1.4aaae4ebb599dp-9 +
              T * (-0x1.683e1c2de8addp-22 +
                   T * (0x1.35ab52f249e03p-35 + T * -0x1.0bfb47e3896cp-49))) +
         0x1.ca56b090d6fd5p+1);
    h[sO2CHO] = 0x1.10865abf22dacp+7 *
                (T * (0x1.cf687884a10fp+2 +
                      T * (0x1.2fa303b93cb4p-9 +
                           T * (-0x1.24f11ac33b5f1p-21 +
                                T * (0x1.1d8cef74001dep-34 +
                                     T * -0x1.b8e406b46d5d2p-49)))) -
                 0x1.243b0c154c986p+14);
    cp[sO2CHO] =
        0x1.10865abf22dacp+7 *
        (0x1.cf687884a10fp+2 +
         T * (0x1.2fa303b93cb4p-8 +
              T * (-0x1.b769a824d90e9p-20 +
                   T * (0x1.1d8cef74001dep-32 + T * -0x1.138e8430c45a3p-46))));
    s[sO2CHO] =
        0x1.10865abf22dacp+7 *
        (0x1.cf687884a10fp+2 * log(T) +
         T * (0x1.2fa303b93cb4p-8 +
              T * (-0x1.b769a824d90e9p-21 +
                   T * (0x1.7cbbe9f00027dp-34 + T * -0x1.138e8430c45a3p-48))) -
         0x1.9fb5d0b1deb4dp+2);
    h[sHOCHO] = 0x1.694b465e7f561p+7 *
                (T * (0x1.27490455d0163p+2 +
                      T * (0x1.a6aef168ba33p-9 +
                           T * (-0x1.99f588e7276bbp-21 +
                                T * (0x1.93b27c1812aa7p-34 +
                                     T * -0x1.3b3b9def78c84p-48)))) -
                 0x1.7335b33333333p+15);
    cp[sHOCHO] =
        0x1.694b465e7f561p+7 *
        (0x1.27490455d0163p+2 +
         T * (0x1.a6aef168ba33p-8 +
              T * (-0x1.337826ad5d90cp-19 +
                   T * (0x1.93b27c1812aa7p-32 + T * -0x1.8a0a856b56fa4p-46))));
    s[sHOCHO] =
        0x1.694b465e7f561p+7 *
        (0x1.27490455d0163p+2 * log(T) +
         T * (0x1.a6aef168ba33p-8 +
              T * (-0x1.337826ad5d90cp-20 +
                   T * (0x1.0d21a8100c71ap-33 + T * -0x1.8a0a856b56fa4p-48))) +
         0x1.b21dd451507fap-1);
    h[sHOCH2O] = 0x1.618d10efa02bbp+7 *
                 (T * (0x1.994b347c088f2p+2 +
                       T * (0x1.e75fa1fc6bad4p-9 +
                            T * (-0x1.c025f44665ca3p-21 +
                                 T * (0x1.a72e020a6d527p-34 +
                                      T * -0x1.3f9de4f37299cp-48)))) -
                  0x1.78bb89374bc6ap+14);
    cp[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.994b347c088f2p+2 +
         T * (0x1.e75fa1fc6bad4p-8 +
              T * (-0x1.501c7734cc57ap-19 +
                   T * (0x1.a72e020a6d527p-32 + T * -0x1.8f855e304f403p-46))));
    s[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.994b347c088f2p+2 * log(T) +
         T * (0x1.e75fa1fc6bad4p-8 +
              T * (-0x1.501c7734cc57ap-20 +
                   T * (0x1.1a1eac06f38c5p-33 + T * -0x1.8f855e304f403p-48))) -
         0x1.a8dfbcb3cffbbp+2);
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
    s[sHOCH2OCO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.6bf5abb377913p+3 * log(T) +
         T * (0x1.0bee9e8151789p-7 +
              T * (-0x1.87f62024c689bp-20 +
                   T * (0x1.56176cd31dbc3p-33 + T * -0x1.f1b2323925201p-48))) -
         0x1.c9a80b673c4f4p+4);
    h[sCH2OH] = 0x1.0bea1f388b9cfp+8 *
                (T * (0x1.45f610fe55051p+2 +
                      T * (0x1.85c8619b35011p-9 +
                           T * (-0x1.718aa63a0f48p-21 +
                                T * (0x1.6326b7cd0c28fp-34 +
                                     T * -0x1.0f1e351552636p-48)))) -
                 0x1.f84315b573eabp+11);
    cp[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.45f610fe55051p+2 +
         T * (0x1.85c8619b35011p-8 +
              T * (-0x1.1527fcab8b76p-19 +
                   T * (0x1.6326b7cd0c28fp-32 + T * -0x1.52e5c25aa6fc3p-46))));
    s[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.45f610fe55051p+2 * log(T) +
         T * (0x1.85c8619b35011p-8 +
              T * (-0x1.1527fcab8b76p-20 +
                   T * (0x1.d988f511658bfp-34 + T * -0x1.52e5c25aa6fc3p-48))) -
         0x1.d8cf6ab6d818ep+0);
    h[sCH3O2H] = 0x1.5a2209942d73ap+7 *
                 (T * (0x1.f0fbfed405fc2p+2 +
                       T * (0x1.1a4bd571592f3p-8 +
                            T * (-0x1.0aa6ef7a10ea9p-20 +
                                 T * (0x1.01a2f1ef778aep-33 +
                                      T * -0x1.8cce416275897p-48)))) -
                  0x1.1de7fe5c91d15p+14);
    cp[sCH3O2H] =
        0x1.5a2209942d73ap+7 *
        (0x1.f0fbfed405fc2p+2 +
         T * (0x1.1a4bd571592f3p-7 +
              T * (-0x1.8ffa6737195fdp-19 +
                   T * (0x1.01a2f1ef778aep-31 + T * -0x1.f001d1bb12ebcp-46))));
    s[sCH3O2H] =
        0x1.5a2209942d73ap+7 *
        (0x1.f0fbfed405fc2p+2 * log(T) +
         T * (0x1.1a4bd571592f3p-7 +
              T * (-0x1.8ffa6737195fdp-20 +
                   T * (0x1.5783ed3f4a0e8p-33 + T * -0x1.f001d1bb12ebcp-48))) -
         0x1.ccc6ca1e9ca34p+3);
    h[sCH3O2] = 0x1.618d10efa02bbp+7 *
                (T * (0x1.9eb37c0e18719p+2 +
                      T * (0x1.e7d9c6efe3df3p-9 +
                           T * (-0x1.c39867846fc2p-21 +
                                T * (0x1.ac5842a43d63p-34 +
                                     T * -0x1.4485a4457a2e5p-48)))) -
                 0x1.86b6aea747d8p+10);
    cp[sCH3O2] =
        0x1.618d10efa02bbp+7 *
        (0x1.9eb37c0e18719p+2 +
         T * (0x1.e7d9c6efe3df3p-8 +
              T * (-0x1.52b24da353d18p-19 +
                   T * (0x1.ac5842a43d63p-32 + T * -0x1.95a70d56d8b9ep-46))));
    s[sCH3O2] =
        0x1.618d10efa02bbp+7 *
        (0x1.9eb37c0e18719p+2 * log(T) +
         T * (0x1.e7d9c6efe3df3p-8 +
              T * (-0x1.52b24da353d18p-20 +
                   T * (0x1.1d902c6d7e42p-33 + T * -0x1.95a70d56d8b9ep-48))) -
         0x1.063b8fd8d0912p+3);
    h[sHO2CHO] = 0x1.0c18906567f93p+7 *
                 (T * (0x1.3c005153da9dap+3 +
                       T * (0x1.3085a24b83deep-9 +
                            T * (-0x1.2b454b0918308p-21 +
                                 T * (0x1.275b0ce68fbbap-34 +
                                      T * -0x1.cc007efa1664dp-49)))) -
                  0x1.29447fcb923a3p+15);
    cp[sHO2CHO] =
        0x1.0c18906567f93p+7 *
        (0x1.3c005153da9dap+3 +
         T * (0x1.3085a24b83deep-8 +
              T * (-0x1.c0e7f08da448cp-20 +
                   T * (0x1.275b0ce68fbbap-32 + T * -0x1.1f804f5c4dffp-46))));
    s[sHO2CHO] =
        0x1.0c18906567f93p+7 *
        (0x1.3c005153da9dap+3 * log(T) +
         T * (0x1.3085a24b83deep-8 +
              T * (-0x1.c0e7f08da448cp-21 +
                   T * (0x1.89cebbde14fa3p-34 + T * -0x1.1f804f5c4dffp-48))) -
         0x1.67e713f077cccp+4);
    h[sCH3OCH2O2H] = 0x1.aa02db5960032p+6 *
                     (T * (0x1.d52cbfe70b3c7p+3 +
                           T * (0x1.8b0033a0792e5p-8 +
                                T * (-0x1.72165456b5d5p-20 +
                                     T * (0x1.61e3c4323963ep-33 +
                                          T * -0x1.0da1a7fc29939p-47)))) -
                      0x1.47433cb923a2ap+15);
    cp[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.d52cbfe70b3c7p+3 +
         T * (0x1.8b0033a0792e5p-7 +
              T * (-0x1.1590bf41085fcp-18 +
                   T * (0x1.61e3c4323963ep-31 + T * -0x1.510a11fb33f87p-45))));
    s[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.d52cbfe70b3c7p+3 * log(T) +
         T * (0x1.8b0033a0792e5p-7 +
              T * (-0x1.1590bf41085fcp-19 +
                   T * (0x1.d7da5aeda1da8p-33 + T * -0x1.510a11fb33f87p-47))) -
         0x1.8108dba68463ap+5);
    h[sCH2] = 0x1.286500f4d652p+9 *
              (T * (0x1.92ba938f3e76ep+1 +
                    T * (0x1.8e072a85d6696p-10 +
                         T * (-0x1.64a6efa038311p-22 +
                              T * (0x1.4aeab96f1c774p-35 +
                                   T * -0x1.ee3864013fb42p-50)))) +
               0x1.67b2856041893p+15);
    cp[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.92ba938f3e76ep+1 +
         T * (0x1.8e072a85d6696p-9 +
              T * (-0x1.0b7d33b82a24dp-20 +
                   T * (0x1.4aeab96f1c774p-33 + T * -0x1.34e33e80c7d09p-47))));
    s[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.92ba938f3e76ep+1 * log(T) +
         T * (0x1.8e072a85d6696p-9 +
              T * (-0x1.0b7d33b82a24dp-21 +
                   T * (0x1.b938f73ed09fp-35 + T * -0x1.34e33e80c7d09p-49))) +
         0x1.2e4c77473447p+2);
    h[sCH] = 0x1.3f584141b0716p+9 *
             (T * (0x1.42ae0f7263cabp+1 +
                   T * (0x1.cec790dd3954dp-11 +
                        T * (-0x1.4a567b34ccf05p-23 +
                             T * (0x1.04c23d53c0effp-36 +
                                  T * -0x1.81ef49d937f0bp-51)))) +
              0x1.1522c4dd2f1aap+16);
    cp[sCH] =
        0x1.3f584141b0716p+9 *
        (0x1.42ae0f7263cabp+1 +
         T * (0x1.cec790dd3954dp-10 +
              T * (-0x1.ef81b8cf33688p-22 +
                   T * (0x1.04c23d53c0effp-34 + T * -0x1.e26b1c4f85ecdp-49))));
    s[sCH] =
        0x1.3f584141b0716p+9 *
        (0x1.42ae0f7263cabp+1 * log(T) +
         T * (0x1.cec790dd3954dp-10 +
              T * (-0x1.ef81b8cf33688p-23 +
                   T * (0x1.5bada71a56954p-36 + T * -0x1.e26b1c4f85ecdp-51))) +
         0x1.d9ee8442198p+2);
    h[sCH3OCHO] = 0x1.14e89ec288e65p+7 *
                  (T * (0x1.9559d8b96a198p+2 +
                        T * (0x1.b9e19fe5049e7p-8 +
                             T * (-0x1.b15978bf1ac48p-20 +
                                  T * (0x1.adc133dc1c597p-33 +
                                       T * -0x1.512b88af43bccp-47)))) -
                   0x1.6ddf4de00d1b7p+15);
    cp[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.9559d8b96a198p+2 +
         T * (0x1.b9e19fe5049e7p-7 +
              T * (-0x1.45031a8f54136p-18 +
                   T * (0x1.adc133dc1c597p-31 + T * -0x1.a5766adb14abep-45))));
    s[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.9559d8b96a198p+2 * log(T) +
         T * (0x1.b9e19fe5049e7p-7 +
              T * (-0x1.45031a8f54136p-19 +
                   T * (0x1.1e80cd3d683bap-32 + T * -0x1.a5766adb14abep-47))) -
         0x1.ba9656f9b6e5bp+2);
    h[sCH2OCHO] = 0x1.19a2d4e8daa2p+7 *
                  (T * (0x1.431282b98343fp+3 +
                        T * (0x1.d7c9020d86866p-9 +
                             T * (-0x1.d0f44d6bbb74bp-21 +
                                  T * (0x1.cbb7bb276a362p-34 +
                                       T * -0x1.6672bbb9b9dd4p-48)))) -
                   0x1.715b9b71758e2p+14);
    cp[sCH2OCHO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.431282b98343fp+3 +
         T * (0x1.d7c9020d86866p-8 +
              T * (-0x1.5cb73a10cc978p-19 +
                   T * (0x1.cbb7bb276a362p-32 + T * -0x1.c00f6aa828549p-46))));
    s[sCH2OCHO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.431282b98343fp+3 * log(T) +
         T * (0x1.d7c9020d86866p-8 +
              T * (-0x1.5cb73a10cc978p-20 +
                   T * (0x1.327a7cc4f1797p-33 + T * -0x1.c00f6aa828549p-48))) -
         0x1.b1d4a771c970fp+4);
    h[sC2H4O1X2] = 0x1.797bdf23f83eep+7 *
                   (T * (0x1.5f47e9383d5e3p+2 +
                         T * (0x1.8abac40bb3119p-8 +
                              T * (-0x1.83c5b5176ec53p-20 +
                                   T * (0x1.80fc1606a7dbcp-33 +
                                        T * -0x1.2e4668b73fa5ap-47)))) -
                    0x1.1ee3669ad42c4p+13);
    cp[sC2H4O1X2] =
        0x1.797bdf23f83eep+7 *
        (0x1.5f47e9383d5e3p+2 +
         T * (0x1.8abac40bb3119p-7 +
              T * (-0x1.22d447d19313ep-18 +
                   T * (0x1.80fc1606a7dbcp-31 + T * -0x1.79d802e50f8fp-45))));
    s[sC2H4O1X2] =
        0x1.797bdf23f83eep+7 *
        (0x1.5f47e9383d5e3p+2 * log(T) +
         T * (0x1.8abac40bb3119p-7 +
              T * (-0x1.22d447d19313ep-19 +
                   T * (0x1.00a80eaf1a928p-32 + T * -0x1.79d802e50f8fp-47))) -
         0x1.c51e12a51e322p+2);
    h[sC2H4] = 0x1.286500f4d652p+8 *
               (T * (0x1.fef431eb8a7d2p+1 +
                     T * (0x1.57850e518f5a5p-8 +
                          T * (-0x1.4c9c520463b8ap-20 +
                               T * (0x1.46e680b348f8ep-33 +
                                    T * -0x1.fda295a4b6b8ap-48)))) +
                0x1.0aca8941c8217p+12);
    cp[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.fef431eb8a7d2p+1 +
         T * (0x1.57850e518f5a5p-7 +
              T * (-0x1.f2ea7b069594fp-19 +
                   T * (0x1.46e680b348f8ep-31 + T * -0x1.3e859d86f2336p-45))));
    s[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.fef431eb8a7d2p+1 * log(T) +
         T * (0x1.57850e518f5a5p-7 +
              T * (-0x1.f2ea7b069594fp-20 +
                   T * (0x1.b3de00ef0bf68p-33 + T * -0x1.3e859d86f2336p-47))) -
         0x1.138a2b5ef5d12p-2);
    h[sC2H3] = 0x1.337122eb1b451p+8 *
               (T * (0x1.099dfc1f1fd5ep+2 +
                     T * (0x1.ee27ca0035204p-9 +
                          T * (-0x1.d6a729b0fd579p-21 +
                               T * (0x1.c95e4914c2cf2p-34 +
                                    T * -0x1.61ab63e725efcp-48)))) +
                0x1.088146a7ef9dbp+15);
    cp[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.099dfc1f1fd5ep+2 +
         T * (0x1.ee27ca0035204p-8 +
              T * (-0x1.60fd5f44be01bp-19 +
                   T * (0x1.c95e4914c2cf2p-32 + T * -0x1.ba163ce0ef6bbp-46))));
    s[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.099dfc1f1fd5ep+2 * log(T) +
         T * (0x1.ee27ca0035204p-8 +
              T * (-0x1.60fd5f44be01bp-20 +
                   T * (0x1.30e9860dd734cp-33 + T * -0x1.ba163ce0ef6bbp-48))) +
         0x1.ba6639f0bc962p+0);
    h[sC2H5] = 0x1.1e1d11b23ceecp+8 *
               (T * (0x1.78d26f389cff8p+2 +
                     T * (0x1.51c314bede34ap-8 +
                          T * (-0x1.3659e279efe9p-20 +
                               T * (0x1.24be9ba4a8e5p-33 +
                                    T * -0x1.b9bb3564a607cp-48)))) +
                0x1.67946631f8a09p+13);
    cp[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.78d26f389cff8p+2 +
         T * (0x1.51c314bede34ap-7 +
              T * (-0x1.d186d3b6e7dd8p-19 +
                   T * (0x1.24be9ba4a8e5p-31 + T * -0x1.1415015ee7c4dp-45))));
    s[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.78d26f389cff8p+2 * log(T) +
         T * (0x1.51c314bede34ap-7 +
              T * (-0x1.d186d3b6e7dd8p-20 +
                   T * (0x1.86537a30e1315p-33 + T * -0x1.1415015ee7c4dp-47))) -
         0x1.0fe3791bcab6cp+3);
    h[sOHY] = 0x1.e8db171241376p+8 *
              (T * (0x1.70fd4bf0995abp+1 +
                    T * (0x1.09cea9d476d0ap-11 +
                         T * (-0x1.45f86caf57da5p-24 +
                              T * (0x1.7e9303f6caa02p-38 +
                                   T * -0x1.d8d1745e5c6dfp-54)))) +
               0x1.88b2p+15);
    cp[sOHY] =
        0x1.e8db171241376p+8 *
        (0x1.70fd4bf0995abp+1 +
         T * (0x1.09cea9d476d0ap-10 +
              T * (-0x1.e8f4a30703c77p-23 +
                   T * (0x1.7e9303f6caa02p-36 + T * -0x1.2782e8baf9c4bp-51))));
    s[sOHY] =
        0x1.e8db171241376p+8 *
        (0x1.70fd4bf0995abp+1 * log(T) +
         T * (0x1.09cea9d476d0ap-10 +
              T * (-0x1.e8f4a30703c77p-24 +
                   T * (0x1.fe195a9e63803p-38 + T * -0x1.2782e8baf9c4bp-53))) +
         0x1.66202539756c9p+2);
    h[sCH3OCO] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.c01c2c50c0e2bp+2 +
                       T * (0x1.4e28c046344a7p-8 +
                            T * (-0x1.4727211776365p-20 +
                                 T * (0x1.44113e38a28dp-33 +
                                      T * -0x1.fc182b4f163a5p-48)))) -
                  0x1.61564fdf3b646p+14);
    cp[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.c01c2c50c0e2bp+2 +
         T * (0x1.4e28c046344a7p-7 +
              T * (-0x1.eabab1a331518p-19 +
                   T * (0x1.44113e38a28dp-31 + T * -0x1.3d8f1b116de47p-45))));
    s[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.c01c2c50c0e2bp+2 * log(T) +
         T * (0x1.4e28c046344a7p-7 +
              T * (-0x1.eabab1a331518p-20 +
                   T * (0x1.b016fda0d8bcp-33 + T * -0x1.3d8f1b116de47p-47))) -
         0x1.21af870699f81p+3);
    h[sHOCH2O2H] = 0x1.03a80e5d6edeap+7 *
                   (T * (0x1.742c18570eda9p+3 +
                         T * (0x1.d4ab8605b40e3p-9 +
                              T * (-0x1.abc514cae8f3cp-21 +
                                   T * (0x1.922be348b4157p-34 +
                                        T * -0x1.2eedefbcafa43p-48)))) -
                    0x1.50c7d930be0dfp+15);
    cp[sHOCH2O2H] =
        0x1.03a80e5d6edeap+7 *
        (0x1.742c18570eda9p+3 +
         T * (0x1.d4ab8605b40e3p-8 +
              T * (-0x1.40d3cf982eb6dp-19 +
                   T * (0x1.922be348b4157p-32 + T * -0x1.7aa96babdb8d3p-46))));
    s[sHOCH2O2H] =
        0x1.03a80e5d6edeap+7 *
        (0x1.742c18570eda9p+3 * log(T) +
         T * (0x1.d4ab8605b40e3p-8 +
              T * (-0x1.40d3cf982eb6dp-20 +
                   T * (0x1.0c1d4230780e5p-33 + T * -0x1.7aa96babdb8d3p-48))) -
         0x1.036bdf8f47304p+5);
    h[sHOCH2O2] = 0x1.07cf090c05134p+7 *
                  (T * (0x1.2174673accafap+3 +
                        T * (0x1.d4ba91f6b4eb1p-9 +
                             T * (-0x1.a8235fa3748d7p-21 +
                                  T * (0x1.8bea7d4b2d6dcp-34 +
                                       T * -0x1.28846f9c6a7adp-48)))) -
                   0x1.85b5f4538ef35p+14);
    cp[sHOCH2O2] =
        0x1.07cf090c05134p+7 *
        (0x1.2174673accafap+3 +
         T * (0x1.d4ba91f6b4eb1p-8 +
              T * (-0x1.3e1a87ba976a1p-19 +
                   T * (0x1.8bea7d4b2d6dcp-32 + T * -0x1.72a58b8385198p-46))));
    s[sHOCH2O2] =
        0x1.07cf090c05134p+7 *
        (0x1.2174673accafap+3 * log(T) +
         T * (0x1.d4ba91f6b4eb1p-8 +
              T * (-0x1.3e1a87ba976a1p-20 +
                   T * (0x1.07f1a8dcc8f3dp-33 + T * -0x1.72a58b8385198p-48))) -
         0x1.16bca2120e1f8p+4);
    h[sC2H3CO] = 0x1.2e0c24d5d505ep+7 *
                 (T * (0x1.2bfd5a1dd8aep+3 +
                       T * (0x1.034acb88dd925p-8 +
                            T * (-0x1.de2b81aef2691p-21 +
                                 T * (0x1.c406b4c1f68abp-34 +
                                      T * -0x1.5585cff0bcbf4p-48)))) +
                  0x1.e26c7d2c7b891p+10);
    cp[sC2H3CO] =
        0x1.2e0c24d5d505ep+7 *
        (0x1.2bfd5a1dd8aep+3 +
         T * (0x1.034acb88dd925p-7 +
              T * (-0x1.66a0a14335cedp-19 +
                   T * (0x1.c406b4c1f68abp-32 + T * -0x1.aae743ecebefp-46))));
    s[sC2H3CO] =
        0x1.2e0c24d5d505ep+7 *
        (0x1.2bfd5a1dd8aep+3 * log(T) +
         T * (0x1.034acb88dd925p-7 +
              T * (-0x1.66a0a14335cedp-20 +
                   T * (0x1.2d59cdd6a45c7p-33 + T * -0x1.aae743ecebefp-48))) -
         0x1.816da5f5c86eap+4);
    h[sC2H3CHO] = 0x1.289dd93acb778p+7 *
                  (T * (0x1.4d6451838052ap+3 +
                        T * (0x1.36f4d024dcb58p-8 +
                             T * (-0x1.26a97da44fdb1p-20 +
                                  T * (0x1.1bd3d69c28ed4p-33 +
                                       T * -0x1.b2a2124b268eap-48)))) -
                   0x1.d398398c7e282p+13);
    cp[sC2H3CHO] =
        0x1.289dd93acb778p+7 *
        (0x1.4d6451838052ap+3 +
         T * (0x1.36f4d024dcb58p-7 +
              T * (-0x1.b9fe3c7677c8ap-19 +
                   T * (0x1.1bd3d69c28ed4p-31 + T * -0x1.0fa54b6ef8192p-45))));
    s[sC2H3CHO] =
        0x1.289dd93acb778p+7 *
        (0x1.4d6451838052ap+3 * log(T) +
         T * (0x1.36f4d024dcb58p-7 +
              T * (-0x1.b9fe3c7677c8ap-20 +
                   T * (0x1.7a6fc8d03691bp-33 + T * -0x1.0fa54b6ef8192p-47))) -
         0x1.eb937b21df4dep+4);
    h[sCH3CO] = 0x1.8252e1709b12p+7 *
                (T * (0x1.5413ee5eedcc2p+2 +
                      T * (0x1.2c9b3d0980196p-8 +
                           T * (-0x1.2940405db0acbp-20 +
                                T * (0x1.28944b6b4c6ccp-33 +
                                     T * -0x1.d3b01afc912e5p-48)))) -
                 0x1.c7a1532617c1cp+11);
    cp[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.5413ee5eedcc2p+2 +
         T * (0x1.2c9b3d0980196p-7 +
              T * (-0x1.bde0608c8903p-19 +
                   T * (0x1.28944b6b4c6ccp-31 + T * -0x1.244e10dddabcfp-45))));
    s[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.5413ee5eedcc2p+2 * log(T) +
         T * (0x1.2c9b3d0980196p-7 +
              T * (-0x1.bde0608c8903p-20 +
                   T * (0x1.8b70648f1091p-33 + T * -0x1.244e10dddabcfp-47))) -
         0x1.acfe55051512bp+0);
    h[sC2H5O2H] = 0x1.0bea1f388b9cfp+7 *
                  (T * (0x1.6760dc189a2acp+3 +
                        T * (0x1.8acbba41f485fp-8 +
                             T * (-0x1.62fcf904a252bp-20 +
                                  T * (0x1.4a44b2b1ccbe7p-33 +
                                       T * -0x1.edd26a7ac84f8p-48)))) -
                   0x1.8377032ca57a8p+14);
    cp[sC2H5O2H] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.6760dc189a2acp+3 +
         T * (0x1.8acbba41f485fp-7 +
              T * (-0x1.0a3dbac379bep-18 +
                   T * (0x1.4a44b2b1ccbe7p-31 + T * -0x1.34a3828cbd31bp-45))));
    s[sC2H5O2H] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.6760dc189a2acp+3 * log(T) +
         T * (0x1.8acbba41f485fp-7 +
              T * (-0x1.0a3dbac379bep-19 +
                   T * (0x1.b85b98ed10fdfp-33 + T * -0x1.34a3828cbd31bp-47))) -
         0x1.047c5c71f0de2p+5);
    h[sC2H5O2] = 0x1.10565da70e7acp+7 *
                 (T * (0x1.1c706dfc31ea3p+3 +
                       T * (0x1.bd19211369f7ep-8 +
                            T * (-0x1.b771aa745f58dp-20 +
                                 T * (0x1.b3986a28a8867p-33 +
                                      T * -0x1.553613222afe5p-47)))) -
                  0x1.d1112e9ccb7d4p+12);
    cp[sC2H5O2] =
        0x1.10565da70e7acp+7 *
        (0x1.1c706dfc31ea3p+3 +
         T * (0x1.bd19211369f7ep-7 +
              T * (-0x1.49953fd74782ap-18 +
                   T * (0x1.b3986a28a8867p-31 + T * -0x1.aa8397eab5bdep-45))));
    s[sC2H5O2] =
        0x1.10565da70e7acp+7 *
        (0x1.1c706dfc31ea3p+3 * log(T) +
         T * (0x1.bd19211369f7ep-7 +
              T * (-0x1.49953fd74782ap-19 +
                   T * (0x1.22659c1b1b045p-32 + T * -0x1.aa8397eab5bdep-47))) -
         0x1.3143844eaeb9cp+4);
    h[sCH2CO] = 0x1.8b966bc75924bp+7 *
                (T * (0x1.56f4d64b7ba97p+2 +
                      T * (0x1.c7e54a9db7235p-9 +
                           T * (-0x1.d9e1fdb377569p-21 +
                                T * (0x1.ff58e3be99b1fp-34 +
                                     T * -0x1.bcccbbb8378e7p-48)))) -
                 0x1.edef0ac5c13fdp+12);
    cp[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.56f4d64b7ba97p+2 +
         T * (0x1.c7e54a9db7235p-8 +
              T * (-0x1.63697e469980fp-19 +
                   T * (0x1.ff58e3be99b1fp-32 + T * -0x1.15fff55322b9p-45))));
    s[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.56f4d64b7ba97p+2 * log(T) +
         T * (0x1.c7e54a9db7235p-8 +
              T * (-0x1.63697e469980fp-20 +
                   T * (0x1.54e5ed29bbcbfp-33 + T * -0x1.15fff55322b9p-47))) -
         0x1.fe1ce95a4c26dp+1);
    h[sCH2CHO] = 0x1.8252e1709b12p+7 *
                 (T * (0x1.a2839e701815fp+2 +
                       T * (0x1.ff5651f2414abp-9 +
                            T * (-0x1.eea9534005b03p-21 +
                                 T * (0x1.e617c77e3fbcep-34 +
                                      T * -0x1.7af50283304bfp-48)))) -
                  0x1.29258ab0c88a4p+10);
    cp[sCH2CHO] =
        0x1.8252e1709b12p+7 *
        (0x1.a2839e701815fp+2 +
         T * (0x1.ff5651f2414abp-8 +
              T * (-0x1.72fefe7004442p-19 +
                   T * (0x1.e617c77e3fbcep-32 + T * -0x1.d9b24323fc5eep-46))));
    s[sCH2CHO] =
        0x1.8252e1709b12p+7 *
        (0x1.a2839e701815fp+2 * log(T) +
         T * (0x1.ff5651f2414abp-8 +
              T * (-0x1.72fefe7004442p-20 +
                   T * (0x1.440fda542a7dfp-33 + T * -0x1.d9b24323fc5eep-48))) -
         0x1.1711ba1712963p+3);
    h[sCH3CO3H] = 0x1.b54dcef92163bp+6 *
                  (T * (0x1.90318c9fb6135p+3 +
                        T * (0x1.36925cb7b19fbp-8 +
                             T * (-0x1.27a3909a8f56bp-20 +
                                  T * (0x1.1dab885da5166p-33 +
                                       T * -0x1.b6726016f941dp-48)))) -
                   0x1.674357318fc5p+15);
    cp[sCH3CO3H] =
        0x1.b54dcef92163bp+6 *
        (0x1.90318c9fb6135p+3 +
         T * (0x1.36925cb7b19fbp-7 +
              T * (-0x1.bb7558e7d7021p-19 +
                   T * (0x1.1dab885da5166p-31 + T * -0x1.12077c0e5bc92p-45))));
    s[sCH3CO3H] =
        0x1.b54dcef92163bp+6 *
        (0x1.90318c9fb6135p+3 * log(T) +
         T * (0x1.36925cb7b19fbp-7 +
              T * (-0x1.bb7558e7d7021p-20 +
                   T * (0x1.7ce4b5d231733p-33 + T * -0x1.12077c0e5bc92p-47))) -
         0x1.2f5b547750997p+5);
    h[sCH3CO3] = 0x1.bb2d881f811c1p+6 *
                 (T * (0x1.68126e2c2d857p+3 +
                       T * (0x1.112bdacff3e17p-8 +
                            T * (-0x1.029b17456022cp-20 +
                                 T * (0x1.f1d6ba605933p-34 +
                                      T * -0x1.7cf99473c586dp-48)))) -
                  0x1.969f6f0068db9p+14);
    cp[sCH3CO3] =
        0x1.bb2d881f811c1p+6 *
        (0x1.68126e2c2d857p+3 +
         T * (0x1.112bdacff3e17p-7 +
              T * (-0x1.83e8a2e810342p-19 +
                   T * (0x1.f1d6ba605933p-32 + T * -0x1.dc37f990b6e88p-46))));
    s[sCH3CO3] =
        0x1.bb2d881f811c1p+6 *
        (0x1.68126e2c2d857p+3 * log(T) +
         T * (0x1.112bdacff3e17p-7 +
              T * (-0x1.83e8a2e810342p-20 +
                   T * (0x1.4be47c403b775p-33 + T * -0x1.dc37f990b6e88p-48))) -
         0x1.da3156d4f8eb4p+4);
    h[sH2CC] = 0x1.3f584141b0716p+8 *
               (T * (0x1.11cb4f1e4b44ap+2 +
                     T * (0x1.37b524c4c3598p-9 +
                          T * (-0x1.23b7c98935fbdp-21 +
                               T * (0x1.17f76d60ce1fp-34 +
                                    T * -0x1.ad1216e577b4cp-49)))) +
                0x1.7979604189375p+15);
    cp[sH2CC] =
        0x1.3f584141b0716p+8 *
        (0x1.11cb4f1e4b44ap+2 +
         T * (0x1.37b524c4c3598p-8 +
              T * (-0x1.b593ae4dd0f9bp-20 +
                   T * (0x1.17f76d60ce1fp-32 + T * -0x1.0c2b4e4f6ad0fp-46))));
    s[sH2CC] =
        0x1.3f584141b0716p+8 *
        (0x1.11cb4f1e4b44ap+2 * log(T) +
         T * (0x1.37b524c4c3598p-8 +
              T * (-0x1.b593ae4dd0f9bp-21 +
                   T * (0x1.7549e72bbd7ebp-34 + T * -0x1.0c2b4e4f6ad0fp-48))) +
         0x1.47cd253747141p-1);
    h[sC2H2] = 0x1.3f584141b0716p+8 *
               (T * (0x1.2a29881969888p+2 +
                     T * (0x1.40135d1f55decp-9 +
                          T * (-0x1.1fd08375c80a6p-21 +
                               T * (0x1.0f8d276de39d5p-34 +
                                    T * -0x1.8f812567807a7p-49)))) +
                0x1.927d9de69ad43p+14);
    cp[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.2a29881969888p+2 +
         T * (0x1.40135d1f55decp-8 +
              T * (-0x1.afb8c530ac0f9p-20 +
                   T * (0x1.0f8d276de39d5p-32 + T * -0x1.f3616ec16099p-47))));
    s[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.2a29881969888p+2 * log(T) +
         T * (0x1.40135d1f55decp-8 +
              T * (-0x1.afb8c530ac0f9p-21 +
                   T * (0x1.6a1189e7da271p-34 + T * -0x1.f3616ec16099p-49))) -
         0x1.ffcafaba9bc7cp+1);
    h[sOCHO] = 0x1.716240449121bp+7 *
               (T * (0x1.093658f7bde73p+2 +
                     T * (0x1.6ed4972d16005p-9 +
                          T * (-0x1.658b9a5a20dc1p-21 +
                               T * (0x1.5ba487eeb3f9p-34 +
                                    T * -0x1.0b7fa843dcc72p-48)))) -
                0x1.0d77f46dc5d64p+14);
    cp[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.093658f7bde73p+2 +
         T * (0x1.6ed4972d16005p-8 +
              T * (-0x1.0c28b3c398a51p-19 +
                   T * (0x1.5ba487eeb3f9p-32 + T * -0x1.4e5f9254d3f8ep-46))));
    s[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.093658f7bde73p+2 * log(T) +
         T * (0x1.6ed4972d16005p-8 +
              T * (-0x1.0c28b3c398a51p-20 +
                   T * (0x1.cf860a939aa15p-34 + T * -0x1.4e5f9254d3f8ep-48))) +
         0x1.44fa72d80eca5p+2);
    h[sSC2H4OH] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.96f0671b4225cp+2 +
                        T * (0x1.977d9aec72321p-8 +
                             T * (-0x1.83874b391e7f7p-20 +
                                  T * (0x1.785316267e31bp-33 +
                                       T * -0x1.22e7f2a3625ecp-47)))) -
                   0x1.251808d8ec95cp+13);
    cp[sSC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.96f0671b4225cp+2 +
         T * (0x1.977d9aec72321p-7 +
              T * (-0x1.22a5786ad6df9p-18 +
                   T * (0x1.785316267e31bp-31 + T * -0x1.6ba1ef4c3af67p-45))));
    s[sSC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.96f0671b4225cp+2 * log(T) +
         T * (0x1.977d9aec72321p-7 +
              T * (-0x1.22a5786ad6df9p-19 +
                   T * (0x1.f5c41d88a8424p-33 + T * -0x1.6ba1ef4c3af67p-47))) -
         0x1.834495dc1fd1p+2);
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
    s[sC2H5OH] =
        0x1.68f6f3733777dp+7 *
        (0x1.a3fef5a964e8bp+2 * log(T) +
         T * (0x1.f236422024d26p-7 +
              T * (-0x1.69b1fd1ef62b7p-19 +
                   T * (0x1.3c049c20e638bp-32 + T * -0x1.cdfa32615ab65p-47))) -
         0x1.2f22fb422b1dcp+3);
    h[sCH3CO2] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.114c92d5b0204p+3 +
                       T * (0x1.10f1031dd0d03p-8 +
                            T * (-0x1.fd87a678ef3dp-21 +
                                 T * (0x1.e5e7727bd2f6ep-34 +
                                      T * -0x1.71791675606e8p-48)))) -
                  0x1.d084456d5cfabp+14);
    cp[sCH3CO2] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.114c92d5b0204p+3 +
         T * (0x1.10f1031dd0d03p-7 +
              T * (-0x1.7e25bcdab36dcp-19 +
                   T * (0x1.e5e7727bd2f6ep-32 + T * -0x1.cdd75c12b88a1p-46))));
    s[sCH3CO2] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.114c92d5b0204p+3 * log(T) +
         T * (0x1.10f1031dd0d03p-7 +
              T * (-0x1.7e25bcdab36dcp-20 +
                   T * (0x1.43efa1a7e1f9fp-33 + T * -0x1.cdd75c12b88a1p-48))) -
         0x1.4636b3354c122p+4);
    h[sCH3CHO] = 0x1.797bdf23f83eep+7 *
                 (T * (0x1.59dcf38b7d772p+2 +
                       T * (0x1.80242581cd52bp-8 +
                            T * (-0x1.7a2a05a17d829p-20 +
                                 T * (0x1.77e1ab96723c6p-33 +
                                      T * -0x1.2753ba556825fp-47)))) -
                  0x1.61047ced91687p+14);
    cp[sCH3CHO] =
        0x1.797bdf23f83eep+7 *
        (0x1.59dcf38b7d772p+2 +
         T * (0x1.80242581cd52bp-7 +
              T * (-0x1.1b9f84391e21fp-18 +
                   T * (0x1.77e1ab96723c6p-31 + T * -0x1.7128a8eac22f7p-45))));
    s[sCH3CHO] =
        0x1.797bdf23f83eep+7 *
        (0x1.59dcf38b7d772p+2 * log(T) +
         T * (0x1.80242581cd52bp-7 +
              T * (-0x1.1b9f84391e21fp-19 +
                   T * (0x1.f52ce4c898508p-33 + T * -0x1.7128a8eac22f7p-47))) -
         0x1.bd8a9519d8186p+1);
    h[sHCCO] = 0x1.954e7e0035abdp+7 *
               (T * (0x1.7a8bf952bcbb6p+2 +
                     T * (0x1.e6d016f9d9afcp-10 +
                          T * (-0x1.d1c74be1e87bdp-22 +
                               T * (0x1.c60a28f283f63p-35 +
                                    T * -0x1.5e2207018734fp-49)))) +
                0x1.2e7e8538ef34dp+14);
    cp[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.7a8bf952bcbb6p+2 +
         T * (0x1.e6d016f9d9afcp-9 +
              T * (-0x1.5d5578e96e5cep-20 +
                   T * (0x1.c60a28f283f63p-33 + T * -0x1.b5aa88c1e9022p-47))));
    s[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.7a8bf952bcbb6p+2 * log(T) +
         T * (0x1.e6d016f9d9afcp-9 +
              T * (-0x1.5d5578e96e5cep-21 +
                   T * (0x1.2eb170a1ad4edp-34 + T * -0x1.b5aa88c1e9022p-49))) -
         0x1.605cf0fc81f13p+2);
    h[sHCOH] = 0x1.14e89ec288e65p+8 *
               (T * (0x1.25fff0bb946b3p+3 +
                     T * (0x1.8e7cf50f06665p-11 +
                          T * (-0x1.c1418f821e0bdp-23 +
                               T * (0x1.e296bb8b39281p-36 +
                                    T * -0x1.8d8f2216f35e4p-50)))) +
                0x1.e85a55bab2181p+12);
    cp[sHCOH] =
        0x1.14e89ec288e65p+8 *
        (0x1.25fff0bb946b3p+3 +
         T * (0x1.8e7cf50f06665p-10 +
              T * (-0x1.50f12ba19688ep-21 +
                   T * (0x1.e296bb8b39281p-34 + T * -0x1.f0f2ea9cb035cp-48))));
    s[sHCOH] =
        0x1.14e89ec288e65p+8 *
        (0x1.25fff0bb946b3p+3 * log(T) +
         T * (0x1.8e7cf50f06665p-10 +
              T * (-0x1.50f12ba19688ep-22 +
                   T * (0x1.41b9d25cd0c56p-35 + T * -0x1.f0f2ea9cb035cp-50))) -
         0x1.b57ea7701bf1dp+4);
    h[sC2H3OH] = 0x1.797bdf23f83eep+7 *
                 (T * (0x1.0a6e70ec26596p+3 +
                       T * (0x1.0741027d558fbp-8 +
                            T * (-0x1.d8517ab011e6dp-21 +
                                 T * (0x1.b60ea65abe297p-34 +
                                      T * -0x1.467e9f9edc3ddp-48)))) -
                  0x1.1e48930be0dedp+14);
    cp[sC2H3OH] =
        0x1.797bdf23f83eep+7 *
        (0x1.0a6e70ec26596p+3 +
         T * (0x1.0741027d558fbp-7 +
              T * (-0x1.623d1c040d6d2p-19 +
                   T * (0x1.b60ea65abe297p-32 + T * -0x1.981e4786934d4p-46))));
    s[sC2H3OH] =
        0x1.797bdf23f83eep+7 *
        (0x1.0a6e70ec26596p+3 * log(T) +
         T * (0x1.0741027d558fbp-7 +
              T * (-0x1.623d1c040d6d2p-20 +
                   T * (0x1.2409c43c7ec65p-33 + T * -0x1.981e4786934d4p-48))) -
         0x1.435417ca2120ep+4);
    h[sC2H5O] = 0x1.710a1c4780b3cp+7 *
                (T * (0x1.ac188be8002fp+2 +
                      T * (0x1.ae1a269435ceep-8 +
                           T * (-0x1.a4e57b7c88065p-20 +
                                T * (0x1.a1096bb8bb1edp-33 +
                                     T * -0x1.473cc369bef99p-47)))) -
                 0x1.289c87fcb923ap+12);
    cp[sC2H5O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.ac188be8002fp+2 +
         T * (0x1.ae1a269435ceep-7 +
              T * (-0x1.3bac1c9d6604cp-18 +
                   T * (0x1.a1096bb8bb1edp-31 + T * -0x1.990bf4442eb7fp-45))));
    s[sC2H5O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.ac188be8002fp+2 * log(T) +
         T * (0x1.ae1a269435ceep-7 +
              T * (-0x1.3bac1c9d6604cp-19 +
                   T * (0x1.160647d07cbf3p-32 + T * -0x1.990bf4442eb7fp-47))) -
         0x1.365917939a7c1p+3);
    h[sC2H5CO] = 0x1.2360aa5b5ab7bp+7 *
                 (T * (0x1.a17d005bc578ap+2 +
                       T * (0x1.f9525c864d3d7p-8 +
                            T * (-0x1.ecef71909d437p-20 +
                                 T * (0x1.e705ea04e7686p-33 +
                                      T * -0x1.7d12f1179cfadp-47)))) -
                  0x1.c1c50fba8826bp+12);
    cp[sC2H5CO] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.a17d005bc578ap+2 +
         T * (0x1.f9525c864d3d7p-7 +
              T * (-0x1.71b3952c75f29p-18 +
                   T * (0x1.e705ea04e7686p-31 + T * -0x1.dc57ad5d84398p-45))));
    s[sC2H5CO] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.a17d005bc578ap+2 * log(T) +
         T * (0x1.f9525c864d3d7p-7 +
              T * (-0x1.71b3952c75f29p-19 +
                   T * (0x1.44ae9c0344f04p-32 + T * -0x1.dc57ad5d84398p-47))) -
         0x1.4cb639cdd70dfp+2);
    h[sC2H5CHO] = 0x1.1e520994b8387p+7 *
                  (T * (0x1.dc36ffdbedd2p+2 +
                        T * (0x1.227dbff29474ap-7 +
                             T * (-0x1.1baeec2f9cb6fp-19 +
                                  T * (0x1.187cd4e50ac2cp-32 +
                                       T * -0x1.b729d7b5363cfp-47)))) -
                   0x1.9656535a85879p+14);
    cp[sC2H5CHO] =
        0x1.1e520994b8387p+7 *
        (0x1.dc36ffdbedd2p+2 +
         T * (0x1.227dbff29474ap-6 +
              T * (-0x1.a98662476b126p-18 +
                   T * (0x1.187cd4e50ac2cp-30 + T * -0x1.127a26d141e61p-44))));
    s[sC2H5CHO] =
        0x1.1e520994b8387p+7 *
        (0x1.dc36ffdbedd2p+2 * log(T) +
         T * (0x1.227dbff29474ap-6 +
              T * (-0x1.a98662476b126p-19 +
                   T * (0x1.75fbc686b903bp-32 + T * -0x1.127a26d141e61p-46))) -
         0x1.cd6ce8cc06d43p+3);
    h[sC2H] = 0x1.4c34d172a6531p+8 *
              (T * (0x1.d4d36f5349ffp+1 +
                    T * (0x1.f5571a442a209p-10 +
                         T * (-0x1.e906da2127efp-22 +
                              T * (0x1.d564813c3f32fp-35 +
                                   T * -0x1.6325fc8659684p-49)))) +
               0x1.0660610624dd3p+16);
    cp[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.d4d36f5349ffp+1 +
         T * (0x1.f5571a442a209p-9 +
              T * (-0x1.6ec52398ddf34p-20 +
                   T * (0x1.d564813c3f32fp-33 + T * -0x1.bbef7ba7efc25p-47))));
    s[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.d4d36f5349ffp+1 * log(T) +
         T * (0x1.f5571a442a209p-9 +
              T * (-0x1.6ec52398ddf34p-21 +
                   T * (0x1.38edab7d7f775p-34 + T * -0x1.bbef7ba7efc25p-49))) +
         0x1.f605fe71b579fp+1);
    h[sCH3COCH3] = 0x1.1e520994b8387p+7 *
                   (T * (0x1.d311efac1fd91p+2 +
                         T * (0x1.1fcbd96a9184ap-7 +
                              T * (-0x1.1a9ba4bb963c7p-19 +
                                   T * (0x1.187219b2ec9ap-32 +
                                        T * -0x1.b833d278a4267p-47)))) -
                    0x1.cd83921ff2e49p+14);
    cp[sCH3COCH3] =
        0x1.1e520994b8387p+7 *
        (0x1.d311efac1fd91p+2 +
         T * (0x1.1fcbd96a9184ap-6 +
              T * (-0x1.a7e97719615abp-18 +
                   T * (0x1.187219b2ec9ap-30 + T * -0x1.1320638b6698p-44))));
    s[sCH3COCH3] =
        0x1.1e520994b8387p+7 *
        (0x1.d311efac1fd91p+2 * log(T) +
         T * (0x1.1fcbd96a9184ap-6 +
              T * (-0x1.a7e97719615abp-19 +
                   T * (0x1.75ed77993b78p-32 + T * -0x1.1320638b6698p-46))) -
         0x1.984b1fb902eb7p+3);
    h[sCH3COCH2] = 0x1.2360aa5b5ab7bp+7 *
                   (T * (0x1.e2d2a60a6b3dep+2 +
                         T * (0x1.d608e55dbb92bp-8 +
                              T * (-0x1.c6e445dc83405p-20 +
                                   T * (0x1.bf0fcc81bf444p-33 +
                                        T * -0x1.5c85fa839ef9bp-47)))) -
                    0x1.d3eb90d5a5b96p+12);
    cp[sCH3COCH2] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.e2d2a60a6b3dep+2 +
         T * (0x1.d608e55dbb92bp-7 +
              T * (-0x1.552b346562704p-18 +
                   T * (0x1.bf0fcc81bf444p-31 + T * -0x1.b3a7792486b81p-45))));
    s[sCH3COCH2] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.e2d2a60a6b3dep+2 * log(T) +
         T * (0x1.d608e55dbb92bp-7 +
              T * (-0x1.552b346562704p-19 +
                   T * (0x1.2a0a88567f82dp-32 + T * -0x1.b3a7792486b81p-47))) -
         0x1.6f5616575a59fp+3);
    h[sCH2Y] = 0x1.286500f4d652p+9 *
               (T * (0x1.91483b830de6fp+1 +
                     T * (0x1.7b9398d367f2ep-10 +
                          T * (-0x1.244c014d14444p-22 +
                               T * (0x1.f37f7d2d47a3cp-36 +
                                    T * -0x1.6ec7d24a57748p-50)))) +
                0x1.8a9019ce075f7p+15);
    cp[sCH2Y] =
        0x1.286500f4d652p+9 *
        (0x1.91483b830de6fp+1 +
         T * (0x1.7b9398d367f2ep-9 +
              T * (-0x1.b67201f39e666p-21 +
                   T * (0x1.f37f7d2d47a3cp-34 + T * -0x1.ca79c6dced519p-48))));
    s[sCH2Y] =
        0x1.286500f4d652p+9 *
        (0x1.91483b830de6fp+1 * log(T) +
         T * (0x1.7b9398d367f2ep-9 +
              T * (-0x1.b67201f39e666p-22 +
                   T * (0x1.4cffa8c8da6d3p-35 + T * -0x1.ca79c6dced519p-50))) +
         0x1.03dc0e93ec868p+2);
    h[sOCH2O2H] = 0x1.07cf090c05134p+7 *
                  (T * (0x1.71463e3d5270ep+3 +
                        T * (0x1.5e273a450d07dp-9 +
                             T * (-0x1.457c28f82f5bep-21 +
                                  T * (0x1.36071cb630d7p-34 +
                                       T * -0x1.d77ff1af47eefp-49)))) -
                   0x1.06defedfa43fep+14);
    cp[sOCH2O2H] =
        0x1.07cf090c05134p+7 *
        (0x1.71463e3d5270ep+3 +
         T * (0x1.5e273a450d07dp-8 +
              T * (-0x1.e83a3d744709dp-20 +
                   T * (0x1.36071cb630d7p-32 + T * -0x1.26aff70d8cf55p-46))));
    s[sOCH2O2H] =
        0x1.07cf090c05134p+7 *
        (0x1.71463e3d5270ep+3 * log(T) +
         T * (0x1.5e273a450d07dp-8 +
              T * (-0x1.e83a3d744709dp-21 +
                   T * (0x1.9d5ed0f2ebc95p-34 + T * -0x1.26aff70d8cf55p-48))) -
         0x1.008f7d58f132ep+5);
    h[sHE] = 0x1.9c055f0ad864p+8 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sHE] =
        0x1.9c055f0ad864p+8 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sHE] = 0x1.9c055f0ad864p+8 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.db81b56eaeabcp-1);
    h[sAR] = 0x1.a0439e3ea879p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a0439e3ea879p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sAR] = 0x1.a0439e3ea879p+7 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.184c97fe63f3ap+2);
  } else if (T >= 299.999999) {
    h[sN2] = 0x1.28bba88de9adap+8 *
             (T * (0x1.c3f7fb23cd39bp+1 +
                   T * (-0x1.0355fcc6aa228p-14 +
                        T * (-0x1.680f9ef2bcab5p-23 +
                             T * (0x1.4eb4b77818ea5p-31 +
                                  T * -0x1.3d3c80b1d591dp-42)))) -
              0x1.05be7b5f1bef5p+10);
    cp[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.c3f7fb23cd39bp+1 +
         T * (-0x1.0355fcc6aa228p-13 +
              T * (-0x1.0e0bb7360d808p-21 +
                   T * (0x1.4eb4b77818ea5p-29 + T * -0x1.8c8ba0de4af64p-40))));
    s[sN2] =
        0x1.28bba88de9adap+8 *
        (0x1.c3f7fb23cd39bp+1 * log(T) +
         T * (-0x1.0355fcc6aa228p-13 +
              T * (-0x1.0e0bb7360d808p-22 +
                   T * (0x1.be4649f5768dcp-31 + T * -0x1.8c8ba0de4af64p-42))) +
         0x1.7bd611c4f96ecp+1);
    h[sCH3OCH3] = 0x1.68f6f3733777dp+7 *
                  (T * (0x1.3594cb94129c1p+1 +
                        T * (0x1.32d6c639a4804p-7 +
                             T * (-0x1.f848090065701p-22 +
                                  T * (-0x1.265f43a2a76b4p-30 +
                                       T * 0x1.341647326f48p-42)))) -
                   0x1.71f4b1a9fbe77p+14);
    cp[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.3594cb94129c1p+1 +
         T * (0x1.32d6c639a4804p-6 +
              T * (-0x1.7a3606c04c141p-20 +
                   T * (-0x1.265f43a2a76b4p-28 + T * 0x1.811bd8ff0b1ap-40))));
    s[sCH3OCH3] =
        0x1.68f6f3733777dp+7 *
        (0x1.3594cb94129c1p+1 * log(T) +
         T * (0x1.32d6c639a4804p-6 +
              T * (-0x1.7a3606c04c141p-21 +
                   T * (-0x1.887f04d8df39bp-30 + T * 0x1.811bd8ff0b1ap-42))) +
         0x1.9c68d49fffe53p+3);
    h[sCH4] = 0x1.0325882935d13p+9 *
              (T * (0x1.498b184c7d064p+2 +
                    T * (-0x1.bfaed90747e6dp-8 +
                         T * (0x1.12d748b042dd1p-16 +
                              T * (-0x1.9ff6fd7e220d5p-27 +
                                   T * 0x1.d4f26e61176cdp-39)))) -
               0x1.4034c95182a99p+13);
    cp[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.498b184c7d064p+2 +
         T * (-0x1.bfaed90747e6dp-7 +
              T * (0x1.9c42ed08644bap-15 +
                   T * (-0x1.9ff6fd7e220d5p-25 + T * 0x1.251784fcaea4p-36))));
    s[sCH4] =
        0x1.0325882935d13p+9 *
        (0x1.498b184c7d064p+2 * log(T) +
         T * (-0x1.bfaed90747e6dp-7 +
              T * (0x1.9c42ed08644bap-16 +
                   T * (-0x1.154f53a96c08ep-26 + T * 0x1.251784fcaea4p-38))) -
         0x1.28dcfe88b194ep+2);
    h[sH2] = 0x1.01c3c6b46b46cp+12 *
             (T * (0x1.2c130ac9b2911p+1 +
                   T * (0x1.058175d02a941p-8 +
                        T * (-0x1.b3b80759749dp-18 +
                             T * (0x1.5a4c582f87e3dp-28 +
                                  T * -0x1.9f3d0ec308a15p-40)))) -
              0x1.caf7b3bfb58d1p+9);
    cp[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.2c130ac9b2911p+1 +
         T * (0x1.058175d02a941p-7 +
              T * (-0x1.46ca05831775cp-16 +
                   T * (0x1.5a4c582f87e3dp-26 + T * -0x1.03862939e564dp-37))));
    s[sH2] =
        0x1.01c3c6b46b46cp+12 *
        (0x1.2c130ac9b2911p+1 * log(T) +
         T * (0x1.058175d02a941p-7 +
              T * (-0x1.46ca05831775cp-17 +
                   T * (0x1.cdbb203f5fda7p-28 + T * -0x1.03862939e564dp-39))) +
         0x1.5db38496161b4p-1);
    h[sO2] = 0x1.03d3adab9f55ap+8 *
             (T * (0x1.e42787ae5fa45p+1 +
                   T * (-0x1.88c9b66c8c0dfp-10 +
                        T * (0x1.b88f92d523df1p-19 +
                             T * (-0x1.4ca5927ada53cp-29 +
                                  T * 0x1.6d361ad59bca4p-41)))) -
              0x1.09fc63497b741p+10);
    cp[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.e42787ae5fa45p+1 +
         T * (-0x1.88c9b66c8c0dfp-9 +
              T * (0x1.4a6bae1fdae75p-17 +
                   T * (-0x1.4ca5927ada53cp-27 + T * 0x1.c883a18b02bcdp-39))));
    s[sO2] =
        0x1.03d3adab9f55ap+8 *
        (0x1.e42787ae5fa45p+1 * log(T) +
         T * (-0x1.88c9b66c8c0dfp-9 +
              T * (0x1.4a6bae1fdae75p-18 +
                   T * (-0x1.bb876df9231a5p-29 + T * 0x1.c883a18b02bcdp-41))) +
         0x1.d42eb7e3dc88dp+1);
    h[sCO] = 0x1.28d6c752d48edp+8 *
             (T * (0x1.ca2e275ab7dc8p+1 +
                   T * (-0x1.40004919b01a8p-12 +
                        T * (0x1.6bee98737258p-22 +
                             T * (0x1.f2a1ba06632ddp-33 +
                                  T * -0x1.97510b5c9478p-43)))) -
              0x1.c040b020c49bap+13);
    cp[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.ca2e275ab7dc8p+1 +
         T * (-0x1.40004919b01a8p-11 +
              T * (0x1.10f2f25695c2p-20 +
                   T * (0x1.f2a1ba06632ddp-31 + T * -0x1.fd254e33b996p-41))));
    s[sCO] =
        0x1.28d6c752d48edp+8 *
        (0x1.ca2e275ab7dc8p+1 * log(T) +
         T * (-0x1.40004919b01a8p-11 +
              T * (0x1.10f2f25695c2p-21 +
                   T * (0x1.4c6bd1599773ep-32 + T * -0x1.fd254e33b996p-43))) +
         0x1.c1138e523dba7p+1);
    h[sCO2] = 0x1.79d818115de7cp+7 *
              (T * (0x1.2dac0c62e4d1ap+1 +
                    T * (0x1.2664580d40011p-8 +
                         T * (-0x1.3ea2be239ad7bp-19 +
                              T * (0x1.51ba95ef37188p-31 +
                                   T * -0x1.01664c0a8030bp-45)))) -
               0x1.79e7f126e978dp+15);
    cp[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.2dac0c62e4d1ap+1 +
         T * (0x1.2664580d40011p-7 +
              T * (-0x1.ddf41d3568439p-18 +
                   T * (0x1.51ba95ef37188p-29 + T * -0x1.41bfdf0d203cdp-43))));
    s[sCO2] =
        0x1.79d818115de7cp+7 *
        (0x1.2dac0c62e4d1ap+1 * log(T) +
         T * (0x1.2664580d40011p-7 +
              T * (-0x1.ddf41d3568439p-19 +
                   T * (0x1.c24e1d3ef420bp-31 + T * -0x1.41bfdf0d203cdp-45))) +
         0x1.3cd43393ab431p+3);
    h[sH2O] = 0x1.cd8113bbf9d1p+8 *
              (T * (0x1.0cb67069f5672p+2 +
                    T * (-0x1.0aea4d67f214ap-10 +
                         T * (0x1.23b713c003c4fp-19 +
                              T * (-0x1.7920a18ab86b7p-30 +
                                   T * 0x1.8f03002988b76p-42)))) -
               0x1.d956e76c8b439p+14);
    cp[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.0cb67069f5672p+2 +
         T * (-0x1.0aea4d67f214ap-9 +
              T * (0x1.b5929da005a77p-18 +
                   T * (-0x1.7920a18ab86b7p-28 + T * 0x1.f2c3c033eae53p-40))));
    s[sH2O] =
        0x1.cd8113bbf9d1p+8 *
        (0x1.0cb67069f5672p+2 * log(T) +
         T * (-0x1.0aea4d67f214ap-9 +
              T * (0x1.b5929da005a77p-19 +
                   T * (-0x1.f6d62cb8f5e49p-30 + T * 0x1.f2c3c033eae53p-42))) -
         0x1.b2b14f17eb2e3p-1);
    h[sH] = 0x1.01c3c6b46b46cp+13 *
            (T * (0x1.4p+1 +
                  T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) +
             0x1.8e06a3d70a3d7p+14);
    cp[sH] =
        0x1.01c3c6b46b46cp+13 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sH] = 0x1.01c3c6b46b46cp+13 *
            (0x1.4p+1 * log(T) +
             T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) -
             0x1.c9673aa1bc7ddp-2);
    h[sO] = 0x1.03d3adab9f55ap+9 *
            (T * (0x1.9589c6bdbf12dp+1 +
                  T * (-0x1.add3ae5787a18p-10 +
                       T * (0x1.2934a583cf1ebp-19 +
                            T * (-0x1.a51e14d5c0503p-30 +
                                 T * 0x1.dbba8a635ed8p-42)))) +
             0x1.c709096bb98c8p+14);
    cp[sO] =
        0x1.03d3adab9f55ap+9 *
        (0x1.9589c6bdbf12dp+1 +
         T * (-0x1.add3ae5787a18p-9 +
              T * (0x1.bdcef845b6aep-18 +
                   T * (-0x1.a51e14d5c0503p-28 + T * 0x1.2954967e1b47p-39))));
    s[sO] = 0x1.03d3adab9f55ap+9 *
            (0x1.9589c6bdbf12dp+1 * log(T) +
             T * (-0x1.add3ae5787a18p-9 +
                  T * (0x1.bdcef845b6aep-19 + T * (-0x1.18beb88e80357p-29 +
                                                   T * 0x1.2954967e1b47p-41))) +
             0x1.06a5c1702251ep+1);
    h[sOH] = 0x1.e8db171241376p+8 *
             (T * (0x1.fef956ee7944ep+1 +
                   T * (-0x1.3ab66c9d93f86p-10 +
                        T * (0x1.9d170931f6c01p-20 +
                             T * (-0x1.0a92f54fc020dp-30 +
                                  T * 0x1.32f6d7b9a676fp-42)))) +
              0x1.a51cbf5d78812p+11);
    cp[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.fef956ee7944ep+1 +
         T * (-0x1.3ab66c9d93f86p-9 +
              T * (0x1.35d146e579101p-18 +
                   T * (-0x1.0a92f54fc020dp-28 + T * 0x1.7fb48da81014bp-40))));
    s[sOH] =
        0x1.e8db171241376p+8 *
        (0x1.fef956ee7944ep+1 * log(T) +
         T * (-0x1.3ab66c9d93f86p-9 +
              T * (0x1.35d146e579101p-19 +
                   T * (-0x1.636e9c6a55811p-30 + T * 0x1.7fb48da81014bp-42))) -
         0x1.a9fa4e98c7eb2p-4);
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
    s[sCH3] =
        0x1.148599bc79368p+9 *
        (0x1.d41e76e38c2bfp+1 * log(T) +
         T * (0x1.16bcc8dd4edc2p-9 +
              T * (0x1.6e4e659e7b1a3p-19 +
                   T * (-0x1.2f31e6b257e5dp-29 + T * 0x1.5b0478af296dap-41))) +
         0x1.ac6cd0e3b2c26p+0);
    h[sCH2O] = 0x1.14e89ec288e65p+8 *
               (T * (0x1.32cc5c0e9ea95p+2 +
                     T * (-0x1.44ad1f91a275ap-8 +
                          T * (0x1.a1708294e14d3p-17 +
                               T * (-0x1.45cdb4a6770dfp-27 +
                                    T * 0x1.72e831d82d17dp-39)))) -
                0x1.c1598ff972474p+13);
    cp[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.32cc5c0e9ea95p+2 +
         T * (-0x1.44ad1f91a275ap-7 +
              T * (0x1.391461efa8f9ep-15 +
                   T * (-0x1.45cdb4a6770dfp-25 + T * 0x1.cfa23e4e385dcp-37))));
    s[sCH2O] =
        0x1.14e89ec288e65p+8 *
        (0x1.32cc5c0e9ea95p+2 * log(T) +
         T * (-0x1.44ad1f91a275ap-7 +
              T * (0x1.391461efa8f9ep-16 +
                   T * (-0x1.b2679b889ebd4p-27 + T * 0x1.cfa23e4e385dcp-39))) +
         0x1.34a1f27267955p-1);
    h[sC2H6] = 0x1.148599bc79368p+8 *
               (T * (0x1.12a6b810273f9p+2 +
                     T * (-0x1.688cad1346485p-9 +
                          T * (0x1.4f3af3d2d5847p-16 +
                               T * (-0x1.3048b64944c4ap-26 +
                                    T * 0x1.7a244643bdb59p-38)))) -
                0x1.6811a5119ce07p+13);
    cp[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.12a6b810273f9p+2 +
         T * (-0x1.688cad1346485p-8 +
              T * (0x1.f6d86dbc4046ap-15 +
                   T * (-0x1.3048b64944c4ap-24 + T * 0x1.d8ad57d4ad22fp-36))));
    s[sC2H6] =
        0x1.148599bc79368p+8 *
        (0x1.12a6b810273f9p+2 * log(T) +
         T * (-0x1.688cad1346485p-8 +
              T * (0x1.f6d86dbc4046ap-16 +
                   T * (-0x1.95b64861b1063p-26 + T * 0x1.d8ad57d4ad22fp-38))) +
         0x1.55595f6ccd07ep+1);
    h[sCH3OCH2] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.188a983469bf1p+1 +
                        T * (0x1.37efc24f833b8p-7 +
                             T * (-0x1.30e039d2fc50bp-19 +
                                  T * (0x1.7bce31e86752bp-36 +
                                       T * 0x1.23bed466343a1p-44)))) -
                   0x1.5b5ee1c58255bp+10);
    cp[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.188a983469bf1p+1 +
         T * (0x1.37efc24f833b8p-6 +
              T * (-0x1.c95056bc7a79p-18 +
                   T * (0x1.7bce31e86752bp-34 + T * 0x1.6cae897fc1489p-42))));
    s[sCH3OCH2] =
        0x1.710a1c4780b3cp+7 *
        (0x1.188a983469bf1p+1 * log(T) +
         T * (0x1.37efc24f833b8p-6 +
              T * (-0x1.c95056bc7a79p-19 +
                   T * (0x1.fa68428b346e4p-36 + T * 0x1.6cae897fc1489p-44))) +
         0x1.03e7bbd0fc067p+4);
    h[sO2CH2OCH2O2H] = 0x1.30f32e834a014p+6 *
                       (T * (-0x1.f6163eece0fa9p-4 +
                             T * (0x1.1020f3b7ea89fp-5 +
                                  T * (-0x1.6e7e375331ae9p-16 +
                                       T * (0x1.0d3e93e205bp-27 +
                                            T * -0x1.45f8e294e4298p-40)))) -
                        0x1.0ad2125460aa6p+15);
    cp[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (-0x1.f6163eece0fa9p-4 +
         T * (0x1.1020f3b7ea89fp-4 +
              T * (-0x1.12dea97e6542fp-14 +
                   T * (0x1.0d3e93e205bp-25 + T * -0x1.97771b3a1d33dp-38))));
    s[sO2CH2OCH2O2H] =
        0x1.30f32e834a014p+6 *
        (-0x1.f6163eece0fa9p-4 * log(T) +
         T * (0x1.1020f3b7ea89fp-4 +
              T * (-0x1.12dea97e6542fp-15 +
                   T * (0x1.66fe1a82b24p-27 + T * -0x1.97771b3a1d33dp-40))) +
         0x1.174a3d7e0fd05p+5);
    h[sCH2OCH2O2H] = 0x1.af956cbfc9bdbp+6 *
                     (T * (0x1.d110ea21e1152p-1 +
                           T * (0x1.6effdfa2262afp-6 +
                                T * (-0x1.9d36a20d91085p-17 +
                                     T * (0x1.fdc8a51fdabacp-29 +
                                          T * -0x1.0585feef5cdefp-41)))) -
                      0x1.d7d9daee631f9p+13);
    cp[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.d110ea21e1152p-1 +
         T * (0x1.6effdfa2262afp-5 +
              T * (-0x1.35e8f98a2cc64p-15 +
                   T * (0x1.fdc8a51fdabacp-27 + T * -0x1.46e77eab3416ap-39))));
    s[sCH2OCH2O2H] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.d110ea21e1152p-1 * log(T) +
         T * (0x1.6effdfa2262afp-5 +
              T * (-0x1.35e8f98a2cc64p-16 +
                   T * (0x1.53db18bfe7273p-28 + T * -0x1.46e77eab3416ap-41))) +
         0x1.aa4595bbbe879p+4);
    h[sHO2] = 0x1.f7c8d6a0c5c07p+7 *
              (T * (0x1.1350a8da2956cp+2 +
                    T * (-0x1.373d07403ad69p-9 +
                         T * (0x1.d94d8de7a2103p-18 +
                              T * (-0x1.a110b2d0d403ap-28 +
                                   T * 0x1.058dbbe993b6ap-39)))) +
               0x1.0804bb6ed6777p+8);
    cp[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.1350a8da2956cp+2 +
         T * (-0x1.373d07403ad69p-8 +
              T * (0x1.62fa2a6db98c2p-16 +
                   T * (-0x1.a110b2d0d403ap-26 + T * 0x1.46f12ae3f8a44p-37))));
    s[sHO2] =
        0x1.f7c8d6a0c5c07p+7 *
        (0x1.1350a8da2956cp+2 * log(T) +
         T * (-0x1.373d07403ad69p-8 +
              T * (0x1.62fa2a6db98c2p-17 +
                   T * (-0x1.160b21e08d57cp-27 + T * 0x1.46f12ae3f8a44p-39))) +
         0x1.dbbb9643a3c3cp+1);
    h[sH2O2] = 0x1.e8db171241376p+7 *
               (T * (0x1.142b7127b57bap+2 +
                     T * (-0x1.bc46d8114239fp-12 +
                          T * (0x1.8a9c1b3e6944cp-18 +
                               T * (-0x1.859365c0a3928p-28 +
                                    T * 0x1.ffb185e179fe5p-40)))) -
                0x1.14aaf98c7e282p+14);
    cp[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.142b7127b57bap+2 +
         T * (-0x1.bc46d8114239fp-11 +
              T * (0x1.27f5146ecef39p-16 +
                   T * (-0x1.859365c0a3928p-26 + T * 0x1.3fcef3acec3efp-37))));
    s[sH2O2] =
        0x1.e8db171241376p+7 *
        (0x1.142b7127b57bap+2 * log(T) +
         T * (-0x1.bc46d8114239fp-11 +
              T * (0x1.27f5146ecef39p-17 +
                   T * (-0x1.03b7992b17b7p-27 + T * 0x1.3fcef3acec3efp-39))) +
         0x1.a309b06d709bdp+1);
    h[sCH3OCH2O2] = 0x1.af956cbfc9bdbp+6 *
                    (T * (0x1.b8942363068a9p+0 +
                          T * (0x1.35c8f46304555p-6 +
                               T * (-0x1.43fd38b7cdb03p-17 +
                                    T * (0x1.914ef37331a14p-29 +
                                         T * -0x1.b32017d7b70cp-42)))) -
                     0x1.40cb705532618p+14);
    cp[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.b8942363068a9p+0 +
         T * (0x1.35c8f46304555p-5 +
              T * (-0x1.e5fbd513b4884p-16 +
                   T * (0x1.914ef37331a14p-27 + T * -0x1.0ff40ee6d2678p-39))));
    s[sCH3OCH2O2] =
        0x1.af956cbfc9bdbp+6 *
        (0x1.b8942363068a9p+0 * log(T) +
         T * (0x1.35c8f46304555p-5 +
              T * (-0x1.e5fbd513b4884p-17 +
                   T * (0x1.0b89f7a221163p-28 + T * -0x1.0ff40ee6d2678p-41))) +
         0x1.5edfc349964cp+4);
    h[sCH3O] = 0x1.0bea1f388b9cfp+8 *
               (T * (0x1.db1c6d4903aap+1 +
                     T * (-0x1.6f9bde8a3adf3p-10 +
                          T * (0x1.a52a46122c5f3p-17 +
                               T * (-0x1.965da70bf7a12p-27 +
                                    T * 0x1.069990277f7e5p-38)))) +
                0x1.43eca57a786c2p+10);
    cp[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.db1c6d4903aap+1 +
         T * (-0x1.6f9bde8a3adf3p-9 +
              T * (0x1.3bdfb48da1476p-15 +
                   T * (-0x1.965da70bf7a12p-25 + T * 0x1.483ff4315f5dep-36))));
    s[sCH3O] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.db1c6d4903aap+1 * log(T) +
         T * (-0x1.6f9bde8a3adf3p-9 +
              T * (0x1.3bdfb48da1476p-16 +
                   T * (-0x1.0ee91a07fa6b7p-26 + T * 0x1.483ff4315f5dep-38))) +
         0x1.a4a257d92fdc9p+2);
    h[sCH3OH] = 0x1.037c7db28a41cp+8 *
                (T * (0x1.6a250944216d1p+2 +
                      T * (-0x1.0b08335c2133cp-7 +
                           T * (0x1.82f5b99b82ae8p-16 +
                                T * (-0x1.45b7fc2b1b3a3p-26 +
                                     T * 0x1.8aaaac3af75ccp-38)))) -
                 0x1.902fe4f765fd9p+14);
    cp[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.6a250944216d1p+2 +
         T * (-0x1.0b08335c2133cp-6 +
              T * (0x1.22384b34a202ep-14 +
                   T * (-0x1.45b7fc2b1b3a3p-24 + T * 0x1.ed555749b533ep-36))));
    s[sCH3OH] =
        0x1.037c7db28a41cp+8 *
        (0x1.6a250944216d1p+2 * log(T) +
         T * (-0x1.0b08335c2133cp-6 +
              T * (0x1.22384b34a202ep-15 +
                   T * (-0x1.b24aa58ecef84p-26 + T * 0x1.ed555749b533ep-38))) -
         0x1.cb6ee783204a4p-1);
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
    s[sCH3OCH2O] =
        0x1.10565da70e7acp+7 *
        (0x1.a1236b2999ac6p+1 * log(T) +
         T * (0x1.6bf6efabb70dp-6 +
              T * (-0x1.053d7af9638p-18 +
                   T * (-0x1.6204ff3b4eee3p-34 + T * 0x1.fccf7c342220dp-44))) +
         0x1.8bc6b66806799p+3);
    h[sOCH2OCHO] = 0x1.bb2d881f811c1p+6 *
                   (T * (0x1.f6c7ef82da3a7p+1 +
                         T * (0x1.304d51e4c68b1p-7 +
                              T * (-0x1.28d2aad5835dfp-23 +
                                   T * (-0x1.dcf46930e5c2dp-30 +
                                        T * 0x1.058f261d3e2e7p-41)))) -
                    0x1.3e062b367a0f9p+15);
    cp[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.f6c7ef82da3a7p+1 +
         T * (0x1.304d51e4c68b1p-6 +
              T * (-0x1.bd3c0040450cep-22 +
                   T * (-0x1.dcf46930e5c2dp-28 + T * 0x1.46f2efa48dbap-39))));
    s[sOCH2OCHO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.f6c7ef82da3a7p+1 * log(T) +
         T * (0x1.304d51e4c68b1p-6 +
              T * (-0x1.bd3c0040450cep-23 +
                   T * (-0x1.3df84620992c9p-29 + T * 0x1.46f2efa48dbap-41))) +
         0x1.468e3cd9a5226p+3);
    h[sHO2CH2OCHO] = 0x1.694b465e7f561p+6 *
                     (T * (0x1.07e6d56a32d75p-1 +
                           T * (0x1.86019a3c8a528p-6 +
                                T * (-0x1.96fd95f9a8e7fp-17 +
                                     T * (0x1.aad58b39fd6efp-29 +
                                          T * -0x1.5399dbe8fe478p-42)))) -
                      0x1.c7420af4f0d84p+15);
    cp[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.07e6d56a32d75p-1 +
         T * (0x1.86019a3c8a528p-5 +
              T * (-0x1.313e307b3eadfp-15 +
                   T * (0x1.aad58b39fd6efp-27 + T * -0x1.a88052e33dd96p-40))));
    s[sHO2CH2OCHO] =
        0x1.694b465e7f561p+6 *
        (0x1.07e6d56a32d75p-1 * log(T) +
         T * (0x1.86019a3c8a528p-5 +
              T * (-0x1.313e307b3eadfp-16 +
                   T * (0x1.1c8e5cd1539f5p-28 + T * -0x1.a88052e33dd96p-42))) +
         0x1.c330351c417f2p+4);
    h[sHCO] = 0x1.1e87151786dcap+8 *
              (T * (0x1.0f33f48eb2b5dp+2 +
                    T * (-0x1.b341f75f5d33cp-10 +
                         T * (0x1.393e0dcd17733p-18 +
                              T * (-0x1.cd3ebf92b0d7p-29 +
                                   T * 0x1.ec7ca35ee47e2p-41)))) +
               0x1.e40d2de00d1b7p+11);
    cp[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.0f33f48eb2b5dp+2 +
         T * (-0x1.b341f75f5d33cp-9 +
              T * (0x1.d5dd14b3a32cdp-17 +
                   T * (-0x1.cd3ebf92b0d7p-27 + T * 0x1.33cde61b4ecedp-38))));
    s[sHCO] =
        0x1.1e87151786dcap+8 *
        (0x1.0f33f48eb2b5dp+2 * log(T) +
         T * (-0x1.b341f75f5d33cp-9 +
              T * (0x1.d5dd14b3a32cdp-18 +
                   T * (-0x1.337f2a61cb3ap-28 + T * 0x1.33cde61b4ecedp-40))) +
         0x1.a777f849a83fap+1);
    h[sO2CHO] = 0x1.10865abf22dacp+7 *
                (T * (0x1.faf4b6e128239p+1 +
                      T * (0x1.5b5928149012bp-8 +
                           T * (-0x1.d6667cec396ccp-20 +
                                T * (0x1.1798c85f5816ap-32 +
                                     T * -0x1.9e50345189b2fp-48)))) -
                 0x1.0f3fc0d1b7176p+14);
    cp[sO2CHO] =
        0x1.10865abf22dacp+7 *
        (0x1.faf4b6e128239p+1 +
         T * (0x1.5b5928149012bp-7 +
              T * (-0x1.60ccddb12b119p-18 +
                   T * (0x1.1798c85f5816ap-30 + T * -0x1.02f220b2f60fdp-45))));
    s[sO2CHO] =
        0x1.10865abf22dacp+7 *
        (0x1.faf4b6e128239p+1 * log(T) +
         T * (0x1.5b5928149012bp-7 +
              T * (-0x1.60ccddb12b119p-19 +
                   T * (0x1.74cbb5d475738p-32 + T * -0x1.02f220b2f60fdp-47))) +
         0x1.78fbe3dbdd0b1p+3);
    h[sHOCHO] = 0x1.694b465e7f561p+7 *
                (T * (0x1.f2fd834dfdb9dp+1 +
                      T * (-0x1.d274d321e1b85p-10 +
                           T * (0x1.8d4a5db43331fp-17 +
                                T * (-0x1.78ab146067505p-27 +
                                     T * 0x1.e18a77a43a465p-39)))) -
                 0x1.6d6537ced9168p+15);
    cp[sHOCHO] =
        0x1.694b465e7f561p+7 *
        (0x1.f2fd834dfdb9dp+1 +
         T * (-0x1.d274d321e1b85p-9 +
              T * (0x1.29f7c64726657p-15 +
                   T * (-0x1.78ab146067505p-25 + T * 0x1.2cf68ac6a46bfp-36))));
    s[sHOCHO] =
        0x1.694b465e7f561p+7 *
        (0x1.f2fd834dfdb9dp+1 * log(T) +
         T * (-0x1.d274d321e1b85p-9 +
              T * (0x1.29f7c64726657p-16 +
                   T * (-0x1.f639708089c07p-27 + T * 0x1.2cf68ac6a46bfp-38))) +
         0x1.d65edbc309d57p+2);
    h[sHOCH2O] = 0x1.618d10efa02bbp+7 *
                 (T * (0x1.07283f191a834p+2 +
                       T * (0x1.ee0b28e595de3p-9 +
                            T * (0x1.51a2c161bff7fp-20 +
                                 T * (-0x1.72393312daf39p-30 +
                                      T * 0x1.47e5d8012317p-42)))) -
                  0x1.6449334d6a162p+14);
    cp[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.07283f191a834p+2 +
         T * (0x1.ee0b28e595de3p-8 +
              T * (0x1.fa7422129ff3ep-19 +
                   T * (-0x1.72393312daf39p-28 + T * 0x1.99df4e016bdccp-40))));
    s[sHOCH2O] =
        0x1.618d10efa02bbp+7 *
        (0x1.07283f191a834p+2 * log(T) +
         T * (0x1.ee0b28e595de3p-8 +
              T * (0x1.fa7422129ff3ep-20 +
                   T * (-0x1.eda1996e7944cp-30 + T * 0x1.99df4e016bdccp-42))) +
         0x1.ddf4e686dd296p+2);
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
    s[sHOCH2OCO] =
        0x1.bb2d881f811c1p+6 *
        (0x1.853c57a9e00dcp+2 * log(T) +
         T * (0x1.a5f2ba8881b83p-7 +
              T * (0x1.125dfa29a9147p-20 +
                   T * (-0x1.1787b9d5157d7p-29 + T * 0x1.fa265e4f865bap-42))) +
         0x1.45308fd54a9fbp+1);
    h[sCH2OH] = 0x1.0bea1f388b9cfp+8 *
                (T * (0x1.1e9d2ec4b982fp+2 +
                      T * (-0x1.621426907e0b2p-11 +
                           T * (0x1.377aece9515d9p-17 +
                                T * (-0x1.396b8e087efebp-27 +
                                     T * 0x1.a0528cc47799p-39)))) -
                 0x1.b597532617c1cp+11);
    cp[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.1e9d2ec4b982fp+2 +
         T * (-0x1.621426907e0b2p-10 +
              T * (0x1.d338635dfa0c5p-16 +
                   T * (-0x1.396b8e087efebp-25 + T * 0x1.043397facabfap-36))));
    s[sCH2OH] =
        0x1.0bea1f388b9cfp+8 *
        (0x1.1e9d2ec4b982fp+2 * log(T) +
         T * (-0x1.621426907e0b2p-10 +
              T * (0x1.d338635dfa0c5p-17 +
                   T * (-0x1.a1e4bd60a9539p-27 + T * 0x1.043397facabfap-38))) +
         0x1.a791bc5586445p+1);
    h[sCH3O2H] = 0x1.5a2209942d73ap+7 *
                 (T * (0x1.73e470edd3d88p+1 +
                       T * (0x1.1eb61c984d4f8p-7 +
                            T * (0x1.d8aa160d6df13p-20 +
                                 T * (-0x1.b25a9c89775aap-28 +
                                      T * 0x1.7a36802f1e85cp-39)))) -
                  0x1.07e5da5119cep+14);
    cp[sCH3O2H] =
        0x1.5a2209942d73ap+7 *
        (0x1.73e470edd3d88p+1 +
         T * (0x1.1eb61c984d4f8p-6 +
              T * (0x1.627f908a1274ep-18 +
                   T * (-0x1.b25a9c89775aap-26 + T * 0x1.d8c4203ae6273p-37))));
    s[sCH3O2H] =
        0x1.5a2209942d73ap+7 *
        (0x1.73e470edd3d88p+1 * log(T) +
         T * (0x1.1eb61c984d4f8p-6 +
              T * (0x1.627f908a1274ep-19 +
                   T * (-0x1.2191bdb0fa3c7p-27 + T * 0x1.d8c4203ae6273p-39))) +
         0x1.6bf96f8d56954p+3);
    h[sCH3O2] = 0x1.618d10efa02bbp+7 *
                (T * (0x1.f93038acefb7bp+0 +
                      T * (0x1.f720a6714bfb3p-8 +
                           T * (-0x1.1d213e380ac7fp-19 +
                                T * (0x1.5fc472c5e98c6p-34 +
                                     T * 0x1.fc5af9258a787p-45)))) +
                 0x1.fc8ec3760bf5dp+7);
    cp[sCH3O2] =
        0x1.618d10efa02bbp+7 *
        (0x1.f93038acefb7bp+0 +
         T * (0x1.f720a6714bfb3p-7 +
              T * (-0x1.abb1dd54102bfp-18 +
                   T * (0x1.5fc472c5e98c6p-32 + T * 0x1.3db8dbb7768b4p-42))));
    s[sCH3O2] =
        0x1.618d10efa02bbp+7 *
        (0x1.f93038acefb7bp+0 * log(T) +
         T * (0x1.f720a6714bfb3p-7 +
              T * (-0x1.abb1dd54102bfp-19 +
                   T * (0x1.d505ee5d3765dp-34 + T * 0x1.3db8dbb7768b4p-44))) +
         0x1.0eb5f3519bd4p+4);
    h[sHO2CHO] = 0x1.0c18906567f93p+7 *
                 (T * (0x1.365ad767049bfp+1 +
                       T * (0x1.67f788ebaa104p-7 +
                            T * (-0x1.79634e9981ff4p-18 +
                                 T * (0x1.adeadc701196cp-30 +
                                      T * -0x1.9a919c152b457p-43)))) -
                  0x1.153599e83e426p+15);
    cp[sHO2CHO] =
        0x1.0c18906567f93p+7 *
        (0x1.365ad767049bfp+1 +
         T * (0x1.67f788ebaa104p-6 +
              T * (-0x1.1b0a7af3217f7p-16 +
                   T * (0x1.adeadc701196cp-28 + T * -0x1.009b018d3b0b6p-40))));
    s[sHO2CHO] =
        0x1.0c18906567f93p+7 *
        (0x1.365ad767049bfp+1 * log(T) +
         T * (0x1.67f788ebaa104p-6 +
              T * (-0x1.1b0a7af3217f7p-17 +
                   T * (0x1.1e9c92f5610f3p-29 + T * -0x1.009b018d3b0b6p-42))) +
         0x1.180b629f31891p+4);
    h[sCH3OCH2O2H] = 0x1.aa02db5960032p+6 *
                     (T * (0x1.798870aff401ap-1 +
                           T * (0x1.7ee6fb7577688p-6 +
                                T * (-0x1.a1406417d4bdp-17 +
                                     T * (0x1.02f617093dcfbp-28 +
                                          T * -0x1.12795f0823b3bp-41)))) -
                      0x1.239c0eb1c432dp+15);
    cp[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.798870aff401ap-1 +
         T * (0x1.7ee6fb7577688p-5 +
              T * (-0x1.38f04b11df8dcp-15 +
                   T * (0x1.02f617093dcfbp-26 + T * -0x1.5717b6ca2ca0ap-39))));
    s[sCH3OCH2O2H] =
        0x1.aa02db5960032p+6 *
        (0x1.798870aff401ap-1 * log(T) +
         T * (0x1.7ee6fb7577688p-5 +
              T * (-0x1.38f04b11df8dcp-16 +
                   T * (0x1.59481eb6fd14fp-28 + T * -0x1.5717b6ca2ca0ap-41))) +
         0x1.9cb14bb78e772p+4);
    h[sCH2] = 0x1.286500f4d652p+9 *
              (T * (0x1.dbd99c6901cc3p+1 +
                    T * (0x1.4df2d3d1e1a52p-11 +
                         T * (0x1.84f5432c1736bp-21 +
                              T * (-0x1.df77ac364c39ap-31 +
                                   T * 0x1.74046d5bc53a4p-42)))) +
               0x1.6660c5f06f694p+15);
    cp[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.dbd99c6901cc3p+1 +
         T * (0x1.4df2d3d1e1a52p-10 +
              T * (0x1.23b7f2611169p-19 +
                   T * (-0x1.df77ac364c39ap-29 + T * 0x1.d10588b2b688dp-40))));
    s[sCH2] =
        0x1.286500f4d652p+9 *
        (0x1.dbd99c6901cc3p+1 * log(T) +
         T * (0x1.4df2d3d1e1a52p-10 +
              T * (0x1.23b7f2611169p-20 +
                   T * (-0x1.3fa51d7988267p-30 + T * 0x1.d10588b2b688dp-42))) +
         0x1.c0c342e04f609p+0);
    h[sCH] = 0x1.3f584141b0716p+9 *
             (T * (0x1.beb06664b8e74p+1 +
                   T * (0x1.54136aa56c34p-13 +
                        T * (-0x1.2e6ecd946abcep-21 +
                             T * (0x1.b2b29c1af28cap-31 +
                                  T * -0x1.3ca4c67e12c19p-42)))) +
              0x1.13d4a56041893p+16);
    cp[sCH] =
        0x1.3f584141b0716p+9 *
        (0x1.beb06664b8e74p+1 +
         T * (0x1.54136aa56c34p-12 +
              T * (-0x1.c5a6345ea01b5p-20 +
                   T * (0x1.b2b29c1af28cap-29 + T * -0x1.8bcdf81d9771fp-40))));
    s[sCH] =
        0x1.3f584141b0716p+9 *
        (0x1.beb06664b8e74p+1 * log(T) +
         T * (0x1.54136aa56c34p-12 +
              T * (-0x1.c5a6345ea01b5p-21 +
                   T * (0x1.21cc6811f7087p-30 + T * -0x1.8bcdf81d9771fp-42))) +
         0x1.0ac9d24689515p+1);
    h[sCH3OCHO] = 0x1.14e89ec288e65p+7 *
                  (T * (0x1.7decabe54b9ecp+2 +
                        T * (-0x1.33644f1b2c678p-8 +
                             T * (0x1.8bbee4c3259c8p-16 +
                                  T * (-0x1.6474036d37483p-26 +
                                       T * 0x1.b93e890e42f8dp-38)))) -
                   0x1.6406a74538ef3p+15);
    cp[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.7decabe54b9ecp+2 +
         T * (-0x1.33644f1b2c678p-7 +
              T * (0x1.28cf2b925c356p-14 +
                   T * (-0x1.6474036d37483p-24 + T * 0x1.13c715a8e9db8p-35))));
    s[sCH3OCHO] =
        0x1.14e89ec288e65p+7 *
        (0x1.7decabe54b9ecp+2 * log(T) +
         T * (-0x1.33644f1b2c678p-7 +
              T * (0x1.28cf2b925c356p-15 +
                   T * (-0x1.db4559e6f4604p-26 + T * 0x1.13c715a8e9db8p-37))) +
         0x1.802cb5da5bc56p-1);
    h[sCH2OCHO] = 0x1.19a2d4e8daa2p+7 *
                  (T * (0x1.27b8753c6d18bp+1 +
                        T * (0x1.27b04f2c7ccf8p-7 +
                             T * (-0x1.e5e73fcac7a2p-21 +
                                  T * (-0x1.3cbda76804f0cp-30 +
                                       T * 0x1.7ee3b3aff0369p-42)))) -
                   0x1.3d0c59e83e426p+14);
    cp[sCH2OCHO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.27b8753c6d18bp+1 +
         T * (0x1.27b04f2c7ccf8p-6 +
              T * (-0x1.6c6d6fd815b98p-19 +
                   T * (-0x1.3cbda76804f0cp-28 + T * 0x1.de9ca09bec443p-40))));
    s[sCH2OCHO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.27b8753c6d18bp+1 * log(T) +
         T * (0x1.27b04f2c7ccf8p-6 +
              T * (-0x1.6c6d6fd815b98p-20 +
                   T * (-0x1.a652348ab141p-30 + T * 0x1.de9ca09bec443p-42))) +
         0x1.127ac4212ca07p+4);
    h[sC2H4O1X2] = 0x1.797bdf23f83eep+7 *
                   (T * (0x1.e128a7bef64a7p+1 +
                         T * (-0x1.355ead45af1cfp-8 +
                              T * (0x1.c11fe0055f15dp-16 +
                                   T * (-0x1.b0f76facdd394p-26 +
                                        T * 0x1.19c18262adcp-37)))) -
                    0x1.d88d075f6fd22p+12);
    cp[sC2H4O1X2] =
        0x1.797bdf23f83eep+7 *
        (0x1.e128a7bef64a7p+1 +
         T * (-0x1.355ead45af1cfp-7 +
              T * (0x1.50d7e80407506p-14 +
                   T * (-0x1.b0f76facdd394p-24 + T * 0x1.6031e2fb593p-35))));
    s[sC2H4O1X2] =
        0x1.797bdf23f83eep+7 *
        (0x1.e128a7bef64a7p+1 * log(T) +
         T * (-0x1.355ead45af1cfp-7 +
              T * (0x1.50d7e80407506p-15 +
                   T * (-0x1.20a4f51de8d0dp-25 + T * 0x1.6031e2fb593p-37))) +
         0x1.f662435696e59p+2);
    h[sC2H4] = 0x1.286500f4d652p+8 *
               (T * (0x1.fac7161413885p+1 +
                     T * (-0x1.f02424e25f739p-9 +
                          T * (0x1.3f521c8e0e71ap-16 +
                               T * (-0x1.2908f186fa417p-26 +
                                    T * 0x1.7bd406e98fe6dp-38)))) +
                0x1.3e1c6a6a0125ap+12);
    cp[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.fac7161413885p+1 +
         T * (-0x1.f02424e25f739p-8 +
              T * (0x1.defb2ad515aa7p-15 +
                   T * (-0x1.2908f186fa417p-24 + T * 0x1.dac908a3f3e08p-36))));
    s[sC2H4] =
        0x1.286500f4d652p+8 *
        (0x1.fac7161413885p+1 * log(T) +
         T * (-0x1.f02424e25f739p-8 +
              T * (0x1.defb2ad515aa7p-16 +
                   T * (-0x1.8c0becb3f8574p-26 + T * 0x1.dac908a3f3e08p-38))) +
         0x1.063a32b68b97dp+2);
    h[sC2H3] = 0x1.337122eb1b451p+8 *
               (T * (0x1.ae9039c97a5cfp+1 +
                     T * (0x1.16acf159a70dp-13 +
                          T * (0x1.38c01e5c6a1b8p-17 +
                               T * (-0x1.4064b2299de8ep-27 +
                                    T * 0x1.aab03bcebf425p-39)))) +
                0x1.0d55eaf4f0d84p+15);
    cp[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.ae9039c97a5cfp+1 +
         T * (0x1.16acf159a70dp-12 +
              T * (0x1.d5202d8a9f294p-16 +
                   T * (-0x1.4064b2299de8ep-25 + T * 0x1.0aae256137897p-36))));
    s[sC2H3] =
        0x1.337122eb1b451p+8 *
        (0x1.ae9039c97a5cfp+1 * log(T) +
         T * (0x1.16acf159a70dp-12 +
              T * (0x1.d5202d8a9f294p-17 +
                   T * (-0x1.ab30ed8cd28bdp-27 + T * 0x1.0aae256137897p-38))) +
         0x1.fa910372fc21p+2);
    h[sC2H5] = 0x1.1e1d11b23ceecp+8 *
               (T * (0x1.53ca13340f513p+0 +
                     T * (0x1.216f366b3760bp-7 +
                          T * (-0x1.131d0d5a72c85p-19 +
                               T * (-0x1.4b1c599ec5983p-34 +
                                    T * 0x1.8b125f7e3e24ep-44)))) +
                0x1.a3a338ef34d6ap+13);
    cp[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.53ca13340f513p+0 +
         T * (0x1.216f366b3760bp-6 +
              T * (-0x1.9cab9407ac2c8p-18 +
                   T * (-0x1.4b1c599ec5983p-32 + T * 0x1.edd6f75dcdae1p-42))));
    s[sC2H5] =
        0x1.1e1d11b23ceecp+8 *
        (0x1.53ca13340f513p+0 * log(T) +
         T * (0x1.216f366b3760bp-6 +
              T * (-0x1.9cab9407ac2c8p-19 +
                   T * (-0x1.b97b222907759p-34 + T * 0x1.edd6f75dcdae1p-44))) +
         0x1.12dcdce548c49p+4);
    h[sOHY] = 0x1.e8db171241376p+8 *
              (T * (0x1.d191eeaa6d267p+1 +
                    T * (0x1.8429f9564b43dp-14 +
                         T * (-0x1.2bf61a52a5d29p-21 +
                              T * (0x1.4818341347ec4p-31 +
                                   T * -0x1.7bb7e8707d504p-43)))) +
               0x1.86ca99999999ap+15);
    cp[sOHY] =
        0x1.e8db171241376p+8 *
        (0x1.d191eeaa6d267p+1 +
         T * (0x1.8429f9564b43dp-13 +
              T * (-0x1.c1f1277bf8bbep-20 +
                   T * (0x1.4818341347ec4p-29 + T * -0x1.daa5e28c9ca44p-41))));
    s[sOHY] =
        0x1.e8db171241376p+8 *
        (0x1.d191eeaa6d267p+1 * log(T) +
         T * (0x1.8429f9564b43dp-13 +
              T * (-0x1.c1f1277bf8bbep-21 +
                   T * (0x1.b5759ac45fe5bp-31 + T * -0x1.daa5e28c9ca44p-43))) +
         0x1.5bde481f53826p+0);
    h[sCH3OCO] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.305c570371fa7p+2 +
                       T * (0x1.ffc7d94d84b24p-9 +
                            T * (0x1.6aff97daa7135p-18 +
                                 T * (-0x1.9e65a042e2d3ap-28 +
                                      T * 0x1.0954b2f820328p-39)))) -
                  0x1.502efb7e90ff9p+14);
    cp[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.305c570371fa7p+2 +
         T * (0x1.ffc7d94d84b24p-8 +
              T * (0x1.103fb1e3fd4e8p-16 +
                   T * (-0x1.9e65a042e2d3ap-26 + T * 0x1.4ba9dfb6283f2p-37))));
    s[sCH3OCO] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.305c570371fa7p+2 * log(T) +
         T * (0x1.ffc7d94d84b24p-8 +
              T * (0x1.103fb1e3fd4e8p-17 +
                   T * (-0x1.1443c02c9737cp-27 + T * 0x1.4ba9dfb6283f2p-39))) +
         0x1.31fb543ef1826p+2);
    h[sHOCH2O2H] = 0x1.03a80e5d6edeap+7 *
                   (T * (0x1.db6f4abb900e8p+0 +
                         T * (0x1.08ba1fab53f1cp-6 +
                              T * (-0x1.2de90dc2eea47p-17 +
                                   T * (0x1.7fc781d37a6fep-29 +
                                        T * -0x1.98372b9f1305ep-42)))) -
                    0x1.38bee4ea4a8c1p+15);
    cp[sHOCH2O2H] =
        0x1.03a80e5d6edeap+7 *
        (0x1.db6f4abb900e8p+0 +
         T * (0x1.08ba1fab53f1cp-5 +
              T * (-0x1.c4dd94a465f6ap-16 +
                   T * (0x1.7fc781d37a6fep-27 + T * -0x1.fe44f686d7c75p-40))));
    s[sHOCH2O2H] =
        0x1.03a80e5d6edeap+7 *
        (0x1.db6f4abb900e8p+0 * log(T) +
         T * (0x1.08ba1fab53f1cp-5 +
              T * (-0x1.c4dd94a465f6ap-17 +
                   T * (0x1.ffb4ad19f8953p-29 + T * -0x1.fe44f686d7c75p-42))) +
         0x1.3177e6dc428b9p+4);
    h[sHOCH2O2] = 0x1.07cf090c05134p+7 *
                  (T * (0x1.6d5d82a78fa68p+1 +
                        T * (0x1.7ed59706fe5d2p-7 +
                             T * (-0x1.a4ceefc88db8bp-18 +
                                  T * (0x1.11bf4f5fd78a7p-29 +
                                       T * -0x1.330679572aa97p-42)))) -
                   0x1.672a7a786c227p+14);
    cp[sHOCH2O2] =
        0x1.07cf090c05134p+7 *
        (0x1.6d5d82a78fa68p+1 +
         T * (0x1.7ed59706fe5d2p-6 +
              T * (-0x1.3b9b33d66a4a8p-16 +
                   T * (0x1.11bf4f5fd78a7p-27 + T * -0x1.7fc817acf553dp-40))));
    s[sHOCH2O2] =
        0x1.07cf090c05134p+7 *
        (0x1.6d5d82a78fa68p+1 * log(T) +
         T * (0x1.7ed59706fe5d2p-6 +
              T * (-0x1.3b9b33d66a4a8p-17 +
                   T * (0x1.6cff147fca0dfp-29 + T * -0x1.7fc817acf553dp-42))) +
         0x1.e589adc8fb86fp+3);
    h[sC2H3CO] = 0x1.2e0c24d5d505ep+7 *
                 (T * (0x1.5cc790cdc316fp+0 +
                       T * (0x1.0245be7cdd334p-6 +
                            T * (-0x1.4fca08a925f03p-17 +
                                 T * (0x1.fd192b8503ba4p-29 +
                                      T * -0x1.443a24752bb71p-41)))) +
                  0x1.0a1b3c01a36e3p+12);
    cp[sC2H3CO] =
        0x1.2e0c24d5d505ep+7 *
        (0x1.5cc790cdc316fp+0 +
         T * (0x1.0245be7cdd334p-5 +
              T * (-0x1.f7af0cfdb8e85p-16 +
                   T * (0x1.fd192b8503ba4p-27 + T * -0x1.9548ad9276a4dp-39))));
    s[sC2H3CO] =
        0x1.2e0c24d5d505ep+7 *
        (0x1.5cc790cdc316fp+0 * log(T) +
         T * (0x1.0245be7cdd334p-5 +
              T * (-0x1.f7af0cfdb8e85p-17 +
                   T * (0x1.53661d0357d18p-28 + T * -0x1.9548ad9276a4dp-41))) +
         0x1.1433d54f524dbp+4);
    h[sC2H3CHO] = 0x1.289dd93acb778p+7 *
                  (T * (0x1.2b5f26ce6d584p-2 +
                        T * (0x1.2242963a85054p-6 +
                             T * (-0x1.49e1752def62bp-17 +
                                  T * (0x1.b8260f891ec45p-29 +
                                       T * -0x1.fd3b33df575aap-42)))) -
                   0x1.6c214467381d8p+13);
    cp[sC2H3CHO] =
        0x1.289dd93acb778p+7 *
        (0x1.2b5f26ce6d584p-2 +
         T * (0x1.2242963a85054p-5 +
              T * (-0x1.eed22fc4e714p-16 +
                   T * (0x1.b8260f891ec45p-27 + T * -0x1.3e45006b9698ap-39))));
    s[sC2H3CHO] =
        0x1.289dd93acb778p+7 *
        (0x1.2b5f26ce6d584p-2 * log(T) +
         T * (0x1.2242963a85054p-5 +
              T * (-0x1.eed22fc4e714p-17 +
                   T * (0x1.256eb5061482ep-28 + T * -0x1.3e45006b9698ap-41))) +
         0x1.6e348b220791cp+4);
    h[sCH3CO] = 0x1.8252e1709b12p+7 *
                (T * (0x1.024bb3c81908ep+2 +
                      T * (0x1.cbf486346d062p-12 +
                           T * (0x1.577c5363b078cp-17 +
                                T * (-0x1.51224ee37e8f9p-27 +
                                     T * 0x1.ae918ec1d842fp-39)))) -
                 0x1.4f425c91d14e4p+11);
    cp[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.024bb3c81908ep+2 +
         T * (0x1.cbf486346d062p-11 +
              T * (0x1.019d3e8ac45a9p-15 +
                   T * (-0x1.51224ee37e8f9p-25 + T * 0x1.0d1af9392729dp-36))));
    s[sCH3CO] =
        0x1.8252e1709b12p+7 *
        (0x1.024bb3c81908ep+2 * log(T) +
         T * (0x1.cbf486346d062p-11 +
              T * (0x1.019d3e8ac45a9p-16 +
                   T * (-0x1.c18313d9fe14cp-27 + T * 0x1.0d1af9392729dp-38))) +
         0x1.f72735ceeee6p+2);
    h[sC2H5O2H] = 0x1.0bea1f388b9cfp+7 *
                  (T * (0x1.92c3240191fb7p+0 +
                        T * (0x1.20ab70fb78836p-6 +
                             T * (-0x1.1b3430f783c1bp-17 +
                                  T * (0x1.48c13c5f35171p-29 +
                                       T * -0x1.4da4ad2c702f6p-42)))) -
                   0x1.505f58e219653p+14);
    cp[sC2H5O2H] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.92c3240191fb7p+0 +
         T * (0x1.20ab70fb78836p-5 +
              T * (-0x1.a8ce497345a28p-16 +
                   T * (0x1.48c13c5f35171p-27 + T * -0x1.a10dd8778c3b3p-40))));
    s[sC2H5O2H] =
        0x1.0bea1f388b9cfp+7 *
        (0x1.92c3240191fb7p+0 * log(T) +
         T * (0x1.20ab70fb78836p-5 +
              T * (-0x1.a8ce497345a28p-17 +
                   T * (0x1.b656fb299c1ecp-29 + T * -0x1.a10dd8778c3b3p-42))) +
         0x1.30c158248443cp+4);
    h[sC2H5O2] = 0x1.10565da70e7acp+7 *
                 (T * (0x1.201046138aa16p+2 +
                       T * (0x1.c2dd6e725be29p-9 +
                            T * (0x1.0929184fd2f82p-16 +
                                 T * (-0x1.2955c64b81442p-26 +
                                      T * 0x1.947913f2f48f4p-38)))) -
                  0x1.5137aa6f3f53p+12);
    cp[sC2H5O2] =
        0x1.10565da70e7acp+7 *
        (0x1.201046138aa16p+2 +
         T * (0x1.c2dd6e725be29p-8 +
              T * (0x1.8dbda477bc743p-15 +
                   T * (-0x1.2955c64b81442p-24 + T * 0x1.f99758efb1b31p-36))));
    s[sC2H5O2] =
        0x1.10565da70e7acp+7 *
        (0x1.201046138aa16p+2 * log(T) +
         T * (0x1.c2dd6e725be29p-8 +
              T * (0x1.8dbda477bc743p-16 +
                   T * (-0x1.8c725dba01b03p-26 + T * 0x1.f99758efb1b31p-38))) +
         0x1.fa8dbb94ec0adp+2);
    h[sCH2CO] = 0x1.8b966bc75924bp+7 *
                (T * (0x1.d0710e8b08315p+0 +
                      T * (0x1.460e40a44ec32p-7 +
                           T * (-0x1.ef4c9a0cb058fp-18 +
                                T * (0x1.f2506fa81c7d3p-29 +
                                     T * -0x1.c1187cc2d5855p-41)))) -
                 0x1.b8df302b40f67p+12);
    cp[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.d0710e8b08315p+0 +
         T * (0x1.460e40a44ec32p-6 +
              T * (-0x1.737973898442bp-16 +
                   T * (0x1.f2506fa81c7d3p-27 + T * -0x1.18af4df9c5735p-38))));
    s[sCH2CO] =
        0x1.8b966bc75924bp+7 *
        (0x1.d0710e8b08315p+0 * log(T) +
         T * (0x1.460e40a44ec32p-6 +
              T * (-0x1.737973898442bp-17 +
                   T * (0x1.4c359fc568537p-28 + T * -0x1.18af4df9c5735p-40))) +
         0x1.b37435fd120efp+3);
    h[sCH2CHO] = 0x1.8252e1709b12p+7 *
                 (T * (0x1.65c36976bc1fp+1 +
                       T * (0x1.4b48624b4cfb2p-8 +
                            T * (0x1.69d47f549311p-18 +
                                 T * (-0x1.0a8c621a386ffp-27 +
                                      T * 0x1.887a53b524be6p-39)))) +
                  0x1.45e3d3c361134p+7);
    cp[sCH2CHO] =
        0x1.8252e1709b12p+7 *
        (0x1.65c36976bc1fp+1 +
         T * (0x1.4b48624b4cfb2p-7 +
              T * (0x1.0f5f5f7f6e4ccp-16 +
                   T * (-0x1.0a8c621a386ffp-25 + T * 0x1.ea98e8a26dedfp-37))));
    s[sCH2CHO] =
        0x1.8252e1709b12p+7 *
        (0x1.65c36976bc1fp+1 * log(T) +
         T * (0x1.4b48624b4cfb2p-7 +
              T * (0x1.0f5f5f7f6e4ccp-17 +
                   T * (-0x1.6365d822f5ea9p-27 + T * 0x1.ea98e8a26dedfp-39))) +
         0x1.8bab5766ef226p+3);
    h[sCH3CO3H] = 0x1.b54dcef92163bp+6 *
                  (T * (0x1.1ee4d80666cc6p+1 +
                        T * (0x1.14dc15ff9684bp-6 +
                             T * (-0x1.1bf7e53d89d96p-17 +
                                  T * (0x1.4c758d68c8dedp-29 +
                                       T * -0x1.501e14c6bf9b8p-42)))) -
                   0x1.4bc790f27bb3p+15);
    cp[sCH3CO3H] =
        0x1.b54dcef92163bp+6 *
        (0x1.1ee4d80666cc6p+1 +
         T * (0x1.14dc15ff9684bp-5 +
              T * (-0x1.a9f3d7dc4ec61p-16 +
                   T * (0x1.4c758d68c8dedp-27 + T * -0x1.a42599f86f826p-40))));
    s[sCH3CO3H] =
        0x1.b54dcef92163bp+6 *
        (0x1.1ee4d80666cc6p+1 * log(T) +
         T * (0x1.14dc15ff9684bp-5 +
              T * (-0x1.a9f3d7dc4ec61p-17 +
                   T * (0x1.bb47673661291p-29 + T * -0x1.a42599f86f826p-42))) +
         0x1.1111aad2a7016p+4);
    h[sCH3CO3] = 0x1.bb2d881f811c1p+6 *
                 (T * (0x1.cd472a8befb7ep+1 +
                       T * (0x1.ba7fe7cc1486ap-7 +
                            T * (-0x1.d1f1ccf156711p-18 +
                                 T * (0x1.243e62fec4abap-29 +
                                      T * -0x1.43e9a30471024p-42)))) -
                  0x1.6df21182a9931p+14);
    cp[sCH3CO3] =
        0x1.bb2d881f811c1p+6 *
        (0x1.cd472a8befb7ep+1 +
         T * (0x1.ba7fe7cc1486ap-6 +
              T * (-0x1.5d7559b500d4dp-16 +
                   T * (0x1.243e62fec4abap-27 + T * -0x1.94e40bc58d42dp-40))));
    s[sCH3CO3] =
        0x1.bb2d881f811c1p+6 *
        (0x1.cd472a8befb7ep+1 * log(T) +
         T * (0x1.ba7fe7cc1486ap-6 +
              T * (-0x1.5d7559b500d4dp-17 +
                   T * (0x1.85a883fe5b8f8p-29 + T * -0x1.94e40bc58d42dp-42))) +
         0x1.66729e17ad9bbp+3);
    h[sH2CC] = 0x1.3f584141b0716p+8 *
               (T * (0x1.a409c6525f486p+1 +
                     T * (0x1.c935e59362e0ep-9 +
                          T * (-0x1.aae7fe09f2e5bp-21 +
                               T * (-0x1.4cb95e2d8a5eap-32 +
                                    T * 0x1.ba34d6846b668p-43)))) +
                0x1.7bdb96872b021p+15);
    cp[sH2CC] =
        0x1.3f584141b0716p+8 *
        (0x1.a409c6525f486p+1 +
         T * (0x1.c935e59362e0ep-8 +
              T * (-0x1.402dfe87762c4p-19 +
                   T * (-0x1.4cb95e2d8a5eap-30 + T * 0x1.14610612c3201p-40))));
    s[sH2CC] =
        0x1.3f584141b0716p+8 *
        (0x1.a409c6525f486p+1 * log(T) +
         T * (0x1.c935e59362e0ep-8 +
              T * (-0x1.402dfe87762c4p-20 +
                   T * (-0x1.bba1d2e76328dp-32 + T * 0x1.14610612c3201p-42))) +
         0x1.7ae7afa722186p+2);
    h[sC2H2] = 0x1.3f584141b0716p+8 *
               (T * (0x1.9e0b436642656p-1 +
                     T * (0x1.7ec18d70c4cddp-7 +
                          T * (-0x1.8d40e0138fbadp-17 +
                               T * (0x1.e14c92ee6f7bp-28 +
                                    T * -0x1.de8cbe063b552p-40)))) +
                0x1.9cf3ec56d5cfbp+14);
    cp[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.9e0b436642656p-1 +
         T * (0x1.7ec18d70c4cddp-6 +
              T * (-0x1.29f0a80eabcc2p-15 +
                   T * (0x1.e14c92ee6f7bp-26 + T * -0x1.2b17f6c3e5153p-37))));
    s[sC2H2] =
        0x1.3f584141b0716p+8 *
        (0x1.9e0b436642656p-1 * log(T) +
         T * (0x1.7ec18d70c4cddp-6 +
              T * (-0x1.29f0a80eabcc2p-16 +
                   T * (0x1.40ddb7499fa75p-27 + T * -0x1.2b17f6c3e5153p-39))) +
         0x1.be11d39ccaa68p+3);
    h[sOCHO] = 0x1.716240449121bp+7 *
               (T * (0x1.2c0c705b87bb2p+2 +
                     T * (-0x1.0fe3f19732925p-9 +
                          T * (0x1.1d495841a0f01p-17 +
                               T * (-0x1.e8b8f1c03c047p-28 +
                                    T * 0x1.25ec6236c551ep-39)))) -
                0x1.096ad0ff97247p+14);
    cp[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.2c0c705b87bb2p+2 +
         T * (-0x1.0fe3f19732925p-8 +
              T * (0x1.abee046271682p-16 +
                   T * (-0x1.e8b8f1c03c047p-26 + T * 0x1.6f677ac476a65p-37))));
    s[sOCHO] =
        0x1.716240449121bp+7 *
        (0x1.2c0c705b87bb2p+2 * log(T) +
         T * (-0x1.0fe3f19732925p-8 +
              T * (0x1.abee046271682p-17 +
                   T * (-0x1.45d0a12ad2adap-27 + T * 0x1.6f677ac476a65p-39))) +
         0x1.123164fcd9dadp+2);
    h[sSC2H4OH] = 0x1.710a1c4780b3cp+7 *
                  (T * (0x1.0e42e33eff195p+2 +
                        T * (0x1.4fa8ac128c597p-9 +
                             T * (0x1.85a9ea206ac6dp-17 +
                                  T * (-0x1.a6938cc0e5c65p-27 +
                                       T * 0x1.1b24133ff87cep-38)))) -
                   0x1.006850abb44e5p+13);
    cp[sSC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.0e42e33eff195p+2 +
         T * (0x1.4fa8ac128c597p-8 +
              T * (0x1.243f6f9850152p-15 +
                   T * (-0x1.a6938cc0e5c65p-25 + T * 0x1.61ed180ff69c1p-36))));
    s[sSC2H4OH] =
        0x1.710a1c4780b3cp+7 *
        (0x1.0e42e33eff195p+2 * log(T) +
         T * (0x1.4fa8ac128c597p-8 +
              T * (0x1.243f6f9850152p-16 +
                   T * (-0x1.19b7b32b43d99p-26 + T * 0x1.61ed180ff69c1p-38))) +
         0x1.008945f9df549p+3);
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
    s[sC2H5OH] =
        0x1.68f6f3733777dp+7 *
        (0x1.36f4decf2dd02p+2 * log(T) +
         T * (-0x1.ea3b5dff2daaep-9 +
              T * (0x1.23bc84b8a9a7fp-15 +
                   T * (-0x1.fbb14e2f93e7dp-26 + T * 0x1.35592d14ec954p-37))) +
         0x1.3351958969a0bp+2);
    h[sCH3CO2] = 0x1.19a2d4e8daa2p+7 *
                 (T * (0x1.5fd92e84f8a2cp+0 +
                       T * (0x1.9826a84d76067p-7 +
                            T * (-0x1.85ec221a8c4b8p-18 +
                                 T * (0x1.ad5be45068e73p-30 +
                                      T * -0x1.999c26db8915fp-43)))) -
                  0x1.a9840f5c28f5cp+14);
    cp[sCH3CO2] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.5fd92e84f8a2cp+0 +
         T * (0x1.9826a84d76067p-6 +
              T * (-0x1.24711993e938ap-16 +
                   T * (0x1.ad5be45068e73p-28 + T * -0x1.0001984935adbp-40))));
    s[sCH3CO2] =
        0x1.19a2d4e8daa2p+7 *
        (0x1.5fd92e84f8a2cp+0 * log(T) +
         T * (0x1.9826a84d76067p-6 +
              T * (-0x1.24711993e938ap-17 +
                   T * (0x1.1e3d42e045ef7p-29 + T * -0x1.0001984935adbp-42))) +
         0x1.223fac8889a6ep+4);
    h[sCH3CHO] = 0x1.797bdf23f83eep+7 *
                 (T * (0x1.2eaf76e6106abp+2 +
                       T * (-0x1.a28ce427d2efep-10 +
                            T * (0x1.09d5a4c9ce534p-16 +
                                 T * (-0x1.ed90d262d7aa3p-27 +
                                      T * 0x1.34a728840b03p-38)))) -
                  0x1.511383126e979p+14);
    cp[sCH3CHO] =
        0x1.797bdf23f83eep+7 *
        (0x1.2eaf76e6106abp+2 +
         T * (-0x1.a28ce427d2efep-9 +
              T * (0x1.8ec0772eb57cep-15 +
                   T * (-0x1.ed90d262d7aa3p-25 + T * 0x1.81d0f2a50dc3bp-36))));
    s[sCH3CHO] =
        0x1.797bdf23f83eep+7 *
        (0x1.2eaf76e6106abp+2 * log(T) +
         T * (-0x1.a28ce427d2efep-9 +
              T * (0x1.8ec0772eb57cep-16 +
                   T * (-0x1.490b36ec8fc6dp-26 + T * 0x1.81d0f2a50dc3bp-38))) +
         0x1.0697d0005df3dp+2);
    h[sHCCO] = 0x1.954e7e0035abdp+7 *
               (T * (0x1.e046c2313d64cp+0 +
                     T * (0x1.6a6c46e67385fp-7 +
                          T * (-0x1.916377b66fae8p-17 +
                               T * (0x1.0656ba9f81848p-27 +
                                    T * -0x1.1d14b397f2b55p-39)))) +
                0x1.3b0d89374bc6ap+14);
    cp[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.e046c2313d64cp+0 +
         T * (0x1.6a6c46e67385fp-6 +
              T * (-0x1.2d0a99c8d3c2ep-15 +
                   T * (0x1.0656ba9f81848p-25 + T * -0x1.6459e07def62ap-37))));
    s[sHCCO] =
        0x1.954e7e0035abdp+7 *
        (0x1.e046c2313d64cp+0 * log(T) +
         T * (0x1.6a6c46e67385fp-6 +
              T * (-0x1.2d0a99c8d3c2ep-16 +
                   T * (0x1.5dc8f8d4acb0bp-27 + T * -0x1.6459e07def62ap-39))) +
         0x1.b64c6c54bcf0bp+3);
    h[sHCOH] = 0x1.14e89ec288e65p+8 *
               (T * (-0x1.692957fd97f0dp+1 +
                     T * (0x1.24b9e3ac865dp-6 +
                          T * (-0x1.a9fc890d5d967p-17 +
                               T * (0x1.3fe63e49d77bp-28 +
                                    T * -0x1.85838ecaa1054p-41)))) +
                0x1.60fd566cf41f2p+13);
    cp[sHCOH] =
        0x1.14e89ec288e65p+8 *
        (-0x1.692957fd97f0dp+1 +
         T * (0x1.24b9e3ac865dp-5 +
              T * (-0x1.3f7d66ca0630dp-15 +
                   T * (0x1.3fe63e49d77bp-26 + T * -0x1.e6e4727d49468p-39))));
    s[sHCOH] =
        0x1.14e89ec288e65p+8 *
        (-0x1.692957fd97f0dp+1 * log(T) +
         T * (0x1.24b9e3ac865dp-5 +
              T * (-0x1.3f7d66ca0630dp-16 +
                   T * (0x1.aa88530d1f4ebp-28 + T * -0x1.e6e4727d49468p-41))) +
         0x1.16ca4aea091dbp+5);
    h[sC2H3OH] = 0x1.797bdf23f83eep+7 *
                 (T * (-0x1.061651fbf6302p-3 +
                       T * (0x1.154dde69b422p-6 +
                            T * (-0x1.71d1f20be878bp-17 +
                                 T * (0x1.1b39a3e362fdfp-28 +
                                      T * -0x1.68371dd57e9ep-41)))) -
                  0x1.f3bba29c779a7p+13);
    cp[sC2H3OH] =
        0x1.797bdf23f83eep+7 *
        (-0x1.061651fbf6302p-3 +
         T * (0x1.154dde69b422p-5 +
              T * (-0x1.155d7588ee5a8p-15 +
                   T * (0x1.1b39a3e362fdfp-26 + T * -0x1.c244e54ade458p-39))));
    s[sC2H3OH] =
        0x1.797bdf23f83eep+7 *
        (-0x1.061651fbf6302p-3 * log(T) +
         T * (0x1.154dde69b422p-5 +
              T * (-0x1.155d7588ee5a8p-16 +
                   T * (0x1.79a22fd9d9529p-28 + T * -0x1.c244e54ade458p-41))) +
         0x1.70b3a6a5f196bp+4);
    h[sC2H5O] = 0x1.710a1c4780b3cp+7 *
                (T * (0x1.13ace174fa7dep+2 +
                      T * (0x1.a46526c2aede2p-9 +
                           T * (0x1.5c4ab3c09e1efp-17 +
                                T * (-0x1.7410ad01e0d3dp-27 +
                                     T * 0x1.e648067b59508p-39)))) -
                 0x1.a95813a92a305p+11);
    cp[sC2H5O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.13ace174fa7dep+2 +
         T * (0x1.a46526c2aede2p-8 +
              T * (0x1.053806d076973p-15 +
                   T * (-0x1.7410ad01e0d3dp-25 + T * 0x1.2fed040d17d25p-36))));
    s[sC2H5O] =
        0x1.710a1c4780b3cp+7 *
        (0x1.13ace174fa7dep+2 * log(T) +
         T * (0x1.a46526c2aede2p-8 +
              T * (0x1.053806d076973p-16 +
                   T * (-0x1.f0163c02811a7p-27 + T * 0x1.2fed040d17d25p-38))) +
         0x1.79c3ee6c59c57p+2);
    h[sC2H5CO] = 0x1.2360aa5b5ab7bp+7 *
                 (T * (0x1.90765bbc69525p+2 +
                       T * (-0x1.2caee39a7c123p-8 +
                            T * (0x1.a9b04699a6cecp-16 +
                                 T * (-0x1.84ea6d44e87ccp-26 +
                                      T * 0x1.e73b0d5ab9fbcp-38)))) -
                  0x1.71c2a32f44913p+12);
    cp[sC2H5CO] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.90765bbc69525p+2 +
         T * (-0x1.2caee39a7c123p-7 +
              T * (0x1.3f4434f33d1b1p-14 +
                   T * (-0x1.84ea6d44e87ccp-24 + T * 0x1.3084e858b43d5p-35))));
    s[sC2H5CO] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.90765bbc69525p+2 * log(T) +
         T * (-0x1.2caee39a7c123p-7 +
              T * (0x1.3f4434f33d1b1p-15 +
                   T * (-0x1.0346f38345a88p-25 + T * 0x1.3084e858b43d5p-37))) +
         0x1.1ddcf87e81654p+1);
    h[sC2H5CHO] = 0x1.1e520994b8387p+7 *
                  (T * (0x1.0fb2f16430d71p+2 +
                        T * (0x1.b5f994f68394dp-9 +
                             T * (0x1.13e50275637bbp-16 +
                                  T * (-0x1.209da7dc3b348p-26 +
                                       T * 0x1.782368940785bp-38)))) -
                   0x1.794d33eab367ap+14);
    cp[sC2H5CHO] =
        0x1.1e520994b8387p+7 *
        (0x1.0fb2f16430d71p+2 +
         T * (0x1.b5f994f68394dp-8 +
              T * (0x1.9dd783b015398p-15 +
                   T * (-0x1.209da7dc3b348p-24 + T * 0x1.d62c42b909671p-36))));
    s[sC2H5CHO] =
        0x1.1e520994b8387p+7 *
        (0x1.0fb2f16430d71p+2 * log(T) +
         T * (0x1.b5f994f68394dp-8 +
              T * (0x1.9dd783b015398p-16 +
                   T * (-0x1.80d23525a446p-26 + T * 0x1.d62c42b909671p-38))) +
         0x1.ba129b0d37202p+2);
    h[sC2H] = 0x1.4c34d172a6531p+8 *
              (T * (0x1.7307d70ef007ep+1 +
                    T * (0x1.b3c6d48ccb044p-8 +
                         T * (-0x1.39feb27d48311p-17 +
                              T * (0x1.f154bdde62346p-28 +
                                   T * -0x1.2e979b8febbeep-39)))) +
               0x1.05f59ae147ae1p+16);
    cp[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.7307d70ef007ep+1 +
         T * (0x1.b3c6d48ccb044p-7 +
              T * (-0x1.d6fe0bbbec499p-16 +
                   T * (0x1.f154bdde62346p-26 + T * -0x1.7a3d8273e6ae9p-37))));
    s[sC2H] =
        0x1.4c34d172a6531p+8 *
        (0x1.7307d70ef007ep+1 * log(T) +
         T * (0x1.b3c6d48ccb044p-7 +
              T * (-0x1.d6fe0bbbec499p-17 +
                   T * (0x1.4b8dd3e996cd9p-27 + T * -0x1.7a3d8273e6ae9p-39))) +
         0x1.8bded81225469p+2);
    h[sCH3COCH3] = 0x1.1e520994b8387p+7 *
                   (T * (0x1.639be172763f2p+2 +
                         T * (-0x1.7410cd0490958p-10 +
                              T * (0x1.8aab3bafe58cfp-16 +
                                   T * (-0x1.792785d3b8461p-26 +
                                        T * 0x1.deeabc6b21955p-38)))) -
                    0x1.b2e2283e425afp+14);
    cp[sCH3COCH3] =
        0x1.1e520994b8387p+7 *
        (0x1.639be172763f2p+2 +
         T * (-0x1.7410cd0490958p-9 +
              T * (0x1.28006cc3ec29bp-14 +
                   T * (-0x1.792785d3b8461p-24 + T * 0x1.2b52b5c2f4fd5p-35))));
    s[sCH3COCH3] =
        0x1.1e520994b8387p+7 *
        (0x1.639be172763f2p+2 * log(T) +
         T * (-0x1.7410cd0490958p-9 +
              T * (0x1.28006cc3ec29bp-15 +
                   T * (-0x1.f6df5d1a4b081p-26 + T * 0x1.2b52b5c2f4fd5p-37))) +
         0x1.28e8b9a7d6f42p+1);
    h[sCH3COCH2] = 0x1.2360aa5b5ab7bp+7 *
                   (T * (0x1.2ceb785e8b7bap+2 +
                         T * (0x1.698824aa40f53p-9 +
                              T * (0x1.de283ec99d568p-17 +
                                   T * (-0x1.fed3b1c44b5a2p-27 +
                                        T * 0x1.52bc0bd9eb88cp-38)))) -
                    0x1.7287474fb54ap+12);
    cp[sCH3COCH2] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.2ceb785e8b7bap+2 +
         T * (0x1.698824aa40f53p-8 +
              T * (0x1.669e2f173600ep-15 +
                   T * (-0x1.fed3b1c44b5a2p-25 + T * 0x1.a76b0ed0666afp-36))));
    s[sCH3COCH2] =
        0x1.2360aa5b5ab7bp+7 *
        (0x1.2ceb785e8b7bap+2 * log(T) +
         T * (0x1.698824aa40f53p-8 +
              T * (0x1.669e2f173600ep-16 +
                   T * (-0x1.548d212d87917p-26 + T * 0x1.a76b0ed0666afp-38))) +
         0x1.c846e023c1a63p+2);
    h[sCH2Y] = 0x1.286500f4d652p+9 *
               (T * (0x1.0c5f3e89a88adp+2 +
                     T * (-0x1.31891ed60ea9ep-10 +
                          T * (0x1.6ced6e76a17d1p-19 +
                               T * (-0x1.c799ba2b7a484p-30 +
                                    T * 0x1.b31f5af2aea77p-42)))) +
                0x1.897c72fec56d6p+15);
    cp[sCH2Y] =
        0x1.286500f4d652p+9 *
        (0x1.0c5f3e89a88adp+2 +
         T * (-0x1.31891ed60ea9ep-9 +
              T * (0x1.11b212d8f91ddp-17 +
                   T * (-0x1.c799ba2b7a484p-28 + T * 0x1.0ff398d7ad28ap-39))));
    s[sCH2Y] =
        0x1.286500f4d652p+9 *
        (0x1.0c5f3e89a88adp+2 * log(T) +
         T * (-0x1.31891ed60ea9ep-9 +
              T * (0x1.11b212d8f91ddp-18 +
                   T * (-0x1.2fbbd17251858p-29 + T * 0x1.0ff398d7ad28ap-41))) -
         0x1.7e53f5a080939p-1);
    h[sOCH2O2H] = 0x1.07cf090c05134p+7 *
                  (T * (0x1.f02fe3f359ff5p+0 +
                        T * (0x1.edebe447c645dp-7 +
                             T * (-0x1.23fba667f67ffp-17 +
                                  T * (0x1.781d2ada0882bp-29 +
                                       T * -0x1.918644c15e53ep-42)))) -
                   0x1.afc54cccccccdp+13);
    cp[sOCH2O2H] =
        0x1.07cf090c05134p+7 *
        (0x1.f02fe3f359ff5p+0 +
         T * (0x1.edebe447c645dp-6 +
              T * (-0x1.b5f9799bf1bfep-16 +
                   T * (0x1.781d2ada0882bp-27 + T * -0x1.f5e7d5f1b5e8dp-40))));
    s[sOCH2O2H] =
        0x1.07cf090c05134p+7 *
        (0x1.f02fe3f359ff5p+0 * log(T) +
         T * (0x1.edebe447c645dp-6 +
              T * (-0x1.b5f9799bf1bfep-17 +
                   T * (0x1.f57c3922b6039p-29 + T * -0x1.f5e7d5f1b5e8dp-42))) +
         0x1.2811343a9a2fcp+4);
    h[sHE] = 0x1.9c055f0ad864p+8 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sHE] =
        0x1.9c055f0ad864p+8 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sHE] = 0x1.9c055f0ad864p+8 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.db81b56eaeabcp-1);
    h[sAR] = 0x1.a0439e3ea879p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a0439e3ea879p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sAR] = 0x1.a0439e3ea879p+7 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.184c97fe63f3ap+2);
  } else {
    ComputeThermoData(h, cp, 300.0, s);
    for (i = 0; i < sEnd; i++) {
      h[i] = (T - 300.) * cp[i] + h[i];
      s[i] = log(T / 300) * cp[i] + s[i];
    }
  }
}

double MAX_C(double X1, double X2) { return ((X1 > X2) ? X1 : X2); }
double MIN_C(double X1, double X2) { return ((X2 > X1) ? X1 : X2); }

double GetPlogRateCoeff(double temp, double pressure, double lgt, double rt_inv,
                        double* PlogP, double* PlogA, double* PlogB,
                        double* PlogE, int np) {
  double kR, kR_l, kR_r;
  int i;
  if (pressure <= PlogP[0])
    kR = PlogA[0] * exp(PlogB[0] * lgt - PlogE[0] * rt_inv);
  else if (pressure >= PlogP[np - 1])
    kR = PlogA[np - 1] * exp(PlogB[np - 1] * lgt - PlogE[np - 1] * rt_inv);
  else {
    /* interpolate */
    for (i = 0; i < np; i++) {
      if (pressure <= PlogP[i])
        break;
    }

    kR_l = PlogA[i - 1] * exp(PlogB[i - 1] * lgt - PlogE[i - 1] * rt_inv);
    kR_r = PlogA[i] * exp(PlogB[i] * lgt - PlogE[i] * rt_inv);
    kR = exp(log(MAX_C(kR_l, 1e-60)) +
             (log(MAX_C(kR_r, 1e-60)) - log(MAX_C(kR_l, 1e-60))) /
                 (log(PlogP[i]) - log(PlogP[i - 1])) *
                 (log(pressure) - log(PlogP[i - 1])));
  }
  return MIN_C(kR, DBL_MAX);
}

const char* AramcoMech_DMEonly_74spec::GetMechanismResource() const {
  return R"MECHANISM(ELEMENTS
H  O  N  AR C  HE 
END
SPECIES
N2                CH3OCH3           CH4               H2                O2                CO                CO2               H2O               
H                 O                 OH                CH3               CH2O              C2H6              CH3OCH2           O2CH2OCH2O2H      
CH2OCH2O2H        HO2               H2O2              CH3OCH2O2         CH3O              CH3OH             CH3OCH2O          OCH2OCHO          
HO2CH2OCHO        HCO               O2CHO             HOCHO             HOCH2O            HOCH2OCO          CH2OH             CH3O2H            
CH3O2             HO2CHO            CH3OCH2O2H        CH2               CH                CH3OCHO           CH2OCHO           C2H4O1            
C2H4              C2H3              C2H5              OHY               CH3OCO            HOCH2O2H          HOCH2O2           C2H3CO            
C2H3CHO           CH3CO             C2H5O2H           C2H5O2            CH2CO             CH2CHO            CH3CO3H           CH3CO3            
H2CC              C2H2              OCHO              SC2H4OH           C2H5OH            CH3CO2            CH3CHO            HCCO              
HCOH              C2H3OH            C2H5O             C2H5CO            C2H5CHO           C2H               CH3COCH3          CH3COCH2          
CH2Y              OCH2O2H           HE                AR                
END
REACTIONS
H+O2<=>O+OH                              1.040e+14   -0.000  15286.00
O+H2<=>H+OH                              5.080e+04    2.670   6292.00
OH+H2<=>H+H2O                            4.380e+13   -0.000   6990.00
O+H2O<=>2OH                              2.970e+06    2.020  13400.00
H2+M<=>2H+M                              4.577e+19   -1.400 104400.00
CH4/2.00/ H2/2.50/ CO/1.90/ CO2/3.80/ H2O/12.00/ C2H6/3.00/ HE/0.83/ 
2O+M<=>O2+M                              6.165e+15   -0.500     -0.00
CH4/2.00/ H2/2.50/ CO/1.90/ CO2/3.80/ H2O/12.00/ C2H6/3.00/ HE/0.83/ AR/0.83/ 
O+H+M<=>OH+M                             4.714e+18   -1.000     -0.00
CH4/2.00/ H2/2.50/ CO/1.50/ CO2/2.00/ H2O/12.00/ C2H6/3.00/ HE/0.75/ AR/0.75/ 
H+OH+M<=>H2O+M                           3.500e+22   -2.000     -0.00
CH4/2.00/ H2/0.73/ H2O/3.65/ C2H6/3.00/ AR/0.38/ 
H+O2(+M)<=>HO2(+M)                       4.650e+12    0.440     -0.00
CH4/2.00/ H2/1.30/ CO/1.90/ CO2/3.80/ H2O/10.00/ C2H6/3.00/ HE/0.00/ AR/0.00/ 
     LOW  /  1.737e+19   -1.230     -0.00 /
     TROE/    0.67    1e-30     1e+30     1e+30 /
H+O2(+AR)<=>HO2(+AR)                     4.650e+12    0.440     -0.00
     LOW  /  6.810e+18   -1.200     -0.00 /
     TROE/     0.7    1e-30     1e+30     1e+30 /
H+O2(+HE)<=>HO2(+HE)                     4.650e+12    0.440     -0.00
     LOW  /  9.192e+18   -1.200     -0.00 /
     TROE/    0.59    1e-30     1e+30     1e+30 /
HO2+H<=>2OH                              7.079e+13    0.000    295.00
H2+O2<=>H+HO2                            5.176e+05    2.433  53502.00
HO2+O<=>OH+O2                            3.250e+13   -0.000     -0.00
HO2+OH<=>H2O+O2                          2.456e+13    0.000   -497.00
2HO2<=>H2O2+O2                           1.300e+11   -0.000  -1630.00
     DUPLICATE
2HO2<=>H2O2+O2                           3.658e+14   -0.000  12000.00
     DUPLICATE
H2O2(+H2O)<=>2OH(+H2O)                   2.000e+12    0.900  48749.00
     LOW  /  1.865e+25   -2.300  48749.00 /
!     FCCHECK/    0.49     1e-30   0.51     1e+30      0         0
     TROE/    0.51    1e-30     1e+30        /
H2O2(+M)<=>2OH(+M)                       2.000e+12    0.900  48749.00
N2/1.50/ H2/3.70/ O2/1.20/ CO/2.80/ CO2/1.60/ H2O/0.00/ H2O2/7.70/ HE/0.65/ 
     LOW  /  2.490e+24   -2.300  48749.00 /
!     FCCHECK/    0.57     1e-30   0.43     1e+30      0         0
     TROE/    0.43    1e-30     1e+30        /
H2O2+H<=>H2O+OH                          2.410e+13   -0.000   3970.00
H2O2+H<=>H2+HO2                          2.150e+10    1.000   6000.00
H2O2+O<=>OH+HO2                          9.550e+06    2.000   3970.00
H2O2+OH<=>H2O+HO2                        1.740e+12   -0.000    318.00
     DUPLICATE
H2O2+OH<=>H2O+HO2                        7.590e+13   -0.000   7269.00
     DUPLICATE
CO+O(+M)<=>CO2(+M)                       1.362e+10   -0.000   2384.00
H2/2.00/ CO/1.75/ CO2/3.60/ H2O/12.00/ HE/0.70/ AR/0.70/ 
     LOW  /  1.173e+24   -2.790   4191.00 /
!     FCCHECK/       0         0      0         0      1         0
     TROE/       1     1.00 10000000.00 10000000.00 /
CO+O2<=>CO2+O                            1.119e+12   -0.000  47700.00
CO+OH<=>CO2+H                            7.015e+04    2.053   -355.70
     DUPLICATE
CO+OH<=>CO2+H                            5.757e+12   -0.664    331.80
     DUPLICATE
CO+HO2<=>CO2+OH                          1.570e+05    2.180  17940.00
HCO+M<=>H+CO+M                           5.700e+11    0.660  14870.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/12.00/ C2H6/3.00/ 
HCO+O2<=>CO+HO2                          7.580e+12   -0.000    410.00
HCO+H<=>CO+H2                            7.340e+13   -0.000     -0.00
HCO+O<=>CO+OH                            3.020e+13   -0.000     -0.00
HCO+O<=>CO2+H                            3.000e+13   -0.000     -0.00
HCO+OH<=>CO+H2O                          1.020e+14   -0.000     -0.00
HCO+HO2=>CO2+H+OH                        3.000e+13   -0.000     -0.00
2HCO=>H2+2CO                             3.000e+12   -0.000     -0.00
HCO+CH3<=>CH4+CO                         2.650e+13   -0.000     -0.00
CH2O+O2<=>HCO+HO2                        8.070e+15   -0.000  53420.00
HCO+O2<=>O2CHO                           1.200e+11   -0.000  -1100.00
CH2O+O2CHO<=>HCO+HO2CHO                  1.990e+12   -0.000  11660.00
OCHO+OH<=>HO2CHO                         2.000e+13   -0.000     -0.00
H+CO2<=>OCHO                             7.500e+13   -0.000  29000.00
2HCO<=>CH2O+CO                           1.800e+13   -0.000     -0.00
OHY+H2O<=>OH+H2O                         5.930e+12    0.500   -860.00
OHY+H2<=>OH+H2                           2.950e+12    0.500   -444.00
OHY+N2<=>OH+N2                           1.080e+11    0.500  -1242.00
OHY+OH<=>2OH                             6.010e+12    0.500   -764.00
OHY+H<=>OH+H                             1.310e+12    0.500   -167.00
OHY+AR<=>OH+AR                           1.690e+12   -0.000   4135.00
OHY<=>OH                                 1.450e+06   -0.000     -0.00
OHY+O2<=>OH+O2                           2.100e+12    0.500   -478.00
OHY+CO2<=>OH+CO2                         2.750e+12    0.500   -968.00
OHY+CO<=>OH+CO                           3.230e+12    0.500   -787.00
OHY+CH4<=>OH+CH4                         3.360e+12    0.500   -635.00
CH+O2<=>CO+OHY                           4.040e+13   -0.000     -0.00
HCO+H(+M)<=>CH2O(+M)                     1.090e+12    0.480   -260.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  1.350e+24   -2.570   1425.00 /
     TROE/  0.7824      271      2755      6570 /
CO+H2(+M)<=>CH2O(+M)                     4.300e+07    1.500  79600.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  5.070e+27   -3.420  84348.00 /
     TROE/   0.932      197      1540     10300 /
CH2O+OH<=>HCO+H2O                        7.820e+07    1.630  -1055.00
CH2O+H<=>HCO+H2                          5.740e+07    1.900   2740.00
CH2O+O<=>HCO+OH                          6.260e+09    1.150   2260.00
CH2O+CH3<=>HCO+CH4                       3.830e+01    3.360   4312.00
CH2O+HO2<=>HCO+H2O2                      1.880e+04    2.700  11520.00
CH2O+OH<=>HOCH2O                         4.500e+15   -1.100     -0.00
HOCH2O<=>HOCHO+H                         1.000e+14   -0.000  14900.00
HOCHO<=>CO+H2O                           2.450e+12   -0.000  60470.00
HOCHO<=>CO2+H2                           2.950e+09   -0.000  48520.00
OCHO+HO2<=>HOCHO+O2                      3.500e+10   -0.000  -3275.00
HOCHO+OH=>H2O+CO2+H                      2.620e+06    2.060    916.00
HOCHO+OH=>H2O+CO+OH                      1.850e+07    1.510   -962.00
HOCHO+H=>H2+CO2+H                        4.240e+06    2.100   4868.00
HOCHO+H=>H2+CO+OH                        6.030e+13   -0.350   2988.00
HOCHO+CH3=>CH4+CO+OH                     3.900e-07    5.800   2200.00
OCHO+H2O2<=>HOCHO+HO2                    2.400e+12   -0.000  10000.00
HOCHO+HO2=>H2O2+CO+OH                    1.000e+12   -0.000  11920.00
HOCHO+O=>CO+2OH                          1.770e+18   -1.900   2975.00
CH2O+OCHO<=>HOCHO+HCO                    5.600e+12   -0.000  13600.00
CH3O(+M)<=>CH2O+H(+M)                    6.800e+13   -0.000  26170.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ 
     LOW  /  1.867e+25   -3.000  24307.00 /
     TROE/     0.9     2500      1300     1e+99 /
CH3O+O2<=>CH2O+HO2                       4.380e-19    9.500  -5501.00
CH2O+CH3O<=>CH3OH+HCO                    6.620e+11   -0.000   2294.00
CH3+CH3OH<=>CH4+CH3O                     1.440e+01    3.100   6935.00
CH3O+CH3<=>CH2O+CH4                      1.200e+13   -0.000     -0.00
CH3O+H<=>CH2O+H2                         2.000e+13   -0.000     -0.00
CH3O+HO2<=>CH2O+H2O2                     3.010e+11   -0.000     -0.00
CH2O+H(+M)<=>CH2OH(+M)                   5.400e+11    0.454   3600.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ 
     LOW  /  1.270e+32   -4.820   6530.00 /
     TROE/  0.7187      103      1291      4160 /
CH2OH+O2<=>CH2O+HO2                      1.510e+15   -1.000      0.00
     DUPLICATE
CH2OH+O2<=>CH2O+HO2                      2.410e+14   -0.000   5017.00
     DUPLICATE
CH2OH+H<=>CH2O+H2                        6.000e+12   -0.000     -0.00
CH2OH+HO2<=>CH2O+H2O2                    1.200e+13   -0.000     -0.00
CH2OH+HCO<=>2CH2O                        1.800e+14   -0.000     -0.00
CH2OH+CH3O<=>CH2O+CH3OH                  2.400e+13   -0.000     -0.00
CH3OH+HCO<=>CH2OH+CH2O                   9.630e+03    2.900  13110.00
OH+CH2OH<=>H2O+CH2O                      2.400e+13   -0.000     -0.00
O+CH2OH<=>OH+CH2O                        4.200e+13   -0.000     -0.00
2CH2OH<=>CH2O+CH3OH                      3.000e+12   -0.000     -0.00
CH2OH+HO2<=>HOCH2O+OH                    1.000e+13   -0.000     -0.00
CH2O+HO2<=>OCH2O2H                       1.500e+11   -0.000  11900.00
OCH2O2H<=>HOCH2O2                        3.000e+11   -0.000   8600.00
HOCH2O2+HO2<=>HOCH2O2H+O2                3.500e+10   -0.000  -3275.00
HOCH2O+OH<=>HOCH2O2H                     1.000e+13   -0.000     -0.00
CH3OH(+M)<=>CH3+OH(+M)                   2.084e+18   -0.615  92540.60

     LOW  /  1.500e+43   -6.995  97992.20 /
     TROE/  -0.4748    35580      1116      9023 /
CH3OH(+M)<=>CH2Y+H2O(+M)                 3.121e+18   -1.017  91712.00

     LOW  /  1.430e+47   -8.227  99417.10 /
     TROE/   2.545     3290     47320     47110 /
CH3OH(+M)<=>CH2OH+H(+M)                  7.896e-03    5.038  84467.40

     LOW  /  3.390e+42   -7.244 105230.30 /
     TROE/  -73.91    37050     41500      5220 /
CH3OH+H<=>CH2OH+H2                       3.070e+05    2.550   5440.00
CH3OH+H<=>CH3O+H2                        1.990e+05    2.560  10300.00
CH3OH+O<=>CH2OH+OH                       3.880e+05    2.500   3080.00
CH3OH+OH<=>CH2OH+H2O                     3.080e+04    2.650   -806.70
CH3OH+OH<=>CH3O+H2O                      1.500e+02    3.030   -763.00
CH3OH+O2<=>CH2OH+HO2                     2.050e+13   -0.000  44900.00
CH3OH+HO2<=>CH2OH+H2O2                   1.080e+04    2.550  10530.00
CH3OH+CH3<=>CH2OH+CH4                    3.190e+01    3.170   7172.00
CH3O+CH3OH<=>CH2OH+CH3OH                 3.000e+11   -0.000   4074.00
2CH3O<=>CH3OH+CH2O                       6.030e+13   -0.000     -0.00
CH3+H(+M)<=>CH4(+M)                      1.270e+16   -0.630    383.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  2.477e+33   -4.760   2440.00 /
     TROE/   0.783       74      2941      6964 /
CH4+H<=>CH3+H2                           6.140e+05    2.500   9587.00
CH4+OH<=>CH3+H2O                         5.830e+04    2.600   2190.00
CH4+O<=>CH3+OH                           1.020e+09    1.500   8600.00
CH4+HO2<=>CH3+H2O2                       1.695e+01    3.740  21010.00
CH4+CH2<=>2CH3                           2.460e+06    2.000   8270.00
CH3+OH<=>CH2Y+H2O                        4.936e+14   -0.669   -445.80
PLOG /     0.010 4.936e+14   -0.669   -445.80 /
PLOG /     0.100 1.207e+15   -0.778   -175.60 /
PLOG /     1.000 5.282e+17   -1.518   1772.00 /
PLOG /    10.000 4.788e+23   -3.155   7003.00 /
PLOG /   100.000 8.433e+19   -1.962   8244.00 /
CH3+OH<=>CH2O+H2                         3.502e+05    1.441  -3244.00
PLOG /     0.010 3.502e+05    1.441  -3244.00 /
PLOG /     0.100 8.854e+05    1.327  -2975.00 /
PLOG /     1.000 1.650e+07    0.973  -2010.00 /
PLOG /    10.000 5.374e+09    0.287    280.00 /
PLOG /   100.000 9.494e+18   -2.199   9769.00 /
CH3+OH<=>CH2OH+H                         1.621e+10    0.965   3214.00
PLOG /     0.010 1.621e+10    0.965   3214.00 /
PLOG /     0.100 1.807e+10    0.950   3247.00 /
PLOG /     1.000 4.686e+10    0.833   3566.00 /
PLOG /    10.000 1.525e+13    0.134   5641.00 /
PLOG /   100.000 3.590e+14   -0.186   8601.00 /
CH3+OH<=>H+CH3O                          1.186e+09    1.016  11940.00
PLOG /     0.010 1.186e+09    1.016  11940.00 /
PLOG /     0.100 1.188e+09    1.016  11940.00 /
PLOG /     1.000 1.230e+09    1.011  11950.00 /
PLOG /    10.000 1.798e+09    0.965  12060.00 /
PLOG /   100.000 5.242e+10    0.551  13070.00 /
CH3+OH<=>HCOH+H2                         8.674e+08    0.787  -3046.00
PLOG /     0.010 8.674e+08    0.787  -3046.00 /
PLOG /     0.100 3.115e+09    0.630  -2669.00 /
PLOG /     1.000 1.557e+11    0.156  -1368.00 /
PLOG /    10.000 1.704e+21   -2.641   6412.00 /
PLOG /   100.000 7.250e+20   -2.402   9639.00 /
HCOH+OH<=>HCO+H2O                        2.000e+13   -0.000     -0.00
HCOH+H<=>CH2O+H                          2.000e+14   -0.000     -0.00
HCOH+O=>CO2+2H                           5.000e+13   -0.000     -0.00
HCOH+O=>CO+OH+H                          3.000e+13   -0.000     -0.00
HCOH+O2=>CO2+H+OH                        5.000e+12   -0.000     -0.00
HCOH+O2<=>CO2+H2O                        3.000e+13   -0.000     -0.00
CH3+HO2<=>CH3O+OH                        1.000e+12    0.269   -687.50
CH3+HO2<=>CH4+O2                         1.160e+05    2.230  -3022.00
CH3+O<=>CH2O+H                           5.540e+13    0.050   -136.00
CH3+O2<=>CH3O+O                          7.546e+12   -0.000  28320.00
CH3+O2<=>CH2O+OH                         2.641e+00    3.283   8105.00
CH3+O2(+M)<=>CH3O2(+M)                   7.812e+09    0.900     -0.00

     LOW  /  6.850e+24   -3.000     -0.00 /
     TROE/     0.6     1000        70      1700 /
CH3O2+CH2O<=>CH3O2H+HCO                  1.990e+12   -0.000  11660.00
CH4+CH3O2<=>CH3+CH3O2H                   9.600e-01    3.770  17810.00
CH3OH+CH3O2<=>CH2OH+CH3O2H               1.810e+12   -0.000  13710.00
CH3O2+CH3<=>2CH3O                        5.080e+12   -0.000  -1411.00
CH3O2+HO2<=>CH3O2H+O2                    2.470e+11    0.000  -1570.00
2CH3O2=>CH2O+CH3OH+O2                    3.110e+14   -1.610  -1051.00
2CH3O2=>O2+2CH3O                         1.400e+16   -1.610   1860.00
CH3O2+H<=>CH3O+OH                        9.600e+13   -0.000     -0.00
CH3O2+O<=>CH3O+O2                        3.600e+13   -0.000     -0.00
CH3O2+OH<=>CH3OH+O2                      6.000e+13   -0.000     -0.00
CH3O2H<=>CH3O+OH                         6.310e+14   -0.000  42300.00
CH2Y+N2<=>CH2+N2                         1.500e+13    0.000    600.00
CH2Y+AR<=>CH2+AR                         9.000e+12    0.000    600.00
CH2Y+H<=>CH+H2                           3.000e+13   -0.000     -0.00
CH2Y+O<=>CO+H2                           1.500e+13   -0.000     -0.00
CH2Y+O<=>HCO+H                           1.500e+13   -0.000     -0.00
CH2Y+OH<=>CH2O+H                         3.000e+13   -0.000     -0.00
CH2Y+H2<=>CH3+H                          7.000e+13   -0.000     -0.00
CH2Y+O2=>H+OH+CO                         2.800e+13   -0.000     -0.00
CH2Y+O2<=>CO+H2O                         1.200e+13   -0.000     -0.00
CH2Y+H2O<=>CH2+H2O                       3.000e+13   -0.000     -0.00
CH2Y+CO<=>CH2+CO                         9.000e+12   -0.000     -0.00
CH2Y+CO2<=>CH2+CO2                       7.000e+12   -0.000     -0.00
CH2Y+CO2<=>CH2O+CO                       1.400e+13   -0.000     -0.00
CH2+H(+M)<=>CH3(+M)                      2.500e+16   -0.800     -0.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  3.200e+27   -3.140   1230.00 /
     TROE/    0.68       78      1995      5590 /
CH2+O2<=>HCO+OH                          1.060e+13   -0.000   1500.00
CH2+O2=>CO2+2H                           2.640e+12   -0.000   1500.00
CH2+O=>CO+2H                             5.000e+13   -0.000     -0.00
CH2+H<=>CH+H2                            1.000e+18   -1.560     -0.00
     DUPLICATE
CH2+H<=>CH+H2                            2.700e+11    0.670  25700.00
     DUPLICATE
CH2+OH<=>CH+H2O                          1.130e+07    2.000   3000.00
CH+O2<=>HCO+O                            3.300e+13   -0.000     -0.00
CH+O<=>CO+H                              5.700e+13   -0.000     -0.00
CH+OH<=>HCO+H                            3.000e+13   -0.000     -0.00
CH+H2O<=>H+CH2O                          1.713e+13   -0.000   -755.00
CH+CO2<=>HCO+CO                          1.700e+12    0.000    685.00
2CH3(+M)<=>C2H6(+M)                      2.277e+15   -0.690    174.90
CO/2.00/ CO2/3.00/ H2O/5.00/ 
     LOW  /  8.054e+31   -3.750    981.60 /
     TROE/       0      570     1e+30     1e+30 /
C2H5+H(+M)<=>C2H6(+M)                    5.210e+17   -0.990   1580.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  1.990e+41   -7.080   6685.00 /
     TROE/   0.842      125      2219      6882 /
C2H6+H<=>C2H5+H2                         1.150e+08    1.900   7530.00
C2H6+O<=>C2H5+OH                         3.550e+06    2.400   5830.00
C2H6+OH<=>C2H5+H2O                       1.480e+07    1.900    950.00
C2H6+O2<=>C2H5+HO2                       6.030e+13   -0.000  51870.00
C2H6+CH3<=>C2H5+CH4                      5.480e-01    4.000   8280.00
C2H6+HO2<=>C2H5+H2O2                     3.460e+01    3.610  16920.00
C2H6+CH3O2<=>C2H5+CH3O2H                 1.940e+01    3.640  17100.00
C2H6+CH3O<=>C2H5+CH3OH                   2.410e+11    0.000   7090.00
C2H6+CH<=>C2H5+CH2                       1.100e+14    0.000   -260.00
CH2Y+C2H6<=>CH3+C2H5                     1.200e+14   -0.000     -0.00
C2H4+H(+M)<=>C2H5(+M)                    9.569e+08    1.463   1355.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ AR/0.70/ 
     LOW  /  1.419e+39   -6.642   5769.00 /
     TROE/  -0.569      299     -9147     152.4 /
H2+CH3O2<=>H+CH3O2H                      1.500e+14   -0.000  26030.00
H2+C2H5O2<=>H+C2H5O2H                    1.500e+14   -0.000  26030.00
2C2H4<=>C2H5+C2H3                        4.820e+14   -0.000  71530.00
CH3+C2H5<=>CH4+C2H4                      1.180e+04    2.450  -2921.00
2CH3<=>H+C2H5                            4.740e+12    0.105  10664.30
PLOG /     0.010 4.740e+12    0.105  10664.30 /
PLOG /     0.100 2.570e+13   -0.096  11406.10 /
PLOG /     1.000 3.100e+14   -0.362  13372.50 /
PLOG /    10.000 2.150e+10    0.885  13532.50 /
PLOG /   100.000 1.032e+02    3.230  11236.10 /
C2H5+H<=>C2H4+H2                         2.000e+12   -0.000     -0.00
C2H5+O<=>CH3CHO+H                        1.100e+14   -0.000     -0.00
C2H5+HO2<=>C2H5O+OH                      1.100e+13   -0.000     -0.00
CH3O2+C2H5<=>CH3O+C2H5O                  8.000e+12   -0.000  -1000.00
C2H5O+O2<=>CH3CHO+HO2                    4.280e+10   -0.000   1097.00
CH3+CH2O<=>C2H5O                         3.000e+11   -0.000   6336.00
CH3CHO+H<=>C2H5O                         4.610e+07    1.710   7090.00
C2H5O2+CH2O<=>C2H5O2H+HCO                1.990e+12   -0.000  11660.00
CH4+C2H5O2<=>CH3+C2H5O2H                 1.810e+11   -0.000  18480.00
CH3OH+C2H5O2<=>CH2OH+C2H5O2H             1.810e+12   -0.000  13710.00
C2H5O2+HO2<=>C2H5O2H+O2                  1.750e+10   -0.000  -3275.00
C2H6+C2H5O2<=>C2H5+C2H5O2H               8.600e+00    3.760  17200.00
C2H5O2H<=>C2H5O+OH                       6.310e+14   -0.000  42300.00
C2H5+O2<=>C2H5O2                         3.398e+53  -13.900   9279.00
PLOG /     0.040 3.398e+53  -13.900   9279.00 /
PLOG /     1.000 9.362e+59  -15.280  14240.00 /
PLOG /    10.000 1.262e+60  -14.910  16240.00 /
C2H5+O2<=>C2H4+HO2                       2.094e+09    0.490   -391.40
PLOG /     0.040 2.094e+09    0.490   -391.40 /
PLOG /     1.000 1.843e+07    1.130   -720.60 /
PLOG /    10.000 7.561e+14   -1.010   4749.00 /
     DUPLICATE
C2H5+O2<=>C2H4+HO2                       6.609e+00    3.510  14160.00
     DUPLICATE
C2H5+O2<=>C2H4O1+OH                      1.303e+03    1.930   -502.70
PLOG /     0.040 1.303e+03    1.930   -502.70 /
PLOG /     1.000 2.438e+02    2.180    -62.50 /
PLOG /    10.000 4.621e+09    0.150   5409.00 /
C2H5+O2<=>CH3CHO+OH                      4.908e-06    4.760    254.30
PLOG /     0.040 4.908e-06    4.760    254.30 /
PLOG /     1.000 6.803e-02    3.570   2643.00 /
PLOG /    10.000 8.265e+02    2.410   5285.00 /
C2H5O2<=>CH3CHO+OH                       1.237e+35   -9.420  36360.00
PLOG /     0.040 1.237e+35   -9.420  36360.00 /
PLOG /     1.000 1.687e+36   -9.220  38700.00 /
PLOG /    10.000 2.520e+41  -10.200  43710.00 /
C2H5O2<=>C2H4+HO2                        1.782e+32   -7.100  32840.00
PLOG /     0.040 1.782e+32   -7.100  32840.00 /
PLOG /     1.000 2.701e+37   -8.470  35840.00 /
PLOG /    10.000 1.980e+38   -8.460  37900.00 /
C2H5O2<=>C2H4O1+OH                       5.778e+45  -11.900   4112.00
PLOG /     0.040 5.778e+45  -11.900   4112.00 /
PLOG /     1.000 1.916e+43  -10.750  42400.00 /
PLOG /    10.000 3.965e+43  -10.460  45580.00 /
C2H4O1<=>CH3+HCO                         3.630e+13   -0.000  57200.00
C2H4O1<=>CH3CHO                          7.407e+12   -0.000  53800.00
CH3CHO(+M)<=>CH3+HCO(+M)                 2.450e+22   -1.740  86355.00

     LOW  /  1.030e+59  -11.300  95912.50 /
     TROE/  0.00249    718.1     6.089      3780 /
CH3CHO(+M)<=>CH4+CO(+M)                  2.720e+21   -1.740  86355.00

     LOW  /  1.144e+58  -11.300  95912.50 /
     TROE/  0.00249    718.1     6.089      3780 /
CH3CHO+H<=>CH3CO+H2                      1.310e+05    2.580   1220.00
CH3CHO+H<=>CH2CHO+H2                     2.720e+03    3.100   5210.00
CH3CHO+O<=>CH3CO+OH                      5.940e+12   -0.000   1868.00
CH3CHO+OH<=>CH3CO+H2O                    3.370e+12    0.000   -619.00
CH3CHO+O2<=>CH3CO+HO2                    3.010e+13   -0.000  39150.00
CH3CHO+CH3<=>CH3CO+CH4                   7.080e-04    4.580   1966.00
CH3CHO+HO2<=>CH3CO+H2O2                  3.010e+12   -0.000  11920.00
CH3O2+CH3CHO<=>CH3O2H+CH3CO              3.010e+12   -0.000  11920.00
CH3CHO+CH3CO3<=>CH3CO+CH3CO3H            3.010e+12   -0.000  11920.00
CH3CHO+OH<=>CH3+HOCHO                    3.000e+15   -1.076     -0.00
CH3CHO+OH<=>CH2CHO+H2O                   1.720e+05    2.400    815.00
CH3CO(+M)<=>CH3+CO(+M)                   1.070e+12    0.630  16900.00

     LOW  /  5.650e+18   -0.970  14600.00 /
     TROE/   0.629 8.73e+09      5.52   7.6e+07 /
CH3CO+H<=>CH2CO+H2                       2.000e+13   -0.000     -0.00
CH3CO+O<=>CH2CO+OH                       2.000e+13   -0.000     -0.00
CH3CO+CH3<=>CH2CO+CH4                    5.000e+13   -0.000     -0.00
CH3CO+O2<=>CH3CO3                        1.200e+11   -0.000  -1100.00
CH3CO3+HO2<=>CH3CO3H+O2                  1.750e+10   -0.000  -3275.00
H2O2+CH3CO3<=>HO2+CH3CO3H                2.410e+12   -0.000   9936.00
CH4+CH3CO3<=>CH3+CH3CO3H                 1.810e+11   -0.000  18480.00
CH2O+CH3CO3<=>HCO+CH3CO3H                1.990e+12   -0.000  11660.00
C2H6+CH3CO3<=>C2H5+CH3CO3H               1.700e+13   -0.000  20460.00
CH3CO3H<=>CH3CO2+OH                      5.010e+14   -0.000  40150.00
CH3CO2+M<=>CH3+CO2+M                     4.400e+15   -0.000  10500.00

CH2CHO(+M)<=>CH2CO+H(+M)                 1.430e+15   -0.150  45600.00

     LOW  /  6.000e+29   -3.800  43423.90 /
     TROE/   0.985      393   9.8e+09     5e+09 /
CH2CHO(+M)<=>CH3+CO(+M)                  2.930e+12    0.290  40300.00

     LOW  /  9.520e+33   -5.070  41300.00 /
     TROE/  7.13e-17     1150  4.99e+09  1.79e+09 /
CH2CHO+O2<=>CH2CO+HO2                    1.880e+05    2.370  23730.00
PLOG /     0.010 1.880e+05    2.370  23730.00 /
PLOG /     0.100 1.880e+05    2.370  27370.00 /
PLOG /     1.000 2.510e+05    2.330  23800.00 /
PLOG /    10.000 7.050e+07    1.630  25290.00 /
CH2CHO+O2=>CH2O+CO+OH                    2.680e+17   -1.840   6530.00
PLOG /     0.010 2.680e+17   -1.840   6530.00 /
PLOG /     0.100 1.520e+20   -2.580   8980.00 /
PLOG /     1.000 1.650e+19   -2.220  10340.00 /
PLOG /    10.000 8.953e+13   -0.600  10120.00 /
CH2+CO(+M)<=>CH2CO(+M)                   8.100e+11   -0.000     -0.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  2.690e+33   -5.110   7095.00 /
     TROE/  0.5907      275      1226      5185 /
CH3CO(+M)<=>CH2CO+H(+M)                  9.413e+07    1.917  44987.20

     LOW  /  1.516e+51  -10.270  55390.00 /
     TROE/  0.6009 8.103e+09     667.7     5e+09 /
CH2CO+H<=>HCCO+H2                        1.401e+15   -0.171   8783.20
CH2CO+H<=>CH3+CO                         7.704e+13   -0.171   4183.20
CH2CO+O<=>CH2+CO2                        1.750e+12    0.000   1350.00
CH2CO+O<=>HCCO+OH                        1.000e+13   -0.000   8000.00
CH2CO+OH<=>HCCO+H2O                      1.000e+13   -0.000   2000.00
CH2CO+OH<=>CH2OH+CO                      2.000e+12   -0.000  -1010.00
CH2CO+CH3<=>C2H5+CO                      4.769e+04    2.312   9468.00
CH2Y+CH2CO<=>C2H4+CO                     1.600e+14   -0.000     -0.00
HCCO+OH=>H2+2CO                          1.000e+14   -0.000     -0.00
HCCO+O=>H+2CO                            8.000e+13   -0.000     -0.00
HCCO+H<=>CH2Y+CO                         1.000e+14   -0.000     -0.00
HCCO+O2=>OH+2CO                          1.910e+11   -0.020   1020.00
HCCO+O2=>CO2+CO+H                        4.780e+12   -0.142   1150.00
CH+CO+M<=>HCCO+M                         7.570e+22   -1.900      0.00

CH+CH2O<=>H+CH2CO                        9.460e+13    0.000   -515.00
CH+HCCO<=>CO+C2H2                        5.000e+13   -0.000     -0.00
C2H3+H(+M)<=>C2H4(+M)                    6.080e+12    0.270    280.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  1.400e+30   -3.860   3320.00 /
     TROE/   0.782    207.5      2663      6095 /
C2H4(+M)<=>H2+H2CC(+M)                   8.000e+12    0.440  88770.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ AR/0.70/ 
     LOW  /  7.000e+50   -9.310  99860.00 /
     TROE/  0.7345      180      1035      5417 /
C2H4+H<=>C2H3+H2                         5.070e+07    1.930  12950.00
C2H4+O<=>CH3+HCO                         7.453e+06    1.880    183.00
C2H4+O<=>CH2CHO+H                        6.098e+06    1.880    183.00
C2H4+OH<=>C2H3+H2O                       2.230e+04    2.745   2215.50
C2H4+OH<=>CH3+CH2O                       5.350e+00    2.920  -1732.70
PLOG /     0.010 5.350e+00    2.920  -1732.70 /
PLOG /     0.025 3.190e+01    2.710  -1172.30 /
PLOG /     0.100 5.550e+02    2.360   -180.80 /
PLOG /     1.000 1.780e+05    1.680   2060.50 /
PLOG /    10.000 2.370e+09    0.560   6006.70 /
PLOG /   100.000 2.760e+13   -0.500  11455.10 /
C2H4+OH<=>CH3CHO+H                       2.370e-07    5.300  -2050.60
PLOG /     0.010 2.370e-07    5.300  -2050.60 /
PLOG /     0.025 8.730e-05    4.570   -618.00 /
PLOG /     0.100 4.030e-01    3.540   1881.70 /
PLOG /     1.000 2.380e-02    3.910   1722.70 /
PLOG /    10.000 8.250e+08    1.010  10507.30 /
PLOG /   100.000 6.800e+09    0.810  13867.30 /
C2H4+OH<=>C2H3OH+H                       1.040e+04    2.600   4121.00
PLOG /     0.010 1.040e+04    2.600   4121.00 /
PLOG /     0.025 1.070e+04    2.600   4129.00 /
PLOG /     0.100 1.520e+04    2.560   4238.30 /
PLOG /     1.000 3.190e+05    2.190   5255.60 /
PLOG /    10.000 1.940e+08    1.430   7828.80 /
PLOG /   100.000 8.550e+10    0.750  11490.80 /
C2H3OH+O2<=>CH2CHO+HO2                   5.310e+11    0.210  39830.00
C2H3OH+O<=>CH2CHO+OH                     1.875e+06    1.900   -860.00
C2H3OH+OH<=>CH2CHO+H2O                   3.330e+09    1.100    540.50
C2H3OH+CH3<=>CH2CHO+CH4                  2.030e-08    5.900   1052.00
C2H3OH+CH3O2<=>CH2CHO+CH3O2H             3.400e+03    2.500   8922.00
C2H3OH+H<=>CH2CHO+H2                     1.480e+03    3.077   7220.00
C2H3OH+HO2<=>CH3CHO+HO2                  1.490e+05    1.670   6810.00
C2H3OH<=>CH3CHO                          7.420e+46  -10.560  67420.00
PLOG /     0.100 7.420e+46  -10.560  67420.00 /
PLOG /     1.000 4.420e+42   -9.090  67069.20 /
PLOG /   100.000 2.900e+27   -4.350  61612.90 /
C2H4+CH3<=>C2H3+CH4                      6.620e+00    3.700   9500.00
C2H4+O2<=>C2H3+HO2                       4.220e+13   -0.000  57623.10
C2H4+CH3O<=>C2H3+CH3OH                   1.200e+11   -0.000   6750.00
C2H4+CH3O2<=>C2H3+CH3O2H                 8.590e+00    3.754  27132.00
C2H4+C2H5O2<=>C2H3+C2H5O2H               8.590e+00    3.754  27132.00
C2H4+CH3CO3<=>C2H3+CH3CO3H               1.130e+13   -0.000  30430.00
C2H4+CH3O2<=>C2H4O1+CH3O                 2.820e+12   -0.000  17110.00
C2H4+C2H5O2<=>C2H4O1+C2H5O               2.820e+12   -0.000  17110.00
C2H4+HO2<=>C2H4O1+OH                     5.575e+11   -0.000  17190.00
CH+CH4<=>C2H4+H                          6.000e+13   -0.000     -0.00
CH2Y+CH3<=>C2H4+H                        2.000e+13   -0.000     -0.00
C2H2+H(+M)<=>C2H3(+M)                    1.710e+10    1.266   2709.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  6.346e+31   -4.664   3780.00 /
!     FCCHECK/   0.212    -10200  0.788     1e-30      0         0
     TROE/   0.788   -10200     1e-30        /
C2H3+O2<=>CH2O+HCO                       1.700e+29   -5.312   6503.10
C2H3+O2<=>CH2CHO+O                       7.000e+14   -0.611   5262.40
C2H3+O2=>H+CO+CH2O                       5.190e+15   -1.260   3312.60
CH3+C2H3<=>CH4+C2H2                      3.920e+11   -0.000     -0.00
C2H3+H<=>C2H2+H2                         9.000e+13   -0.000     -0.00
C2H3+H<=>H2CC+H2                         6.000e+13   -0.000     -0.00
C2H3+OH<=>C2H2+H2O                       3.011e+13   -0.000     -0.00
2C2H3<=>C2H2+C2H4                        9.600e+11   -0.000     -0.00
C2H+H(+M)<=>C2H2(+M)                     1.000e+17   -0.000     -0.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ HE/0.70/ AR/0.70/ 
     LOW  /  3.750e+33   -4.800   1900.00 /
     TROE/   0.646      132      1315      5566 /
C2H+O<=>CH+CO                            5.000e+13   -0.000     -0.00
C2H+OH<=>H+HCCO                          2.000e+13   -0.000     -0.00
C2H+O2<=>HCO+CO                          5.000e+13   -0.000   1500.00
C2H+H2<=>H+C2H2                          4.900e+05    2.500    560.00
C2H2(+M)<=>H2CC(+M)                      8.000e+14   -0.520  50750.00
CH4/2.00/ H2/2.00/ CO/1.50/ CO2/2.00/ H2O/6.00/ C2H6/3.00/ C2H4/2.50/ 
C2H2/2.50/ 
     LOW  /  2.450e+15   -0.640  49700.00 /
!     FCCHECK/       0         0      0         0      1         0
     TROE/       1     1.00 10000000.00 10000000.00 /
C2H2+O<=>CH2+CO                          7.395e+08    1.280   2472.00
C2H2+O<=>HCCO+H                          2.958e+09    1.280   2472.00
C2H2+OH<=>C2H+H2O                        2.632e+06    2.140  17060.00
C2H2+OH<=>CH2CO+H                        1.578e+03    2.560   -844.50
PLOG /     0.010 1.578e+03    2.560   -844.50 /
PLOG /     0.025 1.518e+04    2.280   -292.10 /
PLOG /     0.100 3.017e+05    1.920    598.10 /
PLOG /     1.000 7.528e+06    1.550   2106.00 /
PLOG /    10.000 5.101e+06    1.650   3400.00 /
PLOG /   100.000 1.457e+04    2.450   4477.00 /
C2H2+OH<=>CH3+CO                         4.757e+05    1.680   -329.80
PLOG /     0.010 4.757e+05    1.680   -329.80 /
PLOG /     0.025 4.372e+06    1.400    226.50 /
PLOG /     0.100 7.648e+07    1.050   1115.00 /
PLOG /     1.000 1.277e+09    0.730   2579.00 /
PLOG /    10.000 4.312e+08    0.920   3736.00 /
PLOG /   100.000 8.250e+05    1.770   4697.00 /
C2H2+HCO<=>C2H3+CO                       1.000e+07    2.000   6000.00
H2CC+H<=>C2H2+H                          1.000e+14   -0.000     -0.00
H2CC+OH<=>CH2CO+H                        2.000e+13   -0.000     -0.00
H2CC+O2<=>2HCO                           1.000e+13   -0.000     -0.00
C2H5OH<=>C2H4+H2O                        3.410e+59  -14.200  83672.60
PLOG /     0.001 3.410e+59  -14.200  83672.60 /
PLOG /     0.010 2.620e+57  -13.300  85262.20 /
PLOG /     0.100 1.650e+52  -11.500  84745.60 /
PLOG /     1.000 5.230e+43   -8.900  81506.70 /
PLOG /    10.000 4.590e+32   -5.600  76062.40 /
PLOG /   100.000 3.840e+20   -2.060  69465.50 /
C2H5OH<=>CH3+CH2OH                       1.200e+54  -12.900 100005.70
PLOG /     0.001 1.200e+54  -12.900 100005.70 /
PLOG /     0.010 5.180e+59  -14.000  99906.40 /
PLOG /     0.100 1.620e+66  -15.300 105390.50 /
PLOG /     1.000 5.550e+64  -14.500 106183.00 /
PLOG /    10.000 1.550e+58  -12.300 105768.00 /
PLOG /   100.000 1.780e+47   -8.960 101058.80 /
C2H5OH<=>C2H5+OH                         8.100e+46  -11.300 111053.40
PLOG /     0.001 8.100e+46  -11.300 111053.40 /
PLOG /     0.010 1.860e+56  -13.500 107238.40 /
PLOG /     0.100 4.650e+63  -15.000 109622.80 /
PLOG /     1.000 4.460e+65  -14.900 112345.00 /
PLOG /    10.000 2.790e+61  -13.400 113080.20 /
PLOG /   100.000 6.170e+51  -10.300 109940.70 /
C2H5OH+O2<=>SC2H4OH+HO2                  1.500e+13   -0.000  50150.00
C2H5OH+H<=>SC2H4OH+H2                    8.790e+04    2.680   2910.00
C2H5OH+H<=>C2H5O+H2                      9.450e+02    3.140   8701.10
C2H5OH+OH<=>SC2H4OH+H2O                  7.170e+04    2.540  -1534.00
C2H5OH+OH<=>C2H5O+H2O                    5.810e-03    4.280  -3560.00
C2H5OH+HO2<=>SC2H4OH+H2O2                3.500e-05    5.260   7475.10
C2H5OH+HO2<=>C2H5O+H2O2                  6.470e-07    5.300  10533.10
C2H5OH+CH3O2<=>SC2H4OH+CH3O2H            8.200e+03    2.550  10750.00
C2H5OH+CH3O2<=>C2H5O+CH3O2H              2.500e+12   -0.000  24000.00
C2H5OH+O<=>SC2H4OH+OH                    1.450e+05    2.470    876.00
C2H5OH+O<=>C2H5O+OH                      1.460e-03    4.730   1727.00
C2H5OH+CH3<=>SC2H4OH+CH4                 1.993e+01    3.370   7634.00
C2H5OH+CH3<=>C2H5O+CH4                   2.035e+00    3.570   7721.00
C2H5OH+C2H5<=>SC2H4OH+C2H6               5.000e+10   -0.000  10400.00
SC2H4OH<=>CH3CHO+H                       5.690e+52  -13.380  45049.00
PLOG /     0.001 5.690e+52  -13.380  45049.00 /
PLOG /     0.010 3.290e+56  -14.120  48129.00 /
PLOG /     0.100 8.580e+57  -14.160  50743.00 /
PLOG /     1.000 5.360e+55  -13.150  51886.00 /
PLOG /    10.000 1.660e+48  -10.640  50297.00 /
PLOG /    20.000 8.260e+44   -9.590  49218.00 /
PLOG /    50.000 1.010e+40   -8.060  47439.00 /
PLOG /   100.000 1.100e+36   -6.840  45899.00 /
SC2H4OH<=>C2H3OH+H                       5.400e+46  -11.630  44323.00
PLOG /     0.001 5.400e+46  -11.630  44323.00 /
PLOG /     0.010 1.210e+51  -12.550  47240.00 /
PLOG /     0.100 2.870e+54  -13.150  50702.00 /
PLOG /     1.000 3.790e+53  -12.510  52560.00 /
PLOG /    10.000 6.330e+46  -10.200  51441.00 /
PLOG /    20.000 3.870e+43   -9.170  50440.00 /
PLOG /    50.000 5.080e+38   -7.650  48713.00 /
PLOG /   100.000 5.120e+34   -6.410  47182.00 /
SC2H4OH<=>C2H5O                          5.480e+45  -11.630  44328.00
PLOG /     0.001 5.480e+45  -11.630  44328.00 /
PLOG /     0.010 2.540e+49  -12.370  46445.00 /
PLOG /     0.100 1.650e+54  -13.400  50330.00 /
PLOG /     1.000 1.810e+55  -13.310  53132.00 /
PLOG /    10.000 4.580e+49  -11.320  52714.00 /
PLOG /    20.000 4.110e+46  -10.330  51834.00 /
PLOG /    50.000 6.680e+41   -8.830  50202.00 /
PLOG /   100.000 6.540e+37   -7.580  48697.00 /
SC2H4OH+O2<=>CH3CHO+HO2                  5.260e+17   -1.637    838.00
PLOG /     0.010 5.260e+17   -1.637    838.00 /
PLOG /     0.100 5.260e+17   -1.637    838.00 /
PLOG /     1.000 5.280e+17   -1.638    839.00 /
PLOG /    10.000 1.540e+18   -1.771   1120.00 /
PLOG /   100.000 3.780e+20   -2.429   3090.00 /
SC2H4OH+O2<=>C2H3OH+HO2                  5.120e+02    2.496   -414.00
PLOG /     0.010 5.120e+02    2.496   -414.00 /
PLOG /     0.100 5.330e+02    2.490   -402.00 /
PLOG /     1.000 7.620e+02    2.446   -296.00 /
PLOG /    10.000 8.920e+03    2.146    470.00 /
PLOG /   100.000 4.380e+05    1.699   2330.00 /
CH3COCH3<=>CH3CO+CH3                     2.050e+58  -12.796 100030.10
PLOG /     0.010 2.050e+58  -12.796 100030.10 /
PLOG /     0.100 3.300e+51  -10.574  98221.20 /
PLOG /     1.000 1.310e+42   -7.657  94660.60 /
PLOG /    10.000 2.160e+33   -4.989  90916.50 /
PLOG /   100.000 9.400e+28   -3.669  89022.80 /
CH3COCH3+OH<=>CH3COCH2+H2O               1.250e+05    2.483    445.00
CH3COCH3+H<=>CH3COCH2+H2                 9.800e+05    2.430   5160.00
CH3COCH3+O<=>CH3COCH2+OH                 5.130e+11    0.211   4890.00
CH3COCH3+CH3<=>CH3COCH2+CH4              3.960e+11   -0.000   9784.00
CH3COCH3+CH3O<=>CH3COCH2+CH3OH           4.340e+11   -0.000   6460.00
CH3COCH3+O2<=>CH3COCH2+HO2               6.030e+13   -0.000  48500.00
CH3COCH3+HO2<=>CH3COCH2+H2O2             1.700e+13   -0.000  20460.00
CH3COCH3+CH3O2<=>CH3COCH2+CH3O2H         1.700e+13   -0.000  20460.00
CH2CO+CH3<=>CH3COCH2                     1.760e+04    2.480   6130.00
C2H3+HCO<=>C2H3CHO                       1.810e+13   -0.000     -0.00
C2H3CHO+H<=>C2H3CO+H2                    1.340e+13   -0.000   3300.00
C2H3CHO+O<=>C2H3CO+OH                    5.940e+12   -0.000   1868.00
C2H3CHO+OH<=>C2H3CO+H2O                  9.240e+06    1.500   -962.00
C2H3CHO+O2<=>C2H3CO+HO2                  1.005e+13   -0.000  40700.00
C2H3CHO+HO2<=>C2H3CO+H2O2                3.010e+12   -0.000  11920.00
C2H3CHO+CH3<=>C2H3CO+CH4                 2.608e+06    1.780   5911.00
C2H3CHO+C2H3<=>C2H3CO+C2H4               1.740e+12   -0.000   8440.00
C2H3CHO+CH3O<=>C2H3CO+CH3OH              1.000e+12   -0.000   3300.00
C2H3CHO+CH3O2<=>C2H3CO+CH3O2H            3.010e+12   -0.000  11920.00
C2H3+CO<=>C2H3CO                         1.510e+11   -0.000   4810.00
C2H5+HCO<=>C2H5CHO                       1.810e+13   -0.000     -0.00
C2H5CHO+H<=>C2H5CO+H2                    4.000e+13   -0.000   4200.00
C2H5CHO+O<=>C2H5CO+OH                    5.000e+12   -0.000   1790.00
C2H5CHO+OH<=>C2H5CO+H2O                  2.690e+10    0.760   -340.00
C2H5CHO+CH3<=>C2H5CO+CH4                 2.608e+06    1.780   5911.00
C2H5CHO+HO2<=>C2H5CO+H2O2                2.800e+12   -0.000  13600.00
C2H5CHO+CH3O<=>C2H5CO+CH3OH              1.000e+12   -0.000   3300.00
C2H5CHO+CH3O2<=>C2H5CO+CH3O2H            3.010e+12   -0.000  11920.00
C2H5CHO+C2H5<=>C2H5CO+C2H6               1.000e+12   -0.000   8000.00
C2H5CHO+C2H5O<=>C2H5CO+C2H5OH            6.026e+11   -0.000   3300.00
C2H5CHO+C2H5O2<=>C2H5CO+C2H5O2H          3.010e+12   -0.000  11920.00
C2H5CHO+O2<=>C2H5CO+HO2                  1.005e+13   -0.000  40700.00
C2H5CHO+CH3CO3<=>C2H5CO+CH3CO3H          3.010e+12   -0.000  11920.00
C2H5CHO+C2H3<=>C2H5CO+C2H4               1.700e+12   -0.000   8440.00
C2H5+CO<=>C2H5CO                         1.510e+11   -0.000   4810.00
CH3OCH3(+M)<=>CH3+CH3O(+M)               4.380e+21   -1.570  83890.00

     LOW  /  7.520e+15   -0.000  42790.00 /
!     FCCHECK/   0.546     1e-30  0.454      2510      0         0
     TROE/   0.454    1e-30      2510        /
CH3OCH3+OH<=>CH3OCH2+H2O                 6.324e+06    2.000   -651.70
CH3OCH3+H<=>CH3OCH2+H2                   7.721e+06    2.090   3384.00
CH3OCH3+O<=>CH3OCH2+OH                   7.750e+08    1.360   2250.00
CH3OCH3+HO2<=>CH3OCH2+H2O2               8.670e+02    3.010  12090.00
CH3OCH3+CH3O2<=>CH3OCH2+CH3O2H           3.120e+02    3.120  13190.00
CH3OCH3+CH3<=>CH3OCH2+CH4                1.445e-06    5.730   5700.00
CH3OCH3+O2<=>CH3OCH2+HO2                 4.100e+13   -0.000  44910.00
CH3OCH3+CH3O<=>CH3OCH2+CH3OH             6.020e+11   -0.000   4074.00
CH3OCH3+CH3OCH2O2<=>CH3OCH2+CH3OCH2O2H   5.000e+12   -0.000  17690.00
CH3OCH3+O2CHO<=>CH3OCH2+HO2CHO           4.425e+04    2.600  13910.00
CH3OCH3+OCHO<=>CH3OCH2+HOCHO             1.000e+13   -0.000  17690.00
CH3OCH2<=>CH2O+CH3                       1.600e+13   -0.000  25500.00
CH3OCH2+CH3O<=>CH3OCH3+CH2O              2.410e+13   -0.000     -0.00
CH3OCH2+CH2O<=>CH3OCH3+HCO               5.490e+03    2.800   5862.00
CH3OCH2+CH3CHO<=>CH3OCH3+CH3CO           1.260e+12    0.000   8499.00
CH3OCH2+O2<=>CH3OCH2O2                   2.000e+12   -0.000     -0.00
CH3OCH2O2+CH2O<=>CH3OCH2O2H+HCO          1.000e+12   -0.000  11660.00
CH3OCH2O2+CH3CHO<=>CH3OCH2O2H+CH3CO      2.800e+12   -0.000  13600.00
2CH3OCH2O2=>O2+2CH3OCH2O                 1.547e+23   -4.500     -0.00
CH3OCH2O+OH<=>CH3OCH2O2H                 2.000e+13   -0.000     -0.00
CH3O+CH2O<=>CH3OCH2O                     1.000e+11   -0.000   7960.00
CH3OCH2O+O2<=>CH3OCHO+HO2                5.000e+10   -0.000    500.00
CH3OCHO+H<=>CH3OCH2O                     1.000e+13   -0.000   7838.00
CH3OCH2O2<=>CH2OCH2O2H                   6.000e+10   -0.000  21580.00
CH2OCH2O2H=>OH+2CH2O                     1.500e+13   -0.000  20760.00
CH2OCH2O2H+O2<=>O2CH2OCH2O2H             7.000e+11   -0.000     -0.00
O2CH2OCH2O2H<=>HO2CH2OCHO+OH             4.000e+10    0.000  18580.00
HO2CH2OCHO<=>OCH2OCHO+OH                 2.000e+16   -0.000  40500.00
CH2O+OCHO<=>OCH2OCHO                     1.250e+11   -0.000  11900.00
OCH2OCHO<=>HOCH2OCO                      1.000e+11   -0.000  14000.00
HOCH2O+CO<=>HOCH2OCO                     1.500e+11   -0.000   4800.00
CH2OH+CO2<=>HOCH2OCO                     1.500e+11   -0.000  35720.00
CH2OCHO+H<=>CH3OCHO                      1.000e+14   -0.000     -0.00
CH3OCO+H<=>CH3OCHO                       1.000e+14   -0.000     -0.00
CH3OCHO(+M)<=>CH3OH+CO(+M)               1.000e+14   -0.000  62500.00

     LOW  /  6.143e+60  -12.070  75400.00 /
     TROE/    0.78 8.28e+09     438.9   6.7e+08 /
CH3O+HCO<=>CH3OCHO                       3.000e+13   -0.000     -0.00
CH3+OCHO<=>CH3OCHO                       1.000e+13   -0.000     -0.00
CH3OCHO+O2<=>CH3OCO+HO2                  1.000e+13   -0.000  49700.00
CH3OCHO+O2<=>CH2OCHO+HO2                 2.050e+13   -0.000  52000.00
CH3OCHO+OH<=>CH3OCO+H2O                  1.580e+07    1.800    934.00
CH3OCHO+OH<=>CH2OCHO+H2O                 5.270e+09    0.970   1586.00
CH3OCHO+HO2<=>CH3OCO+H2O2                4.820e+03    2.600  13910.00
CH3OCHO+HO2<=>CH2OCHO+H2O2               2.380e+04    2.550  16490.00
CH3OCHO+O<=>CH3OCO+OH                    2.755e+05    2.450   2830.00
CH3OCHO+O<=>CH2OCHO+OH                   9.800e+05    2.430   4750.00
CH3OCHO+H<=>CH3OCO+H2                    6.500e+05    2.400   4471.00
CH3OCHO+H<=>CH2OCHO+H2                   6.650e+05    2.540   6756.00
CH3OCHO+CH3<=>CH3OCO+CH4                 7.550e-01    3.460   5481.00
CH3OCHO+CH3<=>CH2OCHO+CH4                4.520e-01    3.650   7154.00
CH3OCHO+CH3O<=>CH3OCO+CH3OH              5.480e+11   -0.000   5000.00
CH3OCHO+CH3O<=>CH2OCHO+CH3OH             2.170e+11   -0.000   6458.00
CH3OCHO+CH3O2<=>CH3OCO+CH3O2H            4.820e+03    2.600  13910.00
CH3OCHO+CH3O2<=>CH2OCHO+CH3O2H           2.380e+04    2.550  16490.00
CH3OCHO+HCO<=>CH3OCO+CH2O                5.400e+06    1.900  17010.00
CH3OCHO+HCO<=>CH2OCHO+CH2O               1.025e+05    2.500  18430.00
CH3OCO<=>CH2OCHO                         1.629e+12   -0.180  40670.00
CH3+CO2<=>CH3OCO                         4.760e+07    1.540  34700.00
CH3O+CO<=>CH3OCO                         1.550e+06    2.020   5730.00
CH2O+HCO<=>CH2OCHO                       1.500e+11   -0.000  11900.00
END
)MECHANISM";
}

const char* AramcoMech_DMEonly_74spec::GetThermoResource() const {
  return R"THERMORES(THERMO
   300.000  1000.000  5000.000
N2                000000N   2               G    300.00   5000.00 1000.00      1
 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00                   4
CH3OCH3           000000H   6O   1C   2     G    300.00   5000.00 1000.00      1
 7.36358679E+00 1.38910153E-02-4.74083257E-06 7.34874498E-10-4.25867578E-14    2
-2.61148074E+04-1.61332876E+01 2.41860337E+00 1.87279640E-02-1.40894592E-06    3
-4.28367822E-09 1.36818123E-12-2.36771735E+04 1.28877967E+01                   4
CH4               000000H   4C   1          G    300.00   5000.00 1000.00      1
 1.65326226E+00 1.00263099E-02-3.31661238E-06 5.36483138E-10-3.14696758E-14    2
-1.00095936E+04 9.90506283E+00 5.14911468E+00-1.36622009E-02 4.91453921E-05    3
-4.84246767E-08 1.66603441E-11-1.02465983E+04-4.63848842E+00                   4
H2                000000H   2               G    300.00   5000.00 1000.00      1
 2.93286575E+00 8.26608026E-04-1.46402364E-07 1.54100414E-11-6.88804800E-16    2
-8.13065581E+02-1.02432865E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
O2                000000O   2               G    300.00   5000.00 1000.00      1
 3.66096065E+00 6.56365811E-04-1.41149627E-07 2.05797935E-11-1.29913436E-15    2
-1.21597718E+03 3.41536279E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
CO                000000O   1C   1          G    300.00   5000.00 1000.00      1
 3.04848590E+00 1.35172810E-03-4.85794050E-07 7.88536440E-11-4.69807460E-15    2
-1.42661170E+04 6.01709770E+00 3.57953350E+00-6.10353690E-04 1.01681430E-06    3
 9.07005860E-10-9.04424490E-13-1.43440860E+04 3.50840930E+00                   4
CO2               000000O   2C   1          G    300.00   5000.00 1000.00      1
 4.63651110E+00 2.74145690E-03-9.95897590E-07 1.60386660E-10-9.16198570E-15    2
-4.90249040E+04-1.93489550E+00 2.35681300E+00 8.98412990E-03-7.12206320E-06    3
 2.45730080E-09-1.42885480E-13-4.83719710E+04 9.90090350E+00                   4
H2O               000000H   2O   1          G    300.00   5000.00 1000.00      1
 2.67703890E+00 2.97318160E-03-7.73768890E-07 9.44335140E-11-4.26899910E-15    2
-2.98858940E+04 6.88255000E+00 4.19863520E+00-2.03640170E-03 6.52034160E-06    3
-5.48792690E-09 1.77196800E-12-3.02937260E+04-8.49009010E-01                   4
H                 000000H   1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54736600E+04-4.46682850E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54736600E+04-4.46682850E-01                   4
O                 000000O   1               G    300.00   5000.00 1000.00      1
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
OH                000000H   1O   1          G    300.00   5000.00 1000.00      1
 2.83853033E+00 1.10741289E-03-2.94000209E-07 4.20698729E-11-2.42289890E-15    2
 3.69780808E+03 5.84494652E+00 3.99198424E+00-2.40106655E-03 4.61664033E-06    3
-3.87916306E-09 1.36319502E-12 3.36889836E+03-1.03998477E-01                   4
CH3               000000H   3C   1          G    300.00   5000.00 1000.00      1
 2.97812060E+00 5.79785200E-03-1.97558000E-06 3.07297900E-10-1.79174160E-14    2
 1.65095130E+04 4.72247990E+00 3.65717970E+00 2.12659790E-03 5.45838830E-06    3
-6.61810030E-09 2.46570740E-12 1.64227160E+04 1.67353540E+00                   4
CH2O              000000H   2O   1C   1     G    300.00   5000.00 1000.00      1
 3.16952665E+00 6.19320560E-03-2.25056366E-06 3.65975660E-10-2.20149458E-14    2
-1.45486831E+04 6.04207898E+00 4.79372312E+00-9.90833322E-03 3.73219990E-05    3
-3.79285237E-08 1.31772641E-11-1.43791953E+04 6.02798058E-01                   4
C2H6              000000H   6C   2          G    300.00   5000.00 1000.00      1
 4.04666411E+00 1.53538802E-02-5.47039485E-06 8.77826544E-10-5.23167531E-14    2
-1.24473499E+04-9.68698313E-01 4.29142572E+00-5.50154901E-03 5.99438458E-05    3
-7.08466469E-08 2.68685836E-11-1.15222056E+04 2.66678994E+00                   4
CH3OCH2           000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 6.73134670E+00 1.18022530E-02-4.00104434E-06 6.17227325E-10-3.56435023E-14    2
-3.34933412E+03-9.40782513E+00 2.19172957E+00 1.90390966E-02-6.81450700E-06    3
 8.63577609E-11 3.23902492E-13-1.38948253E+03 1.62440756E+01                   4
O2CH2OCH2O2H      000000H   5O   5C   2     G    300.00   5000.00 1000.00      1
 1.96207800E+01 9.96131702E-03-3.41480475E-06 5.31364852E-10-3.08902100E-14    2
-3.98679641E+04-6.72721405E+01-1.22579809E-01 6.64376755E-02-6.55340609E-05    3
 3.13441786E-08-5.79042960E-12-3.41530358E+04 3.49112501E+01                   4
CH2OCH2O2H        000000H   5O   3C   2     G    300.00   5000.00 1000.00      1
 1.51031233E+01 9.18263324E-03-3.15895361E-06 4.92789411E-10-2.87001298E-14    2
-1.96611854E+04-4.85172618E+01 9.08332173E-01 4.47997444E-02-3.69441578E-05    3
 1.48366604E-08-2.32279459E-12-1.50992319E+04 2.66419885E+01                   4
HO2               000000H   1O   2          G    300.00   5000.00 1000.00      1
 4.17228741E+00 1.88117627E-03-3.46277286E-07 1.94657549E-11 1.76256905E-16    2
 3.10206839E+01 2.95767672E+00 4.30179807E+00-4.74912097E-03 2.11582905E-05    3
-2.42763914E-08 9.29225225E-12 2.64018485E+02 3.71666220E+00                   4
H2O2              000000H   2O   2          G    300.00   5000.00 1000.00      1
 4.57977305E+00 4.05326003E-03-1.29844730E-06 1.98211400E-10-1.13968792E-14    2
-1.80071775E+04 6.64970694E-01 4.31515149E+00-8.47390622E-04 1.76404323E-05    3
-2.26762944E-08 9.08950158E-12-1.77067437E+04 3.27373319E+00                   4
CH3OCH2O2         000000H   5O   3C   2     G    300.00   5000.00 1000.00      1
 1.22024149E+01 1.19225104E-02-4.06624713E-06 6.30079576E-10-3.65065312E-14    2
-2.40327175E+04-3.38665462E+01 1.72101041E+00 3.78155492E-02-2.89668870E-05    3
 1.16796117E-08-1.93234479E-12-2.05308597E+04 2.19296296E+01                   4
CH3O              000000H   3O   1C   1     G    300.00   5000.00 1000.00      1
 4.75779238E+00 7.44142474E-03-2.69705176E-06 4.38090504E-10-2.63537098E-14    2
 3.78111940E+02-1.96680028E+00 3.71180502E+00-2.80463306E-03 3.76550971E-05    3
-4.73072089E-08 1.86588420E-11 1.29569760E+03 6.57240864E+00                   4
CH3OH             000000H   4O   1C   1     G    300.00   5000.00 1000.00      1
 3.52726795E+00 1.03178783E-02-3.62892944E-06 5.77448016E-10-3.42182632E-14    2
-2.60028834E+04 5.16758693E+00 5.65851051E+00-1.62983419E-02 6.91938156E-05    3
-7.58372926E-08 2.80427550E-11-2.56119736E+04-8.97330508E-01                   4
CH3OCH2O          000000H   5O   2C   2     G    300.00   5000.00 1000.00      1
 8.60261845E+00 1.35772195E-02-4.84661602E-06 7.77766193E-10-4.62633624E-14    2
-2.13762444E+04-1.75775023E+01 3.25889339E+00 2.22146359E-02-7.78556340E-06    3
-2.41484158E-10 4.51914496E-13-1.92377212E+04 1.23680069E+01                   4
OCH2OCHO          000000H   3O   3C   2     G    300.00   5000.00 1000.00      1
 1.24378042E+01 7.82259106E-03-2.82874305E-06 4.55808935E-10-2.71391601E-14    2
-4.44526361E+04-3.85337353E+01 3.92797655E+00 1.85731220E-02-4.14656828E-07    3
-6.94060480E-09 2.32311215E-12-4.07070844E+04 1.02048630E+01                   4
HO2CH2OCHO        000000H   4O   4C   2     G    300.00   5000.00 1000.00      1
 1.78528337E+01 7.54042388E-03-2.74895624E-06 4.45343367E-10-2.66154007E-14    2
-6.41380252E+04-6.47357373E+01 5.15432996E-01 4.76081860E-02-3.63877925E-05    3
 1.24225089E-08-1.50813145E-12-5.82730214E+04 2.81992694E+01                   4
HCO               000000H   1O   1C   1     G    300.00   5000.00 1000.00      1
 3.92001542E+00 2.52279324E-03-6.71004164E-07 1.05615948E-10-7.43798261E-15    2
 3.65342928E+03 3.58077056E+00 4.23754610E+00-3.32075257E-03 1.40030264E-05    3
-1.34239995E-08 4.37416208E-12 3.87241185E+03 3.30834869E+00                   4
O2CHO             000000H   1O   3C   1     G    300.00   5000.00 1000.00      1
 7.24075139E+00 4.63312951E-03-1.63693995E-06 2.59706693E-10-1.52964699E-14    2
-1.87027618E+04-6.49547212E+00 3.96059309E+00 1.06002279E-02-5.25713351E-06    3
 1.01716726E-09-2.87487602E-14-1.73599383E+04 1.17807483E+01                   4
HOCHO             000000H   2O   2C   1     G    300.00   5000.00 1000.00      1
 4.61383160E+00 6.44963640E-03-2.29082510E-06 3.67160470E-10-2.18736750E-14    2
-4.75148500E+04 8.47883830E-01 3.89836160E+00-3.55877950E-03 3.55205380E-05    3
-4.38499590E-08 1.71077690E-11-4.67706090E+04 7.34953970E+00                   4
HOCH2O            000000H   3O   2C   1     G    300.00   5000.00 1000.00      1
 6.39521515E+00 7.43673043E-03-2.50422354E-06 3.84879712E-10-2.21778689E-14    2
-2.41108840E+04-6.63865583E+00 4.11183145E+00 7.53850697E-03 3.77337370E-06    3
-5.38746005E-09 1.45615887E-12-2.28023001E+04 7.46807254E+00                   4
HOCH2OCO          000000H   3O   3C   2     G    300.00   5000.00 1000.00      1
 1.13737391E+01 8.17663898E-03-2.92034021E-06 4.66695616E-10-2.76276823E-14    2
-4.65575743E+04-2.86035265E+01 6.08180801E+00 1.28768359E-02 2.04419418E-06    3
-6.10154921E-09 1.79820559E-12-4.39526183E+04 2.54054449E+00                   4
CH2OH             000000H   3O   1C   1     G    300.00   5000.00 1000.00      1
 5.09314370E+00 5.94761260E-03-2.06497460E-06 3.23008173E-10-1.88125902E-14    2
-4.03409640E+03-1.84691493E+00 4.47834367E+00-1.35070310E-03 2.78484980E-05    3
-3.64869060E-08 1.47907450E-11-3.50072890E+03 3.30913500E+00                   4
CH3O2H            000000H   4O   2C   1     G    300.00   5000.00 1000.00      1
 7.76538058E+00 8.61499712E-03-2.98006935E-06 4.68638071E-10-2.75339255E-14    2
-1.82979984E+04-1.43992663E+01 2.90540897E+00 1.74994735E-02 5.28243630E-06    3
-2.52827275E-08 1.34368212E-11-1.68894632E+04 1.13741987E+01                   4
CH3O2             000000H   3O   2C   1     G    300.00   5000.00 1000.00      1
 6.47970487E+00 7.44401080E-03-2.52348555E-06 3.89577296E-10-2.25182399E-14    2
-1.56285441E+03-8.19477074E+00 1.97339205E+00 1.53542340E-02-6.37314891E-06    3
 3.19930565E-10 2.82193915E-13 2.54278835E+02 1.69194215E+01                   4
HO2CHO            000000H   2O   3C   1     G    300.00   5000.00 1000.00      1
 9.87503878E+00 4.64663708E-03-1.67230522E-06 2.68624413E-10-1.59595232E-14    2
-3.80502496E+04-2.24939155E+01 2.42464726E+00 2.19706380E-02-1.68705546E-05    3
 6.25612194E-09-9.11645843E-13-3.54828006E+04 1.75027796E+01                   4
CH3OCH2O2H        000000H   6O   3C   2     G    300.00   5000.00 1000.00      1
 1.46617126E+01 1.20544674E-02-4.13604704E-06 6.43721637E-10-3.74188832E-14    2
-4.18896186E+04-4.81293252E+01 7.37369081E-01 4.67410003E-02-3.73051936E-05    3
 1.50735309E-08-2.43781978E-12-3.73260287E+04 2.57932851E+01                   4
CH2               000000H   2C   1          G    300.00   5000.00 1000.00      1
 3.14631886E+00 3.03671259E-03-9.96474439E-07 1.50483580E-10-8.57335515E-15    2
 4.60412605E+04 4.72341711E+00 3.71757846E+00 1.27391260E-03 2.17347251E-06    3
-3.48858500E-09 1.65208866E-12 4.58723866E+04 1.75297945E+00                   4
CH                000000H   1C   1          G    300.00   5000.00 1000.00      1
 2.52093690E+00 1.76536390E-03-4.61476600E-07 5.92896750E-11-3.34745010E-15    2
 7.09467690E+04 7.40518290E+00 3.48975830E+00 3.24321600E-04-1.68997510E-06    3
 3.16284200E-09-1.40618030E-12 7.06126460E+04 2.08428410E+00                   4
CH3OCHO           000000H   4O   2C   2     G    300.00   5000.00 1000.00      1
 6.33360880E+00 1.34851485E-02-4.84305805E-06 7.81719241E-10-4.67917447E-14    2
-4.68316521E+04-6.91542601E+00 5.96757028E+00-9.38085425E-03 7.07648417E-05    3
-8.29932227E-08 3.13522917E-11-4.55713267E+04 7.50341113E-01                   4
CH2OCHO           000000H   3O   2C   2     G    300.00   5000.00 1000.00      1
 1.00960096E+01 7.19887066E-03-2.59813465E-06 4.18110812E-10-2.48723387E-14    2
-2.36389018E+04-2.71144175E+01 2.31031671E+00 1.80474065E-02-2.71519637E-06    3
-4.60918579E-09 1.70037078E-12-2.02910878E+04 1.71549722E+01                   4
C2H4O1            000000H   4O   1C   2     G    300.00   5000.00 1000.00      1
 5.48876410E+00 1.20461900E-02-4.33369310E-06 7.00283110E-10-4.19490880E-14    2
-9.18042510E+03-7.07996050E+00 3.75905320E+00-9.44121800E-03 8.03097210E-05    3
-1.00807880E-07 4.00399210E-11-7.56081430E+03 7.84974750E+00                   4
C2H4              000000H   4C   2          G    300.00   5000.00 1000.00      1
 3.99182724E+00 1.04833908E-02-3.71721342E-06 5.94628366E-10-3.53630386E-14    2
 4.26865851E+03-2.69081762E-01 3.95920063E+00-7.57051373E-03 5.70989993E-05    3
-6.91588352E-08 2.69884190E-11 5.08977598E+03 4.09730213E+00                   4
C2H3              000000H   3C   2          G    300.00   5000.00 1000.00      1
 4.15026763E+00 7.54021341E-03-2.62997847E-06 4.15974048E-10-2.45407509E-14    2
 3.38566380E+04 1.72812235E+00 3.36377642E+00 2.65765722E-04 2.79620704E-05    3
-3.72986942E-08 1.51590176E-11 3.44749589E+04 7.91510092E+00                   4
C2H5              000000H   5C   2          G    300.00   5000.00 1000.00      1
 5.88784390E+00 1.03076793E-02-3.46844396E-06 5.32499257E-10-3.06512651E-14    2
 1.15065499E+04-8.49651771E+00 1.32730217E+00 1.76656753E-02-6.14926558E-06    3
-3.01143466E-10 4.38617775E-13 1.34284028E+04 1.71789216E+01                   4
OHY               000000H   1O   1          G    300.00   5000.00 1000.00      1
 2.88273000E+00 1.01397430E-03-2.27687700E-07 2.17468300E-11-5.12630500E-16    2
 5.02650000E+04 5.59571200E+00 3.63726600E+00 1.85091000E-04-1.67616460E-06    3
 2.38720200E-09-8.43144200E-13 5.00213000E+04 1.35886050E+00                   4
CH3OCO            000000H   3O   2C   2     G    300.00   5000.00 1000.00      1
 7.00171955E+00 1.01977290E-02-3.65621800E-06 5.89475086E-10-3.52561321E-14    2
-2.26135780E+04-9.05267669E+00 4.75563598E+00 7.80915313E-03 1.62272935E-05    3
-2.41210787E-08 9.42644561E-12-2.15157456E+04 4.78096491E+00                   4
HOCH2O2H          000000H   4O   3C   1     G    300.00   5000.00 1000.00      1
 1.16303827E+01 7.15133688E-03-2.39035030E-06 3.65772791E-10-2.10199524E-14    2
-4.31079242E+04-3.24276725E+01 1.85716693E+00 3.23153132E-02-2.69928902E-05    3
 1.11694484E-08-1.81284103E-12-4.00314471E+04 1.90917729E+01                   4
HOCH2O2           000000H   3O   3C   1     G    300.00   5000.00 1000.00      1
 9.04545938E+00 7.15223373E-03-2.37005676E-06 3.60083481E-10-2.05750228E-14    2
-2.49414886E+04-1.74210530E+01 2.85441621E+00 2.33663535E-02-1.88115990E-05    3
 7.96709515E-09-1.36346618E-12-2.29866196E+04 1.51730565E+01                   4
C2H3CO            000000H   3O   1C   3     G    300.00   5000.00 1000.00      1
 9.37467676E+00 7.91296900E-03-2.67198280E-06 4.11115430E-10-2.36978981E-14    2
 1.92969514E+03-2.40892696E+01 1.36242013E+00 3.15273972E-02-3.00218935E-05    3
 1.48167112E-08-2.87971530E-12 4.25770215E+03 1.72626546E+01                   4
C2H3CHO           000000H   4O   1C   3     G    300.00   5000.00 1000.00      1
 1.04184959E+01 9.48963321E-03-3.29310529E-06 5.16279203E-10-3.01587291E-14    2
-1.49630281E+04-3.07235061E+01 2.92355162E-01 3.54321417E-02-2.94936324E-05    3
 1.28100124E-08-2.26144108E-12-1.16521584E+04 2.28878280E+01                   4
CH3CO             000000H   3O   1C   2     G    300.00   5000.00 1000.00      1
 5.31371650E+00 9.17377930E-03-3.32203860E-06 5.39474560E-10-3.24523680E-14    2
-3.64504140E+03-1.67575580E+00 4.03587050E+00 8.77294870E-04 3.07100100E-05    3
-3.92475650E-08 1.52968690E-11-2.68207380E+03 7.86176820E+00                   4
C2H5O2H           000000H   6O   2C   2     G    300.00   5000.00 1000.00      1
 1.12305737E+01 1.20482120E-02-3.96730201E-06 6.00754632E-10-3.42657803E-14    2
-2.47977531E+04-3.25607232E+01 1.57329011E+00 3.52379996E-02-2.53203993E-05    3
 9.56802476E-09-1.48167375E-12-2.15278368E+04 1.90472032E+01                   4
C2H5O2            000000H   5O   2C   2     G    300.00   5000.00 1000.00      1
 8.88872432E+00 1.35833179E-02-4.91116949E-06 7.92343362E-10-4.73525704E-14    2
-7.44107388E+03-1.90789836E+01 4.50099327E+00 6.87965342E-03 4.74143971E-05    3
-6.92287127E-08 2.87395324E-11-5.39547911E+03 7.91490068E+00                   4
CH2CO             000000H   2O   1C   2     G    300.00   5000.00 1000.00      1
 5.35869367E+00 6.95641586E-03-2.64802637E-06 4.65067592E-10-3.08641820E-14    2
-7.90294013E+03-3.98525731E+00 1.81422511E+00 1.99008590E-02-2.21416008E-05    3
 1.45028521E-08-3.98877068E-12-7.05394926E+03 1.36079359E+01                   4
CH2CHO            000000H   3O   1C   2     G    300.00   5000.00 1000.00      1
 6.53928338E+00 7.80238629E-03-2.76413612E-06 4.42098906E-10-2.62954290E-14    2
-1.18858659E+03-8.72091393E+00 2.79502600E+00 1.01099472E-02 1.61750645E-05    3
-3.10303145E-08 1.39436139E-11 1.62944975E+02 1.23646657E+01                   4
CH3CO3H           000000H   4O   3C   2     G    300.00   5000.00 1000.00      1
 1.25060485E+01 9.47789695E-03-3.30402246E-06 5.19630793E-10-3.04233568E-14    2
-4.59856703E+04-3.79195947E+01 2.24135876E+00 3.37963514E-02-2.53887482E-05    3
 9.67583587E-09-1.49266157E-12-4.24677831E+04 1.70668133E+01                   4
CH3CO3            000000H   3O   3C   2     G    300.00   5000.00 1000.00      1
 1.12522498E+01 8.33652672E-03-2.89014530E-06 4.52781734E-10-2.64354456E-14    2
-2.60238584E+04-2.96370457E+01 3.60373432E+00 2.70080341E-02-2.08293438E-05    3
 8.50541104E-09-1.43846110E-12-2.34205171E+04 1.12014914E+01                   4
H2CC              000000H   2C   2          G    300.00   5000.00 1000.00      1
 4.27803400E+00 4.75628040E-03-1.63010090E-06 2.54628060E-10-1.48863790E-14    2
 4.83166880E+04 6.40237010E-01 3.28154830E+00 6.97647910E-03-2.38552440E-06    3
-1.21044320E-09 9.81895450E-13 4.86217940E+04 5.92039100E+00                   4
C2H2              000000H   2C   2          G    300.00   5000.00 1000.00      1
 4.65878489E+00 4.88396667E-03-1.60828888E-06 2.46974544E-10-1.38605959E-14    2
 2.57594042E+04-3.99838194E+00 8.08679682E-01 2.33615762E-02-3.55172234E-05    3
 2.80152958E-08-8.50075165E-12 2.64289808E+04 1.39396761E+01                   4
OCHO              000000H   1O   2C   1     G    300.00   5000.00 1000.00      1
 4.14394211E+00 5.59738818E-03-1.99794019E-06 3.16179193E-10-1.85614483E-14    2
-1.72459887E+04 5.07778617E+00 4.68825921E+00-4.14871834E-03 2.55066010E-05    3
-2.84473900E-08 1.04422559E-11-1.69867041E+04 4.28426480E+00                   4
SC2H4OH           000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 6.35842302E+00 1.24356276E-02-4.33096839E-06 6.84530381E-10-4.03713238E-14    2
-9.37900432E+03-6.05106112E+00 4.22283250E+00 5.12174798E-03 3.48386522E-05    3
-4.91943637E-08 2.01183723E-11-8.20503939E+03 8.01675700E+00                   4
C2H5OH            000000H   6O   1C   2     G    300.00   5000.00 1000.00      1
 6.56243650E+00 1.52042220E-02-5.38967950E-06 8.62250110E-10-5.12897870E-14    2
-3.15256210E+04-9.47302020E+00 4.85869570E+00-3.74017260E-03 6.95553780E-05    3
-8.86547960E-08 3.51688350E-11-2.99961320E+04 4.80185450E+00                   4
CH3CO2            000000H   3O   2C   2     G    300.00   5000.00 1000.00      1
 8.54059736E+00 8.32951214E-03-2.84722010E-06 4.41927196E-10-2.56373394E-14    2
-2.97290678E+04-2.03883545E+01 1.37440768E+00 2.49115604E-02-1.74308894E-05    3
 6.24799508E-09-9.09516835E-13-2.72330150E+04 1.81405454E+01                   4
CH3CHO            000000H   4O   1C   2     G    300.00   5000.00 1000.00      1
 5.40411080E+00 1.17230590E-02-4.22631370E-06 6.83724510E-10-4.09848630E-14    2
-2.25931220E+04-3.48079170E+00 4.72945950E+00-3.19328580E-03 4.75349210E-05    3
-5.74586110E-08 2.19311120E-11-2.15728780E+04 4.10301590E+00                   4
HCCO              000000H   1O   1C   2     G    300.00   5000.00 1000.00      1
 5.91479333E+00 3.71408730E-03-1.30137010E-06 2.06473345E-10-1.21476759E-14    2
 1.93596301E+04-5.50567269E+00 1.87607969E+00 2.21205418E-02-3.58869325E-05    3
 3.05402541E-08-1.01281069E-11 2.01633840E+04 1.36968290E+01                   4
HCOH              000000H   2O   1C   1     G    300.00   5000.00 1000.00      1
 9.18749272E+00 1.52011152E-03-6.27603516E-07 1.09727989E-10-6.89655128E-15    2
 7.81364593E+03-2.73434214E+01-2.82157421E+00 3.57331702E-02-3.80861580E-05    3
 1.86205951E-08-3.45957838E-12 1.12956672E+04 3.48487757E+01                   4
C2H3OH            000000H   4O   1C   2     G    300.00   5000.00 1000.00      1
 8.32598158E+00 8.03387281E-03-2.63928405E-06 3.98410726E-10-2.26551155E-14    2
-1.83221436E+04-2.02080305E+01-1.27972260E-01 3.38506073E-02-3.30644935E-05    3
 1.64858739E-08-3.19935455E-12-1.59914544E+04 2.30438601E+01                   4
C2H5O             000000H   5O   1C   2     G    300.00   5000.00 1000.00      1
 6.68899820E+00 1.31256760E-02-4.70388400E-06 7.58585520E-10-4.54133060E-14    2
-4.74578320E+03-9.69837550E+00 4.30742680E+00 6.41472050E-03 3.11397140E-05    3
-4.33140830E-08 1.72761840E-11-3.40275240E+03 5.90258370E+00                   4
C2H5CO            000000H   5O   1C   3     G    300.00   5000.00 1000.00      1
 6.52325448E+00 1.54211952E-02-5.50898157E-06 8.85889862E-10-5.28846399E-14    2
-7.19631634E+03-5.19862218E+00 6.25722402E+00-9.17612184E-03 7.61190493E-05    3
-9.05514997E-08 3.46198215E-11-5.91616484E+03 2.23330599E+00                   4
C2H5CHO           000000H   6O   1C   3     G    300.00   5000.00 1000.00      1
 7.44085690E+00 1.77301764E-02-6.34081568E-06 1.02040803E-09-6.09461714E-14    2
-2.60055814E+04-1.44195446E+01 4.24529681E+00 6.68296706E-03 4.93337933E-05    3
-6.71986124E-08 2.67262347E-11-2.41473007E+04 6.90738560E+00                   4
C2H               000000H   1C   2          G    300.00   5000.00 1000.00      1
 3.66270248E+00 3.82492252E-03-1.36632500E-06 2.13455040E-10-1.23216848E-14    2
 6.71683790E+04 3.92205792E+00 2.89867676E+00 1.32988489E-02-2.80733327E-05    3
 2.89484755E-08-1.07502351E-11 6.70616050E+04 6.18547632E+00                   4
CH3COCH3          000000H   6O   1C   3     G    300.00   5000.00 1000.00      1
 7.29796974E+00 1.75656913E-02-6.31678065E-06 1.02025553E-09-6.10903592E-14    2
-2.95368927E+04-1.27591704E+01 5.55638920E+00-2.83863547E-03 7.05722951E-05    3
-8.78130984E-08 3.40290951E-11-2.78325393E+04 2.31960221E+00                   4
CH3COCH2          000000H   5O   1C   3     G    300.00   5000.00 1000.00      1
 7.54410697E+00 1.43443222E-02-5.08381081E-06 8.13200521E-10-4.83673315E-14    2
-7.48672286E+03-1.14792587E+01 4.70187196E+00 5.51653762E-03 4.27505858E-05    3
-5.94680816E-08 2.40685378E-11-5.92845491E+03 7.12932590E+00                   4
CH2Y              000000H   2C   1          G    300.00   5000.00 1000.00      1
 3.13501686E+00 2.89593926E-03-8.16668090E-07 1.13572697E-10-6.36262835E-15    2
 5.05040504E+04 4.06030621E+00 4.19331325E+00-2.33105184E-03 8.15676451E-06    3
-6.62985981E-09 1.93233199E-12 5.03662246E+04-7.46734310E-01                   4
OCH2O2H           000000H   3O   3C   1     G    300.00   5000.00 1000.00      1
 1.15398246E+01 5.34291432E-03-1.81878917E-06 2.81968625E-10-1.63584348E-14    2
-1.68237489E+04-3.20700633E+01 1.93823075E+00 3.01465730E-02-2.61053152E-05    3
 1.09463562E-08-1.78312692E-12-1.38166625E+04 1.85042002E+01                   4
HE                000000HE  1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 9.28723974E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.28723974E-01                   4
AR                000000AR  1               G    300.00   5000.00 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967491E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967491E+00                   4
)THERMORES";
}

} // namespace ideal_gas
} // namespace fub