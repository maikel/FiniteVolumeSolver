#include "fub/solver/euler/mechanism/Burke2012.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace fub::euler {
namespace {
double GetLindRateCoeff(double temp, double pressure, double k0, double kInf,
                        double fc, double conc);

double MAX_C(double X1, double X2);
} // namespace

/// This function computes rates of production cdot in [kmole/(m^3s)].
///
/// The parameters w ( reaction rate ), k ( rate coefficient )
/// and M ( third body concentrations ) are just work space for this
/// function. c contains the concentrations of non steady state species in
/// [kmole/m^3] and is workspace for the steady state concentrations, which are
/// computed in this function. temp is the temperature in [K] and pressure is
/// the pressure in [Pa]. Called functions are 'GetLindRateCoeff',
/// 'ComputeSteadyStates', 'CatchZero' and the functions that evaluate the
/// broadening factors of the Troe formulation of pressure dependant rate
/// coefficients 'Fc*'
void Burke2012::ComputeProductionRates(span<double> cdot, span<double> w,
                                       span<double> k, span<double> c,
                                       span<double> M, double temp,
                                       double pressure) const {
  int nSpec = 11;
  int nSpecIn = 11;
  double kTroe0, kTroeInf, fcTroe;
  double RGAS = 8314.34;
  double lgt = log(temp);
  double rt = RGAS * temp;

  M[mM1] = c[sN2] + c[sH] + c[sO2] + c[sO] + c[sOH] + 2.5 * c[sH2] +
           12 * c[sH2O] + c[sHO2] + c[sH2O2];

  M[mM2] = c[sN2] + c[sH] + c[sO2] + c[sO] + c[sOH] + 2.5 * c[sH2] +
           12 * c[sH2O] + c[sHO2] + c[sH2O2];

  M[mM3] = c[sN2] + 0.75 * c[sAR] + c[sH] + c[sO2] + c[sO] + c[sOH] +
           2.5 * c[sH2] + 12 * c[sH2O] + 0.75 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM4] = 2 * c[sN2] + c[sAR] + c[sH] + 1.5 * c[sO2] + c[sO] + c[sOH] +
           3 * c[sH2] + 1.1 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM5] = c[sN2] + 0.67 * c[sAR] + c[sH] + 0.78 * c[sO2] + c[sO] + c[sOH] +
           2 * c[sH2] + 14 * c[sH2O] + 0.8 * c[sHE] + c[sHO2] + c[sH2O2];

  M[mM6] = 1.5 * c[sN2] + c[sAR] + c[sH] + 1.2 * c[sO2] + c[sO] + c[sOH] +
           3.7 * c[sH2] + 7.5 * c[sH2O] + 0.65 * c[sHE] + c[sHO2] +
           7.7 * c[sH2O2];

  k[r1f] = 1.0400000000E+11 * exp(-63957000 / rt);
  k[r1b] = 2.0679845118E+08 * exp(0.434958 * lgt + 6605168.02 / rt);
  k[r2f] = 3.8180000000E+09 * exp(-33254000 / rt);
  k[r2b] = 2.2393817467E+09 * exp(-0.0350345 * lgt - 27487825.08 / rt);
  k[r3f] = 8.7920000000E+11 * exp(-80207000 / rt);
  k[r3b] = 5.1567952637E+11 * exp(-0.0350345 * lgt - 74440825.08 / rt);
  k[r4f] = 2.1600000000E+05 * exp(1.51 * lgt - 14351000 / rt);
  k[r4b] = 2.2182689717E+06 * exp(1.41638 * lgt - 77060225.89 / rt);
  k[r5f] = 3.3400000000E+01 * exp(2.42 * lgt + 8075000 / rt);
  k[r5b] = 5.8480989232E+02 * exp(2.36141 * lgt - 60400400.81 / rt);
  k[r6f] = 4.5770000000E+16 * exp(-1.4 * lgt - 436726000 / rt);
  k[r6b] = 1.9619373578E+14 * exp(-1.7488 * lgt - 3646272.681 / rt);
  k[r7f] = 5.8400000000E+15 * exp(-1.1 * lgt - 436726000 / rt);
  k[r7b] = 2.5033240484E+13 * exp(-1.4488 * lgt - 3646272.681 / rt);
  k[r8f] = 5.8400000000E+15 * exp(-1.1 * lgt - 436726000 / rt);
  k[r8b] = 2.5033240484E+13 * exp(-1.4488 * lgt - 3646272.681 / rt);
  k[r9f] = 6.1650000000E+09 * exp(-0.5 * lgt);
  k[r9b] = 4.2423561386E+14 * exp(-0.621192 * lgt - 497875720.4 / rt);
  k[r10f] = 1.8860000000E+07 * exp(7481000 / rt);
  k[r10b] = 1.2978237919E+12 * exp(-0.121192 * lgt - 490394720.4 / rt);
  k[r11f] = 1.8860000000E+07 * exp(7481000 / rt);
  k[r11b] = 1.2978237919E+12 * exp(-0.121192 * lgt - 490394720.4 / rt);
  k[r12f] = 4.7140000000E+12 * exp(-1 * lgt);
  k[r12b] = 6.4502650942E+14 * exp(-0.686234 * lgt - 427313552.4 / rt);
  k[r13f] = 6.0640000000E+24 * exp(-3.322 * lgt - 505385000 / rt);
  k[r13b] = 2.5310630492E+21 * exp(-3.57718 * lgt - 9596046.787 / rt);
  k[r14f] = 1.0060000000E+23 * exp(-2.44 * lgt - 502833000 / rt);
  k[r14b] = 4.1989601376E+19 * exp(-2.69518 * lgt - 7044046.787 / rt);
  kTroe0 = 6.3660000000E+14 * exp(-1.72 * lgt - 2196000 / rt);
  kTroeInf = 4.6510000000E+09 * exp(0.44 * lgt);
  fcTroe = 0.5 * exp(-temp / 1e-30) + 0.5 * exp(-temp / 1e+30);
  k[r15f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  kTroe0 = 9.0292117899E+17 * exp(-1.73291 * lgt - 206862356.2 / rt);
  kTroeInf = 6.5967427010E+12 * exp(0.427093 * lgt - 204666356.2 / rt);
  k[r15b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM5]);
  k[r16f] = 2.7500000000E+03 * exp(2.09 * lgt + 6071000 / rt);
  k[r16b] = 4.5231945076E+02 * exp(2.45171 * lgt - 222342371.2 / rt);
  k[r17f] = 7.0790000000E+10 * exp(-1234000 / rt);
  k[r17b] = 1.3579714350E+07 * exp(0.761631 * lgt - 153319028.2 / rt);
  k[r18f] = 2.8500000000E+07 * exp(1 * lgt + 3029000 / rt);
  k[r18b] = 2.7494741434E+06 * exp(1.32667 * lgt - 219618196.2 / rt);
  k[r19f] = 2.8900000000E+10 * exp(2079000 / rt);
  k[r19b] = 4.8816975193E+10 * exp(0.268083 * lgt - 289043597.1 / rt);
  k[r20f] = 4.2000000000E+11 * exp(-50133000 / rt);
  k[r20b] = 6.3574616360E+13 * exp(-0.389273 * lgt - 212211053.7 / rt);
  k[r21f] = 1.3000000000E+08 * exp(6817000 / rt);
  k[r21b] = 1.9677857445E+10 * exp(-0.389273 * lgt - 155261053.7 / rt);
  kTroe0 = 2.4900000000E+21 * exp(-2.3 * lgt - 203970000 / rt);
  kTroeInf = 2.0000000000E+12 * exp(0.9 * lgt - 203966000 / rt);
  fcTroe = 0.57 * exp(-temp / 1e-30) + 0.43 * exp(-temp / 1e+30);
  k[r22f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  kTroe0 = 2.2248483092E+12 * exp(-1.13619 * lgt + 10689381.66 / rt);
  kTroeInf = 1.7870267543E+03 * exp(2.06381 * lgt + 10693381.66 / rt);
  k[r22b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM6]);
  k[r23f] = 2.4100000000E+10 * exp(-16610000 / rt);
  k[r23b] = 5.1591045676E+04 * exp(1.41899 * lgt - 297739571.6 / rt);
  k[r24f] = 4.8200000000E+10 * exp(-33263000 / rt);
  k[r24b] = 5.2375134409E+07 * exp(0.75098 * lgt - 99598317.44 / rt);
  k[r25f] = 9.5500000000E+03 * exp(2 * lgt - 16610000 / rt);
  k[r25b] = 6.0865850326E+00 * exp(2.71595 * lgt - 77179142.52 / rt);
  k[r26f] = 1.7400000000E+09 * exp(-1331000 / rt);
  k[r26b] = 1.9417254097E+07 * exp(0.657356 * lgt - 130375543.3 / rt);
  k[r27f] = 7.5900000000E+10 * exp(-30418000 / rt);
  k[r27b] = 8.4699401492E+08 * exp(0.657356 * lgt - 159462543.3 / rt);

  w[r1f] = k[r1f] * c[sH] * c[sO2];
  w[r1b] = k[r1b] * c[sOH] * c[sO];
  w[r2f] = k[r2f] * c[sO] * c[sH2];
  w[r2b] = k[r2b] * c[sOH] * c[sH];
  w[r3f] = k[r3f] * c[sO] * c[sH2];
  w[r3b] = k[r3b] * c[sOH] * c[sH];
  w[r4f] = k[r4f] * c[sH2] * c[sOH];
  w[r4b] = k[r4b] * c[sH] * c[sH2O];
  w[r5f] = k[r5f] * c[sOH] * c[sOH];
  w[r5b] = k[r5b] * c[sH2O] * c[sO];
  w[r6f] = k[r6f] * c[sH2] * M[mM1];
  w[r6b] = k[r6b] * c[sH] * c[sH] * M[mM1];
  w[r7f] = k[r7f] * c[sH2] * c[sAR];
  w[r7b] = k[r7b] * c[sAR] * c[sH] * c[sH];
  w[r8f] = k[r8f] * c[sH2] * c[sHE];
  w[r8b] = k[r8b] * c[sHE] * c[sH] * c[sH];
  w[r9f] = k[r9f] * c[sO] * c[sO] * M[mM2];
  w[r9b] = k[r9b] * c[sO2] * M[mM2];
  w[r10f] = k[r10f] * c[sO] * c[sO] * c[sAR];
  w[r10b] = k[r10b] * c[sAR] * c[sO2];
  w[r11f] = k[r11f] * c[sO] * c[sO] * c[sHE];
  w[r11b] = k[r11b] * c[sHE] * c[sO2];
  w[r12f] = k[r12f] * c[sO] * c[sH] * M[mM3];
  w[r12b] = k[r12b] * c[sOH] * M[mM3];
  w[r13f] = k[r13f] * c[sH2O] * M[mM4];
  w[r13b] = k[r13b] * c[sOH] * c[sH] * M[mM4];
  w[r14f] = k[r14f] * c[sH2O] * c[sH2O];
  w[r14b] = k[r14b] * c[sH2O] * c[sOH] * c[sH];
  w[r15f] = k[r15f] * c[sH] * c[sO2];
  w[r15b] = k[r15b] * c[sHO2];
  w[r16f] = k[r16f] * c[sHO2] * c[sH];
  w[r16b] = k[r16b] * c[sO2] * c[sH2];
  w[r17f] = k[r17f] * c[sHO2] * c[sH];
  w[r17b] = k[r17b] * c[sOH] * c[sOH];
  w[r18f] = k[r18f] * c[sHO2] * c[sO];
  w[r18b] = k[r18b] * c[sOH] * c[sO2];
  w[r19f] = k[r19f] * c[sHO2] * c[sOH];
  w[r19b] = k[r19b] * c[sO2] * c[sH2O];
  w[r20f] = k[r20f] * c[sHO2] * c[sHO2];
  w[r20b] = k[r20b] * c[sO2] * c[sH2O2];
  w[r21f] = k[r21f] * c[sHO2] * c[sHO2];
  w[r21b] = k[r21b] * c[sO2] * c[sH2O2];
  w[r22f] = k[r22f] * c[sH2O2];
  w[r22b] = k[r22b] * c[sOH] * c[sOH];
  w[r23f] = k[r23f] * c[sH2O2] * c[sH];
  w[r23b] = k[r23b] * c[sOH] * c[sH2O];
  w[r24f] = k[r24f] * c[sH2O2] * c[sH];
  w[r24b] = k[r24b] * c[sH2] * c[sHO2];
  w[r25f] = k[r25f] * c[sH2O2] * c[sO];
  w[r25b] = k[r25b] * c[sHO2] * c[sOH];
  w[r26f] = k[r26f] * c[sH2O2] * c[sOH];
  w[r26b] = k[r26b] * c[sH2O] * c[sHO2];
  w[r27f] = k[r27f] * c[sH2O2] * c[sOH];
  w[r27b] = k[r27b] * c[sH2O] * c[sHO2];

  cdot[sN2] = 0.0;

  cdot[sAR] = -w[r7f] + w[r7f] - w[r7b] + w[r7b] - w[r10f] + w[r10f] - w[r10b] +
              w[r10b];

  cdot[sH] = -w[r1f] + w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] + w[r4f] -
             w[r4b] + 2 * w[r6f] - 2 * w[r6b] + 2 * w[r7f] - 2 * w[r7b] +
             2 * w[r8f] - 2 * w[r8b] - w[r12f] + w[r12b] + w[r13f] - w[r13b] +
             w[r14f] - w[r14b] - w[r15f] + w[r15b] - w[r16f] + w[r16b] -
             w[r17f] + w[r17b] - w[r23f] + w[r23b] - w[r24f] + w[r24b];

  cdot[sO2] = -w[r1f] + w[r1b] + w[r9f] - w[r9b] + w[r10f] - w[r10b] + w[r11f] -
              w[r11b] - w[r15f] + w[r15b] + w[r16f] - w[r16b] + w[r18f] -
              w[r18b] + w[r19f] - w[r19b] + w[r20f] - w[r20b] + w[r21f] -
              w[r21b];

  cdot[sO] = w[r1f] - w[r1b] - w[r2f] + w[r2b] - w[r3f] + w[r3b] + w[r5f] -
             w[r5b] - 2 * w[r9f] + 2 * w[r9b] - 2 * w[r10f] + 2 * w[r10b] -
             2 * w[r11f] + 2 * w[r11b] - w[r12f] + w[r12b] - w[r18f] + w[r18b] -
             w[r25f] + w[r25b];

  cdot[sOH] = w[r1f] - w[r1b] + w[r2f] - w[r2b] + w[r3f] - w[r3b] - w[r4f] +
              w[r4b] - 2 * w[r5f] + 2 * w[r5b] + w[r12f] - w[r12b] + w[r13f] -
              w[r13b] + w[r14f] - w[r14b] + 2 * w[r17f] - 2 * w[r17b] +
              w[r18f] - w[r18b] - w[r19f] + w[r19b] + 2 * w[r22f] -
              2 * w[r22b] + w[r23f] - w[r23b] + w[r25f] - w[r25b] - w[r26f] +
              w[r26b] - w[r27f] + w[r27b];

  cdot[sH2] = -w[r2f] + w[r2b] - w[r3f] + w[r3b] - w[r4f] + w[r4b] - w[r6f] +
              w[r6b] - w[r7f] + w[r7b] - w[r8f] + w[r8b] + w[r16f] - w[r16b] +
              w[r24f] - w[r24b];

  cdot[sH2O] = w[r4f] - w[r4b] + w[r5f] - w[r5b] - w[r13f] + w[r13b] -
               2 * w[r14f] + w[r14f] - w[r14b] + 2 * w[r14b] + w[r19f] -
               w[r19b] + w[r23f] - w[r23b] + w[r26f] - w[r26b] + w[r27f] -
               w[r27b];

  cdot[sHE] = -w[r8f] + w[r8f] - w[r8b] + w[r8b] - w[r11f] + w[r11f] - w[r11b] +
              w[r11b];

  cdot[sHO2] = w[r15f] - w[r15b] - w[r16f] + w[r16b] - w[r17f] + w[r17b] -
               w[r18f] + w[r18b] - w[r19f] + w[r19b] - 2 * w[r20f] +
               2 * w[r20b] - 2 * w[r21f] + 2 * w[r21b] + w[r24f] - w[r24b] +
               w[r25f] - w[r25b] + w[r26f] - w[r26b] + w[r27f] - w[r27b];

  cdot[sH2O2] = w[r20f] - w[r20b] + w[r21f] - w[r21b] - w[r22f] + w[r22b] -
                w[r23f] + w[r23b] - w[r24f] + w[r24b] - w[r25f] + w[r25b] -
                w[r26f] + w[r26b] - w[r27f] + w[r27b];
}

/// This function computes enthalpy 'h' and heat capacity 'cp' as
/// function of temperature 'T' for all non steady state species
/// in units [J/kg] and [J/kg K], respectively.
/// The parameter h and cp should provide workspace of length 11
void Burke2012::ComputeThermoData(span<double> h, span<double> cp,
                                  double T) const {
  int i;
  if (T > 1000.0) {
    h[sN2] = 2.96733125e+02 *
             (T * (2.92664000e+00 +
                   T * (7.43988500e-04 +
                        T * (-1.89492033e-07 +
                             T * (2.52426000e-11 + T * -1.35067020e-15)))) -
              9.22797700e+02);
    cp[sN2] = 2.96733125e+02 *
              (2.92664000e+00 +
               T * (1.48797700e-03 +
                    T * (-5.68476100e-07 +
                         T * (1.00970400e-10 + T * -6.75335100e-15))));
    h[sAR] = 2.08132126e+02 *
             (T * (2.50000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 +
                             T * (0.00000000e+00 + T * 0.00000000e+00)))) -
              7.45375000e+02);
    cp[sAR] = 2.08132126e+02 *
              (2.50000000e+00 +
               T * (0.00000000e+00 +
                    T * (0.00000000e+00 +
                         T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sH] = 8.24847438e+03 *
            (T * (2.50000000e+00 +
                  T * (0.00000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 + T * 0.00000000e+00)))) +
             2.54716300e+04);
    cp[sH] = 8.24847438e+03 *
             (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sO2] = 2.59826943e+02 *
             (T * (3.69757800e+00 +
                   T * (3.06759850e-04 +
                        T * (-4.19614000e-08 +
                             T * (4.43820250e-12 + T * -2.27287000e-16)))) -
              1.23393000e+03);
    cp[sO2] = 2.59826943e+02 *
              (3.69757800e+00 +
               T * (6.13519700e-04 +
                    T * (-1.25884200e-07 +
                         T * (1.77528100e-11 + T * -1.13643500e-15))));
    h[sO] = 5.19653886e+02 *
            (T * (2.54206000e+00 +
                  T * (-1.37753100e-05 +
                       T * (-1.03426767e-09 +
                            T * (1.13776675e-12 + T * -8.73610400e-17)))) +
             2.92308000e+04);
    cp[sO] = 5.19653886e+02 *
             (2.54206000e+00 +
              T * (-2.75506200e-05 +
                   T * (-3.10280300e-09 +
                        T * (4.55106700e-12 + T * -4.36805200e-16))));
    h[sOH] = 4.88855960e+02 *
             (T * (2.86472886e+00 +
                   T * (5.28252240e-04 +
                        T * (-8.63609193e-08 +
                             T * (7.63046685e-12 + T * -2.66391752e-16)))) +
              3.68362875e+03);
    cp[sOH] = 4.88855960e+02 *
              (2.86472886e+00 +
               T * (1.05650448e-03 +
                    T * (-2.59082758e-07 +
                         T * (3.05218674e-11 + T * -1.33195876e-15))));
    h[sH2] = 4.12423719e+03 *
             (T * (2.99142300e+00 +
                   T * (3.50032200e-04 +
                        T * (-1.87794300e-08 +
                             T * (-2.30789450e-12 + T * 3.16550400e-16)))) -
              8.35034000e+02);
    cp[sH2] = 4.12423719e+03 *
              (2.99142300e+00 +
               T * (7.00064400e-04 +
                    T * (-5.63382900e-08 +
                         T * (-9.23157800e-12 + T * 1.58275200e-15))));
    h[sH2O] = 4.61504339e+02 *
              (T * (2.67214600e+00 +
                    T * (1.52814650e-03 +
                         T * (-2.91008667e-07 +
                              T * (3.00249000e-11 + T * -1.27832360e-15)))) -
               2.98992100e+04);
    cp[sH2O] = 4.61504339e+02 *
               (2.67214600e+00 +
                T * (3.05629300e-03 +
                     T * (-8.73026000e-07 +
                          T * (1.20099600e-10 + T * -6.39161800e-15))));
    h[sHE] = 2.07861554e+03 *
             (T * (2.50000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 +
                             T * (0.00000000e+00 + T * 0.00000000e+00)))) -
              7.45375000e+02);
    cp[sHE] = 2.07861554e+03 *
              (2.50000000e+00 +
               T * (0.00000000e+00 +
                    T * (0.00000000e+00 +
                         T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sHO2] = 2.51892334e+02 *
              (T * (4.01721090e+00 +
                    T * (1.11991006e-03 +
                         T * (-2.11219383e-07 +
                              T * (2.85615925e-11 + T * -2.15817070e-15)))) +
               1.11856713e+02);
    cp[sHO2] = 2.51892334e+02 *
               (4.01721090e+00 +
                T * (2.23982013e-03 +
                     T * (-6.33658150e-07 +
                          T * (1.14246370e-10 + T * -1.07908535e-14))));
    h[sH2O2] = 2.44427980e+02 *
               (T * (4.57316700e+00 +
                     T * (2.16806800e-03 +
                          T * (-4.91563000e-07 +
                               T * (5.87226000e-11 + T * -2.86330800e-15)))) -
                1.80069600e+04);
    cp[sH2O2] = 2.44427980e+02 *
                (4.57316700e+00 +
                 T * (4.33613600e-03 +
                      T * (-1.47468900e-06 +
                           T * (2.34890400e-10 + T * -1.43165400e-14))));
  } else if (T >= 299.999999) {
    h[sN2] = 2.96733125e+02 *
             (T * (3.29867700e+00 +
                   T * (7.04120000e-04 +
                        T * (-1.32107400e-06 +
                             T * (1.41037875e-09 + T * -4.88971000e-13)))) -
              1.02090000e+03);
    cp[sN2] = 2.96733125e+02 *
              (3.29867700e+00 +
               T * (1.40824000e-03 +
                    T * (-3.96322200e-06 +
                         T * (5.64151500e-09 + T * -2.44485500e-12))));
    h[sAR] = 2.08132126e+02 *
             (T * (2.50000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 +
                             T * (0.00000000e+00 + T * 0.00000000e+00)))) -
              7.45375000e+02);
    cp[sAR] = 2.08132126e+02 *
              (2.50000000e+00 +
               T * (0.00000000e+00 +
                    T * (0.00000000e+00 +
                         T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sH] = 8.24847438e+03 *
            (T * (2.50000000e+00 +
                  T * (0.00000000e+00 +
                       T * (0.00000000e+00 +
                            T * (0.00000000e+00 + T * 0.00000000e+00)))) +
             2.54716300e+04);
    cp[sH] = 8.24847438e+03 *
             (2.50000000e+00 +
              T * (0.00000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sO2] = 2.59826943e+02 *
             (T * (3.21293600e+00 +
                   T * (5.63743000e-04 +
                        T * (-1.91871667e-07 +
                             T * (3.28469250e-10 + T * -1.75371080e-13)))) -
              1.00524900e+03);
    cp[sO2] = 2.59826943e+02 *
              (3.21293600e+00 +
               T * (1.12748600e-03 +
                    T * (-5.75615000e-07 +
                         T * (1.31387700e-09 + T * -8.76855400e-13))));
    h[sO] = 5.19653886e+02 *
            (T * (2.94642900e+00 +
                  T * (-8.19083000e-04 +
                       T * (8.07010667e-07 +
                            T * (-4.00710750e-10 + T * 7.78139200e-14)))) +
             2.91476400e+04);
    cp[sO] = 5.19653886e+02 *
             (2.94642900e+00 +
              T * (-1.63816600e-03 +
                   T * (2.42103200e-06 +
                        T * (-1.60284300e-09 + T * 3.89069600e-13))));
    h[sOH] = 4.88855960e+02 *
             (T * (4.12530561e+00 +
                   T * (-1.61272470e-03 +
                        T * (2.17588230e-06 +
                             T * (-1.44963411e-09 + T * 4.12474758e-13)))) +
              3.34630913e+03);
    cp[sOH] = 4.88855960e+02 *
              (4.12530561e+00 +
               T * (-3.22544939e-03 +
                    T * (6.52764691e-06 +
                         T * (-5.79853643e-09 + T * 2.06237379e-12))));
    h[sH2] = 4.12423719e+03 *
             (T * (3.29812400e+00 +
                   T * (4.12472100e-04 +
                        T * (-2.71433833e-07 +
                             T * (-2.36885850e-11 + T * 8.26974400e-14)))) -
              1.01252100e+03);
    cp[sH2] = 4.12423719e+03 *
              (3.29812400e+00 +
               T * (8.24944200e-04 +
                    T * (-8.14301500e-07 +
                         T * (-9.47543400e-11 + T * 4.13487200e-13))));
    h[sH2O] = 4.61504339e+02 *
              (T * (3.38684200e+00 +
                    T * (1.73749100e-03 +
                         T * (-2.11823200e-06 +
                              T * (1.74214525e-09 + T * -5.01317600e-13)))) -
               3.02081100e+04);
    cp[sH2O] = 4.61504339e+02 *
               (3.38684200e+00 +
                T * (3.47498200e-03 +
                     T * (-6.35469600e-06 +
                          T * (6.96858100e-09 + T * -2.50658800e-12))));
    h[sHE] = 2.07861554e+03 *
             (T * (2.50000000e+00 +
                   T * (0.00000000e+00 +
                        T * (0.00000000e+00 +
                             T * (0.00000000e+00 + T * 0.00000000e+00)))) -
              7.45375000e+02);
    cp[sHE] = 2.07861554e+03 *
              (2.50000000e+00 +
               T * (0.00000000e+00 +
                    T * (0.00000000e+00 +
                         T * (0.00000000e+00 + T * 0.00000000e+00))));
    h[sHO2] = 2.51892334e+02 *
              (T * (4.30179801e+00 +
                    T * (-2.37456025e-03 +
                         T * (7.05276303e-06 +
                              T * (-6.06909735e-09 + T * 1.85845025e-12)))) +
               2.94808040e+02);
    cp[sHO2] = 2.51892334e+02 *
               (4.30179801e+00 +
                T * (-4.74912051e-03 +
                     T * (2.11582891e-05 +
                          T * (-2.42763894e-08 + T * 9.29225124e-12))));
    h[sH2O2] = 2.44427980e+02 *
               (T * (3.38875400e+00 +
                     T * (3.28461300e-03 +
                          T * (-4.95004333e-08 +
                               T * (-1.15645150e-09 + T * 4.94303000e-13)))) -
                1.76631500e+04);
    cp[sH2O2] = 2.44427980e+02 *
                (3.38875400e+00 +
                 T * (6.56922600e-03 +
                      T * (-1.48501300e-07 +
                           T * (-4.62580600e-09 + T * 2.47151500e-12))));
  } else {
    ComputeThermoData(h, cp, 300.0);
    for (i = 0; i < sEnd; i++) {
      h[i] = (T - 300.) * cp[i] + h[i];
    }
  }
}

namespace {
double GetLindRateCoeff(double temp, double pressure, double k0, double kInf,
                        double fc, double conc) {
  const double R = 8314.34; /* [J / kmole K] */
  double Ntmp;
  double kl;
  double f;
  double cCoeff, dCoeff, log10kNull;

  int iTroe = 1;

  if (conc <= 0.0) {
    conc = pressure / (R * temp);
  }
  Ntmp = 0.75 - 1.27 * log10(fc);
  if (iTroe) {
    cCoeff = -0.4 - 0.67 * log10(fc);
    dCoeff = 0.14;
    k0 *= conc / MAX_C(kInf, 1.0e-60);
    log10kNull = log10(k0);
    f = (log10kNull + cCoeff) / (Ntmp - dCoeff * (log10kNull + cCoeff));
    f = pow(fc, 1.0 / (f * f + 1.0));
    kInf *= f * k0 / (1.0 + k0);
  } else {
    k0 = k0 * conc / kInf;
    kl = k0 / (1.0 + k0);
    f = log10(k0) / Ntmp;
    f = pow(fc, 1.0 / (f * f + 1.0));
    kInf = kInf * f * kl;
  }
  return kInf;
}

double MAX_C(double X1, double X2) { return ((X1 > X2) ? X1 : X2); }
} // namespace
} // namespace fub::euler