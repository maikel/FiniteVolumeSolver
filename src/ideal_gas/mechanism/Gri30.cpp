#include "fub/ideal_gas/mechanism/Gri30.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace fub {
namespace ideal_gas {

static double GetLindRateCoeff(double temp, double pressure, double k0,
                               double kInf, double fc, double conc);

double MAX_C(double X1, double X2);

void Gri30::ComputeProductionRates(span<double> cdot, span<double> w,
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

  int nSpec = 53;
  int nSpecIn = 53;
  double kTroe0, kTroeInf, fcTroe;
  double RGAS = 8314.34;
  double lgt = log(temp);
  double rt = RGAS * temp;

  M[mM1] = c[sAR] + c[sO] + c[sO2] + c[sH] + c[sOH] +
           0x1.3333333333333p+1 * c[sH2] + c[sHO2] + c[sH2O2] + c[sCH] +
           0x1.cp+0 * c[sCO] + c[s3XCH2] + c[sHCO] + c[s1XCH2] + c[sCH3] +
           c[sCH2O] + 0x1p+1 * c[sCH4] + 0x1.ccccccccccccdp+1 * c[sCO2] +
           c[sCH2OH] + c[sCH3O] + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO] +
           c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5] + 0x1.8p+1 * c[sC2H6] +
           0x1.ecccccccccccdp+3 * c[sH2O] + c[sN2] + c[sC] + c[sHCCOH] + c[sN] +
           c[sNO] + c[sN2O] + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH] +
           c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN] + c[sHNCO] + c[sH2CN] +
           c[sHCNN] + c[sHCNO] + c[sNH3] + c[sCH2CHO] + c[sCH3CHO] + c[sC3H8] +
           c[sC3H7];

  M[mM2] = c[sAR] + c[sO] + c[sO2] + c[sH] + c[sOH] + 0x1p+1 * c[sH2] +
           c[sHO2] + c[sH2O2] + c[sCH] + 0x1.8p+0 * c[sCO] + c[s3XCH2] +
           c[sHCO] + c[s1XCH2] + c[sCH3] + c[sCH2O] + 0x1p+1 * c[sCH4] +
           0x1p+1 * c[sCO2] + c[sCH2OH] + c[sCH3O] + c[sCH3OH] + c[sC2H] +
           c[sC2H2] + c[sHCCO] + c[sC2H3] + c[sCH2CO] + c[sC2H4] + c[sC2H5] +
           0x1.8p+1 * c[sC2H6] + 0x1.8p+2 * c[sH2O] + c[sN2] + c[sC] +
           c[sHCCOH] + c[sN] + c[sNO] + c[sN2O] + c[sNO2] + c[sNH] + c[sHNO] +
           c[sNH2] + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN] +
           c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO] + c[sNH3] + c[sCH2CHO] +
           c[sCH3CHO] + c[sC3H8] + c[sC3H7];

  M[mM3] = c[sAR] + c[sO] + 0x1.8p+2 * c[sO2] + c[sH] + c[sOH] +
           0x1p+1 * c[sH2] + c[sHO2] + c[sH2O2] + c[sCH] + 0x1.8p+0 * c[sCO] +
           c[s3XCH2] + c[sHCO] + c[s1XCH2] + c[sCH3] + c[sCH2O] +
           0x1p+1 * c[sCH4] + 0x1.cp+1 * c[sCO2] + c[sCH2OH] + c[sCH3O] +
           c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO] + c[sC2H3] + c[sCH2CO] +
           c[sC2H4] + c[sC2H5] + 0x1.8p+1 * c[sC2H6] + 0x1.8p+2 * c[sH2O] +
           c[sN2] + c[sC] + c[sHCCOH] + c[sN] + c[sNO] + c[sN2O] + c[sNO2] +
           c[sNH] + c[sHNO] + c[sNH2] + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN] +
           c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO] + c[sNH3] +
           c[sCH2CHO] + c[sCH3CHO] + c[sC3H8] + c[sC3H7];

  M[mM4] = c[sAR] + c[sO] + c[sH] + c[sOH] + c[sH2] + c[sHO2] + c[sH2O2] +
           c[sCH] + 0x1.8p-1 * c[sCO] + c[s3XCH2] + c[sHCO] + c[s1XCH2] +
           c[sCH3] + c[sCH2O] + c[sCH4] + 0x1.8p+0 * c[sCO2] + c[sCH2OH] +
           c[sCH3O] + c[sCH3OH] + c[sC2H] + c[sC2H2] + c[sHCCO] + c[sC2H3] +
           c[sCH2CO] + c[sC2H4] + c[sC2H5] + 0x1.8p+0 * c[sC2H6] + c[sC] +
           c[sHCCOH] + c[sN] + c[sNO] + c[sN2O] + c[sNO2] + c[sNH] + c[sHNO] +
           c[sNH2] + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN] +
           c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO] + c[sNH3] + c[sCH2CHO] +
           c[sCH3CHO] + c[sC3H8] + c[sC3H7];

  M[mM5] = c[sAR] + c[sO] + c[sO2] + c[sH] + c[sOH] + c[sHO2] + c[sH2O2] +
           c[sCH] + c[sCO] + c[s3XCH2] + c[sHCO] + c[s1XCH2] + c[sCH3] +
           c[sCH2O] + 0x1p+1 * c[sCH4] + c[sCH2OH] + c[sCH3O] + c[sCH3OH] +
           c[sC2H] + c[sC2H2] + c[sHCCO] + c[sC2H3] + c[sCH2CO] + c[sC2H4] +
           c[sC2H5] + 0x1.8p+1 * c[sC2H6] + c[sN2] + c[sC] + c[sHCCOH] + c[sN] +
           c[sNO] + c[sN2O] + c[sNO2] + c[sNH] + c[sHNO] + c[sNH2] + c[sNNH] +
           c[sCN] + c[sNCO] + c[sHCN] + c[sHOCN] + c[sHNCO] + c[sH2CN] +
           c[sHCNN] + c[sHCNO] + c[sNH3] + c[sCH2CHO] + c[sCH3CHO] + c[sC3H8] +
           c[sC3H7];

  M[mM6] = c[sAR] + c[sO] + c[sO2] + c[sH] + c[sOH] +
           0x1.75c28f5c28f5cp-1 * c[sH2] + c[sHO2] + c[sH2O2] + c[sCH] +
           c[sCO] + c[s3XCH2] + c[sHCO] + c[s1XCH2] + c[sCH3] + c[sCH2O] +
           0x1p+1 * c[sCH4] + c[sCO2] + c[sCH2OH] + c[sCH3O] + c[sCH3OH] +
           c[sC2H] + c[sC2H2] + c[sHCCO] + c[sC2H3] + c[sCH2CO] + c[sC2H4] +
           c[sC2H5] + 0x1.8p+1 * c[sC2H6] + 0x1.d333333333333p+1 * c[sH2O] +
           c[sN2] + c[sC] + c[sHCCOH] + c[sN] + c[sNO] + c[sN2O] + c[sNO2] +
           c[sNH] + c[sHNO] + c[sNH2] + c[sNNH] + c[sCN] + c[sNCO] + c[sHCN] +
           c[sHOCN] + c[sHNCO] + c[sH2CN] + c[sHCNN] + c[sHCNO] + c[sNH3] +
           c[sCH2CHO] + c[sCH3CHO] + c[sC3H8] + c[sC3H7];

  k[r1f] = 0x1.bf08eb0000003p+36 * exp(-0x1p+0 * lgt);
  k[r1b] = 0x1.de35cb62ed243p+52 *
           exp(-0x1.1f22e28e9d4p+0 * lgt - 0x1.dac2b01f5aa41p+28 / rt);
  k[r2f] = 0x1.d1a94a2000003p+38 * exp(-0x1p+0 * lgt);
  k[r2b] = 0x1.f9b7650801d8p+45 *
           exp(-0x1.5ea1cc111fp-1 * lgt - 0x1.958186286074cp+28 / rt);
  k[r3f] = 0x1.359999999999ap+5 *
           exp(0x1.599999999999ap+1 * lgt - 0x1.8fb1ep+24 / rt);
  k[r3b] = 0x1.70435e6268c8p+4 *
           exp(0x1.551ca7aa65509p+1 * lgt - 0x1.1756a2baf28c9p+24 / rt);
  k[r4f] = 0x1.2a05f20000001p+34;
  k[r4b] = 0x1.cc26c642b4b0cp+30 *
           exp(0x1.4edbcce9a7ep-2 * lgt - 0x1.a4a39766f5ab5p+27 / rt);
  k[r5f] = 0x1.2cfp+13 * exp(0x1p+1 * lgt - 0x1.feca800000001p+23 / rt);
  k[r5b] = 0x1.36c2fea809651p+4 *
           exp(0x1.4933f13b56c34p+1 * lgt - 0x1.20bb1a8db7c2ap+26 / rt);
  k[r6f] = 0x1.a8aedf4000001p+35;
  k[r6b] = 0x1.d21fe02ba4021p+41 *
           exp(-0x1.c4180de8138p-5 * lgt - 0x1.608f56c942f31p+29 / rt);
  k[r7f] = 0x1.2a05f20000001p+36;
  k[r7b] = 0x1.c751834c85fa3p+42 *
           exp(-0x1.61f62f98b79p-2 * lgt - 0x1.6d9cdfa458efep+28 / rt);
  k[r8f] = 0x1.bf08eb0000001p+33;
  k[r8b] = 0x1.5f9855d9d7805p+35 *
           exp(0x1.6663a74ea76p-3 * lgt - 0x1.77bea3115d2abp+29 / rt);
  k[r9f] = 0x1.bf08eb0000001p+33;
  k[r9b] = 0x1.f1f46f53d2e95p+39 *
           exp(-0x1.a1be2be72bfp-2 * lgt - 0x1.91300c63f59aap+28 / rt);
  k[r10f] = 0x1.78ffd74000001p+35;
  k[r10b] = 0x1.8733f5670c45p+42 *
            exp(-0x1.1a6113b128ap-2 * lgt - 0x1.1380e4afade44p+28 / rt);
  k[r11f] = 0x1.f20c000000001p+19 * exp(0x1.8p+0 * lgt - 0x1.128cf6p+25 / rt);
  k[r11b] = 0x1.5c3677f2c4f5dp+8 *
            exp(0x1.f9fbb51c59f96p+0 * lgt - 0x1.22718d8c6852ep+24 / rt);
  kTroe0 = 0x1.1f0e540000002p+29 * exp(-0x1.7f17ep+23 / rt);
  kTroeInf = 0x1.12a8800000001p+24 * exp(-0x1.308f3p+23 / rt);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r12f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM3]);
  kTroe0 = 0x1.9a2fcc4021454p+58 *
           exp(-0x1.ed3bba4f794p-1 * lgt - 0x1.04f3d8104f9fep+29 / rt);
  kTroeInf = 0x1.887893cd182dcp+53 *
             exp(-0x1.ed3bba4f794p-1 * lgt - 0x1.03b9b5504f9fep+29 / rt);
  k[r12b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM3]);
  k[r13f] = 0x1.bf08eb0000001p+34;
  k[r13b] = 0x1.7773114d2cbd2p+29 *
            exp(0x1.1884380a6e84p-1 * lgt - 0x1.56c785ea73e1cp+28 / rt);
  k[r14f] = 0x1.bf08eb0000001p+34;
  k[r14b] = 0x1.ee00c754f3017p+51 *
            exp(-0x1.7615b633ec6cp-1 * lgt - 0x1.bf34f0e2b2ad6p+28 / rt);
  k[r15f] = 0x1.229298c000001p+35 * exp(-0x1.c40dp+23 / rt);
  k[r15b] = 0x1.26b460776a09p+25 *
            exp(0x1.a855fb2a7faap-2 * lgt - 0x1.0c98f0ff0201p+26 / rt);
  k[r16f] = 0x1.2a05f20000001p+33;
  k[r16b] = 0x1.097902eba39fbp+27 *
            exp(0x1.64a0b277e358p-1 * lgt - 0x1.206df200a650ap+28 / rt);
  k[r17f] = 0x1.2a05f20000001p+33;
  k[r17b] = 0x1.e03ad3b05e86dp+23 *
            exp(0x1.4661b31bd394p-1 * lgt - 0x1.3ba2a7e0d30c3p+28 / rt);
  k[r18f] = 0x1.8400000000001p+8 * exp(0x1.4p+1 * lgt - 0x1.8bdd68p+23 / rt);
  k[r18b] = 0x1.b9e0f58b91e98p+0 *
            exp(0x1.5685a13c764c6p+1 * lgt - 0x1.ef5ac345a1253p+24 / rt);
  k[r19f] = 0x1.0400000000001p+7 * exp(0x1.4p+1 * lgt - 0x1.3f3e9p+24 / rt);
  k[r19b] = 0x1.475fff740e68ap+2 *
            exp(0x1.5e1561137a3b7p+1 * lgt - 0x1.6abe8285aaf07p+23 / rt);
  k[r20f] = 0x1.74876e8000001p+35;
  k[r20b] = 0x1.647458ce9dd0bp+24 *
            exp(0x1.033cb1ceb8a8p+0 * lgt - 0x1.368967b14a19ep+28 / rt);
  k[r21f] =
      0x1.a5e0000000001p+13 * exp(0x1p+1 * lgt - 0x1.e540880000001p+22 / rt);
  k[r21b] = 0x1.4547588288b15p+15 *
            exp(0x1.d1675210992ep+0 * lgt - 0x1.566cdafbb54d2p+26 / rt);
  k[r22f] = 0x1.46d9833756001p+55 *
            exp(-0x1.68f5c28f5c28fp+0 * lgt - 0x1.ce1ad4p+26 / rt);
  k[r22b] = 0x1.106c8af8917e6p+43 *
            exp(-0x1.b271b512aa0f2p-1 * lgt + 0x1.357affdf68d8cp+23 / rt);
  k[r23f] =
      0x1.b1c0000000001p+12 * exp(0x1p+1 * lgt - 0x1.e540880000001p+22 / rt);
  k[r23b] = 0x1.27e5f8031a2bep-6 *
            exp(0x1.a88a24a33a2e8p+1 * lgt - 0x1.7b4c49ac0fc71p+27 / rt);
  k[r24f] = 0x1.bf08eb0000001p+34;
  k[r24b] = 0x1.cf31dde898758p+46 *
            exp(-0x1.33a0a823b408p-1 * lgt - 0x1.698393d05e103p+28 / rt);
  k[r25f] = 0x1.86a0000000001p+13 *
            exp(0x1.d47ae147ae148p+0 * lgt - 0x1.c17f2p+19 / rt);
  k[r25b] = 0x1.280212bf1cc23p-3 *
            exp(0x1.5ab4699b36cdcp+1 * lgt - 0x1.9d6c71ddbd07dp+26 / rt);
  k[r26f] = 0x1.4dc9380000001p+34;
  k[r26b] = 0x1.28dfd0297f633p+24 *
            exp(0x1.cd2be798f93cp-1 * lgt - 0x1.36a264caeb805p+28 / rt);
  k[r27f] = 0x1.5ec8000000001p+16 *
            exp(0x1.eb851eb851eb8p+0 * lgt - 0x1.6b4d14p+24 / rt);
  k[r27b] = 0x1.4184b5e7a4972p+3 *
            exp(0x1.38d8d1daa3be8p+1 * lgt - 0x1.a125d0c931565p+24 / rt);
  k[r28f] = 0x1.74876e8000001p+36;
  k[r28b] = 0x1.856e2310e9f04p+10 *
            exp(0x1.62c26d5d7414p+0 * lgt - 0x1.96edcd7d86eccp+28 / rt);
  k[r29f] = 0x1.2a05f20000001p+33 * exp(-0x1.feca800000001p+24 / rt);
  k[r29b] = 0x1.163feaeac4e4ap+22 *
            exp(0x1.1fcc6ebad593p-1 * lgt - 0x1.be153351c73ddp+23 / rt);
  k[r30f] = 0x1.a13b860000001p+30 * exp(-0x1.58c8d8p+22 / rt);
  k[r30b] = 0x1.c57e95a577cccp+21 *
            exp(0x1.908c6ee830fp-1 * lgt - 0x1.87befc13cf251p+27 / rt);
  k[r31f] = 0x1.2a05f20000001p+31 * exp(-0x1.7d7f93p+27 / rt);
  k[r31b] = 0x1.8e188ccf4e6dfp+44 *
            exp(-0x1.aef5f5323a64p-1 * lgt - 0x1.c3d8150289482p+27 / rt);
  k[r32f] = 0x1.74876e8000001p+36 * exp(-0x1.3f3e9p+27 / rt);
  k[r32b] = 0x1.e9689a9c974aap+29 *
            exp(0x1.65e8b9037cdd4p-4 * lgt - 0x1.29a84622d8b3bp+21 / rt);
  k[r33f] = 0x1.45f680b000002p+41 * exp(-0x1.b851eb851eb85p-1 * lgt);
  k[r33b] = 0x1.ca8c2b3218a0bp+51 *
            exp(-0x1.be619e0b11a8p-1 * lgt - 0x1.865f74e9cb3eap+27 / rt);
  k[r34f] = 0x1.2eae09c800002p+44 * exp(-0x1.3d70a3d70a3d7p+0 * lgt);
  k[r34b] = 0x1.a9cb4cae853c2p+54 *
            exp(-0x1.40787d1a03c2p+0 * lgt - 0x1.865f74e9cb3f6p+27 / rt);
  k[r35f] = 0x1.47b5899b00002p+43 * exp(-0x1.851eb851eb852p-1 * lgt);
  k[r35b] = 0x1.cd0108af2b494p+53 *
            exp(-0x1.8b2e6ad7de9ep-1 * lgt - 0x1.865f74e9cb3fep+27 / rt);
  k[r36f] = 0x1.7a598c3a00003p+44 * exp(-0x1.3d70a3d70a3d7p+0 * lgt);
  k[r36b] = 0x1.0a1f0fed12dd6p+55 *
            exp(-0x1.40787d1a03b5p+0 * lgt - 0x1.865f74e9cb3ecp+27 / rt);
  k[r38f] = 0x1.81a0316280001p+44 *
            exp(-0x1.5765fd8adab9fp-1 * lgt - 0x1.100328p+26 / rt);
  k[r38b] = 0x1.877e8e8674866p+35 *
            exp(-0x1.df0811faf72f2p-3 * lgt + 0x1.405ff6fa2b679p+20 / rt);
  k[r39f] = 0x1.d1a94a2000003p+39 * exp(-0x1p+0 * lgt);
  k[r39b] = 0x1.a92848d40e53cp+47 *
            exp(-0x1.4cae04544d38p-1 * lgt - 0x1.9d0739fcb14d2p+28 / rt);
  k[r40f] = 0x1.4f46b04000002p+36 * exp(-0x1.3333333333333p-1 * lgt);
  k[r40b] = 0x1.321d013c85e32p+44 *
            exp(-0x1.ff84de1e022p-3 * lgt - 0x1.9d0739fcb14ddp+28 / rt);
  k[r41f] = 0x1.b48eb57e00003p+45 * exp(-0x1.4p+0 * lgt);
  k[r41b] = 0x1.8e95c446cd068p+53 *
            exp(-0x1.ccae04544d24p-1 * lgt - 0x1.9d0739fcb14cfp+28 / rt);
  k[r42f] = 0x1.f438daa060004p+48 * exp(-0x1p+1 * lgt);
  k[r42b] = 0x1.c8b6463bca9d4p+56 *
            exp(-0x1.a657022a268cp+0 * lgt - 0x1.9d0739fcb14cfp+28 / rt);
  k[r43f] = 0x1.38a388a43c002p+54 * exp(-0x1p+1 * lgt);
  k[r43b] = 0x1.8bafdb5867d32p+65 *
            exp(-0x1.c12d7d1f3508p+0 * lgt - 0x1.daddfc03e033cp+28 / rt);
  k[r44f] = 0x1.d942c90000001p+31 * exp(-0x1.56bd9p+21 / rt);
  k[r44b] = 0x1.b04484e0960aap+23 *
            exp(0x1.6358b153a442p-1 * lgt - 0x1.aa35257000cbdp+27 / rt);
  k[r45f] = 0x1.4dc9380000001p+35 * exp(-0x1.10c33p+22 / rt);
  k[r45b] = 0x1.b145db8ebfc95p+32 *
            exp(0x1.72c35c634b68p-2 * lgt - 0x1.bc35188f975cep+27 / rt);
  k[r46f] = 0x1.38eca48000001p+36 * exp(-0x1.445a6p+21 / rt);
  k[r46b] = 0x1.ea830edc45572p+23 *
            exp(0x1.8711df80ef3ep-1 * lgt - 0x1.1f32acf9014ep+27 / rt);
  k[r47f] = 0x1.7a20000000001p+13 * exp(0x1p+1 * lgt - 0x1.4c03b4p+24 / rt);
  k[r47b] = 0x1.4844be3895ee2p+5 *
            exp(0x1.4db0e32a8b3p+1 * lgt - 0x1.51f986defb27ap+26 / rt);
  k[r48f] = 0x1.2a05f20000001p+33 * exp(-0x1.cbb6b8p+23 / rt);
  k[r48b] = 0x1.19185e7c13e78p+16 *
            exp(0x1.44143b207fb6p+0 * lgt - 0x1.19033ffb6e555p+28 / rt);
  k[r49f] = 0x1.3356219000001p+37;
  k[r49b] = 0x1.397b640e65027p+36 *
            exp(0x1.e017c9cd2fc4p-3 * lgt - 0x1.73898965d6547p+26 / rt);
  kTroe0 = 0x1.68d28e3f00282p+66 *
           exp(-0x1.6147ae147ae14p+1 * lgt - 0x1.98a228p+22 / rt);
  kTroeInf = 0x1.176592e000001p+39;
  fcTroe = 0x1.c083126e978d5p-2 * exp(-temp / 0x1.6cp+6) +
           0x1.1fbe76c8b4396p-1 * exp(-temp / 0x1.6ccp+12) +
           0x1p+0 * exp(-0x1.0b4p+13 / temp);
  k[r50f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.1c6d7e2ff3da8p+83 *
           exp(-0x1.76ed83fb043fp+1 * lgt - 0x1.c0fa357d4b006p+28 / rt);
  kTroeInf = 0x1.b87bc99ca4de3p+55 *
             exp(-0x1.5a5d5e6895dcp-3 * lgt - 0x1.ba97acdd4b006p+28 / rt);
  k[r50b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r51f] = 0x1.bf08eb0000001p+34;
  k[r51b] = 0x1.405606466019cp+30 *
            exp(0x1.d769aac8ace2p-3 * lgt - 0x1.72f4c481a3785p+25 / rt);
  kTroe0 = 0x1.0ee6d2930f2a2p+91 *
           exp(-0x1.30a3d70a3d70ap+2 * lgt - 0x1.3794d8p+23 / rt);
  kTroeInf = 0x1.948b11ff00001p+43 *
             exp(-0x1.116872b020c4ap-1 * lgt - 0x1.11c8cp+21 / rt);
  fcTroe = 0x1.bc6a7ef9db22dp-3 * exp(-temp / 0x1.28p+6) +
           0x1.90e5604189375p-1 * exp(-temp / 0x1.6fap+11) +
           0x1p+0 * exp(-0x1.b34p+12 / temp);
  k[r52f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.a4cc63551db24p+109 *
           exp(-0x1.3af6fdd377b2p+2 * lgt - 0x1.af68b2cf99f29p+28 / rt);
  kTroeInf = 0x1.3a3191b6734d5p+62 *
             exp(-0x1.6401a8f9f2cfap-1 * lgt - 0x1.a7cf9d8f99f29p+28 / rt);
  k[r52b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r53f] = 0x1.4244000000001p+19 *
            exp(0x1.9eb851eb851ecp+0 * lgt - 0x1.5a0f9ap+25 / rt);
  k[r53b] = 0x1.7ad839423f38p+8 *
            exp(0x1.10d6f5732408ap+1 * lgt - 0x1.14e90968baf17p+25 / rt);
  kTroe0 = 0x1.12399f4e99b82p+61 *
           exp(-0x1.48f5c28f5c28fp+1 * lgt - 0x1.b22c200000001p+20 / rt);
  kTroeInf = 0x1.03e0520000001p+30 *
             exp(0x1.eb851eb851eb8p-2 * lgt + 0x1.099c3ffffffffp+20 / rt);
  fcTroe = 0x1.bda5119ce075fp-3 * exp(-temp / 0x1.0fp+8) +
           0x1.9096bb98c7e28p-1 * exp(-temp / 0x1.586p+11) +
           0x1p+0 * exp(-0x1.9aap+12 / temp);
  k[r54f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.25a343956f966p+78 *
           exp(-0x1.55a8f4f8f3c8p+1 * lgt - 0x1.622dde089ff5ap+28 / rt);
  kTroeInf = 0x1.1645ece4e8748p+47 *
             exp(0x1.85eb8b6b94f3p-2 * lgt - 0x1.5f7215a89ff5ap+28 / rt);
  k[r54b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r55f] = 0x1.116fb1e000001p+36;
  k[r55b] = 0x1.82228565b4689p+31 *
            exp(0x1.2a77ffc74018p-1 * lgt - 0x1.5e4d39bec4bafp+28 / rt);
  kTroe0 = 0x1.a4352f2e5ea7cp+86 *
           exp(-0x1.347ae147ae148p+2 * lgt - 0x1.a0ef28p+24 / rt);
  kTroeInf = 0x1.017df80000001p+29 *
             exp(0x1.d0e5604189375p-2 * lgt - 0x1.cbb6b8p+23 / rt);
  fcTroe = 0x1.200d1b71758e2p-2 * exp(-temp / 0x1.9cp+6) +
           0x1.6ff972474538fp-1 * exp(-temp / 0x1.42cp+10) +
           0x1p+0 * exp(-0x1.04p+12 / temp);
  k[r56f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.00276e91cd8b9p+100 *
           exp(-0x1.4ce33118ce496p+2 * lgt - 0x1.1e450d4f744b2p+27 / rt);
  kTroeInf = 0x1.39edbd6608e96p+42 *
             exp(0x1.29818cbe17a54p-4 * lgt - 0x1.06e293cf744b2p+27 / rt);
  k[r56b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.d1de3d2d5c716p+80 *
           exp(-0x1.3333333333333p+2 * lgt - 0x1.630014p+24 / rt);
  kTroeInf = 0x1.017df80000001p+29 *
             exp(0x1.d0e5604189375p-2 * lgt - 0x1.4c035p+23 / rt);
  fcTroe = 0x1.ef9db22d0e56p-3 * exp(-temp / 0x1.78p+6) +
           0x1.84189374bc6a8p-1 * exp(-temp / 0x1.84cp+10) +
           0x1p+0 * exp(-0x1.068p+12 / temp);
  k[r57f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.39faad856763cp+97 *
           exp(-0x1.47d3a318d1786p+2 * lgt - 0x1.c03b7e1e35a81p+26 / rt);
  kTroeInf = 0x1.5b14ebe0988afp+45 *
             exp(0x1.0dbcc3cf49c8ap-3 * lgt - 0x1.90fbe31e35a81p+26 / rt);
  k[r57b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r58f] = 0x1.c070000000002p+15 *
            exp(0x1.e666666666666p+0 * lgt - 0x1.5e257ffffffffp+23 / rt);
  k[r58b] = 0x1.7e5d230452b67p+6 *
            exp(0x1.2cbae487b783cp+1 * lgt - 0x1.1df2d0504565cp+26 / rt);
  kTroe0 = 0x1.208545e35d9d9p+85 *
           exp(-0x1.299999999999ap+2 * lgt - 0x1.445a6p+24 / rt);
  kTroeInf = 0x1.f7102e0000002p+29 * exp(0x1p-1 * lgt - 0x1.5f6ccp+18 / rt);
  fcTroe = 0x1.999999999999ap-2 * exp(-temp / 0x1.9p+6) +
           0x1.3333333333333p-1 * exp(-temp / 0x1.5f9p+16) +
           0x1p+0 * exp(-0x1.388p+13 / temp);
  k[r59f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.1321fd6f1072ep+100 *
           exp(-0x1.20b0a3b9f882p+2 * lgt - 0x1.97306b340663fp+28 / rt);
  kTroeInf = 0x1.dfb86274ad164p+44 *
             exp(0x1.4747aefd08bdp-1 * lgt - 0x1.8342a0640663fp+28 / rt);
  k[r59b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r60f] = 0x1.2a05f20000001p+34;
  k[r60b] = 0x1.be5e12d53f49ap+28 *
            exp(0x1.76947a34b4cp-1 * lgt - 0x1.27f3a5d4f72a2p+28 / rt);
  k[r61f] = 0x1.3ab6680000001p+27 *
            exp(0x1.4cccccccccccdp-1 * lgt + 0x1.22212p+20 / rt);
  k[r61b] = 0x1.0e29282ed11c3p+14 *
            exp(0x1.9f4f048ea21edp+0 * lgt - 0x1.795d861f0da0cp+23 / rt);
  k[r62f] = 0x1.e8c2120000001p+34 *
            exp(-0x1.70a3d70a3d70ap-4 * lgt - 0x1.3795500000001p+21 / rt);
  k[r62b] = 0x1.cdf9d13a06c47p+16 *
            exp(0x1.5bee8e79130e1p+0 * lgt - 0x1.81eb4afe1e9b7p+23 / rt);
  kTroe0 = 0x1.66fe4ba326d82p+118 *
           exp(-0x1.dc28f5c28f5c3p+2 * lgt - 0x1.c17f18p+25 / rt);
  kTroeInf = 0x1.21adb70000001p+31 *
             exp(0x1.07ae147ae147bp-1 * lgt - 0x1.98a2p+17 / rt);
  fcTroe = 0x1.3333333333333p-2 * exp(-temp / 0x1.9p+6) +
           0x1.6666666666666p-1 * exp(-temp / 0x1.5f9p+16) +
           0x1p+0 * exp(-0x1.388p+13 / temp);
  k[r63f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.35a2c21c12004p+130 *
           exp(-0x1.d707dfce701ap+2 * lgt - 0x1.d64f5e143320ap+28 / rt);
  kTroeInf = 0x1.f3b3b7c1100f2p+42 *
             exp(0x1.30b6c41bdb593p-1 * lgt - 0x1.9e528f543320ap+28 / rt);
  k[r63b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r64f] = 0x1.4438000000001p+15 *
            exp(0x1.a147ae147ae14p+0 * lgt - 0x1.eb6197fffffffp+22 / rt);
  k[r64b] = 0x1.253fe2e03918bp+12 *
            exp(0x1.92282e6673274p+0 * lgt - 0x1.1711e20165dfcp+25 / rt);
  k[r65f] = 0x1.2a05f20000001p+34;
  k[r65b] = 0x1.93bb0f4e7d51bp+25 *
            exp(0x1.58557ad8a5a8p-1 * lgt - 0x1.43285bb523e3ep+28 / rt);
  k[r66f] =
      0x1.65a0bc0000001p+30 * exp(0x1p-1 * lgt + 0x1.c17f000000001p+18 / rt);
  k[r66b] = 0x1.15ad236d018e5p+14 *
            exp(0x1.69c91e7a33f6p+0 * lgt - 0x1.3d8b1b8929423p+25 / rt);
  k[r67f] = 0x1.e80355e000001p+37 *
            exp(-0x1.d70a3d70a3d71p-3 * lgt - 0x1.1145f8p+22 / rt);
  k[r67b] = 0x1.a13611d1540e5p+16 *
            exp(0x1.28f8048d9a9e4p+0 * lgt - 0x1.48cfebc0ed839p+25 / rt);
  k[r68f] = 0x1.09a0000000001p+14 *
            exp(0x1.0cccccccccccdp+1 * lgt - 0x1.36f19p+24 / rt);
  k[r68b] = 0x1.fca496cbcc218p+6 *
            exp(0x1.27cf5ff87786cp+1 * lgt - 0x1.6c5c6e45575ebp+25 / rt);
  k[r69f] = 0x1.0680000000001p+12 *
            exp(0x1.0cccccccccccdp+1 * lgt - 0x1.36f19p+24 / rt);
  k[r69b] = 0x1.15df202dd9a9p+8 *
            exp(0x1.2f5f1fcf7b6d7p+1 * lgt - 0x1.256d7e87e3085p+24 / rt);
  kTroe0 = 0x1.83bdac6ae9bc5p+91 *
           exp(-0x1.3333333333333p+2 * lgt - 0x1.e540880000001p+22 / rt);
  kTroeInf = 0x1.6bcc41e900001p+46 * exp(-0x1p+0 * lgt);
  fcTroe = 0x1.6a161e4f765fep-2 * exp(-temp / 0x1.08p+7) +
           0x1.4af4f0d844d01p-1 * exp(-temp / 0x1.48cp+10) +
           0x1p+0 * exp(-0x1.5bep+12 / temp);
  k[r70f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.f9388451339f6p+110 *
           exp(-0x1.42f6a6b6d881p+2 * lgt - 0x1.0d248aa3ade0ep+29 / rt);
  kTroeInf = 0x1.da05fa1696238p+65 *
             exp(-0x1.3f0dce0e95374p+0 * lgt - 0x1.095a0993ade0ep+29 / rt);
  k[r70b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.d462db8b40f6dp+114 *
           exp(-0x1.d147ae147ae14p+2 * lgt - 0x1.ccfd48p+24 / rt);
  kTroeInf = 0x1.4dc9380000001p+32 * exp(-0x1.3279dp+23 / rt);
  fcTroe = 0x1.fe90ff9724745p-3 * exp(-temp / 0x1.8ap+6) +
           0x1.805bc01a36e2fp-1 * exp(-temp / 0x1.458p+10) +
           0x1p+0 * exp(-0x1.047p+12 / temp);
  k[r71f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.956bb28aecdfep+122 *
           exp(-0x1.c6478be55c7e4p+2 * lgt - 0x1.519ff3b8c301dp+27 / rt);
  kTroeInf = 0x1.20ea38eb49d8bp+40 *
             exp(0x1.600445e3cc6p-3 * lgt - 0x1.2b27e7b8c301cp+27 / rt);
  k[r71b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.287626ee5219ap+80 *
           exp(-0x1.ee147ae147ae1p+1 * lgt - 0x1.a7f4dp+23 / rt);
  kTroeInf = 0x1.6a65700000001p+32 *
             exp(0x1.147ae147ae148p-2 * lgt - 0x1.1e0adffffffffp+20 / rt);
  fcTroe = 0x1.be76c8b439581p-3 * exp(-temp / 0x1.9fp+7) +
           0x1.90624dd2f1aap-1 * exp(-temp / 0x1.4cep+11) +
           0x1p+0 * exp(-0x1.7cfp+12 / temp);
  k[r72f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.3a4a51f727b98p+99 *
           exp(-0x1.ff70ea811698p+1 * lgt - 0x1.ca27f40d6eb5ep+28 / rt);
  kTroeInf = 0x1.8030b3c3f9fa9p+51 *
             exp(0x1.132ec8926d8ap-3 * lgt - 0x1.be06586d6eb5ep+28 / rt);
  k[r72b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r73f] = 0x1.bf08eb0000001p+34;
  k[r73b] = 0x1.d78a63f5178c6p+34 *
            exp(0x1.6d43a8cb00fp-3 * lgt - 0x1.110714a04fce6p+28 / rt);
  kTroe0 = 0x1.ce3922c2af446p+118 *
           exp(-0x1.e7ae147ae147bp+2 * lgt - 0x1.bd06f4p+24 / rt);
  kTroeInf = 0x1.017df80000001p+29 *
             exp(0x1.d0e5604189375p-2 * lgt - 0x1.d0d1e80000001p+22 / rt);
  fcTroe = 0x1.94af4f0d844dp-6 * exp(-temp / 0x1.a4p+7) +
           0x1.f35a858793dd9p-1 * exp(-temp / 0x1.ecp+9) +
           0x1p+0 * exp(-0x1.116p+12 / temp);
  k[r74f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.a5b133a8b84d8p+129 *
           exp(-0x1.ef71ae271c43cp+2 * lgt - 0x1.5848328a476fbp+27 / rt);
  kTroeInf = 0x1.d5d3a1bc76d01p+39 *
             exp(0x1.54abc57dd9765p-2 * lgt - 0x1.2f2de34a476fcp+27 / rt);
  k[r74b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r75f] = 0x1.4b40000000001p+10 *
            exp(0x1.43d70a3d70a3dp+1 * lgt - 0x1.86c166p+25 / rt);
  k[r75b] = 0x1.1d47eddeaefb8p-1 *
            exp(0x1.8207f8c82c923p+1 * lgt - 0x1.0f7192f429b64p+24 / rt);
  kTroe0 = 0x1.329ba8fa506bbp+117 *
           exp(-0x1.c51eb851eb852p+2 * lgt - 0x1.aad45cp+24 / rt);
  kTroeInf = 0x1.d9d8c3ed90001p+48 *
             exp(-0x1.fae147ae147aep-1 * lgt - 0x1.9386800000001p+22 / rt);
  fcTroe = 0x1.432ca57a786c2p-3 * exp(-temp / 0x1.f4p+6) +
           0x1.af34d6a161e4fp-1 * exp(-temp / 0x1.156p+11) +
           0x1p+0 * exp(-0x1.ae2p+12 / temp);
  k[r76f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.6b49b64a6dfbfp+137 *
           exp(-0x1.d27e13134ca38p+2 * lgt - 0x1.acd1401bcd61fp+28 / rt);
  kTroeInf = 0x1.18b8a5967d05ap+69 *
             exp(-0x1.32ee0edc8eb6fp+0 * lgt - 0x1.9872145bcd61fp+28 / rt);
  k[r76b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r77f] = 0x1.dcd6500000001p+30;
  k[r77b] = 0x1.dd34b74fbc571p+27 *
            exp(0x1.e2dd921b167p-2 * lgt - 0x1.0cb38ff78d972p+28 / rt);
  k[r78f] = 0x1.c138000000001p+16 *
            exp(0x1.e666666666666p+0 * lgt - 0x1.e0c8780000001p+24 / rt);
  k[r78b] = 0x1.5a278f0324df5p+4 *
            exp(0x1.3ac667a0e2712p+1 * lgt - 0x1.477e39071f745p+25 / rt);
  k[r79f] = 0x1.74876e8000001p+36;
  k[r79b] = 0x1.c4127d49747b6p+16 *
            exp(0x1.8f9ef6497844p+0 * lgt - 0x1.11df055df797ep+26 / rt);
  k[r80f] = 0x1.74876e8000001p+35 * exp(-0x1.feca800000001p+24 / rt);
  k[r80b] = 0x1.24681f5a7f7fbp+25 *
            exp(0x1.31c03677a7244p-1 * lgt - 0x1.5765d6edf12fdp+24 / rt);
  k[r81f] = 0x1.50c4288000001p+33 * exp(-0x1.b5bf47fffffffp+23 / rt);
  k[r81b] = 0x1.93d88e932a5cdp+11 *
            exp(0x1.939868cec26ap+0 * lgt - 0x1.11a6214d26a9ap+27 / rt);
  k[r82f] = 0x1.2a05f20000001p+33;
  k[r82b] = 0x1.62bdb2165b37p+28 *
            exp(0x1.a4dec7796eap-2 * lgt - 0x1.dfd95b4ff7eb5p+26 / rt);
  kTroe0 = 0x1.12d86259fd1e9p+72 *
           exp(-0x1.b5c28f5c28f5cp+1 * lgt - 0x1.509a69p+28 / rt);
  kTroeInf = 0x1.4ff0000000001p+15 * exp(0x1.8p+0 * lgt - 0x1.3da61d8p+28 / rt);
  fcTroe = 0x1.16872b020c49cp-4 * exp(-temp / 0x1.8ap+7) +
           0x1.dd2f1a9fbe76dp-1 * exp(-temp / 0x1.81p+10) +
           0x1p+0 * exp(-0x1.41ep+13 / temp);
  k[r83f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.a0cfd9608b786p+93 *
           exp(-0x1.0689e0dbc72d8p+2 * lgt - 0x1.52c8e129db48cp+28 / rt);
  kTroeInf = 0x1.fd75d7f92f165p+36 *
             exp(0x1.a2bb36926a6bp-1 * lgt - 0x1.3fd495a9db48cp+28 / rt);
  k[r83b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r84f] = 0x1.a5e0000000001p+17 *
            exp(0x1.828f5c28f5c29p+0 * lgt - 0x1.b600e8p+23 / rt);
  k[r84b] = 0x1.246761a9043f7p+21 *
            exp(0x1.67b8e133e7aep+0 * lgt - 0x1.2e1b251cbb992p+26 / rt);
  kTroe0 = 0x1.0bc1576c00002p+41 *
           exp(-0x1.ccccccccccccdp-1 * lgt + 0x1.b22c48p+22 / rt);
  kTroeInf = 0x1.13abe64000001p+36 * exp(-0x1.7ae147ae147aep-2 * lgt);
  fcTroe = 0x1.0fc504816f007p-2 * exp(-temp / 0x1.78p+6) +
           0x1.781d7dbf487fdp-1 * exp(-temp / 0x1.b7p+10) +
           0x1p+0 * exp(-0x1.43ep+12 / temp);
  k[r85f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.6749e579ebf5cp+69 *
           exp(-0x1.eba81ea61af3p+0 * lgt - 0x1.92df8150e3bc8p+27 / rt);
  kTroeInf = 0x1.71e93e2e4714fp+64 *
             exp(-0x1.63fa0a2b39ab6p+0 * lgt - 0x1.a070e390e3bc9p+27 / rt);
  k[r85b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r86f] = 0x1.1d9999999999ap+5 *
            exp(0x1.3333333333333p+1 * lgt + 0x1.0d7144p+23 / rt);
  k[r86b] = 0x1.4cd5fcd0cd0fp+9 *
            exp(0x1.2a44e7a7e068p+1 * lgt - 0x1.e7875ddbfdf37p+25 / rt);
  k[r87f] = 0x1.b022388000001p+33 * exp(0x1.feca800000001p+20 / rt);
  k[r87b] = 0x1.84c93f9b47f12p+34 *
            exp(0x1.0769708f1298p-2 * lgt - 0x1.15af770efa92fp+28 / rt);
  k[r88f] = 0x1.dcd6500000001p+30 * exp(-0x1.b4374p+20 / rt);
  k[r88b] = 0x1.1eebece1e36afp+26 *
            exp(0x1.011696c00fd8p-1 * lgt - 0x1.fd247efbb6b93p+26 / rt);
  k[r89f] = 0x1.8289060790001p+50 * exp(-0x1.d572b6p+26 / rt);
  k[r89b] = 0x1.d12baa398f24cp+45 *
            exp(0x1.011696c012d8p-1 * lgt - 0x1.e5e32bfddb67p+27 / rt);
  k[r90f] = 0x1.74876e8000001p+35;
  k[r90b] = 0x1.51024087844acp+43 *
            exp(-0x1.04a75729f6dp-2 * lgt - 0x1.35e0ff86b094dp+29 / rt);
  k[r91f] = 0x1.bf08eb0000001p+34;
  k[r91b] = 0x1.241abe6af109fp+46 *
            exp(-0x1.34c5b8e8ef5cp-1 * lgt - 0x1.6a5727a81203cp+28 / rt);
  k[r92f] = 0x1.2a05f20000001p+34;
  k[r92b] = 0x1.c0ef3e46251b9p+50 *
            exp(-0x1.852615619b38p-1 * lgt - 0x1.38970b64986ffp+28 / rt);
  k[r93f] = 0x1.6120000000001p+13 * exp(0x1p+1 * lgt - 0x1.7f17ep+23 / rt);
  k[r93b] = 0x1.e11a8b350de3cp+12 *
            exp(0x1.18045cbbd24dcp+1 * lgt - 0x1.526bb35f1aaa4p+26 / rt);
  k[r94f] = 0x1.bf08eb0000001p+34;
  k[r94b] = 0x1.eaf9217c946fbp+50 *
            exp(-0x1.a50a1388d558p-1 * lgt - 0x1.5c2a3824351a7p+28 / rt);
  kTroe0 = 0x1.93e5939a08cedp+101 *
           exp(-0x1.7ae147ae147aep+2 * lgt - 0x1.90f87p+23 / rt);
  kTroeInf = 0x1.3d2fafdd8c001p+51 *
             exp(-0x1.6e147ae147ae1p+0 * lgt - 0x1.53ad3p+22 / rt);
  fcTroe = 0x1.2d0e560418937p-1 * exp(-temp / 0x1.86p+7) +
           0x1.a5e353f7ced91p-2 * exp(-temp / 0x1.70cp+12) +
           0x1p+0 * exp(-0x1.8fap+12 / temp);
  k[r95f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.c0abae4489dadp+129 *
           exp(-0x1.b032795882648p+2 * lgt - 0x1.82857b630df7dp+28 / rt);
  kTroeInf = 0x1.60590d2a4d6efp+79 *
             exp(-0x1.21aca0c57faa4p+1 * lgt - 0x1.7b4c6ca30df7dp+28 / rt);
  k[r95b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r96f] = 0x1.b580000000001p+15 *
            exp(0x1.999999999999ap+0 * lgt - 0x1.5a0fccp+24 / rt);
  k[r96b] = 0x1.5f382f8540af7p+10 *
            exp(0x1.01dbe423bbe9cp+1 * lgt - 0x1.af3a5f34a99e6p+25 / rt);
  k[r97f] = 0x1.24db779e20001p+49 *
            exp(-0x1.570a3d70a3d71p+0 * lgt - 0x1.69e54ffffffffp+22 / rt);
  k[r97b] = 0x1.42753ce769d71p+44 *
            exp(-0x1.b9f41f5e51746p-1 * lgt - 0x1.2d5bd37c444abp+21 / rt);
  k[r98f] = 0x1.86a0000000001p+16 *
            exp(0x1.999999999999ap+0 * lgt - 0x1.8e6a880000001p+23 / rt);
  k[r98b] = 0x1.3e4722f08761fp+9 *
            exp(0x1.00dc5bcfa6f4p+1 * lgt - 0x1.069510d119095p+26 / rt);
  k[r99f] = 0x1.73e0000000001p+15 *
            exp(0x1.3a5e353f7ced9p+0 * lgt - 0x1.1e0bp+18 / rt);
  k[r99b] = 0x1.e94cc3ed80508p+37 *
            exp(-0x1.9dd83bf60fp-5 * lgt - 0x1.a2d3b6e0fb2c6p+26 / rt);
  k[r100f] = 0x1.74876e8000001p+35;
  k[r100b] = 0x1.6c9f17275772fp+34 *
             exp(0x1.e99613ba471p-2 * lgt - 0x1.9c23fbc5f3a0dp+28 / rt);
  k[r101f] = 0x1.a2b3800000001p+21 *
             exp(0x1.2e147ae147ae1p+0 * lgt + 0x1.c8a5ep+20 / rt);
  k[r101b] = 0x1.eee3bb6d2a179p+15 *
             exp(0x1.864d629541f9p+0 * lgt - 0x1.e26690ed00f7ap+26 / rt);
  k[r102f] = 0x1.2a05f20000001p+32;
  k[r102b] = 0x1.356134070f78cp+30 *
             exp(0x1.40e7844a988p-1 * lgt - 0x1.65ca67dc260e9p+28 / rt);
  k[r103f] = 0x1.2a05f20000001p+32;
  k[r103b] = 0x1.17d3f2cb80c09p+27 *
             exp(0x1.22a884ee88b4p-1 * lgt - 0x1.80ff1dbc52cabp+28 / rt);
  k[r104f] = 0x1.6800000000001p+10 * exp(0x1p+1 * lgt + 0x1.ad10ap+21 / rt);
  k[r104b] = 0x1.ddccb6665437bp+6 *
             exp(0x1.0d9755b1237b4p+1 * lgt - 0x1.5264563f67401p+26 / rt);
  k[r105f] = 0x1.89c0000000001p+12 * exp(0x1p+1 * lgt - 0x1.7f1808p+22 / rt);
  k[r105b] = 0x1.20e4353b24277p+12 *
             exp(0x1.15271588276f4p+1 * lgt - 0x1.0aeb843eb4547p+26 / rt);
  k[r106f] = 0x1.2a05f20000001p+34;
  k[r106b] = 0x1.13b1399c7476ap+48 *
             exp(-0x1.7cab2bead90ep-1 * lgt - 0x1.9671833bd13d6p+27 / rt);
  k[r107f] =
      0x1.d426c47622577p-23 * exp(0x1.2p+2 * lgt + 0x1.feca800000001p+21 / rt);
  k[r107b] = 0x1.829c4889eae25p-10 *
             exp(0x1.e0c08d5997484p+1 * lgt - 0x1.701278117c6bdp+26 / rt);
  k[r108f] = 0x1.f800000000002p+8 *
             exp(0x1.2666666666666p+1 * lgt - 0x1.aefb0ep+25 / rt);
  k[r108b] = 0x1.5dab23a0a71d8p+26 *
             exp(0x1.251635a1a06d4p+0 * lgt - 0x1.deb3df0612301p+24 / rt);
  k[r109f] = 0x1.0748p+15 * exp(0x1p+1 * lgt - 0x1.bef162p+25 / rt);
  k[r109b] = 0x1.ff7807112dbd8p+6 *
             exp(0x1.3ef02877b072cp+1 * lgt - 0x1.02a3908ea84bbp+15 / rt);
  k[r110f] =
      0x1.034f03b80a36ep-21 * exp(0x1p+2 * lgt + 0x1.fecaa7fffffffp+22 / rt);
  k[r110b] = 0x1.00cc4179329a7p-30 *
             exp(0x1.354660e07c37cp+2 * lgt - 0x1.a6583d95e4dc1p+27 / rt);
  k[r111f] = 0x1.2a05f20000001p+32;
  k[r111b] = 0x1.b3c50f01f3517p+35 *
             exp(0x1.2d1fa24520ap-4 * lgt - 0x1.4eddd6a77eb29p+28 / rt);
  k[r112f] = 0x1.c200000000001p+11 * exp(0x1p+1 * lgt - 0x1.3f3e9p+23 / rt);
  k[r112b] = 0x1.0c9d31a6e84a8p+4 *
             exp(0x1.30c5b110349p+1 * lgt - 0x1.3f7d17b38bf8ep+25 / rt);
  k[r113f] = 0x1.ba80000000001p+11 *
             exp(0x1.0f5c28f5c28f6p+1 * lgt - 0x1.bc637p+21 / rt);
  k[r113b] = 0x1.d8aa95b971021p+2 *
             exp(0x1.49841fe8ea864p+1 * lgt - 0x1.30cb22204b4aap+26 / rt);
  k[r114f] = 0x1.bf08eb0000001p+32 * exp(-0x1.fecaa7fffffffp+22 / rt);
  k[r114b] = 0x1.e667a2ba0d87bp+25 *
             exp(0x1.f826811b143cp-2 * lgt - 0x1.dadd10b06fb95p+25 / rt);
  k[r115f] = 0x1.efe9200000001p+26 * exp(0x1.a04b9p+22 / rt);
  k[r115b] = 0x1.72be35065728fp+32 *
             exp(-0x1.f58779e21c1p-3 * lgt - 0x1.273055a019cc3p+27 / rt);
  k[r116f] = 0x1.8727cda000001p+38 * exp(-0x1.7f17ep+25 / rt);
  k[r116b] = 0x1.246d910e09168p+44 *
             exp(-0x1.f58779e2164p-3 * lgt - 0x1.93f8aa2019d0ap+27 / rt);
  k[r117f] = 0x1.2a05f20000001p+34;
  k[r117b] = 0x1.5fdac5cc275cbp+38 *
             exp(0x1.ebca1f544p-9 * lgt - 0x1.c5a7ad211915dp+28 / rt);
  k[r118f] = 0x1.dcd6500000001p+29;
  k[r118b] = 0x1.0742c8f2256c9p+38 *
             exp(-0x1.32260f0f7d4p-3 * lgt - 0x1.c4f8a33568a54p+27 / rt);
  k[r119f] = 0x1.19a1c74000001p+35;
  k[r119b] = 0x1.1c49123840b02p+39 *
             exp(-0x1.320175cde3bcp-3 * lgt - 0x1.93bb7a2d6dfd2p+26 / rt);
  k[r120f] = 0x1.176592e000001p+37 * exp(-0x1.78b5670000001p+26 / rt);
  k[r120b] = 0x1.201fc934a09c6p+47 *
             exp(-0x1.07880ebd6914p-1 * lgt - 0x1.53ab6674bf75fp+28 / rt);
  k[r121f] = 0x1.5e00000000001p+12 * exp(0x1p+1 * lgt - 0x1.7f17ep+25 / rt);
  k[r121b] = 0x1.57c1493496cc6p+11 *
             exp(0x1.d7ad9c53f2f58p+0 * lgt - 0x1.6582ece29493cp+25 / rt);
  k[r122f] = 0x1.b022388000001p+35 * exp(-0x1.26376p+21 / rt);
  k[r122b] = 0x1.8ce1297245492p+34 *
             exp(0x1.754135dc826p-3 * lgt - 0x1.1466a1eb337cap+29 / rt);
  k[r123f] = 0x1.74876e8000001p+35;
  k[r123b] = 0x1.9bc079007b2e4p+49 *
             exp(-0x1.0281370aec5p+0 * lgt - 0x1.387e4f585dfaep+28 / rt);
  k[r124f] = 0x1.74876e8000001p+35;
  k[r124b] = 0x1.544dc6dcfd6c4p+52 *
             exp(-0x1.1643594c6ea8p+0 * lgt - 0x1.909ab5a26ebb7p+28 / rt);
  k[r125f] = 0x1.f3ef17e000001p+35;
  k[r125b] = 0x1.4ba43d6a39b27p+38 *
             exp(-0x1.5486ff73509p-3 * lgt - 0x1.2515fdb117d5p+28 / rt);
  k[r126f] = 0x1.9254d38000001p+36 * exp(-0x1.8d23f8p+23 / rt);
  k[r126b] = 0x1.995b8fd990ff8p+40 *
             exp(-0x1.2b7cd1b2ca944p-2 * lgt - 0x1.9db3ef683cec6p+20 / rt);
  k[r127f] = 0x1.5457af8000001p+32 * exp(0x1.81a5ap+21 / rt);
  k[r127b] = 0x1.784dd71fc2457p+49 *
             exp(-0x1.e5378850e424p-1 * lgt - 0x1.d9e32499a38bcp+27 / rt);
  k[r128f] = 0x1.2a05f20000001p+35;
  k[r128b] = 0x1.df7eb7b403582p+59 *
             exp(-0x1.5f3509b5b4dp+0 * lgt - 0x1.0586c56e3f014p+29 / rt);
  k[r129f] = 0x1.bf08eb0000001p+34;
  k[r129b] = 0x1.8ae226eb6e882p+50 *
             exp(-0x1.07e8d52c28dp+0 * lgt - 0x1.b8ec06b72904dp+27 / rt);
  k[r130f] = 0x1.bf08eb0000001p+35;
  k[r130b] = 0x1.0d81a4aa08107p+52 *
             exp(-0x1.01551946ddaap+0 * lgt - 0x1.e76489b2d28adp+27 / rt);
  kTroe0 = 0x1.6c901f7d6c461p+74 *
           exp(-0x1.deb851eb851ecp+1 * lgt - 0x1.ee7248p+22 / rt);
  kTroeInf = 0x1.74876e8000001p+35;
  fcTroe = 0x1.b27bb2fec56d6p-2 * exp(-temp / 0x1.dap+7) +
           0x1.26c226809d495p-1 * exp(-temp / 0x1.9dp+10) +
           0x1p+0 * exp(-0x1.3cdp+12 / temp);
  k[r131f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.7ec598cf37b8p+106 *
           exp(-0x1.4b94f468efb9p+2 * lgt - 0x1.31eaa934fefabp+28 / rt);
  kTroeInf = 0x1.8722c3613a564p+67 *
             exp(-0x1.70e32dccb4a68p+0 * lgt - 0x1.2a30e014fefaap+28 / rt);
  k[r131b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r132f] = 0x1.61e70f6000001p+37 * exp(-0x1.f826b20000001p+25 / rt);
  k[r132b] = 0x1.5f81794b3ee43p+26 *
             exp(0x1.59d435556cc4p-1 * lgt - 0x1.40ee92efd33c8p+28 / rt);
  k[r133f] = 0x1.6069972000001p+36 * exp(0x1.070ed00000001p+21 / rt);
  k[r133b] = 0x1.36d7879f7eeafp+57 *
             exp(-0x1.5b20e451b89ep+0 * lgt - 0x1.306a94dc44adcp+28 / rt);
  k[r134f] = 0x1.74876e8000001p+35;
  k[r134b] = 0x1.0926f7f237c48p+40 *
             exp(0x1.03bf6c01328p-3 * lgt - 0x1.398c3c79cc49ep+29 / rt);
  k[r135] = 0x1.2a05f20000001p+32 * exp(-0x1.7f1808p+22 / rt);
  k[r136f] = 0x1.f400000000002p+8 * exp(0x1p+1 * lgt - 0x1.cda09p+24 / rt);
  k[r136b] = 0x1.afafd0888eaf6p+17 *
             exp(0x1.7b0b565d142e8p+0 * lgt - 0x1.d353df04cda0fp+25 / rt);
  k[r137f] = 0x1.74876e8000001p+40 * exp(-0x1.7d4e42p+25 / rt);
  k[r137b] = 0x1.268a5f406e7acp+61 *
             exp(-0x1.1455d549018cp+0 * lgt - 0x1.22c15f768ae51p+29 / rt);
  k[r138f] = 0x1.2a05f20000001p+35;
  k[r138b] = 0x1.2c6dca0e79d3bp+58 *
             exp(-0x1.396b7dd4ed48p+0 * lgt - 0x1.072282bce9ad5p+28 / rt);
  k[r139f] = 0x1.3380000000001p+11 * exp(0x1p+1 * lgt - 0x1.0803e4p+25 / rt);
  k[r139b] = 0x1.38191778cdbc4p+9 *
             exp(0x1.fe00ef57d69e8p+0 * lgt - 0x1.af60ea6d88842p+25 / rt);
  kTroe0 = 0x1.1623b50660ab4p+91 *
           exp(-0x1.470a3d70a3d71p+2 * lgt - 0x1.c501ec0000001p+24 / rt);
  kTroeInf = 0x1.823cf40000001p+29 * exp(0x1p-1 * lgt - 0x1.1ff594p+24 / rt);
  fcTroe = 0x1.a31f8a0902dep-2 * exp(-temp / 0x1.13p+8) +
           0x1.2e703afb7e91p-1 * exp(-temp / 0x1.328p+10) +
           0x1p+0 * exp(-0x1.441p+12 / temp);
  k[r140f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.6daa535e2a3bbp+129 *
           exp(-0x1.b6c342979906p+2 * lgt - 0x1.5bc2b536b7adbp+28 / rt);
  kTroeInf = 0x1.fbc7e1629aa7bp+67 *
             exp(-0x1.3ee4149bd4bbcp+0 * lgt - 0x1.5171efb6b7adap+28 / rt);
  k[r140b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r141f] = 0x1.bf08eb0000001p+34;
  k[r141b] = 0x1.28782531181ccp+35 *
             exp(0x1.ed75669bc74p-3 * lgt - 0x1.6cdcd023e0899p+28 / rt);
  k[r142f] = 0x1.bf08eb0000001p+33 * exp(-0x1.32798p+21 / rt);
  k[r142b] = 0x1.45ee2d53eef6dp+33 *
             exp(-0x1.fe3fe273a4d8p-5 * lgt - 0x1.2fc0fdfce55acp+25 / rt);
  k[r144f] = 0x1.a13b860000001p+34;
  k[r144b] = 0x1.6ce34d0565d4p+19 *
             exp(0x1.0bd5ce6826p-2 * lgt - 0x1.0d34e22f0ed7cp+28 / rt);
  k[r145f] = 0x1.65a0bc0000001p+33;
  k[r145b] = 0x1.8bd77731a313cp+29 *
             exp(0x1.038fecf5a8f8p-1 * lgt - 0x1.74096f197785cp+29 / rt);
  k[r146f] = 0x1.04c533c000001p+36;
  k[r146b] = 0x1.484cb37d4eec7p+44 *
             exp(-0x1.29cd516d12c1p-1 * lgt - 0x1.048e7e80d977dp+26 / rt);
  kTroe0 = 0x1.289c98651e77dp+107 *
           exp(-0x1.970a3d70a3d71p+2 * lgt - 0x1.41cc78p+24 / rt);
  kTroeInf = 0x1.b6605ec820001p+48 *
             exp(-0x1.28f5c28f5c28fp+0 * lgt - 0x1.246d9p+22 / rt);
  fcTroe = 0x1.96d5cfaacd9e8p-2 * exp(-temp / 0x1.ap+7) +
           0x1.3495182a9930cp-1 * exp(-temp / 0x1.ea4p+11) +
           0x1p+0 * exp(-0x1.3e2p+13 / temp);
  k[r147f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.2b3efee68d5ep+140 *
           exp(-0x1.eadf7a8b70768p+2 * lgt - 0x1.8d675cfc156f2p+28 / rt);
  kTroeInf = 0x1.ba451935a52dfp+81 *
             exp(-0x1.3c255b7d47536p+1 * lgt - 0x1.7ddc4bbc156f1p+28 / rt);
  k[r147b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r148f] = 0x1.bf08eb0000001p+34;
  k[r148b] = 0x1.45ee2d53ee30ep+34 *
             exp(-0x1.fe3fe273a4bp-5 * lgt - 0x1.1c9965fce5584p+25 / rt);
  k[r149f] = 0x1.65a0bc0000001p+33 * exp(0x1.2326bp+21 / rt);
  k[r149b] = 0x1.06d959f2c3caep+56 *
             exp(-0x1.495d7ce88a6ep+0 * lgt - 0x1.286f621c8657ap+28 / rt);
  k[r150f] = 0x1.dcd6500000001p+33 * exp(0x1.2326bp+21 / rt);
  k[r150b] = 0x1.60db8a58dd84dp+31 *
             exp(-0x1.1f10fbbc6e84p-4 * lgt - 0x1.b1c4016a6dcaep+25 / rt);
  k[r151f] = 0x1.0c388d0000001p+33;
  k[r151b] = 0x1.871dcffe513e9p+32 *
             exp(-0x1.fe3fe273a548p-5 * lgt - 0x1.1c9965fce5592p+25 / rt);
  k[r152f] = 0x1.a13b860000001p+32;
  k[r152b] = 0x1.3033a1c5cd3e9p+32 *
             exp(-0x1.fe3fe273a4b8p-5 * lgt - 0x1.1c9965fce557dp+25 / rt);
  k[r153f] = 0x1.a13b860000001p+33;
  k[r153b] = 0x1.5c45316b11561p+27 *
             exp(0x1.d31fb56b0a88p-2 * lgt - 0x1.e7799a57eca06p+27 / rt);
  k[r154f] = 0x1.2a05f20000001p+35 * exp(0x1.18ef5ffffffffp+21 / rt);
  k[r154b] = 0x1.211e17e49d6d9p+31 *
             exp(-0x1.700ff6cace4p-6 * lgt - 0x1.275402046929bp+26 / rt);
  k[r155f] = 0x1.093d9c8000001p+35 * exp(-0x1.e687c2p+26 / rt);
  k[r155b] = 0x1.5acf75ede9129p+42 *
             exp(-0x1.e7dc87d09394p-2 * lgt - 0x1.87e06afc15ef4p+23 / rt);
  k[r156f] = 0x1.135f9b0000001p+31 * exp(-0x1.4445de0000001p+26 / rt);
  k[r156b] = 0x1.2218cd62a725ep+29 *
             exp(0x1.49cdbcce251p-3 * lgt - 0x1.1f513238b3b82p+28 / rt);
  k[r157f] = 0x1.8800000000001p+4 *
             exp(0x1.3c28f5c28f5c3p+1 * lgt - 0x1.4abccp+24 / rt);
  k[r157b] = 0x1.217d4e28d8218p+7 *
             exp(0x1.485f0c6fb95acp+1 * lgt - 0x1.743b122a9db6cp+26 / rt);
  kTroe0 = 0x1.05ed2d4c301e4p+118 *
           exp(-0x1.c1eb851eb851fp+2 * lgt - 0x1.60b368p+23 / rt);
  kTroeInf = 0x1.ec95139c40001p+45 *
             exp(-0x1.2e147ae147ae1p+0 * lgt - 0x1.4e0ec00000001p+21 / rt);
  fcTroe = 0x1.8624dd2f1a9fcp-2 * exp(-temp / 0x1.24ccccccccccdp+6) +
           0x1.3ced916872b02p-1 * exp(-temp / 0x1.27p+10) +
           0x1p+0 * exp(-0x1.3878p+13 / temp);
  k[r158f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.6a1457effb6bep+155 *
           exp(-0x1.0d4b37072598p+3 * lgt - 0x1.7a0815808fc55p+28 / rt);
  kTroeInf = 0x1.5477369d882cbp+83 *
             exp(-0x1.48600f4fc9932p+1 * lgt - 0x1.719e97c08fc55p+28 / rt);
  k[r158b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r159f] = 0x1.97b21e0000001p+32 *
             exp(0x1.999999999999ap-4 * lgt - 0x1.5266460000001p+25 / rt);
  k[r159b] = 0x1.dba8052b207ep+49 *
             exp(-0x1.13949f1f2ca82p+0 * lgt - 0x1.cad22930997p+22 / rt);
  k[r160f] = 0x1.8a952f0000001p+34;
  k[r160b] = 0x1.d9feb89ef8536p+40 *
             exp(0x1.24666e8dd92p-4 * lgt - 0x1.66f20bd1ad5fp+28 / rt);
  k[r161f] = 0x1.a8f5c28f5c29p+1 *
             exp(0x1.67ae147ae147bp+1 * lgt - 0x1.7627980000001p+24 / rt);
  k[r161b] = 0x1.343b098dc4415p+3 *
             exp(0x1.5fbaf9520475cp+1 * lgt - 0x1.724b4e9be7f88p+26 / rt);
  k[r162f] = 0x1.d4c0000000001p+14 * exp(0x1.8p+0 * lgt - 0x1.3d5454p+25 / rt);
  k[r162b] = 0x1.7dc68f9f3cec7p+18 *
             exp(0x1.330f8d5c934d8p+0 * lgt - 0x1.29af456e4e42ep+26 / rt);
  k[r163f] = 0x1.3880000000001p+13 * exp(0x1.8p+0 * lgt - 0x1.3d5454p+25 / rt);
  k[r163b] = 0x1.1965573618f0ep+20 *
             exp(0x1.422f0d0a9b118p+0 * lgt - 0x1.79b8dbdb36abfp+25 / rt);
  k[r164f] = 0x1.c600000000001p+7 * exp(0x1p+1 * lgt - 0x1.25b47ap+25 / rt);
  k[r164b] = 0x1.4c9a436aab0e3p+7 *
             exp(0x1.f96c441ab54e2p+0 * lgt - 0x1.af49b84567d6fp+23 / rt);
  k[r165f] = 0x1.7fc0000000002p+12 *
             exp(0x1.bd70a3d70a3d7p+0 * lgt - 0x1.4d9c4cp+25 / rt);
  k[r165b] = 0x1.f7165b214326cp+10 *
             exp(0x1.c9a173b7a6778p+0 * lgt - 0x1.e9dcd99e6497fp+25 / rt);
  k[r166f] =
      0x1.550f7dca70001p+50 * exp(-0x1p+0 * lgt - 0x1.0f5bad0000001p+26 / rt);
  k[r166b] = 0x1.07c1e691a62dep+38 *
             exp(-0x1.88d9fbe470d01p-1 * lgt - 0x1.473ac084dc0eep+22 / rt);
  k[r167f] =
      0x1.5426a92560001p+47 * exp(-0x1p+0 * lgt - 0x1.0f5bad0000001p+26 / rt);
  k[r167b] = 0x1.070dd79e3b38fp+35 *
             exp(-0x1.88d9fbe470ceep-1 * lgt - 0x1.473ac084dc0a6p+22 / rt);
  k[r168f] = 0x1.90d75b4000001p+33 * exp(-0x1.98a2p+20 / rt);
  k[r168b] = 0x1.b4135d88eab69p+31 *
             exp(0x1.c45946566a6p-3 * lgt - 0x1.0c1cb86df2189p+27 / rt);
  k[r169f] = 0x1.0c388d0000001p+34 * exp(-0x1.cbb64p+21 / rt);
  k[r169b] = 0x1.357c89a3ee65cp+31 *
             exp(0x1.7a6598061ea2p-2 * lgt - 0x1.46ce4b34adefp+26 / rt);
  k[r170f] = 0x1.ed734dacb8924p-52 *
             exp(0x1.e666666666666p+2 * lgt + 0x1.c2c5a8p+23 / rt);
  k[r170b] = 0x1.017da88893323p-57 *
             exp(0x1.fa44dffb464f4p+2 * lgt - 0x1.6ceabbb560d91p+26 / rt);
  k[r171f] = 0x1.2a05f20000001p+33 * exp(0x1.81a5ap+21 / rt);
  k[r171b] = 0x1.7a569ea6afa9ap+24 *
             exp(0x1.b157a3c09d8p-1 * lgt - 0x1.2c4e0d1130f69p+29 / rt);
  k[r172f] = 0x1.b159800000001p+25 *
             exp(0x1.ccccccccccccdp-1 * lgt - 0x1.fd01p+22 / rt);
  k[r172b] = 0x1.353878fddf5bep+37 *
             exp(0x1.36be6a07ep-2 * lgt - 0x1.f68374aaa9cf4p+26 / rt);
  k[r173f] = 0x1.4d3d25d880001p+45 *
             exp(-0x1.63d70a3d70a3dp+0 * lgt - 0x1.033a08p+22 / rt);
  k[r173b] = 0x1.03af37dd8ae62p+39 *
             exp(-0x1.9a2eb3d7f82p-1 * lgt - 0x1.602dc72531379p+28 / rt);
  kTroe0 = 0x1.14c1a67dbf7f6p+160 *
           exp(-0x1.299999999999ap+3 * lgt - 0x1.8646e38p+28 / rt);
  kTroeInf =
      0x1.d1a94a2p+42 * exp(0x1.c28f5c28f5c29p-2 * lgt - 0x1.5a42d28p+28 / rt);
  fcTroe = 0x1.0fdf3b645a1cbp-2 * exp(-temp / 0x1.68p+7) +
           0x1.7810624dd2f1bp-1 * exp(-temp / 0x1.02cp+10) +
           0x1p+0 * exp(-0x1.529p+12 / temp);
  k[r174f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.135e19dbf90bbp+141 *
           exp(-0x1.1f8d6f0e79324p+3 * lgt - 0x1.b4cb5525c2544p+27 / rt);
  kTroeInf = 0x1.cf530dce4978cp+23 *
             exp(0x1.820a56c681574p-1 * lgt - 0x1.5cc33325c2544p+27 / rt);
  k[r174b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r175f] = 0x1.908b100000001p+29 * exp(-0x1.eed49p+23 / rt);
  k[r175b] = 0x1.34cf6dd23835ep+29 *
             exp(0x1.c068d6df3048p-4 * lgt - 0x1.094ad3bf07a6p+26 / rt);
  k[r176f] = 0x1.7d78400000001p+31 * exp(-0x1.b4374p+21 / rt);
  k[r176b] = 0x1.94d864f9060c9p-4 *
             exp(0x1.d29469e381f2p+0 * lgt - 0x1.551512068cbcdp+28 / rt);
  k[r177f] = 0x1.2a05f20000001p+33;
  k[r177b] = 0x1.940fd33626c5bp+5 *
             exp(0x1.915b1b4cdb2cp+0 * lgt - 0x1.48e798de9998ap+28 / rt);
  k[r178f] = 0x1.9254d38000001p+34 * exp(-0x1.6aa94p+20 / rt);
  k[r178b] = 0x1.3ff2e44807314p+35 *
             exp(0x1.0099698be71p-3 * lgt - 0x1.2d0f5a763e84dp+28 / rt);
  k[r179f] = 0x1.12a8800000001p+23 * exp(0x1p+0 * lgt - 0x1.9f0488p+24 / rt);
  k[r179b] = 0x1.b269ebfbde5dbp+19 *
             exp(0x1.173f2d44abb3p+0 * lgt - 0x1.2f0e5fd543d35p+27 / rt);
  k[r180f] = 0x1.f4add40000002p+34 * exp(-0x1.894f8p+20 / rt);
  k[r180b] = 0x1.8603bf58c20e7p+40 *
             exp(-0x1.624b3d0588c8p-2 * lgt - 0x1.88c2c1c3382ep+27 / rt);
  k[r181f] = 0x1.4dc9380000001p+30 * exp(-0x1.591a7cp+25 / rt);
  k[r181b] = 0x1.8555cdabc5731p+17 *
             exp(0x1.0010a9d06fb8p+0 * lgt - 0x1.63561e37db59ap+28 / rt);
  k[r182f] = 0x1.b022388000001p+34 * exp(-0x1.718689p+26 / rt);
  k[r182b] = 0x1.f541125bc2ccbp+17 *
             exp(0x1.ee7953c73e4p-1 * lgt - 0x1.cd0d4e587d7bdp+27 / rt);
  k[r183f] = 0x1.686bfd7800001p+38 * exp(-0x1.2d5de7p+26 / rt);
  k[r183b] = 0x1.aacd387ab391dp+16 *
             exp(0x1.6fe2a6567de4p+0 * lgt - 0x1.3e491e80e12b9p+28 / rt);
  k[r184f] = 0x1.dcd6500000001p+30 * exp(-0x1.502a1bp+26 / rt);
  k[r184b] = 0x1.6839d28e54c03p+21 *
             exp(0x1.58b36d2c0cd4p-1 * lgt - 0x1.73d71388c10ap+27 / rt);
  kTroe0 = 0x1.28a0514400001p+39 * exp(-0x1.c40d19p+27 / rt);
  kTroeInf = 0x1.26aba37p+36 * exp(-0x1.bf1a02p+27 / rt);
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r185f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.436f7bd609b05p+10 *
           exp(0x1.1f338c5f0f5p+0 * lgt - 0x1.fbb558c405ecap+25 / rt);
  kTroeInf = 0x1.414d8d962c822p+7 *
             exp(0x1.1f338c5f0f5p+0 * lgt - 0x1.e7e8fcc405ec8p+25 / rt);
  k[r185b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r186f] = 0x1.f7102e0000002p+30 * exp(0x1.ea5bdffffffffp+20 / rt);
  k[r186b] = 0x1.0fd8e8e91a701p+36 *
             exp(-0x1.03f2f1e8c2628p-2 * lgt - 0x1.b481fe37dd572p+24 / rt);
  k[r187f] = 0x1.81a0316280003p+46 * exp(-0x1.68f5c28f5c28fp+0 * lgt);
  k[r187b] = 0x1.20c08c671abf4p+71 *
             exp(-0x1.0e662a6949dfp+1 * lgt - 0x1.25a3602f5da3p+28 / rt);
  k[r188f] = 0x1.d0ea8e0000001p+31 * exp(0x1.ea5c7ffffffffp+19 / rt);
  k[r188b] = 0x1.4c186fdc5f638p+23 *
             exp(0x1.29675f6934eep-1 * lgt - 0x1.6854435ffa013p+27 / rt);
  k[r189f] = 0x1.ebbd028000001p+36 * exp(-0x1.6fc55ffffffffp+20 / rt);
  k[r189b] = 0x1.6499b066dd233p+19 *
             exp(0x1.0485ac3aa81cp+0 * lgt - 0x1.c537ad640b484p+26 / rt);
  k[r190f] = 0x1.2a05f20000001p+35;
  k[r190b] = 0x1.f857fd4044771p+40 *
             exp(-0x1.13904dc025fp-2 * lgt - 0x1.1c3fb353e824bp+28 / rt);
  k[r191f] = 0x1.dcd6500000001p+34 * exp(-0x1.511f2p+20 / rt);
  k[r191b] = 0x1.b3730fadf8ac2p+35 *
             exp(0x1.ca89fafc18e8p-4 * lgt - 0x1.88f9d39a739a6p+26 / rt);
  k[r192f] = 0x1.2a05f20000001p+34;
  k[r192b] = 0x1.61ea16d209ef4p+47 *
             exp(-0x1.78afbeb1e0e6p-1 * lgt - 0x1.236a56604cb7bp+26 / rt);
  k[r193f] = 0x1.e848000000001p+20 * exp(0x1.3333333333333p+0 * lgt);
  k[r193b] = 0x1.350e4ce11d1b4p+25 *
             exp(0x1.350557ede684p+0 * lgt - 0x1.3d882f9b9797fp+27 / rt);
  k[r194f] = 0x1.cd00000000001p+8 * exp(0x1p+1 * lgt - 0x1.9f0488p+24 / rt);
  k[r194b] = 0x1.15e4c331bde96p+13 *
             exp(0x1.b37a1d2d1d7e8p+0 * lgt - 0x1.d89b42119007fp+24 / rt);
  k[r195f] = 0x1.4000000000001p+10 * exp(0x1.8p+0 * lgt - 0x1.98a2p+18 / rt);
  k[r195b] = 0x1.12e35b86db324p+7 *
             exp(0x1.aaede9160424p+0 * lgt - 0x1.aec963b9dbeb8p+27 / rt);
  k[r196f] = 0x1.bf08eb0000001p+33;
  k[r196b] = 0x1.2ccdfa46dc9dbp+40 *
             exp(-0x1.268731f4662p-3 * lgt - 0x1.23f232451355p+29 / rt);
  k[r197f] = 0x1.2a05f20000001p+34 * exp(-0x1.ba2768p+25 / rt);
  k[r197b] = 0x1.fe9f1126bbe44p+43 *
             exp(-0x1.4302c8c7c30cp-1 * lgt - 0x1.09230243912c6p+26 / rt);
  k[r198f] = 0x1.41dd760000001p+34 * exp(-0x1.d70a3d70a3d71p-3 * lgt);
  k[r198b] = 0x1.16084d3436358p+35 *
             exp(-0x1.c7d7aacfbc8p-6 * lgt - 0x1.850c53288a92bp+28 / rt);
  k[r199f] = 0x1.53eec80800001p+38 * exp(-0x1.ccccccccccccdp-2 * lgt);
  k[r199b] = 0x1.eff0c5c699065p+60 *
             exp(-0x1.af53f086daap+0 * lgt - 0x1.24355ccf52d52p+27 / rt);
  k[r200f] = 0x1.65a0bc0000001p+31;
  k[r200b] = 0x1.ef88fd7b8392ep+26 *
             exp(0x1.0012cca58bf4p-2 * lgt - 0x1.47a313a474aeep+25 / rt);
  k[r201f] = 0x1.229298c000001p+35;
  k[r201b] = 0x1.de21442de9f03p+43 *
             exp(-0x1.f14cb0be35d8p-2 * lgt - 0x1.c73be032870f1p+26 / rt);
  k[r202f] = 0x1.2a05f20000001p+35 * exp(-0x1.d21917fffffffp+23 / rt);
  k[r202b] = 0x1.5b2a3a07393c5p+31 *
             exp(0x1.23fa5c1f2fc8p-2 * lgt - 0x1.f856f846fb7dfp+25 / rt);
  k[r203f] = 0x1.5f90000000001p+16 * exp(0x1.8p+0 * lgt + 0x1.d5edep+20 / rt);
  k[r203b] = 0x1.1bd9610615ec6p+16 *
             exp(0x1.ae281c12bd818p+0 * lgt - 0x1.b1eba9c039519p+26 / rt);
  k[r204f] = 0x1.3ab668p+28;
  k[r204b] = 0x1.777d142080db7p+16 *
             exp(0x1.22970786df73p-2 * lgt - 0x1.db97ded5cf17ep+24 / rt);
  k[r205f] = 0x1.e449a94000001p+36 *
             exp(-0x1.c28f5c28f5c29p-4 * lgt - 0x1.3df79cp+24 / rt);
  k[r205b] = 0x1.20e7df9b71f69p+25 *
             exp(0x1.63e660f94628p-3 * lgt - 0x1.8cc7bd6ae7929p+25 / rt);
  k[r206f] = 0x1.2a05f20000001p+32;
  k[r206b] = 0x1.f434c5241f48ep+30 *
             exp(0x1.1677a27af9dcp-2 * lgt - 0x1.c1d270c485209p+27 / rt);
  k[r207f] = 0x1.74876e8000001p+34;
  k[r207b] = 0x1.e2b3e42914d16p+29 *
             exp(0x1.32a9b7b250cp-1 * lgt - 0x1.b33b0415bd666p+28 / rt);
  k[r208f] = 0x1.04c533c000001p+36;
  k[r208b] = 0x1.87295c9a4ece9p+30 *
             exp(0x1.964bcb594b1cp-2 * lgt - 0x1.7175876996a21p+25 / rt);
  k[r209f] = 0x1.74876e8000001p+35;
  k[r209b] = 0x1.95cf476b5a872p+31 *
             exp(0x1.449d7f6f226p-1 * lgt - 0x1.bac0b7ea0e3fp+28 / rt);
  k[r210f] = 0x1.2a05f20000001p+34;
  k[r210b] = 0x1.c207bb5a4aa98p+33 *
             exp(0x1.0ef0898505b8p-1 * lgt - 0x1.f89779f13d24cp+28 / rt);
  k[r211f] = 0x1.74876e8000001p+34;
  k[r211b] = 0x1.59341c3fdb874p+41 *
             exp(0x1.f5926bcceaap-4 * lgt - 0x1.c36589fcf6e4p+28 / rt);
  k[r212f] = 0x1.45f680b000002p+45 *
             exp(-0x1.51eb851eb851fp+0 * lgt - 0x1.79fc6p+21 / rt);
  k[r212b] = 0x1.f0d4278a6690bp+59 *
             exp(-0x1.78b037102ebep+0 * lgt - 0x1.8a20c25916fcfp+27 / rt);
  k[r213f] = 0x1.74876e8000001p+34;
  k[r213b] = 0x1.096f681be100dp+27 *
             exp(0x1.ddcf2fa39bfp-2 * lgt - 0x1.a6ca3b77a9edp+27 / rt);
  k[r214f] = 0x1.ad27480000001p+29 *
             exp(0x1.70a3d70a3d70ap-1 * lgt - 0x1.511f7p+21 / rt);
  k[r214b] = 0x1.011242a220709p+23 *
             exp(0x1.38bf9b4c6e7bp+0 * lgt - 0x1.bb1a20e04b9fdp+27 / rt);
  k[r215f] = 0x1.9640000000001p+13 *
             exp(0x1.e666666666666p+0 * lgt + 0x1.e5406p+21 / rt);
  k[r215b] = 0x1.515636b0cba85p+10 *
             exp(0x1.25fecd9c53eap+1 * lgt - 0x1.14f712d754b54p+28 / rt);
  k[r216f] = 0x1.2a05f20000001p+33 * exp(-0x1.9f0488p+25 / rt);
  k[r216b] = 0x1.130f69ee6d7afp+29 *
             exp(0x1.1de6c573ecc4p-3 * lgt - 0x1.a79f1842d11e9p+25 / rt);
  k[r217f] = 0x1.1ed8ec2000001p+36;
  k[r217b] = 0x1.f79dc7e3bfa25p+37 *
             exp(0x1.7c202dbf558p-5 * lgt - 0x1.369d2463311cp+28 / rt);
  k[r218f] = 0x1.2a05f20000001p+35;
  k[r218b] = 0x1.715da0821e948p+56 *
             exp(-0x1.6d8968d0167ap+0 * lgt - 0x1.fcbd82181f5cp+26 / rt);
  k[r219f] = 0x1.dcd6500000001p+32 * exp(-0x1.dc5053fffffffp+24 / rt);
  k[r219b] = 0x1.193e5dc2afca8p+41 *
             exp(-0x1.787d91a230e1p-1 * lgt - 0x1.c93d7fa211ab9p+25 / rt);
  k[r220f] = 0x1.6df8f70000001p+32 * exp(0x1.c17f4p+20 / rt);
  k[r220b] = 0x1.cc7be0baf117ep+44 *
             exp(-0x1.fb6ed89411bfp-1 * lgt - 0x1.c165ba786d412p+25 / rt);
  k[r221f] = 0x1.2700000000001p+8 *
             exp(0x1.399999999999ap+1 * lgt - 0x1.1e0b080000001p+23 / rt);
  k[r221b] = 0x1.e262947d69c66p+19 *
             exp(0x1.9c1def6d0c36p+0 * lgt - 0x1.88a713edc4695p+26 / rt);
  k[r222f] = 0x1.5e2d62c000001p+34;
  k[r222b] = 0x1.8269c7dae49b2p+20 *
             exp(0x1.20d79afcaed6p+0 * lgt - 0x1.7a45d53ec55cfp+28 / rt);
  k[r223f] = 0x1.9254d38000001p+35;
  k[r223b] = 0x1.06581961a013p+16 *
             exp(0x1.65bbae6cb837p+0 * lgt - 0x1.781887ab74e44p+26 / rt);
  k[r224f] = 0x1.2a05f20000001p+31;
  k[r224b] = 0x1.2ed09b8d0ed54p+10 *
             exp(0x1.a051020a7c96cp-1 * lgt + 0x1.b3bb0e99b179dp+24 / rt);
  k[r225f] = 0x1.2a05f20000001p+34;
  k[r225b] = 0x1.0586032dcc642p+21 *
             exp(0x1.40eac82e2b68p+0 * lgt - 0x1.52f5433a81f1bp+29 / rt);
  k[r226f] = 0x1.dcd6500000001p+30 * exp(-0x1.3f3ea9p+26 / rt);
  k[r226b] = 0x1.5f6eab3e19687p+30 *
             exp(0x1.2572818e3fep-2 * lgt - 0x1.ed41c08009fd7p+28 / rt);
  k[r227f] = 0x1.20b5c27000001p+38 * exp(-0x1.af611bp+27 / rt);
  k[r227b] = 0x1.7896bba80b6ddp+11 *
             exp(0x1.28bb5046a3508p+0 * lgt + 0x1.98cd32cdc94fep+22 / rt);
  k[r228f] = 0x1.599ba503c0001p+47 *
             exp(-0x1.851eb851eb852p+0 * lgt - 0x1.79fc6p+21 / rt);
  k[r228b] = 0x1.48c82839c1c5ep+50 *
             exp(-0x1.5b83c738da8cp+0 * lgt - 0x1.e62992250d458p+27 / rt);
  k[r229f] = 0x1.b0028e44b0001p+51 * exp(-0x1p+1 * lgt - 0x1.98a2p+21 / rt);
  k[r229b] = 0x1.402b9547b2283p+55 *
             exp(-0x1.adcf5faf9f64p+0 * lgt - 0x1.27589205d34aep+29 / rt);
  k[r230f] = 0x1.581b6d300d024p+86 *
             exp(-0x1.a666666666666p+1 * lgt - 0x1.f93485p+28 / rt);
  k[r230b] = 0x1.ccf8359d1f3e9p+66 *
             exp(-0x1.67b0436e3c472p+1 * lgt - 0x1.79ef23eed7c1ep+21 / rt);
  k[r231f] = 0x1.44ccccccccccdp+4 *
             exp(0x1.51eb851eb851fp+1 * lgt - 0x1.3df79cp+24 / rt);
  k[r231b] = 0x1.24d32bcc3e189p+13 *
             exp(0x1.023480aa8c4f2p+1 * lgt - 0x1.927dcdb22f2e5p+25 / rt);
  k[r232f] = 0x1.447ae147ae148p+2 *
             exp(0x1.51eb851eb851fp+1 * lgt - 0x1.3df79cp+24 / rt);
  k[r232b] = 0x1.7d80a3525c535p-9 *
             exp(0x1.b51257e0e8788p+1 * lgt - 0x1.20abb742463c1p+27 / rt);
  k[r233f] = 0x1.dd4b800000001p+21 *
             exp(0x1.947ae147ae148p+0 * lgt - 0x1.a8984ap+26 / rt);
  k[r233b] = 0x1.5b31adeec3a97p+9 *
             exp(0x1.314b2097b6cf6p+1 * lgt - 0x1.2cde3e07c25cep+23 / rt);
  k[r234f] = 0x1.13p+10 * exp(0x1.03d70a3d70a3dp+1 * lgt - 0x1.aad48ep+25 / rt);
  k[r234b] = 0x1.1928cce8e02cbp+22 *
             exp(0x1.56bc741b1222ap+0 * lgt - 0x1.570ebcd1442ccp+24 / rt);
  k[r235] = 0x1.199999999999bp+2 *
            exp(0x1.2147ae147ae14p+1 * lgt - 0x1.98a228p+24 / rt);
  k[r236f] = 0x1.47ae147ae147cp-3 *
             exp(0x1.47ae147ae147bp+1 * lgt - 0x1.1f51e8p+25 / rt);
  k[r236b] = 0x1.160bb7cd41718p-9 *
             exp(0x1.8ad28da86006p+1 * lgt - 0x1.ddb0f1b252244p+26 / rt);
  kTroe0 = 0x1.e5b8fa8fe2ac3p+66 *
           exp(-0x1.b333333333333p+1 * lgt - 0x1.e540880000001p+22 / rt);
  kTroeInf = 0x1.ebbd028000001p+34;
  fcTroe = 0x1p+0 * exp(-0x0p+0 / temp);
  k[r237f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.2775a84ef92b6p+79 *
           exp(-0x1.c2ac1cbfa23bp+1 * lgt - 0x1.a3eee7d5acb0bp+26 / rt);
  kTroeInf = 0x1.2b1e70e97ff24p+47 *
             exp(-0x1.ef1d318de0fap-4 * lgt - 0x1.859adf55acb0bp+26 / rt);
  k[r237b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r238f] = 0x1.bf08eb0000001p+35 * exp(-0x1.98a2p+20 / rt);
  k[r238b] = 0x1.ec1763f6a50a8p+24 *
             exp(0x1.bab2386e30f8p-1 * lgt - 0x1.35a9d12f723f4p+28 / rt);
  k[r239f] = 0x1.d562f6c000001p+35 * exp(-0x1.6f4a5ep+27 / rt);
  k[r239b] = 0x1.866e6e34d25cbp+35 *
             exp(-0x1.48b35289072cp-4 * lgt + 0x1.b2b8edaa8686p+20 / rt);
  k[r240f] = 0x1.7cdc000000001p+21 *
             exp(0x1.c28f5c28f5c29p-1 * lgt - 0x1.4151dp+26 / rt);
  k[r240b] = 0x1.08310f6ff201dp+32 *
             exp(0x1.8d5172fb405cp-3 * lgt - 0x1.34616c9cf08cbp+26 / rt);
  kTroe0 = 0x1.68d28e3f00282p+63 *
           exp(-0x1.947ae147ae148p+1 * lgt - 0x1.79fc6p+21 / rt);
  kTroeInf = 0x1.718c7e0000001p+31 * exp(0x1.3333333333333p-3 * lgt);
  fcTroe = 0x1.54fdf3b645a1dp-2 * exp(-temp / 0x1.d6p+7) +
           0x1.55810624dd2f2p-1 * exp(-temp / 0x1.08ap+11) +
           0x1p+0 * exp(-0x1.1b8p+12 / temp);
  k[r241f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.0b1dfac96ae22p+96 *
           exp(-0x1.34211a05cce4ap+2 * lgt - 0x1.0f79b2e0d58dcp+27 / rt);
  kTroeInf = 0x1.1193c3a54c167p+64 *
             exp(-0x1.81283f2171032p+0 * lgt - 0x1.0991c160d58dbp+27 / rt);
  k[r241b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r242f] = 0x1.2a05f20000001p+33 * exp(-0x1.274d2bp+28 / rt);
  k[r242b] = 0x1.bcfccd14dc0dp+38 *
             exp(-0x1.02cdd5f03d92p-1 * lgt - 0x1.9bde50626e13bp+27 / rt);
  k[r243f] = 0x1.7d78400000001p+26 * exp(-0x1.0362d5p+28 / rt);
  k[r243b] = 0x1.9f47bd5cfcd82p+31 *
             exp(-0x1.22b1d4177896p-1 * lgt - 0x1.9b2ffde1a7656p+27 / rt);
  k[r244f] = 0x1.1b1f3f8000001p+34;
  k[r244b] = 0x1.768dcddd292dfp+34 *
             exp(0x1.70ff011d3e2p-5 * lgt - 0x1.c93325224fd9dp+26 / rt);
  k[r245f] = 0x1.b022388000001p+34;
  k[r245b] = 0x1.f5dab7eb851dp+36 *
             exp(0x1.768f976e496p-4 * lgt - 0x1.a8e9edabc512ap+28 / rt);
  k[r246f] = 0x1.31794b4000001p+35;
  k[r246b] = 0x1.5104fc385036dp+46 *
             exp(-0x1.1f14a5072ep-1 * lgt - 0x1.2868985d7aa42p+28 / rt);
  k[r247f] = 0x1.e2cc310000001p+33;
  k[r247b] = 0x1.e037b73db4855p+53 *
             exp(-0x1.2ef85b6bef46p+0 * lgt - 0x1.46d8d853c088bp+28 / rt);
  k[r248f] = 0x1.6e918d8000001p+34;
  k[r248b] = 0x1.337d5b3736639p+40 *
             exp(-0x1.074034cc5634p-2 * lgt - 0x1.4efe2c8cebdafp+27 / rt);
  k[r249f] = 0x1.19f17fe160001p+48 *
             exp(-0x1.6147ae147ae14p+0 * lgt - 0x1.445a38p+22 / rt);
  k[r249b] = 0x1.3bb3d8afe691ep+70 *
             exp(-0x1.5b6321459a34p+1 * lgt - 0x1.778d31c07c861p+28 / rt);
  k[r250f] = 0x1.0e15635000001p+38 *
             exp(-0x1.6147ae147ae14p-1 * lgt - 0x1.8433affffffffp+21 / rt);
  k[r250b] = 0x1.5c59ee2cd3705p+44 *
             exp(-0x1.fc91b1ff14f4p-1 * lgt - 0x1.2eb6b7b9c18f2p+28 / rt);
  k[r251f] = 0x1.1b1f3f8000001p+35 *
             exp(-0x1.70a3d70a3d70ap-2 * lgt - 0x1.2842300000001p+21 / rt);
  k[r251b] = 0x1.30d23af0b7639p+59 *
             exp(-0x1.ed580e457bf4p+0 * lgt - 0x1.86a84b7deafd5p+26 / rt);
  k[r252f] = 0x1.19f17fe160001p+48 *
             exp(-0x1.6147ae147ae14p+0 * lgt - 0x1.445a38p+22 / rt);
  k[r252b] = 0x1.cc5a84462d5a7p+69 *
             exp(-0x1.635c20cf68c7p+1 * lgt - 0x1.9b205e8019309p+28 / rt);
  k[r253f] = 0x1.0e15635000001p+38 *
             exp(-0x1.6147ae147ae14p-1 * lgt - 0x1.8433affffffffp+21 / rt);
  k[r253b] = 0x1.fbf62c4980fc6p+43 *
             exp(-0x1.0e3ad81327a4p+0 * lgt - 0x1.5249e4795e3a9p+28 / rt);
  k[r254f] = 0x1.1b1f3f8000001p+35 *
             exp(-0x1.70a3d70a3d70ap-2 * lgt - 0x1.2842300000001p+21 / rt);
  k[r254b] = 0x1.bc7c90a70e6d1p+58 *
             exp(-0x1.fd4a0d59190ep+0 * lgt - 0x1.0a7a7f3e2ed3cp+27 / rt);
  k[r255f] = 0x1.65a0bc0000001p+36 * exp(-0x1.cbb686p+26 / rt);
  k[r255b] = 0x1.724c3f7234c5bp+37 *
             exp(0x1.c792cb8924ap-4 * lgt - 0x1.bee2410056c73p+28 / rt);
  k[r256f] = 0x1.dcd6500000001p+29 * exp(-0x1.5b2dbcp+26 / rt);
  k[r256b] = 0x1.da97afd9ed5fbp+31 *
             exp(-0x1.052ca5045d28p-2 * lgt - 0x1.4a46528f0df67p+25 / rt);
  k[r257f] = 0x1.47d3570000001p+34;
  k[r257b] = 0x1.e60a22fb1e4e1p+7 *
             exp(0x1.996de518970cp+0 * lgt - 0x1.1e2ae6710d8f4p+29 / rt);
  k[r258f] = 0x1.dcd6500000001p+30;
  k[r258b] = 0x1.79fc36e06b8fep+21 *
             exp(0x1.e01ae111a37p-1 * lgt - 0x1.d0549e00cde54p+28 / rt);
  k[r259f] = 0x1.65a0bc0000001p+33;
  k[r259b] = 0x1.4076d78b7c988p+3 *
             exp(0x1.7cfdc5996d4dp+0 * lgt - 0x1.409a3a015a1dcp+27 / rt);
  k[r260f] = 0x1.65a0bc0000001p+33;
  k[r260b] = 0x1.3ba90e9b150b4p+12 *
             exp(0x1.0d2bc9135fbp+0 * lgt - 0x1.cb1c8def4e7b2p+27 / rt);
  k[r261f] = 0x1.74876e8000001p+36;
  k[r261b] = 0x1.d37775d51cf4dp+15 *
             exp(0x1.b6586ef0fep+0 * lgt - 0x1.0d72ed7baec34p+28 / rt);
  k[r262f] = 0x1.7ed0000000001p+16 *
             exp(0x1.68f5c28f5c28fp+0 * lgt - 0x1.0f5b94p+25 / rt);
  k[r262b] = 0x1.551293475232p+12 *
             exp(0x1.efd00375fbcbp+0 * lgt - 0x1.80038da54de5ap+27 / rt);
  k[r263f] = 0x1.24f8000000001p+17 *
             exp(0x1.91eb851eb851fp+0 * lgt - 0x1.5f2b6ap+27 / rt);
  k[r263b] = 0x1.d72cbb0f65cafp+3 *
             exp(0x1.51dd6ee5cb6cp+1 * lgt - 0x1.2e1933f27b5b2p+28 / rt);
  k[r264f] = 0x1.13p+11 * exp(0x1.0e147ae147ae1p+1 * lgt - 0x1.6bf05cp+25 / rt);
  k[r264b] = 0x1.1d9479d7a93c6p+4 *
             exp(0x1.424a3fadd21bfp+1 * lgt - 0x1.4834bbe2bc03fp+22 / rt);
  k[r265f] = 0x1.5f90000000001p+14 *
             exp(0x1.b333333333333p+0 * lgt - 0x1.e54038p+23 / rt);
  k[r265b] = 0x1.579d6b774ee9dp-8 *
             exp(0x1.a0aadc07ce72ep+1 * lgt - 0x1.bde88a5d991d6p+24 / rt);
  k[r266f] = 0x1.a400000000002p+6 * exp(0x1.4p+1 * lgt - 0x1.a8987cp+25 / rt);
  k[r266b] = 0x1.6eae0a591f711p+0 *
             exp(0x1.78b2b6bbbedcfp+1 * lgt - 0x1.43b8ac3dbc9ap+24 / rt);
  k[r267f] = 0x1.01d0000000001p+15 * exp(0x1.8p+0 * lgt - 0x1.cbb6b8p+23 / rt);
  k[r267b] = 0x1.380303d8fa7c4p+12 *
             exp(0x1.d68ef2826f078p+0 * lgt - 0x1.5ae7985855689p+25 / rt);
  k[r268f] = 0x1.9c80000000002p+11 * exp(0x1.8p+0 * lgt - 0x1.cbb6b8p+23 / rt);
  k[r268b] = 0x1.093dd587ae97fp+12 *
             exp(0x1.c6d58dbd3c65p+0 * lgt - 0x1.06ff4f3c30b84p+27 / rt);
  k[r269f] = 0x1.576cd9de00001p+43 * exp(-0x1.5214958p+28 / rt);
  k[r269b] = 0x1.ac42a356243f7p+9 *
             exp(0x1.7d781e0e6069bp+0 * lgt + 0x1.b880e9bf07e7cp+23 / rt);
  k[r270f] = 0x1.e8f1c10800001p+40 *
             exp(-0x1.6147ae147ae14p-1 * lgt - 0x1.6bf0bffffffffp+23 / rt);
  k[r270b] = 0x1.fc83f45b0798dp+38 *
             exp(-0x1.d3cd4bf82928p-2 * lgt - 0x1.1e81c06101c75p+28 / rt);
  k[r271f] = 0x1.017df80000001p+28 *
             exp(0x1.70a3d70a3d70ap-3 * lgt - 0x1.0eb837fffffffp+23 / rt);
  k[r271b] = 0x1.3478175749333p+10 *
             exp(0x1.719e916ee74bp+0 * lgt - 0x1.a99507348da18p+27 / rt);
  k[r272f] = 0x1.3ca6512000001p+37 * exp(-0x1.8p-1 * lgt - 0x1.710bc8p+23 / rt);
  k[r272b] = 0x1.41e122d39daa8p+13 *
             exp(0x1.09d308e89cbp+0 * lgt - 0x1.2b5f1f86db585p+28 / rt);
  k[r273] =
      0x1.3880000000001p+14 * exp(0x1p+1 * lgt - 0x1.fecaa7fffffffp+22 / rt);
  k[r274f] = 0x1.0c388d0000001p+33;
  k[r274b] = 0x1.ff0150636be65p+36 *
             exp(-0x1.182214d116bp-4 * lgt - 0x1.8ec8f92d2a9eep+27 / rt);
  k[r275f] = 0x1.1c0daaa800001p+39 *
             exp(-0x1.3d70a3d70a3d7p-2 * lgt - 0x1.28427ffffffffp+20 / rt);
  k[r275b] = 0x1.b874ac3f4eda3p+46 *
             exp(-0x1.d27442f07a6ap-1 * lgt - 0x1.2cfb5e66fba45p+27 / rt);
  k[r276f] = 0x1.b9130a0000001p+31 *
             exp(0x1.3333333333333p-3 * lgt + 0x1.6fc5p+18 / rt);
  k[r276b] = 0x1.00a3c5ca77885p+35 *
             exp(0x1.48cfb4a8d6p-6 * lgt - 0x1.d099fd9ac3f3ap+28 / rt);
  k[r277f] = 0x1.0ep+9 * exp(0x1.3333333333333p+1 * lgt - 0x1.3c8808p+25 / rt);
  k[r277b] = 0x1.6f2ff0b4f0e88p+1 *
             exp(0x1.618189e956731p+1 * lgt - 0x1.197cd78656518p+24 / rt);
  k[r278f] = 0x1.86a0000000001p+15 *
             exp(0x1.999999999999ap+0 * lgt - 0x1.e7ce6ffffffffp+21 / rt);
  k[r278b] = 0x1.7033104285b73p+11 *
             exp(0x1.db5fcc10d171p+0 * lgt - 0x1.5d695afca2459p+25 / rt);
  k[r279f] = 0x1.25c0000000001p+13 *
             exp(0x1.f0a3d70a3d70ap+0 * lgt - 0x1.9c769ffffffffp+24 / rt);
  k[r279b] = 0x1.db2e547819d3ep+4 *
             exp(0x1.2223504c0d74ap+1 * lgt + 0x1.dbbeadf5ba25p+21 / rt);
  k[r280f] = 0x1.2a05f20000001p+33 * exp(-0x1.ca1dbcp+25 / rt);
  k[r280b] = 0x1.0cfae728e8dd6p+24 *
             exp(0x1.15ea2f8c7b3e8p-1 * lgt - 0x1.9b0e21fd465dcp+24 / rt);
  k[r281f] = 0x1.668f272800001p+42 *
             exp(-0x1.810624dd2f1aap-1 * lgt - 0x1.60724p+20 / rt);
  k[r281b] = 0x1.42443462f4603p+46 *
             exp(-0x1.2986cf0405dep+0 * lgt - 0x1.e0dbf17e15512p+27 / rt);
  k[r282f] = 0x1.836e210000001p+31 * exp(0x1.681b8p+21 / rt);
  k[r282b] = 0x1.5fad300a99dd8p+39 *
             exp(-0x1.92c59cb7434p-4 * lgt - 0x1.c59c2a43c83ddp+28 / rt);
  k[r283f] = 0x1.65a0bc0000001p+31 * exp(-0x1.68bf2cp+25 / rt);
  k[r283b] = 0x1.a773da7dad1f3p+14 *
             exp(0x1.dd744fbb975cp-1 * lgt - 0x1.0f0517d2ba9cep+27 / rt);
  k[r284] = 0x1.f62b4c4000002p+34;
  k[r285f] = 0x1.a2c0000000001p+12 *
             exp(0x1.d47ae147ae148p+0 * lgt - 0x1.c17f2p+19 / rt);
  k[r285b] = 0x1.97b037b64d252p+15 *
             exp(0x1.6795b4508cc0cp+0 * lgt - 0x1.bf74e04c30d6cp+25 / rt);
  k[r286f] = 0x1.984ab48000001p+36;
  k[r286b] = 0x1.a3c23df6b6659p+46 *
             exp(-0x1.0696b78d12f8p-1 * lgt - 0x1.2e77eafe57ba3p+28 / rt);
  k[r287f] = 0x1.2309ce5400001p+42 * exp(-0x1.14a01dp+26 / rt);
  k[r287b] = 0x1.05d81740a0223p+43 *
             exp(0x1.0769708f15cp-2 * lgt - 0x1.5cd648cefa97p+28 / rt);
  k[r288] = 0x1.e848000000001p+22 * exp(0x1p-1 * lgt + 0x1.c03838p+22 / rt);
  kTroe0 = 0x1.4e7466500eea2p+65 *
           exp(-0x1.6666666666666p+1 * lgt - 0x1.2d5dfffffffffp+21 / rt);
  kTroeInf = 0x1.d5af420000001p+30 *
             exp(0x1.b851eb851eb85p-2 * lgt + 0x1.79fc6p+20 / rt);
  fcTroe = 0x1.b020c49ba5e35p-2 * exp(-temp / 0x1.e8p+6) +
           0x1.27ef9db22d0e5p-1 * exp(-temp / 0x1.3cep+11) +
           0x1p+0 * exp(-0x1.24a8p+13 / temp);
  k[r289f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.0c3f320be8ab4p+86 *
           exp(-0x1.a17bd683492ep+1 * lgt - 0x1.b226fd0cb33c5p+28 / rt);
  kTroeInf = 0x1.78b4fd4d2021dp+51 *
             exp(-0x1.02ccab0fbc258p-5 * lgt - 0x1.ae5244acb33c5p+28 / rt);
  k[r289b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r290] = 0x1.59b4fa0000001p+32 * exp(-0x1.7f1808p+22 / rt);
  k[r291f] = 0x1.1e1a300000001p+31 * exp(-0x1.7f1808p+22 / rt);
  k[r291b] = 0x1.b58933d2ae7dep+38 *
             exp(-0x1.4b0438aaff28p-2 * lgt - 0x1.f2a4831b3c825p+27 / rt);
  k[r292] = 0x1.74876e8000001p+37 * exp(-0x1.5ed16p+25 / rt);
  k[r293] = 0x1.0429900000001p+26 * exp(0x1p-2 * lgt + 0x1.dd972p+21 / rt);
  k[r294f] = 0x1.20f69c0000001p+28 *
             exp(0x1.28f5c28f5c28fp-2 * lgt - 0x1.6798cccccccccp+15 / rt);
  k[r294b] = 0x1.16ce396c7a85cp+34 *
             exp(-0x1.31edc822596fp-3 * lgt - 0x1.93eb6ae0092afp+24 / rt);
  k[r295f] = 0x1.4e40000000001p+10 *
             exp(0x1.9c28f5c28f5c3p+0 * lgt + 0x1.8849a00000001p+20 / rt);
  k[r295b] = 0x1.0f9dc2b613482p+13 *
             exp(0x1.6d2093c31caf8p+0 * lgt - 0x1.ad3a5bc420fe4p+25 / rt);
  k[r296f] = 0x1.5c17540000001p+32 * exp(-0x1.cdc1600000001p+22 / rt);
  k[r296b] = 0x1.8866b93d7261ep+21 *
             exp(0x1.0c475eef89c1cp-1 * lgt - 0x1.4fad31e6b1ff9p+24 / rt);
  k[r297] = 0x1.5c17540000001p+32 * exp(-0x1.cdc1600000001p+22 / rt);
  k[r298] = 0x1.c086634000001p+34 * exp(-0x1.3875bfp+27 / rt);
  k[r299f] = 0x1.f47d000000002p+20 *
             exp(0x1.28f5c28f5c28fp+0 * lgt - 0x1.331d18p+23 / rt);
  k[r299b] = 0x1.da52a59fe3339p+10 *
             exp(0x1.b81355e589dfcp+0 * lgt - 0x1.ee26a32bbf916p+24 / rt);
  k[r300] = 0x1.f47d000000002p+20 *
            exp(0x1.28f5c28f5c28fp+0 * lgt - 0x1.331d18p+23 / rt);
  k[r301] = 0x1.6583700000001p+24 *
            exp(0x1.75c28f5c28f5cp-1 * lgt + 0x1.1c416p+22 / rt);
  k[r302] = 0x1.66d1e90000001p+31 * exp(-0x1.7ca294p+25 / rt);
  k[r303] = 0x1.5400000000001p+11 *
            exp(0x1.c51eb851eb852p+0 * lgt - 0x1.79fc74p+24 / rt);
  kTroe0 = 0x1.85cee78dff543p+119 *
           exp(-0x1.e851eb851eb85p+2 * lgt - 0x1.ec25d7fffffffp+23 / rt);
  kTroeInf = 0x1.cff66a0000001p+28 *
             exp(0x1.b020c49ba5e35p-2 * lgt + 0x1.c03838p+22 / rt);
  fcTroe = 0x1.11eb851eb851fp-1 * exp(-temp / 0x1.92p+7) +
           0x1.dc28f5c28f5c3p-2 * exp(-temp / 0x1.bb4p+10) +
           0x1p+0 * exp(-0x1.4d5p+12 / temp);
  k[r304f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.844c3653f6cd1p+129 *
           exp(-0x1.e5c5598e57fe8p+2 * lgt - 0x1.33a789ed2d808p+27 / rt);
  kTroeInf = 0x1.ce2a29021f5dp+38 *
             exp(0x1.d8e9e40811805p-2 * lgt - 0x1.06e36aad2d809p+27 / rt);
  k[r304b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r305] = 0x1.176592e000001p+37;
  k[r306] = 0x1.142f200000001p+24;
  k[r307] = 0x1.6694e00000001p+24;
  k[r308f] = 0x1.47d3570000001p+34;
  k[r308b] = 0x1.fe5200dd4b019p+14 *
             exp(0x1.4dd31ee5e0c48p+0 * lgt - 0x1.7b64036f493f5p+25 / rt);
  k[r309f] = 0x1.47d3570000001p+33;
  k[r309b] = 0x1.2c79a2313aad4p+31 *
             exp(0x1.3ddad7eafa18p-2 * lgt - 0x1.1294a3c61a8edp+28 / rt);
  k[r310f] = 0x1.65a0bc0000001p+33;
  k[r310b] = 0x1.c663324b52e6ap+34 *
             exp(0x1.a501d82d832p-3 * lgt - 0x1.506b65cd49738p+28 / rt);
  k[r311f] = 0x1.c086634000001p+34;
  k[r311b] = 0x1.96ad073c6a29fp+28 *
             exp(0x1.53aa02f693edp-2 * lgt - 0x1.13fb98e785d9bp+25 / rt);
  kTroe0 = 0x1.41a99040f6102p+227 *
           exp(-0x1.0d1eb851eb852p+4 * lgt - 0x1.a117c8p+25 / rt);
  kTroeInf = 0x1.190930c000001p+33;
  fcTroe = 0x1.b1d14e3bcd35bp-1 * exp(-temp / 0x1.23p+8) +
           0x1.38bac710cb296p-3 * exp(-temp / 0x1.56cp+11) +
           0x1p+0 * exp(-0x1.e44p+12 / temp);
  k[r312f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.b33488203ecfp+269 *
           exp(-0x1.294cd65c2733ep+4 * lgt - 0x1.9a58ed7546cefp+28 / rt);
  kTroeInf = 0x1.7c3cf48c8921p+75 *
             exp(-0x1.c2e1e0a3baecp+0 * lgt - 0x1.6635f47546cfp+28 / rt);
  k[r312b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r313f] = 0x1.8200000000001p+7 *
             exp(0x1.570a3d70a3d71p+1 * lgt - 0x1.da862p+23 / rt);
  k[r313b] = 0x1.0e5bf12a12932p-6 *
             exp(0x1.9be4e455dee99p+1 * lgt - 0x1.05b6b07deda7bp+24 / rt);
  k[r314f] = 0x1.4a00000000001p+10 *
             exp(0x1.451eb851eb852p+1 * lgt - 0x1.af5ce8p+24 / rt);
  k[r314b] = 0x1.84a2358f8be26p-3 *
             exp(0x1.8e7651265b104p+1 * lgt - 0x1.2015e2e17da06p+25 / rt);
  k[r315f] = 0x1.edc0000000001p+14 *
             exp(0x1.ccccccccccccdp+0 * lgt - 0x1.dd143p+21 / rt);
  k[r315b] = 0x1.93069916c9a11p+5 *
             exp(0x1.2252c1c04eb34p+1 * lgt - 0x1.2a77610d7a60bp+26 / rt);
  k[r316] = 0x1.83126e978d4ffp-2 *
            exp(0x1.5c28f5c28f5c3p+1 * lgt - 0x1.7f1808p+22 / rt);
  k[r317f] = 0x1.d96e9bbf0dc78p-11 *
             exp(0x1.d333333333333p+1 * lgt - 0x1.c8c65ffffffffp+24 / rt);
  k[r317b] = 0x1.da4847fe3bf18p-13 *
             exp(0x1.db0fff8a418d8p+1 * lgt - 0x1.71f12f78c2c0ap+25 / rt);
  kTroe0 = 0x1.e965d885673ap+190 *
           exp(-0x1.d333333333333p+3 * lgt - 0x1.2208a2p+26 / rt);
  kTroeInf = 0x1.3ec0000000001p+11 *
             exp(0x1.999999999999ap+0 * lgt - 0x1.6bf05cp+24 / rt);
  fcTroe = 0x1.9f06f69446738p-1 * exp(-temp / 0x1.15p+8) +
           0x1.83e425aee632p-3 * exp(-temp / 0x1.116p+13) +
           0x1p+0 * exp(-0x1.ed3p+12 / temp);
  k[r318f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.85990f05fba82p+223 *
           exp(-0x1.04283ad1b3b12p+4 * lgt - 0x1.5622f5b3d1d68p+27 / rt);
  kTroeInf = 0x1.fb7f9ba87f719p+43 *
             exp(-0x1.ea0f3d00fbdcp-5 * lgt - 0x1.e5396067a3adp+26 / rt);
  k[r318b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r319f] = 0x1.671e344000001p+36;
  k[r319b] = 0x1.ab0f1ee3cf366p+21 *
             exp(0x1.434367e46c0cp+0 * lgt - 0x1.41453c5ae8b21p+28 / rt);
  kTroe0 = 0x1.cd78333f6505ap+184 *
           exp(-0x1.b170a3d70a3d7p+3 * lgt - 0x1.6a90fep+25 / rt);
  kTroeInf = 0x1.0d30819000001p+35;
  fcTroe = 0x1.5eb851eb851ecp-1 * exp(-temp / 0x1.71p+8) +
           0x1.428f5c28f5c29p-2 * exp(-temp / 0x1.9aap+11) +
           0x1p+0 * exp(-0x1.a0bp+12 / temp);
  k[r320f] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  kTroe0 = 0x1.65c37bd558c9p+205 *
           exp(-0x1.b8916a516ad58p+3 * lgt - 0x1.c14c6be0819c4p+28 / rt);
  kTroeInf = 0x1.a163a262d3cd5p+55 *
             exp(-0x1.c8319e982604p-3 * lgt - 0x1.93fa4c20819c4p+28 / rt);
  k[r320b] = GetLindRateCoeff(temp, pressure, kTroe0, kTroeInf, fcTroe, M[mM2]);
  k[r321f] = 0x1.fb80000000001p+11 *
             exp(0x1.1851eb851eb85p+1 * lgt - 0x1.c69acp+21 / rt);
  k[r321b] = 0x1.22ccf5e51b77cp-10 *
             exp(0x1.dd3fc1ed79c56p+1 * lgt - 0x1.8a8c6959d6793p+25 / rt);
  k[r322f] = 0x1.671e344000001p+34;
  k[r322b] = 0x1.df6c18888e8b3p+25 *
             exp(0x1.21e61d50f4768p-1 * lgt - 0x1.06ba52d213143p+25 / rt);
  k[r323f] = 0x1.8519600000001p+24 *
             exp(0x1.051eb851eb852p-2 * lgt + 0x1.e1ad1p+21 / rt);
  k[r323b] = 0x1.acdf1aa5201fbp+34 *
             exp(0x1.692a708dd8ep-5 * lgt - 0x1.9a0e6f1737f84p+27 / rt);
  k[r324] = 0x1.671e344000001p+34;
  k[r325f] = 0x1.1f25376000001p+34 * exp(-0x1.47ae147ae147bp-2 * lgt);
  k[r325b] = 0x1.7fed34dacc418p+29 *
             exp(0x1.583ddf26dfe8p-5 * lgt - 0x1.545af1ffa64e8p+23 / rt);

  w[r1f] = k[r1f] * c[sO] * c[sO] * M[mM1];
  w[r1b] = k[r1b] * c[sO2] * M[mM1];
  w[r2f] = k[r2f] * c[sO] * c[sH] * M[mM2];
  w[r2b] = k[r2b] * c[sOH] * M[mM2];
  w[r3f] = k[r3f] * c[sO] * c[sH2];
  w[r3b] = k[r3b] * c[sOH] * c[sH];
  w[r4f] = k[r4f] * c[sO] * c[sHO2];
  w[r4b] = k[r4b] * c[sO2] * c[sOH];
  w[r5f] = k[r5f] * c[sO] * c[sH2O2];
  w[r5b] = k[r5b] * c[sHO2] * c[sOH];
  w[r6f] = k[r6f] * c[sO] * c[sCH];
  w[r6b] = k[r6b] * c[sCO] * c[sH];
  w[r7f] = k[r7f] * c[sO] * c[s3XCH2];
  w[r7b] = k[r7b] * c[sHCO] * c[sH];
  w[r8f] = k[r8f] * c[sO] * c[s1XCH2];
  w[r8b] = k[r8b] * c[sCO] * c[sH2];
  w[r9f] = k[r9f] * c[sO] * c[s1XCH2];
  w[r9b] = k[r9b] * c[sHCO] * c[sH];
  w[r10f] = k[r10f] * c[sO] * c[sCH3];
  w[r10b] = k[r10b] * c[sCH2O] * c[sH];
  w[r11f] = k[r11f] * c[sO] * c[sCH4];
  w[r11b] = k[r11b] * c[sCH3] * c[sOH];
  w[r12f] = k[r12f] * c[sO] * c[sCO];
  w[r12b] = k[r12b] * c[sCO2];
  w[r13f] = k[r13f] * c[sO] * c[sHCO];
  w[r13b] = k[r13b] * c[sCO] * c[sOH];
  w[r14f] = k[r14f] * c[sO] * c[sHCO];
  w[r14b] = k[r14b] * c[sCO2] * c[sH];
  w[r15f] = k[r15f] * c[sO] * c[sCH2O];
  w[r15b] = k[r15b] * c[sHCO] * c[sOH];
  w[r16f] = k[r16f] * c[sO] * c[sCH2OH];
  w[r16b] = k[r16b] * c[sCH2O] * c[sOH];
  w[r17f] = k[r17f] * c[sO] * c[sCH3O];
  w[r17b] = k[r17b] * c[sCH2O] * c[sOH];
  w[r18f] = k[r18f] * c[sO] * c[sCH3OH];
  w[r18b] = k[r18b] * c[sCH2OH] * c[sOH];
  w[r19f] = k[r19f] * c[sO] * c[sCH3OH];
  w[r19b] = k[r19b] * c[sCH3O] * c[sOH];
  w[r20f] = k[r20f] * c[sO] * c[sC2H];
  w[r20b] = k[r20b] * c[sCO] * c[sCH];
  w[r21f] = k[r21f] * c[sO] * c[sC2H2];
  w[r21b] = k[r21b] * c[sHCCO] * c[sH];
  w[r22f] = k[r22f] * c[sO] * c[sC2H2];
  w[r22b] = k[r22b] * c[sC2H] * c[sOH];
  w[r23f] = k[r23f] * c[sO] * c[sC2H2];
  w[r23b] = k[r23b] * c[s3XCH2] * c[sCO];
  w[r24f] = k[r24f] * c[sO] * c[sC2H3];
  w[r24b] = k[r24b] * c[sCH2CO] * c[sH];
  w[r25f] = k[r25f] * c[sO] * c[sC2H4];
  w[r25b] = k[r25b] * c[sHCO] * c[sCH3];
  w[r26f] = k[r26f] * c[sO] * c[sC2H5];
  w[r26b] = k[r26b] * c[sCH2O] * c[sCH3];
  w[r27f] = k[r27f] * c[sO] * c[sC2H6];
  w[r27b] = k[r27b] * c[sC2H5] * c[sOH];
  w[r28f] = k[r28f] * c[sO] * c[sHCCO];
  w[r28b] = k[r28b] * c[sCO] * c[sCO] * c[sH];
  w[r29f] = k[r29f] * c[sO] * c[sCH2CO];
  w[r29b] = k[r29b] * c[sHCCO] * c[sOH];
  w[r30f] = k[r30f] * c[sO] * c[sCH2CO];
  w[r30b] = k[r30b] * c[sCO2] * c[s3XCH2];
  w[r31f] = k[r31f] * c[sO2] * c[sCO];
  w[r31b] = k[r31b] * c[sCO2] * c[sO];
  w[r32f] = k[r32f] * c[sO2] * c[sCH2O];
  w[r32b] = k[r32b] * c[sHCO] * c[sHO2];
  w[r33f] = k[r33f] * c[sH] * c[sO2] * M[mM4];
  w[r33b] = k[r33b] * c[sHO2] * M[mM4];
  w[r34f] = k[r34f] * c[sH] * c[sO2] * c[sO2];
  w[r34b] = k[r34b] * c[sO2] * c[sHO2];
  w[r35f] = k[r35f] * c[sH] * c[sO2] * c[sH2O];
  w[r35b] = k[r35b] * c[sH2O] * c[sHO2];
  w[r36f] = k[r36f] * c[sH] * c[sO2] * c[sN2];
  w[r36b] = k[r36b] * c[sN2] * c[sHO2];
  w[r38f] = k[r38f] * c[sH] * c[sO2];
  w[r38b] = k[r38b] * c[sOH] * c[sO];
  w[r39f] = k[r39f] * c[sH] * c[sH] * M[mM5];
  w[r39b] = k[r39b] * c[sH2] * M[mM5];
  w[r40f] = k[r40f] * c[sH] * c[sH] * c[sH2];
  w[r40b] = k[r40b] * c[sH2] * c[sH2];
  w[r41f] = k[r41f] * c[sH] * c[sH] * c[sH2O];
  w[r41b] = k[r41b] * c[sH2O] * c[sH2];
  w[r42f] = k[r42f] * c[sH] * c[sH] * c[sCO2];
  w[r42b] = k[r42b] * c[sCO2] * c[sH2];
  w[r43f] = k[r43f] * c[sH] * c[sOH] * M[mM6];
  w[r43b] = k[r43b] * c[sH2O] * M[mM6];
  w[r44f] = k[r44f] * c[sH] * c[sHO2];
  w[r44b] = k[r44b] * c[sH2O] * c[sO];
  w[r45f] = k[r45f] * c[sH] * c[sHO2];
  w[r45b] = k[r45b] * c[sH2] * c[sO2];
  w[r46f] = k[r46f] * c[sH] * c[sHO2];
  w[r46b] = k[r46b] * c[sOH] * c[sOH];
  w[r47f] = k[r47f] * c[sH] * c[sH2O2];
  w[r47b] = k[r47b] * c[sH2] * c[sHO2];
  w[r48f] = k[r48f] * c[sH] * c[sH2O2];
  w[r48b] = k[r48b] * c[sH2O] * c[sOH];
  w[r49f] = k[r49f] * c[sH] * c[sCH];
  w[r49b] = k[r49b] * c[sH2] * c[sC];
  w[r50f] = k[r50f] * c[sH] * c[s3XCH2];
  w[r50b] = k[r50b] * c[sCH3];
  w[r51f] = k[r51f] * c[sH] * c[s1XCH2];
  w[r51b] = k[r51b] * c[sH2] * c[sCH];
  w[r52f] = k[r52f] * c[sH] * c[sCH3];
  w[r52b] = k[r52b] * c[sCH4];
  w[r53f] = k[r53f] * c[sH] * c[sCH4];
  w[r53b] = k[r53b] * c[sH2] * c[sCH3];
  w[r54f] = k[r54f] * c[sH] * c[sHCO];
  w[r54b] = k[r54b] * c[sCH2O];
  w[r55f] = k[r55f] * c[sH] * c[sHCO];
  w[r55b] = k[r55b] * c[sCO] * c[sH2];
  w[r56f] = k[r56f] * c[sH] * c[sCH2O];
  w[r56b] = k[r56b] * c[sCH2OH];
  w[r57f] = k[r57f] * c[sH] * c[sCH2O];
  w[r57b] = k[r57b] * c[sCH3O];
  w[r58f] = k[r58f] * c[sH] * c[sCH2O];
  w[r58b] = k[r58b] * c[sH2] * c[sHCO];
  w[r59f] = k[r59f] * c[sH] * c[sCH2OH];
  w[r59b] = k[r59b] * c[sCH3OH];
  w[r60f] = k[r60f] * c[sH] * c[sCH2OH];
  w[r60b] = k[r60b] * c[sCH2O] * c[sH2];
  w[r61f] = k[r61f] * c[sH] * c[sCH2OH];
  w[r61b] = k[r61b] * c[sCH3] * c[sOH];
  w[r62f] = k[r62f] * c[sH] * c[sCH2OH];
  w[r62b] = k[r62b] * c[sH2O] * c[s1XCH2];
  w[r63f] = k[r63f] * c[sH] * c[sCH3O];
  w[r63b] = k[r63b] * c[sCH3OH];
  w[r64f] = k[r64f] * c[sH] * c[sCH3O];
  w[r64b] = k[r64b] * c[sCH2OH] * c[sH];
  w[r65f] = k[r65f] * c[sH] * c[sCH3O];
  w[r65b] = k[r65b] * c[sCH2O] * c[sH2];
  w[r66f] = k[r66f] * c[sH] * c[sCH3O];
  w[r66b] = k[r66b] * c[sCH3] * c[sOH];
  w[r67f] = k[r67f] * c[sH] * c[sCH3O];
  w[r67b] = k[r67b] * c[sH2O] * c[s1XCH2];
  w[r68f] = k[r68f] * c[sH] * c[sCH3OH];
  w[r68b] = k[r68b] * c[sH2] * c[sCH2OH];
  w[r69f] = k[r69f] * c[sH] * c[sCH3OH];
  w[r69b] = k[r69b] * c[sH2] * c[sCH3O];
  w[r70f] = k[r70f] * c[sH] * c[sC2H];
  w[r70b] = k[r70b] * c[sC2H2];
  w[r71f] = k[r71f] * c[sH] * c[sC2H2];
  w[r71b] = k[r71b] * c[sC2H3];
  w[r72f] = k[r72f] * c[sH] * c[sC2H3];
  w[r72b] = k[r72b] * c[sC2H4];
  w[r73f] = k[r73f] * c[sH] * c[sC2H3];
  w[r73b] = k[r73b] * c[sC2H2] * c[sH2];
  w[r74f] = k[r74f] * c[sH] * c[sC2H4];
  w[r74b] = k[r74b] * c[sC2H5];
  w[r75f] = k[r75f] * c[sH] * c[sC2H4];
  w[r75b] = k[r75b] * c[sH2] * c[sC2H3];
  w[r76f] = k[r76f] * c[sH] * c[sC2H5];
  w[r76b] = k[r76b] * c[sC2H6];
  w[r77f] = k[r77f] * c[sH] * c[sC2H5];
  w[r77b] = k[r77b] * c[sC2H4] * c[sH2];
  w[r78f] = k[r78f] * c[sH] * c[sC2H6];
  w[r78b] = k[r78b] * c[sH2] * c[sC2H5];
  w[r79f] = k[r79f] * c[sH] * c[sHCCO];
  w[r79b] = k[r79b] * c[sCO] * c[s1XCH2];
  w[r80f] = k[r80f] * c[sH] * c[sCH2CO];
  w[r80b] = k[r80b] * c[sH2] * c[sHCCO];
  w[r81f] = k[r81f] * c[sH] * c[sCH2CO];
  w[r81b] = k[r81b] * c[sCO] * c[sCH3];
  w[r82f] = k[r82f] * c[sH] * c[sHCCOH];
  w[r82b] = k[r82b] * c[sCH2CO] * c[sH];
  w[r83f] = k[r83f] * c[sH2] * c[sCO];
  w[r83b] = k[r83b] * c[sCH2O];
  w[r84f] = k[r84f] * c[sOH] * c[sH2];
  w[r84b] = k[r84b] * c[sH2O] * c[sH];
  w[r85f] = k[r85f] * c[sOH] * c[sOH];
  w[r85b] = k[r85b] * c[sH2O2];
  w[r86f] = k[r86f] * c[sOH] * c[sOH];
  w[r86b] = k[r86b] * c[sH2O] * c[sO];
  w[r87f] = k[r87f] * c[sOH] * c[sHO2];
  w[r87b] = k[r87b] * c[sH2O] * c[sO2];
  w[r88f] = k[r88f] * c[sOH] * c[sH2O2];
  w[r88b] = k[r88b] * c[sH2O] * c[sHO2];
  w[r89f] = k[r89f] * c[sOH] * c[sH2O2];
  w[r89b] = k[r89b] * c[sH2O] * c[sHO2];
  w[r90f] = k[r90f] * c[sOH] * c[sC];
  w[r90b] = k[r90b] * c[sCO] * c[sH];
  w[r91f] = k[r91f] * c[sOH] * c[sCH];
  w[r91b] = k[r91b] * c[sHCO] * c[sH];
  w[r92f] = k[r92f] * c[sOH] * c[s3XCH2];
  w[r92b] = k[r92b] * c[sCH2O] * c[sH];
  w[r93f] = k[r93f] * c[sOH] * c[s3XCH2];
  w[r93b] = k[r93b] * c[sH2O] * c[sCH];
  w[r94f] = k[r94f] * c[sOH] * c[s1XCH2];
  w[r94b] = k[r94b] * c[sCH2O] * c[sH];
  w[r95f] = k[r95f] * c[sOH] * c[sCH3];
  w[r95b] = k[r95b] * c[sCH3OH];
  w[r96f] = k[r96f] * c[sOH] * c[sCH3];
  w[r96b] = k[r96b] * c[sH2O] * c[s3XCH2];
  w[r97f] = k[r97f] * c[sOH] * c[sCH3];
  w[r97b] = k[r97b] * c[sH2O] * c[s1XCH2];
  w[r98f] = k[r98f] * c[sOH] * c[sCH4];
  w[r98b] = k[r98b] * c[sH2O] * c[sCH3];
  w[r99f] = k[r99f] * c[sOH] * c[sCO];
  w[r99b] = k[r99b] * c[sCO2] * c[sH];
  w[r100f] = k[r100f] * c[sOH] * c[sHCO];
  w[r100b] = k[r100b] * c[sCO] * c[sH2O];
  w[r101f] = k[r101f] * c[sOH] * c[sCH2O];
  w[r101b] = k[r101b] * c[sH2O] * c[sHCO];
  w[r102f] = k[r102f] * c[sOH] * c[sCH2OH];
  w[r102b] = k[r102b] * c[sCH2O] * c[sH2O];
  w[r103f] = k[r103f] * c[sOH] * c[sCH3O];
  w[r103b] = k[r103b] * c[sCH2O] * c[sH2O];
  w[r104f] = k[r104f] * c[sOH] * c[sCH3OH];
  w[r104b] = k[r104b] * c[sH2O] * c[sCH2OH];
  w[r105f] = k[r105f] * c[sOH] * c[sCH3OH];
  w[r105b] = k[r105b] * c[sH2O] * c[sCH3O];
  w[r106f] = k[r106f] * c[sOH] * c[sC2H];
  w[r106b] = k[r106b] * c[sHCCO] * c[sH];
  w[r107f] = k[r107f] * c[sOH] * c[sC2H2];
  w[r107b] = k[r107b] * c[sCH2CO] * c[sH];
  w[r108f] = k[r108f] * c[sOH] * c[sC2H2];
  w[r108b] = k[r108b] * c[sHCCOH] * c[sH];
  w[r109f] = k[r109f] * c[sOH] * c[sC2H2];
  w[r109b] = k[r109b] * c[sH2O] * c[sC2H];
  w[r110f] = k[r110f] * c[sOH] * c[sC2H2];
  w[r110b] = k[r110b] * c[sCO] * c[sCH3];
  w[r111f] = k[r111f] * c[sOH] * c[sC2H3];
  w[r111b] = k[r111b] * c[sC2H2] * c[sH2O];
  w[r112f] = k[r112f] * c[sOH] * c[sC2H4];
  w[r112b] = k[r112b] * c[sH2O] * c[sC2H3];
  w[r113f] = k[r113f] * c[sOH] * c[sC2H6];
  w[r113b] = k[r113b] * c[sH2O] * c[sC2H5];
  w[r114f] = k[r114f] * c[sOH] * c[sCH2CO];
  w[r114b] = k[r114b] * c[sH2O] * c[sHCCO];
  w[r115f] = k[r115f] * c[sHO2] * c[sHO2];
  w[r115b] = k[r115b] * c[sH2O2] * c[sO2];
  w[r116f] = k[r116f] * c[sHO2] * c[sHO2];
  w[r116b] = k[r116b] * c[sH2O2] * c[sO2];
  w[r117f] = k[r117f] * c[sHO2] * c[s3XCH2];
  w[r117b] = k[r117b] * c[sCH2O] * c[sOH];
  w[r118f] = k[r118f] * c[sHO2] * c[sCH3];
  w[r118b] = k[r118b] * c[sCH4] * c[sO2];
  w[r119f] = k[r119f] * c[sHO2] * c[sCH3];
  w[r119b] = k[r119b] * c[sCH3O] * c[sOH];
  w[r120f] = k[r120f] * c[sHO2] * c[sCO];
  w[r120b] = k[r120b] * c[sCO2] * c[sOH];
  w[r121f] = k[r121f] * c[sHO2] * c[sCH2O];
  w[r121b] = k[r121b] * c[sH2O2] * c[sHCO];
  w[r122f] = k[r122f] * c[sC] * c[sO2];
  w[r122b] = k[r122b] * c[sCO] * c[sO];
  w[r123f] = k[r123f] * c[sC] * c[s3XCH2];
  w[r123b] = k[r123b] * c[sC2H] * c[sH];
  w[r124f] = k[r124f] * c[sC] * c[sCH3];
  w[r124b] = k[r124b] * c[sC2H2] * c[sH];
  w[r125f] = k[r125f] * c[sCH] * c[sO2];
  w[r125b] = k[r125b] * c[sHCO] * c[sO];
  w[r126f] = k[r126f] * c[sCH] * c[sH2];
  w[r126b] = k[r126b] * c[s3XCH2] * c[sH];
  w[r127f] = k[r127f] * c[sCH] * c[sH2O];
  w[r127b] = k[r127b] * c[sCH2O] * c[sH];
  w[r128f] = k[r128f] * c[sCH] * c[s3XCH2];
  w[r128b] = k[r128b] * c[sC2H2] * c[sH];
  w[r129f] = k[r129f] * c[sCH] * c[sCH3];
  w[r129b] = k[r129b] * c[sC2H3] * c[sH];
  w[r130f] = k[r130f] * c[sCH] * c[sCH4];
  w[r130b] = k[r130b] * c[sC2H4] * c[sH];
  w[r131f] = k[r131f] * c[sCH] * c[sCO];
  w[r131b] = k[r131b] * c[sHCCO];
  w[r132f] = k[r132f] * c[sCH] * c[sCO2];
  w[r132b] = k[r132b] * c[sCO] * c[sHCO];
  w[r133f] = k[r133f] * c[sCH] * c[sCH2O];
  w[r133b] = k[r133b] * c[sCH2CO] * c[sH];
  w[r134f] = k[r134f] * c[sCH] * c[sHCCO];
  w[r134b] = k[r134b] * c[sC2H2] * c[sCO];
  w[r135] = k[r135] * c[s3XCH2] * c[sO2];
  w[r136f] = k[r136f] * c[s3XCH2] * c[sH2];
  w[r136b] = k[r136b] * c[sCH3] * c[sH];
  w[r137f] = k[r137f] * c[s3XCH2] * c[s3XCH2];
  w[r137b] = k[r137b] * c[sC2H2] * c[sH2];
  w[r138f] = k[r138f] * c[s3XCH2] * c[sCH3];
  w[r138b] = k[r138b] * c[sC2H4] * c[sH];
  w[r139f] = k[r139f] * c[s3XCH2] * c[sCH4];
  w[r139b] = k[r139b] * c[sCH3] * c[sCH3];
  w[r140f] = k[r140f] * c[s3XCH2] * c[sCO];
  w[r140b] = k[r140b] * c[sCH2CO];
  w[r141f] = k[r141f] * c[s3XCH2] * c[sHCCO];
  w[r141b] = k[r141b] * c[sCO] * c[sC2H3];
  w[r142f] = k[r142f] * c[s1XCH2] * c[sN2];
  w[r142b] = k[r142b] * c[sN2] * c[s3XCH2];
  w[r144f] = k[r144f] * c[s1XCH2] * c[sO2];
  w[r144b] = k[r144b] * c[sCO] * c[sOH] * c[sH];
  w[r145f] = k[r145f] * c[s1XCH2] * c[sO2];
  w[r145b] = k[r145b] * c[sH2O] * c[sCO];
  w[r146f] = k[r146f] * c[s1XCH2] * c[sH2];
  w[r146b] = k[r146b] * c[sH] * c[sCH3];
  w[r147f] = k[r147f] * c[s1XCH2] * c[sH2O];
  w[r147b] = k[r147b] * c[sCH3OH];
  w[r148f] = k[r148f] * c[s1XCH2] * c[sH2O];
  w[r148b] = k[r148b] * c[sH2O] * c[s3XCH2];
  w[r149f] = k[r149f] * c[s1XCH2] * c[sCH3];
  w[r149b] = k[r149b] * c[sC2H4] * c[sH];
  w[r150f] = k[r150f] * c[s1XCH2] * c[sCH4];
  w[r150b] = k[r150b] * c[sCH3] * c[sCH3];
  w[r151f] = k[r151f] * c[s1XCH2] * c[sCO];
  w[r151b] = k[r151b] * c[sCO] * c[s3XCH2];
  w[r152f] = k[r152f] * c[s1XCH2] * c[sCO2];
  w[r152b] = k[r152b] * c[sCO2] * c[s3XCH2];
  w[r153f] = k[r153f] * c[s1XCH2] * c[sCO2];
  w[r153b] = k[r153b] * c[sCH2O] * c[sCO];
  w[r154f] = k[r154f] * c[s1XCH2] * c[sC2H6];
  w[r154b] = k[r154b] * c[sC2H5] * c[sCH3];
  w[r155f] = k[r155f] * c[sCH3] * c[sO2];
  w[r155b] = k[r155b] * c[sCH3O] * c[sO];
  w[r156f] = k[r156f] * c[sCH3] * c[sO2];
  w[r156b] = k[r156b] * c[sCH2O] * c[sOH];
  w[r157f] = k[r157f] * c[sCH3] * c[sH2O2];
  w[r157b] = k[r157b] * c[sCH4] * c[sHO2];
  w[r158f] = k[r158f] * c[sCH3] * c[sCH3];
  w[r158b] = k[r158b] * c[sC2H6];
  w[r159f] = k[r159f] * c[sCH3] * c[sCH3];
  w[r159b] = k[r159b] * c[sC2H5] * c[sH];
  w[r160f] = k[r160f] * c[sCH3] * c[sHCO];
  w[r160b] = k[r160b] * c[sCO] * c[sCH4];
  w[r161f] = k[r161f] * c[sCH3] * c[sCH2O];
  w[r161b] = k[r161b] * c[sCH4] * c[sHCO];
  w[r162f] = k[r162f] * c[sCH3] * c[sCH3OH];
  w[r162b] = k[r162b] * c[sCH4] * c[sCH2OH];
  w[r163f] = k[r163f] * c[sCH3] * c[sCH3OH];
  w[r163b] = k[r163b] * c[sCH4] * c[sCH3O];
  w[r164f] = k[r164f] * c[sCH3] * c[sC2H4];
  w[r164b] = k[r164b] * c[sCH4] * c[sC2H3];
  w[r165f] = k[r165f] * c[sCH3] * c[sC2H6];
  w[r165b] = k[r165b] * c[sCH4] * c[sC2H5];
  w[r166f] = k[r166f] * c[sHCO] * c[sH2O];
  w[r166b] = k[r166b] * c[sH2O] * c[sCO] * c[sH];
  w[r167f] = k[r167f] * c[sHCO] * M[mM2];
  w[r167b] = k[r167b] * c[sCO] * c[sH] * M[mM2];
  w[r168f] = k[r168f] * c[sHCO] * c[sO2];
  w[r168b] = k[r168b] * c[sCO] * c[sHO2];
  w[r169f] = k[r169f] * c[sCH2OH] * c[sO2];
  w[r169b] = k[r169b] * c[sCH2O] * c[sHO2];
  w[r170f] = k[r170f] * c[sCH3O] * c[sO2];
  w[r170b] = k[r170b] * c[sCH2O] * c[sHO2];
  w[r171f] = k[r171f] * c[sC2H] * c[sO2];
  w[r171b] = k[r171b] * c[sCO] * c[sHCO];
  w[r172f] = k[r172f] * c[sC2H] * c[sH2];
  w[r172b] = k[r172b] * c[sC2H2] * c[sH];
  w[r173f] = k[r173f] * c[sC2H3] * c[sO2];
  w[r173b] = k[r173b] * c[sCH2O] * c[sHCO];
  w[r174f] = k[r174f] * c[sC2H4];
  w[r174b] = k[r174b] * c[sC2H2] * c[sH2];
  w[r175f] = k[r175f] * c[sC2H5] * c[sO2];
  w[r175b] = k[r175b] * c[sC2H4] * c[sHO2];
  w[r176f] = k[r176f] * c[sHCCO] * c[sO2];
  w[r176b] = k[r176b] * c[sCO] * c[sCO] * c[sOH];
  w[r177f] = k[r177f] * c[sHCCO] * c[sHCCO];
  w[r177b] = k[r177b] * c[sC2H2] * c[sCO] * c[sCO];
  w[r178f] = k[r178f] * c[sN] * c[sNO];
  w[r178b] = k[r178b] * c[sO] * c[sN2];
  w[r179f] = k[r179f] * c[sN] * c[sO2];
  w[r179b] = k[r179b] * c[sO] * c[sNO];
  w[r180f] = k[r180f] * c[sN] * c[sOH];
  w[r180b] = k[r180b] * c[sH] * c[sNO];
  w[r181f] = k[r181f] * c[sN2O] * c[sO];
  w[r181b] = k[r181b] * c[sO2] * c[sN2];
  w[r182f] = k[r182f] * c[sN2O] * c[sO];
  w[r182b] = k[r182b] * c[sNO] * c[sNO];
  w[r183f] = k[r183f] * c[sN2O] * c[sH];
  w[r183b] = k[r183b] * c[sOH] * c[sN2];
  w[r184f] = k[r184f] * c[sN2O] * c[sOH];
  w[r184b] = k[r184b] * c[sHO2] * c[sN2];
  w[r185f] = k[r185f] * c[sN2O];
  w[r185b] = k[r185b] * c[sO] * c[sN2];
  w[r186f] = k[r186f] * c[sHO2] * c[sNO];
  w[r186b] = k[r186b] * c[sOH] * c[sNO2];
  w[r187f] = k[r187f] * c[sNO] * c[sO] * M[mM2];
  w[r187b] = k[r187b] * c[sNO2] * M[mM2];
  w[r188f] = k[r188f] * c[sNO2] * c[sO];
  w[r188b] = k[r188b] * c[sO2] * c[sNO];
  w[r189f] = k[r189f] * c[sNO2] * c[sH];
  w[r189b] = k[r189b] * c[sOH] * c[sNO];
  w[r190f] = k[r190f] * c[sNH] * c[sO];
  w[r190b] = k[r190b] * c[sH] * c[sNO];
  w[r191f] = k[r191f] * c[sNH] * c[sH];
  w[r191b] = k[r191b] * c[sH2] * c[sN];
  w[r192f] = k[r192f] * c[sNH] * c[sOH];
  w[r192b] = k[r192b] * c[sH] * c[sHNO];
  w[r193f] = k[r193f] * c[sNH] * c[sOH];
  w[r193b] = k[r193b] * c[sH2O] * c[sN];
  w[r194f] = k[r194f] * c[sNH] * c[sO2];
  w[r194b] = k[r194b] * c[sO] * c[sHNO];
  w[r195f] = k[r195f] * c[sNH] * c[sO2];
  w[r195b] = k[r195b] * c[sOH] * c[sNO];
  w[r196f] = k[r196f] * c[sNH] * c[sN];
  w[r196b] = k[r196b] * c[sH] * c[sN2];
  w[r197f] = k[r197f] * c[sNH] * c[sH2O];
  w[r197b] = k[r197b] * c[sH2] * c[sHNO];
  w[r198f] = k[r198f] * c[sNH] * c[sNO];
  w[r198b] = k[r198b] * c[sOH] * c[sN2];
  w[r199f] = k[r199f] * c[sNH] * c[sNO];
  w[r199b] = k[r199b] * c[sH] * c[sN2O];
  w[r200f] = k[r200f] * c[sNH2] * c[sO];
  w[r200b] = k[r200b] * c[sNH] * c[sOH];
  w[r201f] = k[r201f] * c[sNH2] * c[sO];
  w[r201b] = k[r201b] * c[sHNO] * c[sH];
  w[r202f] = k[r202f] * c[sNH2] * c[sH];
  w[r202b] = k[r202b] * c[sH2] * c[sNH];
  w[r203f] = k[r203f] * c[sNH2] * c[sOH];
  w[r203b] = k[r203b] * c[sH2O] * c[sNH];
  w[r204f] = k[r204f] * c[sNNH];
  w[r204b] = k[r204b] * c[sH] * c[sN2];
  w[r205f] = k[r205f] * c[sNNH] * M[mM2];
  w[r205b] = k[r205b] * c[sH] * c[sN2] * M[mM2];
  w[r206f] = k[r206f] * c[sNNH] * c[sO2];
  w[r206b] = k[r206b] * c[sN2] * c[sHO2];
  w[r207f] = k[r207f] * c[sNNH] * c[sO];
  w[r207b] = k[r207b] * c[sN2] * c[sOH];
  w[r208f] = k[r208f] * c[sNNH] * c[sO];
  w[r208b] = k[r208b] * c[sNO] * c[sNH];
  w[r209f] = k[r209f] * c[sNNH] * c[sH];
  w[r209b] = k[r209b] * c[sN2] * c[sH2];
  w[r210f] = k[r210f] * c[sNNH] * c[sOH];
  w[r210b] = k[r210b] * c[sN2] * c[sH2O];
  w[r211f] = k[r211f] * c[sNNH] * c[sCH3];
  w[r211b] = k[r211b] * c[sN2] * c[sCH4];
  w[r212f] = k[r212f] * c[sH] * c[sNO] * M[mM2];
  w[r212b] = k[r212b] * c[sHNO] * M[mM2];
  w[r213f] = k[r213f] * c[sHNO] * c[sO];
  w[r213b] = k[r213b] * c[sOH] * c[sNO];
  w[r214f] = k[r214f] * c[sHNO] * c[sH];
  w[r214b] = k[r214b] * c[sNO] * c[sH2];
  w[r215f] = k[r215f] * c[sHNO] * c[sOH];
  w[r215b] = k[r215b] * c[sH2O] * c[sNO];
  w[r216f] = k[r216f] * c[sHNO] * c[sO2];
  w[r216b] = k[r216b] * c[sNO] * c[sHO2];
  w[r217f] = k[r217f] * c[sCN] * c[sO];
  w[r217b] = k[r217b] * c[sN] * c[sCO];
  w[r218f] = k[r218f] * c[sCN] * c[sOH];
  w[r218b] = k[r218b] * c[sH] * c[sNCO];
  w[r219f] = k[r219f] * c[sCN] * c[sH2O];
  w[r219b] = k[r219b] * c[sOH] * c[sHCN];
  w[r220f] = k[r220f] * c[sCN] * c[sO2];
  w[r220b] = k[r220b] * c[sO] * c[sNCO];
  w[r221f] = k[r221f] * c[sCN] * c[sH2];
  w[r221b] = k[r221b] * c[sH] * c[sHCN];
  w[r222f] = k[r222f] * c[sNCO] * c[sO];
  w[r222b] = k[r222b] * c[sCO] * c[sNO];
  w[r223f] = k[r223f] * c[sNCO] * c[sH];
  w[r223b] = k[r223b] * c[sCO] * c[sNH];
  w[r224f] = k[r224f] * c[sNCO] * c[sOH];
  w[r224b] = k[r224b] * c[sCO] * c[sH] * c[sNO];
  w[r225f] = k[r225f] * c[sNCO] * c[sN];
  w[r225b] = k[r225b] * c[sCO] * c[sN2];
  w[r226f] = k[r226f] * c[sNCO] * c[sO2];
  w[r226b] = k[r226b] * c[sCO2] * c[sNO];
  w[r227f] = k[r227f] * c[sNCO] * M[mM2];
  w[r227b] = k[r227b] * c[sCO] * c[sN] * M[mM2];
  w[r228f] = k[r228f] * c[sNCO] * c[sNO];
  w[r228b] = k[r228b] * c[sCO] * c[sN2O];
  w[r229f] = k[r229f] * c[sNCO] * c[sNO];
  w[r229b] = k[r229b] * c[sCO2] * c[sN2];
  w[r230f] = k[r230f] * c[sHCN] * M[mM2];
  w[r230b] = k[r230b] * c[sCN] * c[sH] * M[mM2];
  w[r231f] = k[r231f] * c[sHCN] * c[sO];
  w[r231b] = k[r231b] * c[sH] * c[sNCO];
  w[r232f] = k[r232f] * c[sHCN] * c[sO];
  w[r232b] = k[r232b] * c[sCO] * c[sNH];
  w[r233f] = k[r233f] * c[sHCN] * c[sO];
  w[r233b] = k[r233b] * c[sOH] * c[sCN];
  w[r234f] = k[r234f] * c[sHCN] * c[sOH];
  w[r234b] = k[r234b] * c[sH] * c[sHOCN];
  w[r235] = k[r235] * c[sHCN] * c[sOH];
  w[r236f] = k[r236f] * c[sHCN] * c[sOH];
  w[r236b] = k[r236b] * c[sCO] * c[sNH2];
  w[r237f] = k[r237f] * c[sH] * c[sHCN];
  w[r237b] = k[r237b] * c[sH2CN];
  w[r238f] = k[r238f] * c[sH2CN] * c[sN];
  w[r238b] = k[r238b] * c[s3XCH2] * c[sN2];
  w[r239f] = k[r239f] * c[sC] * c[sN2];
  w[r239b] = k[r239b] * c[sN] * c[sCN];
  w[r240f] = k[r240f] * c[sCH] * c[sN2];
  w[r240b] = k[r240b] * c[sN] * c[sHCN];
  w[r241f] = k[r241f] * c[sCH] * c[sN2];
  w[r241b] = k[r241b] * c[sHCNN];
  w[r242f] = k[r242f] * c[s3XCH2] * c[sN2];
  w[r242b] = k[r242b] * c[sNH] * c[sHCN];
  w[r243f] = k[r243f] * c[s1XCH2] * c[sN2];
  w[r243b] = k[r243b] * c[sHCN] * c[sNH];
  w[r244f] = k[r244f] * c[sC] * c[sNO];
  w[r244b] = k[r244b] * c[sO] * c[sCN];
  w[r245f] = k[r245f] * c[sC] * c[sNO];
  w[r245b] = k[r245b] * c[sN] * c[sCO];
  w[r246f] = k[r246f] * c[sCH] * c[sNO];
  w[r246b] = k[r246b] * c[sO] * c[sHCN];
  w[r247f] = k[r247f] * c[sCH] * c[sNO];
  w[r247b] = k[r247b] * c[sNCO] * c[sH];
  w[r248f] = k[r248f] * c[sCH] * c[sNO];
  w[r248b] = k[r248b] * c[sHCO] * c[sN];
  w[r249f] = k[r249f] * c[s3XCH2] * c[sNO];
  w[r249b] = k[r249b] * c[sHNCO] * c[sH];
  w[r250f] = k[r250f] * c[s3XCH2] * c[sNO];
  w[r250b] = k[r250b] * c[sHCN] * c[sOH];
  w[r251f] = k[r251f] * c[s3XCH2] * c[sNO];
  w[r251b] = k[r251b] * c[sHCNO] * c[sH];
  w[r252f] = k[r252f] * c[s1XCH2] * c[sNO];
  w[r252b] = k[r252b] * c[sHNCO] * c[sH];
  w[r253f] = k[r253f] * c[s1XCH2] * c[sNO];
  w[r253b] = k[r253b] * c[sHCN] * c[sOH];
  w[r254f] = k[r254f] * c[s1XCH2] * c[sNO];
  w[r254b] = k[r254b] * c[sHCNO] * c[sH];
  w[r255f] = k[r255f] * c[sCH3] * c[sNO];
  w[r255b] = k[r255b] * c[sH2O] * c[sHCN];
  w[r256f] = k[r256f] * c[sCH3] * c[sNO];
  w[r256b] = k[r256b] * c[sOH] * c[sH2CN];
  w[r257f] = k[r257f] * c[sHCNN] * c[sO];
  w[r257b] = k[r257b] * c[sN2] * c[sH] * c[sCO];
  w[r258f] = k[r258f] * c[sHCNN] * c[sO];
  w[r258b] = k[r258b] * c[sNO] * c[sHCN];
  w[r259f] = k[r259f] * c[sHCNN] * c[sO2];
  w[r259b] = k[r259b] * c[sN2] * c[sHCO] * c[sO];
  w[r260f] = k[r260f] * c[sHCNN] * c[sOH];
  w[r260b] = k[r260b] * c[sN2] * c[sHCO] * c[sH];
  w[r261f] = k[r261f] * c[sHCNN] * c[sH];
  w[r261b] = k[r261b] * c[sN2] * c[s3XCH2];
  w[r262f] = k[r262f] * c[sHNCO] * c[sO];
  w[r262b] = k[r262b] * c[sCO2] * c[sNH];
  w[r263f] = k[r263f] * c[sHNCO] * c[sO];
  w[r263b] = k[r263b] * c[sCO] * c[sHNO];
  w[r264f] = k[r264f] * c[sHNCO] * c[sO];
  w[r264b] = k[r264b] * c[sOH] * c[sNCO];
  w[r265f] = k[r265f] * c[sHNCO] * c[sH];
  w[r265b] = k[r265b] * c[sCO] * c[sNH2];
  w[r266f] = k[r266f] * c[sHNCO] * c[sH];
  w[r266b] = k[r266b] * c[sNCO] * c[sH2];
  w[r267f] = k[r267f] * c[sHNCO] * c[sOH];
  w[r267b] = k[r267b] * c[sH2O] * c[sNCO];
  w[r268f] = k[r268f] * c[sHNCO] * c[sOH];
  w[r268b] = k[r268b] * c[sCO2] * c[sNH2];
  w[r269f] = k[r269f] * c[sHNCO] * M[mM2];
  w[r269b] = k[r269b] * c[sCO] * c[sNH] * M[mM2];
  w[r270f] = k[r270f] * c[sHCNO] * c[sH];
  w[r270b] = k[r270b] * c[sHNCO] * c[sH];
  w[r271f] = k[r271f] * c[sHCNO] * c[sH];
  w[r271b] = k[r271b] * c[sHCN] * c[sOH];
  w[r272f] = k[r272f] * c[sHCNO] * c[sH];
  w[r272b] = k[r272b] * c[sCO] * c[sNH2];
  w[r273] = k[r273] * c[sHOCN] * c[sH];
  w[r274f] = k[r274f] * c[sHCCO] * c[sNO];
  w[r274b] = k[r274b] * c[sCO] * c[sHCNO];
  w[r275f] = k[r275f] * c[sCH3] * c[sN];
  w[r275b] = k[r275b] * c[sH] * c[sH2CN];
  w[r276f] = k[r276f] * c[sCH3] * c[sN];
  w[r276b] = k[r276b] * c[sH2] * c[sHCN];
  w[r277f] = k[r277f] * c[sNH3] * c[sH];
  w[r277b] = k[r277b] * c[sH2] * c[sNH2];
  w[r278f] = k[r278f] * c[sNH3] * c[sOH];
  w[r278b] = k[r278b] * c[sH2O] * c[sNH2];
  w[r279f] = k[r279f] * c[sNH3] * c[sO];
  w[r279b] = k[r279b] * c[sOH] * c[sNH2];
  w[r280f] = k[r280f] * c[sNH] * c[sCO2];
  w[r280b] = k[r280b] * c[sCO] * c[sHNO];
  w[r281f] = k[r281f] * c[sCN] * c[sNO2];
  w[r281b] = k[r281b] * c[sNO] * c[sNCO];
  w[r282f] = k[r282f] * c[sNCO] * c[sNO2];
  w[r282b] = k[r282b] * c[sCO2] * c[sN2O];
  w[r283f] = k[r283f] * c[sN] * c[sCO2];
  w[r283b] = k[r283b] * c[sCO] * c[sNO];
  w[r284] = k[r284] * c[sO] * c[sCH3];
  w[r285f] = k[r285f] * c[sO] * c[sC2H4];
  w[r285b] = k[r285b] * c[sCH2CHO] * c[sH];
  w[r286f] = k[r286f] * c[sO] * c[sC2H5];
  w[r286b] = k[r286b] * c[sCH3CHO] * c[sH];
  w[r287f] = k[r287f] * c[sOH] * c[sHO2];
  w[r287b] = k[r287b] * c[sH2O] * c[sO2];
  w[r288] = k[r288] * c[sOH] * c[sCH3];
  w[r289f] = k[r289f] * c[sCH] * c[sH2];
  w[r289b] = k[r289b] * c[sCH3];
  w[r290] = k[r290] * c[s3XCH2] * c[sO2];
  w[r291f] = k[r291f] * c[s3XCH2] * c[sO2];
  w[r291b] = k[r291b] * c[sCH2O] * c[sO];
  w[r292] = k[r292] * c[s3XCH2] * c[s3XCH2];
  w[r293] = k[r293] * c[s1XCH2] * c[sH2O];
  w[r294f] = k[r294f] * c[sC2H3] * c[sO2];
  w[r294b] = k[r294b] * c[sCH2CHO] * c[sO];
  w[r295f] = k[r295f] * c[sC2H3] * c[sO2];
  w[r295b] = k[r295b] * c[sC2H2] * c[sHO2];
  w[r296f] = k[r296f] * c[sO] * c[sCH3CHO];
  w[r296b] = k[r296b] * c[sCH2CHO] * c[sOH];
  w[r297] = k[r297] * c[sO] * c[sCH3CHO];
  w[r298] = k[r298] * c[sO2] * c[sCH3CHO];
  w[r299f] = k[r299f] * c[sH] * c[sCH3CHO];
  w[r299b] = k[r299b] * c[sH2] * c[sCH2CHO];
  w[r300] = k[r300] * c[sH] * c[sCH3CHO];
  w[r301] = k[r301] * c[sOH] * c[sCH3CHO];
  w[r302] = k[r302] * c[sHO2] * c[sCH3CHO];
  w[r303] = k[r303] * c[sCH3] * c[sCH3CHO];
  w[r304f] = k[r304f] * c[sH] * c[sCH2CO];
  w[r304b] = k[r304b] * c[sCH2CHO];
  w[r305] = k[r305] * c[sO] * c[sCH2CHO];
  w[r306] = k[r306] * c[sO2] * c[sCH2CHO];
  w[r307] = k[r307] * c[sO2] * c[sCH2CHO];
  w[r308f] = k[r308f] * c[sH] * c[sCH2CHO];
  w[r308b] = k[r308b] * c[sHCO] * c[sCH3];
  w[r309f] = k[r309f] * c[sH] * c[sCH2CHO];
  w[r309b] = k[r309b] * c[sH2] * c[sCH2CO];
  w[r310f] = k[r310f] * c[sOH] * c[sCH2CHO];
  w[r310b] = k[r310b] * c[sCH2CO] * c[sH2O];
  w[r311f] = k[r311f] * c[sOH] * c[sCH2CHO];
  w[r311b] = k[r311b] * c[sCH2OH] * c[sHCO];
  w[r312f] = k[r312f] * c[sCH3] * c[sC2H5];
  w[r312b] = k[r312b] * c[sC3H8];
  w[r313f] = k[r313f] * c[sO] * c[sC3H8];
  w[r313b] = k[r313b] * c[sC3H7] * c[sOH];
  w[r314f] = k[r314f] * c[sH] * c[sC3H8];
  w[r314b] = k[r314b] * c[sH2] * c[sC3H7];
  w[r315f] = k[r315f] * c[sOH] * c[sC3H8];
  w[r315b] = k[r315b] * c[sH2O] * c[sC3H7];
  w[r316] = k[r316] * c[sC3H7] * c[sH2O2];
  w[r317f] = k[r317f] * c[sCH3] * c[sC3H8];
  w[r317b] = k[r317b] * c[sCH4] * c[sC3H7];
  w[r318f] = k[r318f] * c[sCH3] * c[sC2H4];
  w[r318b] = k[r318b] * c[sC3H7];
  w[r319f] = k[r319f] * c[sO] * c[sC3H7];
  w[r319b] = k[r319b] * c[sCH2O] * c[sC2H5];
  w[r320f] = k[r320f] * c[sH] * c[sC3H7];
  w[r320b] = k[r320b] * c[sC3H8];
  w[r321f] = k[r321f] * c[sH] * c[sC3H7];
  w[r321b] = k[r321b] * c[sC2H5] * c[sCH3];
  w[r322f] = k[r322f] * c[sOH] * c[sC3H7];
  w[r322b] = k[r322b] * c[sCH2OH] * c[sC2H5];
  w[r323f] = k[r323f] * c[sHO2] * c[sC3H7];
  w[r323b] = k[r323b] * c[sC3H8] * c[sO2];
  w[r324] = k[r324] * c[sHO2] * c[sC3H7];
  w[r325f] = k[r325f] * c[sCH3] * c[sC3H7];
  w[r325b] = k[r325b] * c[sC2H5] * c[sC2H5];

  cdot[sAR] = 0.0;

  cdot[sO] =
      -0x1p+1 * w[r1f] + 0x1p+1 * w[r1b] - w[r2f] + w[r2b] - w[r3f] + w[r3b] -
      w[r4f] + w[r4b] - w[r5f] + w[r5b] - w[r6f] + w[r6b] - w[r7f] + w[r7b] -
      w[r8f] + w[r8b] - w[r9f] + w[r9b] - w[r10f] + w[r10b] - w[r11f] +
      w[r11b] - w[r12f] + w[r12b] - w[r13f] + w[r13b] - w[r14f] + w[r14b] -
      w[r15f] + w[r15b] - w[r16f] + w[r16b] - w[r17f] + w[r17b] - w[r18f] +
      w[r18b] - w[r19f] + w[r19b] - w[r20f] + w[r20b] - w[r21f] + w[r21b] -
      w[r22f] + w[r22b] - w[r23f] + w[r23b] - w[r24f] + w[r24b] - w[r25f] +
      w[r25b] - w[r26f] + w[r26b] - w[r27f] + w[r27b] - w[r28f] + w[r28b] -
      w[r29f] + w[r29b] - w[r30f] + w[r30b] + w[r31f] - w[r31b] + w[r38f] -
      w[r38b] + w[r44f] - w[r44b] + w[r86f] - w[r86b] + w[r122f] - w[r122b] +
      w[r125f] - w[r125b] + w[r155f] - w[r155b] + w[r178f] - w[r178b] +
      w[r179f] - w[r179b] - w[r181f] + w[r181b] - w[r182f] + w[r182b] +
      w[r185f] - w[r185b] - w[r187f] + w[r187b] - w[r188f] + w[r188b] -
      w[r190f] + w[r190b] + w[r194f] - w[r194b] - w[r200f] + w[r200b] -
      w[r201f] + w[r201b] - w[r207f] + w[r207b] - w[r208f] + w[r208b] -
      w[r213f] + w[r213b] - w[r217f] + w[r217b] + w[r220f] - w[r220b] -
      w[r222f] + w[r222b] - w[r231f] + w[r231b] - w[r232f] + w[r232b] -
      w[r233f] + w[r233b] + w[r244f] - w[r244b] + w[r246f] - w[r246b] -
      w[r257f] + w[r257b] - w[r258f] + w[r258b] + w[r259f] - w[r259b] -
      w[r262f] + w[r262b] - w[r263f] + w[r263b] - w[r264f] + w[r264b] -
      w[r279f] + w[r279b] - w[r284] - w[r285f] + w[r285b] - w[r286f] +
      w[r286b] + w[r291f] - w[r291b] + w[r294f] - w[r294b] - w[r296f] +
      w[r296b] - w[r297] - w[r305] - w[r313f] + w[r313b] - w[r319f] + w[r319b];

  cdot[sO2] = w[r1f] - w[r1b] + w[r4f] - w[r4b] - w[r31f] + w[r31b] - w[r32f] +
              w[r32b] - w[r33f] + w[r33b] - 0x1p+1 * w[r34f] + w[r34f] -
              w[r34b] + 0x1p+1 * w[r34b] - w[r35f] + w[r35b] - w[r36f] +
              w[r36b] - w[r38f] + w[r38b] + w[r45f] - w[r45b] + w[r87f] -
              w[r87b] + w[r115f] - w[r115b] + w[r116f] - w[r116b] + w[r118f] -
              w[r118b] - w[r122f] + w[r122b] - w[r125f] + w[r125b] - w[r135] -
              w[r144f] + w[r144b] - w[r145f] + w[r145b] - w[r155f] + w[r155b] -
              w[r156f] + w[r156b] - w[r168f] + w[r168b] - w[r169f] + w[r169b] -
              w[r170f] + w[r170b] - w[r171f] + w[r171b] - w[r173f] + w[r173b] -
              w[r175f] + w[r175b] - w[r176f] + w[r176b] - w[r179f] + w[r179b] +
              w[r181f] - w[r181b] + w[r188f] - w[r188b] - w[r194f] + w[r194b] -
              w[r195f] + w[r195b] - w[r206f] + w[r206b] - w[r216f] + w[r216b] -
              w[r220f] + w[r220b] - w[r226f] + w[r226b] - w[r259f] + w[r259b] +
              w[r287f] - w[r287b] - w[r290] - w[r291f] + w[r291b] - w[r294f] +
              w[r294b] - w[r295f] + w[r295b] - w[r298] - w[r306] - w[r307] +
              w[r323f] - w[r323b];

  cdot[sH] =
      -w[r2f] + w[r2b] + w[r3f] - w[r3b] + w[r6f] - w[r6b] + w[r7f] - w[r7b] +
      w[r9f] - w[r9b] + w[r10f] - w[r10b] + w[r14f] - w[r14b] + w[r21f] -
      w[r21b] + w[r24f] - w[r24b] + w[r28f] - w[r28b] - w[r33f] + w[r33b] -
      w[r34f] + w[r34b] - w[r35f] + w[r35b] - w[r36f] + w[r36b] - w[r38f] +
      w[r38b] - 0x1p+1 * w[r39f] + 0x1p+1 * w[r39b] - 0x1p+1 * w[r40f] +
      0x1p+1 * w[r40b] - 0x1p+1 * w[r41f] + 0x1p+1 * w[r41b] -
      0x1p+1 * w[r42f] + 0x1p+1 * w[r42b] - w[r43f] + w[r43b] - w[r44f] +
      w[r44b] - w[r45f] + w[r45b] - w[r46f] + w[r46b] - w[r47f] + w[r47b] -
      w[r48f] + w[r48b] - w[r49f] + w[r49b] - w[r50f] + w[r50b] - w[r51f] +
      w[r51b] - w[r52f] + w[r52b] - w[r53f] + w[r53b] - w[r54f] + w[r54b] -
      w[r55f] + w[r55b] - w[r56f] + w[r56b] - w[r57f] + w[r57b] - w[r58f] +
      w[r58b] - w[r59f] + w[r59b] - w[r60f] + w[r60b] - w[r61f] + w[r61b] -
      w[r62f] + w[r62b] - w[r63f] + w[r63b] - w[r64f] + w[r64f] - w[r64b] +
      w[r64b] - w[r65f] + w[r65b] - w[r66f] + w[r66b] - w[r67f] + w[r67b] -
      w[r68f] + w[r68b] - w[r69f] + w[r69b] - w[r70f] + w[r70b] - w[r71f] +
      w[r71b] - w[r72f] + w[r72b] - w[r73f] + w[r73b] - w[r74f] + w[r74b] -
      w[r75f] + w[r75b] - w[r76f] + w[r76b] - w[r77f] + w[r77b] - w[r78f] +
      w[r78b] - w[r79f] + w[r79b] - w[r80f] + w[r80b] - w[r81f] + w[r81b] -
      w[r82f] + w[r82f] - w[r82b] + w[r82b] + w[r84f] - w[r84b] + w[r90f] -
      w[r90b] + w[r91f] - w[r91b] + w[r92f] - w[r92b] + w[r94f] - w[r94b] +
      w[r99f] - w[r99b] + w[r106f] - w[r106b] + w[r107f] - w[r107b] + w[r108f] -
      w[r108b] + w[r123f] - w[r123b] + w[r124f] - w[r124b] + w[r126f] -
      w[r126b] + w[r127f] - w[r127b] + w[r128f] - w[r128b] + w[r129f] -
      w[r129b] + w[r130f] - w[r130b] + w[r133f] - w[r133b] + w[r135] +
      w[r136f] - w[r136b] + w[r138f] - w[r138b] + w[r144f] - w[r144b] +
      w[r146f] - w[r146b] + w[r149f] - w[r149b] + w[r159f] - w[r159b] +
      w[r166f] - w[r166b] + w[r167f] - w[r167b] + w[r172f] - w[r172b] +
      w[r180f] - w[r180b] - w[r183f] + w[r183b] - w[r189f] + w[r189b] +
      w[r190f] - w[r190b] - w[r191f] + w[r191b] + w[r192f] - w[r192b] +
      w[r196f] - w[r196b] + w[r199f] - w[r199b] + w[r201f] - w[r201b] -
      w[r202f] + w[r202b] + w[r204f] - w[r204b] + w[r205f] - w[r205b] -
      w[r209f] + w[r209b] - w[r212f] + w[r212b] - w[r214f] + w[r214b] +
      w[r218f] - w[r218b] + w[r221f] - w[r221b] - w[r223f] + w[r223b] +
      w[r224f] - w[r224b] + w[r230f] - w[r230b] + w[r231f] - w[r231b] +
      w[r234f] - w[r234b] + w[r235] - w[r237f] + w[r237b] + w[r247f] -
      w[r247b] + w[r249f] - w[r249b] + w[r251f] - w[r251b] + w[r252f] -
      w[r252b] + w[r254f] - w[r254b] + w[r257f] - w[r257b] + w[r260f] -
      w[r260b] - w[r261f] + w[r261b] - w[r265f] + w[r265b] - w[r266f] +
      w[r266b] - w[r270f] + w[r270f] - w[r270b] + w[r270b] - w[r271f] +
      w[r271b] - w[r272f] + w[r272b] - w[r273] + w[r273] + w[r275f] - w[r275b] -
      w[r277f] + w[r277b] + w[r284] + w[r285f] - w[r285b] + w[r286f] -
      w[r286b] + 0x1p+1 * w[r290] + 0x1p+1 * w[r292] - w[r299f] + w[r299b] -
      w[r300] - w[r304f] + w[r304b] + w[r305] - w[r308f] + w[r308b] - w[r309f] +
      w[r309b] - w[r314f] + w[r314b] - w[r320f] + w[r320b] - w[r321f] +
      w[r321b];

  cdot[sOH] =
      w[r2f] - w[r2b] + w[r3f] - w[r3b] + w[r4f] - w[r4b] + w[r5f] - w[r5b] +
      w[r11f] - w[r11b] + w[r13f] - w[r13b] + w[r15f] - w[r15b] + w[r16f] -
      w[r16b] + w[r17f] - w[r17b] + w[r18f] - w[r18b] + w[r19f] - w[r19b] +
      w[r22f] - w[r22b] + w[r27f] - w[r27b] + w[r29f] - w[r29b] + w[r38f] -
      w[r38b] - w[r43f] + w[r43b] + 0x1p+1 * w[r46f] - 0x1p+1 * w[r46b] +
      w[r48f] - w[r48b] + w[r61f] - w[r61b] + w[r66f] - w[r66b] - w[r84f] +
      w[r84b] - 0x1p+1 * w[r85f] + 0x1p+1 * w[r85b] - 0x1p+1 * w[r86f] +
      0x1p+1 * w[r86b] - w[r87f] + w[r87b] - w[r88f] + w[r88b] - w[r89f] +
      w[r89b] - w[r90f] + w[r90b] - w[r91f] + w[r91b] - w[r92f] + w[r92b] -
      w[r93f] + w[r93b] - w[r94f] + w[r94b] - w[r95f] + w[r95b] - w[r96f] +
      w[r96b] - w[r97f] + w[r97b] - w[r98f] + w[r98b] - w[r99f] + w[r99b] -
      w[r100f] + w[r100b] - w[r101f] + w[r101b] - w[r102f] + w[r102b] -
      w[r103f] + w[r103b] - w[r104f] + w[r104b] - w[r105f] + w[r105b] -
      w[r106f] + w[r106b] - w[r107f] + w[r107b] - w[r108f] + w[r108b] -
      w[r109f] + w[r109b] - w[r110f] + w[r110b] - w[r111f] + w[r111b] -
      w[r112f] + w[r112b] - w[r113f] + w[r113b] - w[r114f] + w[r114b] +
      w[r117f] - w[r117b] + w[r119f] - w[r119b] + w[r120f] - w[r120b] +
      w[r135] + w[r144f] - w[r144b] + w[r156f] - w[r156b] + w[r176f] -
      w[r176b] - w[r180f] + w[r180b] + w[r183f] - w[r183b] - w[r184f] +
      w[r184b] + w[r186f] - w[r186b] + w[r189f] - w[r189b] - w[r192f] +
      w[r192b] - w[r193f] + w[r193b] + w[r195f] - w[r195b] + w[r198f] -
      w[r198b] + w[r200f] - w[r200b] - w[r203f] + w[r203b] + w[r207f] -
      w[r207b] - w[r210f] + w[r210b] + w[r213f] - w[r213b] - w[r215f] +
      w[r215b] - w[r218f] + w[r218b] + w[r219f] - w[r219b] - w[r224f] +
      w[r224b] + w[r233f] - w[r233b] - w[r234f] + w[r234b] - w[r235] -
      w[r236f] + w[r236b] + w[r250f] - w[r250b] + w[r253f] - w[r253b] +
      w[r256f] - w[r256b] - w[r260f] + w[r260b] + w[r264f] - w[r264b] -
      w[r267f] + w[r267b] - w[r268f] + w[r268b] + w[r271f] - w[r271b] -
      w[r278f] + w[r278b] + w[r279f] - w[r279b] - w[r287f] + w[r287b] -
      w[r288] + w[r296f] - w[r296b] + w[r297] - w[r301] + w[r306] + w[r307] -
      w[r310f] + w[r310b] - w[r311f] + w[r311b] + w[r313f] - w[r313b] -
      w[r315f] + w[r315b] - w[r322f] + w[r322b] + w[r324];

  cdot[sH2] =
      -w[r3f] + w[r3b] + w[r8f] - w[r8b] + w[r39f] - w[r39b] - w[r40f] +
      0x1p+1 * w[r40f] - 0x1p+1 * w[r40b] + w[r40b] + w[r41f] - w[r41b] +
      w[r42f] - w[r42b] + w[r45f] - w[r45b] + w[r47f] - w[r47b] + w[r49f] -
      w[r49b] + w[r51f] - w[r51b] + w[r53f] - w[r53b] + w[r55f] - w[r55b] +
      w[r58f] - w[r58b] + w[r60f] - w[r60b] + w[r65f] - w[r65b] + w[r68f] -
      w[r68b] + w[r69f] - w[r69b] + w[r73f] - w[r73b] + w[r75f] - w[r75b] +
      w[r77f] - w[r77b] + w[r78f] - w[r78b] + w[r80f] - w[r80b] - w[r83f] +
      w[r83b] - w[r84f] + w[r84b] - w[r126f] + w[r126b] - w[r136f] + w[r136b] +
      w[r137f] - w[r137b] - w[r146f] + w[r146b] - w[r172f] + w[r172b] +
      w[r174f] - w[r174b] + w[r191f] - w[r191b] + w[r197f] - w[r197b] +
      w[r202f] - w[r202b] + w[r209f] - w[r209b] + w[r214f] - w[r214b] -
      w[r221f] + w[r221b] + w[r266f] - w[r266b] + w[r276f] - w[r276b] +
      w[r277f] - w[r277b] + w[r284] + w[r288] - w[r289f] + w[r289b] + w[r293] +
      w[r299f] - w[r299b] + w[r300] + w[r309f] - w[r309b] + w[r314f] - w[r314b];

  cdot[sHO2] = -w[r4f] + w[r4b] + w[r5f] - w[r5b] + w[r32f] - w[r32b] +
               w[r33f] - w[r33b] + w[r34f] - w[r34b] + w[r35f] - w[r35b] +
               w[r36f] - w[r36b] - w[r44f] + w[r44b] - w[r45f] + w[r45b] -
               w[r46f] + w[r46b] + w[r47f] - w[r47b] - w[r87f] + w[r87b] +
               w[r88f] - w[r88b] + w[r89f] - w[r89b] - 0x1p+1 * w[r115f] +
               0x1p+1 * w[r115b] - 0x1p+1 * w[r116f] + 0x1p+1 * w[r116b] -
               w[r117f] + w[r117b] - w[r118f] + w[r118b] - w[r119f] + w[r119b] -
               w[r120f] + w[r120b] - w[r121f] + w[r121b] + w[r157f] - w[r157b] +
               w[r168f] - w[r168b] + w[r169f] - w[r169b] + w[r170f] - w[r170b] +
               w[r175f] - w[r175b] + w[r184f] - w[r184b] - w[r186f] + w[r186b] +
               w[r206f] - w[r206b] + w[r216f] - w[r216b] - w[r287f] + w[r287b] +
               w[r295f] - w[r295b] + w[r298] - w[r302] + w[r316] - w[r323f] +
               w[r323b] - w[r324];

  cdot[sH2O2] = -w[r5f] + w[r5b] - w[r47f] + w[r47b] - w[r48f] + w[r48b] +
                w[r85f] - w[r85b] - w[r88f] + w[r88b] - w[r89f] + w[r89b] +
                w[r115f] - w[r115b] + w[r116f] - w[r116b] + w[r121f] -
                w[r121b] - w[r157f] + w[r157b] + w[r302] - w[r316];

  cdot[sCH] = -w[r6f] + w[r6b] + w[r20f] - w[r20b] - w[r49f] + w[r49b] +
              w[r51f] - w[r51b] - w[r91f] + w[r91b] + w[r93f] - w[r93b] -
              w[r125f] + w[r125b] - w[r126f] + w[r126b] - w[r127f] + w[r127b] -
              w[r128f] + w[r128b] - w[r129f] + w[r129b] - w[r130f] + w[r130b] -
              w[r131f] + w[r131b] - w[r132f] + w[r132b] - w[r133f] + w[r133b] -
              w[r134f] + w[r134b] - w[r240f] + w[r240b] - w[r241f] + w[r241b] -
              w[r246f] + w[r246b] - w[r247f] + w[r247b] - w[r248f] + w[r248b] -
              w[r289f] + w[r289b];

  cdot[sCO] =
      w[r6f] - w[r6b] + w[r8f] - w[r8b] - w[r12f] + w[r12b] + w[r13f] -
      w[r13b] + w[r20f] - w[r20b] + w[r23f] - w[r23b] + 0x1p+1 * w[r28f] -
      0x1p+1 * w[r28b] - w[r31f] + w[r31b] + w[r55f] - w[r55b] + w[r79f] -
      w[r79b] + w[r81f] - w[r81b] - w[r83f] + w[r83b] + w[r90f] - w[r90b] -
      w[r99f] + w[r99b] + w[r100f] - w[r100b] + w[r110f] - w[r110b] - w[r120f] +
      w[r120b] + w[r122f] - w[r122b] - w[r131f] + w[r131b] + w[r132f] -
      w[r132b] + w[r134f] - w[r134b] + w[r135] - w[r140f] + w[r140b] +
      w[r141f] - w[r141b] + w[r144f] - w[r144b] + w[r145f] - w[r145b] -
      w[r151f] + w[r151f] - w[r151b] + w[r151b] + w[r153f] - w[r153b] +
      w[r160f] - w[r160b] + w[r166f] - w[r166b] + w[r167f] - w[r167b] +
      w[r168f] - w[r168b] + w[r171f] - w[r171b] + 0x1p+1 * w[r176f] -
      0x1p+1 * w[r176b] + 0x1p+1 * w[r177f] - 0x1p+1 * w[r177b] + w[r217f] -
      w[r217b] + w[r222f] - w[r222b] + w[r223f] - w[r223b] + w[r224f] -
      w[r224b] + w[r225f] - w[r225b] + w[r227f] - w[r227b] + w[r228f] -
      w[r228b] + w[r232f] - w[r232b] + w[r236f] - w[r236b] + w[r245f] -
      w[r245b] + w[r257f] - w[r257b] + w[r263f] - w[r263b] + w[r265f] -
      w[r265b] + w[r269f] - w[r269b] + w[r272f] - w[r272b] + w[r274f] -
      w[r274b] + w[r280f] - w[r280b] + w[r283f] - w[r283b] + w[r284] + w[r297] +
      w[r298] + w[r300] + w[r301] + w[r302] + w[r303] + w[r306];

  cdot[s3XCH2] =
      -w[r7f] + w[r7b] + w[r23f] - w[r23b] + w[r30f] - w[r30b] - w[r50f] +
      w[r50b] - w[r92f] + w[r92b] - w[r93f] + w[r93b] + w[r96f] - w[r96b] -
      w[r117f] + w[r117b] - w[r123f] + w[r123b] + w[r126f] - w[r126b] -
      w[r128f] + w[r128b] - w[r135] - w[r136f] + w[r136b] - 0x1p+1 * w[r137f] +
      0x1p+1 * w[r137b] - w[r138f] + w[r138b] - w[r139f] + w[r139b] - w[r140f] +
      w[r140b] - w[r141f] + w[r141b] + w[r142f] - w[r142b] + w[r148f] -
      w[r148b] + w[r151f] - w[r151b] + w[r152f] - w[r152b] + w[r238f] -
      w[r238b] - w[r242f] + w[r242b] - w[r249f] + w[r249b] - w[r250f] +
      w[r250b] - w[r251f] + w[r251b] + w[r261f] - w[r261b] - w[r290] -
      w[r291f] + w[r291b] - 0x1p+1 * w[r292] + w[r305];

  cdot[sHCO] = w[r7f] - w[r7b] + w[r9f] - w[r9b] - w[r13f] + w[r13b] - w[r14f] +
               w[r14b] + w[r15f] - w[r15b] + w[r25f] - w[r25b] + w[r32f] -
               w[r32b] - w[r54f] + w[r54b] - w[r55f] + w[r55b] + w[r58f] -
               w[r58b] + w[r91f] - w[r91b] - w[r100f] + w[r100b] + w[r101f] -
               w[r101b] + w[r121f] - w[r121b] + w[r125f] - w[r125b] + w[r132f] -
               w[r132b] - w[r160f] + w[r160b] + w[r161f] - w[r161b] - w[r166f] +
               w[r166b] - w[r167f] + w[r167b] - w[r168f] + w[r168b] + w[r171f] -
               w[r171b] + w[r173f] - w[r173b] + w[r248f] - w[r248b] + w[r259f] -
               w[r259b] + w[r260f] - w[r260b] + 0x1p+1 * w[r307] + w[r308f] -
               w[r308b] + w[r311f] - w[r311b];

  cdot[s1XCH2] = -w[r8f] + w[r8b] - w[r9f] + w[r9b] - w[r51f] + w[r51b] +
                 w[r62f] - w[r62b] + w[r67f] - w[r67b] + w[r79f] - w[r79b] -
                 w[r94f] + w[r94b] + w[r97f] - w[r97b] - w[r142f] + w[r142b] -
                 w[r144f] + w[r144b] - w[r145f] + w[r145b] - w[r146f] +
                 w[r146b] - w[r147f] + w[r147b] - w[r148f] + w[r148b] -
                 w[r149f] + w[r149b] - w[r150f] + w[r150b] - w[r151f] +
                 w[r151b] - w[r152f] + w[r152b] - w[r153f] + w[r153b] -
                 w[r154f] + w[r154b] - w[r243f] + w[r243b] - w[r252f] +
                 w[r252b] - w[r253f] + w[r253b] - w[r254f] + w[r254b] - w[r293];

  cdot[sCH3] =
      -w[r10f] + w[r10b] + w[r11f] - w[r11b] + w[r25f] - w[r25b] + w[r26f] -
      w[r26b] + w[r50f] - w[r50b] - w[r52f] + w[r52b] + w[r53f] - w[r53b] +
      w[r61f] - w[r61b] + w[r66f] - w[r66b] + w[r81f] - w[r81b] - w[r95f] +
      w[r95b] - w[r96f] + w[r96b] - w[r97f] + w[r97b] + w[r98f] - w[r98b] +
      w[r110f] - w[r110b] - w[r118f] + w[r118b] - w[r119f] + w[r119b] -
      w[r124f] + w[r124b] - w[r129f] + w[r129b] + w[r136f] - w[r136b] -
      w[r138f] + w[r138b] + 0x1p+1 * w[r139f] - 0x1p+1 * w[r139b] + w[r146f] -
      w[r146b] - w[r149f] + w[r149b] + 0x1p+1 * w[r150f] - 0x1p+1 * w[r150b] +
      w[r154f] - w[r154b] - w[r155f] + w[r155b] - w[r156f] + w[r156b] -
      w[r157f] + w[r157b] - 0x1p+1 * w[r158f] + 0x1p+1 * w[r158b] -
      0x1p+1 * w[r159f] + 0x1p+1 * w[r159b] - w[r160f] + w[r160b] - w[r161f] +
      w[r161b] - w[r162f] + w[r162b] - w[r163f] + w[r163b] - w[r164f] +
      w[r164b] - w[r165f] + w[r165b] - w[r211f] + w[r211b] - w[r255f] +
      w[r255b] - w[r256f] + w[r256b] - w[r275f] + w[r275b] - w[r276f] +
      w[r276b] - w[r284] - w[r288] + w[r289f] - w[r289b] + w[r297] + w[r298] +
      w[r300] + w[r301] + w[r302] - w[r303] + w[r303] + w[r308f] - w[r308b] -
      w[r312f] + w[r312b] - w[r317f] + w[r317b] - w[r318f] + w[r318b] +
      w[r321f] - w[r321b] - w[r325f] + w[r325b];

  cdot[sCH2O] =
      w[r10f] - w[r10b] - w[r15f] + w[r15b] + w[r16f] - w[r16b] + w[r17f] -
      w[r17b] + w[r26f] - w[r26b] - w[r32f] + w[r32b] + w[r54f] - w[r54b] -
      w[r56f] + w[r56b] - w[r57f] + w[r57b] - w[r58f] + w[r58b] + w[r60f] -
      w[r60b] + w[r65f] - w[r65b] + w[r83f] - w[r83b] + w[r92f] - w[r92b] +
      w[r94f] - w[r94b] - w[r101f] + w[r101b] + w[r102f] - w[r102b] + w[r103f] -
      w[r103b] + w[r117f] - w[r117b] - w[r121f] + w[r121b] + w[r127f] -
      w[r127b] - w[r133f] + w[r133b] + w[r153f] - w[r153b] + w[r156f] -
      w[r156b] - w[r161f] + w[r161b] + w[r169f] - w[r169b] + w[r170f] -
      w[r170b] + w[r173f] - w[r173b] + w[r288] + w[r291f] - w[r291b] + w[r293] +
      w[r306] + w[r319f] - w[r319b] + w[r324];

  cdot[sCH4] = -w[r11f] + w[r11b] + w[r52f] - w[r52b] - w[r53f] + w[r53b] -
               w[r98f] + w[r98b] + w[r118f] - w[r118b] - w[r130f] + w[r130b] -
               w[r139f] + w[r139b] - w[r150f] + w[r150b] + w[r157f] - w[r157b] +
               w[r160f] - w[r160b] + w[r161f] - w[r161b] + w[r162f] - w[r162b] +
               w[r163f] - w[r163b] + w[r164f] - w[r164b] + w[r165f] - w[r165b] +
               w[r211f] - w[r211b] + w[r303] + w[r317f] - w[r317b];

  cdot[sCO2] = w[r12f] - w[r12b] + w[r14f] - w[r14b] + w[r30f] - w[r30b] +
               w[r31f] - w[r31b] - w[r42f] + w[r42f] - w[r42b] + w[r42b] +
               w[r99f] - w[r99b] + w[r120f] - w[r120b] - w[r132f] + w[r132b] -
               w[r152f] + w[r152f] - w[r152b] + w[r152b] - w[r153f] + w[r153b] +
               w[r226f] - w[r226b] + w[r229f] - w[r229b] + w[r262f] - w[r262b] +
               w[r268f] - w[r268b] - w[r280f] + w[r280b] + w[r282f] - w[r282b] -
               w[r283f] + w[r283b] + w[r290] + w[r305];

  cdot[sCH2OH] = -w[r16f] + w[r16b] + w[r18f] - w[r18b] + w[r56f] - w[r56b] -
                 w[r59f] + w[r59b] - w[r60f] + w[r60b] - w[r61f] + w[r61b] -
                 w[r62f] + w[r62b] + w[r64f] - w[r64b] + w[r68f] - w[r68b] -
                 w[r102f] + w[r102b] + w[r104f] - w[r104b] + w[r162f] -
                 w[r162b] - w[r169f] + w[r169b] + w[r311f] - w[r311b] +
                 w[r322f] - w[r322b];

  cdot[sCH3O] = -w[r17f] + w[r17b] + w[r19f] - w[r19b] + w[r57f] - w[r57b] -
                w[r63f] + w[r63b] - w[r64f] + w[r64b] - w[r65f] + w[r65b] -
                w[r66f] + w[r66b] - w[r67f] + w[r67b] + w[r69f] - w[r69b] -
                w[r103f] + w[r103b] + w[r105f] - w[r105b] + w[r119f] -
                w[r119b] + w[r155f] - w[r155b] + w[r163f] - w[r163b] -
                w[r170f] + w[r170b];

  cdot[sCH3OH] = -w[r18f] + w[r18b] - w[r19f] + w[r19b] + w[r59f] - w[r59b] +
                 w[r63f] - w[r63b] - w[r68f] + w[r68b] - w[r69f] + w[r69b] +
                 w[r95f] - w[r95b] - w[r104f] + w[r104b] - w[r105f] + w[r105b] +
                 w[r147f] - w[r147b] - w[r162f] + w[r162b] - w[r163f] +
                 w[r163b];

  cdot[sC2H] = -w[r20f] + w[r20b] + w[r22f] - w[r22b] - w[r70f] + w[r70b] -
               w[r106f] + w[r106b] + w[r109f] - w[r109b] + w[r123f] - w[r123b] -
               w[r171f] + w[r171b] - w[r172f] + w[r172b];

  cdot[sC2H2] =
      -w[r21f] + w[r21b] - w[r22f] + w[r22b] - w[r23f] + w[r23b] + w[r70f] -
      w[r70b] - w[r71f] + w[r71b] + w[r73f] - w[r73b] - w[r107f] + w[r107b] -
      w[r108f] + w[r108b] - w[r109f] + w[r109b] - w[r110f] + w[r110b] +
      w[r111f] - w[r111b] + w[r124f] - w[r124b] + w[r128f] - w[r128b] +
      w[r134f] - w[r134b] + w[r137f] - w[r137b] + w[r172f] - w[r172b] +
      w[r174f] - w[r174b] + w[r177f] - w[r177b] + w[r292] + w[r295f] - w[r295b];

  cdot[sHCCO] = w[r21f] - w[r21b] - w[r28f] + w[r28b] + w[r29f] - w[r29b] -
                w[r79f] + w[r79b] + w[r80f] - w[r80b] + w[r106f] - w[r106b] +
                w[r114f] - w[r114b] + w[r131f] - w[r131b] - w[r134f] +
                w[r134b] - w[r141f] + w[r141b] - w[r176f] + w[r176b] -
                0x1p+1 * w[r177f] + 0x1p+1 * w[r177b] - w[r274f] + w[r274b];

  cdot[sC2H3] = -w[r24f] + w[r24b] + w[r71f] - w[r71b] - w[r72f] + w[r72b] -
                w[r73f] + w[r73b] + w[r75f] - w[r75b] - w[r111f] + w[r111b] +
                w[r112f] - w[r112b] + w[r129f] - w[r129b] + w[r141f] -
                w[r141b] + w[r164f] - w[r164b] - w[r173f] + w[r173b] -
                w[r294f] + w[r294b] - w[r295f] + w[r295b];

  cdot[sCH2CO] = w[r24f] - w[r24b] - w[r29f] + w[r29b] - w[r30f] + w[r30b] -
                 w[r80f] + w[r80b] - w[r81f] + w[r81b] + w[r82f] - w[r82b] +
                 w[r107f] - w[r107b] - w[r114f] + w[r114b] + w[r133f] -
                 w[r133b] + w[r140f] - w[r140b] - w[r304f] + w[r304b] +
                 w[r309f] - w[r309b] + w[r310f] - w[r310b];

  cdot[sC2H4] = -w[r25f] + w[r25b] + w[r72f] - w[r72b] - w[r74f] + w[r74b] -
                w[r75f] + w[r75b] + w[r77f] - w[r77b] - w[r112f] + w[r112b] +
                w[r130f] - w[r130b] + w[r138f] - w[r138b] + w[r149f] -
                w[r149b] - w[r164f] + w[r164b] - w[r174f] + w[r174b] +
                w[r175f] - w[r175b] - w[r285f] + w[r285b] - w[r318f] + w[r318b];

  cdot[sC2H5] = -w[r26f] + w[r26b] + w[r27f] - w[r27b] + w[r74f] - w[r74b] -
                w[r76f] + w[r76b] - w[r77f] + w[r77b] + w[r78f] - w[r78b] +
                w[r113f] - w[r113b] + w[r154f] - w[r154b] + w[r159f] -
                w[r159b] + w[r165f] - w[r165b] - w[r175f] + w[r175b] -
                w[r286f] + w[r286b] - w[r312f] + w[r312b] + w[r319f] -
                w[r319b] + w[r321f] - w[r321b] + w[r322f] - w[r322b] + w[r324] +
                0x1p+1 * w[r325f] - 0x1p+1 * w[r325b];

  cdot[sC2H6] = -w[r27f] + w[r27b] + w[r76f] - w[r76b] - w[r78f] + w[r78b] -
                w[r113f] + w[r113b] - w[r154f] + w[r154b] + w[r158f] -
                w[r158b] - w[r165f] + w[r165b];

  cdot[sH2O] = -w[r35f] + w[r35f] - w[r35b] + w[r35b] - w[r41f] + w[r41f] -
               w[r41b] + w[r41b] + w[r43f] - w[r43b] + w[r44f] - w[r44b] +
               w[r48f] - w[r48b] + w[r62f] - w[r62b] + w[r67f] - w[r67b] +
               w[r84f] - w[r84b] + w[r86f] - w[r86b] + w[r87f] - w[r87b] +
               w[r88f] - w[r88b] + w[r89f] - w[r89b] + w[r93f] - w[r93b] +
               w[r96f] - w[r96b] + w[r97f] - w[r97b] + w[r98f] - w[r98b] +
               w[r100f] - w[r100b] + w[r101f] - w[r101b] + w[r102f] - w[r102b] +
               w[r103f] - w[r103b] + w[r104f] - w[r104b] + w[r105f] - w[r105b] +
               w[r109f] - w[r109b] + w[r111f] - w[r111b] + w[r112f] - w[r112b] +
               w[r113f] - w[r113b] + w[r114f] - w[r114b] - w[r127f] + w[r127b] +
               w[r145f] - w[r145b] - w[r147f] + w[r147b] - w[r148f] + w[r148f] -
               w[r148b] + w[r148b] - w[r166f] + w[r166f] - w[r166b] + w[r166b] +
               w[r193f] - w[r193b] - w[r197f] + w[r197b] + w[r203f] - w[r203b] +
               w[r210f] - w[r210b] + w[r215f] - w[r215b] - w[r219f] + w[r219b] +
               w[r255f] - w[r255b] + w[r267f] - w[r267b] + w[r278f] - w[r278b] +
               w[r287f] - w[r287b] - w[r293] + w[r301] + w[r310f] - w[r310b] +
               w[r315f] - w[r315b];

  cdot[sN2] = -w[r36f] + w[r36f] - w[r36b] + w[r36b] - w[r142f] + w[r142f] -
              w[r142b] + w[r142b] + w[r178f] - w[r178b] + w[r181f] - w[r181b] +
              w[r183f] - w[r183b] + w[r184f] - w[r184b] + w[r185f] - w[r185b] +
              w[r196f] - w[r196b] + w[r198f] - w[r198b] + w[r204f] - w[r204b] +
              w[r205f] - w[r205b] + w[r206f] - w[r206b] + w[r207f] - w[r207b] +
              w[r209f] - w[r209b] + w[r210f] - w[r210b] + w[r211f] - w[r211b] +
              w[r225f] - w[r225b] + w[r229f] - w[r229b] + w[r238f] - w[r238b] -
              w[r239f] + w[r239b] - w[r240f] + w[r240b] - w[r241f] + w[r241b] -
              w[r242f] + w[r242b] - w[r243f] + w[r243b] + w[r257f] - w[r257b] +
              w[r259f] - w[r259b] + w[r260f] - w[r260b] + w[r261f] - w[r261b];

  cdot[sC] = w[r49f] - w[r49b] - w[r90f] + w[r90b] - w[r122f] + w[r122b] -
             w[r123f] + w[r123b] - w[r124f] + w[r124b] - w[r239f] + w[r239b] -
             w[r244f] + w[r244b] - w[r245f] + w[r245b];

  cdot[sHCCOH] = -w[r82f] + w[r82b] + w[r108f] - w[r108b];

  cdot[sN] = -w[r178f] + w[r178b] - w[r179f] + w[r179b] - w[r180f] + w[r180b] +
             w[r191f] - w[r191b] + w[r193f] - w[r193b] - w[r196f] + w[r196b] +
             w[r217f] - w[r217b] - w[r225f] + w[r225b] + w[r227f] - w[r227b] -
             w[r238f] + w[r238b] + w[r239f] - w[r239b] + w[r240f] - w[r240b] +
             w[r245f] - w[r245b] + w[r248f] - w[r248b] - w[r275f] + w[r275b] -
             w[r276f] + w[r276b] - w[r283f] + w[r283b];

  cdot[sNO] = -w[r178f] + w[r178b] + w[r179f] - w[r179b] + w[r180f] - w[r180b] +
              0x1p+1 * w[r182f] - 0x1p+1 * w[r182b] - w[r186f] + w[r186b] -
              w[r187f] + w[r187b] + w[r188f] - w[r188b] + w[r189f] - w[r189b] +
              w[r190f] - w[r190b] + w[r195f] - w[r195b] - w[r198f] + w[r198b] -
              w[r199f] + w[r199b] + w[r208f] - w[r208b] - w[r212f] + w[r212b] +
              w[r213f] - w[r213b] + w[r214f] - w[r214b] + w[r215f] - w[r215b] +
              w[r216f] - w[r216b] + w[r222f] - w[r222b] + w[r224f] - w[r224b] +
              w[r226f] - w[r226b] - w[r228f] + w[r228b] - w[r229f] + w[r229b] -
              w[r244f] + w[r244b] - w[r245f] + w[r245b] - w[r246f] + w[r246b] -
              w[r247f] + w[r247b] - w[r248f] + w[r248b] - w[r249f] + w[r249b] -
              w[r250f] + w[r250b] - w[r251f] + w[r251b] - w[r252f] + w[r252b] -
              w[r253f] + w[r253b] - w[r254f] + w[r254b] - w[r255f] + w[r255b] -
              w[r256f] + w[r256b] + w[r258f] - w[r258b] - w[r274f] + w[r274b] +
              w[r281f] - w[r281b] + w[r283f] - w[r283b];

  cdot[sN2O] = -w[r181f] + w[r181b] - w[r182f] + w[r182b] - w[r183f] +
               w[r183b] - w[r184f] + w[r184b] - w[r185f] + w[r185b] + w[r199f] -
               w[r199b] + w[r228f] - w[r228b] + w[r282f] - w[r282b];

  cdot[sNO2] = w[r186f] - w[r186b] + w[r187f] - w[r187b] - w[r188f] + w[r188b] -
               w[r189f] + w[r189b] - w[r281f] + w[r281b] - w[r282f] + w[r282b];

  cdot[sNH] = -w[r190f] + w[r190b] - w[r191f] + w[r191b] - w[r192f] + w[r192b] -
              w[r193f] + w[r193b] - w[r194f] + w[r194b] - w[r195f] + w[r195b] -
              w[r196f] + w[r196b] - w[r197f] + w[r197b] - w[r198f] + w[r198b] -
              w[r199f] + w[r199b] + w[r200f] - w[r200b] + w[r202f] - w[r202b] +
              w[r203f] - w[r203b] + w[r208f] - w[r208b] + w[r223f] - w[r223b] +
              w[r232f] - w[r232b] + w[r242f] - w[r242b] + w[r243f] - w[r243b] +
              w[r262f] - w[r262b] + w[r269f] - w[r269b] - w[r280f] + w[r280b];

  cdot[sHNO] = w[r192f] - w[r192b] + w[r194f] - w[r194b] + w[r197f] - w[r197b] +
               w[r201f] - w[r201b] + w[r212f] - w[r212b] - w[r213f] + w[r213b] -
               w[r214f] + w[r214b] - w[r215f] + w[r215b] - w[r216f] + w[r216b] +
               w[r263f] - w[r263b] + w[r280f] - w[r280b];

  cdot[sNH2] = -w[r200f] + w[r200b] - w[r201f] + w[r201b] - w[r202f] +
               w[r202b] - w[r203f] + w[r203b] + w[r236f] - w[r236b] + w[r265f] -
               w[r265b] + w[r268f] - w[r268b] + w[r272f] - w[r272b] + w[r277f] -
               w[r277b] + w[r278f] - w[r278b] + w[r279f] - w[r279b];

  cdot[sNNH] = -w[r204f] + w[r204b] - w[r205f] + w[r205b] - w[r206f] +
               w[r206b] - w[r207f] + w[r207b] - w[r208f] + w[r208b] - w[r209f] +
               w[r209b] - w[r210f] + w[r210b] - w[r211f] + w[r211b];

  cdot[sCN] = -w[r217f] + w[r217b] - w[r218f] + w[r218b] - w[r219f] + w[r219b] -
              w[r220f] + w[r220b] - w[r221f] + w[r221b] + w[r230f] - w[r230b] +
              w[r233f] - w[r233b] + w[r239f] - w[r239b] + w[r244f] - w[r244b] -
              w[r281f] + w[r281b];

  cdot[sNCO] = w[r218f] - w[r218b] + w[r220f] - w[r220b] - w[r222f] + w[r222b] -
               w[r223f] + w[r223b] - w[r224f] + w[r224b] - w[r225f] + w[r225b] -
               w[r226f] + w[r226b] - w[r227f] + w[r227b] - w[r228f] + w[r228b] -
               w[r229f] + w[r229b] + w[r231f] - w[r231b] + w[r247f] - w[r247b] +
               w[r264f] - w[r264b] + w[r266f] - w[r266b] + w[r267f] - w[r267b] +
               w[r281f] - w[r281b] - w[r282f] + w[r282b];

  cdot[sHCN] = w[r219f] - w[r219b] + w[r221f] - w[r221b] - w[r230f] + w[r230b] -
               w[r231f] + w[r231b] - w[r232f] + w[r232b] - w[r233f] + w[r233b] -
               w[r234f] + w[r234b] - w[r235] - w[r236f] + w[r236b] - w[r237f] +
               w[r237b] + w[r240f] - w[r240b] + w[r242f] - w[r242b] + w[r243f] -
               w[r243b] + w[r246f] - w[r246b] + w[r250f] - w[r250b] + w[r253f] -
               w[r253b] + w[r255f] - w[r255b] + w[r258f] - w[r258b] + w[r271f] -
               w[r271b] + w[r276f] - w[r276b];

  cdot[sHOCN] = w[r234f] - w[r234b] - w[r273];

  cdot[sHNCO] = w[r235] + w[r249f] - w[r249b] + w[r252f] - w[r252b] - w[r262f] +
                w[r262b] - w[r263f] + w[r263b] - w[r264f] + w[r264b] -
                w[r265f] + w[r265b] - w[r266f] + w[r266b] - w[r267f] +
                w[r267b] - w[r268f] + w[r268b] - w[r269f] + w[r269b] +
                w[r270f] - w[r270b] + w[r273];

  cdot[sH2CN] = w[r237f] - w[r237b] - w[r238f] + w[r238b] + w[r256f] -
                w[r256b] + w[r275f] - w[r275b];

  cdot[sHCNN] = w[r241f] - w[r241b] - w[r257f] + w[r257b] - w[r258f] +
                w[r258b] - w[r259f] + w[r259b] - w[r260f] + w[r260b] -
                w[r261f] + w[r261b];

  cdot[sHCNO] = w[r251f] - w[r251b] + w[r254f] - w[r254b] - w[r270f] +
                w[r270b] - w[r271f] + w[r271b] - w[r272f] + w[r272b] +
                w[r274f] - w[r274b];

  cdot[sNH3] = -w[r277f] + w[r277b] - w[r278f] + w[r278b] - w[r279f] + w[r279b];

  cdot[sCH2CHO] = w[r285f] - w[r285b] + w[r294f] - w[r294b] + w[r296f] -
                  w[r296b] + w[r299f] - w[r299b] + w[r304f] - w[r304b] -
                  w[r305] - w[r306] - w[r307] - w[r308f] + w[r308b] - w[r309f] +
                  w[r309b] - w[r310f] + w[r310b] - w[r311f] + w[r311b];

  cdot[sCH3CHO] = w[r286f] - w[r286b] - w[r296f] + w[r296b] - w[r297] -
                  w[r298] - w[r299f] + w[r299b] - w[r300] - w[r301] - w[r302] -
                  w[r303];

  cdot[sC3H8] = w[r312f] - w[r312b] - w[r313f] + w[r313b] - w[r314f] +
                w[r314b] - w[r315f] + w[r315b] + w[r316] - w[r317f] + w[r317b] +
                w[r320f] - w[r320b] + w[r323f] - w[r323b];

  cdot[sC3H7] = w[r313f] - w[r313b] + w[r314f] - w[r314b] + w[r315f] -
                w[r315b] - w[r316] + w[r317f] - w[r317b] + w[r318f] - w[r318b] -
                w[r319f] + w[r319b] - w[r320f] + w[r320b] - w[r321f] +
                w[r321b] - w[r322f] + w[r322b] - w[r323f] + w[r323b] - w[r324] -
                w[r325f] + w[r325b];
}

double GetLindRateCoeff(double temp, double pressure, double k0, double kInf,
                        double fc, double conc) {
  const double R = 8314.34; /* [J / kmole K] */
  double Ntmp;
  double kl;
  double f;
  double cCoeff, dCoeff, log10kNull;

  int iTroe = 1;

  if (isinf(conc)) {
    conc = pressure / (R * temp);
  }
  Ntmp = 0.75 - 1.27 * (fc < 1e-200 ? -200 : log10(fc));
  if (iTroe) {
    cCoeff = -0.4 - 0.67 * (fc < 1e-200 ? -200 : log10(fc));
    dCoeff = 0.14;
    k0 *= conc / MAX_C(kInf, 1.0e-60);
    log10kNull = k0 < 1e-200 ? -200 : log10(k0);
    f = (log10kNull + cCoeff) / (Ntmp - dCoeff * (log10kNull + cCoeff));
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

void Gri30::ComputeThermoData(span<double> h, span<double> cp, double T,
                              span<double> s) const {
  /*
          This function computes entropy 's', enthalpy 'h' and heat
          capacity 'cp' as function of temperature 'T' for all non steady
          state species in units [K/kg K], [J/kg] and [J/kg K], respectively.
          The parameter s, h and cp should provide workspace of length 53
          If you compile using cdecl, 's' is optional for backwards
     compatibility: This function will only store entropies if *s == 0, to
     detect omitted arguments (in most cases). On the other hand, if *s == 0
     after invocation, you were using an old mechanism without support for
     entropies. */

  int i;
  if (T > 1000.0) {
    h[sAR] = 0x1.a042152c51ec9p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a042152c51ec9p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sAR] = 0x1.a042152c51ec9p+7 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.176c8b439581p+2);
    h[sO] = 0x1.03d2b851eb852p+9 *
            (T * (0x1.48e2c7b5e1a3dp+1 +
                  T * (-0x1.689a006cdb509p-15 +
                       T * (0x1.e0722ea52787bp-27 +
                            T * (-0x1.607823c010ceap-39 +
                                 T * 0x1.1b3c369230ddbp-52)))) +
             0x1.c88650ff97247p+14);
    cp[sO] =
        0x1.03d2b851eb852p+9 *
        (0x1.48e2c7b5e1a3dp+1 +
         T * (-0x1.689a006cdb509p-14 +
              T * (0x1.6855a2fbdda5cp-25 +
                   T * (-0x1.607823c010ceap-37 + T * 0x1.620b4436bd152p-50))));
    s[sO] =
        0x1.03d2b851eb852p+9 *
        (0x1.48e2c7b5e1a3dp+1 * log(T) +
         T * (-0x1.689a006cdb509p-14 +
              T * (0x1.6855a2fbdda5cp-26 +
                   T * (-0x1.d5f585001668dp-39 + T * 0x1.620b4436bd152p-52))) +
         0x1.32329ab1f280ep+2);
    h[sO2] = 0x1.03d2b851eb852p+8 *
             (T * (0x1.a42a332f575b8p+1 +
                   T * (0x1.84c851ecf574dp-11 +
                        T * (-0x1.0f496e356c6d3p-22 +
                             T * (0x1.cca1706e20348p-35 +
                                  T * -0x1.3852c0eebb014p-48)))) -
              0x1.101d4b48d3ae7p+10);
    cp[sO2] =
        0x1.03d2b851eb852p+8 *
        (0x1.a42a332f575b8p+1 +
         T * (0x1.84c851ecf574dp-10 +
              T * (-0x1.96ee255022a3dp-21 +
                   T * (0x1.cca1706e20348p-33 + T * -0x1.8667712a69c18p-46))));
    s[sO2] =
        0x1.03d2b851eb852p+8 *
        (0x1.a42a332f575b8p+1 * log(T) +
         T * (0x1.84c851ecf574dp-10 +
              T * (-0x1.96ee255022a3dp-22 +
                   T * (0x1.33164af415785p-34 + T * -0x1.8667712a69c18p-48))) +
         0x1.5d01bdd004baap+2);
    h[sH] = 0x1.01c2d34d34d35p+13 *
            (T * (0x1.40000015798eep+1 +
                  T * (-0x1.961a6ec6a1e45p-37 +
                       T * (0x1.840f10457b303p-48 +
                            T * (-0x1.5d647fbe2e672p-60 +
                                 T * 0x1.e1d3b0ee8129dp-74)))) +
             0x1.8e06a3bcd35a8p+14);
    cp[sH] =
        0x1.01c2d34d34d35p+13 *
        (0x1.40000015798eep+1 +
         T * (-0x1.961a6ec6a1e45p-36 +
              T * (0x1.230b4c341c642p-46 +
                   T * (-0x1.5d647fbe2e672p-58 + T * 0x1.2d244e9510ba2p-71))));
    s[sH] =
        0x1.01c2d34d34d35p+13 *
        (0x1.40000015798eep+1 * log(T) +
         T * (-0x1.961a6ec6a1e45p-36 +
              T * (0x1.230b4c341c642p-47 +
                   T * (-0x1.d1db54fd93343p-60 + T * 0x1.2d244e9510ba2p-73))) -
         0x1.c9673eed3f77dp-2);
    h[sOH] = 0x1.e8d94973d693ep+8 *
             (T * (0x1.8be3be406d029p+1 +
                   T * (0x1.1f88fd8e1440ep-12 +
                        T * (0x1.6a395011eed9bp-25 +
                             T * (-0x1.82ca91827ab6dp-36 +
                                  T * 0x1.526b0ac38c4c7p-49)))) +
              0x1.e2550624dd2f2p+11);
    cp[sOH] =
        0x1.e8d94973d693ep+8 *
        (0x1.8be3be406d029p+1 +
         T * (0x1.1f88fd8e1440ep-11 +
              T * (0x1.0faafc0d73234p-23 +
                   T * (-0x1.82ca91827ab6dp-34 + T * 0x1.a705cd746f5f8p-47))));
    s[sOH] =
        0x1.e8d94973d693ep+8 *
        (0x1.8be3be406d029p+1 * log(T) +
         T * (0x1.1f88fd8e1440ep-11 +
              T * (0x1.0faafc0d73234p-24 +
                   T * (-0x1.01dc6101a7249p-35 + T * 0x1.a705cd746f5f8p-49))) +
         0x1.1e82305be85e2p+2);
    h[sH2] = 0x1.01c2d34d34d35p+12 *
             (T * (0x1.ab2bf6fecf7e5p+1 +
                   T * (-0x1.9e6b00cea2ba4p-16 +
                        T * (0x1.65866c28c4ed8p-23 +
                             T * (-0x1.8adee4a4b6baep-35 +
                                  T * 0x1.2099318342be3p-48)))) -
              0x1.db14578e5c4ebp+9);
    cp[sH2] =
        0x1.01c2d34d34d35p+12 *
        (0x1.ab2bf6fecf7e5p+1 +
         T * (-0x1.9e6b00cea2ba4p-15 +
              T * (0x1.0c24d11e93b22p-21 +
                   T * (-0x1.8adee4a4b6baep-33 + T * 0x1.68bf7de4136dcp-46))));
    s[sH2] =
        0x1.01c2d34d34d35p+12 *
        (0x1.ab2bf6fecf7e5p+1 * log(T) +
         T * (-0x1.9e6b00cea2ba4p-15 +
              T * (0x1.0c24d11e93b22p-22 +
                   T * (-0x1.073f431879d1fp-34 + T * 0x1.68bf7de4136dcp-48))) -
         0x1.9a3e342daf0fdp+1);
    h[sHO2] = 0x1.f7c6fae98a1dp+7 *
              (T * (0x1.0119fbbf289f6p+2 +
                    T * (0x1.2593e46a1f9dfp-10 +
                         T * (-0x1.c597158096605p-23 +
                              T * (0x1.f675fa32f618cp-36 +
                                   T * -0x1.370671f4bb474p-49)))) +
               0x1.bf6d462c343b7p+6);
    cp[sHO2] =
        0x1.f7c6fae98a1dp+7 *
        (0x1.0119fbbf289f6p+2 +
         T * (0x1.2593e46a1f9dfp-9 +
              T * (-0x1.5431502070c84p-21 +
                   T * (0x1.f675fa32f618cp-34 + T * -0x1.84c80e71ea191p-47))));
    s[sHO2] =
        0x1.f7c6fae98a1dp+7 *
        (0x1.0119fbbf289f6p+2 * log(T) +
         T * (0x1.2593e46a1f9dfp-9 +
              T * (-0x1.5431502070c84p-22 +
                   T * (0x1.4ef951774ebb3p-35 + T * -0x1.84c80e71ea191p-49))) +
         0x1.e47e3a2d2278p+1);
    h[sH2O2] = 0x1.e8d94973d693ep+7 *
               (T * (0x1.0a8f681d1fcb7p+2 +
                     T * (0x1.41abe4bc5704bp-9 +
                          T * (-0x1.544474230e883p-21 +
                               T * (0x1.981f91177d524p-34 +
                                    T * -0x1.9eeb6a97eb6fcp-48)))) -
                0x1.1717269ad42c4p+14);
    cp[sH2O2] =
        0x1.e8d94973d693ep+7 *
        (0x1.0a8f681d1fcb7p+2 +
         T * (0x1.41abe4bc5704bp-8 +
              T * (-0x1.fe66ae3495cc5p-20 +
                   T * (0x1.981f91177d524p-32 + T * -0x1.0353229ef325dp-45))));
    s[sH2O2] =
        0x1.e8d94973d693ep+7 *
        (0x1.0a8f681d1fcb7p+2 * log(T) +
         T * (0x1.41abe4bc5704bp-8 +
              T * (-0x1.fe66ae3495cc5p-21 +
                   T * (0x1.10150b64fe36dp-33 + T * -0x1.0353229ef325dp-47))) +
         0x1.75449ec074fabp+1);
    h[sCH] = 0x1.3f5713b4542aap+9 *
             (T * (0x1.70718843050d6p+1 +
                   T * (0x1.fd09d40e9c422p-12 +
                        T * (0x1.9d97c8eb05901p-25 +
                             T * (-0x1.1f62b7f017dbdp-35 +
                                  T * 0x1.fb83a7106ba7cp-49)))) +
              0x1.15646fb7e91p+16);
    cp[sCH] =
        0x1.3f5713b4542aap+9 *
        (0x1.70718843050d6p+1 +
         T * (0x1.fd09d40e9c422p-11 +
              T * (0x1.3631d6b0442c1p-23 +
                   T * (-0x1.1f62b7f017dbdp-33 + T * 0x1.3d32486a4348dp-46))));
    s[sCH] =
        0x1.3f5713b4542aap+9 *
        (0x1.70718843050d6p+1 * log(T) +
         T * (0x1.fd09d40e9c422p-11 +
              T * (0x1.3631d6b0442c1p-24 +
                   T * (-0x1.7f2e4a9575251p-35 + T * 0x1.3d32486a4348dp-48))) +
         0x1.5f09e98310ec1p+2);
    h[sCO] = 0x1.28d5af05f0d41p+8 *
             (T * (0x1.5b8b33bac2892p+1 +
                   T * (0x1.0e56efb7c8adbp-10 +
                        T * (-0x1.657e610d71e2bp-22 +
                             T * (0x1.f9e45483ed15bp-35 +
                                  T * -0x1.257cbef78f17p-48)))) -
              0x1.ba3efaacd9e84p+13);
    cp[sCO] =
        0x1.28d5af05f0d41p+8 *
        (0x1.5b8b33bac2892p+1 +
         T * (0x1.0e56efb7c8adbp-9 +
              T * (-0x1.0c1ec8ca156ap-20 +
                   T * (0x1.f9e45483ed15bp-33 + T * -0x1.6edbeeb572dccp-46))));
    s[sCO] =
        0x1.28d5af05f0d41p+8 *
        (0x1.5b8b33bac2892p+1 * log(T) +
         T * (0x1.0e56efb7c8adbp-9 +
              T * (-0x1.0c1ec8ca156ap-21 +
                   T * (0x1.5142e3029e0e7p-34 + T * -0x1.6edbeeb572dccp-48))) +
         0x1.f465612dc25bp+2);
    h[s3XCH2] = 0x1.2863e91362288p+9 *
                (T * (0x1.6fe28bbb5f921p+1 +
                      T * (0x1.df4030068403dp-10 +
                           T * (-0x1.f8480a2ca2df1p-22 +
                                T * (0x1.1e1208519dbd9p-34 +
                                     T * -0x1.0e8b3f67fea4fp-48)))) +
                 0x1.696f353f7ced9p+15);
    cp[s3XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.6fe28bbb5f921p+1 +
         T * (0x1.df4030068403dp-9 +
              T * (-0x1.7a3607a17a275p-20 +
                   T * (0x1.1e1208519dbd9p-32 + T * -0x1.522e0f41fe4e3p-46))));
    s[s3XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.6fe28bbb5f921p+1 * log(T) +
         T * (0x1.df4030068403dp-9 +
              T * (-0x1.7a3607a17a275p-21 +
                   T * (0x1.7d6d606cd2521p-34 + T * -0x1.522e0f41fe4e3p-48))) +
         0x1.8af4d47dc6297p+2);
    h[sHCO] = 0x1.1e86068742b68p+8 *
              (T * (0x1.62d69c2e745cap+1 +
                    T * (0x1.44dbe8babc438p-9 +
                         T * (-0x1.bc9c59a38633fp-21 +
                              T * (0x1.43e5261c638b2p-33 +
                                   T * -0x1.806efc8d2c14ep-47)))) +
               0x1.f57d617c1bda5p+11);
    cp[sHCO] =
        0x1.1e86068742b68p+8 *
        (0x1.62d69c2e745cap+1 +
         T * (0x1.44dbe8babc438p-8 +
              T * (-0x1.4d75433aa4a6fp-19 +
                   T * (0x1.43e5261c638b2p-31 + T * -0x1.e08abbb0771a1p-45))));
    s[sHCO] =
        0x1.1e86068742b68p+8 *
        (0x1.62d69c2e745cap+1 * log(T) +
         T * (0x1.44dbe8babc438p-8 +
              T * (-0x1.4d75433aa4a6fp-20 +
                   T * (0x1.afdc32d084b98p-33 + T * -0x1.e08abbb0771a1p-47))) +
         0x1.398c0aa54a7cdp+3);
    h[s1XCH2] = 0x1.2863e91362288p+9 *
                (T * (0x1.256183d389aa6p+1 +
                      T * (0x1.3120cfb16b35fp-9 +
                           T * (-0x1.680c0914ed793p-21 +
                                T * (0x1.cb7e14e4de1f4p-34 +
                                     T * -0x1.e9953783e1eedp-48)))) +
                 0x1.8ddbffd8adabap+15);
    cp[s1XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.256183d389aa6p+1 +
         T * (0x1.3120cfb16b35fp-8 +
              T * (-0x1.0e0906cfb21aep-19 +
                   T * (0x1.cb7e14e4de1f4p-32 + T * -0x1.31fd42b26d354p-45))));
    s[s1XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.256183d389aa6p+1 * log(T) +
         T * (0x1.3120cfb16b35fp-8 +
              T * (-0x1.0e0906cfb21aep-20 +
                   T * (0x1.32540dede96a3p-33 + T * -0x1.31fd42b26d354p-47))) +
         0x1.140c4d45ae0a1p+3);
    h[sCH3] = 0x1.1484949ef95bfp+9 *
              (T * (0x1.249265f3a4d2ep+1 +
                    T * (0x1.da795f50c1384p-9 +
                         T * (-0x1.0b48fa9d24b6p-20 +
                              T * (0x1.477b292040de5p-33 +
                                   T * -0x1.509ec68687becp-47)))) +
               0x1.061e5652bd3c3p+14);
    cp[sCH3] =
        0x1.1484949ef95bfp+9 *
        (0x1.249265f3a4d2ep+1 +
         T * (0x1.da795f50c1384p-8 +
              T * (-0x1.90ed77ebb711p-19 +
                   T * (0x1.477b292040de5p-31 + T * -0x1.a4c6782829ae7p-45))));
    s[sCH3] =
        0x1.1484949ef95bfp+9 *
        (0x1.249265f3a4d2ep+1 * log(T) +
         T * (0x1.da795f50c1384p-8 +
              T * (-0x1.90ed77ebb711p-20 +
                   T * (0x1.b4a436d5abd31p-33 + T * -0x1.a4c6782829ae7p-47))) +
         0x1.0f5cbf83b907ap+3);
    h[sCH2O] = 0x1.14e79947885d3p+8 *
               (T * (0x1.c2bc95c7fddb1p+0 +
                     T * (0x1.2d77335291c21p-8 +
                          T * (-0x1.8bb9fb0f6d53cp-20 +
                               T * (0x1.14a3f4fe94722p-32 +
                                    T * -0x1.3e714b400a41ap-46)))) -
                0x1.b55ea88ce703bp+13);
    cp[sCH2O] =
        0x1.14e79947885d3p+8 *
        (0x1.c2bc95c7fddb1p+0 +
         T * (0x1.2d77335291c21p-7 +
              T * (-0x1.28cb7c4b91fedp-18 +
                   T * (0x1.14a3f4fe94722p-30 + T * -0x1.8e0d9e100cd2p-44))));
    s[sCH2O] =
        0x1.14e79947885d3p+8 *
        (0x1.c2bc95c7fddb1p+0 * log(T) +
         T * (0x1.2d77335291c21p-7 +
              T * (-0x1.28cb7c4b91fedp-19 +
                   T * (0x1.70da9bfe1b42dp-32 + T * -0x1.8e0d9e100cd2p-46))) +
         0x1.b5009917939a8p+3);
    h[sCH4] = 0x1.03249373f39ddp+9 *
              (T * (0x1.32977b314eac9p-4 +
                    T * (0x1.b6cb6711ca116p-8 +
                         T * (-0x1.007bd5ae50b8fp-19 +
                              T * (0x1.5027b89a2ce85p-32 +
                                   T * -0x1.6ed3f97ad9a02p-46)))) -
               0x1.27e2c1b866e44p+13);
    cp[sCH4] =
        0x1.03249373f39ddp+9 *
        (0x1.32977b314eac9p-4 +
         T * (0x1.b6cb6711ca116p-7 +
              T * (-0x1.80b9c08579157p-18 +
                   T * (0x1.5027b89a2ce85p-30 + T * -0x1.ca88f7d990082p-44))));
    s[sCH4] =
        0x1.03249373f39ddp+9 *
        (0x1.32977b314eac9p-4 * log(T) +
         T * (0x1.b6cb6711ca116p-7 +
              T * (-0x1.80b9c08579157p-19 +
                   T * (0x1.c034f622e68b1p-32 + T * -0x1.ca88f7d990082p-46))) +
         0x1.26ff4128bf3bfp+4);
    h[sCO2] = 0x1.79d6b3468efa8p+7 *
              (T * (0x1.edc1423f95973p+1 +
                    T * (0x1.214cd7e62db66p-9 +
                         T * (-0x1.8c5b3dc37e1cfp-21 +
                              T * (0x1.1fcab1b9912a3p-33 +
                                   T * -0x1.542c2854e86a4p-47)))) -
               0x1.7cee54fdf3b64p+15);
    cp[sCO2] =
        0x1.79d6b3468efa8p+7 *
        (0x1.edc1423f95973p+1 +
         T * (0x1.214cd7e62db66p-8 +
              T * (-0x1.29446e529e95bp-19 +
                   T * (0x1.1fcab1b9912a3p-31 + T * -0x1.a937326a2284cp-45))));
    s[sCO2] =
        0x1.79d6b3468efa8p+7 *
        (0x1.edc1423f95973p+1 * log(T) +
         T * (0x1.214cd7e62db66p-8 +
              T * (-0x1.29446e529e95bp-20 +
                   T * (0x1.7fb8ecf76c384p-33 + T * -0x1.a937326a2284cp-47))) +
         0x1.22c509340641ep+1);
    h[sCH2OH] = 0x1.0be9223bc1817p+8 *
                (T * (0x1.d8a944f2ce3e4p+1 +
                      T * (0x1.1b4df55708be3p-8 +
                           T * (-0x1.4fa2824b870b8p-20 +
                                T * (0x1.b0c96ce8203e2p-33 +
                                     T * -0x1.d35522361598ap-47)))) -
                 0x1.9550335d249e4p+11);
    cp[sCH2OH] =
        0x1.0be9223bc1817p+8 *
        (0x1.d8a944f2ce3e4p+1 +
         T * (0x1.1b4df55708be3p-7 +
              T * (-0x1.f773c3714a914p-19 +
                   T * (0x1.b0c96ce8203e2p-31 + T * -0x1.24153561cd7f6p-44))));
    s[sCH2OH] =
        0x1.0be9223bc1817p+8 *
        (0x1.d8a944f2ce3e4p+1 * log(T) +
         T * (0x1.1b4df55708be3p-7 +
              T * (-0x1.f773c3714a914p-20 +
                   T * (0x1.2086489ac0297p-32 + T * -0x1.24153561cd7f6p-46))) +
         0x1.73de1ecef8203p+2);
    h[sCH3O] = 0x1.0be9223bc1817p+8 *
               (T * (0x1.e2a98aa8650e7p+1 +
                     T * (0x1.01eee717c07fdp-8 +
                          T * (-0x1.db60e105b3b4p-21 +
                               T * (0x1.b1b1dcc557332p-34 +
                                    T * -0x1.3075c5faa9763p-48)))) +
                0x1.ff54801f75105p+6);
    cp[sCH3O] =
        0x1.0be9223bc1817p+8 *
        (0x1.e2a98aa8650e7p+1 +
         T * (0x1.01eee717c07fdp-7 +
              T * (-0x1.6488a8c446c7p-19 +
                   T * (0x1.b1b1dcc557332p-32 + T * -0x1.7c93377953d3bp-46))));
    s[sCH3O] =
        0x1.0be9223bc1817p+8 *
        (0x1.e2a98aa8650e7p+1 * log(T) +
         T * (0x1.01eee717c07fdp-7 +
              T * (-0x1.6488a8c446c7p-20 +
                   T * (0x1.21213dd8e4cccp-33 + T * -0x1.7c93377953d3bp-48))) +
         0x1.76fc504816fp+1);
    h[sCH3OH] = 0x1.037b88ab2ad7p+8 *
                (T * (0x1.ca2a4c2ed7aedp+0 +
                      T * (0x1.cdd39bbea4fb5p-8 +
                           T * (-0x1.1cc401d83630dp-19 +
                                T * (0x1.7bcd41e07455p-32 +
                                     T * -0x1.a5c0fef86370dp-46)))) -
                 0x1.8c7b7fb15b574p+14);
    cp[sCH3OH] =
        0x1.037b88ab2ad7p+8 *
        (0x1.ca2a4c2ed7aedp+0 +
         T * (0x1.cdd39bbea4fb5p-7 +
              T * (-0x1.ab2602c451494p-18 +
                   T * (0x1.7bcd41e07455p-30 + T * -0x1.07989f5b3e268p-43))));
    s[sCH3OH] =
        0x1.037b88ab2ad7p+8 *
        (0x1.ca2a4c2ed7aedp+0 * log(T) +
         T * (0x1.cdd39bbea4fb5p-7 +
              T * (-0x1.ab2602c451494p-19 +
                   T * (0x1.fa6702809b1cp-32 + T * -0x1.07989f5b3e268p-45))) +
         0x1.d0135a1a27c97p+3);
    h[sC2H] = 0x1.4c3397c02c83cp+8 *
              (T * (0x1.957aaf1dba502p+1 +
                    T * (0x1.377101463a6fcp-9 +
                         T * (-0x1.48e6585748e57p-21 +
                              T * (0x1.4e75f1b05c789p-34 +
                                   T * -0x1.fed6b3b564915p-49)))) +
               0x1.063110a3d70a4p+16);
    cp[sC2H] =
        0x1.4c3397c02c83cp+8 *
        (0x1.957aaf1dba502p+1 +
         T * (0x1.377101463a6fcp-8 +
              T * (-0x1.ed598482ed583p-20 +
                   T * (0x1.4e75f1b05c789p-32 + T * -0x1.3f4630515edadp-46))));
    s[sC2H] =
        0x1.4c3397c02c83cp+8 *
        (0x1.957aaf1dba502p+1 * log(T) +
         T * (0x1.377101463a6fcp-8 +
              T * (-0x1.ed598482ed583p-21 +
                   T * (0x1.bdf29795d0a0cp-34 + T * -0x1.3f4630515edadp-48))) +
         0x1.a8b27fe4bcadap+2);
    h[sC2H2] = 0x1.3f5713b4542aap+8 *
               (T * (0x1.0971c7ee6baep+2 +
                     T * (0x1.86b42b3f9ab1ep-9 +
                          T * (-0x1.a8a7da8ef5fe7p-21 +
                               T * (0x1.00f66a3bafde2p-33 +
                                    T * -0x1.044c229f7f34fp-47)))) +
                0x1.953fff2e48e8ap+14);
    cp[sC2H2] =
        0x1.3f5713b4542aap+8 *
        (0x1.0971c7ee6baep+2 +
         T * (0x1.86b42b3f9ab1ep-8 +
              T * (-0x1.3e7de3eb387edp-19 +
                   T * (0x1.00f66a3bafde2p-31 + T * -0x1.455f2b475f023p-45))));
    s[sC2H2] =
        0x1.3f5713b4542aap+8 *
        (0x1.0971c7ee6baep+2 * log(T) +
         T * (0x1.86b42b3f9ab1ep-8 +
              T * (-0x1.3e7de3eb387edp-20 +
                   T * (0x1.569de2fa3fd2dp-33 + T * -0x1.455f2b475f023p-47))) -
         0x1.3af3b599d553bp+0);
    h[sHCCO] = 0x1.954cff46b5264p+7 *
               (T * (0x1.683486198a14cp+2 +
                     T * (0x1.0bbca21f5e9bfp-9 +
                          T * (-0x1.1d28ea5b6f3c1p-21 +
                               T * (0x1.3abf2c56d8d5fp-34 +
                                    T * -0x1.17b243118233fp-48)))) +
                0x1.2dfcdc28f5c29p+14);
    cp[sHCCO] =
        0x1.954cff46b5264p+7 *
        (0x1.683486198a14cp+2 +
         T * (0x1.0bbca21f5e9bfp-8 +
              T * (-0x1.abbd5f8926da2p-20 +
                   T * (0x1.3abf2c56d8d5fp-32 + T * -0x1.5d9ed3d5e2c0fp-46))));
    s[sHCCO] =
        0x1.954cff46b5264p+7 *
        (0x1.683486198a14cp+2 * log(T) +
         T * (0x1.0bbca21f5e9bfp-8 +
              T * (-0x1.abbd5f8926da2p-21 +
                   T * (0x1.a3a99073cbc7fp-34 + T * -0x1.5d9ed3d5e2c0fp-48))) -
         0x1.f712be48a58b4p+1);
    h[sC2H3] = 0x1.3370009b17838p+8 *
               (T * (0x1.8224031487768p+1 +
                     T * (0x1.52803e497ede9p-8 +
                          T * (-0x1.a2d53f39c89fcp-20 +
                               T * (0x1.17b98c3c9686ep-32 +
                                    T * -0x1.36c974e43b3b7p-46)))) +
                0x1.0e69bf6fd21ffp+15);
    cp[sC2H3] =
        0x1.3370009b17838p+8 *
        (0x1.8224031487768p+1 +
         T * (0x1.52803e497ede9p-7 +
              T * (-0x1.3a1fef6b5677dp-18 +
                   T * (0x1.17b98c3c9686ep-30 + T * -0x1.847bd21d4a0a4p-44))));
    s[sC2H3] =
        0x1.3370009b17838p+8 *
        (0x1.8224031487768p+1 * log(T) +
         T * (0x1.52803e497ede9p-7 +
              T * (-0x1.3a1fef6b5677dp-19 +
                   T * (0x1.74f765a61e093p-32 + T * -0x1.847bd21d4a0a4p-46))) +
         0x1.f26383479da37p+2);
    h[sCH2CO] = 0x1.8b94f63b4adbcp+7 *
                (T * (0x1.20b91864fbad3p+2 +
                      T * (0x1.2707a64c0b6f9p-8 +
                           T * (-0x1.75123ec2bd6e9p-20 +
                                T * (0x1.fb9d615c62415p-33 +
                                     T * -0x1.1e5ee26606538p-46)))) -
                 0x1.d7f0d989df117p+12);
    cp[sCH2CO] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.20b91864fbad3p+2 +
         T * (0x1.2707a64c0b6f9p-7 +
              T * (-0x1.17cdaf120e12fp-18 +
                   T * (0x1.fb9d615c62415p-31 + T * -0x1.65f69aff87e86p-44))));
    s[sCH2CO] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.20b91864fbad3p+2 * log(T) +
         T * (0x1.2707a64c0b6f9p-7 +
              T * (-0x1.17cdaf120e12fp-19 +
                   T * (0x1.5268eb92ec2b9p-32 + T * -0x1.65f69aff87e86p-46))) +
         0x1.43b5e7d8ecc0bp-1);
    h[sC2H4] = 0x1.2863e91362288p+8 *
               (T * (0x1.049f4a5d9c3d6p+1 +
                     T * (0x1.dfe6a5720731ep-8 +
                          T * (-0x1.2c3c348d37b5bp-19 +
                               T * (0x1.94aeec0be93f2p-32 +
                                    T * -0x1.c4e760753f615p-46)))) +
                0x1.34be2da122fadp+12);
    cp[sC2H4] =
        0x1.2863e91362288p+8 *
        (0x1.049f4a5d9c3d6p+1 +
         T * (0x1.dfe6a5720731ep-7 +
              T * (-0x1.c25a4ed3d3908p-18 +
                   T * (0x1.94aeec0be93f2p-30 + T * -0x1.1b109c49479cdp-43))));
    s[sC2H4] =
        0x1.2863e91362288p+8 *
        (0x1.049f4a5d9c3d6p+1 * log(T) +
         T * (0x1.dfe6a5720731ep-7 +
              T * (-0x1.c25a4ed3d3908p-19 +
                   T * (0x1.0dc9f2b29b7f7p-31 + T * -0x1.1b109c49479cdp-45))) +
         0x1.49c595d6967a3p+3);
    h[sC2H5] = 0x1.1e1c03861415ep+8 *
               (T * (0x1.f4645cf6d1024p+0 +
                     T * (0x1.1d0972c8defc6p-7 +
                          T * (-0x1.651c9302484bfp-19 +
                               T * (0x1.e1a27cc16cb18p-32 +
                                    T * -0x1.0d91ff12db09fp-45)))) +
                0x1.91cc28f5c28f6p+13);
    cp[sC2H5] =
        0x1.1e1c03861415ep+8 *
        (0x1.f4645cf6d1024p+0 +
         T * (0x1.1d0972c8defc6p-6 +
              T * (-0x1.0bd56e41b638fp-17 +
                   T * (0x1.e1a27cc16cb18p-30 + T * -0x1.50f67ed791cc7p-43))));
    s[sC2H5] =
        0x1.1e1c03861415ep+8 *
        (0x1.f4645cf6d1024p+0 * log(T) +
         T * (0x1.1d0972c8defc6p-6 +
              T * (-0x1.0bd56e41b638fp-18 +
                   T * (0x1.4116fdd648765p-31 + T * -0x1.50f67ed791cc7p-45))) +
         0x1.aecc4304618e9p+3);
    h[sC2H6] = 0x1.1484949ef95bfp+8 *
               (T * (0x1.1266d373affbp+0 +
                     T * (0x1.634a9ae4e575dp-7 +
                          T * (-0x1.c089bdb6616a1p-19 +
                               T * (0x1.304e6c1a9e602p-31 +
                                    T * -0x1.56475dfc0f594p-45)))) -
                0x1.651325460aa65p+13);
    cp[sC2H6] =
        0x1.1484949ef95bfp+8 *
        (0x1.1266d373affbp+0 +
         T * (0x1.634a9ae4e575dp-6 +
              T * (-0x1.50674e48c90f9p-17 +
                   T * (0x1.304e6c1a9e602p-29 + T * -0x1.abd9357b132f9p-43))));
    s[sC2H6] =
        0x1.1484949ef95bfp+8 *
        (0x1.1266d373affbp+0 * log(T) +
         T * (0x1.634a9ae4e575dp-6 +
              T * (-0x1.50674e48c90f9p-18 +
                   T * (0x1.95bde578d32adp-31 + T * -0x1.abd9357b132f9p-45))) +
         0x1.e3b31535f22a5p+3);
    h[sH2O] = 0x1.cd7f5ff1730a8p+8 *
              (T * (0x1.8459ddac6e07ap+1 +
                    T * (0x1.1d553f9364039p-10 +
                         T * (-0x1.d5ca6d748cbc5p-25 +
                              T * (-0x1.aacb906a5f7ccp-36 +
                                   T * 0x1.e4ce6d148b33p-49)))) -
               0x1.d4d1303afb7e9p+14);
    cp[sH2O] =
        0x1.cd7f5ff1730a8p+8 *
        (0x1.8459ddac6e07ap+1 +
         T * (0x1.1d553f9364039p-9 +
              T * (-0x1.6057d217698d4p-23 +
                   T * (-0x1.aacb906a5f7ccp-34 + T * 0x1.2f01042cd6ffep-46))));
    s[sH2O] =
        0x1.cd7f5ff1730a8p+8 *
        (0x1.8459ddac6e07ap+1 * log(T) +
         T * (0x1.1d553f9364039p-9 +
              T * (-0x1.6057d217698d4p-24 +
                   T * (-0x1.1c87b59c3fa88p-35 + T * 0x1.2f01042cd6ffep-48))) +
         0x1.3ddf8fb2900aap+2);
    h[sN2] = 0x1.28ba905aa1e55p+8 *
             (T * (0x1.769c23b7952d2p+1 +
                   T * (0x1.86106ec5d7fc6p-11 +
                        T * (-0x1.96ee5424b7af7p-23 +
                             T * (0x1.bc128a9b8a724p-36 +
                                  T * -0x1.854ddeba4197bp-50)))) -
              0x1.cd661b089a027p+9);
    cp[sN2] =
        0x1.28ba905aa1e55p+8 *
        (0x1.769c23b7952d2p+1 +
         T * (0x1.86106ec5d7fc6p-10 +
              T * (-0x1.3132bf1b89c39p-21 +
                   T * (0x1.bc128a9b8a724p-34 + T * -0x1.e6a15668d1fd9p-48))));
    s[sN2] =
        0x1.28ba905aa1e55p+8 *
        (0x1.769c23b7952d2p+1 * log(T) +
         T * (0x1.86106ec5d7fc6p-10 +
              T * (-0x1.3132bf1b89c39p-22 +
                   T * (0x1.280c5c67b1a18p-35 + T * -0x1.e6a15668d1fd9p-50))) +
         0x1.7ec0f88333fc8p+2);
    h[sC] = 0x1.5a24731aa4fb9p+9 *
            (T * (0x1.3f0fc61badb08p+1 +
                  T * (0x1.928f7063e132p-16 +
                       T * (-0x1.9ecca30d2609fp-26 +
                            T * (0x1.493ad8d0d6b66p-37 +
                                 T * -0x1.18e5865166a27p-50)))) +
             0x1.4dcb4b98c7e28p+16);
    cp[sC] =
        0x1.5a24731aa4fb9p+9 *
        (0x1.3f0fc61badb08p+1 +
         T * (0x1.928f7063e132p-15 +
              T * (-0x1.37197a49dc877p-24 +
                   T * (0x1.493ad8d0d6b66p-35 + T * -0x1.5f1ee7e5c04b1p-48))));
    s[sC] =
        0x1.5a24731aa4fb9p+9 *
        (0x1.3f0fc61badb08p+1 * log(T) +
         T * (0x1.928f7063e132p-15 +
              T * (-0x1.37197a49dc877p-25 +
                   T * (0x1.b6f92116739ddp-37 + T * -0x1.5f1ee7e5c04b1p-50))) +
         0x1.334bd64cfe358p+2);
    h[sHCCOH] = 0x1.8b94f63b4adbcp+7 *
                (T * (0x1.7b200416e5f59p+2 +
                      T * (0x1.bd24e4100a643p-9 +
                           T * (-0x1.cb2d8a19f367cp-21 +
                                T * (0x1.eea583d58eb03p-34 +
                                     T * -0x1.af7b79e42769ap-48)))) +
                 0x1.c60a04189374cp+12);
    cp[sHCCOH] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.7b200416e5f59p+2 +
         T * (0x1.bd24e4100a643p-8 +
              T * (-0x1.58622793768ddp-19 +
                   T * (0x1.eea583d58eb03p-32 + T * -0x1.0dad2c2e98a2p-45))));
    s[sHCCOH] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.7b200416e5f59p+2 * log(T) +
         T * (0x1.bd24e4100a643p-8 +
              T * (-0x1.58622793768ddp-20 +
                   T * (0x1.49c3ad3909cadp-33 + T * -0x1.0dad2c2e98a2p-47))) -
         0x1.e68377ef24e67p+2);
    h[sN] = 0x1.28ba905aa1e55p+9 *
            (T * (0x1.353d9df0406f6p+1 +
                  T * (0x1.6ec5b3e8752afp-14 +
                       T * (-0x1.54cd49949c5cap-25 +
                            T * (0x1.09df73ee980cbp-37 +
                                 T * -0x1.d57e0610cd838p-52)))) +
             0x1.b68b8bc6a7efap+15);
    cp[sN] =
        0x1.28ba905aa1e55p+9 *
        (0x1.353d9df0406f6p+1 +
         T * (0x1.6ec5b3e8752afp-13 +
              T * (-0x1.ff33ee5eea8afp-24 +
                   T * (0x1.09df73ee980cbp-35 + T * -0x1.256ec3ca80723p-49))));
    s[sN] =
        0x1.28ba905aa1e55p+9 *
        (0x1.353d9df0406f6p+1 * log(T) +
         T * (0x1.6ec5b3e8752afp-13 +
              T * (-0x1.ff33ee5eea8afp-25 +
                   T * (0x1.627f453e2010fp-37 + T * -0x1.256ec3ca80723p-51))) +
         0x1.29933424cabaap+2);
    h[sNO] = 0x1.150d649310ce2p+8 *
             (T * (0x1.a15b863893c54p+1 +
                   T * (0x1.383da80ef9d95p-11 +
                        T * (-0x1.3336527a26839p-23 +
                             T * (0x1.317a620fff6a8p-36 +
                                  T * -0x1.d10b26b6fb28dp-51)))) +
              0x1.3607cbfb15b57p+13);
    cp[sNO] =
        0x1.150d649310ce2p+8 *
        (0x1.a15b863893c54p+1 +
         T * (0x1.383da80ef9d95p-10 +
              T * (-0x1.ccd17bb739c56p-22 +
                   T * (0x1.317a620fff6a8p-34 + T * -0x1.22a6f8325cf98p-48))));
    s[sNO] =
        0x1.150d649310ce2p+8 *
        (0x1.a15b863893c54p+1 * log(T) +
         T * (0x1.383da80ef9d95p-10 +
              T * (-0x1.ccd17bb739c56p-23 +
                   T * (0x1.974dd815548ep-36 + T * -0x1.22a6f8325cf98p-50))) +
         0x1.97a2a7cab4c4ap+2);
    h[sN2O] = 0x1.79c0ba18e3f82p+7 *
              (T * (0x1.34ad39f4ee54p+2 +
                    T * (0x1.585455c7d19e1p-10 +
                         T * (-0x1.57104ac5ba269p-22 +
                              T * (0x1.5fdc00313900ap-35 +
                                   T * -0x1.19c077f8f0906p-49)))) +
               0x1.f8967a0f9096cp+12);
    cp[sN2O] =
        0x1.79c0ba18e3f82p+7 *
        (0x1.34ad39f4ee54p+2 +
         T * (0x1.585455c7d19e1p-9 +
              T * (-0x1.014c38144b9cfp-20 +
                   T * (0x1.5fdc00313900ap-33 + T * -0x1.603095f72cb47p-47))));
    s[sN2O] =
        0x1.79c0ba18e3f82p+7 *
        (0x1.34ad39f4ee54p+2 * log(T) +
         T * (0x1.585455c7d19e1p-9 +
              T * (-0x1.014c38144b9cfp-21 +
                   T * (0x1.d5255596f6ab8p-35 + T * -0x1.603095f72cb47p-49))) -
         0x1.19d1fbe0b68eap+1);
    h[sNO2] = 0x1.696a1b0a8416fp+7 *
              (T * (0x1.389fd0147fe9cp+2 +
                    T * (0x1.1cbd801ca8a9cp-10 +
                         T * (-0x1.2860a148121a7p-22 +
                              T * (0x1.5a4a99a75d176p-35 +
                                   T * -0x1.2ef4b74d1a737p-49)))) +
               0x1.218ff212d7732p+11);
    cp[sNO2] =
        0x1.696a1b0a8416fp+7 *
        (0x1.389fd0147fe9cp+2 +
         T * (0x1.1cbd801ca8a9cp-9 +
              T * (-0x1.bc90f1ec1b27ap-21 +
                   T * (0x1.5a4a99a75d176p-33 + T * -0x1.7ab1e52061104p-47))));
    s[sNO2] =
        0x1.696a1b0a8416fp+7 *
        (0x1.389fd0147fe9cp+2 * log(T) +
         T * (0x1.1cbd801ca8a9cp-9 +
              T * (-0x1.bc90f1ec1b27ap-22 +
                   T * (0x1.cdb8ccdf26c9dp-35 + T * -0x1.7ab1e52061104p-49))) -
         0x1.e0f09883efe43p-4);
    h[sNH] = 0x1.14cfff745b811p+9 *
             (T * (0x1.64500bb10e2a6p+1 +
                   T * (0x1.5c9c40c67750dp-11 +
                        T * (-0x1.3011d8266aad2p-23 +
                             T * (0x1.589491ea05ea6p-36 +
                                  T * -0x1.3d4f4cd0ccdbep-50)))) +
              0x1.4911b22d0e56p+15);
    cp[sNH] =
        0x1.14cfff745b811p+9 *
        (0x1.64500bb10e2a6p+1 +
         T * (0x1.5c9c40c67750dp-10 +
              T * (-0x1.c81ac439a003bp-22 +
                   T * (0x1.589491ea05ea6p-34 + T * -0x1.8ca320050012dp-48))));
    s[sNH] =
        0x1.14cfff745b811p+9 *
        (0x1.64500bb10e2a6p+1 * log(T) +
         T * (0x1.5c9c40c67750dp-10 +
              T * (-0x1.c81ac439a003bp-23 +
                   T * (0x1.cb70c28d5d388p-36 + T * -0x1.8ca320050012dp-50))) +
         0x1.6f68f019022f8p+2);
    h[sHNO] = 0x1.0c0c830fdcdfcp+8 *
              (T * (0x1.7d5817ef0a0e7p+1 +
                    T * (0x1.ca04ce1e70287p-10 +
                         T * (-0x1.1923fec713eaep-22 +
                              T * (0x1.f99889ffd471p-37 +
                                   T * -0x1.64af48192a2ap-55)))) +
               0x1.6f34a7ef9db23p+13);
    cp[sHNO] =
        0x1.0c0c830fdcdfcp+8 *
        (0x1.7d5817ef0a0e7p+1 +
         T * (0x1.ca04ce1e70287p-9 +
              T * (-0x1.a5b5fe2a9de05p-21 +
                   T * (0x1.f99889ffd471p-35 + T * -0x1.bddb1a1f74b48p-53))));
    s[sHNO] =
        0x1.0c0c830fdcdfcp+8 *
        (0x1.7d5817ef0a0e7p+1 * log(T) +
         T * (0x1.ca04ce1e70287p-9 +
              T * (-0x1.a5b5fe2a9de05p-22 +
                   T * (0x1.51105bffe2f6p-36 + T * -0x1.bddb1a1f74b48p-55))) +
         0x1.136767ee25e2fp+3);
    h[sNH2] = 0x1.0366cf1665468p+9 *
              (T * (0x1.6ad8d4420c141p+1 +
                    T * (0x1.a46367a785cacp-10 +
                         T * (-0x1.4e423aa696f9bp-22 +
                              T * (0x1.2d54c47ce6c2dp-35 +
                                   T * -0x1.c897a32c4c11fp-50)))) +
               0x1.5a6fd3f7ced91p+14);
    cp[sNH2] =
        0x1.0366cf1665468p+9 *
        (0x1.6ad8d4420c141p+1 +
         T * (0x1.a46367a785cacp-9 +
              T * (-0x1.f56357f9e2768p-21 +
                   T * (0x1.2d54c47ce6c2dp-33 + T * -0x1.1d5ec5fbaf8b3p-47))));
    s[sNH2] =
        0x1.0366cf1665468p+9 *
        (0x1.6ad8d4420c141p+1 * log(T) +
         T * (0x1.a46367a785cacp-9 +
              T * (-0x1.f56357f9e2768p-22 +
                   T * (0x1.91c65b5133ae7p-35 + T * -0x1.1d5ec5fbaf8b3p-49))) +
         0x1.a14e802b338a7p+2);
    h[sNNH] = 0x1.1e6cc1bcc36e4p+8 *
              (T * (0x1.e225021808348p+1 +
                    T * (0x1.7afeea4f6b578p-10 +
                         T * (-0x1.74d34a440fb4p-22 +
                              T * (0x1.725f60051ba0bp-35 +
                                   T * -0x1.22e10c914e65bp-49)))) +
               0x1.bfaac9ba5e354p+14);
    cp[sNNH] =
        0x1.1e6cc1bcc36e4p+8 *
        (0x1.e225021808348p+1 +
         T * (0x1.7afeea4f6b578p-9 +
              T * (-0x1.179e77b30bc7p-20 +
                   T * (0x1.725f60051ba0bp-33 + T * -0x1.6b994fb5a1ff2p-47))));
    s[sNNH] =
        0x1.1e6cc1bcc36e4p+8 *
        (0x1.e225021808348p+1 * log(T) +
         T * (0x1.7afeea4f6b578p-9 +
              T * (-0x1.179e77b30bc7p-21 +
                   T * (0x1.edd48006cf80fp-35 + T * -0x1.6b994fb5a1ff2p-49))) +
         0x1.1e1cc8224320ep+2);
    h[sCN] = 0x1.3f8958be799aep+8 *
             (T * (0x1.df7c49fd7a13cp+1 +
                   T * (0x1.6c7dd42c97ac9p-16 +
                        T * (0x1.a949952e905e8p-24 +
                             T * (-0x1.2def0fddeb7a1p-36 +
                                  T * 0x1.fcd5164ab9e94p-51)))) +
              0x1.92a0604189375p+15);
    cp[sCN] =
        0x1.3f8958be799aep+8 *
        (0x1.df7c49fd7a13cp+1 +
         T * (0x1.6c7dd42c97ac9p-15 +
              T * (0x1.3ef72fe2ec46ep-22 +
                   T * (-0x1.2def0fddeb7a1p-34 + T * 0x1.3e052deeb431cp-48))));
    s[sCN] =
        0x1.3f8958be799aep+8 *
        (0x1.df7c49fd7a13cp+1 * log(T) +
         T * (0x1.6c7dd42c97ac9p-15 +
              T * (0x1.3ef72fe2ec46ep-23 +
                   T * (-0x1.92941527e4a2cp-36 + T * 0x1.3e052deeb431cp-50))) +
         0x1.64b48e11a61abp+1);
    h[sNCO] = 0x1.8bbb85aa76addp+7 *
              (T * (0x1.49bd640e9d51bp+2 +
                    T * (0x1.2e24dfec0aebfp-10 +
                         T * (-0x1.3b153abdf5ee2p-22 +
                              T * (0x1.45373865bd7d7p-35 +
                                   T * -0x1.0639e98f3f73bp-49)))) +
               0x1.b5a0fbe76c8b4p+13);
    cp[sNCO] =
        0x1.8bbb85aa76addp+7 *
        (0x1.49bd640e9d51bp+2 +
         T * (0x1.2e24dfec0aebfp-9 +
              T * (-0x1.d89fd81cf0e53p-21 +
                   T * (0x1.45373865bd7d7p-33 + T * -0x1.47c863f30f509p-47))));
    s[sNCO] =
        0x1.8bbb85aa76addp+7 *
        (0x1.49bd640e9d51bp+2 * log(T) +
         T * (0x1.2e24dfec0aebfp-9 +
              T * (-0x1.d89fd81cf0e53p-22 +
                   T * (0x1.b19ef5dcfca74p-35 + T * -0x1.47c863f30f509p-49))) -
         0x1.45aa821f2990fp+1);
    h[sHCN] = 0x1.339e97ece7b7p+8 *
              (T * (0x1.e6afc62bc8dbap+1 +
                    T * (0x1.9c686e0cffc11p-10 +
                         T * (-0x1.7c8a6ce08d1dp-22 +
                              T * (0x1.6d78ea0eecc37p-35 +
                                   T * -0x1.1a7571994cb7dp-49)))) +
               0x1.c23a560418937p+13);
    cp[sHCN] =
        0x1.339e97ece7b7p+8 *
        (0x1.e6afc62bc8dbap+1 +
         T * (0x1.9c686e0cffc11p-9 +
              T * (-0x1.1d67d1a869d5cp-20 +
                   T * (0x1.6d78ea0eecc37p-33 + T * -0x1.6112cdff9fe5cp-47))));
    s[sHCN] =
        0x1.339e97ece7b7p+8 *
        (0x1.e6afc62bc8dbap+1 * log(T) +
         T * (0x1.9c686e0cffc11p-9 +
              T * (-0x1.1d67d1a869d5cp-21 +
                   T * (0x1.e74be2be91049p-35 + T * -0x1.6112cdff9fe5cp-49))) +
         0x1.93515a65a723cp+0);
    h[sHOCN] = 0x1.82763b11595dp+7 *
               (T * (0x1.79765b05e013dp+2 +
                     T * (0x1.9f38e1a731283p-10 +
                          T * (-0x1.9026cee799af7p-22 +
                               T * (0x1.85c3050860138p-35 +
                                    T * -0x1.2cbcb94ae7d7bp-49)))) -
                0x1.cf5110e022142p+11);
    cp[sHOCN] =
        0x1.82763b11595dp+7 *
        (0x1.79765b05e013dp+2 +
         T * (0x1.9f38e1a731283p-9 +
              T * (-0x1.2c1d1b2db3439p-20 +
                   T * (0x1.85c3050860138p-33 + T * -0x1.77ebe79da1cd9p-47))));
    s[sHOCN] =
        0x1.82763b11595dp+7 *
        (0x1.79765b05e013dp+2 * log(T) +
         T * (0x1.9f38e1a731283p-9 +
              T * (-0x1.2c1d1b2db3439p-21 +
                   T * (0x1.03d758b0400dp-34 + T * -0x1.77ebe79da1cd9p-49))) -
         0x1.8ba09dcf893fbp+2);
    h[sHNCO] = 0x1.82763b11595dp+7 *
               (T * (0x1.8e5538004c811p+2 +
                     T * (0x1.a0a17608fd1e1p-10 +
                          T * (-0x1.877b5818cda63p-22 +
                               T * (0x1.77735a207cf8ap-35 +
                                    T * -0x1.1ecba951ee841p-49)))) -
                0x1.044fbcd35a858p+14);
    cp[sHNCO] =
        0x1.82763b11595dp+7 *
        (0x1.8e5538004c811p+2 +
         T * (0x1.a0a17608fd1e1p-9 +
              T * (-0x1.259c82129a3cap-20 +
                   T * (0x1.77735a207cf8ap-33 + T * -0x1.667e93a66a251p-47))));
    s[sHNCO] =
        0x1.82763b11595dp+7 *
        (0x1.8e5538004c811p+2 * log(T) +
         T * (0x1.a0a17608fd1e1p-9 +
              T * (-0x1.259c82129a3cap-21 +
                   T * (0x1.f499cd80a6a0dp-35 + T * -0x1.667e93a66a251p-49))) -
         0x1.0c3b5eeb9dc85p+3);
    h[sH2CN] = 0x1.288f36628b032p+8 *
               (T * (0x1.4d6bc621b7e0bp+2 +
                     T * (0x1.8530e0556750fp-10 +
                          T * (-0x1.98d27283907c8p-24 +
                               T * (-0x1.67a94795f37d6p-35 +
                                    T * 0x1.b6946fa7f6caap-48)))) +
                0x1.b0746f9db22d1p+14);
    cp[sH2CN] =
        0x1.288f36628b032p+8 *
        (0x1.4d6bc621b7e0bp+2 +
         T * (0x1.8530e0556750fp-9 +
              T * (-0x1.329dd5e2ac5d6p-22 +
                   T * (-0x1.67a94795f37d6p-33 + T * 0x1.121cc5c8fa3eap-45))));
    s[sH2CN] =
        0x1.288f36628b032p+8 *
        (0x1.4d6bc621b7e0bp+2 * log(T) +
         T * (0x1.8530e0556750fp-9 +
              T * (-0x1.329dd5e2ac5d6p-23 +
                   T * (-0x1.df8c5f7299fc8p-35 + T * 0x1.121cc5c8fa3eap-47))) -
         0x1.1c7253da72a7cp+2);
    h[sHCNN] = 0x1.9533b6ca1714dp+7 *
               (T * (0x1.7941b83134557p+2 +
                     T * (0x1.05764fea2d6dcp-9 +
                          T * (-0x1.1e040d7e2a8fp-21 +
                               T * (0x1.4199ba4ec2266p-34 +
                                    T * -0x1.21984e547d5e8p-48)))) +
                0x1.a199e1cac0831p+15);
    cp[sHCNN] =
        0x1.9533b6ca1714dp+7 *
        (0x1.7941b83134557p+2 +
         T * (0x1.05764fea2d6dcp-8 +
              T * (-0x1.ad06143d3fd68p-20 +
                   T * (0x1.4199ba4ec2266p-32 + T * -0x1.69fe61e99cb62p-46))));
    s[sHCNN] =
        0x1.9533b6ca1714dp+7 *
        (0x1.7941b83134557p+2 * log(T) +
         T * (0x1.05764fea2d6dcp-8 +
              T * (-0x1.ad06143d3fd68p-21 +
                   T * (0x1.acccf86902dddp-34 + T * -0x1.69fe61e99cb62p-48))) -
         0x1.46985fddb6292p+2);
    h[sHCNO] = 0x1.82763b11595dp+7 *
               (T * (0x1.a64f89801bef2p+2 +
                     T * (0x1.8cdba5ee94b04p-10 +
                          T * (-0x1.817d26917d8d7p-22 +
                               T * (0x1.797fa9f8139b6p-35 +
                                    T * -0x1.24610f343c403p-49)))) +
                0x1.18b8891d14e3cp+14);
    cp[sHCNO] =
        0x1.82763b11595dp+7 *
        (0x1.a64f89801bef2p+2 +
         T * (0x1.8cdba5ee94b04p-9 +
              T * (-0x1.211ddced1e2a1p-20 +
                   T * (0x1.797fa9f8139b6p-33 + T * -0x1.6d7953014b503p-47))));
    s[sHCNO] =
        0x1.82763b11595dp+7 *
        (0x1.a64f89801bef2p+2 * log(T) +
         T * (0x1.8cdba5ee94b04p-9 +
              T * (-0x1.211ddced1e2a1p-21 +
                   T * (0x1.f754e2a01a248p-35 + T * -0x1.6d7953014b503p-49))) -
         0x1.4a94c4121328p+3);
    h[sNH3] = 0x1.e81a453143a8ep+8 *
              (T * (0x1.5135b9f630736p+1 +
                    T * (0x1.73580035f49d6p-9 +
                         T * (-0x1.3536c58500f2ap-21 +
                              T * (0x1.066c15fb54892p-34 +
                                   T * -0x1.6a8f0ef5ec3d4p-49)))) -
               0x1.990b21ff2e48fp+12);
    cp[sNH3] =
        0x1.e81a453143a8ep+8 *
        (0x1.5135b9f630736p+1 +
         T * (0x1.73580035f49d6p-8 +
              T * (-0x1.cfd22847816bfp-20 +
                   T * (0x1.066c15fb54892p-32 + T * -0x1.c532d2b3674c8p-47))));
    s[sNH3] =
        0x1.e81a453143a8ep+8 *
        (0x1.5135b9f630736p+1 * log(T) +
         T * (0x1.73580035f49d6p-8 +
              T * (-0x1.cfd22847816bfp-21 +
                   T * (0x1.5de572a470b6dp-34 + T * -0x1.c532d2b3674c8p-49))) +
         0x1.a43e2427fd751p+2);
    h[sCH2CHO] = 0x1.825174a3f53f5p+7 *
                 (T * (0x1.7e7160956c0d7p+2 +
                       T * (0x1.0a6c5738986f4p-8 +
                            T * (-0x1.eafd9a9007a33p-21 +
                                 T * (0x1.bf88df5305066p-34 +
                                      T * -0x1.3998db7e70d3fp-48)))) +
                  0x1.ea52617c1bda5p+8);
    cp[sCH2CHO] =
        0x1.825174a3f53f5p+7 *
        (0x1.7e7160956c0d7p+2 +
         T * (0x1.0a6c5738986f4p-7 +
              T * (-0x1.703e33ec05ba6p-19 +
                   T * (0x1.bf88df5305066p-32 + T * -0x1.87ff125e0d08ep-46))));
    s[sCH2CHO] =
        0x1.825174a3f53f5p+7 *
        (0x1.7e7160956c0d7p+2 * log(T) +
         T * (0x1.0a6c5738986f4p-7 +
              T * (-0x1.703e33ec05ba6p-20 +
                   T * (0x1.2a5b3f8cae044p-33 + T * -0x1.87ff125e0d08ep-48))) -
         0x1.42e56473471f8p+2);
    h[sCH3CHO] = 0x1.797a7ab03edafp+7 *
                 (T * (0x1.59dcf38b7d772p+2 +
                       T * (0x1.80242581cd52bp-8 +
                            T * (-0x1.7a2a05a17d829p-20 +
                                 T * (0x1.77e1ab96723c6p-33 +
                                      T * -0x1.2753ba556825fp-47)))) -
                  0x1.61047ced91687p+14);
    cp[sCH3CHO] =
        0x1.797a7ab03edafp+7 *
        (0x1.59dcf38b7d772p+2 +
         T * (0x1.80242581cd52bp-7 +
              T * (-0x1.1b9f84391e21fp-18 +
                   T * (0x1.77e1ab96723c6p-31 + T * -0x1.7128a8eac22f7p-45))));
    s[sCH3CHO] =
        0x1.797a7ab03edafp+7 *
        (0x1.59dcf38b7d772p+2 * log(T) +
         T * (0x1.80242581cd52bp-7 +
              T * (-0x1.1b9f84391e21fp-19 +
                   T * (0x1.f52ce4c898508p-33 + T * -0x1.7128a8eac22f7p-47))) -
         0x1.bd8a9519d8186p+1);
    h[sC3H8] = 0x1.791e6f137eed1p+7 *
               (T * (0x1.e22f4c1de5c41p+2 +
                     T * (0x1.3533e853aabbcp-7 +
                          T * (-0x1.18990808f2843p-19 +
                               T * (0x1.f6e487e665bdbp-33 +
                                    T * -0x1.58b5a95d75eb2p-47)))) -
                0x1.014e10624dd2fp+14);
    cp[sC3H8] =
        0x1.791e6f137eed1p+7 *
        (0x1.e22f4c1de5c41p+2 +
         T * (0x1.3533e853aabbcp-6 +
              T * (-0x1.a4e58c0d6bc65p-18 +
                   T * (0x1.f6e487e665bdbp-31 + T * -0x1.aee313b4d365ep-45))));
    s[sC3H8] =
        0x1.791e6f137eed1p+7 *
        (0x1.e22f4c1de5c41p+2 * log(T) +
         T * (0x1.3533e853aabbcp-6 +
              T * (-0x1.a4e58c0d6bc65p-19 +
                   T * (0x1.4f43054443d3dp-32 + T * -0x1.aee313b4d365ep-47))) -
         0x1.1e470fbeb9e49p+4);
    h[sC3H7] = 0x1.81f10d092f8b8p+7 *
               (T * (0x1.ecf903f7dc451p+2 +
                     T * (0x1.06de43cb39826p-7 +
                          T * (-0x1.d8be5fda28d19p-20 +
                               T * (0x1.a374b81e9310cp-33 +
                                    T * -0x1.1bd9f02a1e858p-47)))) +
                0x1.0353780346dc6p+13);
    cp[sC3H7] =
        0x1.81f10d092f8b8p+7 *
        (0x1.ecf903f7dc451p+2 +
         T * (0x1.06de43cb39826p-6 +
              T * (-0x1.628ec7e39e9d3p-18 +
                   T * (0x1.a374b81e9310cp-31 + T * -0x1.62d06c34a626ep-45))));
    s[sC3H7] =
        0x1.81f10d092f8b8p+7 *
        (0x1.ecf903f7dc451p+2 * log(T) +
         T * (0x1.06de43cb39826p-6 +
              T * (-0x1.628ec7e39e9d3p-19 +
                   T * (0x1.17a32569b7608p-32 + T * -0x1.62d06c34a626ep-47))) -
         0x1.ef5da272862f6p+3);
  } else if (T >= 299.999999) {
    h[sAR] = 0x1.a042152c51ec9p+7 *
             (T * (0x1.4p+1 +
                   T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) -
              0x1.74bp+9);
    cp[sAR] =
        0x1.a042152c51ec9p+7 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sAR] = 0x1.a042152c51ec9p+7 *
             (0x1.4p+1 * log(T) +
              T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
              0x1.176c8b439581p+2);
    h[sO] = 0x1.03d2b851eb852p+9 *
            (T * (0x1.9589c6bdbf12dp+1 +
                  T * (-0x1.add3ae5787a18p-10 +
                       T * (0x1.2934a583cf1ebp-19 +
                            T * (-0x1.a51e14d5c0503p-30 +
                                 T * 0x1.dbba8a635ed8p-42)))) +
             0x1.c709096bb98c8p+14);
    cp[sO] =
        0x1.03d2b851eb852p+9 *
        (0x1.9589c6bdbf12dp+1 +
         T * (-0x1.add3ae5787a18p-9 +
              T * (0x1.bdcef845b6aep-18 +
                   T * (-0x1.a51e14d5c0503p-28 + T * 0x1.2954967e1b47p-39))));
    s[sO] = 0x1.03d2b851eb852p+9 *
            (0x1.9589c6bdbf12dp+1 * log(T) +
             T * (-0x1.add3ae5787a18p-9 +
                  T * (0x1.bdcef845b6aep-19 + T * (-0x1.18beb88e80357p-29 +
                                                   T * 0x1.2954967e1b47p-41))) +
             0x1.06a5c1702251ep+1);
    h[sO2] = 0x1.03d2b851eb852p+8 *
             (T * (0x1.e42787ae5fa45p+1 +
                   T * (-0x1.88c9b66c8c0dfp-10 +
                        T * (0x1.b88f92d523df1p-19 +
                             T * (-0x1.4ca5927ada53cp-29 +
                                  T * 0x1.6d361ad59bca4p-41)))) -
              0x1.09fc63497b741p+10);
    cp[sO2] =
        0x1.03d2b851eb852p+8 *
        (0x1.e42787ae5fa45p+1 +
         T * (-0x1.88c9b66c8c0dfp-9 +
              T * (0x1.4a6bae1fdae75p-17 +
                   T * (-0x1.4ca5927ada53cp-27 + T * 0x1.c883a18b02bcdp-39))));
    s[sO2] =
        0x1.03d2b851eb852p+8 *
        (0x1.e42787ae5fa45p+1 * log(T) +
         T * (-0x1.88c9b66c8c0dfp-9 +
              T * (0x1.4a6bae1fdae75p-18 +
                   T * (-0x1.bb876df9231a5p-29 + T * 0x1.c883a18b02bcdp-41))) +
         0x1.d42eb7e3dc88dp+1);
    h[sH] = 0x1.01c2d34d34d35p+13 *
            (T * (0x1.4p+1 + T * (0x1.8d112bff6d58p-42 +
                                  T * (-0x1.7f85eab1cab54p-51 +
                                       T * (0x1.538a620b2373fp-61 +
                                            T * -0x1.c09fb33011ea8p-73)))) +
             0x1.8e06a3bcd35a8p+14);
    cp[sH] =
        0x1.01c2d34d34d35p+13 *
        (0x1.4p+1 +
         T * (0x1.8d112bff6d58p-41 +
              T * (-0x1.1fa470055807fp-49 +
                   T * (0x1.538a620b2373fp-59 + T * -0x1.1863cffe0b329p-70))));
    s[sH] =
        0x1.01c2d34d34d35p+13 *
        (0x1.4p+1 * log(T) +
         T * (0x1.8d112bff6d58p-41 +
              T * (-0x1.1fa470055807fp-50 +
                   T * (0x1.c4b882b9849a9p-61 + T * -0x1.1863cffe0b329p-72))) -
         0x1.c9673ad546a18p-2);
    h[sOH] = 0x1.e8d94973d693ep+8 *
             (T * (0x1.fefa5c927d1abp+1 +
                   T * (-0x1.3abed86e71c7dp-10 +
                        T * (0x1.9d34c5377fa67p-20 +
                             T * (-0x1.0ab59e9e67441p-30 +
                                  T * 0x1.332bdbce74b0fp-42)))) +
              0x1.c3e293f290abbp+11);
    cp[sOH] =
        0x1.e8d94973d693ep+8 *
        (0x1.fefa5c927d1abp+1 +
         T * (-0x1.3abed86e71c7dp-9 +
              T * (0x1.35e793e99fbcdp-18 +
                   T * (-0x1.0ab59e9e67441p-28 + T * 0x1.7ff6d2c211dd3p-40))));
    s[sOH] =
        0x1.e8d94973d693ep+8 *
        (0x1.fefa5c927d1abp+1 * log(T) +
         T * (-0x1.3abed86e71c7dp-9 +
              T * (0x1.35e793e99fbcdp-19 +
                   T * (-0x1.639cd37ddf057p-30 + T * 0x1.7ff6d2c211dd3p-42))) -
         0x1.a9adbdb54f242p-4);
    h[sH2] = 0x1.01c2d34d34d35p+12 *
             (T * (0x1.2c130ac9b2911p+1 +
                   T * (0x1.058175d02a941p-8 +
                        T * (-0x1.b3b80759749dp-18 +
                             T * (0x1.5a4c582f87e3dp-28 +
                                  T * -0x1.9f3d0ec308a15p-40)))) -
              0x1.caf7b3bfb58d1p+9);
    cp[sH2] =
        0x1.01c2d34d34d35p+12 *
        (0x1.2c130ac9b2911p+1 +
         T * (0x1.058175d02a941p-7 +
              T * (-0x1.46ca05831775cp-16 +
                   T * (0x1.5a4c582f87e3dp-26 + T * -0x1.03862939e564dp-37))));
    s[sH2] =
        0x1.01c2d34d34d35p+12 *
        (0x1.2c130ac9b2911p+1 * log(T) +
         T * (0x1.058175d02a941p-7 +
              T * (-0x1.46ca05831775cp-17 +
                   T * (0x1.cdbb203f5fda7p-28 + T * -0x1.03862939e564dp-39))) +
         0x1.5db38496161b4p-1);
    h[sHO2] = 0x1.f7c6fae98a1dp+7 *
              (T * (0x1.1350a899bcaa1p+2 +
                    T * (-0x1.373d054674594p-9 +
                         T * (0x1.d94d8bda368ddp-18 +
                              T * (-0x1.a110b0905e0fcp-28 +
                                   T * 0x1.058dba0c9e34ap-39)))) +
               0x1.26cedbb59ddc2p+8);
    cp[sHO2] =
        0x1.f7c6fae98a1dp+7 *
        (0x1.1350a899bcaa1p+2 +
         T * (-0x1.373d054674594p-8 +
              T * (0x1.62fa28e3a8ea6p-16 +
                   T * (-0x1.a110b0905e0fcp-26 + T * 0x1.46f1288fc5c1cp-37))));
    s[sHO2] =
        0x1.f7c6fae98a1dp+7 *
        (0x1.1350a899bcaa1p+2 * log(T) +
         T * (-0x1.373d054674594p-8 +
              T * (0x1.62fa28e3a8ea6p-17 +
                   T * (-0x1.160b20603eb53p-27 + T * 0x1.46f1288fc5c1cp-39))) +
         0x1.dbbb985c82b7dp+1);
    h[sH2O2] = 0x1.e8d94973d693ep+7 *
               (T * (0x1.11abd48f63e0ap+2 +
                     T * (-0x1.1c98643a78c6ap-12 +
                          T * (0x1.7652d932e9588p-18 +
                               T * (-0x1.72b101d35ae59p-28 +
                                    T * 0x1.e584c5d1abf9dp-40)))) -
                0x1.149a541205bcp+14);
    cp[sH2O2] =
        0x1.e8d94973d693ep+7 *
        (0x1.11abd48f63e0ap+2 +
         T * (-0x1.1c98643a78c6ap-11 +
              T * (0x1.18be22e62f026p-16 +
                   T * (-0x1.72b101d35ae59p-26 + T * 0x1.2f72fba30b7c2p-37))));
    s[sH2O2] =
        0x1.e8d94973d693ep+7 *
        (0x1.11abd48f63e0ap+2 * log(T) +
         T * (-0x1.1c98643a78c6ap-11 +
              T * (0x1.18be22e62f026p-17 +
                   T * (-0x1.ee4157c479321p-28 + T * 0x1.2f72fba30b7c2p-39))) +
         0x1.b7afbe1e3346dp+1);
    h[sCH] = 0x1.3f5713b4542aap+9 *
             (T * (0x1.beb24fde64a4cp+1 +
                   T * (0x1.5390f0ed10ec5p-13 +
                        T * (-0x1.2e41b3d3fb43bp-21 +
                             T * (0x1.b29b14b88d49cp-31 +
                                  T * -0x1.3c9f9bc981c96p-42)))) +
              0x1.148d4b1c432cap+16);
    cp[sCH] =
        0x1.3f5713b4542aap+9 *
        (0x1.beb24fde64a4cp+1 +
         T * (0x1.5390f0ed10ec5p-12 +
              T * (-0x1.c5628dbdf8e58p-20 +
                   T * (0x1.b29b14b88d49cp-29 + T * -0x1.8bc782bbe23bbp-40))));
    s[sCH] =
        0x1.3f5713b4542aap+9 *
        (0x1.beb24fde64a4cp+1 * log(T) +
         T * (0x1.5390f0ed10ec5p-12 +
              T * (-0x1.c5628dbdf8e58p-21 +
                   T * (0x1.21bcb87b08dbdp-30 + T * -0x1.8bc782bbe23bbp-42))) +
         0x1.0ac0e0048d028p+1);
    h[sCO] = 0x1.28d5af05f0d41p+8 *
             (T * (0x1.ca2e271a4b2fdp+1 +
                   T * (-0x1.400048c1ba1bp-12 +
                        T * (0x1.6bee992797344p-22 +
                             T * (0x1.f2a1bae3bf93ap-33 +
                                  T * -0x1.97510ba095003p-43)))) -
              0x1.c040b020c49bap+13);
    cp[sCO] =
        0x1.28d5af05f0d41p+8 *
        (0x1.ca2e271a4b2fdp+1 +
         T * (-0x1.400048c1ba1bp-11 +
              T * (0x1.10f2f2ddb1673p-20 +
                   T * (0x1.f2a1bae3bf93ap-31 + T * -0x1.fd254e88ba403p-41))));
    s[sCO] =
        0x1.28d5af05f0d41p+8 *
        (0x1.ca2e271a4b2fdp+1 * log(T) +
         T * (-0x1.400048c1ba1bp-11 +
              T * (0x1.10f2f2ddb1673p-21 +
                   T * (0x1.4c6bd1ed2a627p-32 + T * -0x1.fd254e88ba403p-43))) +
         0x1.c1138e274a9cbp+1);
    h[s3XCH2] = 0x1.2863e91362288p+9 *
                (T * (0x1.e19f746480de1p+1 +
                      T * (0x1.fbf7d158743a2p-12 +
                           T * (0x1.f42aa333e052bp-21 +
                                T * (-0x1.08a1f3bb46844p-30 +
                                     T * 0x1.7bf8fa708a0d4p-42)))) +
                 0x1.67681487fcb92p+15);
    cp[s3XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.e19f746480de1p+1 +
         T * (0x1.fbf7d158743a2p-11 +
              T * (0x1.771ffa66e83ep-19 +
                   T * (-0x1.08a1f3bb46844p-28 + T * 0x1.daf7390cac908p-40))));
    s[s3XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.e19f746480de1p+1 * log(T) +
         T * (0x1.fbf7d158743a2p-11 +
              T * (0x1.771ffa66e83ep-20 +
                   T * (-0x1.60d7efa45e05bp-30 + T * 0x1.daf7390cac908p-42))) +
         0x1.9002165ab5584p+0);
    h[sHCO] = 0x1.1e86068742b68p+8 *
              (T * (0x1.0e27e8a748d9cp+2 +
                    T * (-0x1.a9301251f436fp-10 +
                         T * (0x1.34408c6c7a954p-18 +
                              T * (-0x1.c97ac8618dbf3p-29 +
                                   T * 0x1.e8615cf58dbadp-41)))) +
               0x1.dff21426fe719p+11);
    cp[sHCO] =
        0x1.1e86068742b68p+8 *
        (0x1.0e27e8a748d9cp+2 +
         T * (-0x1.a9301251f436fp-9 +
              T * (0x1.ce60d2a2b7dfep-17 +
                   T * (-0x1.c97ac8618dbf3p-27 + T * 0x1.313cda197894cp-38))));
    s[sHCO] =
        0x1.1e86068742b68p+8 *
        (0x1.0e27e8a748d9cp+2 * log(T) +
         T * (-0x1.a9301251f436fp-9 +
              T * (0x1.ce60d2a2b7dfep-18 +
                   T * (-0x1.30fc85965e7f7p-28 + T * 0x1.313cda197894cp-40))) +
         0x1.b27acbb8a5a37p+1);
    h[s1XCH2] = 0x1.2863e91362288p+9 *
                (T * (0x1.0cb5ee035346ap+2 +
                      T * (-0x1.36326518bab75p-10 +
                           T * (0x1.7056247403375p-19 +
                                T * (-0x1.cb9b5a07323d1p-30 +
                                     T * 0x1.b58ed1c91768cp-42)))) +
                 0x1.8a81a1f212d77p+15);
    cp[s1XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.0cb5ee035346ap+2 +
         T * (-0x1.36326518bab75p-9 +
              T * (0x1.14409b5702698p-17 +
                   T * (-0x1.cb9b5a07323d1p-28 + T * 0x1.1179431daea17p-39))));
    s[s1XCH2] =
        0x1.2863e91362288p+9 *
        (0x1.0cb5ee035346ap+2 * log(T) +
         T * (-0x1.36326518bab75p-9 +
              T * (0x1.14409b5702698p-18 +
                   T * (-0x1.3267915a217e1p-29 + T * 0x1.1179431daea17p-41))) -
         0x1.89c9f613ff21ap-1);
    h[sCH3] = 0x1.1484949ef95bfp+9 *
              (T * (0x1.d63835d17324cp+1 +
                    T * (0x1.079458000b018p-10 +
                         T * (0x1.005d9a7f42e15p-19 +
                              T * (-0x1.d82ef9cab9cf9p-30 +
                                   T * 0x1.1e69b21fd8e05p-41)))) +
               0x1.00f3fec56d5dp+14);
    cp[sCH3] =
        0x1.1484949ef95bfp+9 *
        (0x1.d63835d17324cp+1 +
         T * (0x1.079458000b018p-9 +
              T * (0x1.808c67bee451fp-18 +
                   T * (-0x1.d82ef9cab9cf9p-28 + T * 0x1.66041ea7cf186p-39))));
    s[sCH3] =
        0x1.1484949ef95bfp+9 *
        (0x1.d63835d17324cp+1 * log(T) +
         T * (0x1.079458000b018p-9 +
              T * (0x1.808c67bee451fp-19 +
                   T * (-0x1.3ac9fbdc7bdfbp-29 + T * 0x1.66041ea7cf186p-41))) +
         0x1.9ac4ba59ad9b6p+0);
    h[sCH2O] = 0x1.14e79947885d3p+8 *
               (T * (0x1.32cc5c2ed4ffap+2 +
                     T * (-0x1.44ad209405144p-8 +
                          T * (0x1.a17083e6a66a1p-17 +
                               T * (-0x1.45cdb600576d1p-27 +
                                    T * 0x1.72e833dfa2eeep-39)))) -
                0x1.bf27a7525460bp+13);
    cp[sCH2O] =
        0x1.14e79947885d3p+8 *
        (0x1.32cc5c2ed4ffap+2 +
         T * (-0x1.44ad209405144p-7 +
              T * (0x1.391462ecfccf9p-15 +
                   T * (-0x1.45cdb600576d1p-25 + T * 0x1.cfa240d78baa9p-37))));
    s[sCH2O] =
        0x1.14e79947885d3p+8 *
        (0x1.32cc5c2ed4ffap+2 * log(T) +
         T * (-0x1.44ad209405144p-7 +
              T * (0x1.391462ecfccf9p-16 +
                   T * (-0x1.b2679d55c9e6cp-27 + T * 0x1.cfa240d78baa9p-39))) +
         0x1.34a3e47636bep-1);
    h[sCH4] = 0x1.03249373f39ddp+9 *
              (T * (0x1.4997920d33445p+2 +
                    T * (-0x1.bff87b6cd6efap-8 +
                         T * (0x1.1308ea9253b9fp-16 +
                              T * (-0x1.a0641e640b81dp-27 +
                                   T * 0x1.d533a7731d0cap-39)))) -
               0x1.40352e48e8a72p+13);
    cp[sCH4] =
        0x1.03249373f39ddp+9 *
        (0x1.4997920d33445p+2 +
         T * (-0x1.bff87b6cd6efap-7 +
              T * (0x1.9c8d5fdb7d96ep-15 +
                   T * (-0x1.a0641e640b81dp-25 + T * 0x1.254048a7f227ep-36))));
    s[sCH4] =
        0x1.03249373f39ddp+9 *
        (0x1.4997920d33445p+2 * log(T) +
         T * (-0x1.bff87b6cd6efap-7 +
              T * (0x1.9c8d5fdb7d96ep-16 +
                   T * (-0x1.15981442b2569p-26 + T * 0x1.254048a7f227ep-38))) -
         0x1.290b1eed001ep+2);
    h[sCO2] = 0x1.79d6b3468efa8p+7 *
              (T * (0x1.2daac1343d496p+1 +
                    T * (0x1.266842a5bf315p-8 +
                         T * (-0x1.3eb3eab0b4aa9p-19 +
                              T * (0x1.51fd10511d1a8p-31 +
                                   T * -0x1.02ddb83a1a923p-45)))) -
               0x1.79e7f07c84b5ep+15);
    cp[sCO2] =
        0x1.79d6b3468efa8p+7 *
        (0x1.2daac1343d496p+1 +
         T * (0x1.266842a5bf315p-7 +
              T * (-0x1.de0de0090effep-18 +
                   T * (0x1.51fd10511d1a8p-29 + T * -0x1.43952648a136cp-43))));
    s[sCO2] =
        0x1.79d6b3468efa8p+7 *
        (0x1.2daac1343d496p+1 * log(T) +
         T * (0x1.266842a5bf315p-7 +
              T * (-0x1.de0de0090effep-19 +
                   T * (0x1.c2a6c06c26cep-31 + T * -0x1.43952648a136cp-45))) +
         0x1.3cd56b771c6c2p+3);
    h[sCH2OH] = 0x1.0be9223bc1817p+8 *
                (T * (0x1.ee93ebafbbefep+1 +
                      T * (0x1.6ec96e6bebfdfp-9 +
                           T * (0x1.096ce0decb881p-19 +
                                T * (-0x1.672b5483bb261p-29 +
                                     T * 0x1.ebfb3e2ef3b1cp-41)))) -
                 0x1.8f3d3cc8de2acp+11);
    cp[sCH2OH] =
        0x1.0be9223bc1817p+8 *
        (0x1.ee93ebafbbefep+1 +
         T * (0x1.6ec96e6bebfdfp-8 +
              T * (0x1.8e23514e314c1p-18 +
                   T * (-0x1.672b5483bb261p-27 + T * 0x1.337d06dd584f1p-38))));
    s[sCH2OH] =
        0x1.0be9223bc1817p+8 *
        (0x1.ee93ebafbbefep+1 * log(T) +
         T * (0x1.6ec96e6bebfdfp-8 +
              T * (0x1.8e23514e314c1p-19 +
                   T * (-0x1.dee470afa432cp-29 + T * 0x1.337d06dd584f1p-40))) +
         0x1.5e45ffdec7f7p+2);
    h[sCH3O] = 0x1.0be9223bc1817p+8 *
               (T * (0x1.0d9817b95a294p+1 +
                     T * (0x1.d8f25f83733c9p-9 +
                          T * (0x1.ddadaadf4f2cbp-20 +
                               T * (-0x1.fafcbebd8fc7ap-30 +
                                    T * 0x1.d362c52c6841ap-42)))) +
                0x1.e94cf0d844d01p+9);
    cp[sCH3O] =
        0x1.0be9223bc1817p+8 *
        (0x1.0d9817b95a294p+1 +
         T * (0x1.d8f25f83733c9p-8 +
              T * (0x1.664240277b618p-18 +
                   T * (-0x1.fafcbebd8fc7ap-28 + T * 0x1.241dbb3bc129p-39))));
    s[sCH3O] =
        0x1.0be9223bc1817p+8 *
        (0x1.0d9817b95a294p+1 * log(T) +
         T * (0x1.d8f25f83733c9p-8 +
              T * (0x1.664240277b618p-19 +
                   T * (-0x1.51fdd47e5fda7p-29 + T * 0x1.241dbb3bc129p-41))) +
         0x1.a4dea24cc6823p+3);
    h[sCH3OH] = 0x1.037b88ab2ad7p+8 *
                (T * (0x1.6dc90b8ca6163p+2 +
                      T * (-0x1.f316286598c3ep-8 +
                           T * (0x1.6cdf1d364a3a3p-16 +
                                T * (-0x1.314a0b40b1447p-26 +
                                     T * 0x1.6fd23baa36b93p-38)))) -
                 0x1.90ab0ff972474p+14);
    cp[sCH3OH] =
        0x1.037b88ab2ad7p+8 *
        (0x1.6dc90b8ca6163p+2 +
         T * (-0x1.f316286598c3ep-7 +
              T * (0x1.11a755e8b7abap-14 +
                   T * (-0x1.314a0b40b1447p-24 + T * 0x1.cbc6ca94c4677p-36))));
    s[sCH3OH] =
        0x1.037b88ab2ad7p+8 *
        (0x1.6dc90b8ca6163p+2 * log(T) +
         T * (-0x1.f316286598c3ep-7 +
              T * (0x1.11a755e8b7abap-15 +
                   T * (-0x1.970d645641b09p-26 + T * 0x1.cbc6ca94c4677p-38))) -
         0x1.810c94e3d24cfp+0);
    h[sC2H] = 0x1.4c3397c02c83cp+8 *
              (T * (0x1.71e04a987f933p+1 +
                    T * (0x1.b76ae82ebca6cp-8 +
                         T * (-0x1.3e82612c20167p-17 +
                              T * (0x1.fa727902a68efp-28 +
                                   T * -0x1.33bda806adcfp-39)))) +
               0x1.051764a8c154dp+16);
    cp[sC2H] =
        0x1.4c3397c02c83cp+8 *
        (0x1.71e04a987f933p+1 +
         T * (0x1.b76ae82ebca6cp-7 +
              T * (-0x1.ddc391c23021ap-16 +
                   T * (0x1.fa727902a68efp-26 + T * -0x1.80ad12085942bp-37))));
    s[sC2H] =
        0x1.4c3397c02c83cp+8 *
        (0x1.71e04a987f933p+1 * log(T) +
         T * (0x1.b76ae82ebca6cp-7 +
              T * (-0x1.ddc391c23021ap-17 +
                   T * (0x1.51a1a601c45f5p-27 + T * -0x1.80ad12085942bp-39))) +
         0x1.8e450c6411777p+2);
    h[sC2H2] = 0x1.3f5713b4542aap+8 *
               (T * (0x1.9e0b72c73f3bap-1 +
                     T * (0x1.7ec17f28e481ap-7 +
                          T * (-0x1.8d40c15d06efbp-17 +
                               T * (0x1.e14c5845a203fp-28 +
                                    T * -0x1.de8c6d30d264p-40)))) +
                0x1.9cf3ec3c9eeccp+14);
    cp[sC2H2] =
        0x1.3f5713b4542aap+8 *
        (0x1.9e0b72c73f3bap-1 +
         T * (0x1.7ec17f28e481ap-6 +
              T * (-0x1.29f09105c533cp-15 +
                   T * (0x1.e14c5845a203fp-26 + T * -0x1.2b17c43e837e8p-37))));
    s[sC2H2] =
        0x1.3f5713b4542aap+8 *
        (0x1.9e0b72c73f3bap-1 * log(T) +
         T * (0x1.7ec17f28e481ap-6 +
              T * (-0x1.29f09105c533cp-16 +
                   T * (0x1.40dd902e6c02ap-27 + T * -0x1.2b17c43e837e8p-39))) +
         0x1.be12106e0c4d1p+3);
    h[sHCCO] = 0x1.954cff46b5264p+7 *
               (T * (0x1.203868265a06ep+1 +
                     T * (0x1.2142867388492p-7 +
                          T * (-0x1.0967cefa39fe4p-17 +
                               T * (0x1.28cb9772e0568p-28 +
                                    T * -0x1.1d37b00a8be7ep-40)))) +
                0x1.396dcbc6a7efap+14);
    cp[sHCCO] =
        0x1.954cff46b5264p+7 *
        (0x1.203868265a06ep+1 +
         T * (0x1.2142867388492p-6 +
              T * (-0x1.8e1bb67756fd6p-16 +
                   T * (0x1.28cb9772e0568p-26 + T * -0x1.64859c0d2ee1dp-38))));
    s[sHCCO] =
        0x1.954cff46b5264p+7 *
        (0x1.203868265a06ep+1 * log(T) +
         T * (0x1.2142867388492p-6 +
              T * (-0x1.8e1bb67756fd6p-17 +
                   T * (0x1.8bba1f43d5c8bp-28 + T * -0x1.64859c0d2ee1dp-40))) +
         0x1.8fb17efe0ce0cp+3);
    h[sC2H3] = 0x1.3370009b17838p+8 *
               (T * (0x1.9b3219c31fa4ep+1 +
                     T * (0x1.8d17f1df63fcdp-11 +
                          T * (0x1.21ebbad5b9bbdp-17 +
                               T * (-0x1.3339cad4b897fp-27 +
                                    T * 0x1.9e3160f1cd55p-39)))) +
                0x1.1057b18fc5048p+15);
    cp[sC2H3] =
        0x1.3370009b17838p+8 *
        (0x1.9b3219c31fa4ep+1 +
         T * (0x1.8d17f1df63fcdp-10 +
              T * (0x1.b2e198409699cp-16 +
                   T * (-0x1.3339cad4b897fp-25 + T * 0x1.02dedc9720552p-36))));
    s[sC2H3] =
        0x1.3370009b17838p+8 *
        (0x1.9b3219c31fa4ep+1 * log(T) +
         T * (0x1.8d17f1df63fcdp-10 +
              T * (0x1.b2e198409699cp-17 +
                   T * (-0x1.99a263c64b754p-27 + T * 0x1.02dedc9720552p-38))) +
         0x1.10565881a1555p+3);
    h[sCH2CO] = 0x1.8b94f63b4adbcp+7 *
                (T * (0x1.116315790e08dp+1 +
                      T * (0x1.28dc0ec708b6bp-7 +
                           T * (-0x1.851d295f1200dp-18 +
                                T * (0x1.410e7ab1ff68dp-29 +
                                     T * -0x1.c5a468865bdep-42)))) -
                 0x1.b82eb04ab606bp+12);
    cp[sCH2CO] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.116315790e08dp+1 +
         T * (0x1.28dc0ec708b6bp-6 +
              T * (-0x1.23d5df074d80ap-16 +
                   T * (0x1.410e7ab1ff68dp-27 + T * -0x1.1b86c153f96acp-39))));
    s[sCH2CO] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.116315790e08dp+1 * log(T) +
         T * (0x1.28dc0ec708b6bp-6 +
              T * (-0x1.23d5df074d80ap-17 +
                   T * (0x1.ac134e42a9e11p-29 + T * -0x1.1b86c153f96acp-41))) +
         0x1.86e696a26e547p+3);
    h[sC2H4] = 0x1.2863e91362288p+8 *
               (T * (0x1.fac71d356ff96p+1 +
                     T * (-0x1.f0244a6c1abf8p-9 +
                          T * (0x1.3f5227836c76ep-16 +
                               T * (-0x1.2908fcd07cb16p-26 +
                                    T * 0x1.7bd417ca8706dp-38)))) +
                0x1.3e1c6a35935fcp+12);
    cp[sC2H4] =
        0x1.2863e91362288p+8 *
        (0x1.fac71d356ff96p+1 +
         T * (-0x1.f0244a6c1abf8p-8 +
              T * (0x1.defb3b4522b25p-15 +
                   T * (-0x1.2908fcd07cb16p-24 + T * 0x1.dac91dbd28c88p-36))));
    s[sC2H4] =
        0x1.2863e91362288p+8 *
        (0x1.fac71d356ff96p+1 * log(T) +
         T * (-0x1.f0244a6c1abf8p-8 +
              T * (0x1.defb3b4522b25p-16 +
                   T * (-0x1.8c0bfbc0a641dp-26 + T * 0x1.dac91dbd28c88p-38))) +
         0x1.063aaba285a67p+2);
    h[sC2H5] = 0x1.1e1c03861415ep+8 *
               (T * (0x1.139d223a3c91dp+2 +
                     T * (-0x1.125f4e7e42173p-9 +
                          T * (0x1.1605bc90bd661p-16 +
                               T * (-0x1.0152aabf3f38cp-26 +
                                    T * 0x1.44699f0472426p-38)))) +
                0x1.914d03126e979p+13);
    cp[sC2H5] =
        0x1.1e1c03861415ep+8 *
        (0x1.139d223a3c91dp+2 +
         T * (-0x1.125f4e7e42173p-8 +
              T * (0x1.a1089ad91c191p-15 +
                   T * (-0x1.0152aabf3f38cp-24 + T * 0x1.958406c58ed2fp-36))));
    s[sC2H5] =
        0x1.1e1c03861415ep+8 *
        (0x1.139d223a3c91dp+2 * log(T) +
         T * (-0x1.125f4e7e42173p-8 +
              T * (0x1.a1089ad91c191p-16 +
                   T * (-0x1.5718e3a9a9a1p-26 + T * 0x1.958406c58ed2fp-38))) +
         0x1.2d42ea8b4ea83p+2);
    h[sC2H6] = 0x1.1484949ef95bfp+8 *
               (T * (0x1.12a6b4b528ec3p+2 +
                     T * (-0x1.688c91f95b2dfp-9 +
                          T * (0x1.4f3aed97cf26bp-16 +
                               T * (-0x1.3048b11b68ad6p-26 +
                                    T * 0x1.7a244044f8d16p-38)))) -
                0x1.6811a4dd2f1aap+13);
    cp[sC2H6] =
        0x1.1484949ef95bfp+8 *
        (0x1.12a6b4b528ec3p+2 +
         T * (-0x1.688c91f95b2dfp-8 +
              T * (0x1.f6d86463b6bap-15 +
                   T * (-0x1.3048b11b68ad6p-24 + T * 0x1.d8ad50563705bp-36))));
    s[sC2H6] =
        0x1.1484949ef95bfp+8 *
        (0x1.12a6b4b528ec3p+2 * log(T) +
         T * (-0x1.688c91f95b2dfp-8 +
              T * (0x1.f6d86463b6bap-16 +
                   T * (-0x1.95b64179e0e73p-26 + T * 0x1.d8ad50563705bp-38))) +
         0x1.555a7618352bp+1);
    h[sH2O] = 0x1.cd7f5ff1730a8p+8 *
              (T * (0x1.0cb686e536fbfp+2 +
                    T * (-0x1.0aeb63b84c925p-10 +
                         T * (0x1.23b7c52ad124bp-19 +
                              T * (-0x1.7921667072d59p-30 +
                                   T * 0x1.8f03963eb52f4p-42)))) -
               0x1.d956e8240b78p+14);
    cp[sH2O] =
        0x1.cd7f5ff1730a8p+8 *
        (0x1.0cb686e536fbfp+2 +
         T * (-0x1.0aeb63b84c925p-9 +
              T * (0x1.b593a7c039b7p-18 +
                   T * (-0x1.7921667072d59p-28 + T * 0x1.f2c47bce627b1p-40))));
    s[sH2O] =
        0x1.cd7f5ff1730a8p+8 *
        (0x1.0cb686e536fbfp+2 * log(T) +
         T * (-0x1.0aeb63b84c925p-9 +
              T * (0x1.b593a7c039b7p-19 +
                   T * (-0x1.f6d73340991ccp-30 + T * 0x1.f2c47bce627b1p-42))) -
         0x1.b2b4597d38a9bp-1);
    h[sN2] = 0x1.28ba905aa1e55p+8 *
             (T * (0x1.a63b0c4588a05p+1 +
                   T * (0x1.712969da0405p-11 +
                        T * (-0x1.629f839620fc5p-20 +
                             T * (0x1.83ae94da0fe59p-30 +
                                  T * -0x1.13441e6a0b3dfp-41)))) -
              0x1.fe732fec56d5dp+9);
    cp[sN2] =
        0x1.28ba905aa1e55p+8 *
        (0x1.a63b0c4588a05p+1 +
         T * (0x1.712969da0405p-10 +
              T * (-0x1.09f7a2b098bd4p-18 +
                   T * (0x1.83ae94da0fe59p-28 + T * -0x1.581526048e0d6p-39))));
    s[sN2] =
        0x1.28ba905aa1e55p+8 *
        (0x1.a63b0c4588a05p+1 * log(T) +
         T * (0x1.712969da0405p-10 +
              T * (-0x1.09f7a2b098bd4p-19 +
                   T * (0x1.0274633c0a991p-29 + T * -0x1.581526048e0d6p-41))) +
         0x1.f9a5ca29845ddp+1);
    h[sC] = 0x1.5a24731aa4fb9p+9 *
            (T * (0x1.46f15252b32b4p+1 +
                  T * (-0x1.51282024e7f51p-13 +
                       T * (0x1.06a26cec32fb4p-22 +
                            T * (-0x1.928ce62f9cbp-33 +
                                 T * 0x1.e01f5296c2ba4p-45)))) +
             0x1.4dc3e219652bdp+16);
    cp[sC] =
        0x1.5a24731aa4fb9p+9 *
        (0x1.46f15252b32b4p+1 +
         T * (-0x1.51282024e7f51p-12 +
              T * (0x1.89f3a3624c78ep-21 +
                   T * (-0x1.928ce62f9cbp-31 + T * 0x1.2c13939e39b46p-42))));
    s[sC] =
        0x1.5a24731aa4fb9p+9 *
        (0x1.46f15252b32b4p+1 * log(T) +
         T * (-0x1.51282024e7f51p-12 +
              T * (0x1.89f3a3624c78ep-22 +
                   T * (-0x1.0c5deeca68755p-32 + T * 0x1.2c13939e39b46p-44))) +
         0x1.2200f5486bff8p+2);
    h[sHCCOH] = 0x1.8b94f63b4adbcp+7 *
                (T * (0x1.3e0c2d34ec70dp+0 +
                      T * (0x1.fd1641c705f4ap-7 +
                           T * (-0x1.1c77d6cfea29fp-16 +
                                T * (0x1.728b8de31e2a5p-27 +
                                     T * -0x1.8a79cae1a9a06p-39)))) +
                 0x1.f5f9d42c3c9efp+12);
    cp[sHCCOH] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.3e0c2d34ec70dp+0 +
         T * (0x1.fd1641c705f4ap-6 +
              T * (-0x1.aab3c237df3eep-15 +
                   T * (0x1.728b8de31e2a5p-25 + T * -0x1.ed183d9a14087p-37))));
    s[sHCCOH] =
        0x1.8b94f63b4adbcp+7 *
        (0x1.3e0c2d34ec70dp+0 * log(T) +
         T * (0x1.fd1641c705f4ap-6 +
              T * (-0x1.aab3c237df3eep-16 +
                   T * (0x1.ee0f67d97d8dcp-27 + T * -0x1.ed183d9a14087p-39))) +
         0x1.bbfa6bd6e8af8p+3);
    h[sN] = 0x1.28ba905aa1e55p+9 *
            (T * (0x1.4p+1 +
                  T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0)))) +
             0x1.b6514624dd2f2p+15);
    cp[sN] =
        0x1.28ba905aa1e55p+9 *
        (0x1.4p+1 + T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))));
    s[sN] = 0x1.28ba905aa1e55p+9 *
            (0x1.4p+1 * log(T) +
             T * (0x0p+0 + T * (0x0p+0 + T * (0x0p+0 + T * 0x0p+0))) +
             0x1.0c6900093a3b6p+2);
    h[sNO] = 0x1.150d649310ce2p+8 *
             (T * (0x1.0dfb8404dcde9p+2 +
                   T * (-0x1.30051a34f94adp-9 +
                        T * (0x1.edf78adb7f8ccp-19 +
                             T * (-0x1.40c983f44a665p-29 +
                                  T * 0x1.3ba79af6c716ep-41)))) +
              0x1.33a4fbe76c8b4p+13);
    cp[sNO] =
        0x1.150d649310ce2p+8 *
        (0x1.0dfb8404dcde9p+2 +
         T * (-0x1.30051a34f94adp-8 +
              T * (0x1.7279a8249fa99p-17 +
                   T * (-0x1.40c983f44a665p-27 + T * 0x1.8a9181b478dc9p-39))));
    s[sNO] =
        0x1.150d649310ce2p+8 *
        (0x1.0dfb8404dcde9p+2 * log(T) +
         T * (-0x1.30051a34f94adp-8 +
              T * (0x1.7279a8249fa99p-18 +
                   T * (-0x1.abb75a9b0dddcp-29 + T * 0x1.8a9181b478dc9p-41))) +
         0x1.23f2c65b9983dp+1);
    h[sN2O] = 0x1.79c0ba18e3f82p+7 *
              (T * (0x1.20ea4c3994764p+1 +
                    T * (0x1.726eee866c268p-8 +
                         T * (-0x1.31d27d9c2c03bp-18 +
                              T * (0x1.4cab9a1c76912p-29 +
                                   T * -0x1.49f8336f5dcdfp-41)))) +
               0x1.112e31f8a0903p+13);
    cp[sN2O] =
        0x1.79c0ba18e3f82p+7 *
        (0x1.20ea4c3994764p+1 +
         T * (0x1.726eee866c268p-7 +
              T * (-0x1.cabbbc6a42059p-17 +
                   T * (0x1.4cab9a1c76912p-27 + T * -0x1.9c76404b35416p-39))));
    s[sN2O] =
        0x1.79c0ba18e3f82p+7 *
        (0x1.20ea4c3994764p+1 * log(T) +
         T * (0x1.726eee866c268p-7 +
              T * (-0x1.cabbbc6a42059p-18 +
                   T * (0x1.bb8f7825f36c3p-29 + T * -0x1.9c76404b35416p-41))) +
         0x1.584178705425fp+3);
    h[sNO2] = 0x1.696a1b0a8416fp+7 *
              (T * (0x1.f8d603ad33aa5p+1 +
                    T * (-0x1.9f9c56d20d983p-11 +
                         T * (0x1.74a102a6557fdp-18 +
                              T * (-0x1.5fc3e039d8b6cp-28 +
                                   T * 0x1.b913100af35e8p-40)))) +
               0x1.6a13c5d638866p+11);
    cp[sNO2] =
        0x1.696a1b0a8416fp+7 *
        (0x1.f8d603ad33aa5p+1 +
         T * (-0x1.9f9c56d20d983p-10 +
              T * (0x1.1778c1fcc01fep-16 +
                   T * (-0x1.5fc3e039d8b6cp-26 + T * 0x1.13abea06d81b1p-37))));
    s[sNO2] =
        0x1.696a1b0a8416fp+7 *
        (0x1.f8d603ad33aa5p+1 * log(T) +
         T * (-0x1.9f9c56d20d983p-10 +
              T * (0x1.1778c1fcc01fep-17 +
                   T * (-0x1.d5052af7cb9e5p-28 + T * 0x1.13abea06d81b1p-39))) +
         0x1.93f7ac0907e68p+2);
    h[sNH] = 0x1.14cfff745b811p+9 *
             (T * (0x1.bf17a02fb5d03p+1 +
                   T * (0x1.46f005b7511bfp-13 +
                        T * (-0x1.0a79c05240476p-21 +
                             T * (0x1.551317c9ecfb5p-31 +
                                  T * -0x1.d26fb323ff778p-43)))) +
              0x1.4731420c49ba6p+15);
    cp[sNH] =
        0x1.14cfff745b811p+9 *
        (0x1.bf17a02fb5d03p+1 +
         T * (0x1.46f005b7511bfp-12 +
              T * (-0x1.8fb6a07b606b1p-20 +
                   T * (0x1.551317c9ecfb5p-29 + T * -0x1.2385cff67faabp-40))));
    s[sNH] =
        0x1.14cfff745b811p+9 *
        (0x1.bf17a02fb5d03p+1 * log(T) +
         T * (0x1.46f005b7511bfp-12 +
              T * (-0x1.8fb6a07b606b1p-21 +
                   T * (0x1.c6c41fb7e6a47p-31 + T * -0x1.2385cff67faabp-42))) +
         0x1.d92c02bd49a21p+0);
    h[sHNO] = 0x1.0c0c830fdcdfcp+8 *
              (T * (0x1.2224b9f3ac34ap+2 +
                    T * (-0x1.739064067cd64p-9 +
                         T * (0x1.9d3d17dd76251p-18 +
                              T * (-0x1.2669bc89c46ap-28 +
                                   T * 0x1.382e76512416p-40)))) +
               0x1.68e2604189375p+13);
    cp[sHNO] =
        0x1.0c0c830fdcdfcp+8 *
        (0x1.2224b9f3ac34ap+2 +
         T * (-0x1.739064067cd64p-8 +
              T * (0x1.35edd1e6189bdp-16 +
                   T * (-0x1.2669bc89c46ap-26 + T * 0x1.863a13e56d1b7p-38))));
    s[sHNO] =
        0x1.0c0c830fdcdfcp+8 *
        (0x1.2224b9f3ac34ap+2 * log(T) +
         T * (-0x1.739064067cd64p-8 +
              T * (0x1.35edd1e6189bdp-17 +
                   T * (-0x1.888cfb625b38p-28 + T * 0x1.863a13e56d1b7p-40))) +
         0x1.bff5a02aad52bp+0);
    h[sNH2] = 0x1.0366cf1665468p+9 *
              (T * (0x1.0d0e622df2819p+2 +
                    T * (-0x1.140e47f4e9d5p-10 +
                         T * (0x1.3df453ff82b78p-19 +
                              T * (-0x1.819ee607a7abap-30 +
                                   T * 0x1.723647e9362c7p-42)))) +
               0x1.55f7a3d70a3d7p+14);
    cp[sNH2] =
        0x1.0366cf1665468p+9 *
        (0x1.0d0e622df2819p+2 +
         T * (-0x1.140e47f4e9d5p-9 +
              T * (0x1.dcee7dff44134p-18 +
                   T * (-0x1.819ee607a7abap-28 + T * 0x1.cec3d9e383b78p-40))));
    s[sNH2] =
        0x1.0366cf1665468p+9 *
        (0x1.0d0e622df2819p+2 * log(T) +
         T * (-0x1.140e47f4e9d5p-9 +
              T * (0x1.dcee7dff44134p-19 +
                   T * (-0x1.0114995a6fc7cp-29 + T * 0x1.cec3d9e383b78p-42))) -
         0x1.227e4f6644ad8p-3);
    h[sNNH] = 0x1.1e6cc1bcc36e4p+8 *
              (T * (0x1.160f71f86ae05p+2 +
                    T * (-0x1.3dd495d1b5803p-9 +
                         T * (0x1.c0b8f5e3cdf25p-18 +
                              T * (-0x1.7541ffcb2aa7ep-28 +
                                   T * 0x1.bf5facf3a76f5p-40)))) +
               0x1.c1dfe45a1cac1p+14);
    cp[sNNH] =
        0x1.1e6cc1bcc36e4p+8 *
        (0x1.160f71f86ae05p+2 +
         T * (-0x1.3dd495d1b5803p-8 +
              T * (0x1.508ab86ada75cp-16 +
                   T * (-0x1.7541ffcb2aa7ep-26 + T * 0x1.179bcc1848a59p-37))));
    s[sNNH] =
        0x1.1e6cc1bcc36e4p+8 *
        (0x1.160f71f86ae05p+2 * log(T) +
         T * (-0x1.3dd495d1b5803p-8 +
              T * (0x1.508ab86ada75cp-17 +
                   T * (-0x1.f1ad550ee38a8p-28 + T * 0x1.179bcc1848a59p-39))) +
         0x1.7d2d2bb23571dp+1);
    h[sCN] = 0x1.3f8958be799aep+8 *
             (T * (0x1.ce74a8488905dp+1 +
                   T * (-0x1.f4f6d1f6dfc83p-12 +
                        T * (0x1.7fbcad61cb559p-21 +
                             T * (-0x1.5a86901540efep-34 +
                                  T * -0x1.a2351c9ec9918p-44)))) +
              0x1.93f8ae147ae14p+15);
    cp[sCN] =
        0x1.3f8958be799aep+8 *
        (0x1.ce74a8488905dp+1 +
         T * (-0x1.f4f6d1f6dfc83p-11 +
              T * (0x1.1fcd820958803p-19 +
                   T * (-0x1.5a86901540efep-32 + T * -0x1.056131e33dfafp-41))));
    s[sCN] =
        0x1.3f8958be799aep+8 *
        (0x1.ce74a8488905dp+1 * log(T) +
         T * (-0x1.f4f6d1f6dfc83p-11 +
              T * (0x1.1fcd820958803p-20 +
                   T * (-0x1.ce08c01c56953p-34 + T * -0x1.056131e33dfafp-43))) +
         0x1.fd8101f31f46fp+1);
    h[sNCO] = 0x1.8bbb85aa76addp+7 *
              (T * (0x1.69d8de53070e1p+1 +
                    T * (0x1.20871c0410b29p-8 +
                         T * (-0x1.7735f28930a85p-19 +
                              T * (0x1.49f8561e273d6p-30 +
                                   T * -0x1.2bcba65f8f60ep-42)))) +
               0x1.cad3d0e560419p+13);
    cp[sNCO] =
        0x1.8bbb85aa76addp+7 *
        (0x1.69d8de53070e1p+1 +
         T * (0x1.20871c0410b29p-7 +
              T * (-0x1.196875e6e47e4p-17 +
                   T * (0x1.49f8561e273d6p-28 + T * -0x1.76be8ff773391p-40))));
    s[sNCO] =
        0x1.8bbb85aa76addp+7 *
        (0x1.69d8de53070e1p+1 * log(T) +
         T * (0x1.20871c0410b29p-7 +
              T * (-0x1.196875e6e47e4p-18 +
                   T * (0x1.b7f5c8283451dp-30 + T * -0x1.76be8ff773391p-42))) +
         0x1.319d67efd3621p+3);
    h[sHCN] = 0x1.339e97ece7b7p+8 *
              (T * (0x1.212689d784b6bp+1 +
                    T * (0x1.495b5337e06c1p-8 +
                         T * (-0x1.2aac838f4b80ap-18 +
                              T * (0x1.5ac53d990c77ep-29 +
                                   T * -0x1.52c5b7f52c4d1p-41)))) +
               0x1.cbc510624dd2fp+13);
    cp[sHCN] =
        0x1.339e97ece7b7p+8 *
        (0x1.212689d784b6bp+1 +
         T * (0x1.495b5337e06c1p-7 +
              T * (-0x1.c002c556f140fp-17 +
                   T * (0x1.5ac53d990c77ep-27 + T * -0x1.a77725f277605p-39))));
    s[sHCN] =
        0x1.339e97ece7b7p+8 *
        (0x1.212689d784b6bp+1 * log(T) +
         T * (0x1.495b5337e06c1p-7 +
              T * (-0x1.c002c556f140fp-18 +
                   T * (0x1.ce5c522165f53p-29 + T * -0x1.a77725f277605p-41))) +
         0x1.1d537df6a5e43p+3);
    h[sHOCN] = 0x1.82763b11595dp+7 *
               (T * (0x1.e49d454ab7df3p+1 +
                     T * (0x1.c3534e0727423p-9 +
                          T * (-0x1.1fa997b70c679p-20 +
                               T * (0x1.1c54d54eed6dep-33 +
                                    T * 0x1.5808b762828a2p-49)))) -
                0x1.615f7ced91687p+11);
    cp[sHOCN] =
        0x1.82763b11595dp+7 *
        (0x1.e49d454ab7df3p+1 +
         T * (0x1.c3534e0727423p-8 +
              T * (-0x1.af7e6392929b6p-19 +
                   T * (0x1.1c54d54eed6dep-31 + T * 0x1.ae0ae53b232cap-47))));
    s[sHOCN] =
        0x1.82763b11595dp+7 *
        (0x1.e49d454ab7df3p+1 * log(T) +
         T * (0x1.c3534e0727423p-8 +
              T * (-0x1.af7e6392929b6p-20 +
                   T * (0x1.7b1bc713e73d3p-33 + T * 0x1.ae0ae53b232cap-49))) +
         0x1.6881c9aeb534bp+2);
    h[sHNCO] = 0x1.82763b11595dp+7 *
               (T * (0x1.d0c366b210b3dp+1 +
                     T * (0x1.de990c66cfc7ap-9 +
                          T * (-0x1.981c83338d0a3p-21 +
                               T * (-0x1.6b89a95ebf73ep-33 +
                                    T * 0x1.4645e0ad41077p-44)))) -
                0x1.e71ae8a71de6ap+13);
    cp[sHNCO] =
        0x1.82763b11595dp+7 *
        (0x1.d0c366b210b3dp+1 +
         T * (0x1.de990c66cfc7ap-8 +
              T * (-0x1.32156266a9c7ap-19 +
                   T * (-0x1.6b89a95ebf73ep-31 + T * 0x1.97d758d891494p-42))));
    s[sHNCO] =
        0x1.82763b11595dp+7 *
        (0x1.d0c366b210b3dp+1 * log(T) +
         T * (0x1.de990c66cfc7ap-8 +
              T * (-0x1.32156266a9c7ap-20 +
                   T * (-0x1.e4b78c7e549a8p-33 + T * 0x1.97d758d891494p-44))) +
         0x1.8c73f438cc7a4p+2);
    h[sH2CN] = 0x1.288f36628b032p+8 *
               (T * (0x1.6d033a4723abp+1 +
                     T * (0x1.753e27e85841dp-9 +
                          T * (0x1.7f603d8bd0c55p-22 +
                               T * (-0x1.be052b31bb853p-32 +
                                    T * -0x1.a789b938fcdc2p-45)))) +
                0x1.bf7747ae147aep+14);
    cp[sH2CN] =
        0x1.288f36628b032p+8 *
        (0x1.6d033a4723abp+1 +
         T * (0x1.753e27e85841dp-8 +
              T * (0x1.1f882e28dc94p-20 +
                   T * (-0x1.be052b31bb853p-30 + T * -0x1.08b613c39e099p-42))));
    s[sH2CN] =
        0x1.288f36628b032p+8 *
        (0x1.6d033a4723abp+1 * log(T) +
         T * (0x1.753e27e85841dp-8 +
              T * (0x1.1f882e28dc94p-21 +
                   T * (-0x1.2958c7767d037p-31 + T * -0x1.08b613c39e099p-44))) +
         0x1.1fc49df4722d4p+3);
    h[sHCNN] = 0x1.9533b6ca1714dp+7 *
               (T * (0x1.431ce5e9d4449p+1 +
                     T * (0x1.057fb0284029dp-7 +
                          T * (-0x1.a4ea2aafe7ecp-18 +
                               T * (0x1.a0a160b67c399p-29 +
                                    T * -0x1.6c4fcad5a137p-41)))) +
                0x1.a7ebf7ced9168p+15);
    cp[sHCNN] =
        0x1.9533b6ca1714dp+7 *
        (0x1.431ce5e9d4449p+1 +
         T * (0x1.057fb0284029dp-6 +
              T * (-0x1.3bafa003edf1p-16 +
                   T * (0x1.a0a160b67c399p-27 + T * -0x1.c763bd8b0984cp-39))));
    s[sHCNN] =
        0x1.9533b6ca1714dp+7 *
        (0x1.431ce5e9d4449p+1 * log(T) +
         T * (0x1.057fb0284029dp-6 +
              T * (-0x1.3bafa003edf1p-17 +
                   T * (0x1.15c0eb2452d11p-28 + T * -0x1.c763bd8b0984cp-41))) +
         0x1.75a0ba1f4b1eep+3);
    h[sHCNO] = 0x1.82763b11595dp+7 *
               (T * (0x1.52da11437449p+1 +
                     T * (0x1.a1cf3bb2a0b62p-8 +
                          T * (-0x1.d4d76aef56155p-19 +
                               T * (0x1.2f59af8ea3ceep-30 +
                                    T * -0x1.55284761bf8f9p-43)))) +
                0x1.2d8c19ce075f7p+14);
    cp[sHCNO] =
        0x1.82763b11595dp+7 *
        (0x1.52da11437449p+1 +
         T * (0x1.a1cf3bb2a0b62p-7 +
              T * (-0x1.5fa19033809p-17 +
                   T * (0x1.2f59af8ea3ceep-28 + T * -0x1.aa72593a2f737p-41))));
    s[sHCNO] =
        0x1.82763b11595dp+7 *
        (0x1.52da11437449p+1 * log(T) +
         T * (0x1.a1cf3bb2a0b62p-7 +
              T * (-0x1.5fa19033809p-18 +
                   T * (0x1.947794be2fbe8p-30 + T * -0x1.aa72593a2f737p-43))) +
         0x1.57772bb087f2bp+3);
    h[sNH3] = 0x1.e81a453143a8ep+8 *
              (T * (0x1.124e45de30a26p+2 +
                    T * (-0x1.316e99de047ap-9 +
                         T * (0x1.e5d5bcc67c9d3p-18 +
                              T * (-0x1.87da8bbf9ca4bp-28 +
                                   T * 0x1.d135f9b4ce075p-40)))) -
               0x1.a55ba7ef9db23p+12);
    cp[sNH3] =
        0x1.e81a453143a8ep+8 *
        (0x1.124e45de30a26p+2 +
         T * (-0x1.316e99de047ap-8 +
              T * (0x1.6c604d94dd75ep-16 +
                   T * (-0x1.87da8bbf9ca4bp-26 + T * 0x1.22c1bc1100c49p-37))));
    s[sNH3] =
        0x1.e81a453143a8ep+8 *
        (0x1.124e45de30a26p+2 * log(T) +
         T * (-0x1.316e99de047ap-8 +
              T * (0x1.6c604d94dd75ep-17 +
                   T * (-0x1.053c5d2a686ddp-27 + T * 0x1.22c1bc1100c49p-39))) -
         0x1.4030dc15eaf8ep-1);
    h[sCH2CHO] = 0x1.825174a3f53f5p+7 *
                 (T * (0x1.b45c24c404a73p+1 +
                       T * (0x1.5fe1b0115dd4p-8 +
                            T * (0x1.527ee4c6ea8afp-21 +
                                 T * (-0x1.ebef1fbb3686fp-30 +
                                      T * 0x1.42d6bee6fe47p-41)))) +
                  0x1.7c5e809d49518p+10);
    cp[sCH2CHO] =
        0x1.825174a3f53f5p+7 *
        (0x1.b45c24c404a73p+1 +
         T * (0x1.5fe1b0115dd4p-7 +
              T * (0x1.fbbe572a5fd07p-20 +
                   T * (-0x1.ebef1fbb3686fp-28 + T * 0x1.938c6ea0bdd8cp-39))));
    s[sCH2CHO] =
        0x1.825174a3f53f5p+7 *
        (0x1.b45c24c404a73p+1 * log(T) +
         T * (0x1.5fe1b0115dd4p-7 +
              T * (0x1.fbbe572a5fd07p-21 +
                   T * (-0x1.47f4bfd22459fp-29 + T * 0x1.938c6ea0bdd8cp-41))) +
         0x1.31dd82fd75e2p+3);
    h[sCH3CHO] = 0x1.797a7ab03edafp+7 *
                 (T * (0x1.2eaf76e6106abp+2 +
                       T * (-0x1.a28ce427d2efep-10 +
                            T * (0x1.09d5a4c9ce534p-16 +
                                 T * (-0x1.ed90d262d7aa3p-27 +
                                      T * 0x1.34a728840b03p-38)))) -
                  0x1.511383126e979p+14);
    cp[sCH3CHO] =
        0x1.797a7ab03edafp+7 *
        (0x1.2eaf76e6106abp+2 +
         T * (-0x1.a28ce427d2efep-9 +
              T * (0x1.8ec0772eb57cep-15 +
                   T * (-0x1.ed90d262d7aa3p-25 + T * 0x1.81d0f2a50dc3bp-36))));
    s[sCH3CHO] =
        0x1.797a7ab03edafp+7 *
        (0x1.2eaf76e6106abp+2 * log(T) +
         T * (-0x1.a28ce427d2efep-9 +
              T * (0x1.8ec0772eb57cep-16 +
                   T * (-0x1.490b36ec8fc6dp-26 + T * 0x1.81d0f2a50dc3bp-38))) +
         0x1.0697d0005df3dp+2);
    h[sC3H8] = 0x1.791e6f137eed1p+7 *
               (T * (0x1.ddfac3d6032c6p-1 +
                     T * (0x1.b0f0b7a76578fp-7 +
                          T * (0x1.112d354968b8dp-19 +
                               T * (-0x1.79921013a4d49p-28 +
                                    T * 0x1.0bd243305f727p-39)))) -
                0x1.b43428f5c28f6p+13);
    cp[sC3H8] =
        0x1.791e6f137eed1p+7 *
        (0x1.ddfac3d6032c6p-1 +
         T * (0x1.b0f0b7a76578fp-6 +
              T * (0x1.99c3cfee1d154p-18 +
                   T * (-0x1.79921013a4d49p-26 + T * 0x1.4ec6d3fc774f1p-37))));
    s[sC3H8] =
        0x1.791e6f137eed1p+7 *
        (0x1.ddfac3d6032c6p-1 * log(T) +
         T * (0x1.b0f0b7a76578fp-6 +
              T * (0x1.99c3cfee1d154p-19 +
                   T * (-0x1.f76d6ac4dbc61p-28 + T * 0x1.4ec6d3fc774f1p-39))) +
         0x1.333a20578e5c5p+4);
    h[sC3H7] = 0x1.81f10d092f8b8p+7 *
               (T * (0x1.0d327faf0cc86p+0 +
                     T * (0x1.a9da4403baf59p-7 +
                          T * (0x1.a9ed608259079p-21 +
                               T * (-0x1.50e3cbe90afe8p-28 +
                                    T * 0x1.07d55cdad8895p-39)))) +
                0x1.4c3ee76c8b439p+13);
    cp[sC3H7] =
        0x1.81f10d092f8b8p+7 *
        (0x1.0d327faf0cc86p+0 +
         T * (0x1.a9da4403baf59p-6 +
              T * (0x1.3f720861c2c5bp-19 +
                   T * (-0x1.50e3cbe90afe8p-26 + T * 0x1.49cab4118eabap-37))));
    s[sC3H7] =
        0x1.81f10d092f8b8p+7 *
        (0x1.0d327faf0cc86p+0 * log(T) +
         T * (0x1.a9da4403baf59p-6 +
              T * (0x1.3f720861c2c5bp-20 +
                   T * (-0x1.c12fba8c0ea8bp-28 + T * 0x1.49cab4118eabap-39))) +
         0x1.51f6006d0d499p+4);
  } else {
    ComputeThermoData(h, cp, 300.0, s);
    for (i = 0; i < sEnd; i++) {
      h[i] = (T - 300.) * cp[i] + h[i];
      s[i] = log(T / 300) * cp[i] + s[i];
    }
  }
}

double MAX_C(double X1, double X2) { return ((X1 > X2) ? X1 : X2); }

} // namespace ideal_gas
} // namespace fub