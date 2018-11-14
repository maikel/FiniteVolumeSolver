#ifndef _TCdefsHSeen_
#define _TCdefsHSeen_

#ifdef __cplusplus
extern "C" {
#endif

/*! #
 * \file TC_defs.h
 * \brief Definitions of variables names used by the library.
 */
#include "TC_params.h"
#include "copyright.h"
#include <stdio.h>

/* --------------------------------------------
   numbers, lengths, etc.
   -------------------------------------------- */
/** \var static int TC_maxSpecInReac_
 *  \ingroup maxpar
 *  \brief Maximum number of species in a reaction */
extern int TC_maxSpecInReac_;
/** \var extern int TC_maxTbInReac_
 *  \ingroup maxpar
    \brief Max # of third-body efficiencies in a reaction */
extern int TC_maxTbInReac_;
/** \var extern int TC_nNASAinter_
 *  \ingroup maxpar
    \brief # of temperature regions for thermo fits */
extern int TC_nNASAinter_;
/** \var extern int TC_nCpCoef_
 *  \ingroup maxpar
    \brief # of polynomial coefficients for thermo fits */
extern int TC_nCpCoef_;
/** \var extern int TC_nArhPar_
 *  \ingroup maxpar
    \brief # of Arrhenius parameters */
extern int TC_nArhPar_;
/** \var extern int TC_nLtPar_
 *  \ingroup maxpar
    \brief # of parameters for Landau-Teller reactions */
extern int TC_nLtPar_;
/** \var extern int TC_nFallPar_
 *  \ingroup maxpar
    \brief # of parameters for pressure-dependent reactions */
extern int TC_nFallPar_;
/** \var extern int TC_nJanPar_
 *  \ingroup maxpar
    \brief # of parameters for Jannev-Langer fits (JAN) */
extern int TC_nJanPar_;
/** \var extern int TC_maxOrdPar_
 *  \ingroup maxpar
    \brief # of parameters for arbitrary order reactions */
extern int TC_maxOrdPar_;
/** \var extern int TC_nFit1Par_
 *  \ingroup maxpar
    \brief # of parameters for FIT1 fits */
extern int TC_nFit1Par_;

/** \var extern int TC_Nvars_
 *  \ingroup nospec
    \brief # of variables = no. of species + 1 */
extern int TC_Nvars_;
/** \var extern int TC_Nvjac_
 *  \ingroup nospec
    \brief # of lines/cols in the Jacobian = no. of species + 3 */
extern int TC_Nvjac_;
/** \var extern int TC_Nelem_
 *  \ingroup nospec
    \brief # of chemical elements */
extern int TC_Nelem_;
/** \var extern int TC_Nspec_
 *  \ingroup nospec
    \brief # of species */
extern int TC_Nspec_;
/** \var extern int TC_nIonSpec_
 *  \ingroup nospec
    \brief # of ion species */
extern int TC_nIonSpec_;
/** \var extern int TC_electrIndx_
 *  \ingroup nospec
    \brief Index of the electron species */
extern int TC_electrIndx_;
/** \var extern int TC_nIonEspec_
 *  \ingroup nospec
    \brief # of ion species excluding the electron species */
extern int TC_nIonEspec_;
/** \var extern int TC_nNASA9coef_
 *  \ingroup nospec
    \brief # of species with 9-term NASA polynomial fits */
extern int TC_nNASA9coef_;

/** \var extern int TC_Nreac_
 *  \ingroup noreac
    \brief # of reactions */
extern int TC_Nreac_;
/** \var extern int TC_nRevReac_
 *  \ingroup noreac
    \brief # of reactions with REV given */
extern int TC_nRevReac_;
/** \var extern int TC_nFallReac_
 *  \ingroup noreac
    \brief # of pressure-dependent reactions */
extern int TC_nFallReac_;
/** \var extern int TC_nPlogReac_
 *  \ingroup noreac
    \brief # of PLOG reactions */
extern int TC_nPlogReac_;
/** \var extern int TC_nThbReac_
 *  \ingroup noreac
    \brief # of reactions using third-body efficiencies */
extern int TC_nThbReac_;
/** \var extern int TC_nLtReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions */
extern int TC_nLtReac_;
/** \var extern int TC_nRltReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions with RLT given */
extern int TC_nRltReac_;
/** \var extern int TC_nHvReac_
 *  \ingroup noreac
    \brief # of reactions with HV */
extern int TC_nHvReac_;
/** \var extern int TC_nJanReac_
 *  \ingroup noreac
    \brief # of reactions with JAN fits */
extern int TC_nJanReac_;
/** \var extern int TC_nFit1Reac_
 *  \ingroup noreac
    \brief # of reactions with FIT1 fits */
extern int TC_nFit1Reac_;
/** \var extern int TC_nExciReac_
 *  \ingroup noreac
    \brief # of reactions with EXCI */
extern int TC_nExciReac_;
/** \var extern int TC_nMomeReac_
 *  \ingroup noreac
    \brief # of reactions with MOME */
extern int TC_nMomeReac_;
/** \var extern int TC_nXsmiReac_
 *  \ingroup noreac
    \brief # of reactions with XSMI */
extern int TC_nXsmiReac_;
/** \var extern int TC_nTdepReac_
 *  \ingroup noreac
    \brief # of reactions with TDEP */
extern int TC_nTdepReac_;
/** \var extern int TC_nRealNuReac_
 *  \ingroup noreac
    \brief # of reactions with non-int stoichiometric coefficients */
extern int TC_nRealNuReac_;
/** \var extern int TC_nOrdReac_
 *  \ingroup noreac
    \brief # of reactions with arbitrary order */
extern int TC_nOrdReac_;

/* --------------------------------------------
      elements & species -> names, masses, and element count
   -------------------------------------------- */
/** \var extern char TC_sNames_
    \brief species names, name of species i stored at LENGTHOFSPECNAME*i */
extern char* TC_sNames_;
/** \var extern char TC_eNames_
    \brief species names, name of species i stored at LENGTHOFELEMNAME*i */
extern char* TC_eNames_;
/** \var extern double TC_sMass_
    \brief array of species molar masses */
extern double* TC_sMass_;
/** \var extern double TC_eMass_
    \brief array of element molar masses */
extern double* TC_eMass_;
/** \var extern int TC_elemcount_
    \brief no. of atoms of element j in each species i at (i*TC_Nelem_+j)*/
extern int* TC_elemcount_;

/* --------------------------------------------
     species charges, number of temp fits, phase
   -------------------------------------------- */
/** \var extern int TC_sCharge_
    \brief species electrical charges */
extern int* TC_sCharge_;
/** \var extern int TC_sTfit_
    \brief no. of temperature fits for thermodynamic properties */
extern int* TC_sTfit_;
/** \var extern int TC_sPhase_
    \brief species phase id */
extern int* TC_sPhase_;

/* temperature limits for NASA polynomials and their coefficients */
extern double *TC_Tlo_, *TC_Tmi_, *TC_Thi_;
extern double* TC_cppol_;

/* temperature limits for NASA polynomials and their coefficients */
extern int *TC_spec9t_, *TC_spec9nrng_;
extern double *TC_spec9trng_, *TC_spec9coefs_;

/* list of non-electron ion species */
extern int* TC_sNion_;

/* Reaction data */
/* is reaction reversible ? */
extern int* TC_isRev_;

/* no. of reac+prod, no. of reac only, no. of prod only, stoichiom.coef.
 * indicator */
extern int *TC_reacNrp_, *TC_reacNreac_, *TC_reacNprod_, *TC_reacScoef_;

/* stoichiometric coeffs and reactants and product indices */
extern int *TC_reacNuki_, *TC_reacSidx_;
extern double* TC_reacNukiDbl_;

/* real stoichiometric coeffs */
extern int* TC_reacRnu_;
extern double* TC_reacRealNuki_;

/* list of reactions with given reverse Arrhenius parameters */
extern int* TC_reacRev_;

/* Arrhenius parameters, forward and reverse */
extern double *TC_reacArhenFor_, *TC_reacArhenRev_;

/* Pressure dependent reactions */
extern int *TC_reacPfal_, *TC_reacPtype_, *TC_reacPlohi_, *TC_reacPspec_;
extern double* TC_reacPpar_;

/* Third-body data */
extern int *TC_reacTbdy_, *TC_reacTbno_, *TC_specTbdIdx_;
extern double* TC_specTbdEff_;

/* Reactions with arbitrary orders  */
extern int *TC_reacAOrd_, *TC_specAOidx_;
extern double* TC_specAOval_;

/* Radiation wavelength data */
extern int* TC_reacHvIdx_;
extern double* TC_reacHvPar_;

/* PLOG reactions */
extern int* TC_reacPlogIdx_;
extern int* TC_reacPlogPno_;
extern double* TC_reacPlogPars_;

/* Other work arrays */
extern double *TC_kfor, *TC_krev, *TC_kforP, *TC_krevP, *TC_kforPder,
    *TC_krevPder, *TC_ropFor, *TC_ropRev, *TC_rop;
extern double *TC_cpks, *TC_hks, *TC_gk, *TC_gkp;
extern double *TC_PrDer, *TC_Crnd, *TC_CrndDer, *TC_dFfac, *TC_omg, *TC_omgP,
    *TC_jacFull;
extern double *TC_Xconc, *TC_scalIn, *TC_Mconc;
extern int *TC_sigNu, *TC_NuIJ;
extern double *TC_sigRealNu, *TC_RealNuIJ;
extern double* qfr;
extern double* TC_y2x2y;
extern double* TC_entr;

/* variables */
extern double TC_reltol, TC_abstol;
extern double TC_rhoref_, TC_pref_, TC_Tref_, TC_Wref_, TC_Daref_, TC_omgref_,
    TC_cpref_, TC_href_;
extern double TC_timref_;
extern int TC_isInit_, TC_tab_, TC_nonDim_, TC_RVset_;

/* tables */
extern int TC_Ntab_;
extern double TC_delT_, TC_odelT_;
extern double *TC_cptab, *TC_cpPtab, *TC_htab, *TC_gktab, *TC_gkPtab;
extern double *TC_kfortab, *TC_krevtab, *TC_kforPtab, *TC_krevPtab;

/* thermodynamic pressure */
extern double TC_pressure_, TC_prescgs_;

/* density */
extern double TC_rho_, TC_rhocgs_;
extern int TC_rhoset_;

/* universal gas constant */
extern double TC_Runiv_, TC_Rcal_, TC_Rcgs_;

/* Backup Arrhenius factors */
extern int TC_ArhenForChg_, TC_ArhenRevChg_;
extern double *TC_reacArhenForSave_, *TC_reacArhenRevSave_;

/* kc coefficients - depend on pressure only */
extern double* TC_kc_coeff;

/* Backups for reactional removal */
extern int TC_initRemoved;
extern int TC_NreacBackup_, TC_NrevBackup_, TC_NfalBackup_, TC_NthbBackup_;
extern int TC_nRealNuReacBackup_, TC_nOrdReacBackup_;
extern int *TC_isRevBackup_, *TC_reacNrpBackup_, *TC_reacNreacBackup_,
    *TC_reacNprodBackup_;
extern int *TC_reacNukiBackup_, *TC_reacSidxBackup_, *TC_reacScoefBackup_;
extern double *TC_reacArhenForBackup_, *TC_reacArhenRevBackup_,
    *TC_reacNukiDblBackup_;
extern int *TC_reacPfalBackup_, *TC_reacPtypeBackup_, *TC_reacPlohiBackup_,
    *TC_reacPspecBackup_;
extern int *TC_reacTbdyBackup_, *TC_reacTbnoBackup_, *TC_specTbdIdxBackup_;
extern double *TC_reacPparBackup_, *TC_specTbdEffBackup_;
extern int *TC_reacRnuBackup_, *TC_reacAOrdBackup_, *TC_specAOidxBackup_;
extern double *TC_reacRealNukiBackup_, *TC_sigRealNuBackup, *TC_RealNuIJBackup;
extern double *TC_specAOvalBackup_, *TC_kc_coeffBackup;
extern int *TC_sigNuBackup, *TC_NuIJBackup;

/* functions */
int TC_kmodint_(const char* mechfile, int* lmech, const char* thermofile,
                int* lthrm);
int TC_kmodint__(FILE* mechin, FILE* thermoin);

/*
             _____ ____     _       _ _
            |_   _/ ___|   (_)_ __ (_) |_   ___
              | || |       | | '_ \| | __| / __|
              | || |___    | | | | | | |_ | (__
              |_| \____|___|_|_| |_|_|\__(_)___|
                      |_____|


*/

int TC_makeSpace();

int TC_createTables(double delT);

void TC_errorMSG(int msgID, char const* func, int var1, int var2);
void TC_errorINI(FILE* errfile, char const* msg);

/*
        _____ ____          _   _ _
       |_   _/ ___|   _   _| |_(_) |___   ___
         | || |      | | | | __| | / __| / __|
         | || |___   | |_| | |_| | \__ \| (__
         |_| \____|___\__,_|\__|_|_|___(_)___|
                 |_____|
 */
double fastIntPow(double val, int exponent);

/*
          _____ ____
         |_   _/ ___|    _ __ _ __        ___
           | || |       | '__| '__|      / __|
           | || |___    | |  | |     _  | (__
           |_| \____|___|_|  |_|    (_)  \___|
                   |_____|
*/

/* molar rate of progress [kmol/(m3.s)] */
int TC_getRopsLocal(double* scal);

/* molar reaction rates [kmol/(m3.s)] */
int TC_getReacRates(double* scal, int Nvars, double* omega);

int TC_getgk(double t1, double t_1, double tln);
int TC_getgkFcn(double t1, double t_1, double tln);
int TC_getgkTab(double t1);
int TC_getgkFcn9t(double t, int icnt, double* gki);

int TC_getgkp(double t1, double t_1, double tln);
int TC_getgkpFcn(double t1, double t_1, double tln);
int TC_getgkpTab(double t1);
int TC_getgkpFcn9t(double t, int icnt, double* gki);

double TC_getSumNuGk(int i, double* gkLoc);
double TC_getSumRealNuGk(int i, int ir, double* gkLoc);

int TC_get3rdBdyConc(double* concX, double* concM);

int TC_getkForRev(double t1, double t_1, double tln);
int TC_getkForRevFcn(double t_1, double tln);
int TC_getkForRevTab(double t1);

int TC_getkForRevP(double t1, double t_1);
int TC_getkForRevPFcn(double t_1);
int TC_getkForRevPTab(double t1);

int TC_getRateofProg(double* concX);
int TC_getRateofProgDer(double* concX, int ireac, int ispec, double* qfr);
int TC_getCrnd(double t1, double t_1, double tln, double* concX, double* concM);
int TC_getCrndDer(int ireac, int* itbdy, int* ipfal, double t1, double t_1,
                  double tln, double* concX, double* concM);

/*
      _____ ____    _   _
     |_   _/ ___|  | |_| |__   ___ _ __ _ __ ___   ___         ___
       | || |      | __| '_ \ / _ \ '__| '_ ` _ \ / _ \       / __|
       | || |___   | |_| | | |  __/ |  | | | | | | (_) |  _  | (__
       |_| \____|___\__|_| |_|\___|_|  |_| |_| |_|\___/  (_)  \___|
               |_____|
*/
int TC_getCpSpecMsFcn(double t, double* cpi);
int TC_getCpSpecMsTab(double t1, double* cpi);
int TC_getCpSpecMlFcn(double t, double* cpi);
int TC_getCpSpecMlTab(double t1, double* cpi);

int TC_getCpSpecMs1Fcn(double t, int i, double* cpi);
int TC_getCpSpecMl1Fcn(double t, int i, double* cpi);

int TC_getCpMixMsP(double* scal, int Nvars, double* cpmix);
int TC_getCpSpecMsP(double t, int Nspec, double* cpi);

int TC_getCpSpecMsPFcn(double t, double* cpi);
int TC_getCpSpecMsPtab(double t1, double* cpi);

int TC_getCpSpecMs1PFcn(double t, int i, double* cpi);

int TC_getHspecMsFcn(double t, double* hi);
int TC_getHspecMsTab(double t1, double* hi);
int TC_getHspecMlFcn(double t, double* hi);
int TC_getHspecMlTab(double t1, double* hi);

int TC_getUspecMsFcn(double t, double* ui);
int TC_getUspecMsTab(double t1, double* ui);
int TC_getUspecMlFcn(double t, double* ui);
int TC_getUspecMlTab(double t1, double* ui);

int TC_getCpFcn9t(double t, int icnt, double* cpi);
int TC_getCpFcnP9t(double t, int icnt, double* cpi);
int TC_getHspecFcn9t(double t, int icnt, double* hi);

/* other defs */
#define TMAX 3500.0
#define TMIN 290.0
#define DTMID 10.0 /* buffer around middle temperature for NASA polynomials */

#ifdef __cplusplus
}
#endif

#endif