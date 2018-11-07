#include "TC_defs.h"
#include "TC_interface.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* --------------------------------------------
   numbers, lengths, etc.
   -------------------------------------------- */
/** \var static int TC_maxSpecInReac_
 *  \ingroup maxpar
 *  \brief Maximum number of species in a reaction */
int TC_maxSpecInReac_ ;
/** \var int TC_maxTbInReac_
 *  \ingroup maxpar
    \brief Max # of third-body efficiencies in a reaction */
int TC_maxTbInReac_   ;
/** \var int TC_nNASAinter_
 *  \ingroup maxpar
    \brief # of temperature regions for thermo fits */
int TC_nNASAinter_    ;
/** \var int TC_nCpCoef_
 *  \ingroup maxpar
    \brief # of polynomial coefficients for thermo fits */
int TC_nCpCoef_       ;
/** \var int TC_nArhPar_
 *  \ingroup maxpar
    \brief # of Arrhenius parameters */
int TC_nArhPar_       ;
/** \var int TC_nLtPar_
 *  \ingroup maxpar
    \brief # of parameters for Landau-Teller reactions */
int TC_nLtPar_        ;
/** \var int TC_nFallPar_
 *  \ingroup maxpar
    \brief # of parameters for pressure-dependent reactions */
int TC_nFallPar_      ;
/** \var int TC_nJanPar_
 *  \ingroup maxpar
    \brief # of parameters for Jannev-Langer fits (JAN) */
int TC_nJanPar_       ;
/** \var int TC_maxOrdPar_
 *  \ingroup maxpar
    \brief # of parameters for arbitrary order reactions */
int TC_maxOrdPar_     ;
/** \var int TC_nFit1Par_
 *  \ingroup maxpar
    \brief # of parameters for FIT1 fits */
int TC_nFit1Par_      ;

/** \var int TC_Nvars_
 *  \ingroup nospec
    \brief # of variables = no. of species + 1 */
int TC_Nvars_         ;
/** \var int TC_Nvjac_
 *  \ingroup nospec
    \brief # of lines/cols in the Jacobian = no. of species + 3 */
int TC_Nvjac_         ;
/** \var int TC_Nelem_
 *  \ingroup nospec
    \brief # of chemical elements */
int TC_Nelem_         ;
/** \var int TC_Nspec_
 *  \ingroup nospec
    \brief # of species */
int TC_Nspec_         ;
/** \var int TC_nIonSpec_
 *  \ingroup nospec
    \brief # of ion species */
int TC_nIonSpec_      ;
/** \var int TC_electrIndx_
 *  \ingroup nospec
    \brief Index of the electron species */
int TC_electrIndx_    ;
/** \var int TC_nIonEspec_
 *  \ingroup nospec
    \brief # of ion species excluding the electron species */
int TC_nIonEspec_     ;
/** \var int TC_nNASA9coef_
 *  \ingroup nospec
    \brief # of species with 9-term NASA polynomial fits */
int TC_nNASA9coef_    ;

/** \var int TC_Nreac_
 *  \ingroup noreac
    \brief # of reactions */
int TC_Nreac_         ;
/** \var int TC_nRevReac_
 *  \ingroup noreac
    \brief # of reactions with REV given */
int TC_nRevReac_      ;
/** \var int TC_nFallReac_
 *  \ingroup noreac
    \brief # of pressure-dependent reactions */
int TC_nFallReac_     ;
/** \var int TC_nPlogReac_
 *  \ingroup noreac
    \brief # of PLOG reactions */
int TC_nPlogReac_     ;
/** \var int TC_nThbReac_
 *  \ingroup noreac
    \brief # of reactions using third-body efficiencies */
int TC_nThbReac_      ;
/** \var int TC_nLtReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions */
int TC_nLtReac_       ;
/** \var int TC_nRltReac_
 *  \ingroup noreac
    \brief # of Landau-Teller reactions with RLT given */
int TC_nRltReac_      ;
/** \var int TC_nHvReac_
 *  \ingroup noreac
    \brief # of reactions with HV */
int TC_nHvReac_       ;
/** \var int TC_nJanReac_
 *  \ingroup noreac
    \brief # of reactions with JAN fits */
int TC_nJanReac_      ;
/** \var int TC_nFit1Reac_
 *  \ingroup noreac
    \brief # of reactions with FIT1 fits */
int TC_nFit1Reac_     ;
/** \var int TC_nExciReac_
 *  \ingroup noreac
    \brief # of reactions with EXCI */
int TC_nExciReac_     ;
/** \var int TC_nMomeReac_
 *  \ingroup noreac
    \brief # of reactions with MOME */
int TC_nMomeReac_     ;
/** \var int TC_nXsmiReac_
 *  \ingroup noreac
    \brief # of reactions with XSMI */
int TC_nXsmiReac_     ;
/** \var int TC_nTdepReac_
 *  \ingroup noreac
    \brief # of reactions with TDEP */
int TC_nTdepReac_     ;
/** \var int TC_nRealNuReac_
 *  \ingroup noreac
    \brief # of reactions with non-int stoichiometric coefficients */
int TC_nRealNuReac_   ;
/** \var int TC_nOrdReac_
 *  \ingroup noreac
    \brief # of reactions with arbitrary order */
int TC_nOrdReac_      ;

/* --------------------------------------------
      elements & species -> names, masses, and element count
   -------------------------------------------- */
/** \var char TC_sNames_
    \brief species names, name of species i stored at LENGTHOFSPECNAME*i */
char   *TC_sNames_    ;
/** \var char TC_eNames_
    \brief species names, name of species i stored at LENGTHOFELEMNAME*i */
char   *TC_eNames_    ;
/** \var double TC_sMass_
    \brief array of species molar masses */
double *TC_sMass_     ;
/** \var double TC_eMass_
    \brief array of element molar masses */
double *TC_eMass_     ;
/** \var int TC_elemcount_
    \brief no. of atoms of element j in each species i at (i*TC_Nelem_+j)*/
int  *TC_elemcount_ ;

/* --------------------------------------------
     species charges, number of temp fits, phase
   -------------------------------------------- */
/** \var int TC_sCharge_
    \brief species electrical charges */
int *TC_sCharge_ ;
/** \var int TC_sTfit_
    \brief no. of temperature fits for thermodynamic properties */
int *TC_sTfit_   ;
/** \var int TC_sPhase_
    \brief species phase id */
int *TC_sPhase_  ;

/* temperature limits for NASA polynomials and their coefficients */
double *TC_Tlo_,*TC_Tmi_,*TC_Thi_ ;
double *TC_cppol_ ;

/* temperature limits for NASA polynomials and their coefficients */
int *TC_spec9t_, *TC_spec9nrng_ ;
double *TC_spec9trng_, *TC_spec9coefs_ ;

/* list of non-electron ion species */
int *TC_sNion_ ;

/* Reaction data */
/* is reaction reversible ? */
int *TC_isRev_ ;

/* no. of reac+prod, no. of reac only, no. of prod only, stoichiom.coef. indicator */
int *TC_reacNrp_, *TC_reacNreac_, *TC_reacNprod_, *TC_reacScoef_;

/* stoichiometric coeffs and reactants and product indices */
int *TC_reacNuki_, *TC_reacSidx_ ;
double *TC_reacNukiDbl_ ;

/* real stoichiometric coeffs */
int    *TC_reacRnu_ ;
double *TC_reacRealNuki_ ;

/* list of reactions with given reverse Arrhenius parameters */
int *TC_reacRev_ ;

/* Arrhenius parameters, forward and reverse */
double *TC_reacArhenFor_, *TC_reacArhenRev_ ;

/* Pressure dependent reactions */
int *TC_reacPfal_, *TC_reacPtype_, *TC_reacPlohi_, *TC_reacPspec_ ;
double *TC_reacPpar_ ;

/* Third-body data */
int *TC_reacTbdy_, *TC_reacTbno_, *TC_specTbdIdx_ ;
double *TC_specTbdEff_ ;

/* Reactions with arbitrary orders  */
int    *TC_reacAOrd_, *TC_specAOidx_ ;
double *TC_specAOval_ ;

/* Radiation wavelength data */
int    *TC_reacHvIdx_ ;
double *TC_reacHvPar_ ;

/* PLOG reactions */
int    *TC_reacPlogIdx_  ;
int    *TC_reacPlogPno_  ;
double *TC_reacPlogPars_ ;

/* Other work arrays */
double *TC_kfor, *TC_krev, *TC_kforP, *TC_krevP, *TC_kforPder, *TC_krevPder, *TC_ropFor, *TC_ropRev, *TC_rop;
double *TC_cpks, *TC_hks, *TC_gk, *TC_gkp ;
double *TC_PrDer, *TC_Crnd, *TC_CrndDer, *TC_dFfac, *TC_omg, *TC_omgP, *TC_jacFull ;
double *TC_Xconc,*TC_scalIn, *TC_Mconc ;
int    *TC_sigNu, *TC_NuIJ ;
double *TC_sigRealNu, *TC_RealNuIJ ;
double *qfr ;
double *TC_y2x2y ;
double *TC_entr  ;

/* variables */
double TC_reltol,TC_abstol ;
double TC_rhoref_, TC_pref_, TC_Tref_, TC_Wref_, TC_Daref_, TC_omgref_, TC_cpref_, TC_href_ ;
double TC_timref_ ;
int    TC_isInit_, TC_tab_, TC_nonDim_, TC_RVset_ ;

/* tables */
int    TC_Ntab_ ;
double TC_delT_, TC_odelT_ ;
double *TC_cptab,   *TC_cpPtab,  *TC_htab,     *TC_gktab, *TC_gkPtab ;
double *TC_kfortab, *TC_krevtab, *TC_kforPtab, *TC_krevPtab ;

/* thermodynamic pressure */
double TC_pressure_, TC_prescgs_ ;

/* density */
double TC_rho_, TC_rhocgs_ ;
int TC_rhoset_ ;

/* universal gas constant */
double TC_Runiv_, TC_Rcal_, TC_Rcgs_ ;

/* Backup Arrhenius factors */
int  TC_ArhenForChg_, TC_ArhenRevChg_ ;
double *TC_reacArhenForSave_, *TC_reacArhenRevSave_ ;

/* kc coefficients - depend on pressure only */
double *TC_kc_coeff ;

/* Backups for reactional removal */
int TC_initRemoved ;
int TC_NreacBackup_, TC_NrevBackup_, TC_NfalBackup_, TC_NthbBackup_ ;
int TC_nRealNuReacBackup_, TC_nOrdReacBackup_ ;
int *TC_isRevBackup_, *TC_reacNrpBackup_, *TC_reacNreacBackup_, *TC_reacNprodBackup_ ;
int *TC_reacNukiBackup_, *TC_reacSidxBackup_, *TC_reacScoefBackup_ ;
double *TC_reacArhenForBackup_, *TC_reacArhenRevBackup_, *TC_reacNukiDblBackup_ ;
int *TC_reacPfalBackup_, *TC_reacPtypeBackup_, *TC_reacPlohiBackup_, *TC_reacPspecBackup_ ;
int *TC_reacTbdyBackup_, *TC_reacTbnoBackup_, *TC_specTbdIdxBackup_ ;
double *TC_reacPparBackup_, *TC_specTbdEffBackup_;
int *TC_reacRnuBackup_, *TC_reacAOrdBackup_, *TC_specAOidxBackup_ ;
double *TC_reacRealNukiBackup_, *TC_sigRealNuBackup, *TC_RealNuIJBackup;
double *TC_specAOvalBackup_, *TC_kc_coeffBackup ;
int *TC_sigNuBackup, *TC_NuIJBackup ;

/*!
 * \file TC_init.c
 * \brief Initialize chemical library
*/
/**
 * \ingroup init
 * \brief Initialize library
*/
char data_path[256];

int icomp (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
int TC_initChem(const char *mechfile, const char *thermofile, int tab, double delT)
{
  /**
    \param mechfile : name of file containing kinetic model in chemkin format
    \param thermofile : name of file containing coefficients for NASA polynomials
    \param tab : flag to tabulate temperature dependent terms (tab=1 -> create tables)
    \param delT : temperature step size [K] for the tables; used only tab=1
  */

  /* loop variables */
  int i, j, k, ir, ij ;

  int  lmech, lthrm, ierr, itemp ;
  double sml, reacbalance ;
  char charvar4[9] ;

  FILE *chemfile, *echofile, *errfile ;

  /* zero-out variables */
  TC_isInit_   = 0 ;
  TC_pressure_ = 0.0 ;
  TC_prescgs_  = 0.0 ;
  TC_rho_      = 0.0 ;
  TC_rhocgs_   = 0.0 ;
  TC_rhoset_   = 0   ;
  TC_Runiv_    = 0.0 ;
  TC_Rcal_     = 0.0 ;
  TC_Rcgs_     = 0.0 ;

  /* integers */
  TC_Nvars_ = TC_Nvjac_ = 0 ;
  TC_maxSpecInReac_ = TC_maxTbInReac_ = TC_nNASAinter_ = TC_nCpCoef_ = TC_nArhPar_ = TC_nLtPar_ = 0 ;
  TC_nFallPar_ = TC_nJanPar_ = TC_maxOrdPar_ = TC_nFit1Par_  = 0;
  TC_Nelem_ = TC_Nspec_ = TC_Nreac_ = TC_nRevReac_ = TC_nFallReac_ = TC_nPlogReac_ = TC_nThbReac_ = 0 ;
  TC_nLtReac_ = TC_nRltReac_ = TC_nHvReac_  = TC_nIonSpec_ = TC_nJanReac_ = TC_nFit1Reac_ = 0 ;
  TC_nExciReac_  = TC_nMomeReac_ = TC_nXsmiReac_ = TC_nTdepReac_ = TC_nRealNuReac_ = TC_nOrdReac_ = 0 ;
  TC_electrIndx_ = TC_nIonEspec_ = TC_nNASA9coef_ = 0 ;

  /* elements & species -> names, masses, and element count */
  TC_sNames_    = 0 ;
  TC_eNames_    = 0 ;
  TC_sMass_     = 0 ;
  TC_eMass_     = 0 ;
  TC_elemcount_ = 0 ;

  /* species charges, number of temp fits, phase */
  TC_sCharge_ = 0 ;
  TC_sTfit_   = 0 ;
  TC_sPhase_  = 0 ;

  /* temperature limits for NASA polynomials and their coefficients */
  TC_Tlo_   = 0 ;
  TC_Tmi_   = 0 ;
  TC_Thi_   = 0 ;
  TC_cppol_ = 0 ;

  /* list of non-electron ion species */
  TC_sNion_ = 0 ;

  /* is reaction reversible ? */
  TC_isRev_ = 0 ;

  /* no. of reac+prod, no. of reac only, no. of prod only, stoichiom.coef. indicator */
  TC_reacNrp_   = 0 ;
  TC_reacNreac_ = 0 ;
  TC_reacNprod_ = 0 ;
  TC_reacScoef_ = 0 ;

  /* stoichiometric coeffs and reactants and product indices */
  TC_reacNuki_ = 0 ;
  TC_reacSidx_ = 0 ;

  /* real stoichiometric coeffs */
  TC_reacRnu_      = 0 ;
  TC_reacRealNuki_ = 0 ;

  /* list of reactions with given reverse Arrhenius parameters */
  TC_reacRev_ = 0 ;

  /* Arrhenius parameters, forward and reverse */
  TC_reacArhenFor_ = 0 ;
  TC_reacArhenRev_ = 0 ;

  /* Pressure dependent reactions */
  TC_reacPfal_  = 0 ;
  TC_reacPtype_ = 0 ;
  TC_reacPlohi_ = 0 ;
  TC_reacPspec_ = 0 ;
  TC_reacPpar_  = 0 ;

  /* Third-body data */
  TC_reacTbdy_ = 0 ;
  TC_reacTbno_ = 0 ;
  TC_specTbdIdx_ = 0 ;
  TC_specTbdEff_ = 0 ;

  /* Reactions with arbitrary orders */
  TC_reacAOrd_  = 0 ;
  TC_specAOidx_ = 0 ;
  TC_specAOval_ = 0 ;

  /* Other work arrays */
  TC_kfor   = 0 ;
  TC_krev   = 0 ;
  TC_kforP  = 0 ;
  TC_krevP  = 0 ;
  TC_ropFor = 0 ;
  TC_ropRev = 0 ;
  TC_rop    = 0 ;

  TC_cpks = 0 ;
  TC_hks  = 0 ;
  TC_gk   = 0 ;
  TC_gkp  = 0 ;

  TC_PrDer = 0 ;
  TC_Crnd  = 0 ;
  TC_CrndDer = 0 ;
  TC_dFfac   = 0 ;
  TC_omg     = 0 ;
  TC_omgP    = 0 ;
  TC_jacFull = 0 ;
  TC_Xconc   = 0 ;
  TC_scalIn  = 0 ;
  TC_Mconc   = 0 ;
  TC_sigNu   = 0 ;
  TC_NuIJ    = 0 ;
  TC_sigRealNu = 0 ;
  TC_RealNuIJ  = 0 ;

  TC_kc_coeff = 0;

  /* work variables */
  TC_reltol  = 0.0 ;
  TC_abstol  = 0.0 ;
  TC_rhoref_ = 0.0 ;
  TC_pref_   = 0.0 ;
  TC_Tref_   = 0.0 ;
  TC_Wref_   = 0.0 ;
  TC_Daref_  = 0.0 ;
  TC_omgref_ = 0.0 ;
  TC_cpref_  = 0.0 ;
  TC_href_   = 0.0 ;
  TC_timref_ = 0.0 ;
  TC_isInit_ = 0 ;
  TC_tab_    = 0 ;
  TC_nonDim_ = 0 ;
  TC_RVset_  = 0 ;

  /* tables */
  TC_Ntab_  = 0   ;
  TC_delT_  = 0.0 ;
  TC_odelT_ = 0.0 ;

  TC_cptab    = 0 ;
  TC_cpPtab   = 0 ;
  TC_htab     = 0 ;
  TC_gktab    = 0 ;
  TC_gkPtab   = 0 ;
  TC_kfortab  = 0 ;
  TC_krevtab  = 0 ;
  TC_kforPtab = 0 ;
  TC_krevPtab = 0 ;

  /* Backup Arrhenius factors */
  TC_ArhenForChg_       = 0 ;
  TC_reacArhenForSave_  = 0 ;
  TC_ArhenRevChg_       = 0 ;
  TC_reacArhenRevSave_  = 0 ;

  /* Other backup reaction parameters for use by TC_removeReaction */
  TC_initRemoved = 0 ;
  TC_NreacBackup_ = TC_NrevBackup_ = TC_NfalBackup_ = TC_NthbBackup_ = 0 ;
  TC_nRealNuReacBackup_  = TC_nOrdReacBackup_ = 0 ;

  TC_isRevBackup_        = 0;
  TC_reacNrpBackup_      = 0;
  TC_reacNreacBackup_    = 0;
  TC_reacNprodBackup_    = 0;
  TC_reacArhenForBackup_ = 0;
  TC_reacArhenRevBackup_ = 0;
  TC_reacNukiBackup_     = 0;
  TC_reacSidxBackup_     = 0;
  TC_reacScoefBackup_    = 0;
  TC_reacNukiDblBackup_  = 0;
  TC_reacPfalBackup_     = 0;
  TC_reacPtypeBackup_    = 0;
  TC_reacPlohiBackup_    = 0;
  TC_reacPspecBackup_    = 0;
  TC_reacTbdyBackup_     = 0;
  TC_reacTbnoBackup_     = 0;
  TC_specTbdIdxBackup_   = 0;
  TC_reacPparBackup_     = 0;
  TC_specTbdEffBackup_   = 0;
  TC_reacRnuBackup_      = 0;
  TC_reacAOrdBackup_     = 0;
  TC_specAOidxBackup_    = 0;
  TC_reacRealNukiBackup_ = 0;
  TC_sigRealNuBackup     = 0;
  TC_RealNuIJBackup      = 0;
  TC_specAOvalBackup_    = 0;
  TC_kc_coeffBackup      = 0;
  TC_sigNuBackup         = 0;
  TC_NuIJBackup          = 0;

  /* determine machine precision parameters (for numerical Jac) */
  sml = 1.0 ;
  while ( 1.0+sml != 1.0 ) sml *= 0.5 ;
  TC_reltol = sqrt(2.0*sml) ;
  TC_abstol = TC_reltol ;

#ifdef VERBOSE
  printf("TC_initChem() : \n") ;
  printf("  -> Machine epsilon    = %e\n",sml   ) ;
  printf("  -> Absolute tolerance = %e\n",TC_abstol) ;
  printf("  -> Relative tolerance = %e\n",TC_reltol) ;
#endif

  /* Set kinetic model input and thermodynamic files */
  lmech = strlen( mechfile   ) ;
  lthrm = strlen( thermofile ) ;

  /* call the interpreter */
  ierr = TC_kmodint_(mechfile,&lmech,thermofile,&lthrm) ;
  if ( ierr != 0 )
  {
    printf("TC_initChem() : received error from the kinetic model interpreter !!!\n") ;
    exit(ierr) ;
  }

  /* Retrieve things from kmod.list */
  char fn_chem[256] = {0}, fn_echo[256] = {0}, fn_erro[256] = {0};
  strcat(fn_chem, "kmod.list");
  strcat(fn_echo, "kmod.echo");
  strcat(fn_erro, "kmod.err" );
  chemfile = fopen(fn_chem, "r") ;
  echofile = fopen(fn_echo, "w") ;
  errfile  = fopen(fn_erro, "w") ;
  //chemfile = fopen("kmod.list","r") ;
  //echofile = fopen("kmod.echo","w") ;
  //errfile  = fopen("kmod.err" ,"w") ;

  fscanf(chemfile,"%s",charvar4);
  /* printf("%s\n",charvar4) ; */
  if ( strcmp(charvar4,"ERRO") == 0 )
    TC_errorINI( errfile, "kmod.list : Error when interpreting kinetic model  !!!" ) ;

  fscanf(chemfile,"%d",&TC_maxSpecInReac_ ) ;
  fscanf(chemfile,"%d",&TC_maxTbInReac_   ) ;
  fscanf(chemfile,"%d",&TC_nNASAinter_    ) ;
  fscanf(chemfile,"%d",&TC_nCpCoef_       ) ;
  fscanf(chemfile,"%d",&TC_nArhPar_       ) ;
  fscanf(chemfile,"%d",&TC_nLtPar_        ) ;
  fscanf(chemfile,"%d",&TC_nFallPar_      ) ;
  fscanf(chemfile,"%d",&TC_nJanPar_       ) ;
  fscanf(chemfile,"%d",&TC_maxOrdPar_     ) ;
  fscanf(chemfile,"%d",&TC_nFit1Par_      ) ;

  fscanf(chemfile,"%d",&TC_Nelem_         ) ;
  fscanf(chemfile,"%d",&TC_Nspec_         ) ;
  fscanf(chemfile,"%d",&TC_Nreac_         ) ;
  fscanf(chemfile,"%d",&TC_nRevReac_      ) ;
  fscanf(chemfile,"%d",&TC_nFallReac_     ) ;
  fscanf(chemfile,"%d",&TC_nPlogReac_     ) ;
  fscanf(chemfile,"%d",&TC_nThbReac_      ) ;
  fscanf(chemfile,"%d",&TC_nLtReac_       ) ;
  fscanf(chemfile,"%d",&TC_nRltReac_      ) ;
  fscanf(chemfile,"%d",&TC_nHvReac_       ) ;
  fscanf(chemfile,"%d",&TC_nIonSpec_      ) ;
  fscanf(chemfile,"%d",&TC_nJanReac_      ) ;
  fscanf(chemfile,"%d",&TC_nFit1Reac_     ) ;
  fscanf(chemfile,"%d",&TC_nExciReac_     ) ;
  fscanf(chemfile,"%d",&TC_nMomeReac_     ) ;
  fscanf(chemfile,"%d",&TC_nXsmiReac_     ) ;
  fscanf(chemfile,"%d",&TC_nTdepReac_     ) ;
  fscanf(chemfile,"%d",&TC_nRealNuReac_   ) ;
  fscanf(chemfile,"%d",&TC_nOrdReac_      ) ;
  fscanf(chemfile,"%d",&TC_electrIndx_    ) ;
  fscanf(chemfile,"%d",&TC_nIonEspec_     ) ;
  fscanf(chemfile,"%d",&TC_nNASA9coef_    ) ;

  fprintf(echofile,"kmod.list : Max # of species in a reaction                    : %d\n",TC_maxSpecInReac_ ) ;
  fprintf(echofile,"kmod.list : Max # of third-body efficiencies in a reaction    : %d\n",TC_maxTbInReac_   ) ;
  fprintf(echofile,"kmod.list : # of temperature regions for thermo fits          : %d\n",TC_nNASAinter_    ) ;
  fprintf(echofile,"kmod.list : # of polynomial coefficients for thermo fits      : %d\n",TC_nCpCoef_       ) ;
  fprintf(echofile,"kmod.list : # of species with 9 coefficients thermo props     : %d\n",TC_nNASA9coef_    ) ;
  fprintf(echofile,"kmod.list : # of Arrhenius parameters                         : %d\n",TC_nArhPar_   ) ;
  fprintf(echofile,"kmod.list : # of parameters for Landau-Teller reactions       : %d\n",TC_nLtPar_    ) ;
  fprintf(echofile,"kmod.list : # of parameters for pressure-dependent reactions  : %d\n",TC_nFallPar_  ) ;
  fprintf(echofile,"kmod.list : # of parameters for Jannev-Langer fits (JAN)      : %d\n",TC_nJanPar_   ) ;
  fprintf(echofile,"kmod.list : # of parameters for arbitrary order reactions     : %d\n",TC_maxOrdPar_ ) ;
  fprintf(echofile,"kmod.list : # of parameters for FIT1 fits                     : %d\n",TC_nFit1Par_  ) ;
  fprintf(echofile,"kmod.list : # of elements                                     : %d\n",TC_Nelem_  ) ;
  fprintf(echofile,"kmod.list : # of species                                      : %d\n",TC_Nspec_  ) ;
  fprintf(echofile,"kmod.list : # of reactions                                    : %d\n",TC_Nreac_  ) ;
  fprintf(echofile,"kmod.list : # of reactions with REV given                     : %d\n",TC_nRevReac_    ) ;
  fprintf(echofile,"kmod.list : # of pressure-dependent reactions                 : %d\n",TC_nFallReac_   ) ;
  fprintf(echofile,"kmod.list : # of PLOG reactions                               : %d\n",TC_nPlogReac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions using third-body efficiencies      : %d\n",TC_nThbReac_    ) ;
  fprintf(echofile,"kmod.list : # of Landau-Teller reactions                      : %d\n",TC_nLtReac_     ) ;
  fprintf(echofile,"kmod.list : # of Landau-Teller reactions with RLT given       : %d\n",TC_nRltReac_    ) ;
  fprintf(echofile,"kmod.list : # of reactions with HV                            : %d\n",TC_nHvReac_     ) ;
  fprintf(echofile,"kmod.list : # of ion species                                  : %d\n",TC_nIonSpec_    ) ;
  fprintf(echofile,"kmod.list : # of reactions with JAN fits                      : %d\n",TC_nJanReac_    ) ;
  fprintf(echofile,"kmod.list : # of reactions with FIT1 fits                     : %d\n",TC_nFit1Reac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions with EXCI                          : %d\n",TC_nExciReac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions with MOME                          : %d\n",TC_nMomeReac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions with XSMI                          : %d\n",TC_nXsmiReac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions with TDEP                          : %d\n",TC_nTdepReac_   ) ;
  fprintf(echofile,"kmod.list : # of reactions with non-int stoichiometric coeffs : %d\n",TC_nRealNuReac_ ) ;
  fprintf(echofile,"kmod.list : # of reactions with arbitrary order               : %d\n",TC_nOrdReac_    ) ;
  fprintf(echofile,"kmod.list : Index of the electron species                     : %d\n",TC_electrIndx_ ) ;
  fprintf(echofile,"kmod.list : # of ion species excluding the electron species   : %d\n",TC_nIonEspec_  ) ;
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;
  fflush(echofile) ;

  fscanf(chemfile,"%lf", &reacbalance ) ;
  fprintf(echofile,"kmod.list : Tolerance for reaction balance                    : %e\n",reacbalance ) ;
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;
  fflush(echofile) ;

  /* Elements' name and weights */
  TC_eNames_ = (char   *) malloc ( TC_Nelem_*LENGTHOFELEMNAME*sizeof(char  ) ) ;
  TC_eMass_  = (double *) malloc ( TC_Nelem_*sizeof(double) ) ;

  for ( i=0 ; i<TC_Nelem_*LENGTHOFELEMNAME ; i++ ) TC_eNames_[i] = 0 ;
  for ( i=0 ; i<TC_Nelem_                  ; i++ ) TC_eMass_ [i] = 0 ;

  for ( i=0 ; i<TC_Nelem_ ; i++ )
    fscanf(chemfile,"%s", &(TC_eNames_[i*LENGTHOFELEMNAME]) ) ;
  for ( i=0 ; i<TC_Nelem_ ; i++ )
    fscanf(chemfile,"%lf", &(TC_eMass_[i]) ) ;

  fprintf( echofile, "No. \t Element \t Mass\n" ) ;
  for ( i=0 ; i<TC_Nelem_ ; i++ )
    fprintf( echofile, "%-3d\t%-4s\t%f10.7\n",
	     i+1,&(TC_eNames_[i*LENGTHOFELEMNAME]),TC_eMass_[i]) ;
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;
  fflush(echofile) ;

  /* Species' name and weights */
  TC_sNames_ = (char   *) malloc ( TC_Nspec_*LENGTHOFSPECNAME*sizeof(char) ) ;
  TC_sMass_  = (double *) malloc ( TC_Nspec_*sizeof(double) ) ;

  for ( i=0 ; i<TC_Nspec_*LENGTHOFSPECNAME ; i++ ) TC_sNames_[i] = 0 ;
  for ( i=0 ; i<TC_Nspec_                  ; i++ ) TC_sMass_ [i] = 0 ;

  for ( i=0 ; i<TC_Nspec_ ; i++ )
    fscanf(chemfile,"%s", &(TC_sNames_[i*LENGTHOFSPECNAME]) ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ )
    fscanf(chemfile,"%lf", &(TC_sMass_[i]) ) ;

  fprintf( echofile, "No. \t Species \t Mass\n" ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ )
    fprintf( echofile, "%-3d\t%-32s\t%12.7f\n",
	     i+1,&(TC_sNames_[i*LENGTHOFSPECNAME]),TC_sMass_[i]) ;
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;
  fflush(echofile) ;

  /* Species' elemental composition */
  TC_elemcount_ = (int *) malloc ( TC_Nspec_*TC_Nelem_ * sizeof(int) ) ;
  for ( i=0 ; i<TC_Nspec_*TC_Nelem_ ; i++ )
    fscanf(chemfile,"%d", &(TC_elemcount_[i]) ) ;

  fprintf( echofile, "Elemental composition of species\n" ) ;
  fprintf( echofile, "No. \t Species \t\t Element\n\t\t\t\t" ) ;
  for ( i=0 ; i<TC_Nelem_ ; i++ )
    fprintf( echofile, "%s\t",&(TC_eNames_[i*LENGTHOFELEMNAME])) ;
  fprintf(echofile,"\n") ;

  for ( i=0 ; i<TC_Nspec_ ; i++ )
  {
    fprintf( echofile, "%-3d\t%-32s",
	     i+1,&(TC_sNames_[i*LENGTHOFSPECNAME])) ;
    for ( j=0 ; j < TC_Nelem_ ; j++ )
      fprintf( echofile, "%-3d\t",TC_elemcount_[i*TC_Nelem_+j]) ;
    fprintf(echofile,"\n") ;
  }

  /* Species charges, no of tempfits, phase */
  TC_sCharge_ = (int *) malloc ( TC_Nspec_ * sizeof(int) ) ;
  TC_sTfit_   = (int *) malloc ( TC_Nspec_ * sizeof(int) ) ;
  TC_sPhase_  = (int *) malloc ( TC_Nspec_ * sizeof(int) ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ ) fscanf(chemfile,"%d", &(TC_sCharge_[i] ) ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ ) fscanf(chemfile,"%d", &(TC_sTfit_  [i] ) ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ ) fscanf(chemfile,"%d", &(TC_sPhase_ [i] ) ) ;

  /* Range of temperatures for thermo fits */
  TC_Tlo_ = (double *) malloc ( TC_Nspec_ * sizeof(double) ) ;
  TC_Tmi_ = (double *) malloc ( TC_Nspec_ * sizeof(double) ) ;
  TC_Thi_ = (double *) malloc ( TC_Nspec_ * sizeof(double) ) ;
  for ( i=0 ; i< TC_Nspec_ ;i++ )
  {
    fscanf(chemfile,"%lf", &(TC_Tlo_[i] ) ) ;
    fscanf(chemfile,"%lf", &(TC_Tmi_[i] ) ) ;
    fscanf(chemfile,"%lf", &(TC_Thi_[i] ) ) ;
  }

  fprintf( echofile, "Range of temperature for thermodynamic fits\n" ) ;
  fprintf( echofile, "No. \t Species \t\t Tlow \tTmid \tThigh\n") ;
  for ( i=0 ; i<TC_Nspec_ ; i++ )
  {
    fprintf( echofile, "%-3d\t%-32s %12.4f\t%12.4f\t%12.4f\n",
	     i+1,&(TC_sNames_[i*LENGTHOFSPECNAME]),TC_Tlo_[i],TC_Tmi_[i],TC_Thi_[i]) ;
  }
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  /* Polynomial coeffs for thermo fits */
  TC_cppol_ = (double *) malloc ( TC_nNASAinter_*(TC_nCpCoef_+2)*TC_Nspec_*sizeof(double) ) ;
  for ( i=0 ; i< TC_Nspec_ ; i++ )
    for ( j=0 ; j<TC_nNASAinter_ ; j++ )
      for ( k=0 ; k<TC_nCpCoef_+2 ; k++ )
	fscanf( chemfile,"%lf", &(TC_cppol_[i*(TC_nCpCoef_+2)*TC_nNASAinter_+j*(TC_nCpCoef_+2)+k]) ) ;

  fprintf( echofile, "List of coefficients for thermodynamic fits\n" ) ;
  for ( i=0 ; i<TC_Nspec_ ; i++ )
  {
    int indx = i*(TC_nCpCoef_+2)*TC_nNASAinter_ ;
    fprintf( echofile, "%-4d %-32s\n %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	     i+1,&(TC_sNames_[i*LENGTHOFSPECNAME]),TC_cppol_[indx],TC_cppol_[indx+1],TC_cppol_[indx+2],
	     TC_cppol_[indx+3],TC_cppol_[indx+4],TC_cppol_[indx+5],TC_cppol_[indx+6]) ;
    fprintf( echofile, " %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	     TC_cppol_[indx+7],TC_cppol_[indx+8],TC_cppol_[indx+9],
	     TC_cppol_[indx+10],TC_cppol_[indx+11],TC_cppol_[indx+12],TC_cppol_[indx+13]) ;
  }
  fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  /* Polynomial coeffs for 9-term thermo fits */
  if ( TC_nNASA9coef_ > 0 )
  {

    fscanf( chemfile,"%d",  &TC_nNASA9coef_ ) ;

    TC_spec9t_    = (int *) malloc ( TC_nNASA9coef_ * sizeof(int) ) ;
    TC_spec9nrng_ = (int *) malloc ( TC_nNASA9coef_ * sizeof(int) ) ;
    TC_spec9trng_  = (double *) malloc ( TC_nNASA9coef_*NTH9RNGMAX*2 * sizeof(double) ) ;
    TC_spec9coefs_ = (double *) malloc ( TC_nNASA9coef_*NTH9RNGMAX*9 * sizeof(double) ) ;


    for (i=0;i<TC_nNASA9coef_;i++)
    {
      fscanf( chemfile,"%d", &TC_spec9t_[i]   ) ; TC_spec9t_[i] -= 1;
      fscanf( chemfile,"%d", &TC_spec9nrng_[i] ) ;

      for ( j = 0; j<TC_spec9nrng_[i]; j++ )
      {
	fscanf( chemfile,"%lf", &TC_spec9trng_[i*NTH9RNGMAX*2+2*j  ] ) ;
	fscanf( chemfile,"%lf", &TC_spec9trng_[i*NTH9RNGMAX*2+2*j+1] ) ;
        for ( k = 0; k < 9; k++ )
	  fscanf( chemfile,"%lf", &TC_spec9coefs_[i*NTH9RNGMAX*9+k]);

      } /* done loop over temperature ranges */
    } /* done loop over species */
  } /* done if for species with 9-coefficients thermo props */


  /* Ionic species */
  if ( TC_nIonEspec_ > 0 )
  {
    fscanf( chemfile,"%d", &TC_nIonEspec_ ) ;
    TC_sNion_ = (int *) malloc ( TC_nIonEspec_ * sizeof(int) ) ;
    for ( i=0 ; i<TC_nIonEspec_ ; i++ )
      { fscanf( chemfile,"%d", &(TC_sNion_[i]) ) ; TC_sNion_[i] -= 1 ;}
  }

  /* reaction info */
  if ( TC_Nreac_ > 0)
  {

    /* no of reactants and products */
    TC_isRev_   = (int *) malloc ( TC_Nreac_ * sizeof(int) ) ;
    TC_reacNrp_ = (int *) malloc ( TC_Nreac_ * sizeof(int) ) ;
    for ( i=0 ; i<TC_Nreac_ ; i++ )
    {
      fscanf( chemfile,"%d", &itemp ) ;
      if (itemp<0) TC_isRev_[i] = 0 ;
      else         TC_isRev_[i] = 1 ;
      TC_reacNrp_[i] = abs(itemp) ;
    }

    /* no of reactants only */
    TC_reacNreac_ = (int *) malloc ( TC_Nreac_ * sizeof(int) ) ;
    for ( i=0 ; i<TC_Nreac_ ; i++ ) fscanf( chemfile,"%d", &(TC_reacNreac_[i]) ) ;

    /* no of products */
    TC_reacNprod_ = (int *) malloc ( TC_Nreac_ * sizeof(int) ) ;
    for ( i=0 ; i<TC_Nreac_ ; i++ ) TC_reacNprod_[i] = TC_reacNrp_[i]-TC_reacNreac_[i] ;

    /* Stoichiometric coefficients */
    TC_reacNuki_  = (int *) malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(int) ) ;
    TC_reacNukiDbl_ = (double *) malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(double) ) ;
    TC_reacSidx_  = (int *) malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(int) ) ;
    TC_reacScoef_ = (int *) malloc( TC_Nreac_* sizeof(int) ) ;

    for ( i=0 ; i<TC_Nreac_ ; i++ )
    {

      int nusumk = 0 ;
      /* by default reaction has integer stoichiometric coefficients */
      TC_reacScoef_[i] = -1 ;

      nusumk = 0 ;
      for ( j = 0 ; j<TC_maxSpecInReac_ ; j++ )
      {
	      fscanf( chemfile,"%d", &(TC_reacNuki_[i*TC_maxSpecInReac_+j]) ) ;
	      fscanf( chemfile,"%d", &(TC_reacSidx_[i*TC_maxSpecInReac_+j]) ) ;
	      TC_reacSidx_[i*TC_maxSpecInReac_+j] -= 1 ;
	      nusumk += TC_reacNuki_[i*TC_maxSpecInReac_+j] ;
        TC_reacNukiDbl_[i*TC_maxSpecInReac_+j] = TC_reacNuki_[i*TC_maxSpecInReac_+j];
      }
      fscanf( chemfile,"%d", &itemp ) ;
      assert(itemp == nusumk) ;
    }

    /* Arrhenius parameters */
    TC_reacArhenFor_  = (double *) malloc( TC_Nreac_*3* sizeof(double) ) ;
    for ( i=0 ; i< TC_Nreac_ ; i++ )
    {
      fscanf( chemfile,"%lf", &(TC_reacArhenFor_[i*3  ]) ) ;
      fscanf( chemfile,"%lf", &(TC_reacArhenFor_[i*3+1]) ) ;
      fscanf( chemfile,"%lf", &(TC_reacArhenFor_[i*3+2]) ) ;
    }

    fprintf( echofile, "Reaction data : species and Arrhenius pars\n" ) ;
    for ( i=0 ; i<TC_Nreac_ ; i++ )
    {
      int indx = i*TC_maxSpecInReac_ ;
      fprintf(echofile,"%-5d\t%1d\t%2d\t%2d | ",i+1,TC_isRev_[i],TC_reacNreac_[i],TC_reacNprod_[i]);
      for ( j = 0; j<TC_reacNreac_[i] ; j++)
	fprintf(echofile,"%d*%s | ",TC_reacNuki_[indx+j],
		&(TC_sNames_[TC_reacSidx_[indx+j]*LENGTHOFSPECNAME]));

      indx = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
	fprintf(echofile,"%d*%s | ",TC_reacNuki_[indx+j],
		&(TC_sNames_[TC_reacSidx_[indx+j]*LENGTHOFSPECNAME]));

      indx = i*3;
      fprintf(echofile,"%16.8e\t%16.8e\t%16.8e\n",
	      TC_reacArhenFor_[indx],
	      TC_reacArhenFor_[indx+1],
	      TC_reacArhenFor_[indx+2]) ;

    }
    fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  }
#ifdef VERBOSE
  printf("TC_initChem() : Done reading reaction data\n") ;
#endif

  /* Reactions with reversible Arrhenius parameters given */
  if ( TC_nRevReac_ > 0 )
  {
    /* No. of such reactions */
    fscanf( chemfile,"%d", &TC_nRevReac_ ) ;

    /* Their indices */
    TC_reacRev_= (int *) malloc ( TC_nRevReac_ * sizeof(int) ) ;
    for ( i=0 ; i < TC_nRevReac_ ; i++ )
      { fscanf( chemfile,"%d", &(TC_reacRev_[i])) ; TC_reacRev_[i] -= 1 ; }

    /* reverse Arrhenius parameters */
    TC_reacArhenRev_ = (double *) malloc( (TC_nRevReac_*3) * sizeof(double) ) ;
    for ( i=0 ; i< TC_nRevReac_ ; i++ )
    {
      fscanf( chemfile,"%lf", &(TC_reacArhenRev_[i*3  ]) ) ;
      fscanf( chemfile,"%lf", &(TC_reacArhenRev_[i*3+1]) ) ;
      fscanf( chemfile,"%lf", &(TC_reacArhenRev_[i*3+2]) ) ;
    }

    fprintf( echofile, "Reaction data : reverse Arrhenius pars for %d reactions :\n",TC_nRevReac_ ) ;
    for ( i=0 ; i<TC_nRevReac_ ; i++ )
    {
      int indx = i*3;
      fprintf(echofile,"%-4d\t%16.8e\t%16.8e\t%16.8e\n",TC_reacRev_[i]+1,
	      TC_reacArhenRev_[indx],
	      TC_reacArhenRev_[indx+1],
	      TC_reacArhenRev_[indx+2]) ;

    }
    fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  } /* Done if  TC_nRevReac_ > 0 */

#ifdef VERBOSE
  printf("TC_initChem() : Done reading reverse reaction data\n") ;
#endif

  /* Pressure-dependent reactions */
  if ( TC_nFallReac_ > 0 )
  {

    fscanf( chemfile,"%d", &TC_nFallReac_   ) ;
    fscanf( chemfile,"%d", &TC_nFallPar_ ) ;

    TC_reacPfal_  = (int *) malloc ( TC_nFallReac_ * sizeof(int) ) ;
    TC_reacPtype_ = (int *) malloc ( TC_nFallReac_ * sizeof(int) ) ;
    TC_reacPlohi_ = (int *) malloc ( TC_nFallReac_ * sizeof(int) ) ;
    TC_reacPspec_ = (int *) malloc ( TC_nFallReac_ * sizeof(int) ) ;

    for ( i=0 ; i < TC_nFallReac_ ; i++ )
    {
      fscanf( chemfile,"%d", &(TC_reacPfal_ [i])) ; TC_reacPfal_ [i] -= 1 ;
      fscanf( chemfile,"%d", &(TC_reacPtype_[i])) ;
      fscanf( chemfile,"%d", &(TC_reacPlohi_[i])) ;
      fscanf( chemfile,"%d", &(TC_reacPspec_[i])) ; TC_reacPspec_[i] -= 1 ;
    }

    TC_reacPpar_ = (double *) malloc( TC_nFallReac_*TC_nFallPar_ *sizeof(double)) ;
    for ( i=0 ; i < TC_nFallReac_*TC_nFallPar_ ; i++ )  fscanf( chemfile,"%lf", &(TC_reacPpar_[i])) ;

    fprintf( echofile,"Reaction data : Pressure dependencies for %d reactions :\n",TC_nFallReac_) ;
    for ( i=0 ; i<TC_nFallReac_ ; i++ )
    {
      fprintf( echofile,"%-4d\t",TC_reacPfal_[i]+1) ;

      if (TC_reacPtype_[i] == 1) fprintf( echofile,"Lind \t") ;
      if (TC_reacPtype_[i] == 2) fprintf( echofile,"SRI  \t") ;
      if (TC_reacPtype_[i] == 3) fprintf( echofile,"Troe3\t") ;
      if (TC_reacPtype_[i] == 4) fprintf( echofile,"Troe4\t") ;
      if (TC_reacPtype_[i] == 6) fprintf( echofile,"Cheb \t") ;

      if (TC_reacPlohi_[i] == 0) fprintf( echofile,"Low  \t") ;
      if (TC_reacPlohi_[i] == 1) fprintf( echofile,"High \t") ;

      if (TC_reacPspec_[i] <  0) fprintf( echofile,"Mixture") ;
      if (TC_reacPspec_[i] >= 0) fprintf( echofile,"%s\n", &(TC_sNames_[TC_reacPspec_[i]*LENGTHOFSPECNAME])) ;

    }
    fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  } /* Done if fall-off reactions */
#ifdef VERBOSE
  printf("TC_initChem() : Done reading pressure-dependent reaction data\n" );
#endif

  /* Third-body reactions */
  if ( TC_nThbReac_ > 0 )
  {

    fscanf( chemfile,"%d", &TC_nThbReac_  ) ;
    fscanf( chemfile,"%d", &TC_maxTbInReac_ ) ;

    TC_reacTbdy_ = (int *) malloc ( TC_nThbReac_ * sizeof(int) ) ;
    TC_reacTbno_ = (int *) malloc ( TC_nThbReac_ * sizeof(int) ) ;

    for ( i=0 ; i < TC_nThbReac_ ; i++ )
    {
      fscanf( chemfile,"%d", &(TC_reacTbdy_[i])) ; TC_reacTbdy_ [i] -= 1 ;
      fscanf( chemfile,"%d", &(TC_reacTbno_[i])) ;
    }

    TC_specTbdIdx_ = (int    *) malloc ( TC_nThbReac_*TC_maxTbInReac_ * sizeof(int   ) ) ;
    TC_specTbdEff_ = (double *) malloc ( TC_nThbReac_*TC_maxTbInReac_ * sizeof(double) ) ;

    for ( i=0 ; i < TC_nThbReac_ ; i++ )
      for ( j=0 ; j < TC_maxTbInReac_ ; j++ )
      {
	fscanf( chemfile,"%d", &itemp) ;
	TC_specTbdIdx_[i*TC_maxTbInReac_+j] = itemp-1 ;
      }
    for ( i=0 ; i < TC_nThbReac_*TC_maxTbInReac_ ; i++ ) fscanf( chemfile,"%lf", &(TC_specTbdEff_[i])) ;

    fprintf( echofile,"Reaction data : Third body efficiencies :\n") ;
    for ( i=0 ; i<TC_nThbReac_ ; i++ )
    {
      int indx = i*TC_maxTbInReac_;
      fprintf( echofile,"%-4d\t",TC_reacTbdy_[i]+1) ;
      for ( j=0 ; j < TC_reacTbno_[i] ; j++ )
	fprintf( echofile,"%s->%5.2f, ",&(TC_sNames_[TC_specTbdIdx_[indx+j]*LENGTHOFSPECNAME]),
		 TC_specTbdEff_[indx+j]) ;
      fprintf( echofile,"\n") ;

    }
    fprintf(echofile,"-------------------------------------------------------------------------\n") ;

  }
#ifdef VERBOSE
  printf("TC_initChem() : Done reading third-body data\n" );
#endif

  /* Reactions with real stoichiometric coefficients */
  if ( TC_nRealNuReac_ > 0 )
  {

    fscanf( chemfile,"%d", &TC_nRealNuReac_  ) ;

    TC_reacRnu_   = (int *) malloc ( TC_nRealNuReac_ * sizeof(int) ) ;
    for ( i=0 ; i < TC_nRealNuReac_ ; i++ )
    {
      fscanf( chemfile,"%d", &(TC_reacRnu_[i])) ; TC_reacRnu_[i] -= 1 ;
      TC_reacScoef_[TC_reacRnu_[i]] = i ;
    }

    /* Real stoichiometric coefficients */
    TC_reacRealNuki_ = (double *) malloc( TC_nRealNuReac_ * TC_maxSpecInReac_ * sizeof(double) ) ;
    for ( i=0 ; i<TC_nRealNuReac_ ; i++ )
    {
      double dtemp ;
      double nusumk = 0.0 ;
      for ( j = 0 ; j<TC_maxSpecInReac_ ; j++ )
      {
	fscanf( chemfile,"%lf", &(TC_reacRealNuki_[i*TC_maxSpecInReac_+j])) ;
	nusumk += TC_reacRealNuki_[i*TC_maxSpecInReac_+j] ;
      }
      fscanf( chemfile,"%lf", &dtemp) ;
      assert(fabs(dtemp-nusumk)<reacbalance) ;
    }

  }

  /* Arbitrary reaction orders */
  if ( TC_nOrdReac_ > 0 )
  {

    fscanf( chemfile,"%d", &TC_nOrdReac_   ) ;
    fscanf( chemfile,"%d", &TC_maxOrdPar_ ) ;

    TC_reacAOrd_ = (int *) malloc ( TC_nOrdReac_ * sizeof(int) ) ;

    for ( i=0 ; i < TC_nOrdReac_ ; i++ )
    {
      fscanf( chemfile,"%d", &(TC_reacAOrd_[i])) ; TC_reacAOrd_[i] -= 1 ;
    }

    TC_specAOidx_ = (int    *) malloc( TC_nOrdReac_*TC_maxOrdPar_*sizeof(int   )) ;
    TC_specAOval_ = (double *) malloc( TC_nOrdReac_*TC_maxOrdPar_*sizeof(double)) ;
    for ( i=0 ; i < TC_nOrdReac_*TC_maxOrdPar_ ; i++ ) fscanf( chemfile,"%d", &(TC_specAOidx_[i])) ;
    for ( i=0 ; i < TC_nOrdReac_*TC_maxOrdPar_ ; i++ ) fscanf( chemfile,"%lf", &(TC_specAOval_[i])) ;

    printf("Reading arbitrary orders: %d, %d \n",TC_nOrdReac_,TC_maxOrdPar_) ;

    for ( i=0 ; i < TC_nOrdReac_*TC_maxOrdPar_ ; i++ ) printf("%d, ",TC_specAOidx_[i] ) ;
    printf("\n") ;

    for ( i=0 ; i < TC_nOrdReac_*TC_maxOrdPar_ ; i++ ) printf("%e, ",TC_specAOval_[i] ) ;
    printf("\n") ;

  }

  /* Landau-Teller reactions */
  if ( TC_nLtReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Landau-Teller reactions are not supported yet !!!") ;

  /* reverse Landau-Teller reactions */
  if ( TC_nRltReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Landau-Teller reactions are not supported yet !!!") ;

  /* radiation wavelength */
  if ( TC_nHvReac_ > 0 )
  {

    TC_reacHvIdx_ = (int    *) malloc ( TC_nHvReac_ * sizeof(int   ) ) ;
    TC_reacHvPar_ = (double *) malloc ( TC_nHvReac_ * sizeof(double) ) ;

    for ( i=0 ; i < TC_nHvReac_ ; i++ )
    {
      fscanf( chemfile,"%d", &(TC_reacHvIdx_[i])) ; TC_reacHvIdx_[i] -= 1 ;
    }

    for ( i=0 ; i < TC_nHvReac_ ; i++ ) fscanf( chemfile,"%lf", &(TC_reacHvPar_[i])) ;
  }
  //TC_errorINI(errfile,"kmod.list : Error: radiation wavelength not supported yet !!!") ;

  /* JAN fits */
  if ( TC_nJanReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: JAN fits are not supported yet !!!") ;

  /* FIT1 fits */
  if ( TC_nFit1Reac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: FIT1 fits are not supported yet !!!") ;

  /* Excitation reactions */
  if ( TC_nExciReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Excitation reactions are not supported yet !!!") ;

  /* Plasma Momentum transfer collision */
  if ( TC_nMomeReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Plasma Momentum transfer collision is not supported yet !!!") ;

  /* Ion Momentum transfer collision */
  if ( TC_nXsmiReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Ion Momentum transfer collision is not supported yet !!!") ;

  /* Species temperature dependency */
  if ( TC_nTdepReac_ > 0 ) TC_errorINI(errfile,"kmod.list : Error: Custom species temperature dependency is not supported yet !!!") ;

  /* Reactions with PLOG formulation */
  if ( TC_nPlogReac_ > 0 ) {
    fscanf( chemfile,"%d", &TC_nPlogReac_  ) ;

    TC_reacPlogIdx_ = (int *) malloc ( TC_nPlogReac_     * sizeof(int) ) ;
    TC_reacPlogPno_ = (int *) malloc ( (TC_nPlogReac_+1) * sizeof(int) ) ;

    /* Their indices and no of PLOG intervals */
    TC_reacPlogPno_[0] = 0;
    for ( i = 0; i < TC_nPlogReac_; i++) {
      fscanf( chemfile,"%d", &(TC_reacPlogIdx_[i  ]) ) ; TC_reacPlogIdx_[i] -= 1 ;
      fscanf( chemfile,"%d", &(TC_reacPlogPno_[i+1]) ) ; TC_reacPlogPno_[i+1] += TC_reacPlogPno_[i] ;
    }

    /* Plog parameters */
    TC_reacPlogPars_ = (double *) malloc( TC_reacPlogPno_[TC_nPlogReac_] * 4 * sizeof(double) ) ;
    for ( i=0 ; i<TC_reacPlogPno_[TC_nPlogReac_] ; i++ ) {
      fscanf( chemfile,"%lf", &(TC_reacPlogPars_[4*i  ])) ; TC_reacPlogPars_[4*i  ] = log(TC_reacPlogPars_[4*i  ]);
      fscanf( chemfile,"%lf", &(TC_reacPlogPars_[4*i+1])) ; TC_reacPlogPars_[4*i+1] = log(TC_reacPlogPars_[4*i+1]);
      fscanf( chemfile,"%lf", &(TC_reacPlogPars_[4*i+2])) ;
      fscanf( chemfile,"%lf", &(TC_reacPlogPars_[4*i+3])) ;
    }

  } /* Done if nPlogReac > 0 */

  fclose(chemfile) ;
  fclose(errfile ) ;
  fclose(echofile) ;

  /* no. of variables (temperature + species) */
  TC_Nvars_ = TC_Nspec_ + 1 ;

  /* no. of variables int the Jacobian(rho,P,T,species) */
  TC_Nvjac_ = TC_Nspec_ + 3 ;

  /* default thermodynamic pressure */
  TC_pressure_ = ATMPA ;
  TC_prescgs_  = TC_pressure_*10.0 ;

  /* universal gas constant */
  TC_Runiv_ = RUNIV*1.0e3 ;
  TC_Rcal_  = TC_Runiv_/(CALJO*1.0e3) ;
  TC_Rcgs_  = TC_Runiv_*1.e4 ;

  /* allocate work arrays */
  TC_makeSpace() ;

  /* populate some arrays */
  /* sum(nu) for each reaction */
  for ( i = 0 ; i < TC_Nreac_ ; i++ )
  {
    int indx ;
    TC_sigNu[i] = 0 ;
    /* reactants */
    indx = i*TC_maxSpecInReac_ ;
    for ( j = 0; j<TC_reacNreac_[i] ; j++) TC_sigNu[i] += TC_reacNuki_[indx+j] ;
    /* products */
    indx = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( j = 0; j<TC_reacNprod_[i] ; j++) TC_sigNu[i] += TC_reacNuki_[indx+j] ;

  }
  if ( TC_nRealNuReac_ > 0 )
    for ( ir = 0 ; ir < TC_nRealNuReac_ ; ir++ )
    {
      int indx, irR ;
      irR = TC_reacRnu_[ir] ;
      TC_sigRealNu[ir] = 0.0 ;
      /* reactants */
      indx = ir*TC_maxSpecInReac_ ;
      for ( j = 0; j<TC_reacNreac_[irR] ; j++) TC_sigRealNu[ir] += TC_reacRealNuki_[indx+j] ;
      /* products */
      indx = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[irR] ; j++) TC_sigRealNu[ir] += TC_reacRealNuki_[indx+j] ;

    }

  /* NuIJ=NuII''-NuIJ' */
  for ( ij = 0 ; ij < TC_Nreac_*TC_Nspec_ ; ij++ ) TC_NuIJ[ij] = 0 ;
  for ( j  = 0 ; j  < TC_Nreac_        ; j++ )
  {
    /* reactants */
    int indx = j*TC_maxSpecInReac_, kspec ;
    for ( i = 0; i<TC_reacNreac_[j] ; i++)
    {
      kspec = TC_reacSidx_[indx+i] ;
      TC_NuIJ[j*TC_Nspec_+kspec] += TC_reacNuki_[indx+i] ;
    }
    /* products */
    indx = j*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( i = 0; i<TC_reacNprod_[j] ; i++)
    {
      kspec = TC_reacSidx_[indx+i] ;
      TC_NuIJ[j*TC_Nspec_+kspec] += TC_reacNuki_[indx+i] ;
    }
  }

  if ( TC_nRealNuReac_ > 0 )
  {
    for ( ij = 0 ; ij < TC_nRealNuReac_*TC_Nspec_ ; ij++ ) TC_RealNuIJ[ij] = 0.0 ;
    for ( j  = 0 ; j  < TC_nRealNuReac_        ; j++ )
    {
      /* reactants */
      int jR=TC_reacRnu_[j];
      int indxR = j *TC_maxSpecInReac_, kspec ;
      int indx  = jR*TC_maxSpecInReac_ ;
      for ( i = 0; i<TC_reacNreac_[jR] ; i++)
      {
        kspec = TC_reacSidx_[indx+i] ;
        TC_RealNuIJ[j*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR+i] ;
      }
      printf("got here...........\n");fflush(stdout);
      /* products */
      indxR = j *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indx  = jR*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( i = 0; i<TC_reacNprod_[jR] ; i++)
      {
        kspec = TC_reacSidx_[indx+i] ;
        TC_RealNuIJ[j*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR+i] ;
      }
    }
  }

  TC_setThermoPres( TC_pressure_ );

  /* Store coefficients for kc */
  TC_kc_coeff = (double*) malloc( TC_Nreac_ * sizeof(double));
  for ( i = 0 ; i<TC_Nreac_ ; i++ )
  {
    int ir = TC_reacScoef_[i] ;
    if ( ir == -1 )
    {
      //TC_kc_coeff[i] = fastIntPow(TC_prescgs_/TC_Rcgs_,TC_sigNu[i]) ;
      TC_kc_coeff[i] = fastIntPow(ATMPA*10.0/TC_Rcgs_,TC_sigNu[i]) ;
    }
    else
    {
      //TC_kc_coeff[i] = pow(TC_prescgs_/TC_Rcgs_,TC_sigRealNu[ir]);
      TC_kc_coeff[i] = pow(ATMPA*10.0/TC_Rcgs_,TC_sigRealNu[ir]);
    }
  }

  if ( tab == 1 ) {
    if ( TC_nPlogReac_ > 0 ) {
      TC_errorINI(errfile,"Error: Cannot create tables when PLOG keyword is present !!!") ;
    }
    TC_createTables( delT ) ;
  } /* done if tables are wanted */

  /* done */
  TC_isInit_ = 1 ;

  return ( 0 ) ;

}
/*
               __            _
     _ __ ___ / _|_   ____ _| |___
    | '__/ _ \ |_\ \ / / _` | / __|
    | | |  __/  _|\ V / (_| | \__ \
    |_|  \___|_|   \_/ \__,_|_|___/

*/
/**
 * \ingroup init
 * \brief Set reference values to the library.
*/
void TC_setRefVal(double rhoref, double pref, double Tref, double Wref, double Daref,
		  double omgref, double cpref, double href, double timref)
{
/**
  \param rhoref : density \f$[kg/m^3]\f$
  \param pref : pressure \f$[N/m^2]\f$
  \param Tref : temperature \f$[K]\f$
  \param Wref : molecular weight \f$[kg/kmol]\f$
  \param Daref : Damkohler number []
  \param omgref : molar reaction rate \f$[kmol/(m^3\cdot s)]\f$
  \param cpref : specific heat at constant pressure \f$[J/(kg\cdot K)]\f$
  \param href : specific enthalphy \f$[J/kg]\f$
  \param timref : time \f$[s]\f$
*/

  TC_rhoref_ = rhoref ;
  TC_pref_   = pref   ;
  TC_Tref_   = Tref   ;
  TC_Wref_   = Wref   ;
  TC_Daref_  = Daref  ;
  TC_omgref_ = omgref ;
  TC_cpref_  = cpref  ;
  TC_href_   = href   ;
  TC_timref_ = timref ;

  TC_RVset_  = 1 ;

  return ;

}

/*
    ____  _              ___   _             ____  _
   |  _ \(_)_ __ ___    / / \ | | ___  _ __ |  _ \(_)_ __ ___
   | | | | | '_ ` _ \  / /|  \| |/ _ \| '_ \| | | | | '_ ` _ \
   | |_| | | | | | | |/ / | |\  | (_) | | | | |_| | | | | | | |
   |____/|_|_| |_| |_/_/  |_| \_|\___/|_| |_|____/|_|_| |_| |_|

*/
/**
 * \ingroup init
 *  \brief Set library to function in non-dimensional mode
 */
void TC_setNonDim()
{
  if ( TC_RVset_ == 0 )
  {
    printf("TC_setNonDim() : Cannot set non-dimensional flag before reference values\n");
    printf("                 -> Abort !!!\n") ;
    exit(1) ;
  }

  TC_nonDim_ = 1 ;

  return ;

}
/**
 * \ingroup init
 * \brief Set library to function in dimensional mode (default)
*/
void TC_setDim()
{
  TC_nonDim_ = 0 ;
  return ;
}
/**
 * \ingroup init
 * \brief Send thermodynamic pressure to the library.
 */
void TC_setThermoPres(double pressure)
{
/**
  \param pressure : thermodynamic pressure \f$[N/m^2]\f$
*/
  TC_pressure_ = pressure ;

  /* Need dimensional pressure -> convert the input if necessary */
  if ( TC_nonDim_ == 1 ) TC_pressure_ *= TC_pref_ ;

  TC_prescgs_  = TC_pressure_*10.0 ;

  return ;

}

/**
 * \ingroup init
 * \brief Retrieve thermodynamic pressure from the library.
 */
double TC_getThermoPres() {
  return ( TC_pressure_ ) ;
}

/**
 * \ingroup init
 * \brief Send density to the library.
 */
void TC_setDens(double density) {
/**
  \param density : mixture density \f$[kg/m^3]\f$
*/
  TC_rho_ = density ;

  /* Need dimensional density -> convert the input if necessary */
  if ( TC_nonDim_ == 1 ) TC_rho_ *= TC_rhoref_ ;

  TC_rhocgs_ = TC_rho_ * 0.001 ;

  TC_rhoset_ = 1 ;

  return ;

}
/**
 * \ingroup init
 * \brief Retrieve density from the library.
 */
double TC_getDens() { return ( TC_rho_ ) ; }



/*------------------------------------------------------------------------
                         _                  _            _
      ___  ___ _ __ ___ (_)      _ __  _ __(_)_   ____ _| |_ ___
     / __|/ _ \ '_ ` _ \| |_____| '_ \| '__| \ \ / / _` | __/ _ \
     \__ \  __/ | | | | | |_____| |_) | |  | |\ V / (_| | ||  __/
     |___/\___|_| |_| |_|_|     | .__/|_|  |_| \_/ \__,_|\__\___|
                                |_|

------------------------------------------------------------------------ */
/*
                       _        ____
       _ __ ___   __ _| | _____/ ___| _ __   __ _  ___ ___
      | '_ ` _ \ / _` | |/ / _ \___ \| '_ \ / _` |/ __/ _ \
      | | | | | | (_| |   <  __/___) | |_) | (_| | (_|  __/
      |_| |_| |_|\__,_|_|\_\___|____/| .__/ \__,_|\___\___|
                                     |_|
*/
/**
 * \ingroup init
 * \brief Allocate internal work arrays for the library. Should not be called by external functions.
*/
int TC_makeSpace()
{

  /* allocate space based on Nspec_ */
  TC_y2x2y   = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_cpks    = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_entr    = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_hks     = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_Xconc   = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_gk      = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_gkp     = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_omg     = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_omgP    = (double *) malloc( TC_Nspec_ * sizeof(double) ) ;
  TC_CrndDer = (double *) malloc((TC_Nspec_+1) * sizeof(double)) ;
  TC_PrDer   = (double *) malloc((TC_Nspec_+1) * sizeof(double)) ;
  TC_dFfac   = (double *) malloc((TC_Nspec_+1) * sizeof(double)) ;
  TC_scalIn  = (double *) malloc((TC_Nspec_+1) * sizeof(double)) ;

  TC_jacFull = (double *) malloc( (TC_Nvjac_ * TC_Nvjac_) * sizeof(double)) ;

  /* allocate space based on Nreac_ */
  TC_Mconc  = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_kfor   = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_krev   = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_kforP  = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_krevP  = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_kforPder = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_krevPder = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_ropFor = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_ropRev = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_rop    = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;
  TC_Crnd   = (double *) malloc( TC_Nreac_ * sizeof(double) ) ;

  /* allocate more space */
  TC_sigNu = (int *) malloc( TC_Nreac_ * sizeof(int) ) ;
  TC_NuIJ  = (int *) malloc( TC_Nreac_ * TC_Nspec_ * sizeof(int) ) ;

  if ( TC_nRealNuReac_ > 0 )
  {
    TC_sigRealNu = (double *) malloc( TC_nRealNuReac_ * sizeof(double) ) ;
    TC_RealNuIJ  = (double *) malloc( TC_nRealNuReac_ * TC_Nspec_ * sizeof(double) ) ;
  }

  /* work space */
  qfr = (double *) malloc( 2 * sizeof(double) ) ;

  return ( 0 ) ;

}
/*
                         _      _____     _     _
      ___ _ __ ___  __ _| |_ __|_   _|_ _| |__ | | ___  ___
     / __| '__/ _ \/ _` | __/ _ \| |/ _` | '_ \| |/ _ \/ __|
    | (__| | |  __/ (_| | ||  __/| | (_| | |_) | |  __/\__ \
     \___|_|  \___|\__,_|\__\___||_|\__,_|_.__/|_|\___||___/

*/
/**
 * \ingroup init
 * \brief Create tables for temperature dependent terms. Should not be called by external functions
*/
int TC_createTables( double delT )
{
/**
  \param delT : temperature step size [K]
*/

  int i, j, itabS, itabR ;
  int spectabsize, reactabsize ;
  double temper ;

  TC_tab_   = 1        ;
  TC_delT_  = delT     ;
  TC_odelT_ = 1.0/delT ;

  /* number of elements */
  TC_Ntab_ = 1+(int) ( (TMAX-TMIN)*TC_odelT_ );

  /* tables for species variables */
  spectabsize = TC_Nspec_*TC_Ntab_ ;
  TC_cptab  = (double *) malloc( spectabsize * sizeof(double) ) ;
  TC_cpPtab = (double *) malloc( spectabsize * sizeof(double) ) ;
  TC_htab   = (double *) malloc( spectabsize * sizeof(double) ) ;
  TC_gktab  = (double *) malloc( spectabsize * sizeof(double) ) ;
  TC_gkPtab = (double *) malloc( spectabsize * sizeof(double) ) ;

  /* tables for reaction variables */
  reactabsize = TC_Nreac_*TC_Ntab_ ;
  TC_kfortab  = (double *) malloc( reactabsize * sizeof(double) ) ;
  TC_krevtab  = (double *) malloc( reactabsize * sizeof(double) ) ;
  TC_kforPtab = (double *) malloc( reactabsize * sizeof(double) ) ;
  TC_krevPtab = (double *) malloc( reactabsize * sizeof(double) ) ;

  /* create tables with thermo properties */
  temper = TMIN ;
  for ( i=0,itabS=0,itabR=0 ; i<TC_Ntab_; i++,itabS+=TC_Nspec_,itabR+=TC_Nreac_)
  {
    double t_1, tln ;

    /* tables for cp,dcp/dT, and h */
    TC_getCpSpecMsFcn ( temper, &(TC_cptab [itabS]) ) ;
    TC_getCpSpecMsPFcn( temper, &(TC_cpPtab[itabS]) ) ;
    TC_getHspecMsFcn  ( temper, &(TC_htab  [itabS]) ) ;

    /* tables for gk and dgk/dT */
    t_1 = 1.0/temper ;
    tln = log(temper) ;
    TC_getgkFcn (temper,t_1,tln) ;
    TC_getgkpFcn(temper,t_1,tln) ;
    for ( j = 0 ; j < TC_Nspec_ ; j++ ) TC_gktab [itabS+j] = TC_gk [j] ;
    for ( j = 0 ; j < TC_Nspec_ ; j++ ) TC_gkPtab[itabS+j] = TC_gkp[j] ;

    /* tables for kfor, krev, dkfor/dT, and dkrev/dT */
    TC_getkForRevFcn (t_1,tln) ;
    TC_getkForRevPFcn(t_1    ) ;
    for ( j = 0 ; j < TC_Nreac_ ; j++ ) TC_kfortab [itabR+j] = TC_kfor [j] ;
    for ( j = 0 ; j < TC_Nreac_ ; j++ ) TC_krevtab [itabR+j] = TC_krev [j] ;
    for ( j = 0 ; j < TC_Nreac_ ; j++ ) TC_kforPtab[itabR+j] = TC_kforP[j] ;
    for ( j = 0 ; j < TC_Nreac_ ; j++ ) TC_krevPtab[itabR+j] = TC_krevP[j] ;

    temper += TC_delT_ ;

  }

#ifdef DEBUG_ADD
  /* Fix the discontinuity in the derivative of cp at Tmi */
  for ( i=0,itabS=0; i<TC_Ntab_; i++,itabS+=TC_Nspec_)
  {
    temper = TMIN+( (double) i ) * TC_delT_ ;

    for (j=0; j<TC_Nspec_; j++)
    {
      if ( (temper>TC_Tmi_[j]-DTMID) && (temper<TC_Tmi_[j]+DTMID) )
      {
	double cp1, cp2 ;
	TC_getCpSpecMs1PFcn( TC_Tmi_[j]-DTMID, j, &cp1 ) ;
	TC_getCpSpecMs1PFcn( TC_Tmi_[j]-DTMID, j, &cp2 ) ;
        TC_cpPtab[itabS+j] = cp1+(cp2-cp1)/(2.0*DTMID)*(temper-TC_Tmi_[j]+DTMID);
      }
    }
  }
#endif

#ifdef DEBUG_OUT
  FILE *fout1, *fout2, *fout3, *fout4, *fout5 ;
  fout1 = fopen("cpTab.dat" , "w") ;
  fout2 = fopen("cpPTab.dat", "w") ;
  fout3 = fopen("hTab.dat"  , "w") ;
  fout4 = fopen("gkTab.dat" , "w") ;
  fout5 = fopen("gkPTab.dat", "w") ;
  for ( i = 0, itabS = 0 ; i < TC_Ntab_ ; i++, itabS+=TC_Nspec_ )
  {
    temper = TMIN+( (double) i ) * TC_delT_ ;
    fprintf(fout1,"%5d %25.16 ", i, temper ) ;
    fprintf(fout2,"%5d %25.16 ", i, temper ) ;
    fprintf(fout3,"%5d %25.16 ", i, temper ) ;
    fprintf(fout4,"%5d %25.16 ", i, temper ) ;
    fprintf(fout5,"%5d %25.16 ", i, temper ) ;
    for (j=0; j<TC_Nspec_; j++)
    {
      fprintf(fout1,"%25.16 ", TC_cptab [itabS+j] ) ;
      fprintf(fout2,"%25.16 ", TC_cpPtab[itabS+j] ) ;
      fprintf(fout3,"%25.16 ", TC_htab  [itabS+j] ) ;
      fprintf(fout4,"%25.16 ", TC_gktab [itabS+j] ) ;
      fprintf(fout5,"%25.16 ", TC_gkPtab[itabS+j] ) ;
    }
    fprintf(fout1,"\n") ;
    fprintf(fout2,"\n") ;
    fprintf(fout3,"\n") ;
    fprintf(fout4,"\n") ;
    fprintf(fout5,"\n") ;
  }
  fclose(fout1) ;
  fclose(fout2) ;
  fclose(fout3) ;
  fclose(fout4) ;
  fclose(fout5) ;

  fout1 = fopen("kforTab.dat" , "w") ;
  fout2 = fopen("krevTab.dat" , "w") ;
  fout3 = fopen("kforPTab.dat", "w") ;
  fout4 = fopen("krevPTab.dat", "w") ;
  for ( i = 0, itabR = 0 ; i < TC_Ntab_ ; i++, itabR+=TC_Nreac_ )
  {
    temper = TMIN+( (double) i ) * TC_delT_ ;
    fprintf(fout1,"%5d %25.16 ", i, temper ) ;
    fprintf(fout2,"%5d %25.16 ", i, temper ) ;
    fprintf(fout3,"%5d %25.16 ", i, temper ) ;
    fprintf(fout4,"%5d %25.16 ", i, temper ) ;
    for (j=0; j<TC_Nreac_; j++)
    {
      fprintf(fout1,"%25.16 ", TC_kfortab [itabR+j] ) ;
      fprintf(fout2,"%25.16 ", TC_krevtab [itabR+j] ) ;
      fprintf(fout3,"%25.16 ", TC_kforPtab[itabR+j] ) ;
      fprintf(fout4,"%25.16 ", TC_krevPtab[itabR+j] ) ;
    }
    fprintf(fout1,"\n") ;
    fprintf(fout2,"\n") ;
    fprintf(fout3,"\n") ;
    fprintf(fout4,"\n") ;
  }
  fclose(fout1) ;
  fclose(fout2) ;
  fclose(fout3) ;
  fclose(fout4) ;
#endif

  return ( 0 ) ;

}
/*
                          __  __ ____   ____
  ___ _ __ _ __ ___  _ __|  \/  / ___| / ___|
 / _ \ '__| '__/ _ \| '__| |\/| \___ \| |  _
|  __/ |  | | | (_) | |  | |  | |___) | |_| |
 \___|_|  |_|  \___/|_|  |_|  |_|____/ \____|

*/
/**
 * \brief Outputs error messages for the library, then exits the execution
*/
void TC_errorMSG(int msgID, char const* func, int var1, int var2)
{
/**
  \param msgID : message ID
  \param func : name of function calling this error function
  \param var1 : value \#1 to be printed in the message
  \param var2 : value \#2 to be printed in the message
*/

  if ( msgID == 10 )
  {
    fprintf( stderr,"%s() : disagreement in Nspec : %d vs %d\n",
	     func, var1, var2 ) ;
    fprintf( stderr,"              -> Abort !!!\n") ;
    fflush( stderr ) ;
    exit(-1) ;
  }
  else if ( msgID == 20 )
  {
    fprintf( stderr,"%s() : disagreement in Nvars : %d vs %d\n",
	     func, var1, var2 ) ;
    fprintf( stderr,"              -> Abort !!!\n") ;
    fflush( stderr ) ;
    exit(-1) ;
  }
  else if ( msgID == 30 )
  {
    fprintf( stderr,"%s() : disagreement in Nreac : %d vs %d\n",
	     func, var1, var2 ) ;
    fprintf( stderr,"              -> Abort !!!\n") ;
    fflush( stderr ) ;
    exit(-1) ;
  }
  else
  {
    fprintf( stderr,"Unknown message if in %s() : %d vs %d\n",
	     func, var1, var2 ) ;
    fprintf( stderr,"              -> Abort !!!\n") ;
    fflush( stderr ) ;
    exit(-1) ;
  }

  return ;

}
/**
 * \brief Outputs error messages generated by TC_initChem
*/
void TC_errorINI(FILE *errfile, char const* msg)
{
/**
  \param errfile : file pointer for output
  \param msg : error message
*/


  fprintf(errfile,"%s\n",msg );
  fprintf(errfile,"-> Abort !!!\n") ;
  fflush(errfile) ;
  fclose(errfile) ;
  exit(1) ;

  return ;

}
