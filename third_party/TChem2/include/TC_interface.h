/*! \file TC_interface.h
    \brief Header file to be included in user's code. Contains function definitions.
*/ 
#ifndef _TCinterfHSeen_
#define _TCinterfHSeen_

#include "copyright.h"

/* This is a C library, not C++ */
#ifdef __cplusplus
extern "C"
{
#endif

/*
     __       __                        _       __            ____              
    / /______/ /_  ___  ____ ___       (_)___  / /____  _____/ __/___ _________ 
   / __/ ___/ __ \/ _ \/ __ `__ \     / / __ \/ __/ _ \/ ___/ /_/ __ `/ ___/ _ \
  / /_/ /__/ / / /  __/ / / / / /    / / / / / /_/  __/ /  / __/ /_/ / /__/  __/
  \__/\___/_/ /_/\___/_/ /_/ /_/____/_/_/ /_/\__/\___/_/  /_/  \__,_/\___/\___/ 
                              /_____/                                           
*/
#define TCSMALL 1.e-30

int TC_getReacStr(int ireac, int lenstr,char *rstr);
/*

       _____ ____         _                 
      |_   _/ ___|    ___| |__   __ _   ___ 
        | || |       / __| '_ \ / _` | / __|
        | || |___   | (__| | | | (_| || (__ 
        |_| \____|___\___|_| |_|\__, (_)___|
                |_____|         |___/       

*/
int TC_chgArhenFor( int ireac, int ipos, double newval ) ;
int TC_chgArhenRev( int ireac, int ipos, double newval ) ;
int TC_chgArhenPresDepFor( int ireac, int ipos, double newval) ;
int TC_chgArhenForBack( int ireac, int ipos) ;
int TC_chgArhenRevBack( int ireac, int ipos) ;

/*
       _____ ____                _         
      |_   _/ ___|     __ _  ___| |_   ___ 
        | || |        / _` |/ _ \ __| / __|
        | || |___    | (_| |  __/ |_ | (__ 
        |_| \____|____\__, |\___|\__(_)___|
                |_____|___/                
*/
int TC_getArhenRev( int ireac, int ipos, double* val ) ;
int TC_getArhenFor( int ireac, int ipos, double* val ) ;
int TC_getArhenPresDepFor( int ireac, int ipos, double *val ) ;

double TC_getRUniv();

/*
         _____ ____             _ _ _
        |_   _/ ___|    ___  __| (_) |_   ___
          | || |       / _ \/ _` | | __| / __|
          | || |___   |  __/ (_| | | |_ | (__
          |_| \____|___\___|\__,_|_|\__(_)___|
                  |_____|
*/
/* Free memory and reset variables to allow TC_initchem to be called again */
void TC_reset() ;

/* Remove numRemoveReacs reactions in reacArr (which must be in ascending order). 
   Set revOnly to 1 to remove reverse only, any other number to remove forwards 
   and reverse. */
void TC_removeReaction(int* reacArr, int numRemoveReacs, int revOnly) ;

/* Restores mechainsm to original state before calling TC_removeReaction.
   Any changes from TC_chgArhenFor() and TC_chgArhenRev() are also reset. */
void TC_restoreReactions() ;

/*
             _____ ____     _       _ _         
            |_   _/ ___|   (_)_ __ (_) |_   ___ 
              | || |       | | '_ \| | __| / __|
              | || |___    | | | | | | |_ | (__ 
              |_| \____|___|_|_| |_|_|\__(_)___|
                      |_____|                   
      
  
*/
/* functions controlling the 7/9 term NASA polynomial settings */
void TC_set7term_  () ;
void TC_unset7term_() ;
void TC_set9term_  () ;
void TC_unset9term_() ;
void TC_setthf9_(char *thf9)  ;
void TC_unsetthf9_ () ;

/* initializes library */
int TC_initChem(const char *mechfile, const char *thermofile, int tab, double delT) ;
int TC_initChemFromFile(FILE *mechin, FILE *thermoin, int tab, double delT);

/* set reference values, can be called only after TC_setNonDim() */
void TC_setRefVal(double rhoref, double pref, double Tref, double Wref, 
		  double Daref, double omgref, double cpref, double href,
		  double timref) ;  

/* set non-dimensional/dimensional option; default is dimensional */
void TC_setNonDim() ;  
void TC_setDim() ;

/* thermodynamic pressure */
void   TC_setThermoPres( double pressure ) ;
double TC_getThermoPres() ;

/* density */
void TC_setDens( double density ) ;
double TC_getDens() ;

/* 
          _____ ____                               
         |_   _/ ___|    ___ _ __   ___  ___ 
           | || |       / __| '_ \ / _ \/ __|
           | || |___    \__ \ |_) |  __/ (__ 
           |_| \____|___|___/ .__/ \___|\___|
                   |_____|  |_|                    

*/
int    TC_getNspec() ;      /* get no. of species   */
int    TC_getNelem() ;      /* get no. of elements  */
int    TC_getNvars() ;      /* get no. of variables */
double TC_getThermoPres() ; /* get thermo pressure  */

int    TC_getEnameLen(); /* get maximum length of element names */
int    TC_getSnameLen(); /* get maximum length of species names */
int    TC_getEnames(int Nelem,char *enames) ; /* get element names */
int    TC_getSnames(int Nspec,char *snames) ; /* get species names */

int    TC_getSpos  (const char *sname, const int slen)  ; /* get species position  */
int    TC_getSmass (int Nspec,double *Wi)   ; /* get species molar weights */

/* Returns list of element counts */
int TC_getECount(int *elemcnt) ;
 
/*
    _____ ____     _        __       
   |_   _/ ___|   (_)_ __  / _| ___  
     | || |       | | '_ \| |_ / _ \ 
     | || |___    | | | | |  _| (_) |
     |_| \____|___|_|_| |_|_|  \___/ 
             |_____|                 

*/
/* no. of reactions */
int TC_getNreac() ;       
/* no. of reactions with low/high pressure limit rate constants */
int TC_getNPresDepReac();
/* no. of reactions with reaction rates enhanced by third-body efficiencies */
int TC_getNThbReac() ;
  /* no. of reactions with reaction rates in PLOG format */
int TC_getNplogReac() ;


/* Returns 1 if a reaction is reversible, 0 otherwise */
int TC_getReacRev(int ireac) ;
/* Returns 1 if a reaction employs third-body efficiencies, 0 otherwise */
int TC_getReacThb(int ireac) ;
/* Returns 1 if a reaction is of "LOW" type and 0 otherwise */
int TC_getReacLow(int ireac) ; 
/* Returns 1 if a reaction is of "HIGH" type and 0 otherwise */
int TC_getReacHIGH(int ireac) ;
/* Returns 1 if a reaction is of "PLOG" type and 0 otherwise */
int TC_getReacPLOG(int ireac) ;

/* Returns type for pressure dependent reaction
   1-Lindeman/2-SRI/3-Troe3/4-Troe4/6-Chebyshev and 0 otherwise (not pressure dependent */
int TC_getReacPtype(int ireac) ;
/* Returns the number of third-body efficiencies for a certain
 reaction. Also returns the list of species indices and the 
 corresponding third-body efficiencies */
int TC_getThbReacEff(int ireac, int *specidx, double *effs) ;

/*  Returns the index of reaction in the list of reactions with third-body efficiencies */
int TC_getReacThbIdx(int ireac);
/* Returns the index of reaction in the list of fall-off reactions */
int TC_getReacFallIdx(int ireac);
/*  Returns the index of reaction in the list of reactions with real stoichiometric coeffs */
int TC_getReacRealStIdx(int ireac);

/* stoichiometric coefficients */
int TC_getStoiCoef    ( int Nspec, int Nreac, double *stoicoef ) ;
int TC_getStoiCoefFR  ( int Nspec, int Nreac, double *stoicoef ) ;
int TC_getStoiCoefReac( int Nspec, int Nreac, int ireac, int idx,   double *stoicoef ) ;
int TC_getStoiCoefIJ  ( int Nspec, int Nreac, int ireac, int kspec, double *stoicoef ) ;

/* Returns reaction string */
int TC_getReacStr(int ireac, int lenstr, char *rstr) ;

/*
          _____ ____                          
         |_   _/ ___|    _ __ _ __  
           | || |       | '__| '__| 
           | || |___    | |  | |    
           |_| \____|___|_|  |_|    
                   |_____|                    
*/

int TC_getkForRevPder(double t_1, double tln);

/* forward and reverse rate constants */
int TC_getkFR(double t, double *kF, double *kR);

/* T, Y_i -> molar reaction rates [kmol/(m3.s)]*/
int TCDND_getTY2RRml(double *scal, int Nvars, double *omega) ;
int TC_getTY2RRml(double *scal, int Nvars, double *omega) ;

/* T, Y_i -> mass reaction rates [kg/(m3.s)]*/
int TCDND_getTY2RRms(double *scal, int Nvars, double *omega) ;
int TC_getTY2RRms(double *scal, int Nvars, double *omega) ;

/* T, C_i [kmol/m3] -> molar reaction rates [kmol/(m3.s)]*/
int TCDND_getTXC2RRml(double *scal, int Nvars, double *omega) ;
int TC_getTXC2RRml(double *scal, int Nvars, double *omega) ;

/* T, C_i [kmol/m3] -> mass reaction rates [kg/(m3.s)]*/
int TCDND_getTXC2RRms(double *scal, int Nvars, double *omega) ;
int TC_getTXC2RRms(double *scal, int Nvars, double *omega) ;

/* rate-of-progress given temperature and species mass fractions */
int TC_getRops(double *scal, int Nvars, double *datarop) ;

/* forward & reverse rate-of-progress given temperature and 
   species mass fractions */
int TC_getRfrb(double *scal, int Nvars, double *datarop) ;

/*
      _____ ____    _   _                                          
     |_   _/ ___|  | |_| |__   ___ _ __ _ __ ___   ___         ___ 
       | || |      | __| '_ \ / _ \ '__| '_ ` _ \ / _ \       / __|
       | || |___   | |_| | | |  __/ |  | | | | | | (_) |  _  | (__ 
       |_| \____|___\__|_| |_|\___|_|  |_| |_| |_|\___/  (_)  \___|
               |_____|                                             
*/
int TC_getEntr0SpecMsFcn(double t,double *s0i);
int TC_getCpSpecMsFcn(double t,double *cpi);
int TC_getEntr0SpecMlFcn(double t,double *s0i);
int TC_getMl2EntrMixMl(double *scal, int Nvars, double *smix);
int TC_getMs2EntrMixMs(double *scal, int Nvars, double *smix);

/* T, Y_i -> rho [kg/m3] */
int TCDND_getRhoMixMs(double *scal,int Nvars,double *rhomix) ;
int TC_getRhoMixMs   (double *scal,int Nvars,double *rhomix) ;

/* T, X_i -> rho [kg/m3] */
int TCDND_getRhoMixMl(double *scal,int Nvars,double *rhomix) ;
int TC_getRhoMixMl   (double *scal,int Nvars,double *rhomix) ;

/* rho[kg/m3], Y_i -> T [K] */
int TCDND_getTmixMs(double *scal,int Nvars,double *Tmix) ;
int TC_getTmixMs   (double *scal,int Nvars,double *Tmix) ;

/* rho[kg/m3], X_i -> T [K] */
int TCDND_getTmixMl(double *scal,int Nvars,double *Tmix) ;
int TC_getTmixMl   (double *scal,int Nvars,double *Tmix) ;

/* T, Y_i -> cp [J/(kg.K)] */
int TCDND_getMs2CpMixMs(double *scal,int Nvars,double *cpmix) ;
int TC_getMs2CpMixMs   (double *scal,int Nvars,double *cpmix) ;

/* T, Y_i -> cv [J/(kg.K)] */
int TCDND_getMs2CvMixMs(double *scal,int Nvars,double *cvmix) ;
int TC_getMs2CvMixMs   (double *scal,int Nvars,double *cvmix) ;

/* T, X_i -> Cp [J/(kmol.K)] */
int TCDND_getMl2CpMixMl(double *scal,int Nvars,double *cpmix) ;
int TC_getMl2CpMixMl   (double *scal,int Nvars,double *cpmix) ;

/* T, X_i -> Cv [J/(kmol.K)] */
int TC_getMl2CvMixMl(double *scal,int Nvars,double *cvmix) ;

/* T -> cp_i [J/(kg.K)] */
int TCDND_getCpSpecMs(double t,int Nspec,double *cpi) ;
int TC_getCpSpecMs   (double t,int Nspec,double *cpi) ;

/* T -> cv_i [J/(kg.K)] */
int TC_getCvSpecMs   (double t,int Nspec,double *cvi) ;

/* T -> Cp_i [J/(kmol.K)] */
int TCDND_getCpSpecMl(double t,int Nspec,double *cpi) ;
int TC_getCpSpecMl   (double t,int Nspec,double *cpi) ;

/* T -> Cv_i [J/(kmol.K)] */
int TC_getCvSpecMl   (double t,int Nspec,double *cvi) ;

/* T, Y_i -> h [J/kg] */
int TCDND_getMs2HmixMs(double *scal,int Nvars,double *hmix) ; 
int TC_getMs2HmixMs   (double *scal,int Nvars,double *hmix) ; 

/* T, X_i -> h [J/kmol] */
int TCDND_getMl2HmixMl(double *scal,int Nvars,double *hmix) ;
int TC_getMl2HmixMl   (double *scal,int Nvars,double *hmix) ;

/* T -> h_i [J/kg] */
int TCDND_getHspecMs(double t,int Nspec,double *hi) ;
int TC_getHspecMs   (double t,int Nspec,double *hi) ;

/* T -> u_i [J/kg] */
int TC_getUspecMs   (double t,int Nspec,double *ui) ;

/* T -> H_i [J/kmol] */
int TCDND_getHspecMl(double t,int Nspec,double *hi) ;
int TC_getHspecMl   (double t,int Nspec,double *hi) ;

/* T -> U_i [J/kmol] */
int TC_getUspecMl   (double t,int Nspec,double *ui) ;

/*
      _____ ____               _                          
     |_   _/ ___|    _ __ ___ | |_ __ ___  ___        ___ 
       | || |       | '_ ` _ \| | '_ ` _ \/ __|      / __|
       | || |___    | | | | | | | | | | | \__ \  _  | (__ 
       |_| \____|___|_| |_| |_|_|_| |_| |_|___/ (_)  \___|
               |_____|                                    
  
*/
/* T, Y_i -> C_i [kmol/m3] */
int TCDND_getMs2Cc(double *scal,int Nvars,double *concX) ;
int TC_getMs2Cc   (double *scal,int Nvars,double *concX) ;

/* X_i -> Y_i */
int TCDND_getMl2Ms(double *Xspec,int Nspec,double *Yspec) ;
int TC_getMl2Ms   (const double *Xspec,int Nspec,double *Yspec) ;

/* Y_i -> X_i */
int TCDND_getMs2Ml(double *Yspec,int Nspec,double *Xspec) ;
int TC_getMs2Ml   (const double *Yspec,int Nspec,double *Xspec) ;

/* Y_i -> Wmix [kg/kmol] */
int TCDND_getMs2Wmix(double *Yspec,int Nspec,double *Wmix) ;
int TC_getMs2Wmix   (const double *Yspec,int Nspec,double *Wmix) ;

/* X_i -> Wmix [kg/kmol] */
int TCDND_getMl2Wmix(double *Xspec,int Nspec,double *Wmix) ;
int TC_getMl2Wmix   (const double *Xspec,int Nspec,double *Wmix) ;

/*
   _____ ____                        
  |_   _/ ___|    ___ _ __ ___   ___ 
    | || |       / __| '__/ __| / __|
    | || |___    \__ \ | | (__ | (__ 
    |_| \____|___|___/_|  \___(_)___|
            |_____|                  

*/

/* non-conservative src term given temperature and species mass fractions */
int TCDND_getSrc(double *scal,int Nvars,double *omega) ;
int TC_getSrc   (double *scal,int Nvars,double *omega) ;
int TC_getSrcCV (double *scal,int Nvars,double *omega) ;

/* Smat */
int TC_getSmat  (double *scal, int Nvars, double *Smat) ;
int TC_getSmatCV(double *scal, int Nvars, double *Smat) ;

/* conservative src term given density and species mass fractions */
int TCDND_getSrcCons(double *scal,int Nvars,double *omega) ;
int TC_getSrcCons   (double *scal,int Nvars,double *omega) ;

/* jacobian for non-cons system: T, Y_1, Y_2,..., Y_{N-1} */
int TCDND_getJacTYNm1anl(double *scal, int Nspec, double *jac) ;
int TC_getJacTYNm1anl   (double *scal, int Nspec, double *jac) ;
int TCDND_getJacTYNm1   (double *scal, int Nspec, double *jac, unsigned int useJacAnl) ;
int TC_getJacTYNm1      (double *scal, int Nspec, double *jac, unsigned int useJacAnl) ;

/* jacobian for non-cons system: T, Y_1, Y_2,..., Y_{N} */
int TCDND_getJacTYNanl(double *scal, int Nspec, double *jac) ;
int TC_getJacTYNanl   (double *scal, int Nspec, double *jac) ;
int TCDND_getJacTYN   (double *scal, int Nspec, double *jac, unsigned int useJacAnl) ;
int TC_getJacTYN      (double *scal, int Nspec, double *jac, unsigned int useJacAnl) ;

/* jacobian for cons system: rho, P, T, Y_1, Y_2,..., Y_{N} */
int TC_getJacRPTYN   (double *scal, int Nspec, double *jac, unsigned int useJacAnl) ;
int TC_getJacRPTYNanl(double *scal, int Nspec, double *jac) ;
int TC_getJacRPTYNnum(double *scal, int Nspec, double *jac) ;

/* jacobian for constant volume */
int TC_getJacCVTYNanl(double *scal, int Nspec, double *jac);
int TC_getJacCVTYNnum(double *scal, int Nspec, double *jac);
  
/* Eliminate species and reactions from a kinetic model */
void TC_reduce(char *mechIn, char *slist, char *mechOut);

void TC_splitFR(char *minp, int *rlist) ;
int  TC_parseFG(char *minp);

#ifdef __cplusplus
} /* End of extern "C" */
#endif

#endif /* end if #ifndef _TCinterfHSeen_  */
