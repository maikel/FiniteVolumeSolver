#include "TC_defs.h"
#include "TC_interface.h"

#include <string.h>

int tcsplitfr_(char *mechfile, int *lmech, int *rlist) {
  memset(&mechfile  [*lmech], 0, 1) ;
  TC_splitFR(mechfile,rlist);
  return (0);
}

/*
    _____ ____     _        __
   |_   _/ ___|   (_)_ __  / _| ___
     | || |       | | '_ \| |_ / _ \
     | || |___    | | | | |  _| (_) |
     |_| \____|___|_|_| |_|_|  \___/
             |_____|

*/
int tcgetnreac_() { return ( TC_getNreac() ) ; }
int tcgetnpresdepreac_() { return ( TC_getNPresDepReac() ) ; }
int tcgetnthbreac_() { return ( TC_getNThbReac() ) ; }

int tcgetstoicoef_( int *Nspec, int *Nreac, double *stoicoef ) {
  return (TC_getStoiCoef( *Nspec, *Nreac, stoicoef ) ) ;
}
int tcgetstoicoeffr_( int *Nspec, int *Nreac, double *stoicoef ) {
  return ( TC_getStoiCoefFR( *Nspec, *Nreac, stoicoef ) ) ;
}
int tcgetstoicoefreac_( int *Nspec, int *Nreac, int *ireac, int *idx, double *stoicoef ) {
  return ( TC_getStoiCoefReac( *Nspec, *Nreac, *ireac, *idx, stoicoef ) ) ;
}
int tcgetstoicoefij_( int *Nspec, int *Nreac, int *ireac, int *kspec, double *stoicoef ) {
  return ( TC_getStoiCoefIJ( *Nspec, *Nreac, *ireac, *kspec, stoicoef ) );
}

int tcgetreacrev_(int *ireac)  { return ( TC_getReacRev(*ireac)  ) ; }
int tcgetreacthb_(int *ireac)  { return ( TC_getReacThb(*ireac)  ) ; }
int tcgetreaclow_(int *ireac)  { return ( TC_getReacLow(*ireac)  ) ; }
int tcgetreachigh_(int *ireac) { return ( TC_getReacHIGH(*ireac) ) ; }

int tcgetreacptype_(int *ireac) { return ( TC_getReacPtype(*ireac) ) ; }
int tcgetthbreaceff_(int *ireac,int *specidx, double *effs) {
  return ( TC_getThbReacEff(*ireac,specidx,effs) ) ;
}

int tcgetteacthbidx_   (int *ireac) { return ( TC_getReacThbIdx(*ireac)    ) ; }
int tcgetreacfallidx_  (int *ireac) { return ( TC_getReacFallIdx(*ireac)   ) ; }
int tcgetreacrealstidx_(int *ireac) { return ( TC_getReacRealStIdx(*ireac) ) ; }

int tcgetreacstr_(int *ireac,int *lenstr, char *rstr) {
  int i, ans = 0 ;
  ans = TC_getReacStr(*ireac,*lenstr,rstr)  ;
  for ( i=0; i<*lenstr; i++)
    if ( rstr[i] == 0 ) rstr[i] = 32;
  return ( ans ) ;
}

/*

       _____ ____         _
      |_   _/ ___|    ___| |__   __ _
        | || |       / __| '_ \ / _` |
        | || |___   | (__| | | | (_| |
        |_| \____|___\___|_| |_|\__, |
                |_____|         |___/

*/
int tcchgarhenfor_( int *ireac, int *ipos, double *newval )
{
  int ans = 0 ;
  return ( ans ) ;
  ans =  TC_chgArhenFor( *ireac, *ipos, *newval ) ;
}

int tcchgarhenforback_( int *ireac, int *ipos)
{
  int ans = 0 ;
  ans =  TC_chgArhenForBack( *ireac, *ipos) ;
  return ( ans ) ;
}

int tcchgarhenrev_( int *ireac, int *ipos, double *newval )
{
  int ans = 0 ;
  ans =  TC_chgArhenRev( *ireac, *ipos, *newval ) ;
  return ( ans ) ;
}

int tcchgarhenrevback_( int *ireac, int *ipos)
{
  int ans = 0 ;
  ans =  TC_chgArhenRevBack( *ireac, *ipos) ;
  return ( ans ) ;
}
int tcchgarhenpresdepfor_( int *ireac, int *ipos, double *newval )
{
  int ans = 0 ;
  return ( ans ) ;
  ans =  TC_chgArhenPresDepFor( *ireac, *ipos, *newval ) ;
}

/*
         _____ ____             _ _ _
        |_   _/ ___|    ___  __| (_) |_
          | || |       / _ \/ _` | | __|
          | || |___   |  __/ (_| | | |_
          |_| \____|___\___|\__,_|_|\__|
                  |_____|
*/
void tcreset_()
{
  TC_reset() ;
  return ;
}

/*
             _____ ____     _       _ _
            |_   _/ ___|   (_)_ __ (_) |_   ___
              | || |       | | '_ \| | __| / __|
              | || |___    | | | | | | |_ | (__
              |_| \____|___|_|_| |_|_|\__(_)___|
                      |_____|


*/
int tcinitchem_(char *mechfile,int *lmech,char *thermofile,int *lthrm,
                int *tab, double *delT)
{
  int ans = 0 ;
  /* add a NULL at the ends of mechfile and thermofile */
  memset(&mechfile  [*lmech], 0, 1) ;
  memset(&thermofile[*lthrm], 0, 1) ;
  ans = TC_initChem(mechfile,thermofile,*tab,*delT);
  return ( ans ) ;
}

void tcsetrefval_(double *rhoref, double *pref, double *Tref, double *Wref,
                  double *Daref, double *omgref, double *cpref, double *href,
                  double *timref)
{
  TC_setRefVal(*rhoref, *pref, *Tref, *Wref, *Daref,
	       *omgref, *cpref, *href, *timref);
  return ;
}

void tcsetnondim_()
{
  TC_setNonDim()  ;
  return ;
}

void tcsetdim_()
{
  TC_setDim()  ;
  return ;
}

void tcsetthermopres_(double *pressure)
{
  TC_setThermoPres(*pressure) ;
  return ;
}

void tcsetdens_(double *density)
{
  TC_setDens(*density) ;
  return ;
}

/*
      _____ ____               _
     |_   _/ ___|    _ __ ___ | |_ __ ___  ___        ___
       | || |       | '_ ` _ \| | '_ ` _ \/ __|      / __|
       | || |___    | | | | | | | | | | | \__ \  _  | (__
       |_| \____|___|_| |_| |_|_|_| |_| |_|___/ (_)  \___|
               |_____|

*/
int tcdndgetms2cc_(double *scal,int *Nvars,double *concX)
{
  int ans = 0 ;
  ans = TCDND_getMs2Cc( scal, (*Nvars), concX) ;
  return ( ans ) ;
}

int tcgetms2cc_(double *scal,int *Nvars,double *concX)
{
  int ans = 0 ;
  ans = TC_getMs2Cc( scal, (*Nvars), concX) ;
  return ( ans ) ;
}

int tcdndgetml2ms_(double *Xspec,int *Nspec,double *Yspec)
{
  int ans = 0 ;
  ans = TCDND_getMl2Ms( Xspec, (*Nspec), Yspec) ;
  return ( ans ) ;
}

int tcgetml2ms_(double *Xspec,int *Nspec,double *Yspec)
{
  int ans = 0 ;
  ans = TC_getMl2Ms( Xspec, (*Nspec), Yspec) ;
  return ( ans ) ;
}

int tcdndgetms2ml_(double *Yspec,int *Nspec,double *Xspec)
{
  int ans = 0 ;
  ans = TCDND_getMs2Ml( Yspec, (*Nspec), Xspec) ;
  return ( ans ) ;
}

int tcgetms2ml_(double *Yspec,int *Nspec,double *Xspec)
{
  int ans = 0 ;
  ans = TC_getMs2Ml( Yspec, (*Nspec), Xspec) ;
  return ( ans ) ;
}

int tcdndgetms2wmix_(double *Yspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TCDND_getMs2Wmix(Yspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcgetms2wmix_(double *Yspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TC_getMs2Wmix(Yspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcdndgetml2wmix_(double *Xspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TCDND_getMl2Wmix(Xspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

int tcgetml2wmix_(double *Xspec,int *Nspec,double *Wmix)
{
  int ans = 0 ;
  ans = TC_getMl2Wmix(Xspec,(*Nspec),Wmix) ;
  return ( ans ) ;
}

/*
          _____ ____
         |_   _/ ___|    _ __ _ __        ___
           | || |       | '__| '__|      / __|
           | || |___    | |  | |     _  | (__
           |_| \____|___|_|  |_|    (_)  \___|
                   |_____|
*/

int tcgetarhenfor_( int *ireac, int *ipos, double *val )
{
  int ans = 0 ;
  ans = TC_getArhenFor( *ireac, *ipos, val ) ;
  return ( ans );
}

int tcgetarhenrev_( int *ireac, int *ipos, double *val )
{
  int ans = 0 ;
  ans = TC_getArhenRev( *ireac, *ipos, val ) ;
  return ( ans );
}

int tcdndgetty2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTY2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetty2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTY2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgetty2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTY2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetty2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTY2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgettxc2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTXC2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgettxc2rrml_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTXC2RRml( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcdndgettxc2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TCDND_getTXC2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgettxc2rrms_(double *scal, int *Nvars, double *omega)
{
  int ans = 0 ;
  ans = TC_getTXC2RRms( scal, *Nvars, omega ) ;
  return ( ans );
}

int tcgetrops_(double *scal, int *Nvars, double *datarop)
{
  int ans = 0 ;
  ans = TC_getRops( scal, *Nvars, datarop) ;
  return ( ans );
}

int tcgetrfrb_(double *scal, int *Nvars, double *dataRfrb)
{
  int ans = 0 ;
  ans = TC_getRfrb(scal, *Nvars, dataRfrb) ;
  return ( ans );
}
/*
          _____ ____
         |_   _/ ___|    ___ _ __   ___  ___   ___
           | || |       / __| '_ \ / _ \/ __| / __|
           | || |___    \__ \ |_) |  __/ (__ | (__
           |_| \____|___|___/ .__/ \___|\___(_)___|
                   |_____|  |_|

*/
int tcgetnvars_() { return ( TC_getNvars() ) ; }

int tcgetnelem_() { return ( TC_getNelem() ) ; }

int tcgetnspec_() { return ( TC_getNspec() ) ; }

int tcgetenamelen_() { return ( TC_getEnameLen() ) ; }

int tcgetsnamelen_() { return ( TC_getSnameLen() ) ; }

int tcgetenames_(int *Nelem,char *enames)
{
  int i, ans = 0 ;
  ans = TC_getEnames( (*Nelem), enames);
  for ( i=0; i<LENGTHOFELEMNAME * TC_Nelem_; i++)
    if ( enames[i] == 0 ) enames[i] = 32;
  return ( ans ) ;
}

int tcgetsnames_(int *Nspec,char *snames)
{
  int i, ans = 0 ;
  ans = TC_getSnames( (*Nspec), snames);
  for ( i=0; i<LENGTHOFSPECNAME * TC_Nspec_; i++)
    if ( snames[i] == 0 ) snames[i] = 32;
  return ( ans ) ;
}

int tcgetspos_(const char *sname, const int *slen) {
  return ( TC_getSpos( sname, *slen) ) ;
}
int tcgetsmass_(int *Nspec,double *Wi) { return ( TC_getSmass( *Nspec, Wi) ); }

int tcgetecount_(int *elemcnt) { return ( TC_getECount(elemcnt) ); }


/*
   _____ ____
  |_   _/ ___|    ___ _ __ ___   ___
    | || |       / __| '__/ __| / __|
    | || |___    \__ \ | | (__ | (__
    |_| \____|___|___/_|  \___(_)___|
            |_____|

*/
int tcdndgetsrc_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TCDND_getSrc(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsrc_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TC_getSrc(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsrccv_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TC_getSrcCV(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsmat_(double *scal,int *Nvars,double *Smat)
{
  int ans = 0 ;
  ans = TC_getSmat(scal,*Nvars,Smat) ;
  return ( ans );
}

int tcgetsmatcv_(double *scal,int *Nvars,double *Smat)
{
  int ans = 0 ;
  ans = TC_getSmatCV(scal,*Nvars,Smat) ;
  return ( ans );
}

int tcdndgetsrccons_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TCDND_getSrcCons(scal,*Nvars,omega) ;
  return ( ans );
}

int tcgetsrccons_(double *scal,int *Nvars,double *omega)
{
  int ans = 0 ;
  ans = TC_getSrcCons(scal,*Nvars,omega) ;
  return ( ans );
}

int tcdndgetjactynm1anl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TCDND_getJacTYNm1anl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjactynm1anl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacTYNm1anl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcdndgetjactynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TCDND_getJacTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjactynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjaccvtynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacCVTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcdndgetjactynm1_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TCDND_getJacTYNm1(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjactynm1_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacTYNm1(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcdndgetjactyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TCDND_getJacTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjactyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjacrptyn_(double *scal, int *Nspec, double *jac, int *useJacAnl)
{
  int ans = 0 ;
  unsigned int jactype = (unsigned int) *useJacAnl;
  ans = TC_getJacRPTYN(scal,(*Nspec),jac,jactype) ;
  return ( ans ) ;
}

int tcgetjacrptynanl_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacRPTYNanl(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

int tcgetjacrptynnum_(double *scal, int *Nspec, double *jac)
{
  int ans = 0 ;
  ans = TC_getJacRPTYNnum(scal,(*Nspec),jac) ;
  return ( ans ) ;
}

/*
      _____ ____    _   _
     |_   _/ ___|  | |_| |__   ___ _ __ _ __ ___   ___         ___
       | || |      | __| '_ \ / _ \ '__| '_ ` _ \ / _ \       / __|
       | || |___   | |_| | | |  __/ |  | | | | | | (_) |  _  | (__
       |_| \____|___\__|_| |_|\___|_|  |_| |_| |_|\___/  (_)  \___|
               |_____|
*/

/* T, Y_i -> rho [?] */
int tcdndgetrhomixms_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TCDND_getRhoMixMs( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

/* T, Y_i -> rho [kg/m3] */
int tcgetrhomixms_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TC_getRhoMixMs( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

/* T, X_i -> rho [?] */
int tcdndgetrhomixml_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TCDND_getRhoMixMl( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

/* T, X_i -> rho [kg/m3] */
int tcgetrhomixml_(double *scal,int *Nvars,double *rhomix)
{
  int ans = 0 ;
  ans = TC_getRhoMixMl( scal, (*Nvars), rhomix) ;
  return ( ans ) ;
}

/* rho[kg/m3], Y_i -> T [?] */
int tcdndgettmixms_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TCDND_getTmixMs( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

/* rho[kg/m3], Y_i -> T [K] */
int tcgettmixms_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TC_getTmixMs( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

/* rho[kg/m3], X_i -> T [?] */
int tcdndgettmixml_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TCDND_getTmixMl( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

/* rho[kg/m3], X_i -> T [K] */
int tcgettmixml_(double *scal,int *Nvars,double *Tmix)
{
  int ans = 0 ;
  ans = TC_getTmixMl( scal, (*Nvars), Tmix) ;
  return ( ans ) ;
}

/* T, Y_i -> cp [?] */
int tcdndgetms2cpmixms_(double *scal,int *Nvars,double *cpmix)
{
  int ans = 0 ;
  ans = TCDND_getMs2CpMixMs( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

/* T, Y_i -> cp [J/(kg.K)] */
int tcgetms2cpmixms_(double *scal,int *Nvars,double *cpmix)
{
  int ans = 0 ;
  ans = TC_getMs2CpMixMs( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

/* T, Y_i -> S [J/(kg.K)] */
int tcgetms2entrmixms_(double *scal,int *Nvars,double *smix)
{
  int ans = 0 ;
  ans = TC_getMs2EntrMixMs( scal, (*Nvars), smix) ;
  return ( ans ) ;
}

/* T, Y_i -> cv [?] */
int tcdndgetms2cvmixms_(double *scal,int *Nvars,double *cvmix)
{
  int ans = 0 ;
  ans = TCDND_getMs2CvMixMs( scal, (*Nvars), cvmix) ;
  return ( ans ) ;
}

/* T, Y_i -> cv [J/(kg.K)] */
int tcgetms2cvmixms_(double *scal,int *Nvars,double *cvmix)
{
  int ans = 0 ;
  ans = TC_getMs2CvMixMs( scal, (*Nvars), cvmix) ;
  return ( ans ) ;
}

/* T, X_i -> Cp [?] */
int tcdndgetml2cpmixml_(double *scal,int *Nvars,double *cpmix)
{
  int ans = 0 ;
  ans = TCDND_getMl2CpMixMl( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

/* T, X_i -> Cp [J/(kmol.K)] */
int tcgetml2cpmixml_(double *scal,int *Nvars,double *cpmix)
{
  int ans = 0 ;
  ans = TC_getMl2CpMixMl( scal, (*Nvars), cpmix) ;
  return ( ans ) ;
}

/* T, X_i -> S [J/(kmol.K)] */
int tcgetml2entrmixml_(double *scal,int *Nvars,double *smix)
{
  int ans = 0 ;
  ans = TC_getMl2EntrMixMl( scal, (*Nvars), smix) ;
  return ( ans ) ;
}

/* T, X_i -> Cv [J/(kmol.K)] */
int tcgetml2cvmixml_(double *scal,int *Nvars,double *cvmix)
{
  int ans = 0 ;
  ans = TC_getMl2CvMixMl( scal, (*Nvars), cvmix) ;
  return ( ans ) ;
}

/* T -> cp_i [?] */
int tcdndgetcpspecms_(double *t,int *Nspec,double *cpi)
{
  int ans = 0 ;
  ans = TCDND_getCpSpecMs( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

/* T -> cp_i [J/(kg.K)] */
int tcgetcpspecms_(double *t,int *Nspec,double *cpi)
{
  int ans = 0 ;
  ans = TC_getCpSpecMs( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

/* T -> s_i^0 [J/(kg.K)] */
int tcgetentr0specms_(double *t, double *si0)
{
  int ans = 0 ;
  ans = TC_getEntr0SpecMsFcn( *t, si0) ;
  return ( ans ) ;
}

/* T -> cv_i [J/(kg.K)] */
int tcgetcvspecms_(double *t,int *Nspec,double *cvi)
{
  int ans = 0 ;
  ans = TC_getCvSpecMs( *t, (*Nspec), cvi) ;
  return ( ans ) ;
}

/* T -> Cp_i [?] */
int tcdndgetcpspecml_(double *t,int *Nspec,double *cpi)
{
  int ans = 0 ;
  ans = TCDND_getCpSpecMl( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

/* T -> Cp_i [J/(kmol.K)] */
int tcgetcpspecml_(double *t,int *Nspec,double *cpi)
{
  int ans = 0 ;
  ans = TC_getCpSpecMl( *t, (*Nspec), cpi) ;
  return ( ans ) ;
}

/* T -> S_i^0 [J/(kmol.K)] */
int tcgetentr0specml_(double *t, double *si0)
{
  int ans = 0 ;
  ans = TC_getEntr0SpecMlFcn( *t, si0) ;
  return ( ans ) ;
}

/* T -> Cv_i [J/(kmol.K)] */
int tcgetcvspecml_(double *t,int *Nspec,double *cvi)
{
  int ans = 0 ;
  ans = TC_getCvSpecMl( *t, (*Nspec), cvi) ;
  return ( ans ) ;
}

/* T, Y_i -> h [?] */
int tcdndgetms2hmixms_(double *scal,int *Nvars,double *hmix)
{
  int ans = 0 ;
  ans = TCDND_getMs2HmixMs( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

/* T, Y_i -> h [J/kg] */
int tcgetms2hmixms_(double *scal,int *Nvars,double *hmix)
{
  int ans = 0 ;
  ans = TC_getMs2HmixMs( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

/* T, X_i -> h [?] */
int tcdndgetml2hmixml_(double *scal,int *Nvars,double *hmix)
{
  int ans = 0 ;
  ans = TCDND_getMl2HmixMl( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

/* T, X_i -> h [J/kmol] */
int tcgetml2hmixml_(double *scal,int *Nvars,double *hmix)
{
  int ans = 0 ;
  ans = TC_getMl2HmixMl( scal, (*Nvars), hmix) ;
  return ( ans ) ;
}

/* T -> h_i [?] */
int tcdndgethspecms_(double *t,int *Nspec,double *hi)
{
  int ans = 0 ;
  ans = TCDND_getHspecMs( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

/* T -> h_i [J/kg] */
int tcgethspecms_(double *t,int *Nspec,double *hi)
{
  int ans = 0 ;
  ans = TC_getHspecMs( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

/* T -> u_i [J/kg] */
int tcgetuspecms_(double *t,int *Nspec,double *ui)
{
  int ans = 0 ;
  ans = TC_getUspecMs( *t, (*Nspec), ui) ;
  return ( ans ) ;
}

/* T -> H_i [?] */
int tcdndgethspecml_(double *t,int *Nspec,double *hi)
{
  int ans = 0 ;
  ans = TCDND_getHspecMl( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

/* T -> H_i [J/kmol] */
int tcgethspecml_(double *t,int *Nspec,double *hi)
{
  int ans = 0 ;
  ans = TC_getHspecMl( *t, (*Nspec), hi) ;
  return ( ans ) ;
}

/* T -> U_i [J/kmol] */
int tcgetuspecml_(double *t,int *Nspec,double *ui)
{
  int ans = 0 ;
  ans = TC_getUspecMl( *t, (*Nspec), ui) ;
  return ( ans ) ;
}
