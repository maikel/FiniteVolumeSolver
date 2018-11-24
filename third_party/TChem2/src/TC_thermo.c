#include "TC_defs.h"
#include "TC_interface.h"

#include <math.h>

double TC_getRUniv() {
  return TC_Runiv_;
}

/*! \file TC_thermo.c
    \brief Equation of state and thermodynamic functions
*/
/*
                  _   ____  _           __  __ _
        __ _  ___| |_|  _ \| |__   ___ |  \/  (_)_  __
       / _` |/ _ \ __| |_) | '_ \ / _ \| |\/| | \ \/ /
      | (_| |  __/ |_|  _ <| | | | (_) | |  | | |>  <
       \__, |\___|\__|_| \_\_| |_|\___/|_|  |_|_/_/\_\
       |___/

*/
/**
 * \ingroup eqstate
 * \brief Computes density based on temperature and species mass fractions using the equation of state.
 * Input temperature is normalized, output density also normalized before exit.
*/
int TCDND_getRhoMixMs(double *scal,int Nvars,double *rhomix)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>),
                 temperature T [K], mass fractions Y []
   \param Nvars : no. of variables = N<sub>spec</sub>+1
   \return rhomix : pointer to mixture density [kg/m<sup>3</sup>]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getRhoMixMs", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getRhoMixMs(scal, Nvars, rhomix) ;

  if ( TC_nonDim_ == 1 )
  {
    *rhomix /= TC_rhoref_ ;
    scal[0] /= TC_Tref_   ;
  }

  return ( ans ) ;

}
/**
 * \ingroup eqstate
 * \brief Computes density based on temperature and species mass fractions using the equation of state.
 */
int TC_getRhoMixMs(double *scal,int Nvars,double *rhomix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>),
                temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return rhomix : pointer to mixture density [kg/m<sup>3</sup>]
*/

  double temperature, *Yspec, sumYoW ;
  int i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getRhoMixMs", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Yspec       = &scal[1] ;

  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i]/TC_sMass_[i] ;
  /*
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) printf("TC_getRhoMixMs(): %d %e %e\n", i, Yspec[i],TC_sMass_[i]) ;
  */
  *rhomix = TC_pressure_/(TC_Runiv_*sumYoW*temperature) ;

  return ( 0 ) ;

}

/**
 * \ingroup eqstate
 * \brief Computes density based on temperature and species mole fractions using the equation of state.
 * Input temperature is normalized, output density also normalized before exit.
 */
int TCDND_getRhoMixMl(double *scal,int Nvars,double *rhomix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>),
                temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return rhomix : pointer to mixture density [kg/m<sup>3</sup>]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getRhoMixMl", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getRhoMixMl(scal, Nvars, rhomix) ;

  if ( TC_nonDim_ == 1 )
  {
    *rhomix /= TC_rhoref_ ;
    scal[0] /= TC_Tref_   ;
  }

  return ( ans ) ;

}

/**
 * \ingroup eqstate
 * \brief Computes density based on temperature and species mole fractions using the equation of state.
 */
int TC_getRhoMixMl(double *scal,int Nvars,double *rhomix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>),
                temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return rhomix : pointer to mixture density [kg/m<sup>3</sup>]
*/

  double temperature, *Xspec, Wmix ;
  int i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getRhoMixMl", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Xspec       = &scal[1] ;

  Wmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) Wmix += Xspec[i]*TC_sMass_[i] ;

  *rhomix = (TC_pressure_*Wmix)/(TC_Runiv_*temperature) ;

  return ( 0 ) ;

}
/*
                      _  _____ __  __ _
            __ _  ___| ||_   _|  \/  (_)_  __
           / _` |/ _ \ __|| | | |\/| | \ \/ /
          | (_| |  __/ |_ | | | |  | | |>  <
           \__, |\___|\__||_| |_|  |_|_/_/\_\
           |___/

*/
/**
 * \ingroup eqstate
 * \brief Computes temperature based on density and species mass fractions using the equation of state.
 *        Input density is normalized, output temperature also normalized before exit.
 */
int TCDND_getTmixMs(double *scal,int Nvars,double *Tmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (rho,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>),
                density rho [kg/m<sup>3</sup>], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return Tmix : pointer to temperature [K]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTmixMs", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_rhoref_ ;

  ans = TC_getTmixMs(scal, Nvars, Tmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *Tmix   /= TC_Tref_   ;
    scal[0] /= TC_rhoref_ ;
  }

  return ( ans ) ;

}
/**
 * \ingroup eqstate
 *  \brief Computes temperature based on density and species mass fractions using the equation of state.
 */
int TC_getTmixMs(double *scal,int Nvars,double *Tmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (rho,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), density rho [kg/m<sup>3</sup>], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return Tmix : pointer to temperature [K]
*/

  double rho, *Yspec, sumYoW ;
  int i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTmixMs", Nvars, TC_Nvars_ ) ;

  rho   =  scal[0] ;
  Yspec = &scal[1] ;

  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i]/TC_sMass_[i] ;

  *Tmix = TC_pressure_/(TC_Runiv_*sumYoW*rho) ;

  return ( 0 ) ;

}
/**
 * \ingroup eqstate
 * \brief Computes temperature based on density and species mole fractions using the equation of state.
 *        Input density is normalized, output temperature also normalized before exit.
 */
int TCDND_getTmixMl(double *scal,int Nvars,double *Tmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (rho,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), density rho [kg/m<sup>3</sup>], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return Tmix : pointer to temperature [K]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTmixMl", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_rhoref_ ;

  ans = TC_getTmixMl(scal, Nvars, Tmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *Tmix   /= TC_Tref_   ;
    scal[0] /= TC_rhoref_ ;
  }

  return ( ans ) ;

}
/**
 * \ingroup eqstate
 * \brief Computes temperature based on density and species mole fractions using the equation of state.
 */
int TC_getTmixMl(double *scal,int Nvars,double *Tmix)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (rho,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), density rho [kg/m<sup>3</sup>], mole fractions X []
   \param Nvars : no. of variables = N<sub>spec</sub>+1
   \return Tmix : pointer to temperature [K]
*/

  double rho, *Xspec, Wmix ;
  int i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTmixMl", Nvars, TC_Nvars_ ) ;

  rho   =  scal[0] ;
  Xspec = &scal[1] ;

  Wmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) Wmix += Xspec[i]*TC_sMass_[i] ;

  *Tmix = ( TC_pressure_ * Wmix ) / ( TC_Runiv_ * rho ) ;

  return ( 0 ) ;

}
/*
                _   __  __     ____   ____       __  __ _
      __ _  ___| |_|  \/  |___|___ \ / ___|_ __ |  \/  (_)_  _
     / _` |/ _ \ __| |\/| / __| __) | |   | '_ \| |\/| | \ \/ /
    | (_| |  __/ |_| |  | \__ \/ __/| |___| |_) | |  | | |>  <|
     \__, |\___|\__|_|  |_|___/_____|\____| .__/|_|  |_|_/_/\_\
     |___/                                |_|

*/
/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant pressure based on temperature and
   species mass fractions. Input temperature is normalized, output specific heat is also normalized before exit.
*/
int TCDND_getMs2CpMixMs(double *scal,int Nvars,double *cpmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cpmix : pointer to mixture specific heat at constant pressure [J/(kg.K)]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getMs2CpMixMs", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMs2CpMixMs(scal,Nvars,cpmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *cpmix  /= TC_cpref_ ;
    scal[0] /= TC_Tref_  ;
  }

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant pressure based on temperature and
         species mass fractions.
*/
int TC_getMs2CpMixMs(double *scal,int Nvars,double *cpmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>),
                temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cpmix : pointer to mixture specific heat at constant pressure [J/(kg.K)]
*/

  int ans, i ;
  double temperature, *Yspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMs2CpMixMs", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Yspec       = &scal[1] ;

  /* compute species cp's */
  ans = TC_getCpSpecMs(temperature,TC_Nspec_,TC_cpks) ;

  *cpmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *cpmix += Yspec[i]*TC_cpks[i] ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant volume based on temperature and
   species mass fractions. Input temperature is normalized, output specific heat is also normalized before exit.
*/
int TCDND_getMs2CvMixMs(double *scal,int Nvars,double *cvmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cpmix : pointer to mixture specific heat at constant volume [J/(kg.K)]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getMs2CvMixMs", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMs2CvMixMs(scal,Nvars,cvmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *cvmix  /= TC_cpref_ ;
    scal[0] /= TC_Tref_  ;
  }

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant volume based on temperature and
         species mass fractions.
*/
int TC_getMs2CvMixMs(double *scal,int Nvars,double *cvmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>),
                temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cvmix : pointer to mixture specific heat at constant volume [J/(kg.K)]
*/

  int ans, i ;
  double *Yspec, sumYoW, cpmix ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMs2CvMixMs", Nvars, TC_Nvars_ ) ;

  /* Compute molar weight */
  Yspec = &scal[1] ;
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i]/TC_sMass_[i] ;

  /* get cp first */
  ans = TC_getMs2CpMixMs(scal, Nvars, &cpmix) ;

  *cvmix = cpmix - sumYoW*TC_Runiv_ ;

  return ( ans ) ;

}
/*
                  _   __  __ _ ____   ____       __  __ _
        __ _  ___| |_|  \/  | |___ \ / ___|_ __ |  \/  (_)_  __
       / _` |/ _ \ __| |\/| | | __) | |   | '_ \| |\/| | \ \/ /
      | (_| |  __/ |_| |  | | |/ __/| |___| |_) | |  | | |>  <
       \__, |\___|\__|_|  |_|_|_____|\____| .__/|_|  |_|_/_/\_\
       |___/                              |_|

*/
/**
 * \ingroup thermo
 * \brief Computes mixture heat capacity at constant pressure based on temperature and
         species mole fractions. Input temperature is normalized, output specific heat
         is also normalized before exit.
*/
int TCDND_getMl2CpMixMl(double *scal,int Nvars,double *cpmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cpmix : pointer to mixture heat capacity at constant pressure [J/(kmol.K)]
*/

  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getMl2CpMixMl", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMl2CpMixMl(scal,Nvars,cpmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *cpmix  /= ( TC_cpref_ * TC_Wref_ ) ;
    scal[0] /= TC_Tref_  ;
  }

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant pressure based on temperature and
         species mole fractions.
*/
int TC_getMl2CpMixMl(double *scal,int Nvars,double *cpmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cpmix : pointer to mixture specific heat at constant pressure [J/(kmol.K)]
*/

  int ans, i ;
  double temperature, *Xspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMl2CpMixMl", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Xspec       = &scal[1] ;

  /* compute species cp's */
  ans = TC_getCpSpecMl(temperature,TC_Nspec_,TC_cpks) ;

  *cpmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *cpmix += Xspec[i]*TC_cpks[i] ;

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture specific heat at constant volume based on temperature and
         species mole fractions.
*/
int TC_getMl2CvMixMl(double *scal,int Nvars,double *cvmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return cvmix : pointer to mixture specific heat at constant volume [J/(kmol.K)]
*/

  int ans, i ;
  double temperature, *Xspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMl2CvMixMl", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Xspec       = &scal[1] ;

  /* compute species cp's */
  ans = TC_getCpSpecMl(temperature,TC_Nspec_,TC_cpks) ;

  *cvmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *cvmix += Xspec[i]*TC_cpks[i] ;
  *cvmix -= TC_Runiv_ ;

  return ( ans ) ;

}
/*
                   _    ____      ____
         __ _  ___| |_ / ___|_ __/ ___| _ __   ___  ___
        / _` |/ _ \ __| |   | '_ \___ \| '_ \ / _ \/ __|
       | (_| |  __/ |_| |___| |_) |__) | |_) |  __/ (__
        \__, |\___|\__|\____| .__/____/| .__/ \___|\___|
        |___/               |_|        |_|

*/
/**
 * \ingroup thermo
 * \brief Computes species specific heat at constant pressure based on temperature.
         Input temperature is normalized, output specific heats are also normalized before exit.
*/
int TCDND_getCpSpecMs(double t,int Nspec,double *cpi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return cpi : array with species specific heats at constant pressure [J/(kg.K)]
*/

  int ans, i ;
  double tdim = t ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TCDND_getCpSpecMs", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) tdim *= TC_Tref_ ;

  ans = TC_getCpSpecMs(tdim,Nspec,cpi)  ;

  if ( TC_nonDim_ == 1 )
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i]  /= TC_cpref_ ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species specific heat at constant pressure based on temperature.
*/
int TC_getCpSpecMs(double t,int Nspec,double *cpi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return cpi : array with species specific heats at constant pressure [J/(kg.K)]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getCpSpecMs", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getCpSpecMsTab(t,cpi) ;
  else
    ans = TC_getCpSpecMsFcn(t,cpi) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species specific heat at constant volume based on temperature.
*/
int TC_getCvSpecMs(double t,int Nspec,double *cvi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return cpi : array with species specific heats at constant volume [J/(kg.K)]
*/

  int i, ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getCvSpecMs", Nspec, TC_Nspec_ ) ;

  ans = TC_getCvSpecMl(t,Nspec,cvi) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cvi[i] /= TC_sMass_[i] ;

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes species heat capacities at constant pressure based on temperature.
         Input temperature is normalized, output heat capacities are also normalized before exit.
*/
int TCDND_getCpSpecMl(double t,int Nspec,double *cpi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return cpi : array with species heat capacities at constant pressure [J/(kmol.K)]
*/

  int ans, i ;
  double tdim = t ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TCDND_getCpSpecMl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) tdim *= TC_Tref_ ;

  ans = TC_getCpSpecMl(tdim,Nspec,cpi)  ;

  if ( TC_nonDim_ == 1 )
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i]  /= ( TC_cpref_ * TC_Wref_ ) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species heat capacities at constant pressure based on temperature.
*/
int TC_getCpSpecMl(double t,int Nspec,double *cpi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return cpi : array with species heat capacities at constant pressure [J/(kmol.K)]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getCpSpecMl", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getCpSpecMlTab(t,cpi) ;
  else
    ans = TC_getCpSpecMlFcn(t,cpi) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species heat capacities at constant volume based on temperature.
*/
int TC_getCvSpecMl(double t,int Nspec,double *cvi)
{
/**
  \param t     : temperature T [K]
  \param Nspec : no. of species
  \return cvi  : array with species heat capacities at constant pressure [J/(kmol.K)]
*/

  int i, ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getCpSpecMl", Nspec, TC_Nspec_ ) ;

  ans = TC_getCpSpecMl(t,Nspec,cvi) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cvi[i] -= TC_Runiv_ ;

  return ( ans ) ;

}

/*
                  _   __  __     ____  _   _           _
        __ _  ___| |_|  \/  |___|___ \| | | |_ __ ___ (_)_  __
       / _` |/ _ \ __| |\/| / __| __) | |_| | '_ ` _ \| \ \/ /
      | (_| |  __/ |_| |  | \__ \/ __/|  _  | | | | | | |>  <
       \__, |\___|\__|_|  |_|___/_____|_| |_|_| |_| |_|_/_/\_\
       |___/

*/
/**
 * \ingroup thermo
 * \brief Computes mixture specific enthalpy based on temperature and species mass fractions.
         Input temperature is normalized, output enthalpy is normalized before exit.
*/
int TCDND_getMs2HmixMs(double *scal,int Nvars,double *hmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K], mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return hmix : pointer to mixture specific enthalpy [J/kg]
*/
  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getMs2HmixMs", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMs2HmixMs(scal,Nvars,hmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *hmix   /= ( TC_cpref_ * TC_Tref_ ) ;
    scal[0] /= TC_Tref_  ;
  }

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture specific enthalpy based on temperature and species mass fractions.
*/
int TC_getMs2HmixMs(double *scal,int Nvars,double *hmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,
                Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K],
                mass fractions Y []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return hmix : pointer to mixture specific enthalpy [J/kg]
*/

  int ans, i ;
  double temperature, *Yspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMs2HmixMs", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Yspec       = &scal[1] ;

  /* compute species enthalpies */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;

  *hmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *hmix += Yspec[i]*TC_hks[i] ;

  return ( ans ) ;

}
/*
                   _   __  __ _ ____  _   _           _
         __ _  ___| |_|  \/  | |___ \| | | |_ __ ___ (_)_  __
        / _` |/ _ \ __| |\/| | | __) | |_| | '_ ` _ \| \ \/ /
       | (_| |  __/ |_| |  | | |/ __/|  _  | | | | | | |>  <
        \__, |\___|\__|_|  |_|_|_____|_| |_|_| |_| |_|_/_/\_\
        |___/

*/
/**
 * \ingroup thermo
 * \brief Computes mixture molar enthalpy based on temperature and species mole fractions.
         Input temperature is normalized, output enthalpy is normalized before exit.
 */
int TCDND_getMl2HmixMl(double *scal,int Nvars,double *hmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return hmix : pointer to mixture molar enthalpy [J/kmol]
*/
  int ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getMl2HmixMl", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMl2HmixMl(scal,Nvars,hmix) ;

  if ( TC_nonDim_ == 1 )
  {
    *hmix   /= ( TC_cpref_ * TC_Wref_ * TC_Tref_ ) ;
    scal[0] /= TC_Tref_  ;
  }

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes mixture molar enthalpy based on temperature and species mole fractions.
*/
int TC_getMl2HmixMl(double *scal,int Nvars,double *hmix)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,X<sub>1</sub>,X<sub>2</sub>,...,X<sub>Nspec</sub>), temperature T [K], mole fractions X []
  \param Nvars : no. of variables = N<sub>spec</sub>+1
  \return hmix : pointer to mixture molar enthalpy [J/kmol]
*/

  int ans, i ;
  double temperature, *Xspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMl2HmixMl", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Xspec       = &scal[1] ;

  /* compute species enthalpies */
  ans = TC_getHspecMl(temperature,TC_Nspec_,TC_hks) ;

  *hmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *hmix += Xspec[i]*TC_hks[i] ;

  return ( ans ) ;

}

/*
                   _   _   _
         __ _  ___| |_| | | |___ _ __   ___  ___
        / _` |/ _ \ __| |_| / __| '_ \ / _ \/ __|
       | (_| |  __/ |_|  _  \__ \ |_) |  __/ (__
        \__, |\___|\__|_| |_|___/ .__/ \___|\___|
        |___/                   |_|

*/
/**
 * \ingroup thermo
 * \brief Computes species specific enthalpies based on temperature.
         Input temperature is normalized, output enthalpies are also normalized before exit.
*/
int TCDND_getHspecMs(double t,int Nspec,double *hi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species specific enthalpies [J/kg]
*/

  int ans, i ;
  double tdim = t ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TCDND_getHspecMs", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) tdim *= TC_Tref_ ;

  ans = TC_getHspecMs(tdim,Nspec,hi)  ;

  if ( TC_nonDim_ == 1 )
    for ( i = 0 ; i < TC_Nspec_ ; i++ )
      hi[i]  /= ( TC_cpref_ * TC_Tref_ ) ;

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes species specific enthalpies based on temperature.
*/
int TC_getHspecMs(double t,int Nspec,double *hi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species specific enthalpies [J/kg]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getHspecMs", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getHspecMsTab(t,hi) ;
  else
    ans = TC_getHspecMsFcn(t,hi) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species specific internal energies based on temperature.
*/
int TC_getUspecMs(double t,int Nspec,double *ui)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species specific internal energies [J/kg]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getUspecMs", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getUspecMsTab(t,ui) ;
  else
    ans = TC_getUspecMsFcn(t,ui) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species molar enthalpies based on temperature.
 *        Input temperature is normalized, output enthalpies are also normalized before exit.
*/
int TCDND_getHspecMl(double t,int Nspec,double *hi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species molar enthalpies [J/kmol]
*/

  int ans, i ;
  double tdim = t ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TCDND_getHspecMl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) tdim *= TC_Tref_ ;

  ans = TC_getHspecMl(tdim,Nspec,hi)  ;

  if ( TC_nonDim_ == 1 )
    for ( i = 0 ; i < TC_Nspec_ ; i++ )
      hi[i]  /= ( TC_cpref_ * TC_Wref_ * TC_Tref_ ) ;

  return ( ans ) ;

}
/**
 * \ingroup thermo
 * \brief Computes species molar enthalpies based on temperature.
*/
int TC_getHspecMl(double t,int Nspec,double *hi)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species molar enthalpies [J/kmol]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getHspecMl", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getHspecMlTab(t,hi) ;
  else
    ans = TC_getHspecMlFcn(t,hi) ;

  return ( ans ) ;

}

/**
 * \ingroup thermo
 * \brief Computes species molar internal energies based on temperature.
*/
int TC_getUspecMl(double t,int Nspec,double *ui)
{
/**
  \param t : temperature T [K]
  \param Nspec : no. of species
  \return hi : array with species molar internal energies [J/kmol]
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getUspecMl", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getUspecMlTab(t,ui) ;
  else
    ans = TC_getUspecMlFcn(t,ui) ;

  return ( ans ) ;

}
/*------------------------------------------------------------------------
                         _                  _            _
      ___  ___ _ __ ___ (_)      _ __  _ __(_)_   ____ _| |_ ___
     / __|/ _ \ '_ ` _ \| |_____| '_ \| '__| \ \ / / _` | __/ _ \
     \__ \  __/ | | | | | |_____| |_) | |  | |\ V / (_| | ||  __/
     |___/\___|_| |_| |_|_|     | .__/|_|  |_| \_/ \__,_|\__\___|
                                |_|

------------------------------------------------------------------------ */
/*
              _    ____      ____
    __ _  ___| |_ / ___|_ __/ ___| _ __   ___  ___
   / _` |/ _ \ __| |   | '_ \___ \| '_ \ / _ \/ __|
  | (_| |  __/ |_| |___| |_) |__) | |_) |  __/ (__
   \__, |\___|\__|\____| .__/____/| .__/ \___|\___|
   |___/               |_|        |_|

     _____             _______     _
    |  ___|__ _ __    / /_   _|_ _| |__
    | |_ / __| '_ \  / /  | |/ _` | '_ \
    |  _| (__| | | |/ /   | | (_| | |_) |
    |_|  \___|_| |_/_/    |_|\__,_|_.__/


*/
int TC_getCpSpecMsFcn(double t,double *cpi)
{

  int ans, i ;

  ans = TC_getCpSpecMlFcn(t,cpi) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i] /= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getCpSpecMs1Fcn(double t, int i, double *cpi)
{

  int ans ;

  ans = TC_getCpSpecMl1Fcn(t,i,cpi) ;

  *cpi /= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getCpSpecMlFcn(double t,double *cpi)
{
  int i, ipol, i9t ;

  i9t = 0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    int icpst ;

    if ( i9t < TC_nNASA9coef_ )
    {
      if ( i == TC_spec9t_[i9t] )
      {
	TC_getCpFcn9t(t, i9t, &cpi[i]) ; i9t++ ;
        continue ;
      }
    }

    ipol = 0 ;
    if (t > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;
    cpi[i] = TC_cppol_[icpst]+t*(TC_cppol_[icpst+1]+
			         t*(TC_cppol_[icpst+2]+
				    t*(TC_cppol_[icpst+3]+
				       t*TC_cppol_[icpst+4]
				       )
				    )
				 ) ;
  }

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i] *= TC_Runiv_ ;

  return ( 0 ) ;

}

int TC_getEntr0SpecMlFcn(double t,double *s0i)
{
  int i, ipol ;

  double one2 = 1.0/2.0 ;
  double one3 = 1.0/3.0 ;
  double one4 = 1.0/4.0 ;
  double tln = log(t);

  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    int icpst ;
    ipol = 0 ;
    if (t > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;
    s0i[i] = TC_cppol_[icpst]*tln+
           t*(TC_cppol_[icpst+1]+
              t*(one2*TC_cppol_[icpst+2]+
				         t*(one3*TC_cppol_[icpst+3]+
				            t*one4*TC_cppol_[icpst+4]
		 		           )
				        )
				     ) + TC_cppol_[icpst+6];
  }
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) s0i[i] *= TC_Runiv_ ;

  return ( 0 ) ;

}

int TC_getEntr0SpecMsFcn(double t,double *s0i)
{
  int i, ans ;

  ans = TC_getEntr0SpecMlFcn(t,s0i) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    s0i[i] /= TC_sMass_[i] ;

  return ( 0 ) ;

}

int TC_getMl2EntrMixMl(double *scal, int Nvars, double *smix)
{
  int i, ans ;

  double temperature = scal[0];
  double *Xspec      = &scal[1];

  double lprat = log(TC_pressure_/ATMPA) ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMl2EntrMixMl", Nvars, TC_Nvars_ ) ;

  ans = TC_getEntr0SpecMlFcn(temperature,TC_entr) ;
  *smix = 0.0;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    (*smix) += (TC_entr[i] - TC_Runiv_ * (log(Xspec[i]) + lprat)) * Xspec[i];

  return ( 0 ) ;

}

int TC_getMs2EntrMixMs(double *scal, int Nvars, double *smix)
{
  int i, ans ;
  double Wmix, temperature, *Yspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getMs2EntrMixMs", Nvars, TC_Nvars_ ) ;

  temperature = scal[0];
  Yspec       = &scal[1];

  double lprat = log(TC_pressure_/ATMPA) ;

  ans = TC_getMs2Ml(Yspec,Nvars-1,TC_y2x2y) ;
  ans = TC_getEntr0SpecMlFcn(temperature,TC_entr) ;
  *smix = 0.0;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    (*smix) += (TC_entr[i] - TC_Runiv_ * (log(TC_y2x2y[i]) + lprat)) * TC_y2x2y[i];

  Wmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    Wmix += TC_y2x2y[i]*TC_sMass_[i] ;

  *smix /= Wmix ;

  return ( 0 ) ;

}

int TC_getCpSpecMl1Fcn(double t, int i, double *cpi)
{
  int ipol, icpst, i9t;

  if ( TC_nNASA9coef_>0 )
  {
    for ( i9t=0; i9t<TC_nNASA9coef_; i9t++ )
      if ( i == TC_spec9t_[i9t] )
      {
        TC_getCpFcn9t(t, i9t, cpi) ;
        return ( 0 ) ;
      }
  }

  ipol = 0 ;
  if (t > TC_Tmi_[i]) ipol = 1 ;
  icpst = i*7*2+ipol*7 ;
  *cpi = TC_cppol_[icpst]+
         t*(TC_cppol_[icpst+1]+
			      t*(TC_cppol_[icpst+2]+
			         t*(TC_cppol_[icpst+3]+
				          t*TC_cppol_[icpst+4]
				         )
				      )
			    ) ;

  *cpi *= TC_Runiv_ ;

  return ( 0 ) ;

}

int TC_getCpSpecMsTab(double t1,double *cpi)
{

  int i ;
  double trat = (t1-TMIN)*TC_odelT_ ;
  int iTab = (int) trat ;

  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nspec_ ;

  for ( i = 0 ; i < TC_Nspec_ ; i++,iTab++ )
    cpi[i] = TC_cptab[iTab]+trat*(TC_cptab[iTab+TC_Nspec_]-TC_cptab[iTab]) ;

  return ( 0 ) ;

}

int TC_getCpSpecMlTab(double t1,double *cpi)
{

  int ans, i ;

  ans = TC_getCpSpecMsFcn(t1,cpi) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i] *= TC_sMass_[i] ;

  return ( ans ) ;

}
/*
                 _    ____       __  __ _      __  __     ____
       __ _  ___| |_ / ___|_ __ |  \/  (_)_  _|  \/  |___|  _ \
      / _` |/ _ \ __| |   | '_ \| |\/| | \ \/ / |\/| / __| |_) |
     | (_| |  __/ |_| |___| |_) | |  | | |>  <| |  | \__ \  __/
      \__, |\___|\__|\____| .__/|_|  |_|_/_/\_\_|  |_|___/_|
      |___/               |_|

*/
int TC_getCpMixMsP(double *scal,int Nvars,double *cpmix)
{

  int ans, i ;
  double temperature, *Yspec ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getCpMixMsP", Nvars, TC_Nvars_ ) ;

  temperature =  scal[0] ;
  Yspec       = &scal[1] ;

  /* compute species cp's */
  ans = TC_getCpSpecMsP(temperature,TC_Nspec_,TC_cpks) ;

  *cpmix = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) *cpmix += Yspec[i]*TC_cpks[i] ;

  return ( ans ) ;

}

/*
                 _    ____      ____                  __  __     ____
       __ _  ___| |_ / ___|_ __/ ___| _ __   ___  ___|  \/  |___|  _ \
      / _` |/ _ \ __| |   | '_ \___ \| '_ \ / _ \/ __| |\/| / __| |_) |
     | (_| |  __/ |_| |___| |_) |__) | |_) |  __/ (__| |  | \__ \  __/
      \__, |\___|\__|\____| .__/____/| .__/ \___|\___|_|  |_|___/_|
      |___/               |_|        |_|

*/
int TC_getCpSpecMsP(double t,int Nspec,double *cpi)
{

  int ans ;
  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getCpSpecMsP", Nspec, TC_Nspec_ ) ;

  if ( TC_tab_ )
    ans = TC_getCpSpecMsPtab(t, cpi) ;
  else
    ans = TC_getCpSpecMsPFcn(t, cpi) ;

  return ( ans ) ;

}

int TC_getCpSpecMsPFcn(double t,double *cpi)
{

  int i, ipol, i9t ;

  i9t = 0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    int icpst ;

    if ( i9t < TC_nNASA9coef_ )
    {
      if ( i == TC_spec9t_[i9t] )
      {
	TC_getCpFcnP9t(t, i9t, &cpi[i]) ; i9t++ ;
        continue ;
      }
    }

    ipol = 0 ;
    if (t > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;
    cpi[i] = TC_cppol_[icpst+1]+t*(2.0*TC_cppol_[icpst+2]+
			           t*(3.0*TC_cppol_[icpst+3]+
				      t*4.0*TC_cppol_[icpst+4]
				      )
				   ) ;
  }
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) cpi[i] *= TC_Runiv_/TC_sMass_[i] ;

  return ( 0 ) ;

}

int TC_getCpSpecMs1PFcn(double t, int i, double *cpi)
{

  int ipol, icpst, i9t ;

  if ( TC_nNASA9coef_>0 )
  {
    for ( i9t=0; i9t<TC_nNASA9coef_; i9t++ )
      if ( i == TC_spec9t_[i9t] )
      {
        TC_getCpFcnP9t(t, i9t, cpi) ;
        return ( 0 ) ;
      }
  }

  ipol = 0 ;
  if (t > TC_Tmi_[i]) ipol = 1 ;
  icpst = i*7*2+ipol*7 ;
  *cpi = TC_cppol_[icpst+1]+t*(2.0*TC_cppol_[icpst+2]+
			           t*(3.0*TC_cppol_[icpst+3]+
				      t*4.0*TC_cppol_[icpst+4]
				      )
				   ) ;
  *cpi *= TC_Runiv_/TC_sMass_[i] ;

  return ( 0 ) ;

}

int TC_getCpSpecMsPtab(double t1,double *cpi)
{

  int i ;
  double trat = (t1-TMIN)*TC_odelT_ ;
  int iTab = (int) trat ;

  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nspec_ ;

  for ( i = 0 ; i < TC_Nspec_ ; i++,iTab++ )
    cpi[i] = TC_cpPtab[iTab]+trat*(TC_cpPtab[iTab+TC_Nspec_]-TC_cpPtab[iTab]) ;

  return ( 0 ) ;

}
/*
                _   _   _                          _____             _______     _
      __ _  ___| |_| | | |___ _ __   ___  ___     |  ___|__ _ __    / /_   _|_ _| |__
     / _` |/ _ \ __| |_| / __| '_ \ / _ \/ __|    | |_ / __| '_ \  / /  | |/ _` | '_ \
    | (_| |  __/ |_|  _  \__ \ |_) |  __/ (__     |  _| (__| | | |/ /   | | (_| | |_) |
     \__, |\___|\__|_| |_|___/ .__/ \___|\___|    |_|  \___|_| |_/_/    |_|\__,_|_.__/
     |___/                   |_|

*/
int TC_getHspecMsFcn(double t,double *hi)
{

  int ans, i ;

  ans = TC_getHspecMlFcn(t,hi) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ )
    hi[i] /= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getUspecMsFcn(double t,double *ui)
{

  int ans, i ;

  ans = TC_getUspecMlFcn(t,ui) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) ui[i] /= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getHspecMlFcn(double t,double *hi)
{

  int i, ipol, i9t ;
  double one2, one3, one4, one5 ;

  one2 = 1.0/2.0 ;
  one3 = 1.0/3.0 ;
  one4 = 1.0/4.0 ;
  one5 = 1.0/5.0 ;
  i9t = 0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    int icpst ;

    if ( i9t < TC_nNASA9coef_ )
    {
      if ( i == TC_spec9t_[i9t] )
      {
	TC_getHspecFcn9t(t, i9t, &hi[i]) ; i9t++ ;
        continue ;
      }
    }

    ipol = 0 ;
    if (t > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;
    hi[i] = t*(TC_cppol_[icpst]+
	       t*(one2*TC_cppol_[icpst+1]+
		  t*(one3*TC_cppol_[icpst+2]+
		     t*(one4*TC_cppol_[icpst+3]+
			t*(one5*TC_cppol_[icpst+4]
			   )
			)
		     )
		  )
	       ) + TC_cppol_[icpst+5];
  }

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) hi[i] *= TC_Runiv_ ;

  return ( 0 ) ;

}

int TC_getUspecMlFcn(double t,double *ui)
{

  int ans, i ;
  ans = TC_getHspecMlFcn(t, ui) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) ui[i] -= TC_Runiv_*t ;

  return ( ans ) ;

}

int TC_getHspecMsTab(double t1,double *hi)
{

  int i ;
  double trat = (t1-TMIN)*TC_odelT_ ;
  int iTab = (int) trat ;

  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nspec_ ;

  for ( i = 0 ; i < TC_Nspec_ ; i++,iTab++ )
    hi[i] = TC_htab[iTab]+trat*(TC_htab[iTab+TC_Nspec_]-TC_htab[iTab]) ;

  return ( 0 ) ;

}

int TC_getUspecMsTab(double t1,double *ui)
{

  int ans, i ;
  ans = TC_getUspecMlTab(t1, ui);
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) ui[i] /= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getHspecMlTab(double t1,double *hi)
{

  int ans, i ;

  ans = TC_getHspecMsTab(t1,hi) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) hi[i] *= TC_sMass_[i] ;

  return ( ans ) ;

}

int TC_getUspecMlTab(double t1,double *ui)
{

  int ans, i ;
  ans = TC_getHspecMlTab(t1, ui);
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) ui[i] -= TC_Runiv_*t1 ; ;

  return ( ans ) ;

}

/*
 */
int TC_getCpFcn9t(double t, int icnt, double *cpi)
{
  int ist, irng ;

  ist = icnt*NTH9RNGMAX;
  if ( ( t <  TC_spec9trng_[ist*2] ) ||
       ( t >  TC_spec9trng_[ist*2+2*TC_spec9nrng_[icnt]-1]) )
  {
    printf("Error: temperature outside the range for species %d\n",TC_spec9t_[icnt]) ;
    fflush(stdout);
    return (-1) ;
  }

  /* determine range */
  irng = 0;
  while ( t>TC_spec9trng_[ist*2+2*irng+1] ) irng++ ;

  ist = (ist+irng)*9 ;
  (*cpi) = (TC_spec9coefs_[ist]/t+TC_spec9coefs_[ist+1])/t ;
  (*cpi)+= TC_spec9coefs_[ist+2]+t*(TC_spec9coefs_[ist+3]+
                                    t*(TC_spec9coefs_[ist+4]+
                                       t*(TC_spec9coefs_[ist+5]+
                                          t*TC_spec9coefs_[ist+6]
                                         )
                                      )
                                  ) ;

  return ( 0 ) ;

}

int TC_getCpFcnP9t(double t, int icnt, double *cpi)
{
  int ist, irng ;

  ist = icnt*NTH9RNGMAX;
  if ( ( t <  TC_spec9trng_[ist*2] ) ||
       ( t >  TC_spec9trng_[ist*2+2*TC_spec9nrng_[icnt]-1]) )
  {
    printf("Error: temperature outside the range for species %d\n",TC_spec9t_[icnt]) ;
    fflush(stdout);
    return (-1) ;
  }

  /* determine range */
  irng = 0;
  while ( t>TC_spec9trng_[ist*2+2*irng+1] ) irng++ ;

  ist = (ist+irng)*9 ;
  (*cpi) = (-2.0*TC_spec9coefs_[ist]/t-TC_spec9coefs_[ist+1])/(t*t) ;
  (*cpi)+= TC_spec9coefs_[ist+3]+t*(TC_spec9coefs_[ist+4]/2.0 +
                                    t*(TC_spec9coefs_[ist+5]/3.0 +
                                       t*TC_spec9coefs_[ist+6]/4.0
                                      )
				   ) ;

  return ( 0 ) ;

}

int TC_getHspecFcn9t(double t, int icnt, double *hi)
{
  int ist, irng ;

  ist = icnt*NTH9RNGMAX;
  if ( ( t <  TC_spec9trng_[ist*2] ) ||
       ( t >  TC_spec9trng_[ist*2+2*TC_spec9nrng_[icnt]-1]) )
  {
    printf("Error: temperature outside the range for species %d\n",TC_spec9t_[icnt]) ;
    fflush(stdout);
    return (-1) ;
  }

  /* determine range */
  irng = 0;
  while ( t>TC_spec9trng_[ist*2+2*irng+1] ) irng++ ;

  ist = (ist+irng)*9 ;
  (*hi) = -TC_spec9coefs_[ist]/t+TC_spec9coefs_[ist+1]*log(t) ;
  (*hi)+= t*(TC_spec9coefs_[ist+2] +
             t*(TC_spec9coefs_[ist+3]/2.0 +
                t*(TC_spec9coefs_[ist+4]/3.0 +
                   t*(TC_spec9coefs_[ist+5]/4.0 +
                      t*TC_spec9coefs_[ist+6]/5.0
		     )
		  )
	       )
	    ) ;
  (*hi)+= TC_spec9coefs_[ist+7] ;

  return ( 0 ) ;

}
