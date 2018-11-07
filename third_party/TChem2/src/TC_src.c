#include "TC_defs.h"
#include "TC_interface.h"

#include <math.h>
#include <stdlib.h>


/*! \file TC_src.c
 *    \brief Source term and Jacobian functions
 */ 

/*
      _____              _                          
     /  ___|            | |                         
     \ `--. _ __ ___    | |_ ___ _ __ _ __ ___  ___ 
      `--. \ '__/ __|   | __/ _ \ '__| '_ ` _ \/ __|
     /\__/ / | | (__ _  | ||  __/ |  | | | | | \__ \
     \____/|_|  \___(_)  \__\___|_|  |_| |_| |_|___/
                                                          
*/
/*
                          _   ____           
                __ _  ___| |_/ ___| _ __ ___ 
               / _` |/ _ \ __\___ \| '__/ __|
              | (_| |  __/ |_ ___) | | | (__ 
               \__, |\___|\__|____/|_|  \___|
               |___/                         

*/
/**
 * \ingroup srcs
 * \brief Returns dimensional/non-dimensional source term for 
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TCDND_getSrc(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return omega : array of N<sub>spec</sub>+1 source terms (possibly normalized) for temperature and 
                   species mass fractions equations: 
		   omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TCDND_getSrc", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getSrc(scal, Nvars, omega) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    omega[0] *= (TC_timref_/TC_Tref_) ;
    for ( i = 1 ; i < TC_Nspec_+1 ; i++ ) omega[i] *= TC_timref_ ;
  }

  return ( ans ) ;
 
}
/**
 * \ingroup srcs
 * \brief Returns source term for 
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSrc(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return omega : array of N<sub>spec</sub>+1 source terms for temperature and species
                   mass fractions: omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int i, ans ;
  double temperature, *Yspec, *omegaspec, rhomix, orho, cpmix ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TC_getSrc", Nvars, TC_Nvars_ ) ;

  ans         = 0 ;
  temperature =  scal[0] ;
  Yspec       = &scal[1] ;
  omegaspec   = &omega[1] ;

  /* get species molar reaction rates */
  ans = TC_getTY2RRms(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    printf("TC_getSrc() : Yspec[%-3d] = %e, omega[%-3d] = %e\n",i+1,Yspec[i],i+1,omega[i+1]) ;
#endif

  /* get density, cpmix */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  printf("TC_getSrc() : (rhomix,cpmix)=(%e,%e)\n",rhomix,cpmix) ;
#endif

  /* get species enthalpies */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  
  /* transform reaction rate to source term (*Wi/rho) */
  orho = 1.0/rhomix ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= orho ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) omegaspec[i] = 0.0 ;
#endif
  
  /* compute source term for temperature */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] -= omegaspec[i]*TC_hks[i] ;
  omega[0] /= cpmix ;

#ifdef DEBUGMSG
  printf("TC_getSrc() : omega[(%-3d] = %e\n",0,omega[0]) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) printf("TC_getSrc() : omega[(%-3d] = %e\n",i+1,omega[i+1]) ;
  exit(1) ;
#endif

  return ( ans ) ;
 
}

/**
 * \ingroup srcs
 * \brief Returns \f$S\f$ term for 
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSmat(double *scal,int Nvars, double *Smat)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return Smat : array of \f$(N_{spec}+1)\times 2N_{reac}\f$ holding
                  the S components in column major format
*/

  int i, j, kspec, indx, indxR, ans ;
  double temperature, rhomix, orho, orhocp, cpmix ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TC_getSrc", Nvars, TC_Nvars_ ) ;

  ans         = 0 ;
  temperature = scal[0] ;

  /* clean Smat */
  for ( i = 0 ; i < (TC_Nspec_+1)*2*TC_Nreac_ ; i++ ) Smat[i] = 0.0 ;

  /* get density, cpmix */
  ans = TC_getRhoMixMs  (scal,Nvars,&rhomix) ; 
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ; 
  if ( ans != 0 ) return ( ans ) ;
  orho   = 1.0/rhomix ;
  orhocp = orho/cpmix ;

  /* get species enthalpies and multiply by the weights */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i<TC_Nspec_; i++) TC_hks[i] *= TC_sMass_[i] ;
  
  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         -= TC_reacNukiDbl_[indx]*TC_hks[kspec] ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         -= TC_reacNukiDbl_[indx]*TC_hks[kspec] ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

  } /* done loop over reactions */
  
  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;

    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      int i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;

      for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)]         -= TC_reacRealNuki_[indxR]*TC_hks[kspec] ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)        ] -= TC_reacRealNuki_[indxR]*TC_hks[kspec] ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  for ( i = 0 ; i< TC_Nreac_; i++)
    Smat[i*(TC_Nspec_+1)] *= orhocp;

  for ( i = 0 ; i< TC_Nreac_*(TC_Nspec_+1); i++)
    Smat[i+TC_Nreac_*(TC_Nspec_+1)] = -Smat[i];

  return ( ans ) ;
 
}

/**
 * \ingroup srcs
 * \brief Returns source term for constant volume ignition system
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSrcCV(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return omega : array of N<sub>spec</sub>+1 source terms for temperature and species
                   mass fractions: omega[0] : [(K/s)], omega[1...N<sub>spec</sub>] : [(1/s)]
*/

  int i, ans ;
  double temperature, *Yspec, *omegaspec, rhomix, orho, cpmix, wmix, sumYoW, gamma ;
  double psum ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TC_getSrcCV", Nvars, TC_Nvars_ ) ;

  if ( TC_rhoset_ == 0 ) {
    printf("TC_getSrcCV() - density was not set -> Abort !\n") ;
    exit(1);
  }

  /* re-scale pressure to get appropriate density */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  TC_setThermoPres(TC_pressure_ * TC_rho_/rhomix) ;
  rhomix = TC_rho_; // bugfix based on MV observation 2015/01/19

  // Checking equivalence in formulations
  /* wmix = 0.0 ; */
  /* for ( i = 0 ; i < TC_Nspec_ ; i++ ) wmix += scal[i+1]/TC_sMass_[i] ; */
  /* TC_pressure_ = TC_rho_*TC_Runiv_*scal[0]*wmix; */
  /* TC_prescgs_  = TC_pressure_*10.0 ; */
  /* rhomix = TC_rho_; */

  ans         = 0 ;
  temperature =  scal[0]  ;
  Yspec       = &scal[1]  ;
  omegaspec   = &omega[1] ;

  /* get species molar reaction rates */
  ans = TC_getTY2RRms(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef DEBUGMSG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    printf("TC_getSrc() : Yspec[%-3d] = %e, omega[%-3d] = %e\n",i+1,Yspec[i],i+1,omega[i+1]) ;
#endif

  /* get cpmix, gamma, and gamma factors */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i] / TC_sMass_[i] ;
  wmix = 1.0 / sumYoW ;

  /* gamma */
  gamma = cpmix / ( cpmix - TC_Runiv_/wmix );

#ifdef DEBUGMSG
  printf("TC_getSrc() : (cpmix,gamma)=(%e,%e)\n",cpmix,gamma) ;
#endif

  /* get species enthalpies */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  
  /* transform reaction rate to source term (*Wi/rho) */
  orho = 1.0/rhomix ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= orho ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) omegaspec[i] = 0.0 ;
#endif
  
  /* compute source term for temperature */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] += omegaspec[i]*TC_hks[i] ;
  omega[0] = gamma * ( -omega[0] / cpmix ) ;

  /* add second term */
  psum = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) psum += omegaspec[i]/TC_sMass_[i] ;
  omega[0] += (gamma-1) * temperature * wmix * psum ;
  

#ifdef DEBUGMSG
  printf("TC_getSrc() : omega[(%-3d] = %e\n",0,omega[0]) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) printf("TC_getSrcCV() : omega[(%-3d] = %e\n",i+1,omega[i+1]) ;
  exit(1) ;
#endif

  return ( ans ) ;
 
}

/**
 * \ingroup srcs
 * \brief Returns \f$S\f$ matrix for constant volume ignition system
 * \f[\frac{\partial T}{\partial t}=\omega_0,\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 *  based on temperature T and species mass fractions Y's
 */

int TC_getSmatCV(double *scal,int Nvars,double *Smat)
{

/**
   \param scal : array of Nspec+1 doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$ N_{vars}= N_{spec}+1\f$ 
   \return Smat : array of \f$(N_{spec}+1)\times 2N_{reac}\f$ holding
                  the S components in column major format
*/

  int i, j, kspec, indx, indxR, ans ;
  double temperature, *Yspec, rhomix, orho, cpmix, wmix, sumYoW, gamma, c1g, c2g ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TC_getSmatCV", Nvars, TC_Nvars_ ) ;

  if ( TC_rhoset_ == 0 ) {
    printf("TC_getSmatCV() - density was not set -> Abort !\n") ;
    exit(1);
  }

  /* clean Smat */
  for ( i = 0 ; i < (TC_Nspec_+1)*2*TC_Nreac_ ; i++ ) Smat[i] = 0.0 ;

  /* re-scale pressure to get appropriate density */
  ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;
  TC_setThermoPres(TC_pressure_ * TC_rho_/rhomix) ;

  ans         = 0 ;
  temperature =  scal[0]  ;
  Yspec       = &scal[1]  ;

  /* get cpmix, gamma, and gamma factors */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix ) ;
  if ( ans != 0 ) return ( ans ) ;
  sumYoW = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumYoW += Yspec[i] / TC_sMass_[i] ;
  wmix = 1.0 / sumYoW ;

  orho   = 1.0/TC_rho_ ;

  /* gamma */
  gamma = cpmix / ( cpmix - TC_Runiv_/wmix );
  c1g = gamma*orho/cpmix ;
  c2g = (gamma-1.0)*temperature*wmix*orho ;

  /* get species enthalpies and multiply by the weights */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i<TC_Nspec_; i++) TC_hks[i] *= TC_sMass_[i] ;
  
  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         += TC_reacNukiDbl_[indx]*(-TC_hks[kspec]*c1g+c2g) ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      Smat[i*(TC_Nspec_+1)]         += TC_reacNukiDbl_[indx]*(-TC_hks[kspec]*c1g+c2g) ;
      Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacNukiDbl_[indx]*TC_sMass_[kspec]*orho ;
    }

  } /* done loop over reactions */
  
  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;

    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      int i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;

      for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)]         += TC_reacRealNuki_[indxR]*(-TC_hks[kspec]*c1g+c2g) ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
        Smat[i*(TC_Nspec_+1)        ] += TC_reacRealNuki_[indxR]*(-TC_hks[kspec]*c1g+c2g) ;
        Smat[i*(TC_Nspec_+1)+1+kspec] += TC_reacRealNuki_[indxR]*TC_sMass_[kspec]*orho ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  for ( i = 0 ; i< TC_Nreac_*(TC_Nspec_+1); i++)
    Smat[i+TC_Nreac_*(TC_Nspec_+1)] = -Smat[i];

  return ( ans ) ;
 
}

/*
                _   ____            ____                
      __ _  ___| |_/ ___| _ __ ___ / ___|___  _ __  ___ 
     / _` |/ _ \ __\___ \| '__/ __| |   / _ \| '_ \/ __|
    | (_| |  __/ |_ ___) | | | (__| |__| (_) | | | \__ \
     \__, |\___|\__|____/|_|  \___|\____\___/|_| |_|___/
     |___/                                              

*/
/**
 * \ingroup srcs
 * \brief Returns source term (dimensional/non-dimensional) for 
 * \f[\frac{\partial \rho}{\partial t}=\omega_0,\rho\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 * based on \f$\rho\f$ and Y's
 */
int TCDND_getSrcCons(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(\rho,Y_1,Y_2,...,Y_N)\f$
	        density \f$[kg/m^3]\f$, mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub>+1 source terms for density and species mf
                   conservative formulation: 
		   omega[0] : \f$[kg/(m^3\cdot s)]\f$, 
		   omega[1...N<sub>spec</sub>] : \f$[kg/(m^3\cdot s)]\f$
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TCDND_getSrcCons", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getSrcCons(scal, Nvars, omega) ;

  if ( TC_nonDim_ == 1 ) 
  {

    scal[0] /= TC_Tref_ ;
    omega[0] *= (TC_timref_/TC_rhoref_) ;
    for ( i = 1 ; i < TC_Nspec_+1 ; i++ ) omega[i] *= TC_timref_ ;

  }

  return ( ans ) ;

}
/**
 * \ingroup srcs
 * \brief Returns source term for 
 * \f[\frac{\partial \rho}{\partial t}=\omega_0,\rho\frac{\partial Y_i}{\partial t}=\omega_i,\f]
 * based on \f$\rho\f$ and Y's
 */
int TC_getSrcCons(double *scal,int Nvars,double *omega)
{

/**
   \param scal : array of Nspec+1 doubles \f$(\rho,Y_1,Y_2,...,Y_N)\f$
	        density \f$[kg/m^3]\f$, mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub>+1 source terms for density and species mf
                   conservative formulation: 
		   omega[0] : \f$[kg/(m^3\cdot s)]\f$, 
		   omega[1...N<sub>spec</sub>] : \f$[kg/(m^3\cdot s)]\f$
*/

  int ans, i ;
  double temperature, rhomix, *Yspec, *omegaspec, cpmix, Wmix, sumOmega ;

  if ( Nvars != TC_Nvars_ ) 
    TC_errorMSG( 20, "TC_getSrcCons", Nvars, TC_Nvars_ ) ;

  ans       = 0 ;
  rhomix    =  scal[0] ;
  Yspec     = &scal[1] ;
  omegaspec = &omega[1] ;

  /* get temperature, cpmix, and mixture molecular weight */
  ans = TC_getTmixMs(scal,Nvars,&temperature) ;
  if ( ans != 0 ) return ( ans ) ;
  scal[0] = temperature ;

  /* get species molar reaction rates [kmol/(m3.s)]*/
  ans = TC_getTY2RRml(scal,Nvars,omegaspec) ;
  if ( ans != 0 ) return ( ans ) ;

  /* get cpmix [J/(kg.K)] and mixture molecular weight [kg/kmol] */
  ans = TC_getMs2CpMixMs(scal,Nvars,&cpmix) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix(Yspec,TC_Nspec_,&Wmix) ;
  if ( ans != 0 ) return ( ans ) ;

  /* get species enthalpies [J/kg] */
  ans = TC_getHspecMs(temperature,TC_Nspec_,TC_hks) ;
  if ( ans != 0 ) return ( ans ) ;

#ifdef NONNEG
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) 
    if ( ( Yspec[i]<0 ) && ( omegaspec[i] < 0 ) ) omegaspec[i] = 0.0 ;
#endif

  /* sum molar reaction rates (for density source term) */
  sumOmega = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) sumOmega += omegaspec[i] ;
  
  /* transform molar reaction rate to mass reaction rate(*Wi) [kg/(m3.s)]*/
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omegaspec[i] *= TC_sMass_[i] ;

  /* compute source term for density */
  omega[0] = 0.0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[0] += omegaspec[i]*TC_hks[i] ;
  omega[0] /= (cpmix*temperature) ;
  omega[0] -= Wmix*sumOmega ;

  return ( ans ) ;

}
/*
           _                 _     _                 
          | | __ _  ___ ___ | |__ (_) __ _ _ __  ___ 
       _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \/ __|
      | |_| | (_| | (_| (_) | |_) | | (_| | | | \__ \
       \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|___/

 */
/*
              _      _           _______   ___   _           _             _ 
    __ _  ___| |_   | | __ _  __|_   _\ \ / / \ | |_ __ ___ / | __ _ _ __ | |
   / _` |/ _ \ __|  | |/ _` |/ __|| |  \ V /|  \| | '_ ` _ \| |/ _` | '_ \| |
  | (_| |  __/ || |_| | (_| | (__ | |   | | | |\  | | | | | | | (_| | | | | |
   \__, |\___|\__\___/ \__,_|\___||_|   |_| |_| \_|_| |_| |_|_|\__,_|_| |_|_|
   |___/                                                                     

*/
/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian (dimensional/non-dimensional) for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on temperature T and species mass fractions Y's
 */
int TCDND_getJacTYNm1anl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNm1anl( scal, Nspec, jac) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<Nspec*Nspec; i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element*/
    for ( i=Nspec ; i < Nspec*Nspec ; i+=Nspec ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done with TCDND_getJacTYNm1anl */

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on T and Y's
 */
int TC_getJacTYNm1anl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec-1 ; k++ )
  {
    double drhodYk = - rhomix * Wmix * ( 1.0 / TC_sMass_[k] - 1.0 / TC_sMass_[Nspec-1] ) ;
    double dYndYk  = - 1.0 ;
    for ( i = 0,iJacR=Nspec*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
      jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk 
                 + TC_jacFull[iJacF+3+Nspec-1]*dYndYk ;
  }

  return ( ans ) ;

} /* done with TC_getJacTYNm1anl */

/*
               _      _           _______   ___   _             _ 
     __ _  ___| |_   | | __ _  __|_   _\ \ / / \ | | __ _ _ __ | |
    / _` |/ _ \ __|  | |/ _` |/ __|| |  \ V /|  \| |/ _` | '_ \| |
   | (_| |  __/ || |_| | (_| | (__ | |   | | | |\  | (_| | | | | |
    \__, |\___|\__\___/ \__,_|\___||_|   |_| |_| \_|\__,_|_| |_|_|
    |___/                                                         
*/
/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian (dimensional/non-dimensional) for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TCDND_getJacTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNanl( scal, Nspec, jac) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<(Nspec+1)*(Nspec+1); i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec+1 ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec+1 ; i < (Nspec+1)*(Nspec+1) ; i+=(Nspec+1) ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done with TCDND_getJacTYNanl */

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1anl", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */

  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec+1 ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec; k++ )
  {
    double drhodYk = - rhomix * Wmix / TC_sMass_[k]  ;
    for ( i = 0,iJacR=(Nspec+1)*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec+1 ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk ;
  }

  return ( ans ) ;

} /* done with TC_getJacTYNanl */
/*
                     _      _           _______   ___   _           _ 
           __ _  ___| |_   | | __ _  __|_   _\ \ / / \ | |_ __ ___ / |
          / _` |/ _ \ __|  | |/ _` |/ __|| |  \ V /|  \| | '_ ` _ \| |
         | (_| |  __/ || |_| | (_| | (__ | |   | | | |\  | | | | | | |
          \__, |\___|\__\___/ \__,_|\___||_|   |_| |_| \_|_| |_| |_|_|
          |___/                                                       

*/
/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian (dimensional/non-dimensional) 
 * for the system \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on temperature T and species mass fractions Y's
*/
int TCDND_getJacTYNm1(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYNm1", Nspec, TC_Nspec_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYNm1( scal, Nspec, jac, useJacAnl) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<Nspec*Nspec; i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec ; i < Nspec*Nspec ; i+=Nspec ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done TCDND_getJacTYNm1 */

/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_{N-1})\f$ 
 * based on T and Y's
*/
int TC_getJacTYNm1(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYNm1", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0,iJacR=0,iJacF = TC_Nvjac_*2 ; 
	i < Nspec ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k<Nspec-1 ; k++)
  {
    double drhodYk = - rhomix * Wmix * ( 1.0 / TC_sMass_[k] - 1.0 / TC_sMass_[Nspec-1] ) ;
    double dYndYk  = - 1.0 ;
    for ( i = 0,iJacR=Nspec*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec ; 
	  i++,iJacR++,iJacF+=TC_Nvjac_ )
      jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk 
      + TC_jacFull[iJacF+3+Nspec-1]*dYndYk ;
  }
 
  return ( ans ) ;

} /* done TC_getJacTYNm1 */

/*
                 _      _           _______   ___   _ 
       __ _  ___| |_   | | __ _  __|_   _\ \ / / \ | |
      / _` |/ _ \ __|  | |/ _` |/ __|| |  \ V /|  \| |
     | (_| |  __/ || |_| | (_| | (__ | |   | | | |\  |
      \__, |\___|\__\___/ \__,_|\___||_|   |_| |_| \_|
      |___/                                           
*/
/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
 */
int TCDND_getJacTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TCDND_getJacTYN", Nspec, TC_Nspec_ ) ;
	     
  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getJacTYN( scal, Nspec, jac, useJacAnl) ;

  if ( TC_nonDim_ == 1 ) 
  {
    scal[0] /= TC_Tref_ ;
    for ( i=0; i<(Nspec+1)*(Nspec+1); i++ ) jac[i] *= TC_timref_ ;
    /* first column, except first element */
    for ( i=1 ; i < Nspec+1 ; i++ ) jac[i] *= TC_Tref_ ;
    /* first line, except first element */
    for ( i=Nspec+1 ; i < (Nspec+1)*(Nspec+1) ; i+=(Nspec+1) ) jac[i] /= TC_Tref_ ;
  }

  return ( ans ) ;

} /* done TCDND_getJacTYN */
/**
 * \ingroup jacs
 * \brief Computes (analytical or numerical) Jacobian for the system 
 * \f$(T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
 */
int TC_getJacTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans, i, k, iJacR,iJacF ;
  double rhomix, Wmix, drhodT, *Yspec ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacTYN", Nspec, TC_Nspec_ ) ;

  ans   = 0 ;
  Yspec = &scal[1] ;

  /* first, get full Jacobian */
  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, TC_jacFull) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, TC_jacFull) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reduced Jacobian */
  ans = TC_getRhoMixMs(  scal, TC_Nvars_, &rhomix ) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getMs2Wmix ( Yspec, TC_Nspec_, &Wmix   ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* first column */
  drhodT = -rhomix/scal[0] ;
  for ( i = 0, iJacR=0, iJacF = TC_Nvjac_*2 ; 
	i < Nspec+1 ; 
	i++,iJacR++,iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+2]+TC_jacFull[iJacF]*drhodT ;

  /* rest of the columns */
  for ( k = 0 ; k < Nspec; k++ )
  {
    double drhodYk = - rhomix * Wmix / TC_sMass_[k] ;
    for ( i = 0,iJacR=(Nspec+1)*(k+1),iJacF = TC_Nvjac_*2 ; i < Nspec+1 ; 
	  i++, iJacR++, iJacF+=TC_Nvjac_ )
    jac[iJacR] = TC_jacFull[iJacF+3+k] + TC_jacFull[iJacF]*drhodYk ;
  }

  /* /\* skip rho and pressure lines and columns *\/   */
  /* for (iJacR = 0, iJacF=2*TC_Nvjac_+2, k = 0;       */
  /*      k<Nspec+1;                                   */
  /*      k++, iJacF +=2 )                             */
  /*   for ( i = 0; i<Nspec+1; i++, iJacR++, iJacF++ ) */
  /*     jac[ iJacR ] = TC_jacFull[ iJacF ] ;          */
	     
  return ( ans ) ;

} /* done TC_getJacTYN */

/*
                     _      _            ____  ____ _______   ___   _ 
           __ _  ___| |_   | | __ _  ___|  _ \|  _ \_   _\ \ / / \ | |
          / _` |/ _ \ __|  | |/ _` |/ __| |_) | |_) || |  \ V /|  \| |
         | (_| |  __/ || |_| | (_| | (__|  _ <|  __/ | |   | | | |\  |
          \__, |\___|\__\___/ \__,_|\___|_| \_\_|    |_|   |_| |_| \_|
          |___/                                                       

*/
/**
 * \ingroup jacs
 * \brief Computes (analytical) Jacobian for the system 
 * \f$(\rho,P,T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacRPTYN(double *scal, int Nspec, double *jac, unsigned int useJacAnl)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \param useJacAnl : flag for Jacobian type (1-analytical,other values-numerical) 
 *  \return jac  : Jacobian array 
 */

  int ans ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacRPTYN", Nspec, TC_Nspec_ ) ;

  if ( useJacAnl )
    ans = TC_getJacRPTYNanl(scal, Nspec, jac) ;
  else
    ans = TC_getJacRPTYNnum(scal, Nspec, jac) ;

  return ( ans ) ;

}

/**
 * \ingroup jacs
 * \brief Computes analytical Jacobian for the system 
 * \f$(\rho,P,T,Y_1,Y_2,\ldots,Y_N)\f$ 
 * based on T and Y's
*/
int TC_getJacRPTYNanl(double *scal, int Nspec, double *jac)
{
/**
 *   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_{N_{spec}})\f$
 *                 temperature T [K], mass fractions Y [] 
 *  \param Nspec : no. of species \f$N_{spec}\f$ 
 *  \return jac  : Jacobian array 
 */

  int ans, i, j, k, kspec, ir, indx, indxR, iJac, itbdy, ipfal, iord ;
  double *Yspec, t1, t_1, tln ;
  double rhomix,cpmix,cpmixP ;
  double sum1, sum2, sum3 ;

  if ( Nspec != TC_Nspec_ ) 
    TC_errorMSG( 10, "TC_getJacRPTYNanl", Nspec, TC_Nspec_ ) ;

  ans = 0 ;

  /* scal array contains T and all species mass fractions */
  /* clear data */
  for ( i = 0 ; i<TC_Nspec_; i++) TC_omg[i] = 0.0 ;

  /* Compute molar concentrations (kmol/m3) and transform them to moles/cm3 */
  ans = TC_getMs2Cc(scal,TC_Nspec_+1,TC_Xconc) ;
  if ( ans != 0 ) return ( ans ) ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_Xconc[i] *= 1.e-3 ;

  /* set species array */
  Yspec = &scal[1] ;

  /* compute functions and powers of temperature */
  t1  = scal[0] ;
  t_1 = 1.0/t1  ;
  tln = log(t1) ;

  /* get 3rd-body concentrations */
  ans = TC_get3rdBdyConc(TC_Xconc,TC_Mconc) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgk (t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getgkp(t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute forward and reverse rate constants */
  ans = TC_getkForRev (t1,t_1,tln) ;
  if ( ans != 0 ) return ( ans ) ;
  ans = TC_getkForRevP(t1,t_1) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute rate-of-progress */
  ans = TC_getRateofProg(TC_Xconc) ;
  if ( ans != 0 ) return ( ans ) ;

  /* compute pressure dependent factors */
  ans = TC_getCrnd( t1, t_1, tln, TC_Xconc, TC_Mconc ) ;
  if ( ans != 0 ) return ( ans ) ;

  /* assemble reaction rates */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {

    TC_rop[i] *= TC_Crnd[i] ;

    indx  = i*TC_maxSpecInReac_ ;
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

    indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      int kspec = TC_reacSidx_[indx] ;
      TC_omg[kspec] += ((double) TC_reacNuki_[indx])*TC_rop[i] ;
      indx++ ;
    }

  } /* done loop over the number of reactions */

  /* check for reactions with real stoichiometric coefficients */
  if ( TC_nRealNuReac_ > 0 )
    for ( ir = 0 ; ir < TC_nRealNuReac_; ir++ ) {

      i = TC_reacRnu_[ir] ;

      indx  = i *TC_maxSpecInReac_ ;
      indxR = ir*TC_maxSpecInReac_ ;      
      for ( j = 0 ; j<TC_reacNreac_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++) {
        int kspec = TC_reacSidx_[indx] ;
        TC_omg[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
	      indx++  ;
	      indxR++ ;
      }
    } /* done loop over the number of reactions */

  /* transform from mole/(cm3.s) to kg/(m3.s) */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_omg[i] *= (1.e3*TC_sMass_[i]) ;

  /* assemble jacobian terms, F[3+i,2...] */
  for ( i = 0 ; i < TC_Nvjac_* TC_Nvjac_ ; i++) jac[i] = 0.0 ;
  itbdy = 0 ; 
  ipfal = 0 ; 
  iord  = 0 ;
  for (j = 0 ; j < TC_Nreac_; j++) {

    unsigned int arbord ;
    /* compute {\partial Crnd}/{\partial T,Y1,Y2,...,YN} */
    TC_getCrndDer( j, &itbdy, &ipfal, t1, t_1, tln, TC_Xconc, TC_Mconc ) ;

    arbord = 0 ;
    if ( iord < TC_nOrdReac_ )
      if (TC_reacAOrd_[iord] == j) arbord = 1 ;

    if ( TC_reacScoef_[j] == -1 ) {

      /* reaction has integer stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
      
        /* skip if species i is not involved in reaction j */
        int indxIJ = j*TC_Nspec_+i ;
        if ( TC_NuIJ[indxIJ] == 0 ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_NuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;
      
        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {
          indx = j*TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx  = j*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            indx++ ;
          } /* done loop over species k in products (2nd term) */
        } /* done section for non-arbitrary order reactions */
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_NuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
	      } /* done section for arbitrary order reaction */
      } /* done loop over species (counter i) for integer stoichiometric coefficients */
    } /* end if integer coeffs */
    else {

      /* reaction has real stoichiometric coefficients */
      for ( i = 0 ; i < TC_Nspec_ ; i++) {
        int jr = TC_reacScoef_[j] ;
        /* skip if species i is not involved in reaction j */
        int indxIJ = jr*TC_Nspec_+i ;
        if ( TC_RealNuIJ[indxIJ] < TCSMALL ) continue ;
      
        /* add 1st term for the T,Y1,Y2,...,YN - derivatives */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        for ( k = -1 ; k < TC_Nspec_ ; k++,iJac++)
          jac[iJac] += TC_RealNuIJ[indxIJ]*TC_CrndDer[k+1]*(TC_ropFor[j]-TC_ropRev[j]) ;

        /* add 2nd term for T-derivative F[3+i,2] */
        iJac = TC_Nvjac_*(i+3)+2 ; 
        jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*(TC_ropFor[j]*TC_kforP[j]-TC_ropRev[j]*TC_krevP[j]) ;
      
        /* add 2nd term for species derivatives F[3+i,3+k] */
        if ( !arbord ) {

          indx = j *TC_maxSpecInReac_ ;
          for ( kspec = 0; kspec<TC_reacNreac_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            indx++ ;
          } /* done loop over species k in reactants (2nd term) */

          indx = j *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          for ( kspec = 0; kspec<TC_reacNprod_[j] ; kspec++) {
            k = TC_reacSidx_[indx] ;
            iJac = TC_Nvjac_*(i+3)+3+k ; 
            TC_getRateofProgDer(TC_Xconc, j, k, qfr);
            jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            indx++ ;	  
          } /* done loop over species k in products (2nd term) */
        }
        else {
          /* arbitrary order reaction */
          indx = iord*TC_maxOrdPar_ ;
          for ( kspec = 0; kspec<TC_maxOrdPar_ ; kspec++) {
            if (TC_specAOidx_[indx]<0) {
              k = -TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] += TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[0];
            }
            else if (TC_specAOidx_[indx]>0) {
              k = TC_specAOidx_[indx]-1 ;
              iJac = TC_Nvjac_*(i+3)+3+k ; 
              TC_getRateofProgDer(TC_Xconc, j, k, qfr);
              jac[iJac] -= TC_RealNuIJ[indxIJ]*TC_Crnd[j]*qfr[1];
            }
            indx++ ;
          }
        } /* done section for arbitrary order reaction */
      } /* done loop over species i */
    } /* done if real stoichiometric coeffs */

    if ( arbord ) iord++ ;

  } /* done loop over the number of reactions */
  
  /* get density, cpmix, species cp, species enthalpies */
  TC_getRhoMixMs  (scal,TC_Nvars_,&rhomix) ;
  TC_getCpMixMsP  (scal,TC_Nvars_,&cpmixP) ; 
  TC_getMs2CpMixMs(scal,TC_Nvars_,&cpmix ) ;
  TC_getHspecMs   (t1  ,TC_Nspec_,TC_hks ) ;

  /* Multiply by appropriate masses */
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {
    /* F[3+i,2] */
    iJac = TC_Nvjac_ * ( i + 3 ) + 2 ; 
    jac[iJac] *= TC_sMass_[i] / rhomix ;

    /* F[3+i,3+k] */
    iJac++ ;
    for ( k = 0 ; k < TC_Nspec_ ; k++, iJac++ )
      jac[iJac] *= TC_sMass_[i] ;
  }

  for ( k = 0 ; k < TC_Nspec_ ; k++ )
  {    
    double oMassk ;
    iJac = TC_Nvjac_*3 + 3 + k ;
    oMassk = 1.0 / TC_sMass_[k] ;
    /* F[3+i,3+k] */
    for ( i = 0 ; i < TC_Nspec_ ; i++, iJac += TC_Nvjac_ )
      jac[iJac] *= oMassk ;
  }

  /* transform T-derivatives to SI */
  for ( i=0, iJac=TC_Nvjac_*3+2 ; i < TC_Nspec_ ; i++, iJac+=TC_Nvjac_)
    jac[iJac] *= 1.e3 ;

  /* compute F[3+i,0] */
  iJac = TC_Nvjac_*3 ;
  for ( i=0, iJac=TC_Nvjac_*3 ; i < TC_Nspec_ ; i++, iJac+=TC_Nvjac_)
  {
    for ( k = 0 ; k<TC_Nspec_ ; k++ )
      jac[iJac] += Yspec[k]*jac[iJac+3+k] ;

    jac[iJac] = ( jac[iJac] - TC_omg[i] / rhomix ) / rhomix ;
  }

  /* compute F[2,0] */
  iJac = TC_Nvjac_*2 ;
  for ( k = 0 ; k < TC_Nspec_ ; k++ )
    jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;
  jac[iJac] /= cpmix ;


  /* compute F[2,2] */
  sum1 = 0.0; 
  sum2 = 0.0; 
  sum3 = 0.0 ;
  for ( k = 0, iJac=TC_Nvjac_*3+2 ; k < TC_Nspec_ ; k++,iJac+=TC_Nvjac_ )
  {
    sum1 += TC_hks [k]*TC_omg[k] ;
    sum2 += TC_cpks[k]*TC_omg[k] ;
    sum3 += TC_hks [k]*jac[iJac] ;
  }
  iJac = TC_Nvjac_*2+2 ;
  jac[iJac] = (sum1*cpmixP/cpmix-sum2)/(rhomix*cpmix)-sum3/cpmix ;

  /* compute F[2,3+i] */
  sum1 /= (rhomix*cpmix*cpmix) ;
  for ( i = 0,iJac=TC_Nvjac_*2+3 ; i < TC_Nspec_ ; i++,iJac++ )
  {
    for ( k = 0 ; k < Nspec ; k++)
      jac[iJac] -= TC_hks[k]*jac[iJac+TC_Nvjac_*(k+1)] ;

    jac[iJac] = jac[iJac]/cpmix+sum1*TC_cpks[i] ;

  }

  return ( ans ) ;

}
