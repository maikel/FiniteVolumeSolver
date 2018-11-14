#include "TC_defs.h"
#include "TC_interface.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


/*! \file TC_rr.c
 *  \brief Reaction rate functions
 */ 
/* #define DEBUGMSG */

/*
    ______                _   _              ______      _            
    | ___ \              | | (_)             | ___ \    | |           
    | |_/ /___  __ _  ___| |_ _  ___  _ __   | |_/ /__ _| |_ ___  ___ 
    |    // _ \/ _` |/ __| __| |/ _ \| '_ \  |    // _` | __/ _ \/ __|
    | |\ \  __/ (_| | (__| |_| | (_) | | | | | |\ \ (_| | ||  __/\__ \
    \_| \_\___|\__,_|\___|\__|_|\___/|_| |_| \_| \_\__,_|\__\___||___/
*/

/** 
   \brief Returns forward and reverse rate constants for all reactions
   at a specified temperature
*/
int TC_getkFR(double t, double *kF, double *kR)
{

  int i, ans ;

  double t_1, tln;
  t_1 = 1.0/t;
  tln = log(t) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgk(t,t_1,tln) ;

  /* compute forward and reverse rate constants */
  ans = TC_getkForRev(t,t_1,tln) ;

  for ( i = 0 ; i<TC_Nreac_; i++) kF[i] = TC_kfor[i] ;
  for ( i = 0 ; i<TC_Nreac_; i++) kR[i] = TC_krev[i] ;
  
  return (ans) ;

}

/**
 \brief Returns forward and reverse rate constants for all reactions
 at a specified temperature
 */
int TC_getkFRp(double t, double *kF, double *kR)
{
  
  int i, ans ;
  
  double t_1, tln;
  t_1 = 1.0/t;
  tln = log(t) ;
  
  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgkp(t,t_1,tln) ;
  
  /* compute forward and reverse rate constants */
  ans = TC_getkForRevPFcn(t_1) ;
  
  for ( i = 0 ; i<TC_Nreac_; i++) kF[i] = TC_kforP[i] ;
  for ( i = 0 ; i<TC_Nreac_; i++) kR[i] = TC_krevP[i] ;
  
  return (ans) ;
  
}

/**
 \brief Returns forward and reverse rate constants for all reactions
 at a specified temperature
 */
int TC_getkFRder(double t, double *kF, double *kR)
{
  
  int i, ans ;
  
  double t_1, tln;
  t_1 = 1.0/t;
  tln = log(t) ;
  
  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  ans = TC_getgk(t,t_1,tln) ;
  
  /* compute forward and reverse rate constants */
  ans = TC_getkForRevPder(t_1,tln) ;
  
  for ( i = 0 ; i<TC_Nreac_; i++) kF[i] = TC_kforPder[i] ;
  for ( i = 0 ; i<TC_Nreac_; i++) kR[i] = TC_krevPder[i] ;
  
  return (ans) ;
  
}

/**
   \brief Return current value of the Arrhenius parameters for forward rate constants.
   Return -1 if no data available, otherwise return 0 and store value in val
*/

int TC_getArhenFor( int ireac, int ipos, double *val )
{
/**
 * \param ireac : reaction index
 * \param ipos : index of Arrhenius parameter (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param *val : value of Arrhenius parameter
 */
  if ( ( ireac < 0 ) || ( ireac >= TC_Nreac_ ))
  {
    printf("TC_getArhenFor() : reaction index is wrong -> %d\n", ireac) ;
    return ( -1 ) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    printf("TC_getArhenFor() : index of Arrhenius parameter is wrong -> %d\n", ipos) ;
    return ( -1 ) ;
  }

  *val = TC_reacArhenFor_[ireac*3+ipos];

  return ( 0 ) ;

} /* Done with TC_getArhenFor */

/** 
   \brief Return current value of the Arrhenius parameters for the low/high pressure limit rate constants for Pressure-dependent reactions.
   Return -1 if no data available, otherwise return 0 and store value in val.
*/
int TC_getArhenPresDepFor( int ireac, int ipos, double *val )
{
/**
 * \param ireac : reaction index
 * \param ipos : index of Arrhenius parameter (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param *val : value of Arrhenius parameter
 */
  if ( ( ireac < 0 ) || ( ireac >= TC_nFallReac_ ))
  {
    printf("TC_getArhenresDepFor() : reaction index is wrong -> %d\n", ireac) ;
    return ( -1 ) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= TC_nFallPar_ ))
  {
    printf("TC_getArhenresDepFor() : index of Arrhenius parameter is wrong -> %d\n", ipos) ;
    return ( -1 ) ;
  }

  *val = TC_reacPpar_[ireac*TC_nFallPar_ + ipos];

  return ( 0 ) ;

} /* Done with TC_TC_getArhenPresDepFor */

/** 
   \brief Return current value of the Arrhenius parameters for reverse rate constants.
   Return -1 if no data available, otherwise return 0 and store value in val
*/
int TC_getArhenRev( int ireac, int ipos, double* val )
{
/**
 * \param ireac : reaction index
 * \param ipos : index of of Arrhenius parameter (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param *val : value of Arrhenius parameter
 */

  int i;
  int irevindx = -1 ;

  for ( i = 0 ; i < TC_nRevReac_; i++ )
    if ( TC_reacRev_[i] == ireac ) { irevindx = i; break ;}

  if ( irevindx < 0 )
  {
    return ( -1 ) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    return ( -1 ) ;
  }

  *val = TC_reacArhenRev_[irevindx*3+ipos];

  return ( 0 ) ;

} /* Done with TC_getArhenRev */

/*
                  _  _______   ______  ____  ____            _ 
        __ _  ___| ||_   _\ \ / /___ \|  _ \|  _ \ _ __ ___ | |
       / _` |/ _ \ __|| |  \ V /  __) | |_) | |_) | '_ ` _ \| |
      | (_| |  __/ |_ | |   | |  / __/|  _ <|  _ <| | | | | | |
       \__, |\___|\__||_|   |_| |_____|_| \_\_| \_\_| |_| |_|_|
       |___/                                                   

*/
/**
 * \brief Returns non-dimensional molar reaction rates, 
 * \f$\dot{\omega}_i\cdot t_{ref}\frac{W_{ref}}{\rho_{ref}}\f$, 
 * based on T and Y's
*/
int TCDND_getTY2RRml(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,Y_1,Y_2,...,Y_N)
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> molar reaction rates [kmol/(m3.s)] (but non-dimensional)
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTY2RRml", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getTY2RRml( scal, Nvars, omega ) ;

  if ( TC_nonDim_ == 1 ) 
  {

    /* non-dimensionalize mass reaction rates [kmol/(m3.s)] */
    double makeND = TC_timref_ * TC_Wref_ / TC_rhoref_ ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= makeND ;

    scal[0] /= TC_Tref_ ;

  }

  return ( ans ) ; 

}
/**
  \brief Returns molar reaction rates, \f$\dot{\omega}_i\f$, based on T and Y's
*/
int TC_getTY2RRml(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of \f$N_{spec}+1\f$ doubles \f$(T,Y_1,Y_2,...,Y_N)\f$:
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables \f$N_{vars}=N_{spec}+1\f$ 
   \return omega : array of \f$N_{spec}\f$ molar reaction rates \f$\dot{\omega}_i\f$ \f$\left[kmol/(m^3\cdot s)\right]\f$
*/
  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTY2RRml", Nvars, TC_Nvars_ ) ;

  /* compute molar concentrations */ 
  ans = TC_getMs2Cc(scal,Nvars,TC_Xconc) ;

  TC_scalIn[0] = scal[0] ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_scalIn[i+1] = TC_Xconc[i] ;

  /* get molar reaction rates */
  ans = TC_getReacRates(TC_scalIn, Nvars, omega) ;

  return ( ans ) ; 

}
/*
                _  _______   ______  ____  ____                
      __ _  ___| ||_   _\ \ / /___ \|  _ \|  _ \ _ __ ___  ___ 
     / _` |/ _ \ __|| |  \ V /  __) | |_) | |_) | '_ ` _ \/ __|
    | (_| |  __/ |_ | |   | |  / __/|  _ <|  _ <| | | | | \__ \
     \__, |\___|\__||_|   |_| |_____|_| \_\_| \_\_| |_| |_|___/
     |___/                                                     

*/
/**
  \brief Returns non-dimensional mass reaction rates based on T and Y's
*/
int TCDND_getTY2RRms(double *scal, int Nvars, double *omega)
{

/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,Y_1,Y_2,...,Y_N)
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> mass reaction rates  [kg/(m3.s)] (but non-dimensional)
*/
  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTY2RRms", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getTY2RRms( scal, Nvars, omega ) ;

  if ( TC_nonDim_ == 1 ) 
  {
    /* non-dimensionalize mass reaction rates [kg/(m3.s)] */
    double makeND = TC_timref_/TC_rhoref_ ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= makeND ;

    scal[0] /= TC_Tref_ ;

  }

  return ( ans ) ; 
}
/**
  \brief Returns mass reaction rates based on T and Y's
*/
int TC_getTY2RRms(double *scal, int Nvars, double *omega)
{
	
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,Y_1,Y_2,...,Y_N)
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> mass reaction rates  [kg/(m3.s)] 
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTY2RRms", Nvars, TC_Nvars_ ) ;

  /* compute molar concentrations */
  ans = TC_getMs2Cc(scal,Nvars,TC_Xconc) ;

  TC_scalIn[0] = scal[0] ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_scalIn[i+1] = TC_Xconc[i] ;

  /* get molar reaction rates */
  ans = TC_getReacRates(TC_scalIn, Nvars, omega) ;

  /* transform molar reaction rates to mass reaction rates */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= TC_sMass_[i] ;
  
  return ( ans ) ; 
}
/*
                  _  _______  ______ ____  ____  ____            _ 
        __ _  ___| ||_   _\ \/ / ___|___ \|  _ \|  _ \ _ __ ___ | |
       / _` |/ _ \ __|| |  \  / |     __) | |_) | |_) | '_ ` _ \| |
      | (_| |  __/ |_ | |  /  \ |___ / __/|  _ <|  _ <| | | | | | |
       \__, |\___|\__||_| /_/\_\____|_____|_| \_\_| \_\_| |_| |_|_|
       |___/                                                       

*/
/**
  \brief Returns non-dimensional molar reaction rates based on temperature T and molar concentrations XC
*/
int TCDND_getTXC2RRml(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,XC_1,XC_2,...,XC_N)
	        temperature T [K], molar concentrations XC [kmol/m3] (but non-dimensional)
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> (molar) reaction rates [kmol/(m3.s)]  (but non-dimensional)
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTXC2RRml", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) 
  {
    double makeND = TC_rhoref_ / TC_Wref_ ;
    scal[0] *= TC_Tref_ ;
    for ( i = 1 ; i < TC_Nvars_ ; i++ ) scal[i] *= makeND ;
  }

  ans = TC_getTXC2RRml( scal, Nvars, omega ) ;

  if ( TC_nonDim_ == 1 ) 
  {
    /* non-dimensionalize mass reaction rates [kmol/(m3.s)] */
    double makeND = TC_timref_ * TC_Wref_ / TC_rhoref_ ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= makeND ;

    scal[0] /= TC_Tref_ ;
    makeND = TC_Wref_ / TC_rhoref_ ;
    for ( i = 1 ; i < TC_Nvars_ ; i++ ) scal[i] *= makeND ;

  }

  return ( ans ) ; 

}
/**
  \brief Returns molar reaction rates based on temperature T and molar concentrations XC
*/
int TC_getTXC2RRml(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,XC_1,XC_2,...,XC_N)
	        temperature T [K], molar concentrations XC [kmol/m3] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> (molar) reaction rates [kmol/(m3.s)] 
*/	

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTXC2RRml", Nvars, TC_Nvars_ ) ;

  for ( i = 0 ; i < TC_Nspec_+1 ; i++ ) TC_scalIn[i] = scal[i] ;

  /* get molar reaction rates */
  ans = TC_getReacRates(TC_scalIn, Nvars, omega) ;

  return ( ans ) ; 
}
/*
                _  _______  ______ ____  ____  ____                
      __ _  ___| ||_   _\ \/ / ___|___ \|  _ \|  _ \ _ __ ___  ___ 
     / _` |/ _ \ __|| |  \  / |     __) | |_) | |_) | '_ ` _ \/ __|
    | (_| |  __/ |_ | |  /  \ |___ / __/|  _ <|  _ <| | | | | \__ \
     \__, |\___|\__||_| /_/\_\____|_____|_| \_\_| \_\_| |_| |_|___/
     |___/                                                         

*/
/**
  \brief Returns non-dimensional mass reaction rates based on T and molar concentrations
*/
int TCDND_getTXC2RRms(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of N<sub>spec</sub>+1 doubles (T,XC_1,XC_2,...,XC_N)
	        temperature T [K], molar concentrations XC [kmol/m3]  (but non-dimensional)
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> (mass) reaction rates [kg/(m3.s)]  (but non-dimensional)
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TCDND_getTXC2RRms", Nvars, TC_Nvars_ ) ;

  if ( TC_nonDim_ == 1 ) 
  {
    double makeND = TC_rhoref_ / TC_Wref_ ;
    scal[0] *= TC_Tref_ ;
    for ( i = 1 ; i < TC_Nvars_ ; i++ ) scal[i] *= makeND ;
  }

  ans = TC_getTXC2RRms( scal, Nvars, omega ) ;

  if ( TC_nonDim_ == 1 ) 
  {
    /* non-dimensionalize mass reaction rates [kg/(m3.s)] */
    double makeND = TC_timref_/TC_rhoref_ ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= makeND ;

    scal[0] /= TC_Tref_ ;
    makeND = TC_Wref_ / TC_rhoref_ ;
    for ( i = 1 ; i < TC_Nvars_ ; i++ ) scal[i] *= makeND ;

  }

  return ( ans ) ; 
}
/**
  \brief Returns mass reaction rates based on T and molar concentrations
*/
int TC_getTXC2RRms(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of Nspec+1 doubles (T,XC_1,XC_2,...,XC_N)
	        temperature T [K], molar concentrations XC [kmol/m3] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return omega : array of N<sub>spec</sub> (mass) reaction rates [kg/(m3.s)] 
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getTXC2RRms", Nvars, TC_Nvars_ ) ;

  for ( i = 0 ; i < TC_Nspec_+1 ; i++ ) TC_scalIn[i] = scal[i] ;

  /* get molar reaction rates */
  ans = TC_getReacRates(TC_scalIn, Nvars, omega) ;

  /* transform molar reaction rates to mass reaction rates */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= TC_sMass_[i] ;
  
  return ( ans ) ; 

}
/**
  \brief Returns rate-of-progress variables based on temperature T and species mass fractions Y's
*/
int TC_getRops(double *scal, int Nvars, double *datarop) 
{
/**
   \param scal : array of N<sub>reac</sub>+1 doubles (T,Y_1,Y_2,...,Y_N)
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return datarop : array of N<sub>reac</sub> rate-of-progress variables  [kmol/(m3.s)]
*/


  int i, ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getRops", Nvars, TC_Nvars_ ) ;

  /* compute molar concentrations */ 
  ans = TC_getMs2Cc(scal,Nvars,TC_Xconc) ;

  TC_scalIn[0] = scal[0] ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_scalIn[i+1] = TC_Xconc[i] ;

  TC_getRopsLocal( TC_scalIn ) ;  

  for ( i = 0 ; i < TC_Nreac_ ; i++ ) datarop[i] = TC_rop[i] ;
  
  return ( ans ) ;

}

/**
  \brief Returns forward and reverse rate-of-progress variables based on T and Y's
*/
int TC_getRfrb(double *scal, int Nvars, double *dataRfrb) 
{
/**
   \param scal : array of Nspec+1 doubles (T,Y_1,Y_2,...,Y_N)
	        temperature T [K], mass fractions Y [] 
   \param Nvars : no. of variables = N<sub>spec</sub>+1 
   \return datarop : array of N<sub>reac</sub> forward rate-of-progress variables  
                     and N<sub>reac</sub> reverse rate-of-progress variables [kmol/(m3.s)]
*/

  int i, ans ;

  if ( Nvars != TC_Nvars_ )
    TC_errorMSG( 20, "TC_getRfrb", Nvars, TC_Nvars_ ) ;

  /* compute molar concentrations */ 
  ans = TC_getMs2Cc(scal,Nvars,TC_Xconc) ;

  TC_scalIn[0] = scal[0] ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) TC_scalIn[i+1] = TC_Xconc[i] ;

  TC_getRopsLocal( TC_scalIn ) ;  

  for ( i = 0 ; i < TC_Nreac_ ; i++ ) dataRfrb[i]            = TC_ropFor[i] ;
  for ( i = 0 ; i < TC_Nreac_ ; i++ ) dataRfrb[i+TC_Nreac_ ] = TC_ropRev[i] ;

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
int TC_getRopsLocal(double *scal)
{
  /*
    Input :
         scal  : array of Nspec+1 doubles (T,XC_1,XC_2,...,XC_N)
	         temperature T [K], molar concentrations XC [kmol/m3]
    Output -> stored in TC_rop [kmol/(m3.s)] 
  */ 

  /* clear data */
  int i ;
  double t1, t_1, tln ;

  /* Transform molar concentrations (kmol/m3) to (moles/cm3) */
  for ( i = 0 ; i<TC_Nspec_; i++) TC_Xconc[i] = scal[i+1]*1.e-3 ;

  /* compute functions and powers of temperature */
  t1  = scal[0] ;
  t_1 = 1.0/t1  ;
  tln = log(t1) ;

  /* get 3rd-body concentrations */
  TC_get3rdBdyConc(TC_Xconc, TC_Mconc) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  TC_getgk(t1,t_1,tln) ;

  /* compute forward and reverse rate constants */
  TC_getkForRev(t1,t_1,tln) ;

  /* compute rate-of-progress */
  TC_getRateofProg(TC_Xconc) ;

  /* compute pressure dependent factors */
  TC_getCrnd( t1, t_1, tln, TC_Xconc, TC_Mconc ) ;

  /* assemble final rate of progress */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    TC_rop   [i] *= TC_Crnd[i] ;
    TC_ropFor[i] *= TC_Crnd[i] ;
    TC_ropRev[i] *= TC_Crnd[i] ;

  }

  /* transform from mole/(cm3.s) to kmol/(m3.s) */
  for ( i = 0 ; i < TC_Nreac_ ; i++ ) TC_rop   [i] *= 1.e3 ;
  for ( i = 0 ; i < TC_Nreac_ ; i++ ) TC_ropFor[i] *= 1.e3 ;
  for ( i = 0 ; i < TC_Nreac_ ; i++ ) TC_ropRev[i] *= 1.e3 ;

  return ( 0 ) ;

}
/*
                 _   ____                 ____       _            
       __ _  ___| |_|  _ \ ___  __ _  ___|  _ \ __ _| |_ ___  ___ 
      / _` |/ _ \ __| |_) / _ \/ _` |/ __| |_) / _` | __/ _ \/ __|
     | (_| |  __/ |_|  _ <  __/ (_| | (__|  _ < (_| | ||  __/\__ \
      \__, |\___|\__|_| \_\___|\__,_|\___|_| \_\__,_|\__\___||___/
      |___/                                                       
                                                                                       
*/
/**
  \brief Returns molar reaction rates, \f$\dot{\omega}_i\f$, based on T and molar 
         concentractions XC's (semi-private function)
*/
int TC_getReacRates(double *scal, int Nvars, double *omega)
{
/**
   \param scal : array of \f$N_{spec}+1\f$ doubles \f$((T,XC_1,XC_2,...,XC_N)\f$:
	        temperature T [K], molar concentrations XC \f$[kmol/m^3]\f$ 
   \param Nvars : no. of variables \f$N_{vars}=N_{spec}+1\f$ 
   \return omega : array of \f$N_{spec}\f$ molar reaction rates \f$\dot{\omega}_i\f$ \f$\left[kmol/(m^3\cdot s)\right]\f$
*/

  /* clear data */
  int i, j ;
  double t1, t_1, tln ;
  int indx, indxR, kspec ;

  for ( i = 0 ; i<TC_Nspec_; i++) omega[i] = 0.0 ;

  /* Transform molar concentrations (kmol/m3) to (moles/cm3) */
  for ( i = 0 ; i<TC_Nspec_; i++)  TC_Xconc[i] = scal[i+1]*1.e-3 ;

  /* compute functions and powers of temperature */
  t1  = scal[0] ;
  t_1 = 1.0/t1  ;
  tln = log(t1) ;

  /* get 3rd-body concentrations */
  TC_get3rdBdyConc(TC_Xconc,TC_Mconc) ;

  /* compute (-ln(T)+dS/R-dH/RT) for each species */
  TC_getgk(t1,t_1,tln) ;

  /* compute forward and reverse rate constants */
  TC_getkForRev(t1,t_1,tln) ;

  /* compute rate-of-progress */
  TC_getRateofProg(TC_Xconc) ;

  /* compute pressure dependent factors */
  TC_getCrnd( t1, t_1, tln, TC_Xconc, TC_Mconc ) ;

  /* assemble reaction rates */

  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    TC_rop[i] *= TC_Crnd[i] ;

    indx = i * TC_maxSpecInReac_;
    for ( j = 0; j<TC_reacNreac_[i] ; j++, indx++ )
    {
      kspec = TC_reacSidx_[indx] ;
      omega[kspec] += TC_reacNukiDbl_[indx] * TC_rop[i] ;
    }

    indx = i * TC_maxSpecInReac_ + TC_maxSpecInReac_/2;
    for ( j = 0; j<TC_reacNprod_[i] ; j++, indx++ )
    {
      kspec = TC_reacSidx_[indx] ;
      omega[kspec] += TC_reacNukiDbl_[indx] * TC_rop[i] ;
    }

  } /* done loop over reactions */

  /* check for reactions with real stoichiometric coefficients */
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
        omega[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
        kspec = TC_reacSidx_[indx] ;
        omega[kspec] += TC_reacRealNuki_[indxR]*TC_rop[i] ;
        indx++  ;
        indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  /* transform from mole/(cm3.s) to kmol/(m3.s) */
  for ( i = 0 ; i < TC_Nspec_ ; i++ ) omega[i] *= 1.e3 ;

  return ( 0 ) ;

}

/*
          _         _____           ____    ______            _   _    ______ _____ 
         | |_ __   |_   _|    _    / ___|  / /  _ \          | | | |  / /  _ \_   _|
  _____  | | '_ \    | |    _| |_  \___ \ / /| |_) |  _____  | |_| | / /| |_) || |  
 |_____| | | | | |   | |   |_   _|  ___) / / |  _ <  |_____| |  _  |/ / |  _ < | |  
         |_|_| |_|   |_|     |_|   |____/_/  |_| \_\         |_| |_/_/  |_| \_\|_|  

                                                                       
*/

double TC_getgi(double t1,int ispec)
{
  
  int ipol, icpst ;
  double t2, t3, t4, gi ;
  t2  = t1*t1 ;
  t3  = t2*t1 ;
  t4  = t3*t1 ;

  t2 = t2 /  6.0 ;
  t3 = t3 / 12.0 ;
  t4 = t4 / 20.0 ;

  ipol = 0 ; if (t1 > TC_Tmi_[ispec]) ipol = 1 ;
  icpst = ispec*7*2+ipol*7 ;

  gi = TC_cppol_[icpst+6]-TC_cppol_[icpst]
     +(TC_cppol_[icpst]-1.0)*log(t1)
     +TC_cppol_[icpst+1]*t1*0.5
     +TC_cppol_[icpst+2]*t2
     +TC_cppol_[icpst+3]*t3
     +TC_cppol_[icpst+4]*t4
     -TC_cppol_[icpst+5]/t1; 

  return ( gi ) ;

}

int TC_getgk(double t1,double t_1, double tln)
{
  
  if ( TC_tab_ )
    TC_getgkTab(t1) ;
  else
    TC_getgkFcn(t1,t_1,tln) ;

  /* done computing (gk=-lnT+S/R-H/RT) */
  
  return ( 0 ) ;

}

int TC_getgkFcn(double t1,double t_1, double tln)
{
  
  int i, i9t ;
  double t2, t3, t4 ;
  t2  = t1*t1 ;
  t3  = t2*t1 ;
  t4  = t3*t1 ;

  t2 = t2 /  6.0 ;
  t3 = t3 / 12.0 ;
  t4 = t4 / 20.0 ;

  i9t = 0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {

    int ipol, icpst ;

    if ( i9t < TC_nNASA9coef_ )
    {
      if ( i == TC_spec9t_[i9t] )
      {
	TC_getgkFcn9t(t1, i9t, &TC_gk[i]) ; i9t++ ;
        continue ;
      }
    }

    ipol = 0 ; if (t1 > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;

    TC_gk[i] = TC_cppol_[icpst+6]-TC_cppol_[icpst]
             +(TC_cppol_[icpst]-1.0)*tln
              +TC_cppol_[icpst+1]*t1*0.5
              +TC_cppol_[icpst+2]*t2
              +TC_cppol_[icpst+3]*t3
              +TC_cppol_[icpst+4]*t4
              -TC_cppol_[icpst+5]*t_1
     ;

  } 

  /* done computing (gk=-lnT+S/R-H/RT) */
  
  return ( 0 ) ;

}

int TC_getgkFcn9t(double t, int icnt, double *gki) 
{
  int ist, irng ;
  double pfac[4] = {0.5,1.0/6.0,1.0/12.0,1.0/20.0} ;
  
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
  (*gki) = TC_spec9coefs_[ist+8]-TC_spec9coefs_[ist+7]/t-log(t) ;
  (*gki) += TC_spec9coefs_[ist]*0.5/(t*t);
  (*gki) -= TC_spec9coefs_[ist+1]*(log(t)+1.0)/t;
  (*gki) += TC_spec9coefs_[ist+2]*(log(t)-1.0);
  (*gki) += t*(TC_spec9coefs_[ist+3]*pfac[0] + 
	       t*(TC_spec9coefs_[ist+4]*pfac[1] + 
                  t*(TC_spec9coefs_[ist+5]*pfac[2] + 
		     t*TC_spec9coefs_[ist+6]*pfac[3]
                    )
		 )
	       ) ; 

  return ( 0 ) ;
  
}

int TC_getgkTab(double t1)
{

  int i ;
  double trat = (t1-TMIN)*TC_odelT_ ;

  int iTab = (int) trat ;
  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nspec_ ;

  for ( i = 0 ; i < TC_Nspec_ ; i++,iTab++ )
    TC_gk[i] = TC_gktab[iTab]+trat*(TC_gktab[iTab+TC_Nspec_]-TC_gktab[iTab]) ;

  /* done computing (gk=-lnT+S/R-H/RT) */
  return ( 0 ) ;

}

int TC_getgkp(double t1,double t_1, double tln)
{
  
  if ( TC_tab_ )
    TC_getgkpTab(t1) ;
  else
    TC_getgkpFcn(t1,t_1,tln) ;

  /* done computing gkp=d(gk)/dT */
  return ( 0 ) ;

}

int TC_getgkpFcn(double t1,double t_1, double tln)
{
  
  int i, i9t ;
  double t2, t3, one3 ;
  t2  = t1*t1 ;
  t3  = t2*t1 ;
  one3 = 1.0/3.0 ;

  t2 = t2 / 4.0 ;
  t3 = t3 / 5.0 ;

  i9t = 0 ;
  for ( i = 0 ; i < TC_Nspec_ ; i++ )
  {

    int ipol, icpst ;

    if ( i9t < TC_nNASA9coef_ )
    {
      if ( i == TC_spec9t_[i9t] )
      {
	TC_getgkpFcn9t(t1, i9t, &TC_gkp[i]) ; i9t++ ;
        continue ;
      }
    }
    
    ipol = 0 ; if (t1 > TC_Tmi_[i]) ipol = 1 ;
    icpst = i*7*2+ipol*7 ;

    TC_gkp[i] = t_1*(TC_cppol_[icpst]-1.0+TC_cppol_[icpst+5]*t_1)
                    +TC_cppol_[icpst+1]*0.5
                    +TC_cppol_[icpst+2]*t1*one3
                    +TC_cppol_[icpst+3]*t2
                    +TC_cppol_[icpst+4]*t3
               ;

  } 

  /* done computing gkp=d(gk)/dT */
  return ( 0 ) ;

}

int TC_getgkpFcn9t(double t, int icnt, double *gki) 
{
  int ist, irng ;
  double pfac[4] = {0.5,1.0/3.0,1.0/4.0,1.0/5.0} ;
  
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
  (*gki) = (TC_spec9coefs_[ist+7]/t-1.0)/t ;
  (*gki) -= TC_spec9coefs_[ist]/(t*t*t);
  (*gki) += TC_spec9coefs_[ist+1]*log(t)/(t*t);
  (*gki) += TC_spec9coefs_[ist+2]/t;
  (*gki) += TC_spec9coefs_[ist+3]*pfac[0] + 
             t*(TC_spec9coefs_[ist+4]*pfac[1] + 
               t*(TC_spec9coefs_[ist+5]*pfac[2] + 
                  t*TC_spec9coefs_[ist+6]*pfac[3]
                 )
		) ;
  return ( 0 ) ;
  
}

int TC_getgkpTab(double t1)
{

  int i ;

  double trat = (t1-TMIN)*TC_odelT_ ;

  int iTab = (int) trat ;
  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nspec_ ;

  for ( i = 0 ; i < TC_Nspec_ ; i++,iTab++ )
    TC_gkp[i] = TC_gkPtab[iTab]+trat*(TC_gkPtab[iTab+TC_Nspec_]-TC_gkPtab[iTab]) ;

  /* done computing gkp=d(gk)/dT */
  return ( 0 ) ;

}

/*
      ______        
      \  ___)       
       \ \     _  __           __ _ 
        > >   | |/ /          / _` |
       / /__  | / /  _    _  | (_| |  _    
      /_____) |__/  | | _(_)  \__, | | | __
                    | |/ / |  |___/  | |/ /
                    |   <| |         |   < 
                    |_|\_\_|         |_|\_\
                          
 */
double TC_getSumNuGk(int i, double *gkLoc)
{
  
  int  j, indx, kspec ;
  double sumNuGk ;

  sumNuGk = 0.0 ;
  indx = i*TC_maxSpecInReac_ ;
  for ( j = 0; j < TC_reacNreac_[i] ; j++, indx++)
  {
    kspec = TC_reacSidx_[indx] ;
    sumNuGk += TC_reacNukiDbl_[indx]*gkLoc[kspec] ;
    
  } /* done for loop over reactants */

  indx = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
  for ( j = 0; j< TC_reacNprod_[i] ; j++, indx++)
  {
    kspec = TC_reacSidx_[indx] ;
    sumNuGk += TC_reacNukiDbl_[indx]*gkLoc[kspec] ;
    
  } /* done loop over products */

  return ( sumNuGk ) ;
  
}

double TC_getSumRealNuGk(int i,int ir, double *gkLoc)
{
  
  int j, indx, indxR, kspec ;
  double sumNuGk ;

  sumNuGk = 0.0 ;
  indx  = i *TC_maxSpecInReac_ ;
  indxR = ir*TC_maxSpecInReac_ ;
  for ( j = 0; j < TC_reacNreac_[i] ; j++)
  {
    kspec = TC_reacSidx_[indx] ;
    sumNuGk += TC_reacRealNuki_[indxR]*gkLoc[kspec] ;
    indx++  ;
    indxR++ ;
  } /* done for loop over reactants */

  indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
  indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
  for ( j = 0; j < TC_reacNprod_[i] ; j++)
  {
    kspec = TC_reacSidx_[indx] ;
    sumNuGk += TC_reacRealNuki_[indxR]*gkLoc[kspec] ;
    indx++  ;
    indxR++ ;    
  } /* done loop over products */

  return ( sumNuGk ) ;
  
}

/*
     _____         _       _               _                              
    |___ / _ __ __| |     | |__   ___   __| |_   _    ___ ___  _ __   ___ 
      |_ \| '__/ _` |_____| '_ \ / _ \ / _` | | | |  / __/ _ \| '_ \ / __|
     ___) | | | (_| |_____| |_) | (_) | (_| | |_| | | (_| (_) | | | | (__ 
    |____/|_|  \__,_|     |_.__/ \___/ \__,_|\__, |  \___\___/|_| |_|\___|
                                             |___/                        

*/
int TC_get3rdBdyConc(double *concX, double *concM)
{

  int i ;
  for ( i = 0 ; i < TC_Nreac_ ; i++ ) concM[i] = 1.0 ;

  if ( TC_nThbReac_ > 0 )
  {

    int j, indx, kspec, ireac ;
    double concSum ;

    concSum = 0.0 ;
    for ( j = 0 ; j<TC_Nspec_ ; j++) concSum += concX[j] ;
    
    for ( i = 0 ; i<TC_nThbReac_ ; i++ )
    {

      ireac = TC_reacTbdy_[i] ;
      concM[ireac] = concSum ;

      for ( j = 0 ; j<TC_reacTbno_[i] ; j++ )
      {
        indx  = i*TC_maxTbInReac_+j ;
        kspec = TC_specTbdIdx_[indx] ;
        concM[ireac] += (TC_specTbdEff_[indx]-1.0)*concX[kspec] ;

      } /* done loop over efficiencies */

    } /* done loop over third-body reaction list */

  } /* done if there are third-body reactions */

  return ( 0 ) ;

}

/*
       _     __                ___      _                   
      | | __/ _| ___  _ __    ( _ )    | | ___ __ _____   __
      | |/ / |_ / _ \| '__|   / _ \/\  | |/ / '__/ _ \ \ / /
      |   <|  _| (_) | |     | (_>  <  |   <| | |  __/\ V / 
      |_|\_\_|  \___/|_|      \___/\/  |_|\_\_|  \___| \_/  

*/
int TC_getkForRev(double t1, double t_1, double tln)
{
  
  if ( TC_tab_ )
    TC_getkForRevTab( t1 ) ;
  else
    TC_getkForRevFcn( t_1, tln ) ;

#ifdef DEBUGMSG
  {
    int i ;
    FILE *fdbg = fopen("dbg_kfr.dat","w") ;
    for ( i = 0 ; i<TC_Nreac_; i++ )
    {
      fprintf(fdbg,"%d : %20.12e %20.12e\n",i,TC_kfor[i],TC_krev[i]) ;
    }
    fclose(fdbg);
    exit(1) ;
  }
#endif

  return ( 0 ) ;

}

int TC_getkForRevFcn(double t_1, double tln)
{

  int i, j, irev, iplog, indx, plogtest ;
  double logP, ki, ki1, *api, *api1 ;
  
  if ( TC_nPlogReac_ > 0 ) logP = log(TC_pressure_/ATMPA) ;
  
  irev  = 0 ;
  iplog = 0 ;
  for ( i = 0 ; i<TC_Nreac_; i++ )
  {

    plogtest = 0 ;
    if ( iplog < TC_nPlogReac_ )
      plogtest = (i == TC_reacPlogIdx_[iplog]);

    if (plogtest) {
      indx = 4*TC_reacPlogPno_[iplog];
      api  = &(TC_reacPlogPars_[indx]) ;
      
      if (logP <= (*api)) {
        TC_kfor[i] = exp(*(api+1) + *(api+2) * tln - *(api+3) * t_1 ) ;
      } /* Done with branch p<pmin */
      
      else {
        indx = 4*(TC_reacPlogPno_[iplog+1]-1);
        api  = &(TC_reacPlogPars_[indx]) ;

        if (logP >= (*api)) {
          TC_kfor[i] = exp( *(api +1) + *(api +2) * tln - *(api +3) * t_1 ) ;
        } /* Done with branch p>pmax */
        
        else {
          
          for (j=TC_reacPlogPno_[iplog]; j<TC_reacPlogPno_[iplog+1]-1;j++) {
            indx = 4*j;
            if (logP<TC_reacPlogPars_[indx+4]) {
              api  = &(TC_reacPlogPars_[indx  ]);
              api1 = &(TC_reacPlogPars_[indx+4]);
              ki  = *(api +1) + *(api +2) * tln - *(api +3) * t_1 ;
              ki1 = *(api1+1) + *(api1+2) * tln - *(api1+3) * t_1 ;
              TC_kfor[i] = exp( ki+(logP - (*api))/( (*api1) - (*api))*(ki1-ki) );
              break ;
            }
            
          } /* Done loop over all intervals */
          
        } /* Done with branch pmin<p<pmax */
        
      } /* Done with branch p>pmin */
      
      iplog++;
      
    } /* Done if reaction has a PLOG form */
    else {
      indx = i*3 ;
      TC_kfor[i] = TC_reacArhenFor_[indx]
                  *exp( TC_reacArhenFor_[indx+1] * tln - TC_reacArhenFor_[indx+2] * t_1 ) ;
    }
    TC_krev[i] = 0.0 ;

    /* is reaction reversible ? */
    if ( TC_isRev_[i] )
    {
      /* are reverse Arhenius parameters given ? */
      if ( irev < TC_nRevReac_)
      {
        if ( TC_reacRev_[irev] == i )
        {
          /* yes, reverse Arhenius parameters are given */
          indx = irev*3 ;
          TC_krev[i] = TC_reacArhenRev_[indx]
                      *exp( TC_reacArhenRev_[indx+1]*tln-TC_reacArhenRev_[indx+2]*t_1 ) ; 
          irev++ ;
        } /* done if section for reverse Arhenius parameters */
        else
        {
          /* no, need to compute equilibrium constant */
          double sumNuGk, kc ;
          int ir = TC_reacScoef_[i] ;
          if ( ir == -1 )
          {
            sumNuGk = TC_getSumNuGk(i,TC_gk) ;
          }
          else
          {
            sumNuGk = TC_getSumRealNuGk(i,ir,TC_gk) ;
          }
          kc         = TC_kc_coeff[i]*exp(sumNuGk) ;
          TC_krev[i] = TC_kfor[i]/kc ;

        } /* done if section for equilibrium constant */

      }
      else
      {
        /* no, need to compute equilibrium constant */
        double sumNuGk, kc ;
        int ir = TC_reacScoef_[i] ;
        if ( ir == -1 )
        {
          sumNuGk = TC_getSumNuGk(i,TC_gk) ;
        }
        else
        {
          sumNuGk = TC_getSumRealNuGk(i,ir,TC_gk) ;
        }
        kc         = TC_kc_coeff[i]*exp(sumNuGk) ;
        TC_krev[i] = TC_kfor[i]/kc ;

      } /* done if section for equilibrium constant */

    } /* done if reaction is reversible */

  } /* done computing kforward and kreverse rate constants */

  return ( 0 ) ;

}

int TC_getkForRevTab(double t1)
{
  int i ;
  double trat = (t1-TMIN)*TC_odelT_ ;

  int iTab = (int) trat ;
  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nreac_ ;

  for ( i = 0 ; i < TC_Nreac_ ; i++,iTab++ )
    TC_kfor[i] = TC_kfortab[iTab]+trat*(TC_kfortab[iTab+TC_Nreac_]-TC_kfortab[iTab]) ;

  iTab -= TC_Nreac_ ;
  for ( i = 0 ; i < TC_Nreac_ ; i++,iTab++ )
  {
    TC_krev[i] = 0.0 ;
    if ( TC_isRev_[i] )    
      TC_krev[i] = TC_krevtab[iTab]+trat*(TC_krevtab[iTab+TC_Nreac_]-TC_krevtab[iTab]) ;
  }

  return ( 0 ) ;

}

int TC_getkForRevP(double t1, double t_1)
{
  
  if ( TC_tab_ )
    TC_getkForRevPTab(t1) ;
  else
    TC_getkForRevPFcn(t_1) ;

  return ( 0 ) ;

}


int TC_getkForRevPFcn(double t_1)
{

  int indx, i, j, irev, iplog, plogtest ;
  double logP, ki, ki1, *api, *api1 ;

  if ( TC_nPlogReac_ > 0 ) logP = log(TC_pressure_/ATMPA) ;

  irev  = 0 ;
  iplog = 0 ;
  for ( i = 0 ; i<TC_Nreac_; i++ )
  {

    plogtest = 0 ;
    if ( iplog < TC_nPlogReac_ )
      plogtest = (i == TC_reacPlogIdx_[iplog]);
    
    if (plogtest) {
      indx = 4*TC_reacPlogPno_[iplog];
      api  = &(TC_reacPlogPars_[indx]) ;
      
      if (logP <= (*api)) {
        TC_kforP[i] = t_1 * ( *(api+2) + *(api+3) * t_1 )  ;
      } /* Done with branch p<pmin */
      
      else {
        indx = 4*(TC_reacPlogPno_[iplog+1]-1);
        api  = &(TC_reacPlogPars_[indx]) ;
        
        if (logP >= (*api)) {
          TC_kforP[i] = t_1 * ( *(api +2) + *(api +3) * t_1 ) ;
        } /* Done with branch p>pmax */
        
        else {
          
          for (j=TC_reacPlogPno_[iplog]; j<TC_reacPlogPno_[iplog+1]-1;j++) {
            indx = 4*j;
            if (logP<TC_reacPlogPars_[indx+4]) {
              api  = &(TC_reacPlogPars_[indx  ]);
              api1 = &(TC_reacPlogPars_[indx+4]);
              ki  = t_1 * ( *(api +2) + *(api +3) * t_1 ) ;
              ki1 = t_1 * ( *(api1+2) + *(api1+3) * t_1 ) ;
              TC_kforP[i] = ki+(logP - (*api))/( (*api1) - (*api))*(ki1-ki);
              break ;
            }
            
          } /* Done loop over all intervals */
          
        } /* Done with branch pmin<p<pmax */
        
      } /* Done with branch p>pmin */
      
      iplog++;
      
    } /* Done if reaction has a PLOG form */
    else {
      indx = i*3 ;
      TC_kforP[i] = t_1*(TC_reacArhenFor_[indx+1]+TC_reacArhenFor_[indx+2] * t_1) ;
    }
    TC_krevP[i] = 0.0 ;

    /* is reaction reversible ? */
    if ( TC_isRev_[i] )
    {
      /* are reverse Arhenius parameters given ? */
      if (irev < TC_nRevReac_)
      {
        if (TC_reacRev_[irev] == i)
        {
          /* yes, reverse Arhenius parameters are given */
          indx = irev*3 ;
          TC_krevP[i] = t_1*(TC_reacArhenRev_[indx+1]+TC_reacArhenRev_[indx+2]*t_1) ;
          irev++ ;
        } /* done if section for reverse Arhenius parameters */
        else
        { 
          /* no, need to compute krevP through derivatives of kc */
          double sumNuGkp ;
          if ( TC_reacScoef_[i] == -1)
            sumNuGkp = TC_getSumNuGk(i,TC_gkp) ;
          else
          {
            int ir = TC_reacScoef_[i] ;
            sumNuGkp = TC_getSumRealNuGk(i,ir,TC_gkp) ;
          }
          TC_krevP[i] = TC_kforP[i] - sumNuGkp ;

        } /* done if section for equilibrium constant */
  
      }
      else
      {
        /* no, need to compute krevP through derivatives of kc */
        double sumNuGkp ;
        if ( TC_reacScoef_[i] == -1)
          sumNuGkp = TC_getSumNuGk(i,TC_gkp) ;
        else
        {
          int ir = TC_reacScoef_[i] ;
          sumNuGkp = TC_getSumRealNuGk(i,ir,TC_gkp) ;
        }
        TC_krevP[i] = TC_kforP[i] - sumNuGkp ;
        //printf("xxxx:%d: %e\n",i+1,sumNuGkp);

      } /* done if section for equilibrium constant */

    } /* done if reaction is reversible */

  } /* done computing kforward and kreverse rate constants */

  return ( 0 ) ;

}

int TC_getkForRevPTab(double t1)
{

  int i ;

  double trat = (t1-TMIN)*TC_odelT_ ;

  int iTab = (int) trat ;
  iTab = MIN(MAX(iTab,0),TC_Ntab_-2) ;

  trat -= (double) iTab ;
  iTab *= TC_Nreac_ ;

  for ( i = 0 ; i < TC_Nreac_ ; i++,iTab++ )
    TC_kforP[i] = TC_kforPtab[iTab]+trat*(TC_kforPtab[iTab+TC_Nreac_]-TC_kforPtab[iTab]) ;

  iTab -= TC_Nreac_ ;
  for ( i = 0 ; i < TC_Nreac_ ; i++,iTab++ )
  {
    TC_krevP[i] = 0.0 ;
    if ( TC_isRev_[i] )    
      TC_krevP[i] = TC_krevPtab[iTab]+trat*(TC_krevPtab[iTab+TC_Nreac_]-TC_krevPtab[iTab]) ;
  }

  return ( 0 ) ;

}

int TC_getkForRevPder(double t_1, double tln)
{
  
  int indx, i, j, irev, iplog, plogtest ;
  double logP, ki, ki1, sumNuGk, kc, *api, *api1 ;
  
  if ( TC_nPlogReac_ > 0 ) logP = log(TC_pressure_/ATMPA) ;
  
  irev  = 0 ;
  iplog = 0 ;
  for ( i = 0 ; i<TC_Nreac_; i++ )
  {
    
    plogtest = 0 ;
    if ( iplog < TC_nPlogReac_ )
      plogtest = (i == TC_reacPlogIdx_[iplog]);
    
    if (plogtest) {
      
      indx = 4*TC_reacPlogPno_[iplog];
      api  = &(TC_reacPlogPars_[indx]) ;
      
      if (logP <= (*api)) 
        TC_kforPder[i] = 0.0  ;
        /* Done with branch p<pmin */
      
      else {
        
        indx = 4*(TC_reacPlogPno_[iplog+1]-1);
        api  = &(TC_reacPlogPars_[indx]) ;
        
        if (logP >= (*api)) 
          TC_kforPder[i] = 0.0 ;
          /* Done with branch p>pmax */
        
        else {
          
          for (j=TC_reacPlogPno_[iplog]; j<TC_reacPlogPno_[iplog+1]-1;j++) {
            indx = 4*j;
            if (logP<TC_reacPlogPars_[indx+4]) {
              api  = &(TC_reacPlogPars_[indx  ]);
              api1 = &(TC_reacPlogPars_[indx+4]);
              ki  = *(api +1) + *(api +2) * tln - *(api +3) * t_1 ;
              ki1 = *(api1+1) + *(api1+2) * tln - *(api1+3) * t_1 ;
              TC_kforPder[i] = (ki1-ki)/((*api1) - (*api))/TC_pressure_;
              break ;
            }
            
          } /* Done loop over all intervals */
          
        } /* Done with branch pmin<p<pmax */
        
      } /* Done with branch p>pmin */
      
      iplog++;
      
    } /* Done if reaction has a PLOG form */
    else {
      TC_kforPder[i] = 0.0;
    }
    TC_krevPder[i] = 0.0 ;
    
    /* is reaction reversible ? */
    if ( TC_isRev_[i] )
    {
      /* are reverse Arhenius parameters given ? */
      if (irev < TC_nRevReac_)
      {
        if (TC_reacRev_[irev] == i)
        {
          /* yes, reverse Arhenius parameters are given */
	  if (plogtest) {
            printf("Error: Reverse Arhnenius factors cannot be present for PLOG\n");
            exit(1);
	  }
          else
            TC_krevPder[i]=0.0; /* Currently given reverse rates do
                                   not depend on pressure */
        } /* done if section for reverse Arhenius parameters */
        else
        {
          /* no, need to compute krevP through derivatives of kc */
          if ( TC_reacScoef_[i] == -1)
            sumNuGk = TC_getSumNuGk(i,TC_gkp) ;
          else
          {
            int ir = TC_reacScoef_[i] ;
            sumNuGk = TC_getSumRealNuGk(i,ir,TC_gkp) ;
          }
          kc = TC_kc_coeff[i]*exp(sumNuGk) ;
          TC_krevPder[i] = TC_kforPder[i]/kc ;
          
        } /* done if section for equilibrium constant */
        
      }
      else
      {
        /* no, need to compute krevP through derivatives of kc */
        if ( TC_reacScoef_[i] == -1)
          sumNuGk = TC_getSumNuGk(i,TC_gkp) ;
        else
        {
          int ir = TC_reacScoef_[i] ;
          sumNuGk = TC_getSumRealNuGk(i,ir,TC_gkp) ;
        }
        kc = TC_kc_coeff[i]*exp(sumNuGk) ;
        TC_krevPder[i] = TC_kforPder[i]/kc ;
        
      } /* done if section for equilibrium constant */
      
    } /* done if reaction is reversible */
    
  } /* done computing kforward and kreverse rate constants */
  
  return ( 0 ) ;
  
}

/*
                   _                    __
         _ __ __ _| |_ ___        ___  / _|      _ __  _ __ ___   __ _ 
        | '__/ _` | __/ _ \_____ / _ \| |_ _____| '_ \| '__/ _ \ / _` |
        | | | (_| | ||  __/_____| (_) |  _|_____| |_) | | | (_) | (_| |
        |_|  \__,_|\__\___|      \___/|_|       | .__/|_|  \___/ \__, |
                                                |_|              |___/ 


 */
int TC_getRateofProg(double *concX)
{

  int i, j, indx, indxR, kspec, irnu, iord ;

  irnu = 0 ;
  iord = 0 ;
  for ( i = 0 ; i<TC_Nreac_; i++)
  {

    TC_ropFor[i] = TC_kfor[i] ;
    TC_ropRev[i] = TC_krev[i] ;
    
    /* compute forward rop */
    indx = i*TC_maxSpecInReac_ ;
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      int niup ;
      kspec = TC_reacSidx_[indx] ;
      niup  = abs(TC_reacNuki_[indx]) ;
      TC_ropFor[i] *= fastIntPow(concX[kspec],niup) ;
#ifdef NONNEG
      if (concX[kspec]<0) TC_ropFor[i] = 0.0 ;
#endif
      indx++ ;
    }

    if ( TC_isRev_[i] )
    {
      /* compute reverse rop */
      indx     = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
        int nius ;
        kspec = TC_reacSidx_[indx] ;
        nius  = TC_reacNuki_[indx] ;
        TC_ropRev[i] *= fastIntPow(concX[kspec],nius) ;
#ifdef NONNEG
        if (concX[kspec]<0) TC_ropRev[i] = 0.0 ;
#endif
        indx++ ;
      }
    }

    /* check for real stoichiometric coefficients */
    if ( irnu < TC_nRealNuReac_ )
      if ( TC_reacRnu_[irnu] == i )
      {
        TC_ropFor[i] = TC_kfor[i] ;
        TC_ropRev[i] = TC_krev[i] ;
    
        /* compute forward rop */
        indx  = i   *TC_maxSpecInReac_ ;
        indxR = irnu*TC_maxSpecInReac_ ;
        for ( j = 0; j<TC_reacNreac_[i] ; j++)
        {
    
          double niup ;
          kspec = TC_reacSidx_[indx] ;
          niup  = fabs(TC_reacRealNuki_[indxR]) ;
          TC_ropFor[i] *= pow(concX[kspec],niup) ;
#ifdef NONNEG
          if (concX[kspec]<0) TC_ropFor[i] = 0.0 ;
#endif
          indx++ ;
          indxR++ ;
        }

        if ( TC_isRev_[i] )
        {
          indx  = i   *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          indxR = irnu*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
          /* compute reverse rop */
          for ( j = 0; j<TC_reacNprod_[i] ; j++)
          {
            double nius ;
            kspec = TC_reacSidx_[indx] ;
            nius  = TC_reacRealNuki_[indxR] ;
            TC_ropRev[i] *= pow(concX[kspec],nius) ;
#ifdef NONNEG
            if (concX[kspec]<0) TC_ropRev[i] = 0.0 ;
#endif
            indx++ ;
            indxR++ ;
          }
        }

        irnu++ ;

      } /* done if real stoichiometric coefficients */

    /* check for arbitrary order reaction */
    if ( iord < TC_nOrdReac_ )
      if (TC_reacAOrd_[iord] == i)
      {
        /* found arbitrary order -> need to recompute rop's */
        TC_ropFor[i] = TC_kfor[i] ;
        TC_ropRev[i] = TC_krev[i] ;
        for ( j = 0; j<TC_maxOrdPar_ ; j++)
        {

          indx  = iord*TC_maxOrdPar_+j ;
          if (TC_specAOidx_[indx]<0) 
          {
            kspec = -TC_specAOidx_[indx]-1 ;
            TC_ropFor[i] *= pow(concX[kspec],TC_specAOval_[indx]) ;
#ifdef NONNEG
            if (concX[kspec]<0) TC_ropFor[i] = 0.0 ;
#endif
          }
          else if (TC_specAOidx_[indx]>0) 
          {
            kspec = TC_specAOidx_[indx]-1 ;
            TC_ropRev[i] *= pow(concX[kspec],TC_specAOval_[indx]) ;
#ifdef NONNEG
            if (concX[kspec]<0) TC_ropRev[i] = 0.0 ;
#endif
          }
        }

        iord++ ;

      } /* done if arbitrary order reaction */

    TC_rop[i] = TC_ropFor[i] - TC_ropRev[i] ;

  } /* done loop over all reactions */

  return ( 0 ) ;

}

int TC_getRateofProgDer(double *concX, int ireac, int ispec, double *qfr)
{

  int j, indx, indxR, kspec, irnu, iord ;

  qfr[0] = TC_kfor[ireac] ;
  qfr[1] = 0.0            ;
 
  /* compute forward qfor */
  indx = ireac*TC_maxSpecInReac_ ;
  for ( j = 0; j < TC_reacNreac_[ireac] ; j++)
  {
    int niup ;
    kspec = TC_reacSidx_[indx] ;
    niup  = abs(TC_reacNuki_[indx]) ;
    if ( ispec == kspec )
    {  
      if ( niup != 1 ) qfr[0] *= ((double) niup)*fastIntPow(concX[kspec],niup-1) ;
    }
    else
      qfr[0] *= fastIntPow(concX[kspec],niup) ;
    indx++ ;
  } /* done loop over reactants */
    
  if ( TC_isRev_[ireac] ) {
    /* compute reverse qrev */
    qfr[1] = TC_krev[ireac] ;
    indx     = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
    for ( j = 0; j<TC_reacNprod_[ireac] ; j++)
    {
      int nius ;
      kspec = TC_reacSidx_[indx] ;
      nius  = TC_reacNuki_[indx] ;
      if ( ispec == kspec ) {
        if ( nius != 1 ) qfr[1] *= ((double) nius)*fastIntPow(concX[kspec],nius-1) ;
      }
      else
        qfr[1] *= fastIntPow(concX[kspec],nius) ;
      indx++ ;
    }
  }

  /* Check for real stoichiometric coefficients */
  for ( irnu = 0 ; irnu < TC_nRealNuReac_ ; irnu++ )
    if ( TC_reacRnu_[irnu] == ireac ) {

      /* compute forward qfor */
      qfr[0] = TC_kfor[ireac] ;
      qfr[1] = 0.0            ;
      
      indx  = ireac*TC_maxSpecInReac_ ;
      indxR = irnu *TC_maxSpecInReac_ ;
      for ( j = 0; j<TC_reacNreac_[ireac] ; j++) {
    
        double niup ;
        kspec = TC_reacSidx_[indx] ;
        niup  = fabs(TC_reacRealNuki_[indxR]) ;
        if ( concX[kspec] > 0.0 ) {
	  if ( ispec == kspec ) {
	    if ( niup != 1.0 )
	      qfr[0] *= niup*pow(concX[kspec],niup-1) ;
	  }
	  else
	    qfr[0] *= pow(concX[kspec],niup) ;
	}
        else {
	  if ( ispec == kspec ) {
	    if ( niup != 1.0 ) 
	      qfr[0] *= niup*pow(1.e-20,niup-1) ;
	  }
	  else
	    qfr[0] *= pow(1.e-20,niup) ;
	}
	  
        indx++ ;
        indxR++ ;
      } /* Done loop over all reactants */

      if ( TC_isRev_[ireac] )
      {
        /* compute reverse qrev */
        qfr[1] = TC_krev[ireac] ;
        indx  = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
        indxR = irnu *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
        for ( j = 0; j<TC_reacNprod_[ireac] ; j++)
        {
          double nius ;
          kspec = TC_reacSidx_[indx] ;
          nius  = TC_reacRealNuki_[indxR] ;
	  if ( concX[kspec] > 0.0 ) {
	    if ( ispec == kspec ) {
	      if ( nius != 1.0 ) 
		qfr[1] *= nius*pow(concX[kspec],nius-1) ;
	    }
	    else
	      qfr[1] *= pow(concX[kspec],nius) ;
	  }
          else {
	    if ( ispec == kspec ) {
	      if ( nius != 1.0 ) 
		    qfr[1] *= nius*pow(1.e-20,nius-1) ;
	    }
	    else
	      qfr[1] *= pow(1.e-20,nius) ;
	  }
          indx++ ;
          indxR++ ;
        } /* Done loop over all products */
      } /* Done if real coef reac is rev */
    } /* Done loop over all reactions with real coeffs */
    
  
  /* check for arbitrary order reaction */
  for ( iord = 0 ; iord < TC_nOrdReac_; iord++ )
    if (TC_reacAOrd_[iord] == ireac)
    {
      /* found arbitrary order -> need to recompute qfor and qrev */
      qfr[0] = TC_kfor[ireac] ;
      qfr[1] = TC_krev[ireac] ;
      for ( j = 0; j<TC_maxOrdPar_ ; j++)
      {
        indx  = iord*TC_maxOrdPar_+j ;
        if (TC_specAOidx_[indx]<0) 
        {
          kspec = -TC_specAOidx_[indx]-1 ;
         
          if ( ispec == kspec )
	  {
	    if ( concX[kspec] > 0.0 ) {
	      if (TC_specAOval_[indx] != 1.0)
		qfr[0] *= TC_specAOval_[indx]*pow(concX[kspec],TC_specAOval_[indx]-1.0);
	    } else {
	      if (TC_specAOval_[indx] != 1.0)
		qfr[0] *= TC_specAOval_[indx]*pow(1.e-20,TC_specAOval_[indx]-1.0);
	    }
	  }
	  else {
	    if ( concX[kspec] > 0.0 ) 
	      qfr[0] *= pow(concX[kspec],TC_specAOval_[indx]);
	    else
	      qfr[0] *= pow(1.e-20,TC_specAOval_[indx]);
	  }

        }
        else if (TC_specAOidx_[indx]>0) 
        {
          kspec = TC_specAOidx_[indx]-1 ;
          if ( ispec == kspec )
	  {
	    if ( concX[kspec] > 0.0 ) {
	      if (TC_specAOval_[indx] != 1.0)
		qfr[1] *= TC_specAOval_[indx]*pow(concX[kspec],TC_specAOval_[indx]-1.0) ;
	    } else {
	      if (TC_specAOval_[indx] != 1.0)
		qfr[1] *= TC_specAOval_[indx]*pow(1.e-20,TC_specAOval_[indx]-1.0) ;
	    }
	  }
          else {
	    if ( concX[kspec] > 0.0 ) {
	      qfr[1] *= pow(concX[kspec],TC_specAOval_[indx]) ;
	    } else {
	      qfr[1] *= pow(1.e-20,TC_specAOval_[indx]) ;
	    }
	  }
        }
      }
    } /* done if arbitrary order reaction */

  return ( 0 ) ;

}
/*
         ____                           _    __            _             
        / ___|_ __ ___  _   _ _ __   __| |  / _| __ _  ___| |_ ___  _ __ 
       | |   | '__/ _ \| | | | '_ \ / _` | | |_ / _` |/ __| __/ _ \| '__|
       | |___| | | (_) | |_| | | | | (_| | |  _| (_| | (__| || (_) | |   
        \____|_|  \___/ \__,_|_| |_|\__,_| |_|  \__,_|\___|\__\___/|_|   

*/
int TC_getCrnd(double t1,double t_1,double tln,double *concX, double *concM)
{

  int i, ipfal ;

  ipfal = 0 ;
  for ( i = 0 ; i<TC_Nreac_; i++)
  {

    TC_Crnd[i] = concM[i] ;

    if ( ipfal < TC_nFallReac_ )
      if (TC_reacPfal_[ipfal] == i)
      {

        double *Arh0, k0, kinf, Pr, logPr ;
        int    indx0 = ipfal * TC_nFallPar_ ;
        if (TC_reacPlohi_[ipfal] == 0) {
          /* LOW reaction */
          Arh0 = &(TC_reacPpar_[indx0]) ;
          k0    = Arh0[0]*exp(Arh0[1]*tln-Arh0[2]*t_1) ;
          Pr    = k0/TC_kfor[i] ;
          if ( TC_reacPspec_[ipfal] >= 0 )
            Pr *= concX[TC_reacPspec_[ipfal]] ;
          else 
            Pr *= concM[i] ;
          TC_Crnd[i] = Pr / ( 1.0 + Pr ) ; /* At least Lindemann form */
	} else {
          /* HIGH reaction */
          Arh0 = &(TC_reacPpar_[indx0]) ;
          kinf = Arh0[0]*exp(Arh0[1]*tln-Arh0[2]*t_1) ;
          Pr   = TC_kfor[i]/kinf ;
          if ( TC_reacPspec_[ipfal] >= 0 )
            Pr *= concX[TC_reacPspec_[ipfal]] ;
          else 
            Pr *= concM[i] ;
          TC_Crnd[i] = 1.0 / ( 1.0 + Pr ) ; /* At least Lindemann form */
	}

        logPr = log10(MAX(Pr,TCSMALL)) ;

        if ( TC_reacPtype_[ipfal] == 2 )
        {

          /*
                  _______  ____
                 / __/ _ \/  _/
                _\ \/ , _// /  
               /___/_/|_/___/  
                                     
          */
          double *Psri, Xpres, Ffac ;

          indx0 += 3 ;
          Psri  = &(TC_reacPpar_[indx0]) ;
          Xpres = 1.0/(1.0 + logPr*logPr) ;
          Ffac  = pow(Psri[0]*exp(-Psri[1]*t_1)+exp(-t1/Psri[2]),Xpres)
            *Psri[3]*pow(t1,Psri[4]) ;

          TC_Crnd[i] *= Ffac ;

        } /* done with SRI form */

        else if ( TC_reacPtype_[ipfal] >= 3 )
        {

          /*
                  ______            
                 /_  __/______  ___ 
                  / / / __/ _ \/ -_)
                 /_/ /_/  \___/\__/ 
                                             
          */
          double *Ptroe, Fc, logFc, Atroe, Btroe, logFfac, Ffac, Atroe_Btroe ;
          indx0 += 3;
          Ptroe = &(TC_reacPpar_[indx0]) ;

	  Fc = 0.0 ;
	  if ( fabs(1.0-Ptroe[0]) > 0.0 ) Fc += (1.0-Ptroe[0])*exp(-t1/Ptroe[1]);
	  if ( fabs(    Ptroe[0]) > 0.0 ) Fc +=      Ptroe[0] *exp(-t1/Ptroe[2]) ;

          Fc = (1.0-Ptroe[0])*exp(-t1/Ptroe[1])+Ptroe[0]*exp(-t1/Ptroe[2]) ;

          if ( TC_reacPtype_[ipfal] == 4 ) Fc += exp(-Ptroe[3]*t_1) ;

          logFc   = log10(Fc) ;
          Atroe   = logPr-0.40-0.67*logFc ;
          Btroe   = 0.75 - 1.27 * logFc - 0.14 * Atroe;
          Atroe_Btroe = Atroe / Btroe;
          logFfac = logFc / ( 1.0 + Atroe_Btroe*Atroe_Btroe ) ;
          Ffac    = pow(10.0,logFfac) ;

          TC_Crnd[i] *= Ffac ;

        } /* done with Troe form */

        ipfal++ ;

      } /* done if pressure dependent reaction */

  }

  return ( 0 ) ;

}

/*
       ____                           _       _           _           
      / ___|_ __ ___  _   _ _ __   __| |   __| | ___ _ __(_)_   _____ 
     | |   | '__/ _ \| | | | '_ \ / _` |  / _` |/ _ \ '__| \ \ / / __|
     | |___| | | (_) | |_| | | | | (_| | | (_| |  __/ |  | |\ V /\__ \
      \____|_|  \___/ \__,_|_| |_|\__,_|  \__,_|\___|_|  |_| \_/ |___/

*/
int TC_getCrndDer(int ireac, int *itbdy, int *ipfal, 
      double t1,double t_1,double tln,double *concX, double *concM)
{
  int i, j ;
  for ( i = 0 ; i < TC_Nspec_+1 ; i++) TC_CrndDer[i] = 0.0 ;

  /* Check if reaction involves third-body */
  if ( (*itbdy) < TC_nThbReac_ )
  {

    if ( ireac == TC_reacTbdy_[*itbdy] )
    {
      for ( i = 1 ; i < TC_Nspec_+1 ; i++) TC_CrndDer[i] = 1.0 ;

      for ( j = 0 ; j < TC_reacTbno_[*itbdy] ; j++ )
      {
        int indx, kspec ;
        indx  = (*itbdy)*TC_maxTbInReac_+j ;
        kspec = TC_specTbdIdx_[indx] ;
        TC_CrndDer[kspec+1] = TC_specTbdEff_[indx]  ;

      } /* done loop over efficiencies */
      
      (*itbdy) += 1 ;

    } /* done if reaction involves third-body */

  } /* done if there are third-body reaction left */

  /* Check if reaction is pressure dependent */
  if ( (*ipfal) < TC_nFallReac_ )
    if ( ireac == TC_reacPfal_[ (*ipfal) ] )
    {
      /* Compute Pr and its derivatives */
      int indx0, indxI ;
      double *Arh0, *ArhI, k0, kinf, Pr ;
      indx0 = (*ipfal) * TC_nFallPar_ ;
      indxI = ( ireac) * 3          ;
      Arh0 = &(TC_reacPpar_[indx0]) ;
      ArhI = &(TC_reacArhenFor_[indxI]) ;
      if (TC_reacPlohi_[*ipfal] == 0) {
        /* LOW reaction */
        k0   = Arh0[0]*exp(Arh0[1]*tln-Arh0[2]*t_1) ;
        Pr   = k0/TC_kfor[ireac] ;
      } else {
        /* HIGH reaction */
        kinf = Arh0[0]*exp(Arh0[1]*tln-Arh0[2]*t_1) ;
        Pr   = TC_kfor[ireac]/kinf ;
      }
      for ( i = 1 ; i < TC_Nspec_+1 ; i++) TC_PrDer[i] = 0.0 ;
      
      if ( TC_reacPspec_[(*ipfal)] >= 0 )
      {
        TC_PrDer[TC_reacPspec_[(*ipfal)]+1] = Pr ;
        Pr *= concX[TC_reacPspec_[(*ipfal)]] ;
      }
      else 
      {
        for ( i = 1 ; i<TC_Nspec_+1 ; i++ ) TC_PrDer[i] = Pr*TC_CrndDer[i] ;
        Pr *= concM[ireac] ;
      }

      if (TC_reacPlohi_[*ipfal] == 0) {
        /* LOW reaction */
        TC_PrDer[0] = Pr * t_1 * ( Arh0[1] - ArhI[1] + t_1 * ( Arh0[2]-ArhI[2] ) ) ;
      } else {
        /* HIGH reaction */
        TC_PrDer[0] = Pr * t_1 * ( ArhI[1] - Arh0[1] + t_1 * ( ArhI[2]-Arh0[2] ) ) ;
      }
      
      if ( TC_reacPtype_[(*ipfal)] == 1 )
      {

        /*
            __   _         __                       
           / /  (_)__  ___/ /__ __ _  ___ ____  ___ 
          / /__/ / _ \/ _  / -_)  ' \/ _ `/ _ \/ _ \
         /____/_/_//_/\_,_/\__/_/_/_/\_,_/_//_/_//_/

        */
        double Prfac = 1.0/((1.0+Pr)*(1.0+Pr)) ;
        for ( i = 0 ; i<TC_Nspec_+1 ; i++ ) TC_CrndDer[i] = Prfac*TC_PrDer[i] ;

      } /* done with Lindemann form */

      else if ( TC_reacPtype_[(*ipfal)] == 2 )
      {

        /*
              _______  ____
            / __/ _ \/  _/
           _\ \/ , _// /  
          /___/_/|_/___/  
                                 
        */
        double *Psri, logPr, Xp, dXp, abcS, Ffac, Prfac ;
  
        indx0 += 3;
        Psri  = &(TC_reacPpar_[indx0]) ;
        logPr = log(Pr)/log(10.0) ;
        Xp    = 1.0/(1.0 + logPr*logPr) ;
        dXp   = -Xp*Xp*2.0*logPr / ( Pr * log(10.0) ) ;
        abcS  = Psri[0]*exp(-Psri[1]*t_1)+exp(-t1/Psri[2]) ;
        Ffac  = pow(abcS,Xp)*Psri[3]*pow(t1,Psri[4]) ;
        TC_dFfac[0] = Ffac*(Psri[4]*t_1+dXp*TC_PrDer[0]*log(abcS)+
                     Xp*(Psri[0]*Psri[1]*t_1*t_1*exp(-Psri[1]*t_1)
                    -exp(-t1/Psri[2])/Psri[2])/abcS) ;
        abcS = log(abcS)*dXp ;
        for ( i = 1 ; i<TC_Nspec_+1 ; i++ ) 
          TC_dFfac[i] = Ffac*abcS*TC_PrDer[i] ;

        Ffac = Ffac/((1.0+Pr)*(1.0+Pr)) ;
        Prfac = Pr/(1.0+Pr) ;
        for ( i = 0 ; i<TC_Nspec_+1 ; i++ ) 
          TC_CrndDer[i] = Ffac*TC_PrDer[i]+Prfac*TC_dFfac[i] ;

      } /* done with SRI form */

      else if ( TC_reacPtype_[(*ipfal)] >= 3 )
      {

        /*
          ______            
         /_  __/______  ___ 
          / / / __/ _ \/ -_)
         /_/ /_/  \___/\__/ 
                                             
        */
        double *Ptroe, Fc1, Fc2, Fc, FcDer ;
        double logFc, logPr, Atroe, Btroe, oABtroe, logFfac, Ffac, Atroe_Btroe;
        double oPr, oFc, Afc, Bfc, Apr, Bpr, Gfac, Prfac ;

        indx0 += 3 ;
        Ptroe = &(TC_reacPpar_[indx0]) ;
	Fc1 = Fc2 = FcDer = 0.0 ;
	if ( 1.0-Ptroe[0] > 0.0 ) Fc1 = (1.0-Ptroe[0])*exp(-t1/Ptroe[1]) ;
        if (     Ptroe[0] > 0.0 ) Fc2 =      Ptroe[0] *exp(-t1/Ptroe[2]) ;
        Fc  =  Fc1 + Fc2 ;
	if ( 1.0-Ptroe[0] > 0.0 ) FcDer -= Fc1/Ptroe[1] ;
        if (     Ptroe[0] > 0.0 ) FcDer -= Fc2/Ptroe[2] ;
  
        if ( TC_reacPtype_[(*ipfal)] == 4 ) 
        {
          double Fc3 = exp(-Ptroe[3]*t_1) ;
          Fc    += Fc3 ;
          FcDer += Fc3*Ptroe[3]*t_1*t_1 ;
        }

        logFc   = log(Fc) / log(10.0) ;

	if ( Pr > 0.0 )
	{
          logPr   = log(Pr) / log(10.0) ;
	  Atroe   = logPr-0.40-0.67*logFc ;
	  Btroe   = 0.75 - 1.27 * logFc - 0.14 * Atroe;
	  Atroe_Btroe = Atroe / Btroe;
	}
	else
	  Atroe_Btroe = -1.0/0.14;

        oABtroe = 1.0/( 1.0 + Atroe_Btroe*Atroe_Btroe ) ;
        logFfac = logFc * oABtroe ; 
        Ffac    = pow(10.0,logFfac) ;

	if ( Pr > 0.0 )
	{
	  oPr = 1.0/(Pr*log(10.0)) ;
	  oFc = 1.0/(Fc*log(10.0)) ;
	  Afc = -0.67  *oFc ;
	  Bfc = -1.1762*oFc ;
	  Apr =  1.0   *oPr ;
	  Bpr = -0.14  *oPr ;
	  Gfac = Ffac*log(Fc)*2.0*Atroe/(Btroe*Btroe*Btroe)*oABtroe*oABtroe  ;
	  
	  /* dF/dPr */
	  TC_dFfac[0] = -Gfac * ( Apr*Btroe - Bpr*Atroe ) ;
	  /* dF/dFc */
	  TC_dFfac[1] = Ffac/Fc*oABtroe-Gfac * ( Afc*Btroe - Bfc*Atroe ) ;
	}
	else
	{
	  TC_dFfac[0] = 0.0;
	  TC_dFfac[1] = 0.0;
	}

        Prfac = Pr/(1.0+Pr) ;
        Ffac = Ffac/((1.0+Pr)*(1.0+Pr))+Prfac*TC_dFfac[0] ;
  
        TC_CrndDer[0] = Ffac*TC_PrDer[0]+Prfac*TC_dFfac[1]*FcDer ;
        for ( i = 1 ; i<TC_Nspec_+1 ; i++ ) 
          TC_CrndDer[i] = Ffac*TC_PrDer[i] ;    

      } /* done with Troe form */

      (*ipfal) += 1 ;

    } /* done if reaction is pressure dependent */

  return ( 0 ) ;

}
