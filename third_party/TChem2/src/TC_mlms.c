#include "TC_defs.h"
#include "TC_interface.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*! \file TC_mlms.c
    \brief Mass fractions - Mole fractions - Molar concetrations - Molecular weight
*/ 

/*
                 _   __  __     ____   ____     
       __ _  ___| |_|  \/  |___|___ \ / ___|___ 
      / _` |/ _ \ __| |\/| / __| __) | |   / __|
     | (_| |  __/ |_| |  | \__ \/ __/| |__| (__ 
      \__, |\___|\__|_|  |_|___/_____|\____\___|
      |___/                                     

*/
/**
 * \brief Computes molar concentrations based on temperature and species mass fractions. 
 *        Input temperature is normalized, output concentrations are also normalized before exit.
 *        \f[ \overline{[X_k]}=[X_k]\cdot\frac{W_{ref}}{\rho_{ref}}=Y_k\cdot\frac{\rho}{W_k}\cdot\frac{W_{ref}}{\rho_{ref}} \f]
 */
int TCDND_getMs2Cc(double *scal,int Nvars,double *concX)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K], mass fractions Y [] 
  \param Nvars : no. of variables = N<sub>spec</sub>+1 
  \return concX : array of doubles containing species molar concentrations [kmol/m<sup>3</sup>]
*/

  int ans, i ;

  if ( Nvars != TC_Nvars_ )
  {
    printf("TCDND_getMs2Cc() : disagreement in Nvars : %d vs %d\n",Nvars,TC_Nvars_) ;
    printf("                   -> Abort !!!\n") ;
    exit(1) ;
  }

  if ( TC_nonDim_ == 1 ) scal[0] *= TC_Tref_ ;

  ans = TC_getMs2Cc(scal,Nvars,concX) ;

  if ( TC_nonDim_ == 1 ) 
  {
    double makeND = TC_Wref_ / TC_rhoref_ ;
    scal[0] /= TC_Tref_ ;
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) concX[i] *= makeND ;
  }

  return ( ans ) ;

}
/**
 * \brief Computes molar concentrations based on temperature and species mass fractions. 
 * \f[ [X_k]=Y_k\cdot\rho/W_k \f]
 */
int TC_getMs2Cc(double *scal,int Nvars,double *concX)
{
/**
  \param scal : array of N<sub>spec</sub>+1 doubles (T,Y<sub>1</sub>,Y<sub>2</sub>,...,Y<sub>Nspec</sub>), temperature T [K], mass fractions Y [] 
  \param Nvars : no. of variables = N<sub>spec</sub>+1 
  \return concX : array of doubles containing species molar concentrations [kmol/m<sup>3</sup>]
*/

  int ans=0, i ;
  double rhomix ;
  double *Yspec ;

  if ( Nvars != TC_Nvars_ )
  {
    printf("TC_getMs2Cc() : disagreement in Nvars : %d vs %d\n",Nvars,TC_Nvars_) ;
    printf("                -> Abort !!!\n") ;
    exit(1) ;
  }

  /* compute density */
  if ( TC_rhoset_ == 1 )
    rhomix = TC_rho_;
  else
    ans = TC_getRhoMixMs(scal,Nvars,&rhomix) ;

  /* Compute molar mass */
  Yspec = &scal[1] ;
  for (i = 0 ; i < TC_Nspec_ ; i++ ) concX[i] = Yspec[i]*rhomix/TC_sMass_[i] ;

  return ( ans ) ;

}
/*
                _   __  __ _ ____  __  __     
      __ _  ___| |_|  \/  | |___ \|  \/  |___ 
     / _` |/ _ \ __| |\/| | | __) | |\/| / __|
    | (_| |  __/ |_| |  | | |/ __/| |  | \__ \
     \__, |\___|\__|_|  |_|_|_____|_|  |_|___/
     |___/                                    

*/
/**
 * \brief Transforms mole fractions to mass fractions (same as TC_getMl2Ms()).
 * \f[ Y_k=X_k\cdot W_k/W_{mix} \f]
 */
int TCDND_getMl2Ms(double *Xspec,int Nspec,double *Yspec)
{
/**
  \param Xspec : array of N<sub>spec</sub> mole fractions 
  \param Nspec : no. of species
  \return Yspec : array of N<sub>spec</sub> mass fractions
*/

  int ans ;
  if ( Nspec != TC_Nspec_ )
  {
    printf("TCDND_getMl2Ms() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                   -> Abort !!!\n") ;
    exit(1) ;
  }

  ans = TC_getMl2Ms(Xspec, Nspec, Yspec) ;

  return ( 0 ) ;

}
/**
 *  \brief Transforms mole fractions to mass fractions. 
 *  \f[ Y_k=X_k\cdot W_k/W_{mix} \f]
 */
int TC_getMl2Ms(double *Xspec,int Nspec,double *Yspec)
{
/**
  \param Xspec : array of N<sub>spec</sub> mole fractions 
  \param Nspec : no. of species
  \return Yspec : array of N<sub>spec</sub> mass fractions
*/

  int ans, i ;
  double Wmix ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TC_getMl2Ms() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                -> Abort !!!\n") ;
    exit(1) ;
  }

  /* Compute molar mass */
  ans = TC_getMl2Wmix(Xspec, Nspec, &Wmix) ;

  for (i = 0 ; i < TC_Nspec_ ; i++ ) Yspec[i] = Xspec[i]*TC_sMass_[i]/Wmix ;
  /*
  for (int i = 0 ; i < TC_Nspec_ ; i++ )
    std::cout<<"get_mole2mass() "<<i<<","<<Yspec[i]<<", "<<Xspec[i]<<", "<<TC_sMass_[i]
  	     <<std::endl<<std::flush; 
  */

  return ( ans ) ;

}
/*
                _   __  __     ____  __  __ _ 
      __ _  ___| |_|  \/  |___|___ \|  \/  | |
     / _` |/ _ \ __| |\/| / __| __) | |\/| | |
    | (_| |  __/ |_| |  | \__ \/ __/| |  | | |
     \__, |\___|\__|_|  |_|___/_____|_|  |_|_|
     |___/                                    
*/
/**
 *  \brief Transforms mass fractions to mole fractions (same as TC_getMs2Ml())
 *  \f[ X_k=Y_k\cdot W_{mix}/W_k \f]
 */
int TCDND_getMs2Ml(double *Yspec,int Nspec,double *Xspec)
{
/**
  \param Yspec : array of N<sub>spec</sub> mass fractions
  \param Nspec : no. of species
  \return Xspec : array of N<sub>spec</sub> mole fractions 
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TCDND_getMs2Ml() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                   -> Abort !!!\n") ;
    exit(1) ;
  }

  ans = TC_getMs2Ml(Yspec,Nspec,Xspec) ;

  return ( ans ) ;

}
/**
 *  \brief Transforms mass fractions to mole fractions. 
 *  \f[ X_k=Y_k\cdot W_{mix}/W_k \f]
 */
int TC_getMs2Ml(double *Yspec,int Nspec,double *Xspec)
{
/**
  \param Yspec : array of N<sub>spec</sub> mass fractions
  \param Nspec : no. of species
  \return Xspec : array of N<sub>spec</sub> mole fractions 
*/

  int ans, i ;
  double Wmix ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TC_getMs2Ml() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                -> Abort !!!\n") ;
    exit(1) ;
  }

  /* Compute molar mass */
  ans = TC_getMs2Wmix(Yspec,Nspec,&Wmix) ;

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) Xspec[i] = Yspec[i]*Wmix/TC_sMass_[i] ;

  return ( ans ) ;

}

/*
            _   __  __     ______        __         _      
  __ _  ___| |_|  \/  |___|___ \ \      / / __ ___ (_)_  __
 / _` |/ _ \ __| |\/| / __| __) \ \ /\ / / '_ ` _ \| \ \/ /
| (_| |  __/ |_| |  | \__ \/ __/ \ V  V /| | | | | | |>  < 
 \__, |\___|\__|_|  |_|___/_____| \_/\_/ |_| |_| |_|_/_/\_\
 |___/                                                     

*/
/**
 *  \brief Computes mixture molecular weight based on species mass fractions. 
 *         The molecular weight ([kg/kmol]=[g/mol]) is normalized before output.
 *         \f[ \bar{W}_{mix}=\frac{W_{mix}}{W_{ref}}=\frac{1}{W_{ref}}\left(\sum_{k=1}^{N_{spec}}Y_k/W_k\right)^{-1} \f]
*/
int TCDND_getMs2Wmix(double *Yspec,int Nspec,double *Wmix)
{
/**
  \param Yspec : array of N<sub>spec</sub> mass fractions
  \param Nspec : no. of species
  \return Wmix : pointer to normalized mixture molecular weight
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TCDND_getMs2Wmix() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                     -> Abort !!!\n") ;
    exit(1) ;
  }

  ans = TC_getMs2Wmix( Yspec, Nspec, Wmix) ;

  if ( TC_nonDim_ == 1 ) (*Wmix) /= TC_Wref_ ;

  return ( ans ) ;

}

/**
 *  \brief Computes mixture molecular weight based on species mass fractions. 
 *         \f[ W_{mix}=\left(\sum_{k=1}^{N_{spec}}Y_k/W_k\right)^{-1} \f]
*/
int TC_getMs2Wmix(double *Yspec,int Nspec,double *Wmix)
{
/**
  \param Yspec : array of N<sub>spec</sub> mass fractions
  \param Nspec : no. of species
  \return Wmix : pointer to mixture molecular weight [kg/kmol]=[g/mol]
*/

  int i ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TC_getMs2Wmix() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                  -> Abort !!!\n") ;
    exit(1) ;
  }

  (*Wmix) = 0.0 ;
  for (i = 0 ; i < TC_Nspec_ ; i++ ) (*Wmix) += Yspec[i]/TC_sMass_[i] ;

  (*Wmix) = 1.0/(*Wmix) ;

  return ( 0 ) ;

}
/*
             _   __  __ _ ______        __         _      
   __ _  ___| |_|  \/  | |___ \ \      / / __ ___ (_)_  __
  / _` |/ _ \ __| |\/| | | __) \ \ /\ / / '_ ` _ \| \ \/ /
 | (_| |  __/ |_| |  | | |/ __/ \ V  V /| | | | | | |>  < 
  \__, |\___|\__|_|  |_|_|_____| \_/\_/ |_| |_| |_|_/_/\_\
  |___/                                                   

*/
/**
 *  \brief Computes mixture molecular weight based on species mole fractions. 
 *         The molecular weight ([kg/kmol]=[g/mol]) is normalized before output.
 *         \f[ \bar{W}_{mix}=\frac{W_{mix}}{W_{ref}}=\frac{1}{W_{ref}}\sum_{k=1}^{N_{spec}}X_kW_k\f]
 */
int TCDND_getMl2Wmix(double *Xspec,int Nspec,double *Wmix)
{
/**
  \param Xspec : array of N<sub>spec</sub> mole fractions
  \param Nspec : no. of species
  \return Wmix : pointer to normalized mixture molecular weight 
*/

  int ans ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TCDND_getMl2Wmix() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                     -> Abort !!!\n") ;
    exit(1) ;
  }

  ans = TC_getMs2Wmix( Xspec, Nspec, Wmix) ;

  if ( TC_nonDim_ == 1 ) (*Wmix) /= TC_Wref_ ;

  return ( ans ) ;

}
/**
 *  \brief Computes mixture molecular weight based on species mole fractions. 
 *         \f[ W_{mix}=\sum_{k=1}^{N_{spec}}X_kW_k\f]
 */
int TC_getMl2Wmix(double *Xspec,int Nspec,double *Wmix)
{
/**
  \param Xspec : array of N<sub>spec</sub> mole fractions
  \param Nspec : no. of species
  \return Wmix : pointer to mixture molecular weight [kg/kmol]=[g/mol]
*/

  int i ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("TC_getMl2Wmix() : disagreement in Nspec : %d vs %d\n",Nspec,TC_Nspec_) ;
    printf("                  -> Abort !!!\n") ;
    exit(1) ;
  }

  (*Wmix) = 0.0 ;
  for (i = 0 ; i < TC_Nspec_ ; i++ ) (*Wmix) += Xspec[i]*TC_sMass_[i] ;

  return ( 0 ) ;

}
