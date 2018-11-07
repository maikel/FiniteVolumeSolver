#include "TC_defs.h"
#include "TC_interface.h"

#include <stdlib.h>
#include <string.h>

/*! \file TC_chg.c
    \brief Functions for changing Arrhenius rate factors
*/ 
/**
 * \brief Change parameters for forward rate constants
 */
int TC_chgArhenFor( int ireac, int ipos, double newval )
{
/**
 * \param ireac : reaction index
 * \param ipos : index of parameter to be changed (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param newval : new parameter value
 */

  if ( TC_ArhenForChg_ == 0 )
  {

    /* Save original factors */
    TC_reacArhenForSave_ = (double *) malloc( TC_Nreac_*3*sizeof(double) ) ;
    memcpy(TC_reacArhenForSave_, TC_reacArhenFor_, 3*TC_Nreac_*sizeof(double) );
    TC_ArhenForChg_ = 1 ;

  }

  if ( ( ireac < 0 ) || ( ireac >= TC_Nreac_ ))
  {
    printf("TC_chgArhenFor() : reaction index is wrong -> %d\n", ireac) ;
    exit(1) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    printf("TC_chgArhenFor() : index of Arrhenius parameter is wrong -> %d\n", ipos) ;
    exit(1) ;
  }

  TC_reacArhenFor_[ireac*3+ipos] = newval ;

  return ( 0 ) ;

}

/**
 * \brief Change parameters for forward rate constants of low/high pressure limits for Pressure-dependent reactions
 */
int TC_chgArhenPresDepFor( int ireac, int ipos, double newval)
{
/**
 * \param ireac : Pressure-dependent reaction index
 * \param ipos : index of parameter to be changed (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param newval : new parameter value
 */

  if ( ( ireac < 0 ) || ( ireac >= TC_nFallReac_ ))
  {
    printf("TC_chgArhenPresDepFor() : reaction index is wrong -> %d\n", ireac) ;
    exit(1) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= TC_nFallPar_ ))
  {
    printf("TC_chgArhenPresDepFor() : index of Arrhenius parameter is wrong -> %d\n", ipos) ;
    exit(1) ;
  }

  TC_reacPpar_[ireac*TC_nFallPar_ + ipos] = newval ;

  return ( 0 ) ;

} /* Done with TC_chgArhenPresDepFor */

/**
 * \brief Reverse changes for forward rate constants' parameters
 */
int TC_chgArhenForBack( int ireac, int ipos)
{
/**
 * \param ireac : reaction index
 * \param ipos : index of parameter to be changed (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 */

  if ( TC_ArhenForChg_ == 0 )
  {

    printf("TC_chgArhenForBack() : Arrhenius paremeters were not backed up !!!\n") ;
    exit(1) ;

  }

  if ( ( ireac < 0 ) || ( ireac >= TC_Nreac_ ))
  {
    printf("TC_chgArhenForBack() : reaction index is wrong -> %d\n", ireac) ;
    exit(1) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    printf("TC_chgArhenForBack() : index of Arrhenius parameter is wrong -> %d\n", ipos) ;
    exit(1) ;
  }

  TC_reacArhenFor_[ireac*3+ipos] = TC_reacArhenForSave_[ireac*3+ipos] ;

  return ( 0 ) ;

}

/**
 * \brief Change parameters for reverse rate constants
 */
int TC_chgArhenRev( int ireac, int ipos, double newval )
{
/**
 * \param ireac : reaction index
 * \param ipos : index of parameter to be changed (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 * \param newval : new parameter value
 */

  int i, irevindx ;

  if ( TC_ArhenRevChg_ == 0 )
  {

    /* Save original factors */
    TC_reacArhenRevSave_ = (double *) malloc( TC_nRevReac_*3* sizeof(double) ) ;
    memcpy(TC_reacArhenRevSave_, TC_reacArhenRev_, 3*TC_nRevReac_*sizeof(double) );

    TC_ArhenRevChg_ = 1 ;

  }

  irevindx = -1 ;
  for ( i = 0 ; i < TC_nRevReac_; i++ )
    if ( TC_reacRev_[i] == ireac ) { irevindx = i; break ;}

  if ( irevindx < 0 )
  {
    printf("TC_chgArhenRev() : reaction %d does not have reversible Arrh. params. given -> \n", 
	   ireac) ;
    exit(1) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    printf("TC_chgArhenRev() : index of Arrhenius parameter %d is wrong\n", ipos) ;
    exit(1) ;
  }

  TC_reacArhenRev_[irevindx*3+ipos] = newval ;

  return ( 0 ) ;

}

/**
 * \brief Reverse changes for reverse rate constants' parameters
 */
int TC_chgArhenRevBack( int ireac, int ipos)
{
/**
 * \param ireac : reaction index
 * \param ipos : index of parameter to be changed (0) pre-exponential factor
                 (1) temperature exponent, (2) activation energy
 */

  int i, irevindx ;

  if ( TC_ArhenRevChg_ == 0 )
  {

    printf("TC_chgArhenRevBack() : Arrhenius paremeters were not saved !!!\n") ;
    exit(1) ;

  }

  irevindx = -1 ;
  for ( i = 0 ; i < TC_nRevReac_; i++ )
    if ( TC_reacRev_[i] == ireac ) { irevindx = i; break ;}

  if ( irevindx < 0 )
  {
    printf("TC_chgArhenRevBack() : reaction %d does not have reversible Arrh. params. given -> \n", 
	   ireac) ;
    exit(1) ;
  }
  if ( ( ipos < 0 ) || ( ipos >= 3 ))
  {
    printf("TC_chgArhenRevBack() : index of Arrhenius parameter %d is wrong\n", ipos) ;
    exit(1) ;
  }

  TC_reacArhenRev_[irevindx*3+ipos] = TC_reacArhenRevSave_[irevindx*3+ipos] ;

  return ( 0 ) ;

}
