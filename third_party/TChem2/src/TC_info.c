#include "TC_defs.h"

#include <stdlib.h>
#include <string.h>

/*! \file TC_info.c
 *  \brief Various functions providing info about the kinetic model
 */ 

/**
 * \brief Returns number of reactions \f$N_{reac}\f$
*/
int TC_getNreac() { return ( TC_Nreac_ ) ; }   

/**
 * \brief Returns number of reactions with low/high pressure limit
 * rate constants
*/
int TC_getNPresDepReac() { return ( TC_nFallReac_ ) ; }

/**
 * \brief Returns number of reactions with reaction rates enhanced 
 * by third-body efficiencies
*/
int TC_getNThbReac() { return ( TC_nThbReac_ ) ; }

/**
 * \brief Returns number of reactions with reaction rates in PLOG format
 */
int TC_getNplogReac() { return ( TC_nPlogReac_ ) ; }


/**
 * \brief Returns 1 if a reaction is reversible, 0 otherwise
*/
int TC_getReacRev(int ireac) { return (TC_isRev_[ireac]) ; }

/**
 * \brief Returns 1 if a reaction employs third-body efficiencies, 0 otherwise
*/
int TC_getReacThb(int ireac) { 

  int i, ithb ;

  if ( TC_nThbReac_ <= 0 ) return (0);

  ithb = 0;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacTbdy_[i] == ireac) {
      ithb = 1; /* found match */
      break ;
    }
  }
  return (ithb) ;

}

/**
 * \brief Returns 1 if a reaction is of "LOW" type, 0 otherwise
*/
int TC_getReacLow(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, ipfal;

  if ( TC_nFallReac_ <= 0 ) return (0);

  ipfal = 0;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacPfal_[i] == ireac) {
      if (TC_reacPlohi_[i] == 0) ipfal = 1; /* found match, check if LOW */
      break ;
    }
  }
  return (ipfal) ;

}

/**
 * \brief Returns 1 if a reaction is of "HIGH" type, 0 otherwise
*/
int TC_getReacHIGH(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, ipfal;

  if ( TC_nFallReac_ <= 0 ) return (0);

  ipfal = 0;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacPfal_[i] == ireac) {
      if (TC_reacPlohi_[i] == 1) ipfal = 1; /* found match, check if HIGH */
      break ;
    }
  }

  return (ipfal) ;

}

/**
 * \brief Returns type for pressure dependent reaction
   (1) Lindeman/(2) SRI/(3) Troe3/(4) Troe4/(6) Chebyshev, or 0 otherwise
*/
int TC_getReacPtype(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, iptype;

  if ( TC_nFallReac_ <= 0 ) return (0);

  iptype = 0;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacPfal_[i] == ireac) {
      iptype = TC_reacPtype_[i]; /* found match */
      break ;
    }
  }
  return (iptype) ;

}

/**
 * \brief Returns 1 if a reaction is of "PLOG" type, 0 otherwise
 */
int TC_getReacPLOG(int ireac) {
  /**
   \param ireac : reaction index
   */
  
  int i, iplog;
  
  if ( TC_nPlogReac_ <= 0 ) return (0);
  
  iplog = 0;
  for ( i=0; i<TC_nPlogReac_; i++) {
    if (TC_reacPlogIdx_[i] == ireac) {
      iplog = 1; /* found match */
      break ;
    }
  }
  
  return (iplog) ;
  
}

/**
 * \brief Returns the number of third-body efficiencies for a certain
 * reaction. Also returns the list of species indices and the
 * corresponding third-body efficiencies
*/
int TC_getThbReacEff(int ireac, int *specidx, double *effs) { 
/**
   \param ireac    : reaction index
   \return specidx : array with species indices 
   \return effs    : array with species efficiencies 
*/

  int i, j, nthb;

  if ( TC_nThbReac_ <= 0 ) return (0);

  nthb = 0;
  for ( i=0; i<TC_nThbReac_; i++) {
    if (TC_reacTbdy_[i] == ireac) {
      nthb = TC_reacTbno_[i]; /* found match */
      int indx = i*TC_maxTbInReac_;
      for ( j=0 ; j < nthb ; j++ ) {
        specidx[j] = TC_specTbdIdx_[indx+j];
        effs[j]    = TC_specTbdEff_[indx+j];
      }
      break ;
    }
  }
  return (nthb) ;

}

/**
 * \brief Returns the index of reaction in the list of reactions with third-body efficiencies
*/
int TC_getReacThbIdx(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, ithb;

  if ( TC_nThbReac_ <= 0 ) return (-1);

  ithb = -1;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacTbdy_[i] == ireac) {
      ithb = i; /* found match */
      break ;
    }
  }
  return (ithb) ;

}

/**
 * \brief Returns the index of reaction in the list of fall-off reactions
*/
int TC_getReacFallIdx(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, ipfal;

  if ( TC_nFallReac_ <= 0 ) return (-1);

  ipfal = -1;
  for ( i=0; i<TC_nFallReac_; i++) {
    if (TC_reacPfal_[i] == ireac) {
      ipfal = i;
      break ;
    }
  }
  return (ipfal) ;

}

/**
 * \brief Returns the index of reaction in the list of reactions with real stoichiometric coeffs
*/
int TC_getReacRealStIdx(int ireac) { 
/**
   \param ireac : reaction index
*/

  int i, irst;

  if (  TC_nRealNuReac_ <= 0 ) return (-1);

  irst = -1;
  for ( i=0; i< TC_nRealNuReac_; i++) {
    if (TC_reacRnu_[i] == ireac) {
      irst = i;
      break ;
    }
  }
  return (irst) ;

}

/**
 * \brief Returns stoichiometric coefficients' matrix. 
 * The stoichiometric coefficient for species "j" in reaction "i"
 * is stored at position \f$i\cdot N_{spec}+j\f$. It assumes that 
 * stoicoef was dimentioned to at least \f$N_{reac}\cdot N_{spec}\f$
 */
int TC_getStoiCoef( int Nspec, int Nreac, double *stoicoef )
{
  /**
   \param Nspec : no. of species 
   \param Nreac : no. of reactions 
   \return stoicoef : array of stoichiometric coefficients
  */

  int i, j, indx, indxR, kspec ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getStoiCoef", Nspec, TC_Nspec_ ) ;
  if ( Nreac != TC_Nreac_ )
    TC_errorMSG( 30, "TC_getStoiCoef", Nreac, TC_Nreac_ ) ;

  for ( i = 0 ; i < TC_Nspec_*TC_Nreac_ ; i++ ) stoicoef[i] = 0.0 ;

  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[i*TC_Nspec_+kspec] += (double) TC_reacNuki_[indx] ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[i*TC_Nspec_+kspec] += (double) TC_reacNuki_[indx] ;
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
	stoicoef[i*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
	stoicoef[i*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  return ( 0 ) ;

} /* done with TC_getStoiCoef */

/**
 * \brief Returns stoichiometric coefficients' matrix, saving forward
 * and reverse components separately. It assumes that 
 * stoicoef was dimentioned to at least \f$2*N_{reac}\cdot N_{spec}\f$
 */
int TC_getStoiCoefFR( int Nspec, int Nreac, double *stoicoef )
{
  /**
   \param Nspec : no. of species 
   \param Nreac : no. of reactions 
   \return stoicoef : 1D array of stoichiometric coefficients, stored
                      in a column major format
  */

  int i, j, indx, indxR, kspec ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getStoiCoef", Nspec, TC_Nspec_ ) ;
  if ( Nreac != TC_Nreac_ )
    TC_errorMSG( 30, "TC_getStoiCoef", Nreac, TC_Nreac_ ) ;

  for ( i = 0 ; i < TC_Nspec_*TC_Nreac_*2 ; i++ ) stoicoef[i] = 0.0 ;

  /* assemble matrix based on integer stoichiometric coefficients */
  for ( i = 0 ; i<TC_Nreac_; i++)
  {
    for ( j = 0; j<TC_reacNreac_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[i*TC_Nspec_+kspec] += (double) TC_reacNuki_[indx] ;
    }

    for ( j = 0; j<TC_reacNprod_[i] ; j++)
    {
      indx  = i*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[(i+TC_Nreac_)*TC_Nspec_+kspec] += (double) TC_reacNuki_[indx] ;
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
	stoicoef[i*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }

      indx  = i *TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
	kspec = TC_reacSidx_[indx] ;
	stoicoef[(i+TC_Nreac_)*TC_Nspec_+kspec] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  return ( 0 ) ;

} /* done with TC_getStoiCoefFR */

/**
 * \brief Returns stoichiometric coefficients' array for reaction 'ireac' 
 * for either reactants (idx=0) or products (idx=1) 
 * The stoichiometric coefficient for species "j" in reaction "ireac"
 * is stored at position \f$j\f$. It assumes that 
 * stoicoef was dimentioned to at least \f$N_{spec}\f$
 */
int TC_getStoiCoefReac( int Nspec, int Nreac, int ireac, int idx, double *stoicoef )
{
/**
   \param Nspec : no. of species 
   \param Nreac : no. of reactions 
   \param ireac : reaction index
   \param idx : 0-reactants, 1-products
   \return stoicoef : array of stoichiometric coefficients
*/

  int i, j, indx, indxR, kspec ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getStoiCoef", Nspec, TC_Nspec_ ) ;
  if ( Nreac != TC_Nreac_ )
    TC_errorMSG( 30, "TC_getStoiCoef", Nreac, TC_Nreac_ ) ;
  if ( (ireac < 0) || (ireac >= Nreac) )
  {
    printf("TC_getStoiCoefReac() : Error, illegal reaction number : %d -> terminate\n",ireac) ;
    exit(-1);
  }
  if ( (idx != 0) && (idx != 1) )
  {
    printf("TC_getStoiCoefReac() : Error, illegal idx value : %d -> terminate\n",idx) ;
    exit(-1);
  }

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) stoicoef[i] = 0.0 ;

  /* assemble matrix based on integer stoichiometric coefficients */
  if ( idx == 0 )
  {
    for ( j = 0; j<TC_reacNreac_[ireac] ; j++)
    {
      indx  = ireac*TC_maxSpecInReac_+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[kspec] += (double) TC_reacNuki_[indx] ;
    }
  }
  else 
  {
    for ( j = 0; j<TC_reacNprod_[ireac] ; j++)
    {
      indx  = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
      kspec = TC_reacSidx_[indx] ;
      stoicoef[kspec] += (double) TC_reacNuki_[indx] ;
    }

  }

  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;
    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      /* skip if this is not the reaction I was looking for */
      if ( TC_reacRnu_[ir] != ireac ) continue ;

      if ( idx == 0 )
      {
	indx  = TC_reacRnu_[ir] * TC_maxSpecInReac_ ;
	indxR = ir              * TC_maxSpecInReac_ ;
	for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
	{
	  kspec = TC_reacSidx_[indx] ;
	  stoicoef[kspec] += TC_reacRealNuki_[indxR] ;
	  indx++  ;
	  indxR++ ;
	}
      }
      else
      {
	indx  = TC_reacRnu_[ir] * TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
	indxR = ir              * TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
	for ( j = 0; j<TC_reacNprod_[i] ; j++)
	{
	  kspec = TC_reacSidx_[indx] ;
	  stoicoef[kspec] += TC_reacRealNuki_[indxR] ;
	  indx++  ;
	  indxR++ ;
	}
      } /* done if idx = 0,1 */

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  return ( 0 ) ;

} /* done with TC_getStoiCoefReac */

/**
 * \brief Returns stoichiometric coefficients' array for species
 * 'kspec' in reaction 'ireac' 
 * The stoichiometric coefficient for the reactant side is stored in
 * stoicoef[0], and for product side in stoicoef[1]. It assumes that 
 * stoicoef was dimentioned to at least \f$2\f$
 */
int TC_getStoiCoefIJ( int Nspec, int Nreac, int ireac, int kspec, double *stoicoef )
{
/**
   \param Nspec : no. of species 
   \param Nreac : no. of reactions 
   \param ireac : reaction index
   \param kspec : species index
   \return stoicoef : array of stoichiometric coefficients
*/

  int i, j, indx, indxR ;

  if ( Nspec != TC_Nspec_ )
    TC_errorMSG( 10, "TC_getStoiCoefIJ", Nspec, TC_Nspec_ ) ;
  if ( Nreac != TC_Nreac_ )
    TC_errorMSG( 30, "TC_getStoiCoefIJ", Nreac, TC_Nreac_ ) ;
  if ( (ireac < 0) || (ireac >= Nreac) )
  {
    printf("TC_getStoiCoefIJ() : Error, illegal reaction index : %d -> terminate\n",ireac) ;
    exit(-1);
  }
  if ( (kspec < 0) || (kspec >= Nspec) )
  {
    printf("TC_getStoiCoefIJ() : Error, illegal specie index : %d -> terminate\n",kspec) ;
    exit(-1);
  }

  for ( i = 0 ; i < 2 ; i++ ) stoicoef[i] = 0.0 ;

  /* assemble matrix based on integer stoichiometric coefficients */
  for ( j = 0; j<TC_reacNreac_[ireac] ; j++)
  {
    indx  = ireac*TC_maxSpecInReac_+j ;
    if ( kspec == TC_reacSidx_[indx]) 
      stoicoef[0] += (double) TC_reacNuki_[indx] ;
  }
  for ( j = 0; j<TC_reacNprod_[ireac] ; j++)
  {
    indx  = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j ;
    if ( kspec == TC_reacSidx_[indx] )
      stoicoef[1] += (double) TC_reacNuki_[indx] ;
  }

  /* add contributions from real stoichiometric coefficients if any */
  if ( TC_nRealNuReac_ > 0 )
  {
    int ir ;
    for ( ir = 0 ; ir<TC_nRealNuReac_; ir++)
    {

      /* skip if this is not the reaction I was looking for */
      if ( TC_reacRnu_[ir] != ireac ) continue ;

      indx  = TC_reacRnu_[ir] * TC_maxSpecInReac_ ;
      indxR = ir              * TC_maxSpecInReac_ ;
      for ( j = 0 ; j<TC_reacNreac_[i] ; j++)
      {
        if ( kspec == TC_reacSidx_[indx] )
	  stoicoef[0] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }
      indx  = TC_reacRnu_[ir] * TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      indxR = ir              * TC_maxSpecInReac_+TC_maxSpecInReac_/2 ;
      for ( j = 0; j<TC_reacNprod_[i] ; j++)
      {
        if ( kspec == TC_reacSidx_[indx] )
	  stoicoef[1] += TC_reacRealNuki_[indxR] ;
	indx++  ;
	indxR++ ;
      }

    } /* done loop over the number of reactions */

  } /* done if TC_nRealNuReac_ > 0 */

  return ( 0 ) ;

} /* done with TC_getStoiCoefIJ */

/**
 * \brief Returns reaction string
*/
int TC_getReacStr(int ireac, int lenstr, char *rstr) {

  int ir, j, indx, kspec, ifall ;
  char tmpstr[100];

  memset(rstr, 0, lenstr) ; 
  memset(tmpstr, 0, 100) ; 

  ir    = TC_getReacRealStIdx(ireac);
  ifall = TC_getReacFallIdx(ireac);

  /* Output reactant string */
  for ( j = 0 ; j<TC_reacNreac_[ireac] ; j++) {
    if ( ir >= 0 ) {
      indx = ir*TC_maxSpecInReac_+j;
      sprintf(tmpstr,"%f",-TC_reacRealNuki_[indx]) ; 
    }
    else {
      indx = ireac*TC_maxSpecInReac_+j;
      if ( -TC_reacNuki_[indx] != 1) {
        sprintf(tmpstr,"%d",-TC_reacNuki_[indx]);
	strcat(rstr, tmpstr) ;
        memset(tmpstr, 0, 100) ; 
      }
    }
    indx = ireac*TC_maxSpecInReac_+j;
    kspec = TC_reacSidx_[indx] ;
    sprintf(tmpstr,"%s",&TC_sNames_[kspec*LENGTHOFSPECNAME]);
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
    if (j != TC_reacNreac_[ireac]-1) 
      sprintf(tmpstr,"%s","+") ;  
      strcat(rstr, tmpstr) ;
      memset(tmpstr, 0, 100) ; 
  }

  /* Check for third-body */
  if ( (ifall<0) && (TC_getReacThb(ireac)==1)) {
    sprintf(tmpstr,"%s","+M") ;
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
  }

  /* Check for fall-off */
  if ( ifall>=0 ) {
    if (TC_reacPspec_[ifall]>=0)
      sprintf(tmpstr,"%s%s%s","(+",&TC_sNames_[TC_reacPspec_[ifall]*LENGTHOFSPECNAME],")");
    else 
      sprintf(tmpstr,"(+M)");
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
  }

  /* Output delimiter */
  if ( TC_isRev_[ireac])
    sprintf(tmpstr,"%s","<=>") ;
  else
    sprintf(tmpstr,"%s","=>") ;
  strcat(rstr, tmpstr) ;
  memset(tmpstr, 0, 100) ; 

  /* Output product string */
  for ( j = 0 ; j<TC_reacNprod_[ireac] ; j++) {
    if ( ir >= 0 ) {
      indx = ir*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j;
      sprintf(tmpstr,"%f",-TC_reacRealNuki_[indx]) ; 
    }
    else {
      indx = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j;
      if ( TC_reacNuki_[indx] != 1) {
        sprintf(tmpstr,"%d",TC_reacNuki_[indx]);
	strcat(rstr, tmpstr) ;
        memset(tmpstr, 0, 100) ; 
      }
    }
    indx = ireac*TC_maxSpecInReac_+TC_maxSpecInReac_/2+j;
    kspec = TC_reacSidx_[indx] ;
    sprintf(tmpstr,"%s",&TC_sNames_[kspec*LENGTHOFSPECNAME]);
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
    if (j != TC_reacNprod_[ireac]-1) 
      sprintf(tmpstr,"%s","+") ;  
      strcat(rstr, tmpstr) ;
      memset(tmpstr, 0, 100) ;
 
  } /* Done loop over products*/

  /* Check for third-body */
  if ( (ifall<0) && (TC_getReacThb(ireac)==1)) {
    sprintf(tmpstr,"%s","+M") ;
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
  }

  /* Check for fall-off */
  if ( ifall>=0 ) {
    if (TC_reacPspec_[ifall]>=0)
      sprintf(tmpstr,"%s%s%s","(+",&TC_sNames_[TC_reacPspec_[ifall]*LENGTHOFSPECNAME],")");
    else 
      sprintf(tmpstr,"(+M)");
    strcat(rstr, tmpstr) ;
    memset(tmpstr, 0, 100) ; 
  }

  return ( 0 ) ;

}
