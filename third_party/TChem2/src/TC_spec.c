#include "TC_defs.h"

#include <stdlib.h>
#include <string.h>

/*! \file TC_spec.c
    \brief Species info
*/ 
/*
	            _   _   _      
	  __ _  ___| |_| \ | |_/\__
	 / _` |/ _ \ __|  \| \    /
	| (_| |  __/ |_| |\  /_  _\
	 \__, |\___|\__|_| \_| \/  
	 |___/                     
*/
/**
  \brief Returns no. of species N<sub>spec</sub>
*/
int TC_getNspec() { return ( TC_Nspec_ ) ; }

/**
  \brief Returns no. of elements N<sub>elem</sub>
*/
int TC_getNelem() { return ( TC_Nelem_ ) ; }

/**
  \brief Returns no. of variables (N<sub>spec</sub>+1)
*/
int TC_getNvars() { return ( TC_Nvars_ ) ; }

/*
                 _   ____                                  
       __ _  ___| |_/ ___| _ __   __ _ _ __ ___   ___  ___ 
      / _` |/ _ \ __\___ \| '_ \ / _` | '_ ` _ \ / _ \/ __|
     | (_| |  __/ |_ ___) | | | | (_| | | | | | |  __/\__ \
      \__, |\___|\__|____/|_| |_|\__,_|_| |_| |_|\___||___/
      |___/                                                
*/
/**
  \brief Returns species names 
*/
int TC_getSnames(int Nspec,char *snames)
{
/**
  \param Nspec : no. of species 
  \return snames: array of characters containing species names, each name 
                  is LENGTHOFSPECNAME characters
*/

  if ( Nspec != TC_Nspec_ )
  {
    printf("Error : TC_getSnames - disagreement in the number of species : %d vs %d\n",
	   Nspec, TC_Nspec_ ) ;
    exit( -1 ) ;
  }

  memcpy( snames, TC_sNames_, LENGTHOFSPECNAME * TC_Nspec_ *sizeof(char)) ;
  
  return ( 0 ) ;
  
}

/**
  \brief Returns element names 
*/
int TC_getEnames(int Nelem,char *enames)
{
/**
  \param Nspec : no. of species 
  \return snames: array of characters containing species names, each name 
                  is LENGTHOFSPECNAME characters
*/

  if ( Nelem != TC_Nelem_ )
  {
    printf("Error : TC_getEnames - disagreement in the number of elements : %d vs %d\n",
	   Nelem, TC_Nelem_ ) ;
    exit( -1 ) ;
  }

  memcpy( enames, TC_eNames_, LENGTHOFELEMNAME * TC_Nelem_ *sizeof(char)) ;
  
  return ( 0 ) ;
  
}

/**
  \brief Returns length of species names 
*/
int TC_getSnameLen()
{
  return ( LENGTHOFSPECNAME ) ;
}
/**
  \brief Returns length of element names 
*/
int TC_getEnameLen()
{
  return ( LENGTHOFELEMNAME ) ;
}
/*
                    _   ____                  
          __ _  ___| |_/ ___| _ __   ___  ___ 
         / _` |/ _ \ __\___ \| '_ \ / _ \/ __|
        | (_| |  __/ |_ ___) | |_) | (_) \__ \
         \__, |\___|\__|____/| .__/ \___/|___/
         |___/               |_|              

*/
/**
  \brief Returns position a species in the list of species 
*/
int TC_getSpos(const char *sname, const int slen)
{
/**
  \param sname : string containing the name of the species
  \param slen : length of species "sname" name 
  \return position of species sname in the list of species, 0...(N<sub>spec</sub>-1)
*/
  int i ;

  for ( i = 0 ; i < TC_Nspec_ ; i ++ )
  {
    /* printf("%s  %s\n",&TC_sNames_[i*LENGTHOFSPECNAME],sname) ; */
    if ( strncmp(&TC_sNames_[i*LENGTHOFSPECNAME],sname,LENGTHOFSPECNAME) == 0 ) return ( i ) ;
  }
  return ( -1 ) ;
  
}

/*
              _   ____                          
    __ _  ___| |_/ ___| _ __ ___   __ _ ___ ___ 
   / _` |/ _ \ __\___ \| '_ ` _ \ / _` / __/ __|
  | (_| |  __/ |_ ___) | | | | | | (_| \__ \__ \
   \__, |\___|\__|____/|_| |_| |_|\__,_|___/___/
   |___/                                        

*/
/**
  \brief Returns species molar weights. 
*/
int TC_getSmass(int Nspec,double *Wi)
{
/**
  \param Nspec : no. of species 
  \return Wi : array of species molar weights [kg/kmol]=[g/mol]
*/

  int i ;

  if ( Nspec != TC_Nspec_ )
  {
    printf("Error in TC_getSmass() : disagreement in Nspec : %d vs %d\n",
	   Nspec, TC_Nspec_ ) ;
    exit( -1 ) ;
  }

  for ( i = 0 ; i < TC_Nspec_ ; i++ ) Wi[i] = TC_sMass_[i] ;

  if ( TC_nonDim_ == 1 ) 
    for ( i = 0 ; i < TC_Nspec_ ; i++ ) Wi[i] /= TC_Wref_ ;

  return ( 0 ) ;
  
}

/**
  \brief Returns list of element counts 
*/
int TC_getECount(int *elemcnt)
{
/**
  \return elemcnt : no. of atoms of element j in each species i at (i*TC_Nelem_+j)*
*/

  int i ;
  for ( i = 0 ; i < TC_Nspec_*TC_Nelem_ ; i++ ) elemcnt[i] = TC_elemcount_[i];

  return ( 0 ) ;
  
}
