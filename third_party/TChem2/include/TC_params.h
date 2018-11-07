#ifndef TCParamsHSeen
#define TCParamsHSeen

/*! 
 * \file TC_params.h
 * \brief Definitions of parameters and constants
*/ 
#include "copyright.h"

/* -------------------------------------------------------------------------
                           Macros
   ------------------------------------------------------------------------- */
/** 
 * \def MAX
 * Maximum of two expressions 
 */
#define MAX(A,B) ( ((A) > (B)) ? (A) : (B) )
/** 
 * \def MIN
 * Minimum of two expressions 
 */
#define MIN(A,B) ( ((A) < (B)) ? (A) : (B) )


/* -------------------------------------------------------------------------
                           Parameter definitions
   ------------------------------------------------------------------------- */
/** 
 * \def LENGTHOFELEMNAME
 * Maximum number of characters for element names 
 */
#define LENGTHOFELEMNAME 3

/** 
 * \def LENGTHOFSPECNAME
 *  Maximum number of characters for species names
 */
#define LENGTHOFSPECNAME 18

/**
 * \def NUMBEROFELEMINSPEC
 * Maximum number of (different) elements that compose a species
 */
#define NUMBEROFELEMINSPEC 5

/**
 * \def NTHRDBMAX
 * Maximum number of third body efficiencies 
 */
#define NTHRDBMAX 10

/**
 * \def NPLOGMAX
 * Maximum number of interpolation ranges for PLOG 
 */
#define NPLOGMAX 10

/**
 * \def NTH9RNGMAX
 * Maximum number of temperature ranges for 9-coefficients NASA polynomials
 */
#define NTH9RNGMAX 5

/**
 * \def NSPECREACMAX
 * Maximum number of reactant or product species in a reaction 
 */
#define NSPECREACMAX 6

/** 
 * \def REACBALANCE
 * Treshold for checking reaction balance with 
 *        real stoichiometric coefficients
*/
#define REACBALANCE 1.e-4

/** 
 * \def RUNIV
 * Universal gas constant \f$J/(mol\cdot K)\f$
*/
#define RUNIV 8.314472

/** 
 * \def NAVOG 
 * Avogadro's number
*/
#define NAVOG 6.02214179E23

/** 
 * \def ATMPA
 * Standard atmospheric pressure \f$[Pa]\f$
 */
#define ATMPA 101325.0

/** 
 * \def CALJO
 * Conversion from calories to Joule
 */
#define CALJO 4.184

/** 
 * \def KBOLT
 * Boltzmann's constant \f$(k_B)\f$ \f$[JK^{-1}]\f$
 */
#define KBOLT 1.3806504E-23

/** 
 * \def EVOLT
 * electron volt (eV) unit \f$[J]\f$
 */
#define EVOLT 1.60217653E-19

#endif
