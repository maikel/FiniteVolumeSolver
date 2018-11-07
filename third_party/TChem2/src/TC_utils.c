#include <math.h>

/*! \file TC_utils.c
 *  \brief Various utilities used by other functions
 */ 
/** 
 * \brief Much faster version of pow for small integers 
 *  (considering that reactions are generally no more than 3rd order).
 *  The C version of pow requires a double exponent; pow is the slowest part of the 
 *  code when used instead of fastIntPow. 
 */
double fastIntPow(double val, int exponent)
{
  switch(exponent)
  {
    case 0:
      return ( 1.0 );
    case 1:
      return ( val );
    case 2:
      return ( val*val );
    case 3:
      return ( val*val*val );
    case 4:
      return ( val*val*val*val );
    case 5:
      return ( val*val*val*val*val );
    case -1:
      return ( 1.0/val );
    case -2:
      return ( 1.0/(val*val) );
    case -3:
      return ( 1.0/(val*val*val) );
    default:
      return ( pow(val,exponent) );
  }
}

