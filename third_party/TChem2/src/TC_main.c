#include "../include/TC_interface.h"
#include "../include/TC_defs.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "assert.h"

/*! \mainpage Thermo-chemical software library
 * \authors Cosmin Safta, Habib Najm, Omar Knio
 
The TChem toolkit is a software library that enables numerical simulations 
using complex chemistry and facilitates the analysis of detailed kinetic
models. The toolkit provide capabilities for thermodynamic properties 
based on NASA polynomials, kinetic model reaction rates (both for individual
reactions and for species). It incorporates methods that can selectively 
modify reaction parameters for sensitivity analysis and uncertainty
quantification. The library contains several functions that provide 
analytically computed Jacobian matrices necessary for the efficient 
time advancement and analysis of detailed kinetic models. 

The software is released subject the terms of the 
GNU Lesser General Public License 
(http://www.gnu.org/copyleft/lesser.html)
Please read the license information before downloading the software. 

Copyright 2011 Sandia Corporation. 
Under the terms of Contract DE-AC04-94AL85000 with 
Sandia Corporation, the U.S. Government retains certain rights 
in this software.
*/
/**
 * \defgroup eqstate Equation of state
 * \defgroup thermo  Thermodynamic properties
 * \defgroup init    Initialization
 * \defgroup jacs    Jacobians
 * \defgroup srcs    Source terms
 * \defgroup maxpar  Max. no. of parameters
 * \defgroup noreac  No. of reactions
 * \defgroup nospec  No. of species, variables, etc.
 */
#include "TC_chg.c"
#include "TC_edit.c"
#include "TC_init.c"
#include "TC_mlms.c"
#include "TC_rr.c"
#include "TC_spec.c"
#include "TC_src.c"
#include "TC_thermo.c"
#include "TC_utils.c"
#include "TC_info.c"
#include "TC_jac.c"

#include "TC_for.c"

