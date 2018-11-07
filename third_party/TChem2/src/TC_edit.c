#include "TC_defs.h"

#include <stdlib.h>
#include <string.h>


/** 
  \brief Frees all memory and sets variables to 0 so that TC_initChem can be called
         again without a memory leak.  Not designed for use with tables.   
*/
void TC_reset()
{
  /* For straightforward error checking, arrays are freed in the same order as 
  allocation in TC_init.c */

  /* From TC_initChem() */
  free(TC_eNames_);       TC_eNames_       = 0;
  free(TC_eMass_);        TC_eMass_        = 0;
  free(TC_sNames_);       TC_sNames_       = 0;
  free(TC_sMass_);        TC_sMass_        = 0;
  free(TC_elemcount_);    TC_elemcount_    = 0;
  free(TC_sCharge_);      TC_sCharge_      = 0;
  free(TC_sTfit_);        TC_sTfit_        = 0;
  free(TC_sPhase_);       TC_sPhase_       = 0;
  free(TC_Tlo_);          TC_Tlo_          = 0;
  free(TC_Tmi_);          TC_Tmi_          = 0;
  free(TC_Thi_);          TC_Thi_          = 0;
  free(TC_cppol_);        TC_cppol_        = 0;
  free(TC_sNion_);        TC_sNion_        = 0;
  free(TC_isRev_);        TC_isRev_        = 0;
  free(TC_reacNrp_);      TC_reacNrp_      = 0;
  free(TC_reacNreac_);    TC_reacNreac_    = 0;
  free(TC_reacNprod_);    TC_reacNprod_    = 0;
  free(TC_reacNuki_);     TC_reacNuki_     = 0;
  free(TC_reacNukiDbl_);  TC_reacNukiDbl_  = 0;
  free(TC_reacSidx_);     TC_reacSidx_     = 0;
  free(TC_reacScoef_);    TC_reacScoef_    = 0;
  free(TC_reacArhenFor_); TC_reacArhenFor_ = 0;
  free(TC_reacRev_);      TC_reacRev_      = 0;
  free(TC_reacArhenRev_); TC_reacArhenRev_ = 0;
  free(TC_reacPfal_);     TC_reacPfal_     = 0;
  free(TC_reacPtype_);    TC_reacPtype_    = 0;
  free(TC_reacPlohi_);    TC_reacPlohi_    = 0;
  free(TC_reacPspec_);    TC_reacPspec_    = 0;
  free(TC_reacPpar_);     TC_reacPpar_     = 0;
  free(TC_reacTbdy_);     TC_reacTbdy_     = 0;
  free(TC_reacTbno_);     TC_reacTbno_     = 0;
  free(TC_specTbdIdx_);   TC_specTbdIdx_   = 0;
  free(TC_specTbdEff_);   TC_specTbdEff_   = 0;
  free(TC_reacRnu_);      TC_reacRnu_      = 0;
  free(TC_reacRealNuki_); TC_reacRealNuki_ = 0;
  free(TC_reacAOrd_);     TC_reacAOrd_     = 0;
  free(TC_specAOidx_);    TC_specAOidx_    = 0;
  free(TC_specAOval_);    TC_specAOval_    = 0;
  free(TC_kc_coeff);      TC_kc_coeff      = 0;

  /* From TC_makeSpace() */
  free(TC_cpks);          TC_cpks    = 0;
  free(TC_hks);           TC_hks     = 0;
  free(TC_Xconc);         TC_Xconc   = 0;
  free(TC_gk);            TC_gk      = 0;
  free(TC_gkp);           TC_gkp     = 0;
  free(TC_omg);           TC_omg     = 0;
  free(TC_omgP);          TC_omgP    = 0;
  free(TC_CrndDer);       TC_CrndDer = 0;
  free(TC_PrDer);         TC_PrDer   = 0;
  free(TC_dFfac);         TC_dFfac   = 0;
  free(TC_scalIn);        TC_scalIn  = 0;
  free(TC_jacFull);       TC_jacFull = 0;

  free(TC_Mconc);         TC_Mconc  = 0;
  free(TC_kfor);          TC_kfor   = 0;
  free(TC_krev);          TC_krev   = 0;
  free(TC_kforP);         TC_kforP  = 0;
  free(TC_krevP);         TC_krevP  = 0;
  free(TC_ropFor);        TC_ropFor = 0;
  free(TC_ropRev);        TC_ropRev = 0;
  free(TC_rop);           TC_rop    = 0;
  free(TC_Crnd);          TC_Crnd   = 0;

  free(TC_sigNu);         TC_sigNu     = 0;
  free(TC_NuIJ);          TC_NuIJ      = 0;
  free(TC_sigRealNu);     TC_sigRealNu = 0;
  free(TC_RealNuIJ);      TC_RealNuIJ  = 0;

  free(qfr);              qfr = 0;

  /* From TC_createTables() */
  free(TC_cptab);         TC_cptab  = 0;
  free(TC_cpPtab);        TC_cpPtab = 0;
  free(TC_htab);          TC_htab   = 0;
  free(TC_gktab);         TC_gktab  = 0;
  free(TC_gkPtab);        TC_gkPtab = 0;

  free(TC_kfortab);       TC_kfortab  = 0;
  free(TC_krevtab);       TC_krevtab  = 0;
  free(TC_kforPtab);      TC_kforPtab = 0;
  free(TC_krevPtab);      TC_krevPtab = 0;

  /* Now free and reset Arrhenius backup parameters from TC_chg.c */
  free(TC_reacArhenRevSave_); TC_reacArhenRevSave_ = 0;
  free(TC_reacArhenForSave_); TC_reacArhenForSave_ = 0;
  TC_ArhenRevChg_ = 0;
  TC_ArhenForChg_ = 0;

  /* Now free reaction backup data */
  free(TC_isRevBackup_);         TC_isRevBackup_        = 0;
  free(TC_reacNrpBackup_);       TC_reacNrpBackup_      = 0;
  free(TC_reacNreacBackup_);     TC_reacNreacBackup_    = 0;
  free(TC_reacNprodBackup_);     TC_reacNprodBackup_    = 0;
  free(TC_reacArhenForBackup_);  TC_reacArhenForBackup_ = 0;
  free(TC_reacArhenRevBackup_);  TC_reacArhenRevBackup_ = 0;
  free(TC_reacNukiBackup_);      TC_reacNukiBackup_     = 0;
  free(TC_reacSidxBackup_);      TC_reacSidxBackup_     = 0;
  free(TC_reacScoefBackup_);     TC_reacScoefBackup_    = 0;
  free(TC_reacNukiDblBackup_);   TC_reacNukiDblBackup_  = 0;
  free(TC_reacPfalBackup_);      TC_reacPfalBackup_     = 0;
  free(TC_reacPtypeBackup_);     TC_reacPtypeBackup_    = 0;
  free(TC_reacPlohiBackup_);     TC_reacPlohiBackup_    = 0;
  free(TC_reacPspecBackup_);     TC_reacPspecBackup_    = 0;
  free(TC_reacTbdyBackup_);      TC_reacTbdyBackup_     = 0;
  free(TC_reacTbnoBackup_);      TC_reacTbnoBackup_     = 0;
  free(TC_specTbdIdxBackup_);    TC_specTbdIdxBackup_   = 0;
  free(TC_reacPparBackup_);      TC_reacPparBackup_     = 0;
  free(TC_specTbdEffBackup_);    TC_specTbdEffBackup_   = 0;
  free(TC_reacRnuBackup_);       TC_reacRnuBackup_      = 0;
  free(TC_reacAOrdBackup_);      TC_reacAOrdBackup_     = 0;
  free(TC_specAOidxBackup_);     TC_specAOidxBackup_    = 0;
  free(TC_reacRealNukiBackup_);  TC_reacRealNukiBackup_ = 0;
  free(TC_sigRealNuBackup);      TC_sigRealNuBackup     = 0;
  free(TC_RealNuIJBackup);       TC_RealNuIJBackup      = 0;
  free(TC_specAOvalBackup_);     TC_specAOvalBackup_    = 0;
  free(TC_kc_coeffBackup);       TC_kc_coeffBackup      = 0;
  free(TC_sigNuBackup);          TC_sigNuBackup         = 0;
  free(TC_NuIJBackup);           TC_NuIJBackup          = 0;

  /* Reset tchem */
  TC_isInit_     = 0 ;
  TC_tab_        = 0 ;
  TC_initRemoved = 0 ;

  return ;

} /* end of TC_reset() */


/** 
  \brief Removes a reaction from the mechanism.  Not designed for use with tables.

  \param reacArr        : 0-based reaction indices in ascending order.
  \param numRemoveReacs : length of rearArr array.
  \param revOnly        : set to 1 to remove reverse only, 0 to remove forward and 
                          reverse.
*/
void TC_removeReaction(int* reacArr, int numRemoveReacs, int revOnly)
{  
  int i, j, ireac, ireacLast, loopMin, finalIdx, destIdx, curIdx ;
  int ispec,inext; /* Indices for special reaction types: third body, etc */
  int numMoveItems;
  int numPfal, numthb;

  if( TC_tab_ == 1 ) /* Not designed for use with tables. */
  {
    fprintf(stderr,"Reaction removal with tables has not been implemented!\n");
    exit(1) ;
  }
  ireacLast = reacArr[0] - 1 ;
  for( i=0 ; i < numRemoveReacs ; ++i )
  {
    ireac = reacArr[i];
    if ( ( ireac < 0 ) || ( ireac >= TC_Nreac_ ) )
    {
      fprintf(stderr,"When removing reaction: index is wrong -> %d\n", ireac) ;
      exit(1) ;
    }
    else if( ireac <= ireacLast )
    {
      fprintf(stderr,
        "When removing reaction: repeated reaction or not in ascending order -> %d,%d\n", 
        ireacLast, ireac) ;
      exit(1) ;
    }
    ireacLast = ireac ;
  }

  /* Passed error checks, so initiate backup before removing reaction */
  if ( TC_initRemoved == 0)
  {
    TC_initRemoved = 1 ;

    TC_NreacBackup_ = TC_Nreac_ ;

    /* If reaction is reversible */
    TC_isRevBackup_ = (int*) malloc( TC_Nreac_ * sizeof(int) ) ;
    memcpy(TC_isRevBackup_, TC_isRev_, TC_Nreac_*sizeof(int) ) ;

    /* Number of reactants and products */
    TC_reacNrpBackup_ = (int*) malloc( TC_Nreac_ * sizeof(int) ) ;
    memcpy(TC_reacNrpBackup_, TC_reacNrp_, TC_Nreac_*sizeof(int) ) ;

    /* Number of reactants */
    TC_reacNreacBackup_ = (int*) malloc ( TC_Nreac_ * sizeof(int) ) ;
    memcpy(TC_reacNreacBackup_, TC_reacNreac_, TC_Nreac_*sizeof(int) ) ;

    /* Number of products */
    TC_reacNprodBackup_ = (int*) malloc ( TC_Nreac_ * sizeof(int) ) ;
    memcpy(TC_reacNprodBackup_, TC_reacNprod_, TC_Nreac_*sizeof(int) ) ;

    /* Stoichiometric coefficients */
    TC_reacNukiBackup_    = (int*)    malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(int) ) ;
    TC_reacNukiDblBackup_ = (double*) malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(double) ) ;
    TC_reacSidxBackup_    = (int*)    malloc( TC_Nreac_*TC_maxSpecInReac_* sizeof(int) ) ;
    memcpy(TC_reacNukiBackup_,    TC_reacNuki_,    TC_maxSpecInReac_*TC_Nreac_*sizeof(int) ) ;
    memcpy(TC_reacNukiDblBackup_, TC_reacNukiDbl_, TC_maxSpecInReac_*TC_Nreac_*sizeof(double) );
    memcpy(TC_reacSidxBackup_,    TC_reacSidx_,    TC_maxSpecInReac_*TC_Nreac_*sizeof(int) ) ;
    /* TC_reacScoef_ determines if stoichiometric coefficients are int or real */
    TC_reacScoefBackup_ = (int*) malloc( TC_Nreac_* sizeof(int) ) ;
    memcpy(TC_reacScoefBackup_, TC_reacScoef_, TC_Nreac_*sizeof(int) ) ;

    /* Arrhenius parameters */
    TC_reacArhenForBackup_ = (double*) malloc( TC_Nreac_*3*sizeof(double) ) ;
    if( TC_ArhenForChg_ == 1)
    {
      memcpy(TC_reacArhenForBackup_, TC_reacArhenForSave_, 3*TC_Nreac_*sizeof(double) );
    }
    else
    {
      memcpy(TC_reacArhenForBackup_, TC_reacArhenFor_, 3*TC_Nreac_*sizeof(double) );
    }

    if (TC_nRevReac_ > 0 ) /* No. reactions with reverse Arrhenius parameters */
    {
      TC_NrevBackup_ = TC_nRevReac_;

      TC_reacArhenRevBackup_ = (double*) malloc( TC_nRevReac_*3* sizeof(double) );
      if( TC_ArhenRevChg_ == 1)
      {
        memcpy(TC_reacArhenRevBackup_, TC_reacArhenRevSave_, 3*TC_nRevReac_*sizeof(double) );
      }
      else
      {
        memcpy(TC_reacArhenRevBackup_, TC_reacArhenRev_, 3*TC_nRevReac_*sizeof(double) ); 
      }
    }

    /* Pressure-dependent reactions */
    if ( TC_nFallReac_ > 0 )
    {
      TC_NfalBackup_ = TC_nFallReac_;

      TC_reacPfalBackup_  = (int*) malloc( TC_nFallReac_ * sizeof(int) ) ;
      TC_reacPtypeBackup_ = (int*) malloc( TC_nFallReac_ * sizeof(int) ) ;
      TC_reacPlohiBackup_ = (int*) malloc( TC_nFallReac_ * sizeof(int) ) ;
      TC_reacPspecBackup_ = (int*) malloc( TC_nFallReac_ * sizeof(int) ) ;
      TC_reacPparBackup_  = (double*) malloc( TC_nFallReac_*TC_nFallPar_*sizeof(double)) ;

      memcpy(TC_reacPfalBackup_,  TC_reacPfal_,  TC_nFallReac_*sizeof(int) );
      memcpy(TC_reacPtypeBackup_, TC_reacPtype_, TC_nFallReac_*sizeof(int) );
      memcpy(TC_reacPlohiBackup_, TC_reacPlohi_, TC_nFallReac_*sizeof(int) );
      memcpy(TC_reacPspecBackup_, TC_reacPspec_, TC_nFallReac_*sizeof(int) );
      memcpy(TC_reacPparBackup_,  TC_reacPpar_,  TC_nFallReac_*TC_nFallPar_*sizeof(double) );
    }

    /* Third-body reactions */
    if ( TC_nThbReac_ > 0 )
    {
      TC_NthbBackup_ = TC_nThbReac_ ;

      TC_reacTbdyBackup_   = (int*)    malloc( TC_nThbReac_ * sizeof(int) ) ;
      TC_reacTbnoBackup_   = (int*)    malloc( TC_nThbReac_ * sizeof(int) ) ;
      TC_specTbdIdxBackup_ = (int*)    malloc( TC_nThbReac_*TC_maxTbInReac_ * sizeof(int) ) ;
      TC_specTbdEffBackup_ = (double*) malloc( TC_nThbReac_*TC_maxTbInReac_ * sizeof(double) ) ;

      memcpy(TC_reacTbdyBackup_,   TC_reacTbdy_,   TC_nThbReac_ * sizeof(int) );
      memcpy(TC_reacTbnoBackup_,   TC_reacTbno_,   TC_nThbReac_ * sizeof(int) );
      memcpy(TC_specTbdIdxBackup_, TC_specTbdIdx_, TC_nThbReac_ * TC_maxTbInReac_ * sizeof(int) );
      memcpy(TC_specTbdEffBackup_, TC_specTbdEff_, TC_nThbReac_ * TC_maxTbInReac_ * sizeof(double) );
    }

    /* Reactions with real stoichiometric coefficients */
    if ( TC_nRealNuReac_ > 0 ) 
    {
      TC_nRealNuReacBackup_ =  TC_nRealNuReac_ ;

      TC_reacRnuBackup_      = (int*)    malloc( TC_nRealNuReac_ * sizeof(int) ) ;
      TC_reacRealNukiBackup_ = (double*) malloc( TC_nRealNuReac_ * TC_maxSpecInReac_ * sizeof(double) ) ;
      TC_sigRealNuBackup     = (double*) malloc( TC_nRealNuReac_ * sizeof(double) ) ;
      TC_RealNuIJBackup      = (double*) malloc( TC_nRealNuReac_ * TC_Nspec_ * sizeof(double) ) ;

      memcpy(TC_reacRnuBackup_,      TC_reacRnu_,      TC_nRealNuReac_ * sizeof(int) );
      memcpy(TC_reacRealNukiBackup_, TC_reacRealNuki_, TC_nRealNuReac_ * TC_maxSpecInReac_ * sizeof(double) );
      memcpy(TC_sigRealNuBackup,     TC_sigRealNu,     TC_nRealNuReac_ * sizeof(double) );
      memcpy(TC_RealNuIJBackup,      TC_RealNuIJ,      TC_nRealNuReac_ * TC_Nspec_ * sizeof(double) );
    }

    /* Arbitrary reaction orders */
    if ( TC_nOrdReac_ > 0 ) 
    {
      TC_nOrdReacBackup_ = TC_nOrdReac_ ;

      TC_reacAOrdBackup_  = (int*)    malloc( TC_nOrdReac_ * sizeof(int) ) ;
      TC_specAOidxBackup_ = (int*)    malloc( TC_nOrdReac_*TC_maxOrdPar_*sizeof(int)) ;
      TC_specAOvalBackup_ = (double*) malloc( TC_nOrdReac_*TC_maxOrdPar_*sizeof(double)) ;

      memcpy(TC_reacAOrdBackup_,  TC_reacAOrd_,  TC_nOrdReac_*sizeof(int) );
      memcpy(TC_specAOidxBackup_, TC_specAOidx_, TC_nOrdReac_*TC_maxOrdPar_*sizeof(int) );
      memcpy(TC_specAOvalBackup_, TC_specAOval_, TC_nOrdReac_*TC_maxOrdPar_*sizeof(double) );
    }

    /* sum(nu) and TC_kc_coeff for each reaction */
    TC_sigNuBackup    = (int*)    malloc( TC_Nreac_ * sizeof(int) ) ;
    TC_kc_coeffBackup = (double*) malloc( TC_Nreac_ * sizeof(double));
    memcpy(TC_sigNuBackup,    TC_sigNu,    TC_Nreac_ * sizeof(int) );
    memcpy(TC_kc_coeffBackup, TC_kc_coeff, TC_Nreac_ * sizeof(double) );

    /* NuIJ=NuII''-NuIJ' */
    TC_NuIJBackup  = (int*) malloc( TC_Nreac_ * TC_Nspec_ * sizeof(int) ) ; 
    memcpy(TC_NuIJBackup, TC_NuIJ, TC_Nreac_ * TC_Nspec_ * sizeof(int) );
  }

  /*------------------------------------------------------------------------
    Backup complete - now remove reactions
  ------------------------------------------------------------------------ */

  if( revOnly == 1) /* Remove reverse only */
  {
      for( j=0 ; j < numRemoveReacs ; ++j )
      {
        ireac = reacArr[j];
      if( TC_isRev_[ireac] == 1)
      {
        TC_isRev_[ireac] = 0;

        /* Check for reverse Arrhenius parameters and remove */
        ispec = -1;
        loopMin = (ireac+1 < TC_nRevReac_) ? ireac+1 : TC_nRevReac_;
        for ( i = 0 ; i < loopMin ; i++ )
        {
          if ( TC_reacRev_[i] == ireac ) { ispec = i; break ;}
        }
        if(ispec > -1)
        {
          memmove(&TC_reacRev_[ispec],&TC_reacRev_[ispec+1],3*(TC_nRevReac_-ispec-1)*sizeof(int));
          memmove(&TC_reacArhenRev_[ispec*3],&TC_reacArhenRev_[ispec*3+3],
            3*(TC_nRevReac_-ispec-1)*sizeof(double));

          /* TC_chg.c reverse array */
          if ( TC_ArhenRevChg_ == 1 )
          {
            memmove(&TC_reacArhenRevSave_[ispec*3],&TC_reacArhenRevSave_[ispec*3+3],
              3*(TC_nRevReac_-ispec-1)*sizeof(double));
          }

          --TC_nRevReac_;
        }
        /* else reaction is not reversible already */
      }
      return;
    }
  } /* End of remove reverse only */

  /* If we have reached this point, then remove the reactions completely.
  Remove one item from each reaction data array per reaction.
  Arrays are listed here in same order as declaration in TC_init.c */

  destIdx = reacArr[0] ; /* Overwrite this location on first iteration */
  finalIdx = destIdx ;
  for( j=0 ; j < numRemoveReacs ; ++j )
  {
    ireac = finalIdx + 1 ; /* Skip reaction finalIdx - ireac denotes first reaction to keep after destIdx */
    curIdx = reacArr[j] - j;  /* Reaction being removed at this step in new reaction indices (i.e. subtract off number of reactions already removed) */

    /* Remove reactions as far as one before finalIdx */
    if( j == numRemoveReacs - 1 )
      finalIdx = TC_Nreac_ ; 
    else
      finalIdx = reacArr[j+1] ;

    numMoveItems = finalIdx - ireac ;
    /* Still continue if numMoveItems is 0 - need to deal with corrections to indices of 
       special case reactions (e.g. third body, arbitrary order, etc. ) */

    /* Use memmove rather than memcpy due to overlapping memory regions */

    /* If reaction is reversible */
    memmove(&TC_isRev_[destIdx], &TC_isRev_[ireac], numMoveItems * sizeof(int) );

    /* Number of reactants and products */
    memmove(&TC_reacNrp_[destIdx], &TC_reacNrp_[ireac], numMoveItems * sizeof(int) );

    /* Number of reactants */
    memmove(&TC_reacNreac_[destIdx], &TC_reacNreac_[ireac], numMoveItems * sizeof(int) );

    /* Number of products */
    memmove(&TC_reacNprod_[destIdx], &TC_reacNprod_[ireac], numMoveItems * sizeof(int) );

    /* Stoichiometric coefficients */
    memmove(&TC_reacNuki_[destIdx*TC_maxSpecInReac_], &TC_reacNuki_[ireac*TC_maxSpecInReac_], 
      numMoveItems * TC_maxSpecInReac_ * sizeof(int) );
    memmove(&TC_reacNukiDbl_[destIdx*TC_maxSpecInReac_], &TC_reacNukiDbl_[ireac*TC_maxSpecInReac_], 
      numMoveItems * TC_maxSpecInReac_ * sizeof(double) );
    memmove(&TC_reacSidx_[destIdx*TC_maxSpecInReac_], &TC_reacSidx_[ireac*TC_maxSpecInReac_], 
      numMoveItems * TC_maxSpecInReac_ * sizeof(int) );
    /* TC_reacScoef_ determines if stoichiometric coefficients are int or real.
       Real values are edited again later, so TC_reacScoef_ is treated after that. */

    /* Arrhenius parameters */
    memmove(&TC_reacArhenFor_[3*destIdx], &TC_reacArhenFor_[3*ireac],
      3 * numMoveItems * sizeof(double) );
    

    if (TC_nRevReac_ > 0 ) /* No. reactions with reverse Arrhenius parameters */
    {
      /* Reverse Arrhenius parameter index, if applicable */
      ispec = -1 ;
      inext = -1 ;
      for ( i = 0 ; i < TC_nRevReac_ ; i++ )
      {
        if ( TC_reacRev_[i] == curIdx ) 
        { 
          ispec = i; 
          break ; 
        }
        else if( TC_reacRev_[i] > curIdx ) /* Moved too far */
        {
          inext = i ; /* Index of first reaction with index greater than curIdx */
          break;
        }
      }
      if( ispec > -1 )
      {
        /* Note: this may be somewhat inefficient, as it requires iterating over all
           remaining reversible reactions (even though some may be removed later), but 
           the alternative is to check to see if any of the other reactions to be removed 
           are also reversible. This is likely to be expensive. */

        /* Move reaction index array down by 1 and subtract 1 from each index */
        for ( i = ispec ; i < TC_nRevReac_-1 ; i++ )
        {
          TC_reacRev_[i] = TC_reacRev_[i+1] - 1 ;
        }
        memmove(&TC_reacArhenRev_[ispec*3], &TC_reacArhenRev_[ispec*3+3],
          3*(TC_nRevReac_-ispec-1)*sizeof(double) ) ;

        /* TC_chg.c reverse array */
        if ( TC_ArhenRevChg_ == 1 )
        {
          memmove(&TC_reacArhenRevSave_[ispec*3], &TC_reacArhenRevSave_[ispec*3+3],
            3*(TC_nRevReac_-ispec-1)*sizeof(double) ) ;
        }

        --TC_nRevReac_;
      }
      else if( inext > -1 )
      {
        /* Array remains the same size, but must subtract 1 from each index for
        reactions above the reaction that has been removed */
        for ( i = inext ; i < TC_nRevReac_ ; ++i )
        {
          TC_reacRev_[i] -= 1 ;
        }     
      }
    }

    /* Pressure-dependent reactions */
    if ( TC_nFallReac_ > 0 )
    {
      /* Pressure-dependent index, if applicable */
      ispec = -1;
      inext = -1;
      for ( i = 0 ; i < TC_nFallReac_; i++ )
      {
        if ( TC_reacPfal_[i] == curIdx )
        { 
          ispec = i ;
          break ;
        }
        else if( TC_reacPfal_[i] > curIdx ) /* Moved too far */
        {
          inext = i ; /* Index of first reaction with index greater than curIdx */
          break ;
        }
      }
      if( ispec > -1)
      {
        /* Move reaction number array down by 1 and subtract 1 from each value */
        for ( i = ispec ; i < TC_nFallReac_-1 ; i++ )
        {
          TC_reacPfal_[i] = TC_reacPfal_[i+1] - 1;
        }

        numPfal=TC_nFallReac_-ispec-1;
        memmove(&TC_reacPtype_[ispec],&TC_reacPtype_[ispec+1],
          numPfal*sizeof(int) );
        memmove(&TC_reacPlohi_[ispec],&TC_reacPlohi_[ispec+1],
          numPfal*sizeof(int) );
        memmove(&TC_reacPspec_[ispec],&TC_reacPspec_[ispec+1],
          numPfal*sizeof(int) );
        memmove(&TC_reacPpar_[ispec*TC_nFallPar_], &TC_reacPpar_[(ispec+1)*TC_nFallPar_],
          numPfal*TC_nFallPar_*sizeof(double) );

        --TC_nFallReac_;
      }
      else if( inext > -1)
      {
        /* Array remains the same size, but must subtract 1 from each index for
        reactions above the reaction that has been removed */
        for ( i = inext ; i < TC_nFallReac_ ; ++i )
        {
          TC_reacPfal_[i] -= 1 ;
        }     
      }
    }

    /* Third-body reactions */
    if ( TC_nThbReac_ > 0 )
    {
      /* Third-body index, if applicable */
      ispec = -1;
      inext = -1;
      for ( i = 0 ; i < TC_nThbReac_ ; i++ )
      {
        if ( TC_reacTbdy_[i] == curIdx )
        { 
          ispec = i ; 
          break ;
        }
        else if( TC_reacTbdy_[i] > curIdx ) /* Moved too far */
        {
          inext = i ; /* Index of first reaction with index greater than curIdx */
          break ;
        }
      }
      if( ispec > -1 )
      {
        /* Move reaction number array down by 1 and subtract 1 from each value */
        for ( i = ispec ; i < TC_nThbReac_-1 ; i++ )
        {
          TC_reacTbdy_[i] = TC_reacTbdy_[i+1] - 1;
        }
        numthb = TC_nThbReac_ - ispec - 1 ;
        memmove(&TC_reacTbno_[ispec],&TC_reacTbno_[ispec+1],
          numthb*sizeof(int) );  
        memmove(&TC_specTbdIdx_[ispec*TC_maxTbInReac_], &TC_specTbdIdx_[(ispec+1)*TC_maxTbInReac_],
          numthb*TC_maxTbInReac_*sizeof(int) );
        memmove(&TC_specTbdEff_[ispec*TC_maxTbInReac_], &TC_specTbdEff_[(ispec+1)*TC_maxTbInReac_],
          numthb*TC_maxTbInReac_*sizeof(double) );

        --TC_nThbReac_;
      }
      else if( inext > -1 )
      {
        /* Array remains the same size, but must subtract 1 from each index for
        reactions above the reaction that has been removed */
        for ( i = inext ; i < TC_nThbReac_ ; ++i )
        {
          TC_reacTbdy_[i] -= 1 ;
        }     
      }
    }

    /* Reactions with real stoichiometric coefficients */
    if ( TC_nRealNuReac_ > 0 ) 
    {
      /* Unlike other cases here, ispec is stored explicitly.
         TC_reacScoef_ was only rearranged as far as finalIdx on previous iterations,
         so add j to compensate (i.e. the relevant entry is still numbered by original 
         reaction indices) */
      ispec = TC_reacScoef_[curIdx+j] ; 
      if( ispec > -1)
      {    
        /* Move reaction number array down by 1 and subtract 1 from each value. */
        /* Also set TC_reacScoef_ to new values. */
        for ( i = ispec ; i < TC_nRealNuReac_-1 ; i++ )
        {
           /* TC_reacScoef will be moved to new position below, so do not subtract 1 here. Add
              j because relevant entries of TC_reacScoef are still numbered by original
              reaction indices. */
          TC_reacScoef_[TC_reacRnu_[i+1]+j] = i ;
          TC_reacRnu_[i] = TC_reacRnu_[i+1] - 1;
        }
        memmove(&TC_reacRealNuki_[ispec*TC_maxSpecInReac_], &TC_reacRealNuki_[(ispec+1)*TC_maxSpecInReac_],
          (TC_nRealNuReac_-ispec-1)*TC_maxSpecInReac_*sizeof(double) );

        /* sum(real nu) for each reaction */
        memmove(&TC_sigRealNu[ispec], &TC_sigRealNu[ispec+1],
          (TC_nRealNuReac_-ispec-1)*sizeof(double) );

        /* NuIJ=NuII''-NuIJ' for real reactions */
        memmove(&TC_RealNuIJ[ispec*TC_Nspec_], &TC_RealNuIJ[(ispec+1)*TC_Nspec_],
          (TC_nRealNuReac_-ispec-1)*TC_Nspec_*sizeof(double) );

        --TC_nRealNuReac_;
      }
      else
      {
        /* Find first reaction with index greater than ireac. Unlike other cases, 
        this did not occur above */
        inext = -1;
        for ( i = 0 ; i < TC_nRealNuReac_; ++i )
        {
          if( TC_reacRnu_[i] > ireac )
          {
            inext = i ; /* Index of first reaction with index greater than ireac */
            break ;
          }
        }
        if( inext > -1 )
        {
          for ( i = inext ; i < TC_nRealNuReac_ ; ++i )
          {
            TC_reacRnu_[i] -= 1 ;
            /* Note: TC_reacScoef_ is moved down below and TC_nRealNuReac has not changed, 
            so no further changes needed. */
          }
        }
      }
    }

    /* Now safe to deal with TC_reacScoef_ */
    memmove(&TC_reacScoef_[destIdx], &TC_reacScoef_[ireac], numMoveItems * sizeof(int));

    /* Arbitrary reaction orders */
    if ( TC_nOrdReac_ > 0 ) 
    {
      ispec = -1;
      inext = -1;
      /* Arbitrary order reaction index, if applicable */
      for ( i = 0 ; i < TC_nOrdReac_; i++ )
      {
        if ( TC_reacAOrd_[i] == curIdx )
        { 
          ispec = i ; 
          break ;
        }
        else if( TC_reacAOrd_[i] > curIdx ) /* Moved too far */
        {
          inext = i ; /* Index of first reaction with index greater than curIdx */
          break ;
        }
      }
      if( ispec > -1)
      {   
        /* Move reaction number array down by 1 and subtract 1 from each value */
        for ( i = ispec ; i < TC_nOrdReac_-1 ; i++ )
        {
          TC_reacAOrd_[i] = TC_reacAOrd_[i+1] - 1;
        }     
        memmove(&TC_specAOidx_[ispec*TC_maxOrdPar_], &TC_specAOidx_[(ispec+1)*TC_maxOrdPar_],
          (TC_nOrdReac_-ispec-1)*TC_maxOrdPar_*sizeof(int) );
        memmove(&TC_specAOval_[ispec*TC_maxOrdPar_], &TC_specAOval_[(ispec+1)*TC_maxOrdPar_],
          (TC_nOrdReac_-ispec-1)*TC_maxOrdPar_*sizeof(double) );

        --TC_nOrdReac_;
      }
      else if( inext > -1)
      {
        /* Array remains the same size, but must subtract 1 from each index for
        reactions above the reaction that has been removed */
        for ( i = inext ; i < TC_nOrdReac_ ; ++i )
        {
          TC_reacAOrd_[i] -= 1 ;
        }     
      }
    }

    /* sum(nu) and TC_kc_coeff for each reaction */
    memmove(&TC_sigNu[destIdx], &TC_sigNu[ireac], numMoveItems * sizeof(int) );
    memmove(&TC_kc_coeff[destIdx], &TC_kc_coeff[ireac], numMoveItems * sizeof(double) );

    /* NuIJ=NuII''-NuIJ' */
    memmove(&TC_NuIJ[destIdx*TC_Nspec_], &TC_NuIJ[ireac*TC_Nspec_],
      numMoveItems * TC_Nspec_ * sizeof(int) );

    /* TC_chg.c forward array */
    if( TC_ArhenForChg_ == 1)
    {
      memmove(&TC_reacArhenForSave_[3*destIdx], &TC_reacArhenForSave_[3*ireac],
        3*numMoveItems*sizeof(double) );    
    }

    destIdx += numMoveItems ;
  }

  TC_Nreac_ -= numRemoveReacs ;

  return ;

} /* end of TC_removeReaction() */


/** 
  \brief Restores reaction mechanism to original state. Any changes from TC_chgArhenFor and TC_chgArhenRev are also reset.
*/
void TC_restoreReactions()
{   
  if ( TC_initRemoved == 0 ) return;

  if( TC_tab_ == 1 ) /* Not designed for use with tables. */
  {
    fprintf(stderr,"Reaction removal with tables has not been implemented!\n");
    exit(1) ;
  }

  TC_Nreac_ = TC_NreacBackup_;

  /* If reaction is reversible */
  memcpy(TC_isRev_, TC_isRevBackup_, TC_Nreac_*sizeof(int) );

  /* Number of reactants and products */
  memcpy(TC_reacNrp_, TC_reacNrpBackup_, TC_Nreac_*sizeof(int) );

  /* Number of reactants */
  memcpy(TC_reacNreac_, TC_reacNreacBackup_, TC_Nreac_*sizeof(int) );

  /* Number of products */
  memcpy(TC_reacNprod_, TC_reacNprodBackup_, TC_Nreac_*sizeof(int) );

  /* Stoichiometric coefficients */
  memcpy(TC_reacNuki_,    TC_reacNukiBackup_, TC_maxSpecInReac_*TC_Nreac_*sizeof(int) );
  memcpy(TC_reacNukiDbl_, TC_reacNukiDblBackup_, TC_maxSpecInReac_*TC_Nreac_*sizeof(double) );
  memcpy(TC_reacSidx_,    TC_reacSidxBackup_, TC_maxSpecInReac_*TC_Nreac_*sizeof(int) );

  /* TC_reacScoef_ determines if stoichiometric coefficients are int or real */
  memcpy(TC_reacScoef_, TC_reacScoefBackup_, TC_Nreac_*sizeof(int) );

  /* Arrhenius parameters */
  memcpy(TC_reacArhenFor_, TC_reacArhenForBackup_, 3*TC_Nreac_*sizeof(double) );

  if( TC_ArhenForChg_ == 1)
  {
    free(TC_reacArhenForSave_);
    TC_reacArhenForSave_ = 0;
    TC_ArhenForChg_ = 0;
  }

  if (TC_NrevBackup_ > 0 ) /* No. reactions with reverse Arrhenius parameters */
  {
    TC_nRevReac_ = TC_NrevBackup_;    
    memcpy(TC_reacArhenRev_, TC_reacArhenRevBackup_, 3*TC_nRevReac_*sizeof(double) ) ;
  }

  if( TC_ArhenRevChg_ == 1)
  {
    free(TC_reacArhenRevSave_);
    TC_reacArhenRevSave_ = 0;
    TC_ArhenRevChg_ = 0;
  }

  /* Pressure-dependent reactions */
  if ( TC_NfalBackup_ > 0 )
  {
    TC_nFallReac_ = TC_NfalBackup_;

    memcpy(TC_reacPfal_,  TC_reacPfalBackup_,  TC_nFallReac_*sizeof(int) );
    memcpy(TC_reacPtype_, TC_reacPtypeBackup_, TC_nFallReac_*sizeof(int) );
    memcpy(TC_reacPlohi_, TC_reacPlohiBackup_, TC_nFallReac_*sizeof(int) );
    memcpy(TC_reacPspec_, TC_reacPspecBackup_, TC_nFallReac_*sizeof(int) );
    memcpy(TC_reacPpar_,  TC_reacPparBackup_,  TC_nFallReac_*TC_nFallPar_*sizeof(double) );
  }

  /* Third-body reactions */
  if ( TC_NthbBackup_ > 0 )
  {
    TC_nThbReac_ = TC_NthbBackup_ ;

    memcpy(TC_reacTbdy_,   TC_reacTbdyBackup_,   TC_nThbReac_ * sizeof(int) );
    memcpy(TC_reacTbno_,   TC_reacTbnoBackup_,   TC_nThbReac_ * sizeof(int) );
    memcpy(TC_specTbdIdx_, TC_specTbdIdxBackup_, TC_nThbReac_ * TC_maxTbInReac_ * sizeof(int) );
    memcpy(TC_specTbdEff_, TC_specTbdEffBackup_, TC_nThbReac_ * TC_maxTbInReac_ * sizeof(double) );
  }

  /* Reactions with real stoichiometric coefficients */
  if ( TC_nRealNuReacBackup_ > 0 ) 
  {
    TC_nRealNuReac_ =  TC_nRealNuReacBackup_ ;

    memcpy(TC_reacRnu_,      TC_reacRnuBackup_,      TC_nRealNuReac_ * sizeof(int) );
    memcpy(TC_reacRealNuki_, TC_reacRealNukiBackup_, TC_nRealNuReac_ * TC_maxSpecInReac_ * sizeof(double) );
    memcpy(TC_sigRealNu,     TC_sigRealNuBackup,     TC_nRealNuReac_ * sizeof(double) );
    memcpy(TC_RealNuIJ,      TC_RealNuIJBackup,      TC_nRealNuReac_ * TC_Nspec_ * sizeof(double) );
  }

  /* Arbitrary reaction orders */
  if ( TC_nOrdReacBackup_ > 0 ) 
  {
    TC_nOrdReac_ = TC_nOrdReacBackup_ ;

    memcpy(TC_reacAOrd_,  TC_reacAOrdBackup_,  TC_nOrdReac_ * sizeof(int) );
    memcpy(TC_specAOidx_, TC_specAOidxBackup_, TC_nOrdReac_ * TC_maxOrdPar_ * sizeof(int) );
    memcpy(TC_specAOval_, TC_specAOvalBackup_, TC_nOrdReac_ * TC_maxOrdPar_ * sizeof(double) );
  }

  /* sum(nu) and TC_kc_coeff for each reaction */
  memcpy(TC_sigNu,    TC_sigNuBackup, TC_Nreac_ * sizeof(int) );
  memcpy(TC_kc_coeff, TC_kc_coeffBackup, TC_Nreac_ * sizeof(double) );

  /* NuIJ=NuII''-NuIJ' */
  memcpy(TC_NuIJ, TC_NuIJBackup, TC_Nreac_ * TC_Nspec_ * sizeof(int) );

  return ;

} /* end of TC_restoreReactions() */
