#ifndef DUNE_ONE_STEP_PATTERNENGINE_RCOUPLING_HH
#define DUNE_ONE_STEP_PATTERNENGINE_RCOUPLING_HH

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/residualcoupling/gridoperator/onestep/enginebase_rcoupling.hh>

namespace Dune{
  namespace PDELab{

    /**
       \brief The local assembler engine for OneStep sub triangulations which
       creates the matrix pattern

       \tparam LA The local udg assembler

    */
    template<typename OSLA>
    class OneStepLocalPatternAssemblerEngineCoupling
      : public OneStepLocalAssemblerEngineBaseCoupling<OSLA,
                                               typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
                                               typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine,
                                               typename OSLA::LocalAssemblerCouplingDT0::LocalPatternAssemblerEngine,
                                               typename OSLA::LocalAssemblerCouplingDT1::LocalPatternAssemblerEngine
                                               >
    {

      typedef OneStepLocalAssemblerEngineBaseCoupling<OSLA,
                                              typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
                                              typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine,
                                              typename OSLA::LocalAssemblerCouplingDT0::LocalPatternAssemblerEngine,
                                              typename OSLA::LocalAssemblerCouplingDT1::LocalPatternAssemblerEngine
                                              > BaseT;

      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::laec0;
      using BaseT::laec1;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;
      using BaseT::setLocalAssemblerEngineCouplingDT0;
      using BaseT::setLocalAssemblerEngineCouplingDT1;

    public:
      //! The type of the wrapping local assembler
      typedef OSLA LocalAssembler;

      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;
      typedef typename OSLA::LocalAssemblerCouplingDT0 LocalAssemblerCouplingDT0;
      typedef typename OSLA::LocalAssemblerCouplingDT1 LocalAssemblerCouplingDT1;

      //! The type of the matrix pattern container
      typedef typename LocalAssembler::Traits::MatrixPattern Pattern;
      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor

         \param [in] la_ The local assembler object which
         creates this engine
      */
      OneStepLocalPatternAssemblerEngineCoupling(const LocalAssembler & la_)
        : BaseT(la_),
          invalid_pattern(static_cast<Pattern*>(0)), pattern(invalid_pattern)
      {}

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){

        // Set pointer to global pattern
        pattern = &pattern_;

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineDT1(la.la1.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineCouplingDT0(la.lac0.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineCouplingDT1(la.lac1.localPatternAssemblerEngine(pattern_));
      }


      //! @name Notification functions
      //! @{
      void preAssembly(){
        implicit = la.osp_method->implicit();

        lae0->preAssembly();
        lae1->preAssembly();
        laec0->preAssembly();
        laec1->preAssembly();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
        laec0->postAssembly(gfsu,gfsv);
        laec1->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

    }; // End of class OneStepLocalPatternAssemblerEngineCoupling

  }
}

#endif
