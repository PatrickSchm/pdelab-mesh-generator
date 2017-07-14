#ifndef DUNE_PDELAB_ONESTEP_PRESTAGEENGINE_HH
#define DUNE_PDELAB_ONESTEP_PRESTAGEENGINE_HH

#include <dune/pdelab/residualcoupling/gridoperator/onestep/enginebase_rcoupling.hh>
#include <cmath>
#include <vector>

namespace Dune
{
namespace PDELab
{

/**
 \brief The local assembler engine for one step methods which
 assembles the constant part of the residual vector

 \tparam LA The local one step assembler

 */
template<typename OSLA, typename GOSL>
class OneStepLocalPreStageAssemblerEngineCoupling: public OneStepLocalAssemblerEngineBaseCoupling <
    OSLA,
    typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
    typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine,
    typename OSLA::LocalAssemblerCouplingDT0::LocalResidualAssemblerEngine,
    typename OSLA::LocalAssemblerCouplingDT1::LocalResidualAssemblerEngine >
{

    typedef OneStepLocalAssemblerEngineBaseCoupling<OSLA,
            typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
            typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine,
            typename OSLA::LocalAssemblerCouplingDT0::LocalResidualAssemblerEngine,
            typename OSLA::LocalAssemblerCouplingDT1::LocalResidualAssemblerEngine> BaseT;

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

    //! Types of the subordinate assemblers and engines
    //! @{
    typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
    typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;
    typedef typename OSLA::LocalAssemblerCouplingDT0 LocalAssemblerCouplingDT0;
    typedef typename OSLA::LocalAssemblerCouplingDT1 LocalAssemblerCouplingDT1;

    typedef typename LocalAssemblerDT0::LocalResidualAssemblerEngine ResidualEngineDT0;
    typedef typename LocalAssemblerDT1::LocalResidualAssemblerEngine ResidualEngineDT1;
    typedef typename LocalAssemblerCouplingDT0::LocalResidualAssemblerEngine ResidualEngineCouplingDT0;
    typedef typename LocalAssemblerCouplingDT1::LocalResidualAssemblerEngine ResidualEngineCouplingDT1;

    //! @}

    //! The type of the residual vector
    typedef typename OSLA::Traits::Residual Residual;
    typedef typename Residual::ElementType ResidualElement;

    //! The type of the solution vector
    typedef typename OSLA::Traits::Solution Solution;
    typedef typename Solution::ElementType SolutionElement;

    typedef typename GOSL::Traits::Domain SolutionSl;

    //! The type for real numbers
    typedef typename OSLA::Real Real;

    //! The type of the solution container
    typedef std::vector<Solution*> Solutions;
    typedef std::vector<SolutionSl*> SolutionsSl;

    /**
     \brief Constructor

     \param [in] la_ The local assembler object which
     creates this engine
     */
    OneStepLocalPreStageAssemblerEngineCoupling(LocalAssembler & la_) :
        BaseT(la_), invalid_residual(static_cast<Residual*>(0)), invalid_solutions(
            static_cast<Solutions*>(0)), invalid_solutionsSl(
                static_cast<SolutionsSl*>(0)), const_residual_0(
                    invalid_residual), const_residual_1(
                        invalid_residual), solutions(invalid_solutions), solutionsSl(
                            invalid_solutionsSl)
    {
    }

    //! Query methods for the global grid assembler
    //! @{
    bool requireSkeleton() const
    {
        return lae0->requireSkeleton() || lae1->requireSkeleton();
    }

    //! @}

    //! Set current solution vector. Must be called before
    //! setConstResidual()! Should be called prior to assembling.
    void setSolutions(const Solutions & solutions_)
    {
        solutions = &solutions_;
    }

    void setSolutionsOld(const Solutions & solutionsOld_)
    {
        solutionsOld = &solutionsOld_;
    }

    void setSolutionsSl(const SolutionsSl & solutionsSl_)
    {
        solutionsSl = &solutionsSl_;
    }

    //! Set current const residual vector. Should be called prior to
    //! assembling.
    void setConstResiduals(Residual & const_residual_0_,
                           Residual & const_residual_1_)
    {
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;


        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(
            la.la0.localResidualAssemblerEngine(*const_residual_0,
                                                *((*solutions)[0])));
        setLocalAssemblerEngineDT1(
            la.la1.localResidualAssemblerEngine(*const_residual_1,
                                                *((*solutions)[0])));
        setLocalAssemblerEngineCouplingDT0(
            la.lac0.localResidualAssemblerEngineCoupling(
                *const_residual_0, *((*solutions)[0])),
            *((*solutionsSl)[0]));
        setLocalAssemblerEngineCouplingDT1(
            la.lac1.localResidualAssemblerEngineCoupling(
                *const_residual_1, *((*solutions)[0])),
            *((*solutionsSl)[0]));
    }

    //! Set current const residual vector. Should be called prior to
    //! assembling.
    void setConstResidual(Residual & const_residual_)
    {
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;


        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(
            la.la0.localResidualAssemblerEngine(*const_residual_0,
                                                *((*solutions)[0])));
        setLocalAssemblerEngineDT1(
            la.la1.localResidualAssemblerEngine(*const_residual_1,
                                                *((*solutions)[0])));
        setLocalAssemblerEngineCouplingDT0(
            la.lac0.localResidualAssemblerEngineCoupling(
                *const_residual_0, *((*solutions)[0]), *((*solutions)[0]),
                *((*solutionsSl)[0])));
        setLocalAssemblerEngineCouplingDT1(
            la.lac1.localResidualAssemblerEngineCoupling(
                *const_residual_1, *((*solutions)[0]), *((*solutions)[0]),
                *((*solutionsSl)[0])));
//                std::cout <<  ((*solutionsSl)[0]) << std::endl;
    }

    //! Methods for loading of the local function's
    //! coefficients. These methods are blocked. The loading of the
    //! coefficients is done in each assemble call.
    //!@{
    template<typename LFSU>
    void loadCoefficientsLFSUInside(const LFSU & lfsu_s)
    {
    }

    template<typename LFSU>
    void loadCoefficientsLFSUInsideOld(const LFSU & lfsu_s)
    {
    }
    template<typename LFSU>
    void loadCoefficientsLFSUOutside(const LFSU & lfsu_n)
    {
    }
    template<typename LFSU>
    void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c)
    {
    }
    //! @}

    //! Method setting time for la1 local assembler.
    //! This function must be called for explicit methods
    //! before jacobian_engine->assemble.. was called
    void setTimeInLastStage()
    {
        la.la1.setTime(la.time + la.osp_method->d(la.stage) * la.dt);
    }

    //! Notifier functions, called immediately before and after assembling
    //! @{
    void preAssembly()
    {
        lae0->preAssembly();
        lae1->preAssembly();
        laec0->preAssembly();
        laec1->preAssembly();

        *const_residual_0 = 0.0;
        *const_residual_1 = 0.0;


        // Extract the coefficients of the time step scheme
        a.resize(la.stage);
        b.resize(la.stage);
        d.resize(la.stage);
        do0.resize(la.stage);
        do1.resize(la.stage);
        doc0.resize(la.stage);
        doc1.resize(la.stage);
        for (int i = 0; i < la.stage; ++i)
        {
            a[i] = la.osp_method->a(la.stage, i);
            b[i] = la.osp_method->b(la.stage, i);
            d[i] = la.osp_method->d(i);
            do0[i] = (std::abs(b[i]) > 1E-6);
            do1[i] = (std::abs(a[i]) > 1E-6);
            doc0[i] = (std::abs(b[i]) > 1E-6);
            doc1[i] = (std::abs(a[i]) > 1E-6);
        }

        // prepare local operators for stage
        la.la0.preStage(la.time + la.osp_method->d(la.stage) * la.dt,
                        la.stage);
        la.la1.preStage(la.time + la.osp_method->d(la.stage) * la.dt,
                        la.stage);
        la.lac0.preStage(la.time + la.osp_method->d(la.stage) * la.dt,
                         la.stage);
        la.lac1.preStage(la.time + la.osp_method->d(la.stage) * la.dt,
                         la.stage);
    }

    template<typename GFSU, typename GFSV>
    void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
    {
        lae0->postAssembly(gfsu, gfsv);
        lae1->postAssembly(gfsu, gfsv);
        laec0->postAssembly(gfsu, gfsv);
        laec1->postAssembly(gfsu, gfsv);
    }
    //! @}

    //! @ Assembling methods
    //! @{

    template<typename EG, typename LFSU, typename LFSV>
    void assembleUVVolume(const EG & eg, const LFSU & lfsu,
                          const LFSV & lfsv, int iMat)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu);
            lae1->loadCoefficientsLFSUInside(lfsu);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleUVVolume(eg, lfsu, lfsv, iMat);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleUVVolume(eg, lfsu, lfsv, iMat);
            }
        }
    }


    template<typename IG, typename EGC, typename EGP, typename LFSU, typename LFSV, typename LFSUC, typename LFSVC>
    void assembleCBoundary(const IG& ig, const EGC& egC, const EGP& egP,
                           const LFSU& lfsu, const LFSV& lfsv, const LFSU& lfsu_pair, const LFSV& lfsv_pair, const LFSUC& lfsuC,
                           const LFSVC& lfsvC) {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.lac0.setTime(la.time + d[s] * la.dt);
            la.lac1.setTime(la.time + d[s] * la.dt);

            laec0->setSolution(*((*solutions)[s]));
            laec1->setSolution(*((*solutions)[s]));

            laec0->loadCoefficientsLFSUInside(lfsu);
            laec1->loadCoefficientsLFSUInside(lfsu);

            laec0->setSolutionOld(*((*solutionsOld)[0]));
            laec1->setSolutionOld(*((*solutionsOld)[0]));

            laec0->loadCoefficientsLFSUInsideOld(lfsu_pair);
            laec1->loadCoefficientsLFSUInsideOld(lfsu_pair);

            laec0->setSolutionSl(*((*solutionsSl)[0]));
            laec1->setSolutionSl(*((*solutionsSl)[0]));

            laec0->loadCoefficientsLFSUInsideSl(lfsuC);
            laec1->loadCoefficientsLFSUInsideSl(lfsuC);
            if (do0[s])
            {
                la.lac0.setWeight(b[s] * la.dt_factor0);
                laec0->assembleCBoundary(ig, egC, egP, lfsu, lfsv, lfsu_pair, lfsv_pair,
                                         lfsuC, lfsvC);
            }
            if (do1[s])
            {
                la.lac1.setWeight(a[s] * la.dt_factor1);
                laec1->assembleCBoundary(ig, egC, egP, lfsu, lfsv, lfsu_pair, lfsv_pair,
                                         lfsuC, lfsvC);
            }
        }
    }

    template<typename IG, typename EGC, typename LFSU, typename LFSV, typename LFSUC, typename LFSVC, typename FH>
    void assembleCBoundaryReversed(const IG& ig, const EGC& egC, const LFSU& lfsu,
                                   const LFSV& lfsv, const LFSUC& lfsuC,
                                   const LFSVC& lfsvC, std::vector<FH>& fractureHeight) {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.lac0.setTime(la.time + d[s] * la.dt);
            la.lac1.setTime(la.time + d[s] * la.dt);
            laec0->setSolution(*((*solutions)[s]));
            laec0->setSolutionSl(*((*solutionsSl)[s]));

            laec1->setSolution(*((*solutions)[s]));
            laec1->setSolutionSl(*((*solutionsSl)[s]));

            laec0->loadCoefficientsLFSUInside(lfsu);
            laec1->loadCoefficientsLFSUInside(lfsu);

            laec0->loadCoefficientsLFSUInsideSl(lfsuC);
            laec1->loadCoefficientsLFSUInsideSl(lfsuC);
            if (do0[s])
            {
                la.lac0.setWeight(b[s] * la.dt_factor0);
                laec0->assembleCBoundaryReversed(ig, egC, lfsu, lfsv,
                                                 lfsuC, lfsvC, fractureHeight);
            }
            if (do1[s])
            {
                la.lac1.setWeight(a[s] * la.dt_factor1);
                laec1->assembleCBoundaryReversed(ig, egC, lfsu, lfsv,
                                                 lfsuC, lfsvC, fractureHeight);
            }
        }
    }

    template<typename EG, typename LFSV>
    void assembleVVolume(const EG & eg, const LFSV & lfsv)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVVolume(eg, lfsv);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVVolume(eg, lfsv);
            }

        }
    }

    template<typename IG, typename LFSU_S, typename LFSV_S,
             typename LFSU_N, typename LFSV_N>
    void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s,
                            const LFSV_S & lfsv_s, const LFSU_N & lfsu_n,
                            const LFSV_N & lfsv_n)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);
            lae0->loadCoefficientsLFSUOutside(lfsu_n);
            lae1->loadCoefficientsLFSUOutside(lfsu_n);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->setSolution(*((*solutions)[s]));
                lae0->assembleUVSkeleton(ig, lfsu_s, lfsv_s, lfsu_n,
                                         lfsv_n);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->setSolution(*((*solutions)[s]));
                lae1->assembleUVSkeleton(ig, lfsu_s, lfsv_s, lfsu_n,
                                         lfsv_n);
            }
        }
    }

    template<typename IG, typename LFSV_S, typename LFSV_N>
    void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s,
                           const LFSV_N & lfsv_n)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVSkeleton(ig, lfsv_s, lfsv_n);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVSkeleton(ig, lfsv_s, lfsv_n);
            }
        }
    }

    template<typename IG, typename LFSU_S, typename LFSV_S>
    void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s,
                            const LFSV_S & lfsv_s)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleUVBoundary(ig, lfsu_s, lfsv_s);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleUVBoundary(ig, lfsu_s, lfsv_s);
            }
        }
    }

    template<typename IG, typename LFSV_S>
    void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVBoundary(ig, lfsv_s);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVBoundary(ig, lfsv_s);
            }
        }
    }

    template<typename IG, typename LFSU_S, typename LFSV_S>
    void assembleUVProcessor(const IG & ig, const LFSU_S & lfsu_s,
                             const LFSV_S & lfsv_s)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleUVProcessor(ig, lfsu_s, lfsv_s);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleUVProcessor(ig, lfsu_s, lfsv_s);
            }
        }
    }

    template<typename IG, typename LFSV_S>
    void assembleVProcessor(const IG & ig, const LFSV_S & lfsv_s)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVProcessor(ig, lfsv_s);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVProcessor(ig, lfsv_s);
            }
        }
    }

    template<typename IG, typename LFSU_S, typename LFSV_S,
             typename LFSU_N, typename LFSV_N, typename LFSU_C,
             typename LFSV_C>
    void assembleUVEnrichedCoupling(const IG & ig,
                                    const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                    const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                    const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu_s);
            lae1->loadCoefficientsLFSUInside(lfsu_s);

            lae0->loadCoefficientsLFSUOutside(lfsu_n);
            lae1->loadCoefficientsLFSUOutside(lfsu_n);

            lae0->loadCoefficientsLFSUCoupling(lfsu_c);
            lae1->loadCoefficientsLFSUCoupling(lfsu_c);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleUVEnrichedCoupling(ig, lfsu_s, lfsv_s,
                                                 lfsu_n, lfsv_n, lfsu_c, lfsv_c);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleUVEnrichedCoupling(ig, lfsu_s, lfsv_s,
                                                 lfsu_n, lfsv_n, lfsu_c, lfsv_c);
            }
        }
    }

    template<typename IG, typename LFSV_S, typename LFSV_N,
             typename LFSV_C>
    void assembleVEnrichedCoupling(const IG & ig, const LFSV_S & lfsv_s,
                                   const LFSV_N & lfsv_n, const LFSV_C & lfsv_c)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVEnrichedCoupling(ig, lfsv_s, lfsv_n,
                                                lfsv_c);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVEnrichedCoupling(ig, lfsv_s, lfsv_n,
                                                lfsv_c);
            }

        }
    }

    template<typename EG, typename LFSU, typename LFSV>
    void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu,
                                      const LFSV & lfsv)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            lae0->setSolution(*((*solutions)[s]));
            lae1->setSolution(*((*solutions)[s]));

            lae0->loadCoefficientsLFSUInside(lfsu);
            lae1->loadCoefficientsLFSUInside(lfsu);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleUVVolumePostSkeleton(eg, lfsu, lfsv);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleUVVolumePostSkeleton(eg, lfsu, lfsv);
            }

        }
    }

    template<typename EG, typename LFSV>
    void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
    {
        for (int s = 0; s < la.stage; ++s)
        {
            // Reset the time in the local assembler
            la.la0.setTime(la.time + d[s] * la.dt);
            la.la1.setTime(la.time + d[s] * la.dt);

            if (do0[s])
            {
                la.la0.setWeight(b[s] * la.dt_factor0);
                lae0->assembleVVolumePostSkeleton(eg, lfsv);
            }

            if (do1[s])
            {
                la.la1.setWeight(a[s] * la.dt_factor1);
                lae1->assembleVVolumePostSkeleton(eg, lfsv);
            }
        }
    }
//! @}

private:

//! Default value indicating an invalid residual pointer
    Residual * const invalid_residual;

//! Default value indicating an invalid solution pointer
    Solutions * const invalid_solutions;
    SolutionsSl * const invalid_solutionsSl;

//! Pointer to the current constant part residual vector in
//! which to assemble the residual corresponding to the operator
//! representing the time derivative of order zero and one.
//! @{
    Residual * const_residual_0;
    Residual * const_residual_1;

//! @}

//! Pointer to the current residual vector in which to assemble
    const Solutions * solutions;
    const Solutions * solutionsOld;
    const SolutionsSl * solutionsSl;

//! Coefficients of time stepping scheme
    std::vector<Real> a;
    std::vector<Real> b;
    std::vector<Real> d;
    std::vector<bool> do0;
    std::vector<bool> do1;
    std::vector<bool> doc0;
    std::vector<bool> doc1;

};
// End of class OneStepLocalPreStageAssemblerEngineCoupling

}
}
#endif
