#ifndef DUNE_ONE_STEP_RESIDUALENGINE_HH
#define DUNE_ONE_STEP_RESIDUALENGINE_HH

#include <dune/pdelab/residualcoupling/gridoperator/onestep/enginebase_rcoupling.hh>

#include <cmath>

namespace Dune
{
namespace PDELab
{

/**
 \brief The local assembler engine for one step methods which
 assembles the residual vector
 \tparam LA The local one step assembler
 */
template<typename OSLA, typename GOSL>
class OneStepLocalResidualAssemblerEngineCoupling: public OneStepLocalAssemblerEngineBaseCoupling <
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
    typedef OSLA OneStepLocalAssembler;

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

    //! The type of the solution vector
    typedef typename OSLA::Traits::Solution Solution;
    typedef typename GOSL::Traits::Domain SolutionSl;

    //! The type for real numbers
    typedef typename OSLA::Real Real;

    typedef OSLA LocalAssembler;

    /**
     \brief Constructor
     \param [in] local_assembler_ The local assembler object which
     creates this engine
     */
    OneStepLocalResidualAssemblerEngineCoupling(
        const LocalAssembler & local_assembler_) :
        BaseT(local_assembler_), invalid_residual(
            static_cast<Residual*>(0)), invalid_solution(
                static_cast<Solution*>(0)), invalid_solutionSl(
                    static_cast<SolutionSl*>(0)), residual_0(
                        invalid_residual), residual_1(invalid_residual), residual_c0(
                            invalid_residual), residual_c1(invalid_residual), const_residual_0(
                                invalid_residual), const_residual_1(
                                    invalid_residual), const_residual_c0(
                                        invalid_residual), const_residual_c1(
                                            invalid_residual), solution(invalid_solution), solutionSl(
                                                invalid_solutionSl)
    {
    }

    //! Set current solution vector. Must be called before
    //! setResidual(). Should be called prior to assembling.
    void setSolution(const Solution & solution_)
    {
        solution = &solution_;
    }

        void setSolutionOld(const Solution & solutionOld_)
    {
        solutionOld = &solutionOld_;
    }

    void setSolutionSl(const SolutionSl & solutionSl_)
    {
        solutionSl = &solutionSl_;
//                Dune::printvector(std::cout, solutionSl_.base(), "onestep solutionSl ----------",
//                        "row");
    }

    //! Set current const residual vector. Must be called before
    //! setResidual(). Should be called prior to assembling.
    void setConstResidual(const Residual &const_residual_)
    {
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;
        const_residual_c0 = &const_residual_;
        const_residual_c1 = &const_residual_;
    }

    //! Set current const residual vector. Should be called prior to
    //! assembling.
    void setResidual(Residual & residual_)
    {
        residual_0 = &residual_;
        residual_1 = &residual_;
        residual_c0 = &residual_;
        residual_c1 = &residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(
            la.la0.localResidualAssemblerEngine(*residual_0,
                                                *solution));
        setLocalAssemblerEngineDT1(
            la.la1.localResidualAssemblerEngine(*residual_1,
                                                *solution));
        setLocalAssemblerEngineCouplingDT0(
            la.lac0.localResidualAssemblerEngineCoupling(
                *residual_0, *solution, *solutionOld, *solutionSl));
        setLocalAssemblerEngineCouplingDT1(
            la.lac1.localResidualAssemblerEngineCoupling(
                *residual_1, *solution, *solutionOld, *solutionSl));
    }

    //! Set current const residual vectors. Must be called before
    //! setResidual(). Should be called prior to assembling. Here,
    //! separate vectors are used for the operators corresponding to
    //! the time dervatives of order one and zero.
    void setConstResiduals(const Residual &const_residual_0_,
                           const Residual &const_residual_1_,
                           const Residual &const_residual_c0_,
                           const Residual &const_residual_c1_)
    {
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;
        const_residual_c0 = &const_residual_c0_;
        const_residual_c1 = &const_residual_c1_;
    }

    //! Set current const residual vectors. Should be called prior
    //! to assembling. Here, separate vectors are used for the
    //! operators corresponding to the time dervatives of order one
    //! and zero.
    void setResiduals(Residual & residual_0_, Residual & residual_1_,
                      Residual & residual_c0_, Residual & residual_c1_)
    {
        residual_0 = &residual_0_;
        residual_1 = &residual_1_;
        residual_c0 = &residual_c0_;
        residual_c1 = &residual_c1_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(
            la.la0.localResidualAssemblerEngine(*residual_0,
                                                *solution));
        setLocalAssemblerEngineDT1(
            la.la1.localResidualAssemblerEngine(*residual_1,
                                                *solution));
        setLocalAssemblerEngineCouplingDT0(
            la.lac0.localResidualAssemblerEngineCouplingDT0(
                *residual_0, *solution, *solutionSl));
        setLocalAssemblerEngineCouplingDT1(
            la.lac1.localResidualAssemblerEngineCouplingDT1(
                *residual_1, *solution, *solutionSl));
    }

    //! When multiple engines are combined in one assembling
    //! procedure, this method allows to reset the weights which may
    //! have been changed by the other engines.
    void setWeights()
    {
        la.la0.setWeight(b_rr * la.dt_factor0);
        la.la1.setWeight(la.dt_factor1);

        la.lac0.setWeight(b_rr * la.dt_factor0);
        la.lac1.setWeight(la.dt_factor1);
    }

    //! Notifier functions, called immediately before and after assembling
    //! @{
    void preAssembly()
    {
        lae0->preAssembly();
        lae1->preAssembly();
        laec0->preAssembly();
        laec1->preAssembly();

        // Extract the coefficients of the time step scheme
        b_rr = la.osp_method->b(la.stage, la.stage);
        d_r = la.osp_method->d(la.stage);
        implicit = std::abs(b_rr) > 1e-6;

        // prepare local operators for stage
        la.la0.setTime(la.time + d_r * la.dt);
        la.la1.setTime(la.time + d_r * la.dt);
        la.lac0.setTime(la.time + d_r * la.dt);
        la.lac1.setTime(la.time + d_r * la.dt);

        setWeights();
    }

    template<typename GFSU, typename GFSV>
    void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
    {

        // Update residual vectors with constant part
        assert(const_residual_0 != invalid_residual);
        assert(const_residual_1 != invalid_residual);
        // assert(const_residual_c0 != invalid_residual);
        // assert(const_residual_c1 != invalid_residual);

        *residual_0 += *const_residual_0;
        if (residual_0 != residual_1)
        {
            assert(const_residual_0 != const_residual_1);
            *residual_1 += *const_residual_1;
        }

        // if (residual_0 != residual_c0)
        // {
        //     assert(const_residual_0 != const_residual_0);
        //     *residual_c0 += *const_residual_0;
        // }

        // if (residual_0 != residual_c1)
        // {
        //     assert(const_residual_0 != const_residual_1);
        //     *residual_c1 += *const_residual_1;
        // }


        lae0->postAssembly(gfsu, gfsv);
        lae1->postAssembly(gfsu, gfsv);
        laec0->postAssembly(gfsu, gfsv);
        laec1->postAssembly(gfsu, gfsv);
    }
    //! @}

private:

    //! Default value indicating an invalid residual pointer
    Residual * const invalid_residual;

    //! Default value indicating an invalid solution pointer
    Solution * const invalid_solution;
    SolutionSl * const invalid_solutionSl;

    //! Pointer to the current constant part residual vector in
    //! which to assemble the residual corresponding to the operator
    //! representing the time derivative of order zero and one.
    //! @{
    Residual * residual_0;
    Residual * residual_1;
    Residual * residual_c0;
    Residual * residual_c1;
    //! @}

    //! Pointer to the current constant part residual vectors in
    //! which to assemble the residual corresponding to the operator
    //! representing the time derivative of order zero and one.
    //! @{
    const Residual * const_residual_0;
    const Residual * const_residual_1;
    const Residual * const_residual_c0;
    const Residual * const_residual_c1;
    //! @}

    //! Pointer to the current residual vector in which to assemble
    const Solution * solution;
    const Solution * solutionOld;
    const SolutionSl * solutionSl;

    //! Coefficients of time stepping scheme
    Real b_rr, d_r;

};
// End of class OneStepLocalResidualAssemblerEngineCoupling

}
}

#endif