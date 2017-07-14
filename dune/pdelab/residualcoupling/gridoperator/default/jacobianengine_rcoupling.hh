#ifndef DUNE_PDELAB_DEFAULT_JACOBIANENGINE_RCOUPLING_HH
#define DUNE_PDELAB_DEFAULT_JACOBIANENGINE_RCOUPLING_HH

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/diagonallocalmatrix.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/residualcoupling/gridoperator/common/localassemblerenginebase_rcoupling.hh>
#include <dune/pdelab/residualcoupling/localoperator/callswitch_rcoupling.hh>
// #include <dune/pdelab/localoperator/callswitch.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune
{
namespace PDELab
{

/**
 \brief The local assembler engine for DUNE grids which
 assembles the jacobian matrix

 \tparam LA The local assembler

 */
template<typename LA, typename LASL>
class DefaultLocalJacobianAssemblerEngineRCoupling: public LocalAssemblerEngineBaseRCoupling
{
public:

    template<typename TrialConstraintsContainer,
             typename TestConstraintsContainer>
    bool needsConstraintsCaching(const TrialConstraintsContainer& cu,
                                 const TestConstraintsContainer& cv)
    {
        return cu.containsNonDirichletConstraints()
               || cv.containsNonDirichletConstraints();
    }

    //! The type of the wrapping local assembler
    typedef LA LocalAssembler;

    //! The type of the local operator
    typedef typename LA::LocalOperator LOP;

    //! The local function spaces
    typedef typename LA::LFSU LFSU;
    typedef typename LA::LFSUCache LFSUCache;
    typedef typename LFSU::Traits::GridFunctionSpace GFSU;
    typedef typename LA::LFSV LFSV;
    typedef typename LA::LFSVCache LFSVCache;
    typedef typename LFSV::Traits::GridFunctionSpace GFSV;

    //! The type of the jacobian matrix
    typedef typename LA::Traits::Jacobian Jacobian;
    typedef typename Jacobian::ElementType JacobianElement;
    typedef typename Jacobian::template LocalView<LFSVCache, LFSUCache> JacobianView;

    //! The type of the solution vector
    typedef typename LA::Traits::Solution Solution;
    typedef typename Solution::ElementType SolutionElement;
    typedef typename Solution::template ConstLocalView<LFSUCache> SolutionView;

    typedef typename LA::Traits::SolutionSl SolutionSl;
    typedef typename SolutionSl::ElementType SolutionElementSl;

    typedef typename LASL::LFSU LFSUSL;
    typedef typename LASL::LFSUCache LFSUCacheSL;
    typedef typename LFSUSL::Traits::GridFunctionSpace GFSUSL;
    typedef typename LASL::LFSV LFSVSL;
    typedef typename LASL::LFSVCache LFSVCacheSL;
    typedef typename LFSVSL::Traits::GridFunctionSpace GFSVSL;

    typedef typename SolutionSl::template ConstLocalView<LFSUCacheSL> SolutionViewSl;

    typedef typename LASL::Traits::Jacobian JacobianSl;
    typedef typename JacobianSl::ElementType JacobianElementSl;
    typedef typename JacobianSl::template LocalView<LFSVCache, LFSUCache> JacobianViewSl;

    /**
     \brief Constructor

     \param [in] local_assembler_ The local assembler object which
     creates this engine
     */
    DefaultLocalJacobianAssemblerEngineRCoupling(
        const LocalAssembler & local_assembler_) :
        local_assembler(local_assembler_), lop(
            local_assembler_.lop), al_view(al, 1.0), al_sn_view(
                al_sn, 1.0), al_ns_view(al_ns, 1.0), al_nn_view(
                    al_nn, 1.0)
    {
    }

    //! Query methods for the global grid assembler
    //! @{
    bool requireSkeleton() const
    {
        return local_assembler.doAlphaSkeleton();
    }
    bool requireSkeletonTwoSided() const
    {
        return local_assembler.doSkeletonTwoSided();
    }
    bool requireUVVolume() const
    {
        return local_assembler.doAlphaVolume();
    }
    bool requireCoupling() const
    {
        return local_assembler.doBoundaryCoupling();
    }
    bool requireUVSkeleton() const
    {
        return local_assembler.doAlphaSkeleton();
    }
    bool requireUVBoundary() const
    {
        return local_assembler.doAlphaBoundary();
    }
    bool requireUVVolumePostSkeleton() const
    {
        return local_assembler.doAlphaVolumePostSkeleton();
    }
    //! @}

    //! Public access to the wrapping local assembler
    const LocalAssembler & localAssembler() const
    {
        return local_assembler;
    }

    //! Trial space constraints
    const typename LocalAssembler::Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
    {
        return localAssembler().trialConstraints();
    }

    //! Test space constraints
    const typename LocalAssembler::Traits::TestGridFunctionSpaceConstraints& testConstraints() const
    {
        return localAssembler().testConstraints();
    }

    //! Set current residual vector. Should be called prior to
    //! assembling.
    void setJacobian(Jacobian & jacobian_)
    {
        global_a_ss_view.attach(jacobian_);
        global_a_sn_view.attach(jacobian_);
        global_a_ns_view.attach(jacobian_);
        global_a_nn_view.attach(jacobian_);
    }

    //! Set current solution vector. Should be called prior to
    //! assembling.
    void setSolution(const Solution & solution_)
    {
        global_s_s_view.attach(solution_);
        global_s_s_view_pair.attach(solution_);

        global_s_n_view.attach(solution_);
    }

    void setSolutionOld(const Solution & solutionOld_)
    {
        global_s_s_view_old.attach(solutionOld_);
        global_s_n_view_old.attach(solutionOld_);
        global_s_s_view_pair_old.attach(solutionOld_);
    }

    void setSolutionSl(const SolutionSl & solutionSl_)
    {
        global_s_s_viewSl.attach(solutionSl_);
        global_s_n_viewSl.attach(solutionSl_);
    }

    //! Called immediately after binding of local function space in
    //! global assembler.
    //! @{
    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache,
                     const LFSVC & lfsv_cache)
    {
        global_s_s_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());

        global_a_ss_view.bind(lfsv_cache, lfsu_cache);
        al.assign(lfsv_cache.size(), lfsu_cache.size(), 0.0);
    }

    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVOld(const EG & eg, const LFSUC & lfsu_cache,
                        const LFSVC & lfsv_cache)
    {
        global_s_s_view_old.bind(lfsu_cache);
        xlOld.resize(lfsu_cache.size());
    }

    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVPair(const EG & eg, const LFSUC & lfsu_cache,
                         const LFSVC & lfsv_cache)
    {
        global_s_s_view_pair.bind(lfsu_cache);
        xlPair.resize(lfsu_cache.size());
    }

        template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVPairOld(const EG & eg, const LFSUC & lfsu_cache,
                         const LFSVC & lfsv_cache)
    {
        global_s_s_view_pair_old.bind(lfsu_cache);
        xlPairOld.resize(lfsu_cache.size());
    }

    template<typename EGSL, typename LFSUCSL, typename LFSVCSL>
    void onBindLFSUVSL(const EGSL & egSl, const LFSUCSL & lfsu_cacheSl,
                       const LFSVCSL & lfsv_cacheSl)
    {
        global_s_s_viewSl.bind(lfsu_cacheSl);
        xlSl.resize(lfsu_cacheSl.size());
    }


    template<typename IG, typename LFSUC, typename LFSVC>
    void onBindLFSUVOutside(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache, const LFSUC & lfsu_n_cache,
                            const LFSVC & lfsv_n_cache)
    {
        global_s_n_view.bind(lfsu_n_cache);
        xn.resize(lfsu_n_cache.size());
        global_a_sn_view.bind(lfsv_s_cache, lfsu_n_cache);
        al_sn.assign(lfsv_s_cache.size(), lfsu_n_cache.size(), 0.0);
        global_a_ns_view.bind(lfsv_n_cache, lfsu_s_cache);
        al_ns.assign(lfsv_n_cache.size(), lfsu_s_cache.size(), 0.0);
        global_a_nn_view.bind(lfsv_n_cache, lfsu_n_cache);
        al_nn.assign(lfsv_n_cache.size(), lfsu_n_cache.size(), 0.0);
    }

    //! @}

    //! Called when the local function space is about to be rebound or
    //! discarded
    //! @{
    template<typename EG, typename LFSUC, typename LFSVC>
    void onUnbindLFSUV(const EG & eg, const LFSUC & lfsu_cache,
                       const LFSVC & lfsv_cache)
    {
        local_assembler.scatter_jacobian(al, global_a_ss_view, false);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void onUnbindLFSUVOutside(const IG & ig, const LFSUC & lfsu_s_cache,
                              const LFSVC & lfsv_s_cache, const LFSUC & lfsu_n_cache,
                              const LFSVC & lfsv_n_cache)
    {
        local_assembler.scatter_jacobian(al_sn, global_a_sn_view,
                                         false);
        local_assembler.scatter_jacobian(al_ns, global_a_ns_view,
                                         false);
        local_assembler.scatter_jacobian(al_nn, global_a_nn_view,
                                         false);
    }

    //! @}

    //! Methods for loading of the local function's coefficients
    //! @{
    template<typename LFSUC>
    void loadCoefficientsLFSUInside(const LFSUC & lfsu_cache)
    {
        global_s_s_view.read(xl);
    }

    template<typename LFSUC>
    void loadCoefficientsLFSUInsideOld(const LFSUC & lfsu_cacheOld)
    {
        global_s_s_view_old.read(xlOld);
    }
    template<typename LFSUC>
    void loadCoefficientsLFSUInsidePair(const LFSUC & lfsu_cachePair)
    {
        global_s_s_view_pair.read(xlPair);
    }
    template<typename LFSUC>
    void loadCoefficientsLFSUInsidePairOld(const LFSUC & lfsu_cachePair)
    {
        global_s_s_view_pair_old.read(xlPairOld);
    }
    template<typename LFSUCSL>
    void loadCoefficientsLFSUInsideSl(const LFSUCSL & lfsu_cacheSl)
    {
        global_s_s_viewSl.read(xlSl);
    }

    template<typename LFSUC>
    void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache)
    {
        global_s_n_view.read(xn);
    }
    template<typename LFSUC>
    void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "No coupling lfsu_cache available for ");
    }
    //! @}

    //! Notifier functions, called immediately before and after assembling
    //! @{
    void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
    {
        Jacobian& jacobian = global_a_ss_view.container();
        global_s_s_view.detach();
        global_s_n_view.detach();
        global_a_ss_view.detach();
        global_a_sn_view.detach();
        global_a_ns_view.detach();
        global_a_nn_view.detach();

        if (local_assembler.doPostProcessing)
        {
            local_assembler.handle_dirichlet_constraints(gfsv,
                                                         jacobian);
        }
    }
    //! @}

    //! Assembling methods
    //! @{

    /** Assemble on a given cell without function spaces.

     \return If true, the assembling for this cell is assumed to
     be complete and the assembler continues with the next grid
     cell.
     */
    template<typename EG>
    bool assembleCell(const EG & eg)
    {
        return LocalAssembler::isNonOverlapping
               && eg.entity().partitionType() != Dune::InteriorEntity;
    }

    template<typename EG, typename LFSUC, typename LFSVC>
    void assembleUVVolume(const EG & eg, const LFSUC & lfsu_cache,
                          const LFSVC & lfsv_cache, int iMat)
    {
        al_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaVolume>::jacobian_volume(
            lop, eg, lfsu_cache.localFunctionSpace(), xl,
            lfsv_cache.localFunctionSpace(), iMat, al_view);
    }

    template<typename IG, typename EGC, typename EGP, typename LFSU, typename LFSV, typename LFSUC, typename LFSVC>
    void assembleCBoundary(const IG& ig, const EGC& egC, const EGP& egP,
                           const LFSU& lfsu_cache, const LFSV& lfsv_cache, const LFSU& lfsu_cache_pair, const LFSV& lfsv_cache_pair, const LFSUC& lfsuC_cache,
                           const LFSVC& lfsvC_cache)
    {
        al_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitchRCoupling<LOP,
             LOP::doBoundaryCoupling>::jacobian_boundary_coupling(lop, ig, egC, lfsu_cache.localFunctionSpace(),
                                                                  xl, xlOld, lfsv_cache.localFunctionSpace(),
                                                                  lfsuC_cache.localFunctionSpace(), xlSl,
                                                                  lfsvC_cache.localFunctionSpace(), al_view) ;
    }

    template<typename EG, typename IGC, typename LFSU, typename LFSV, typename LFSUC, typename LFSVC, typename FH>
    void assembleCBoundaryReversed(const EG& eg, const IGC& igC, const LFSU& lfsu_cache,
                                   const LFSV& lfsv_cache, const LFSUC& lfsuC_cache,
                                   const LFSVC& lfsvC_cache, std::vector<FH> fractureHeight)
    {
        al_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitchRCoupling<LOP,
             LOP::doBoundaryCouplingReversed>::jacobian_boundary_coupling_reversed(lop, eg, igC, lfsu_cache.localFunctionSpace(),
                                                                                   xl, lfsv_cache.localFunctionSpace(),
                                                                                   lfsuC_cache.localFunctionSpace(), xlSl,
                                                                                   lfsvC_cache.localFunctionSpace(), fractureHeight, al_view) ;
    }


    template<typename IG, typename LFSUC, typename LFSVC>
    void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache, const LFSUC & lfsu_n_cache,
                            const LFSVC & lfsv_n_cache)
    {
        al_view.setWeight(local_assembler.weight);
        al_sn_view.setWeight(local_assembler.weight);
        al_ns_view.setWeight(local_assembler.weight);
        al_nn_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaSkeleton>::jacobian_skeleton(
            lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
            lfsv_s_cache.localFunctionSpace(),
            lfsu_n_cache.localFunctionSpace(), xn,
            lfsv_n_cache.localFunctionSpace(), al_view, al_sn_view,
            al_ns_view, al_nn_view);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache)
    {
        al_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaBoundary>::jacobian_boundary(
            lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
            lfsv_s_cache.localFunctionSpace(), al_view);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    static void assembleUVEnrichedCoupling(const IG & ig,
                                           const LFSUC & lfsu_s_cache, const LFSVC & lfsv_s_cache,
                                           const LFSUC & lfsu_n_cache, const LFSVC & lfsv_n_cache,
                                           const LFSUC & lfsu_coupling_cache,
                                           const LFSVC & lfsv_coupling_cache)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Assembling of coupling spaces is not implemented for ");
    }

    template<typename IG, typename LFSVC>
    static void assembleVEnrichedCoupling(const IG & ig,
                                          const LFSVC & lfsv_s_cache, const LFSVC & lfsv_n_cache,
                                          const LFSVC & lfsv_coupling_cache)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Assembling of coupling spaces is not implemented for ");
    }

    template<typename EG, typename LFSUC, typename LFSVC>
    void assembleUVVolumePostSkeleton(const EG & eg,
                                      const LFSUC & lfsu_cache, const LFSVC & lfsv_cache)
    {
        al_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,
             LOP::doAlphaVolumePostSkeleton>::jacobian_volume_post_skeleton(
                 lop, eg, lfsu_cache.localFunctionSpace(), xl,
                 lfsv_cache.localFunctionSpace(), al_view);
    }

//! @}

private:
//! Reference to the wrapping local assembler object which
//! constructed this engine
    const LocalAssembler & local_assembler;

//! Reference to the local operator
    const LOP & lop;

//! Pointer to the current solution vector for which to assemble
    SolutionView global_s_s_view;
    SolutionView global_s_n_view;

    SolutionView global_s_s_view_pair;
    SolutionView global_s_s_view_pair_old;


    SolutionView global_s_s_view_old;
    SolutionView global_s_n_view_old;

    SolutionViewSl global_s_s_viewSl;
    SolutionViewSl global_s_n_viewSl;

//! Pointer to the current residual vector in which to assemble
    JacobianView global_a_ss_view;
    JacobianView global_a_sn_view;
    JacobianView global_a_ns_view;
    JacobianView global_a_nn_view;

    JacobianViewSl global_a_ss_viewSl;
    JacobianViewSl global_a_sn_viewSl;
    JacobianViewSl global_a_ns_viewSl;
    JacobianViewSl global_a_nn_viewSl;

//! The local vectors and matrices as required for assembling
//! @{
    typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
    typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;

    typedef Dune::PDELab::LocalVector<SolutionElement,
            LocalTrialSpaceTag> SolutionVector;
    typedef typename std::conditional <
    std::is_base_of<lop::DiagonalJacobian, LOP>::value,
        Dune::PDELab::DiagonalLocalMatrix<JacobianElement>,
        Dune::PDELab::LocalMatrix<JacobianElement> >::type JacobianMatrix;

    typedef Dune::PDELab::LocalVector<SolutionElementSl,
            LocalTrialSpaceTag> SolutionVectorSl;

    SolutionVector xl;
    SolutionVector xlPair;
    SolutionVector xlPairOld;
    SolutionVector xn;

    SolutionVector xlOld;
    SolutionVector xnOld;


    SolutionVectorSl xlSl;
    SolutionVectorSl xnSl;

    JacobianMatrix al;
    JacobianMatrix al_sn;
    JacobianMatrix al_ns;
    JacobianMatrix al_nn;

    typename JacobianMatrix::WeightedAccumulationView al_view;
    typename JacobianMatrix::WeightedAccumulationView al_sn_view;
    typename JacobianMatrix::WeightedAccumulationView al_ns_view;
    typename JacobianMatrix::WeightedAccumulationView al_nn_view;


//! @}

};
// End of class DefaultLocalJacobianAssemblerEngine

}
}
#endif
