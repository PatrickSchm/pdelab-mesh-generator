#ifndef DUNE_PDELAB_DEFAULT_RESIDUALENGINE_COUPLING_RCOUPLING_HH
#define DUNE_PDELAB_DEFAULT_RESIDUALENGINE_COUPLING_RCOUPLING_HH

#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/residualcoupling/gridoperator/common/localassemblerenginebase_rcoupling.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/residualcoupling/localoperator/callswitch_rcoupling.hh>

namespace Dune
{
namespace PDELab
{

/**
 \brief The local assembler engine for DUNE grids which
 assembles the residual vector

 \tparam LA The local assembler

 */
template<typename LA, typename LASL>
class DefaultLocalResidualAssemblerEngineRCoupling: public LocalAssemblerEngineBaseRCoupling
{
public:

    template<typename TrialConstraintsContainer,
             typename TestConstraintsContainer>
    bool needsConstraintsCaching(const TrialConstraintsContainer& cu,
                                 const TestConstraintsContainer& cv) const
    {
        return false;
    }

    //! The type of the wrapping local assembler
    typedef LA LocalAssembler;

    //! The type of the local operator
    typedef typename LA::LocalOperator LOP;
//            typedef typename LASL::LocalOperator LOPSL;

//! The type of the residual vector
    typedef typename LA::Traits::Residual Residual;
    typedef typename Residual::ElementType ResidualElement;

    //! The type of the solution vector
    typedef typename LA::Traits::Solution Solution;
    typedef typename Solution::ElementType SolutionElement;

    typedef typename LA::Traits::SolutionSl SolutionSl;
    typedef typename SolutionSl::ElementType SolutionElementSl;

    //! The local function spaces
    typedef typename LA::LFSU LFSU;
    typedef typename LA::LFSUCache LFSUCache;
    typedef typename LFSU::Traits::GridFunctionSpace GFSU;
    typedef typename LA::LFSV LFSV;
    typedef typename LA::LFSVCache LFSVCache;
    typedef typename LFSV::Traits::GridFunctionSpace GFSV;

    typedef typename LASL::LFSU LFSUSL;
    typedef typename LASL::LFSUCache LFSUCacheSL;
    typedef typename LFSUSL::Traits::GridFunctionSpace GFSUSL;
    typedef typename LASL::LFSV LFSVSL;
    typedef typename LASL::LFSVCache LFSVCacheSL;
    typedef typename LFSVSL::Traits::GridFunctionSpace GFSVSL;

    typedef typename Solution::template ConstLocalView<LFSUCache> SolutionView;
    typedef typename SolutionSl::template ConstLocalView<LFSUCacheSL> SolutionViewSl;
    typedef typename Residual::template LocalView<LFSVCache> ResidualView;

    /**
     \brief Constructor

     \param [in] local_assembler_ The local assembler object which
     creates this engine
     */
    DefaultLocalResidualAssemblerEngineRCoupling(
        const LocalAssembler & local_assembler_) :
        local_assembler(local_assembler_), lop(
            local_assembler_.lop), rl_view(rl, 1.0), rn_view(rn,
                                                             1.0)
    {
    }

    //! Query methods for the global grid assembler
    //! @{
    bool requireSkeleton() const
    {
        return (local_assembler.doAlphaSkeleton()
                || local_assembler.doLambdaSkeleton());
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
    bool requireCouplingReversed() const
    {
        return local_assembler.doBoundaryCouplingReversed();
    }
    bool requireVVolume() const
    {
        return local_assembler.doLambdaVolume();
    }
    bool requireUVSkeleton() const
    {
        return local_assembler.doAlphaSkeleton();
    }
    bool requireVSkeleton() const
    {
        return local_assembler.doLambdaSkeleton();
    }
    bool requireUVBoundary() const
    {
        return local_assembler.doAlphaBoundary();
    }
    bool requireVBoundary() const
    {
        return local_assembler.doLambdaBoundary();
    }
    bool requireUVVolumePostSkeleton() const
    {
        return local_assembler.doAlphaVolumePostSkeleton();
    }
    bool requireVVolumePostSkeleton() const
    {
        return local_assembler.doLambdaVolumePostSkeleton();
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
    void setResidual(Residual & residual_)
    {
        global_rl_view.attach(residual_);
        global_rn_view.attach(residual_);
    }

    //! Set current solution vector. Should be called prior to
    //! assembling.
    void setSolution(const Solution & solution_)
    {
        global_sl_view.attach(solution_);
        global_sn_view.attach(solution_);
        global_sl_view_pair.attach(solution_);
    }

    void setSolutionOld(const Solution & solutionOld_)
    {
        global_sl_view_old.attach(solutionOld_);
        global_sn_view_old.attach(solutionOld_);
        global_sl_view_pair_old.attach(solutionOld_);
    }

    void setSolutionSl(const SolutionSl & solutionSl_)
    {
        global_sl_viewSl.attach(solutionSl_);
        global_sn_view_sl.attach(solutionSl_);
    }

    //! Called immediately after binding of local function space in
    //! global assembler.
    //! @{
    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUV(const EG & eg, const LFSUC & lfsu_cache,
                     const LFSVC & lfsv_cache)
    {
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
    }

    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVPair(const EG & eg, const LFSUC & lfsu_cache,
                         const LFSVC & lfsv_cache)
    {
        global_sl_view_pair.bind(lfsu_cache);
        xlPair.resize(lfsu_cache.size());
    }
    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVPairOld(const EG & eg, const LFSUC & lfsu_cache,
                            const LFSVC & lfsv_cache)
    {
        global_sl_view_pair_old.bind(lfsu_cache);
        xlPairOld.resize(lfsu_cache.size());
    }
    template<typename EG, typename LFSUC, typename LFSVC>
    void onBindLFSUVOld(const EG & eg, const LFSUC & lfsu_cache,
                        const LFSVC & lfsv_cache)
    {
        global_sl_view_old.bind(lfsu_cache);
        xlOld.resize(lfsu_cache.size());
    }

    template<typename EGSL, typename LFSUCSL, typename LFSVCSL>
    void onBindLFSUVSL(const EGSL & egSl, const LFSUCSL & lfsu_cacheSl,
                       const LFSVCSL & lfsv_cacheSl)
    {
        global_sl_viewSl.bind(lfsu_cacheSl);
        xlSl.resize(lfsu_cacheSl.size());
    }

    template<typename EG, typename LFSVC>
    void onBindLFSV(const EG & eg, const LFSVC & lfsv_cache)
    {
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(), 0.0);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void onBindLFSUVInside(const IG & ig, const LFSUC & lfsu_cache,
                           const LFSVC & lfsv_cache)
    {
        global_sl_view.bind(lfsu_cache);
        xl.resize(lfsu_cache.size());
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void onBindLFSUVOutside(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache, const LFSUC & lfsu_n_cache,
                            const LFSVC & lfsv_n_cache)
    {
        global_sn_view.bind(lfsu_n_cache);
        xn.resize(lfsu_n_cache.size());
    }

    template<typename IG, typename LFSVC>
    void onBindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
    {
        global_rl_view.bind(lfsv_cache);
        rl.assign(lfsv_cache.size(), 0.0);
    }

    template<typename IG, typename LFSVC>
    void onBindLFSVOutside(const IG & ig, const LFSVC & lfsv_s_cache,
                           const LFSVC & lfsv_n_cache)
    {
        global_rn_view.bind(lfsv_n_cache);
        rn.assign(lfsv_n_cache.size(), 0.0);
    }

    //! @}

    //! Called when the local function space is about to be rebound or
    //! discarded
    //! @{
    template<typename EG, typename LFSVC>
    void onUnbindLFSV(const EG & eg, const LFSVC & lfsv_cache)
    {
        global_rl_view.add(rl);
        global_rl_view.commit();
    }

    template<typename IG, typename LFSVC>
    void onUnbindLFSVInside(const IG & ig, const LFSVC & lfsv_cache)
    {
        global_rl_view.add(rl);
        global_rl_view.commit();
    }

    template<typename IG, typename LFSVC>
    void onUnbindLFSVOutside(const IG & ig, const LFSVC & lfsv_s_cache,
                             const LFSVC & lfsv_n_cache)
    {
        global_rn_view.add(rn);
        global_rn_view.commit();
    }
    //! @}

    //! Methods for loading of the local function's coefficients
    //! @{
    template<typename LFSUC>
    void loadCoefficientsLFSUInside(const LFSUC & lfsu_s_cache)
    {
        global_sl_view.read(xl);
    }

    template<typename LFSUC>
    void loadCoefficientsLFSUInsideOld(const LFSUC & lfsu_s_cache)
    {
        global_sl_view_old.read(xlOld);
    }

    template<typename LFSUC>
    void loadCoefficientsLFSUInsidePair(const LFSUC & lfsu_s_cache_pair)
    {
        global_sl_view_pair.read(xlPair);
    }

    template<typename LFSUC>
    void loadCoefficientsLFSUInsidePairOld(const LFSUC & lfsu_s_cache_pair)
    {
        global_sl_view_pair_old.read(xlPairOld);
    }

    template<typename LFSUCSL>
    void loadCoefficientsLFSUInsideSl(const LFSUCSL & lfsu_s_cacheSl)
    {
        global_sl_viewSl.read(xlSl);
    }
    template<typename LFSUC>
    void loadCoefficientsLFSUOutside(const LFSUC & lfsu_n_cache)
    {
        global_sn_view.read(xn);
    }
    template<typename LFSUC>
    void loadCoefficientsLFSUCoupling(const LFSUC & lfsu_c_cache)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "No coupling lfsu available for ");
    }
    //! @}

    //! Notifier functions, called immediately before and after assembling
    //! @{

    void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
    {
        if (local_assembler.doPostProcessing)
        {
            Dune::PDELab::constrain_residual(
                *(local_assembler.pconstraintsv),
                global_rl_view.container());
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
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaVolume>::alpha_volume(
            lop, eg, lfsu_cache.localFunctionSpace(), xl,
            lfsv_cache.localFunctionSpace(), rl_view, iMat);
    }

    template<typename IG, typename EGC, typename EGP, typename LFSUC, typename LFSVC, typename LFSUCC, typename LFSVCC>
    void assembleCBoundary(const IG& ig, const EGC& egC, const EGP& egP,
                           const LFSUC& lfsu_s_cache, const LFSVC& lfsv_s_cache, const LFSUC& lfsu_pair_cache, const LFSVC& lfsv_pair_cache, const LFSUCC& lfsuC_cache,
                           const LFSVCC& lfsvC_cache)
    {
        rl_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitchRCoupling<LOP,
             LOP::doBoundaryCoupling>::boundary_coupling(lop, ig, egC, egP,
                                                         lfsu_s_cache.localFunctionSpace(), lfsu_pair_cache.localFunctionSpace(),
                                                         xl, xlPair, xlOld, xlPairOld, lfsv_s_cache.localFunctionSpace(), lfsv_pair_cache.localFunctionSpace(),
                                                         lfsuC_cache.localFunctionSpace(), xlSl,
                                                         lfsvC_cache.localFunctionSpace(), rl_view);
    }

    template<typename EG, typename IGC, typename LFSUC, typename LFSVC, typename LFSUCC, typename LFSVCC, typename FH>
    void assembleCBoundaryReversed(const EG& eg, const IGC& igC, const LFSUC& lfsu_cache,
                                   const LFSVC& lfsv_cache, const LFSUCC& lfsuC_cache,
                                   const LFSVCC& lfsvC_cache, std::vector<FH>& fractureHeight)
    {
        rl_view.setWeight(local_assembler.weight);

        Dune::PDELab::LocalAssemblerCallSwitchRCoupling<LOP,
             LOP::doBoundaryCouplingReversed>::boundary_coupling_reversed(lop, eg, igC,
                                                                          lfsu_cache.localFunctionSpace(), xl,
                                                                          lfsv_cache.localFunctionSpace(), rl_view,
                                                                          lfsuC_cache.localFunctionSpace(), xlSl,
                                                                          lfsvC_cache.localFunctionSpace(), fractureHeight);
    }

    template<typename EG, typename LFSVC>
    void assembleVVolume(const EG & eg, const LFSVC & lfsv_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doLambdaVolume>::lambda_volume(
            lop, eg, lfsv_cache.localFunctionSpace(), rl_view);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void assembleUVSkeleton(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache, const LFSUC & lfsu_n_cache,
                            const LFSVC & lfsv_n_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaSkeleton>::alpha_skeleton(
            lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
            lfsv_s_cache.localFunctionSpace(),
            lfsu_n_cache.localFunctionSpace(), xn,
            lfsv_n_cache.localFunctionSpace(), rl_view, rn_view);
    }

    template<typename IG, typename LFSVC>
    void assembleVSkeleton(const IG & ig, const LFSVC & lfsv_s_cache,
                           const LFSVC & lfsv_n_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        rn_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,
             LOP::doLambdaSkeleton>::lambda_skeleton(lop, ig,
                                                     lfsv_s_cache.localFunctionSpace(),
                                                     lfsv_n_cache.localFunctionSpace(), rl_view, rn_view);
    }

    template<typename IG, typename LFSUC, typename LFSVC>
    void assembleUVBoundary(const IG & ig, const LFSUC & lfsu_s_cache,
                            const LFSVC & lfsv_s_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP, LOP::doAlphaBoundary>::alpha_boundary(
            lop, ig, lfsu_s_cache.localFunctionSpace(), xl,
            lfsv_s_cache.localFunctionSpace(), rl_view);
    }

    template<typename IG, typename LFSVC>
    void assembleVBoundary(const IG & ig, const LFSVC & lfsv_s_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,
             LOP::doLambdaBoundary>::lambda_boundary(lop, ig,
                                                     lfsv_s_cache.localFunctionSpace(), rl_view);
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
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,
             LOP::doAlphaVolumePostSkeleton>::alpha_volume_post_skeleton(
                 lop, eg, lfsu_cache.localFunctionSpace(), xl,
                 lfsv_cache.localFunctionSpace(), rl_view);
    }

    template<typename EG, typename LFSVC>
    void assembleVVolumePostSkeleton(const EG & eg,
                                     const LFSVC & lfsv_cache)
    {
        rl_view.setWeight(local_assembler.weight);
        Dune::PDELab::LocalAssemblerCallSwitch<LOP,
             LOP::doLambdaVolumePostSkeleton>::lambda_volume_post_skeleton(
                 lop, eg, lfsv_cache.localFunctionSpace(), rl_view);
    }

    //! @}

private:
    //! Reference to the wrapping local assembler object which
    //! constructed this engine
    const LocalAssembler & local_assembler;

    //! Reference to the local operator
    const LOP & lop;
//            const LOPSL & lopSl;

    //! Pointer to the current residual vector in which to assemble
    ResidualView global_rl_view;
    ResidualView global_rn_view;

    //! Pointer to the current residual vector in which to assemble
    SolutionView global_sl_view;
    SolutionView global_sl_view_pair;
    SolutionView global_sl_view_pair_old;
    SolutionView global_sn_view;
    SolutionView global_sl_view_old;
    SolutionView global_sn_view_old;
    SolutionViewSl global_sl_viewSl;
    SolutionViewSl global_sn_view_sl;

    //! The local vectors and matrices as required for assembling
    //! @{
    typedef Dune::PDELab::TrialSpaceTag LocalTrialSpaceTag;
    typedef Dune::PDELab::TestSpaceTag LocalTestSpaceTag;

    typedef Dune::PDELab::LocalVector<SolutionElement,
            LocalTrialSpaceTag> SolutionVector;

    typedef Dune::PDELab::LocalVector<SolutionElementSl,
            LocalTrialSpaceTag> SolutionVectorSl;
    typedef Dune::PDELab::LocalVector<ResidualElement, LocalTestSpaceTag> ResidualVector;

    //! Inside local coefficients
    SolutionVector xl;
    SolutionVector xlPair;
    SolutionVector xlPairOld;
    SolutionVector xlOld;
    SolutionVectorSl xlSl;
    //! Outside local coefficients
    SolutionVector xn;
    //! Inside local residual
    ResidualVector rl;
    //! Outside local residual
    ResidualVector rn;
    //! Inside local residual weighted view
    typename ResidualVector::WeightedAccumulationView rl_view;
    //! Outside local residual weighted view
    typename ResidualVector::WeightedAccumulationView rn_view;
    //! @}

};
// End of class DefaultLocalResidualAssemblerEngine

}
}
#endif
