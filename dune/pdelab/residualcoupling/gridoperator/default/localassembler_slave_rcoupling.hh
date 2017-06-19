#ifndef DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_SLAVE_HH
#define DUNE_PDELAB_DEFAULT_LOCAL_ASSEMBLER_SLAVE_HH

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/residualcoupling/gridoperator/default/residualengine_normal_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/patternengine_normal_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/jacobianengine_slave_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/jacobianapplyengine_rcoupling.hh>
//#include <dune/pdelab/residualcoupling/gridoperator/default/nonlinearjacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune
{
namespace PDELab
{

/**
 \brief The local assembler for DUNE grids

 \tparam GFSU GridFunctionSpace for ansatz functions
 \tparam GFSV GridFunctionSpace for test functions
 \tparam X The solution vector representation type
 \tparam R The residual vector representation type
 \tparam A The jacobian matrix representation type
 \tparam B The matrix backend
 \tparam P The matrix pattern representation type
 \tparam CU   Constraints maps for the individual dofs (trial space)
 \tparam CV   Constraints maps for the individual dofs (test space)

 */
template<typename GO, typename LOP, bool nonoverlapping_mode = false>
class DefaultLocalAssemblerSlave: public Dune::PDELab::LocalAssemblerBase <
    typename GO::Traits::MatrixBackendSl,
    typename GO::Traits::TrialGridFunctionSpaceConstraintsSl,
    typename GO::Traits::TestGridFunctionSpaceConstraintsSl >
{

    // The GridOperator has to be a friend to modify the do{Pre,Post}Processing flags
    template<typename, typename, typename, typename, typename, typename,
             typename, typename, typename, typename, typename, typename,
             typename, typename, typename, typename, typename, typename,
             typename, typename, typename, int, bool>
    friend class GridOperatorCoupling;
    template<typename, typename, typename, typename, typename, typename,
             typename, typename, typename, int>
    friend class GridOperatorNormalRCoupling;
public:

    //! The traits class
    typedef Dune::PDELab::LocalAssemblerTraitsCoupling<GO> Traits;

    //! The local operators type for real numbers e.g. time
    typedef typename Traits::ResidualSl::ElementType RangeField;
    typedef RangeField Real;

    typedef typename Traits::TrialGridFunctionSpaceSl GFSU;
    typedef typename Traits::TestGridFunctionSpaceSl GFSV;

    typedef typename Traits::TrialGridFunctionSpaceConstraintsSl CU;
    typedef typename Traits::TestGridFunctionSpaceConstraintsSl CV;

    //! The base class of this local assembler
    typedef Dune::PDELab::LocalAssemblerBase <
    typename Traits::MatrixBackendSl, CU, CV > Base;

    //! The local operator
    typedef LOP LocalOperator;

    static const bool isNonOverlapping = nonoverlapping_mode;

    //! The local function spaces
    //! @{
    // Types of local function spaces
    typedef Dune::PDELab::LocalFunctionSpace<GFSU,
            Dune::PDELab::TrialSpaceTag> LFSU;
    typedef Dune::PDELab::LocalFunctionSpace<GFSV,
            Dune::PDELab::TestSpaceTag> LFSV;
    typedef LFSIndexCache<LFSU, CU> LFSUCache;
    typedef LFSIndexCache<LFSV, CV> LFSVCache;

    //! @}

    //! The local assembler engines
    //! @{
    typedef DefaultLocalPatternAssemblerEngineNormalRCoupling <
    DefaultLocalAssemblerSlave > LocalPatternAssemblerEngine;
    typedef DefaultLocalResidualAssemblerEngineNormalRCoupling <
    DefaultLocalAssemblerSlave > LocalResidualAssemblerEngine;
    typedef DefaultLocalJacobianAssemblerEngineSlave <
    DefaultLocalAssemblerSlave > LocalJacobianAssemblerEngine;
    typedef DefaultLocalJacobianApplyAssemblerEngineRCoupling <
    DefaultLocalAssemblerSlave > LocalJacobianApplyAssemblerEngine;
//          typedef DefaultLocalNonlinearJacobianApplyAssemblerEngine<
//                  DefaultLocalAssemblerSlave> LocalNonlinearJacobianApplyAssemblerEngine;

    friend class DefaultLocalPatternAssemblerEngineNormalRCoupling <
        DefaultLocalAssemblerSlave > ;
    friend class DefaultLocalResidualAssemblerEngineNormalRCoupling <
        DefaultLocalAssemblerSlave > ;
    friend class DefaultLocalJacobianAssemblerEngineSlave <
        DefaultLocalAssemblerSlave > ;
    friend class DefaultLocalJacobianApplyAssemblerEngineRCoupling <
        DefaultLocalAssemblerSlave > ;
//          friend class DefaultLocalNonlinearJacobianApplyAssemblerEngine<
//                  DefaultLocalAssemblerSlave> ;
    //! @}

    //! Constructor with empty constraints
    DefaultLocalAssemblerSlave(LOP & lop_,
                               shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger) :
        lop(lop_), weight(1.0), doPreProcessing(true), doPostProcessing(
            true), pattern_engine(*this, border_dof_exchanger), residual_engine(
                *this), jacobian_engine(*this), jacobian_apply_engine(
                    *this), _reconstruct_border_entries(
                        isNonOverlapping)
    {
    }

    //! Constructor for non trivial constraints
    DefaultLocalAssemblerSlave(LOP & lop_, const CU& cu_, const CV& cv_,
                               shared_ptr<typename GO::BorderDOFExchanger> border_dof_exchanger) :
        Base(cu_, cv_), lop(lop_), weight(1.0), doPreProcessing(
            true), doPostProcessing(true), pattern_engine(*this,
                                                          border_dof_exchanger), residual_engine(*this), jacobian_engine(
                                                              *this), jacobian_apply_engine(*this), _reconstruct_border_entries(
                                                                  isNonOverlapping)
    {
    }

    //! Notifies the local assembler about the current time of
    //! assembling. Should be called before assembling if the local
    //! operator has time dependencies.
    void setTime(Real time_)
    {
        lop.setTime(time_);
    }

    //! Notifies the assembler about the current weight of assembling.
    void setWeight(RangeField weight_)
    {
        weight = weight_;
    }

    //! Time stepping interface
    //! @{
    void preStage(Real time_, int r_)
    {
        lop.preStage(time_, r_);
    }
    void preStep(Real time_, Real dt_, std::size_t stages_)
    {
        lop.preStep(time_, dt_, stages_);
    }
    void postStep()
    {
        lop.postStep();
    }
    void postStage()
    {
        lop.postStage();
    }
    Real suggestTimestep(Real dt) const
    {
        return lop.suggestTimestep(dt);
    }
    //! @}

    bool reconstructBorderEntries() const
    {
        return _reconstruct_border_entries;
    }

    //! Access methods which provid "ready to use" engines
    //! @{

    //! Returns a reference to the requested engine. This engine is
    //! completely configured and ready to use.
    LocalPatternAssemblerEngine & localPatternAssemblerEngine(
        typename Traits::MatrixPattern & p)
    {
        pattern_engine.setPattern(p);
        return pattern_engine;
    }

    //! Returns a reference to the requested engine. This engine is
    //! completely configured and ready to use.
    LocalResidualAssemblerEngine & localResidualAssemblerEngine(
        typename Traits::Residual & r,
        const typename Traits::Solution & x)
    {
        residual_engine.setResidual(r);
        residual_engine.setSolution(x);
        return residual_engine;
    }

    //! Returns a reference to the requested engine. This engine is
    //! completely configured and ready to use.
    LocalJacobianAssemblerEngine & localJacobianAssemblerEngine(
        typename Traits::Jacobian & a,
        const typename Traits::Solution & x)
    {
        jacobian_engine.setJacobian(a);
        jacobian_engine.setSolution(x);
        return jacobian_engine;
    }

    //! Returns a reference to the requested engine. This engine is
    //! completely configured and ready to use.
    LocalJacobianApplyAssemblerEngine & localJacobianApplyAssemblerEngine(
        typename Traits::Residual & r,
        const typename Traits::Solution & x)
    {
        jacobian_apply_engine.setResidual(r);
        jacobian_apply_engine.setSolution(x);
        return jacobian_apply_engine;
    }

    //! Returns a reference to the requested engine. This engine is
    //! completely configured and ready to use.
//          LocalNonlinearJacobianApplyAssemblerEngine & localNonlinearJacobianApplyAssemblerEngine(
//                  typename Traits::Residual & r,
//                  const typename Traits::Solution & x,
//                  const typename Traits::Solution & z)
//          {
//              nonlinear_jacobian_apply_engine.setResidual(r);
//              nonlinear_jacobian_apply_engine.setSolution(x);
//              nonlinear_jacobian_apply_engine.setUpdate(z);
//              return nonlinear_jacobian_apply_engine;
//          }

    //! @}

    //! \brief Query methods for the assembler engines. Theses methods
    //! do not belong to the assembler interface, but simplify the
    //! implementations of query methods in the engines;
    //! @{
    static bool doAlphaVolume()
    {
        return LOP::doAlphaVolume;
    }

    static bool doBoundaryCoupling()
    {
        return LOP::doBoundaryCoupling;
    }

    static bool doBoundaryCouplingReversed()
    {
        return LOP::doBoundaryCouplingReversed;
    }
    static bool doLambdaVolume()
    {
        return LOP::doLambdaVolume;
    }
    static bool doAlphaSkeleton()
    {
        return LOP::doAlphaSkeleton;
    }
    static bool doLambdaSkeleton()
    {
        return LOP::doLambdaSkeleton;
    }
    static bool doAlphaBoundary()
    {
        return LOP::doAlphaBoundary;
    }
    static bool doLambdaBoundary()
    {
        return LOP::doLambdaBoundary;
    }
    static bool doAlphaVolumePostSkeleton()
    {
        return LOP::doAlphaVolumePostSkeleton;
    }
    static bool doLambdaVolumePostSkeleton()
    {
        return LOP::doLambdaVolumePostSkeleton;
    }
    static bool doSkeletonTwoSided()
    {
        return LOP::doSkeletonTwoSided;
    }
    static bool doPatternVolume()
    {
        return LOP::doPatternVolume;
    }

    static bool doPatternCoupling()
    {
        return LOP::doPatternCoupling;
    }
    static bool doPatternSkeleton()
    {
        return LOP::doPatternSkeleton;
    }
    static bool doPatternBoundary()
    {
        return LOP::doPatternBoundary;
    }
    static bool doPatternVolumePostSkeleton()
    {
        return LOP::doPatternVolumePostSkeleton;
    }
    //! @}

    //! This method allows to set the behavior with regard to any
    //! preprocessing within the engines. It is called by the
    //! setupGridOperators() method of the GridOperator and should
    //! not be called directly.
    void preProcessing(bool v)
    {
        doPreProcessing = v;
    }

    //! This method allows to set the behavior with regard to any
    //! postprocessing within the engines. It is called by the
    //! setupGridOperators() method of the GridOperator and should
    //! not be called directly.
    void postProcessing(bool v)
    {
        doPostProcessing = v;
    }

private:

    //! The local operator
    LOP & lop;

    //! The current weight of assembling
    RangeField weight;

    //! Indicates whether this local operator has to perform pre
    //! processing
    bool doPreProcessing;

    //! Indicates whether this local operator has to perform post
    //! processing
    bool doPostProcessing;

    //! The engine member objects
    //! @{
    LocalPatternAssemblerEngine pattern_engine;
    LocalResidualAssemblerEngine residual_engine;
    LocalJacobianAssemblerEngine jacobian_engine;
    LocalJacobianApplyAssemblerEngine jacobian_apply_engine;
//          LocalNonlinearJacobianApplyAssemblerEngine nonlinear_jacobian_apply_engine;
    //! @}

    bool _reconstruct_border_entries;

};

}
}
#endif
