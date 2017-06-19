#ifndef DUNE_PDELAB_GRIDOPERATOR_COUPLING_RCOUPLING_HH
#define DUNE_PDELAB_GRIDOPERATOR_COUPLING_RCOUPLING_HH

#include <dune/common/tupleutility.hh>

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/residualcoupling/gridoperator/common/gridoperatorutilities_coupling_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/assembler_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/localassembler_rcoupling.hh>
#include <dune/pdelab/residualcoupling/gridoperator/default/localassembler_slave_rcoupling.hh>

namespace Dune
{
namespace PDELab
{

/**
 \brief Standard grid operator implementation

 \tparam GFSU GridFunctionSpace for ansatz functions
 \tparam GFSV GridFunctionSpace for test functions
 \tparam MB The matrix backend to be used for representation of the jacobian
 \tparam DF The domain field type of the operator
 \tparam RF The range field type of the operator
 \tparam JF The jacobian field type
 \tparam CU   Constraints maps for the individual dofs (trial space)
 \tparam CV   Constraints maps for the individual dofs (test space)

 */
template < typename GFSUMA, typename GFSVMA, typename GFSUSL,
           typename GFSVSL, typename FP, typename MaterialFlags, typename InPara,
           typename LOP, typename LOPSL, typename MB,
           typename MBSL, typename DF, typename DFSL, typename RF,
           typename RFSL, typename JF, typename JFSL,
           typename CU = Dune::PDELab::EmptyTransformation,
           typename CV = Dune::PDELab::EmptyTransformation,
           typename CUSL = Dune::PDELab::EmptyTransformation,
           typename CVSL = Dune::PDELab::EmptyTransformation,
           int nonoverlapping_mode = -1, bool coupling_assembly = true >
class GridOperatorCoupling: public impl::warn_on_deprecated_nonoverlapping_mode_parameter <
    nonoverlapping_mode >
{
public:

    static_assert(nonoverlapping_mode == -1 ||
                  nonoverlapping_mode == 0 ||
                  nonoverlapping_mode == 1,
                  "invalid value for nonoverlapping_mode! This parameter is also deprecated in PDELab 2.4, so please remove it from your typedefs!");

    //! The global assembler type
    typedef DefaultAssemblerCoupling<GFSUMA, GFSVMA, GFSUSL, GFSVSL, CU,
            CV, CUSL, CVSL, FP, MaterialFlags, InPara> Assembler;

    //! The type of the domain (solution).
    using Domain = Dune::PDELab::Backend::Vector<GFSUMA, DF>;
    //! The type of the range (residual).
    using DomainSl = Dune::PDELab::Backend::Vector<GFSUSL, DFSL>;
    //! The type of the range (residual).
    using Range = Dune::PDELab::Backend::Vector<GFSVMA, RF>;
    //! The type of the jacobian.
    using RangeSl = Dune::PDELab::Backend::Vector<GFSVSL, RF>;
    //! The type of the jacobian.
    using Jacobian = Dune::PDELab::Backend::Matrix<MB, Domain, Range, JF>;

    //! The sparsity pattern container for the jacobian matrix
    typedef typename MB::template Pattern<Jacobian, GFSVMA, GFSUMA> Pattern;

    //! The local assembler type
    typedef DefaultLocalAssemblerSlave<GridOperatorCoupling, LOPSL,
            GFSUSL::Traits::EntitySet::Partitions::partitionIterator()
            == InteriorBorder_Partition> LocalAssemblerSlave;

    typedef DefaultLocalAssemblerCoupling<GridOperatorCoupling,
            LocalAssemblerSlave, LOP,
            GFSUMA::Traits::EntitySet::Partitions::partitionIterator()
            == InteriorBorder_Partition> LocalAssembler;

    // Fix this as soon as the default Partitions are constexpr
    typedef typename std::conditional <
    GFSUMA::Traits::EntitySet::Partitions::partitionIterator()
    == InteriorBorder_Partition,
    NonOverlappingBorderDOFExchanger<GridOperatorCoupling>,
    OverlappingBorderDOFExchanger<GridOperatorCoupling> >::type BorderDOFExchanger;

    typedef typename std::conditional <
    GFSUSL::Traits::EntitySet::Partitions::partitionIterator()
    == InteriorBorder_Partition,
    NonOverlappingBorderDOFExchanger<GridOperatorCoupling>,
    OverlappingBorderDOFExchanger<GridOperatorCoupling> >::type BorderDOFExchangerSl;

    //! The grid operator traits
    typedef Dune::PDELab::GridOperatorTraitsRCoupling<GFSUMA, GFSVMA,
            GFSUSL, GFSVSL, MB, MBSL, DF, DFSL, RF, RFSL, JF, JFSL, CU,
            CV, CUSL, CVSL, Assembler, LocalAssembler,
            LocalAssemblerSlave> Traits;

    template<typename MFT>
    struct MatrixContainer
    {
        typedef typename Traits::Jacobian Type;
    };

    //! Constructor for non trivial constraints
    GridOperatorCoupling(const GFSUMA & gfsuMa_, const GFSUSL & gfsuSl_,
                         const CU & cu_, const CUSL & cuSl_, const GFSVMA & gfsvMa_,
                         const GFSVSL & gfsvSl_, const CV & cv_, const CVSL & cvSl_,
                         const FP &couplingPairs_,
                         const MaterialFlags & materialFlags_, const InPara & indicesParallel_, LOP & lop_, LOPSL & lopSl_,
                         const MB& mb_ = MB()) :
        global_assembler(gfsuMa_, gfsvMa_, gfsuSl_, gfsvSl_, couplingPairs_,
                         materialFlags_, indicesParallel_, cu_, cv_, cuSl_, cvSl_), dof_exchanger(
                             std::make_shared < BorderDOFExchanger > (*this)), local_assembler(
                                 lop_, cu_, cv_, dof_exchanger), backend(mb_)
    {
    }

    //! Constructor for empty constraints
    GridOperatorCoupling(const GFSUMA & gfsuMa_, const GFSVMA & gfsvMa_,
                         const GFSUSL & gfsuSl_, const GFSVSL & gfsvSl_,
                         const FP &couplingPairs_,
                         const MaterialFlags & materialFlags_, const InPara & indicesParallel_, LOP & lop_, LOPSL & lopSl_,
                         const MB& mb_ = MB()) :
        global_assembler(gfsuMa_, gfsvMa_, gfsuSl_, gfsvSl_, couplingPairs_,
                         materialFlags_, indicesParallel_), dof_exchanger(
                             std::make_shared < BorderDOFExchanger > (*this)), local_assembler(
                                 lop_, dof_exchanger), backend(mb_)
    {
    }

    //! Get the trial grid function space
    const GFSUMA& trialGridFunctionSpace() const
    {
        return global_assembler.trialGridFunctionSpaceMa();
    }

    //! Get the trial grid function space
    const GFSUSL& trialGridFunctionSpaceSl() const
    {
        return global_assembler.trialGridFunctionSpaceSl();
    }

    //! Get the test grid function space
    const GFSVMA& testGridFunctionSpace() const
    {
        return global_assembler.testGridFunctionSpaceMa();
    }

    //! Get the test grid function space
    const GFSVSL& testGridFunctionSpaceSl() const
    {
        return global_assembler.testGridFunctionSpaceSl();
    }

    //! Get dimension of space u
    typename GFSUMA::Traits::SizeType globalSizeU() const
    {
        return trialGridFunctionSpace().globalSize();
    }

    //! Get dimension of space v
    typename GFSVMA::Traits::SizeType globalSizeV() const
    {
        return testGridFunctionSpace().globalSize();
    }

    Assembler & assembler()
    {
        return global_assembler;
    }

    const Assembler & assembler() const
    {
        return global_assembler;
    }

    LocalAssembler & localAssembler() const
    {
        return local_assembler;
    }

//! Visitor which is called in the method setupGridOperators for
//! each tuple element.
    template<typename GridOperatorTuple>
    struct SetupGridOperator
    {
        SetupGridOperator() :
            index(0), size(
                Dune::tuple_size < GridOperatorTuple > ::value)
        {
        }

        template<typename T>
        void visit(T& elem)
        {
            elem.localAssembler().doPreProcessing = index == 0;
            elem.localAssembler().doPostProcessing = index == size - 1;
            ++index;
        }

        int index;
        const int size;
    };

    //! Method to set up a number of grid operators which are used
    //! in a joint assembling. It is assumed that all operators are
    //! specializations of the same template type
    template<typename GridOperatorTuple>
    static void setupGridOperators(GridOperatorTuple tuple)
    {
        Dune::ForEachValue<GridOperatorTuple> forEach(tuple);
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        forEach.apply(setup_visitor);
    }

    //! Interpolate the constrained dofs from given function
    template<typename F, typename X>
    void interpolate(const X& xold, F& f, X& x) const
    {
        // Interpolate f into grid function space and set corresponding coefficients
        Dune::PDELab::interpolate(f,
                                  global_assembler.trialGridFunctionSpace(), x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(
            local_assembler.trialConstraints(), xold, x);
    }

    //! Fill pattern of jacobian matrix
    void fill_pattern(Pattern & p) const
    {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        PatternEngine & pattern_engine =
            local_assembler.localPatternAssemblerEngine(p);
        global_assembler.assemble(pattern_engine);
    }

    //! Assemble residual
    void residual(const Domain & x, const DomainSl & xSl,
                  Range & r) const
    {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine =
            local_assembler.localResidualAssemblerEngineCoupling(r,
                                                                 x, xSl);
        global_assembler.assemble(residual_engine);
    }

    //! Assembler jacobian
    void jacobian(const Domain & x, Jacobian & a) const
    {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine =
            local_assembler.localJacobianAssemblerEngine(a, x);
        global_assembler.assemble(jacobian_engine);
    }

    //! Apply jacobian matrix without explicitly assembling it
    void jacobian_apply(const Domain & z, Range & r) const
    {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        JacobianApplyEngine & jacobian_apply_engine =
            local_assembler.localJacobianApplyAssemblerEngine(r, z);
        global_assembler.assemble(jacobian_apply_engine);

    }

    //! Apply jacobian matrix without explicitly assembling it
    void nonlinear_jacobian_apply(const Domain & x, const Domain & z,
                                  Range & r) const
    {
        global_assembler.assemble_coupling(
            local_assembler.localNonlinearJacobianApplyAssemblerEngine(
                r, x, z));
    }

    void make_consistent(Jacobian& a) const
    {
        dof_exchanger->accumulateBorderEntries(*this, a);
        dof_exchangerSl->accumulateBorderEntries(*this, a);
    }

    void update()
    {
        // the DOF exchanger has matrix information, so we need to update it
        dof_exchanger->update(*this);
        dof_exchangerSl->update(*this);
    }

    //! Get the matrix backend for this grid operator.
    const typename Traits::MatrixBackend& matrixBackend() const
    {
        return backend;
    }

private:
    Assembler global_assembler;
    shared_ptr<BorderDOFExchanger> dof_exchanger;
    shared_ptr<BorderDOFExchangerSl> dof_exchangerSl;

    mutable LocalAssembler local_assembler;

    MB backend;

};

}
}
#endif
