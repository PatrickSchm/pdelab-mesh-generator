// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLERUTILITIES_RCOUPLING_HH
#define DUNE_PDELAB_ASSEMBLERUTILITIES_RCOUPLING_HH

#include <algorithm>
#include <tuple>

#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>

namespace Dune
{
    namespace PDELab
    {

        /** Traits of the local assembler

         \tparam GO The grid operator

         */
        template<typename GO>
        struct LocalAssemblerTraitsCoupling
        {
            //! The trial grid function space.
            typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
            typedef typename GO::Traits::TrialGridFunctionSpaceSl TrialGridFunctionSpaceSl;

            //! The test grid function space.
            typedef typename GO::Traits::TestGridFunctionSpace TestGridFunctionSpace;
            typedef typename GO::Traits::TestGridFunctionSpaceSl TestGridFunctionSpaceSl;

            //! The type of the trial grid function space constraints.
            typedef typename GO::Traits::TrialGridFunctionSpaceConstraints TrialGridFunctionSpaceConstraints;
            typedef typename GO::Traits::TrialGridFunctionSpaceConstraintsSl TrialGridFunctionSpaceConstraintsSl;

            //! The type of the test grid function space constraints.
            typedef typename GO::Traits::TestGridFunctionSpaceConstraints TestGridFunctionSpaceConstraints;
            typedef typename GO::Traits::TestGridFunctionSpaceConstraintsSl TestGridFunctionSpaceConstraintsSl;

            //! The matrix backend of the grid operator.
            typedef typename GO::Traits::MatrixBackend MatrixBackend;
            typedef typename GO::Traits::MatrixBackendSl MatrixBackendSl;

            //! The field type of the domain (solution).
            typedef typename GO::Traits::DomainField DomainField;
            typedef typename GO::Traits::DomainFieldSl DomainFieldSl;

            //! The type of the domain (solution).
            typedef typename GO::Traits::Domain Solution;
            typedef typename GO::Traits::DomainSl SolutionSl;

            //! The field type of the range (residual).
            typedef typename GO::Traits::RangeField RangeField;

            //! The type of the range (residual).
            typedef typename GO::Traits::Range Residual;
            typedef typename GO::Traits::RangeSl ResidualSl;

            //! The field type of the jacobian.
            typedef typename GO::Traits::JacobianField JacobianField;
            typedef typename GO::Traits::JacobianFieldSl JacobianFieldSl;

            //! The type of the jacobian.
            typedef typename GO::Traits::Jacobian Jacobian;
            typedef typename GO::Traits::JacobianSl JacobianSl;

            //! The matrix pattern
            typedef typename MatrixBackend::template Pattern<Jacobian,
                    TestGridFunctionSpace, TrialGridFunctionSpace> MatrixPattern;
            typedef typename MatrixBackendSl::template Pattern<JacobianSl,
                    TestGridFunctionSpaceSl, TrialGridFunctionSpaceSl> MatrixPatternSl;

            //! The helper class to exchange data on the processor boundary
            typedef typename GO::BorderDOFExchanger BorderDOFExchanger;

            //! The helper class to exchange data on the processor boundary
            typedef typename GO::BorderDOFExchanger BorderDOFExchangerSl;
        };

        // ********************************************************************************
        // default local pattern implementation
        // ********************************************************************************

        /**
         \brief Entry in sparsity pattern

         The sparsity pattern of a linear operator is described by by connecting
         degrees of freedom in one element with degrees of freedom in the
         same element (intra) or an intersecting element (inter).

         This numbering is with respect to the depth-first canonical order of the
         degrees of freedom of an entity.

         \nosubgrouping
         */
    }
} // namespace Dune

#endif //DUNE_PDELAB_ASSEMBLERUTILITIES_HH
