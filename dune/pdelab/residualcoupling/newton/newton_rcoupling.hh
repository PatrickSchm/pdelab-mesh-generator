// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_NEWTON_NEWTON_RCOUPLING_HH
#define DUNE_PDELAB_NEWTON_NEWTON_RCOUPLING_HH

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>

#include <math.h>

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/solver.hh>

namespace Dune
{
namespace PDELab
{





// Exception classes used in NewtonSolver
class NewtonError : public Exception {};
class NewtonDefectError : public NewtonError {};
class NewtonLinearSolverError : public NewtonError {};
class NewtonLineSearchError : public NewtonError {};
class NewtonNotConverged : public NewtonError {};



template<typename GFS, typename GFSSL, typename PNew, typename POld, typename UOld, typename FP>
void updatePNewton(const GFS & gfs, GFSSL & gfsSl, FP & fracturePairs, UOld & uold, UOld & unew, POld & pold, PNew & pnew) {

  const typename GFS::Ordering& ordering = gfs.ordering();
  const typename GFSSL::Ordering& orderingP = gfsSl.ordering();

  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFSU;
  LFSU lfsu(gfs);
  LFSU lfsuPair(gfs);

  typedef Dune::PDELab::LocalFunctionSpace<GFSSL> LFSUSL;
  LFSUSL lfsuSl(gfsSl);
  // using LFSUSL = Dune::TypeTree::Child<LFSUSL, Dune::TypeTree::Indices::_0>;
  typedef typename LFSUSL::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType URangeC;

  using EntitySet = typename GFS::Traits::EntitySet;
  using Element = typename EntitySet::Element;
  using Intersection = typename EntitySet::Intersection;
  using EntitySetC = typename GFSSL::Traits::EntitySet;
  using ElementC = typename EntitySetC::Element;

  const unsigned int usize = lfsu.size();

  typedef typename GFSSL::Traits::SizeType SizeType;
  typedef typename GFS::Traits::GridView GV;
  auto entitySet = gfs.entitySet();
  auto entitySetSl = gfsSl.entitySet();
  auto& indexSet = entitySet.indexSet();
  auto& indexSetSl = entitySetSl.indexSet();


  for (const auto & element : elements(entitySet)) {
    auto ids = indexSet.uniqueIndex(element);
    auto idsindex = indexSet.index(element);
    // std::cout << "------------------------- ELEMENT: " << idsindex << std::endl;

    auto fractureElements = fracturePairs.bulkFractureCoupledElements[idsindex];
    auto pElement = fracturePairs.bulkCoupledElements[idsindex];

    lfsu.bind(element);

    if (fractureElements.size() > 0 && pElement.size() > 0) {

      lfsuPair.bind(pElement[0]);

      Dune::PDELab::ElementGeometry<Element> egPair(pElement[0]);

      auto cBegin = fracturePairs.bulkFractureCoupledElements[idsindex].begin();
      auto cEnd = fracturePairs.bulkFractureCoupledElements[idsindex].end();
      for (auto & cElement = cBegin; cElement != cEnd; cElement++) {
        Dune::PDELab::ElementGeometry<ElementC> egC(*cElement);

        lfsuSl.bind(*cElement);

        unsigned int intersection_index = 0;
        for (const auto& intersection : intersections(entitySet, element))
        {

          if (intersection.boundary())
          {
            // if (idsindex == 5652)
            Dune::PDELab::ElementGeometry<Element> eg(element);
            Dune::PDELab::IntersectionGeometry<Intersection> ig(intersection, intersection_index);

            auto geo = eg.geometry();
            auto geoPair = egPair.geometry();
            auto geoC = egC.geometry();
            for (int j = 0; j < geoC.corners(); j++) {
              const typename GFSSL::Ordering::Traits::DOFIndex& diP = lfsuSl.dofIndex(j);
              typename GFSSL::Ordering::Traits::ContainerIndex ciP;

              orderingP.mapIndex(diP.view(), ciP);
              auto localC = geo.local(geoC.corner(j));
              auto igLocalC = ig.geometry().local(geoC.corner(j));
              auto localCPair = geoPair.local(geoC.corner(j));
              std::vector < URangeC > phi(usize);
              std::vector < URangeC > phiPair(usize);
              lfsu.child(0).finiteElement().localBasis().evaluateFunction(localC,
                                                                          phi);
              lfsu.child(0).finiteElement().localBasis().evaluateFunction(localCPair,
                                                                          phiPair);
              double ux_old = 0.0;
              double ux_new = 0.0;
              double ux_pair_old = 0.0;
              double ux_pair_new = 0.0;

              double uy_old = 0.0;
              double uy_new = 0.0;
              double uy_pair_old = 0.0;
              double uy_pair_new = 0.0;

              auto normal = ig.unitOuterNormal(igLocalC);

              for (unsigned i = 0; i < lfsu.child(0).size(); ++i)
              {
                const typename GFS::Ordering::Traits::DOFIndex& di0x = lfsu.child(0).dofIndex(i);
                typename GFS::Ordering::Traits::ContainerIndex cix;
                ordering.mapIndex(di0x.view(), cix);

                const typename GFS::Ordering::Traits::DOFIndex& di0Pairx = lfsuPair.child(0).dofIndex(i);
                typename GFS::Ordering::Traits::ContainerIndex ciPairx;
                ordering.mapIndex(di0Pairx.view(), ciPairx);


                ux_old += unew[cix] * phi[i];
                ux_new += uold[cix] * phi[i];

                ux_pair_old += unew[ciPairx] * phiPair[i];
                ux_pair_new += uold[ciPairx] * phiPair[i];
                const typename GFS::Ordering::Traits::DOFIndex& di0y = lfsu.child(1).dofIndex(i);
                typename GFS::Ordering::Traits::ContainerIndex ciy;
                ordering.mapIndex(di0y.view(), ciy);

                const typename GFS::Ordering::Traits::DOFIndex& di0Pairy = lfsuPair.child(1).dofIndex(i);
                typename GFS::Ordering::Traits::ContainerIndex ciPairy;
                ordering.mapIndex(di0Pairy.view(), ciPairy);

                uy_old += unew[ciy] * phi[i];
                uy_new += uold[ciy] * phi[i];

                uy_pair_old += unew[ciPairy] * phiPair[i];
                uy_pair_new += uold[ciPairy] * phiPair[i];
              }

              // double unewAvg = fabs(uy_new) + fabs(uy_pair_new);
              // double unewAvg = (ux_new - ux_pair_new) * normal[0] - (uy_new - uy_pair_new) * normal[1] ;

              // std::cout << "NORMAL X: " << normal[0]  << std::endl;
              // std::cout << "NORMAL Y: " << normal[1]  << std::endl;
              // std::cout << "DFX: " << dfx << std::endl;
              // std::cout << "DFY: " << dfy << std::endl;

              // std::cout << "OPEN: " << (uy_new - uy_old) - (uy_pair_new - uy_pair_old) << std::endl;
              // std::cout << "NEW: " << unewAvg + 0.005 << std::endl;
              if (fabs(normal[0]) < 1.e-4)
                normal[0] = 0;
              if (fabs(normal[1]) < 1.e-4)
                normal[1] = 0;
              double dfx =  ((ux_new - ux_old) * normal[0] - (ux_pair_new - ux_pair_old) * normal[0]);
              double dfy =  ((uy_new - uy_old) * normal[1] - (uy_pair_new - uy_pair_old) * normal[1]);
// std::cout << "Coordinate: " << localC << std::endl;
              auto global = geoC.corner(j);
              pold[ciP] = pnew[ciP] + (dfx + dfy) / (1.e-4) * 2.2e9 / 2.0;

              if (global[1] <= 5.0 + 1.e-6 && global [0] <= 2.275 + 1.e-6 && global [0] >= 2.275 - 1.e-6) {
                std::cout << "Coordinate: " << geoC.corner(j) << std::endl;
                std::cout << "U NEW: " << uy_new  << std::endl;
                std::cout << "U OLD: " << uy_old  << std::endl;
                std::cout << "UPNEW: " << uy_pair_new  << std::endl;
                std::cout << "UPOLD: " << uy_pair_old << std::endl;
                std::cout << "UPDATE ADD X: " << ((ux_new - ux_pair_new) * normal[0] - (ux_old - ux_pair_old) * normal[0]) / (0.005) * 2.2e7 / 2 << std::endl;
                std::cout << "UPDATE ADD Y: " << ((uy_new - uy_pair_new) * normal[1] - (uy_old - uy_pair_old) * normal[1]) / (0.005) * 2.2e7 / 2 << std::endl;
                std::cout << "P: " << pold[ciP]  << std::endl;

              }
// std::cout << "UPDATE ADD y: " <<-((uy_new - uy_pair_new) * normal[1] - (uy_old - uy_pair_old) * normal[1])/ (0.005) * 2.2e7 << std::endl;


            }
          }
          ++intersection_index;
        }
      }
    }
  }
  pnew = pold;
  // unew = uold;

}



// Status information of Newton's method
template<class RFType>
struct NewtonResult : LinearSolverResult<RFType>
{
  RFType first_defect;       // the first defect
  RFType defect;             // the final defect
  double assembler_time;     // Cumulative time for matrix assembly
  double linear_solver_time; // Cumulative time for linear solver
  int linear_solver_iterations; // Total number of linear iterations

  NewtonResult() :
    first_defect(0.0), defect(0.0), assembler_time(0.0), linear_solver_time(0.0),
    linear_solver_iterations(0) {}
};

template<class GOS, class GOSSL, class FP, class TrlV, class TrlVSl, class TstV>
class NewtonBase
{
  typedef GOS GridOperator;
  typedef GOSSL GridOperatorSl;
  typedef TrlV TrialVector;
  typedef TrlVSl TrialVectorSl;
  typedef TstV TestVector;

  typedef typename TestVector::ElementType RFType;
  typedef typename GOS::Traits::Jacobian Matrix;


public:
  // export result type
  typedef NewtonResult<RFType> Result;

  void setVerbosityLevel(unsigned int verbosity_level)
  {
    if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank() > 0)
      verbosity_level_ = 0;
    else
      verbosity_level_ = verbosity_level;
  }

  //! Set whether the jacobian matrix should be kept across calls to apply().
  void setKeepMatrix(bool b)
  {
    keep_matrix_ = b;
  }

  //! Return whether the jacobian matrix is kept across calls to apply().
  bool keepMatrix() const
  {
    return keep_matrix_;
  }

  //! Discard the stored Jacobian matrix.
  void discardMatrix()
  {
    if (A_)
      A_.reset();
  }

protected:
  const GridOperator& gridoperator_;
  const GridOperatorSl &gridoperatorSl_;
  const FP *fracturePairs_;
  TrialVector *u_;
  TrialVector *uOld_;
  TrialVectorSl *uSl_;

  std::shared_ptr<TrialVector> z_;
  std::shared_ptr<TestVector> r_;
  std::shared_ptr<Matrix> A_;
  Result res_;
  unsigned int verbosity_level_;
  RFType prev_defect_;
  RFType linear_reduction_;
  bool reassembled_;
  RFType reduction_;
  RFType abs_limit_;
  bool keep_matrix_;

  NewtonBase(const GridOperator& go, const GridOperatorSl& goSl, const FP& fracturePairs, TrialVector& u)
    : gridoperator_(go)
    , gridoperatorSl_(goSl)
    , fracturePairs_(&fracturePairs)
    , u_(&u)
    , verbosity_level_(1)
    , keep_matrix_(true)
  {
    if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank() > 0)
      verbosity_level_ = 0;
  }

  NewtonBase(const GridOperator& go)
    : gridoperator_(go)
    , u_(0)
    , verbosity_level_(1)
    , keep_matrix_(true)
  {
    if (gridoperator_.trialGridFunctionSpace().gridView().comm().rank() > 0)
      verbosity_level_ = 0;
  }

  virtual ~NewtonBase() { }

  virtual bool terminate() = 0;
  virtual void prepare_step(Matrix& A, TestVector& r) = 0;
  virtual void line_search(TrialVector& z, TestVector& r) = 0;
  virtual void defect(TestVector& r) = 0;
}; // end class NewtonBase

template<class GOS, class GOSSL, class FP, class S, class TrlV, class TrlVSl, class TstV>
class NewtonSolver : public virtual NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
{
  typedef S Solver;
  typedef GOS GridOperator;
  typedef TrlV TrialVector;
  typedef TrlVSl TrialVectorSl;
  typedef TstV TestVector;
  const FP *fracturePairs_;

  typedef typename TestVector::ElementType RFType;
  typedef typename GOS::Traits::Jacobian Matrix;

public:
  typedef NewtonResult<RFType> Result;

  NewtonSolver(const GridOperator& go, const GOSSL& goSl, const FP& fracturePairs, TrialVector& u_, TrialVector& uOld_, TrialVectorSl& uSl, Solver& solver)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs, u_)
    , fracturePairs_(&fracturePairs)
    , solver_(solver)
    , result_valid_(false)
  {}

  NewtonSolver(const GridOperator& go, Solver& solver)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , solver_(solver)
    , result_valid_(false)
  {}

  void apply();

  void apply(TrialVector& u_, GOSSL &goSl, FP& fracturePairs_, TrialVector& uOld_, TrialVectorSl& uSl_);

  const Result& result() const
  {
    if (!result_valid_)
      DUNE_THROW(NewtonError,
                 "NewtonSolver::result() called before NewtonSolver::solve()");
    return this->res_;
  }

protected:
  virtual void defect(TestVector& r)
  {
    r = 0.0;                                        // TODO: vector interface

    this->gridoperator_.residual(*this->u_, *this->uOld_, *this->uSl_, r);
    this->res_.defect = this->solver_.norm(r);                    // TODO: solver interface
    if (!std::isfinite(this->res_.defect))
      DUNE_THROW(NewtonDefectError,
                 "NewtonSolver::defect(): Non-linear defect is NaN or Inf");
  }


private:
  void linearSolve(Matrix& A, TrialVector& z, TestVector& r) const
  {
    if (this->verbosity_level_ >= 4)
      std::cout << "      Solving linear system..." << std::endl;
    z = 0.0;                                        // TODO: vector interface
    this->solver_.apply(A, z, r, this->linear_reduction_);        // TODO: solver interface

    ios_base_all_saver restorer(std::cout); // store old ios flags

    if (!this->solver_.result().converged)                 // TODO: solver interface
      DUNE_THROW(NewtonLinearSolverError,
                 "NewtonSolver::linearSolve(): Linear solver did not converge "
                 "in " << this->solver_.result().iterations << " iterations");
    if (this->verbosity_level_ >= 4)
      std::cout << "          linear solver iterations:     "
                << std::setw(12) << solver_.result().iterations << std::endl
                << "          linear defect reduction:      "
                << std::setw(12) << std::setprecision(4) << std::scientific
                << solver_.result().reduction << std::endl;
  }

  Solver& solver_;
  bool result_valid_;
}; // end class NewtonSolver

template<class GOS, class GOSSL, class FP, class S, class TrlV, class TrlVSl, class TstV>
void NewtonSolver<GOS, GOSSL, FP, S, TrlV, TrlVSl, TstV>::apply(TrialVector& u, GOSSL& goSl, FP& fracturePairs, TrialVector& uOld, TrialVectorSl& uSl)
{
  this->u_ = &u;
  this->uOld_ = &uOld;
  this->uSl_ = &uSl;
  this->fracturePairs_ = &fracturePairs;
  apply();
}

template<class GOS, class GOSSL, class FP, class S, class TrlV, class TrlVSl, class TstV>
void NewtonSolver<GOS, GOSSL, FP, S, TrlV, TrlVSl, TstV>::apply()
{
  this->res_.iterations = 0;
  this->res_.converged = false;
  this->res_.reduction = 1.0;
  this->res_.conv_rate = 1.0;
  this->res_.elapsed = 0.0;
  this->res_.assembler_time = 0.0;
  this->res_.linear_solver_time = 0.0;
  this->res_.linear_solver_iterations = 0;
  result_valid_ = true;
  Timer timer;

  try
  {
    if (!this->r_) {
      // std::cout << "=== Setting up residual vector ..." << std::endl;
      this->r_ = std::make_shared<TestVector>(this->gridoperator_.testGridFunctionSpace());
    }
    // residual calculation in member function "defect":
    //--------------------------------------------------
    // - set residual vector to zero
    // - calculate new residual
    // - store norm of residual in "this->res_.defect"
    this->defect(*this->r_);
    this->res_.first_defect = this->res_.defect;
    this->prev_defect_ = this->res_.defect;

    if (this->verbosity_level_ >= 2)
    {
      // store old ios flags
      ios_base_all_saver restorer(std::cout);
      std::cout << "  Initial defect: "
                << std::setw(12) << std::setprecision(4) << std::scientific
                << this->res_.defect << std::endl;
    }

    if (!this->A_) {
      // std::cout << "==== Setting up jacobian matrix ... " << std::endl;
      this->A_ = std::make_shared<Matrix>(this->gridoperator_);
    }
    if (!this->z_) {
      // std::cout << "==== Setting up correction vector ... " << std::endl;
      this->z_ = std::make_shared<TrialVector>(this->gridoperator_.trialGridFunctionSpace());
    }

    while (!this->terminate())
    {
      if (this->verbosity_level_ >= 3)
        std::cout << "  Newton iteration " << this->res_.iterations
                  << " --------------------------------" << std::endl;

      Timer assembler_timer;
      try
      {
        // jacobian calculation in member function "prepare_step"
        //-------------------------------------------------------
        // - if above reassemble threshold
        //   - set jacobian to zero
        //   - calculate new jacobian
        // - set linear reduction
        this->prepare_step(*this->A_, *this->r_);
      }
      catch (...)
      {
        this->res_.assembler_time += assembler_timer.elapsed();
        throw;
      }
      double assembler_time = assembler_timer.elapsed();
      this->res_.assembler_time += assembler_time;
      if (this->verbosity_level_ >= 3)
        std::cout << "      matrix assembly time:             "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << assembler_time << std::endl;

      Timer linear_solver_timer;
      try
      {
        // solution of linear system in member function "linearSolve"
        //-----------------------------------------------------------
        // - set initial guess for correction z to zero
        // - call linear solver
        this->linearSolve(*this->A_, *this->z_, *this->r_);
      }
      catch (...)
      {
        this->res_.linear_solver_time += linear_solver_timer.elapsed();
        this->res_.linear_solver_iterations += this->solver_.result().iterations;
        throw;
      }
      double linear_solver_time = linear_solver_timer.elapsed();
      this->res_.linear_solver_time += linear_solver_time;
      this->res_.linear_solver_iterations += this->solver_.result().iterations;

      try
      {
        // line search with correction z
        // the undamped version is also integrated in here
        // this->uOld_ = &prev_u;

        this->line_search(*this->z_, *this->r_);
        auto gfs = this->gridoperator_.testGridFunctionSpace();
        auto gfsSl = this->gridoperatorSl_.testGridFunctionSpace();
        auto pold = *this->uSl_;
        updatePNewton(gfs, gfsSl, *this->fracturePairs_, *this->u_, *this->uOld_, pold, *this->uSl_);
        this->uOld_ = this->u_;
      }
      catch (NewtonLineSearchError)
      {
        if (this->reassembled_)
          throw;
        if (this->verbosity_level_ >= 3)
          std::cout << "      line search failed - trying again with reassembled matrix" << std::endl;
        continue;
      }

      this->res_.reduction = this->res_.defect / this->res_.first_defect;
      this->res_.iterations++;
      this->res_.conv_rate = std::pow(this->res_.reduction, 1.0 / this->res_.iterations);

      // store old ios flags
      ios_base_all_saver restorer(std::cout);

      if (this->verbosity_level_ >= 3)
        std::cout << "      linear solver time:               "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << linear_solver_time << std::endl
                  << "      defect reduction (this iteration):"
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.defect / this->prev_defect_ << std::endl
                  << "      defect reduction (total):         "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.reduction << std::endl
                  << "      new defect:                       "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.defect << std::endl;
      if (this->verbosity_level_ == 2)
        std::cout << "  Newton iteration " << std::setw(2) << this->res_.iterations
                  << ".  New defect: "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.defect
                  << ".  Reduction (this): "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.defect / this->prev_defect_
                  << ".  Reduction (total): "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << this->res_.reduction << std::endl;
    } // end while
  } // end try
  catch (...)
  {
    this->res_.elapsed = timer.elapsed();
    throw;
  }
  this->res_.elapsed = timer.elapsed();

  ios_base_all_saver restorer(std::cout); // store old ios flags

  if (this->verbosity_level_ == 1)
    std::cout << "  Newton converged after " << std::setw(2) << this->res_.iterations
              << " iterations.  Reduction: "
              << std::setw(12) << std::setprecision(4) << std::scientific
              << this->res_.reduction
              << "   (" << std::setprecision(4) << this->res_.elapsed << "s)"
              << std::endl;

  if (!this->keep_matrix_)
    this->A_.reset();
} // end apply

template<class GOS, class GOSSL, class FP, class TrlV, class TrlVSl, class TstV>
class NewtonTerminate : public virtual NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
{
  typedef GOS GridOperator;
  typedef TrlV TrialVector;

  typedef typename TstV::ElementType RFType;

public:
  NewtonTerminate(const GridOperator& go, const GOSSL goSl, const FP fracturePairs_, TrialVector& u_)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , maxit_(40)
    , force_iteration_(false)
  {
    this->reduction_ = 1e-8;
    this->abs_limit_ = 1e-12;
  }

  NewtonTerminate(const GridOperator& go)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , maxit_(40)
    , force_iteration_(false)
  {
    this->reduction_ = 1e-8;
    this->abs_limit_ = 1e-12;
  }

  void setReduction(RFType reduction)
  {
    this->reduction_ = reduction;
  }

  void setMaxIterations(unsigned int maxit)
  {
    maxit_ = maxit;
  }

  void setForceIteration(bool force_iteration)
  {
    force_iteration_ = force_iteration;
  }

  void setAbsoluteLimit(RFType abs_limit_)
  {
    this->abs_limit_ = abs_limit_;
  }

  virtual bool terminate()
  {
    if (force_iteration_ && this->res_.iterations == 0)
      return false;
    this->res_.converged = this->res_.defect < this->abs_limit_
                           || this->res_.defect < this->res_.first_defect * this->reduction_;

    if (this->res_.iterations >= maxit_ && !this->res_.converged)
      DUNE_THROW(NewtonNotConverged,
                 "NewtonTerminate::terminate(): Maximum iteration count reached");
    return this->res_.converged;
  }

private:
  unsigned int maxit_;
  bool force_iteration_;
}; // end class NewtonTerminate

template<class GOS, class GOSSL, class FP, class TrlV, class TrlVSl, class TstV>
class NewtonPrepareStep : public virtual NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
{
  typedef GOS GridOperator;
  typedef TrlV TrialVector;

  typedef typename TstV::ElementType RFType;
  typedef typename GOS::Traits::Jacobian Matrix;

public:
  NewtonPrepareStep(const GridOperator& go, const GOSSL& goSl, const FP fracturePairs_, TrialVector& u_)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , min_linear_reduction_(1e-3)
    , fixed_linear_reduction_(0.0)
    , reassemble_threshold_(0.0)
  {}

  NewtonPrepareStep(const GridOperator& go)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , min_linear_reduction_(1e-3)
    , fixed_linear_reduction_(0.0)
    , reassemble_threshold_(0.0)
  {}

  /**\brief set the minimal reduction in the linear solver

     \note with min_linear_reduction > 0, the linear reduction will be
     determined as mininum of the min_linear_reduction and the
     linear_reduction needed to achieve second order
     Newton convergence. */
  void setMinLinearReduction(RFType min_linear_reduction)
  {
    min_linear_reduction_ = min_linear_reduction;
  }

  /**\brief set a fixed reduction in the linear solver (overwrites setMinLinearReduction)

     \note with fixed_linear_reduction > 0, the linear reduction
     rate will always be fixed to min_linear_reduction. */
  void setFixedLinearReduction(bool fixed_linear_reduction)
  {
    fixed_linear_reduction_ = fixed_linear_reduction;
  }

  /**\brief set a threshold, when the linear operator is reassembled

     We allow to keep the linear operator over several newton
     iterations. If the reduction in the newton drops below a
     given threshold the linear operator is reassembled to ensure
     convergence.
   */
  void setReassembleThreshold(RFType reassemble_threshold)
  {
    reassemble_threshold_ = reassemble_threshold;
  }

  virtual void prepare_step(Matrix& A, TstV& )
  {
    this->reassembled_ = false;
    if (this->res_.defect / this->prev_defect_ > reassemble_threshold_)
    {
      if (this->verbosity_level_ >= 3)
        std::cout << "      Reassembling matrix..." << std::endl;
      A = 0.0;                                    // TODO: Matrix interface
      this->gridoperator_.jacobian(*this->u_, *this->uOld_, *this->uSl_, A);
      this->reassembled_ = true;
    }

    if (fixed_linear_reduction_ == true)
      this->linear_reduction_ = min_linear_reduction_;
    else {
      // determine maximum defect, where Newton is converged.
      RFType stop_defect =
        std::max(this->res_.first_defect * this->reduction_,
                 this->abs_limit_);

      /*
        To achieve second order convergence of newton
        we need a linear reduction of at least
        current_defect^2/prev_defect^2.
        For the last newton step a linear reduction of
        1/10*end_defect/current_defect
        is sufficient for convergence.
      */
      if ( stop_defect / (10 * this->res_.defect) >
           this->res_.defect * this->res_.defect / (this->prev_defect_ * this->prev_defect_) )
        this->linear_reduction_ =
          stop_defect / (10 * this->res_.defect);
      else
        this->linear_reduction_ =
          std::min(min_linear_reduction_, this->res_.defect * this->res_.defect / (this->prev_defect_ * this->prev_defect_));
    }

    this->prev_defect_ = this->res_.defect;

    ios_base_all_saver restorer(std::cout); // store old ios flags

    if (this->verbosity_level_ >= 3)
      std::cout << "      requested linear reduction:       "
                << std::setw(12) << std::setprecision(4) << std::scientific
                << this->linear_reduction_ << std::endl;
  }

private:
  RFType min_linear_reduction_;
  bool fixed_linear_reduction_;
  RFType reassemble_threshold_;
}; // end class NewtonPrepareStep

template<class GOS, class GOSSL, class FP, class TrlV, class TrlVSl, class TstV>
class NewtonLineSearch : public virtual NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
{
  typedef GOS GridOperator;
  typedef TrlV TrialVector;
  typedef TstV TestVector;

  typedef typename TestVector::ElementType RFType;

public:
  enum Strategy {
    /** \brief don't do any line search or damping */
    noLineSearch,
    /** \brief perform a linear search for the optimal damping parameter with multiples of damping

     the strategy was described in <a href="http://dx.doi.org/10.1007/BF01406516">[Hackbusch and Reusken, 1989]</a> */
    hackbuschReusken,
    /** \brief same as hackbuschReusken, but doesn't fail if the best update is still not good enough */
    hackbuschReuskenAcceptBest
  };

  NewtonLineSearch(const GridOperator& go, const GOSSL& goSl, const FP & fracturePairs_, TrialVector& u_)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , strategy_(hackbuschReusken)
    , maxit_(10)
    , damping_factor_(0.5)
  {}

  NewtonLineSearch(const GridOperator& go)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , strategy_(hackbuschReusken)
    , maxit_(10)
    , damping_factor_(0.5)
  {}

  void setLineSearchStrategy(Strategy strategy)
  {
    strategy_ = strategy;
  }

  void setLineSearchStrategy(std::string strategy)
  {
    strategy_ = strategyFromName(strategy);
  }

  void setLineSearchMaxIterations(unsigned int maxit)
  {
    maxit_ = maxit;
  }

  void setLineSearchDampingFactor(RFType damping_factor)
  {
    damping_factor_ = damping_factor;
  }

  virtual void line_search(TrialVector& z, TestVector& r)
  {
    if (strategy_ == noLineSearch)
    {
      this->u_->axpy(-1.0, z);                     // TODO: vector interface
      this->defect(r);
      return;
    }

    if (this->verbosity_level_ >= 4)
      std::cout << "      Performing line search..." << std::endl;
    RFType lambda = 1.0;
    RFType best_lambda = 0.0;
    RFType best_defect = this->res_.defect;
    TrialVector prev_u(*this->u_);  // TODO: vector interface
    unsigned int i = 0;
    ios_base_all_saver restorer(std::cout); // store old ios flags

    while (1)
    {
      if (this->verbosity_level_ >= 4)
        std::cout << "          trying line search damping factor:   "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << lambda
                  << std::endl;


      this->u_->axpy(-lambda, z);                  // TODO: vector interface

      // std::cout << "NORM: " << r.two_norm() << std::endl;
      try {
        this->defect(r);
      }
      catch (NewtonDefectError)
      {
        if (this->verbosity_level_ >= 4)
          std::cout << "          NaNs detected" << std::endl;
      }       // ignore NaNs and try again with lower lambda

      if (this->res_.defect <= (1.0 - lambda / 4) * this->prev_defect_)
      {
        if (this->verbosity_level_ >= 4)
          std::cout << "          line search converged" << std::endl;
        break;
      }

      if (this->res_.defect < best_defect)
      {
        best_defect = this->res_.defect;
        best_lambda = lambda;
      }
std::cout << "DEFECT: " << this->res_.defect <<std::endl;
      std::cout << "LAMBDA: " << lambda << std::endl;
      std::cout << "BEST LAMBDA: " << best_lambda << std::endl;
      if (++i >= maxit_)
      {
        if (this->verbosity_level_ >= 4)
          std::cout << "          max line search iterations exceeded" << std::endl;
        switch (strategy_)
        {
        case hackbuschReusken:
          *this->u_ = prev_u;
          this->defect(r);
          DUNE_THROW(NewtonLineSearchError,
                     "NewtonLineSearch::line_search(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve enough");
        case hackbuschReuskenAcceptBest:
          if (best_lambda == 0.0)
          {
            *this->u_ = prev_u;
            this->defect(r);
            // DUNE_THROW(NewtonLineSearchError,
            //            "NewtonLineSearch::line_search(): line search failed, "
            //            "max iteration count reached, "
            //            "defect did not improve in any of the iterations");
          }
          if (best_lambda != lambda)
          {
            *this->u_ = prev_u;
            this->u_->axpy(-best_lambda, z);
            this->defect(r);
          }
          break;
        case noLineSearch:
          break;
        }
        break;
      }

      lambda *= damping_factor_;
      *this->u_ = prev_u;
    }
    if (this->verbosity_level_ >= 4)
      std::cout << "          line search damping factor:   "
                << std::setw(12) << std::setprecision(4) << std::scientific
                << lambda << std::endl;
  } // end line_search

protected:
  /** helper function to get the different strategies from their name */
  Strategy strategyFromName(const std::string & s) {
    if (s == "noLineSearch")
      return noLineSearch;
    if (s == "hackbuschReusken")
      return hackbuschReusken;
    if (s == "hackbuschReuskenAcceptBest")
      return hackbuschReuskenAcceptBest;
    DUNE_THROW(Exception, "unknown line search strategy" << s);
  }

private:
  Strategy strategy_;
  unsigned int maxit_;
  RFType damping_factor_;
}; // end class NewtonLineSearch

template<class GOS, class GOSSL, class FP, class S, class TrlV, class TrlVSl, class TstV = TrlV>
class Newton : public NewtonSolver<GOS, GOSSL, FP, S, TrlV, TrlVSl, TstV>
  , public NewtonTerminate<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
  , public NewtonLineSearch<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
  , public NewtonPrepareStep<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>
{
  typedef GOS GridOperator;
  typedef GOSSL GridOperatorSl;
  typedef S Solver;
  typedef TrlV TrialVector;
  typedef TrlVSl TrialVectorSl;

public:
  Newton(const GridOperator& go, const GridOperatorSl& goSl, const FP & fracturePairs_, TrialVector& u_, TrialVector& uOld_, TrialVectorSl& uSl_, Solver& solver_)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , NewtonSolver<GOS, GOSSL, FP, S, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_,  u_, uOld_, uSl_, solver_)
    , NewtonTerminate<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , NewtonLineSearch<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
    , NewtonPrepareStep<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go, goSl, fracturePairs_, u_)
  {}
  Newton(const GridOperator& go, Solver& solver_)
    : NewtonBase<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , NewtonSolver<GOS, GOSSL, FP, S, TrlV, TrlVSl, TstV>(go, solver_)
    , NewtonTerminate<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , NewtonLineSearch<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
    , NewtonPrepareStep<GOS, GOSSL, FP, TrlV, TrlVSl, TstV>(go)
  {}
  //! interpret a parameter tree as a set of options for the newton solver
  /**

     example configuration:

     \code
     [NewtonParameters]

     ReassembleThreshold = 0.1
     LineSearchMaxIterations = 10
     MaxIterations = 7
     AbsoluteLimit = 1e-6
     Reduction = 1e-4
     LinearReduction = 1e-3
     LineSearchDamping  = 0.9
     \endcode

     and invocation in the code:
     \code
     newton.setParameters(param.sub("NewtonParameters"));
     \endcode
  */
  void setParameters(Dune::ParameterTree & param)
  {
    typedef typename TstV::ElementType RFType;
    if (param.hasKey("VerbosityLevel"))
      this->setVerbosityLevel(
        param.get<unsigned int>("VerbosityLevel"));
    if (param.hasKey("Reduction"))
      this->setReduction(
        param.get<RFType>("Reduction"));
    if (param.hasKey("MaxIterations"))
      this->setMaxIterations(
        param.get<unsigned int>("MaxIterations"));
    if (param.hasKey("ForceIteration"))
      this->setForceIteration(
        param.get<bool>("ForceIteration"));
    if (param.hasKey("AbsoluteLimit"))
      this->setAbsoluteLimit(
        param.get<RFType>("AbsoluteLimit"));
    if (param.hasKey("MinLinearReduction"))
      this->setMinLinearReduction(
        param.get<RFType>("MinLinearReduction"));
    if (param.hasKey("FixedLinearReduction"))
      this->setFixedLinearReduction(
        param.get<bool>("FixedLinearReduction"));
    if (param.hasKey("ReassembleThreshold"))
      this->setReassembleThreshold(
        param.get<RFType>("ReassembleThreshold"));
    if (param.hasKey("LineSearchStrategy"))
      this->setLineSearchStrategy(
        param.get<std::string>("LineSearchStrategy"));
    if (param.hasKey("LineSearchMaxIterations"))
      this->setLineSearchMaxIterations(
        param.get<unsigned int>("LineSearchMaxIterations"));
    if (param.hasKey("LineSearchDampingFactor"))
      this->setLineSearchDampingFactor(
        param.get<RFType>("LineSearchDampingFactor"));
    if (param.hasKey("KeepMatrix"))
      this->setKeepMatrix(
        param.get<bool>("KeepMatrix"));
  }
}; // end class Newton
} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_NEWTON_NEWTON_RCOUPLING_HH
