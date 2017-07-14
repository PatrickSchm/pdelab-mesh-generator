// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_RCOUPLING_HH
#define DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_RCOUPLING_HH

#include <dune/common/concept.hh>

namespace Dune
{
namespace PDELab
{

#ifndef DOXYGEN

namespace impl
{

// ********************************************************************************
// concept checks that test whether a local operator provides a given apply method
// these are used to emit better error messages for the two variants of
// apply methods
// ********************************************************************************


template<typename ... Args>
struct HasJacobianApplyCoupling
{
  template<typename LO>
  auto require(LO&& lo) -> decltype(
    lo.jacobian_apply_coupling(std::declval<Args>()...)
  );
};

} // namespace impl

#endif // DOXYGEN

// compile time switching of function call
template<typename LA, bool doIt>
struct LocalAssemblerCallSwitchRCoupling
{
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_volume(const LA& la, const LFSU& lfsu,
                             const LFSV& lfsv, LocalPattern& pattern)
  {
  }

  template<typename LFSU, typename LFSV, typename LFSUC,
           typename LFSVC, typename LocalPattern>
  static void pattern_coupling(const LA& la, const LFSU& lfsu,
                               const LFSV& lfsv, const LFSUC& lfsuC,
                               const LFSVC& lfsvC, LocalPattern& pattern_coupling)
  {
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_volume_post_skeleton(const LA& la,
                                           const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
  {
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_skeleton(const LA& la, const LFSU& lfsu_s,
                               const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                               LocalPattern& pattern_sn, LocalPattern& pattern_ns)
  {
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_boundary(const LA& la, const LFSU& lfsu_s,
                               const LFSV& lfsv_s, LocalPattern& pattern_ss)
  {
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_volume(const LA& la, const EG& eg,
                           const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r, int iMat)
  {
  }

  template<typename IG, typename EGP, typename EGC, typename LFSU, typename LFSV, typename X, typename R, typename LFSUC,
           typename XC, typename LFSVC>
  static void boundary_coupling(const LA& la, const IG& ig,
                                const EGC& egC, const EGP& egPair, const LFSU& lfsu, const LFSU& lfsu_pair, const X& x, const X& xPair, const X& xOld, const LFSV& lfsv, const LFSV& lfsv_pair, const LFSUC& lfsuC,
                                const XC& xC, const LFSVC& lfsvC,
                                R& r)
  {
  }

  template<typename EG, typename IGC, typename LFSU,
           typename X, typename LFSV, typename LFSUC, typename XC,
           typename LFSVC, typename FH, typename R>
  static void boundary_coupling_reversed(const LA& la, const EG& eg,
                                         const IGC& igC, const LFSU& lfsu,
                                         const X& x, const LFSV& lfsv, const LFSUC& lfsuC,
                                         const XC& xC, const LFSVC& lfsvC, std::vector<FH> fractureHeight,
                                         R& r)
  {
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_volume_post_skeleton(const LA& la, const EG& eg,
                                         const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_skeleton(const LA& la, const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                             R& r_s, R& r_n)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_boundary(const LA& la, const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             R& r_s)
  {
  }

  template<typename EG, typename LFSV, typename R>
  static void lambda_volume(const LA& la, const EG& eg,
                            const LFSV& lfsv, R& r)
  {
  }
  template<typename EG, typename LFSV, typename R>
  static void lambda_volume_post_skeleton(const LA& la, const EG& eg,
                                          const LFSV& lfsv, R& r)
  {
  }
  template<typename IG, typename LFSV, typename R>
  static void lambda_skeleton(const LA& la, const IG& ig,
                              const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n)
  {
  }
  template<typename IG, typename LFSV, typename R>
  static void lambda_boundary(const LA& la, const IG& ig,
                              const LFSV& lfsv, R& r)
  {
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_volume(const LA& la, const EG& eg,
                                    const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
  {
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_coupling(const LA& la, const EG& eg,
                                      const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
  {
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_volume_post_skeleton(const LA& la,
                                                  const EG& eg, const LFSU& lfsu, const X& x,
                                                  const LFSV& lfsv, Y& y)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_skeleton(const LA& la, const IG& ig,
                                      const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                      const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                      Y& y_s, Y& y_n)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_boundary(const LA& la, const IG& ig,
                                      const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                      Y& y_s)
  {
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_volume(const LA& la,
                                              const EG& eg, const LFSU& lfsu, const X& x, const X& z,
                                              const LFSV& lfsv, Y& y)
  {
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_volume_post_skeleton(
    const LA& la, const EG& eg, const LFSU& lfsu, const X& x,
    const X& z, const LFSV& lfsv, Y& y)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_skeleton(const LA& la,
                                                const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                                const X& z_s, const LFSV& lfsv_s, const LFSU& lfsu_n,
                                                const X& x_n, const X& z_n, const LFSV& lfsv_n, Y& y_s,
                                                Y& y_n)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_boundary(const LA& la,
                                                const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                                const X& z_s, const LFSV& lfsv_s, Y& y_s)
  {
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_volume(const LA& la, const EG& eg,
                              const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
  {
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_coupling(const LA& la, const EG& eg,
                                const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
  {
  }

  // volume integral depending on test and ansatz functions
  template< typename EG, typename IGC, typename LFSU,
            typename X, typename LFSV, typename LFSUC,
            typename XC, typename LFSVC, typename FH, typename M>
  static void jacobian_boundary_coupling_reversed(const LA& la, const EG& eg, const IGC& igC, const LFSU& lfsu,
                                                  const X& x, const X& xOld, const LFSV& lfsv, const LFSUC& lfsuC,
                                                  const XC& xC, const LFSVC& lfsvC, const std::vector<FH>& fractureHeight,  const M&mat)
  {
  }

  // volume integral depending on test and ansatz functions
  template< typename IG, typename EGC, typename LFSU, typename X, typename LFSV, typename LFSUC,
            typename XC, typename LFSVC, typename M>
  static void jacobian_boundary_coupling(const LA& la, const IG& ig,
                                         const EGC& egC, const LFSU& lfsu, const X& x, const X& xOld, const LFSV& lfsv, const LFSUC& lfsuC,
                                         const XC& xC, const LFSVC& lfsvC, M& mat)
  {
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_volume_post_skeleton(const LA & la,
                                            const EG & eg, const LFSU & lfsu, const X & x,
                                            const LFSV & lfsv, M & mat)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_skeleton(const LA & la, const IG & ig,
                                const LFSU & lfsu_s, const X & x_s, const LFSV & lfsv_s,
                                const LFSU & lfsu_n, const X & x_n, const LFSV & lfsv_n,
                                M & mat_ss, M & mat_sn, M & mat_ns, M & mat_nn)
  {
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_boundary(const LA & la, const IG & ig,
                                const LFSU & lfsu_s, const X & x_s, const LFSV & lfsv_s,
                                M & mat_ss)
  {
  }
};
template<typename LA>
struct LocalAssemblerCallSwitchRCoupling<LA, true>
{
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_volume(const LA& la, const LFSU& lfsu,
                             const LFSV& lfsv, LocalPattern& pattern)
  {
    la.pattern_volume(lfsu, lfsv, pattern);
  }

  template<typename LFSU, typename LFSV, typename LFSUC,
           typename LFSVC, typename LocalPattern>
  static void pattern_coupling(const LA& la, const LFSU& lfsu,
                               const LFSV& lfsv, const LFSUC& lfsuC,
                               const LFSVC& lfsvC, LocalPattern& pattern_coupling)
  {
    la.pattern_coupling(lfsu, lfsv, lfsuC, lfsvC,
                        pattern_coupling);
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_volume_post_skeleton(const LA& la,
                                           const LFSU& lfsu, const LFSV& lfsv, LocalPattern& pattern)
  {
    la.pattern_volume_post_skeleton(lfsu, lfsv, pattern);
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_skeleton(const LA& la, const LFSU& lfsu_s,
                               const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
                               LocalPattern& pattern_sn, LocalPattern& pattern_ns)
  {
    la.pattern_skeleton(lfsu_s, lfsv_s, lfsu_n, lfsv_n, pattern_sn,
                        pattern_ns);
  }
  template<typename LFSU, typename LFSV, typename LocalPattern>
  static void pattern_boundary(const LA& la, const LFSU& lfsu_s,
                               const LFSV& lfsv_s, LocalPattern& pattern_ss)
  {
    la.pattern_boundary(lfsu_s, lfsv_s, pattern_ss);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_volume(const LA& la, const EG& eg,
                           const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r, int iMat)
  {
    la.alpha_volume(eg, lfsu, x, lfsv, r, iMat);
  }

  template<typename IG, typename EGP, typename EGC, typename LFSU, typename LFSV, typename X, typename R, typename LFSUC,
           typename XC, typename LFSVC>
  static void boundary_coupling(const LA& la, const IG& ig,
                                const EGC& egC, const EGP& egPair, const LFSU& lfsu, const LFSU& lfsu_pair,
                                const X& x, const X& xPair, const X& xOld,const X& xPairOld, const LFSV& lfsv,
                                const LFSV& lfsv_pair, const LFSUC& lfsuC,
                                const XC& xC, const LFSVC& lfsvC, R& r)
  {
    la.boundary_coupling(ig, egC, egPair, lfsu, lfsu_pair, x, xPair, xPairOld, xOld, lfsv, lfsv_pair, r, lfsuC,
                         xC, lfsvC);
  }

  template<typename EG, typename IGC, typename LFSU,
           typename X, typename LFSV, typename R, typename LFSUC,
           typename XC, typename LFSVC, typename FH>
  static void boundary_coupling_reversed(const LA& la, const EG& eg,
                                         const IGC& igC, const LFSU& lfsu,
                                         const X& x, const LFSV& lfsv, R& r, const LFSUC& lfsuC,
                                         const XC& xC, const LFSVC& lfsvC, std::vector<FH>& fractureHeight)
  {
    la.boundary_coupling_reversed(eg, igC, lfsu, x, lfsv, r, lfsuC,
                                  xC, lfsvC, fractureHeight);
  }
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_volume_post_skeleton(const LA& la, const EG& eg,
                                         const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
  {
    la.alpha_volume_post_skeleton(eg, lfsu, x, lfsv, r);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_skeleton(const LA& la, const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                             R& r_s, R& r_n)
  {
    la.alpha_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n,
                      r_s, r_n);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename R>
  static void alpha_boundary(const LA& la, const IG& ig,
                             const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                             R& r_s)
  {
    la.alpha_boundary(ig, lfsu_s, x_s, lfsv_s, r_s);
  }

  template<typename EG, typename LFSV, typename R>
  static void lambda_volume(const LA& la, const EG& eg,
                            const LFSV& lfsv, R& r)
  {
    la.lambda_volume(eg, lfsv, r);
  }
  template<typename EG, typename LFSV, typename R>
  static void lambda_volume_post_skeleton(const LA& la, const EG& eg,
                                          const LFSV& lfsv, R& r)
  {
    la.lambda_volume_post_skeleton(eg, lfsv, r);
  }
  template<typename IG, typename LFSV, typename R>
  static void lambda_skeleton(const LA& la, const IG& ig,
                              const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n)
  {
    la.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
  }
  template<typename IG, typename LFSV, typename R>
  static void lambda_boundary(const LA& la, const IG& ig,
                              const LFSV& lfsv, R& r)
  {
    la.lambda_boundary(ig, lfsv, r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_volume(const LA& la, const EG& eg,
                                    const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
  {
    static_assert(
      models<impl::HasJacobianApplyVolume<EG, LFSU, X, LFSV, Y&>, LA>(),
      "Your local operator does not implement jacobian_apply_volume() for linear problems (without explicit Jacobian evaluation point)"
    );
    la.jacobian_apply_volume(eg, lfsu, x, lfsv, y);
  }

  // template<typename EG, typename EGC, typename LFSU, typename X,
  //          typename LFSV, typename Y>
  // static void jacobian_apply_coupling(const LA& la, const EG& eg,
  //                                     const EGC& egC, const LFSU& lfsu, const X& x,
  //                                     const LFSV& lfsv, Y& y)
  // {
  //   static_assert(
  //     models<impl::HasJacobianApplyCoupling<EG, EGC, LFSU, X, LFSV, Y&>, LA>(),
  //     "Your local operator does not implement jacobian_apply_volume() for linear problems (without explicit Jacobian evaluation point)"
  //   );
  //   la.jacobian_apply_volume(eg, lfsu, x, lfsv, y);
  // }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_volume_post_skeleton(const LA& la,
                                                  const EG& eg, const LFSU& lfsu, const X& x,
                                                  const LFSV& lfsv, Y& y)
  {
    static_assert(
      models<impl::HasJacobianApplyVolumePostSkeleton<EG, LFSU, X, LFSV, Y&>, LA>(),
      "Your local operator does not implement jacobian_apply_volume_post_skeleton() for linear problems (without explicit Jacobian evaluation point)"
    );
    la.jacobian_apply_volume_post_skeleton(eg, lfsu, x, lfsv, y);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_skeleton(const LA& la, const IG& ig,
                                      const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                      const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                      Y& y_s, Y& y_n)
  {
    static_assert(
      models<impl::HasJacobianApplySkeleton<IG, LFSU, X, LFSV, LFSU, X, LFSV, Y&, Y&>, LA>(),
      "Your local operator does not implement jacobian_apply_skeleton() for linear problems (without explicit Jacobian evaluation point)"
    );
    la.jacobian_apply_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n,
                               lfsv_n, y_s, y_n);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void jacobian_apply_boundary(const LA& la, const IG& ig,
                                      const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                      Y& y_s)
  {
    static_assert(
      models<impl::HasJacobianApplyBoundary<IG, LFSU, X, LFSV, Y&>, LA>(),
      "Your local operator does not implement jacobian_apply_boundary() for linear problems (without explicit Jacobian evaluation point)"
    );
    la.jacobian_apply_boundary(ig, lfsu_s, x_s, lfsv_s, y_s);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static auto nonlinear_jacobian_apply_volume(const LA& la,
                                              const EG& eg, const LFSU& lfsu, const X& x, const X& z,
                                              const LFSV& lfsv,
                                              Y& y)
  -> typename std::enable_if<models<impl::HasJacobianApplyVolume<EG, LFSU, X, X, LFSV, Y&>, LA>()>::type
  {
    la.jacobian_apply_volume(eg, lfsu, x, z, lfsv, y);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static auto nonlinear_jacobian_apply_volume(const LA& la,
                                              const EG& eg, const LFSU& lfsu, const X& x, const X& z,
                                              const LFSV& lfsv,
                                              Y& y)
  -> typename std::enable_if<not models<impl::HasJacobianApplyVolume<EG, LFSU, X, X, LFSV, Y&>, LA>()>::type
  {
    static_assert(AlwaysFalse<EG>::value, "Your local operator does not implement jacobian_apply_volume() for nonlinear problems (Jacobian evaluation point given as argument)");
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_volume_post_skeleton(
    const LA& la, const EG& eg, const LFSU& lfsu, const X& x,
    const X& z, const LFSV& lfsv, Y& y)
  {
    la.jacobian_apply_volume_post_skeleton(eg, lfsu, x, z, lfsv, y);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_skeleton(const LA& la,
                                                const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                                const X& z_s, const LFSV& lfsv_s, const LFSU& lfsu_n,
                                                const X& x_n, const X& z_n, const LFSV& lfsv_n, Y& y_s,
                                                Y& y_n)
  {
    la.jacobian_apply_skeleton(ig, lfsu_s, x_s, z_s, lfsv_s, lfsu_n,
                               x_n, z_n, lfsv_n, y_s, y_n);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  static void nonlinear_jacobian_apply_boundary(const LA& la,
                                                const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                                const X& z_s, const LFSV& lfsv_s, Y& y_s)
  {
    la.jacobian_apply_boundary(ig, lfsu_s, x_s, z_s, lfsv_s, y_s);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_volume(const LA& la, const EG& eg,
                              const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat)
  {
    la.jacobian_volume(eg, lfsu, x, lfsv, mat);
  }
  // volume integral depending on test and ansatz functions
  template< typename EG, typename IGC, typename LFSU,
            typename X, typename LFSV, typename LFSUC,
            typename XC, typename LFSVC, typename FH, typename M>
  static void jacobian_boundary_coupling_reversed(const LA& la, const EG& eg, const IGC& igC, const LFSU& lfsu,
                                                  const X& x, const LFSV& lfsv, const LFSUC& lfsuC,
                                                  const XC& xC, const LFSVC& lfsvC, std::vector<FH>& fractureHeight, M& mat)
  {
    la.jacobian_boundary_coupling_reversed(eg, igC, lfsu, x, lfsv, lfsuC, xC, lfsvC, fractureHeight,  mat);
  }
  // volume integral depending on test and ansatz functions
  template< typename IG, typename EGC, typename LFSU,
            typename X, typename LFSV, typename LFSUC,
            typename XC, typename LFSVC, typename M>
  static void jacobian_boundary_coupling(const LA& la, const IG& ig, const EGC& egC, const LFSU& lfsu,
                                         const X& x, const X& xOld, const LFSV& lfsv, const LFSUC& lfsuC,
                                         const XC& xC, const LFSVC& lfsvC, M& mat)
  {
    la.jacobian_boundary_coupling(ig, egC, lfsu, x, xOld, lfsv, lfsuC, xC, lfsvC, mat);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_volume_post_skeleton(const LA& la,
                                            const EG& eg, const LFSU& lfsu, const X& x,
                                            const LFSV& lfsv, M & mat)
  {
    la.jacobian_volume_post_skeleton(eg, lfsu, x, lfsv, mat);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_skeleton(const LA& la, const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                M & mat_ss, M & mat_sn, M & mat_ns, M & mat_nn)
  {
    la.jacobian_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n,
                         lfsv_n, mat_ss, mat_sn, mat_ns, mat_nn);
  }
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename M>
  static void jacobian_boundary(const LA& la, const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                M & mat_ss)
  {
    la.jacobian_boundary(ig, lfsu_s, x_s, lfsv_s, mat_ss);
  }

};

} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_CALLSWITCH_HH
