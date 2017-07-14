#ifndef DUNE_PDELAB_DEFAULT_ASSEMBLER_RCOUPLING_HH
#define DUNE_PDELAB_DEFAULT_ASSEMBLER_RCOUPLING_HH

#include <dune/common/typetraits.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/residualcoupling/gridoperator/common/assemblerutilities_rcoupling.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/residualcoupling/submesh/index_determiner.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/common/iteratorrange.hh>

extern bool print;

namespace Dune
{
namespace PDELab
{

/**
 \brief The assembler for standard DUNE grid

 * \tparam GFSU GridFunctionSpace for ansatz functions
 * \tparam GFSV GridFunctionSpace for test functions
 * \tparam nonoverlapping_mode Indicates whether assembling is done for overlap cells
 */
template<typename ...> class DefaultAssemblerCoupling;

template<typename GFSUMA, typename GFSVMA, typename GFSUC,
         typename GFSVC, typename CU, typename CV, typename CUC,
         typename CVC, typename FP, typename MaterialFlags, typename InPara>
class DefaultAssemblerCoupling<GFSUMA, GFSVMA, GFSUC, GFSVC, CU, CV,
          CUC, CVC, FP, MaterialFlags, InPara>
{
public:

    //! Types related to current grid view
    //! @{
    using EntitySet = typename GFSUMA::Traits::EntitySet;
    using EntitySetC = typename GFSUC::Traits::EntitySet;
    using Element = typename EntitySet::Element;
    using ElementC = typename EntitySetC::Element;
    using Intersection = typename EntitySet::Intersection;
    using IntersectionC = typename EntitySetC::Intersection;

    //! @}
    typedef typename GFSUMA::Traits::GridViewType GV0;
    typedef typename GFSUC::Traits::GridViewType GV1;

    typedef typename GV0::Traits::template Codim<0>::Iterator ElementIteratorMa;
    typedef typename GV1::Traits::template Codim<0>::Iterator ElementIteratorC;
    //! Grid function spaces
    //! @{
    typedef GFSUMA TrialGridFunctionSpace;
    typedef GFSVMA TestGridFunctionSpace;
    typedef GFSUC TrialGridFunctionSpaceC;
    typedef GFSVC TestGridFunctionSpaceC;
    //! @}

    //! Size type as used in grid function space
    typedef typename GFSUMA::Traits::SizeType SizeType;
    typedef typename GFSUC::Traits::SizeType SizeTypeC;

    //! Static check on whether this is a Galerkin method
    static const bool isGalerkinMethod =
        std::is_same<GFSUMA, GFSVMA>::value;

    DefaultAssemblerCoupling(const GFSUMA& gfsuMa_,
                             const GFSVMA& gfsvMa_, const GFSUC& gfsuC_,
                             const GFSVC& gfsvC_, const FP& fracturePairs_, const MaterialFlags & materialFlags_,
                             const InPara & indicesParallel_, const CU& cu_,
                             const CV& cv_, const CUC& cuC_, const CVC& cvC_) :
        gfsuMa(gfsuMa_), gfsvMa(gfsvMa_), gfsuC(gfsuC_), gfsvC(
            gfsvC_), fracturePairs(fracturePairs_), materialFlags(
                materialFlags_), indicesParallel(indicesParallel_), cu(cu_), cv(cv_), cuC(cuC_), cvC(
                    cvC_), lfsuMa(gfsuMa_), lfsvMa(gfsvMa_), lfsuPair(gfsuMa_), lfsvPair(gfsvMa_), lfsunMa(
                        gfsuMa_), lfsvnMa(gfsvMa_), lfsuC(gfsuC_), lfsvC(
                            gfsvC_), lfsunC(gfsuC_), lfsvnC(gfsvC_)
    {
    }

    DefaultAssemblerCoupling(const GFSUMA& gfsuMa_,
                             const GFSVMA& gfsvMa_, const CU& cu_, const CV& cv_) :
        gfsuMa(gfsuMa_), gfsvMa(gfsvMa_), gfsuC(), gfsvC(), cu(
            cu_), cv(cv_), cuC(), cvC(), fracturePairs(), materialFlags(), indicesParallel(), lfsuMa(
                gfsuMa_), lfsvMa(gfsvMa_), lfsuPair(gfsuMa_), lfsvPair(gfsvMa_), lfsunMa(gfsuMa_), lfsvnMa(
                    gfsvMa_)
    {
    }

    //! Get the trial grid function space
    const GFSUMA& trialGridFunctionSpaceMa() const
    {
        return gfsuMa;
    }

    //! Get the test grid function space
    const GFSVMA& testGridFunctionSpaceMa() const
    {
        return gfsvMa;
    }

    //! Get the trial grid function space
    const GFSUC& trialGridFunctionSpaceSl() const
    {
        return gfsuC;
    }

    //! Get the test grid function space
    const GFSVMA& testGridFunctionSpaceSl() const
    {
        return gfsvC;
    }

    // Assembler (const GFSU& gfsu_, const GFSV& gfsv_)
    //   : gfsu(gfsu_), gfsv(gfsv_), lfsu(gfsu_), lfsv(gfsv_),
    //     lfsun(gfsu_), lfsvn(gfsv_),
    //     sub_triangulation(ST(gfsu_.gridview(),Dune::PDELab::NoSubTriangulationImp()))
    // { }

    template<class LocalAssemblerEngine>
    void assemble(LocalAssemblerEngine & assembler_engine) const
    {
        using lowElement = typename GV1::Traits::template Codim<0>::Entity;
        typedef std::map<int, std::vector<lowElement>> BFMAP;
        typedef LFSIndexCache<LFSUMA, CU> LFSUMACache;
        typedef LFSIndexCache<LFSVMA, CV> LFSVMACache;

        typedef LFSIndexCache<LFSUC, CUC> LFSUCCache;
        typedef LFSIndexCache<LFSVC, CVC> LFSVCCache;

        const bool needs_constraints_caching =
            assembler_engine.needsConstraintsCaching(cu, cv);

        const bool needs_constraints_cachingC =
            assembler_engine.needsConstraintsCaching(cuC, cvC);

        LFSUMACache lfsu_cacheMa(lfsuMa, cu, needs_constraints_caching);
        LFSVMACache lfsv_cacheMa(lfsvMa, cv, needs_constraints_caching);

        LFSUMACache lfsu_cacheMaOld(lfsuMa, cu, needs_constraints_caching);
        LFSVMACache lfsv_cacheMaOld(lfsvMa, cv, needs_constraints_caching);

        LFSUMACache lfsu_cachePair(lfsuPair, cu, needs_constraints_caching);
        LFSVMACache lfsv_cachePair(lfsvPair, cv, needs_constraints_caching);

        LFSUMACache lfsu_cachePairOld(lfsuPair, cu, needs_constraints_caching);
        LFSVMACache lfsv_cachePairOld(lfsvPair, cv, needs_constraints_caching);


        LFSUMACache lfsun_cacheMa(lfsunMa, cu,
                                  needs_constraints_caching);
        LFSVMACache lfsvn_cacheMa(lfsvnMa, cv,
                                  needs_constraints_caching);

        LFSUCCache lfsu_cacheC(lfsuC, cuC,
                               needs_constraints_cachingC);
        LFSVCCache lfsv_cacheC(lfsvC, cvC,
                               needs_constraints_cachingC);
        LFSUCCache lfsun_cacheC(lfsunC, cuC,
                                needs_constraints_cachingC);
        LFSVCCache lfsvn_cacheC(lfsvnC, cvC,
                                needs_constraints_cachingC);

        // Notify assembler engine about oncoming assembly
        assembler_engine.preAssembly();

        // Extract integration requirements from the local assembler
        const bool require_uv_skeleton =
            assembler_engine.requireUVSkeleton();
        const bool require_v_skeleton =
            assembler_engine.requireVSkeleton();
        const bool require_uv_boundary =
            assembler_engine.requireUVBoundary();
        const bool require_v_boundary =
            assembler_engine.requireVBoundary();
        const bool require_uv_processor =
            assembler_engine.requireUVBoundary();
        const bool require_v_processor =
            assembler_engine.requireVBoundary();
        const bool require_uv_post_skeleton =
            assembler_engine.requireUVVolumePostSkeleton();
        const bool require_v_post_skeleton =
            assembler_engine.requireVVolumePostSkeleton();
        const bool require_skeleton_two_sided =
            assembler_engine.requireSkeletonTwoSided();
        const bool require_v_coupling =
            assembler_engine.requireCoupling();
        bool reverse = false;
        if (GV0::dimension < GV1::dimension)
            reverse = true;

        auto entity_set = gfsuMa.entitySet();
        auto& index_set = entity_set.indexSet();

        auto entity_setC = gfsuC.entitySet();
        auto& index_setC = entity_setC.indexSet();

// Traverse grid view
        for (const auto& element : elements(entity_set))
        {
            assembler_engine.requireUVSkeleton();
            // Compute unique id
            auto ids = index_set.uniqueIndex(element);
            int idsindex = index_set.index(element);
            // std::cout << "ELEMENT: " << idsindex << std::endl;
            // int iMat = materialFlags.getMaterial(idsindex);

            int iMat = 0;

            ElementGeometry<Element> eg(element);

            if (assembler_engine.assembleCell(eg))
                continue;

            // Bind local test function space to element
            lfsvMa.bind( element );
            lfsv_cacheMa.update();
            lfsv_cacheMaOld.update();
            // Notify assembler engine about bind

            // Notify assembler engine about bind
            assembler_engine.onBindLFSV(eg, lfsv_cacheMa);

            // Volume integration
            assembler_engine.assembleVVolume(eg, lfsv_cacheMa);

            // Bind local trial function space to element
            lfsuMa.bind(element);
            lfsu_cacheMa.update();
            lfsu_cacheMaOld.update();


            assembler_engine.onBindLFSUV(eg, lfsu_cacheMa, lfsv_cacheMa);
            assembler_engine.onBindLFSUVOld(eg, lfsu_cacheMaOld, lfsv_cacheMaOld);

            // Load coefficients of local functions
            assembler_engine.loadCoefficientsLFSUInside(lfsu_cacheMa);
            assembler_engine.loadCoefficientsLFSUInsideOld(lfsu_cacheMaOld);

            // Volume integration
            assembler_engine.assembleUVVolume(eg, lfsu_cacheMa, lfsv_cacheMa, iMat);

            auto fractureElements = fracturePairs.bulkFractureCoupledElements[idsindex];
            auto pElement = fracturePairs.bulkCoupledElements[idsindex];

            if (fractureElements.size() > 0 && pElement.size() > 0) {

                lfsuPair.bind(pElement[0]);
                lfsu_cachePair.update();
                lfsu_cachePairOld.update();
                lfsvPair.bind(pElement[0]);
                lfsv_cachePair.update();
                lfsv_cachePairOld.update();


                ElementGeometry<Element> egPair(pElement[0]);
                assembler_engine.onBindLFSUVPair(egPair, lfsu_cachePair, lfsv_cachePair);
                assembler_engine.onBindLFSUVPairOld(egPair, lfsu_cachePairOld, lfsv_cachePairOld);
                assembler_engine.loadCoefficientsLFSUInsidePair(lfsu_cachePair);
                assembler_engine.loadCoefficientsLFSUInsidePairOld(lfsu_cachePairOld);

                auto cBegin = fracturePairs.bulkFractureCoupledElements[idsindex].begin();
                auto cEnd = fracturePairs.bulkFractureCoupledElements[idsindex].end();
                for (auto & cElement = cBegin; cElement != cEnd; cElement++) {
                    ElementGeometry<ElementC> egC(*cElement);

                    lfsuC.bind(*cElement);
                    lfsu_cacheC.update();
                    lfsvC.bind(*cElement);
                    lfsv_cacheC.update();

                    assembler_engine.onBindLFSUVSL(egC, lfsu_cacheC, lfsv_cacheC);

                    // Load coefficients of local functions
                    assembler_engine.loadCoefficientsLFSUInsideSl(lfsu_cacheC);
                    unsigned int intersection_index = 0;

                    for (const auto& intersection : intersections(entity_set, element))
                    {
                        IntersectionGeometry<Intersection> ig(intersection, intersection_index);
                        // Boundary integration
                        if (intersection.boundary()) {
                            // std::cout << "------------------------- ELEMENT: " << idsindex << std::endl;
                            assembler_engine.assembleCBoundary(ig, egC, egPair , lfsu_cacheMa, lfsv_cacheMa, lfsu_cachePair, lfsv_cachePair, lfsu_cacheC, lfsv_cacheC);
                        }
                        ++intersection_index;
                    }

                    assembler_engine.onUnbindLFSUVC(egC,
                                                    lfsu_cacheC, lfsv_cacheC);
                    // Notify assembler engine about unbinds
                    assembler_engine.onUnbindLFSVC(egC,
                                                   lfsv_cacheC);
                }
            }

            // Skip if no intersection iterator is needed
            if (require_uv_skeleton || require_v_skeleton ||
                    require_uv_boundary || require_v_boundary ||
                    require_uv_processor || require_v_processor)
            {
                // Traverse intersections
                unsigned int intersection_index = 0;
                for (const auto& intersection : intersections(entity_set, element))
                {

                    IntersectionGeometry<Intersection> ig(intersection, intersection_index);

                    auto intersection_data = classifyIntersection(entity_set, intersection);
                    auto intersection_type = std::get<0>(intersection_data);
                    auto& outside_element = std::get<1>(intersection_data);

                    switch (intersection_type)
                    {
                    case IntersectionType::skeleton:
                    // the specific ordering of the if-statements in the old code caused periodic
                    // boundary intersection to be handled the same as skeleton intersections
                    case IntersectionType::periodic:
                        if (require_uv_skeleton || require_v_skeleton)
                        {
                            // compute unique id for neighbor

                            auto idn = index_set.uniqueIndex(outside_element);

                            // Visit face if id is bigger
                            bool visit_face = ids > idn || require_skeleton_two_sided;

                            // unique vist of intersection
                            if (visit_face)
                            {
                                // Bind local test space to neighbor element
                                lfsvnMa.bind(outside_element);
                                lfsvn_cacheMa.update();

                                // Notify assembler engine about binds
                                assembler_engine.onBindLFSVOutside(ig, lfsv_cacheMa, lfsvn_cacheMa);

                                // Skeleton integration
                                assembler_engine.assembleVSkeleton(ig, lfsv_cacheMa, lfsvn_cacheMa);

                                if (require_uv_skeleton) {

                                    // Bind local trial space to neighbor element
                                    lfsunMa.bind(outside_element);
                                    lfsun_cacheMa.update();

                                    // Notify assembler engine about binds
                                    assembler_engine.onBindLFSUVOutside(ig,
                                                                        lfsu_cacheMa, lfsv_cacheMa,
                                                                        lfsun_cacheMa, lfsvn_cacheMa);

                                    // Load coefficients of local functions
                                    assembler_engine.loadCoefficientsLFSUOutside(lfsun_cacheMa);

                                    // Skeleton integration
                                    assembler_engine.assembleUVSkeleton(ig, lfsu_cacheMa, lfsv_cacheMa, lfsun_cacheMa, lfsvn_cacheMa);

                                    // Notify assembler engine about unbinds
                                    assembler_engine.onUnbindLFSUVOutside(ig,
                                                                          lfsu_cacheMa, lfsv_cacheMa,
                                                                          lfsun_cacheMa, lfsvn_cacheMa);
                                }

                                // Notify assembler engine about unbinds
                                assembler_engine.onUnbindLFSVOutside(ig, lfsv_cacheMa, lfsvn_cacheMa);
                            }
                        }
                        break;

                    case IntersectionType::boundary:
                        if (require_uv_boundary || require_v_boundary )
                        {

                            // Boundary integration
                            assembler_engine.assembleVBoundary(ig, lfsv_cacheMa);

                            if (require_uv_boundary) {
                                // Boundary integration
                                assembler_engine.assembleUVBoundary(ig, lfsu_cacheMa, lfsv_cacheMa);
                            }
                        }
                        break;

                    case IntersectionType::processor:
                        if (require_uv_processor || require_v_processor )
                        {

                            // Processor integration
                            assembler_engine.assembleVProcessor(ig, lfsv_cacheMa);

                            if (require_uv_processor) {
                                // Processor integration
                                assembler_engine.assembleUVProcessor(ig, lfsu_cacheMa, lfsv_cacheMa);
                            }
                        }
                        break;
                    } // switch

                    ++intersection_index;
                } // iit
            } // do skeleton

            if (require_uv_post_skeleton || require_v_post_skeleton) {
                // Volume integration
                assembler_engine.assembleVVolumePostSkeleton(eg, lfsv_cacheMa);

                if (require_uv_post_skeleton) {
                    // Volume integration
                    assembler_engine.assembleUVVolumePostSkeleton(eg, lfsu_cacheMa, lfsv_cacheMa);
                }
            }

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSUV(eg, lfsu_cacheMa, lfsv_cacheMa);

            // Notify assembler engine about unbinds
            assembler_engine.onUnbindLFSV(eg, lfsv_cacheMa);

        } // it
        assembler_engine.postAssembly(gfsuMa, gfsvMa);
    }
private:

    /* global function spaces */
    const GFSUMA& gfsuMa;
    const GFSVMA& gfsvMa;
    const GFSUC& gfsuC;
    const GFSVC& gfsvC;
    const FP & fracturePairs;
    const MaterialFlags& materialFlags;
    const InPara& indicesParallel;
    typename std::conditional <
    std::is_same<CU, EmptyTransformation>::value, const CU,
        const CU& >::type cu;
    typename std::conditional <
    std::is_same<CV, EmptyTransformation>::value, const CV,
        const CV& >::type cv;
    typename std::conditional <
    std::is_same<CUC, EmptyTransformation>::value, const CUC,
        const CUC& >::type cuC;
    typename std::conditional <
    std::is_same<CVC, EmptyTransformation>::value, const CVC,
        const CVC& >::type cvC;

    /* local function spaces */
    typedef LocalFunctionSpace<GFSUMA, TrialSpaceTag> LFSUMA;
    typedef LocalFunctionSpace<GFSVMA, TestSpaceTag> LFSVMA;
    typedef LocalFunctionSpace<GFSUC, TrialSpaceTag> LFSUC;
    typedef LocalFunctionSpace<GFSVC, TestSpaceTag> LFSVC;
// local function spaces in local cell
    mutable LFSUMA lfsuMa;
    mutable LFSVMA lfsvMa;
    // local function spaces in local cell
    mutable LFSUMA lfsuPair;
    mutable LFSVMA lfsvPair;
// local function spaces in neighbor
    mutable LFSUMA lfsunMa;
    mutable LFSVMA lfsvnMa;

// local function spaces in local cell
    mutable LFSUC lfsuC;
    mutable LFSVC lfsvC;
// local function spaces in neighbor
    mutable LFSUC lfsunC;
    mutable LFSVC lfsvnC;

};
}
}
#endif
