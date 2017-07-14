namespace Dune {
namespace FractureCoupling {

template<typename Grid0, typename Grid1>
class Pairs {

    using highElement = typename Grid0::LeafGridView::template Codim<0>::Entity;
    using lowElement = typename  Grid1::LeafGridView::template Codim<0>::Entity;
public:

    std::vector<std::vector<highElement>> bulkCoupledElements;
    std::vector<std::vector<lowElement>> bulkFractureCoupledElements;
    std::vector<std::vector<highElement>> fractureBulkCoupledElements;
private:
    const Grid0& grid0_;
    const Grid1& grid1_;

public:
    Pairs(const Grid0& grid0, const Grid1& grid1): grid0_(grid0), grid1_(grid1) {}

    void setBulkPairs(std::vector<int> boundary_id_to_physical_entity, std::map<int, int> bulkPairs) {
        auto gv0 = grid0_.leafGridView();
        bulkCoupledElements.resize(gv0.size(0));
        const auto & indexSetHigh = gv0.indexSet();

        for (const auto &element : elements(gv0)) {
            for (const auto& intersection : intersections(gv0, element))
            {
                if (intersection.boundary()) {

                    int physicalTag = boundary_id_to_physical_entity[intersection.boundarySegmentIndex() - 1];
                    if (physicalTag >= 10000) {

                        int idx = indexSetHigh.index(element);
                        bool pairFound = false;
                        std::vector<highElement> element_dummy;

                        for (const auto &elementPair : elements(gv0)) {
                            int idxPair = indexSetHigh.index(elementPair);

                            for (const auto& intersectionPair : intersections(gv0, elementPair))
                            {
                                if (intersectionPair.boundary()) {
                                    if (boundary_id_to_physical_entity[intersectionPair.boundarySegmentIndex() - 1] == bulkPairs[physicalTag]) {
                                        auto activeElementHigh = gv0.template begin<0>();
                                        std::advance(activeElementHigh, idxPair);
                                        element_dummy.push_back(*activeElementHigh);


                                        pairFound = true;
                                        break;
                                    }
                                }
                            }
                            if (pairFound)
                                break;
                        }
                        bulkCoupledElements[idx] = element_dummy;
                    }
                }
            }
        }
    }

    void setBulkFracturePairs(std::vector<int> boundary_id_to_physical_entity, std::vector<int> element_id_to_physical_entity_low, std::map<int, std::vector<int>> bulkFracturePairs) {
        auto gv0 = grid0_.leafGridView();
        auto gv1 = grid1_.leafGridView();

        const auto & indexSetHigh = gv0.indexSet();
        const auto & indexSetLow = gv1.indexSet();
        bulkFractureCoupledElements.resize(gv0.size(0));
        fractureBulkCoupledElements.resize(gv1.size(0));
        for (const auto &element : elements(gv0)) {
            for (const auto& intersection : intersections(gv0, element))
            {
                if (intersection.boundary()) {
                    int idx = indexSetHigh.index(element);

                    int physicalTag = boundary_id_to_physical_entity[intersection.boundarySegmentIndex() - 1];
                    if (physicalTag >= 10000) {
                        bool pairFound = false;
                        auto bulkFractureElements = bulkFracturePairs[physicalTag];
                        std::vector<lowElement> element_dummy;
                        for (const auto &elementLow : elements(gv1)) {
                            int idxLow = indexSetLow.index(elementLow);
                            std::vector<highElement> element_dummy_high;

                            int physicalTagLow = element_id_to_physical_entity_low[idxLow];
                            if (std::find(bulkFractureElements.begin(), bulkFractureElements.end(), physicalTagLow) != bulkFractureElements.end()) {
                                auto activeElementLow = gv1.template begin<0>();
                                std::advance(activeElementLow, idxLow);

                                element_dummy.push_back(*activeElementLow);
                                element_dummy_high.push_back(element);
                                pairFound = true;
                            }
                            if (pairFound)
                                fractureBulkCoupledElements[idxLow] = element_dummy_high;

                        }
                        bulkFractureCoupledElements[idx] = element_dummy;
                    }
                }
            }
        }
    }
};
}
}