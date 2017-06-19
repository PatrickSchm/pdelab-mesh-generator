
#ifndef MATERIALREADER_HH
#define MATERIALREADER_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace Dune {
namespace PDELab {

template<typename Grid>
class MaterialFlags {
private:
    int nEle;
    std::string activate = "$Elements";
    std::string deactivate = "$EndElements";
    std::vector<int> materialFlags;
    bool elementFlag;
public:
    const int getMaterial(int iEle) const
    {
        return this->materialFlags[iEle];
    }

    // void setMaterial(std::vector<int> input, typename Dune::GridFactory<Grid> factory, typename Grid::LeafGridView gv)
    // {
    //     typedef typename Grid::LeafGridView  GV;
    //     typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    //     ElementIterator it = gv.template begin<0>();
    //     ElementIterator itEnd = gv.template end<0>();

    //     for (it; it != itEnd; ++it)
    //     {
    //         materialFlags.push_back(input[factory.insertionIndex(it)]);
    //     }
    // }
};
}
}
#endif