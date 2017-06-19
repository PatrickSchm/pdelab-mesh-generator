#include <iostream>
#include <fstream>

using namespace std;

namespace Dune
{
namespace SubIntegration
{
template<class SubMesh, int DimHigh, int DimLow>
class SubIntegrationOutput
{
private:
    const SubMesh& subMesh;

public:
//            using SubIntg = typename Dune::SubIntegration::Integrationpoints;

    SubIntegrationOutput(const SubMesh& subMesh_) :
        subMesh(subMesh_)
    {
    }

    void MatlabGauss() const
    {
        for (int iel = 0; iel < subMesh.getSubEleSize(); iel++)
        {
            auto activeElement = subMesh.getSubEle(iel);
            if (activeElement->getId() == 82 || activeElement->getId() == 3222)
            {
                for (int isub = 0;
                        isub < activeElement->getLowDimIdSize(); isub++)
                {
                std::cout << "print" << std::endl;
                    ofstream gaussfile;
                    std::string filePath;
                    std::string file =
                        "/scratch/local/schmidt/documents/tracksubmeshes/gauss_element_";
                    std::ostringstream stringStream;
                    stringStream << file << activeElement->getId()
                                 << "_subelement_" << isub << ".m";
                    filePath = stringStream.str();
                    gaussfile.open(filePath);

                    ofstream coordfile;
                    std::string filePathcoord;
                    std::string filecoord =
                        "/scratch/local/schmidt/documents/tracksubmeshes/coord_element_";
                    std::ostringstream stringStreamcoord;
                    stringStreamcoord << filecoord
                                      << activeElement->getId() << "_subelement_"
                                      << isub << ".m";
                    filePathcoord = stringStreamcoord.str();
                    coordfile.open(filePathcoord);

                    std::vector<int>* intPointIndices = activeElement->getIntPointIndices(isub);

                    for (int iint = 0;
                            iint < activeElement->getIntPointSize();
                            iint++)
                    {

                        std::vector<double>*  intCoord =
                            activeElement->getIntPointPosition(
                                iint);
                        int intGauss = 0.0;
                        if (std::find(intPointIndices->begin(), intPointIndices->end(), iint) != intPointIndices ->end())
                        {
                            int iintgauss = find(
                            (*intPointIndices).begin(),
                            (*intPointIndices).end(),
                            iint)
                        - (*intPointIndices).begin();

                            intGauss =
                                activeElement->getIntPointGauss(isub,
                                                                iintgauss);
                        }
                        if (iint == 0)
                        {
                            coordfile << "Coord" << std::endl
                                      << (*intCoord)[0] << " " << (*intCoord)[1]
                                      << " " << (*intCoord)[2]
                                      << std::endl;

                            gaussfile << "Gauss" << std::endl << intGauss
                                      << std::endl ;
                        } else if (iint
                                   == activeElement->getIntPointSize() - 1)
                        {
                            coordfile << (*intCoord)[0] << " "
                                      << (*intCoord)[1] << " " << (*intCoord)[2]
                                      << std::endl;

                            gaussfile << intGauss << "];" << std::endl;
                        } else
                        {
                            coordfile << (*intCoord)[0] << " "
                                      << (*intCoord)[1] << " " << (*intCoord)[2]
                                      << std::endl;

                            gaussfile << intGauss << std::endl ;
                        }

                    }
                    gaussfile.close();
                    coordfile.close();
                }
            }
        }
    }
};
}
}
