// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_PATTERN_RCOUPLING_HH
#define DUNE_PDELAB_PATTERN_RCOUPLING_HH

#include<dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

namespace Dune
{
    namespace PDELab
    {
        class FullCouplingPattern
        {
        public:

            // define sparsity pattern of operator representation
            template<typename LFSU, typename LFSV,typename LFSUSL, typename LFSVSL, typename LocalPattern>
                      void pattern_coupling(const LFSU& lfsu, const LFSV& lfsv,const LFSUSL& lfsuSl, const LFSVSL& lfsvSl,
                              LocalPattern& pattern_coupling) const
            {
                for (size_t i = 0; i < lfsv.size(); ++i)
                    for (size_t j = 0; j < lfsu.size(); ++j)
                        pattern_coupling.addLink(lfsv, i, lfsu, j);
                // for (size_t i = 0; i < lfsvSl.size(); ++i)
                //     for (size_t j = 0; j < lfsuSl.size(); ++j)
                //         pattern_coupling.addLink(lfsvSl, i, lfsuSl, j);
            }
        };
    //! \} group GridFunctionSpace
    }// namespace PDELab
} // namespace Dune

#endif
