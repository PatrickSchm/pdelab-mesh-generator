exclude_from_headercheck(assembler_rcoupling.hh)

install(FILES     assembler_rcoupling.hh
        assemblerutilities_rcoupling.hh
        gridoperatorutilities_coupling_rcoupling.hh
        localassemblerenginebase_rcoupling.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/residualcoupling/gridoperator/common)
