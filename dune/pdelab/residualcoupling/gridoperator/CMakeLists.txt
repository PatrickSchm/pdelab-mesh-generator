add_subdirectory(common)
add_subdirectory(default)
add_subdirectory(onestep)

install(FILES    gridoperator_coupling_rcoupling.hh
        onestep_rcoupling.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/residualcoupling/gridoperator)
