cmake_minimum_required(VERSION 3.22)

project(Poisson LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

ifem_add_library(
  NAME
    CommonPoisson
  SOURCES
    Poisson.C
    PoissonSolutions.C
    PoissonSource.C
    SIMPoisson.C
  HEADERS
    Poisson.h
    PoissonSolutions.h
    PoissonSource.h
    SIMPoisson.h
  LIBRARIES
    IFEMAppCommon
    IFEM
)

ifem_add_application(
  NAME
    Poisson
  SOURCES
    main_Poisson.C
  LIBRARIES
    CommonPoisson
)

enable_testing()

# Regression tests
ifem_add_hdf5_test(
  TARGET
    Poisson
  TEST_FILES
    MPI/Square-2-LR.hreg
  PARALLEL
    2
)

ifem_add_hdf5_test(
  TARGET
    Poisson
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    interior_layer_adap_res.hreg
    interior_layer_adap_rec.hreg
)

ifem_add_hdf5_test(
  TARGET
    Poisson
  TEST_FILES
    Square.hreg
    Square-eig1.hreg
    Square-eig4.hreg
)

ifem_add_vtf_test(
  TARGET
    Poisson
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    round_stair_3patch_3D_LR_adap.vreg
    interior_layer_adap_res.vreg
)

ifem_add_vtf_test(
  TARGET
    Poisson
  TEST_FILES
    Square.vreg
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    Zoltan_FOUND
  TEST_FILES
    MPI/Cube.reg
    MPI/Square.reg
    MPI/Square-eig4.reg
    MPI/CGL2-petsc.reg
  PARALLEL
    4
)

ifem_add_regression_test(
  TARGET
    Poisson
  TEST_FILES
    MPI/Square-2-LR.reg
  PARALLEL
    2
)

ifem_add_regression_test(
  TARGET
    Poisson
  TEST_FILES
    MPI/Cube-3-LR.reg
  PARALLEL
    3
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    SLEPc_FOUND
    Zoltan_FOUND
  TEST_FILES
    MPI/Square-eig13-slepc.reg
  PARALLEL
    3
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    LRSpline_FOUND
    Zoltan_FOUND
  TEST_FILES
    MPI/interior_layer_adap_dorfel.reg
    MPI/Waterfall3D_adap.reg
  PARALLEL
    4
)

ifem_add_regression_test(
  TARGET
    Poisson
  TEST_FILES
    Cube.reg
    CGL2.reg
    Cube-lag.reg
    Line.reg
    Lshape.reg
    Square.reg
    Square-robin.reg
    Square-lag.reg
    Cube-robin.reg
    Waterfall3D.reg
    exact_p1_2D.reg
    exact_p2_2D.reg
    exact_p3_2D.reg
    Sinus2D.reg
    Annulus_sep_proj.reg
    Square-neumann.reg
    Square-neumann-lag.reg
    Square-eig4.reg
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    PETSc_FOUND
  TEST_FILES
    Square-eig4-petsc.reg
    CGL2-petsc.reg
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    SLEPc_FOUND
  TEST_FILES
    Square-eig13-slepc.reg
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    interior_layer_adap.reg
    interior_layer_adap_rec.reg
    interior_layer_adap_res.reg
    interior_layer_adap_dorfel.reg
    interior_layer_adap_symmetrized.reg
    CubeLR.reg
    CGL2-LR.reg
    LshapeLR.reg
    Waterfall3DLR.reg
    Waterfall3D_adap.reg
    exact_p1_2D_LR.reg
    exact_p2_2D_LR.reg
    exact_p3_2D_LR.reg
    Square-2-LR.reg
    Square-2-LR-r.reg
    Backstep-3-LR.reg
    Cube-2-LR.reg
    Cube-3-LR.reg
    round_stair_3patch_3D_LR_adap.reg
    lshape-r2-3patch-LR-adap.reg
    Annulus_sep_proj_LR.reg
)

ifem_add_regression_test(
  TARGET
    Poisson
  CONDITIONS
    HDF5_FOUND
  TEST_FILES
    Cube-mat.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    Poisson
  SOURCES
    Test/TestSIMPoisson.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
  LIBRARIES
    CommonPoisson
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET Poisson)

# Installation
include(GNUInstallDirs)
install(TARGETS Poisson DESTINATION ${CMAKE_INSTALL_BINDIR})
