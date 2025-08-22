cmake_minimum_required(VERSION 3.5)

project(Poisson)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      ${PROJECT_SOURCE_DIR}/../../cmake/Scripts)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

set(Poisson_SRCS Poisson.C
                 PoissonSolutions.C
                 PoissonSource.C
                 SIMPoisson.C
                 Poisson.h
                 PoissonSolutions.h
                 PoissonSource.h
                 SIMPoisson.h)

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonPoisson STATIC ${Poisson_SRCS})
list(APPEND CHECK_SOURCES ${Poisson_SRCS} main_Poisson.C)

enable_testing()
include(IFEMTesting)

# Regression tests
if(MPI_FOUND)
  ifem_add_test(MPI/Cube.reg Poisson 4)
  ifem_add_test(MPI/Cube-subdivision.reg Poisson 4)
  ifem_add_test(MPI/Square.reg Poisson 4)
  ifem_add_test(MPI/Square-2-LR.reg Poisson 2)
  ifem_add_test(MPI/Square-eig4.reg Poisson 3)
  ifem_add_test(MPI/Cube-3-LR.reg Poisson 3)
  ifem_add_test(MPI/CGL2-petsc.reg Poisson 4)
  ifem_add_hdf5_test(MPI/Square-2-LR.hreg Poisson 2)
  if(SLEPC_FOUND)
    ifem_add_test(MPI/Square-eig13-slepc.reg Poisson 3)
  endif()
  if(LRSpline_FOUND OR LRSPLINE_FOUND)
    ifem_add_test(MPI/interior_layer_adap_dorfel.reg Poisson 4)
    ifem_add_test(MPI/Waterfall3D_adap.reg Poisson 4)
  endif()
endif()
IF(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Cube.reg Poisson)
  ifem_add_test(CGL2.reg Poisson)
  ifem_add_test(Cube-subdivision.reg Poisson)
  ifem_add_test(Cube-lag.reg Poisson)
  ifem_add_test(Line.reg Poisson)
  ifem_add_test(Line-subdivision.reg Poisson)
  ifem_add_test(Lshape.reg Poisson)
  ifem_add_test(Square.reg Poisson)
  ifem_add_test(Square-robin.reg Poisson)
  ifem_add_test(Square-lag.reg Poisson)
  ifem_add_test(Cube-robin.reg Poisson)
  ifem_add_test(Square-subdivision-nxisnotny.reg Poisson)
  ifem_add_test(Square-subdivision-pxisnotpy.reg Poisson)
  ifem_add_test(Square-subdivision-nonhom.reg Poisson)
  ifem_add_test(Waterfall3D.reg Poisson)
  ifem_add_test(exact_p1_2D.reg Poisson)
  ifem_add_test(exact_p2_2D.reg Poisson)
  ifem_add_test(exact_p3_2D.reg Poisson)
  ifem_add_test(Sinus2D.reg Poisson)
  ifem_add_test(Annulus_sep_proj.reg Poisson)
  ifem_add_test(Square-neumann.reg Poisson)
  ifem_add_test(Square-neumann-lag.reg Poisson)
  ifem_add_test(Square-eig4.reg Poisson)
  if(PETSC_FOUND)
    ifem_add_test(Square-eig4-petsc.reg Poisson)
    ifem_add_test(CGL2-petsc.reg Poisson)
  endif()
  if(SLEPC_FOUND)
    ifem_add_test(Square-eig13-slepc.reg Poisson)
  endif()
  if(LRSpline_FOUND OR LRSPLINE_FOUND)
    ifem_add_test(interior_layer_adap.reg Poisson)
    ifem_add_test(interior_layer_adap_rec.reg Poisson)
    ifem_add_test(interior_layer_adap_res.reg Poisson)
    ifem_add_test(interior_layer_adap_dorfel.reg Poisson)
    ifem_add_test(interior_layer_adap_symmetrized.reg Poisson)
    ifem_add_test(CubeLR.reg Poisson)
    ifem_add_test(CGL2-LR.reg Poisson)
    ifem_add_test(LshapeLR.reg Poisson)
    ifem_add_test(Waterfall3DLR.reg Poisson)
    ifem_add_test(Waterfall3D_adap.reg Poisson)
    ifem_add_test(exact_p1_2D_LR.reg Poisson)
    ifem_add_test(exact_p2_2D_LR.reg Poisson)
    ifem_add_test(exact_p3_2D_LR.reg Poisson)
    ifem_add_test(Square-2-LR.reg Poisson)
    ifem_add_test(Square-2-LR-r.reg Poisson)
    ifem_add_test(Backstep-3-LR.reg Poisson)
    ifem_add_test(Cube-2-LR.reg Poisson)
    ifem_add_test(Cube-3-LR.reg Poisson)
    ifem_add_test(round_stair_3patch_3D_LR_adap.reg Poisson)
    ifem_add_test(lshape-r2-3patch-LR-adap.reg Poisson)
    ifem_add_vtf_test(round_stair_3patch_3D_LR_adap.vreg Poisson)
    ifem_add_vtf_test(interior_layer_adap_res.vreg Poisson)
    ifem_add_hdf5_test(interior_layer_adap_res.hreg Poisson)
    ifem_add_hdf5_test(interior_layer_adap_rec.hreg Poisson)
    ifem_add_test(Annulus_sep_proj_LR.reg Poisson)
  endif()
  ifem_add_vtf_test(Square.vreg Poisson)
  ifem_add_hdf5_test(Square.hreg Poisson)
  ifem_add_hdf5_test(Square-eig1.hreg Poisson)
  ifem_add_hdf5_test(Square-eig4.hreg Poisson)
  if(HDF5_FOUND)
    ifem_add_test(Cube-mat.reg Poisson)
  endif()
endif()
list(APPEND TEST_APPS Poisson)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  Poisson 0
                  CommonPoisson IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  add_check_target()
endif()

add_executable(Poisson main_Poisson.C)
target_link_libraries(Poisson CommonPoisson IFEMAppCommon ${IFEM_LIBRARIES})

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(Poisson Poisson)

# Installation
install(TARGETS Poisson DESTINATION bin COMPONENT bin)
