project(Poisson)

cmake_minimum_required(VERSION 2.6)

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
                 PoissonArgs.C
                 PoissonSolutions.C
                 SIMPoisson1D.C
                 SIMPoisson2D.C
                 SIMPoisson3D.C)

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonPoisson STATIC ${Poisson_SRCS})
list(APPEND CHECK_SOURCES ${Poisson_SRCS} main_Poisson.C)

enable_testing()
include(IFEMTesting)

# Regression tests
if(MPI_FOUND)
  ifem_add_test(MPI/Cube-subdivision.reg Poisson 4)
endif()
IF(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Cube.reg Poisson)
  ifem_add_test(Cube-subdivision.reg Poisson)
  ifem_add_test(Line.reg Poisson)
  ifem_add_test(Line-subdivision.reg Poisson)
  ifem_add_test(Lshape.reg Poisson)
  ifem_add_test(Square.reg Poisson)
  ifem_add_test(Square-subdivision-nxisnotny.reg Poisson)
  ifem_add_test(Square-subdivision-pxisnotpy.reg Poisson)
  ifem_add_test(Square-subdivision-nonhom.reg Poisson)
  ifem_add_test(Waterfall3D.reg Poisson)
  ifem_add_test(exact_p1_2D.reg Poisson)
  ifem_add_test(exact_p2_2D.reg Poisson)
  ifem_add_test(exact_p3_2D.reg Poisson)
  if(LRSpline_FOUND OR LRSPLINE_FOUND)
    ifem_add_test(interior_layer_adap.reg Poisson)
    ifem_add_test(CubeLR.reg Poisson)
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
  endif()
endif()
list(APPEND TEST_APPS Poisson)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  Poisson
                  CommonPoisson IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
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
