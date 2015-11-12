PROJECT(Poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

IF(NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
ENDIF(NOT IFEM_CONFIGURED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

set(Poisson_SRCS Poisson.C
                 PoissonSolutions.C
                 SIMPoisson1D.C
                 SIMPoisson2D.C
                 SIMPoisson3D.C)

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR})

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonPoisson STATIC ${Poisson_SRCS})

ADD_EXECUTABLE(Poisson main_Poisson.C)
TARGET_LINK_LIBRARIES(Poisson CommonPoisson CommonIFEM ${IFEM_LIBRARIES})

# Installation
INSTALL(TARGETS Poisson DESTINATION bin COMPONENT bin)

ENABLE_TESTING()

# Regression tests
IF(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
ELSE(IFEM_USE_PARALLEL_PETSC)
  FILE(GLOB POISSON_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                       ${PROJECT_SOURCE_DIR}/Test/*.reg)
  FOREACH(TESTFILE ${POISSON_TESTFILES})
    IFEM_ADD_TEST(${TESTFILE} Poisson)
  ENDFOREACH()
ENDIF(IFEM_USE_PARALLEL_PETSC)
list(APPEND TEST_APPS Poisson)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  Poisson
                  CommonPoisson CommonIFEM ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
