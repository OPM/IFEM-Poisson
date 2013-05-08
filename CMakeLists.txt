PROJECT(Poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

# To allow standalone compilation out of the main tree
IF(NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
ENDIF(NOT IFEM_CONFIGURED)

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../Common/scripts/doregtests.sh.in
  doregtests.sh @ONLY)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../Common/scripts/check-patchseries.sh.in
               check-patchseries.sh)

FILE(GLOB Poisson_SRCS *.C)

INCLUDE_DIRECTORIES(${IFEM_INCLUDES})

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

ADD_EXECUTABLE(Poisson ${Poisson_SRCS})
TARGET_LINK_LIBRARIES(Poisson ${IFEM_LIBRARIES})

# Installation
INSTALL(TARGETS Poisson DESTINATION bin COMPONENT bin)

ENABLE_TESTING()

# Regression tests
IF(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
ELSE(IFEM_USE_PARALLEL_PETSC)
  FILE(GLOB POISSON_TESTFILES "${PROJECT_SOURCE_DIR}/Test/*.reg")
  FOREACH(TESTFILE ${POISSON_TESTFILES})
    ADD_TEST(${TESTFILE} regtest.sh "${EXECUTABLE_OUTPUT_PATH}/Poisson" "${TESTFILE}")
  ENDFOREACH()
ENDIF(IFEM_USE_PARALLEL_PETSC)
