PROJECT(Poisson)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# To allow standalone compilation out of the main tree
IF(NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

  IF(NOT WIN32)
    # Emit position-independent code, suitable for dynamic linking
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    # Enable all warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  ENDIF(NOT WIN32)
ENDIF(NOT IFEM_CONFIGURED)

FILE(GLOB Poisson_SRCS *.C)
FILE(GLOB Poisson_HDRS *.h;CMakeLists.txt)
FILE(GLOB Poisson_TSTS Test/*.xinp;Test/*.g2)
ADD_EXECUTABLE(Poisson ${Poisson_SRCS})
TARGET_LINK_LIBRARIES(Poisson IFEM ${IFEM_DEPLIBS})

# Installation
INSTALL(TARGETS Poisson DESTINATION bin COMPONENT bin)
INSTALL(FILES ${Poisson_SRCS} DESTINATION share/doc/IFEM/Apps/Poisson COMPONENT examples)
INSTALL(FILES ${Poisson_HDRS} DESTINATION share/doc/IFEM/Apps/Poisson COMPONENT examples)
INSTALL(FILES ${Poisson_TSTS} DESTINATION share/doc/IFEM/Apps/Poisson/Test COMPONENT examples)
