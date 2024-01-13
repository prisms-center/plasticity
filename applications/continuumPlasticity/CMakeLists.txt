##
#  CMake script for the PRISMS-Plasticity continuum applications:
##

# Set the name of the project and target:
SET(TARGET "main")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc
  # You can specify additional files here!
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# Find deal.II installation
FIND_PACKAGE(deal.II 8.3.0 REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

# Check to make sure deal.II is configured with p4est
IF(NOT ${DEAL_II_WITH_P4EST})
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II was not installed with p4est. ***\n\n"
    “The p4est library is a mandatory prerequisite for PRISMS-PF. Please consult the \n”
    “user guide to confirm that deal.II and p4est were installed and configured correctly.”
    )
ENDIF()

#
# Are all dependencies fullfilled?
#
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_P4EST = ON
One or all of these are OFF in your installation but are required for this tutorial step."
    )
ENDIF()
DEAL_II_INITIALIZE_CACHED_VARIABLES()

INCLUDE_DIRECTORIES(${DEAL_II_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
INCLUDE_DIRECTORIES(/../../include)
SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING
	"Choose the type of build, options are: Debug, Release"
	)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-deprecated-declarations -Wno-reorder -Wno-unused-variable -Wno-extra")

PROJECT(${TARGET})

ADD_LIBRARY(PLibrary SHARED ${CMAKE_CURRENT_SOURCE_DIR}/../../src/materialModels/continuumPlasticity/models/PLibrary.cc)
DEAL_II_SETUP_TARGET(PLibrary)

##ADD_LIBRARY(BVPLibrary STATIC ${EBVPfunctions})
##DEAL_II_SETUP_TARGET(BVPLibrary)


ADD_EXECUTABLE(${TARGET} ${TARGET}.cc)
DEAL_II_SETUP_TARGET(${TARGET})

TARGET_LINK_LIBRARIES(${TARGET} ${DEAL_II_LIBRARIES} PLibrary)
TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../libprisms_cp_debug.a)

ADD_CUSTOM_TARGET(run
  COMMAND ${TARGET} -pc_type jacobi
  COMMENT "Run ${TARGET} with ${CMAKE_BUILD_TYPE} configuration"
  )

#Debug or release
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
