# - Basic CMake setup
# Check version meets our requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(MLResurrectionModule)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
add_library(MLResurrectionModule SHARED MLResurrection.h MLResurrection.cpp lib/FLSimulateCommandLine.h lib/FLSimulateErrors.h lib/FLSimulateArgs.h lib/FLSimulateArgs.cc)

# Build an executable from our sources
#add_executable(MLResurrectionModule MLResurrectionModule.cpp MLResurrectionModule.h  FLSimulateArgs.h  FLSimulateCommandLine.h  FLSimulateErrors.h FLSimulateArgs.cc)

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
# - Explicit link to ROOT Physics library which is not linked
#   by flreconstruct
target_link_libraries(MLResurrectionModule
  PUBLIC
    Falaise::FalaiseModule
    Bayeux::Bayeux_mctools_geant4		
    ${ROOT_Physics_LIBRARY}
    )

# Configure example pipeline script for use from the build dir
configure_file("MLResurrectionModule.conf.in" "MLResurrectionModule.conf" @ONLY)

file(COPY official-2.0.0_noCAT.conf DESTINATION ${PROJECT_BINARY_DIR})
#file(COPY official-2.0.0_from_CAT.conf DESTINATION ${PROJECT_BINARY_DIR})
file(COPY list_shorts.txt DESTINATION ${PROJECT_BINARY_DIR})

# Add a basic test: simulated a brio, reconstruct up to calibration and runs the module on it
enable_testing()

# - Simulate
add_test(NAME testMLResurrectionModule_simulate
  COMMAND Falaise::flsimulate -o test-sim.brio
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

# - Reconstruct up to calibration
add_test(NAME testMLResurrectionModule_reconstruct
	COMMAND Falaise::flreconstruct -i test-sim.brio -p official-2.0.0_noCAT.conf -o test-cal.brio
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

# - Run Module
add_test(NAME testMLResurrectionModule_module
  COMMAND Falaise::flreconstruct -i test-cal.brio -p MLResurrectionModule.conf
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
