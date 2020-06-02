# - Basic CMake setup
# Check version meets our requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
<<<<<<< HEAD
project(LazarusModule)
=======
project(MLResurrectionModule)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
>>>>>>> Working version of CNN module for dead cells

find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
<<<<<<< HEAD
add_library(LazarusModule SHARED LazarusModule.h LazarusModule.cpp lib/FLSimulateCommandLine.h lib/FLSimulateErrors.h lib/FLSimulateArgs.h lib/FLSimulateArgs.cc)

# Build an executable from our sources
#add_executable(LazarusModule LazarusModule.cpp LazarusModule.h  FLSimulateArgs.h  FLSimulateCommandLine.h  FLSimulateErrors.h FLSimulateArgs.cc)
=======
add_library(MLResurrectionModule SHARED MLResurrection.h MLResurrection.cpp lib/FLSimulateCommandLine.h lib/FLSimulateErrors.h lib/FLSimulateArgs.h lib/FLSimulateArgs.cc)

# Build an executable from our sources
#add_executable(MLResurrectionModule MLResurrectionModule.cpp MLResurrectionModule.h  FLSimulateArgs.h  FLSimulateCommandLine.h  FLSimulateErrors.h FLSimulateArgs.cc)
>>>>>>> Working version of CNN module for dead cells

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
# - Explicit link to ROOT Physics library which is not linked
#   by flreconstruct
<<<<<<< HEAD
target_link_libraries(LazarusModule
=======
target_link_libraries(MLResurrectionModule
>>>>>>> Working version of CNN module for dead cells
  PUBLIC
    Falaise::FalaiseModule
    Bayeux::Bayeux_mctools_geant4		
    ${ROOT_Physics_LIBRARY}
    )

# Configure example pipeline script for use from the build dir
<<<<<<< HEAD
configure_file("LazarusModule.conf.in" "LazarusModule.conf" @ONLY)

file(COPY official-2.0.0_only_calibration.conf DESTINATION ${PROJECT_BINARY_DIR})
file(COPY official-2.0.0_from_CAT.conf DESTINATION ${PROJECT_BINARY_DIR})
file(COPY random_dead_cells_100.txt DESTINATION ${PROJECT_BINARY_DIR})
=======
configure_file("MLResurrectionModule.conf.in" "MLResurrectionModule.conf" @ONLY)

file(COPY official-2.0.0_noCAT.conf DESTINATION ${PROJECT_BINARY_DIR})
#file(COPY official-2.0.0_from_CAT.conf DESTINATION ${PROJECT_BINARY_DIR})
file(COPY list_shorts.txt DESTINATION ${PROJECT_BINARY_DIR})
>>>>>>> Working version of CNN module for dead cells

# Add a basic test: simulated a brio, reconstruct up to calibration and runs the module on it
enable_testing()

# - Simulate
<<<<<<< HEAD
add_test(NAME testLazarusModule_simulate
=======
add_test(NAME testMLResurrectionModule_simulate
>>>>>>> Working version of CNN module for dead cells
  COMMAND Falaise::flsimulate -o test-sim.brio
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

# - Reconstruct up to calibration
<<<<<<< HEAD
add_test(NAME testLazarusModule_reconstruct
  COMMAND Falaise::flreconstruct -i test-sim.brio -p official-2.0.0_only_calibration.conf -o test-cal.brio
=======
add_test(NAME testMLResurrectionModule_reconstruct
	COMMAND Falaise::flreconstruct -i test-sim.brio -p official-2.0.0_noCAT.conf -o test-cal.brio
>>>>>>> Working version of CNN module for dead cells
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

# - Run Module
<<<<<<< HEAD
add_test(NAME testLazarusModule_module
  COMMAND Falaise::flreconstruct -i test-cal.brio -p LazarusModule.conf
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

=======
add_test(NAME testMLResurrectionModule_module
  COMMAND Falaise::flreconstruct -i test-cal.brio -p MLResurrectionModule.conf
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
>>>>>>> Working version of CNN module for dead cells
