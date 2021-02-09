#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
message(STATUS
            "CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

# Check if we have a fortran compiler
find_program(FORT NAMES $ENV{FC}
                        gfortran
                        ifort)
if(NOT FORT)
  project(FuelPinSF C CXX)
  message(STATUS
          "No FORTRAN compiler")
else()
  project(FuelPinSF C CXX Fortran)
  message(STATUS
          "Found FORTRAN compiler")
endif()

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

if(NOT NO_MPI_BUILD)
  find_package(G4mpi REQUIRED)
  add_definitions(-DTOOLS_USE_NATIVE_MPI) 
  add_definitions(-DIS_MPI_COMPILED)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/SaG4n/include
                    ${Geant4_INCLUDE_DIR}
                    ${G4mpi_INCLUDE_DIR}
                    ${PROJECT_BINARY_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/SaG4n/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                   ${PROJECT_SOURCE_DIR}/SaG4n/include/*.hh
		  ${PROJECT_BINARY_DIR}/Config.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
if(NOT NO_MPI_BUILD)
  add_executable(FuelPinSF-mpi FuelPinSF-mpi.cc ${sources} ${headers})
  target_link_libraries(FuelPinSF-mpi ${Geant4_LIBRARIES} ${G4mpi_LIBRARIES})
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build rdecay01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(rdecay01_SCRIPTS
    macros/Cf252.mac macros/U238.mac macros/Cm244.mac macros/singleDecay.mac macros/timeWindow.mac
    macros/timeWindowBiased.mac macros/vis.mac
  )

foreach(_script ${rdecay01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#set environment variable to point to the directory of input files
#set(ENV{FUEL_SIM_INPUT_DIR} ${PROJECT_BINARY_DIR}/input)
#set($ENV{FUEL_SIM_INPUT_DIR} input)
set(FUEL_SIM_INPUT_DIR ${PROJECT_SOURCE_DIR}/input)
#add_definitions( -DFUEL_SIM_INPUT_DIR=${PROJECT_BINARY_DIR}/input )
configure_file(Config.h.in Config.h)
#target_include_directories(FuelPinSF PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           )

#----------------------------------------------------------------------------
# Copy all input data to the build directory.
# This is so that we can run the executable directly because it
# relies on these data being in the current working directory.
#
set(input_files
    input/actinides.txt
    input/activities.txt
    input/activity-nat-Pu.txt
    input/atomic-mass-eval-2016.txt
    input/DummyActivity.txt
    input/elements.txt
    input/MetaStables.txt
    input/SKB-TR-10-13-PWR-5years-nuclide-vector.plt
    input/SKB-TR-10-13-PWR-5years-nuclide-vector-truncated.plt
    input/SKB-TR-10-13-PWR_gamma_spectrum.txt
    input/SKB-TR-10-13-PWR_MCNP_matls.inp
    input/SKB-TR-10-13-PWR_MCNP_matls_new.inp
    input/SKB-TR-10-13-PWR_neutron_spectrum.txt
    input/test_mat.plt
  )

foreach(_script ${input_files})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
if(NOT NO_MPI_BUILD)
install(TARGETS FuelPinSF-mpi DESTINATION bin)

# Check if we override built-in version with new version (libFission.*) of the library
if(USE_FISSION_NEW)
  find_library(LFISSION Fission
               PATHS ${CMAKE_SOURCE_DIR}/../lib ${FISSION_DIR}/lib
               DOC "fission library"
               NO_DEFAULT_PATH
              )
  if(NOT LFISSION)
    message(FATAL_ERROR
            "fission library not found")
  else()
    target_link_libraries(FuelPinSF-mpi 
                          ${LFISSION} 
                         )
    configure_file(${LFISSION} .
                  COPYONLY
                  )
    message(STATUS ${LFISSION})
  endif()
  message(STATUS
          "Overriding fission library libFission.* with version in ${FISSION_DIR}")
  include_directories(${CMAKE_SOURCE_DIR}/../include ${FISSION_DIR}/include)
  add_definitions(-DFISSION_NEW)

  # Check if we use FREYA 
  if(USE_FREYA)
    if(NOT FORT)
      message(FATAL_ERROR
              "No Fortran compiler found, so FREYA is not included")
    endif()
    message(STATUS
            "Using FREYA")
    add_definitions(-DUSEFREYA)
  endif()
endif()

if(DEBUG_ALPHA_N)
message(STATUS
            "DEBUGING : Compiling in (alpha,n) mode only !!!!!")
  add_definitions(-DDEBUG_ALPHA_N)
ENDIF()

endif()
#------------------------------------------------------------------------------------------

add_executable(FuelPinSF FuelPinSF.cc ${sources} ${headers})
#target_link_libraries(FuelPinSF ${Geant4_LIBRARIES} ${G4mpi_LIBRARIES})
target_link_libraries(FuelPinSF ${Geant4_LIBRARIES})

target_compile_definitions(FuelPinSF PUBLIC -DNOT_USING_MPI)



# Check if we override built-in version with new version (libFission.*) of the library
if(USE_FISSION_NEW)
  find_library(LFISSION Fission
               PATHS ${CMAKE_SOURCE_DIR}/../lib ${FISSION_DIR}/lib
               DOC "fission library"
               NO_DEFAULT_PATH
              )
  if(NOT LFISSION)
    message(FATAL_ERROR
            "fission library not found")
  else()
    target_link_libraries(FuelPinSF 
                          ${LFISSION} 
                         )
    configure_file(${LFISSION} .
                  COPYONLY
                  )
    message(STATUS ${LFISSION})
  endif()
  message(STATUS
          "Overriding fission library libFission.* with version in ${FISSION_DIR}")
  include_directories(${CMAKE_SOURCE_DIR}/../include ${FISSION_DIR}/include)
  add_definitions(-DFISSION_NEW)

  # Check if we use FREYA 
  if(USE_FREYA)
    if(NOT FORT)
      message(FATAL_ERROR
              "No Fortran compiler found, so FREYA is not included")
    endif()
    message(STATUS
            "Using FREYA")
    add_definitions(-DUSEFREYA)
  endif()
endif()

if(DEBUG_ALPHA_N)
message(STATUS
            "DEBUGING : Compiling in (alpha,n) mode only !!!!!")
  add_definitions(-DDEBUG_ALPHA_N)
ENDIF()