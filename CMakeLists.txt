cmake_minimum_required(VERSION 3.1)
project(NekRS)
enable_language(C CXX Fortran)
include(FortranCInterface)

set(CMAKE_VERBOSE_MAKEFILE on)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/occa/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/ins")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/elliptic")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/parAlmond")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/libparanumal")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/gatherScatter")

###############################################################################
# Check compiler/MPI vendor, version and Fortran/C compatibility              #
###############################################################################
macro(check_compiler_id compiler)
  if(NOT "${compiler}" STREQUAL "GNU")
    message(FATAL_ERROR "NekRS only supports GNU compilers")
  endif()
endmacro()

macro(check_compiler_version compiler_version)
  if(${compiler_version} LESS 6.2)
    message(FATAL_ERROR "NekRS only supports GNU compiler versions >= 6.2")
  endif()
endmacro()

check_compiler_id("${CMAKE_CXX_COMPILER_ID}")
check_compiler_id("${CMAKE_C_COMPILER_ID}")
check_compiler_id("${CMAKE_Fortran_COMPILER_ID}")

check_compiler_version("${CMAKE_CXX_COMPILER_VERSION}")
check_compiler_version("${CMAKE_C_COMPILER_VERSION}")
check_compiler_version("${CMAKE_Fortran_COMPILER_VERSION}")

## Sanity checks for MPI
set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
find_package(MPI REQUIRED)
if(NOT MPI_FOUND)
  message(FATAL_ERROR "MPI is not found for one of the languages: C/CXX/Fortran")
endif()
#if (MPI_VERSION VERSION_LESS "3.1")
#  message(FATAL_ERROR "NekRS requires an MPI version >= 3.1")
#endif()

## Check Fortran/C interface
FortranCInterface_VERIFY(CXX QUIET)
if (NOT FortranCInterface_VERIFIED_CXX)
    message(WARNING "C and Fortran may not work properly")
endif()

## Use c++14
set (CMAKE_CXX_STANDARD 14)

###############################################################################
# Setup cmake variables                                                       #
###############################################################################
set(BIN_THIRD_PARTY_DIR "${CMAKE_BINARY_DIR}/3rd_party" CACHE STRING "")
set(SRC_THIRD_PARTY_DIR "${CMAKE_SOURCE_DIR}/3rd_party" CACHE STRING "")

set(LIBPDIR "${BIN_THIRD_PARTY_DIR}/libparanumal/")
set(ELLIPTICDIR "${LIBPDIR}/solvers/elliptic/")
set(INSDIR "${LIBPDIR}/solvers/ins/")
set(PARALMONDDIR "${LIBPDIR}/libs/parAlmond/")
set(OGSDIR "${LIBPDIR}/libs/gatherScatter/")
set(NEKDIR "${BIN_THIRD_PARTY_DIR}/nek5000/")
set(OCCADIR "${BIN_THIRD_PARTY_DIR}/occa/")

set(NEKINTERFACEDIR "${CMAKE_SOURCE_DIR}/src/nekInterface/")

###############################################################################
# Check if installed already                                                  #
###############################################################################
if (EXISTS "${CMAKE_INSTALL_PREIFX}/occa")
  message(WARNING "Found OCCA already installed at ${CMAKE_INSTALL_PREIFX}/.")
endif()
if (EXISTS "${CMAKE_INSTALL_PREIFX}/libparanumal")
  message(WARNING "Found libparanumal already installed at ${CMAKE_INSTALL_PREIFX}/.")
endif()

###############################################################################
# Set compiler flags                                                          #
###############################################################################
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

if (NOT $ENV{NEKRS_LIBP_DEFINES} STREQUAL "")
  set(LIBP_DEFINES "$ENV{NEKRS_LIBP_DEFINES}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LIBP_DEFINES}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBP_DEFINES}")
endif()

###############################################################################
# Start building nekRS                                                        #
###############################################################################
file(MAKE_DIRECTORY ${BIN_THIRD_PARTY_DIR})

set(NEKRS_INCLUDES)
set(THIRD_PARTY_INCLUDES)

macro(set_static_lib target_name builddir)
   set(${target_name}_lib
     "${builddir}/${CMAKE_STATIC_LIBRARY_PREFIX}${target_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endmacro()

macro(set_shared_lib target_name builddir)
   set(${target_name}_lib
     "${builddir}/${CMAKE_SHARED_LIBRARY_PREFIX}${target_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
endmacro()

###############################################################################
# Build OCCA                                                                  #
###############################################################################

set_shared_lib(occa "${OCCADIR}/lib")
file(REMOVE ${occa_lib})
add_custom_target(occa_build DEPENDS ${occa_lib})
add_custom_command(OUTPUT ${occa_lib}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_CXX=${CMAKE_CXX_COMPILER}
    NEKRS_CXXFLAGS=${CMAKE_CXX_FLAGS}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
    ${SRC_THIRD_PARTY_DIR}/occa.build
)
add_library(occa SHARED IMPORTED)
set_target_properties(occa PROPERTIES IMPORTED_LOCATION ${occa_lib})
set(THIRD_PARTY_INCLUDES ${OCCADIR}/include ${OCCADIR})

###############################################################################
# Build Nek5000                                                               #
###############################################################################
set(NEK_BUILD_LIBRARIES)
set(NEK_INCLUDES)
macro(FIND_NEK_DEPS target_name builddir include_)
  set_static_lib(${target_name} ${builddir})
  set(NEK_BUILD_LIBRARIES ${NEK_BUILD_LIBRARIES} ${${target_name}_lib})
  set(NEK_INCLUDES ${NEK_INCLUDES} ${include_})
endmacro()

FIND_NEK_DEPS(gs ${NEKDIR}/3rd_party/gslib/gslib/lib
  ${NEKDIR}/3rd_party/gslib/gslib/include)
FIND_NEK_DEPS(parRSB ${NEKDIR}/3rd_party/parRSB/parRSB/lib
  ${NEKDIR}/3rd_party/parRSB/parRSB/include)
FIND_NEK_DEPS(blasLapack ${NEKDIR}/3rd_party/blasLapack
  ${NEKDIR}/3rd_party/blasLapack)

file(REMOVE ${NEK_BUILD_LIBRARIES})
add_custom_target(nek_build DEPENDS ${NEK_BUILD_LIBRARIES})
add_custom_command(OUTPUT ${NEK_BUILD_LIBRARIES}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_FC=${CMAKE_Fortran_COMPILER}
    NEKRS_FFLAGS=${CMAKE_Fortran_FLAGS}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_NEK5000_PPLIST=${NEK5000_PPLIST}
    ${SRC_THIRD_PARTY_DIR}/nek5000.build
)

###############################################################################
# Set include dirs                                                            #
###############################################################################
set(THIRD_PARTY_INCLUDES ${THIRD_PARTY_INCLUDES}
  ${LIBPDIR}/include
  ${PARALMONDDIR}
  ${ELLIPTICDIR}
  ${OGSDIR}
  ${INSDIR}
)

set(NEKRS_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/udf
  ${NEKINTERFACEDIR}
)

###############################################################################
# Build libP                                                                  #
###############################################################################
macro(FIND_DEPS target_name builddir)
  set_shared_lib(${target_name} ${builddir})
  add_library(${target_name} SHARED IMPORTED)
  set_target_properties(${target_name} PROPERTIES IMPORTED_LOCATION
    ${${target_name}_lib})
endmacro()

set(LIBP_BUILD_LIBRARIES)
macro(FIND_LIBP_DEPS target_name builddir)
  FIND_DEPS(${target_name} ${builddir})
  set(LIBP_BUILD_LIBRARIES ${LIBP_BUILD_LIBRARIES} ${${target_name}_lib})
endmacro()

FIND_LIBP_DEPS(P "${LIBPDIR}")
FIND_LIBP_DEPS(ins "${INSDIR}")
FIND_LIBP_DEPS(elliptic "${ELLIPTICDIR}")
FIND_LIBP_DEPS(parAlmond "${PARALMONDDIR}")
FIND_LIBP_DEPS(ogs "${OGSDIR}")

file(REMOVE ${LIBP_BUILD_LIBRARIES})
add_custom_target(libparanumal_build DEPENDS ${LIBP_BUILD_LIBRARIES})
add_custom_command(OUTPUT ${LIBP_BUILD_LIBRARIES}
  WORKING_DIRECTORY ${BIN_THIRD_PARTY_DIR}
  COMMAND NEKRS_CC=${CMAKE_C_COMPILER}
    NEKRS_CFLAGS=${CMAKE_C_FLAGS}
    NEKRS_CXX=${CMAKE_CXX_COMPILER}
    NEKRS_CXXFLAGS=${CMAKE_CXX_FLAGS}
    NEKRS_FC=${CMAKE_Fortran_COMPILER}
    NEKRS_FFLAGS=${CMAKE_Fortran_FLAGS}
    NEKRS_LIBP_DIR=${CMAKE_INSTALL_PREFIX}
    NEKRS_LIBP_DEFINES=${LIBP_DEFINES}
    NEKRS_SRC_THIRD_PARTY_DIR=${SRC_THIRD_PARTY_DIR}
    NEKRS_BIN_THIRD_PARTY_DIR=${BIN_THIRD_PARTY_DIR}
    NEKRS_OCCA_DIR=${BIN_THIRD_PARTY_DIR}/occa
    NEKRS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
    ${SRC_THIRD_PARTY_DIR}/libparanumal/libparanumal.build
  DEPENDS occa_build
)

###############################################################################
# Add nekrs target                                                            #
###############################################################################
set(SRC src/nekrs.cpp
    src/timestepper.cpp
    src/udf/udf.cpp
    src/meshPhysicalNodesHex3D.cpp
    src/occaDeviceConfig.cpp
    src/report.cpp
    src/nekInterface/nekInterfaceAdapter.cpp
    src/parReader.cpp
    ## To get fortran flags
    src/dummy.f
)

add_executable(nekrs ${SRC})
set_target_properties(nekrs PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(nekrs ${CMAKE_DL_LIBS} ins elliptic parAlmond P ogs occa)
target_include_directories(nekrs PRIVATE ${THIRD_PARTY_INCLUDES} ${NEKRS_INCLUDES})
add_dependencies(nekrs occa_build libparanumal_build nek_build)

#################################################################################
### Install                                                                     #
#################################################################################
# Install dependencies and examples
configure_file(install_script.in.cmake install_script.cmake @ONLY)
install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/install_script.cmake)

#install nekRS
install(TARGETS nekrs
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  PRIVATE_HEADER DESTINATION include
)
