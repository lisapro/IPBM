cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Create main project.
project(brom2 Fortran)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  endif()
endif()

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

find_package(NetCDF REQUIRED)
include_directories("${NetCDF_INCLUDE_DIRS}")
if (NetCDF_STATIC_MSVC_BUILD)
  # On Windows with a statically-compiled NetCDF library - compile all code against static runtime.
  # This MUST be done before any targets are added.
  add_compile_options("/libs:static")
endif()

find_path(FABM_BASE src/fabm.F90 DOC "Path to FABM source directory.")
mark_as_advanced(CLEAR FABM_BASE)
set(FABM_FORCED_HOST brom)
add_subdirectory(${FABM_BASE}/src fabm)
set(FABM_LIBRARIES fabm)
set(FABM_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fabm/modules)

add_executable(brom2 WIN32
               transport.F90
               diff_mod.F90
               variables_mod.F90
               input_mod.F90
               output_mod.F90
               types_mod.F90
               list_mod.F90
               item_mod.F90
               ../include/brom.h
               ../include/parameters.h)

set_property(TARGET brom2 APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
set_property(TARGET brom2 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
set_property(TARGET brom2 PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set_property(TARGET brom2 APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")

target_link_libraries(brom2 ${FABM_LIBRARIES})
target_link_libraries(brom2 "${NetCDF_LIBRARIES}")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    target_compile_options(brom2 PRIVATE /nothreads /libs:static /libs:qwin)
  endif()
endif()

if (NetCDF_STATIC_MSVC_BUILD)
  set_property(TARGET brom2 PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmt\"")
  set_property(TARGET brom2 PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmtd\"")
endif()
