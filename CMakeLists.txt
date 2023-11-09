##################################################################
# CMAPLE cmake build definition
# Copyright (c) 2022 Nhan Ly-Trong, Chris Bielow, Nicola De Maio, Bui Quang Minh 
##################################################################

# Windows example usages:
#------------------------
# cmake -G "Visual Studio 12" <source_dir>          (32-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" <source_dir>    (64-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" <source_dir>        (64-bit version, compiled with ICC)
# cmake -G "MinGW Makefiles" <source_dir> (TDM-GCC)
# cmake -G "Unix Makefiles" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_MAKE_PROGRAM=mingw32-make <source_dir> (TDM-GCC and clang)

# Linux example usages:
#----------------------
# cmake <source_dir>                           
# cmake  -DCMAKE_CXX_FLAGS="-g -fno-omit-frame-pointer -march=native" <source_dir>   ## some custom compile flags
#
# To compile with CLANG on Linux:CMAPLE_
# export CC=/usr/bin/clang
# export CXX=/usr/bin/clang++
# Best practices for setting up CMAKE for diffrent compiler can be found here:
# http://stackoverflow.com/questions/7031126/switching-between-gcc-and-clang-llvm-using-cmake
#
# Mac OSX example usages:
#------------------------
#
# Assuming clang3.7 and clang++3.7 are the newly built compilers, then:
# cmake -DCMAKE_C_COMPILER=clang3.7 -DCMAKE_CXX_COMPILER=clang++3.7 <source_dir>   (OpenMP version)
#
# Xcode project example usages:
#------------------------------
#
# To generate Xcode project without OpenMP:
# cmake -G Xcode <source_dir>
#
#

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

set(GCC_MIN_VERSION "5")


##############################################################
# Clang-tidy
##############################################################
message("Setup clang-tidy")
#set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-header-filter=.;-checks=performance-*")

project(cmaple)
add_definitions(-DCMAPLE)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE) # Enable IPO (LTO) by default

# The version number.
set (cmaple_VERSION_MAJOR 1)
set (cmaple_VERSION_MINOR 0)
set (cmaple_VERSION_PATCH ".0")

option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("Build mode   : Release")
endif()

if (CMAKE_GENERATOR MATCHES "Xcode")
    set(CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
endif()

include_directories("${PROJECT_SOURCE_DIR}")
include_directories(after system "/usr/local/include")

##################################################################
# Setup compiler, currently supported GCC, CLANG, MSVC, and ICC
##################################################################


##################################################################
# Setup compiler flags
##################################################################

## enable 'SSE/AVX' on x86-64, 'neon' on arm to achive faster computations (mainly the Matrix::dotProduct())
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86") 
    add_compile_options(-msse -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2 -mavx) # needed for simde instructions
  elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    add_compile_options(-neon)  
  endif() 
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86") 
    add_compile_options(/arch:AVX2) ## not really needed, but may be beneficial in other places when in Release mode (i.e. with optimizations turned on)
  endif()
endif()


##################################################################
# check existence of a few basic functions
##################################################################
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (getrusage HAVE_GETRUSAGE)
check_function_exists (GlobalMemoryStatusEx HAVE_GLOBALMEMORYSTATUSEX)
check_function_exists (strndup HAVE_STRNDUP)
check_function_exists (strtok_r HAVE_STRTOK_R)

find_package(Backtrace)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/cmaple_config.h.in"
  "${PROJECT_BINARY_DIR}/cmaple_config.h"
  )

# add the binary tree to the search path for include files
# so that we will find cmaple_config.h
include_directories("${PROJECT_BINARY_DIR}")

#zlib will be detected for appearance
#include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7")


if (NOT IQTREE_FLAGS MATCHES "nozlib")
    find_package(ZLIB)
endif()

if(ZLIB_FOUND)
  message ("Using system zlib")
  include_directories(${ZLIB_INCLUDE_DIRS})
else(ZLIB_FOUND)
  message ("Using own zlib-1.2.7")
  include_directories("${PROJECT_BINARY_DIR}/libraries/zlib-1.2.7" "${PROJECT_SOURCE_DIR}/libraries/zlib-1.2.7")
  if (NOT TARGET zlibstatic)
    add_subdirectory(libraries/zlib-1.2.7)
  endif()
endif(ZLIB_FOUND)

##################################################################
# fetch GoogleTest framework when CMake runs
##################################################################
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

##################################################################
# subdirectories containing necessary libraries for the build
##################################################################
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
if (NOT TARGET ncl)
	add_subdirectory(libraries/ncl)
endif() 
if (NOT TARGET nclextra)
	add_subdirectory(libraries/nclextra)
endif() 
add_subdirectory(utils)
add_subdirectory(libraries/SIMDe) 
include_directories(${SIMDE_INCLUDE_DIRECTORY})
add_subdirectory(model)
add_subdirectory(alignment)
add_subdirectory(tree)
add_subdirectory(maple)
add_subdirectory(unittest)


##################################################################
# the main executable
##################################################################
add_executable(cmaple main/main.cpp)
add_executable(cmaple-aa main/main.cpp)

if(Backtrace_FOUND)
  include_directories(${Backtrace_INCLUDE_DIR})
  target_link_libraries(cmaple ${Backtrace_LIBRARY})
  target_link_libraries(cmaple-aa ${Backtrace_LIBRARY})
endif(Backtrace_FOUND)

##################################################################
# setup linking flags
##################################################################

target_link_libraries(cmaple maple cmaple_tree cmaple_alignment cmaple_model ncl nclextra cmaple_utils)
target_link_libraries(cmaple-aa maple-aa cmaple_tree-aa cmaple_alignment-aa cmaple_model-aa ncl nclextra cmaple_utils)

##############################################################
# add the install targets
##############################################################
install (TARGETS cmaple DESTINATION bin)
install (TARGETS cmaple-aa DESTINATION bin)
#install (FILES "${PROJECT_SOURCE_DIR}/example/input.fa" DESTINATION .)
#install (FILES "${PROJECT_SOURCE_DIR}/example/input.fa.ref" DESTINATION .)
#install (FILES "${PROJECT_SOURCE_DIR}/example/output.maple" DESTINATION .)
#install (FILES "${PROJECT_BINARY_DIR}/cmaple_config.h" DESTINATION "include")
#install(DIRECTORY "${PROJECT_SOURCE_DIR}/" # source directory
#        DESTINATION "include/cmaple" # target directory
#        FILES_MATCHING # install only matched files
#        PATTERN "*.h" # select header files
#)

if (WIN32)
    install (FILES "${BINARY_DIR}/cmaple${EXE_SUFFIX}-click.exe" DESTINATION bin)
    install (FILES "${BINARY_DIR}/cmaple-aa${EXE_SUFFIX}-click.exe" DESTINATION bin)
endif()

##############################################################
# build a CPack driven installer package
##############################################################
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${cmaple_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${cmaple_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${cmaple_VERSION_PATCH}")
if(WIN32 OR APPLE)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

#set(CPACK_SOURCE_PACKAGE_FILE_NAME
#  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build.*/;/debug.*/;/examples/;/test_scripts/;/manual/;/.bzr/;~$;/\\\\.svn/;/\\\\.git/;/pllrepo/;${CPACK_SOURCE_IGNORE_FILES}")


set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}${CPACK_PACKAGE_VERSION_PATCH}-${SYSTEM_NAME}")

if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()

include (CPack)

#add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
