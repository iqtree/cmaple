##################################################################
# CMAPLE cmake build definition
# Copyright (c) 2022 Bui Quang Minh, Nhan Ly-Trong
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
# cmake -DCMAPLE_FLAGS="m32" <source_dir>  (32-bit version)
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
# Compile 32-bit version: cmake -DCMAPLE_FLAGS=m32 ....
# Compile static version: cmake -DCMAPLE_FLAGS=static ....
#

#NOTE: Static linking with clang windows: make a symlink libgcc_eh.a to libgcc.a (administrator required)
# C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0>mklink libgcc_eh.a libgcc.a
#
#

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

set(GCC_MIN_VERSION "5")

project(cmaple)
add_definitions(-DCMAPLE)

# The version number.
set (cmaple_VERSION_MAJOR 1)
set (cmaple_VERSION_MINOR 0)
set (cmaple_VERSION_PATCH ".0")

option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)

message("CMAPLE flags : ${CMAPLE_FLAGS}")

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

set (__ARM_NEON "FALSE")
set (NEON 0)
if (APPLE)
	execute_process(
      COMMAND uname -m
      RESULT_VARIABLE result
      OUTPUT_VARIABLE OSX_NATIVE_ARCHITECTURE
      OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	if(OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
		message("APPLE SILICON (M1) DETECTED!")
		set (__ARM_NEON "TRUE")
		set (NEON 1)
	endif()
endif()

##################################################################
# Detect target platforms
##################################################################
if (WIN32)
    message("Target OS     : Windows")
    # build as static binary to run on most machines
    if (CMAPLE_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (APPLE)
    message("Target OS     : Mac OS X")
	if(OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
		add_definitions("--target=arm64-apple-macos10.5")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --target=arm64-apple-macos12.0.1")
	else()
		# to be compatible back to Mac OS X 10.7
    	if (CMAPLE_FLAGS MATCHES "oldmac")
        	add_definitions("-mmacosx-version-min=10.5")
        	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.5")
    	else()
        	add_definitions("--target=x86_64-apple-macos10.7")
        	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --target=x86_64-apple-macos10.7")
    	endif()
	endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (UNIX)
    message("Target OS     : Unix")
    # build as static binary to run on most machines
    if (NOT CMAPLE_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
else()
    # Note that CMAPLE has NOT been tested on other platforms
    message("Target OS     : Unknown and untested yet")
endif()

##################################################################
# Setup compiler, currently supported GCC, CLANG, MSVC, and ICC
##################################################################

set(GCC "FALSE")   #  GNU compiler
set(CLANG "FALSE") # Clang compiler
set(ICC "FALSE")   # Intel compiler
set(VCC "FALSE")   # MS Visual C Compiler, note that it is different from MSVC variable
set(CLANG_UNDER_VS "FALSE") #Clang compiler, used from inside Visual Studio
# using C++11 standard
# disable AVX for NEON
if (__ARM_NEON)
    set(IQTREE_FLAGS "${IQTREE_FLAGS} novx")
elseif (CMAKE_CXX_COMPILER MATCHES "VISUAL STUDIO")
    set(CLANG_UNDER_VS "TRUE")
    #it won't recognize the -std=c++11 parameter.
    #Todo: don't hard-code this; figure out some way it can be passed in (though ideally, not the whole shebang).
    include_directories("C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\Llvm\\lib\\clang\\10.0.0\\include")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Compiler      : GNU Compiler (gcc)")
    set(GCC "TRUE")
#    set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -fms-extensions -Wno-deprecated")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -ffunction-sections -fdata-sections -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -g -ffunction-sections -fdata-sections -Wall")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()
    # require at least gcc ${GCC_MIN_VERSION}
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN_VERSION)
        message(FATAL_ERROR "GCC version must be at least ${GCC_MIN_VERSION}!")
    endif()
    if (WIN32)
        # disable AVX on Windows due to memory alignment
        set(CMAPLE_FLAGS "${CMAPLE_FLAGS} novx")
        message("WARNING: AVX is disabled on Windows as GCC does not properly suport memory alignment")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiler      : Clang")
    set(CLANG "TRUE")
#    set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -Wno-nested-anon-types")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections -Wall")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()

    # use libc++ per default in MacOS
    if (APPLE)
        SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()

    #remove -rdynamic for Clang under Linux
    if (UNIX AND CMAPLE_FLAGS MATCHES "static")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(VCC "TRUE")
    message("Compiler      : MS Visual C++ Compiler")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("Compiler      : Intel C++ Compiler (icc)")
    set(ICC "TRUE")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
else()
    message("Compiler      : Unknown and untested yet")
endif()
message("Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

set(EXE_SUFFIX "")

if (MSVC)
    # MS Visual Studio environment
    message("Exporting MS Visual Studio projects...")
    if (CLANG_UNDER_VS)
        #see https://clang.llvm.org/docs/UsersManual.html#clang-cl
        #note .GX is how you say -fexceptions
        add_definitions(/D_UWIN)
        set(CMAKE_C_FLAGS_RELEASE "/O2 /GX")
        set(CMAKE_C_FLAGS_DEBUG   "/D_UWIN /GX")
    else()
        add_definitions(/MP) # enable multi-processor compilation
    endif()
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_definitions(/Ot /Oi)
        if (VCC)
            add_definitions(/O2)
        elseif (ICC)
            add_definitions(/O3)
        endif()
    endif()
endif()

# enable link time optimization
if (CMAPLE_FLAGS MATCHES "lto")
    #if (CLANG)
    #    set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto=thin")
    #else()
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto")
    #endif()
endif()

##################################################################
# detect 32 or 64 bit binary
##################################################################
set (BINARY32 "FALSE")
if(CMAKE_SIZEOF_VOID_P EQUAL 4 OR CMAPLE_FLAGS MATCHES "m32")
    set(BINARY32 "TRUE")
    message("Target binary : 32-bit")
    if (CMAKE_GENERATOR MATCHES "Win64")
        error("Both 32-bit and 64-bit mode cannot be specified")
    endif()
    if (GCC OR CLANG)
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -m32")
    endif()
    add_definitions(-DBINARY32)
    if (WIN32)
        add_definitions(-DWIN32)
    endif()
else()
    message("Target binary : 64-bit")
    if (CLANG_UNDER_VS)
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -m64")
    endif()
    if (WIN32)
        add_definitions(-DWIN64)
    endif()
endif()

##################################################################
# Setup compiler flags
##################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O2 -fno-omit-frame-pointer -g")
set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O2 -fno-omit-frame-pointer -g")

if(CLANG AND CMAPLE_FLAGS MATCHES "static")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -Wl,--allow-multiple-definition")
endif()

if (CMAPLE_FLAGS MATCHES "libcxx")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()


if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Profile")
    message("C flags       : ${CMAKE_C_FLAGS_PROFILE} ")
    message("CXX flags     : ${CMAKE_CXX_FLAGS_PROFILE} ")
endif()

message("LINKER flags  : ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")


if (GCC)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
elseif (CLANG AND NOT CLANG_UNDER_VS)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
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
  include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7" "${PROJECT_SOURCE_DIR}/zlib-1.2.7")
  add_subdirectory(zlib-1.2.7)
endif(ZLIB_FOUND)

##################################################################
# subdirectories containing necessary libraries for the build
##################################################################

add_subdirectory(main)
add_subdirectory(utils)
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
add_subdirectory(model)
add_subdirectory(alignment)
add_subdirectory(tree)

##################################################################
# the main executable
##################################################################
add_executable(cmaple main/main.cpp)

if(Backtrace_FOUND)
  include_directories(${Backtrace_INCLUDE_DIR})
  target_link_libraries(cmaple ${Backtrace_LIBRARY})
endif(Backtrace_FOUND)

##################################################################
# setup linking flags
##################################################################

# link special lib for WIN32
if (WIN32)
    set(PLATFORM_LIB "ws2_32")
else()
    set(PLATFORM_LIB "m")
endif()

if (CMAPLE_FLAGS MATCHES "libcxx")
    set(STD_LIB "c++abi")
endif()

# basic linking librararies
target_link_libraries(cmaple utils model alignment tree main ${PLATFORM_LIB} ${STD_LIB})

# setup the executable name
##################################################################
set_target_properties(cmaple PROPERTIES OUTPUT_NAME "cmaple${EXE_SUFFIX}")

# strip the release build
if (NOT CMAPLE_FLAGS MATCHES "nostrip" AND CMAKE_BUILD_TYPE STREQUAL "Release" AND (GCC OR CLANG) AND NOT APPLE) # strip is not necessary for MSVC
    if (WIN32)
        ADD_CUSTOM_COMMAND(TARGET cmaple POST_BUILD COMMAND strip $<TARGET_FILE:cmaple>)
    elseif (NOT APPLE)
        ADD_CUSTOM_COMMAND(TARGET cmaple POST_BUILD COMMAND strip $<TARGET_FILE:cmaple>)
    endif()
endif()

if (MSVC)
    set (BINARY_DIR "${PROJECT_BINARY_DIR}/Release")
else()
    set (BINARY_DIR "${PROJECT_BINARY_DIR}")
endif()

if (WIN32)
    if (MSVC)
        if (CLANG_UNDER_VS)
            ADD_CUSTOM_COMMAND(TARGET cmaple POST_BUILD COMMAND copy "cmaple${EXE_SUFFIX}.exe" "cmaple${EXE_SUFFIX}-click.exe")
        else()
            ADD_CUSTOM_COMMAND(TARGET cmaple POST_BUILD COMMAND copy "Release\\cmaple${EXE_SUFFIX}.exe" "Release\\cmaple${EXE_SUFFIX}-click.exe")
        endif()
    else()
        ADD_CUSTOM_COMMAND(TARGET cmaple POST_BUILD COMMAND copy "cmaple${EXE_SUFFIX}.exe" "cmaple${EXE_SUFFIX}-click.exe")
    endif()
endif()

##############################################################
# add the install targets
##############################################################
install (TARGETS cmaple DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/example/input.fa" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/example/input.fa.ref" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/example/output.diff" DESTINATION .)

if (WIN32)
    install (FILES "${BINARY_DIR}/cmaple${EXE_SUFFIX}-click.exe" DESTINATION bin)
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

set (SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    if (CMAPLE_FLAGS MATCHES "oldmac")
        set (SYSTEM_NAME "MacOS10.5")
    else()
        set (SYSTEM_NAME "MacOSX")
    endif()
endif()

if (BINARY32)
    set (SYSTEM_NAME "${SYSTEM_NAME}32")
endif()

if (CMAPLE_FLAGS MATCHES "KNL")
    set (SYSTEM_NAME "${SYSTEM_NAME}KNL")
endif()

set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}${CPACK_PACKAGE_VERSION_PATCH}-${SYSTEM_NAME}")

if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()

include (CPack)

#add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
