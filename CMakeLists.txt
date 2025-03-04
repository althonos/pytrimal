cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# Platform information
message(STATUS "Building for processor: ${CMAKE_SYSTEM_PROCESSOR}")

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindAVX2.cmake")
include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")

# --- Setup include folders ----------------------------------------------------

include("src/scripts/cmake/CythonExtension.cmake")
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Compile code -------------------------------------------------------------

add_subdirectory("src")
