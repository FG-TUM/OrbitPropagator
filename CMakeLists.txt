cmake_minimum_required(VERSION 3.10)

project(OrbitPropagator)
set(CMAKE_CXX_STANDARD 17)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

enable_testing()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (PROJECT_BINARY debris_sim)
add_subdirectory(src)

option(DebrisSim_Tests "Build tests" ON)
if(DebrisSim_Tests)
    add_subdirectory(lib/googletest)
    add_subdirectory(test)
endif()

option(DebrisSim_Thesis "Build the thesis document" OFF)
if(DebrisSim_Thesis)
    add_subdirectory(Thesis)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(doxygen)
