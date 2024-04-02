# cmake minimum version
cmake_minimum_required(VERSION 3.0)
# Project name
project(TestParser)
# C++ 11 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif ()
message(STATUS "Build type: '${CMAKE_BUILD_TYPE}'")


# find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

# include library from /home/ming/mcnptools
include_directories(/home/ming/mcnptools/include)
link_directories(/home/ming/mcnptools/lib)

# include(FetchContent)
# FetchContent_Declare(
#   mcnptools
#   GIT_REPOSITORY https://github.com/fm140905/mcnptools.git
# )
# FetchContent_MakeAvailable(mcnptools)

# add source subdirectory
add_subdirectory(source)
