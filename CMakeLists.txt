cmake_minimum_required(VERSION 3.2)
project(vuptool VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)

find_package(savvy REQUIRED)

add_executable(vuptool src/main.cpp src/inv_norm.cpp)
target_link_libraries(vuptool savvy)

target_compile_definitions(vuptool PUBLIC -DVUPTOOL_VERSION="${PROJECT_VERSION}")
