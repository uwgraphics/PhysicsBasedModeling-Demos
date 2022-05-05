cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 14)

project (Collision-Tests)

set (CMAKE_CXX_FLAGS "-fopenmp -Wno-deprecated ${CMAKE_CXX_FLAGS}")

add_definitions(-DBUILD_COMPONENT_SRC_PREFIX="" -DBUILD_OPTLEVEL_DEV)

set(USD_ROOT_DIR /nobackup/sifakis/opt)
set(USD_LIB_DIR ${USD_ROOT_DIR}/lib/)
set(USD_INC_DIR ${USD_ROOT_DIR}/include/)

add_subdirectory(tests)
