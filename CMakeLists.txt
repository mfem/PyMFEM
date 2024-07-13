cmake_minimum_required(VERSION 3.15...3.26)
# project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)
project(pymfem_dev)
include(ExternalProject)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../mfem-build mfem)
message(STATUS "SKBUILD_DATA_DIR : ${SKBUILD_DATA_DIR}")

# (Step 1) Add mfem
ExternalProject_Add(
    mfem
    GIT_REPOSITORY https://github.com/mfem/mfem.git
    CMAKE_ARGS -DBUILD_SHARED_LIBS=1 -DMFEM_ENABLE_EXAMPLES=1 -DMFEM_ENABLE_MINIAPPS=1 -DCMAKE_SHARED_LINKER_FLAGS= -DMFEM_USE_ZLIB=1 -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CXX_COMPILER=c++ -DMFEM_USE_EXCEPTIONS=1# -DCMAKE_INSTALL_PREFIX=/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser -DCMAKE_INSTALL_RPATH=/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/lib -DCMAKE_BUILD_WITH_INSTALL_RPATH=1
)
# cmake .. -DBUILD_SHARED_LIBS=1 -DMFEM_ENABLE_EXAMPLES=1 -DMFEM_ENABLE_MINIAPPS=1 -DCMAKE_SHARED_LINKER_FLAGS= -DMFEM_USE_ZLIB=1 -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CXX_COMPILER=c++ -DMFEM_USE_EXCEPTIONS=1

# /home/jgl/.local/bin/swig -Wall -c++ -python -fastproxy -olddefs -keyword -I/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/include -I/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/include/mfem -I/home/jgl/git/WORK/pymfem-temp/external/mfem lininteg.i