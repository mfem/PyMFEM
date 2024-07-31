cmake_minimum_required(VERSION 3.15...3.26)
# project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)
project(pymfem LANGUAGES CXX)

find_package(
    Python
    COMPONENTS Interpreter Development.Module NumPy
    REQUIRED
)
find_package(nanobind CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # fPIC
# set(CMAKE_CXX_FLAGS -fPIC)
# set(MFEM_USE_ZLIB= 1)


# (Step 1) Add mfem (fetchcontent)
# include(FetchContent)
# FetchContent_Declare(
#     mfem
#     GIT_REPOSITORY https://github.com/mfem/mfem.git
#     GIT_TAG        master
# )
# FetchContent_MakeAvailable(mfem)

# (Step 1) Add mfem (pre-built)
include_directories("../mfem-shared/build")
link_directories("../mfem-shared/build")

# (Step 2) Create bindings
nanobind_add_module(
    pymfem
    NB_STATIC
    bindings/pymfem.cpp
    bindings/mesh.cpp
)

target_link_libraries(pymfem PRIVATE mfem rt)

install(TARGETS pymfem LIBRARY DESTINATION .)

# target_include_directories(mfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mfem)

# # (Step 1) Add mfem
# include(ExternalProject)
# ExternalProject_Add(
#     mfem
#     GIT_REPOSITORY https://github.com/mfem/mfem.git
#     GIT_TAG        master
#     SOURCE_DIR     "${CMAKE_BINARY_DIR}/mfem/"
#     BINARY_DIR     "${CMAKE_BINARY_DIR}/mfem/"
#     CMAKE_ARGS
#         -DBUILD_SHARED_LIBS=1
#         -DMFEM_ENABLE_EXAMPLES=1
#         -DMFEM_ENABLE_MINIAPPS=1
#         -DMFEM_USE_ZLIB=1
#         -DMFEM_USE_EXCEPTIONS=1
#         -DCMAKE_SHARED_LINKER_FLAGS=
#         -DCMAKE_CXX_FLAGS=-std=c++11
#         -DCMAKE_CXX_COMPILER=c++
#         -DCMAKE_MAKE_PROGRAM=/usr/bin/make
# )

# # cmake .. -DBUILD_SHARED_LIBS=1 -DMFEM_ENABLE_EXAMPLES=1 -DMFEM_ENABLE_MINIAPPS=1
# # -DCMAKE_SHARED_LINKER_FLAGS= -DMFEM_USE_ZLIB=1 -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_CXX_COMPILER=c++ -DMFEM_USE_EXCEPTIONS=1
# # -DCMAKE_INSTALL_PREFIX=/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser
# # -DCMAKE_INSTALL_RPATH=/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/lib
# # -DCMAKE_BUILD_WITH_INSTALL_RPATH=1

# # /home/jgl/.local/bin/swig -Wall -c++ -python -fastproxy -olddefs -keyword -I/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/include -I/home/jgl/.local/lib/python3.8/site-packages/mfem/external/ser/include/mfem -I/home/jgl/git/WORK/pymfem-temp/external/mfem lininteg.i