cmake_minimum_required (VERSION 3.5)

# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "g++")


project(gplspec)
# get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS  "-march=native")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2b -fconcepts-diagnostics-depth=2")
# set(CMAKE_CXX_FLAGS "-O2")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/FindLAPACKE.cmake")

# Locate  FFTW3


# Exact Eigen package used by the verified baseline; the source shim also
# satisfies nested GSHTrans FetchContent calls without fetching another copy.
include(FetchContent)
find_package(Eigen3 3.4.0 EXACT CONFIG REQUIRED)
FetchContent_Declare(
  Eigen3
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/eigen-system"
)

find_package(FFTW)

# Declare the complete dependency closure before populating any project. This
# ensures these immutable revisions win over transitive branch declarations.
FetchContent_Declare(
  NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG 888126b44a979fda5dcc9d0b11d80cb8b3f61b52
)
FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY https://github.com/da380/FFTWpp.git
  GIT_TAG 06cc1fb04c4398407839e637cb5414b179b0e9c3
)
FetchContent_Declare(
  GaussQuad
  GIT_REPOSITORY https://github.com/da380/GaussQuad.git
  GIT_TAG fbe37c7eef93695317dbf9c48715110bf64c800f
)
FetchContent_Declare(
  Interpolation
  GIT_REPOSITORY https://github.com/da380/Interpolation.git
  GIT_TAG 1557aad571ceac1899c2d3795696aa4047f41d8a
)
FetchContent_Declare(
  PlanetaryModel
  GIT_REPOSITORY https://github.com/da380/PlanetaryModel.git
  GIT_TAG a801b9e9c2205c3b32daa745fd523d596a9d2ede
)
FetchContent_Declare(
  GSHTRANS
  GIT_REPOSITORY https://github.com/adcm2/GSHTrans.git
  GIT_TAG f0a0e24a3579ceb7c99f8782328cae14bef195de
)
FetchContent_Declare(
  TomographyModels
  GIT_REPOSITORY https://github.com/adcm2/TomographyModels.git
  GIT_TAG 97869aeea3f3d87901705db59b5de975960bd9e1
)

# Populate in the dependency order required by the interface targets.
FetchContent_MakeAvailable(Eigen3)
FetchContent_MakeAvailable(NumericConcepts)
FetchContent_MakeAvailable(FFTWpp)
FetchContent_MakeAvailable(GaussQuad)
FetchContent_MakeAvailable(Interpolation)
FetchContent_MakeAvailable(PlanetaryModel)
FetchContent_MakeAvailable(GSHTRANS)
FetchContent_MakeAvailable(TomographyModels)

# Set up the library
add_library (gplspec INTERFACE)
target_link_libraries(gplspec INTERFACE Eigen3::Eigen)
# target_link_libraries(Gravitational_Field INTERFACE FFTWpp)
target_link_libraries(gplspec INTERFACE GaussQuad)
target_link_libraries(gplspec INTERFACE GSHTrans)
# target_link_libraries(gplspec INTERFACE PlanetaryModel)
# target_link_libraries(gplspec INTERFACE GSHTRANS)
# target_link_libraries(Gravitational_Field INTERFACE TomographyModels)
# target_link_libraries(Gravitational_Field INTERFACE Interpolation)
target_include_directories (gplspec INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)



# blas and lapack
find_package(OpenMP REQUIRED)





# optionally add in the examples and tests
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # We're in the root, define additional targets for developers.
    option(MY_PROJECT_BUILD_EXAMPLES   "whether or not examples should be built" ON)
    option(GPLSPEC_BUILD_BASELINE_HARNESS "build deterministic stage-00 baseline harness" OFF)
    if(GPLSPEC_BUILD_BASELINE_HARNESS)
      add_subdirectory(tests)
    endif()
    # option(MY_PROJECT_BUILD_TESTS      "whether or not tests should be built" ON)


    
    # include(FetchContent)
    # FetchContent_Declare(
    #   Interpolation
    #   GIT_REPOSITORY https://github.com/da380/Interpolation.git
    #   GIT_TAG main
    #   )
    # FetchContent_MakeAvailable(Interpolation)

    if(MY_PROJECT_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
    # if(MY_PROJECT_BUILD_TESTS)
    #     enable_testing()
    #     add_subdirectory(tests)
    # endif()
endif()

# add_subdirectory(Gravitational_Field)
# add_subdirectory(examples)
