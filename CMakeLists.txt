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
# set(CMAKE_CXX_FLAGS "-Ofast")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/FindLAPACKE.cmake")

# Locate  FFTW3


# add in Eigen3
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_MODULE_PATH "/usr/local/include/eigen/cmake/;${CMAKE_MODULE_PATH}")
find_package(FFTW)


# find_package (FFTWpp)
include(FetchContent)
FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY git@github.com:da380/FFTWpp
  GIT_TAG main
)
FetchContent_MakeAvailable(FFTWpp)

include(FetchContent)
FetchContent_Declare(
  GaussQuad
  GIT_REPOSITORY git@github.com:da380/GaussQuad
  GIT_TAG main
)
FetchContent_MakeAvailable(GaussQuad)

include(FetchContent)
FetchContent_Declare(
  PlanetaryModel
  GIT_REPOSITORY git@github.com:da380/PlanetaryModel
  GIT_TAG newdev
)
FetchContent_MakeAvailable(PlanetaryModel)

include(FetchContent)
FetchContent_Declare(
  GSHTRANS
  GIT_REPOSITORY git@github.com:da380/GSHTrans
  GIT_TAG main
)
FetchContent_MakeAvailable(GSHTRANS)


include(FetchContent)
FetchContent_Declare(
  TomographyModels
  GIT_REPOSITORY git@github.com:adcm2/TomographyModels
  GIT_TAG main
)
FetchContent_MakeAvailable(TomographyModels)

include(FetchContent)
FetchContent_Declare(
  Interpolation
  GIT_REPOSITORY git@github.com:da380/Interpolation.git
  GIT_TAG develop
)
FetchContent_MakeAvailable(Interpolation)

# Set up the library
add_library (gplspec INTERFACE)
target_link_libraries(gplspec INTERFACE Eigen3::Eigen)
# target_link_libraries(Gravitational_Field INTERFACE FFTWpp)
# target_link_libraries(Gravitational_Field INTERFACE GaussQuad)
# target_link_libraries(Gravitational_Field INTERFACE PlanetaryModel)
# target_link_libraries(Gravitational_Field INTERFACE GSHTRANS)
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
