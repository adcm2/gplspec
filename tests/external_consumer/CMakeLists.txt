cmake_minimum_required(VERSION 3.20)
project(gplspec_external_consumer LANGUAGES CXX)

if(NOT GPLSPEC_SOURCE_DIR OR NOT GPLSPEC_EIGEN_INCLUDE_DIR)
  message(FATAL_ERROR "GPLSPEC_SOURCE_DIR and GPLSPEC_EIGEN_INCLUDE_DIR are required")
endif()

add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
set_property(TARGET Eigen3::Eigen PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES "${GPLSPEC_EIGEN_INCLUDE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${GPLSPEC_SOURCE_DIR}/cmake")
add_subdirectory("${GPLSPEC_SOURCE_DIR}" gplspec)

add_executable(external_consumer main.cpp peer.cpp)
target_compile_features(external_consumer PRIVATE cxx_std_23)
target_link_libraries(external_consumer PRIVATE gplspec)
