add_executable(stage00_reference stage00_reference.cpp)
target_link_libraries(stage00_reference PRIVATE gplspec PlanetaryModel TomographyModels)

# Public umbrella headers must compile independently and remain ODR-safe when
# included by more than one translation unit.
foreach(header IN ITEMS Core GeneralModels SimpleModels All)
  string(TOLOWER "${header}" header_lower)
  add_library(stage02_${header_lower}_smoke OBJECT stage02_${header_lower}_smoke.cpp)
  target_link_libraries(stage02_${header_lower}_smoke PRIVATE gplspec)
endforeach()

add_executable(stage02_header_link stage02_header_main.cpp stage02_header_peer.cpp)
target_link_libraries(stage02_header_link PRIVATE gplspec)


if(BUILD_TESTING)
  find_package(PythonInterp 3 REQUIRED)
  add_test(NAME stage02_patch_application
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_stage02_patch_application.py)
endif()

add_executable(stage03_utilities stage03_utilities.cpp)
target_link_libraries(stage03_utilities PRIVATE gplspec)
if(BUILD_TESTING)
  add_test(NAME stage03_utilities COMMAND stage03_utilities)
endif()

add_executable(stage04_construction_storage stage04_construction_storage.cpp)
target_link_libraries(stage04_construction_storage PRIVATE gplspec)
if(BUILD_TESTING)
  add_test(NAME stage04_construction_storage COMMAND stage04_construction_storage)
endif()

add_executable(stage04_mapping_geometry stage04_mapping_geometry.cpp)
target_link_libraries(stage04_mapping_geometry PRIVATE gplspec)
if(BUILD_TESTING)
  add_test(NAME stage04_mapping_geometry COMMAND stage04_mapping_geometry)
endif()

add_executable(stage04_density_sampling stage04_density_sampling.cpp)
target_link_libraries(stage04_density_sampling PRIVATE gplspec)
if(BUILD_TESTING)
  add_test(NAME stage04_density_sampling COMMAND stage04_density_sampling)
endif()

add_executable(stage04_discretization_setup stage04_discretization_setup.cpp)
target_link_libraries(stage04_discretization_setup PRIVATE gplspec)
if(BUILD_TESTING)
  add_test(NAME stage04_discretization_setup COMMAND stage04_discretization_setup)
endif()
