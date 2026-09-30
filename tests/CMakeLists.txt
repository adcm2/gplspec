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
