#ifndef GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H
#define GPLSPEC_MODEL_CONSTRUCTION_UTILITIES_H

#include <cstddef>
#include <vector>

namespace GPLSpec::model_detail {

// Allocate an element/node/spatial scalar field with the requested initial value.
inline std::vector<std::vector<std::vector<double>>> InitializeScalarFieldStorage(
    std::size_t element_count, std::size_t node_count,
    std::size_t spatial_count, double initial_value) {
   return std::vector<std::vector<std::vector<double>>>(
       element_count,
       std::vector<std::vector<double>>(
           node_count, std::vector<double>(spatial_count, initial_value)));
}

}  // namespace GPLSpec::model_detail

#endif
