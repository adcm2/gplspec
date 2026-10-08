#include "../gplspec/src/GeneralModels/Model_Construction_Utilities.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using ScalarField = std::vector<std::vector<std::vector<double>>>;

void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

ScalarField original_storage(std::size_t elements, std::size_t nodes,
                             std::size_t spatial, double initial) {
   return ScalarField(elements,
                      std::vector<std::vector<double>>(
                          nodes, std::vector<double>(spatial, initial)));
}

void check_case(std::size_t elements, std::size_t nodes, std::size_t spatial,
                double initial) {
   const auto expected = original_storage(elements, nodes, spatial, initial);
   const auto actual = GPLSpec::model_detail::InitializeScalarFieldStorage(
       elements, nodes, spatial, initial);
   require(actual == expected, "nested scalar storage differs from original initializer");
}
}  // namespace

int main() {
   check_case(0, 4, 6, 0.0);
   check_case(3, 0, 6, 0.0);
   check_case(2, 4, 0, 0.0);
   check_case(1, 3, 5, 0.0);
   check_case(2, 4, 7, 1.0);
   check_case(3, 2, 4, -2.5);
}
