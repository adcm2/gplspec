#include "../gplspec/src/GeneralModels/Model_Construction_Utilities.h"
#include <GSHTrans/All>

#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

template <class Range>
std::vector<double> as_vector(const Range &range) {
   return {range.begin(), range.end()};
}

using Grid = GSHTrans::GaussLegendreGrid<double, GSHTrans::All, GSHTrans::All>;

void check_case(int requested_order, int l_max) {
   GaussQuad::Quadrature1D<double> expected_quadrature;
   Grid expected_grid;
   int expected_polynomial_order = -1;
   expected_quadrature =
       GaussQuad::GaussLobattoLegendreQuadrature1D<double>(requested_order + 1);
   expected_grid = Grid(l_max, 2);
   expected_polynomial_order = expected_quadrature.N() - 1;

   GaussQuad::Quadrature1D<double> actual_quadrature;
   Grid actual_grid;
   int actual_polynomial_order = -1;
   GPLSpec::model_detail::InitializeModelDiscretization(
       actual_quadrature, actual_grid, actual_polynomial_order,
       requested_order, l_max);

   require(actual_polynomial_order == expected_polynomial_order,
           "polynomial order differs from original setup");
   require(actual_quadrature.N() == expected_quadrature.N(),
           "quadrature node count differs from original setup");
   require(as_vector(actual_quadrature.Points()) ==
               as_vector(expected_quadrature.Points()),
           "quadrature points differ from original setup");
   require(as_vector(actual_quadrature.Weights()) ==
               as_vector(expected_quadrature.Weights()),
           "quadrature weights differ from original setup");
   require(actual_grid.NumberOfLongitudes() ==
               expected_grid.NumberOfLongitudes() &&
               actual_grid.NumberOfCoLatitudes() ==
                   expected_grid.NumberOfCoLatitudes(),
           "transform-grid extents differ from original setup");
   require(as_vector(actual_grid.Longitudes()) ==
               as_vector(expected_grid.Longitudes()) &&
               as_vector(actual_grid.CoLatitudes()) ==
                   as_vector(expected_grid.CoLatitudes()),
           "transform-grid coordinates differ from original setup");
}
}  // namespace

int main() {
   check_case(2, 2);
   check_case(3, 3);
   check_case(4, 4);
}
