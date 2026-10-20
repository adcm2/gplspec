#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#ifndef GPLSPEC_SOURCE_DIR
#error "GPLSPEC_SOURCE_DIR must identify the fixture tree"
#endif

namespace {

void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

bool close(double actual, double expected, double scale) {
   constexpr double rtol = 1e-12;
   constexpr double atol_factor = 1e-14;
   return std::abs(actual - expected) <= rtol * std::abs(expected) +
                                         atol_factor * scale;
}

struct SingleLayerModel {
   double LengthNorm() const { return 1.0; }
   double MassNorm() const { return 1.0; }
   double TimeNorm() const { return 1.0; }
   int NumberOfLayers() const { return 1; }
   double LowerRadius(int) const { return 0.0; }
   double UpperRadius(int) const { return 1.0; }
   double OuterRadius() const { return 1.0; }
   auto Density(int) const { return [](double) { return 1.0; }; }
};

}  // namespace

int main() {
   constexpr int poly_order = 3;
   constexpr int lmax = 2;
   constexpr double max_step = 0.5;
   constexpr double outer_radius = 1.2;
   const std::string prem_path =
       std::string(GPLSPEC_SOURCE_DIR) + "/modeldata/PREM500";
   const auto prem = GPLSpec::ModelInput::ReadPremDeck(prem_path);
   GeneralEarthModels::RadialEarthModel<> input(prem);
   const auto q =
       GaussQuad::GaussLobattoLegendreQuadrature1D<double>(poly_order + 1);
   GeneralEarthModels::spherical_1D radial(input, q, max_step, outer_radius);
   auto density = GeneralEarthModels::Density3D::OneDimensionalPlanetFromModel(
       input, poly_order, lmax, max_step, outer_radius);

   const auto &mesh = radial.Node_Information();
   const auto &density_mesh = density.Node_Information();
   const double radius_scale = outer_radius;
   require(mesh.NumberOfElements() == density_mesh.NumberOfElements(),
           "radial owners disagree on element count");
   require(mesh.PolynomialOrder() == static_cast<std::size_t>(poly_order) &&
               density_mesh.PolynomialOrder() == mesh.PolynomialOrder(),
           "radial owners disagree on polynomial order");

   std::size_t interface_element = mesh.NumberOfElements();
   std::size_t exterior_element = mesh.NumberOfElements();
   for (std::size_t element = 0; element < mesh.NumberOfElements(); ++element) {
      const int layer = mesh.LayerNumber(element);
      require(layer == density_mesh.LayerNumber(element),
              "radial owners disagree on the element material layer");
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         const double radius = mesh.NodeRadius(element, node);
         require(close(radius, density_mesh.NodeRadius(element, node),
                       radius_scale),
                 "radial owners disagree on an element quadrature radius");
         const double expected_density =
             layer < input.NumberOfLayers()
                 ? input.Density(layer)(radius)
                 : 0.0;
         const double expected_vpv =
             layer < input.NumberOfLayers() ? input.VPV(layer)(radius) : 0.0;
         const double expected_vsv =
             layer < input.NumberOfLayers() ? input.VSV(layer)(radius) : 0.0;
         const double density_scale = std::max(1.0, std::abs(expected_density));
         const double velocity_scale = std::max(1.0, std::abs(expected_vpv));

         require(close(radial.Density(element, node), expected_density,
                       density_scale),
                 "radial density samples disagree with the input model");
         require(close(radial.VPV(element, node), expected_vpv,
                       velocity_scale) &&
                     close(radial.VSV(element, node), expected_vsv,
                           std::max(1.0, std::abs(expected_vsv))),
                 "radial elastic samples disagree with the input model");

         const auto density_samples = density.DensityAtRadialNode(element, node);
         for (const auto sample : density_samples) {
            require(close(sample, expected_density, density_scale),
                    "Density3D radial density differs from the input model");
         }
         if (layer == static_cast<int>(input.NumberOfLayers())) {
            exterior_element = element;
            require(expected_density == 0.0 && expected_vpv == 0.0 &&
                        std::all_of(density_samples.begin(),
                                    density_samples.end(),
                                    [](double value) { return value == 0.0; }),
                    "exterior material did not remain zero valued");
         }
      }
      if (element > 0 && mesh.LayerNumber(element - 1) != layer) {
         interface_element = element;
      }
   }

   require(interface_element < mesh.NumberOfElements(),
           "PREM fixture has no material interface");
   require(exterior_element < mesh.NumberOfElements(),
           "radial fixture has no exterior element");
   const std::size_t last_node = mesh.PolynomialOrder();
   require(mesh.GlobalNodeIndex(interface_element - 1, last_node) ==
               mesh.GlobalNodeIndex(interface_element, 0),
           "one-sided material samples do not share a global interface DOF");
   const double left_density = radial.Density(interface_element - 1, last_node);
   const double right_density = radial.Density(interface_element, 0);
   require(left_density != right_density,
           "PREM interface did not preserve distinct one-sided densities");

   auto matrix_owner = radial.SpectralElementInformation();
   const auto weak_form = matrix_owner.fullmatrix<double>(mesh, 0, 0);
   require(weak_form.rows() ==
               static_cast<Eigen::Index>(mesh.NumberOfGlobalRadialNodes()),
           "weak-form matrix dimensions disagree with canonical radial mesh");

   const SingleLayerModel single_layer;
   const Radial_Tools::RadialMesh one_element_mesh(single_layer, q, 2.0, 1.2,
                                                   false);
   require(one_element_mesh.NumberOfElements() == 1,
           "single-layer fixture did not produce one radial element");
   const auto one_element_matrix =
       matrix_owner.fullmatrix<double>(one_element_mesh, 0, 0);
   require(one_element_matrix.rows() ==
                   static_cast<Eigen::Index>(
                       one_element_mesh.NumberOfGlobalRadialNodes()) &&
               one_element_matrix.norm() > 0.0,
           "weak-form assembly rejected or emptied a one-element mesh");
}
