#include <gplspec/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
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

template <class Mesh>
void check_radial_layout(const Mesh &mesh) {
   const std::size_t order = mesh.PolynomialOrder();
   const std::size_t elements = mesh.NumberOfElements();
   require(mesh.NumberOfGlobalRadialNodes() == elements * order + 1,
           "radial global-node count changed");
   for (std::size_t element = 0; element < elements; ++element) {
      for (std::size_t node = 0; node < mesh.NumNodesWithinElement(); ++node) {
         require(mesh.GlobalNodeIndex(element, node) == element * order + node,
                 "radial element-local node has the wrong global index");
      }
      if (element + 1 < elements) {
         require(mesh.GlobalNodeIndex(element, order) ==
                     mesh.GlobalNodeIndex(element + 1, 0),
                 "neighboring element endpoints do not share a global node");
      }
   }
}

double radial_weight(const Radial_Tools::RadialMesh &mesh,
                     const GaussQuad::Quadrature1D<double> &q,
                     std::size_t element, std::size_t node) {
   const double lower = mesh.ElementLowerRadius(element);
   const double upper = mesh.ElementUpperRadius(element);
   const double radius =
       0.5 * (upper + lower) + 0.5 * (upper - lower) * q.X(node);
   return q.W(node) * 0.5 * (upper - lower) * radius * radius;
}

}  // namespace

int main() {
   using Complex = std::complex<double>;
   constexpr int lmax = 2;
   constexpr int stride = (lmax + 1) * (lmax + 1);
   constexpr double pi = 3.1415926535;

   const auto quadrature =
       GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
   GeneralEarthModels::spherical_1D radial_model(
       std::string(GPLSPEC_SOURCE_DIR) + "/modeldata/PREM500", quadrature,
       0.5, 1.2);
   const auto &radial_mesh = radial_model.Node_Information();
   check_radial_layout(radial_mesh);

   const Eigen::VectorXcd monopole = Gravity_Tools::FindForce(radial_model, lmax);
   Eigen::VectorXcd expected_monopole = Eigen::VectorXcd::Zero(monopole.size());
   const double one_dimensional_scale =
       4.0 * pi * radial_model.GravitationalConstant() * std::sqrt(4.0 * pi);
   for (std::size_t element = 0; element < radial_mesh.NumberOfElements();
        ++element) {
      for (std::size_t node = 0; node < radial_mesh.NumNodesWithinElement();
           ++node) {
         const auto radial = radial_mesh.GlobalNodeIndex(element, node);
         expected_monopole(radial * stride) +=
             one_dimensional_scale *
             radial_weight(radial_mesh, quadrature, element, node) *
             radial_model.Density(element, node);
      }
   }
   require(monopole.size() == expected_monopole.size() &&
               (monopole - expected_monopole).norm() <=
                   1e-12 * expected_monopole.norm() +
                       1e-14 * expected_monopole.norm(),
           "1-D monopole source disagrees with independent radial quadrature");
   for (Eigen::Index index = 0; index < monopole.size(); ++index) {
      if (index % stride != 0) {
         require(monopole(index) == Complex(0.0, 0.0),
                 "1-D source populated a non-monopole harmonic");
      }
   }

   auto layered = GeneralEarthModels::Density3D::SphericalLayeredPlanet(
       {0.5, 1.0}, {2.0, 1.0}, 1.0, 1.0, 1.0, 0.25, 1.2);
   const auto &layered_mesh = layered.Node_Information();
   check_radial_layout(layered_mesh);
   const Eigen::VectorXcd layered_source = Gravity_Tools::FindForce(layered);
   require(layered_source.size() ==
               static_cast<Eigen::Index>(layered_mesh.NumberOfGlobalRadialNodes() *
                                         stride),
           "3-D source shape differs from the canonical radial/harmonic layout");

   std::size_t interface_element = layered_mesh.NumberOfElements();
   for (std::size_t element = 1; element < layered_mesh.NumberOfElements();
        ++element) {
      if (layered_mesh.LayerNumber(element - 1) !=
          layered_mesh.LayerNumber(element)) {
         interface_element = element;
         break;
      }
   }
   require(interface_element < layered_mesh.NumberOfElements(),
           "layered source fixture has no material interface");
   const std::size_t order = layered_mesh.PolynomialOrder();
   const std::size_t interface_dof =
       layered_mesh.GlobalNodeIndex(interface_element - 1, order);
   require(interface_dof == layered_mesh.GlobalNodeIndex(interface_element, 0),
           "material interface is not one shared global radial degree of freedom");
   const double left_density = layered.DensityAtRadialNode(interface_element - 1,
                                                           order)[0];
   const double right_density = layered.DensityAtRadialNode(interface_element, 0)[0];
   require(std::abs(left_density - right_density) > 1e-12,
           "layered source fixture lost its one-sided density values");

   const auto layered_q = layered.q();
   const double interface_expected =
       4.0 * std::numbers::pi_v<double> * layered.GravitationalConstant() *
       std::sqrt(4.0 * std::numbers::pi_v<double>) *
       (radial_weight(layered_mesh, layered_q, interface_element - 1, order) *
            left_density +
        radial_weight(layered_mesh, layered_q, interface_element, 0) *
            right_density);
   const Eigen::Index interface_source_index = interface_dof * stride;
   require(std::abs(layered_source(interface_source_index).real() -
                    interface_expected) <=
               1e-12 * std::abs(interface_expected) +
                   1e-14 * std::abs(interface_expected),
           "3-D source did not accumulate both one-sided interface contributions");

   double nonmonopole_norm = 0.0;
   for (std::size_t radial = 0;
        radial < layered_mesh.NumberOfGlobalRadialNodes(); ++radial) {
      const Eigen::Index base = radial * stride;
      require(layered_source(base + lmax * lmax) == Complex(0.0, 0.0) &&
                  layered_source(base + lmax * lmax + 2 * lmax) ==
                      Complex(0.0, 0.0),
              "3-D source endpoint harmonic slots are not zero");
      for (int harmonic = 1; harmonic < stride; ++harmonic) {
         nonmonopole_norm += std::norm(layered_source(base + harmonic));
      }
   }
   for (int harmonic = 1; harmonic < stride; ++harmonic) {
      require(layered_source(harmonic) == Complex(0.0, 0.0),
              "3-D source origin retained a non-monopole coefficient");
   }
   const double source_scale = layered_source.norm();
   require(std::sqrt(nonmonopole_norm) <=
               1e-11 * source_scale + 1e-13 * source_scale,
           "spherical layered source has nonphysical non-monopole content");
}
