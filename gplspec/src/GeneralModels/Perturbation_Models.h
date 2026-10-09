#ifndef PERTURBATION_MODEL_3D_H
#define PERTURBATION_MODEL_3D_H

#include "../Radial_Tools.h"
#include "Density_Model_Return.h"
#include "Earth_Density_Models_3D.h"
#include "Mapping_Perturbation_Utilities.h"
#include "../Model_Input.h"
#include <GSHTrans/All>
#include <PlanetaryModel/All>
#include <algorithm>
#include <complex>
#include <vector>

namespace GeneralEarthModels {

class MappingPerturbation {
 public:
   using TensorField = MappingPerturbationDetail::TensorField;

   explicit MappingPerturbation(const Density3D &model)
       : _delta_tensor(ZeroTensor(model)) {}

   template <class Map>
      requires PlanetaryModel::RadialMappingClass<Map>
   MappingPerturbation(const Density3D &model, const Map &mapping)
       : _delta_tensor(BuildRadialMapping(model, mapping)) {}

   MappingPerturbation(const Density3D &model,
                       const GPLSpec::ModelInput::PhobosHarmonics &harmonics,
                       int lmin, int lmax)
       : _delta_tensor(BuildFileMapping(model, harmonics, lmin, lmax)) {}

   const TensorField &delta_laplace_tensor() const { return _delta_tensor; }

 private:
   static TensorField ZeroTensor(const Density3D &model) {
      const auto &mesh = model.Node_InformationP();
      const auto &grid = model.GSH_Grid();
      const auto spatial = grid.NumberOfCoLatitudes() *
                           grid.NumberOfLongitudes();
      return TensorField(
          mesh.NumberOfElements(),
          std::vector<std::vector<Eigen::Matrix3cd>>(
              mesh.NumNodesWithinElement(),
              std::vector<Eigen::Matrix3cd>(spatial,
                                             Eigen::Matrix3cd::Zero())));
   }

   template <class Map>
      requires PlanetaryModel::RadialMappingClass<Map>
   static TensorField BuildRadialMapping(const Density3D &model,
                                         const Map &mapping) {
      const auto &mesh = model.Node_InformationP();
      const auto &grid = model.GSH_Grid();
      const std::size_t elements = mesh.NumberOfElements();
      const std::size_t nodes = mesh.NumNodesWithinElement();
      const std::size_t spatial =
          grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();
      MappingPerturbationDetail::SpatialDisplacement displacement(
          elements,
          std::vector<std::vector<Eigen::Vector3cd>>(
              nodes, std::vector<Eigen::Vector3cd>(
                         spatial, Eigen::Vector3cd::Zero())));
      const double planet_radius = mesh.PlanetRadius();
      const double outer_radius = mesh.OuterRadius();

      for (std::size_t element = 0; element < elements; ++element) {
         const int layer = mesh.LayerNumber(element);
         const auto radial_map = mapping.RadialMapping(layer);
         for (std::size_t node = 0; node < nodes; ++node) {
            const double radius = mesh.NodeRadius(element, node);
            const bool exterior = radius > planet_radius;
            const double evaluation_radius =
                exterior ? planet_radius : radius;
            const double taper =
                exterior ? (outer_radius - radius) /
                              (outer_radius - planet_radius)
                         : 1.0;
            std::size_t point = 0;
            for (double theta : grid.CoLatitudes()) {
               for (double phi : grid.Longitudes()) {
                  displacement[element][node][point++](1) =
                      radial_map(evaluation_radius, theta, phi) * taper;
               }
            }
         }
      }
      return MappingPerturbationDetail::BuildPerturbationTensor(model,
                                                                 displacement);
   }

   static TensorField BuildFileMapping(
       const Density3D &model,
       const GPLSpec::ModelInput::PhobosHarmonics &harmonics, int lmin,
       int requested_lmax) {
      const auto &mesh = model.Node_InformationP();
      const auto &grid = model.GSH_Grid();
      const std::size_t elements = mesh.NumberOfElements();
      const std::size_t nodes = mesh.NumNodesWithinElement();
      const std::size_t spatial =
          grid.NumberOfCoLatitudes() * grid.NumberOfLongitudes();
      const int lmax = grid.MaxDegree();
      const std::size_t nonnegative_count =
          GSHTrans::GSHIndices<GSHTrans::NonNegative>(lmax, lmax, 0).Size();
      const std::size_t model_lmax = harmonics.InferredMaximumDegree();
      const int maximum_degree = std::min(
          {lmax, static_cast<int>(model_lmax), requested_lmax});

      std::vector<std::complex<double>> coefficients(nonnegative_count,
                                                     std::complex<double>{});
      std::size_t coefficient = 0;
      for (int degree = 0; degree < lmin; ++degree) {
         coefficient += static_cast<std::size_t>(degree + 1);
      }
      const std::complex<double> imaginary{0.0, 1.0};
      for (int degree = lmin; degree <= maximum_degree; ++degree) {
         for (int order = 0; order <= degree; ++order) {
            const double normalization = order == 0 ? 1.0 : 0.5;
            coefficients[coefficient] =
                normalization *
                (harmonics.rows[coefficient].cosine -
                 imaginary * harmonics.rows[coefficient].sine);
            ++coefficient;
         }
      }

      std::vector<double> surface_displacement(spatial, 0.0);
      grid.InverseTransformation(lmax, 0, coefficients,
                                 surface_displacement);
      MappingPerturbationDetail::SpatialDisplacement displacement(
          elements,
          std::vector<std::vector<Eigen::Vector3cd>>(
              nodes, std::vector<Eigen::Vector3cd>(
                         spatial, Eigen::Vector3cd::Zero())));
      const double planet_radius = mesh.PlanetRadius();
      const double outer_radius = mesh.OuterRadius();
      for (std::size_t element = 0; element < elements; ++element) {
         for (std::size_t node = 0; node < nodes; ++node) {
            const double radius = mesh.NodeRadius(element, node);
            const double taper =
                radius > planet_radius
                    ? (outer_radius - radius) / (outer_radius - planet_radius)
                    : radius / planet_radius;
            for (std::size_t point = 0; point < spatial; ++point) {
               displacement[element][node][point](1) =
                   surface_displacement[point] / model.LengthNorm() * taper;
            }
         }
      }
      return MappingPerturbationDetail::BuildPerturbationTensor(model,
                                                                 displacement);
   }

   TensorField _delta_tensor;
};

}  // namespace GeneralEarthModels

#endif
