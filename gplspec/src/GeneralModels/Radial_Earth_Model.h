#ifndef GPLSPEC_RADIAL_EARTH_MODEL_H
#define GPLSPEC_RADIAL_EARTH_MODEL_H

#include "../Model_Input.h"
#include <Interpolation/All>
#include <PlanetaryModel/All>
#include <cassert>
#include <concepts>
#include <stdexcept>
#include <string>
#include <vector>

namespace GeneralEarthModels {

template <typename FLOAT = double, typename INTEGRAL = int>
   requires std::same_as<FLOAT, double>
class RadialEarthModel {
 public:
   using size_type = INTEGRAL;
   using value_type = FLOAT;
   // using InterpA = Interpolation::Akima<std::vector<double>::iterator,
   //                                      std::vector<double>::iterator>;
   // template <typename NFLOAT = FLOAT>
   // using InterpA = Interpolation::CubicSpline<std::vector<double>::iterator,
   //                                            std::vector<double>::iterator>;
   using myvector = std::vector<FLOAT>;
   using myiter = myvector::iterator;
   using InterpA = Interpolation::CubicSpline<myiter, myiter>;
   // using InterpA = Interpolation::CubicSpline<std::vector<FLOAT>::iterator,
   //                                            std::vector<FLOAT>::iterator>;
   // This model stores double-valued radial profiles. Returned spline objects
   // borrow its layer vectors, and derived A/C/N/L/F/Kappa/Mu functions capture
   // this owner. Discard returned functions before moving or destroying it.
   // Copying is disabled because internal splines refer to owned storage.
   RadialEarthModel(const RadialEarthModel &) = delete;
   RadialEarthModel &operator=(const RadialEarthModel &) = delete;
   RadialEarthModel(RadialEarthModel &&) = default;
   RadialEarthModel &operator=(RadialEarthModel &&) = default;

   // constructors
   explicit RadialEarthModel(const GPLSpec::ModelInput::PremDeck &);
   template <template <typename> class ParameterModel>
   RadialEarthModel(const GPLSpec::ModelInput::PremDeck &,
              const ParameterModel<FLOAT> &);

   // norms
   FLOAT LengthNorm() const { return length_norm; };
   FLOAT MassNorm() const { return mass_norm; };
   FLOAT TimeNorm() const { return time_norm; }
   FLOAT DensityNorm() const { return density_norm; };
   FLOAT InertiaNorm() const { return inertia_norm; };
   FLOAT VelocityNorm() const { return velocity_norm; };
   FLOAT AccelerationNorm() const { return acceleration_norm; };
   FLOAT ForceNorm() const { return force_norm; };
   FLOAT StressNorm() const { return stress_norm; };
   FLOAT GravitationalConstant() const { return gravitational_constant; };

   // Geometry of model
   int NumberOfLayers() const { return _numlayers; };
   int LayerLowerIndex(int i) const { return _vec_indices[i][0]; };
   int LayerUpperIndex(int i) const { return _vec_indices[i][1]; };
   int LayerIndexDifference(int i) const {
      return _vec_indices[i][1] - _vec_indices[i][0];
   };
   auto LayerRadii() const { return layered_radii; };
   auto LayerRadii(int i) const { return layered_radii[i]; };
   auto LayerRadiiNumber(int i) const { return layered_radii[i].size(); };

   FLOAT LowerRadius(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside the number of layers in the model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return vec_layers[i];
   }
   FLOAT UpperRadius(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside the number of layers in the model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return vec_layers[i + 1];
   }
   FLOAT OuterRadius() const { return vec_layers[_numlayers]; }

   // Isotropy/fluid/solid etc
   bool IsIsotropic() const { return _isisotropic; };

   // Solid or fluid
   bool IsSolid(INTEGRAL i) const { return _issolid[i]; }
   bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

   // Density
   InterpA Density(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_rho[i];
      // return func_rhoc[i];
   };

   // Velocities
   // InterpA VP(INTEGRAL i) const { return func_vpv[i]; };
   InterpA VPV(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };

      return func_vpv[i];
   };
   InterpA VPH(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vph[i];
   };
   // InterpA VS(INTEGRAL i) const { return vec_s_velocity[i]; };
   InterpA VSV(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsv[i];
   };
   InterpA VSH(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsh[i];
   };

   ///////////////////////////////////////////////////////////////
   ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
   ///////////////////////////////////////////////////////////////
   InterpA VS(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vsv[i];
   };
   InterpA VP(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_vpv[i];
   };
   ///////////////////////////////////////////////////////////////
   ////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!!! ////////////////
   ///////////////////////////////////////////////////////////////

   // Returning eta:
   auto Eta(INTEGRAL i) const {
      if (i < 0) {
         throw std::invalid_argument("Negative layer index");
      } else if (i > _numlayers - 1) {
         assert("Outside model");
         throw std::invalid_argument(
             "Layer index greater than number of layers");
      };
      return func_eta[i];
   }

   // returning A, C, N, L, kappa, mu
   auto A(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
      };
      return aret;
   };
   auto C(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
      };
      return aret;
   };
   auto N(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
      };
      return aret;
   };
   auto L(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
      };
      return aret;
   };
   auto F(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x));
      };
      return aret;
   };
   auto Kappa(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
      };
      return aret;
   };
   auto Mu(INTEGRAL i) const {
      auto aret = [i, this](FLOAT x) {
         return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
                 2.0 * F(i)(x)) /
                15.0;
      };
      return aret;
   };

 private:
   using vecdb = std::vector<double>;
   using vvecdb = std::vector<vecdb>;
   vecdb vec_radius, vec_rho, vec_vpv, vec_vsv, vec_qkappa, vec_qshear, vec_vph,
       vec_vsh, vec_eta, vec_layers;
   vvecdb layered_radii = vvecdb(1, vecdb());
   vvecdb layered_rho = vvecdb(1, vecdb());
   vvecdb layered_vpv = vvecdb(1, vecdb());
   vvecdb layered_vsv = vvecdb(1, vecdb());
   vvecdb layered_qkappa = vvecdb(1, vecdb());
   vvecdb layered_qshear = vvecdb(1, vecdb());
   vvecdb layered_vph = vvecdb(1, vecdb());
   vvecdb layered_vsh = vvecdb(1, vecdb());
   vvecdb layered_eta = vvecdb(1, vecdb());

   bool _isisotropic = true;
   std::vector<bool> _issolid = std::vector<bool>(1, true);
   std::vector<std::vector<int>> _vec_indices;
   // density
   // Interpolation::CubicSpline<std::vector<double>::iterator,
   //                            std::vector<double>::iterator>
   //     checkval;
   std::vector<InterpA> func_rho, func_vpv, func_vsv, func_qkappa, func_qshear,
       func_vph, func_vsh, func_eta;
   // std::vector<InterpC> func_rhoc;
   // vvecdb layered_radii = vvecdb(1, vecdb());
   // std::vector<std::vector<double>> layered_vpv, layered_vsv,
   //     layered_qkappa, layered_qshear, layered_vph, layered_vsh, layered_eta;

   // model information
   std::string modeltitle;
   int ifanis, ifdeck, numnodes, nic, noc, _numlayers;
   double tref;

   FLOAT length_norm, mass_norm, time_norm, density_norm, inertia_norm,
       velocity_norm, acceleration_norm, force_norm, stress_norm,
       gravitational_constant;

   ///////////////////////////////////////////////
   /////////////// ??????????????? ///////////////
   ///////////////////////////////////////////////

   // find layers of model
   std::vector<double> findlayers(const std::vector<double> &vec_sorted) {
      std::vector<double> vec_bounds;
      auto i1 = vec_sorted.begin();
      while (i1 != vec_sorted.end()) {
         vec_bounds.push_back(*i1);
         i1 = std::adjacent_find(++i1, vec_sorted.end());
      }
      vec_bounds.push_back(*(--i1));
      return vec_bounds;
   };

   // find indices
   std::vector<std::size_t>
   layerindices(const std::vector<double> &vec_sorted) {
      std::vector<std::size_t> vec_indices;
      auto i1 = vec_sorted.begin();
      while (i1 != vec_sorted.end()) {
         vec_indices.push_back(std::distance(vec_sorted.begin(), i1));
         i1 = std::adjacent_find(++i1, vec_sorted.end());
      }
      vec_indices.push_back(
          std::distance(vec_sorted.begin(), vec_sorted.end()) - 1);
      return vec_indices;
   }
};

template <typename FLOAT, typename INTEGRAL>
   requires std::same_as<FLOAT, double>
RadialEarthModel<FLOAT, INTEGRAL>::RadialEarthModel(
    const GPLSpec::ModelInput::PremDeck &deck)
    : RadialEarthModel(deck, EarthModels::EarthConstants<FLOAT>()) {}

// Construction from owned parsed PREM values. Existing material-layer, scaling,
// and interpolation operations are retained; parsing lives in ModelInput.
template <typename FLOAT, typename INTEGRAL>
   requires std::same_as<FLOAT, double>
template <template <typename> class ParameterModel>
RadialEarthModel<FLOAT, INTEGRAL>::RadialEarthModel(
    const GPLSpec::ModelInput::PremDeck &deck,
    const ParameterModel<FLOAT> &ModelConstants)
    : modeltitle(deck.title), ifanis(deck.anisotropy_flag),
      ifdeck(deck.deck_flag), numnodes(deck.node_count),
      nic(deck.inner_core_node_count), noc(deck.outer_core_node_count),
      tref(deck.reference_time), length_norm(ModelConstants.LengthNorm()),
      mass_norm(ModelConstants.MassNorm()), time_norm(ModelConstants.TimeNorm()),
      density_norm(ModelConstants.DensityNorm()),
      velocity_norm(ModelConstants.VelocityNorm()),
      acceleration_norm(ModelConstants.AccelerationNorm()),
      force_norm(ModelConstants.ForceNorm()),
      stress_norm(ModelConstants.StressNorm()),
      inertia_norm(ModelConstants.InertiaNorm()),
      gravitational_constant(ModelConstants.GravitationalConstant()) {
   GPLSpec::ModelInput::ValidatePremDeck(deck);
   int laynum = 0;
   int idxouter = 0;
   int idxinner = 0;

   for (const auto &node : deck.nodes) {
      const double radius = node.radius;
      const double rho = node.density;
      const double vpv = node.vpv;
      const double vsv = node.vsv;
      const double qkappa = node.qkappa;
      const double qshear = node.qshear;
      const double vph = node.vph;
      const double vsh = node.vsh;
      const double eta = node.eta;

      if (idxinner > 0 && (radius / this->length_norm ==
                           layered_radii[laynum][idxinner - 1])) {
         layered_radii.push_back({radius / this->length_norm});
         layered_rho.push_back({rho / this->density_norm});
         layered_vpv.push_back({vpv / this->velocity_norm});
         layered_vsv.push_back({vsv / this->velocity_norm});
         layered_qkappa.push_back({qkappa});
         layered_qshear.push_back({qshear});
         layered_vph.push_back({vph / this->velocity_norm});
         layered_vsh.push_back({vsh / this->velocity_norm});
         layered_eta.push_back({eta});

         if (_isisotropic && vpv != vsv) {
            _isisotropic = false;
         }
         if (vsv == 0.0 && vsh == 0) {
            _issolid.push_back(false);
         } else {
            _issolid.push_back(true);
         }
         idxinner = 0;
         ++laynum;
      } else {
         layered_radii[laynum].push_back(radius / this->length_norm);
         layered_rho[laynum].push_back(rho / this->density_norm);
         layered_vpv[laynum].push_back(vpv / this->velocity_norm);
         layered_vsv[laynum].push_back(vsv / this->velocity_norm);
         layered_qkappa[laynum].push_back(qkappa);
         layered_qshear[laynum].push_back(qshear);
         layered_vph[laynum].push_back(vph / this->velocity_norm);
         layered_vsh[laynum].push_back(vsh / this->velocity_norm);
         layered_eta[laynum].push_back(eta);

         if (idxouter == 0 && vsv == 0.0 && vsh == 0.0) {
            _issolid[laynum] = false;
         }
      }
      ++idxinner;
      ++idxouter;
   }

   _numlayers = layered_radii.size();
   vec_layers.reserve(_numlayers + 1);
   vec_layers.push_back(0.0);
   _vec_indices =
       std::vector<std::vector<int>>(_numlayers, std::vector<int>(2, 0));
   _vec_indices[0][0] = 0;
   for (int idx = 0; idx < _numlayers; ++idx) {
      vec_layers.push_back(layered_radii[idx].back());
      if (idx != 0) {
         _vec_indices[idx][0] = _vec_indices[idx - 1][1] + 1;
      }
      _vec_indices[idx][1] =
          _vec_indices[idx][0] + layered_radii[idx].size() - 1;
   }

   for (int idx = 0; idx < _numlayers; ++idx) {
      auto it1 = layered_radii[idx].begin();
      auto it2 = layered_radii[idx].end();
      auto it_rho = layered_rho[idx].begin();
      auto it_vpv = layered_vpv[idx].begin();
      auto it_vsv = layered_vsv[idx].begin();
      auto it_qkappa = layered_qkappa[idx].begin();
      auto it_qshear = layered_qshear[idx].begin();
      auto it_vph = layered_vph[idx].begin();
      auto it_vsh = layered_vsh[idx].begin();
      auto it_eta = layered_eta[idx].begin();
      func_rho.push_back(InterpA(it1, it2, it_rho));
      func_vpv.push_back(InterpA(it1, it2, it_vpv));
      func_vsv.push_back(InterpA(it1, it2, it_vsv));
      func_qkappa.push_back(InterpA(it1, it2, it_qkappa));
      func_qshear.push_back(InterpA(it1, it2, it_qshear));
      func_vph.push_back(InterpA(it1, it2, it_vph));
      func_vsh.push_back(InterpA(it1, it2, it_vsh));
      func_eta.push_back(InterpA(it1, it2, it_eta));
   }
}


}  // namespace GeneralEarthModels

#endif  // GPLSPEC_RADIAL_EARTH_MODEL_H
