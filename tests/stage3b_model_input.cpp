#include <gplspec/GeneralModels>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef GPLSPEC_SOURCE_DIR
#error "GPLSPEC_SOURCE_DIR must identify the fixture tree"
#endif

namespace {

template <typename Scalar>
concept RadialModelScalarSupported = requires {
   typename GeneralEarthModels::RadialEarthModel<Scalar>;
};

static_assert(RadialModelScalarSupported<double>);
static_assert(!RadialModelScalarSupported<float>);

void require(bool condition, const std::string &message) {
   if (!condition) {
      throw std::runtime_error(message);
   }
}

bool close(double actual, double expected, double tolerance = 1e-12) {
   return std::isfinite(actual) &&
          std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected));
}

}  // namespace

int main() {
   const std::string root = GPLSPEC_SOURCE_DIR;
   const auto prem = GPLSpec::ModelInput::ReadPremDeck(root + "/modeldata/PREM500");
   require(prem.title == " prem500             ", "PREM title changed");
   require(prem.anisotropy_flag == 1 && close(prem.reference_time, 1.0) &&
               prem.deck_flag == 1,
           "PREM header fields changed");
   require(prem.node_count == 504 && prem.inner_core_node_count == 96 &&
               prem.outer_core_node_count == 273 && prem.nodes.size() == 504,
           "PREM node metadata changed");
   require(prem.additional_count_values == std::vector<int>{2},
           "PREM trailing count field was not retained");
   require(close(prem.nodes.front().radius, 0.0) &&
               close(prem.nodes.front().density, 13088.50) &&
               close(prem.nodes.front().vpv, 11262.20) &&
               close(prem.nodes.front().vsv, 3667.80),
           "PREM center row differs from the independent fixture values");
   require(close(prem.nodes[1].radius, 12858.0) &&
               close(prem.nodes[1].density, 13088.46),
           "PREM second row differs from the independent fixture values");
   bool found_one_sided_boundary = false;
   for (std::size_t i = 1; i < prem.nodes.size(); ++i) {
      if (prem.nodes[i - 1].radius == prem.nodes[i].radius) {
         found_one_sided_boundary = true;
         break;
      }
   }
   require(found_one_sided_boundary,
           "PREM reader lost duplicated-radius material interfaces");

   const auto prem200 = GPLSpec::ModelInput::ReadPremDeck(root + "/modeldata/prem.200");
   require(prem200.node_count == 220 && prem200.nodes.size() == 220 &&
               prem200.additional_header_values == std::vector<double>{2.0},
           "prem.200 extra header field or row count changed");
   require(close(prem200.nodes.front().density, 13088.50) &&
               close(prem200.nodes[1].radius, 31321.0),
           "prem.200 known raw values changed");

   using RadialModel = GeneralEarthModels::RadialEarthModel<>;
   static_assert(!std::is_copy_constructible_v<RadialModel>);
   static_assert(!std::is_copy_assignable_v<RadialModel>);
   static_assert(std::is_move_constructible_v<RadialModel>);
   const RadialModel earth(prem);
   const EarthModels::EarthConstants<double> units;
   int expected_layers = 1;
   std::size_t first_boundary_index = prem.nodes.size();
   for (std::size_t i = 1; i < prem.nodes.size(); ++i) {
      if (prem.nodes[i - 1].radius == prem.nodes[i].radius) {
         ++expected_layers;
         if (first_boundary_index == prem.nodes.size()) {
            first_boundary_index = i;
         }
      }
   }
   require(earth.NumberOfLayers() == expected_layers,
           "PREM model construction changed one-sided layer count");
   require(close(earth.OuterRadius(),
                 prem.nodes.back().radius / units.LengthNorm()),
           "PREM model radius normalization changed");
   require(close(earth.LayerRadii(0).front(),
                 prem.nodes.front().radius / units.LengthNorm()) &&
               close(earth.LayerRadii(0).back(),
                     prem.nodes[first_boundary_index - 1].radius /
                         units.LengthNorm()),
           "PREM first-layer boundaries changed");
   require(close(earth.Density(0)(earth.LayerRadii(0).front()),
                 prem.nodes.front().density / units.DensityNorm()),
           "PREM density scale changed during model construction");

   auto invalid_deck = prem;
   invalid_deck.node_count += 1;
   bool rejected_inconsistent_deck = false;
   try {
      (void)GeneralEarthModels::RadialEarthModel<>(invalid_deck);
   } catch (const std::invalid_argument &) {
      rejected_inconsistent_deck = true;
   }
   require(rejected_inconsistent_deck,
           "PREM model construction accepted a mismatched row count");

   auto short_lived_deck =
       GPLSpec::ModelInput::ReadPremDeck(root + "/modeldata/PREM500");
   const GeneralEarthModels::RadialEarthModel<> owned_earth(short_lived_deck);
   const auto owned_earth_layers = owned_earth.NumberOfLayers();
   const auto owned_earth_density =
       owned_earth.Density(0)(owned_earth.LayerRadii(0).front());
   short_lived_deck.nodes.clear();
   require(owned_earth.NumberOfLayers() == owned_earth_layers &&
               close(owned_earth.Density(0)(owned_earth.LayerRadii(0).front()),
                     owned_earth_density),
           "radial model retained a borrow into the parsed PREM deck");
   const auto moved_earth = [&root]() {
      auto deck = GPLSpec::ModelInput::ReadPremDeck(root + "/modeldata/PREM500");
      RadialModel source(deck);
      return RadialModel(std::move(source));
   }();
   require(moved_earth.NumberOfLayers() == expected_layers &&
               close(moved_earth.Density(0)(moved_earth.LayerRadii(0).front()),
                     owned_earth_density),
           "moved radial model lost its owned interpolation storage");

   auto phobos =
       GPLSpec::ModelInput::ReadPhobosHarmonics(root + "/modeldata/TABLEA1.DAT");
   require(phobos.rows.size() == 1081 && phobos.InferredMaximumDegree() == 45,
           "Phobos TABLEA1 triangular coefficient extent changed");
   require(phobos.rows.front().degree == 0 && phobos.rows.front().order == 0 &&
               close(phobos.rows.front().cosine, 38970.54) &&
               close(phobos.rows.front().cosine_sigma, 1.685397) &&
               close(phobos.rows.front().sine, 0.0) &&
               close(phobos.rows.front().sine_sigma, 0.0),
           "Phobos degree-zero row differs from independent fixture values");
   require(phobos.rows[2].degree == 1 && phobos.rows[2].order == 1 &&
               close(phobos.rows[2].cosine, -1.539084e-3) &&
               close(phobos.rows[2].sine, -4.209417e-5),
           "Phobos degree-one complex row differs from independent fixture values");
   require(phobos.rows.back().degree == 45 && phobos.rows.back().order == 45,
           "Phobos maximum-degree row changed");
   bool rejected_empty_phobos = false;
   try {
      (void)GPLSpec::ModelInput::PhobosHarmonics{}.InferredMaximumDegree();
   } catch (const std::logic_error &) {
      rejected_empty_phobos = true;
   }
   require(rejected_empty_phobos,
           "Phobos maximum degree accepted an empty row list");

   const auto malformed_path =
       std::filesystem::temp_directory_path() / "gplspec-stage3b-input-contract.tmp";
   const auto write_input_fixture = [&malformed_path](const std::string &contents) {
      std::ofstream file(malformed_path);
      file << contents;
      file.close();
      require(static_cast<bool>(file), "failed to create temporary parser fixture");
   };
   const std::string tiny_prem_rows =
       "2 0 0\n0 1 2 3 4 5 6 7 8\n1 1 2 3 4 5 6 7 8\n";
   const auto prem_error_contains = [&malformed_path](const std::string &expected) {
      try {
         (void)GPLSpec::ModelInput::ReadPremDeck(malformed_path.string());
      } catch (const std::runtime_error &error) {
         return std::string(error.what()).find(expected) != std::string::npos;
      }
      return false;
   };
   write_input_fixture(std::string("tiny PREM\n1 1.0 1 3.5\n2 0 0 4\n") +
                       "0 1 2 3 4 5 6 7 8\n1 1 2 3 4 5 6 7 8\n");
   const auto valid_tiny_prem =
       GPLSpec::ModelInput::ReadPremDeck(malformed_path.string());
   require(valid_tiny_prem.additional_header_values == std::vector<double>{3.5} &&
               valid_tiny_prem.additional_count_values == std::vector<int>{4} &&
               valid_tiny_prem.nodes.size() == 2,
           "PREM reader rejected a complete valid tiny deck");

   write_input_fixture("tiny PREM\n1 1.0 1 trailing\n" + tiny_prem_rows);
   require(prem_error_contains("malformed optional numeric field at line 2"),
           "PREM reader accepted a nonnumeric optional header field");
   write_input_fixture("tiny PREM\n1 1.0 1 2e\n" + tiny_prem_rows);
   require(prem_error_contains("malformed optional numeric field at line 2"),
           "PREM reader accepted an incomplete optional numeric field");
   write_input_fixture(std::string("tiny PREM\n1 1.0 1\n2 0 0 trailing\n") +
                       "0 1 2 3 4 5 6 7 8\n1 1 2 3 4 5 6 7 8\n");
   require(prem_error_contains("malformed optional numeric field at line 3"),
           "PREM reader accepted a nonnumeric optional count field");
   write_input_fixture("0 0 1 0 0 0\n1 1 1 0 0 0\n1 0 1 0 0 0\n");
   bool rejected_phobos_order = false;
   try {
      (void)GPLSpec::ModelInput::ReadPhobosHarmonics(malformed_path.string());
   } catch (const std::runtime_error &error) {
      rejected_phobos_order =
          std::string(error.what()).find("triangular degree/order order") !=
          std::string::npos;
   }
   require(rejected_phobos_order,
           "Phobos reader accepted rows outside triangular degree/order order");
   std::filesystem::remove(malformed_path);

   GPLSpec::ModelInput::PhobosHarmonics nonfinite_phobos;
   nonfinite_phobos.rows.push_back(
       {0, 0, std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0});
   bool rejected_nonfinite_phobos = false;
   try {
      (void)nonfinite_phobos.InferredMaximumDegree();
   } catch (const std::invalid_argument &) {
      rejected_nonfinite_phobos = true;
   }
   require(rejected_nonfinite_phobos,
           "Phobos triangular contract accepted a non-finite coefficient");

   const double length_norm = 10000.0;
   const double time_norm = 3600.0;
   const double mass_norm = std::pow(length_norm, 3.0) * 1342.0;
   const double density_scale = 1860.0 / 1342.0;
   const double expected_phobos_radius =
       phobos.rows.front().cosine / std::sqrt(4.0 * std::acos(-1.0));
   GeneralEarthModels::Density3D owned_phobos(
       density_scale, phobos, 2, 2, length_norm, time_norm, mass_norm, 0.5,
       1.5);
   const double owned_phobos_density = owned_phobos.Density_Point(0, 0, 0);
   require(std::isfinite(owned_phobos_density) && owned_phobos_density > 0.0,
           "Phobos model construction returned invalid owned density");
   require(close(owned_phobos.Node_Information().PlanetRadius(),
                 expected_phobos_radius / length_norm),
           "Phobos model geometry changed its A00 radius convention");
   phobos.rows.clear();
   require(close(owned_phobos.Node_Information().PlanetRadius(),
                 expected_phobos_radius / length_norm) &&
               close(owned_phobos.Density_Point(0, 0, 0),
                     owned_phobos_density),
           "Phobos model retained a borrow into parsed harmonics");
}
