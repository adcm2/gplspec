#include <gplspec/GeneralModels>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef GPLSPEC_SOURCE_DIR
#error "GPLSPEC_SOURCE_DIR must identify the fixture tree"
#endif

namespace {

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

   const TestTools::EarthModel<> earth(prem);
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
      (void)TestTools::EarthModel<>(invalid_deck);
   } catch (const std::invalid_argument &) {
      rejected_inconsistent_deck = true;
   }
   require(rejected_inconsistent_deck,
           "PREM model construction accepted a mismatched row count");

   const auto phobos =
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
}
