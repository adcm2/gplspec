#ifndef GPLSPEC_MODEL_INPUT_H
#define GPLSPEC_MODEL_INPUT_H

#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace GPLSpec::ModelInput {

struct PremNode {
   double radius = 0.0;
   double density = 0.0;
   double vpv = 0.0;
   double vsv = 0.0;
   double qkappa = 0.0;
   double qshear = 0.0;
   double vph = 0.0;
   double vsh = 0.0;
   double eta = 0.0;
};

struct PremDeck {
   std::string title;
   int anisotropy_flag = 0;
   double reference_time = 0.0;
   int deck_flag = 0;
   std::vector<double> additional_header_values;
   int node_count = 0;
   int inner_core_node_count = 0;
   int outer_core_node_count = 0;
   std::vector<int> additional_count_values;
   std::vector<PremNode> nodes;
};

struct PhobosHarmonicRow {
   int degree = 0;
   int order = 0;
   double cosine = 0.0;
   double cosine_sigma = 0.0;
   double sine = 0.0;
   double sine_sigma = 0.0;
};

struct PhobosHarmonics {
   std::vector<PhobosHarmonicRow> rows;

   // TABLEA1's historical constructor infers the maximum degree from its row
   // count using this triangular-layout relationship.
   std::size_t InferredMaximumDegree() const {
      if (rows.empty()) {
         throw std::logic_error("cannot infer Phobos degree from no rows");
      }
      const auto maximum_degree = static_cast<std::size_t>(
          (std::sqrt(1.0 + 8.0 * rows.size()) - 3.0) / 2.0);
      if ((maximum_degree + 1) * (maximum_degree + 2) / 2 != rows.size()) {
         throw std::invalid_argument(
             "Phobos harmonic row count is not triangular");
      }

      std::size_t row_index = 0;
      for (std::size_t degree = 0; degree <= maximum_degree; ++degree) {
         for (std::size_t order = 0; order <= degree; ++order, ++row_index) {
            const auto &row = rows[row_index];
            if (row.degree != static_cast<int>(degree) ||
                row.order != static_cast<int>(order)) {
               throw std::invalid_argument(
                   "Phobos harmonic rows must be in triangular degree/order order");
            }
            if (!std::isfinite(row.cosine) ||
                !std::isfinite(row.cosine_sigma) ||
                !std::isfinite(row.sine) ||
                !std::isfinite(row.sine_sigma)) {
               throw std::invalid_argument(
                   "Phobos harmonic row contains a non-finite value");
            }
         }
      }
      return maximum_degree;
   }
};

inline void ValidatePremDeck(const PremDeck &deck) {
   if (deck.node_count < 2 ||
       deck.nodes.size() != static_cast<std::size_t>(deck.node_count)) {
      throw std::invalid_argument(
          "PREM deck node count must match at least two stored nodes");
   }
   if (deck.inner_core_node_count < 0 || deck.outer_core_node_count < 0 ||
       !std::isfinite(deck.reference_time)) {
      throw std::invalid_argument("PREM deck metadata is outside its domain");
   }
   std::size_t nodes_in_layer = 1;
   const auto finite_node = [](const PremNode &node) {
      return std::isfinite(node.radius) && std::isfinite(node.density) &&
             std::isfinite(node.vpv) && std::isfinite(node.vsv) &&
             std::isfinite(node.qkappa) && std::isfinite(node.qshear) &&
             std::isfinite(node.vph) && std::isfinite(node.vsh) &&
             std::isfinite(node.eta);
   };
   if (!finite_node(deck.nodes.front())) {
      throw std::invalid_argument("PREM deck contains a non-finite node");
   }
   for (std::size_t i = 1; i < deck.nodes.size(); ++i) {
      if (!finite_node(deck.nodes[i]) ||
          deck.nodes[i].radius < deck.nodes[i - 1].radius) {
         throw std::invalid_argument(
             "PREM deck nodes must be finite and radius ordered");
      }
      if (deck.nodes[i].radius == deck.nodes[i - 1].radius) {
         if (nodes_in_layer < 2) {
            throw std::invalid_argument(
                "PREM interpolation layers require at least two nodes");
         }
         nodes_in_layer = 1;
      } else {
         ++nodes_in_layer;
      }
   }
   if (nodes_in_layer < 2) {
      throw std::invalid_argument(
          "PREM interpolation layers require at least two nodes");
   }
}

namespace detail {

template <class... Values>
void ReadExactFields(std::istream &input, const std::string &path,
                     std::size_t line, Values &...values) {
   if (!(input >> ... >> values)) {
      throw std::runtime_error(path + ": malformed input at line " +
                               std::to_string(line));
   }
}

inline bool IsBlank(const std::string &line) {
   return line.find_first_not_of(" \t\r") == std::string::npos;
}

template <class Value>
void ReadOptionalNumericFields(std::istream &input, std::vector<Value> &values,
                               const std::string &path, std::size_t line) {
   std::string token;
   while (input >> token) {
      std::istringstream field(token);
      Value value{};
      if (!(field >> value)) {
         throw std::runtime_error(path + ": malformed optional numeric field at line " +
                                  std::to_string(line));
      }
      char extra = '\0';
      if (field >> extra) {
         throw std::runtime_error(path + ": malformed optional numeric field at line " +
                                  std::to_string(line));
      }
      if constexpr (std::is_floating_point_v<Value>) {
         if (!std::isfinite(value)) {
            throw std::runtime_error(path + ": non-finite optional numeric field at line " +
                                     std::to_string(line));
         }
      }
      values.push_back(value);
   }
}

inline void RequireNoTrailingFields(std::istream &input,
                                    const std::string &path,
                                    std::size_t line) {
   std::string extra;
   if (input >> extra) {
      throw std::runtime_error(path + ": unexpected extra field at line " +
                               std::to_string(line));
   }
}

inline std::string NextNonblankLine(std::ifstream &file,
                                    const std::string &path,
                                    std::size_t &line_number) {
   std::string line;
   while (std::getline(file, line)) {
      ++line_number;
      if (!IsBlank(line)) {
         return line;
      }
   }
   throw std::runtime_error(path + ": unexpected end of file");
}

}  // namespace detail

inline PremDeck ReadPremDeck(const std::string &path) {
   std::ifstream file(path);
   if (!file) {
      throw std::runtime_error("cannot open PREM model: " + path);
   }

   PremDeck deck{};
   std::size_t line_number = 0;
   deck.title = detail::NextNonblankLine(file, path, line_number);

   std::string line = detail::NextNonblankLine(file, path, line_number);
   std::istringstream header(line);
   detail::ReadExactFields(header, path, line_number, deck.anisotropy_flag,
                           deck.reference_time, deck.deck_flag);
   detail::ReadOptionalNumericFields(header, deck.additional_header_values, path,
                                    line_number);

   line = detail::NextNonblankLine(file, path, line_number);
   std::istringstream counts(line);
   detail::ReadExactFields(counts, path, line_number, deck.node_count,
                           deck.inner_core_node_count,
                           deck.outer_core_node_count);
   detail::ReadOptionalNumericFields(counts, deck.additional_count_values, path,
                                    line_number);
   if (deck.node_count < 0) {
      throw std::runtime_error(path + ": negative PREM node count");
   }

   deck.nodes.reserve(static_cast<std::size_t>(deck.node_count));
   for (int index = 0; index < deck.node_count; ++index) {
      line = detail::NextNonblankLine(file, path, line_number);
      std::istringstream node_line(line);
      PremNode node{};
      detail::ReadExactFields(node_line, path, line_number, node.radius,
                              node.density, node.vpv, node.vsv, node.qkappa,
                              node.qshear, node.vph, node.vsh, node.eta);
      detail::RequireNoTrailingFields(node_line, path, line_number);
      deck.nodes.push_back(node);
   }
   ValidatePremDeck(deck);
   return deck;
}

inline PhobosHarmonics ReadPhobosHarmonics(const std::string &path) {
   std::ifstream file(path);
   if (!file) {
      throw std::runtime_error("cannot open Phobos harmonic model: " + path);
   }

   PhobosHarmonics harmonics;
   std::string line;
   std::size_t line_number = 0;
   while (std::getline(file, line)) {
      ++line_number;
      if (detail::IsBlank(line)) {
         continue;
      }
      std::istringstream row_line(line);
      PhobosHarmonicRow row{};
      detail::ReadExactFields(row_line, path, line_number, row.degree,
                              row.order, row.cosine, row.cosine_sigma,
                              row.sine, row.sine_sigma);
      detail::RequireNoTrailingFields(row_line, path, line_number);
      harmonics.rows.push_back(row);
   }
   if (harmonics.rows.empty()) {
      throw std::runtime_error(path + ": empty Phobos harmonic model");
   }
   try {
      (void)harmonics.InferredMaximumDegree();
   } catch (const std::invalid_argument &error) {
      throw std::runtime_error(path + ": " + error.what());
   }
   return harmonics;
}

}  // namespace GPLSpec::ModelInput

#endif  // GPLSPEC_MODEL_INPUT_H
